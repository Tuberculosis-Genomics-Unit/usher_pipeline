import argparse
import shutil
import subprocess
import sys
import tomllib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from Bio import Phylo

REPO = Path(__file__).resolve().parent
SCRIPTS = REPO / "scripts"


def resolve(path: str) -> Path:
    return (REPO / path).resolve()


def run(cmd: list, **kwargs) -> None:
    subprocess.run([str(c) for c in cmd], check=True, **kwargs)


def read_samples(multifasta: Path) -> list[str]:
    with multifasta.open() as f:
        return [line[1:].split()[0] for line in f if line.startswith(">")]


def sample_to_diff(
    sample: str,
    ref_fasta: Path,
    snp_table: Path,
    bed_mask: Path,
    split_dir: Path,
    vcf_dir: Path,
    diff_dir: Path,
) -> None:
    aligned = split_dir / f"{sample}_aligned.fasta"
    vcf = vcf_dir / f"{sample}.vcf"
    # The scripts print every sequence/diff to stdout; silence it.
    run(
        [sys.executable, SCRIPTS / "align_ref_and_rescued_fasta.py", ref_fasta, split_dir / f"{sample}.fasta", snp_table],
        stdout=subprocess.DEVNULL,
    )
    run(["faToVcf", aligned, vcf])
    aligned.unlink()
    run(
        [sys.executable, SCRIPTS / "vcf_to_diff.py", "-v", vcf, "-d", diff_dir, "-smf", bed_mask],
        stdout=subprocess.DEVNULL,
    )


def prune_ref(nwk: Path) -> None:
    tree = Phylo.read(nwk, "newick")
    tree.prune("ref")
    # Keep UShER's integer branch lengths (mutation counts); Biopython defaults to %1.5f.
    Phylo.write(tree, nwk, "newick", format_branch_length="%g")


def map_mutations(final_pb: Path, mapping_table: Path, work_dir: Path) -> None:
    jsonl_gz = work_dir / "taxonium_tree.jsonl.gz"
    jsonl = work_dir / "taxonium_tree.jsonl"
    json_file = work_dir / "taxonium_tree.json"
    run(["usher_to_taxonium", "-i", final_pb, "-o", jsonl_gz])
    run(["gunzip", "-f", jsonl_gz])
    with json_file.open("w") as dst:
        run(["jq", "-s", ".", jsonl], stdout=dst)
    run(["Rscript", SCRIPTS / "taxonium_phylogeny_mutations_to_dataframe.R", json_file, mapping_table])
    # The .pb still contains ref, so drop its row here.
    lines = mapping_table.read_text().splitlines(keepends=True)
    mapping_table.write_text("".join(l for l in lines if not l.startswith("ref\t")))
    jsonl.unlink()
    json_file.unlink()


def main(args: argparse.Namespace) -> None:
    with args.config.open("rb") as f:
        config = tomllib.load(f)
    inp, data, out = config["input"], config["data"], config["output"]

    # Command-line paths resolve from the cwd; config paths from the repo.
    multifasta = args.multifasta.resolve() if args.multifasta else resolve(inp["multifasta"])
    snp_table = args.snp_table.resolve() if args.snp_table else resolve(inp["snp_table"])
    ref_fasta = resolve(data["ref_fasta"])
    bed_mask = resolve(data["bed_mask"])
    output_dir = args.output.resolve() if args.output else resolve(out["dir"])
    keep_intermediates = args.keep_intermediates or out["keep_intermediates"]
    threads = args.threads or config["resources"]["threads"]
    mapping = args.mapping or config["mapping"]["enabled"]
    n_steps = 5 if mapping else 4

    final_diff = output_dir / "final.diff"
    final_pb = output_dir / "usher_tree.pb"
    final_nwk = output_dir / "usher_tree.nwk"
    mapping_table = output_dir / "usher_mapped_mutations.tsv"
    usher_dir = output_dir / "usher_output"
    split_dir = output_dir / "split"
    vcf_dir = output_dir / "vcf"
    diff_dir = output_dir / "diff"
    for d in (vcf_dir, diff_dir, usher_dir):
        d.mkdir(parents=True, exist_ok=True)

    samples = read_samples(multifasta)
    print(f"[1/{n_steps}] Splitting multifasta ({len(samples)} samples)", flush=True)
    run(["bash", SCRIPTS / "split_multifasta.sh", multifasta, split_dir])

    print(f"[2/{n_steps}] Building per-sample .diff files ({threads} parallel jobs)", flush=True)
    with ThreadPoolExecutor(max_workers=threads) as pool:
        list(pool.map(
            lambda s: sample_to_diff(s, ref_fasta, snp_table, bed_mask, split_dir, vcf_dir, diff_dir),
            samples,
        ))
    shutil.rmtree(split_dir)

    print(f"[3/{n_steps}] Merging .diff files", flush=True)
    # Concatenate only this run's samples, so stale .diff files in diff_dir are ignored.
    with final_diff.open("wb") as dst:
        for name in sorted(f"{s}.diff" for s in samples):
            with (diff_dir / name).open("rb") as f:
                shutil.copyfileobj(f, dst)
    if not keep_intermediates:
        shutil.rmtree(vcf_dir)
        shutil.rmtree(diff_dir)

    print(f"[4/{n_steps}] Building tree with usher-sampled", flush=True)
    (usher_dir / "base_tree.nwk").write_text("();\n")
    (usher_dir / "ref.diff").write_text(">ref\n")
    run(
        ["usher-sampled", "-t", "base_tree.nwk", "--diff", "ref.diff", "--ref", ref_fasta, "-o", "ref.pb"],
        cwd=usher_dir,
    )
    run(
        ["usher-sampled", "-i", "ref.pb", "--diff", final_diff, "--ref", ref_fasta,
         "-o", final_pb, "-u", "-d", ".", "-T", threads],
        cwd=usher_dir,
    )
    shutil.move(usher_dir / "uncondensed-final-tree.nh", final_nwk)
    prune_ref(final_nwk)
    for name in ("base_tree.nwk", "ref.diff", "ref.pb"):
        (usher_dir / name).unlink()

    if mapping:
        print(f"[5/{n_steps}] Mapping mutations to tree nodes", flush=True)
        map_mutations(final_pb, mapping_table, output_dir)

    print(f"Done. Results in {output_dir}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Values given here override config.toml.")
    parser.add_argument("config", nargs="?", type=Path, default=REPO / "config.toml")
    parser.add_argument("--multifasta", type=Path)
    parser.add_argument("--snp-table", type=Path)
    parser.add_argument("-o", "--output", type=Path)
    parser.add_argument("--threads", "-t", type=int)
    parser.add_argument("--mapping", action="store_true")
    parser.add_argument("--keep-intermediates", action="store_true")
    main(parser.parse_args())
