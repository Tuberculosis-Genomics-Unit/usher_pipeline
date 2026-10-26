# UShER phylogeny pipeline

Builds a phylogenetic tree with [UShER](https://github.com/yatisht/usher) from the multifasta and SNP table produced by [ThePipeline](https://github.com/Tuberculosis-Genomics-Unit/ThePipeline_v3). Optionally, it also writes a table with the mutations on every branch of the tree.

Starting from `.annoF` files (no rescue mode) is not supported (deprecated, but can be added in the future).

<br>

## Repository contents

```
usher_pipeline/
├── pipeline.py            # the pipeline: runs every step in order
├── config.toml            # settings: input files, output paths, threads
├── setup.sh               # installs the software environment
├── data/
│   ├── MTB_ancestor_reference.fasta   # reference genome (MTB ancestor)
│   └── R00000039_repregions.bed       # regions of the genome to mask
├── scripts/               # scripts called by the pipeline
└── envs/                  # exact software versions used by setup.sh
```

<br>

## What you need

- A Linux server.
- `conda`. To check, run `conda --version`. If you get `command not found`, see [Installing conda](#installing-conda).
- About 2 GB of free space in your home directory, for the software environment.
- Your input files: a **multifasta** and a **SNP table**, both from [`ThePipeline consensus`](https://github.com/Tuberculosis-Genomics-Unit/ThePipeline_v3) (see [Step 3](#step-3--check-your-input-files)).

## Step 1 – Get the code

```shell
# SSH (requires an SSH key registered with GitHub)
git clone git@github.com:Tuberculosis-Genomics-Unit/usher_pipeline.git
# or HTTPS
git clone https://github.com/Tuberculosis-Genomics-Unit/usher_pipeline.git

cd usher_pipeline
```

All the following commands are run from inside this `usher_pipeline` folder.

## Step 2 – Install the software environment

Only once per user and server:

```shell
bash setup.sh
```

This shouldn't take long. It installs UShER, R, taxonium, jq and everything else the pipeline needs into your own environment at `~/.conda/envs/usher-env`. The installed versions are exactly the same on every server. Each user gets their own copy, even when `conda` is a shared installation on the server.

Running `bash setup.sh` again deletes the environment and installs it again from scratch.

## Step 3 – Check your input files

| File | Where it comes from |
|---|---|
| Multifasta | [`ThePipeline consensus`](https://github.com/Tuberculosis-Genomics-Unit/ThePipeline_v3) |
| SNP table | [`ThePipeline consensus`](https://github.com/Tuberculosis-Genomics-Unit/ThePipeline_v3) |
| Reference | Already in `data/` |
| Mask (`.bed`) | Already in `data/` |

## Step 4 – Set parameters

### In `config.toml`

Open `config.toml` in any text editor. Normally you only change three lines: your two input files and the output folder.

```toml
[input]
multifasta = "data/aln.fasta"                       # <- your multifasta
snp_table = "data/snp_table.txt"                    # <- your SNP table

[data]
ref_fasta = "data/MTB_ancestor_reference.fasta"     # already provided, don't change
bed_mask = "data/R00000039_repregions.bed"          # already provided, don't change

[output]
dir = "output"                                      # <- folder where all results go
keep_intermediates = false                          # true = keep the per-sample files

[mapping]
enabled = false                                     # true = also make the mutation table

[resources]
threads = 8                                         # CPU cores to use
```

You only choose the output **folder**. The file names inside it are always the same (see [Step 6](#step-6--results)). The folder is created if it doesn't exist. If it already holds results from an earlier run, they are overwritten.

Rules for paths:

- You can use **absolute paths** (`/data/myproject/aln.fasta`) or **relative paths** (`data/aln.fasta`). Relative paths are read from the `usher_pipeline` folder, wherever you run the command from.
- The output folder path **must not contain dots** (for example `/data/run.2024/` will fail), and must not start with `./`. The scripts in `scripts/` read the sample name by cutting the path at the first dot.

### On the command line (without editing `config.toml`)

Some settings can be given on the command line. They take priority over `config.toml`:

```shell
python pipeline.py \
    --multifasta aln.fasta \
    --snp-table snp_table.txt \
    -o results \
    --threads 4 \
    --mapping
```

| Option | What it does |
|---|---|
| `--multifasta FILE` | Multifasta to use |
| `--snp-table FILE` | SNP table to use |
| `-o FOLDER`, `--output FOLDER` | Folder where all results go |
| `--threads N` | Number of CPU cores |
| `--mapping` | Also make the mutation table |
| `--keep-intermediates` | Keep the per-sample files (`vcf/`, `diff/`) instead of deleting them |
| `another_config.toml` | Use a different config file instead of `config.toml` |

Run `python pipeline.py --help` to see all options.

Give each dataset its own output folder (`-o`). Running twice with the same folder overwrites the first run's results.

## Step 5 – Run the pipeline

```shell
conda activate ~/.conda/envs/usher-env
python pipeline.py
```

The pipeline prints its progress:

```
[1/4] Splitting multifasta (250 samples)
[2/4] Building per-sample .diff files (8 parallel jobs)
[3/4] Merging .diff files
[4/4] Building tree with usher-sampled
Done. Results in /home/<you>/usher_pipeline/output
```

With the mutation table turned on (`--mapping`), there is a fifth step, `[5/5] Mapping mutations to tree nodes`.

### Long runs: keep it running after you log out

```shell
nohup python pipeline.py > pipeline.log 2>&1 &
```

## Step 6 – Results

Everything is inside the output folder (`output/` by default), always with these names:

```
output/
├── usher_tree.nwk
├── usher_tree.pb
├── usher_mapped_mutations.tsv
├── final.diff
├── usher_output/
├── vcf/          (only with keep_intermediates = true)
└── diff/         (only with keep_intermediates = true)
```

| File | What it is |
|---|---|
| `usher_tree.nwk` | **The final tree** in Newick format. Open it with FigTree, iTOL, R (`ape`), etc. |
| `usher_tree.pb` | The same tree in UShER's format, with the mutations on every branch. Used by UShER tools and taxonium. |
| `usher_mapped_mutations.tsv` | **The mutation table** (only with `--mapping` or `[mapping] enabled = true`). See below. |
| `final.diff` | Variants of all samples, in the format UShER reads. |
| `usher_output/` | Extra files written by UShER (for example `placement_stats.tsv`, with how confidently each sample was placed). |
| `vcf/`, `diff/` | One file per sample from the intermediate steps, plus one log file per sample in `diff/`. Deleted at the end unless `keep_intermediates = true` (or `--keep-intermediates`). |

> [!WARNING]
> Just in case you need to use it: `usher_tree.pb` contains an extra sample called `ref`. It is a "template", used by UShER to build the tree, not one of your samples. The pipeline removes it from `usher_tree.nwk` and from `usher_mapped_mutations.tsv` to keep everything clean, but not from `usher_tree.pb`.

### The mutation table

Tab-separated, one row per sample and per internal node of the tree:

```
Sample_or_node	Mutations
S3	A1315C,C1926A,T4352A
node_8	T630A,T7381A
S1	T1164A,C1640T
```

- Each row lists the mutations on the **branch leading to that node**, not all the mutations accumulated from the root.
- Mutations are written as `<reference base><position><new base>`, for example `A1315C`.
- The names of internal nodes (`node_8`, …) are the same as in `usher_tree.nwk`, so the table can be joined to the tree.

## Troubleshooting

| Message | What to do |
|---|---|
| `conda: command not found` | conda is not installed or not loaded. See [Installing conda](#installing-conda). |
| `CondaError: Run 'conda init' before 'conda activate'` | Run `conda init bash`, then close and reopen the terminal. |
| `FileNotFoundError: ... aln.fasta` (or another input) | A path in `config.toml` is wrong. Check it with `ls <path>`. Remember that relative paths are read from the `usher_pipeline` folder. |
| `FileNotFoundError: 'faToVcf'`, `'usher-sampled'` or `'Rscript'` | The environment is not active. Run `conda activate ~/.conda/envs/usher-env`. |
| `Couldn't open ..._aligned.fasta` | A sample name or the output path contains a dot, or the path starts with `./`. Rename it; see the path rules in [Step 4](#step-4--set-parameters). The failed run may also leave stray `*_aligned.fasta` files in the folder above the dotted one; delete them. |
| `IndexError: list index out of range` | The multifasta has sequences split over several lines. Each sequence must be on a single line. |

### Installing conda

If the server has no `conda`, install Miniforge in your home directory (no admin rights needed):

```shell
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

Accept the defaults, answer `yes` when it asks to initialise conda, then close and reopen the terminal.
