# taxonium_phylogeny_mutations_to_dataframe.R
# Script para convertir una filogenia en json a una tabla con cada nodo y muestra y sus mutaciones

library(jsonlite)
library(stringr)
library(dplyr)
library(progress)

args <- commandArgs(trailingOnly = TRUE) # 1: input .json, 2: output .tsv
tree_json <- fromJSON(args[1])
mutations_codes <- tree_json[[2]][[1]] # Save the mutation info in a dataframe

df_mutations_per_node <- as.data.frame(matrix(ncol = 2, nrow = 1))
colnames(df_mutations_per_node) <- c("Sample_or_node", "Mutations")
node_names <- tree_json$name
node_names <- node_names[-1]
node_ids <- tree_json$node_id
node_ids <- node_ids[-1]
vector_final <- c()

total_iterations <- length(node_ids)
pb <- progress_bar$new(
  format = "  Progress: [:bar] :percent Completed | Current: :current/:total | Elapsed: :elapsed | ETA: :eta",
  total = total_iterations,
  clear = FALSE
)

for(i in 1:length(node_ids)){
  pb$tick()
  if(node_names[i] == ""){
    vector_final <- c(vector_final, node_ids[i])
  } else{
    vector_final <- c(vector_final, node_names[i])
  }
}


count <- 1
pb <- progress_bar$new(
  format ="  Progress: [:bar] :percent Completed | Current: :current/:total | Elapsed: :elapsed | ETA: :eta",
  total = length(node_names),
  clear = FALSE
)

for(name in vector_final){
  pb$tick()
  #name <- vector_final[i]
  
  if(!grepl("node", vector_final[i])){
    codes <- tree_json[tree_json$name == name,]$mutations[[2]]} #else{
    #codes <- tree_json[tree_json$name == as.numeric(name),]
  #}
    
  vector_mutations <- c()
  
  for(code in codes){
    mutation <- paste(mutations_codes[mutations_codes$mutation_id == code, c(2,3,4)], collapse = "") # Extract the mutation using the code
    vector_mutations <- c(vector_mutations, mutation)
  }
  
  df_mutations_per_node[count,1] <- name
  df_mutations_per_node[count,2] <- paste(vector_mutations, collapse = ",")
  count <- count+1
}

write.table(df_mutations_per_node, args[2], sep = "\t", quote = FALSE, row.names = FALSE)