merge_to_phy <- function(all_merges, all_leaves) {
  all_edge_name = unique(c(names(all_merges), unlist(all_merges)))
  all_edge_len = rep(1, length(all_edge_name))
  names(all_edge_len) = all_edge_name

  newick_str = map(all_leaves, function(leaf) {
    paste0(leaf, ":", all_edge_len[leaf])
  })
  names(newick_str) = all_leaves

  merge_newick <- function(newick_str1, newick_str2, node_label, edge_len) {
    paste0("(", newick_str1, ",", newick_str2, ")", node_label, ":", edge_len)
  }
  for (merge_node in names(all_merges)) {
    if (merge_node == "root") {
      newick_str[[merge_node]] = paste0(newick_str[[all_merges[[merge_node]][1]]], ";")
    } else {
      assertthat::assert_that(all_merges[[merge_node]][1] %in% names(newick_str))
      assertthat::assert_that(all_merges[[merge_node]][2] %in% names(newick_str))
      newick_str[[merge_node]] =
        merge_newick(newick_str[[all_merges[[merge_node]][1]]],
                     newick_str[[all_merges[[merge_node]][2]]],
                     merge_node,
                     all_edge_len[merge_node])
    }
  }
  tr = read.tree(text = newick_str[["root"]])
  tr
}

library(ape)
library(qfm)
library(tidyverse)
library(igraph)
library(ggraph)
setwd('/home/ssrikan2/data-kreza1/smriti/qfm2')
devtools::load_all()
global_step_size = 0.01
global_target_time = 9.0

cell_state_meta = read_csv("./metadata/simple_test 1.csv")
cell_state_col = cell_state_meta$color
names(cell_state_col) = cell_state_meta$name

cell_state_param_list = make_cell_state_param(cell_state_meta)
all_state_bias = map(cell_state_param_list, "downstream_state_bias")
all_state_merges = map(all_state_bias, names)
all_node_merges = all_state_merges[!map_lgl(all_state_merges, is.null)]
all_leaves = names(all_state_merges)[map_lgl(all_state_merges, is.null)]
all_merges = c(rev(all_node_merges), list("root" = "TypeS"))

phy = merge_to_phy(all_merges, all_leaves)
qfm_spec = make_qfm_spec(phy, cell_state_param_list)

count_graph = initialize_count_graph(qfm_spec)
tt = 0.
while (length(count_graph$active_nodes) > 0) {
  tt = tt + global_step_size
  tt = round(tt, digits = 3)
  print(paste0('time: ', format(tt, digits = 3)))
  print(paste0('number of active nodes: ', length(count_graph$active_nodes)))
  for (node_x in count_graph$active_nodes) {
    # message(node_x)
    run_step(node_x)
  }
  merge_step()
}

count_graph_original = set_fixed_sample_size(count_graph, sample_size_val = 5000)
count_graph = generate_count_graph_sample_size(count_graph_original)
count_graph = generate_phylogeny(count_graph)

job_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
setwd('/home/ssrikan2/data-kreza1/smriti/MF_Signal_Simulation')
save(count_graph, file = paste0('./three_cell_cg_ss5000/count_graph_', job_id, '.rda'))
