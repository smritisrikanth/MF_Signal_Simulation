library(ape)
library(qfm)
library(tidyverse)
library(igraph)
library(ggraph)
setwd('/home/ssrikan2/data-kreza1/smriti/qfm2')
devtools::load_all()
setwd('/home/ssrikan2/data-kreza1/smriti/MF_Signal_Simulation')

job_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

global_step_size = 0.01
global_target_time = 9.0

cell_state_meta = read_csv("./signal_record_test.csv")
cell_state_col = cell_state_meta$color
names(cell_state_col) = cell_state_meta$name

cell_state_param_list = make_cell_state_param(cell_state_meta)
qfm_spec = list(
  founder_size = 1,
  root_id = "Default",
  node_id = "Default",
  tip_id = "Default",
  merge = list("Default" = NULL),
  cell_state_params = cell_state_param_list
)

global_step_size = 0.01
global_target_time = 9.0
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

count_graph = set_fixed_sample_size(count_graph, sample_size_val = 5000)

count_graph = generate_count_graph_sample_size(count_graph)
count_graph = generate_phylogeny(count_graph)

tr = phylo_edges_to_tr(count_graph$phylo_edges)

save(count_graph, tr, list = c('count_graph', 'tr'), file = paste0('./output_2/count_graph_', job_id, '.rda'))
