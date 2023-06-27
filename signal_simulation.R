library(ape)
library(qfm)
library(tidyverse)
library(igraph)
library(ggraph)
setwd('~/Documents/qfm2')
devtools::load_all()
setwd('~/Documents/R/QFM')
devtools::load_all()
global_step_size = 0.01
global_target_time = 9.0

cell_state_meta = read_csv("./metadata/signal_record_test.csv")
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

count_graph = set_fixed_sample_size(count_graph, sample_size_val = 2000)
count_graph = generate_count_graph_sample_size(count_graph)
count_graph = generate_phylogeny(count_graph)

tr = phylo_edges_to_tr(count_graph$phylo_edges)

# define a signal curve
signal_func = function(x_vec) {
  1e-9 * ((1/6)* pmin(pmax(0, x_vec), 4) + (-2/15) * pmax(0, x_vec-4))
}
x_vec = seq(0, global_target_time, by = 0.05)
y_vec = signal_func(x_vec)
plot(x_vec, y_vec, type = "l")

# somatic with signal
tr_dd = list_dd_and_tips_mod2(tr)$dd
tr_node_time = count_graph$phylo_edges$length
tr_node_start = count_graph$phylo_edges$in_time
tr_node_end = count_graph$phylo_edges$out_time
names(tr_node_time) = names(tr_node_start) = names(tr_node_end) = count_graph$phylo_edges$out_node

tr_tips = list_dd_and_tips_mod2(tr)$tips
make_mut <- function() {
  paste0(sample(c(letters, 1:9), size = 20, replace = T), collapse = "")
}



cell_mut_tb = bind_rows(map(names(tr_tips), function(node_id) {
  # constant rate
  # Lambda = mut_rate * tr_node_time[[node_id]]
  # signal
  Lambda = integrate(signal_func, tr_node_start[node_id], tr_node_end[node_id])$value
  n_mut = rbinom(1, 10^9, 1-exp(-Lambda))
  # this version value is number of occurrences along the edge
  # tibble(cell = tr_tips[[node_id]],
  #        mut = make_mut(),
  #        value = n_mut)
  # this version generates one row for each mut that occurred
  as_tibble(
    expand.grid(cell = tr_tips[[node_id]],
                mut = map_chr(1:n_mut, function(i) {
                  make_mut()
                }))
  )
}))

# occurrences <- table(cell_mut_tb$mut)
# subset_mut <- occurrences[occurrences > 1]
# cell_mut_tb <- subset(cell_mut_tb, mut %in% names(subset_mut))
cell_mut_tb$value = 1
m <- spread(cell_mut_tb, mut, value, fill = 0)
# m$type <- m$cell
# m <- m[,c(1,16994,2:16993)]
chr_mat = as.matrix(m[-1])
mf_vec = colMeans(chr_mat)