library(ape)
library(qfm)
library(tidyverse)
library(igraph)
library(ggraph)
setwd('/home/ssrikan2/data-kreza1/smriti/qfm2')
devtools::load_all()
setwd('/home/ssrikan2/data-kreza1/smriti/MF_Signal_Simulation')

job_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
load('param_tb.rda')

load(paste0('./output/',param_tb$filename[job_id]))

# define a signal curve
signal_func = function(x_vec) {
  val = 1e-9 * ((1/6)* pmin(pmax(0, x_vec), 4) + (-4/15) * pmax(0, x_vec-4))
  pmax(val, 0)
}


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

cell_mut_tb$value = 1
m <- spread(cell_mut_tb, mut, value, fill = 0)
# m$type <- m$cell
# m <- m[,c(1,16994,2:16993)]
chr_mat = as.matrix(m[-1])
mf_vec = colMeans(chr_mat)

save(tibble(sequence = names(mf_vec), mosaic_fraction = mf_vec), file = paste0('./output2/', param_tb$filename[job_id], '_', param_tb$num_sim[job_id], '.rda'))


