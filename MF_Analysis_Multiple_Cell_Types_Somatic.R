library(ape)
library(qfm)
library(tidyverse)
library(igraph)
library(ggraph)
setwd('~/Documents/qfm2')
devtools::load_all()

mut_p = readRDS("./metadata//mut_p_marc1.rds")
mut_p$mut_rate_func = map(1:length(mut_p$mut_rate), function(i) {
  signal_func_hgRNA
})
mut_p$mut_rate = NULL

setwd("~/Documents/R/MosaicFractionAnalysis/R")
devtools::load_all()


setwd('~/Documents/R/MF_Signal_Simulation/')

load('phy.rda')

signal_func_hgRNA = function(x_vec) {
  val = 5e-1 * ((1/6)* pmin(pmax(0, x_vec-1), 4) + (-1/6) * pmax(0, x_vec-5))
  out_vec = pmax(val, 0)
  #out_vec[x_vec < 1] = 0
  out_vec
}
signal_func = function(x_vec) {
  val = 1e-9 * ((1/6)* pmin(pmax(0, x_vec-1), 4) + (-1/6) * pmax(0, x_vec-5))
  out_vec = pmax(val, 0) #+ 2*10^(-9)/3
  #out_vec[x_vec < 1] = 0
  out_vec
}

raw_data = load_simulated_data('./one_cell_somatic_mut_ss5000_015_2_signal4/', 100, mat_name = 'mf_table')

raw_data = annotate_mut_node(raw_data, phy)
nested_data = nest_raw(raw_data, c('ID', 'sample'))
breaks <- c(-Inf, seq(-13.5, -0.5, by = 1), 1)
nested_data = append_density(nested_data, breaks)

#correct cum_mf
node_order_list = list('TypeS'=1, 'TypeT'=2, 'TypeC'=3, 'TypeA'=4, 'TypeB'=5)
nested_data = increment_cum_mf(nested_data, node_order_list)

#cumulative adjustment, 1/MF normalization, and averaging across node type
nested_data = normalize_density_list(nested_data, apply_cum_adjustment = FALSE)
nested_data$nested$condition = nested_data$nested$ID
nested_density_raw = average_density_by(nested_data, c('condition'), raw = TRUE)
nested_density = average_density_by(nested_data, c('condition'), normalized = TRUE)


nested_density = convert_mf_to_time(nested_density)
nested_density_raw = convert_mf_to_time(nested_density_raw)


doubling_time = 0.65
time_vec = seq(0,length(breaks)*doubling_time, doubling_time)
time_df = align_time_scale(nested_density, time_vec)
log2mf_df = align_counts(nested_density)
# time_df_raw = align_time_scale(nested_density_raw, time_vec)
# time_df_comb = left_join(time_df, time_df_raw, by = 'time_vec')
time_df_norm = time_df %>% select(-time_vec) %>%
  mutate(truth = signal_func(time_vec)) %>%
  mutate(across(everything(), ~ ./sum(.)))
# time_df_norm$truth = signal_func(time_vec)
# time_df_norm = time_df_norm/colSums(time_df_norm)
time_df_norm$time_vec = time_vec

plot_tb = gather(data = time_df_norm %>% select(time_vec, average, average_nonzero, truth),
                 key = 'condition', value = 'count', -time_vec)
ggplot(data = plot_tb) +
  geom_line(aes(x = time_vec, y = count, group = condition, color = condition)) +
  ylab('Mosaic Fraction Density') +
  xlab('Time (days)')

sqrt(sum((time_df_norm$truth - time_df_norm$average_nonzero)^2))
sqrt(sum((time_df_norm$truth - time_df_norm$average)^2))

plot_tb = nested_density_raw$nested_density %>% select(condition, avg_density)
plot_tb$table = map(plot_tb$avg_density, function(density) {
  tibble(time_bin = density$table$bin, count = density$table$count)
})
plot_tb = select(plot_tb, condition, table) %>% unnest(cols = table)
ggplot(data = plot_tb) +
  geom_line(aes(x = time_bin, y = count, group = condition, color = condition))
