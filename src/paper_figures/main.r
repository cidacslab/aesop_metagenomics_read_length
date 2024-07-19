
################################# COMMON SETUP #################################

source("src/paper_figures/0-setup/load_libraries.r")
source("src/paper_figures/0-setup/load_metrics_report_file.r")
source("src/paper_figures/0-setup/load_metrics_reports_all.r")

df_metrics <- load_metrics_reports_all()

results_folder_root <- "results/paper_figures/"
results_folder <- results_folder_root
dir.create(results_folder, recursive = TRUE)

################################## STATISTICS ##################################

source("src/paper_figures/3-statistics/statistic_analysis.r")

################################### FIGURE 1 ###################################

source("src/paper_figures/2-plots/main_plot_boxplots.r")


################################### FIGURE 2 ###################################

source("src/paper_figures/2-plots/main_plot_scatterplot.r")


################################### FIGURE 3 ###################################

source("src/paper_figures/2-plots/main_plot_heatmap.r")


################################################################################