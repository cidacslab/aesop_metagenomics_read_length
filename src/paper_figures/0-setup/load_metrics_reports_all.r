
################################## LOAD DATA ###################################
load_metrics_reports_all <- function() {

    reports_folder_path <- "results/mocks_throat_based/performance_metrics")
    df_mocks <- load_metrics_report_file(reports_folder_path)

    reports_folder_path <- "results/mocks_throat/performance_metrics"
    df_new_mocks <- load_metrics_report_file(reports_folder_path)

    df_metrics <- bind_rows(df_mocks, df_new_mocks)

    df_groups <- df_metrics %>%
    group_by(sample_group, sample_name) %>%
    summarise(n = n())

    df_group_cat <- df_groups %>%
    group_by(sample_group) %>%
    summarise(n_samples = n())
    View(df_group_cat)

    return(df_metrics)
}
