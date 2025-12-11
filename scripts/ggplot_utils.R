library(ggplot2)
library(stringr)

draw_heatmap <- function(
    df,
    size = 1 - 6e-4,
    score = "BTVR",
    fill = "mean",
    scale_fill_limits = NULL,
    labels = c("angl" = "anglemania genes", "hvg" = "highly variable genes"),
    add_geom_text = FALSE,
    ...
    ) {
    tmp_df <- df %>%
            group_by(gene_selection, batch.facLoc, de.facLoc) %>%
            summarise(
                # Summarize the chosen column
                across(
                    all_of(score),
                    list(
                        mean   = ~mean(.x, na.rm = TRUE),
                        median = ~median(.x, na.rm = TRUE),
                        sd     = ~sd(.x, na.rm = TRUE)
                    ),
                    # If you only have one variable, this will produce columns named "mean", "median", "sd"
                    .names = "{.fn}"
                ),
                n = n(),
                .groups = 'drop'
            )
    p <- ggplot(tmp_df, aes(x = as.factor(batch.facLoc), y = as.factor(de.facLoc), fill = !!sym(fill))) +
        geom_tile(width = size, height = size, ...) +
        coord_equal() +
        facet_wrap(~gene_selection, labeller = labeller(gene_selection = labels)) + # Facet by gene_selection
        scale_fill_gradient(low = "grey", high = "purple", limits = scale_fill_limits) +
        labs(
            title = sprintf("Heatmap of %s by Batch and Biological Variability", stringr::str_to_upper(fill)),
            x = "Strength of Technical Noise (batch.facLoc)",
            y = "Biological Variability (de.facLoc)",
            fill = sprintf("%s %s Score", stringr::str_to_upper(fill), stringr::str_to_upper(score))
        ) +
        theme_minimal() +
            theme(
                text = element_text(size = 15, face = "bold"),
                axis.title = element_text(size = 20),
                # increase facet label size
                strip.text = element_text(size = 20)
            )

    if (add_geom_text) {
        p <- p + geom_text(aes(label = round(!!sym(fill), 2)), color = "black", size = 4)
    }
    p
}

draw_heatmap2 <- function(
    df,
    size = 1 - 6e-4,
    score = "BTVR",
    fill = "mean",
    scale_fill_limits = NULL,
    labels = c("angl" = "anglemania genes", "hvg" = "highly variable genes"),
    add_geom_text = FALSE,
    ...
    ) {
    tmp_df <- df %>%
            filter(
                de.facLoc != 0.001,
                batch.facLoc != 0.001
            ) %>%
                mutate(
                    de.facLoc = case_when(
                        de.facLoc == 0 ~ "none",
                        de.facLoc == 0.01 ~ "small",
                        de.facLoc == 0.05 ~ "mild",
                        de.facLoc == 0.1 ~ "moderate",
                        de.facLoc == 0.3 ~ "strong",
                        de.facLoc == 0.6 ~ "very strong"
                    ),
                    de.facLoc = factor(de.facLoc, levels = c("none", "small", "mild", "moderate", "strong", "very strong")),
                    batch.facLoc = case_when(
                        batch.facLoc == 0 ~ "none",
                        batch.facLoc == 0.01 ~ "small",
                        batch.facLoc == 0.05 ~ "mild",
                        batch.facLoc == 0.1 ~ "moderate",
                        batch.facLoc == 0.3 ~ "strong",
                        batch.facLoc == 0.6 ~ "very strong"
                    ),
                    batch.facLoc = factor(batch.facLoc, levels = c("none", "small", "mild", "moderate", "strong", "very strong"))
                ) %>%
                group_by(gene_selection, batch.facLoc, de.facLoc) %>%
                    summarise(
                        # Summarize the chosen column
                        across(
                            all_of(score),
                            list(
                                mean   = ~ mean(.x, na.rm = TRUE),
                                median = ~ median(.x, na.rm = TRUE),
                                sd     = ~ sd(.x, na.rm = TRUE)
                            ),
                            # If you only have one variable, this will produce columns named "mean", "median", "sd"
                            .names = "{.fn}"
                        ),
                        n = n(),
                        .groups = "drop"
                    )
                
    p <- ggplot(tmp_df, aes(x = as.factor(batch.facLoc), y = as.factor(de.facLoc), fill = !!sym(fill))) +
        geom_tile(width = size, height = size, ...) +
        coord_equal() +
        facet_wrap(~gene_selection, labeller = labeller(gene_selection = labels)) + # Facet by gene_selection
        scale_fill_gradient(low = "grey", high = "purple", limits = scale_fill_limits) +
        labs(
            title = sprintf("Heatmap of %s by Batch and Biological Variability", stringr::str_to_upper(fill)),
            x = "Strength of Technical Noise (batch.facLoc)",
            y = "Biological Variability (de.facLoc)",
            fill = sprintf("%s %s Score", stringr::str_to_upper(fill), stringr::str_to_upper(score))
        ) +
        theme_minimal() +
            theme(
                text = element_text(size = 15, face = "bold"),
                axis.title = element_text(size = 20),
                # increase facet label size
                strip.text = element_text(size = 20),
                # rotate x-axis labels
                axis.text.x = element_text(angle = 30, hjust = 1)
            )

    if (add_geom_text) {
        p <- p + geom_text(aes(label = round(!!sym(fill), 2)), color = "black", size = 4)
    }
    p
}

draw_heatmap_thresholds_vis <- function(
    df,
    size = 1 - 6e-4,
    score = "BTVR",
    fill = "mean",
    scale_fill_limits = NULL,
    labels = c("angl" = "anglemania genes"),
    add_geom_text = FALSE,
    ...
){
    tmp_df <- df %>%
        group_by(zscore_mean_threshold, zscore_snr_threshold) %>%
        summarise(
            # Summarize the chosen column
            across(
                all_of(score),
                list(
                    mean = ~ mean(.x, na.rm = TRUE),
                    median = ~ median(.x, na.rm = TRUE),
                    sd = ~ sd(.x, na.rm = TRUE)
                ),
                # If you only have one variable, this will produce columns named "mean", "median", "sd"
                .names = "{.fn}"
            ),
            n = n(),
            .groups = "drop"
        )
    
        fill_legend_text <- if(fill == "n"){
            "Number of Genes"
        } else {
            sprintf("%s %s Score", stringr::str_to_upper(fill), stringr::str_to_upper(score))
        }

        p <- ggplot(
            tmp_df,
            aes(
                x = as.factor(zscore_mean_threshold),
                y = as.factor(zscore_snr_threshold),
                fill = !!sym(fill)
            )
        ) +
            geom_tile(width = size, height = size, ...) +
            coord_equal() +
            scale_fill_gradient(
                low = "grey",
                high = "purple",
                limits = scale_fill_limits
            ) +
            labs(
                title = sprintf(
                    "Heatmap of %s by Zscore Thresholds",
                    stringr::str_to_upper(fill)
                ),
                x = "Zscore Mean Threshold",
                y = "Zscore SNR Threshold",
                fill = fill_legend_text
            ) +
            theme_minimal() +
            theme(
                text = element_text(size = 15, face = "bold"),
                axis.title = element_text(size = 20),
                # increase facet label size
                strip.text = element_text(size = 20),
                # rotate x-axis labels
                axis.text.x = element_text(angle = 30, hjust = 1)
            )

        if (add_geom_text) {
            p <- p +
                geom_text(
                    aes(label = round(!!sym(fill), 2)),
                    color = "black",
                    size = 4
                )
        }
        p
}

draw_heatmap_dropout <- function(
    df,
    size = 1 - 6e-4,
    score = "BTVR",
    fill = "mean",
    scale_fill_limits = NULL,
    labels = c("angl" = "anglemania genes", "hvg" = "highly variable genes"),
    add_geom_text = FALSE,
    ...
) {
    tmp_df <- df %>%
        filter(
            de.facLoc != 0.001,
            batch.facLoc != 0.001
        ) %>%
        mutate(
            de.facLoc = case_when(
                de.facLoc == 0 ~ "none",
                de.facLoc == 0.01 ~ "small",
                de.facLoc == 0.05 ~ "mild",
                de.facLoc == 0.1 ~ "moderate",
                de.facLoc == 0.3 ~ "strong",
                de.facLoc == 0.6 ~ "very strong"
            ),
            de.facLoc = factor(
                de.facLoc,
                levels = c(
                    "none",
                    "small",
                    "mild",
                    "moderate",
                    "strong",
                    "very strong"
                )
            ),
            batch.facLoc = case_when(
                batch.facLoc == 0 ~ "none",
                batch.facLoc == 0.01 ~ "small",
                batch.facLoc == 0.05 ~ "mild",
                batch.facLoc == 0.1 ~ "moderate",
                batch.facLoc == 0.3 ~ "strong",
                batch.facLoc == 0.6 ~ "very strong"
            ),
            batch.facLoc = factor(
                batch.facLoc,
                levels = c(
                    "none",
                    "small",
                    "mild",
                    "moderate",
                    "strong",
                    "very strong"
                )
            )
        ) %>%
        group_by(gene_selection, batch.facLoc, de.facLoc) %>%
        summarise(
            # Summarize the chosen column
            across(
                all_of(score),
                list(
                    mean = ~ mean(.x, na.rm = TRUE),
                    median = ~ median(.x, na.rm = TRUE),
                    sd = ~ sd(.x, na.rm = TRUE)
                ),
                # If you only have one variable, this will produce columns named "mean", "median", "sd"
                .names = "{.fn}"
            ),
            n = n(),
            .groups = "drop"
        )

    p <- ggplot(
        tmp_df,
        aes(
            x = as.factor(batch.facLoc),
            y = as.factor(de.facLoc),
            fill = !!sym(fill)
        )
    ) +
        geom_tile(width = size, height = size, ...) +
        coord_equal() +
        facet_wrap(
            ~gene_selection,
            labeller = labeller(gene_selection = labels)
        ) + # Facet by gene_selection
        scale_fill_gradient(
            low = "grey",
            high = "purple",
            limits = scale_fill_limits
        ) +
        labs(
            title = sprintf(
                "Heatmap of %s by Batch and Biological Variability",
                stringr::str_to_upper(fill)
            ),
            x = "Strength of Technical Noise (batch.facLoc)",
            y = "Biological Variability (de.facLoc)",
            fill = sprintf(
                "%s %s Score",
                stringr::str_to_upper(fill),
                stringr::str_to_upper(score)
            )
        ) +
        theme_minimal() +
        theme(
            text = element_text(size = 15, face = "bold"),
            axis.title = element_text(size = 20),
            # increase facet label size
            strip.text = element_text(size = 20),
            # rotate x-axis labels
            axis.text.x = element_text(angle = 30, hjust = 1)
        )

    if (add_geom_text) {
        p <- p +
            geom_text(
                aes(label = round(!!sym(fill), 2)),
                color = "black",
                size = 4
            )
    }
    p
}

# Define a custom theme with large, easily readable text
theme_big_text <- function(base_size = 15, base_family = "", ...) {
    theme_minimal(base_size = base_size, base_family = base_family) +
        theme(
            plot.title = element_text(
                size = base_size * 1.4,
                face = "bold",
                hjust = 0.5
            ),
            plot.subtitle = element_text(
                size = base_size * 1.2,
                hjust = 0.5
            ),
            axis.title = element_text(
                size = base_size * 1.2,
                face = "bold"
            ),
            axis.text.y = element_text(
                size = base_size,
                color = "black"
            ),
            legend.title = element_text(
                size = base_size * 1.2,
                face = "bold"
            ),
            legend.text = element_text(
                size = base_size
            ),
            strip.text = element_text(
                size = base_size * 1.1,
                face = "bold"
            ),
            axis.text.x = element_text(
                size = base_size,
                color = "black",
                ...
            ),
            # Optionally, you can increase margins for better spacing:
            plot.margin = margin(15, 15, 15, 15)
        )
}

today <- function(){
    format(Sys.Date(), "%Y%m%d")
}
