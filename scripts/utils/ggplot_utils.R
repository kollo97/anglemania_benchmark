big_text_theme <- function() {
    # increase font size of all text
    theme_bw() +
    theme(
        text = element_text(size = 30, face = "bold"),
        aspect.ratio = 1,
        plot.title = element_text(size = 20),
        legend.position = c(0.6, 0.15),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.box = "horizontal"
    )
}