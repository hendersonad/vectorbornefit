#' Background theme for ggplots for Zika paper
#'
#' Background ggplot theme for all plots for Zika paper
#' @param base_size default = 11
#' @param base_family default = ""
#' @export

theme_zika_fiji <- function (base_size = 11, base_family = ""){
  half_line <- base_size/2
  theme(title = element_text(size = 12, face = "bold"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line( size=.1, color=rgb(0,0,0,0.3) ) ,
        axis.ticks = element_blank(),
        line = element_line(colour = "gray70"),
        rect = element_rect(colour = "black", fill = "white")
  )
}

