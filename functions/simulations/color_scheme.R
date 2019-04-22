colors_basis <- gg_color_hue(n=3)
colors_simulation_adjust.batch <- c("Original" = "black",
                                    "ComBat corrected" = colors_basis[3],
                                    "MMUPHin corrected" = colors_basis[1])
colors_simulation_unsupervised <- c("Original" = "black",
                                    "MMUPHin corrected" = colors_basis[1])
