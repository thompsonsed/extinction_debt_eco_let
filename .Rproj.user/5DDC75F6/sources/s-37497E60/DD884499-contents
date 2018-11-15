library(tidyverse)
plot_colours <- c("#A53A2C",
                  "#BD665E",
                  "#D18D89",
                  "#E0B2AF",
                  "#ECD3D1",
                  "#F4EEED",
                  "#EFEFF4",
                  "#D8D7EB",
                  "#BCB9DF",
                  "#9E98D2",
                  "#7D75C5",
                  "#5C4DBC"
)
percent_cover_names = c("0.1"="10% habitat cover", "0.2"="20% habitat cover",
                        "0.4"="40% habitat cover")
size_names = as_labeller(c("50" = "Area == 50^2",
                           "500" = "Area == 500^2", "5000" = "Area == 5000^2"), default=label_parsed)
sigma_names = as_labeller(setNames(sapply(c(2, 4,8, 16, 32), FUN = function(x){paste("sigma == ",x)}), 
                                   sapply(c(2, 4, 8, 16, 32), FUN=function(x) as.character(x))),
                          default=label_parsed)
a_max_names <- as_labeller(c("100" = "A[max] == 100^2",
                             "1000" = "A[max] == 1000^2", "10000" = "A[max] == 10000^2"),
                           default=label_parsed)
