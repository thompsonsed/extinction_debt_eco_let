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
                        "0.4"="40% habitat cover", "0.5"="50% habitat cover",
                        "0.7"="70% habitat cover", "0.9"="90% habitat cover")
percent_cover_names2 = c("0.2"="High habitat\nloss", "0.4"="Medium habitat\nloss",
                        "0.8"="Low habitat\nloss")
size_names = as_labeller(c("50" = "Area == 50^2",
                           "500" = "Area == 500^2", "5000" = "Area == 5000^2"), default=label_parsed)
sigma_names = as_labeller(setNames(sapply(c(2, 4,8, 16, 32), FUN = function(x){paste("sigma == ",x)}), 
                                   sapply(c(2, 4, 8, 16, 32), FUN=function(x) as.character(x))),
                          default=label_parsed)
a_max_names <- as_labeller(c("100" = "A[max] == 100^2",
                             "1000" = "A[max] == 1000^2", "10000" = "A[max] == 10000^2"),
                           default=label_parsed)
a_max_names_large <- as_labeller(c("10000" = "Local scale",
                                   "1e+06" = "Intermediate scale",
                                   "1e+08" = "Regional scale"))


#' Formats a label to 10^# notation
#'
#' @param x the numeric vector to format
#'
#' @return 10^# for each # in x
scientific_10 <- function(x) {
  p <-  parse(text=gsub("e[\\+]*", " %*% 10^", scales::scientific_format()(x)))
  p[x ==0] <- 0
  p[x == 1] <- 1
  return(p)
}

#' Formats a label to 10^# notation without a 1x multiplier
#'
#' @param x the numeric vector to format
#'
#' @return 10^# for each # in x
scientific_10_powers <- function(x) {
  p <-  parse(text=gsub("[1]e[\\+]*", "10^", scales::scientific_format()(x)))
  p[x ==0] <- 0
  p[x == 1] <- 1
  return(p)
}
