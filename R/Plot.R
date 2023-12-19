#' @title Plot the observed and expected number in each condition
#'
#' 
#' @param data A data frame with the following columns:
#' \itemize{
#'   \item{condition: }{The condition names.}
#'   \item{observe: }{The observed number in each condition.}
#'   \item{expect: }{The expected number in each condition.}
#'   \item{z_score: }{The z-score of the observed number.}
#'   \item{log2FC: }{The log2 fold change of the observed number.}
#'   \item{upper.p: }{The upper p-value of the observed number.}
#'   \item{lower.p: }{The lower p-value of the observed number.}
#' }
#' @param width The width of the bar.
#' @param observe.col The color of the observed bar.
#' @param expect.col The color of the expected bar.
#' @return A ggplot object.
#' 
#' @export
ObsExpBarPlot <- function(data, width=0.75, observe.col="red", expect.col="grey"){

  # check input data
  if (!is.data.frame(data)) {
    stop("The input data is not a data frame.")
  }

  # check column names
  check_col_chr <- c(
    "condition", "observe", "expect", "z_score", "log2FC", "upper.p", "lower.p")

  if (!all(colnames(data)%in%check_col_chr)) {
    stop("The column names of the input data frame are not correct.")
  }

  ## transform data
  plot_dat <- lapply(c("observe", "expect"), function(x){
  
    lapply(1:length(data$condition), function(y){
      data.frame(
        condition = data[y,"condition"],
        number    = data[y, x],
        type      = x)
    }) %>% Reduce(rbind, .)
  
  }) %>% Reduce(rbind, .)

  plot_dat$condition <- factor(plot_dat$condition, levels=data$condition)
  plot_dat$type <- factor(plot_dat$type, levels=c("observe", "expect"))

  ggplot2::ggplot(plot_dat, ggplot2::aes(x=condition, y=number, fill=type)) +
    ggplot2::scale_fill_manual(values = c(observe.col, expect.col)) +
    ggplot2::geom_bar(stat="identity", position=ggplot2::position_dodge(), color="black", width=width)

}

