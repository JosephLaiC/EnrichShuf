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

