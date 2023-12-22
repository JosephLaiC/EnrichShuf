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

#' @title Plot the observed and expected number in each condition
#'
#'
#' @param dataList A list of data frames with the following columns:
#' \itemize{
#'   \item{condition: }{The condition names.}
#'   \item{observe: }{The observed number in each condition.}
#'   \item{expect: }{The expected number in each condition.}
#'   \item{z_score: }{The z-score of the observed number.}
#'   \item{log2FC: }{The log2 fold change of the observed number.}
#'   \item{upper.p: }{The upper p-value of the observed number.}
#'   \item{lower.p: }{The lower p-value of the observed number.}
#' }
#' @param color A vector of colors for each condition.
#' @param xlab The x-axis label.
#' @param ylab The y-axis label.
#' @param hline The y-axis value of the horizontal line.
#' @param hline.type The type of the horizontal line.
#' @param hline.color The color of the horizontal line.
#' @param line.width The width of the line.
#' @return A ggplot object.
#'
#' @export
ObsExpCurvePlot <- function(
  dataList, color=NULL, 
  xlab="LogFC (expected/observe)",
  ylab="Distance of feature to nearest genomic regions",
  hline=0, hline.type="solid", hline.color="black",
  line.width=1
) {
  
  if (!typeof(dataList)=="list") {
    stop("The input data is not a list.")
  }

  if (!all(sapply(dataList, is.data.frame))) {
    stop("The input data is not a list of data frames.")
  }

  if (is.null(names(dataList))) {
    stop("The input data list does not have names.")
  }

  ## check column names
  check_col_chr <- c(
    "condition", "observe", "expect", "z_score", "log2FC", "upper.p", "lower.p")
  
  check_bool <- sapply(
    dataList, function(x){
      all(check_col_chr %in% colnames(x))
    }
  )

  if (!all(check_bool)) {
    stop("The column names of the input data frame are not correct.")
  }

  ## check condition names are the same
  if (length(dataList) > 1) {
    check_bool <- sapply(
      dataList, function(x){
        all(dataList[[1]]$condition %in% x$condition)
      }
    )
  } else {
    check_bool <- TRUE
  }

  if (!all(check_bool)) {
    stop("The condition names of the input data frames are not the same.")
  }
  
  condition_chr <- colnames(dataList[[1]])

  # check condition name not contain "intersect"
  if ("intersect" %in% condition_chr) {
    stop("The condition name 'intersect' is not allowed.")
  }

  distance_num <- sapply(
    stringr::str_split(dataList[[1]][,"condition"], "_"),
    function(x){
      as.numeric(x[2])
    }
  )

  ## transform data
  data_tbl <- lapply(names(dataList), function(x) {
    
    data.frame(
      distance = distance_num,
      log2FC   = pull(dataList[[x]], log2FC),
      type     = x)
  
  }) %>% Reduce(rbind, .)

  plot_dat <- ggplot(
    data = data_tbl, aes(x = distance, y = log2FC, color = type)
  ) + 
    geom_line(linewidth = line.width) +
    geom_hline(
      yintercept = hline, linetype = hline.type, color = hline.color
    )
  
  if (!is.null(color)) {

    if (!length(color)==length(names(dataList))) {
      stop("The length of the color vector is not the same as the number of conditions.")
    }

    plot_dat <- plot_dat + scale_color_manual(values = setNames(color, names(dataList)))

  }

  plot_dat <- plot_dat + labs(x = xlab, y = ylab)
  
  print(plot_dat)

}


