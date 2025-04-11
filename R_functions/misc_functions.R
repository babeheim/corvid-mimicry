
na_to_999 <- function(x) {
  x[which(is.na(x))] <- 999
  return(x)
}

dir_init <- function (path, verbose = FALSE) {
    if (substr(path, 1, 2) != "./") 
        stop("path argument must be formatted\n    with \"./\" at beginning")
    contents <- dir(path, recursive = TRUE)
    if (verbose) {
        if (length(contents) == 0) 
            print(paste("folder ", path, " created.", sep = ""))
        if (length(contents) > 0) 
            print(paste("folder ", path, " wiped of ", length(contents), 
                " files/folders.", sep = ""))
    }
    if (dir.exists(path)) 
        unlink(path, recursive = TRUE)
    dir.create(path)
}

write_pandoc_table <- function(df, path, style = "rmarkdown", trim_ws = TRUE, ...) {
  df |>
  pandoc.table(style = style, ...) |>
  capture.output() -> out
  if (trim_ws) {
    out <- out[which(out != "")]
  }
  writeLines(out, path)
}

pr_hidden_single <- function(R, post) {
  (1 - post$rho)^R * post$q / ((1-post$rho)^R * post$q + (1 - post$q))
}

pr_hidden_double <- function(R1, R2, post) {
  (1 - post$rho1)^R1 * (1 - post$rho2)^R2 * post$q / ((1-post$rho1)^R1 * (1-post$rho2)^R2 * post$q + (1 - post$q))
}

pr_hidden_single_true <- function(R, q, rho) {
  (1 - rho)^R * q / ((1 - rho)^R * q + (1 - q))
}

pr_hidden_double_true <- function(R1, R2, q, rho1, rho2) {
  (1 - rho1)^R1 * (1 - rho2)^R2 * q / ((1-rho1)^R1 * (1-rho2)^R2 * q + (1 - q))
}

warnifnot <- function(test) {
  if (!test) {
    warning(deparse(substitute(test)))
  }
}

