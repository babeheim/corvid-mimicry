

gradient_maker <- function(start=NA, stop=NA, cols=c("darkorange", "white", "darkcyan"), vis=FALSE, n_bins=1000){
    if(is.na(start) | is.na(stop)) stop("need to specify start and stop points on a numerical scale")
    colfunc <- colorRampPalette(cols)
    color.list <- colfunc(n_bins)
    color.locations <- seq(start, stop, length=n_bins)
    names(color.locations) <- color.list
    if(vis==TRUE) plot(color.locations, rep(1, n_bins), col=color.list, pch="|", ylim=c(0.9, 1.1), cex=5)
    return(color.locations)
}

data_gradient <- function(data, colors=c("darkorange", "white", "darkcyan"), my.start=NA, my.stop=NA){
    if(is.na(my.start)) my.start <- min(data, na.rm=TRUE)
    if(is.na(my.stop)) my.stop <- max(data, na.rm=TRUE)
    my.gradient <- gradient_maker(start=my.start, stop=my.stop, cols=colors)
    if(any(data > max(my.gradient), na.rm=T) | any(data < min(my.gradient), na.rm=T)) warning("data is not within gradient range")
    data.colors <- rep(NA, length(data))
    for(i in 1:length(data)){
        if(!is.na(data[i])) data.colors[i] <- names(my.gradient)[which.min(abs(data[i]-my.gradient))]
    }
    data.colors
}

col_alpha <- function(acol, alpha = 0.2) {
  acol <- col2rgb(acol)
  acol.red <- acol["red", ] / 255
  acol.green <- acol["green", ] / 255
  acol.blue <- acol["blue", ] / 255
  acol <- mapply(
    function(red, green, blue, alphas) {
      rgb(red, green, blue, alphas)
    },
    acol.red, acol.green, acol.blue, alpha
  )
  return(as.character(acol))
}


data_gradient2 <- function(data, lower_colors = c("darkorange", "white"), middle_color = NA, upper_colors = c("white", "darkcyan"), my.start = NA, middle = 0, my.stop = NA, n_bins = 1000) {
    if(is.na(my.start)) my.start <- min(data, na.rm=TRUE)
    if(is.na(my.stop)) my.stop <- max(data, na.rm=TRUE)
    lower_gradient <- gradient_maker(start=my.start, stop=middle, cols=lower_colors, n_bins = n_bins)
    upper_gradient <- gradient_maker(start=middle, stop=my.stop, cols=upper_colors, n_bins = n_bins)
    if (is.na(middle_color)) {
      my.gradient <- c(lower_gradient, upper_gradient)
    } else {
      names(middle) <- middle_color
      my.gradient <- c(lower_gradient, middle, upper_gradient)
    }
    if(any(data > max(my.gradient), na.rm=T) | any(data < min(my.gradient), na.rm=T)) warning("data is not within gradient range")
    data.colors <- rep(NA, length(data))
    for(i in 1:length(data)){
        if(!is.na(data[i])) {
          if (data[i] < middle) {
            data.colors[i] <- names(lower_gradient)[which.min(abs(data[i]-lower_gradient))]
          } else {
            data.colors[i] <- names(upper_gradient)[which.min(abs(data[i]-upper_gradient))]

          }
        }
    }
    data.colors
}
