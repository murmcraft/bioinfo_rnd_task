lineplot <- function(df, xlim, ylim, title, xlab, ylab, cols, labels) {
  p <- ggplot(df, aes(x = x, 
                      y = value,
                      color = variable,
                      fill = variable)) +
    geom_line() +
    scale_x_continuous(expand = c(0,0),
                       limits = xlim) + 
    scale_y_continuous(expand = c(0,0),
                       limits = ylim) +
    scale_color_manual(values = cols, 
                       labels = labels) +
    labs(title = title, x = xlab, y = ylab) + 
    plottheme
  return(p)
}

barplot <- function(df, xaxis, yaxis, title, xlab, ylab, cols, labels) {
  p <- ggplot(df, aes(x = x,
                      y = value,
                      color = "black",
                      fill = variable)) +
    geom_bar(stat = "identity") +
    xaxis +
    yaxis +
    scale_color_manual(values = "black", 
                       labels = labels) +
    scale_fill_manual(values = cols, 
                       labels = labels) +
    labs(title = title, x = xlab, y = ylab) + 
    plottheme
  return(p)
}

histogram <- function(df, xaxis, yaxis, title, xlab, ylab, cols, labels, bins) {
  p <- ggplot(df, aes(x = x,
                      color = "black")) +
    geom_histogram(bins = bins) +
    xaxis +
    yaxis +
    scale_color_manual(values = "black", 
                       labels = labels) +
    scale_fill_manual(values = cols, 
                      labels = labels) +
    labs(title = title, x = xlab, y = ylab) + 
    plottheme
  return(p)
}
densityplot <- function(df, title, xlab, ylab, cols, labels, alpha, adjust) {
  p <- ggplot(df, 
              aes(x = value,
                  color = "black",
                  fill = variable)) +
    geom_density(alpha = alpha, adjust = adjust) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = cols, labels = labels) +
    scale_color_manual(values = "black", guide = "none") +
    labs(title = title, x = xlab, y = ylab) + 
    plottheme
  return(p)
}