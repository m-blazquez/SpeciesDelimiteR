
plot_alongside_phylogeny <- function(integrated_results, phylogeny, outgroup = NULL, support_threshold = NULL, add_legends_for = NULL) {

  # Ensure support threshold is specified
  if (is.null(support_threshold)) {
    stop("ERROR: You must specify a support threshold value.")
  }

  # Remove outgroup sequences (if provided)
  if (!is.null(outgroup) && length(outgroup) > 0) {
    for (i in seq_along(outgroup)) {
      if (outgroup[i] %in% phylogeny$tip.label) {
        phylogeny <- ape::drop.tip(phylogeny, outgroup[i])
      } else {
        stop(paste("ERROR: Outgroup", outgroup[i], "not found in the phylogeny."))
      }
    }
  } else {
    message("No outgroup provided or outgroup vector is empty. Proceeding with the full tree.")
  }

  # Ensure proper structure of input
  if (!"Specimen" %in% colnames(integrated_results)) {
    stop("ERROR: The input dataframe must contain a 'Specimen' column.")
  }

  if (!all(sapply(integrated_results[, -1], is.factor))) {
    stop("ERROR: All columns except 'Specimen' must be factors.")
  }

  # Convert Specimen in row names and remove it as a variable
  row.names(integrated_results) <- integrated_results$Specimen
  integrated_results$Specimen <- NULL

  # Ensure tree order and reorder dataframe to fit the tip labels
  phylogeny <- ape::reorder.phylo(phylogeny, "cladewise")
  integrated_results <- integrated_results[phylogeny$tip.label, ]

  # Set larger margins
  par(mar = c(5, 4, 4, 18))  # Increased the right margin slightly

  # Initial tree plotting
  phytools::plotTree(phylogeny, plot = FALSE)
  obj <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)

  phytools::plotTree(phylogeny, lwd = 1,
                     ylim = c(0, obj$y.lim[2] * 1.3),
                     xlim = c(0, obj$x.lim[2] * 2.8),  # Extended xlim further for additional space
                     ftype = "off")

  obj <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  h <- max(obj$xx)
  fsize <- 1

  # Add support values and modify supported branches
  if (!is.null(phylogeny$node.label)) {
    for (i in 1:length(phylogeny$node.label)) {
      support <- as.numeric(phylogeny$node.label[i])
      if (!is.na(support) && support >= support_threshold) {
        edgenum <- which(phylogeny$edge[, 2] == (ape::Ntip(phylogeny) + i))
        x_end <- obj$xx[phylogeny$edge[edgenum, 2]]
        y_pos <- obj$yy[phylogeny$edge[edgenum, 2]]
        x_start <- x_end - (x_end - obj$xx[phylogeny$edge[edgenum, 1]])

        segments(x_start, y_pos, x_end, y_pos, lwd = 3, col = "black")
        text(mean(c(x_start, x_end)), y_pos + 1,
             labels = round(support, 2), cex = 0.8, col = "navy", font = 2)
      }
    }
  } else {
    warning("No node labels detected. Support values will not be displayed on the tree.")
  }

  # Draw dotted lines and tip labels
  for (i in 1:ape::Ntip(phylogeny)) {
    lines(c(obj$xx[i], h), rep(obj$yy[i], 2), lty = "dotted")
    text(h, obj$yy[i], phylogeny$tip.label[i], cex = fsize, pos = 4, font = 3, offset = 0.1)
  }

  # Assign unique colors per variable
  s <- max(fsize * strwidth(phylogeny$tip.label))
  start.x <- 1.2 * h + s
  cols <- list()

  # Loop through each specified variable and plot the associated colored columns
  for (i in 1:ncol(integrated_results)) {
    col_name <- gsub("_", " ", colnames(integrated_results)[i])
    text(start.x, max(obj$yy) + 2, col_name, pos = 4, srt = 90, cex = 1, offset = 1)

    # Create color palette for each variable
    num_levels <- length(levels(integrated_results[[i]]))
    cols[[i]] <- setNames(randomcoloR::distinctColorPalette(num_levels), levels(integrated_results[[i]]))

    # Draw the colored boxes
    for (j in 1:nrow(integrated_results)) {
      xy <- c(start.x, obj$yy[j])
      y <- c(xy[2] - 0.5, xy[2] + 0.5, xy[2] + 0.5, xy[2] - 0.5)
      asp <- (par()$usr[2] - par()$usr[1]) / (par()$usr[4] - par()$usr[3]) * par()$pin[2] / par()$pin[1]
      x <- c(xy[1] - 0.5 * asp, xy[1] - 0.5 * asp, xy[1] + 0.5 * asp, xy[1] + 0.5 * asp)
      polygon(x, y, col = cols[[i]][as.character(integrated_results[[i]][j])])
    }

    # Increment x position for the next column
    start.x <- start.x + 2 * asp
  }

  # Once all columns are plotted, draw legends
  if (!is.null(add_legends_for)) {
    legend_x <- start.x + 0.5 * asp  # Positioned closer to columns
    legend_y <- max(obj$yy)

    for (i in seq_along(add_legends_for)) {
      col_name <- add_legends_for[i]
      idx <- which(gsub("_", " ", colnames(integrated_results)) == col_name)
      if (length(idx) > 0) {
        legend(legend_x, legend_y, legend = names(cols[[idx]]),
               fill = cols[[idx]], title = col_name, cex = 0.8, bty = "n", xpd = TRUE)
        legend_y <- legend_y - 10  # Adjust spacing between legends as needed
      }
    }
  }
}



