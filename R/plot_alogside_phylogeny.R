
plot_alongside_phylogeny <- function(integrated_results, phylogeny, outgroup = NULL, support_threshold = NULL, add_legends_for = NULL) {

  # Ensure support threshold is specified
  if (is.null(support_threshold)) {

    stop("ERROR: You must specify a support threshold value.")

  }

  # Remove outgroup sequences (if provided)
  if (!is.null(outgroup) && length(outgroup) > 0) {

    for (i in seq_along(outgroup)) {

      if (outgroup[i] %in% phylogeny$tip.label) {

        phylogeny <- drop.tip(phylogeny, outgroup[i])

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
  phylogeny <- reorder(phylogeny, "cladewise")
  integrated_results <- integrated_results[phylogeny$tip.label, ]

  # Initial tree plotting
  plotTree(phylogeny, plot = FALSE) # CHECK FROM WHICH PACKAGE THE FUNTION plotTree CAME AND ADD IT TO THE DESCRIPTION!!!
  obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)

  plotTree(phylogeny, lwd = 1,
           ylim = c(0, obj$y.lim[2] * 1.05),
           xlim = c(0, obj$x.lim[2] * 1.5),
           ftype = "off")

  obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  h <- max(obj$xx)
  fsize <- 1

  # Add support values and modify supported branches
  if (!is.null(phylogeny$node.label)) {

    for (i in 1:length(phylogeny$node.label)) {

      support <- as.numeric(phylogeny$node.label[i])

      if (!is.na(support) && support >= support_threshold) {

        edgenum <- which(phylogeny$edge[, 2] == (Ntip(phylogeny) + i))
        lines(obj$xx[phylogeny$edge[edgenum, ]],
              obj$yy[phylogeny$edge[edgenum, ]],
              lwd = 2, col = "black")
        text(mean(obj$xx[phylogeny$edge[edgenum, ]]),
             mean(obj$yy[phylogeny$edge[edgenum, ]]),
             labels = support, cex = 0.7, col = "red")

      }

    }

  }

  # Draw dotted lines and tip labels
  for (i in 1:Ntip(phylogeny)) {

    lines(c(obj$xx[i], h), rep(obj$yy[i], 2), lty = "dotted")
    text(h, obj$yy[i], phylogeny$tip.label[i], cex = fsize, pos = 4, font = 3, offset = 0.1)

  }

  # Assign unique colors per variable
  s <- max(fsize * strwidth(phylogeny$tip.label))
  start.x <- 1.05 * h + s
  cols <- list()

  for (i in 1:ncol(integrated_results)) {

    text(start.x, max(obj$yy) + 1, paste(colnames(integrated_results)[i]), pos = 4, srt = 60, cex = 0.8, offset = 0)
    num_levels <- length(levels(integrated_results[[i]]))
    cols[[i]] <- setNames(distinctColorPalette(num_levels), levels(integrated_results[[i]]))

    for (j in 1:nrow(integrated_results)) {

      xy <- c(start.x, obj$yy[j])
      y <- c(xy[2] - 0.5, xy[2] + 0.5, xy[2] + 0.5, xy[2] - 0.5)
      asp <- (par()$usr[2] - par()$usr[1]) / (par()$usr[4] - par()$usr[3]) * par()$pin[2] / par()$pin[1]
      x <- c(xy[1] - 0.5 * asp, xy[1] - 0.5 * asp, xy[1] + 0.5 * asp, xy[1] + 0.5 * asp)
      polygon(x, y, col = cols[[i]][as.character(integrated_results[[i]][j])])

    }

    start.x <- start.x + 2 * asp

  }

  # Add legends if specified
  if (!is.null(add_legends_for)) {

    legend_x <- start.x + 2 * asp
    legend_y <- max(obj$yy)

    for (var in add_legends_for) {

      if (var %in% colnames(integrated_results)) {

        legend(legend_x, legend_y, legend = names(cols[[var]]),
               fill = cols[[var]], title = var, cex = 0.8, bty = "n")
        legend_y <- legend_y - length(cols[[var]]) * 0.5

      }

    }

  }

}

