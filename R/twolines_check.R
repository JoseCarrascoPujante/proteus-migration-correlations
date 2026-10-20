# Two-lines test of U-shaped relations (Simonsohn 2018, Adv Methods Pract Psychol Sci 1:538,
# doi:10.1177/2515245918805755) for every metric pair, both directions, three scenarios, corrected as one
# family like the RESET and Lind-Mehlum tests of dependence_rule.m.
# twolines() is the author's implementation in his CRAN package statuser (0.4.0 was used,
# doi:10.32614/CRAN.package.statuser). Usage: Rscript R/twolines_check.R <repository root>

.libPaths(c(file.path(commandArgs(trailingOnly = TRUE)[1], "R", "library"), .libPaths()))
library(statuser)

projectRoot <- commandArgs(trailingOnly = TRUE)[1]
derived <- file.path(projectRoot, "data", "derived")

dep <- read.csv(file.path(derived, "pairwiseDependence.csv"), stringsAsFactors = FALSE)
dep <- dep[order(dep$Rank), ]

metricFiles <- list.files(derived, pattern = "^perTrajectory_.+\\.csv$")
metricOf <- sub("^perTrajectory_(.+)\\.csv$", "\\1", metricFiles)
readMetric <- function(name) {
  file <- metricFiles[metricOf == name]
  stopifnot(length(file) == 1)
  read.csv(file.path(derived, file), check.names = FALSE)
}
scenarios <- c(NoStimuli = "No_stimuli", Galvanotaxis = "Galvanotaxis", Chemotaxis = "Chemotaxis")

rows <- list()
for (i in seq_len(nrow(dep))) {
  ab <- strsplit(dep$Pair[i], " vs ")[[1]]
  valuesA <- readMetric(ab[1])
  valuesB <- readMetric(ab[2])
  for (sc in names(scenarios)) {
    for (d in 1:2) {
      pred <- if (d == 1) ab[1] else ab[2]
      resp <- if (d == 1) ab[2] else ab[1]
      tlData <- data.frame(
        tlPredictor = if (d == 1) valuesA[[scenarios[sc]]] else valuesB[[scenarios[sc]]],
        tlResponse  = if (d == 1) valuesB[[scenarios[sc]]] else valuesA[[scenarios[sc]]])
      res <- twolines(tlResponse ~ tlPredictor, graph = 0, data = tlData, quiet = TRUE)
      rows[[length(rows) + 1]] <- data.frame(
        Pair = dep$Pair[i], Condition = sc, Direction = c("AxB", "BxA")[d],
        Fit = paste(resp, "on", pred), Breakpoint = unname(res$xc),
        b1 = res$b1, z1 = res$z1, p1 = res$p1, b2 = res$b2, z2 = res$z2, p2 = res$p2,
        Ushape_uncorrected = res$u.sig, Note = trimws(res$msg))
    }
  }
}
out <- do.call(rbind, rows)
out$pU <- ifelse(out$b1 * out$b2 < 0, pmax(out$p1, out$p2), 1)
out$qU <- p.adjust(out$pU, method = "BH")
write.csv(out, file.path(derived, "twolines_check.csv"), row.names = FALSE)
cat(sprintf("Two-lines check: %d fits, %d U-shaped by the author's criterion, %d at BH q < 0.05 -> %s\n",
            nrow(out), sum(out$Ushape_uncorrected), sum(out$qU < 0.05), file.path(derived, "twolines_check.csv")))
