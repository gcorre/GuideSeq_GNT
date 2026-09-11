rbo_ext_fast2 <- function(L, S, p = 0.9) {
  
  lenL <- length(L)
  lenS <- length(S)
  
  if (length(S) > length(L)) {
    tmp <- L
    L <- S
    S <- tmp
    print("I changed list order as S was shorter than L")
  }
  
  l <- max(lenL, lenS)
  s <- min(lenL, lenS)
  
  # validations
  if (!is.numeric(p) || length(p) != 1 || is.na(p) || p <= 0 || p >= 1)
    stop("p doit être strictement compris entre 0 et 1.")
  
  # listes non vides
  if (length(L) == 0 || length(S) == 0)
    stop("Les listes ne peuvent pas être vides.")
  
  # unicité
  
  if (anyDuplicated(L))
    stop("L contient des éléments dupliqués.")
  
  if (anyDuplicated(S))
    stop("S contient des éléments dupliqués.")
  
  pd <- cumprod(rep(p, l))
  
  seenL <- new.env(hash = TRUE, parent = emptyenv())
  seenS <- new.env(hash = TRUE, parent = emptyenv())
  
  overlap <- 0L
  X <- integer(l)
  
  for (d in seq_len(l)) {
    
    if (d <= lenL) {
      x <- as.character(L[d])
      
      if (!exists(x, seenL, inherits = FALSE)) {
        assign(x, TRUE, seenL)
        if (exists(x, seenS, inherits = FALSE))
          overlap <- overlap + 1L
      }
    }
    
    if (d <= lenS) {
      x <- as.character(S[d])
      
      if (!exists(x, seenS, inherits = FALSE)) {
        assign(x, TRUE, seenS)
        if (exists(x, seenL, inherits = FALSE))
          overlap <- overlap + 1L
      }
    }
    
    X[d] <- overlap
  }
  
  d <- seq_len(l)
  
  term1 <- sum((X / d) * pd)
  
  term2 <- if (l > s) {
    dd <- (s + 1):l
    sum((X[s] * (dd - s) / (s * dd)) * pd[dd])
  } else 0
  
  term3 <- ((X[l] - X[s]) / l + X[s] / s) * pd[l]
  
  score <- ((1 - p) / p) * (term1 + term2) + term3
  structure(
    list(
      rbo = score,
      p = p,
      overlap_depth = X,
      agreement_depth = X / seq_len(l),
      weighted_agreement = (X / seq_len(l)) * pd,
      overlap_final = X[l],
      length_L = ifelse(lenL<lenS,lenS,lenL),
      length_S = ifelse(lenL<lenS,lenL,lenS)
    ),
    class = "rbo_ext"
  )
}


plot.rbo_ext <- function(
    x,
    type = c("all", "overlap", "agreement", "weighted"),
    ...
) {
  
  type <- match.arg(type)
  
  d <- seq_along(x$overlap_depth)
  
  if (type == "overlap") {
    
    plot(
      d,
      x$overlap_depth,
      type = "l",
      lwd = 2,
      col = "steelblue",
      xlab = "Depth",
      ylab = "Cumulative overlap",
      main = sprintf("RBO overlap (score = %.4f)", x$rbo)
    )
    abline(v = min(x$length_L, x$length_S),
           col = "grey",
           lty = 2)
    abline(coef = c(0,1), lty = 2, col = "grey")
    
    return(invisible(x))
  }
  
  if (type == "agreement") {
    
    plot(
      d,
      x$agreement_depth,
      type = "l",
      lwd = 2,
      col = "darkgreen",
      ylim = c(0, 1),
      xlab = "Depth",
      ylab = "Agreement (Xd/d)",
      main = sprintf("Agreement by depth (p = %.2f)", x$p)
    )
    abline(v = min(x$length_L, x$length_S),
           col = "grey",
           lty = 2)
    return(invisible(x))
  }
  
  if (type == "weighted") {
    
    plot(
      d,
      x$weighted_agreement,
      type = "l",
      lwd = 2,
      col = "firebrick",
      xlab = "Depth",
      ylab = expression((X[d]/d) %.% p^d),
      main = "Weighted contribution to RBO"
    )
    abline(v = min(x$length_L, x$length_S),
           col = "grey",
           lty = 2)
    return(invisible(x))
  }
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  par(mfrow = c(3, 1))
  
  plot(
    d,
    x$overlap_depth,
    type = "l",
    lwd = 2,
    col = "steelblue",
    xlab = "Depth",
    ylab = expression(X[d]),
    main = "Cumulative overlap"
  )
  abline(v = min(x$length_L, x$length_S),
         col = "grey",
         lty = 2)
  abline(coef = c(0,1), lty = 2, col = "grey")
  
  plot(
    d,
    x$agreement_depth,
    type = "l",
    lwd = 2,
    col = "darkgreen",
    ylim = c(0, 1),
    xlab = "Depth",
    ylab = expression(X[d]/d),
    main = "Agreement at depth"
  )
  abline(v = min(x$length_L, x$length_S),
         col = "grey",
         lty = 2)
  plot(
    d,
    x$weighted_agreement,
    type = "l",
    lwd = 2,
    col = "firebrick",
    xlab = "Depth",
    ylab = expression((X[d]/d) %.% p^d),
    main = sprintf(
      "Weighted contribution (RBO = %.4f)",
      x$rbo
    )
  )
  abline(v = min(x$length_L, x$length_S),
         col = "grey",
         lty = 2)
  invisible(x)
}



summary.rbo_ext <- function(object, ...) {
  
  X <- object$overlap_depth
  d <- seq_along(X)
  
  overlap_ratio <- X / d
  
  res <- list(
    rbo = object$rbo,
    p = object$p,
    length_L = object$length_L,
    length_S = object$length_S,
    overlap_final = object$overlap_final,
    max_possible_overlap = min(object$length_L,
                               object$length_S),
    overlap_ratio = overlap_ratio,
    overlap_depth = X
  )
  
  class(res) <- "summary.rbo_ext"
  
  res
}


print.summary.rbo_ext <- function(x,
                                  depths = c(1, 5, 10, 20, 50, 100),
                                  digits = 4,
                                  ...) {
  
  cat("Summary of Rank-Biased Overlap (Extrapolated)\n")
  cat("=============================================\n\n")
  
  cat(sprintf("RBO score          : %.*f\n", digits, x$rbo))
  cat(sprintf("Persistence (p)    : %.3f\n", x$p))
  cat(sprintf("Length L           : %d\n", x$length_L))
  cat(sprintf("Length S           : %d\n", x$length_S))
  cat(sprintf("Final overlap      : %d\n", x$overlap_final))
  cat(sprintf("Maximum possible   : %d\n\n",
              x$max_possible_overlap))
  
  depths <- depths[depths <= length(x$overlap_ratio)]
  
  if (length(depths) > 0) {
    
    tab <- data.frame(
      Depth = depths,
      Overlap = x$overlap_depth[depths],
      Agreement = round(x$overlap_ratio[depths], digits)
    )
    
    cat("Overlap by depth\n")
    print(tab, row.names = FALSE)
  }
  
  invisible(x)
}

print.rbo_ext <- function(x, digits = 4, ...) {
  
  cat("Rank-Biased Overlap (Extrapolated)\n")
  cat("---------------------------------\n")
  cat(sprintf("RBO           : %.*f\n", digits, x$rbo))
  cat(sprintf("p             : %.3f\n", x$p))
  cat(sprintf("Length L      : %d\n", x$length_L))
  cat(sprintf("Length S      : %d\n", x$length_S))
  cat(sprintf("Overlap       : %d\n", x$overlap_final))
  
  invisible(x)
}


#---------------------------------------------##
# Test ----
#---------------------------------------------##
l1 <- names(sort(list_grna$VEGFAs3$`GUIDE-seq_v1`,decreasing = T))
l2 <- (names(sort(list_grna$VEGFAs3$`GUIDE-seq_v2`,decreasing = T)))



res <- rbo_ext_fast2(l1, l2, p = 0.8)
print(res)

plot(res)

