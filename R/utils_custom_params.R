getCustomExpressions <- function(parTable) {
  custom <- parTable$op == ":="

  if (!any(custom))
    return(list())

  idx <- which(custom)
  names <- getParNamesFromParTable(parTable)
  expressions <- emptyNamedList(names[idx])

  for (i in idx) {
    nm <- names[[i]]
    expr <- parse(text = parTable[i, "rhs"])
    pls_stopif(length(expr) != 1L,
      sprintf("Custom parameter (%s) must contain exactly one expression.", nm)
    )
    validateCustomExpression(expr[[1L]], nm)
    expressions[[nm]] <- expr
  }

  expressions
}


evalCustomExpressions <- function(pars, labels, expressions) {
  keep <- intersect(names(pars), names(labels))
  envir <- customExpressionEnv(stats::setNames(pars[keep], nm = labels[keep]))

  names <- names(expressions)
  out <- stats::setNames(
    rep(NA_real_, length(expressions)),
    nm = names
  )

  for (i in seq_along(expressions)) {
    expr <- expressions[[i]]

    value <- tryCatch(eval(expr, envir = envir), error = function(e) {
      pls_msg_warn(
        sprintf("Failed to evaluate expression for parameter (%s)", names[[i]]),
        "Message:", conditionMessage(e)
      )
      NA_real_
    })

    nm <- names[[i]]
    label <- fillna(labels[nm], nm)
    assign(unname(label), value, envir = envir)
    out[[i]] <- value
  }

  out
}


CUSTOM_EXPRESSION_FUNCTIONS <- c(
  "+", "-", "*", "/", "^", "(", "sqrt", "abs", "log", "log10", "log2",
  "exp", "sin", "cos", "tan", "asin", "acos", "atan", "sinh", "cosh",
  "tanh", "min", "max", "sum", "prod", "mean"
)


validateCustomExpression <- function(expr, name) {
  if (!is.call(expr))
    return(invisible(TRUE))

  fun <- expr[[1L]]
  pls_stopif(!is.symbol(fun) || !as.character(fun) %in% CUSTOM_EXPRESSION_FUNCTIONS,
    sprintf("Unsupported expression in custom parameter (%s).", name),
    "Custom parameters only support arithmetic over parameter labels and",
    "the supported math functions:",
    paste0(CUSTOM_EXPRESSION_FUNCTIONS, collapse = ", ")
  )

  for (arg in as.list(expr)[-1L])
    validateCustomExpression(arg, name)

  invisible(TRUE)
}


customExpressionEnv <- function(values) {
  envir <- new.env(parent = emptyenv())

  for (fn in CUSTOM_EXPRESSION_FUNCTIONS)
    assign(fn, get(fn, envir = baseenv()), envir = envir)

  values <- values[nzchar(names(values))]
  list2env(as.list(values), envir = envir)
}
