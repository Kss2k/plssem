estimatePLS_Step0_5 <- function(model) {
  max.iter.0_5 <- model$status$max.iter.0_5

  model <- estimatePLS_Step0(model)

  for (i in seq_len(max.iter.0_5)) {
    model <- model |>
      estimatePLS_Step1() |>
      estimatePLS_Step2() |>
      estimatePLS_Step3() |>
      estimatePLS_Step4() |>
      estimatePLS_Step5()

    if (model$status$convergence) {
      break

    } else if (i >= max.iter.0_5) {
      warning("Convergence reached. Stopping.")
      break
    }
  }

  model$status$iterations <- model$status$iterations + i

  model
}


