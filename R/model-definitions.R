#' @importFrom forecast Arima forecast
NULL

supported_models <- list(
  a = StMoMo::StMoMo("log", TRUE, NULL, NULL),
  p = StMoMo::StMoMo("log", FALSE, "1", NULL),
  c = StMoMo::StMoMo("log", FALSE, NULL, "1"),
  ac = StMoMo::StMoMo("log", TRUE, NULL, "1"),
  ap = StMoMo::StMoMo("log", TRUE, "1", NULL),
  pc = StMoMo::StMoMo("log", FALSE, "1", "1"),
  apc = StMoMo::apc(),
  cbd = StMoMo::cbd(link = "log"),
  m6 = StMoMo::m6(link = "log"),
  m7 = StMoMo::m7(link = "log")
)

model_components <- list(
  a = c(age = TRUE, period = FALSE, cohort = FALSE),
  p = c(age = FALSE, period = TRUE, cohort = FALSE),
  c = c(age = FALSE, period = FALSE, cohort = TRUE),
  ac = c(age = TRUE, period = FALSE, cohort = TRUE),
  ap = c(age = TRUE, period = TRUE, cohort = FALSE),
  pc = c(age = FALSE, period = TRUE, cohort = TRUE),
  apc = c(age = TRUE, period = TRUE, cohort = TRUE),
  cbd = c(age = TRUE, period = TRUE, cohort = FALSE),
  m6 = c(age = TRUE, period = TRUE, cohort = TRUE),
  m7 = c(age = TRUE, period = TRUE, cohort = TRUE)
)
