auto_bi_fixture <- function() {
  data("AutoBI", package = "ChainLadder", envir = environment())
  ChainLadder::incr2cum(AutoBI$AutoBIPaid)
}

fit_age_fixture <- function() {
  set.seed(42)
  suppressWarnings(clmplus(AggregateDataPP(auto_bi_fixture()), "a"))
}
