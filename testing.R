## function to replace itself with an expression of... itself
test_me <- function(x, y, ...) {
  rlang::expr(!!match.call())
}

`+` <- function(x, y) {
  rlang::expr(!!match.call())
}

`+` <- function(x, y) {
  rlang::expr(!!sys.call())
}

`[` <- function(...) {
  rlang::expr(!!sys.call())
}