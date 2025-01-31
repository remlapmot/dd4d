#' Creates a bayesian network object from a list of nodes
#'
#' Converts list to data frame which is a bit easier to work with, and embellishes with some useful columns.
#' The function performs a few checks on the list, for instance to make sure the graph is acyclic
#' and that variables used in the expressions are defined elsewhere or already known.

#' The known_variables argument is for passing a character vector of variables names
#' for variables that are already defined externally in a
#' given dataset, which can be passed to bn_simulate

#' whilst variable_formula is the variable name itself, this is to help with the bn_simulate function
#' it doesn't actually lead to self-dependence (eg var depends on var)

#'
#' @param list of node objects, created by `bn_node`.
#' @param known_variables character vector of variables that will be provided by an external dataset
#'
#' @return data.frame
#' @export
#' @examples
#' library(dplyr)
#' set.seed(314159)
#' pop_n <- 1000
#' index_date <- as.Date("2020-12-08")
#' pfizer_name <- paste("COVID-19 mRNA Vaccine Comirnaty 30micrograms/0.3ml dose conc",
#'   "for susp for inj MDV (Pfizer)")
#' az_name <- "COVID-19 Vaccine Vaxzevria 0.5ml inj multidose vials (AstraZeneca)"
#' # define each variable as a node in a bayesian network
#' sim_list <- lst(
#'
#'   registered = bn_node(
#'     ~ rbernoulli(n = ..n, p = 0.99)
#'   ),
#'   age = bn_node(
#'     ~ as.integer(rnorm(..n, mean= 55, sd=20)),
#'   ),
#'   sex = bn_node(
#'     ~ rfactor(n = ..n, levels = c("female", "male", "intersex", "unknown"),
#'               p = c(0.50, 0.49, 0, 0.01)),
#'     missing_rate = ~0.01
#'   ),
#'   diabetes = bn_node(
#'     ~ rbernoulli(n = ..n, p = 0.10)
#'   ),
#'   vaccine_date1 = bn_node(
#'     ~ runif(n = ..n, 0, 150),
#'     missing_rate = ~0.2
#'   ),
#'   vaccine_date2 = bn_node(
#'     ~ vaccine_date1 + runif(n = ..n, 3*7, 16*7),
#'     missing_rate = ~0.1,
#'     needs = "vaccine_date1"
#'   ),
#'   vaccine_product1 = bn_node(
#'     ~ rcat(n = ..n, c(pfizer_name, az_name), c(0.5, 0.5)),
#'     needs = "vaccine_date1"
#'   ),
#'   vaccine_product2 = bn_node(
#'     ~ if_else(runif(..n) < 0.95, vaccine_product1, c(az_name)),
#'     needs = "vaccine_date2"
#'   ),
#'   death_date = bn_node(
#'     ~ runif(n = ..n,  0, 1000),
#'     missing_rate = ~0.95
#'   )
#' )
#'
#' bn <- bn_create(sim_list,known_variables = c("pfizer_name","az_name"))
bn_create <- function(list, known_variables=NULL){


  stopifnot("'list' must be a list where each element is an object of class 'node'" = all(sapply(list, function(x){"node" %in% class(x)})))

  df <- tibble::enframe(list, name="variable", value="list")
  # df <- tidyr::unnest_wider(df, "list") # this no longer works after update to tidyr, even with ptype specification, so use tidyr::hoist instead
  df <- tidyr::hoist(df, .col="list", variable_formula=1, missing_rate=2, keep=3, needs=4)

  # this bit is needed because if there are no nodes with "needs" specified, this variable does not
  # get unnested, so needs to be created explicitly
  # otherwise if at least one "need" is specified, then for all nodes without needs, `need` is converted from
  # character() to NULL, so need to undo.
  if(!("needs" %in% names(df))){
    df$needs = list(character())
  } else{
    df$needs = purrr::map(df$needs, ~{
      if(is.null(.)){
        character()
      } else
      if(length(.)==1 & all(is.na(.))){
        character()
      } else
        .
      })
  }

  df <- dplyr::mutate(df,
    parents = purrr::map(variable_formula, function(x) {parents <- all.vars(x); parents[parents!="..n"] }),
    missing_formula = purrr::map(missing_rate, ~{
      rhs <- deparse1(rlang::f_rhs(.))
      fun <- stats::as.formula(paste0("~rbernoulli(n=..n, p=", rhs, ")"))
      fun
    }),
    missing_parents = purrr::map(missing_formula, function(x) {parents <- all.vars(x); parents[parents!="..n"] }),
    known=FALSE
  )

  if(!is.null(known_variables)){

    df_available <-
      tibble::tibble(
        variable = known_variables,
        variable_formula = list(~stats::as.formula(paste0("~", variable))),
        missing_rate = list(~0),
        keep = TRUE,
        parents = list(character()),
        missing_formula = list(stats::as.formula("~rbernoulli(n=..n, p=0)")),
        missing_parents = list(character()),
        known = TRUE,
        needs = list(character())
      )

    df <- dplyr::bind_rows(df_available, df)
  }

  df$in_order <- seq_len(nrow(df))

  dagitty <- bn2dagitty(df)


  parents_check <- purrr::map(df$variable, ~ dagitty::parents(dagitty, .)) # should be same as 'parents' above
  stopifnot("mismatch between dependencies and parents" = all(purrr::map2_lgl(df$parents, parents_check, ~all(.x %in% .y) & all(.x %in% .y))))

  df$children <- purrr::map(df$variable, ~ dagitty::children(dagitty, .))

  stopifnot("graph is not acyclic" = dagitty::isAcyclic(dagitty))

  if(!all(purrr::simplify(unique(rlang::flatten(df$parents))) %in% df$variable)){
    print(
      purrr::simplify(unique(rlang::flatten(df$parents)))[!(purrr::simplify(unique(rlang::flatten(df$parents))) %in% df$variable)]
    )
    stop("not all dependencies are defined")
  }

  stopifnot("variable names are not unique" = length(df$variable) == length(unique(df$variable)))

  df
}
