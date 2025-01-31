#' Simulate data from bn_df object
#'
#' @param bn_df initialised bn_df object, with simulation instructions. Created with `bn_create`
#' @param known_df data.frame. Optional data.frame containing upstream variables used for simulation.
#' @param pop_size integer. The size of the dataset to be created.
#' @param keep_all logical. Keep all simulated variables or only keep those specified by `keep`
#' @param .id character. Name of id column placed at the start of the dataset. If NULL (default) then no id column is created.
#'
#' @return tbl
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
#'
#' # simulate the dataset
#' dummy_data <- bn_simulate(bn, pop_size = pop_n, keep_all = FALSE, .id = "patient_id")
bn_simulate <- function(bn_df, known_df=NULL, pop_size, keep_all=FALSE, .id=NULL){

  if (is.null(known_df)) {
    ..n <- pop_size
  } else {
    ..n <- nrow(known_df)
  }

  stopifnot(".id must be NULL or a length 1 character giving the name of the ID variable" = ( length(.id)==1 & is.character(.id) | is.null(.id)))

  dagitty <- bn2dagitty(bn_df)
  rbernoulli <- purrr::rbernoulli

  # creates the simulation function for each variable
  bn_df1 <-
    bn_df %>%
    dplyr::mutate(
      variable_expr = purrr::map(variable_formula, ~rlang::quo_squash(rlang::f_rhs(.))),
      # bn_fun = purrr::pmap(tibble::lst(variable, variable_formula), function(variable, variable_formula){
      #   function(tib){
      #     row_num <- seq_len(nrow(tib))
      #     x <- purrr::simplify(purrr::map(row_num,  ~eval(rlang::f_rhs(variable_formula), tib[.,])))
      #     tib1 <- tib
      #     tib1[variable] <- x
      #     tib1
      #   }
      # }),
    )

  #reorder based on dependencies so that simulation will create variables in the right order
  bn_ordered <-
    dagitty %>%
    dagitty::topologicalOrdering() %>%
    tibble::enframe(name="variable", value='topological_order') %>%
    tidyr::unnest(topological_order) %>%
    dplyr::left_join(bn_df1, ., by='variable') %>%
    dplyr::arrange(topological_order)


  # simulate complete dataset (with a patient ID variable in the initiated dataset)
  if (is.null(known_df)) {
    tbl0 <- tibble::tibble(.id = seq_len(pop_size))
  } else {
    tbl0 <- known_df
  }

  # variables to simulate
  bn_ordered_unknown <- bn_ordered %>% dplyr::filter(!known)

  named_expr <- rlang::set_names(bn_ordered_unknown$variable_expr, bn_ordered_unknown$variable)

  # simulate variables
  tblsim_complete <-
    tbl0 %>%
    dplyr::mutate(
      !!!named_expr
    )

  # create list of formulae that determining missingness
  missing_formula <- rlang::set_names(bn_ordered_unknown$missing_formula, bn_ordered_unknown$variable)

  # introduce NAs according to formulae in `missing_formula`
  tblsim_missing1 <- purrr::pmap_df(
    tibble::lst(variable = tblsim_complete[bn_ordered_unknown$variable], missing_formula, simdat=list(tblsim_complete)),
    function(variable, missing_formula, simdat){

      mask <- eval(rlang::f_rhs(missing_formula), simdat)
      if(!inherits(variable, "factor")){
        NA_type_ <- NA
        mode(NA_type_) <- typeof(variable)
        dplyr::if_else(!mask, variable, NA_type_)
      } else {
        dplyr::if_else(!mask, variable, factor(NA))
      }
    }
  )


  # create list of formulae that determining missingness
  needs <- rlang::set_names(bn_ordered_unknown$needs, bn_ordered_unknown$variable)
  # add all upstream needs
  for (i in seq_along(needs)) {
    if (length(needs[[i]]) > 0) {
      for(j in seq_along(needs[[i]])) {
        needs[[i]] <- c(needs[[needs[[i]][j]]], needs[[i]])
      }
      needs[[i]] <- unique(needs[[i]])
    }
  }

  # introduce NAs according to formulae in `missing_formula`
  tblsim_missing2 <- purrr::pmap_df(
    tibble::lst(variable = tblsim_missing1[bn_ordered_unknown$variable], needs, simdat=list(tblsim_missing1)),
    function(variable, needs, simdat){


      if(length(needs)!=0){
        mask <- simdat %>%
          dplyr::select(tidyselect::all_of(needs)) %>%
          dplyr::mutate(dplyr::across(tidyselect::all_of(needs), ~!is.na(.))) %>%
          dplyr::rowwise() %>%
          dplyr::mutate(
            need_satisfied=!all(dplyr::c_across(tidyselect::all_of(needs)))
          ) %>% purrr::pluck("need_satisfied")
      }
      else{
        mask <- rep(FALSE, nrow(simdat))
      }

      if(!inherits(variable, "factor")){
        NA_type_ <- NA
        mode(NA_type_) <- typeof(variable)
        dplyr::if_else(!mask, variable, NA_type_)
      } else {
        dplyr::if_else(!mask, variable, factor(NA))
      }
    }
  )

  # combine known and simulated variables
  tblsim <- dplyr::bind_cols(tbl0, tblsim_missing2)

  # choose which variables to return
  returnvars <-
    bn_df1 %>%
    dplyr::filter(keep | keep_all, known==FALSE) %>%
    purrr::pluck("variable")
  tblout <-
    tblsim %>%
    dplyr::select(
      names(tbl0),
      tidyselect::all_of(returnvars)
    )

  # rename id variable
  if(is.null(.id)){
    tblout$.id <- NULL
  } else {
    names(tblout)[names(tblout) == '.id'] <- .id
  }

  tblout
}
