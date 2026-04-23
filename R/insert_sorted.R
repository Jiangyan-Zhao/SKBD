#' Insert a new dose into the working grid and keep all trial vectors aligned
#'
#' @description
#' `insert_sorted()` inserts a new dose value (`newdose`) into the current working
#' dose grid (`dose_set_work`), then **sorts** the updated grid and **reorders**
#' all associated vectors accordingly:
#' \itemize{
#'   \item `tox_prob_work`: (simulation-only) true toxicity probabilities aligned with doses
#'   \item `n_dlt_work`: number of DLTs observed at each dose
#'   \item `n_treated_work`: number treated at each dose
#'   \item `is_eliminated`: dose elimination indicators
#' }
#' It returns the updated vectors and the index of the inserted dose in the
#' sorted grid (`d_idx`).
#'
#' @details
#' This helper is designed for adaptive dose insertion routines (e.g., ADM/SKBD
#' insertion) where the dose grid can expand during the trial. By centralizing
#' the "insert + sort + reorder" operation, it prevents common bugs such as
#' mismatched vector lengths, inconsistent ordering across vectors, or incorrect
#' dose indexing after insertion.
#'
#' **Important:** This function does not check whether `newdose` already exists
#' in `dose_set_work`. Callers should guard against duplicates (e.g., using a
#' tolerance check) before calling `insert_sorted()`.
#'
#' @param dose_set_work Numeric vector. Current working dose grid (typically standardized to `[0, 1]`).
#'   Must be strictly increasing and aligned with all other vectors.
#' @param tox_prob_work Numeric vector. Current true toxicity probabilities aligned with `dose_set_work`
#'   (simulation-only). Must have the same length as `dose_set_work`.
#' @param n_dlt_work Integer/numeric vector. Current number of observed DLTs aligned with `dose_set_work`.
#'   Must have the same length as `dose_set_work`.
#' @param n_treated_work Integer/numeric vector. Current number treated aligned with `dose_set_work`.
#'   Must have the same length as `dose_set_work`.
#' @param is_eliminated Logical vector. Dose elimination indicators aligned with `dose_set_work`.
#'   Must have the same length as `dose_set_work`.
#' @param newdose Numeric scalar. The newly inserted dose value.
#' @param tox_new Numeric scalar. The (simulation-only) true toxicity probability at `newdose`.
#'
#' @return A list with components:
#' \describe{
#'   \item{dose_set_work}{Sorted working dose grid after insertion.}
#'   \item{tox_prob_work}{Reordered toxicity probabilities aligned with the new grid.}
#'   \item{n_dlt_work}{Reordered DLT counts aligned with the new grid (0 at the inserted dose).}
#'   \item{n_treated_work}{Reordered treated counts aligned with the new grid (0 at the inserted dose).}
#'   \item{is_eliminated}{Reordered elimination indicators aligned with the new grid (FALSE at the inserted dose).}
#'   \item{d_idx}{Integer index of `newdose` in the sorted `dose_set_work`.}
#' }
#'
#' @examples
#' dose_set_work = c(0.25, 0.50, 0.75)
#' tox_prob_work = c(0.10, 0.25, 0.40)
#' n_dlt_work = c(0, 1, 2)
#' n_treated_work = c(3, 6, 6)
#' is_eliminated = c(FALSE, FALSE, FALSE)
#'
#' newdose = 0.375
#' tox_new = 0.20
#'
#' out = insert_sorted(
#'   dose_set_work, tox_prob_work, n_dlt_work, n_treated_work,
#'   is_eliminated, newdose, tox_new
#' )
#' out$dose_set_work
#' out$d_idx
#'
#' @noRd
insert_sorted = function(
    dose_set_work, tox_prob_work, n_dlt_work, n_treated_work,
    is_eliminated, newdose, tox_new
) {
  # Combine the existing grid with the new dose, then sort and keep the permutation
  tmp = sort(c(dose_set_work, newdose), index.return = TRUE)
  dose_new = tmp$x
  ix = tmp$ix
  
  # Reorder all aligned vectors using the same permutation.
  # The inserted dose receives default values:
  #   n_dlt = 0, n_treated = 0, is_eliminated = FALSE, tox_prob = tox_new.
  out = list(
    dose_set_work = dose_new,
    tox_prob_work = c(tox_prob_work, tox_new)[ix],
    n_dlt_work = c(n_dlt_work, 0)[ix],
    n_treated_work = c(n_treated_work, 0)[ix],
    is_eliminated = c(is_eliminated, FALSE)[ix]
  )
  
  # Locate the inserted dose index in the new sorted grid.
  # Use nearest match to be robust to floating-point representations.
  out$d_idx = which.min(abs(out$dose_set_work - newdose))
  
  return(out)
}
