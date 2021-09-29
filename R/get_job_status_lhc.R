#' Title
#'
#' @param slr_job
#'
#' @return
#' @keywords internal
#' NOT EXPORTED @export
#'
#' @examples
get_job_status_lhc <- function (slr_job)
{
  if (!(class(slr_job) == "slurm_job"))
    stop("input must be a slurm_job")
  stat <- suppressWarnings(system(paste("squeue -n", slr_job$jobname),
                                  intern = TRUE))
  if (length(stat) > 1) {
    res = "Job running or in queue."
    completed = FALSE
  }
  else {
    res = "Job completed or stopped."
    completed = TRUE
  }
  return(completed)
}
