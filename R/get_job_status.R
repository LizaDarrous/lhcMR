#' Title
#'
#' @param slr_job
#'
#' @return
# @export
#'
#' @examples
get_job_status <- function (slr_job)
{
  if (!(class(slr_job) == "slurm_job"))
    stop("input must be a slurm_job")
  stat <- suppressWarnings(system(paste("squeue -n", slr_job$jobname),
                                  intern = TRUE))
  if (length(stat) > 1) {
    res = "Job running or in queue."
  }
  else {
    res = "Job completed or stopped."
    # tmpdir <- paste0("rslurm", slr_job$jobname)
    #  out_files <- file.path(tmpdir, paste0("slurm_", 0:(slr_job$nodes -
    #                                                      1), ".out"))
    #for (outf in out_files) {
    #  cat(paste("\n----", outf, "----\n\n"))
    #  cat(paste(readLines(outf), collapse = "\n"))
  }
  return(res)
}
