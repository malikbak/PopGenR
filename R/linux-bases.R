#' Create Fixed-Width BED Windows with bedtools
#'
#' Create genomic windows from a FASTA index or genome file by calling
#' `bedtools makewindows`.
#'
#' @param Ffile Path to a genome file accepted by `bedtools makewindows -g`.
#' @param Window Window size in base pairs.
#' @param output Output BED file path. Defaults to `test.bed`.
#' @param bedtools Path or command name for bedtools.
#'
#' @return Invisibly returns `output` after bedtools completes.
#'
#' @examples
#' \dontrun{
#' Mbed("genome.fai", Window = 40000, output = "windows.bed")
#' }
#'
#' @export
Mbed <- function(Ffile, Window, output = "test.bed", bedtools = "bedtools") {
  if (!file.exists(Ffile)) {
    stop("`Ffile` does not exist: ", Ffile, call. = FALSE)
  }

  status <- system2(
    bedtools,
    c("makewindows", "-g", Ffile, "-w", Window),
    stdout = output
  )

  if (!identical(status, 0L)) {
    stop("bedtools makewindows failed with status ", status, ".", call. = FALSE)
  }

  invisible(output)
}
