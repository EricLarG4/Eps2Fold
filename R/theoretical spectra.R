# Dependencies----
library('data.table')
library('stringr')

#' Generate theoretical spectra for nucleotide sequences
#'
#' @param sequences Character vector of nucleotide sequences.
#' @param wl_range Integer vector of wavelengths (nm) over which spectra are computed.
#'
#' @return data.table with columns `seq`, `wl`, and `eps`.
#' @export
#'
#' @examples
#' # generate_spectra(c("ATCG", "GCTA"), wl_range = 220:310)
generate_spectra <- function(sequences, wl_range = 215:310) {
  # sanity-check inputs
  stopifnot(length(sequences) > 0)

  # assign ID to each sequence to retain ordering after expansions
  seq_dt <- data.table(seq = sequences)[, seq_id := .I]

  # for every sequence, extract terminal bases and adjacent dinucleotide motifs
  seq_pairs <- seq_dt[
    ,
    {
      n <- str_length(seq)
      if (n == 0) {
        stop("Empty sequence provided.")
      }

      pair_positions <- if (n > 1) paste0(1:(n - 1), "-", 2:n) else character()
      pair_motifs <- if (n > 1) str_sub(seq, 1:(n - 1), 2:n) else character()

      data.table(
        position = c(1, n, pair_positions),
        nn = c(
          str_sub(seq, 1, 1),      # first base contribution
          str_sub(seq, n, n),      # last base contribution
          pair_motifs              # internal dinucleotide contributions
        )
      )
    },
    by = .(seq_id, seq)
  ]

  # attach epsilon contributions (scaled) for each motif
  seq_pairs[nn.260, on = "nn", contrib := epsij * 1000]

  # ensure every motif found a match in the epsilon lookup table
  if (seq_pairs[is.na(contrib), .N] > 0) {
    stop("Some motifs are missing from the epsilon database.")
  }

  # restrict epsilon ratios to requested wavelength range
  nn_param_subset <- nn.param[wl %in% wl_range]

  # join spectra ratios to motifs, sum contributions per sequence/wavelength
  nn_param_subset[
    seq_pairs,
    on = "nn",
    allow.cartesian = TRUE
  ][,
    .(eps = sum(ratio * contrib)),
    by = .(seq, wl)
  ][order(seq, wl)]
}
