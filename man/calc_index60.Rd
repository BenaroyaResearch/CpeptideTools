% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_index60.R
\name{calc_index60}
\alias{calc_index60}
\title{Calculate the T1D Diagnostic Index60 value (Index60), a composite measure of pre-diagnosis T1D status}
\usage{
calc_index60(
  input_data,
  cpeptide_fasting_colname = "pep0",
  cpeptide_60min_colname = "pep60",
  glucose_60min_colname = "glu60")
}
\arguments{
\item{input_data}{data frame containing the C-peptide and glucose values}

\item{cpeptide_fasting_colname}{character or numeric, the column containing the fasting C-peptide values in ng/mL. Defaults to "cpeptide_fasting".}

\item{cpeptide_60min_colname}{character or numeric, the column containing the OGTT 60-minute C-peptide values in ng/mL. Defaults to "cpeptide_60min".}

\item{glucose_60min_colname}{character or numeric, the column containing the OGTT 60-minute glucose values in mg/dL. Defaults to "glucose_60min".}
}
\value{
a vector containing the Index60 values
}
\description{
This function calculates the T1D Diagnostic Index60 value (Index60) of
Sosenko et al. 2015 (PMID 25519451), a composite measure of pre-diagnosis
T1D status from an oral glucose tolerance test (OGTT) that has been found
to predict onset of T1D. Published work indicates that Index60 values above
2.00 or 1.00 may predict onset of clinical T1D. The function takes an input
data frame with fasting C-peptide levels, 60-minute C-peptide levels from
OGTT, and 60-minute glucose levels from OGTT.
NOTE: Values for glucose and C-peptide must be in the standardized units used
in Index60, specifically ng/mL for C-peptide and mg/dL for glucose.
}
\details{
The equation used is:
"Index60 = 0.3695 * log (fasting C-peptide) + 0.0165 * (60-minute glucose) - 0.3644 * (60-minute C-peptide)"
}
