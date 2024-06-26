% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load-extra-mcmc-init-nage.R
\name{load_extra_mcmc_init_nage}
\alias{load_extra_mcmc_init_nage}
\title{Load the inital numbers-at-age as defined in Report_mce_xxxx.sso files}
\usage{
load_extra_mcmc_init_nage(reps, verbose = TRUE, ...)
}
\arguments{
\item{reps}{A list of vectors representing one element per report file}

\item{verbose}{Logical. Show messages}

\item{...}{Absorbs arguments meant for other functions}

\item{progress_n}{Report every time this many list items are processed.
Consider how many posteriors there are, this should be a fairly large
proportion of that (around 1/8th) or there will be too much output and
it will run really slow}
}
\value{
A list of outputs
}
\description{
Load the inital numbers-at-age as defined in Report_mce_xxxx.sso files
}
