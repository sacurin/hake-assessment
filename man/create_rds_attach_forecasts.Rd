% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create-rds-attach-forecasts.R
\name{create_rds_attach_forecasts}
\alias{create_rds_attach_forecasts}
\title{Load catch levels and forecasting and attach to as RDS file}
\usage{
create_rds_attach_forecasts(model_path = NULL, verbose = TRUE, ...)
}
\arguments{
\item{model_path}{Directory name of model to be loaded}

\item{verbose}{Logical. If \code{TRUE}, write more output to the console}

\item{...}{Passed to \code{\link[=check_catch_levels]{check_catch_levels()}}}
}
\value{
Nothing
}
\description{
Load catch levels and forecasting and attach to as RDS file
}
\details{
Load forecasting and attach it to an already-created RDS file as the list
element \code{forecasts} and overwrite the RDS file
}
