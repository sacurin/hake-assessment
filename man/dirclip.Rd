% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-dirclip.R
\name{dirclip}
\alias{dirclip}
\title{Copy the current working directory to the system clipboard}
\usage{
dirclip()
}
\value{
\code{TRUE} or \code{FALSE}. If \code{TRUE}, a character string representing the
current working directory is now in the system clipboard. If \code{FALSE}, the
writing to clipboard was not possible
}
\description{
Copy the current working directory to the system clipboard
}
\details{
An error will be thrown if writing to the clipboard fails
}
