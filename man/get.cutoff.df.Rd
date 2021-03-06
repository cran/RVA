% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_cutoff.R
\name{get.cutoff.df}
\alias{get.cutoff.df}
\title{Create ggplot object for number of differntially expressed genes with
different FDR and fold change cutoff.}
\usage{
get.cutoff.df(datin, pvalues, FCs, FCflag = "logFC", FDRflag = "adj.P.Val")
}
\arguments{
\item{datin}{Dataframe from plot_cutoff_single.}

\item{pvalues}{A set of p-values for FDR cutoff to be checked.}

\item{FCs}{A set of fold change cutoff to be checked.}

\item{FCflag}{The column name of the log2FC in the summary statistics table.}

\item{FDRflag}{The column name of the False Discovery Rate (FDR) in the summary statistics table.}
}
\description{
This function processes dataframe from plot_cutoff_single function
and produces a ggplot object which depicts the number of differntially expressed
genes with different FDR and fold change cutoff.
}
