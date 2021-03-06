% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validate.input.R
\name{validate.pvalflag}
\alias{validate.pvalflag}
\title{Validate pval flag}
\usage{
validate.pvalflag(data, value)
}
\arguments{
\item{data}{A list of summary statistics table (data.frame) from limma or DEseq2, where rownames are gene id.}

\item{value}{P value flag.}
}
\description{
To ensure p value flags are the same accross datasets.
}
\references{
Xingpeng Li, Tatiana Gelaf Romer & Siddhartha Pachhai RVA - RNAseq Visualization Automation tool.
}
