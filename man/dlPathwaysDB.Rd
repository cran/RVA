% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_pathway.R
\name{dlPathwaysDB}
\alias{dlPathwaysDB}
\title{DL Pathways DB}
\usage{
dlPathwaysDB(pathway.db, customized.pathways = NULL, ...)
}
\arguments{
\item{pathway.db}{The databse to be used for encrichment analysis. Can be one of the following, "rWikiPathways", "KEGG", "REACTOME", "Hallmark","rWikiPathways_aug_2020"}

\item{customized.pathways}{the user provided pathway added for analysis.}

\item{...}{pass over parameters}
}
\value{
Returns a dataframe.
}
\description{
Download gene database for enrichment.
}
\references{
Xingpeng Li & Siddhartha Pachhai RVA - RNAseq Visualization Automation tool.
}
