% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/categories_extraction.R
\name{extractEpochs}
\alias{extractEpochs}
\title{extract epochs in a simple list}
\usage{
extractEpochs(experiment, signal, sync, series, diary, category, by, ...)
}
\arguments{
\item{experiment}{}

\item{signal}{}

\item{sync}{}

\item{series}{}

\item{diary}{the Diary}

\item{category}{the Category of the Diary to be extracted}

\item{by}{currently unused. In future will be used to split the data by experimental group, participant, or any other relevant condition}

\item{...}{currently unused.}
}
\description{
This function extracts the selected epochs from every session of a "DyadExperiment" object and puts them in
 a simple list. Categories must be created with epochSeries() beforehand.
}
