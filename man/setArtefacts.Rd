% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generic_filters.R
\name{setArtefacts}
\alias{setArtefacts}
\title{setArtefacts
setArtefacts legge una tabella o lista con le componenti chiamate esattamente: 'dyad', 'session', 'start' ed 'end'
salva le informazioni corrispondenti in una tabella $artefacts nell'oggetto DyadSignal e mette a FALSE le epoche
corrispondenti nello stream 'valid' dell'oggetto.
'end' può anche essere una stringa che contiene le parole 'fine' o 'end'}
\usage{
setArtefacts(x, startEnd, signal)
}
\arguments{
\item{x}{a Dyad... object}

\item{startEnd}{a data.frame (or list) with the following components: 'dyad', 'session', 'start', 'end'.}

\item{signal}{string specifying the name of the signal}
}
\description{
setArtefacts
setArtefacts legge una tabella o lista con le componenti chiamate esattamente: 'dyad', 'session', 'start' ed 'end'
salva le informazioni corrispondenti in una tabella $artefacts nell'oggetto DyadSignal e mette a FALSE le epoche
corrispondenti nello stream 'valid' dell'oggetto.
'end' può anche essere una stringa che contiene le parole 'fine' o 'end'
}
