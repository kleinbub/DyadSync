% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sync_AMICo2_RC.R
\name{peakCor}
\alias{peakCor}
\title{this similarity index rescales the v-p-v triangle according to the 
highest peak observed by that participant during the session
then calculates the root mean square of the absolute differences
ora x e y sono normalizzati su una scala da 0 a 1 dove 1 è il picco
più alto di quel soggetto.
Dimostrazione: max(x)*xmax == pv1$amp[i]
YET ACTUALLY
the maximum different situation is that of a parabola with apex (M,1)
compared to a flat zero line or parabola with apex (M,0)
so the maximum possible difference is the first parabola area
Area = 2/3 * 2M, where 2M equals the length of x}
\usage{
peakCor(x, y, xmax, ymax, diff)
}
\description{
In this new crazy formula, I calculate the maximum using the parabola
and normalize for that maximum to get back to a 0-1 value
this still works Best:
1 - the Root-mean-square deviation of the 2 peaks
max_p = 1
}
