% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signals_PPG.R
\name{PPG2beats}
\alias{PPG2beats}
\title{Extract beats from PPG signal}
\usage{
PPG2beats(
  x,
  r = 0.65,
  f = 0.75,
  buffer = 8,
  FIR_args = list(type = "low", cut = 3, NAsub = 0),
  correctionRangeSeconds = 1/1.5,
  minPeakAmplitude = 0.01,
  plotName = NULL,
  ...
)
}
\arguments{
\item{x}{a rats time-series representing the original PPG signal}

\item{r}{the rejection factor, expressed as a weight of the buffer's median R-R interval}

\item{f}{the fill factor, expressed as a weight of the buffer's median R-R interval}

\item{buffer}{integer. the rolling buffer, expressed in number of beats}

\item{FIR_args}{a named list of arguments to be passed to \link[DyadSignal]{FIR}) filtering function.
'cut' and 'type' must be present, representing . If empty, the filtering is disabled.}

\item{correctionRangeSeconds}{the half range in which the real maximum/minimum
value should be searched, in seconds.
around the derivative shift. Should be less then the periodicity of the
signal.  0.5 is good for skin conductance.}

\item{minPeakAmplitude}{the minimum delta from valley to peak for detection.
Skin conductance traditionally uses 0.05uS, or 0.03uS}

\item{plotName}{an optional filename to save a svg representation of the beats extraction}
}
\value{
a list of:
\describe{
  \item{bool}{a logical vector of the same length of x with TRUE
   value corresponding to a match
   }
  \item{samples}{the position of beats relatively to x's index}
  \item{x}{temporal coordinates of the beats (if x has a frequency attribute)}
  \item{y}{The value of x correspoding to the detected beats}
  \item{type}{a character vector defining for each beats if its a detected peak
              ('p') or an interpolated one ('p_filled')}
  \item{amp}{the through-to-peak amplitudes of the detected beats}
  \item{index}{The absolute position of the beats along the signal. The values
     are independent from the 'mode' argument, allowing interaction between
     different calls of the function.}
  \item{reject_samples}{the position of rejected beats relatively to x's index}
  
  
}
}
\description{
This functions employs a low pass filter to attenuate non-systolic peaks,
then rejects peaks that are too fast compared to the median rate of a buffer.
Finally it interpolates too large gaps.
}
