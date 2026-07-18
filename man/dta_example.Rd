\name{dta_example}
\alias{dta_example}
\docType{data}
\title{Simulated Diagnostic Test Accuracy Studies}
\description{
A simulated dataset of 18 diagnostic-accuracy studies in three subgroups. It
contains study-level two-by-two table counts and identifiers suitable for
demonstrating the main \pkg{DTAtoolkit} workflows.
}
\usage{dta_example}
\format{
A data frame with 18 rows and 9 variables:
\describe{
  \item{study}{Display name for the study.}
  \item{study_id}{Stable study identifier.}
  \item{unique_row_id}{Unique identifier for each dataset row.}
  \item{subgroup}{Simulated test-platform subgroup.}
  \item{setting}{Simulated clinical setting.}
  \item{tp}{True-positive count.}
  \item{fn}{False-negative count.}
  \item{tn}{True-negative count.}
  \item{fp}{False-positive count.}
}
}
\details{
All values are simulated and do not represent real studies, patients, tests,
or clinical outcomes.
}
\examples{
data(dta_example)
head(dta_example)

\dontrun{
dta_forest(
  dta_example,
  subgroup = subgroup,
  study = study,
  tp = tp, tn = tn, fp = fp, fn = fn,
  type = "combined",
  pairwise = "never"
)
}
}
\keyword{datasets}
