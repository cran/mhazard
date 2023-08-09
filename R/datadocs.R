#' Aplastic anemia clinical trial data.
#'
#' Failure time data for 64 patients having severe aplastic anemia. Patients
#' were randomly assigned to receive either cyclosporine plus methotrexate
#' (CSPMTX) or methotrexate alone (MTX) with time from assignment to the
#' diagnosis of stage 2 or greater acute graft versus host disease (A-GVHD)
#' as a key outcome.
#'
#' @usage data(anemia)
#' @format A data frame with 64 rows and 5 variables:
#' \describe{
#'  \item{Treatment:}{Treatment group, either cyclosporine and methotrexate
#'  (CSPMTX) or methotrexate only (MTX).}
#'  \item{Time:}{Time in days to severe (stage >=2) acute graft versus host
#'  disease, death, or last contact.}
#'  \item{LAF:}{Patient assigned to laminar air flow isolation room at
#'  baseline. LAF = 1 signifies yes, LAF = 0 signifies no.}
#'  \item{Age:}{Age of patient in years.}
#'  \item{Censored:}{Boolean for whether or not time to severe A-GVHD is
#'  right-censored; i.e., the patient died without severe A-GVHD or was
#'  without severe A-GVHD at last contact.}
#' }
#'
#' @docType data
#' @source Table 1.1, page 14 of Prentice and Zhao 2019, taken from
#' Table 1.2 on page 3 in The Statistical Analysis of Failure Time Data,
#' Second Edition, by J. D. Kalbfleisch and R. L. Prentice. New York: Wiley
#' and Sons, 2002.
#' @references
#' Prentice, R., Zhao, S. "The statistical analysis of multivariate
#' failure time data: A marginal modeling approach", CRC Press (2019).
"anemia"

#' Bladder tumor recurrence data.
#'
#' Data from a randomized trial of patient recurrences of superficial
#' bladder tumor recurrence as conducted by the Veterans Administration
#' Cooperative Urological Group. These data were used to compare the
#' frequency of recurrences among 48 patients assigned to placebo, among
#' whom there were a total of 87 post-randomization recurrences, and 38
#' patients assigned to treatment with a drug called thiotepa, among whom
#' there were 45 recurrences during the follow-up period, which averaged
#' about 31 months.
#'
#' @usage data(bladder)
#' @format A data frame with 86 rows and 13 variables:
#' \describe{
#'  \item{Group:}{Placebo vs treatment (thiotepa).}
#'  \item{Number:}{Number of tumors present at baseline, with 8 denoting 8
#'  or more.}
#'  \item{Size:}{Size of largest tumor at baseline, in centimeters.}
#'  \item{CensorTime:}{Censor time, in months.}
#'  \item{T1:}{For i = 1,...,9, T_i is the recurrence time, in months.
#'  If there is no recurrence, T_1, ..., T_9 will be NA. If there is only one,
#'  T_1 will be filled and T_2, ..., T_9 will be NA, and so on.}
#'  \item{T2:}{See T1.}
#'  \item{T3:}{See T1.}
#'  \item{T4:}{See T1.}
#'  \item{T5:}{See T1.}
#'  \item{T6:}{See T1.}
#'  \item{T7:}{See T1.}
#'  \item{T8:}{See T1.}
#'  \item{T9:}{See T1.}
#' }
#'
#' @docType data
#' @source Table 1.2, page 18 of Prentice and Zhao 2019, taken from
#' Table 9.2 on page 292 in The Statistical Analysis of Failure Time Data,
#' Second Edition, by J. D. Kalbfleisch and R. L. Prentice. New York: Wiley
#' and Sons, 2002.
#' @references
#' Prentice, R., Zhao, S. "The statistical analysis of multivariate
#' failure time data: A marginal modeling approach", CRC Press (2019).
"bladder"


#' Simulated Clayton-Oakes data.
#'
#' Data simulated under Clayton-Oakes model (Equation 1.7, page 7 of Prentice
#' and Zhao) at various theta values and n = 50, with each observed time having
#' a 1/3 probability of being censored. Variate T_i is right-censored by C_i
#' for i = 1,2 (see pages 52 and 53 of Prentice and Zhao 2019), and
#' min(T_i,C_i) = S_i, with delta_i = 1 if T_i = S_i, and 0 otherwise
#' (i.e., if T_i = S_i then the data was not censored).
#'
#' @usage data(clayton_sim)
#' @format A data frame with 150 rows and 5 variables:
#' \describe{
#'  \item{Theta:}{Parameter value in Clayton-Oakes model.}
#'  \item{S1:}{min(T_1, C_1), where T_1 is time of interest and C_1 is a
#'  censoring variable. We see T_1 if it is not censored (i.e. C_1 > T_1),
#'  and we see C_1 otherwise.}
#'  \item{d1:}{I(T_1=S_1), i.e. 1 if T_1 = S_1 and 0 otherwise.}
#'  \item{S2:}{See S1.}
#'  \item{d2:}{See d1.}
#' }
#'
#' @docType data
#' @source Table 3.4, page 70 of Prentice and Zhao 2019.
#' @references
#' Prentice, R., Zhao, S. "The statistical analysis of multivariate
#' failure time data: A marginal modeling approach", CRC Press (2019).
"clayton_sim"
