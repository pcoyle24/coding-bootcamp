clear
cls
cd /Users/philipcoyle/Documents/School/University_of_Wisconsin/SecondYear/Summer_2020/CodingBootcamp/ProblemSets/PS4
log using PS4_Q1.log, replace
import delimited "lwage.csv"

* Housekeeping
rename v1 lwage
rename v2 college
rename v3 exp
gen exp2 = exp^2

* Regression
reg lwage college exp exp2, r

log close


