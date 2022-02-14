#generate test parameters to use with the model

setwd("~/GitHub/tb-natural-history")

params <- list("p_m"=0.034, #smear progression (from subclinical)
               "p_s"=0.092, #symptom progression (from smear-negative)
               "r_m"=0.015, #smear regression (from subclinical)
               "r_s"=0.35, #symptom regression (from smear-negative)
               "c_sp"=0.37, #spontaneous resolution (from subclinical smear-negative)
               "c_tx"=0.058, #treatment (from symptomatic smear-negative)
               "m_tb"=0.0018, #TB mortality (from symptomatic smear-positive)
               "a_p_m"=1.2, #RR of smear progression if symptomatic (vs. subclinical)
               "a_r_m"=0.077, #RR of smear regression if symptomatic (vs. subclinical)
               "a_p_s"=1.2, #RR of symptom progression if smear-positive (vs. smear-negative)
               "a_r_s"=0.077, #RR of symptom regression if smear-positive (vs. smear-negative)
               "a_m"=9.5, #RR of TB mortality if symptomatic smear-positive (vs. symptomatic smear-negative)
               "a_tx"=2.3, #RR of treatment if symptomatic smear-positive (vs. symptomatic smear-negative)
               "p_c"=0, #return from the spontaneously resolved state
               "m_ac"=0.0004 #non-TB mortality
               )

save(params, file="data/params_test.Rda")
