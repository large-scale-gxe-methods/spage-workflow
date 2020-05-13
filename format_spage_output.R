library(readr)
library(dplyr)


args <- commandArgs(trailingOnly=T)
names(args) <- "resfile"

read_csv(args["resfile"],
	 col_names=c("rsID", "MAC", "MAC.case", "CHR", "POS", "REF", "ALT", 
		     "MAF", "missing", "P_Value_Marginal",
		     "P_Value_Interaction_SPA", "P_Value_Interaction_Norm", "P_Value_Interaction_Firth",
		     "Stat", "Var", "z") %>%
    select(rsID, CHR, POS, REF, ALT, MAF, P_Value_Marginal, P_Value_Interaction=P_Value_Interaction_Norm) %>%
    write_delim("res_fmt.txt", delim=" ")
