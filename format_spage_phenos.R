library(readr)
library(dplyr)


args <- commandArgs(trailingOnly=T)
names(args) <- c("phenofile", "sample_id_header", "outcome", "exposure", "covar_names", "delimiter",
"missing", "samplefile")

phenos <- read_delim(args["phenofile"], delim=args["delimiter"], na=args["missing"])

if (args["samplefile"] != "") {  # Optional sample file ensures correct ordering of IDs
    ids <- read_delim(args["samplefile"], delim=" ", col_names=F, skip=2) %>%
	select(1) %>%
	setNames(args["sample_id_header"])
    phenos <- left_join(ids, phenos, by=unname(args["sample_id_header"]))
}

covars <- if (args["covar_names"] == "") character() else strsplit(args["covar_names"], " ")[[1]]
output_cols <- c(args["sample_id_header"], args["outcome"], args["exposure"], covars)	

phenos %>%
    select(unname(output_cols)) %>%
    write_csv("spage_phenotypes.csv")
