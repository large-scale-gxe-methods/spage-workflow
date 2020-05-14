task process_phenos {
	
	File phenofile
	File? samplefile
	String sample_id_header
	String outcome
	String exposure
	String covar_names
	String delimiter
	String missing
	Int ppmem

	command {
		Rscript /format_spage_phenos.R ${phenofile} ${sample_id_header} ${outcome} ${exposure} "${covar_names}" "${delimiter}" ${missing} "${samplefile}"
	}

	runtime {
		docker: "quay.io/large-scale-gxe-methods/spage-workflow"
		memory: ppmem + "GB"
	}

        output {
                File pheno_fmt = "spage_phenotypes.csv"
	}
}

task run_tests {

	File genofile
	Float maf
	Int mac
	File? samplefile
	File phenofile
	String sample_id_header
	String outcome
	String exposure_names
	String? covar_names
	String missing
	Int memory
	Int disk
	Int monitoring_freq

	command {
		$BGENIX -g ${genofile} -index 
		$BGENIX -g ${genofile} -list > variants.txt

		dstat -c -d -m --nocolor ${monitoring_freq} > system_resource_usage.log &
		atop -x -P PRM ${monitoring_freq} | grep '(R)' > process_resource_usage.log &

		Rscript /SPAGE.R \
			--bgen ${genofile} \
			--bgen-bgi ${genofile}.bgi \
			--variant-name-file variants.txt \
			--pheno-file ${phenofile} \
			--pheno-name ${outcome} \
			--environmental-factors ${exposure_names} \
			${"--covar-names " + covar_names} \
			--sampleid-name ${sample_id_header} \
			--delimiter , \
			--min-maf ${maf} \
			--minMAC ${mac} \
			--out spage_res
	}

	runtime {
		docker: "quay.io/large-scale-gxe-methods/spage-workflow"
		memory: "${memory} GB"
		disks: "local-disk ${disk} HDD"
	}

	output {
		File out = "spage_res"
		File system_resource_usage = "system_resource_usage.log"
		File process_resource_usage = "process_resource_usage.log"
	}
}

task standardize_output {

	File resfile

	command {
		Rscript /format_spage_output.R ${resfile}
	}

	runtime {
		docker: "quay.io/large-scale-gxe-methods/spage-workflow"
		memory: "2 GB"
	}

        output {
                File res_fmt = "res_fmt.txt"
	}
}

task cat_results {

	Array[File] results_array

	command {
		head -1 ${results_array[0]} > all_results.txt && \
			for res in ${sep=" " results_array}; do tail -n +2 $res >> all_results.txt; done
	}
	
	runtime {
		docker: "ubuntu:latest"
		disks: "local-disk 10 HDD"
	}
	output {
		File all_results = "all_results.txt"
	}
}


workflow run_SPAGE {

	Array[File] genofiles
	Float? maf = 0.005
	Int? mac = 10
	File? samplefile
	File phenofile
	String? sample_id_header = "sampleID"
	String outcome
	String exposure_names
	String? covar_names
	String? delimiter = ","
	String? missing = "NA"
	Int? memory = 10
	Int? disk = 50
	Int? monitoring_freq = 1

	Int ppmem = 2 * ceil(size(phenofile, "GB")) + 1

	call process_phenos {
		input:
			phenofile = phenofile,
			samplefile = samplefile,
			sample_id_header = sample_id_header,
			outcome = outcome,
			exposure = exposure_names,
			covar_names = covar_names,
			delimiter = delimiter,
			missing = missing,
			ppmem = ppmem
	}

	scatter (i in range(length(genofiles))) {
		call run_tests {
			input:
				genofile = genofiles[i],
				maf = maf,
				mac = mac,
				samplefile = samplefile,
				phenofile = process_phenos.pheno_fmt,
				sample_id_header = sample_id_header,
				outcome = outcome,
				exposure_names = exposure_names,
				covar_names = covar_names,
				missing = missing,
				memory = memory,
				disk = disk,
				monitoring_freq = monitoring_freq
		}
	}

	scatter (resfile in run_tests.out) {
		call standardize_output {
			input:
				resfile = resfile
		}
	}	

	call cat_results {
		input:
			results_array = standardize_output.res_fmt
	}

	output {
		File results = cat_results.all_results
		Array[File] system_resource_usage = run_tests.system_resource_usage
		Array[File] process_resource_usage = run_tests.process_resource_usage
	}

	parameter_meta {
		genofiles: "Array of genotype filepaths in .bgen format."
		maf: "Minimum minor allele frequency threshold for pre-filtering variants as a fraction (default is 0.005)."
		samplefile: "Optional .sample file accompanying the .bgen file. Required for proper function if .bgen does not store sample identifiers. NOT YET IMPLEMENTED -- phenotype file currently must exactly match genotype file sample set and ordering."
		phenofile: "Phenotype filepath."	
		sample_id_header: "Optional column header name of sample ID in phenotype file."
		outcome: "Column header name of phenotype data in phenotype file."
		exposure_names: "Column header name(s) of the exposures for genotype interaction testing (space-delimited)."
		covar_names: "Column header name(s) of any covariates for which only main effects should be included selected covariates in the pheno data file (space-delimited)."
		delimiter: "Delimiter used in the phenotype file."
		missing: "Missing value key of phenotype file."
		memory: "Requested memory (in GB)."
		disk: "Requested disk space (in GB)."
		monitoring_freq: "Delay between each output for process monitoring (in seconds). Default is 1 second."
	}

        meta {
                author: "Kenny Westerman"
                email: "kewesterman@mgh.harvard.edu"
                description: "Run interaction tests for binary traits using SPAGE and return a table of summary statistics for K-DF interaction and (K+1)-DF joint tests. NOTE: This workflow only includes the basic interaction testing step calling the SPAGE R function. It does not yet include input phenotype alignment to a sample file or output processing to match a pre-specified formatting and set of fields."
        }
}

