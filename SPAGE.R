options(stringsAsFactors = F)
library(data.table)
library(SPAGE)
library(SAIGE)

args <- commandArgs(trailingOnly = TRUE)

hh <- paste(unlist(args),collapse=' ')

listoptions <- unlist(strsplit(hh,'--'))[-1]
options.args <- sapply(listoptions,function(x){
  unlist(strsplit(x, ' '))[-1]
}, simplify=FALSE)
options.names <- sapply(listoptions,function(x){
  option <-  unlist(strsplit(x, ' '))[1]
})
names(options.args) <- unlist(options.names)


spage.options <- c('bgen', 'bgen-bgi', 'variant-name-file', 'pheno-file', 'pheno-name',
             'environmental-factors', 'covar-names', 'sampleid-name', 'delimiter',
             'cutoff', 'impute-method', 'missing-cutoff', 'min-maf', 'Firth-cutoff',
             'BetaG-cutoff', 'BetaG-SPA', 'G-Model', 'minMAC', 'out', 'help')

if(!all(names(options.args) %in% spage.options)){
   stop(paste0('Option(s) ', paste(names(options.args)[!(names(options.args) %in% spage.options)]), ' invalid'))
}

if('help' %in% names(options.args)){
   for(i in 1:length(spage.options)){
       print(spage.options[i])
   }
   quit()
}

if(length(options.args$bgen) == 0){
   stop('No bgen file specified')
}
if(length(options.args$`bgen-bgi`) == 0){
  stop('No bgen index file specified')
}
if(length(options.args$`variant-name-file`) == 0){
  stop('No file with variant names specified')
}
if(length(options.args$`pheno-file`) == 0){
  stop('No phenotype file specified')
} else if(length(options.args$`pheno-file`) > 1){
  stop('Only one phenotype file can be specified')
}
if(length(options.args$`pheno-name`) > 1){
   stop('Only one --pheno-name can be specified.')
}else if(length(options.args$`pheno-name`) == 0){
   stop('--pheno-name is not specified.')
}
if(length(options.args$`environmental-factors`) == 0){
   stop('No environmental factor specified')
}
if(any(options.args$`covar-names` %in% options.args$`environmental-factors`)){
  stop(paste0('Covariates ', 
              paste0(options.args$`covar-names`[which(options.args$`covar-names` %in% options.args$`environmental-factors`)], collapse = ' '),
              ' also specified as environmental factor'))
}
if(length(options.args$`sampleid-name`) == 0){
  stop('The name of the sample identifier column in the phenotype file must be specified')
}
if(length(options.args$delimiter) == 0){
   options.args$delimiter[1] = ','
}
if(length(options.args$cutoff) == 0){
   options.args$cutoff[1] = 2
}
if(length(options.args$`impute-method`) == 0){
  options.args$`impute-method`[1] = "none"
}
if(length(options.args$`missing-cutoff`) == 0){
  options.args$`missing-cutoff`[1] = 0.15
}
if(length(options.args$`min-maf`) == 0){
  options.args$`min-maf`[1] = 0
}
if(length(options.args$`Firth-cutoff`) == 0){
  options.args$`Firth-cutoff`[1] = 0
}
if(length(options.args$`BetaG-cutoff`) == 0){
  options.args$`BetaG-cutoff`[1] = 0.15
}
if(length(options.args$`BetaG-SPA`) == 0){
  options.args$`BetaG-SPA`[1] = F
}else{
  if(!(toupper(options.args$`BetaG.SPA`[1]) %in% c('T', 'F', "TRUE", "FALSE"))){
     stop('--BetaG.SPA must be specified as T, F, TRUE, or FALSE')
  }
}
if(length(options.args$`G-Model`) == 0){
  options.args$`G-Model`[1] = "Add"
}
if(length(options.args$minMAC) == 0){
  options.args$minMAC[1] = 20
}
if(length(options.args$out) == 0){
   options.args$out[1] = 'spage.out'
}


### Read in phenotype data
data <- fread(options.args$`pheno-file`[1], sep = options.args$delimiter[1])
data <- as.data.frame(data)








### To query BGEN file
ids_to_include <- as.character(fread(options.args$`variant-name-file`[1])[,'rsid'])
head(ids_to_include)
ranges_to_include = data.frame(chromosome = NULL, start = NULL, end = NULL)
ranges_to_exclude = data.frame(chromosome = NULL, start = NULL, end = NULL)
ids_to_exclude    = as.character(vector())

Mtest = setgenoTest_bgenDosage(options.args$bgen[1],
                               options.args$`bgen-bgi`[1],
                               ranges_to_exclude = ranges_to_exclude,
                               ranges_to_include = ranges_to_include,
                               ids_to_exclude= ids_to_exclude,
                               ids_to_include=ids_to_include)

if(Mtest == 0){
   stop("Number of variants to be tested is 0")
}
SetSampleIdx(1:nrow(data), nrow(data))







### Create NULL object
if(length(options.args$`covar-names`) == 0){
   null.formula = formula(paste0(options.args$`pheno-name`, " ~ ",
                                 paste0(paste0(options.args$`environmental-factors`, collapse = " + "))))
   print('Null model fitting complete.')
   print(null.formula)
}else{
   null.formula = formula(paste0(options.args$`pheno-name`, " ~ ",
                                 paste0(paste0(options.args$`environmental-factors`, collapse = " + "),
                                        " + ",
                                        paste0(options.args$`covar-names`, collapse = " + "))))
   print('Null model fitting complete.')
   print(null.formula)
}

obj.null = SPAGE_Null_Model(null.formula, subjectID = data[,options.args$`sampleid-name`[1]], data=data, out_type="D")





### Get environmental factors
Envn.mtx <- as.matrix(data[,options.args$`environmental-factors`[1]])
rm(data)










print("Start Analyzing...")
print(paste0('Using SPAGE.one.SNP parameters...'))
cat(sprintf("Cutoff: %d", as.numeric(options.args$cutoff[1])), "\n")
cat(sprintf("impute-method: '%s'", options.args$`impute-method`[1]), "\n")
cat("missing.cutoff: ", as.numeric(options.args$`missing-cutoff`[1]), "\n")
cat("min.maf: ", as.numeric(options.args$`min-maf`[1]), "\n")
cat("Firth.cutoff: ", as.numeric(options.args$`Firth-cutoff`[1]), "\n")
cat("BetaG.cutoff: ", as.numeric(options.args$`BetaG-cutoff`[1]), "\n")
cat("BetaG.SPA: ", as.logical(options.args$`BetaG-SPA`[1]), "\n")
cat(sprintf("G.Model: '%s'", options.args$`G-Model`[1]), "\n")



OUT = NULL
idx = 1;
minMAC  = as.numeric(options.args$minMAC[1])
nSNPs.out = length(ids_to_include)
print(Sys.time())
while(idx <= Mtest){
       Gx = getDosage_bgen_withquery()

       g = Gx$dosages
       g1 = round(g)
       g1.case=g1[obj.null$y==1]
       MAC = min(sum(g1),sum(2-g1))
       MAC.case = min(sum(g1.case),sum(2-g1.case))
       AF = Gx$variants$AF
      if(AF < 0.5){
         MAF = AF
      }else{
         MAF = 1-AF
      }

     # rowHeader=as.vector(unlist(Gx$variants))
     SNPID = Gx$variants$SNPID
     rsid = Gx$variants$rsid
     CHR = Gx$variants$chromosome
     POS = Gx$variants$position
     REF = Gx$variants$allele0
     ALT = Gx$variants$allele1


     # if(MAC<minMAC | markerInfo < minInfo | MAC.case < minMAC.case){}
     if(MAC<minMAC){}
     else{
     out.tmp = SPAGE.one.SNP(g,
                              obj.null = obj.null,
                              Envn.mtx = Envn.mtx,
                              Cutoff = as.numeric(options.args$cutoff[1]),
                              impute.method = options.args$`impute-method`[1],
                              missing.cutoff = as.numeric(options.args$`missing-cutoff`[1]),
                              min.maf = as.numeric(options.args$`min-maf`[1]),
                              Firth.cutoff = as.numeric(options.args$`Firth-cutoff`[1]),
                              BetaG.cutoff = as.numeric(options.args$`BetaG-cutoff`[1]),
                              BetaG.SPA    = as.numeric(options.args$`BetaG-SPA`[1]),
                              G.Model = options.args$`G-Model`[1])
      out = c(rsid, round(MAC,2), MAC.case, CHR, POS, REF, ALT, out.tmp)
      OUT = rbind(OUT, out)
    }

#     ### summarize all results
    if(idx %% nSNPs.out == 0 | idx == Mtest){
       OUT = as.data.frame(OUT)
       fwrite(OUT, options.args$out[1], quote=FALSE, row.names=FALSE, col.names=FALSE, append = TRUE)
       OUT = NULL
    }
    idx=idx+1
 }
 print("Finished")
 print(Sys.time())

#########################################################


