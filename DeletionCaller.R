## file = path to rds object output from mgatk
#cell_heteroplasmy_file = path to heteroplasmy csv, cell_heteroplasmic_df.tsv.gz
# mtDNA_file = path to mtDNA reference , chrM_refAllele.txt

DeletionDetection = function(sample, o_call = 'tenx_call', min_readcount = 1, file, cell_heteroplasmy_file, mtDNA_file){
  
  object = readRDS(file)
  variants = fread(cell_heteroplasmy_file)
  variants = as.data.frame(variants)
  rownames(variants) = variants$V1
  variants = variants[, -1]
  
  coverage = rowSums(object@assays@data$coverage)
  
  coverage_df = object@assays@data$coverage
  colnames(coverage_df) =  object@colData$sample
  
  mtdna = read.csv(mtDNA_file, header = F, sep = '\t')
  
  diff = coverage[2:length(coverage)]-coverage[1:(length(coverage)-1)]
  zscores = (diff-mean(diff))/sd(diff)
  potential = which(abs(zscores)>3.5)
  
  opp_sign = c()
  closeness = c()
  for(i in 1:(length(potential)-1)){
    closeness[i] = potential[(i+1)]-potential[i]
    opp_sign[i] = (zscores[potential[i]] <0) & (zscores[potential[(i+1)]]> 0)
  }
  evaluation = which(closeness ==1 &opp_sign)
  if(length(evaluation)> 0){
    print('potential deletions found')
    potential = potential[evaluation]
    
    for(i in 1:length(potential)){
      position = potential[i]+1
      ref_allele = toupper(mtdna[position, 2])
      assays = paste(ref_allele, c('_counts_fw', '_counts_rev'), sep = '')
      
      ### estimate cell coverage 
      
      est_cov = coverage_df[c((position-5):(position-1), (position+1):position+5), rownames(variants)]
      est_cov = colMeans(est_cov)
      
      ### estimate ref_coverage list 
      
      ref_allele_coverage = rep(0, length(est_cov))
      for(k in 1:length(assays)){
        x = object@assays@data[names(object@assays@data)==assays[k]][[1]]
        colnames(x) =  object@colData$sample
        ref_allele_coverage =ref_allele_coverage +  x[position, rownames(variants)] 
        names(object@assays@data)
      }
      
      
      
      del_coverage = est_cov-ref_allele_coverage 
      name = paste(c(position, ref_allele, '>NULL'), sep = '', collapse = '')
      del_coverage = del_coverage/est_cov
      
      del_coverage[(!est_cov>min_readcount) | est_cov == 0] = NA
      del_coverage[del_coverage<0] = 1
      
      variants[name] = del_coverage
    }
    print('done with adding deletions')
    
    
  }
  
  print('adding full data')
  FullData = data.frame(row.names = object@colData$sample)
  for(i in 1:ncol(variants)){
    name = colnames(variants)[i]
    position = strsplit(name, '')[[1]]
    loc = which(position =='>')
    ref = position[loc-1]
    alt = position[loc+1]
    position = as.numeric(paste(position[1:(loc-2)], sep = '', collapse = ''))
    
    if(alt =='N'){
      ##this is a deletion 
      assays = paste(ref, c('_counts_fw', '_counts_rev'), sep = '')
      
      ### estimate cell coverage 
      
      est_cov = coverage_df[c((position-5):(position-1), (position+1):position+5), rownames(FullData)]
      est_cov = colMeans(est_cov)
      ### estimate ref_coverage list 
      
      ref_allele_coverage = rep(0, length(est_cov))
      for(k in 1:length(assays)){
        x = object@assays@data[names(object@assays@data)==assays[k]][[1]]
        colnames(x) =  object@colData$sample
        ref_allele_coverage =ref_allele_coverage +  x[position, rownames(FullData)] 
        
      }
      
      del_coverage = est_cov-ref_allele_coverage 
      del_coverage = del_coverage/est_cov
      
      del_coverage[(!est_cov>min_readcount) | est_cov == 0] = NA
      del_coverage[del_coverage<0] = 1
      FullData[name] = del_coverage
    } else{
      assays = paste(alt, c('_counts_fw', '_counts_rev'), sep = '')
      
      ### estimate cell coverage 
      est_cov = coverage_df[position, rownames(FullData)]
      
      ### estimate ref_coverage list 
      
      ref_allele_coverage = rep(0, length(est_cov))
      for(k in 1:length(assays)){
        x = object@assays@data[names(object@assays@data)==assays[k]][[1]]
        colnames(x) =  object@colData$sample
        ref_allele_coverage =ref_allele_coverage +  x[position, rownames(FullData)] 
        
      }
      
      ref_allele_coverage = ref_allele_coverage/est_cov
      
      ref_allele_coverage[(!est_cov>min_readcount) | est_cov == 0] = NA
      FullData[name] = ref_allele_coverage
    }
    
  }
  
  return(list(FullData = FullData, variants = variants))
  
  
} 


# Run the original function
results <- DeletionDetection(
  sample = "samplename",
  o_call = "tenx_call",
  min_readcount = 1,
  file = "~/path/to/regular_mgatk_call.rds",
  cell_heteroplasmy_file = "~/path/to/regular_mgatk_call.cell_heteroplasmic_df.tsv.gz",
  mtDNA_file = "~/path/to/chrM_refAllele.txt"
)

