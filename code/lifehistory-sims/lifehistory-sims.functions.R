### Demographic modeling  ######
#### Pre-formatting ######

## eDem_anctree ###
## Simulates many ancestral trees, without mutations, to be posteriorly sampled
## for an ancestral pop
eDem_anctree <- function(params,N,parallel=FALSE,cores=1) {
  require(reticulate)
  require(tidyverse)
  require(parallel)
  require(tictoc)
  ms <- import('msprime')
  tic()
  if (parallel) {
    message('Simulating ',N,' ancestral trees')
    message('Pop size: ',params['pop_size'])
    message('Sample size: ',params['samples'])
    message('Sequence Length: ',params['seq_length'])
    message('Recombination Rate: ',params['rec_rate'])
    message('Saving to: ',params['path'])
    message('Running in parallel with ',cores,' cores.')
    ts <- mclapply(1:N,function(x){return(ms$sim_ancestry(samples=as.numeric(params['samples']),
                                                          recombination_rate=as.numeric(params['rec_rate']),
                                                          sequence_length=as.numeric(params['seq_length']),
                                                          population_size=as.numeric(params['pop_size'])))},
                   mc.cores = cores)
    mclapply(1:N,function(x){ts[[x]]$dump(paste0(params['path'],'/anctree_',x,'.tre'))},mc.cores = cores)
  } else {
    message('Simulating ',N,' ancestral trees')
    message('Pop size: ',params['pop_size'])
    message('Sample size: ',params['samples'])
    message('Sequence Length: ',params['seq_length'])
    message('Recombination Rate: ',params['rec_rate'])
    message('Saving to: ',params['path'])
    ts <- lapply(1:N,function(x){return(ms$sim_ancestry(samples=as.numeric(params['samples']),
                                                        recombination_rate=as.numeric(params['rec_rate']),
                                                        sequence_length=as.numeric(params['seq_length']),
                                                        population_size=as.numeric(params['pop_size'])))})
    lapply(1:N,function(x){ts[[x]]$dump(paste0(params['path'],'/anctree_',x,'.tre'))})
  }
  toc()
  invisible()
}

## eDem_ancpop ###
## Function to make one ancestral pop by:
## 1) Drawing random SNPs from the previously simulated trees or simulating
##    from scratch if it is a sequence and not SNP
## 2) Adding mutations to the trees
## 3) Saving ancPop with a random seed for identifier
eDem_ancpop <- function(G,snp=T,mut_rate,ancTree_path='./',save_to_temp=T,outpath='/.',
                        samples=NULL,rec_rate=NULL,seq_length=NULL,pop_size=NULL) {
  require(tidyverse)
  require(reticulate)
  require(pegas)
  require(vcfR)
  base <- import_builtins()
  tskit <- import('tskit')
  msprime <- import('msprime')
  ## Check if SNP is T or F
  ## If T, draw from previously simulated SNPS, add mutations, save individual VCFs,
  ##   Load VCFs back with pegas, merge and output unified VCF. Use temporary folder
  ## If F, call `msprime` and simulate sequence. Add mutations, save VCF to temporary folder.
  ## Use seed number to identify the ancestral pop in the temp folder
  ## Return temporary folder and seed number to parent function
  if (snp) {
    message('Creating ancestral pop - merging SNPs:')
    files <- list.files(ancTree_path,full.names = T)
    if (length(files) >= G) {
      randomtrees <- files[sample(1:length(files),G)]
    } else {
      stop('Requested number of SNPs higher than the number of available simulated SNPs!')
    }
    dir <- tempdir()
    file.remove(list.files(dir,full.names = T))
    
    ## Getting simulated trees and adding mutations
    message('Adding mutations...')
    trees <- lapply(1:G,function(x){tskit$load(randomtrees[x])}) %>%
      lapply(msprime$sim_mutations,rate=mut_rate,keep=TRUE)
    
    # Saving to VCF
    message('Merging SNPs...')
    lapply(1:G,function(x){
      f <- base$open(paste0(dir,'/snp',x,'.vcf'),'w')
      trees[[x]]$write_vcf(f,position_transform='legacy')
      f$close()
    }
    )
    vcfs <- list.files(dir,pattern = '.vcf',full.names = T) %>% 
      lapply(vcfR::read.vcfR,verbose=FALSE)
    # Saving first vcf for further editing
    vcf <- vcfs[[1]]
    # Extracting random SNP per vcf
    data <- lapply(vcfs,function(x){
      if (nrow(x@gt>0)) {
        # If there's at least one row in the GT matrix, we sample one of those rows randomly.
        var <- sample(1:nrow(x@gt),1)
        fix <- x@fix[var,]
        gt <- x@gt[var,]
      } else {
        # If not, we return NULL
        fix <- gt <- NULL
      }
      return(list(fix,gt))
    })
    
    # Merging SNP data
    fix <- lapply(data,function(x){return(x[[1]])}) %>% do.call(what = rbind)
    gt <- lapply(data,function(x){return(x[[2]])}) %>% do.call(what = rbind)
    
    # Changing positions by randomly selecting positions from 1:G
    # and sorting
    # ID will be a sequence vector
    fix[,2] <- sample(1:G,nrow(fix)) %>% sort() %>% as.character()
    fix[,3] <- as.character(seq(1:nrow(fix)))
    
    # Editing final VCF
    vcf@meta[4] <- paste0('##contig=<ID=1,length=',G,'>')
    vcf@fix <- fix
    vcf@gt <- gt
    message('Saving to VCF file...')
    seed = sample(1:1000000,1)
    if (save_to_temp) {
      ancPop <- paste0(dir,'/ancPop_',seed,'.vcf.gz')
    } else {
      ancPop <- paste0(outpath,'/ancPop_',seed,'.vcf.gz')
    }
    vcfR::write.vcf(vcf,ancPop)
    R.utils::gunzip(ancPop)
    # Removing gz so SLIM can find the file
    ancPop <- gsub('.gz$','',ancPop)
  } else {
    message('Simulating ancestral pop - sequence:')
    message('Pop size: ',pop_size)
    message('Sample size: ',samples)
    message('Sequence Length: ',seq_length)
    message('Recombination Rate: ',rec_rate)
    tree <- ms$sim_ancestry(samples=as.numeric(samples),recombination_rate=as.numeric(rec_rate),
                            sequence_length=as.numeric(seq_length),population_size=as.numeric(pop_size))
    tree <- msprime$sim_mutations(tree,rate=mut_rate,keep=TRUE)
    seed = sample(1:1000000,1)
    if (save_to_temp) {
      ancPop <- paste0(dir,'/ancPop_',seed,'.vcf')
    } else {
      ancPop <- paste0(outpath,'/ancPop_',seed,'.vcf')
    }
    f <- base$open(ancPop,'w')
    trees$write_vcf(f,position_transform='legacy')
    f$close()
  }
  return(ancPop)
}

makeAnc <- function(x,path,name,snps=T,samples=500,rec=0,len=1000,size=10000,mut=1e-5, ncores = 1) {
  library(reticulate)
  ms <- import('msprime')
  require(fs)
  source_python('https://raw.githubusercontent.com/Rilquer/cloud/main/lifehistory-sims.functions.py')
  name <- name
  file <- paste0(path,'/',name,'_rep',x,'.vcf.gz')
  tdir <- tempdir()
  dir.create(tdir)
  dir <- paste0(tdir,'/rep',x)
  if (dir_exists(dir)) {
    dir_delete(dir)
  }
  if (snps) {
    # IN THE FUTUREl, gotta add the option to simulate the alternative to SNPs, e.g., mtDNA haplotypes.
    ms_info=data.frame(ms_path = paste0(dir,'/snp',1:(len+(len*0.5))),
                       samples=samples,rec_rate=rec,
                       seq_length=1,
                       pop_size=size,
                       mut_rate=mut)
    
    # Simulating several SNPs, filtering out tri-allelic and
    # merging them into one ancestral VCF
    br <- FALSE
    while (br == FALSE) {
      dir.create(dir)
      for (i in 1:nrow(ms_info)) {
        ms=demBurnin(ms_info[i,]) # Parallelize this if possible
      }
      files <- list.files(dir,full.names = T, pattern = '.vcf') %>%
        mclapply(read.vcfR, verbose = FALSE, mc.cores = ncores)
      tri <- sapply(1:length(files),function(y){length(grep(3,files[[y]]@gt))})
      vcf <- do.call(rbind,files[which(tri==0)])
      if (nrow(vcf) >= len) {
        br <- TRUE
        if (nrow(vcf) > len) {
          vcf <- vcf[sample(1:nrow(vcf),len)]
        }
      }
      dir_delete(dir)
    }
  }
  vcfR::write.vcf(vcf,file = file)
  R.utils::gunzip(file, overwrite = TRUE, remove = TRUE)
}

### demAncestral ###
demAncestral <- function(path,name,reps=1,snps=T,samples=500,rec=0,len=1000,size=10000,mut=1e-5, ncores = 1) {
  require(reticulate)
  require(vcfR)
  require(tictoc)
  require(parallel)
  message('Simulating ',reps,' ancestral pops:')
  message('N of SNPs: ',len)
  message('Pop size: ',size)
  message('Sample size: ',samples)
  message('Mutation rate: ',mut)
  tic()
  mclapply(1:reps,makeAnc,path = path,name = name,snps=snps,samples=samples,rec=rec,len=len,size=size,mut=mut, ncores = ncores)
  toc()
}

### demRender ###
# Function to create settings for simulations from list of parameters.
# This function will create the lists from which the scripts will be rendered.
# It will draw priors as needed, so the output of this function will be the
# actual parameters used for each simulation.

# Arguments will be the list of parameteres that can be set, both for
# slim and msprime

# This function mostly implement slim_script_render, and it is needed to transform
# the more user-friendly dataframe format of inputting the values to the list format.
# As of Dec 10, 2023, the dataframe format can't handle any string defining the
# path for file reading or writing.
# We also use this function to retrieve ngen before adding the final block
# to set the final generation separately.

## Params - Dataframe with parameters to be set. Number of rows correspond to models to be render. One value per parameter
## Ngen - how many generations the simulation should run by
## script - an R script creating all the slimr_blocks and merging them with slimr_script.
## The value of ngen is expected to be used here.
eDem_render <- function(params,script = 'script.R',parallel=FALSE,ncores=1) {
  # Add code to create ancestral pop
  # Call eDem_ancPop on each line of params
  # Add to params a column with the location of the ancestral pop for each
  nticks <<- unique(params$nticks)
  sim <- source(script,local = TRUE)
  require(tictoc)
  require(future)
  if (parallel == TRUE) {
    tic()
    plan(multisession, workers = ncores)
    message('Rendering ',nrow(params),' models with slimr::slim_script_render...')
    rendered <- slimr::slim_script_render(sim[[1]], template = apply(params,1,as.list), parallel=TRUE)
    message('Done')
    toc()
  } else {
    tic()
    rendered <- slimr::slim_script_render(sim[[1]], template = apply(params,1,as.list))
    toc()
  }
  return(rendered)
}

### demTest ###
# Quick test of of the script with a few subset parameters

### demRun ###
# Function to render and run each simulation.
# Arguments will be the list object with all sims plus number of cores to run
# in parallel.
# This function will also run msprime to create the initial genetic settings

eDem_run <- function(script,nrep=1,parallel=TRUE,ncores=1) {
  require(tictoc)
  require(future)
  require(parallel)
  script <- rep(script,nrep)
  if (parallel == TRUE) {
    tic()
    #plan(multisession, workers = ncores)
    #sr <- slim_run(script , parallel = TRUE, throw_error = FALSE)
    sr <- mclapply(script,slim_run, capture_output = TRUE, show_output = FALSE, throw_error = FALSE, mc.cores = ncores)
    toc()
  } else {
    tic()
    sr <- slim_run(script, capture_output = TRUE, show_output = FALSE, throw_error = FALSE)
    toc()
  }
  return(sr)
}

####  Post-processing   ####

### Basic functions to retract specific data from one slim run ###
# In the functions below, y = one element inside the output of demRun (i.e., one individual slimr output)

getParams <- function(y) {
  y$output_data <- y$output_data %>% rename(tick = generation) %>% 
    mutate_at(c("tick"),as.numeric)
  stop_tick <- max(y$output_data$tick)
  params_temp <- y$output_data %>% filter(tick==1) %>% select(name,data) %>%
    pivot_wider(id_cols = everything(),names_from = name,values_from = data) %>% 
    select(!(matches('_p[[:digit:]]$'))) %>% select(!(matches('_all$'))) %>% 
    mutate_at(vars(-c("name","ev_type","m","rprob")),as.numeric) %>%
    mutate_at(vars(c("m","rprob")),str_split,pattern=' ') %>% 
    mutate(stop_tick = stop_tick)
  return(params_temp)
}

getData <- function(y, final = FALSE) {
  y$output_data <- y$output_data %>% rename(tick = generation) %>% 
    mutate_at(c("tick"),as.numeric)
  if (final==TRUE) {
    last_tick <- max(y$output_data$tick)
  } else {
    last_tick = y$output_data$tick
  }
  data_temp <- y$output_data %>% filter(tick %in% last_tick) %>% 
    select(tick,name,data) %>%
    pivot_wider(id_cols = everything(),names_from = name,values_from = data) %>%
    select(c(tick,matches('_p[[:digit:]]$'),matches('_all$'))) %>%
    mutate_at(vars(matches('^pot[[:alpha:]]'),matches('^act[[:alpha:]]'),matches('^off[[:alpha:]]'),matches('^Rep[[:alpha:]]')),
              str_split,pattern=' ') %>% 
    mutate_if(is.character,as.numeric)
  return(data_temp)
}

getVCF <- function(y) {
  require(tidyverse)
  # vcf_p1 <- y$output_data %>% filter(name=='p1_VCF') %>% select(data)
  # # We check if VCF output exists. It won't, in case there was an error
  # # If no VCF exists, returns NULL.
  # if (nrow(vcf_p1)!=0) {
  #   vcf_p1 <- vcf_p1 %>% as.character() %>%
  #     gsub(pattern = ".*##fileformat*", replacement = "##fileformat")
  # } else {
  #   vcf_p1 <- NULL
  # }
  # vcf_p2 <- y$output_data %>% filter(name=='p2_VCF') %>% select(data)
  # if (nrow(vcf_p2)!=0) {
  #   vcf_p2 <- vcf_p2 %>% as.character() %>%
  #     gsub(pattern = ".*##fileformat*", replacement = "##fileformat")
  # } else {
  #   vcf_p2 <- NULL
  # }
  vcf_p3 <- y$output_data %>% filter(name=='vcf') %>% select(data)
  if (nrow(vcf_p3)!=0) {
    vcf_p3 <- vcf_p3 %>% as.character() %>%
      gsub(pattern = ".*##fileformat*", replacement = "##fileformat")
  } else {
    vcf_p3 <- NULL
  }
  #return(list(vcf_p1,vcf_p2,vcf_p3))
  return(list(vcf_p3))
}

### Functions implementing the basic funcitons above to retrieve info from all runs in a slim_run output ###
# In the functions below, x = output from demRun

### demResults ###
# Return `slimr` output separately: model params, data for each cycle and VCF data.

demResults <- function(x, parallel=FALSE,ncores = 1) {
  require(tidyverse)
  require(tictoc)
  tic()
  getResults <- function(y) {
    out <- vector('list',3)
  
    # Retrieving sim_info
    out[[1]] <- getParams(y)
  
    # Retrieving gen_info
    out[[2]] <- getData(y,final=FALSE)
    
    # Retrieving VCF
    #out[[3]] <- getVCF(y)
  
    names(out) <- c('sim_data','tick_data')
    #names(out) <- c('sim_data','gen_data','VCF')
    return(out)
  }
  if (parallel) {
    temp <- mclapply(x,getResults, mc.cores = ncores)
  } else {
    temp <- lapply(x,getResults)
  }
  # Checking if x is a list of several slimr sims (where names is null), or just one slimr sim.
  #if (is.null(names(x))) {
  #  temp <- mclapply(x,getResults,ncores = ncores)
  #} else {
  #  temp <- getResults(x)
  #}
  return(temp)
  toc()
}

### demCheckFailed ###
# Check all runs and return dataframe of runs that suceeded, runs that failed
# and their respective model parameters and the error message
demCheckFailed <- function(x) {
  getErPar <- function(y) {
    if (y$exit_status==1) {
      params_temp <- getParams(y) %>%
        mutate(exit_status = y$exit_status,error_message = y$error[1])
      return(params_temp)
    }
  }
  getErData <- function(y) {
    if (y$exit_status==1) {
      data_temp <- getData(y,final=FALSE)
      return(data_temp)
    }
  }
  if (is.null(names(x))) {
    message('Retrieving simulation parameters...')
    erPar <- lapply(x,getErPar) %>% do.call(what = rbind.data.frame)
    message('Retrieving simulation data...')
    erData <- lapply(x,getErData)
    erData[sapply(erData, is.null)] <- NULL # Removing NULL elements
    temp <- list(erPar,erData)
    names(temp) <- c('Parameter','Sim_data')
  } else {
    message('Retrieving simulation parameters...')
    erPar <- getErPar(x)
    message('Retrieving simulation data...')
    erData <- getErData(x)
    temp <- list(erPar,erData)
    names(temp) <- c('Parameter','Sim_data')
  }
  return(temp)
}

## demSimData ##
## Getting params and data at final cycle from all runs
## Only gets those that did NOT return exit status = 1, i.e., no error was found.
demSimData <- function(x) {
  require(tidyverse)
  require(tictoc)
  tic()
  getInfo <- function(y) {
    if (y$exit_status!=1) {
      params_temp <- getParams(y)
      data_temp <- getData(y,final=TRUE) %>% select(-generation)
      return(tibble(params_temp,data_temp))
    }
  }
  if (is.null(names(x))) {
    temp <- lapply(x,getInfo) %>% do.call(what = rbind.data.frame)
  } else {
    temp <- getInfo(x)
  }
  return(temp)
  toc()
}

## demParams ##
## Getting parameters of all runs
demParams <- function(x) {
  require(tidyverse)
  require(tictoc)
  tic()
  if (is.null(names(x))) {
    temp <- lapply(x,getParams) %>% do.call(what = rbind.data.frame)
  } else {
    temp <- getParams(x)
  }
  return(temp)
  toc()
}

#######     Calculating pop gen sum stats    #######

## demFormat ##
# Format VCF from sims output into `vcfR` objects for stats calculation
# x = output from demRun
demFormat <- function(x,ncores=0) {
  require(vcfR)
  require(tidyverse)
  require(parallel)
  message('Retrieving VCFs from simulations...')
  if (is.null(names(x))) {
    if (ncores!=0) {
      tic()
      vcfs <- mclapply(x,getVCF,mc.cores = ncores)
      toc()
    } else {
      tic()
      vcfs <- mclapply(x,getVCF)
      toc()
    }
  } else {
    vcfs <- list(getVCF(x))
  }
  dir <- tempdir()
  # Checking non-null VCFs
  id <- sapply(1:length(vcfs),function(y){ifelse(!(is.null(vcfs[[y]][[3]])),return(y),NA)})
  id <- id[which(!(is.na(id)))]
  
  if (ncores!=0) {
    tic()
    message('Converting to `vcfR` object...')
    temp <- mclapply(id,function(y){write_lines(vcfs[[y]][[3]],file = paste0(dir,'/vcf',y))},mc.cores = ncores)
    data <- mclapply(id,function(y){read.vcfR(file = paste0(dir,'/vcf',y),verbose = FALSE)},mc.cores = ncores)
    toc()
  } else {
    tic()
    temp <- mclapply(id,function(y){write_lines(vcfs[[y]][[3]],file = paste0(dir,'/vcf',y))})
    data <- mclapply(id,function(y){read.vcfR(file = paste0(dir,'/vcf',y),verbose = FALSE)})
    toc()
  }

  if (length(vcfs)==1) {
    data <- data[[1]]
  }
  return(data)
}

# General function to bin vectors
bin_vec <- function(z,n=10, sort=FALSE,max=FALSE) {
  if (sort) {
    z=sort(z, decreasing = T)
  }
  splitted <- split(z,  cut(seq_along(z), n, labels = FALSE))
  if (max) {
    return(sapply(splitted,max))
  } else {
    return(sapply(splitted,sum))  
  }
}

# Get stats from genlight or DNABin object
# x = either vcfR or DNABin objects. Single element or a list.
demStats <- function(x, sequence=FALSE, ncores=0) {
  require(vcfR)
  require(pegas)
  require(ape)
  require(tidyverse)
  require(hillR)
  require(tictoc)
  calcStats <- function(y) {
    if (sequence) {
      if (class(y)=='vcfR') {
        y <- vcfR2DNAbin(y)
      }
      # Calculating on raw sequences
      total_pi <- pegas::nuc.div(y, pairwise.deletion=TRUE)
      hd <- pegas::hap.div(y)
      r2 <- pegas::R2.test(y,B=0,quiet=TRUE,plot=FALSE)
      tD <- pegas::tajima.test(y)
      tD <- tD[[1]] %>% as.numeric()
      
      sfs_vec <- as.vector(pegas::site.spectrum(y)) %>% bin_vec()
      sfs <- sum(sapply(1:(length(sfs_vec)-1),function(x){return(sfs_vec[x]-sfs_vec[x+1])}))
      #sfs_vec <- sfs_vec %>% as_tibble_row()
      #colnames(sfs_vec) <- paste0('sfs',1:10)
      
      d <- ape::dist.dna(y, model = 'N',pairwise.deletion = TRUE) %>% as.vector()
      #d_vec <- bin_vec(d, sort=TRUE) %>% as_tibble_row()
      #colnames(d_vec) <- paste0('d',1:10)
      max_d <- max(d)
      avg_d <- mean(d)
      pct75 <- as.numeric(quantile(d)[4])
      d_abvpct75 <- length(which(d>pct75))/length(d)
      h1 <- hillR::hill_taxa(round(d,digits = 2),q=1)
      h4 <- hillR::hill_taxa(round(d,digits = 2),q=4)
      hdiff <- h1-h4
      
      # Calculating on haplotype object
      y <- haplotype(y)
      total_pi_hap <- pegas::nuc.div(y, pairwise.deletion=TRUE)
      hd_hap <- pegas::hap.div(y)
      r2_hap <- pegas::R2.test(y,B=0,quiet=TRUE,plot=FALSE)
      tD_hap <- pegas::tajima.test(y)
      tD_hap <- tD_hap[[1]] %>% as.numeric()
      sfs_vec_hap <- as.vector(pegas::site.spectrum(y)) %>% bin_vec()
      sfs_hap <- sum(sapply(1:(length(sfs_vec_hap)-1),function(x){return(sfs_vec_hap[x]-sfs_vec_hap[x+1])}))
      #sfs_vec_hap <- sfs_vec_hap %>% as_tibble_row()
      #colnames(sfs_vec_hap) <- paste0('sfs_hap',1:10)
      d_hap <- ape::dist.dna(y, model = 'N',pairwise.deletion = TRUE) %>% as.vector()
      #d_hap_vec <- bin_vec(d_hap, sort=TRUE) %>% as_tibble_row()
      #colnames(d_hap_vec) <- paste0('d_hap',1:10)
      max_d_hap <- max(d_hap)
      avg_d_hap <- mean(d_hap)
      pct75_hap <- as.numeric(quantile(d_hap)[4])
      d_abvpct75_hap <- length(which(d_hap>pct75_hap))/length(d_hap)
      h1_hap <- hillR::hill_taxa(round(d_hap,digits = 2),q=1)
      h4_hap <- hillR::hill_taxa(round(d_hap,digits = 2),q=4)
      hdiff_hap <- h1_hap-h4_hap
      
      ss <- tibble(total_pi,hd,r2,tD,max_d,avg_d,pct75,d_abvpct75,h1,h4,hdiff,sfs,
                   total_pi_hap,hd_hap,r2_hap,tD_hap,max_d_hap,avg_d_hap,pct75_hap,d_abvpct75_hap,
             h1_hap,h4_hap,hdiff_hap,sfs_hap)
      #ss <- tibble(ss,sfs_vec,d_vec,ss_hap,sfs_vec_hap,d_hap_vec)
      return(ss)
    } else {
      y_loci <- vcfR2loci(y)
      hs <- mean(heterozygosity(y_loci))
      d <- ape::dist.gene(y_loci, method = 'pairwise', pairwise.deletion = TRUE) %>% as.vector()
      d_vec <- bin_vec(d, sort=TRUE,max=TRUE) %>% as_tibble_row()
      colnames(d_vec) <- paste0('d',1:10)
      max_d <- max(d)
      avg_d <- mean(d)
      pct75 <- as.numeric(quantile(d)[4])
      d_abvpct75 <- length(which(d>pct75))/length(d)
      h1 <- hillR::hill_taxa(round(d,digits = 2),q=1)
      h4 <- hillR::hill_taxa(round(d,digits = 2),q=4)
      hdiff <- h1-h4
      
      y_dna <- vcfR2DNAbin(y)
      sfs_vec <- as.vector(pegas::site.spectrum(y_dna)) %>% bin_vec()
      sfs <- sum(sapply(1:(length(sfs_vec)-1),function(x){return(sfs_vec[x]-sfs_vec[x+1])}))
      sfs_vec <- sfs_vec %>% as_tibble_row()
      colnames(sfs_vec) <- paste0('sfs',1:10)
      
      ss <- tibble(hs,max_d,avg_d,pct75,d_abvpct75,h1,h4,hdiff,sfs)
      ss <- tibble(ss,sfs_vec,d_vec)
      return(ss)
    }
  }
  message('Calculating stats...')
  if (class(x)=='list') {
    if (ncores!=0) {
      tic()
      stats <- mclapply(x,calcStats,mc.cores = ncores) %>% do.call(what = rbind.data.frame)
      toc()
    } else {
      tic()
      stats <- mclapply(x,calcStats) %>% do.call(what = rbind.data.frame)
      toc()
    }
  } else {
    tic()
    stats <- calcStats(x)
    toc()
  }
  return(stats)
}


### demPlot ###
### Plotting parameters with stats
###

demPlot <- function(dataSS, params, stats, file) {
  require(ggthemes)
  param_labs <- c('Annual Mortality','Age of first reproduction','Generation Time','Adult Longevity','Number of Mates')
  names(param_labs) <- c('am','ar','gt','l','n_mates')
  
  stats_labs <- c("Effective Size","Heterozygosity","Differentiation","SFS Shape")
  names(stats_labs) <- stats
  
  plotPerParam <- function(y) {
    dataSS[[y]] %>%
      pivot_longer(y,
                   names_to = 'param', values_to = 'param_value') %>%
      pivot_longer(stats,
                   names_to = 'stat', values_to = 'stat_value') %>%
      mutate(stat = fct_relevel(stat,stats)) %>%
      ggplot(aes(x=param_value,y=stat_value,group=param_value))+
      geom_boxplot()+
      ggh4x::facet_grid2(param~stat,scales='free',independent = 'all',
                         labeller = labeller(param = param_labs, stat = stats_labs))+
      theme(strip.text = element_text(size = 16))+
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=16,face="bold"))+
      xlab('Parameter')+ylab('Summary Statistic')
    ggsave(paste0(file,'_',y,'.png'),width = 10,height = 3)
           
  }
  temp <- lapply(params,plotPerParam)
}
  
## Plotting for am

  


### demInfer ###


