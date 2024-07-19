####  Recapitation - reticulate   ####

recap <- function(x) {
  # Modified from SLiM manual
  msprime <- import("msprime")
  tskit <- import('tskit')
  pyslim <- import('pyslim')
  ts <- tskit$load(paste0(x$outpath,x$outfile,'.trees'))
  if (x$ev_type=='none') {
    rts <- pyslim$recapitate(ts, ancestral_Ne = as.numeric(x$Ne_recap), recombination_rate = 0)  
  }
  if (x$ev_type=='split') {
    demography <- msprime$Demography$from_tree_sequence(ts)
    for (pop in demography$populations) {
      # must set their effective population sizes
      pop$initial_size = as.numeric(x$Ne_recap)
    }
    demography$add_migration_rate_change(time=ts$metadata$SLiM$tick,rate=as.numeric(x$mig_recap),
                                         source="p1", dest="p2")
    demography$add_migration_rate_change(time=ts$metadata$SLiM$tick,rate=as.numeric(x$mig_recap),
                                         source="p2", dest="p1")
    rts <- pyslim$recapitate(ts, demography = demography, recombination_rate = 0)
  }
  
  rts <- msprime$sim_mutations(rts, rate = as.numeric(x$mut_rate))
  return(rts)
}

piCalc <- function(rts) {
  tskit <- import('tskit')
  return(rts$diversity())
}

tajDCalc <- function(rts) {
  tskit <- import('tskit')
  return(rts$Tajimas_D())
}

segCalc <- function(rts) {
  tskit <- import('tskit')
  return(rts$segregating_sites())
}

afsCalc <- function(rts) {
  tskit <- import('tskit')
  return(rts$allele_frequency_spectrum())
}