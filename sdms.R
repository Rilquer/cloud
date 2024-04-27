addCustomOcc <- function(orig,custom) {
  if (is.data.frame(orig)) {
    orig <- list(orig)
  }
  names <- sapply(1:length(orig),function(x){a=x$scientific_name[1]
  
  add <- function(y) {
    orig[[y]] <- orig[[y]] %>% add_row(custom %>% filter(binomial))
    last_ID <- occs[['Sclerurus_cearensis']]$occID[which(is.na(occs[['Sclerurus_cearensis']]$occID))[1]-1]
    final <- length(which(is.na(occs[['Sclerurus_cearensis']]$occID)))+last_ID
    occs[['Sclerurus_cearensis']] <- occs[['Sclerurus_cearensis']] %>%
      mutate(occID = c(na.omit(occID),seq(last_ID+1,final)))
    
  lapply(names,function(x){})
  
  
  }
  
  /