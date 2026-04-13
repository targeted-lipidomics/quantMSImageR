combine_rds = function(rds_names){

  combined = readRDS(rds_names[1])

  for(ind in 2:length(rds_names)){
    combined = combine_MSIs(combined, readRDS(rds_names[ind]))
  }

  return(combined)
}
