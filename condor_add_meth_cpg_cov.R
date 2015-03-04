# nathan dot lazar at gmail dot com

# Makes files and condor script to run add_meth_cpg_cov in parallel

condor_add_meth_cpg_cov <- function(feat.gr, all.bs, outdir) {

  # If drive doesn't exist, create it
  if(!file.exists(outdir)) dir.create(outdir)

  #Save data needed by parallel jobs
  save(feat.gr, all.bs, file=paste0(outdir, '/feat_gr_and_all_bs.dat'))

  # Make the submit file for condor
  writeLines(c('ID=$(Cluster).$(Process)',
    paste0('dir=', outdir),
    'should_transfer_files = IF_NEEDED',
    'when_to_transfer_output = ON_EXIT',
    'executable=/mnt/lustre1/users/lazar/GIBBONS/gibbon_meth/wrap_add_meth_cpg_cov.R',
    'arguments=$(dir)/feat_gr_and_all_bs.dat $$(Cpus) $(dir)/',
    'output=$(dir)/add.$(Process).txt',
    'error=$(dir)/add.stderr.$(ID)',
    'log=$(dir)/add.log.$(ID)',
    'request_cpus=24',
    'request_memory=2 GB',
    'request_disk=2 GB',
    'queue 1'),
    paste0(outdir, '/add.submit'))

  # Run HTCondor script
  system(paste0('condor_submit ', outdir, '/add.submit'))

  # Pause until the condor script is done
  while(!file.exists(paste0(outdir, '/add.0.txt')))
    Sys.sleep(1)
  while(length(readLines(paste0(outdir, '/add.0.txt')))==0)
    Sys.sleep(5)

  # Read in feat.gr (which is saved by wrap_add_meth_cpg_cov.R)
  load(paste0(outdir, '/feat.dat'))

  feat.gr
}
  
