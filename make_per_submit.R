#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# Makes HTCondor submit script for running permutations
# in parallel

# Usage: make_per_submit(drive,
#   executable,
#   arguments, (as a list)
#   output_dir,
#   req_cores,
#   req_mem,
#   req_disk,
#   num_jobs,
#   submit_file)

# Example: make_per_submit('/mnt/lustre1/users/lazar/GIBBONS',
#   '/gibbon_meth/wrap_par_rand.R',
#   c('$(dir)/VOK_GENOME/par_permute.dat', 'all', '1000', '$$(Cpus)'),
#   '/mnt/lustre1/users/lazar/GIBBONS/VOK_GENOME', 16,  '2 GB,'
#   '2 GB,', 63, 'condor.submit')

make_per_submit <- function(drive, 
  executable, arguments=c(), output_dir='.',
  req_cores=16, req_mem='2 GB', req_disk='2 GB',
  num_jobs=1, submit_file) {

  log_dir <- paste0(output_dir, '/logs')
  if (!file.exists(log_dir))
      dir.create(log_dir)

writeLines(c('ID=$(Cluster).$(Process)',
  paste0('dir=', drive),
  'should_transfer_files = IF_NEEDED',
  'when_to_transfer_output = ON_EXIT',
  'Requirements= (machine!="exanode-0-4.local")',
  paste0('executable=', executable),
  paste0('arguments=', paste(arguments, collapse=' ')),
  paste0('output=', output_dir, '/permute.$(Process).txt'),
  paste0('error=', log_dir, '/permute.stderr.$(ID)'),
  paste0('log=', log_dir, '/permute.log.$(ID)'), 
  paste0('request_cpus=', as.character(req_cores)),
  paste0('request_memory=', req_mem),
  paste0('request_disk=', req_disk),
  paste('queue', as.character(num_jobs))), 
submit_file)
}


