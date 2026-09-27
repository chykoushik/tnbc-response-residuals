for (n in c('renv_snapshot_preflight','renv_available_packages','renv_download_method')) {cat('\n',n,'\n');f<-get0(n,envir=asNamespace('renv'));if(!is.null(f))print(f)}
