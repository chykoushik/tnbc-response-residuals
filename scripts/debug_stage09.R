source('R/paths.R');source('R/checksums.R');source('R/prediction_folds.R');source('R/sc_source_audit.R');source('R/sc_annotation.R');source('R/sc_pathway_scoring.R')
options(error=function(){traceback(5);quit(status=1)})
score09()
