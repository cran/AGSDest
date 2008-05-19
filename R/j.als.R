`j.als` <-
function(u,GSD=NULL){
           #if(is.null(GSD)) stop("missing input: specify GSD");
           if(is.null(GSD$als))  GSD$als <- comp.als(GSD=GSD);
           sum(ifelse(u>GSD$als[1:length(GSD$als)],1,0))+1
         }

