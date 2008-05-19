`j.alab` <-
function(h,GSD=NULL){
           if(is.null(GSD)) stop("missing input: specify GSD");
           if(is.null(GSD$alab))  GSD$alab  <-comp.alab(GSD=GSD);
           length(GSD$alab[1:(length(GSD$alab)-3)])+1-sum(ifelse(h>=GSD$alab[1:(length(GSD$alab)-3)],1,0))
         }

