`as.GST` <-
function(GSD,GSDo){
              #if(GSDo$z<GSD$b[GSDo[1][[1]]])print("Stopping rule not met")
              #else{
                GST<-list(GSD=GSD,GSDo=GSDo)
                class(GST)<-"GSTobj"
                return(GST)
                #}
              }

