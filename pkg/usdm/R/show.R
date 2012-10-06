# Author: Babak Naimi, naimi.b@gmail.com
# Date :  Sep. 2012
# Version 1.0
# Licence GPL v3

setMethod ('show' , 'VIF',
           function ( object ) {
             if (length(object@excluded) > 0) {
               cat (length(object@excluded),'variables from the',length(object@variables), 'input variables have collinearity problem:','\n','\n')
               cat (object@excluded,'\n')
             } else cat ('No variable from the',length(object@variables), 'input variables has collinearity problem.','\n')
             cat('\n')
             if (length(object@excluded) > 0) cat('excluding the collinear variables:','\n')
             mx <- .minCor(object@corMatrix)
             cat ('min correlation (',mx[1],'~',mx[2],'): ',object@corMatrix[mx[1],mx[2]], '\n')
             mx <- .maxCor(object@corMatrix)
             cat ('max correlation (',mx[1],'~',mx[2],'): ',object@corMatrix[mx[1],mx[2]], '\n')
             cat ('\n')
             cat('---------- VIFs of the remained variables --------','\n')
             print(object@results)
           }
           )

setMethod ('show' , 'speciesLISA', 
           function(object) {
             cat('class                             :' , class(object), '\n')
             cat('LISA statistic                    :' , object@statistic, '\n\n')
             cat('number of species observations    : ' , nrow(object@LISAs), '\n')
             cat('number of predictor variables     : ' , ncol(object@LISAs), '\n')
             cat('min, mean, max of aggregated LISA : ' , min(object@LISA),',' ,mean(object@LISA),',',max(object@LISA), '\n')
           }
)
