

###################################################################
###Preprocess: normalization, thesholding, and log2 transformation
################################################################### 
hem.preproc <- function (x, cond, data.type = "MAS5")
{

    #Refine!!!
    #Exclude cases with missing values
    x <- as.matrix(na.exclude(x))
 

    #IQR normalization
    if (data.type == "MAS4" || data.type == "MAS5") {
        x <- quant.normalize(x, percent = 50)
    }

    #Tresholding
    if (data.type == "MAS4" || data.type == "dChip") {
        if (length(x[x < 1]) != 0) {
            x[x < 1] <- 1
        }
    }
    if (data.type == "MAS5") {
        if (length(x[x < 0.1]) != 0) {
            x[x < 0.1] <- 0.1
        }
    }
   
    #log2 transformation
    x <- logb(x, 2)

    #LOWESS normalization with rank invariant genes
#    x0 <- dat.rank.invar(x, cond=cond)
#    for(j in 2:ncol(x)){
#       x[,j] <- lowess.normalize.rank.invar(x[,1], x[,j], x0[,1], x0[,j])$y 
#    }

    #LOWESS normalization
#    for(j in 2:ncol(x)){
#       x[,j] <- lowess.normalize(x[,1], x[,j])$y 
#    }


    return(x)
}


###IQR normaization
quant.normalize <- function (x, percent = 50)
{
    quartile <- apply(x, 2, quant.norm, percent = percent)
    max.quartile <- max(quartile)
    ratio <- (quartile/max.quartile)
    ratio.vect <- rep(ratio, nrow(x))
    adjusted <- matrix(ratio.vect, nrow = nrow(x), byrow = TRUE)
    normalized <- data.frame(x/adjusted)
    return(normalized)
}

quant.norm <- function (x, percent = 50)
{
    low <- 0.5 * (100 - percent)/100
    high <- 0.5 * (100 + percent)/100
    difference <- as.vector(diff(quantile(x, probs = c(low, high),
        na.rm = TRUE)))
    return(difference)
}
