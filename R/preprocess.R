

###################################################################
###Preprocess: normalization, thesholding, and log2 transformation
################################################################### 
hem.preproc <- function (x, data.type = "MAS5")
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
