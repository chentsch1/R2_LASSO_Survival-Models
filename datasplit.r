
################################################################################
########### Function: Random Data Splitting for Training and Test ##############
################################################################################


# data = a data set, which should be splitted
# ratio = splitting ratio
# seed = optional seed


datasplit <- function(data, ratio=0.33333, seed=NA) {
	n <- nrow(data)
	testsize <- (n*ratio)
	rownumb <- seq(1:nrow(data))
	data <- as.data.frame(cbind(rownumb,data))
	setClass("folders", representation = representation(train="data.frame", 
		test="data.frame"))
	datafold <- new("folders")
	if(!is.na(seed)){set.seed(seed)}
    s <- sample(data$rownumb,testsize)
		o <- sort(s)
		datafold@test <- data[o,2:ncol(data)]
		datafold@train <- data[-o,2:ncol(data)]
	return(datafold)}


