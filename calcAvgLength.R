args=commandArgs(trailingOnly=T)
lengthVector = c()
for (i in 2:length(args)){
	tryCatch((
		print(paste("Reading file number " , i-1))
		data=read.table(args[i], header=T)
		lengthVector = c(lengthVector,data[4])
	), error=function(e){})
}
paste("The avergage length/width of all peaks is:", mean(as.numeric(unlist(lengthVector))) )
#"The avergage length/width of all peaks is: 396.190491215771"
