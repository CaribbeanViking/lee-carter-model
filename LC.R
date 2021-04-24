## MORTALITY RATES DATA PREPARATION ##
ageLimit <- 90
ageGroups <- ageLimit + 1
timePeriods <- 98
tHor <- 52
mrRaw <- read.csv(file.choose(), header=FALSE) 							# Import dataset
colnames(mrRaw) <- c("Year", "Age", "Females", "Males", "Combined") 	# Name columns
mrReducedAge <- subset(mrRaw, Age >= 0 & Age <= ageLimit)				# Reduce to age interval
mrVec <- t(mrReducedAge[-c(1:4)])										# Remove all columns except "Combined"
mrMatrix <- matrix(mrVec, nrow = ageGroups, ncol = timePeriods)			# Fold into matrix
colnames(mrMatrix) <- 1921:2018											# Assign appropriate column names
rownames(mrMatrix) <- 0:ageLimit 										# Assign appropriate row names
mode(mrMatrix) = "numeric"												# Make matrix numeric

## BASIC LEE-CARTER METHOD ##
## Fitting procedure ##
a_x = rowMeans(log(mrMatrix))											# Create vector of means
M = log(mrMatrix) - a_x													# Define matrix M (centralize)
svdM <- svd(M)															# Perform SVD on M (sing.vec. = 1)
b_x <- -svdM$u[,1]														# Define b_x parameter
k_t <- -svdM$d[1]*svdM$v[,1]											# Define k_t parameter
Mhat <- matrix(nrow = ageGroups, ncol = timePeriods)					# Initiate M-hat matrix
for (i in 1:timePeriods){												# Fit data
	Mhat[,i] <- a_x + b_x*k_t[i]
}
colnames(Mhat) <- 1921:2018												# Appropriate naming
rownames(Mhat) <- 0:ageLimit

## Forecasting procedure ##
dhat <- (k_t[timePeriods] - k_t[1])/(timePeriods - 1)					# Define drift parameter
total = matrix(nrow = timePeriods - 1, ncol = 1)						# Initiate temp vector
for (i in 1:timePeriods-1) {											# Calculate see
	total[i] <- (k_t[i+1] - k_t[i] - dhat)^2
}
see = sqrt(1/(timePeriods-2)*sum(total))								#
kForc_t <- matrix(nrow = tHor, ncol = 1)								# Initiate forcast k_t vector
for (i in 1:tHor) {														# Calculate forecasted k_t
	kForc_t[i] = k_t[timePeriods] + dhat*i + sqrt(i)*rnorm(1,0,see^2)	# Random walk
}
Mtilde <- matrix(nrow = ageGroups, ncol = tHor)							# Initiate M-tilde matrix
for (i in 1:tHor) {														# Fit data
	Mtilde[,i] <- a_x + b_x*kForc_t[i]
}
colnames(Mtilde) <- 2019:2070											# Appropriate naming
rownames(Mtilde) <- 0:ageLimit
Mstar <- cbind(Mhat,Mtilde)												# Concatenate to complete model

## Proportion variance explained procedure ##
num <- matrix(nrow = timePeriods, ncol = 1)								# Initiate temp vectors
den <- matrix(nrow = timePeriods, ncol = 1)
etaSqrd <- matrix(nrow = ageGroups, ncol = 1)							# Initiate eta^2 vector
for (i in 1:ageGroups){													# Calculate prop. variance expl.
	for (j in 1:timePeriods) {												
		num[j] <- (mrMatrix[i,j] - exp(Mhat[i,j]))^2
		den[j] <- (mrMatrix[i,j] - exp(a_x[i]))^2
	}
etaSqrd[i] = 1 - (sum(num)/sum(den))
}

## Confidence intervals procedure ##
upper <- matrix(nrow = ageGroups, ncol = tHor)							# Initiate interval vectors
lower <- matrix(nrow = ageGroups, ncol = tHor)

for (i in 1:ageGroups) {	
	for (j in 1:tHor) {													# Create confidence intervals
	upper[i,j] <- Mtilde[i,j]*exp(-1.96*b_x[i]*see)
	lower[i,j] <- Mtilde[i,j]*exp(1.96*b_x[i]*see)
	}
}