library(MVN)
library(XML)

# Import the data
data <- xmlParse("samples_4_BF.xml")
xml.data <- xmlToList(data)

N <- as.numeric(xml.data$.attrs[["N"]])
beta <- as.numeric(xml.data$.attrs[["beta"]])
n.success <- as.numeric(xml.data$.attrs[["nSuccess"]])
n.failure <- as.numeric(xml.data$.attrs[["nFailure"]])
samples.success <- matrix(unlist(lapply(xml.data$Success, function(entry) as.numeric(unlist(strsplit(entry, "\\s+"))))), ncol = 4 * N, byrow = TRUE)
samples.failure <- matrix(unlist(lapply(xml.data$Failure, function(entry) as.numeric(unlist(strsplit(entry, "\\s+"))))), ncol = 4 * N, byrow = TRUE)
samples.complete <- rbind(samples.success, samples.failure)

# MVN Tests (Korkmaz, Selcuk, Dincer Goksuluk, and Gokmen Zararsiz. "MVN: An R Package for Assessing Multivariate Normality." A peer-reviewed, open-access publication of the R Foundation for Statistical Computing (2014): 151.)
## Mardia's test
pdf(paste("qqplotMardia_", N, "_BF.pdf", sep = ""))
result.failure.mardia <- mardiaTest(samples.failure, qqplot = TRUE)
dev.off()
#result.complete.mardia <- mardiaTest(samples.complete, qqplot = TRUE) # For comparison

## Henze-Zirkler's MVN test
pdf(paste("qqplotHz_", N, "_BF.pdf", sep = ""))
result.failure.hz <- hzTest(samples.failure, qqplot = TRUE)
dev.off()
#result.complete.hz <- hzTest(samples.complete, qqplot = TRUE)  # For comparison

## Royston's MVN test
if (n.success >= 2000 && n.failure >= 2000) {
    pdf(paste("qqplotRoyston_", N, "_BF.pdf"))
    result.failure.royston <- roystonTest(samples.failure[1:2000,], qqplot = TRUE)
    dev.off()
    #result.complete.royston <- roystonTest(samples.complete[1:2000, ], qqplot = TRUE)  # For comparison
}
pdf(paste("qqplotRoyston_", N, "_BF.pdf"))
result.failure.royston <- roystonTest(samples.failure, qqplot = TRUE)
dev.off()

# Log-likelihood test on the parameters based on Wilks's theorem
paramTest <- function(data, meanTarget, covTarget, sl = 0.05) {
    n <- nrow(data)
    p <- ncol(data)
    
    mu <- colMeans(data) # ML estimation of the mean
    Sigma <- cov(data) * (n - 1) / n # ML estimation of the covariance
    
    data.minusMeanTarget <- data - matrix(rep(1, n) %x% meanTarget, ncol = p, byrow = TRUE) # Subtract meanTarget from all samples
    W <- -n * (log(det(Sigma)) - log(det(covTarget))) - n * p + sum(diag(solve(covTarget, t(data.minusMeanTarget) %*% data.minusMeanTarget)))# The test statistic -2 * log(Lambda) 
    
    result <- list();
    result$name <- "Log-likelihood mean and covariance test"
    result$p.value <- 1 - pchisq(W, p * (p + 1) / 2 + p)
    if (result$p.value < sl) result$Result <- "Mean and covariance of the samples do not match the target." # Reject
    else result$Result <- "Mean and covariance of the samples match the target." # Accept
    
    return (result)
}

meanTarget <- -rep(0, 4 * N)
covTarget <- beta / 2 * diag(4 * N)
result.failure.param <- paramTest(samples.failure, meanTarget, covTarget)
#result.complete.param <- paramTest(samples.complete, meanTarget, covTarget)

# Output
cat("N =", N, ", n.failure =", n.failure, ", n.success =", n.success, ", beta =", beta, "\n")
print("MVN distribution tests ********************************************")
print(result.failure.mardia)
print(result.failure.hz)
print(result.failure.royston)
print("Parameter tests ***************************************************")
print(result.failure.param)
