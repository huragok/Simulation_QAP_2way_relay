library(MVN)
library(XML)

# Import the data
data <- xmlParse("samples_4.xml")
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
pdf(paste("qqplotMardia_", N, ".pdf", sep = ""))
result.failure.mardia <- mardiaTest(samples.failure, qqplot = TRUE)
dev.off()
#result.complete.mardia <- mardiaTest(samples.complete, qqplot = TRUE) # For comparison

## Henze-Zirkler's MVN test
pdf(paste("qqplotHz_", N, ".pdf", sep = ""))
result.failure.hz <- hzTest(samples.failure, qqplot = TRUE)
dev.off()
#result.complete.hz <- hzTest(samples.complete, qqplot = TRUE)  # For comparison

## Royston's MVN test
if (n.success >= 2000 && n.failure >= 2000) {
    pdf(paste("qqplotRoyston_", N, ".pdf"))
    result.failure.royston <- roystonTest(samples.failure[1:2000,], qqplot = TRUE)
    dev.off()
    #result.complete.royston <- roystonTest(samples.complete[1:2000, ], qqplot = TRUE)  # For comparison
}

# Log-likelihood test on the parameters based on Wilks's theorem
## Likelyhood ratio test to test the covaiance matrix under unknown mean and MVN assumption (Pinto, Letícia Pereira, and Sueli Aparecida Mingoti. "ON HYPOTHESIS TESTS FOR COVARIANCE MATRICES UNDER MULTIVARIATE NORMALITY." Pesquisa Operacional 35.1 (2015): 123-142.)
covTest <- function(data, covTarget, sl = 0.05) {
    n <- nrow(data)
    p <- ncol(data)
    
    mu <- colMeans(data)
    A <- (n - 1) * cov(data)
    
    matTmp <- solve(covTarget, A)
    W <- -p * n + p * n * log(n) - n * log(det(matTmp)) + sum(diag(matTmp)) # The test statistic 
    
    result <- list();
    result$name <- "Log-likelihood covariance test"
    result$p.value <- 1 - pchisq(W, p * (p + 1) / 2)
    if (result$p.value < sl) result$Result <- "Covariance matrix of the samples does not match the target covariance matrix." # Reject
    else result$Result <- "Covariance matrix of the samples matches the target covariance matrix." # Accept
    
    return (result)
}

covTarget <- beta / 2 * diag(4 * N)
result.failure.cov <- covTest(samples.failure, covTarget)
#result.complete.cov <- covTest(samples.complete, covTarget)

## Mean test under unknown covariance matrix and MVN assumption (Ruey S. Tsay. "Inference about sample mean", lecture notes, Business 41912, The University of Chicago Booth School of Business, Spring Quarter 2010.)
meanTest <- function(data, meanTarget, sl = 0.05) {
    n <- nrow(data)
    p <- ncol(data)
    
    SigmaEst <- cov(data)
    SigmaTarget <- t(data) %*% data / (n - 1)
    W <- -n * log(det(SigmaEst) / det(SigmaTarget))
    
    result <- list();
    result$name <- "Log-likelihood mean test"
    result$p.value <- 1 - pchisq(W, p)
    if (result$p.value < sl) result$Result <- "Mean vector of the samples does not match the target mean vector." # Reject
    else result$Result <- "Mean vector of the samples matches the target mean vector." # Accept
    
    return (result)
}

meanTarget <- rep(0, 4 * N)
result.failure.mean <- meanTest(samples.failure, meanTarget)
#result.complete.mean <- meanTest(samples.complete, meanTarget)

# Output
cat("N =", N, ", n.failure =", n.failure, ", n.success =", n.success, ", beta =", beta, "\n")
print("MVN distribution tests ********************************************")
print(result.failure.mardia)
print(result.failure.hz)
print(result.failure.royston)
print("Parameter tests ***************************************************")
print(result.failure.cov)
print(result.failure.mean)
