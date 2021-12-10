library(boot)
library("rjson")
library("purrr")
library("dplyr")
library("stringr")
library("ggplot2")
library("stats")

createInformationMatrix <- function(dataset)
{
  productInformation <- data.frame(Product_ID=integer(),
                                   Title=character(), 
                                   WebShop=character(),
                                   modelID=character(),
                                   stringsAsFactors=FALSE) 
  
  k <- 1
  for(i in 1:length(dataset))
  {
    for(j in 1:length(dataset[[i]]))
    {
      productInformation <- rbind(productInformation, c(k, dataset[[i]][[j]]$title, dataset[[i]][[j]]$shop, dataset[[i]][[j]]$modelID))
      k <- k+1
    }
  }
  colnames(productInformation) <- c("Product ID", "Title", "Webshop", "ModelID")
  return(productInformation)
}
cleanData <- function(productInformation)
{
  productInformation[,2]<-gsub('"','inch',productInformation$Title)
  productInformation[,2]<-gsub('Inch','inch',productInformation$Title)
  productInformation[,2]<-gsub('inches','inch',productInformation$Title)
  productInformation[,2]<-gsub('-inch','inch',productInformation$Title)
  productInformation[,2]<-gsub(' inch','inch',productInformation$Title)
  
  productInformation[,2]<-gsub('Hertz','hz',productInformation$Title)
  productInformation[,2]<-gsub('hertz','hz',productInformation$Title)
  productInformation[,2]<-gsub('Hz','hz',productInformation$Title)
  productInformation[,2]<-gsub('HZ','hz',productInformation$Title)
  productInformation[,2]<-gsub('hz','hz',productInformation$Title)
  productInformation[,2]<-gsub('-hz','hz',productInformation$Title)
  
  productInformation[,2] <- tolower(productInformation[,2])
  productInformation[,2]<-gsub('[^a-z0-9 ]',"",productInformation$Title)
  
  productInformation[,2]<-gsub('  ',' ',productInformation$Title)
  
  return(productInformation)
}
calculateDuplicates <- function(productInformation)
{
  size <- nrow(productInformation)
  duplicates <- matrix(0 ,nrow = size, ncol = size)
  
  for(i in 1:(size-1))
  {
    for(j in (i+1):size)
    {
      if(productInformation[i,4] == productInformation[j,4])
      {
        duplicates[i,j] = 1
      }
    }
  }
  return(duplicates)
}
createModelWords <- function(productInformation)
{
  modelWords <- data.frame(Product_ID=character(), 
                           stringsAsFactors=FALSE)
  for(i in 1:nrow(productInformation))
  {
    temp <- unlist(strsplit(productInformation[i,2], " "))
    for(j in temp)
    {
      if(grepl("([0-9].*[a-z])", j))
      {
        modelWords <- rbind(modelWords, j)
      }
    }
  }
  colnames(modelWords) <- c("Model Word")
  modelWords <- modelWords %>% distinct()
  return(modelWords)
}
createBinaryVector <- function(productInformation, modelWords)
{
  binaryVector <- matrix(0 ,nrow = nrow(modelWords), ncol = nrow(productInformation))
  for(i in 1:nrow(productInformation))
  {
    for(j in 1:nrow(modelWords))
    {
      if(grepl(modelWords[j,1], productInformation[i,2], fixed = TRUE))
      {
        binaryVector[j,i] = 1
      }
    }
  }
  return(binaryVector)
}
createPermutations <- function(aantalMinhashes, aantalModelWords)
{
  permutations <- matrix(0, nrow = aantalMinhashes, ncol = aantalModelWords)
  for(i in 1:aantalMinhashes)
  {
    permutations[i,] <- sample(aantalModelWords)
  }
  return(permutations)
}
createSignatureMatrix <- function(permutations, binaryVector)
{
  signatureMatrix <- matrix(0, nrow = nrow(permutations), ncol = ncol(binaryVector))
  for(i in 1:nrow(permutations))
  {
    for(k in 1:ncol(binaryVector))
    {
      for(j in 1:ncol(permutations))
      {
        if(binaryVector[match(j,permutations[i,]), k] == 1) 
        {
          signatureMatrix[i,k] = j
          break
        }
      }
    }
  }
  return(signatureMatrix)
}
createPotentialDuplicates <- function(signatureMatrix, numberOfBands, numberOfRows)
{
  potentialDuplicates <- matrix(0 ,nrow = ncol(signatureMatrix), ncol = ncol(signatureMatrix))
  
  for(i in 1:numberOfBands)
  {
    tempBands <- matrix(0 ,nrow = 0, ncol = 2)
    for(j in 1:ncol(signatureMatrix))
    {
      tempSigna <- toString(signatureMatrix[(((i-1)*numberOfRows) + 1):(i*numberOfRows),j])
      if(tempSigna %in% tempBands[,1])
      {
        index <- match(tempSigna, tempBands[,1])
        tempBands[index,2] <- paste(tempBands[index,2], toString(j))       
      }
      else
      {
        tempBands <- rbind(tempBands, c(tempSigna, toString(j)))
      }
    }
    for(k in 1:nrow(tempBands))
    {
      tempDuplicates <- unlist(strsplit(tempBands[k,2], " "))
      if(length(tempDuplicates) >= 2)
      {
        for(m in 1:(length(tempDuplicates)-1))
        {
          for(n in (m+1):length(tempDuplicates))
          {
            potentialDuplicates[as.numeric(tempDuplicates[m]), as.numeric(tempDuplicates[n])] <- 1
          }
        }
      }
    }
  }
  return(potentialDuplicates)
}
createEvaluation <- function(duplicates, bootstrapPotentialDuplicates, bootstrapLSH)
{
  numberOfProducts <- nrow(bootstrapPotentialDuplicates)
  candidatePairs <- 0
  for(i in 1:numberOfProducts)
  {
    candidatePairs <- candidatePairs + sum(bootstrapPotentialDuplicates[i,])
  }
  
  TP <- 0
  FN <- 0
  FP <- 0
  
  for(i in 1:numberOfProducts)
  {
    for(j in 1:numberOfProducts)
    {
      if(bootstrapPotentialDuplicates[i,j] == 1 & duplicates[bootstrapLSH[i],bootstrapLSH[j]] == 1)
      {
        TP <- TP + 1
      }
      else if(bootstrapPotentialDuplicates[i,j] == 1)
      {
        FP <- FP + 1
      }
      else if(duplicates[bootstrapLSH[i],bootstrapLSH[j]] == 1)
      {
        FN <- FN + 1
      }
    }
  }
  
  pairQuality <- TP/candidatePairs
  pairCompleteness <- TP/bootstraptAantalDuplicates
  
  f1Measure <- (2*pairQuality*pairCompleteness)/(pairQuality + pairCompleteness)
  fractionOfComparisons <- candidatePairs/((numberOfProducts*(numberOfProducts-1))/2)
  return(c(pairCompleteness, pairQuality, f1Measure, fractionOfComparisons))
}
createBands <- function(numberOfMinhashes)
{
  totalBands <- vector()
  for(i in 1:minhashes)
  {
    if(minhashes%%i == 0)
    {
      totalBands <- rbind(totalBands,i)
    }
  }
  return(totalBands)
}
createModelWordsKeyValuePairs <- function(dataset)
{
  productInformationKeyValuePairs <- vector()
  
  for(i in 1:length(dataset))
  {
    for(j in 1:length(dataset[[i]]))
    {
      temp <- ""
      for(k in 1:length(dataset[[i]][[j]]$featuresMap))
      {
        temp <- paste(temp, dataset[[i]][[j]]$featuresMap[[k]], sep = " ")
      }
      productInformationKeyValuePairs <- rbind(productInformationKeyValuePairs, temp)
    }
  }
  colnames(productInformationKeyValuePairs) <- c("Value Model Words")
  return(productInformationKeyValuePairs)
}
cleanDataKeyValuePairs <- function(modelWordsKeyValuePairs)
{
  modelWordsKeyValuePairs<-gsub('"','inch',modelWordsKeyValuePairs)
  modelWordsKeyValuePairs<-gsub('Inch','inch',modelWordsKeyValuePairs)
  modelWordsKeyValuePairs<-gsub('inches','inch',modelWordsKeyValuePairs)
  modelWordsKeyValuePairs<-gsub('-inch','inch',modelWordsKeyValuePairs)
  modelWordsKeyValuePairs<-gsub(' inch','inch',modelWordsKeyValuePairs)
  
  modelWordsKeyValuePairs<-gsub('Hertz','hz',modelWordsKeyValuePairs)
  modelWordsKeyValuePairs<-gsub('hertz','hz',modelWordsKeyValuePairs)
  modelWordsKeyValuePairs<-gsub('Hz','hz',modelWordsKeyValuePairs)
  modelWordsKeyValuePairs<-gsub('HZ','hz',modelWordsKeyValuePairs)
  modelWordsKeyValuePairs<-gsub('hz','hz',modelWordsKeyValuePairs)
  modelWordsKeyValuePairs<-gsub('-hz','hz',modelWordsKeyValuePairs)
  
  modelWordsKeyValuePairs <- tolower(modelWordsKeyValuePairs)
  modelWordsKeyValuePairs<-gsub('[^a-z0-9 ]',"",modelWordsKeyValuePairs)
  modelWordsKeyValuePairs<-gsub("[a-z]","",modelWordsKeyValuePairs)
  modelWordsKeyValuePairs<-gsub('\\b\\w{1}\\b',"",modelWordsKeyValuePairs)
  modelWordsKeyValuePairs<-gsub('\\b\\w{9,}\\b',"",modelWordsKeyValuePairs)
  
  modelWordsKeyValuePairs<-gsub('\\s+',' ',modelWordsKeyValuePairs, perl=T)
}
createKeyValueModelWords <- function(productInformationKeyValuePairs)
{
  modelWordsKeyValue <- data.frame(Product_ID=character(), 
                           stringsAsFactors=FALSE)
  for(i in 1:nrow(productInformationKeyValuePairs))
  {
    temp <- unlist(strsplit(productInformationKeyValuePairs[i], " "))
    for(j in temp)
    {
      if(grepl("[0-9]", j))
      {
        modelWordsKeyValue <- rbind(modelWordsKeyValue, j)
      }
    }
  }
  colnames(modelWordsKeyValue) <- c("Model Word")
  modelWordsKeyValue <- modelWordsKeyValue %>% distinct()
  return(modelWordsKeyValue)
}
createKeyValueBinaryVector <- function(productInformationKeyValuePairs, modelWordsKeyValuePairs)
{
  binaryVectorKeyValue <- matrix(0 ,nrow = nrow(modelWordsKeyValuePairs), ncol = nrow(productInformationKeyValuePairs))
  for(i in 1:nrow(productInformationKeyValuePairs))
  {
    for(j in 1:nrow(modelWordsKeyValuePairs))
    {
      if(grepl(modelWordsKeyValuePairs[j,1], productInformationKeyValuePairs[i,1], fixed = TRUE))
      {
        binaryVectorKeyValue[j,i] = 1
      }
    }
  }
  return(binaryVectorKeyValue)
}
Jaccard <- function (x, y) 
{
  M.11 <- sum(x == 1 & y == 1)
  M.10 <- sum(x == 1 & y == 0)
  M.01 <- sum(x == 0 & y == 1)
  return (M.11 / (M.11 + M.10 + M.01))
}
dataset <- fromJSON(file = "C:\\Users\\denni\\Desktop\\Master\\Blok 2 - Computer Science for Business Analytics\\Paper\\TVs-all-merged.json")

productInformation <- createInformationMatrix(dataset)
productInformation <- cleanData(productInformation)
duplicates <- calculateDuplicates(productInformation)
aantalDuplicates <- length(productInformation[,1]) - length(dataset)
modelWords <- createModelWords(productInformation)
aantalModelWords <- nrow(modelWords)
binaryVector <- createBinaryVector(productInformation, modelWords)
minhashes <- 480
permutations <- createPermutations(minhashes, nrow(modelWords))
signatureMatrix <- createSignatureMatrix(permutations, binaryVector)
totalBands <- createBands(minhashes)
totalRows <- minhashes/totalBands
productInformationKeyValuePairs <- createModelWordsKeyValuePairs(dataset)
productInformationKeyValuePairs <- cleanDataKeyValuePairs(productInformationKeyValuePairs)
modelWordsKeyValuePairs <- createKeyValueModelWords(productInformationKeyValuePairs)
binaryVectorKeyValue <- createKeyValueBinaryVector(productInformationKeyValuePairs, modelWordsKeyValuePairs)

evaluationLSH <- matrix(0, nrow=0, ncol=4)
evaluationMSM <- matrix(0, nrow=0, ncol=4)
colnames(evaluationLSH) <- c("pairCompleteness", "pairQuality", "f1Measure", "fractionOfComparisons")
colnames(evaluationMSM) <- c("pairCompleteness", "pairQuality", "f1Measure", "fractionOfComparisons")
numberOfBootstraps <- 5

for(i in 1:numberOfBootstraps)
{
  set.seed(i)
  bootstrap <- sample(nrow(productInformation), replace=T)
  
  bootstraptAantalDuplicates <- 0
  for(j in 1:(length(bootstrap)-1))
  {
    for(k in (j+1):length(bootstrap))
    {
      if(duplicates[bootstrap[j],bootstrap[k]] == 1 || duplicates[bootstrap[k],bootstrap[j]] == 1)
      {
        bootstraptAantalDuplicates <- bootstraptAantalDuplicates + 1
      }
    }
  }
  
  bootstrapSignatureMatrix <- matrix(0, nrow= nrow(signatureMatrix), ncol=length(bootstrap))
  for(j in 1:length(bootstrap))
  {
    bootstrapSignatureMatrix[,j] <- signatureMatrix[,bootstrap[j]]
  }
  
  for(b in 1:length(totalBands))
  {
    bootstrapPotentialDuplicates <- createPotentialDuplicates(bootstrapSignatureMatrix, totalBands[b], totalRows[b])
    
    numberOfProducts <- nrow(bootstrapPotentialDuplicates)
    candidatePairs <- 0
    for(i in 1:numberOfProducts)
    {
      candidatePairs <- candidatePairs + sum(bootstrapPotentialDuplicates[i,])
    }
    
    TP <- 0
    FN <- 0
    FP <- 0
    
    for(i in 1:numberOfProducts)
    {
      for(j in 1:numberOfProducts)
      {
        if(bootstrapPotentialDuplicates[i,j] == 1 & duplicates[bootstrap[i],bootstrap[j]] == 1)
        {
          TP <- TP + 1
        }
        else if(bootstrapPotentialDuplicates[i,j] == 1)
        {
          FP <- FP + 1
        }
        else if(duplicates[bootstrap[i],bootstrap[j]] == 1)
        {
          FN <- FN + 1
        }
      }
    }
    
    pairQuality <- TP/candidatePairs
    pairCompleteness <- TP/bootstraptAantalDuplicates
    
    f1Measure <- (2*pairQuality*pairCompleteness)/(pairQuality + pairCompleteness)
    fractionOfComparisons <- candidatePairs/((numberOfProducts*(numberOfProducts-1))/2)
    
    evaluationLSH <- rbind(evaluationLSH, c(pairCompleteness, pairQuality, f1Measure, fractionOfComparisons))
    
    dissimilarityMatrix = matrix(data = Inf, nrow = nrow(bootstrapPotentialDuplicates), ncol = ncol(bootstrapPotentialDuplicates))
    
    for(i in 1:nrow(bootstrapPotentialDuplicates))
    {
      for(j in 1:ncol(bootstrapPotentialDuplicates))
      {
        indexi <- bootstrap[i]
        indexj <- bootstrap[j]
        if(bootstrapPotentialDuplicates[i,j] == 1 && productInformation[indexi,3] != productInformation[indexj,3])
        {
          titleWords <- Jaccard(binaryVector[,indexi], binaryVector[,indexj])
          keyValueWords <- Jaccard(binaryVectorKeyValue[,indexi], binaryVectorKeyValue[,indexj])
          
          dissimilarityMatrix[i,j] = 1- 2/((1/titleWords)+(1/keyValueWords))
        }
      }
    }
    dissimilarityMatrix <- replace(dissimilarityMatrix, is.na(dissimilarityMatrix), 1)
    
    duplicatesMSM <- matrix(0, nrow = nrow(dissimilarityMatrix), ncol=ncol(dissimilarityMatrix))
    for(i in 1:nrow(dissimilarityMatrix))
    {
      for(j in 1:ncol(dissimilarityMatrix))
      {
        if(dissimilarityMatrix[i,j]<=0.8)
        {
          duplicatesMSM[bootstrap[i],bootstrap[j]] <- 1
        }
      }
    }
    
    TPMSM <- 0
    FNMSM <- 0
    FPMSM <- 0
    som <- 0
    
    for(i in 1:numberOfProducts)
    {
      for(j in 1:numberOfProducts)
      {
        if(duplicates[i,j] == 1)
        {
          som <- som + 1
        }
        if(duplicatesMSM[i,j] && (duplicates[i,j] == 1 || duplicates[j,i] == 1))
        {
          TPMSM <- TPMSM + 1
        }
        else if(duplicatesMSM[i,j] == 1)
        {
          FPMSM <- FPMSM + 1
        }
        else if(duplicates[i,j] == 1)
        {
          FNMSM <- FNMSM + 1
        }
      }
    }
    
    pairQualityMSM <- TPMSM/candidatePairs
    pairCompletenessMSM <- TPMSM/bootstraptAantalDuplicates
    
    f1MeasureMSM <- (2*pairQualityMSM*pairCompletenessMSM)/(pairQualityMSM + pairCompletenessMSM)
    f1starMeasureMSM <- TPMSM / (TPMSM + 0.5*(FPMSM+FNMSM))
    fractionOfComparisonsMSM <- candidatePairs/((numberOfProducts*(numberOfProducts-1))/2)
    
    evaluationMSM <- rbind(evaluationMSM, c(pairCompletenessMSM, pairQualityMSM, f1starMeasureMSM, fractionOfComparisonsMSM))
  }
}