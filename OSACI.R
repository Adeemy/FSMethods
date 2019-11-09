##################################################################################################
## Feature Selection Using OSACI Algorithms: OSACI-Init, OSACI-Update, and OSACI-Init_Update
##
## Paper: "New feature selection methods based on opposition-based learning and 
##         self-adaptive cohort intelligence for predicting patient no-shows"
## Authors: Mohammed Aladeemy, Linda Adwan, Amy Booth, Mohammad T. Khasawneh, and Srikanth Poranki
## Corresponding author: Mohammed Aladeemy (mohammed.aladeemy@gmail.com)
## DOI: https://doi.org/10.1016/j.asoc.2019.105866
##################################################################################################

## Load required packages
library(pROC)
library(caret)
library(doParallel)

# Training and testing sets (Source: UCI repository
# (https://archive.ics.uci.edu/ml/datasets/Parkinson+Dataset+with+replicated+acoustic+features+)
Train <- read.table("Parkinson Train Set.txt", header = TRUE)
Test  <- read.table("Parkinson Test Set.txt", header = TRUE)

# Split features and class of training and testing sets
TrainClass    <- Train[ ,ncol(Train)]  # Training class
TrainFeatures <- Train[ ,-ncol(Train)] # Training set without class

TestClass    <- Test[ ,ncol(Test)]    # Testing class
TestFeatures <- Test[ ,-ncol(Test)]   # Testing set without class

# Global training parameters for "train" function 
TrainControl <- trainControl(method          = "repeatedcv",    
                             number          = 5,               
                             repeats         = 3,		            
                             classProbs      = TRUE,
                             allowParallel   = TRUE,
                             search          = "grid",
                             verboseIter     = FALSE,
                             returnData      = FALSE,
                             seeds           = NULL,
                             summaryFunction = twoClassSummary
)

# Fitness function
ObjFun <- function(Solutions, d, lambda) {
  
  ObjFunValue <- data.frame(ObjFunValue = rep(NA, (nrow(Solutions))))
  Solutions   <- foreach(j = seq_len(nrow(Solutions)), .inorder = FALSE, .combine = 'rbind',
                         .verbose = FALSE) %dopar% {
                           
                           # Repair any solution with no selected features by selecting one feature randomly
                           if (sum(Solutions[j, ]) == 0) {
                             Solutions[j, sample(seq_len(d), size = 1)] <- 1
                           }  
                           
                           # Training set with selected feature subset
                           SelF <- as.data.frame(TrainFeatures[ ,which(Solutions[j, ] == 1)])
                           
                           # Train Decision Tree (DT) on training set with selected feature subset
                           DT_model <- train(SelF, TrainClass, method = "rpart", metric = "ROC", 
                                             parms = list(split = "information"),
                                             tuneLength = 10, trControl = TrainControl)
                           
                           # Fitness value (behavior)
                           ObjFunValue[j, ] <- ((lambda * max(DT_model$results$ROC)) + ((1 - lambda) * (d - ncol(SelF))/d))
                           return(cbind(t(Solutions[j, ]), ObjFunValue[j, ]))
                         }
  return(Solutions)
}

######################################################################################
## OSACI algorithms

OSACI <- function(OSACIVariant, L, S, d, Q, k, EPSILON, TauMax, lambda) {
  
  ## Check input parameters
  
  if (S < 2 || missing(S)) {
    stop("No. of candidates must be greater than 2")
  }
  
  if ((S %% 2) != 0) {
    stop("No. of candidates must be even")
  }
  
  if (Q < 2 || missing(Q)) {
    warning("No. of quality variations must be greater than 2")
  }
  
  if (L < 1 || missing(L)) {
    stop("The max. no. of learning attempts must be greater than 0")
  }
  
  if (is.null(TauMax) == TRUE || missing(TauMax)) {
    TauMax <- Inf
  }
  
  if (is.null(EPSILON) == TRUE || missing(EPSILON)) {
    EPSILON <- 0
  }
  
  # Allocate matrices to save "for" loops outputs
  Candidates <- matrix(NA, nrow = S, ncol = d)
  BestObjVal <- matrix(NA, nrow = L, ncol = 2)
  Results    <- matrix(NA, nrow = L, ncol = d + 1)
  Behaviors  <- matrix(NA, nrow = Q * S, ncol = d)
  
  ## Generate initial candidates
  
  # OBL initialization (Algorithm 1)
  if (OSACIVariant %in% c("OSACI-Init", "OSACI-Init_Update")) { 
    
    # Generate S/2 random Candidates
    for (i in seq_len(S/2)) {
      Candidates[i, ] <- round(runif(d, 0, 1))
    }
    
    # Genereate opposite candidates of the randomly generated (S/2) solutions
    for (j in seq_len(S/2)) {
      Candidates[(S/2) + j, ] <- 1 - Candidates[j, ]
    }
  }
  
  # Random initialization
  if (OSACIVariant %in% c("OSACI-Update")) { 
    
    for (i in seq_len(S)) {
      Candidates[i, ] <- round(runif(d, 0, 1))
    }
  }
  
  # Evaluate initial Candidates  
  Candidates <- ObjFun(Candidates, d, lambda)
  
  # Initiate a counters
  Tau <- 0
  l   <- 1
  
  ## Procedure
  repeat { # While l < L
    
    ## Update the cohort using OBL Update (Algorithm 2)
    
    if (OSACIVariant %in% c("OSACI-Update", "OSACI-Init_Update")) {
      
      # Sort candidates in descending order based on their fitness values
      Candidates <- Candidates[order(Candidates[ ,(d + 1)], decreasing = TRUE), ]
      
      # Allocate memory for S/2 opposite candidates
      OBLCandidates <- Candidates[seq_len(S/2), -(d + 1)]
      
      # Allocate memory for S/2 mutant candidates
      MutantCandidates <- Candidates[seq(from = (S/2) + 1, to = S), -(d + 1)]
      
      for(j in seq_len(S)) {   # Generate S/2 Type-I opposite candidates 
        if (j <= S/2){
          OBLCandidates[j, ] <- 1 - OBLCandidates[j, ]
        }
        
        if (j > S/2) {   # Generate S/2 mutant candidates 
          jj <- sample(seq_len(d), size = (1/S) * d)
          MutantCandidates[j-(S/2), jj] <- 1 - MutantCandidates[j-(S/2), jj]
        }
      }
      
      # Evaluate opposite candidates
      OBLCandidates <- ObjFun(OBLCandidates, d, lambda)
      
      # Select the fittest S/2 candidates from S/2 original candidates and their opposites
      OBLCandidates              <- rbind(Candidates[seq_len(S/2), ], OBLCandidates)
      Candidates[seq_len(S/2), ] <- OBLCandidates[order(OBLCandidates[ ,ncol(OBLCandidates)], decreasing = TRUE), ][seq_len(S/2), ]
      
      # Evaluate mutant candidates
      Candidates[seq(from = (S/2) + 1, to = S), ] <- ObjFun(MutantCandidates, d, lambda)
    }
    
    ## Update the cohort in OSACI-Init
    
    if (OSACIVariant %in% c("OSACI-Init")) {
      
      # Sort candidates in ascending order based on their fitness values to keep the fittest candidate in last row
      Candidates <- Candidates[order(Candidates[ ,(d + 1)], decreasing = FALSE), ]
      
      # Allocate memory for a set of S - 1 mutant candidates
      MutantCandidates <- Candidates[seq_len(S - 1), -(d + 1)]
      
      # Select S - 1 candidates using tournament selection
      for (s in seq_len(S - 1)) {
        randsel               <- sample(seq_len(S), size = k)
        MutantCandidates[s, ] <- Candidates[randsel[which.max(Candidates[randsel, (d + 1)])], -(d + 1)]
      }
      
      # Mutation with rate of 1/S
      for (j in seq_len(S - 1)) {
        jj                      <- sample(seq_len(d), size = (1/S) * d)
        MutantCandidates[j, jj] <- 1 - MutantCandidates[j, jj]
      }
      
      # Evaluate mutant candidates
      Candidates[-S, ] <- ObjFun(MutantCandidates, d, lambda)
    }
    
    # Allocate memory for selected candidates
    SelectedCandidates <- Candidates[order(Candidates[ ,(d + 1)], decreasing = FALSE), ]
    
    ## Each candidate selects a candidate using tournament selection
    for (s in seq_len(S)) {
      randsel                 <- sample(seq_len(S), size = k)
      SelectedCandidates[s, ] <- Candidates[randsel[which.max(Candidates[randsel, (d + 1)])], ]
    }
    
    # Each candidate samples Q solutions (qualities) from its neighborhood
    for (s in seq_len(S)) {
      
      lo <- (s - 1) * Q
      
      # Update mutation rate
      ms_l <- max(1, ceiling(d * (SelectedCandidates[which.max(SelectedCandidates[ ,(d + 1)]) ,(d + 1)] - SelectedCandidates[s ,(d + 1)])))
      
      for (kk in seq_len(Q)) {
        jj                        <- sample(seq_len(d), size = ms_l, replace = FALSE)
        Behaviors[(lo + kk), jj]  <- 1 - SelectedCandidates[s, jj]
        Behaviors[(lo + kk), -jj] <- SelectedCandidates[s, -c(jj, (d + 1))]
      }
    }
    
    # Evaluate solutions sampled from each candidate's neighborhood
    Behaviors <- ObjFun(Behaviors, d, lambda)
    
    # Update current candidates
    kk <- 0
    for (s in seq_len(S)) {
      lo <- (s - 1) * Q
      kk <- kk + Q
      
      # The fittest solution sampled from neighborhood of candidate s
      BestSampled_Q <- Behaviors[which.max(Behaviors[(lo + 1):kk, (d + 1)]) + lo, ]

      Candidates[s, ] <- BestSampled_Q
      
      # if (s < S) {
      #   Candidates[s, ] <- BestSampled_Q
      # } else {
      #   Candidates[s, BestSampled_Q[(d + 1)] >= Candidates[s, (d + 1)]] <- BestSampled_Q
      # }
    }
    
    # Resize Behaviors matrix
    Behaviors <- Behaviors[ ,-ncol(Behaviors)]
    
    ## Save the best solution obtained in the current learning attempt
    
    # Max. and min. objective values to evaluate cohort saturation
    BestObjVal[l, 1] <- Candidates[which.max(Candidates[ ,(d + 1)]), (d + 1)]
    BestObjVal[l, 2] <- Candidates[which.min(Candidates[ ,(d + 1)]), (d + 1)]
    
    # Save the fittest solution in the current learning attempt
    if (l == 1) {
      Results[l, ] <- Candidates[which.max(Candidates[ ,(d + 1)]), ]
    } 
    
    if (Results[(l - 1), (d + 1)] <= BestObjVal[l, 1] && l > 1) {
      Results[l, ] <- Candidates[which.max(Candidates[ ,(d + 1)]), ]
    } 
    
    if (Results[(l - 1), (d + 1)] > BestObjVal[l, 1] && l > 1) {
      Results[l, ] <- Results[(l - 1), ]
    }
    
    ## Cohort saturation reached?
    
    # Update Tau
    if (l > 1 && abs(BestObjVal[l, 1] - BestObjVal[(l - 1), 1]) <= EPSILON) {
      Tau <- Tau + 1
    } else {
      Tau <- 0
    }
    
    ## Stopping criteria
    
    # Terminate if cohort saturation is reached
    if (Tau == TauMax) break()
    
    # Terminate if max. no. of learning attempts is reached
    if (l == L) break()         
    
    # Update counter 
    l <- l + 1
    
  } # End of "repeat" loop
  
  return(Results)
} # End of OSACI function
######################################################################################

## Create a cluster for parallel computing

Numclust <- detectCores(logical = TRUE) - 1

clust <- makeCluster(Numclust, type = "SOCK")

registerDoParallel(clust)

varlist <- ls(envir = parent.frame(), all.names = TRUE)

parallel::clusterExport(clust, varlist = varlist, envir = parent.frame())

invisible(lapply(.packages(), function(pkg)
  parallel::clusterCall(clust, library, package = pkg, character.only = TRUE))[1])

# Start timer for CPU time
StartTime <- proc.time()
Results   <- OSACI(OSACIVariant = "OSACI-Init",        # "OSACI-Init", "OSACI-Update", or "OSACI-Init_Update"
                   L = 50,                             # Max. no. of learning attempts
                   S = 10,                             # No. of candidates
                   d = ncol(TrainFeatures),            # No. of features 
                   Q = 10,                             # No. of quality variations (solutions sampled from neighborhood)
                   k = 2,                              # Tournament size
                   EPSILON = NULL,                     # Convergence tolerence
                   TauMax = NULL,                      # No. of successive learning attempts for cohort saturation
                   lambda = 0.9                        # Trade-off factor in fitness function (lambda*AUC + (1-lambda)*Dim. Reduction)
)

# CPU time
RunTime <- (proc.time() - StartTime)[3]

## Shut down workers of the cluster
stopCluster(clust)

# Best solution
Sol <- Results[which.max(Results[ ,ncol(Train)]), ]

# Selected feature subset
SelectedFeatures <- which(Sol[-c((length(Sol) - 1), (length(Sol)))] == 1)

# Train DT model on training set using selected feature subset
DTmodel <- train(as.data.frame(TrainFeatures[ ,SelectedFeatures]), TrainClass,
                 method = "rpart", metric = "ROC", parms = list(split = "information"),
                 tuneLength = 10, trControl = TrainControl)

## Performance measures on testing set
TestPred <- confusionMatrix(data = predict(DTmodel, TestFeatures[ ,SelectedFeatures]),
                            reference = TestClass, positive = levels(TestClass)[1])

# Sensitivity (true positive rate)
SensTest <- round(TestPred$byClass[1], 3)[1]

# Specificity (true negative rate)
SpecTest <- round(TestPred$byClass[2], 3)[1]

# AUC
ROCTest <- round(roc(as.factor(TestClass), as.numeric(predict(DTmodel, Test[ ,SelectedFeatures])))$auc[1], 3)
