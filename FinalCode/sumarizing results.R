# sumarizing results
library(ggplot2)


# Sim_1 0.5% noise and 30% zero
  #(note: load data by hand)

tucker_result <- cbind(RESULT_BenchmarCV[,1],
                       RESULT_rdCV[, 1],
                       RESULT_BIC[,1],
                       RESULT_IS[,1],
                       RESULT_BoLasso[,1], 
                       RESULT_StabS[,1])
colnames(tucker_result) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stability selection")

boxplot(tucker_result, ylab = "Tucker congruence", main = "0.5% noise and 30% zeros", ylim = c(0,1))

PL <- cbind(RESULT_BenchmarCV[,2],
            RESULT_rdCV[, 2],
            RESULT_BIC[,2],
            RESULT_IS[,2],
            RESULT_BoLasso[,2],
            RESULT_StabS[,2])
colnames(PL) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stability selection")
boxplot(PL, ylab = "Proportion of loadings corrected selected", main = "0.5% noise and 30% zeros", ylim = c(0,1))


# Sim_2 0.5% noise and 50% zero
#(note: load data by hand)

tucker_result <- cbind(RESULT_BenchmarCV[,1],
                       RESULT_rdCV[, 1],
                       RESULT_BIC[,1],
                       RESULT_IS[,1],
                       RESULT_BoLasso[,1],
                       RESULT_StabS[,1])
colnames(tucker_result) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stability selection")

boxplot(tucker_result, ylab = "Tucker congruence", main = "0.5% noise and 50% zeros", ylim = c(0,1))

PL <- cbind(RESULT_BenchmarCV[,2],
            RESULT_rdCV[, 2],
            RESULT_BIC[,2],
            RESULT_IS[,2],
            RESULT_BoLasso[,2],
            RESULT_StabS[,2])
colnames(PL) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stability selection")
boxplot(PL, ylab = "Proportion of loadings corrected selected", main = "0.5% noise and 50% zeros", ylim = c(0,1))

# Sim_3 5% noise and 30% zero
#(note: load data by hand)

tucker_result <- cbind(RESULT_BenchmarCV[,1],
                       RESULT_rdCV[, 1],
                       RESULT_BIC[,1],
                       RESULT_IS[,1],
                       RESULT_BoLasso[,1],
                       RESULT_StabS[,1])
colnames(tucker_result) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stability selection")

boxplot(tucker_result, ylab = "Tucker congruence", main = "5% noise and 30% zeros", ylim = c(0,1))

PL <- cbind(RESULT_BenchmarCV[,2],
            RESULT_rdCV[, 2],
            RESULT_BIC[,2],
            RESULT_IS[,2],
            RESULT_BoLasso[,2],
            RESULT_StabS[,2])
colnames(PL) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stability selection")
boxplot(PL, ylab = "Proportion of loadings corrected selected", main = "5% noise and 30% zeros", ylim = c(0,1))


# Sim_4 5% noise and 50% zero
#(note: load data by hand)

tucker_result <- cbind(RESULT_BenchmarCV[,1],
                       RESULT_rdCV[, 1],
                       RESULT_BIC[,1],
                       RESULT_IS[,1],
                       RESULT_BoLasso[,1],
                       RESULT_StabS[,1])
colnames(tucker_result) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stability selection")

boxplot(tucker_result, ylab = "Tucker congruence", main = "5% noise and 50% zeros", ylim = c(0,1))

PL <- cbind(RESULT_BenchmarCV[,2],
            RESULT_rdCV[, 2],
            RESULT_BIC[,2],
            RESULT_IS[,2],
            RESULT_BoLasso[,2],
            RESULT_StabS[,2])
colnames(PL) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stability selection")
boxplot(PL, ylab = "Proportion of loadings corrected selected", main = "5% noise and 50% zeros", ylim = c(0,1))



# Sim_5 30% noise and 30% zero
#(note: load data by hand)

tucker_result <- cbind(RESULT_BenchmarCV[,1],
                       RESULT_rdCV[, 1],
                       RESULT_BIC[,1],
                       RESULT_IS[,1],
                       RESULT_BoLasso[,1],
                       RESULT_StabS[,1])
colnames(tucker_result) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stability selection")

boxplot(tucker_result, ylab = "Tucker congruence", main = "30% noise and 30% zeros", ylim = c(0,1))

PL <- cbind(RESULT_BenchmarCV[,2],
            RESULT_rdCV[, 2],
            RESULT_BIC[,2],
            RESULT_IS[,2],
            RESULT_BoLasso[,2],
            RESULT_StabS[,2])
colnames(PL) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stability selection")
boxplot(PL, ylab = "Proportion of loadings corrected selected", main = "30% noise and 30% zeros", ylim = c(0,1))


# Sim_6 30% noise and 50% zero
#(note: load data by hand)

tucker_result <- cbind(RESULT_BenchmarCV[,1],
                       RESULT_rdCV[, 1],
                       RESULT_BIC[,1],
                       RESULT_IS[,1],
                       RESULT_BoLasso[,1],
                       RESULT_StabS[,1])
colnames(tucker_result) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stability selection")

boxplot(tucker_result, ylab = "Tucker congruence", main = "30% noise and 50% zeros", ylim = c(0,1))

PL <- cbind(RESULT_BenchmarCV[,2],
            RESULT_rdCV[, 2],
            RESULT_BIC[,2],
            RESULT_IS[,2],
            RESULT_BoLasso[,2],
            RESULT_StabS[,2])
colnames(PL) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stability selection")
boxplot(PL, ylab = "Proportion of loadings corrected selected", main = "30% noise and 50% zeros", ylim = c(0,1))


