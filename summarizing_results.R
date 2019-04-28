################################################################################
###################   sumarizing results of the simulation study    ############
###################               (for revision)                    ############
################################################################################

library(ggplot2)
library(reshape2)
library(gridExtra)

################ PART 1: Boxplots (2 data blocks) ################################

# I. summurizing data; condition: I:20, J1:40, J2:10

# Sim_1 0.5% noise and 30% zero ######
#(note: load data by hand)
tucker_result_Sim1 <- cbind(RESULT_BenchmarCV[,1],
                       RESULT_rdCV[, 1],
                       RESULT_BIC[,1],
                       RESULT_IS[,1],
                       RESULT_BoLasso[,1], 
                       RESULT_StabS[,1], 
                       "0.5% noise, 30% zeros")
colnames(tucker_result_Sim1) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection", "condition")

PL_Sim1 <- cbind(RESULT_BenchmarCV[,2],
            RESULT_rdCV[, 2],
            RESULT_BIC[,2],
            RESULT_IS[,2],
            RESULT_BoLasso[,2],
            RESULT_StabS[,2], 
            "0.5% noise, 30% zeros")
colnames(PL_Sim1) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection", "condition")
save(file = "sim1.RData", tucker_result_Sim1, PL_Sim1)
#######

# Sim_2 0.5% noise and 50% zero #####################
#(note: load data by hand)
tucker_result_Sim2 <- cbind(RESULT_BenchmarCV[,1],
                            RESULT_rdCV[, 1],
                            RESULT_BIC[,1],
                            RESULT_IS[,1],
                            RESULT_BoLasso[,1], 
                            RESULT_StabS[,1], 
                            "0.5% noise, 50% zeros")
colnames(tucker_result_Sim2) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection", "condition")

PL_Sim2 <- cbind(RESULT_BenchmarCV[,2],
                 RESULT_rdCV[, 2],
                 RESULT_BIC[,2],
                 RESULT_IS[,2],
                 RESULT_BoLasso[,2],
                 RESULT_StabS[,2], 
                 "0.5% noise, 50% zeros")
colnames(PL_Sim2) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection", "condition")
save(file = "sim2.RData", tucker_result_Sim2, PL_Sim2)
###########

# Sim_3 30% noise and 30% zero #####################
#(note: load data by hand)
tucker_result_Sim3 <- cbind(RESULT_BenchmarCV[,1],
                            RESULT_rdCV[, 1],
                            RESULT_BIC[,1],
                            RESULT_IS[,1],
                            RESULT_BoLasso[,1], 
                            RESULT_StabS[,1], 
                            "30% noise, 30% zeros")
colnames(tucker_result_Sim3) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection", "condition")

PL_Sim3 <- cbind(RESULT_BenchmarCV[,2],
                 RESULT_rdCV[, 2],
                 RESULT_BIC[,2],
                 RESULT_IS[,2],
                 RESULT_BoLasso[,2],
                 RESULT_StabS[,2], 
                 "30% noise, 30% zeros")
colnames(PL_Sim3) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection", "condition")
save(file = "sim3.RData", tucker_result_Sim3, PL_Sim3)
###########

# Sim_4 30% noise and 50% zero #####################
#(note: load data by hand)
tucker_result_Sim4 <- cbind(RESULT_BenchmarCV[,1],
                            RESULT_rdCV[, 1],
                            RESULT_BIC[,1],
                            RESULT_IS[,1],
                            RESULT_BoLasso[,1], 
                            RESULT_StabS[,1], 
                            "30% noise, 50% zeros")
colnames(tucker_result_Sim4) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection", "condition")

PL_Sim4 <- cbind(RESULT_BenchmarCV[,2],
                 RESULT_rdCV[, 2],
                 RESULT_BIC[,2],
                 RESULT_IS[,2],
                 RESULT_BoLasso[,2],
                 RESULT_StabS[,2], 
                 "30% noise, 50% zeros")
colnames(PL_Sim4) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection", "condition")
save(file = "sim4.RData", tucker_result_Sim4, PL_Sim4)
###########

load("sim1.RData")
load("sim2.RData")
load("sim3.RData")
load("sim4.RData")

PL <- rbind(PL_Sim1, PL_Sim2, PL_Sim3, PL_Sim4)
PL_final<- data.frame(apply(PL[, 1:6], 2, as.numeric))
PL_final$condition <- PL[, 7]
colnames(PL_final)[c(5, 6)] <- c("BL", "SS")
dat_temp <- melt(PL_final,id.vars="condition", measure.vars=c("CV", "RdCV", "BIC", "IS", "BL", "SS"))


p <- ggplot(dat_temp, aes(x = variable, y = value)) +
  geom_boxplot()+
  scale_y_continuous(name = "Proportion of loadings correctedly selected") +
  scale_x_discrete(name = "Variable selection methods") +
  ggtitle("I=20, J1=40, J2=10") +    #do not forget to manually change this.
  theme_bw() +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 14, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 12))+
  facet_grid(. ~ condition)
p

Tucker_results <- rbind(tucker_result_Sim1, tucker_result_Sim2, tucker_result_Sim3, tucker_result_Sim4)
Tucker_final<- data.frame(apply(Tucker_results[, 1:6], 2, as.numeric))
Tucker_final$condition <- Tucker_results[, 7]
colnames(Tucker_final)[c(5, 6)] <- c("BL", "SS")
dat_temp <- melt(Tucker_final,id.vars="condition", measure.vars=c("CV", "RdCV", "BIC", "IS", "BL", "SS"))

p_tucker <- ggplot(dat_temp, aes(x = variable, y = value)) +
  geom_boxplot()+
  scale_y_continuous(name = "Tucker congruence") +
  scale_x_discrete(name = "Variable selection methods") +
  ggtitle("I=20, J1=40, J2=10") +    #do not forget to manually change this.
  theme_bw() +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 14, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 12))+
  facet_grid(. ~ condition)
p_tucker
