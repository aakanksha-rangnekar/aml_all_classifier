library(tidyverse)
library(caret)

###### Data Loading and cleanup ---------------------------------------------------
source("src/load_data.R")

#### Data EDA --------------------------------------------------------------------
train %>% select("Gene Accession Number") %>% unique() ## Will act as our feature list 
train %>% select("Gene Description") %>% unique()

all_data <- train %>% full_join(independent, by = c("Gene Accession Number", "Gene Description") ) 
gene_description = all_data %>% select(c("Gene Accession Number", "Gene Description"))
all_data <- all_data %>% select(-c("Gene Description")) %>% 
  as.data.frame("Gene Accession Number") %>% 
  column_to_rownames("Gene Accession Number") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('patient') %>% 
  as_tibble() %>% 
  mutate(patient = as.numeric(patient)) %>% 
  full_join(actual, by =  'patient')   %>% 
  mutate(patient = ifelse(as.numeric(patient) %in% seq(1:9), paste0("Patient_0", patient), paste0("Patient_", patient) ))
all_data_train <- all_data[1:38,]
all_data_test <-all_data[39:72,]


###### Summarise the data -----------------------------------------------------------
all_data_train %>% 
  summary()
bar_plot_all_vs_aml <- actual %>% 
  ggplot() + 
  geom_bar(aes(cancer)) + 
  ggtitle('AML and ALL samples')

#ggsave('~/R/Bioinformatics/aml_all_classifier/plots/AML_vs_ALL_samples.jpg', bar_plot_all_vs_aml)

actual %>% 
  count(cancer) %>% 
  rowwise() %>%
  mutate(pct = n* 100 /nrow(actual))

#### Datasplitting ----------------------------------------------------- 
set.seed(100)
trainIndex <- createDataPartition(all_data_train$cancer, p = .75, 
                                  list = FALSE, 
                                  times = 1)  # Train and test datasets wil be split in 75, 25 respectively  

testIndex = setdiff(all_data_train$patient %>% str_replace("Patient_0|Patient_", ""), trainIndex)
train_data <- all_data_train %>% filter(str_replace(patient, "Patient_0|Patient_", "")  %in%trainIndex) # Train data index
test_data <- all_data_train %>% filter(str_replace(patient, "Patient_0|Patient_", "") %in% testIndex) # Test data index

train_data  %>% count(cancer) %>% 
  rowwise() %>%
  mutate(pct = n* 100 /nrow(train_data))

test_data  %>% count(cancer) %>% 
  rowwise() %>%
  mutate(pct = n* 100 /nrow(test_data))

response_train_data = train_data$cancer
response_test_data = test_data$cancer


train_data <- train_data %>% 
  select(-cancer) %>% 
  column_to_rownames('patient')


test_data <- test_data %>% 
  select(-cancer) %>% 
  column_to_rownames('patient')

summary_train_data <- t(train_data) %>% summary()


ihs <- function(x) {
  y <- log(x + sqrt(x^2 + 1))
  return(y)
}


cancer_data = all_data_train %>% select(c(patient, cancer))
box_plot_before_normalization  <- t(train_data) %>%  # EDA before normalization
  as.data.frame() %>%
  rownames_to_column('Genes') %>%
  gather(key = patient, value =value, -Genes) %>%
  mutate(ihs_count = ihs(value), 
         patient = (patient)) %>%
  as_tibble() %>%
  left_join(cancer_data, by = "patient") %>%
  group_by(patient, cancer) %>%
  ggplot() +
  aes((patient), ihs_count, color = cancer) +
  geom_boxplot() +
  ggtitle('Before normalization') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6))
ggsave("~//R/Bioinformatics/aml_all_classifier/plots/box_plots_before_normalization.jpeg",  box_plot_before_normalization)

summary_data_train_minimum <- summary_train_data[1,] %>% str_replace_all("Min.   :", "") %>% str_replace_all(" ", "") %>% as.numeric()
summary_data_train_minimum <- summary_train_data[1,] %>% str_replace_all("Min.   :", "") %>% str_replace_all(" ", "") %>% as.numeric()

summary_test_data <- test_data %>% summary()
#### Normalization ----------------------------------------------------- 

#normalized_train_data <-  apply(train_data, 2,  function(x){x/rowSums(train_data)})
#normalized_test_data <-  apply(test_data, 2,  function(x){x/rowSums(test_data)})  

box_plot_after_normalization  <- t(normalized_train_data) %>%
  as.data.frame() %>%
  rownames_to_column('Genes') %>%
  gather(key = Sample, value = value, -Genes) %>%
  ggplot() +
  aes(Sample, (value)) +
  geom_boxplot() +
  ggtitle('After library size normalization')

norm_functions  <- c("corr", 'center','scale', 'pca')

transform_mat <- preProcess(as.matrix(train_data), norm_functions, thresh = 0.85) 
norm_train_data  <- predict(transform_mat, as.matrix(train_data))
norm_test_data  <- predict(transform_mat, as.matrix(test_data))
box_plot_after_preprocess_train_data  <-t(norm_train_data) %>% 
  as.data.frame() %>%
  rownames_to_column('Genes') %>%
  gather(key = patient, value =value, -Genes) %>%
  mutate(ihs_count = ihs(value), 
         patient = (patient)) %>%
  as_tibble() %>%
  left_join(cancer_data, by = "patient") %>%
  group_by(patient, cancer) %>%
  ggplot() +
  aes((patient), ihs_count, color = cancer) +
  geom_boxplot() +
  ggtitle('After normalization and dimensionality reduction') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6))

ggsave("~/R/Bioinformatics/aml_all_classifier/plots/box_plots_after_normalization.jpeg",  box_plot_after_normalization)

fitControl <- trainControl(## 5-fold CV
  method = "repeatedcv",
  number = 3,
  ## repeated ten times
  repeats = 10)


metric = "Accuracy"
set.seed(100)
svmFit1 <- train(norm_train_data, response_train_data, 
                 method = "svmLinear", 
                 trControl = fitControl,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = FALSE,
                 metric = metric)


set.seed(100)
knnfit<-train(norm_train_data,response_train_data,method="knn", trControl = fitControl,metric=metric)
set.seed(100)
rffit<-train(norm_train_data,response_train_data,method="rf",trControl=fitControl,metric=metric,
             verbose=F)  



result<-resamples(list(rf=rffit,svm=svmFit1,knn=knnfit))

jpeg("~/R/Bioinformatics/aml_all_classifier/plots/dot_plot.jpg")
dotplot(result,main="Model Accuracy Results")
dev.off()

pred_test <- predict(rffit, norm_test_data)  
test_confusion_matrix  <- confusionMatrix(pred_test, as.factor(response_test_data))


table <- data.frame(test_confusion_matrix$table)

plotTable <- table %>%
  mutate(goodbad = ifelse(table$Prediction == table$Reference, "good", "bad")) %>%
  mutate(Reference = factor(Reference, levels = c("AML", "ALL")), 
         Prediction = factor(Prediction, levels = c("AML", "ALL"))) %>% 
  group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq)) 




# fill alpha relative to sensitivity/specificity by proportional outcomes within reference groups (see dplyr code above as well as original confusion matrix for comparison)
test_confusion_plot <- ggplot(data = plotTable, mapping = aes(x = Reference, y = Prediction, fill = goodbad, alpha = prop)) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1) +
  scale_fill_manual(values = c(good = "green", bad = "red")) +
  theme_bw() +
  xlim(rev(levels(table$Reference)))


ggsave("~/R/Bioinformatics/aml_all_classifier/plots/test_confusion_matrix_plot.jpg",test_confusion_plot)


validation_data <- all_data_test %>% select(-cancer) %>% column_to_rownames("patient") 
validation_norm <- predict(transform_mat, as.matrix(validation_data))
validation_pred <- predict(rffit, as.matrix(validation_norm))
confusion_val <- confusionMatrix(validation_pred, as.factor(all_data_test$cancer)) 



########### Plotting confusion matrix 
library(ggplot2)
library(dplyr)

table <- data.frame(confusion_val$table)

plotTable <- table %>%
  mutate(goodbad = ifelse(table$Prediction == table$Reference, "good", "bad")) %>%
  group_by(Reference) %>%
  mutate(prop = Freq/sum(Freq))
validation_confusion_matrix <- ggplot(data = plotTable, mapping = aes(x = Reference, y = Prediction, fill = goodbad, alpha = prop)) +
  geom_tile() +
  geom_text(aes(label = Freq), vjust = .5, fontface  = "bold", alpha = 1) +
  scale_fill_manual(values = c(good = "green", bad = "red")) +
  theme_bw() +
  xlim(rev(levels(table$Reference)))
ggsave("~/R/Bioinformatics/aml_all_classifier/plots/validation_confusion_matrix_plot.jpg",validation_confusion_matrix)


# library(caret)
# library(mlbench)
# data(Sonar)
# ctrl <- trainControl(method="cv", 
#                      summaryFunction=twoClassSummary, 
#                      classProbs=T,
#                      savePredictions = T)
# rfFit <- train(Class ~ ., data=Sonar, 
#                method="rf", preProc=c("center", "scale"), 
#                trControl=ctrl)

library(pROC)
# Select a parameter setting
selectedIndices <- rfFit$pred$mtry == 2
# Plot:
plot.roc((rfFit$pred$obs[selectedIndices]),
         rfFit$pred$M[selectedIndices])


AML<- all_data_train %>% select(c(1:5, "cancer")) %>% filter(cancer == "AML") 
ALL<- all_data_train %>% select(c(1:5, "cancer")) %>% filter(cancer == "ALL") 

subset_data <- bind_rows(AML, ALL)
write_csv(all_data_train, "data/all_data_train.csv")


pred_test_svm <- predict(svmFit1, norm_test_data)  
test_confusion_matrix_svm  <- confusionMatrix(pred_test_svm, as.factor(response_test_data))

pred_val_svm <- predict(svmFit1, validation_norm)  
test_confusion_matrix_svm  <- confusionMatrix(pred_val_svm, as.factor(all_data_test$cancer))

