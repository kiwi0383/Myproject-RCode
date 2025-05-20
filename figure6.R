#------------------------------------------------------------------------------
#Construction and evaluation of the ICI response predictive model.
#------------------------------------------------------------------------------

# =============================================================================
# 1. Data Preparation
# =============================================================================
# Set working directory
new_directory <- "E:/project1/data/immunotherapy_cohorts"
setwd(new_directory)

# Load datasets
melanoma1 <- read.xlsx('./1melanoma/melanoma1.xlsx')
melanoma2 <- read.xlsx('./2melanoma/melanoma2.xlsx')
melanoma3 <- read.xlsx('./3melanoma/melanoma3.xlsx')
blca1 <- read.xlsx('./blca/GSE176307/blca-GSE176307.xlsx')
blca2 <- read.xlsx('./blca/blca-imvigor.xlsx')
glioblastoma <- read.xlsx('./Glioblastoma/Glioblastoma.xlsx')
lusc1 <- read.xlsx('./lusc/GSE126044.xlsx')
lusc2 <- read.xlsx('./lusc/GSE135222.xlsx')
lusc3 <- read.xlsx('./lusc/GSE166449.xlsx')
rcc <- read.xlsx('./rcc/rcc.xlsx')
stad <- read.xlsx('./stad.xlsx')

# Combine all datasets
combined_df <- rbindlist(list(melanoma1, melanoma2, melanoma3, blca1, blca2, 
                              glioblastoma, lusc1, lusc2, lusc3, rcc, stad), 
                         fill = TRUE)

# Get common columns across all datasets
common_cols <- Reduce(intersect, list(
  names(melanoma1), names(melanoma2), names(melanoma3),
  names(blca1), names(blca2), names(glioblastoma),
  names(lusc1), names(lusc2), names(lusc3), names(rcc), names(stad)
))

# Subset combined data with common columns
df_all <- combined_df[, common_cols, with = FALSE]

# =============================================================================
# 1.1 Explore relationship between response and risk score
# =============================================================================
# Load coefficients for risk score calculation
coef_df <- read.csv('E:/project1/results/fig1/multicox.csv')
genes <- coef_df$term
coef_values <- coef_df$estimate
names(coef_values) <- coef_df$term

# Calculate risk score
df_all <- df_all %>%
  mutate(
    riskscore = CHEK2 * coef_values[1] + 
      LIPT1 * coef_values[2] + 
      TUFM * coef_values[3] +
      NDUFA10 * coef_values[4] +
      AGK * coef_values[5] + 
      PNPLA2 * coef_values[6] +
      GFM1 * coef_values[7],
    risklevel = ifelse(riskscore > median(riskscore), 'high', 'low')
  )

# =============================================================================
# 1.2 Batch Effect Removal
# =============================================================================
library(sva)

# Convert data types and clean
df_all <- type.convert(df_all)
df_all <- df_all %>%
  filter(across(everything(), ~ .x != 0))

# Recode response variable
df_all$response <- case_when(
  df_all$response == 'N' ~ 'NR',
  TRUE ~ df_all$response
)

# Prepare for ComBat
rownames(df_all) <- df_all$id
batch <- as.factor(df_all$batch)
group <- as.factor(df_all$response)
norm_data <- t(df_all[, 1:7])  # Assuming first 7 columns are gene expression

# Run ComBat
mod <- model.matrix(~as.factor(response), data = df_all) 
expr_combat <- ComBat(
  dat = norm_data, 
  batch = batch,
  mod = mod,
  par.prior = TRUE
)

# Reformat post-ComBat data
df_combat <- expr_combat %>%
  t() %>%
  as.data.frame() %>%
  mutate(
    response = factor(df_all$response),
    batch = df_all$batch
  )

# =============================================================================
# 1.3 Data Splitting
# =============================================================================
# Training cohorts
train_cohorts <- c('GSE126044', 'rcc', 'GSE166449', 
                   'blca-GSE176307', 'melanoma1', 'melanoma3')

# Split data
train <- df_combat[df_combat$batch %in% train_cohorts, ]
validation <- df_combat[!df_combat$batch %in% train_cohorts, ]

# Further split training set
set.seed(0212)
train_index <- sample(nrow(train), 0.7 * nrow(train))
data_train <- train[train_index, ]
data_test <- train[-train_index, ]

# =============================================================================
# 2. Machine Learning Modeling
# =============================================================================
# Define cross-validation parameters
ctrl <- trainControl(
  method = "repeatedcv",
  number = 5,        # 5-fold CV
  repeats = 3,       # 3 repeats
  summaryFunction = twoClassSummary,
  classProbs = TRUE
)

# Model tuning grids
tune_grids <- list(
  nb = expand.grid(
    fL = c(0, 0.5, 1, 1.5, 2.0), 
    usekernel = TRUE, 
    adjust = c(0.5, 0.75, 1, 1.25, 1.5)
  ),
  svmRadial = expand.grid(
    sigma = c(0.0005, 0.001, 0.005, 0.01, 0.05),
    C = c(1, 3, 5, 10, 20)
  ),
  rf = expand.grid(mtry = c(2, 42, 83, 124, 165, 205, 246, 287, 328, 369)),
  kknn = expand.grid(
    kmax = c(5, 7, 9, 11, 13), 
    distance = 2, 
    kernel = 'optimal'
  ),
  adaboost = expand.grid(
    nIter = c(50, 100, 150, 200, 250),
    method = c('Adaboost.M1', 'Real adaboost')
  )
)

# 2.1 Naive Bayes Model
nb_model <- train(
  response ~ .,
  data = data_train,
  method = 'nb',
  metric = "ROC",
  trControl = ctrl,
  tuneGrid = tune_grids$nb
)

model.tune <- train(response ~ .,
                    data = data.train,
                    method = 'nb',
                    metric="ROC",
                    trControl=ctrl)

nb.pred  <- predict(model.tune,newdata=data.test,type = "prob")
roc.nb=roc(data.test$response,nb.pred[,2],smooth=F)
plot(roc.nb,print.auc=TRUE)
a=confusionMatrix(as.factor(nb.pred),data.test$response,positive = 'R') 

results_df1=data.frame()
results_df2=data.frame()
results_df1=data.frame(nb=a[["byClass"]])
results_df2=data.frame(nb=a[["overall"]])

#在独立测试集上
nb.pred  <- predict(model.tune,newdata=validation[,1:7],type = "prob")
threshould<-0.5
class<-ifelse(nb.pred[,2] >threshould,1,0)
print(table(validation$response,class))

roc.nb.v=roc(validation$response,nb.pred[,2],smooth=F)
plot(roc.nb.v,print.auc=TRUE)

#2.2 svmRadial
model.tune <- train(response ~ .,
                    data = data.train,
                    method = 'svmRadial',
                    metric="ROC",
                    trControl=ctrl)

svmRadial.pred <- predict(model.tune,data.test,type="prob")
roc.svmRadial=roc(data.test$response,svmRadial.pred[,2],smooth=F)
plot(roc.svmRadial,print.auc=TRUE)

svmRadial.pred <- predict(model.tune,validation,type="prob")
roc.svmRadial.v=roc(validation$response,svmRadial.pred[,2],smooth=F)
plot(roc.svmRadial.v,print.auc=TRUE)

svmRadial.pred <- predict(model.tune,data.test,type="raw")
a=confusionMatrix(as.factor(svmRadial.pred),data.test$response,positive = 'R') 
results_df$svmRadial=a[["byClass"]]
results_df2$svmRadial=a[["overall"]]


#2.3 svmLinear
model.tune <- train(response ~ .,
                    data = data.train,
                    method = 'svmLinear',
                    metric="ROC",
                    trControl=ctrl)

svmLinear.pred <- predict(model.tune,data.test,type="prob")
roc.svmLinear=roc(data.test$response,svmLinear.pred[,2],smooth=F)
plot(roc.svmLinear,print.auc=TRUE)

svmLinear.pred <- predict(model.tune,validation,type="prob")
roc.svmLinear.v=roc(validation$response,svmLinear.pred[,2],smooth=F)
plot(roc.svmLinear.v,print.auc=TRUE)

svmLinear.pred <- predict(model.tune,data.test,type="raw")
a=confusionMatrix(as.factor(svmLinear.pred),data.test$response,positive = 'R') 
results_df$svmLinear=a[["byClass"]]
results_df2$svmLinear=a[["overall"]]

#2.4 svmPoly
model.tune <- train(response ~ .,
                    data = data.train,
                    method = 'svmPoly',
                    metric="ROC",
                    trControl=ctrl)

svmPoly.pred <- predict(model.tune,data.test,type="prob")
roc.svmPoly=roc(data.test$response,svmPoly.pred[,2],smooth=F)
plot(roc.svmPoly,print.auc=TRUE)

svmPoly.pred <- predict(model.tune,validation,type="prob")
roc.svmPoly.v=roc(validation$response,svmPoly.pred[,2],smooth=F)
plot(roc.svmPoly.v,print.auc=TRUE)

svmPoly.pred <- predict(model.tune,data.test,type="raw")
a=confusionMatrix(as.factor(svmPoly.pred),data.test$response,positive = 'R') 
results_df$svmPoly=a[["byClass"]]
results_df2$svmPoly=a[["overall"]]

#2.5 randomForest
library(randomForest)
set.seed(1214)
fit.forest <- randomForest(response~.,
                           data = data.train,importance=TRUE)
forest.pred <- predict(fit.forest,data.test,type="prob")
roc.rf=roc(data.test$response,forest.pred[,2],smooth=F)
plot(roc.rf,print.auc=TRUE)

forest.pred <- predict(fit.forest,validation,type="prob")
roc.rf.v=roc(validation$response,forest.pred[,2],smooth=F)
plot(roc.rf.v,print.auc=TRUE)

#混淆矩阵
forest.pred <- predict(fit.forest,data.test,type="response")
a=confusionMatrix(as.factor(forest.pred),data.test$response,positive = 'R',mode='everything') 
draw_confusion_matrix(a)
results_df$rf=a[["byClass"]]
results_df2$rf=a[["overall"]]


#2.6 AdaBoost
model.tune <- train(response ~ .,
                    data = data.train,
                    method = 'adaboost',
                    metric="ROC",
                    trControl=ctrl)

ada.pred <- predict(model.tune,data.test,type="prob")
roc.ada=roc(data.test$response,ada.pred[,2],smooth=F)
plot(roc.ada,print.auc=TRUE)

ada.pred <- predict(model.tune,validation,type="prob")
roc.ada.v=roc(validation$response,ada.pred[,2],smooth=F)
plot(roc.ada.v,print.auc=TRUE)

ada.pred<- predict(model.tune,data.test,type="raw")
a=confusionMatrix(as.factor(ada.pred),data.test$response,positive = 'R') 
results_df$adaboost=a[["byClass"]]
results_df2$adaboost=a[["overall"]]

#2.7 xgbTree
model.tune <- train(response ~ .,
                    data = data.train,
                    method = 'xgbTree',
                    metric="ROC",
                    trControl=ctrl)

xgb.pred <- predict(model.tune,data.test,type="prob")
roc.xgb=roc(data.test$response,xgb.pred[,2],smooth=F)
plot(roc.xgb,print.auc=TRUE)

xgb.pred <- predict(model.tune,validation,type="prob")
roc.xgb.v=roc(validation$response,xgb.pred[,2],smooth=F)
plot(roc.xgb.v,print.auc=TRUE)

xgb.pred<- predict(model.tune,data.test,type="raw")
a=confusionMatrix(as.factor(xgb.pred),data.test$response,positive = 'R') 
results_df$xgbTree=a[["byClass"]]
results_df2$xgbTree=a[["overall"]]

#2.8 KNN
model.tune <- train(response ~ ,,
                    data = data.train,
                    method = 'knn',
                    metric="ROC",
                    trControl=ctrl)

knn.pred <- predict(model.tune,data.test,type="prob")
roc.knn=roc(data.test$response,knn.pred[,2],smooth=F)
plot(roc.knn,print.auc=TRUE)

knn.pred <- predict(model.tune,validation,type="prob")
roc.knn.v=roc(validation$response,knn.pred[,2],smooth=F)
plot(roc.knn.v,print.auc=TRUE)

knn.pred<- predict(model.tune,data.test,type="raw")
a=confusionMatrix(as.factor(knn.pred),data.test$response,positive = 'R') 
results_df$knn=a[["byClass"]]
results_df2$knn=a[["overall"]]

#2.9 rpart
model.tune <- train(response ~ .
                    data = data.train,
                    method = 'rpart',
                    metric="ROC",
                    trControl=ctrl)

rpart.pred <- predict(model.tune,data.test,type="prob")
roc.rpart=roc(data.test$response,rpart.pred[,2],smooth=F)
plot(roc.rpart,print.auc=TRUE)


rpart.pred <- predict(model.tune,validation,type="prob")
roc.rpart.v=roc(validation$response,rpart.pred[,2],smooth=F)
plot(roc.rpart.v,print.auc=TRUE)

rpart.pred<- predict(model.tune,data.test,type="raw")
a=confusionMatrix(as.factor(rpart.pred),data.test$response,positive = 'R') 
results_df$rpart=a[["byClass"]]
results_df2$rpart=a[["overall"]]

#===========================================================================
# 3. Model Performance Summary
#===========================================================================

# Combine performance metrics
results_ml <- rbind(results_df, results_df2)
write.xlsx(results_ml, 'results/machine_learning_performance.xlsx')

# Reorder columns and add AUC values
results_ml <- results_ml[, c(5,7,4,8,2,3,9,1,6)]
roc_values <- c(0.955, 0.920, 0.864, 0.964, 0.916, 0.903, 0.951, 0.901, 0.943)
results_ml[nrow(results_ml) + 1, ] <- roc_values
rownames(results_ml)[19] <- 'AUC'

# Define color palette
heatmap_colors <- c("#2171B5", "#6BAED6", "#C6DBEF", "#F7FBFF")

# Create performance heatmap
library(ComplexHeatmap)
pheatmap(results_ml[c(1:7,12,13,19),],
         border_color = "white",
         number_format = "%.3f",
         color = rev(heatmap_colors), 
         cellwidth = 35, 
         cellheight = 30,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_color = "white",
         fontsize_number = 10,
         name = "Performance Value")

#===========================================================================
# 4. ROC Curve Analysis
#===========================================================================

# Define color palette for ROC curves
roc_colors <- RColorBrewer::brewer.pal(n = 10, name = 'Paired')

# Validation set ROC curves
validation_rocs <- list(
  'AdaBoost (0.955)' = roc.ada,
  'KNN (0.920)' = roc.knn,
  'Random Forest (0.964)' = roc.rf,
  'Decision Tree (0.864)' = roc.rpart,
  'SVM Linear (0.916)' = roc.svm,
  'SVM Polynomial (0.903)' = roc.svm.polynomial,
  'SVM Radial (0.951)' = roc.svmRadial,
  'Naive Bayes (0.901)' = roc.nb,
  'XGBoost (0.943)' = roc.xgb
)

# Plot validation ROC curves
g <- ggroc(validation_rocs, legacy.axes = TRUE, size = 0.6) +
  scale_colour_manual(values = c("#68d0ce", "#FDBB84", "#f4a4c9", "#e09a2a",
                                 "#9E9AC8", "#de6a73", "#DE77AE", "#7FBC41", "#2171B5")) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed") +
  theme_minimal() +
  theme_bw(base_size = 12) +
  theme(
    legend.position = c(.7,.4),
    legend.key.width = unit(2, "cm"),
    panel.grid.major = element_blank()
  ) +
  guides(color = guide_legend(title = "Model")) +
  ggtitle("Validation Dataset Performance") +
  theme(plot.title = element_text(hjust = 0.5))

# Test set ROC curves
test_rocs <- list(
  'AdaBoost (0.746)' = roc.ada.v,
  'KNN (0.751)' = roc.knn.v,
  'Random Forest (0.773)' = roc.rf.v,
  'Decision Tree (0.706)' = roc.rpart.v,
  'SVM Linear (0.751)' = roc.svm.linear.v,
  'SVM Polynomial (0.759)' = roc.svm.polynomial.v,
  'SVM Radial (0.772)' = roc.svmRadial.v,
  'Naive Bayes (0.747)' = roc.nb.v,
  'XGBoost (0.749)' = roc.xgb.v
)

# Plot test ROC curves
g_test <- ggroc(test_rocs, legacy.axes = TRUE, size = 0.6) +
  scale_colour_manual(values = c("#68d0ce", "#FDBB84", "#f4a4c9", "#e09a2a",
                                 "#9E9AC8", "#de6a73", "#DE77AE", "#7FBC41", "#2171B5")) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), linetype = "dashed") +
  theme_minimal() +
  theme_bw(base_size = 12) +
  theme(
    legend.position = c(.7,.4),
    legend.key.width = unit(2, "cm"),
    panel.grid.major = element_blank()
  ) +
  guides(color = guide_legend(title = "Model")) +
  ggtitle("Independent Test Dataset Performance") +
  theme(plot.title = element_text(hjust = 0.5))

#===========================================================================
# 5. Survival Analysis
#===========================================================================

# Load and prepare data
df1 <- read.csv('data/immunotherapy_cohorts/combined_11_cohorts.csv')
df1$response <- ifelse(df1$response == 'R', 'R', 'NR')
df1$response <- factor(df1$response)

# Split into training and validation sets
train_cohorts <- c('GSE126044', 'rcc', 'GSE166449', 
                   'blca-GSE176307', 'melanoma1', 'melanoma3')
train <- df1[df1$batch %in% train_cohorts, ]
validation <- df1[!df1$batch %in% train_cohorts, ]

# Train Random Forest model
library(randomForest)
set.seed(1214)
rf_model <- randomForest(
  response ~ CHEK2 + LIPT1 + NDUFA10 + AGK + PNPLA2 + GFM1 + TUFM,
  data = train[,3:10],
  importance = TRUE
)

# Make predictions
validation$prediction <- predict(rf_model, validation, type = "response")
validation$prediction_prob <- predict(rf_model, validation, type = "prob")[,2]

# IMVIGOR cohort survival analysis
imvigor_clinical <- data.frame(
  id = rownames(phenoData),
  OS = phenoData$censOS,
  OS.time = phenoData$os * 30
)

osdata <- merge(validation, imvigor_clinical[, c('id', 'OS', 'OS.time')])

# Kaplan-Meier plot
fit <- survfit(Surv(OS.time, OS) ~ prediction, data = osdata)
ggsurvplot(fit, data = osdata, pval = TRUE)

# Optimal cutoff analysis
res.cut <- surv_cutpoint(
  osdata,
  time = "OS.time",
  event = "OS",
  variables = c("prediction_prob")
)
res.cat <- surv_categorize(res.cut)
fit_prob <- survfit(Surv(OS.time, OS) ~ prediction_prob, data = res.cat)

# Survival plot with probability cutoff
ggsurvplot(
  fit_prob,
  data = res.cat,
  censor.shape = "|",
  censor.size = 4,
  conf.int = TRUE,
  conf.int.style = "ribbon",
  conf.int.alpha = 0.2,
  pval = TRUE,
  ggtheme = theme_survminer(),
  palette = c("#b63d3d", "#557aa4"),
  surv.median.line = "hv",
  legend = "top",
  xlab = "OS Time (days)",
  ylab = "Survival Probability",
  title = "IMVIGOR Cohort Survival Analysis",
  break.x.by = 2000,
  break.y.by = 0.2
)

#===========================================================================
# 6. Feature Comparison Analysis
#===========================================================================

# Load and combine all cohorts
combined_df <- rbindlist(list(
  melanoma1, melanoma2, melanoma3,
  blca1, blca2, Glioblastoma,
  lusc1, lusc2, lusc3, rcc, stad
), fill = TRUE)

# Get common columns
common_cols <- Reduce(intersect, list(
  names(melanoma1), names(melanoma2), names(melanoma3),
  names(blca1), names(blca2), names(Glioblastoma),
  names(lusc1), names(lusc2), names(lusc3),
  names(rcc), names(stad)
))
dfall <- combined_df[, ..common_cols]

# Batch effect correction
dfall <- type.convert(dfall)
dfall$response <- case_when(dfall$response == 'N' ~ 'NR', TRUE ~ dfall$response)
dfall <- column_to_rownames(dfall, var = 'id')

combat_data <- ComBat(
  dat = t(dfall[, 1:(ncol(dfall)-2)]),
  batch = dfall$batch,
  mod = model.matrix(~ as.factor(response), data = dfall),
  par.prior = TRUE
)

dfbatch <- t(combat_data) %>%
  as.data.frame() %>%
  mutate(
    response = factor(dfall$response),
    batch = dfall$batch
  )

# Compare with other gene signatures
gene_signatures <- read.xlsx("data/immunotherapy_cohorts/gene_signatures.xlsx")
gene_signatures$lmgene <- dfcoef$term
gene_signatures$lmgene[8:203] <- NA

signature_results <- data.frame()

for (sig_name in colnames(gene_signatures)[c(1:3,5:11)]) {
  sig_genes <- na.omit(gene_signatures[[sig_name]])
  
  # Prepare data
  sig_data <- dfbatch[, colnames(dfbatch) %in% sig_genes] %>%
    mutate(response = dfbatch$response, batch = dfbatch$batch)
  
  # Split data
  train_sig <- sig_data[sig_data$batch %in% train_cohorts, ]
  test_sig <- sig_data[!sig_data$batch %in% train_cohorts, ]
  
  # Train model
  rf <- randomForest(response ~ ., data = train_sig)
  
  # Calculate performance
  performance_metrics <- numeric()
  
  # Validation set performance
  val_pred <- predict(rf, test_sig, type = "prob")
  val_roc <- roc(test_sig$response, val_pred[, 2])
  performance_metrics['Validation AUC'] <- val_roc$auc
  
  # Test set performance
  test_pred <- predict(rf, validation, type = "prob")
  test_roc <- roc(validation$response, test_pred[, 2])
  performance_metrics['Test AUC'] <- test_roc$auc
  
  # Cohort-specific performance
  for (cohort in unique(sig_data$batch)) {
    cohort_pred <- predict(rf, sig_data[sig_data$batch == cohort,], type = "prob")
    cohort_roc <- roc(sig_data[sig_data$batch == cohort,]$response, cohort_pred[, 2])
    performance_metrics[cohort] <- cohort_roc$auc
  }
  
  signature_results <- rbind(signature_results, performance_metrics)
}

# Add LMgene performance
signature_results$lmgene <- c(0.964, 0.773, 0.969, 0.627, 0.893, 0.928, 0.949)

# Create comparison heatmap
heatmap_colors <- c("#3b374c", "#44598e", "#64a0c0", "#7ec4b7", "#deebcd")
pheatmap(t(signature_results[, c(12, 2:11)]),
         border_color = "white",
         number_format = "%.3f",
         color = rev(heatmap_colors),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         name = "AUC Value")