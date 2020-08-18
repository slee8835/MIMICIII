library(ggplot2)
library(caret)
library(rpart.plot)
library(pROC)
library(randomForest)
library("pwr")
library(corrplot)
library("rgl")
library("pscl")
library("car")

setwd('C:/USFMSHI/HS 614/MIMIC')

#read in our dataset
heart <- read.csv("data_heartattack_final1.csv", header = T)

################################################### EDA ##########################################################
summary(heart)

#in total, we have 1835 patients in this dataset
#drop X
#icd9_list
#drop height (1456 NAs) 
#drop weight (1293 NAs)
#drop total_cholesterol  (961 NAs)
#drop hdl_cholesterol (977 NAs)
#drop ldl_cholesterol  (1018 NAs)
#drop triglycerides  (935 NAs)
#drop crp (1757 NAs)   
#drop hba1c  (1170 NAs) 
#drop mean_bp (1546 NAs) 
#drop heart_rate (1293 NAs)    
#drop max_resp_rate  (1293 NAs)         
#drop procedure_codes 
heart <- heart[ , -c(1, 10, 11, 12, 14, 15, 16, 17, 20, 23, 25, 27, 28, 29)]

#we don't need hadm_id, subject_id, diagnosis, discharge location as well
heart <- heart[ , -c(1, 2, 5, 6)]
summary(heart)

#remove rows with more than 3NAs
heart <- heart[-which(rowSums(is.na(heart)) >= 3),]

#ethnicity
summary(heart$ethnicity)  #there are 24 ethnicity types in this dataset 

#regroup ethnicity types: asian, black, hispanic, white, others 
lev_e <- c("ASIAN", "BLACK", "HISPANIC", "WHITE", "OTHERS")
heart$ethnicities <- factor(rep(NA, nrow(heart)), order = F, levels = lev_e)
heart$ethnicities[heart$ethnicity == "PATIENT DECLINED TO ANSWER" | 
                  heart$ethnicity == "UNABLE TO OBTAIN"|
                  heart$ethnicity == "UNKNOWN/NOT SPECIFIED"] <- NA
heart$ethnicities[heart$ethnicity == "ASIAN"| 
                  heart$ethnicity == "ASIAN - ASIAN INDIAN" |
                  heart$ethnicity == "ASIAN - CAMBODIAN" |
                  heart$ethnicity == "ASIAN - CHINESE" |
                  heart$ethnicity == "ASIAN - VIETNAMESE"] <- "ASIAN"

heart$ethnicities[heart$ethnicity == "BLACK/AFRICAN AMERICAN"| 
                    heart$ethnicity == "BLACK/CAPE VERDEAN" |
                    heart$ethnicity == "BLACK/HAITIAN"] <- "BLACK"

heart$ethnicities[heart$ethnicity == "HISPANIC/LATINO - CUBAN" |
                    heart$ethnicity == "HISPANIC OR LATINO" |
                    heart$ethnicity == "HISPANIC/LATINO - CUBAN"] <- "HISPANIC"

heart$ethnicities[heart$ethnicity == "PORTUGUESE" |
                    heart$ethnicity == "WHITE" |
                    heart$ethnicity == "WHITE - BRAZILIAN" |
                    heart$ethnicity == "WHITE - OTHER EUROPEAN" |
                    heart$ethnicity == "WHITE - RUSSIAN"] <- "WHITE"

heart$ethnicities[heart$ethnicity == "MIDDLE EASTERN" |
                    heart$ethnicity == "MULTI RACE ETHNICITY" |
                    heart$ethnicity == "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER" |
                    heart$ethnicity == "OTHER"] <- "OTHERS"

#make sure new factor has the correct count for each level
summary(heart$ethnicities) #345NAs
#remove ethniity column
heart <- heart[ ,-1]
#remove rows with NAs
heart <- heart[-which(is.na(heart$ethnicities)),]

#expire_flag
#change to factor 
class(heart$expire_flag) #right now is int
heart$expire_flag <- as.factor(heart$expire_flag)
levels(heart$expire_flag) <- c("no", "yes")
summary(heart$expire_flag) #make sure there are no NAs

#age 
summary(heart$age) #we want premature heart attack, making sure maximum age is 65

#iculos (in fractional days)
summary(heart$iculos)
#2NA's
heart[which(is.na(heart$iculos)), ]

#impute with mean 
icu_mean <- mean(heart$iculos, na.rm = TRUE)
heart$iculos[is.na(heart$iculos)] <- icu_mean

#red_blood_cells
class(heart$red_blood_cells) #factor right now, change to numeric

heart$red_blood_cells[heart$red_blood_cells == 'UNABLE TO REPORT'] <- NA
heart$red_blood_cells[heart$red_blood_cells == 'UNABLE TO REPORT DUE MCHC >38.0'] <- NA
heart$red_blood_cells <- as.numeric(as.character(heart$red_blood_cells))
summary(heart$red_blood_cells) 

heart[which(is.na(heart$red_blood_cells)), ] #impute with mean

rbc_mean <- mean(heart$red_blood_cells, na.rm = TRUE)
heart$red_blood_cells[is.na(heart$red_blood_cells)] <- rbc_mean
summary(heart$red_blood_cells)

#hb
summary(heart$hb)  #371 NAs
#drop hb column
heart <- heart[,-6]

#creatinine kinase need to be changed to numeric
heart$creatinine_kinase[heart$creatinine_kinase == "ERROR"] <- NA
heart$creatinine_kinase <- as.numeric(as.character(heart$creatinine_kinase))
summary(heart$creatinine_kinase)  #66 NAs, min at 8 and max at 96400

heart[which(is.na(heart$creatinine_kinase)),]  
#impute with median
kinase_median <- median(heart$creatinine_kinase, na.rm = TRUE)
heart$creatinine_kinase[is.na(heart$creatinine_kinase)] <- kinase_median

#ck_isoenzyme
heart$ck_isoenzyme  #NA, NotDone,Greater than 500  
heart$ck_isoenzyme[heart$ck_isoenzyme == 'NotDone'] <- NA
#substitue greater than 500 to NA
heart$ck_isoenzyme[heart$ck_isoenzyme == "GREATER THAN 500"] <- NA
#change to numeric
heart$ck_isoenzyme <- as.numeric(as.character(heart$ck_isoenzyme))
summary(heart$ck_isoenzyme)  #411 NAs
heart <- heart[,-7]

#glucose
summary(heart$glucose)  #more than 445 NAs
#drop glucose
heart <- heart[,-7]

#heart rhythm
summary(heart$heart_rhythm)  #236 NAs 
heart$heart_rhythm <- droplevels(heart$heart_rhythm)

#remove rows with NA
heart <- heart[-which(is.na(heart$heart_rhythm)), ]

#data cleaning complete
summary(heart)

##correlations 
#create new dataframe that holds only numeric variables 
heart_num <- heart[, -c(1, 2, 7, 8)]

#change all factor variables to numeric 
heart_num$gender <- as.numeric(heart$gender)
heart_num$expire_flag <- as.numeric(heart$expire_flag)
heart_num$heart_rhythm <- as.numeric(heart$heart_rhythm)
heart_num$ethnicities <- as.numeric(heart$ethnicities)

#use spearman method when calculate correlation since our dataframe now has ranks 
heart_cor <- cor(heart_num, use = "pairwise.complete.obs", method = "spearman")
heart_cor
corrplot(heart_cor)

############################################## Data Visualizations ###############################################
#ethnicities
g1 <- ggplot(heart, aes(ethnicities))
g1 + geom_bar(aes(fill = ethnicities), width = 0.75) +
  theme(axis.text.x = element_text(angle = 75, vjust=0.6)) +
  labs(title = "Ethnicity Distribution in Heart Dataset")

#gender
g2 <- ggplot(heart, aes(gender))
g2 + geom_bar(aes(fill = gender), width = 0.75)  + 
  labs(title="Gender Distribution in Heart Dataset") 

#expire flag
g3 <- ggplot(heart, aes(expire_flag))
g3 + geom_bar(aes(fill = gender), width = 0.75)  + 
  labs(title="Expire_flag Distribution in Heart Dataset colored by Gender") 

#age
g4 <- ggplot(heart, aes(age))
g4 + geom_histogram(aes(col = gender), bins = 30) + 
  labs(title="Age Distribution in Heart Dataset colored by Gender") 

#iculos
g5 <- ggplot(heart, aes(iculos))
g5 + geom_histogram(aes(fill = expire_flag), bins = 30) + 
  labs(title="ICU Length of Stay Distribution in Heart Dataset colored by Expire_flag") 

#red blood cells
g6 <- ggplot(heart, aes(red_blood_cells))
g6 + geom_histogram(aes(fill = gender), bins = 30) + 
  labs(title="ICU Length of Stay Distribution in Heart Dataset colored by Gender") 

#creatine kinase
g7 <- ggplot(heart, aes(creatinine_kinase))
g7 + geom_histogram(aes(fill = gender)) + 
  labs(title="Creatine kinase Distribution in Heart Dataset colored by Gender") 

#heart_rhythm 
g8 <- ggplot(heart, aes(heart_rhythm))
g8 + geom_bar(aes(fill = heart_rhythm)) +
  labs(title = "Heart rhythm distribution types")

#heart rhythm and expire flag
g8 + geom_bar(aes(fill = expire_flag))

#iculos and heart heart rhythm and gender 
g9 <- ggplot(heart, aes(heart_rhythm, iculos))
g9 + geom_boxplot(aes(fill = gender))

#iculos and age 
g10 <- ggplot(heart, aes(age, creatinine_kinase))
g10 + geom_point(aes(col = gender)) + ylim(0, 10000)

#expire flag and ethnicity
g3 + geom_bar(aes(fill = ethnicities), width = 0.75)  + 
  labs(title="Expire_flag Distribution in Heart Dataset colored by Ethnicities") 

################################################### Modeling ####################################################
#Split your data into train and test sets.
set.seed(3233)
intrain <- createDataPartition(y = heart$expire_flag, p= 0.7, list = FALSE)
training <- heart[intrain,]
testing <- heart[-intrain,]
dim(intrain); dim(training); dim(testing)

# check into preprocessing
anyNA(heart)
summary(heart)

#Fit and evaluate a logistic regression model
glmfit_all <- glm(expire_flag ~ ., family = binomial(),data = training)
summary(glmfit_all)
pR2(glmfit_all)

fit_null <- glm(expire_flag ~ 1, data = training, family = binomial())
fit_step = step(fit_null, scope=list(lower=fit_null, upper=glmfit_all),direction="both")
summary(fit_step)
pR2(fit_step)

p <- predict(fit_step, testing, type = "response")
p
#confusion matrix
expire.pred = rep("no", dim(testing)[1])
expire.pred[p > 0.5] = "yes"

confusionMatrix(factor(testing$expire_flag), factor(expire.pred), positive = "yes")  #66%

#Fit and evaluate an SVM classifier, trying linear, poly and RBF kernels
#linear
trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
grid <- expand.grid(C = seq(0, 2, length =20))
set.seed(3233)

#Training peformed using either formula interface or matrix
svm_Linear_Grid <- train(expire_flag ~., data = training, method = "svmLinear",
                         trControl=trctrl,
                         preProcess = c("center", "scale"),
                         tuneGrid = grid,
                         tuneLength = 10) # param grid with 10 values for each tuned param
# By default, even with tuneLength=x, svmLinear does not tune C
svm_Linear_Grid 
plot(svm_Linear_Grid)

# standardizes test data the same way as the training data 
test_pred <- predict(svm_Linear_Grid, newdata = testing)
test_pred

# To designate positive class as not the reference level:
confusionMatrix(factor(testing$expire_flag), test_pred, positive="yes") #65%

#poly
svm_poly_Grid <- train(expire_flag ~., data = training, method = "svmPoly",
                         trControl=trctrl,                 
                         preProcess = c("scale","center"),
                         tuneLength = 4)
svm_poly_Grid
plot(svm_poly_Grid)

test_pred <- predict(svm_poly_Grid, newdata = testing)
test_pred

# To designate positive class as not the reference level:
confusionMatrix(factor(testing$expire_flag), test_pred, positive="yes") #65%

#rbf
trctrl_rbf <- trainControl(method = "repeatedcv", number = 10, repeats = 10)
svmGrid <- expand.grid(sigma = c(0,0.01, 0.02, 0.025, 0.03, 0.04,
                                 0.05, 0.06, 0.07,0.08, 0.09, 0.1, 0.25, 0.5, 0.75,0.9),
                       C = c(0,0.01, 0.05, 0.1, 0.25, 0.5, 0.75,
                             1, 1.5, 2,5))
set.seed(3233)
svm_Radial <- train(expire_flag ~., data = training, method = "svmRadial",
                    trControl=trctrl_rbf,
                    preProcess = c("center", "scale"),
                    tuneGrid = svmGrid,
                    tuneLength = 10) #param grid 10 values C selected by caret
svm_Radial
plot(svm_Radial)

# standardizes test data the same way as the training data 
test_pred <- predict(svm_Radial, newdata = testing)
test_pred

confusionMatrix(testing$expire_flag, test_pred, positive = "yes")  #68% prediction accuracy


##8. Fit and evaluate a decision tree. Be sure to tune hyperparameters.
trctrl <- trainControl(summaryFunction=twoClassSummary,classProbs = TRUE,# Use AUC to pick the best model
                       savePredictions = T, method = "cv", number = 5)
set.seed(3233)
dtree_fit <- train(expire_flag ~., data = training, method = "rpart",
                   # selects splits, default is gini (impurity)
                   parms = list(split = "information"), #information gain
                   trControl=trctrl,  # same as above
                   tuneLength = 10)

dtree_fit
#NB: complexity param: cp =0 had highest cv accuracy
# cp complexity parameter will be tuned if you leave it out
test_pred <- predict(dtree_fit, newdata = testing)
test_pred

confusionMatrix(factor(testing$expire_flag), test_pred, positive="yes")

#Alt with default control
set.seed(3233)

dtree_fit2 <- train(expire_flag ~., data = training, method = "rpart")
dtree_fit2
#NB: not cv, selected cp = 0.035

test_pred2 <- predict(dtree_fit2, newdata = testing)
test_pred2

confusionMatrix(factor(testing$expire_flag), test_pred2, positive="yes")

prp(dtree_fit$finalModel, box.palette = "Reds")
prp(dtree_fit2$finalModel, box.palette = "Reds")

#Random Forest
#using random search 
trctrl_rf <- trainControl(method="repeatedcv", number=10, repeats=5, search="random")
set.seed(3233)
mtry <- sqrt(ncol(training))
tunegrid <- expand.grid(.mtry=mtry)
rf_gridsearch <- train(expire_flag ~ ., data = training, method="rf", 
                       trControl=trctrl_rf,
                       preProcess = c("center", "scale"), 
                       tuneLength = 15)
print(rf_gridsearch) 
plot(rf_gridsearch)

# standardizes test data the same way as the training data 
test_pred <- predict(rf_gridsearch, newdata = testing)
test_pred

confusionMatrix(factor(testing$expire_flag), test_pred, positive="yes")
