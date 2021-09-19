#####################################################
#                                                   #
#           BACHELOR THESIS                         #
#           MARKUS CSERVENKA - 11706187             #
#                                                   #
#           VIENNA UNIVERSITY OF                    #
#           ECONOMICS AND BUSINESS                  #
#           JULY 2020                               #
#                                                   #
#####################################################

#video-playlist that covers relevant areas https://www.youtube.com/playlist?list=PLzH6n4zXuckpfMu_4Ff8E7Z1behQks5ba

#set directory
setwd("data")

###
###BASIC PARAMETERS
###
#frequency used while measuring data in Hz (default: 50)
samplingFreq = 50
#aggregating original data; value must be lower than samplingFreq; in case aggTime=0 -> no aggregation (default: 10)
aggFreq = 10
#higher local outlier factors are considered outliers and therefore removed (default: 4)
lofMax = 4
#cutoff frequency for low pass filter in Hz (default: 4)
cutoffFreq = 4
#creating features for certain time window in seconds (default: 1)
windowSize = 1


#Column names of different sensors
accNames = c("Time", "AccX", "AccY", "AccZ") #accelerometer
gyrNames = c("Time", "GyrX", "GyrY", "GyrZ") #gyroscope
ligNames = c("Time", "Light") #light
lacNames = c("Time", "LinAccX", "LinAccY", "LinAccZ") #linear acceleration
magNames = c("Time", "MagX", "MagY", "MagZ") #magnetometer
proNames = c("Time", "Proximity") #proximity


#FUNCTIONS
#Assign column names to the data frames in list
renameCols = function(li) { #li = list that contains data frames
  names(li$Accelerometer) = accNames
  names(li$Gyroscope) = gyrNames
  names(li$Light) = ligNames
  names(li$LinearAcceleration) = lacNames
  names(li$Magnetometer) = magNames
  names(li$Proximity) = proNames
  return(li)
}
#Define degree of roundness to receive discrete timesteps (0, 0.5, 0.25, 0.02)
setRoundRatio = function(x, freq) { #x = value #freq = frequency of data in hz
  round((x * freq), 0) / freq
}
#Round Time for each data frame in list
roundTime = function(li, hz) { #li = list that contains data frames #hz = frequency of data in hz
  for (index in seq(1, length(li), 1)) {
    li[[index]][,1] = sapply(li[[index]][,1], setRoundRatio, freq = hz)
  }
  return(li)
}
#Get and arrange data properly
getData = function(folder, hz) { #folder = folder that contains csv files #hz = frequency of data in hz
  files = list.files(folder, pattern="*.csv", full.names=TRUE) #get files
  frames = lapply(files, read.csv) #read files
  names(frames) = c("Accelerometer","Gyroscope","Light","LinearAcceleration", "Magnetometer", "Proximity")
  frames = renameCols(frames) #rename columns
  final = roundTime(frames, hz) #create discrete timesteps
  return(final)
}
#Join data frames
joinDataFrames = function(li, act) { #li = list of data frames #act = activity performed
  frame = li[[1]] #get first frame in list
  for (index in seq(2, length(li), 1)) { #join all data frames in list
    frame = merge(frame, li[[index]], by = "Time", all = TRUE)
  }
  frame["Activity"] = as.factor(act) #declare column activity as factor
  return(frame)
}
#Missing Data Imputation based on Mean of previous and next Value
avgImpValue = function(frame, colVector) { #frame = data frame #colVector = vector of column numbers
  for (col in colVector) { #for each column
    nas = which(is.na(frame[,col])) #get indices with na value in column
    for (index in nas) { #for every na value
      if (index != 1 & index != nrow(frame)) { #do not take first or last instance
        prevVal = frame[index - 1,col] #previous value
        nextVal = frame[index + 1,col] #next value
        if (!any(is.na(c(prevVal, nextVal)))) { #previous and next value are not na
          frame[index,col] = mean(c(prevVal, nextVal))
        } else if (!is.na(prevVal) & is.na(nextVal) & !is.na(frame[index + 2,col])) { #next value is na
          step = (frame[index + 2,col] - prevVal) / 3 #calc 1/3 of range between prevVal and nextVal+1
          frame[index,col] = prevVal + step #impute value
          frame[index + 1, col] = prevVal + step * 2 #impute next value as well
        }
      }
    }
    
  }
  return(frame)
}
#Aggregate instances based on timeperiod
aggregateTimeSteps = function(sourceFrame, targetFrame, aggFreq, smplFreq) { 
  if (aggFreq != 0) { #in case aggFreq = 0 do nothing
    time = 1 / aggFreq #frequency to time
    timeperiod = smplFreq * time #instances aggregated (originally measured in 50Hz)
    timepoint = time #assigning each instance to a descret timepoint
    for (index in seq(1, nrow(sourceFrame), timeperiod)) { #looping without overlapping
      act = as.character(modus(sourceFrame$Activity[index:(index+timeperiod)])) #activity
      means = colMeans(sourceFrame[index:(index+timeperiod),2:13]) #means of numeric values
      f = c(timepoint, means, act) #new instance
      targetFrame[nrow(targetFrame) + 1,] = f #add instance to target frame
      timepoint = timepoint + time #next timepoint
    }
    return(targetFrame)
  }
}
#Equivalent of colSums or colMeans for Standard Deviation/ Max/ Min
colSd = function(x, na.rm = TRUE) { apply(X=x, MARGIN = 2, FUN = sd, na.rm = na.rm) } #x = column
colMax = function(x, na.rm = TRUE) { apply(X=x, MARGIN = 2, FUN = max, na.rm = na.rm) } #x = column
colMin = function(x, na.rm = TRUE) { apply(X=x, MARGIN = 2, FUN = min, na.rm = na.rm) } #x = column
#Get Mode of Vector
modus = function(vec) { #vec = vector
  uVals = unique(vec)
  uVals[which.max(tabulate(match(vec, uVals)))]
}
#Create Features in time and frequency domain
createFeatures = function(sourceFrame, numCols, actCol, windowSize, samplingRate) { #sourceFrame = data frame that provides data #targetFrame = data frame that reveives calculated data
  #instances = number of instances that should be taken into account (time window size)
  freqs = (0:(windowSize-1))*samplingRate/windowSize
  tf = data.frame() #create data frame
  tf[1,"Timepoint"] = NA #define timepoint variable
  tf[1,"Activity"] = NA #defnie activity variable
  for (col in numCols) { #for each numeric column create new column with defined features
    tf[1,paste(names(sourceFrame)[col],"Mean",sep="")] = NA
    tf[1,paste(names(sourceFrame)[col],"Std",sep="")] = NA
    tf[1,paste(names(sourceFrame)[col],"Max",sep="")] = NA
    tf[1,paste(names(sourceFrame)[col],"Min",sep="")] = NA
    tf[1,paste(names(sourceFrame)[col],"MaxFreq",sep="")] = NA
    tf[1,paste(names(sourceFrame)[col],"FreqWeighted",sep="")] = NA
    tf[1,paste(names(sourceFrame)[col],"Pse",sep="")] = NA
  }
  instance = 1 #starting instance for target data frame
  for (index in seq((windowSize + 1), nrow(sourceFrame), (windowSize / 2))) { #50% overlap
    print(paste("Instance: ",index)) #print current status of process
    for (col in numCols) { #calculate values for each new column
      dft = fft(sourceFrame[(index-windowSize):(index-1),col]) #fast fourier transformation
      tf[instance,paste(names(sourceFrame)[col],"Mean",sep="")] = mean(sourceFrame[(index - windowSize):(index - 1),col])
      tf[instance,paste(names(sourceFrame)[col],"Std",sep="")] = sd(sourceFrame[(index - windowSize):(index - 1),col])
      tf[instance,paste(names(sourceFrame)[col],"Max",sep="")] = max(sourceFrame[(index - windowSize):(index - 1),col])
      tf[instance,paste(names(sourceFrame)[col],"Min",sep="")] = min(sourceFrame[(index - windowSize):(index - 1),col])
      tf[instance,paste(names(sourceFrame)[col],"MaxFreq",sep="")] = freqs[which.max(Re(dft))]
      tf[instance,paste(names(sourceFrame)[col],"FreqWeighted",sep="")] = sum(freqs * Re(dft)) / sum(Re(dft))
      pse = Re(dft)^2/windowSize
      pse = pse/sum(pse)
      tf[instance,paste(names(sourceFrame)[col],"Pse",sep="")] = -sum(log(pse)*pse)
    }
    tf$Timepoint[instance] = instance #assign timepoint to current instance value
    tf$Activity[instance] = as.character(modus(sourceFrame$Activity[(index - windowSize):(index - 1)]))
    instance = instance + 1 #next instance
  }
  tf$Activity = factor(tf$Activity) #declare activity as factor
  return(tf)
}
#Get selected attributes of decision tree
getTreeAttributes = function(tr) { #tr = decision tree
  intnodes = sort(unique(where(tr))) #sort nodes in decison tree
  diffnodes = seq(1:intnodes[length(intnodes)]) #create sequence
  primenodes = setdiff(diffnodes,intnodes) #get nodes
  attributes = vector() #create empty attributes vector
  for (i in primenodes){
    temp <- unlist(nodes(tr, i)[[1]][[5]]) #get current node
    attributes <- append(attributes, as.character(temp[length(temp)])) #add node to attributes
  }
  return(unique(attributes)) #return unique values
}
#Calculate Accuracy of Classification models based on result table
accuracy <- function(tbl){sum(diag(tbl)/(sum(rowSums(tbl)))) * 100} #tbl = table containing results of classification 


###
#Import Data
###
#Get data from activities
walkFrames = getData('walk', samplingFreq) #WALK - measurement 0.02 sec resp. 50Hz
cycleFrames = getData('cycle', samplingFreq) #CYCLE - measurement 0.02 sec resp. 50Hz
runFrames = getData('run', samplingFreq) #RUN - measurement 0.02 sec resp. 50Hz
sitFrames = getData('sit', samplingFreq) #SIT - measurement 0.02 sec resp. 50Hz
standFrames = getData('stand', samplingFreq) #STAND - measurement 0.02 sec resp. 50Hz
lieFrames = getData('lie', samplingFreq) #LIE - measurement 0.02 sec resp. 50Hz
testFrames = getData("test", samplingFreq) #TEST - measurement 0.02 sec resp. 50Hz
rm(accNames, gyrNames, lacNames, ligNames, magNames, proNames) #deleting useless variables

#Joining all data frames based on the activity type
source_walk = joinDataFrames(walkFrames, "walk")
source_cycle = joinDataFrames(cycleFrames, "cycle")
source_run = joinDataFrames(runFrames, "run")
source_sit = joinDataFrames(sitFrames, "sit")
source_stand = joinDataFrames(standFrames, "stand")
source_lie = joinDataFrames(lieFrames, "lie")
testSource = joinDataFrames(testFrames, "lie")
#set activities for test data manually
activities = c("walk", "stand", "sit", "stand", "walk", "run", "walk", "stand", "cycle", "sit", "lie")
timesteps = c(20.05, 24.75, 35.7, 38, 50.35, 82.75, 89.2, 92.3, 170.65, 193.9, 209.3)
i = length(activities) - 1 #assign i for while loop
testSource$Activity = "lie" #activity column to character (with last level measured)
while (i >= 1) { #assign activities to instances
  testSource$Activity[which(testSource$Time < timesteps[i])] = activities[i]
  i = i - 1
}
testSource$Activity = factor(testSource$Activity) #declare activity as factor in test frame
rm(i, activities, timesteps) #delete useless variables
str(testSource)

#Deleting first 9 instances because sensors did not deliver data 
source_walk = source_walk[-(1:9),]
source_cycle = source_cycle[-(1:9),]
source_run = source_run[-(1:9),]
source_sit = source_sit[-(1:9),]
source_lie = source_lie[-(1:9),]
source_stand = source_stand[-(1:9),]
#Delete list of data frames
rm(walkFrames, cycleFrames, runFrames, sitFrames, standFrames, lieFrames, testFrames)

#Bind all source data frames for training set
source_base = rbind(source_cycle, source_lie, source_run, source_sit, source_stand, source_walk)
rm(source_cycle, source_lie, source_run, source_sit, source_stand, source_walk)


###
###Missing Values
###
#Handling Missing Values for training data
sapply(source_base, function(y) sum(is.na(y))/length(y)*100) #share of missing values in % by column
all_base = avgImpValue(source_base, c(2:7,9:14)) #imputate values for Acc, Gyr, Mag
sapply(all_base, function(y) sum(is.na(y))/length(y)*100) #share of missing values in % by column
all_base[,c(8,15)] = NULL #Remove Light and Proximity due to amount of missing values
all_base = all_base[complete.cases(all_base),] #remove left over missing values

#Handling Missing Values for test data
sapply(testSource, function(y) sum(is.na(y))/length(y)*100) #share of missing values in % by column
test = avgImpValue(testSource, c(2:7,9:14)) #imputate values for Acc, Gyr, Mag
test[,c(8,15)] = NULL #Remove Light and Proximity due to amount of missing values
test = test[complete.cases(test),] #remove left over missing values


###
###Aggregating values
###
#Define Column Names for new target frame
cols = c("Time", "AccX", "AccY", "AccZ", "GyrX", "GyrY", "GyrZ", "LinAccX", "LinAccY", "LinAccZ", "MagX", "MagY", "MagZ", "Activity")
#Aggregate data for training data
agg_base = data.frame(matrix(ncol = 14, nrow = 0)) #create target frame
colnames(agg_base) = cols #assign names to columns
agg_base = aggregateTimeSteps(all_base, agg_base, aggFreq, samplingFreq) #aggregate data
sapply(agg_base, function(y) sum(is.na(y))/length(y)*100) #check missing values
agg_base = agg_base[complete.cases(agg_base),] #in case there are missing values remove them
agg_base$Activity = factor(agg_base$Activity) #declare values as factor
agg_base[1:13] = sapply(agg_base[1:13], as.numeric) #declare values as numeric
str(agg_base)

#Aggregate data for test data
testAgg = data.frame(matrix(ncol = 14, nrow = 0)) #create target frame
colnames(testAgg) = cols #assign names to columns
testAgg = aggregateTimeSteps(test, testAgg, aggFreq, samplingFreq) #aggregate data
testAgg = testAgg[complete.cases(testAgg),] #in case there are missing values remove them
testAgg$Activity = factor(testAgg$Activity) #declare values as factor
testAgg[1:13] = sapply(testAgg[1:13], as.numeric) #declare values as numeric
str(testAgg)


###
###Outlier Detection
###
#Multivariate Analysis - Local Outlier Factor
if (!match("dbscan", rownames(installed.packages()))) { install.packages("dbscan") }
library(dbscan) #load dbscan package
lofNeighbours = 5
#Calculate and add LOF values
lofsBase = (lof(agg_base[,2:13], lofNeighbours)) #generate lof values
sum(lofsBase>lofMax) #number of Outlier Instances
agg_base["LOF"] = lofsBase #create lof column and add values
#Visualize lof values for training data
par(mfrow=c(1,1))
plot(agg_base$Time, agg_base$LOF, ylab = "LOF - Value", xlab = "Time (Sec.)", main = "Local Outlier Factor")
abline(lofMax, 0, col = "red") #line shows max lof value
base = agg_base #create new data frame named base
base[which(base$LOF > lofMax), 2:13] = NA #set outliers to na
#Imputate missing values (got rid of outliers)
base = avgImpValue(base, c(2:13)) #imputate missing values
base = base[complete.cases(base),] #remove left over missing values
rm(lofMax, lofNeighbours, lofsBase) #delete useless variables


###
###Data Transformation
###
#Butterworth Low Pass Filter
if (!match("signal", rownames(installed.packages()))) { install.packages("signal") }
library(signal) #load signal package
#Nyquist Frequency is highest Frequency we can observe in data (50Hz resp. 10Hz after aggregating data)
#Nyquist Frequency: fNQ = f0 / 2 -> 10 / 2 = 5
#Therefore measurements have to be taken at least twice per frequency #f0 >= 2 fs (Sample Frequency)
#Parameter w = fcutoff / fNQ -> 4 / 5 = 0.8
#Sources: https://www.youtube.com/watch?v=sgYkOkrlQ_E, https://www.r-bloggers.com/butterworth-filter-overshoot/
#         https://stackoverflow.com/questions/46697836/butterworth-filter-using-scipy

if(aggFreq != 0) {w = cutoffFreq / (aggFreq / 2)} else {w = cutoffFreq / (samplingFreq / 2)} #define w paramter
bf <- butter(1, w, type="low") #Create Butterworth Filter
#Low Pass Filter for training data
base_before = base #for visualization purposes save unfiltered data
base[,2:13] = data.frame(sapply(base[,2:13], function(x) filtfilt(bf, x))) #apply filter
#Visualization
par(mfrow=c(2,1))
plot(base_before$Time[100:700], base_before$GyrX[100:700], type = "l", main = "Data before LPF", ylab = "GyrX", xlab = "Time (Sec)")
plot(base$Time[100:700], base$GyrX[100:700], type = "l", main = "Data after LPF", ylab = "GyrX", xlab = "Time (Sec.)")
rm(base_before) #delete useless variables

#Low Pass Filter for test data
testAgg[,2:13] = data.frame(sapply(testAgg[2:13], function(x) filtfilt(bf, x))) #apply filter
rm(bf) #delete useless variables
test = testAgg #move on with data frame named test

#Principal Component Analysis for training data
pca = prcomp(base[,2:13], scale = TRUE) #principal component analysis
head(pca$x) #show values of pcas
summary(pca)
#Visual Representation
par(mfrow=c(1,1))
varExp = pca$sdev^2 / sum(pca$sdev^2) #explained variance by principal component
plot(varExp, type = "b", ylab = "Explained Variance in %", xlab = "Principal Components", main = "Principal Component Analysis")
axis(side = 1, at=seq(1, length(varExp), 1))
#Add principal components to base data frame
base = cbind(base, pca$x[,1:9])
#Correlation between attributes and principal components
cor(base[,2:13], base[,16:24])
rm(pca, varExp) #delete useless variables

#Principal Component Analysis for test data
pca = prcomp(test[,2:13], scale = TRUE) #principal component analysis
#Add principal components to base data frame
test = cbind(test, pca$x[,1:9])
rm(pca) #delete useless variables


###
###Create Features with overlapped time windows (50%)
### 
#Link for FFT: http://www.di.fc.ul.pt/~jpn/r/fourier/fourier.html
if(aggFreq != 0) {ws = windowSize * aggFreq} else {ws = windowSize * samplingFreq} #instances to take into account (10hz * 1sec = 10 instances)
if(aggFreq != 0) {t = 1 / aggFreq} else {t = 1 / samplingFreq} #define timesteps based on frequency
#Create features for training data
baseFinal = createFeatures(base, c(2:7,16:24), 14, ws, t)
#Create features for test data
testFinal = createFeatures(test, c(2:7,15:23), 14, ws, t)
str(testFinal)
rm(t, ws) #delete useless variables


###
###Machine Learning Preperation
###
set.seed(1234) #set seed for getting same results when executing
#Random assignment of instances to train and test data frames
smpl = sample.int(n = nrow(baseFinal), size = floor(0.7*nrow(baseFinal)), replace = F)
train = baseFinal[smpl,-1] #create train data frame but without Timepoint
valid = baseFinal[-smpl,-1] #create validation data frame but without Timepoint
rm(smpl) #delete useless variables
summary(train)
summary(valid)


###
###Classification Techniques
###
#Decision Tree
if (!match("party", rownames(installed.packages()))) { install.packages("party") }
library(party) #load party package
#Create Decison Tree
#since decision tree is used for feature selection as well, all features are included
form = Activity ~ . #define formula
cltree = ctree(form, train) #create decision tree based on all variables
#Classification Results with decision tree model on train data frame
tb.dc = table(predict(cltree), train$Activity) #in-sample
tb.dc
accuracy(tb.dc)
plot(cltree) #plot decision tree

#Feature Selection
feats = getTreeAttributes(cltree) #create vector of selected features
formFeats = as.formula(paste("Activity~", paste(feats, collapse = "+"), sep = "")) #create formula with selected features

#Classification Results with decision tree model on validation data
predictDt = predict(cltree, newdata = valid)
tb.dt = table(predictDt, valid$Activity) #save confusion matrix
tb.dt #show result table (col = actual, row = predicted)
accuracy(tb.dt) #accuracy of dt on validation data

#K-Nearest-Neighbor
if (!match("class", rownames(installed.packages()))) { install.packages("class") }
library(class) #load class package
#Apply KNN Method (k = 5 Neighbours)
predictKnn = knn(train[,feats], valid[,feats], train$Activity, 5) #apply knn on validation data
tb.knn = table(predictKnn, valid$Activity) #save confusion matrix
tb.knn #show result table (col = actual, row = predicted)
accuracy(tb.knn) #accuracy of knn on validation data

#Support Vector Machines
if (!match("e1071", rownames(installed.packages()))) { install.packages("e1071") }
library(e1071) #load e1071 package
#Create svm model
svmfit = svm(formFeats, data = train, kernel = "linear", cost = 0.01, scale = FALSE)
print(svmfit) #print amount of support vectors
predictSvm = predict(svmfit, valid[,feats], type = "class") #apply svm on validation data
tb.svm = table(predictSvm, valid$Activity) #save confusion matrix
tb.svm #show result table (col = actual, row = predicted)
accuracy(tb.svm) #accuracy of svm on validation data


###
###Test Methods on test data
###
#Decision Tree
predictDtTest = predict(cltree, newdata = testFinal) #apply decision tree on test data
tb.dt.test = table(predictDtTest, testFinal$Activity) #save confusion matrix
tb.dt.test #show result table (col = actual, row = predicted)
accuracy(tb.dt.test) #accuracy of dt on test data

#K-Nearest-Neighbor (k = 5 Neighbours)
predictKnnTest = knn(train[,feats], testFinal[,feats], train$Activity, 5) #apply knn on test data
tb.knn.test = table(predictKnnTest, testFinal$Activity) #save confusion matrix
tb.knn.test #show result table (col = actual, row = predicted)
accuracy(tb.knn.test) #accuracy of knn on test data

#Support Vector Machines
predictSvmTest = predict(svmfit, testFinal[,feats], type = "class") #apply svm on test data
tb.svm.test = table(predictSvmTest, testFinal$Activity) #save confusion matrix
tb.svm.test #show result table (col = actual, row = predicted)
accuracy(tb.svm.test) #accuracy of svm on test data