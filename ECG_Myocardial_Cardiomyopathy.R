# Penalised B-Spline Function

pspline_fn <- function(splinebasis, st, en, y, lambda, 
                       argvals=c(1:30000), plot=FALSE, derv = FALSE){
  #  using a smoothing parameter that is set here to lambda
  lambda = 10^(lambda)
  fdParobj = fdPar(fdobj=splinebasis, Lfdobj=2,
                   lambda=lambda )
  
  ys = smooth.basis(argvals=argvals, y=y, 
                    fdParobj=fdParobj)
  xfd = ys$fd
  if (plot){
    #plot data points and then fit
    plot(y[st:en], type="l",
         main="Original timeseries along with P-spline")
    lines(c(st:en), eval.fd(c(st:en), xfd), 
          col = "red")
  }
  

# Get ECG Feature Function
  get_ecg_feature <- function(heartbeat, plot=FALSE){
    #Comment out below heartbeat definition
    heartbeat <- data.new[26,]
    plot=TRUE
    #Fit spline on each heartbeat
    xfd = pspline_fn(create.bspline.basis(c(1,200),31), 1, 200, heartbeat, 
                     -4.5, argvals=c(1:200), plot=FALSE,derv = TRUE)
    
    temp_crv <- eval.fd(c(1:200), xfd)
    temp_crv <- rescale(temp_crv)
    
    #The median is the baseline of heart beat
    heartbeat_baseln <- median(temp_crv)
    
    #1-27 S wave
    #45-125 T wave
    #126-175 P wave
    #176-200 Q wave
    
    #Use derivates to find keypoints
    derv1 <- round(eval.fd(c(1:200), xfd, Lfd=1), digits = 3)
    derv2 <- round(eval.fd(c(1:200), xfd, Lfd=2), digits = 3)
    
    check1 <- abs(derv1) <= 0.003
    check2 <- abs(derv2) >= 0.001
    
    key_points <- which((check1) & (check2))
    
    #1-25 S wave
    # Check if any element in the list is between 1 and 25 (inclusive)
    #S_wave_condition - slope close to 0, positive second order derivative
    S_wave_condition <- ((abs(derv1) <= 0.003) & (derv2 > 0))
    S_wave_idxs <- which(S_wave_condition)
    S_wave_loc_condition <- (S_wave_idxs >= 1 & S_wave_idxs <= 27)
    between_1_27 <- any(S_wave_loc_condition)
    # Extract points between 1 and 27, where 1st derivative is close to 0, 
    # and 2nd derivative is positive if it exists, otherwise NA
    S_wave_idxs <- if(between_1_27) S_wave_idxs[S_wave_loc_condition] else NA
    #Check all candidate S-wave points.
    #If single candidate, choose it
    #If it is NA (no candidate), then S_idx is NA
    #If multiple candidates, choose point which is minimum
    # if (length(S_wave_idxs) == 1){
    #   S_wave_idx <- S_wave_idxs
    # } else if (any(is.na(S_wave_idxs))) {
    #   S_wave_idx <- NA
    # } else {
    #   S_wave_idx  <-  S_wave_idxs[which.min(temp_crv[S_wave_idxs])]
    # }
    # Alternative 2 - Select the minimum among first 27 points
    S_wave_idx <- c(1:27)[which.min(temp_crv[1:27])]
    
    #45-125 T wave
    #T wave can be inverted or normal both
    T_wave_condition <- ((check1) & (check2))
    T_wave_idxs <- which(T_wave_condition)
    T_wave_loc_condition <-   (T_wave_idxs >= 45 & T_wave_idxs < 125)
    between_45_125 <- any(T_wave_loc_condition)
    # Extract number between 45 and 125 if it exists, otherwise NA
    T_wave_idxs <- if (between_45_125) T_wave_idxs[T_wave_loc_condition] else NA
    #Check all candidate T-wave points.
    #If single candidate, choose it
    #If it is NA (no candidate), then T_idx is NA
    #If multiple candidates, choose point which is minimum
    if (length(T_wave_idxs) == 1){
      T_wave_idx <- T_wave_idxs
    } else if (any(is.na(T_wave_idxs))) {
      T_wave_idx <- NA
    } else{
      T_wave_idx  <-  T_wave_idxs[which.max(abs(temp_crv[T_wave_idxs] - heartbeat_baseln))]
    }
    #Check if T-wave is inverted, second derivative is positive.
    #Store info as T-wave_inversion feature
    if (is.na(T_wave_idx)){
      T_wave_inversion <- NA
    } else if (derv2[T_wave_idx] > 0){
      T_wave_inversion <- 1
    } else {
      T_wave_inversion <- 0 
    }
    
    
    #126-180 P wave
    #P_wave_condition - slope close to 0, negative second order derivative
    P_wave_condition <- ((abs(derv1) <= 0.003) & (derv2 < 0))
    P_wave_idxs <- which(P_wave_condition)
    P_wave_loc_condition <- (P_wave_idxs >= 126 & P_wave_idxs < 180)
    between_126_175 <- any(P_wave_loc_condition)
    P_wave_idxs <- if (between_126_175) P_wave_idxs[P_wave_loc_condition] else NA
    #Check all candidate P-wave points.
    #If single candidate, choose it
    #If it is NA (no candidate), then P_idx is NA
    #If multiple candidates, choose point which is maximum
    if (length(P_wave_idxs) == 1){
      P_wave_idx <- P_wave_idxs
    } else if (any(is.na(P_wave_idxs))) {
      P_wave_idx <- NA
    } else {
      P_wave_idx  <-  P_wave_idxs[which.max(temp_crv[P_wave_idxs])]
    }
    
    #176-200 Q wave
    #Q_wave_condition - slope close to 0, positive second order derivative
    Q_wave_condition <- ((abs(derv1) <= 0.003) & (derv2 > 0))
    Q_wave_idxs <- which(Q_wave_condition)
    Q_wave_loc_condition <- (Q_wave_idxs >= 176 & Q_wave_idxs <= 200)
    between_176_200 <- any(Q_wave_loc_condition)
    # Extract the first number between 176 and 200 if it exists, otherwise NA
    Q_wave_idxs <- if(between_176_200) Q_wave_idxs[Q_wave_loc_condition] else NA
    #Alternative 1
    #Check all candidate Q-wave points.
    #If single candidate, choose it
    #If it is NA (no candidate), then Q_idx is NA
    #If multiple candidates, choose point which is minimum
    # if (length(Q_wave_idxs) == 1){
    #   Q_wave_idx <- Q_wave_idxs
    # } else if (any(is.na(Q_wave_idxs))) {
    #   Q_wave_idx <- NA
    # } else {
    #   Q_wave_idx  <-  Q_wave_idxs[which.min(temp_crv[Q_wave_idxs])]
    # }
    #Alternative 2, just take minimum point between 176-200
    Q_wave_idx <- c(176:200)[which.min(temp_crv[176:200])]
    
    r_s_t_p_q_idx <- c(1, S_wave_idx, T_wave_idx, P_wave_idx, Q_wave_idx)
    #PR interval + TP interval
    #Median from Mid point between P and T, till Q
    heartbeat_baseln_new <- median(temp_crv)
    
    #Plot heartbeat with all features
    if (plot){
      plot(heartbeat, type="l", main="Heartbeat with Spline fit on heartbeat",
           lwd=2)
      lines(temp_crv, col="red", lwd=2, alpha=0.8)
      lines(rep(heartbeat_baseln_new, 200), col="green")
      points(key_points, temp_crv[key_points])
      points(r_s_t_p_q_idx, temp_crv[r_s_t_p_q_idx], col="blue", pch=19)
    }
    
    
    # If any of the detected points is NA, that means not an ideal sample, 
    # don't return features, else return all the features.
    if (any(is.na(r_s_t_p_q_idx))) {
      return (NA)
    } else{
      #Features
      #All amplitudes RSTPQ, calculated in one go below
      # s_wave_amplitude
      # r_amplitude
      all_wave_amplitudes <- (temp_crv[r_s_t_p_q_idx] - heartbeat_baseln_new)
      #S-R + R-Q, R is at 1 and 200
      r_duration <- S_wave_idx + (200 - Q_wave_idx)
      # s_duration <- TBD
      # qrs_duration same as r_duration?
      #qt_interval = R-Q + t_index
      qt_interval = (200 - Q_wave_idx) + T_wave_idx
      # heart_rate <- TBD will have to process at patient level
      #Height different beween Q and S points
      q_s_ht_diff <- temp_crv[Q_wave_idx] - temp_crv[S_wave_idx]
      #Height different beween S and T points
      s_t_ht_diff <- temp_crv[T_wave_idx] - temp_crv[S_wave_idx]
      #Height different beween R and T points
      r_t_ht_diff <- temp_crv[1] - temp_crv[T_wave_idx]
      #R_s height ratio
      r_s_ht_ratio <- abs(all_wave_amplitudes[1]) / abs(all_wave_amplitudes[2])
      #R_slope
      r_slope <- mean(derv1[2:3])
      #S_wave_2nd_derv
      S_wave_2nd_derv <- derv2[S_wave_idx]
      #Q_wave_2nd_derv
      Q_wave_2nd_derv <- derv2[Q_wave_idx]
      #T_wave_2nd_derv
      T_wave_2nd_derv <- derv2[T_wave_idx]
      #T_wave_slope_bf
      T_wave_slope_bf <- mean(derv1[(T_wave_idx-5):(T_wave_idx-3)])
      #T_wave_slope_af
      T_wave_slope_af <- mean(derv1[(T_wave_idx+3):(T_wave_idx+5)])
      #P_wave_2nd_derv
      P_wave_2nd_derv <- derv2[P_wave_idx]
      #P_wave_slope_bf
      P_wave_slope_bf <- mean(derv1[(P_wave_idx-5):(P_wave_idx-3)])
      #P_wave_slope_af
      P_wave_slope_af <- mean(derv1[(P_wave_idx+3):(P_wave_idx+5)])
      #pq_interval
      pq_interval <- Q_wave_idx - P_wave_idx
      
      #TODO
      # Hermite coefficients
      # wavelet features
      #STelevation
      #Assuming S wave ends in 5 points, T wave starts 10 points before
      st_seg_points <- c((S_wave_idx+5),max((T_wave_idx-25), (S_wave_idx+10)))
      #Choice to be made whether to use raw data or spline curve.
      st_elevation <- median(temp_crv[st_seg_points[1]:st_seg_points[2]]) - heartbeat_baseln_new
      if (plot){
        # Plot the ST segment points for debugging
        points(st_seg_points, temp_crv[st_seg_points], col="chocolate", pch=17)
        lines(rep(median(temp_crv[st_seg_points[1]:st_seg_points[2]]),100), col="cyan3") 
      }
      return (c(all_wave_amplitudes,r_duration, qt_interval, q_s_ht_diff, 
                s_t_ht_diff,r_t_ht_diff, r_s_ht_ratio, r_slope, S_wave_2nd_derv, 
                Q_wave_2nd_derv, T_wave_2nd_derv, T_wave_slope_bf, T_wave_slope_af, 
                P_wave_2nd_derv, P_wave_slope_bf, P_wave_slope_af, pq_interval,
                st_elevation, T_wave_inversion))
    }
  }
  
  
#Preprocess Function
  pre_process <- function(data, p=200, plot=FALSE, y_label=NA){
    #empty matrix to store results
    #24 features
    data.n_p<-matrix(NA,dim(data)[1],24);
    
    for (i in 1:dim(data)[1]){
      #consider an individual 
      x <- data[i,]
      
      #Comment out 5 rows below
      x <- X.train[13,]
      y.predict.lda[90]
      plot <- TRUE
      i <- 1
      y_label <- y.train
      
      #Fit bspline to remove non-linear trend.
      # Subtract spline fit from raw data to detrend
      nbasis <- 31
      sp_basis = create.bspline.basis(c(1,30000),nbasis)
      xspline_trend = pspline_fn(sp_basis,1, 30000, x,-4.5, plot = plot)
      x = x-xspline_trend
      
      #plot
      if (plot){
        plot(x,type="l", main=paste("After Bspline for Patient - ", i))
      }
      #Remove outliers
      x = rm_outlier(x)
      
      #rescale between 0-1
      x = rescale(x)
      
      #Drop first and last 2000 points.
      # Since Spline fitting at end extremes is not great
      x = x[2000:27999]
      
      #plot
      if (plot){
        plot(x[1:10000],type="l", main=paste("After Bspline & Outlier RM for Patient - ", i))
      }
      
      # Findeach"bigpeak"perheartbeatcycle
      #Option 1 - Findpeaks function
      peaks1 <- findpeaks(x, MinPeakHeight=0.65, MinPeakDistance=450)
      inds <- peaks1$loc
      if (plot){
        points(inds, x[inds], col="red")
      }
      
      #chop the signal up and represent each heart beat as a vector
      m<-length(inds)-1;
      data.length <- matrix(NA,m,1)
      
      #Calculate length of slices to find unusual long & short slices
      for(m_ in 1:m){
        data.length[m_,] <- length(inds[m_]:inds[m_+1])
      }
      
      #This doesn't work for patient 114
      outliers = ?boxplot.stats(data.length)$out
      if (identical(outliers, integer(0))){
        data.length.f <- data.length
      } else {
        #Calculate variance in heartbeat (R-R) length
        data.length.f <- data.length[-which(data.length %in% outliers)]
      }
      heart_bt_var <- var(rescale(data.length.f))
      
      #New_m is matrix without outlier slices
      new_m = (m - length(outliers))
      data.new <- matrix(NA, new_m,p);
      j = 0
      for(m_ in 1:m){
        # If the slice is not an outlier only then store in data.new
        if (!length(inds[m_]:inds[m_+1]) %in% outliers){
          j = j+1
          t.old<-inds[m_]:inds[m_+1]
          x.old<-x[inds[m_]:inds[m_+1]]
          #Rescale slice from 0-1
          x.old <- rescale(x.old)
          t.new<-?seq(inds[m_],inds[m_+1],length=p)
          x.new<-approx(t.old,x.old,t.new)$y
          data.new[j,]<-x.new
        }
      }
      
      #PLOTS ONLY FOR VISUALIZATION
      #Plot with ECG data with peaks highlighted in red circle
      if (plot){
        plot(x, type="l", main = paste("R-wave peaks for Patient - ", i))
        points(inds,x[inds],col=2)
        
        #Fit P-spline of individual heart beat to smoothen data
        # temp <- pspline_fn(nbasis,1, 30000, data.new[11,],-4.5, plot = plot)
        
        #Plot 1st, 2nd and 11th heartbeat and average heartbeat
        plot(data.new[10,],type="l", col="black", lwd=1, ylim = c(0,1.1), 
             main=paste("Average Heartbeat for Patient - ", i))
        lines(data.new[11,], col="black", lwd=1)
        lines(colMeans(data.new), col="red", lwd=3)
        legend(x = "topright",          # Position
               legend = c("Single beat", "Average"),  # Legend texts
               lty = c(1, 1),           # Line types
               col = c("black", "red"),           # Line colors
               lwd = c(1,3))
      }
      
      # Take the mean of all the heart beats to get a 1x200 vector
      # data.new.mu  = colMeans(data.new)
      
      #store this vector for the patient i
      # data.n_p[i,] = data.new.mu
      final_features <- NA
      for (m_ in 1:new_m){
        features_ecg <- get_ecg_feature(data.new[m_,], plot=plot)
        if (any(is.na(features_ecg))){
          print("Discard heartbeat")
        }else{
          final_features <- rbind(final_features, 
                                  data.frame(t(features_ecg)), make.row.names=FALSE)
        }
      }
      
      colnames(final_features) <- c("R_amp", "S_amp", "T_amp", "P_amp", "Q_amp",
                                    "R_duration",
                                    "QT_interval", "QS_ht_diff", "ST_ht_diff", 
                                    "RT_ht_diff", "r_s_ht_ratio", "r_slope", "S_wave_2nd_derv", 
                                    "Q_wave_2nd_derv", "T_wave_2nd_derv", "T_wave_slope_bf", "T_wave_slope_af", 
                                    "P_wave_2nd_derv", "P_wave_slope_bf", "P_wave_slope_af", "pq_interval",
                                    "ST_elev", "T_inversion")
      
      # Average over valid heartbeats and store this vector for the patient i
      data.n_p[i,] = c(colMeans(final_features[-1,-23]), Mode(final_features[,23]),
                       heart_bt_var)    
      # #Append all heart beats
      # final_features <- cbind(final_features, heart_bt_var[1], y_label[i])
      # data.n_p<- rbind(data.n_p, final_features[-1,], make.row.names=FALSE)
      dim(final_features)
      #End patient loop  
    }
    # Return the pre-processed data
    colnames(data.n_p) <- c("R_amp", "S_amp", "T_amp", "P_amp", "Q_amp",
                            "R_duration",
                            "QT_interval", "QS_ht_diff", "ST_ht_diff", 
                            "RT_ht_diff", "r_s_ht_ratio", "r_slope", "S_wave_2nd_derv", 
                            "Q_wave_2nd_derv", "T_wave_2nd_derv", "T_wave_slope_bf", "T_wave_slope_af", 
                            "P_wave_2nd_derv", "P_wave_slope_bf", "P_wave_slope_af", "pq_interval",
                            "ST_elev", "T_inversion", "heart_bt_var")
    return (data.n_p)
    # return (data)
  }

#Modelling and Cross Validation

  X1.n_p_full = pre_process(X.train[1:115,], plot=FALSE, y_label = y.train)
  X1.n_p_full <- data.frame(X1.n_p_full)
  
  X1.n_p_full_export <- cbind(X1.n_p_full,y.train)
  X1.n_p_full_export <- round(X1.n_p_full_export, digits = 4)
  write.csv(X1.n_p_full_export, "train_data_ecg_features.csv" )
  
  X1.n_p <- X1.n_p_full
  
  #As factor the y labels
  y.train <- as.factor(y.train)
  
  #Cross-val start
  folds = createFolds(y.train, k = 4)
  cv = lapply(folds, function(x) { # start of function
    # in the next two lines we will separate the Training set into it's 10 pieces
    X1.n_p.tr = X1.n_p[-x, ] # training fold =  training set minus (-) it's sub test fold
    X1.n_p.val = X1.n_p[x, ] # here we describe the test fold individually
    
    y.train.tr = y.train[-x] # training fold =  training set minus (-) it's sub test fold
    y.train.val = y.train[x]
    
    
    #Naives Bayes
    nb<-NaiveBayes(x=as.data.frame(X1.n_p.tr),grouping=as.factor(y.train.tr)
                   ,usekernel=TRUE)
    #Naive bayes prediction and accuracy
    dim(X1.n_p.tr)
    y.predict.nb <- predict(nb, newdata = as.data.frame(X1.n_p.val))$class
    misClassError.nb <- mean(y.train.val != y.predict.nb)
    # print(paste('Naive Bayes Accuracy =', 1-misClassError.nb)) 
    cm <- confusionMatrix(y.predict.nb, as.factor(y.train.val))
    mean_sensitivity.nb <- mean(cm$byClass[,1])
    
    #Knn
    y.predict.knn <- knn(train = X1.n_p.tr, 
                         test = X1.n_p.val, 
                         cl = as.factor(y.train.tr), 
                         k = 3) 
    misClassError.knn <- mean(y.train.val != y.predict.knn)
    cm <- confusionMatrix(y.predict.knn, as.factor(y.train.val))
    mean_sensitivity.knn <- mean(cm$byClass[,1])
    # print(paste('KNN Accuracy =', 1-misClassError.knn))
    
    #Random Forrest
    rf <- randomForest(x = X1.n_p.tr, y = as.factor(y.train.tr), ntree = 300, 
                       replace = FALSE, proximity = TRUE) 
    y.predict.rf <- predict(rf, X1.n_p.val)
    misClassError.rf <- mean(y.train.val != y.predict.rf)
    cm <- confusionMatrix(y.predict.rf, as.factor(y.train.val))
    mean_sensitivity.rf <- mean(cm$byClass[,1])
    
    #SVM
    svm_ecg_rad <- svm(x = X1.n_p.tr, y = as.factor(y.train.tr), type = "C-classification", 
                       kernel = "radial", cost = 10, gamma = 0.01, class.weights = "inverse")
    y.predict.svm.rad<- predict(svm_ecg_rad, newdata = X1.n_p.val, type = "decision")
    misClassError.svm <- mean(y.train.val != y.predict.svm.rad)
    cm <- confusionMatrix(y.predict.svm.rad, as.factor(y.train.val))
    mean_sensitivity.svm <- mean(cm$byClass[,1])
    
    #LDA
    lda.fit <- lda(X1.n_p.tr[c(-8,-9,-10)], as.factor(y.train.tr))
    y.predict.lda <- predict(lda.fit,X1.n_p.val[c(-8,-9,-10)])$class
    misClassError.lda <- mean(y.train.val != y.predict.lda)
    cm <- confusionMatrix(y.predict.lda, as.factor(y.train.val))
    mean_sensitivity.lda <- mean(cm$byClass[,1])
    
    #QDA
    # qda.fit <- qda(X1.n_p.tr[c(-8,-9,-10)], as.factor(y.train.tr))
    # y.predict.qda <- predict(qda.fit,X1.n_p.val[c(-8,-9,-10)])$class
    # misClassError.qda <- mean(y.train.val != y.predict.qda)
    
    results <- setNames(c((1-misClassError.knn),
                          (1-misClassError.nb),
                          (1-misClassError.rf),
                          (1-misClassError.svm),
                          (1-misClassError.lda)),
                        c("KNN", "NaiveBayes", "Random_Forest",
                          "SVM", "LDA"))
    
    # results <- setNames(c(mean_sensitivity.knn,mean_sensitivity.nb,
    #                       mean_sensitivity.rf, mean_sensitivity.svm,
    #                       mean_sensitivity.lda),
    #                     c("KNN", "NaiveBayes", "Random_Forest",
    #                       "SVM", "LDA"))
    return(results)
  })
  rowMeans(data.frame(cv))
  write.csv(rowMeans(data.frame(cv)), file ="temp_results.csv" )
  
  #Fit final model on complete training dataset
  rf <- randomForest(x = X1.n_p, y = as.factor(y.train), ntree = 300, 
                     replace = FALSE, proximity = TRUE) 
  
  
  #Preprocess test data
  X1.test.n_p = pre_process(X.test, plot = FALSE)
  
  rf_pred <- predict(rf, X1.test.n_p)
  write.csv(rf_pred,file="ECG_predictions_group_B_week_5.csv",row.names=FALSE)