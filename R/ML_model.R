

ML_model <- function(df, balance = T, sample_size = 100000, train_size = 0.75,
    subset_celltype = F, celltype_col = celltype_col, celltype_name = celltype_name,
    alg = c('all', 'RF', 'XGB'), class_col = class_col, seed = seed, 
    label = label, allowParallel = TRUE, free_cores = 4){
    # Subset G and P
    df <- subset(df, CANARY %in% c('G', 'P'))

    # subset a cell type
    if(subset_celltype){
        df <- subset(df, )
        df <- df[df[,celltype_col] %in% celltype_name,]
    }
    

    # Sample balanced classes
    if(balance){
        k_G <- sample(which(df$CANARY == 'G'), sample_size/2)
        k_P <- sample(which(df$CANARY == 'P'), sample_size/2)
        df <- df[c(k_G,k_P),]
    }else{
        df <- dplyr::sample_n(df, size=sample_size)
    }
    
    # Training and test sets
    train_i <- caret::createDataPartition(y=CANARY, p=train_size, list=FALSE)
    df_train <- df[train_i,]
    df_test <- df[-train_i,]

    # Train model
    TrainModel(df_train, df_test, alg = alg, class_col = class_col, seed = seed, 
    label = label, allowParallel = allowParallel, free_cores = free_cores)


}

TrainModel <- function(TrainSet, TestSet, alg = c('all', 'RF', 'XGB'), 
    class_col = class_col, seed = 40, label = label, allowParallel = TRUE, 
    free_cores = 4){
    # name_0 == 'cells' |  name_1 == 'debris'/'beads'/'dead'
    # Set seed
    set.seed(seed)    
    
    alg <- match.arg(alg, c('all', 'RF', 'XGB'))

    # Ask for channels for training features   
    print(as.matrix(colnames(TrainSet)))
    prompt <- "Enter the column INDICES of the training features (separated by 
              single space only, no comas allowed) \n"
    features <- as.numeric(strsplit(readline(prompt), " ")[[1]])
    features <- colnames(TrainSet)[features]
    features <- features[order(features)]

    TrainSet <- t_asinh(TrainSet[c(class_col, features)])
    TestSet <- t_asinh(TestSet[c(class_col, features)])

    colnames(TrainSet) <- c(class_col, features)

    colnames(TestSet) <- colnames(TrainSet)


    if(alg == 'RF' || alg == 'all'){
        print('Running Random Forest')
        # Initiating parallelization
        if(allowParallel){
        cluster <- parallel::makeCluster(parallel::detectCores() - free_cores) 
        doParallel::registerDoParallel(cluster)
        }
        # Computing training control parameters
        fitControl <- caret::trainControl(
            method = 'repeatedcv',
            number = 10,
            repeats = 3,
            search = 'grid',
            savePredictions = 'final',
            classProbs = T,
            summaryFunction = twoClassSummary,
            allowParallel = allowParallel,
            returnData = FALSE)

        # Setting tune grid
        tunegrid <- expand.grid(.mtry=c(1:length(features)))

        # Fitting the model
        model_rf <- caret::train(TrainSet[,-1], TrainSet[,1],
            method ='rf',
            trControl = fitControl,
            metric = 'ROC',
            tuneGrid=tunegrid)

        # Assesing feature importance
        ftimp_rf <- caret::varImp(model_rf)

        # Test
        # Prediction on test dataset
        pred_rf <- predict(model_rf, TestSet)
        
        # Compute confusion matrix
        conf_rf <- caret::confusionMatrix(
            reference = as.factor(TestSet[,1]),
            data = pred_rf,
            mode = 'everything',
            positive = name_1)

        #save RData
        save(model_rf, ftimp_rf, pred_rf, conf_rf, file = paste0(label, '_RFmodel.RData'))
        print('Random Forest completed')
        
        # Finalizing parallelization
        if(allowParallel){
        parallel::stopCluster(cluster)
        foreach::registerDoSEQ()
        }
    }
    

    if(alg == 'XGB' || alg == 'all'){
        print('Running XGBoost')
        # Rearrangement
        X_train = xgboost::xgb.DMatrix(as.matrix(TrainSet[,-1]))
        y_train = TrainSet[,1]

        # Computing training control parameters
        xgb_trcontrol = caret::trainControl(
            method = 'repeatedcv',
            number = 10,
            repeats = 3,
            search = 'grid',
            savePredictions = 'final',
            summaryFunction = twoClassSummary,
            allowParallel = allowParallel,
            classProbs = TRUE,
            verboseIter = FALSE,
            returnData = FALSE)

        # Specifying grid space
        xgbGrid <- expand.grid(nrounds = c(100,200),  # this is n_estimators in the python code above
                               max_depth = c(10, 15, 20, 25),
                               colsample_bytree = seq(0.5, 0.9, length.out = 5),
                               ## The values below are default values in the sklearn-api.
                               eta = 0.1,
                               gamma=0,
                               min_child_weight = 1,
                               subsample = 1)

        # Fitting the model
        model_xgb = caret::train(
            X_train, y_train,
            trControl = xgb_trcontrol,
            tuneGrid = xgbGrid,
            method = "xgbTree",
            metric = 'ROC',
            nthread = detectCores() - free_cores)

        # Assesing feature importance
        ftimp_xgb <- caret::varImp(model_xgb)
        
        # Test
        X_test = xgboost::xgb.DMatrix(as.matrix(TestSet[,-1]))
        y_test = TestSet[,1]

        # Prediction on test dataset
        pred_xgb <- predict(model_xgb, X_test)
        # Compute confusion matrix
        conf_xgb <- caret::confusionMatrix(
            reference = as.factor(y_test),
            data = pred_xgb,
            mode = 'everything',
            positive = name_1)

        # Save RData
        save(model_xgb, ftimp_xgb, pred_xgb, conf_xgb, file = paste0(label, '_XGBmodel.RData'))
        print('XGBoost completed')
    }

}