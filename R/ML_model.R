equal_sampling <- function(df, group_col = 'CANARY', ptID_col = 'pt_ID', 
    groups = c('G', 'P')){
    
    # Subset groups
    df <- df[which(df[,group_col] %in% groups),]

    ptids <- unique(df[,ptID_col])

    # Get min number of cells among samples
    k <- min(table(df[,ptID_col]))

    idx <- c()
    for(i in ptids){
        idx <- c(idx, sample(which(df[,ptID_col] == i), k))
    }
    df <- df[idx,]

    df
}



ML_model <- function(df, 
                        balance = T, 
                        train_size = 0.75,
                        subset_celltype = F, 
                        celltype_col = celltype_col, 
                        celltype_name = celltype_name,
                        ask_features = T, 
                        ft_idxs, 
                        alg = c('all', 'RF', 'XGB'), 
                        class_col = class_col, 
                        seed = seed, 
                        label = label, 
                        allowParallel = TRUE, 
                        workers = 4){

    # subset a cell type
    if(subset_celltype){
        df <- df[df[,celltype_col] %in% celltype_name,]
    }

    if(ask_features){
        # Ask for channels for training features  
        print(as.matrix(colnames(df)))
        prompt <- "Enter the column INDICES of the training features (separated by 
                  single space only, no comas allowed) \n"
        ft_idxs <- as.numeric(strsplit(readline(prompt), " ")[[1]])
    }

    ft_idxs <- colnames(df)[ft_idxs]
    ft_idxs <- ft_idxs[order(ft_idxs)]
    df <- df[c(class_col, ft_idxs)]

    # Verifying class_col name
    class_idx <- grep(class_col, colnames(df))
    colnames(df)[class_idx] <- 'CANARY'
    df$CANARY <- as.factor(df$CANARY)
    
    # Training and test sets
    set.seed(seed+55)
    train_i <- caret::createDataPartition(df$CANARY, p=train_size, list=FALSE)
    df_train <- df[train_i,]
    df_test <- df[-train_i,]
    message('Partition done!')

    # Balance classes with SMOTE
    set.seed(seed+20)
    df_train <- DMwR::SMOTE(CANARY ~ ., data  = df_train)
    message('SMOTE done!')

    # Train model
    TrainModel(df_train, df_test, alg = alg, class_col = 'CANARY', seed = seed, 
    label = label, allowParallel = allowParallel, workers = workers)


}

TrainModel <- function(TrainSet, TestSet, alg = c('all', 'RF', 'XGB'), 
    class_col = class_col, seed = 40, label = label, allowParallel = TRUE, 
    workers = 4){
    library(caret)

    # Set seed
    set.seed(seed)    
    
    alg <- match.arg(alg, c('all', 'RF', 'XGB'))

    TrainSet <- denoisingCTF::t_asinh(TrainSet)
    TestSet <- denoisingCTF::t_asinh(TestSet)
    message('Arcsin transformation done!')


    if(alg == 'RF' || alg == 'all'){
        print('Running Random Forest')
        # Initiating parallelization
        if(allowParallel){
        cluster <- parallel::makeCluster(workers) 
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
        tunegrid <- expand.grid(.mtry=c(1:ncol(TrainSet)-1))

        # Fitting the model
        model_rf <- caret::train(TrainSet[,-1], TrainSet[,1],
            method ='rf',
            trControl = fitControl,
            metric = 'ROC',
            tuneGrid=tunegrid)
        message('Model training done!')

        # Assesing feature importance
        ftimp_rf <- caret::varImp(model_rf)

        # Test
        # Prediction on test dataset
        pred_rf <- predict(model_rf, TestSet)
        
        # Compute confusion matrix
        conf_rf <- caret::confusionMatrix(
            reference = as.factor(TestSet[,1]),
            data = pred_rf,
            mode = 'everything')

        message('Testing done!')

        #save RData
        save(model_rf, ftimp_rf, pred_rf, conf_rf, TrainSet, TestSet,
            file = paste0(label, '_RFmodel.RData'))
        message('Random Forest completed')
        
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
            nthread = workers)

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
            mode = 'everything')

        # Save RData
        save(model_xgb, ftimp_xgb, pred_xgb, conf_xgb, 
            TrainSet, TestSet, file = paste0(label, '_XGBmodel.RData'))
        print('XGBoost completed')
    }

}