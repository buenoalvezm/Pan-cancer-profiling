
# Function to balance the number of samples in the control group
generate_balanced_controls <- function(cancer, metadata, set_data, n_control_groups) {
  
  if(cancer %in% c("BRC","CVX","ENDC","OVC")) {
    
    female_samples <- 
      metadata %>% 
      filter(Sex == "Female") %>% 
      pull(Sample)
    
    set_data <- 
      set_data %>% 
      filter(Sample %in% female_samples)
  } else if (cancer == "PRC") {
    
    male_samples <- 
      metadata %>% 
      filter(Sex == "Male") %>% 
      pull(Sample)
    
    set_data <- 
      set_data %>% 
      filter(Sample %in% male_samples)
    
  } else {
    set_data <- 
      set_data
  }
  
  control_pool <- 
    set_data %>% 
    filter(GROUP != cancer)
  
  cancer_samples_train <- 
    set_data %>% 
    filter(GROUP == cancer,
           set == "train") %>% 
    pull(Sample) %>% 
    length()
  
  n_train <- 
    cancer_samples_train/n_control_groups 
  
  set.seed(213)
  train_set <- 
    control_pool %>%
    filter(set == "train") %>% 
    group_by(GROUP) %>% 
    sample_n(size = ceiling(n_train)) %>% 
    bind_rows(cancers_split %>% 
                filter(GROUP == cancer,
                       set == "train")) %>% 
    mutate(GROUP = case_when(GROUP == cancer ~ GROUP,
                             T ~ "Other"))
  
  cancer_samples_test <- 
    set_data %>% 
    filter(GROUP == cancer,
           set == "test") %>% 
    pull(Sample) %>% 
    length()
  
  n_test <- 
    cancer_samples_test/n_control_groups
  
  set.seed(213)
  test_set <- 
    control_pool %>%
    filter(set == "test") %>% 
    group_by(GROUP) %>% 
    sample_n(size = ceiling(n_test)) %>% 
    bind_rows(cancers_split %>% 
                filter(GROUP == cancer,
                       set == "test")) %>% 
    mutate(GROUP = case_when(GROUP == cancer ~ GROUP,
                             T ~ "Other"))
  
  train_set %>%
    bind_rows(test_set) 
  
}

# Function to predict disease using the caret package 
disease_classifier <- function(cancer_type, 
                               method, 
                               data, 
                               split_data,
                               metadata, 
                               select_proteins = NULL, 
                               tune_param = 10, 
                               CV = 0) {
  
  # Prepare the expression data
  exp_data <-
    data %>%
    inner_join(split_data, by = "Sample") %>% 
    mutate(Cancer_list=ifelse(GROUP == cancer_type, "1_Cancer", "0_Control")) %>%
    mutate(Cancer_list=as.factor(Cancer_list)) %>%
    arrange(Sample) %>%
    select(Sample,Assay,NPX,Cancer_list) %>% 
    spread(Assay,NPX,-1)%>%
    column_to_rownames("Sample")
  
  # If a set of protein is indicated, select only the given proteins
  if(length(select_proteins)>0) {
    exp_data <-
      exp_data %>%
      as_tibble(rownames = "Sample") %>%
      select(Sample,select_proteins, Cancer_list) %>%
      column_to_rownames("Sample")
  } else {
    exp_data <- exp_data
  }
  
  # Data imputation
  preproc <- preProcess(exp_data,
                        method = "knnImpute")
  
  imputed_data <- predict(preproc,
                          exp_data) %>%
    mutate(Cancer_list = as.factor(Cancer_list))
  
  # Select test and train sets
  train_samples <- 
    split_data %>% 
    filter(set == "train") %>% 
    pull(Sample)
  
  test_samples <- 
    split_data %>% 
    filter(set == "test") %>% 
    pull(Sample)
  
  train_data <- 
    imputed_data %>% 
    rownames_to_column("Sample") %>% 
    filter(Sample %in% train_samples) %>% 
    column_to_rownames("Sample")
  
  test_data <- imputed_data %>% 
    rownames_to_column("Sample") %>% 
    filter(Sample %in% test_samples) %>% 
    column_to_rownames("Sample")
  
  # Fit model
  t <- Sys.time()
  if(CV > 0) {
    set.seed(213)
    ctrl  <- trainControl(method  = "cv", number  = CV)
    model_res <- caret::train(Cancer_list ~ .,
                              data = train_data,
                              method = method,
                              trControl = ctrl,
                              tuneLength = tune_param)
  } else {
    set.seed(213)
    model_res <- caret::train(Cancer_list ~ .,
                              data = train_data,
                              method = method,
                              tuneLength = tune_param)
  }
  
  train_performance <- getTrainPerf(model_res)
  
  speed <- difftime(Sys.time(), t, units="secs") %>% as.numeric()
  
  pred <- predict(model_res, test_data)
  
  
  # Retrieve most important proteins (variable importance function)
  if (!method %in% c("knn", "lda", "svmLinear")) {
    important_prots <- varImp(model_res)
    
    prot_res <-
      important_prots$importance %>%
      as_tibble(rownames = "Protein") %>%
      arrange(-Overall)
    
    prot_plot <-
      prot_res %>%
      mutate(Protein = factor(Protein, levels = rev(prot_res$Protein))) %>%
      head(20) %>%
      ggplot(aes(x = Overall, y = Protein, fill = Overall)) +
      geom_bar(stat = "identity") +
      theme_bw()
  } else {
    prot_res <- ""
    prot_plot <- ggplot() + theme_void()
  }
  
  
  # Generate ROC cruve
  confusion_matrix <- caret::confusionMatrix(table(test_data$Cancer_list,pred))
  
  cm <- table(test_data$Cancer_list,pred) 
  
  cm_plot <-
    cm %>% 
    as.data.frame() %>% 
    rename(True_class = Var1,
           Predicted_class = pred) %>%
    mutate(color = case_when(True_class == Predicted_class ~ "#C5E0B3",
                             TRUE ~ "#FB8875")) %>% 
    ggplot(aes(Predicted_class,True_class, label = Freq)) +
    geom_tile(aes(fill = color), colour = "white") +
    geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
    scale_fill_identity() +
    theme_void()+
    theme_bw() + 
    theme(legend.position = "none") +
    ggtitle(paste(cancer_type, " - confusion matrix"))
  
  
  # Generate ROC curves
  
  prob <- predict(model_res, test_data, type = "prob")
  
  p <- test_data %>% 
    rownames_to_column("Sample") %>% 
    select(Sample, True_class = Cancer_list) %>% 
    mutate(Predicted_class = prob$`1_Cancer`) 
  
  probabilities_test <- 
    test_data %>% 
    rownames_to_column("Sample") %>% 
    select(Sample, True_class = Cancer_list) %>% 
    mutate(Probability_control = prob$`0_Control`,
           Probability_cancer = prob$`1_Cancer`) 
  
  prediction.probabilities <- prob$`1_Cancer`
  
  predicted.classes <- pred
  observed.classes <- test_data$Cancer_list
  
  accuracy <- mean(observed.classes == predicted.classes)
  error <- mean(observed.classes != predicted.classes)
  
  roc_res <- pROC::roc(observed.classes, prediction.probabilities)
  df<-data.frame(y=unlist(roc_res[2]), x=unlist(roc_res[3])) %>% 
    mutate(GROUP = cancer_type)
  
  roc_plot <- df %>% 
    mutate(x = 1-x) %>% 
    arrange(y) %>% 
    ggplot(aes(x, y)) + 
    geom_line(size = 1) + 
    geom_abline(intercept=0, slope=1, linetype="dashed") +
    xlab("False Positive rate (1-Specificity)") + 
    ylab("True Positive rate (Sensitivity)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 8)) +
    annotate("text", x = .75, y = .25, label = paste("AUC",roc_res[[9]] %>% round(3))) +
    ggtitle("ROC curve")
  
  # Generate plot showing the predicted probabilities (using stage 1 patients in the test set)  
  prob_plot <- 
    data.frame(response = roc_res$response, predicted = roc_res$predictor) %>% 
    ggplot(aes(response,predicted, color = response)) +
    geom_point(alpha = 0.4) +
    themes$main
  
  return(list(model_res = model_res,
              split_data = split_data,
              speed = speed,
              train_performance = train_performance,
              confusion_matrix = confusion_matrix,
              prot_res = prot_res,
              prot_plot = prot_plot,
              accuracy = accuracy,
              error = error,
              roc_res = roc_res,
              roc_plot = roc_plot,
              confusion_matrix = confusion_matrix,
              cm_plot = cm_plot,
              prob_plot = prob_plot,
              probabilities_test = probabilities_test
  ))
  
}

# Function to predict stage 1 colorectal & lung cancer using the caret package 
disease_classifier_stage1 <- function(cancer_type, 
                                      method, 
                                      data, 
                                      split_data,
                                      metadata, 
                                      metadata_stage,
                                      select_proteins = NULL, 
                                      tune_param = 10, 
                                      CV = 0) {
  
  # Prepare the expression data
  exp_data <-
    data %>%
    inner_join(split_data, by = "Sample") %>% 
    mutate(Cancer_list=ifelse(GROUP == cancer_type, "1_Cancer", "0_Control")) %>%
    mutate(Cancer_list=as.factor(Cancer_list)) %>%
    arrange(Sample) %>%
    select(Sample,Assay,NPX,Cancer_list) %>% 
    spread(Assay,NPX,-1)%>%
    column_to_rownames("Sample")
  
  # If a set of protein is indicated, select only the given proteins
  if(length(select_proteins)>0) {
    exp_data <-
      exp_data %>%
      as_tibble(rownames = "Sample") %>%
      select(Sample,select_proteins, Cancer_list) %>%
      column_to_rownames("Sample")
  } else {
    exp_data <- exp_data
  }
  
  # Data imputation
  preproc <- preProcess(exp_data,
                        method = "knnImpute")
  
  imputed_data <- predict(preproc,
                          exp_data) %>%
    mutate(Cancer_list = as.factor(Cancer_list))
  
  # Select test and train sets
  train_samples <- 
    split_data %>% 
    filter(set == "train") %>% 
    pull(Sample)
  
  test_samples <- 
    split_data %>% 
    filter(set == "test") %>% 
    pull(Sample)
  
  train_data <- 
    imputed_data %>% 
    rownames_to_column("Sample") %>% 
    filter(Sample %in% train_samples) %>% 
    column_to_rownames("Sample")
  
  test_data <- imputed_data %>% 
    rownames_to_column("Sample") %>% 
    filter(Sample %in% test_samples) %>% 
    column_to_rownames("Sample")
  
  stage1_samples <- 
    split_data %>% 
    select(Sample) %>% 
    left_join(metadata_stage %>% 
                select(Sample,Stage),
              by = "Sample") %>% 
    filter(Stage == "1") %>% 
    pull(Sample)
  
  test_data_stage1 <- 
    imputed_data %>% 
    rownames_to_column("Sample") %>% 
    filter(Sample %in% test_samples,
           Sample %in% stage1_samples) %>% 
    bind_rows(imputed_data %>% 
                rownames_to_column("Sample") %>% 
                filter(Sample %in% test_samples,
                       Cancer_list == "0_Control")) %>% 
    column_to_rownames("Sample")
  
  # Fit model
  t <- Sys.time()
  if(CV > 0) {
    set.seed(213)
    ctrl  <- trainControl(method  = "cv", number  = CV)
    model_res <- caret::train(Cancer_list ~ .,
                              data = train_data,
                              method = method,
                              trControl = ctrl,
                              tuneLength = tune_param)
  } else {
    set.seed(213)
    model_res <- caret::train(Cancer_list ~ .,
                              data = train_data,
                              method = method,
                              tuneLength = tune_param)
  }
  
  train_performance <- getTrainPerf(model_res)
  speed <- difftime(Sys.time(), t, units="secs") %>% as.numeric()
  
  # Retrieve most important proteins (using the varImp function in caret)
  if (!method %in% c("knn", "lda", "svmLinear")) {
    important_prots <- varImp(model_res)
    
    prot_res <-
      important_prots$importance %>%
      as_tibble(rownames = "Protein") %>%
      arrange(-Overall)
    
    prot_plot <-
      prot_res %>%
      mutate(Protein = factor(Protein, levels = rev(prot_res$Protein))) %>%
      head(20) %>%
      ggplot(aes(x = Overall, y = Protein, fill = Overall)) +
      geom_bar(stat = "identity") +
      theme_bw()
  } else {
    prot_res <- ""
    prot_plot <- ggplot() + theme_void()
  }
  
  # Prediction (using all test samples)
  pred <- predict(model_res, test_data)
  
  # Generate confusion matrix (using all test samples)
  confusion_matrix <- caret::confusionMatrix(table(test_data$Cancer_list,pred))
  
  cm <- table(test_data$Cancer_list,pred) 
  
  cm_plot <-
    cm %>% 
    as.data.frame() %>% 
    rename(True_class = Var1,
           Predicted_class = pred) %>%
    mutate(color = case_when(True_class == Predicted_class ~ "#C5E0B3",
                             TRUE ~ "#FB8875")) %>% 
    ggplot(aes(Predicted_class,True_class, label = Freq)) +
    geom_tile(aes(fill = color), colour = "white") +
    geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
    scale_fill_identity() +
    theme_void()+
    theme_bw() + 
    theme(legend.position = "none") +
    ggtitle(paste(cancer_type, " - confusion matrix"))
  
  
  # Generate ROC curves (using all test samples)
  prob <- predict(model_res, test_data, type = "prob")
  
  p <- 
    test_data %>% 
    rownames_to_column("Sample") %>% 
    select(Sample, True_class = Cancer_list) %>% 
    mutate(Predicted_class = prob$`1_Cancer`) 
  
  probabilities_test <- 
    test_data %>% 
    rownames_to_column("Sample") %>% 
    select(Sample, True_class = Cancer_list) %>% 
    mutate(Probability_control = prob$`0_Control`,
           Probability_cancer = prob$`1_Cancer`) 
  
  prediction.probabilities <- prob$`1_Cancer`
  
  predicted.classes <- pred
  observed.classes <- test_data$Cancer_list
  
  accuracy <- mean(observed.classes == predicted.classes)
  error <- mean(observed.classes != predicted.classes)
  
  roc_res <- pROC::roc(observed.classes, prediction.probabilities)
  df<-data.frame(y=unlist(roc_res[2]), x=unlist(roc_res[3])) %>% 
    mutate(GROUP = cancer_type)
  
  roc_plot <- df %>% 
    mutate(x = 1-x) %>% 
    arrange(y) %>% 
    ggplot(aes(x, y)) + 
    geom_line(size = 1) + 
    geom_abline(intercept=0, slope=1, linetype="dashed") +
    xlab("False Positive rate (1-Specificity)") + 
    ylab("True Positive rate (Sensitivity)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 8)) +
    annotate("text", x = .75, y = .25, label = paste("AUC",roc_res[[9]] %>% round(3))) +
    ggtitle("ROC curve")
  
  # Generate plot showing the predicted probabilities (using all test samples)
  prob_plot <- 
    data.frame(response = roc_res$response, predicted = roc_res$predictor) %>% 
    ggplot(aes(response,predicted, color = response)) +
    geom_point(alpha = 0.4) +
    themes$main
  
  
  # Prediction (using stage 1 patients in the test set)  
  pred_s1 <- predict(model_res, test_data_stage1)
  
  # Generate confusion matrix ((using stage 1 patients in the test set)  
  confusion_matrix_s1 <- caret::confusionMatrix(table(test_data_stage1$Cancer_list,pred_s1))
  
  cm_s1 <- table(test_data_stage1$Cancer_list,pred_s1) 
  
  cm_plot_s1 <-
    cm_s1 %>% 
    as.data.frame() %>% 
    rename(True_class = Var1,
           Predicted_class = pred_s1) %>%
    mutate(color = case_when(True_class == Predicted_class ~ "#C5E0B3",
                             TRUE ~ "#FB8875")) %>% 
    ggplot(aes(Predicted_class,True_class, label = Freq)) +
    geom_tile(aes(fill = color), colour = "white") +
    geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
    scale_fill_identity() +
    theme_void()+
    theme_bw() + 
    theme(legend.position = "none") +
    ggtitle(paste(cancer_type, " - confusion matrix"))
  
  # Generate ROC curves (using stage 1 patients in the test set)  
  prob_s1 <- predict(model_res, test_data_stage1, type = "prob")
  
  p_s1 <- test_data_stage1 %>% 
    rownames_to_column("Sample") %>% 
    select(Sample, True_class = Cancer_list) %>% 
    mutate(Predicted_class = prob_s1$`1_Cancer`) 
  
  probabilities_test_s1 <- 
    test_data_stage1 %>% 
    rownames_to_column("Sample") %>% 
    select(Sample, True_class = Cancer_list) %>% 
    mutate(Probability_control = prob_s1$`0_Control`,
           Probability_cancer = prob_s1$`1_Cancer`) 
  
  prediction.probabilities_s1 <- prob_s1$`1_Cancer`
  
  predicted.classes_s1 <- pred_s1
  observed.classes_s1 <- test_data_stage1$Cancer_list
  
  accuracy_s1 <- mean(observed.classes_s1 == predicted.classes_s1)
  error_s1 <- mean(observed.classes_s1 != predicted.classes_s1)
  
  roc_res_s1 <- pROC::roc(observed.classes_s1, prediction.probabilities_s1)
  df_s1<-data.frame(y=unlist(roc_res_s1[2]), x=unlist(roc_res_s1[3])) %>% 
    mutate(GROUP = cancer_type)
  
  roc_plot_s1 <- df_s1 %>% 
    mutate(x = 1-x) %>% 
    arrange(y) %>% 
    ggplot(aes(x, y)) + 
    geom_line(size = 1) + 
    geom_abline(intercept=0, slope=1, linetype="dashed") +
    xlab("False Positive rate (1-Specificity)") + 
    ylab("True Positive rate (Sensitivity)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title = element_text(size = 8)) +
    annotate("text", x = .75, y = .25, label = paste("AUC",roc_res_s1[[9]] %>% round(3))) +
    ggtitle("ROC curve")
  
  # Generate plot showing the predicted probabilities (using stage 1 patients in the test set)  
  prob_plot_s1 <- 
    data.frame(response = roc_res_s1$response, predicted = roc_res_s1$predictor) %>% 
    ggplot(aes(response,predicted, color = response)) +
    geom_point(alpha = 0.4) +
    themes$main
  
  
  return(list(model_res = model_res,
              split_data = split_data,
              speed = speed,
              train_performance = train_performance,
              confusion_matrix = confusion_matrix,
              prot_res = prot_res,
              prot_plot = prot_plot,
              accuracy = accuracy,
              error = error,
              roc_res = roc_res,
              roc_plot = roc_plot,
              confusion_matrix = confusion_matrix,
              cm_plot = cm_plot,
              prob_plot = prob_plot,
              probabilities_test = probabilities_test,
              accuracy_s1 = accuracy_s1,
              error_s1 = error_s1,
              roc_res_s1 = roc_res_s1,
              roc_plot_s1 = roc_plot_s1,
              confusion_matrix_s1 = confusion_matrix_s1,
              cm_plot_s1 = cm_plot_s1,
              prob_plot_s1 = prob_plot_s1,
              probabilities_test_s1 = probabilities_test_s1
  ))
  
}

# Function to predict multiple cancers simultaneously
multiclassification <- function(data, 
                                split_data, 
                                metadata = metadata,
                                palette,
                                select_proteins = NULL,
                                method = "glmnet", 
                                groups, 
                                tune_param = 10, 
                                CV = 0) {
  
  # Prepare the expression data
  exp_data <-
    data %>%
    inner_join(split_data, by = "Sample") %>% 
    filter(GROUP %in% groups) %>% 
    mutate(GROUP=factor(GROUP, levels = groups)) %>%
    arrange(Sample) %>%
    select(Sample,Assay,NPX,GROUP) %>% 
    spread(Assay,NPX,-1,-4)%>%
    column_to_rownames("Sample")
  
  # If a set of protein is indicated, select only the given proteins
  if(length(select_proteins)>0) {
    exp_data <-
      exp_data %>%
      as_tibble(rownames = "Sample") %>%
      select(Sample,select_proteins, GROUP) %>%
      column_to_rownames("Sample")
  } else {
    exp_data <- exp_data
  }
  
  # Data imputation
  preproc <- preProcess(exp_data,
                        method = "knnImpute")
  
  imputed_data <- predict(preproc,
                          exp_data) %>%
    mutate(GROUP = as.factor(GROUP))
  
  # Select test and train sets
  train_samples <- 
    split_data %>% 
    filter(set == "train") %>% 
    pull(Sample)
  
  test_samples <- 
    split_data %>% 
    filter(set == "test") %>% 
    pull(Sample)
  
  train_data <- 
    imputed_data %>% 
    rownames_to_column("Sample") %>% 
    filter(Sample %in% train_samples) %>% 
    column_to_rownames("Sample")
  
  test_data <- imputed_data %>% 
    rownames_to_column("Sample") %>% 
    filter(Sample %in% test_samples) %>% 
    column_to_rownames("Sample")
  
  
  # Fit model
  t <- Sys.time()
  if(CV > 0) {
    set.seed(213)
    ctrl  <- trainControl(method  = "cv", number  = CV)
    model_res <- caret::train(GROUP ~ .,
                              data = train_data,
                              method = method,
                              trControl = ctrl,
                              tuneLength = tune_param)
  } else {
    set.seed(213)
    model_res <- caret::train(GROUP ~ .,
                              data = train_data,
                              method = method,
                              tuneLength = tune_param)
  }
  
  speed <- difftime(Sys.time(), t, units="secs") %>% as.numeric()
  
  # Prediction
  pred <- predict(model_res, test_data)
  
  # Retrieve most important proteins (using the varImp function in caret)
  important_prots <- varImp(model_res)
  
  prot_res <-
    important_prots$importance %>%
    as_tibble(rownames = "Protein")
  
  pal <- palette[groups]
  
  prot_plot <-
    important_prots$importance  %>%
    as_tibble(rownames = "Protein") %>% 
    gather(Cancer, Overall, -Protein) %>% 
    filter(Overall > 0) %>% 
    mutate(Cancer = factor(Cancer, levels = groups)) %>% 
    ggplot(aes(x = Cancer, y = Overall, color = Cancer, label = Protein)) +
    geom_point(alpha = 0.5, show.legend = F) +
    geom_text_repel(size = 3, show.legend = F) +
    scale_color_manual(values = pal) +
    themes$simple +
    ggtitle("Protein importance")
  
  # Generate confusion matrix
  confusion_matrix <- confusionMatrix(table(test_data$GROUP,pred))
  
  cm <- table(test_data$GROUP,pred) 
  
  cm_plot <- 
    cm %>% 
    as.data.frame() %>% 
    rename(True_class = Var1,
           Predicted_class = pred) %>%
    mutate(color = case_when(True_class == Predicted_class ~
                               pal[as.character(True_class)],
                             TRUE ~ paste("grey",90-2*Freq,sep=""))) %>% 
    mutate(color = case_when(color == "grey90" ~ "white",
                             TRUE ~ color)) %>% 
    ggplot(aes(Predicted_class,True_class, label = Freq)) +
    geom_tile(aes(fill = color), alpha = 0.8, colour = "white") +
    geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 0.5, hjust = 0.5) +
    scale_fill_identity()+
    theme_void()+
    themes$simple + 
    theme(legend.position = "none") +
    ggtitle("Confusion matrix") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # Generate ROC/PR curves
  prob <- predict(model_res, test_data, type = "prob")
  
  dat <-
    test_data %>% 
    rownames_to_column("Sample") %>% 
    select(Sample, GROUP) %>% 
    mutate(value = 1) %>% 
    spread(GROUP,value, fill= 0) 
  
  true_dat <- 
    dat %>% 
    set_names(paste(names(dat), "_true", sep = "")) %>%
    rename(Sample = `Sample_true`)
  
  dat_prob <- 
    prob %>% 
    rownames_to_column("Sample")
  
  prob_data <- 
    dat_prob %>% 
    set_names(paste(names(dat_prob), "_pred_glmnet", sep = "")) %>% 
    rename(Sample = Sample_pred_glmnet)
  
  final_df <- 
    true_dat %>% 
    left_join(prob_data, by = "Sample") %>% 
    select(-Sample)
  
  roc_res <- multi_roc(final_df, force_diag=T)
  pr_res <- multi_pr(final_df, force_diag=T)
  
  plot_roc_df <- plot_roc_data(roc_res)
  plot_pr_df <- plot_pr_data(pr_res)
  
  roc_dat <- 
    plot_roc_df %>%
    filter(!Group %in% c("Macro","Micro")) %>% 
    mutate(Performance = paste(Group, ": ", round(AUC, 4), sep = "")) %>% 
    arrange(-AUC)
  
  auc_pal <-
    roc_dat   %>% 
    select(Performance, Group) %>% 
    distinct() %>% 
    left_join(enframe(pal[groups], "Group", "color"), by = "Group") %>% 
    mutate(Performance = factor(Performance, levels = unique(roc_dat$Performance))) %>% 
    arrange(Performance) %>% 
    select(-Group) %>% 
    deframe()
  
  roc_plot <- 
    roc_dat %>% 
    mutate(Performance = factor(Performance, levels = unique(roc_dat$Performance))) %>% 
    arrange(Performance) %>% 
    ggplot(aes(x = 1-Specificity, y=Sensitivity)) +
    geom_path(aes(color = Performance), size=1) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                 colour='grey', linetype = 'dotdash') +
    themes$main + 
    scale_color_manual(values = auc_pal) 
  
  pr_dat <- 
    plot_pr_df %>%
    filter(!Group %in% c("Macro","Micro")) %>% 
    mutate(Performance = paste(Group, ": ", round(AUC, 4), sep = "")) %>% 
    arrange(-AUC)
  
  auc_pr_pal <-
    pr_dat   %>% 
    select(Performance, Group) %>% 
    distinct() %>% 
    left_join(enframe(pal[groups], "Group", "color"), by = "Group") %>% 
    mutate(Performance = factor(Performance, levels = unique(pr_dat$Performance))) %>% 
    arrange(Performance) %>% 
    select(-Group) %>% 
    deframe()
  
  pr_plot <- 
    pr_dat %>% 
    mutate(Performance = factor(Performance, levels = unique(pr_dat$Performance))) %>% 
    arrange(Performance) %>% 
    ggplot(aes(x = Recall, y=Precision)) +
    geom_path(aes(color = Performance), size=1) +
    themes$main + 
    scale_color_manual(values = auc_pr_pal)
  
  # Heatmap summarizing model probabilities for each sample
  annotation_r <-
    rownames(prob) %>% 
    enframe(value = "Sample") %>% 
    left_join(metadata %>% 
                select(Sample,GROUP),
              by = "Sample") %>% 
    select(-name) %>% 
    column_to_rownames("Sample")
  
  
  annotation_pal <- pal[groups]
  
  sample_levels <- 
    prob %>% 
    as_tibble(rownames = "Sample") %>% 
    gather(pred_cancer, prob, -Sample)  %>% 
    left_join(metadata %>% 
                select(Sample,GROUP),
              by = "Sample") %>% 
    group_by(Sample) %>% 
    top_n(1, prob) %>% 
    ungroup() %>% 
    mutate(GROUP = factor(GROUP, groups)) %>% 
    arrange(GROUP,-prob) %>% 
    select(Sample,GROUP)
  
  
  prob_heatmap <-
    prob %>% 
    as_tibble(rownames = "Sample") %>% 
    mutate(Sample = factor(Sample, levels = sample_levels$Sample)) %>% 
    arrange(Sample) %>% 
    column_to_rownames("Sample") %>% 
    pheatmap(clustering_method = "ward.D2",
             show_rownames = F,
             cluster_rows = F,
             annotation_row = annotation_r,
             color = colorRampPalette(c("white","grey30"))(100),
             annotation_colors = list(GROUP = annotation_pal)) %>% 
    as.ggplot()
  
  return(list(model_res = model_res,
              speed = speed,
              confusion_matrix = confusion_matrix,
              prot_res = prot_res,
              prot_plot = prot_plot,
              roc_res = roc_res,
              roc_df = roc_dat,
              roc_plot = roc_plot,
              pr_res = pr_res, 
              pr_df = pr_dat,
              pr_plotv= pr_plot,
              confusion_matrix = confusion_matrix,
              cm_plot = cm_plot,
              prob_data = prob_data, 
              pr_plot = pr_plot,
              prob_heatmap = prob_heatmap))
}
