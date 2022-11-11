
# Create directory to save results
savepath <- 
  function(savename) { 
    dir.create("results/", showWarnings = FALSE)
    result_folder <- paste0("results/", Sys.Date())
    dir.create(result_folder, showWarnings = FALSE)
    
    savename <-
      paste0(result_folder, "/", savename)
    
    
    return(savename)
    
  }

# Create directory to save intermediate data
savepath_data <- 
  function(folder, savename) { 
    result_folder <- paste0("data/processed/", folder)
    dir.create(result_folder, showWarnings = FALSE)
    
    savename <-
      paste0(result_folder, "/", savename)
    
    
    return(savename)
    
  }

# Compute PCA from wide data
do_pca <- 
  function(wide_data, npcs = NULL) {
    suppressMessages(require(tidyverse))
    suppressMessages(require(pcaMethods))
    
    if(is.null(npcs)) {
      npcs <- min(dim(wide_data))
    }
    
    wide_data %>% 
      t() %>% 
      pca(nPcs = npcs) 
  }

# Obtain PCA scores from PCA result
get_pca_scores <- 
  function(pca_res, 
           use_R2cum_PCselection = F,
           R2cum_lim = 0.8,
           use_sDev_PCselection = F) {
    suppressMessages(require(tidyverse))
    
    pc_lim <- ncol(pca_res@scores)
    
    if(use_R2cum_PCselection) {
      pc_lim <- 
        which(pca_res@R2cum > R2cum_lim)[1]
    } else if(use_sDev_PCselection) {
      pc_lim <- 
        rev(which(pca_res@sDev >= 1))[1]
    }
    
    pca_res %>% 
      scores() %>% 
      {.[,1:pc_lim]} 
  }

# Compute UMAP from wide data
do_umap <- 
  function(wide_data, 
           seed = 42, 
           n_neighbors = 15,
           n_components = 2, 
           ...) {
    suppressMessages(require(tidyverse))
    suppressMessages(require(uwot))
    
    set.seed(seed)
    
    umap_res <- 
      uwot::umap(wide_data, 
           n_neighbors = n_neighbors,
           n_components = n_components,
           ...) %>% 
      as.data.frame()
    
    rownames(umap_res) <- rownames(wide_data)
    colnames(umap_res) <- paste0("UMAP", 1:ncol(umap_res))
    
    umap_res
  }

# Adjust p-values (correct for multiple hypothesis testing)
DE_adjust <-
  function(volcano.result) {
    volcano.result %>% 
      as_tibble() %>% 
      group_by(Class) %>% 
      do({
        volcano.data <- .
        
        volcano.data %>%
          mutate(p.val=as.numeric(p.value)) %>%
          mutate(p.adjusted=p.adjust(p.val, method = 'fdr', n = nrow(volcano.data))) %>%
          mutate(sig=ifelse(p.adjusted<0.05 & difference<0, "significant down",
                            ifelse(p.adjusted<0.05 & difference>0, "significant up", 'not significant')))
        
      })
  }

# Function to impute missing values in the data
impute_values <- 
  function(data, ID, wide_data = F) {
    
    if(wide_data == F) {
      data_wide <- 
        data %>% 
        select(ID, Assay, NPX) %>% 
        spread(Assay,NPX) 
      
    } else {
      data_wide <- 
        data
    }
    
    data_imputed <- 
      data_wide %>% 
      column_to_rownames(ID) %>% 
      as.matrix() %>% 
      t() %>% 
      impute.knn() 
    
    final_data <- 
      data_imputed$data %>% 
      t() %>% 
      as_tibble(rownames = ID)
    
    return(final_data)
    
  }

# Generate volcano plots from differential expression results, includes multiple hypothesis testing correction
cancer_volcano <- function(cancer,volcano.result, title) {
  volcano.result.sig<-volcano.result[0,]
  
  volcano.data<-volcano.result %>%
    as_tibble() %>%
    filter(Class==cancer)
  
  volcano.data<-
    volcano.data %>%
    mutate(p.val=as.numeric(p.value)) %>%
    mutate(p.adjusted=p.adjust(p.val, method = 'fdr', n = nrow(volcano.data))) %>%
    mutate(sig=ifelse(p.adjusted<0.05 & difference<0, "significant down",
                      ifelse(p.adjusted<0.05 & difference>0, "significant up", 'not significant')))
  
  volcano.result.sig<-rbind(volcano.result.sig, volcano.data)
  
  top.sig.down<-volcano.data %>%
    filter(p.adjusted<0.05 & difference<0) %>%
    arrange(p.adjusted) %>%
    pull(OlinkID)
  
  top.sig.up<-volcano.data %>%
    filter(p.adjusted<0.05 & difference>0) %>%
    arrange(p.adjusted) %>%
    pull(OlinkID)
  
  top.sig.prot<-c(top.sig.up[1:40], top.sig.down[1:10])
  
  tab<- volcano.data %>%
    mutate(sig.label=ifelse(OlinkID %in% top.sig.prot, "top sig",0))
  
  num.sig.up<-length(top.sig.up)
  num.sig.down<-length(top.sig.down)
  
  volcano.data %>%
    ggplot(aes(x=difference, y=-log10(p.adjusted), color=sig, label=Assay)) +
    geom_point(size = 1, alpha = 0.4)+ 
    geom_text_repel(data=subset(tab, sig.label=="top sig"))+
    geom_hline(yintercept=-log10(0.05), linetype='dashed')+
    geom_vline(xintercept=0, linetype='dashed')+
    scale_color_manual(values=c("grey","blue","red"))+
    ggtitle(label = paste0(cancer, title),
            subtitle = paste0("Num significant up = ",num.sig.up, "\nNum significant down = ",num.sig.down))+
    themes$main +
    theme(legend.position="none", plot.subtitle = element_text(size=10, face="italic"))
}

