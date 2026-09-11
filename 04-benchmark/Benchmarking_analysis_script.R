# CORRE Guillaume @ GENETHON
# 2026-08-14
# Benchmark analysis of GENETHOFF, GUIDE-seq (v1.0.2 / v2) and iGUIDE-seq (v1.2.0)

# R 4.6.1

#---------------------------------------------
# Load libraries ----
#---------------------------------------------

source("rbo_functions.R")

library(tidyverse)
library(GenomicRanges)
library(patchwork)
library(ggprism)
library(UpSetR)
library(ComplexUpset)
library(ggvenn)
library(pheatmap)
library(gridExtra)
library(grid)

slop_window = 10 # use to cluster cleavage sites across pipelines

#---------------------------------------------##
# Run time analysis----
#---------------------------------------------##

runtime <- read.delim("runtime_benchmark.csv", sep =";", check.names = F)
runtime <- runtime %>% mutate(tools = factor(tools,levels = c("GENETHOFF" ,"GUIDE-seq_v1" ,"GUIDE-seq_v2" , "iGUIDE-seq")))


runtime_raw <- ggplot(runtime, aes(total_threads, as.numeric(real_min), 
                                   fill = tools)) +
  geom_point(size = 3, pch = 21, col = "black")+
  ggprism::theme_prism(base_size = 12)+
  stat_summary(
    fun = mean,
    geom = "line",
    aes(group = tools, col = tools),
    linetype = 2,
    linewidth = 1, show.legend = F
  )+
  labs(x = "Provided threads", y = "Runtime (min)")+
  scale_x_continuous(breaks = c(6,12,24))+
  scale_y_continuous()

scalability <- ggplot(runtime, aes(total_threads,`user+sys`/as.numeric(real_min), fill = tools)) +
  geom_point(size = 3, pch = 21, col = "black")+
  geom_abline(lty = 2, col = "grey") + 
  ggprism::theme_prism(base_size = 12)+
  stat_summary(
    fun = mean,
    geom = "line",
    aes(group = tools, col = tools),
    linetype = 2,
    linewidth = 1, show.legend = F
  )+
  labs(x = "Provided threads", y = "Estimated threads")+
  scale_x_continuous(breaks = c(6,12,24))+
  scale_y_continuous(breaks = c(0,6,12,18,24,30))
  

#---------------------------------------------##
# Load each pipeline results table ----
#---------------------------------------------##

  #---------------------------------------------##
  ## GUIDESEQ v1----
  #---------------------------------------------##
  
  # list files of interest 

  files <- list.files("GUIDE-seq_v1.0.2/identified/", full.names = T)
  names(files) <- str_replace(str_remove(basename(files),"_identifiedOfftargets.txt"),"-","_")
  files <- grep(files, pattern = "control", value = T, invert = T)
  
  # read them all 
  guideseq <- lapply(files, read.delim)
  
  guideseq_v1 <- lapply(guideseq, function(x){
    
    x %>% 
      mutate(Chromosome = as.character(Chromosome))%>% 
      select(Chromosome, Position, bi.sum.mi) 
  })
  
  rm(files);
  

  # merge replicates for each gRNA
  targets <- names(guideseq_v1)
  targets_unique <- str_extract(targets,"^[A-Za-z0-9]+") %>% unique()
  
  guideseq_v1_aggregated <- list()
  
  for(i in targets_unique){
    cat("processing",i,"\n")
    
    guideseq_v1_aggregated[[i]] <- bind_rows(guideseq_v1[str_starts(names(guideseq_v1),i)])%>% 
      group_by(Chromosome, Position) %>% 
      summarise(counts = sum(bi.sum.mi))
  }
  

  guideseq_v1 <- bind_rows(guideseq_v1_aggregated,.id = "gRNA")
  
  
  # clean objects
  rm(targets); rm(guideseq); rm(guideseq_v1_aggregated)
  
  #---------------------------------------------##
  ## GUIDESEQ v2 ----
  #---------------------------------------------##
  
  # list files of interest 
    
  files <- list.files("GUIDE-seq_v2/identified/", full.names = T, pattern = "_identifiedOfftargets.txt$")
  files <- grep(files,pattern = "Control",value = T, invert = T)
  names(files) <- str_replace(str_remove(basename(files),"_identifiedOfftargets.txt"),"-","_")

  # read them all 
  guideseq <- lapply(files, read.delim)

  guideseq_v2 <- lapply(guideseq, function(x){
    
    x %>% 
      mutate(Chromosome = as.character(X.BED_Chromosome))%>% 
      select(Chromosome, Position, bi.sum.mi)  
    
  })

  rm(files); rm(guideseq)

  
  # merge replicates for each gRNA
  targets <- names(guideseq_v2)
  targets_unique <- str_extract(targets,"^[A-Za-z0-9]+") %>% unique()
  
  guideseq_v2_aggregated <- list()
  
  for(i in targets_unique){
    
    guideseq_v2_aggregated[[i]] <- bind_rows(guideseq_v2[str_starts(names(guideseq_v2),i)])%>%
      group_by(Chromosome, Position) %>% 
      summarise(counts = sum(bi.sum.mi))
  }
  
  guideseq_v2 <- bind_rows(guideseq_v2_aggregated,.id = "gRNA")
  
  
  # clean objects
  rm(targets); rm(guideseq); rm(guideseq_v2_aggregated)
  
  


  
  #---------------------------------------------##
  ## GENETHOFF ----
  #---------------------------------------------##
  path <- "GENETHOFF/results/paper_review.xlsx"
  
  genethoff <- path %>% 
    readxl::excel_sheets() %>% 
    set_names() %>% 
    map(readxl::read_excel, path = path)
  
  rm(path)
  
  
  # remove sites without gRNA match
  
  genethoff_sub <- lapply(genethoff, function(x){
    x %>%
      filter(!is.na(Alignment),
             abs(relative_distance) < 15,
             PAM_indel_count <=1) %>%
      select(chromosome, N_UMI_cluster,cut_modal_position)  
    
  }
  )
  
  
  # merge libraries for each gRNA
  targets <- names(genethoff_sub)
  targets_unique <- str_extract(targets,"^[A-Za-z0-9]+") %>% unique()
  targets_unique <- grep(targets_unique,pattern = "Mock",ignore.case = T,invert = T, value = T)
  
  genethoff_aggregated <- list()
  
  for(i in targets_unique){
  
  genethoff_aggregated[[i]] <- bind_rows(genethoff_sub[str_starts(names(genethoff_sub),i)])%>% 
    group_by(Chromosome = chromosome, Position = cut_modal_position) %>% 
    summarise(counts = sum(N_UMI_cluster))
  }
  
  
  genethoff <- bind_rows(genethoff_aggregated,.id = "gRNA")

  rm(targets); rm(genethoff_aggregated); rm(genethoff_sub)

  
  #---------------------------------------------##
  ## iGUIDE-seq ----
  #---------------------------------------------##
  path <- "iGUIDE-seq_v1.2.0/manual_data_extraction.xlsx"
  
  iguide <- path %>% 
    readxl::excel_sheets() %>% 
    set_names() %>% 
    map(readxl::read_excel, path = path)
  
  rm(path)
  
  
  iguide_aggregated <- lapply(iguide, function(x){
    
    x %>% select("Edit Site", counts="Abund.") %>%
      separate(convert = T,"Edit Site", into = c("Chromosome","strand","Position"), sep =":") %>% 
      select(-strand)
  })
  
  iguideseq <- bind_rows(iguide_aggregated,.id = "gRNA")

  rm(iguide_aggregated);rm(iguide)
  
  
#---------------------------------------------##
# Aggregate the 4 pipelines results ----
#---------------------------------------------##
  
  all_gRNAs_pipelines <- bind_rows("GUIDE-seq_v1" = guideseq_v1,
                   "GUIDE-seq_v2"= guideseq_v2,
                   "iGUIDE-seq" = iguideseq,
                   "GENETHOFF"= genethoff, .id = "workflow")
  
  

#---------------------------------------------##
# Cluster cleavage sites across pipelines ----
#---------------------------------------------##
  
    #---------------------------------------------##
    ## Reformat columns ----
    #---------------------------------------------##
  
  all_gRNAs_pipelines <- all_gRNAs_pipelines %>%
    mutate(start = Position ,
           end = Position , 
           Chromosome = str_remove_all(Chromosome,"chr"),
           gRNA = case_when(gRNA=="TRAC5"~ "TRAC",
                            TRUE ~ gRNA),
           workflow = factor(workflow))
  
  
    #---------------------------------------------##
    ## Convert to gRange object for clustering ----
    #---------------------------------------------##
  all_gRNAs_pipelines_grange <- makeGRangesFromDataFrame(df = all_gRNAs_pipelines, 
                                                         keep.extra.columns = T,
                                                         seqnames.field = "Chromosome",
                                                         start.field = "start",
                                                         end.field = "end")

  
    #---------------------------------------------##
    ## Cluster all cleavage sites ----
    #---------------------------------------------##
      
    # find overlapping regions in a window
    window = 10
    
    hits <- findOverlaps(all_gRNAs_pipelines_grange, reduce(all_gRNAs_pipelines_grange, min.gapwidth = window))
    cluster_id <- subjectHits(hits)[order(queryHits(hits))]

    # annotate cleavage sites with cluster ID
    all_gRNAs_pipelines$cluster <- cluster_id


    
    #---------------------------------------------##
    # Get relative abundance of cleavage sites per gRNA and cluster ----
    #---------------------------------------------##
    
    # calculate relative abundance of cleavage site per gRNA and method
    all_gRNAs_pipelines <- all_gRNAs_pipelines %>% 
      group_by(workflow,gRNA, Chromosome, cluster) %>% 
      summarise(start = min(start),
                end = max(end),
                counts = sum(counts)) %>% 
      group_by(workflow,gRNA) %>% 
      mutate(prop = counts / sum(counts) * 100) %>% 
      unite(col = "Position",Chromosome,start,end,remove = F,sep = "_")
    
    
    
    
    #---------------------------------------------##
    # Annotate cluster if they are the on target expected site ----
    #---------------------------------------------##
    
    all_gRNAs_pipelines <- all_gRNAs_pipelines %>% 
      mutate(OT = case_when(gRNA=="B2M" & Position %in% c("15_44711567_44711568" ,"15_44711569_44711569","15_44711568_44711568") ~ T,
                                     gRNA=="TRAC" & Position %in% c("14_22547664_22547664","14_22547663_22547664") ~ T,
                                     gRNA=="VEGFAs2" & Position %in% c("6_43770825_43770825","6_43770824_43770824") ~ T,
                                     gRNA=="VEGFAs3" & Position %in%c("6_43769732_43769732", "6_43769733_43769733") ~ T,
                                     TRUE ~ FALSE)
             ) %>% 
      mutate(workflow = factor(workflow, levels = c("GENETHOFF" ,"GUIDE-seq_v1" ,"GUIDE-seq_v2" , "iGUIDE-seq")))

    
    #---------------------------------------------##
    # Pivot table to get pipeline in columns ----
    #---------------------------------------------##
    
    ## add method OT positions
    all_gRNAs_pipelines_wide <- all_gRNAs_pipelines %>% 
      pivot_wider(names_from = "workflow", values_from = c(Position,counts, prop),names_glue ="{workflow}_{.value}", id_cols = c("gRNA","cluster") )
    

#---------------------------------------------##
# Save tables  ----
#---------------------------------------------##

write.table(all_gRNAs_pipelines, paste("complete_OT_table_all_methods_",window,"bp.csv",sep=""), sep=";",row.names = F, quote = F)

write.table(all_gRNAs_pipelines_wide, paste("complete_OT_table_all_methods_wide_",window,"bp.csv",sep=""), sep=";",row.names = F, quote = F)





#---------------------------------------------##
# Make some plots ----
#---------------------------------------------##

  #---------------------------------------------##
  ## Reload datasets ----
  #---------------------------------------------##

  #reload tables if necessary 
  all_gRNAs_pipelines <- read.delim("complete_OT_table_all_methods_10bp.csv", header = T, sep =";")
  all_gRNAs_pipelines_wide <- read.delim("complete_OT_table_all_methods_wide_10bp.csv", header = T, sep =";")
  


  #---------------------------------------------##
  ## Rank-abundance plot ----
  #---------------------------------------------##

rank_ab <- ggplot(all_gRNAs_pipelines %>%
                    filter(counts >0) %>% 
         group_by(gRNA,workflow) %>% 
         mutate(rank = row_number(-counts)), 
       aes(rank,counts, col = workflow)) + 
  geom_line(show.legend = F) +
  geom_point(show.legend = F)+
  facet_wrap(~gRNA, nrow = 2, scale = "free_x") +
  scale_y_log10(breaks = c(1,10,100,1000,10000,100000))+
  scale_x_log10(breaks = c(1,10,100,1000), limits = c(1,NA))+
  ggprism::theme_prism(base_size = 12, border = T) +
  geom_point(data = . %>% filter(OT==T),col = "black", shape = 21, size = 1,stroke = 2,show.legend = F)+
  labs(x = "Rank (descending abundance)", y = "Abundance" ) 

rank_ab


  #---------------------------------------------##
  ## Rank of top cleavage sites ----
  #---------------------------------------------##
  
  # keep site that are in the 4 pipelines
  keep <- all_gRNAs_pipelines %>%
  ungroup %>% 
    count(gRNA,cluster) %>%
    filter(n>=4)

  ggplot(data = all_gRNAs_pipelines %>% 
           group_by(gRNA,workflow) %>% 
           mutate(rank = row_number(-counts)) %>% # calculate rank per gRNA and pipeline
           semi_join(keep) %>%                    # keep sites present in 4 pipelines
           filter(#rank<= 25,
                  gRNA %in% c("VEGFAs2","VEGFAs3")), 
         aes(workflow,rank, col = factor(cluster))) + 
    geom_point(show.legend = c(fill=F),pch = 21,col = "black", aes(size = prop, fill = factor(cluster)), alpha = 0.6) + 
    geom_line(aes(group = cluster), show.legend = F)+
    facet_wrap(~gRNA, scale="free")+
    ggprism::theme_prism(base_size = 12,border = T,axis_text_angle = 45) +
    labs(x = NULL, y = "Cleavage site Rank")
  


  
  #---------------------------------------------##
  # Venn diagrams and upset plots ----
  #---------------------------------------------##


# set parameters

margins <- c(0,0,0,0)
text_size = 4


  venn_plot <- list()
  upset_plot <- list()

  for(i in c("B2M","TRAC","VEGFAs2", "VEGFAs3")){
    
    cat("processing",i,"\n")
    
    df <- all_gRNAs_pipelines %>% filter(gRNA == i )
    
    # make a list
    lst <- split(df$cluster,f=df$workflow)
    
    # make the venn
    venn_plot[[i]] <- ggvenn::ggvenn(lst, 
                                     set_name_color = NA,
                                     fill_alpha = 0.8,
                                     fill_color = scale_fill_hue()$palette(4),
                                     set_name_size =3,
                                     text_size = text_size,
                                     show_percentage = FALSE) +
      ggtitle(i) + 
      ggprism::theme_prism(base_size = 12)+
      theme(plot.margin = unit(margins,units = "cm"), 
            axis.line = element_blank(),
            axis.title = element_blank(), 
            axis.text = element_blank(),
            axis.ticks = element_blank())
    
    
    
    # make the upset plot
    df <- fromList(lst)
    
    upset_plot[[i]]  <- ComplexUpset::upset(
      df,
      intersect = colnames(df),
      min_size=1,
      width_ratio=0.2,name = NULL,height_ratio = 0.6, 
      stripes='white',
      base_annotations=list(
        'Cleavage sites'= ComplexUpset::intersection_size(
          counts=TRUE,
          fill = "black")
      ),
      set_sizes = upset_set_size(
        geom = ggplot2::geom_bar(fill = "black")
      )+
        labs(y="Cleavage sites")
    ) 
  }

  #---------------------------------------------##
  ## Make the final montage ----
  #---------------------------------------------##

  venn_plots <- cowplot::as_grob((venn_plot$B2M + venn_plot$TRAC) / (venn_plot$VEGFAs2 + venn_plot$VEGFAs3))

x11();
(runtime_raw / scalability |plot_spacer()| rank_ab | plot_spacer()| venn_plots) +  
  plot_annotation(tag_levels = 'A')+
  plot_layout(guides = "collect",widths = c(1,0.25, 2,0.25,2))&
  theme(legend.position = "bottom")


ggsave(filename = "Figure2.svg",device = "svg",dpi = 300,width = 16 ,height = 6,units = "in")
ggsave(filename = "Figure2.eps",device = cairo_ps,dpi = 300,width = 16 ,height = 6,units = "in")
#---------------------------------------------##
# Get agreement statistics between pipelines----
#---------------------------------------------##

{
  #---------------------------------------------##
  ## Detection concordance using jaccard index ----
  #---------------------------------------------##

# make functions
jaccard <- function(x,y){
  length(intersect(names(x), names(y))) /
  length(union(names(x), names(y)))
}

jaccard_topk <- function(x, y, k = 20) {
  sx <- names(sort(x,decreasing = T))[1:min(k, length(x))]
  sy <- names(sort(y,decreasing = T))[1:min(k, length(y))]
  
  length(intersect(sx, sy)) /
    length(union(sx, sy))
}



# make a list of named abundances
list_grna <- split(all_gRNAs_pipelines,f = all_gRNAs_pipelines$gRNA)

list_grna <- lapply(list_grna, function(x){

    split(setNames(x$counts,x$cluster),x$workflow)
})


## Jaccard distance matrix
#---------------------------------------------##
jaccard_all <- lapply(list_grna, function(x){
  n <- length(x)
  
  jmat <- matrix(
    NA,
    nrow = n,
    ncol = n,
    dimnames = list(names(x), names(x))
  )
  
  for(i in seq_len(n)){
    for(j in seq_len(n)){
      jmat[i, j] <- jaccard(x[[i]], x[[j]])
    }
  }
  jmat
})


# Create pheatmaps and capture the grobs
ph_list_jaccard_all <- lapply(names(jaccard_all), function(n) {
  mat <- jaccard_all[[n]]
  
  pheatmap(cluster_rows = F, cluster_cols = F,
           mat,
           color = viridisLite::cividis(100),
           breaks = seq(0, 1, length.out = 101),legend = n=="TRAC",number_format = "%.3f",
           main = n,display_numbers = T, number_color = "white",
           silent = TRUE, na_col = "black",border_color = "white"
  )$gtable
})

# Arrange in a 2x2 matrix
jaccard_plot <- grid.arrange(
  grobs = ph_list_jaccard_all,
  ncol = 2,
  bottom = grid::textGrob(
    "Figure S1: Jaccard coefficient matrix",
    gp = grid::gpar(fontsize = 12, fontface = "plain"),
    x = unit(0.02, "npc"), 
    just = "left"
  )
)





## Jaccard top n distance matrix
#---------------------------------------------##

top_sites = 25

jaccard_top <- lapply(list_grna, function(a){
  n <- length(a)
  
  jmat <- matrix(
    NA,
    nrow = n,
    ncol = n,
    dimnames = list(names(a), names(a))
  )
  
  for(i in seq_len(n)){
    for(j in seq_len(n)){
      jmat[i, j] <- jaccard_topk(a[[i]], a[[j]],k=top_sites)
    }
  }
  jmat
})




# Create pheatmaps and capture the grobs
ph_list_jaccard_top <- lapply(names(jaccard_top), function(n) {
  mat <- jaccard_top[[n]]
  
  pheatmap(cluster_rows = F, cluster_cols = F,
           mat,
           color = viridisLite::cividis(100),
           breaks = seq(0, 1, length.out = 101),number_format = "%.3f",
           main = n,display_numbers = T, legend = n=="TRAC", number_color = "black",
           silent = TRUE, na_col = "black",border_color = "white"
  )$gtable
})

# Arrange in a 2x2 matrix
jaccard_top_plot <- grid.arrange(
  grobs = ph_list_jaccard_top,
  ncol = 2,top = grid::textGrob(
    paste("Jaccard distance of top",top_sites,"sites"),
    gp = gpar(fontsize = 16, fontface = "bold")
  )
)



  #---------------------------------------------##
  ## correlation on common sites ----
  #---------------------------------------------##

corr_shared <- function(x,y,type = "spearman"){
  common <- intersect(names(x), names(y))
  
  cor(x[common],
      
      y[common], use = "complete.obs",
      
      method=type)
}


correlation_shared <- lapply(list_grna, function(a){
  n <- length(a)
  
  jmat <- matrix(
    NA,
    nrow = n,
    ncol = n,
    dimnames = list(names(a), names(a))
  )
  
  for(i in seq_len(n)){
    for(j in seq_len(n)){
      jmat[i, j] <- corr_shared(a[[i]], a[[j]],"kendall")
    }
  }
  jmat
})

# Create pheatmaps and capture the grobs
ph_list_correlation <- lapply(names(correlation_shared), function(n) {
  mat <- correlation_shared[[n]]
  
  pheatmap(cluster_rows = F, cluster_cols = F,
           mat,
           color = viridisLite::cividis(100),
           breaks = seq(0, 1, length.out = 101),number_format = "%.3f",legend = n=="TRAC",
           main = n,display_numbers = T, number_color = "black",
           silent = TRUE, na_col = "black",border_color = "white"
  )$gtable
})

# Arrange in a 2x2 matrix
corr_plot <- grid.arrange(
  grobs = ph_list_correlation,
  ncol = 2,
  bottom = grid::textGrob(
    "Figure S2: Kendall correlation matrix (shared sites only).\nWhen a single site is shared, no correlation can be calculated (black cells)",
    gp = grid::gpar(fontsize = 12, fontface = "plain"),
    x = unit(0.02, "npc"), 
    just = "left"
  )
)




  #---------------------------------------------##
  ## RBO metric ----

p = 0.8
top <- Inf # number of top ranked site to evaluate, Inf for all, an integer otherwise

rbo_scores <- lapply(list_grna, function(a,k=top){
  n <- length(a)
  
  jmat <- matrix(
    NA,
    nrow = n,
    ncol = n,
    dimnames = list(names(a), names(a))
  )
  
  for(i in seq_len(n)){
    for(j in seq_len(n)){
      ni <- length(a[[i]])
      nj <- length(a[[j]])
      # If k is set, take the top k sites per list, otherwise take all elements of each list
      if(k!=Inf){
        nij <- min(ni,nj,k)
        ni <- nj <- nij
      }
      
      jmat[i, j] <- rbo_ext_fast2(L = names(sort(a[[i]],decreasing = T))[1:ni],S =  names(sort(a[[j]], decreasing = T))[1:nj],p =  p)$rbo
    }
  }
  jmat
})


# Create pheatmaps and capture the grobs
ph_list_rbo <- lapply(names(rbo_scores), function(n) {
  mat <- rbo_scores[[n]]
  
  pheatmap(cluster_rows = F, cluster_cols = F,
           mat,
           color = viridisLite::cividis(100),number_format = "%.3f",
           breaks = seq(0, 1, length.out = 101),legend = n=="TRAC",
           main = n,display_numbers = T, number_color = "black", 
           silent = TRUE, na_col = "black",border_color = "white"
  )$gtable
})


# Arrange in a 2x2 matrix
rbo_plot <- grid.arrange(
  grobs = ph_list_rbo,
  ncol = 2,
  bottom = grid::textGrob(
    "Figure S3: Rank Biased Overlap (RBO) matrix (Persistence parameter p=0.8)",
    gp = grid::gpar(fontsize = 12, fontface = "plain"),
    x = unit(0.02, "npc"), 
    just = "left"
  )
)


#---------------------------------------------##
## Calculate Kendall W concordance score on ranks ----
#---------------------------------------------##

library(irr)

for(grna in c("B2M","TRAC","VEGFAs2","VEGFAs3")){
  x_df <- all_gRNAs_pipelines_wide %>% ungroup %>% 
    filter(gRNA == grna) %>% 
    select(cluster,ends_with("prop")) %>%
    column_to_rownames("cluster") %>%
    mutate(across(everything(), ~replace_na(.x, 0))) %>% 
    mutate(across(everything(), ~dense_rank(-.x))) %>%         # calculate the rank
    filter(if_any(everything(), ~ . <= 20))                    # keep site that represent more than x % of total abundance
  if(nrow(x_df)>1){
    x = kendall(x_df,correct = T)
  }else {
    x = 1
  }
  cat("####################\n")
  print(x)
  cat("####################\n")
  }
}
#---------------------------------------------##
# generate a pdf supplementary file ----
#---------------------------------------------##


pdf("supplemental_figures1-3.pdf",paper = "a4")
plot(jaccard_plot)
plot(corr_plot)
plot(rbo_plot)
dev.off()
