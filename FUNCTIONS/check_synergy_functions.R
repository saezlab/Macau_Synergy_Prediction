
################################################################ AZ functions ################################################################

check_SYNERGY_AZ <- function(inverse = 1) {  # t,high,low,genomics=c(),exclude=c()
  
  ## see which drug combo act on those targets simultaneously
  target_1 <- rownames(SYN_MATRIX) [ grep(t[1], rownames(SYN_MATRIX)) ] 
  target_2 <- rownames(SYN_MATRIX) [ grep(t[2], rownames(SYN_MATRIX)) ] 
  selected_drug_combo <- intersect(target_1,target_2); selected_drug_combo
  
  common <- intersect(rownames(progeny) , cell_names )
  
  ## compute the metrics: Delta Pathway activity and SNP CNV
  metrics1 <- as.matrix(progeny[ common , high ]) ; metrics1 <- rowMeans(metrics1,na.rm=T)
  metrics2 <- as.matrix(progeny[ common ,low]   ) ; metrics2 <- rowMeans(metrics2,na.rm=T)
  Additional_Parameter <- Additional_Parameter[common, ]
  if(length(genomics)>0) { 
    metrics <-  metrics1 - metrics2 + rowSums(as.matrix(Additional_Parameter[ ,genomics]))
  } else {metrics <-  metrics1 - metrics2}
  names(metrics) <- common
  
  metrics <- metrics*inverse
  
  ## selecte drug combo to compute the Average Synergy: at least 1 or 2 observation 
  if(length(selected_drug_combo)>1) {  
    x <- SYN_MATRIX[selected_drug_combo, rownames(progeny) ] ; minimum_combo <- 1 
  } else if (length(selected_drug_combo)==1 ) {
    x <- SYN_MATRIX[selected_drug_combo, rownames(progeny) ] ; minimum_combo <- 0
    x <- t(matrix(x)); rownames(x)<-selected_drug_combo; colnames(x)<-rownames(progeny); x<-rbind(x,rep(NA,length(x)))
  }
  if(length(exclude)>0) {x <- x[-which(rownames(x) %in% exclude), ]}
  
  to_remove <- colnames(x)[which( colSums(!is.na(x)) <= minimum_combo )]  ## value to remove 
  x <- as.matrix( x[ ,colnames(x) %not in% to_remove])  
  selected_drug_combo <- rownames(x) ; selected_drug_combo <- target_combo[rownames(target_combo) %in% selected_drug_combo, ] ; selected_drug_combo <- selected_drug_combo[ ,-which(colSums(selected_drug_combo)==0)]
  print(selected_drug_combo)
  
  x_vector <- colMeans(x, na.rm=T) ## compute the Average Synergy
  
  a <- metrics[names(metrics)  %not in% to_remove] ; b <- x_vector[names(x_vector) %in% common]
  df <- data.frame(cbind(a,b)) ;  cor.test(df$a,df$b)
  return(list(df,selected_drug_combo))
}

check_SYNERGY_AZ_ALL_DRUG <- function(inverse = 1) {  # t,high,low,genomics=c(),exclude=c()
  
  ## see which drug combo act on those targets simultaneously
  target_1 <- rownames(SYN_MATRIX) [ grep(t[1], rownames(SYN_MATRIX)) ] 
  target_2 <- rownames(SYN_MATRIX) [ grep(t[2], rownames(SYN_MATRIX)) ] 
  selected_drug_combo <- intersect(target_1,target_2); selected_drug_combo
  
  common <- intersect(rownames(progeny) , cell_names )
  
  ## compute the metrics: Delta Pathway activity and SNP CNV
  metrics1 <- as.matrix(progeny[ common , high ]) ; metrics1 <- rowMeans(metrics1,na.rm=T)
  metrics2 <- as.matrix(progeny[ common ,low]   ) ; metrics2 <- rowMeans(metrics2,na.rm=T)
  Additional_Parameter <- Additional_Parameter[common, ]
  if(length(genomics)>0) { metrics <-  metrics1 - metrics2 + 1*Additional_Parameter[ ,genomics] } else {metrics <-  metrics1 - metrics2}
  names(metrics) <- common
  
  metrics <- metrics*inverse
  
  ## selecte drug combo to compute the Average Synergy: at least 1 or 2 observation 
  if(length(selected_drug_combo)>1) {  
    x <- SYN_MATRIX[selected_drug_combo, rownames(progeny) ] ; minimum_combo <- 1 
  } else if (length(selected_drug_combo)==1 ) {
    x <- SYN_MATRIX[selected_drug_combo, rownames(progeny) ] ; minimum_combo <- 0
    x <- t(matrix(x)); rownames(x)<-selected_drug_combo; colnames(x)<-rownames(progeny); x<-rbind(x,rep(NA,length(x)))
  }
  # x <- t(as.matrix(x[1, ])) ; rownames(x)[1] <- selected_drug_combo[1] ; x<-rbind(x,rep(NA,length(colnames(x)))); minimum_combo <- 0 
  if(length(exclude)>0) {x <- x[-which(rownames(x) %in% exclude), ]}
  
  to_remove <- colnames(x)[which( colSums(!is.na(x)) <= minimum_combo )]  ## value to remove 
  x <- as.matrix( x[ ,colnames(x) %not in% to_remove])  
  selected_drug_combo <- rownames(x) ; selected_drug_combo <- target_combo[rownames(target_combo) %in% selected_drug_combo, ] ; selected_drug_combo <- selected_drug_combo[ ,-which(colSums(selected_drug_combo)==0)]
  print(selected_drug_combo)
  
  x_list <- list()
  for(i in 1:length(rownames(x))) {
    x_temp <- x[i,complete.cases(x[i, ])]
   # if(length(x_temp) >= 10) { x_list[[i]] <- x[i, ] }
    x_list[[i]] <- x_temp
  }
  names(x_list) <- rownames(x)
  x_list <- x_list[lapply(x_list,length)>0]
  
  
  df_list <- list()
  for(i in 1:length(x_list)) {
    a <- metrics[names(metrics)  %not in% to_remove] ; b <- x_list[[i]][names(x_list[[i]]) %in% common]
    common <- intersect(names(a),names(b))
    a <- a[common] ; b <- b[common]
    df <- data.frame(cbind(a,b))
    df_list[[i]] <- df
  }
  names(df_list) <- names(x_list)
  return(list(df_list,selected_drug_combo))
}


scatter_plot_AZ <- function(title,df,text_size=33,title_size=3.5,dot_size=12,coeff_x=0.5,coeff_y=0) {  ## ~~italic("p")~"="~p
  # equation, correlation and p value
  out <- cor.test(df$a,df$b) ; r <- out$estimate ; p <- out$p.value
  lm_eqn <- function(df){
    m <- lm(b ~ a, df);
    eq <- substitute(~~italic("r")~"="~r*"",
                     list(a = format(coef(m)[1], digits = 2), 
                          b = format(coef(m)[2], digits = 2), 
                          r = format(r, digits = 2),
                          p = format(p, digits=2)))
    as.character(as.expression(eq));                 
  }
  
  g <- ggplot(df, aes(a, b, color = b)) + 
    geom_point(shape = 16, size = dot_size, show.legend = FALSE, alpha = .8 ) +  geom_smooth(method=lm,se=F,show.legend=F) + 
    labs(x = "Delta Pathway Activity", y="Average Synergy" ) + ggtitle(title) + 
    theme(legend.position="bottom",axis.text=element_text(size= text_size) , axis.title= element_text(size=text_size), plot.title = element_text(size =rel(title_size), hjust = 0.5 ),
          panel.background = element_rect(fill='white'),panel.grid.major = element_line(colour = "grey90") ) + 
    geom_text(x = min(df$a) + coeff_x*(max(df$a)-min(df$a)), y = min(df$b) + coeff_y, label = lm_eqn(df), parse = TRUE,show.legend=F,color="black",size = 20) +  ##  0.02*(max(df$b)-min(df$b))
    scale_color_gradient(low = "#8B3A3A" , high = "#8B3A3A" ) #  
  g
}


############################################################ NOVARTIS functions ###############################################################

check_SYNERGY_NOVARTIS <- function(inverse = 1) {  # t,high,low, genomics=c() ; exclude=c()
  
  selected_drug_combo <- target_combo[ ,t] ; selected_drug_combo <- selected_drug_combo[which(rowSums(selected_drug_combo)==2), ]
  selected_drug_combo <- rownames(selected_drug_combo)
  
  common <- colnames(SYN_MATRIX_tissue)
  
  ## compute the metrics: Delta Pathway activity and SNP CNV
  metrics1 <- as.matrix(progeny[ common , high ]) ; metrics1 <- rowMeans(metrics1,na.rm=T)
  metrics2 <- as.matrix(progeny[ common ,low]   ) ; metrics2 <- rowMeans(metrics2,na.rm=T)
  if(length(genomics)>0) { metrics <-  metrics1 - metrics2 } else {metrics <-  metrics1 - metrics2}
  names(metrics) <- common
  
  metrics <- metrics*inverse
  
  ## selecte drug combo to compute the Average Synergy: at least 1 or 2 observation 
  if(length(selected_drug_combo)>1) {  
    x <- SYN_MATRIX[selected_drug_combo, rownames(progeny) ] ; minimum_combo <- 1 
  } else if (length(selected_drug_combo)==1 ) {
    x <- SYN_MATRIX[selected_drug_combo, rownames(progeny) ] ; minimum_combo <- 0
    x <- t(matrix(x)); rownames(x)<-selected_drug_combo; colnames(x)<-rownames(progeny); x<-rbind(x,rep(NA,length(x)))
  }
  # x <- t(as.matrix(x[1, ])) ; rownames(x)[1] <- selected_drug_combo[1] ; x<-rbind(x,rep(NA,length(colnames(x)))); minimum_combo <- 0 
  if(length(exclude)>0) {x <- x[-which(rownames(x) %in% exclude), ]}
  
  to_remove <- colnames(x)[which( colSums(!is.na(x)) <= minimum_combo )]  ## value to remove 
  x <- as.matrix( x[ ,colnames(x) %not in% to_remove])  
  selected_drug_combo <- rownames(x) ; selected_drug_combo <- target_combo[rownames(target_combo) %in% selected_drug_combo, ] ; selected_drug_combo <- selected_drug_combo[ ,-which(colSums(selected_drug_combo)==0)]
  print(selected_drug_combo)
  
  x_vector <- colMeans(x, na.rm=T) ## compute the Average Synergy
  
  a <- metrics[names(metrics)  %not in% to_remove] ; b <- x_vector[names(x_vector) %in% common]
  df <- data.frame(cbind(a,b)) ;  cor.test(df$a,df$b)
  return(list(df,selected_drug_combo))
}

scatter_plot_NOVARTIS <- function(title,df,text_size=33,title_size=3,dot_size=5) {
  # equation, correlation and p value
  out <- cor.test(df$a,df$b) ; r <- out$estimate ; p <- out$p.value
  lm_eqn <- function(df){
    m <- lm(b ~ a, df);
    eq <- substitute(~~italic("r")~"="~r*","~~italic("p")~"="~p,
                     list(a = format(coef(m)[1], digits = 2), 
                          b = format(coef(m)[2], digits = 2), 
                          r = format(r, digits = 2),
                          p = format(p, digits=2)))
    as.character(as.expression(eq));                 
  }
  
  g <- ggplot(df, aes(a, b, color = b)) + 
    geom_point(shape = 16, size = dot_size, show.legend = FALSE, alpha = .8 ) +  geom_smooth(method=lm,se=F,show.legend=F) + 
    labs(x = "Delta Pathway Activity", y="Average Synergy" ) + ggtitle(title) + 
    theme(legend.position="bottom",axis.text=element_text(size= text_size) , axis.title= element_text(size=text_size), plot.title = element_text(size =rel(title_size), hjust = 0.5 ),
          panel.background = element_rect(fill='white'),panel.grid.major = element_line(colour = "grey90") ) + 
    geom_text(x = min(df$a) + 0.5*(max(df$a)-min(df$a)), y = min(df$b) + 0.02*(max(df$b)-min(df$b)), label = lm_eqn(df), parse = TRUE,show.legend=F,color="black",size = 20) + 
    scale_color_gradient(low = "#0091ff", high = "#f0650e" ) #  
  g
}

############################################################# SANGER functions ################################################################

# synergy_score="XMID" ; inverse = 1 ; switch_anchor=F
check_SYNERGY_SANGER <- function(synergy_data,synergy_score,inverse = 1, switch_anchor=F) {
  ## see which drug combo act on those targets simultaneously
  drug_set <- target[ , union( grep(t[1], colnames(target)) , grep(t[2], colnames(target)) ) ]
  # drug_set <- drug_set[ ,- which(colSums(drug_set) == 0) ] 
  drug_set <- drug_set[ -which(rowSums(drug_set)==0) , ]
  drug_set <- drug_set[ rownames(drug_set) %in% c("Dabrafenib","BMS-754807"), ]
  drug_combo_to_test <- t(combn(x=rownames(drug_set), m=2, FUN = NULL, simplify = TRUE ))
  
  x <- synergy_data
  
  if(switch_anchor==T) { 
    drug_combo_to_test <- paste(drug_combo_to_test[1,2],drug_combo_to_test[1,1],sep="/")
  } else {
    drug_combo_to_test <- paste(drug_combo_to_test[drug,1],drug_combo_to_test[drug,2],sep="/")
  }
  
  ## SELECT CONCENTRATION
  if(conc=="max") {x <- x[ which(x$ANCHOR_CONC == max(x$ANCHOR_CONC)) , ]} else { x <- x[ which(x$ANCHOR_CONC == min(x$ANCHOR_CONC)) , ] }
  x <- x[ ,c("CELL_LINE_NAME",paste0("SYNERGY_DELTA_",synergy_score)) ] ; x_aggr <- aggregate(x[ ,2], list(x$CELL_LINE_NAME), mean ) ; b <- x_aggr$x ; names(b) <- x_aggr$Group.1
  
  progeny <- progeny11 [ unique(x$CELL_LINE_NAME), ] ; progeny <- progeny[complete.cases(progeny) , ]
  common <- rownames(progeny)
  
  ## compute the metrics: Delta Pathway activity and SNP CNV
  metrics1 <- as.matrix(progeny[ common , high ]) ; metrics1 <- rowMeans(metrics1,na.rm=T)
  metrics2 <- as.matrix(progeny[ common ,low]   ) ; metrics2 <- rowMeans(metrics2,na.rm=T)
  Additional_Parameter <- Additional_Parameter[common, ]
  if(length(genomics>0)) { metrics <-  metrics1 - metrics2 + Additional_Parameter[ ,genomics] } else {metrics <-  metrics1 - metrics2  }
  names(metrics) <- common
  
  #   progeny <- as.matrix(progeny)
  #   progeny_high_low <- progeny[ ,c(high,low)]
  #   metrics <- progeny_high_low * pathway[c(high,low)]
  # #  metrics <- progeny_high_low %*% pathway[c(high,low)]
  #   metrics <- rowMeans(metrics[ ,high]) - metrics[ ,low]
  
  metrics <- metrics*inverse
  
  a <- metrics ; common_cell <- intersect(names(a), names(b))
  a <- a[common_cell] ; b <- b[common_cell]
  df <- data.frame(cbind(a,b)) ; df <-df[complete.cases(df), ]  ; df <- df[order(-df$a), ]
  return(list(df,drug_combo_to_test))
}


scatter_plot_SANGER <- function(title,df,synergy_score,inverse=1,switch_anchor=F,text_size=33,title_size=3.5,dot_size=12, coeff_x=0.5,coeff_y=0) {
  # equation, correlation and p value
  out <- cor.test(df$a,df$b) ; r <- out$estimate ; p <- out$p.value
  lm_eqn <- function(df){
    m <- lm(b ~ a, df);
    eq <- substitute(~~italic("r")~"="~r*"",
                     list(a = format(coef(m)[1], digits = 2), 
                          b = format(coef(m)[2], digits = 2), 
                          r = format(r, digits = 2),
                          p = format(p, digits=2)))
    as.character(as.expression(eq));                 
  }
  
  g <- ggplot(df, aes(a, b, color = b)) + 
    geom_point(shape = 16, size = dot_size, show.legend = FALSE, alpha = .8 ) +  geom_smooth(method=lm,se=F,show.legend=F) + 
    labs(x = "Delta Pathway Activity", y=paste0(check_SYNERGY_SANGER(combi_reduce,synergy_score,switch_anchor)[[2]]," Delta ",synergy_score)) + ggtitle(title) + 
    theme(legend.position="bottom",axis.text=element_text(size= text_size) , axis.title= element_text(size=text_size), plot.title = element_text(size=rel(title_size), hjust=0.5),
          panel.background = element_rect(fill='white'),panel.grid.major = element_line(colour = "grey90") ) + 
    geom_text(x = min(df$a) + coeff_x*(max(df$a)-min(df$a)), y = min(df$b) + coeff_y*(max(df$b)-min(df$b)), label = lm_eqn(df), parse = TRUE,show.legend=F,color="black",size = 21) + 
    scale_color_gradient(low = "#8B3A3A" , high = "#8B3A3A" ) #  
  g
}


########################################################## SANGER - MIKE functions #############################################################

check_SYNERGY_SANGER_MIKE <- function(synergy_score="bliss", switch_anchor=F) {
  ## see which drug combo act on those targets simultaneously
  drug_set <- target[ , union( grep(t[1], colnames(target)) , grep(t[2], colnames(target)) ) ]
  # drug_set <- drug_set[ ,- which(colSums(drug_set) == 0) ] 
  drug_set <- drug_set[ -which(rowSums(drug_set)==0) , ]
  drug_set <- drug_set[ rownames(drug_set) %in% names(table(combi_reduce$ANCHOR_NAME)), ]
  drug_combo_to_test <- t(combn(x=rownames(drug_set), m=2, FUN = NULL, simplify = TRUE ))
  
  x <- combi_reduce
  
  if(switch_anchor==T) { 
    x <- x[ which(x$ANCHOR_NAME==drug_combo_to_test[drug,2]) , ]
    x <- x[ which(x$LIBRARY_NAME==drug_combo_to_test[drug,1]) , ]
    drug_combo_to_test <- paste(drug_combo_to_test[1,2],drug_combo_to_test[1,1],sep="/")
  } else {
    x <- x[ which(x$ANCHOR_NAME==drug_combo_to_test[drug,1]) , ]
    x <- x[ which(x$LIBRARY_NAME==drug_combo_to_test[drug,2]) , ]
    drug_combo_to_test <- paste(drug_combo_to_test[1,1],drug_combo_to_test[1,2],sep="/")
  }
  
  ## SELECT LIBRARY CONCENTRATION
  if(conc=="max") {x <- x[ which(x$ANCHOR_CONC == max(x$ANCHOR_CONC)) , ]} else { x <- x[ which(x$ANCHOR_CONC == min(x$ANCHOR_CONC)) , ] }
  
  LIB_CONC <- as.numeric(names(table(x$LIBRARY_CONC)))
  
  df_list <- list()
  for(i in 1:length(LIB_CONC)) {
    
    x_subset <- x[ which(x$LIBRARY_CONC == LIB_CONC[i] ) , ]  
    x_subset <- x_subset[ ,c("CELL_LINE_NAME",synergy_score) ] ; # x_subset[ ,synergy_score] <- as.numeric(x_subset[ ,synergy_score]) 
    x_aggr <- aggregate(x_subset[ ,2], list(x_subset$CELL_LINE_NAME), mean ) ; b <- x_aggr[ ,2] ; names(b) <- x_aggr$Group.1
    
    progeny <- progeny [ unique(x_subset$CELL_LINE_NAME), ] ; progeny <- progeny[complete.cases(progeny) , ]
    common <- rownames(progeny)
    
    ## compute the metrics: Delta Pathway activity and SNP CNV
    metrics1 <- as.matrix(progeny[ common , high ]) ; metrics1 <- rowMeans(metrics1,na.rm=T)
    metrics2 <- as.matrix(progeny[ common ,low]   ) ; metrics2 <- rowMeans(metrics2,na.rm=T)
    Additional_Parameter <- Additional_Parameter[common, ]
    if(length(genomics>0)) { metrics <-  metrics1 - metrics2 + Additional_Parameter[ ,genomics] } else {metrics <-  metrics1 - metrics2  }
    names(metrics) <- common
    
    a <- metrics ; common_cell <- intersect(names(a), names(b))
    b <- b[common_cell]
    df <- data.frame(cbind(a,b)) ; df <-df[complete.cases(df), ]  ; df <- df[order(-df$a), ]
    df_list[[i]] <- df
  }
  return(list(df_list,drug_combo_to_test,LIB_CONC))
}

scatter_plot_SANGER_MIKE <- function(title,df,synergy_score="bliss",switch_anchor=F,text_size=26,title_size=2.5) {
  out <- cor.test(df$a,df$b) ; r <- out$estimate ; p <- out$p.value
  lm_eqn <- function(df){
    m <- lm(b ~ a, df);
    eq <- substitute(italic("y") == a + b %.% italic(x)*","~~italic("r")~"="~r*","~~italic("p")~"="~p,
                     list(a = format(coef(m)[1], digits = 2), 
                          b = format(coef(m)[2], digits = 2), 
                          r = format(r, digits = 2),
                          p = format(p, digits=2)))
    as.character(as.expression(eq));                 
  }
  
  g <- ggplot(df, aes(a, b, color = b)) + 
    geom_point(shape = 16, size = 5, show.legend = FALSE, alpha = .8 ) +  geom_smooth(method=lm,se=F,show.legend=F) + 
    labs(x = "Delta Pathway Activity", y=paste0(check_SYNERGY_SANGER_MIKE()[[2]]," ",synergy_score)) + ggtitle(title) + 
    theme(legend.position="bottom",axis.text=element_text(size= text_size) , axis.title= element_text(size=text_size), plot.title = element_text(size=rel(title_size), hjust=0.5),
          panel.background = element_rect(fill='white'),panel.grid.major = element_line(colour = "grey90") ) + 
    geom_text(x = mean(df$a), y = mean(df$b), label = lm_eqn(df), parse = TRUE,show.legend=F,color="black",size = 10) + 
    scale_color_gradient(low = "#8B3A3A" , high = "#8B3A3A" ) #  
  g
}

