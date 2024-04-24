## Load package
# library(tableone) 
# library(broom)
library(dplyr)

# function to get the length in Python style
len <- length

# simplified functions
fdr <- function(x) {p.adjust(x, "fdr")}
spearman.cor <- function(x, y) cor.test(x, y, method = "spearman")
pearson.cor <- function(x, y) cor.test(x, y, method = "pearson")

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv <- uniqv[uniqv!=""]
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
# https://www.tutorialspoint.com/r/r_mean_median_mode.htm#:~:text=The%20mode%20is%20the%20value,a%20data%20set%20in%20R.


permanova <- function(profile, item, pheno_table, perm = 999) {
  set.seed(1)
  perm <- adonis2(profile ~ pheno_table[, item], 
                 na.action=na.exclude,
                 permutations = perm) # changed from 9999 to save time
  return(perm)
}

adonis.pair<-function(dist.mat,Factor,nper=1000,corr.method="fdr"){
  require(vegan)
  as.factor(Factor)
  comb.fact<-combn(levels(Factor),2)
  pv<-NULL
  R2<-NULL
  SS<-NULL
  MeanSqs<-NULL
  F.Model<-NULL
  for (i in 1:dim(comb.fact)[2]){
    model.temp<-adonis(as.dist(as.matrix(dist.mat)[Factor==comb.fact[1,i] | Factor==comb.fact[2,i],Factor==comb.fact[1,i] | Factor==comb.fact[2,i]])~Factor[Factor==comb.fact[1,i] | Factor==comb.fact[2,i]],permutations=nper)
    pv<-c(pv,model.temp$aov.tab[[6]][1])
    R2<-c(R2,model.temp$aov.tab$R2[1])
    SS<-c(SS,model.temp$aov.tab[[2]][1])
    MeanSqs<-c(MeanSqs,model.temp$aov.tab[[3]][1])
    F.Model<-c(F.Model,model.temp$aov.tab[[4]][1])
  }
  pv.corr<-p.adjust(pv,method=corr.method)
  data.frame(combination=paste(comb.fact[1,],comb.fact[2,],sep=" <-> "),SumsOfSqs=SS,MeanSqs=MeanSqs,F.Model=F.Model,R2=R2,P.value=pv,P.value.corrected=pv.corr)}
# https://rdrr.io/github/GuillemSalazar/EcolUtils/src/R/EcolUtils_functions.R

hash_lineage <- function(x_, hash_) { # hash_ <- hash()
  a <- strsplit(x_, ";?\\S__")[[1]][-1]
  for (i in 2:length(a)) {
    hash_[a[i]] <- a[1:i]
  }
  # return(hash_)
}



inbetween <- function(x, interval) {
  res <- rep(F, length(x))
  for (i in 1:(length(interval)/2)) {
    loc = ifelse(i==1, 1, i+1)
    res <- res | (x >= interval[loc] & x <= interval[loc+1])
  }
  return(res)
}

## median functions
median_iqr<-function(x, digits=2){
  x=na.omit(x)
  median<-round(median(x),digits)
  quantile<-round(quantile(x,probs=c(0.25,0.75)),digits)
  a<-paste0(median," (", quantile[1],", ", quantile[2],")")
  return(a)
}

median_iqr_group<-function(x, group, prefix=""){
  dat <- na.omit(data.frame(x=x, group=group))
  x <- dat$x; group=dat$group
  a <- aggregate(x, by=list(group), median_iqr)
  b <- data.frame(t(a[,2]))
  colnames(b) <- paste0(prefix, a[,1])
  return(b)
}

median_iqr_group_wilcoxon <- function(x, group, prefix="", paired=F, exact=T){
  dat <- na.omit(data.frame(x=x, group=group))
  x <- dat$x; group=dat$group
  stats_ <- median_iqr_group(x, group, prefix=prefix)
  p <- format_p(wilcox.test(x~group, paired=paired, exact=exact)$p.value)
  res <- cbind(stats_, wilcox.p=p)
  return(res)
}

median_iqr_group_wilcoxon.raw_p <- function(x, group, prefix="", exact=T, paired=F){
  dat <- na.omit(data.frame(x=x, group=group))
  x <- dat$x; group=dat$group
  stats_ <- median_iqr_group(x, group, prefix=prefix)
  p <- wilcox.test(x~group, paired=paired, exact=exact)$p.value
  res <- cbind(stats_, wilcox.p=p)
  return(res)
}

median_iqr_group_wilcoxon.raw_p.fdr.matrix <- function(mat, groups) {
  # mat: samples in rows, features in columns
  results <- c()
  for (i in colnames(mat)) {
    x <- as.vector(mat[,i])
    result <- median_iqr_group_wilcoxon.raw_p(x, groups)
    results <- rbind(results, cbind(feature=i, result))
  }
  results$wilcox.fdr <- p.adjust(results$wilcox.p, "fdr")
  return(results)
}

median_iqr_group_dunn <- function(x, group, padjust="none", prefix=""){
  dat <- na.omit(data.frame(x=x, group=group))
  x <- dat$x; group=dat$group
  stats_ <- median_iqr_group(x, group, prefix=prefix)
  res <- dunn.test::dunn.test(x, group, method = padjust)
  p <- sapply(res$P.adjusted, format_p)
  p.df <- data.frame(t(p)); colnames(p.df) <- res$comparisons
  res <- cbind(stats_, p.df)
  return(res)
}

paired_wilcox_test <- function(df) {
  # df: subject, time, var
  rownames(df) <- paste(df[, 1], df[, 2], sep="-")
  time.points <- unique(df[, 2])
  all_pairs <- utils::combn(time.points, 2)
  
  result <- median_iqr_group(df[,3], df[,2])
  
  temp.header <- c()
  temp.result <- c()
  # time point pairs
  for (i in 1:ncol(all_pairs)) {
    
    pair <- all_pairs[, i]
    t1 <- pair[1]
    t2 <- pair[2]
    
    # subjects with both time points
    df.temp <- df[df[, 2] %in% pair, ]
    count.tb <- table(df.temp[, 1])
    subject.both <- names(count.tb)[count.tb==2]
    var.t1 <- df.temp[paste(subject.both, t1, sep="-"), 3]
    var.t2 <- df.temp[paste(subject.both, t2, sep="-"), 3]
    
    temp.result <- c(temp.result, wilcox.test(var.t1, var.t2, paired = T)$p.value)
    temp.header <- c(temp.header, paste(pair, collapse="-"))
  }
  result[, paste(temp.header, "p")] = sapply(temp.result, format_p)
  result[, paste(temp.header, "FDR")] = sapply(p.adjust(temp.result, "fdr"), format_p)
  
  result[, "kruskal"] = format_p(kruskal.test(df[,3]~df[,2])$p.value)
  
  return(result)
}

## mean functions
mean_sd_group_t_test_p <- function(x, group, prefix="", digits=2){
  dat <- na.omit(data.frame(x=x, group=group))
  x <- dat$x; group=dat$group
  
  stats_ <- mean_sd_group(x, group, prefix)
  p <- t.test(x~group)$p.value
  res <- cbind(stats_, t.test.p=p)
  return(res)
}

mean_sd <- function(x){
  x = na.omit(x)
  mean_ <- round(mean(x), 2)
  sd_ <- round(sd(x), 2)
  res <- paste0(mean_, " (", sd_, ")")
  return(res)
}

mean_sd_group <- function(x, group, prefix="") {
  dat <- na.omit(data.frame(x=x, group=group))
  x <- dat$x; group=dat$group
  
  a <- aggregate(x, by=list(group), mean_sd)
  b <- data.frame(t(a[,2]))
  colnames(b) <- paste0(prefix, a[,1])
  return(b)
}

## anova
anova_func <- function(my_data, x, group){
  res.aov <- aov(my_data[, x] ~ my_data[, group])
  p <- tidy(res.aov)$p.value[1]
  # return(round(p, 3))
  return(p)
}


# format proportions; fisher test
n.prop <- function(x, group=1, x.lab="", group.lab="") {
  if(length(unique(group))==1) {
    res <- data.frame(name=names(table(x)), stats=paste0(table(x), " (", round(prop.table(table(x))*100, 1), ")"))
    return(res)
  }
  df_ <- na.omit(data.frame(cbind(x=as.character(x), group=as.character(group))))
  tb <- table(df_)
  res <- data.frame(matrix(paste0(data.frame(tb)[, 3], " (", round(data.frame(prop.table(tb, 2))[, 3]*100, 1), ")"), 
                           nrow=length(na.omit(unique(x)))))
  rownames(res) <- paste0(x.lab, rownames(as.matrix(tb)))
  colnames(res) <- paste0(group.lab, colnames(as.matrix(tb)))
  return(res)
}

n.prop.fisher.raw_p <- function(x, group=1, x.lab="", group.lab="", to_return="all", simulate.p.value = F) {
  if(length(unique(group))==1) {
    res <- data.frame(name=names(table(x)), stats=paste0(table(x), " (", round(prop.table(table(x))*100, 1), ")"))
    return(res)
  }
  df_ <- na.omit(data.frame(cbind(x=as.character(x), group=as.character(group))))
  tb <- table(df_)
  res <- data.frame(matrix(paste0(data.frame(tb)[, 3], " (", round(data.frame(prop.table(tb, 2))[, 3]*100, 1), ")"), 
                           nrow=length(na.omit(unique(x)))))
  rownames(res) <- paste0(x.lab, rownames(as.matrix(tb)))
  colnames(res) <- paste0(group.lab, colnames(as.matrix(tb)))
  
  p <- c()
  if (nrow(tb)==1) {
    p <- 1
  } else {
    p <- fisher.test(tb, simulate.p.value=simulate.p.value)$p.value
  }
  
  res$fisher.p <- p
  
  if (to_return!="all") {
    res <- res[to_return, ]
  }
  
  return(res)
}

# https://rdrr.io/cran/RVAideMemoire/src/R/fisher.multcomp.R
fisher.multcomp <-
function(tab.cont,p.method="fdr") {
  if (is.matrix(tab.cont)) {tab.cont <- as.table(tab.cont)}
  call <- match.call()
  dname <- if(length(call$tab.cont)==1) {
    call$tab.cont
  } else {
    paste(call$tab.cont[1],"(",paste(call$tab.cont[-1],collapse=","),")",sep="")
  }
  if (!is.table(tab.cont)) {stop("'",dname,"' is not a \"table\" object",sep="")}
  if (is.null(colnames(tab.cont)) | any(0%in%nchar(colnames(tab.cont)))) {
    colnames(tab.cont) <- paste0("col",1:ncol(tab.cont))
  }
  if (is.null(rownames(tab.cont)) | any(0%in%nchar(rownames(tab.cont)))) {
    rownames(tab.cont) <- paste0("row",1:nrow(tab.cont))
  }
  colonnes <- combn(colnames(tab.cont),2)
  lignes <- combn(rownames(tab.cont),2)
  colonnes2 <- apply(colonnes,2,function(x) paste(x,collapse=":"))
  lignes2 <- apply(lignes,2,function(x) paste(x,collapse=":"))
  if (ncol(tab.cont)>2) {
    p.no.adjust <- matrix(0,nrow=ncol(lignes),ncol=ncol(colonnes),dimnames=list(lignes2,colonnes2))
    for (i in 1:ncol(colonnes)) {
	for (j in 1:ncol(lignes)) {
	  tab <- tab.cont[lignes[,j],colonnes[,i]]
	  p.no.adjust[j,i] <- fisher.test(tab)$p.value
	}
    }
    p.adjust <- matrix(p.adjust(p.no.adjust,method=p.method),nrow=nrow(p.no.adjust),ncol=ncol(p.no.adjust),
	dimnames=dimnames(p.no.adjust))
  } else {
    p.no.adjust <- matrix(NA,nrow=nrow(tab.cont),ncol=nrow(tab.cont),dimnames=list(rownames(tab.cont),rownames(tab.cont)))
    for (i in 1:ncol(lignes)) {
	tab.temp <- tab.cont[lignes[,i],]
	p.no.adjust[lignes[2,i],lignes[1,i]] <- fisher.test(tab.temp)$p.value
    }
    p.adjust <- matrix(p.adjust(p.no.adjust,method=p.method),nrow=nrow(p.no.adjust),ncol=ncol(p.no.adjust),
	dimnames = dimnames(p.no.adjust))
    p.adjust <- p.adjust[-1,-ncol(p.adjust)]
  }
  call <- match.call()
  dname <- if(length(call$tab.cont)==1) {call$tab.cont} else {paste(call$tab.cont[1],"(",paste(call$tab.cont[-1],collapse=","),")",sep="")}
  result<-list(method="Fisher's exact test for count data",data.name=dname,p.adjust.method=p.method,p.value=p.adjust)
  class(result) <- "RV.multcomp"
  return(result)
}
# get beta (95%CI) from regression models
glm_OR <- function(mod, digits=2) {
  smr <- summary(mod)
  res <- cbind(round(cbind(OR=exp(coefficients(mod)), exp(confint.default(mod))), digits=digits), 
        p=sapply(summary(mod)$coefficients[, ncol(smr$coefficients)], format_p))
  return(res)
}

glm_beta <- function(mod, digits=2) {
smr<-summary(mod)
res <- cbind(round(cbind(coefficients(mod), confint.default(mod)), digits = digits), 
      p=sapply(smr$coefficients[,ncol(smr$coefficients)], format_p))

return(res)
}

effect_size_p <- function(x) {
  coef_ <- coef(summary(x))
  p_ <- sapply(coef_[,4], format_p)
  
  est_ <- round(coef_[,1], 2)
  confint_ <- round(confint(x), 2)
  
  res <- c()
  title <- c()
  for (i in 1:length(p_)) {
    res <- c(res, paste0(est_[i], " (", confint_[i, 1], ", ", confint_[i, 2], ")"), p_[i])
    title <- c(title, names(p_)[i], "p")
  }
  res <- 
    return(cbind(res, title))
}

format_effect_size <- function(glm_res, prefix="") {
  res <- cbind(paste0(glm_res[, 1], " (", glm_res[, 2], ", ", glm_res[, 3], ")"),
               p=glm_res[, 4])
  rownames(res) <- rownames(glm_res)
  colnames(res)[1] <- prefix
  return(res)
}

# fold change
fc <- function(x, y, method="mean", pseudo=1e-5) { # binary only
  cats <- unique(y)
  fc_ = 1
  if(method=="mean") {
    fc_ <- (mean(x[y == cats[1]], na.rm = T) + pseudo) / (mean(x[y == cats[2]], na.rm = T) + pseudo)
  }
  if(method=="median"){
    fc_ <- (median(x[y == cats[1]], na.rm = T) + pseudo) / (median(x[y == cats[2]], na.rm = T) + pseudo)
  }
  
  res <- c(fc_, log2(fc_))
  names(res) <- c(paste0(cats[1], "/", cats[2], " FC"),
                  "log2FC")
  return(res)
}

# Get VIP of variables in a PLSDA model

# https://rdrr.io/cran/RVAideMemoire/src/R/PLSDA.VIP.R
PLSDA.VIP <- function(model,graph=FALSE) {
  if (packageVersion("mixOmics")<"5.0.2") {
    stop(paste("you must update 'mixOmics' to version >= 5.0.2 (actual: ",
               packageVersion("mixOmics"),")",sep=""))
  }
  VIP <- mixOmics::vip(model)
  tab <- as.data.frame(VIP[order(VIP[,ncol(VIP)],decreasing=TRUE),ncol(VIP)])
  colnames(tab) <- "VIP"
  if (graph) {
    opar <- par()
    on.exit(suppressWarnings(par(opar)))
    par(mar=c(5,8,2,2),las=1)
    g <- barplot(rev(tab$VIP),horiz=TRUE,xlab=paste("VIP (",ncol(VIP),ifelse(ncol(VIP)>1," axes)"," axis)"),
                                                    sep=""))
    mtext(rev(rownames(tab)),side=2,line=1,at=g,cex=0.7)
    abline(h=g,lty=3,col="grey40")
    abline(v=1,lty=2,col="red")
  }
  result <- list(tab=tab,sup1=rownames(tab)[which(tab$VIP>1)])
  class(result) <- "PLSDA.VIP"
  return(result)
}

extract_mediation_summary <- function (x) { 
    clp <- 100 * x$conf.level
    isLinear.y <- ((class(x$model.y)[1] %in% c("lm", "rq")) || 
                       (inherits(x$model.y, "glm") && x$model.y$family$family == 
                            "gaussian" && x$model.y$family$link == "identity") || 
                       (inherits(x$model.y, "survreg") && x$model.y$dist == 
                            "gaussian"))
    
    printone <- !x$INT && isLinear.y
    
    if (printone) {
        
        smat <- c(x$d1, x$d1.ci, x$d1.p)
        smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
        smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
        smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
        
        rownames(smat) <- c("ACME", "ADE", "Total Effect", "Prop. Mediated")
        
    } else {
        smat <- c(x$d0, x$d0.ci, x$d0.p)
        smat <- rbind(smat, c(x$d1, x$d1.ci, x$d1.p))
        smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
        smat <- rbind(smat, c(x$z1, x$z1.ci, x$z1.p))
        smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
        smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
        smat <- rbind(smat, c(x$n1, x$n1.ci, x$n1.p))
        smat <- rbind(smat, c(x$d.avg, x$d.avg.ci, x$d.avg.p))
        smat <- rbind(smat, c(x$z.avg, x$z.avg.ci, x$z.avg.p))
        smat <- rbind(smat, c(x$n.avg, x$n.avg.ci, x$n.avg.p))
        
        rownames(smat) <- c("ACME (control)", "ACME (treated)", 
                            "ADE (control)", "ADE (treated)", "Total Effect", 
                            "Prop. Mediated (control)", "Prop. Mediated (treated)", 
                            "ACME (average)", "ADE (average)", "Prop. Mediated (average)")
        
    }
    
    colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep = ""), 
                        paste(clp, "% CI Upper", sep = ""), "p-value")
    smat
    
}

# format p values
format_p <- function(p, digits=3) {
  if (is.na(p)) {return(NA)}
  if (p >= 0.001) {
    return(round(p, digits))
  } else {
    return("<0.001")
  }
}

format_p_list <- function(x) {
	return(sapply(x, format_p))
}

# return the last n characters
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
# https://stackoverflow.com/questions/7963898/extracting-the-last-n-characters-from-a-string-in-r

# remove the last n characters
rmstrRight <- function(x, n){
  substr(x, 1, nchar(x)-n)
}
## text processing for metaphlan profiles
get_species_name <- function(x) {
  if (is.character(x)) {
    return(gsub("_", " ", gsub(".*s__", "", x)))
  }
  if (is.vector(x)) {
    return(sapply(x, function(y) {gsub("_", " ", gsub(".*s__", "", y))}))
  }
  return(NULL)
}

get_phylum_name <- function(x) {
  if (is.character(x)) {
    return(gsub("\\|.*", "", gsub(".*p__", "", x)))
  }
  if (is.vector(x)) {
    return(sapply(x, function(y) {gsub("\\|.*", "", gsub(".*p__", "", x))}))
  }
  return(NULL)
}

# matrix / data.frame
get_shared_count_prop <- function(sample_list, mat_1, mat_2, details=F) {
  richness <- total_ab <- c()
  if(details) {samples <- shared_species <- c()}
  for (f in sample_list) {
    g1 <- rownames(mat_1)[mat_1[, f]>0]
    g2 <- rownames(mat_2)[mat_2[, f]>0]
    idx <- g1 %in% g2
    richness <- c(richness, sum(idx))
    total_ab <- c(total_ab, sum(mat_1[idx, f]))
    
    if(details) {
      samples <- c(samples, 
                   rep(f, sum(idx)))
      shared_species <- c(shared_species,
                          g1[idx])
    }
  }
  result <- data.frame(richness=richness, abundance=total_ab)
  rownames(result) <- sample_list
  if (details) {
    detailed_res=data.frame(sample=samples, shared_species=shared_species)
    return(list(result=result, details=detailed_res))
  } else {
    return(result)
  }
}
		   
# color functions/code
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

# color-blind friendly colors
eight_colors= 
  c("#E69F00", # orange
  "#0072B2", # blue
  "#000000", # black
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#56B4E9", # sky blue
  "#D55E00", # vermilion
  "#CC79A7") # reddish purple

four_colors= 
  c("#E69F00", # orange
    "#F0E442", # yellow
  "#0072B2", # blue
  "#56B4E9" # sky blue
)
		   
two_colors=
  c("#E69F00", # orange
  "#0072B2") # blue

# twenty_colors <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', 
# '#f58231', '#911eb4', '#46f0f0', '#f032e6', 
# '#bcf60c', '#fabebe', '#008080', '#e6beff', 
# '#9a6324', '#fffac8', '#800000', '#aaffc3', 
# '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
# List of 20 Simple, Distinct Colors
# https://sashamaps.net/docs/resources/20-colors/


# library(rcartocolor)
# nColor <- 12
# scales::show_col(carto_pal(nColor, "Safe"))


# https://www.jianshu.com/p/71fc7e2561c4
# https://www.sci666.com.cn/60850.html


# I/O
write.tsv <- function(df, file, sep = "\t", col.names = NA, ...) {
  write.table(df, file, sep = sep, col.names = col.names, quote = F, ...)
}

read.tsv <- function(file, row.names=1, header=T, check.names = F, quote=c(""), as.is=T, ...) {
  read.table(file, sep = "\t", row.names = row.names, header = header, check.names = check.names, quote=quote, as.is=as.is, ...)
}

# pdf to png
pdf_to_png <- function(pdf_, dpi=250) {
  require(pdftools)
  png_ <- gsub("pdf$", "png", pdf_)
  pdf_convert(pdf_, filenames = png_, 
              dpi=dpi)
}

myCol3=c("#009E73", "#CC79A7", "#56B4E9")
myCol3_2=c("#0072B2", "#E69F00", "#009E73")
names(myCol3) = names(myCol3_2) = c("D0", "M1", "M6")

myCol3_3=c("#0072B2", "#E69F00", "#009E73")
names(myCol3_3) = c("Baseline", "1m p.v.", "6m p.v.")

m1_m6_color = c("#E69F00", "#009E73")
names(m1_m6_color) <- c("")
