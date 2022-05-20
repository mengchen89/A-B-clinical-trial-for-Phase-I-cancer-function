threep3ABCDE.desc<-function (truep, A, B, C, D, E, dose = NULL)
{
  
  truep<-c(0.06,0.15)
  A=3; B=3; C=1; D=1; E=1
  dose = NULL
  if (!is.null(dose) & length(dose) != length(truep))
    stop("Length of 'dose' must be the same as the length of 'truep'.")#stop when length of dose not equal to length of truep except there's no dose
  path.mat<-prob<-ssize<-mtd<-dlt.no<-NULL
  exp<-0
  doses <- length(truep)
  mcohort <- 2 * doses
  mcplus1 <- mcohort + 1
  pmat <- as.data.frame(matrix(NA, nrow = 1, ncol = 3 * mcplus1 +
                                 1))#create a matrix
  colnames(pmat) <- c("stop", "desc", paste(c("d", "tox", "ssize"),
                                            rep(1:mcohort, each = 3)), paste("d", mcplus1))#giving names
                                                    #stop: whether the whole process stops ;desc: whether the whole process experiences descendent;
                                                    #d1: the dose level of first experiment; tox1: the number of dlt in first experiment; ssize:cohort number 
  
  pmat[1, 1:3] <- c(0, 0, 1)
  pmat <- pmat[rep(seq_len(nrow(pmat)), rep(A+1, nrow(pmat))),
               ]#copy 4 times; copy each original row of pmat for the times equals A+1; when there is one original row, there will be A+1 rows; when there is two 
                #original rows, there will be 2(A+1) rows in the end
  pmat[, "tox 1"] <- 0:A
  pmat[,"ssize 1"] <- A
  pmat[pmat[, "tox 1"] <= D, "d 2"] <- 1  #first step judging dlt in first experiment
  pmat[pmat[, "tox 1"] < C, "d 2"] <- 2
  pmat[pmat[, "tox 1"] > D, "stop"] <- 1
  pmat[pmat[, "tox 1"] > D, "desc"] <- 1
  stopped.pmat <- pmat[pmat$stop == 1, -2]#select matrix whihc is stopped and deleting the second column (desc) then
  
  if(dim(stopped.pmat)[1]>0){ #to ensure there is stopped path
    path.mat<-stopped.pmat # Edit from published threep3 - keeps record of all paths generated
    dose.mat <- stopped.pmat[, grep("d", names(stopped.pmat))] #select dose level from each experiment
    tox.mat <- stopped.pmat[, grep("tox", names(stopped.pmat))] #select dlt number from each experiment
    prob <- apply(matrix(dbinom(as.matrix(tox.mat), A, truep[as.matrix(dose.mat)]),
                         nrow = nrow(dose.mat)), 1, prod, na.rm = T) #calculate the prob #unclear yet
    ssize <- A * apply(!is.na(stopped.pmat[, grep("d", names(stopped.pmat))]),
                       1, sum)  #how much people used
    last.cohort <- apply(!is.na(stopped.pmat[, grep("d", names(stopped.pmat))]),
                         1, sum)#how many experiment of each possibility
    last.drug.column <- paste("d", last.cohort) #d1,d1; last step of experiment
    last.drug <- sapply(1:nrow(stopped.pmat), function(j) {
      stopped.pmat[j, last.drug.column[j]]
    })#j=1 to 2; the dose level at the last experiment
    previous.drug.column <- paste("d", last.cohort - 1)
    previous.drug <- sapply(1:nrow(stopped.pmat), function(j) {
      ifelse(previous.drug.column[j] == "d 0", 0, stopped.pmat[j,
                                                               previous.drug.column[j]])#
    })#the dlt before last experiment
    
    last.tox.column <- paste("tox", last.cohort)
    last.tox <- sapply(1:nrow(stopped.pmat), function(j) {
      stopped.pmat[j, last.tox.column[j]]
    })# to select the dlts of the last cohort
    mtd <- rep(NA, nrow(stopped.pmat))
    mtd[last.tox <=D & previous.drug == last.drug] <- last.drug[last.tox <=
                                                                  D & previous.drug == last.drug] - 1 # if last dlts<=D and two experiments already, set mtd
    mtd[last.tox <= D & previous.drug != last.drug] <- last.drug[last.tox <=
                                                                   D & previous.drug != last.drug]
    mtd[last.tox < C] <- last.drug[last.tox < C]
    mtd[last.tox > D] <- last.drug[last.tox > D] - 1  #work out mtd based on last.drug 
    exp <- sapply(1:doses, function(j) {
      sum(A * (stopped.pmat[, grep("d", names(stopped.pmat))] ==
                 j) * prob/ssize, na.rm = T)
    })
    
    dlt.no <- apply(stopped.pmat[, grep("tox", names(stopped.pmat))],
                    1, sum, na.rm = T) #whole dlt number
  }
  
  for (i in 3:mcplus1) {
    cat(paste(round(100 * i/mcplus1), "% complete\n", sep = ""))
    dd <- as.character(paste("d", i))
    td <- as.character(paste("tox", i))
    sd <- as.character(paste("ssize", i))
    dc <- as.character(paste("d", i - 1))
    tc <- as.character(paste("tox", i - 1))
    sc <- as.character(paste("ssize", i - 1))
    db <- as.character(paste("d", i - 2))
    tb <- as.character(paste("tox", i - 2))
    sb <- as.character(paste("ssize", i - 2))
    pmat1 <- pmat[rep(which(pmat[, "stop"] == 0 & pmat[, tb] >= C & pmat[,tb]<=D & pmat[, "desc"] == 0 & pmat[,dc]==pmat[,db]), each = B+1), ]
    #copy the middle way B+1 times; process the middle
    pmat2 <- pmat[rep(which(pmat[, "stop"] == 0 & pmat[, "desc"] == 0 & pmat[,dc]!=pmat[,db]), each = A+1),]
    #copy the way which is already increased, each A+1 times;
    if(dim(pmat1)[1]>0){
      pmat1[, tc] <- 0:B #tox=0 to B
      pmat1[, sc] <- B  #cohort=B
      pmat1[pmat1[, tc] + pmat1[,tb] <=E & pmat1[,dc]+1 <=doses & pmat1[,"desc"]==0,dd]<-pmat1[pmat1[, tc] + pmat1[,tb] <=E & pmat1[,dc]+1 <=doses & pmat1[,"desc"]==0,dc]+1
      #process 'are there more than E dlts' left way; ensure there is higher level; no descendent; then next dose level +1
      pmat1[pmat1[, tc] + pmat1[,tb] >E & pmat1[,dc] -1>=1 & pmat1[,"desc"]==0,dd]<-pmat1[pmat1[, tc] + pmat1[,tb] >E & pmat1[,dc] -1 >=1 & pmat1[,"desc"]==0,dc]-1
      #process 'are there more than E dlts' right way; esure there is lower level; no descendent; then dose level -1
      pmat1[pmat1[, tc] + pmat1[,tb] >E & pmat1[,dc] -1>=1 & pmat1[,"desc"]==0,"desc"]<-1
      #set desc as 1 because it's descendent; then next dose level = current level because need another cohort at same lavel
    }else{pmat1<-NULL}
    if(dim(pmat2)[1]>0){
      pmat2[,tc] <- 0:A
      pmat2[, sc] <- A
      pmat2[pmat2[,tc] < C & pmat2[,dc]+1<=doses & pmat2[,"desc"]==0, dd]<-pmat2[pmat2[,tc] < C & pmat2[,dc]+1<=doses & pmat2[,"desc"]==0, dc]+1
      #process 'how many dlts at dj', when dlts <C ; ensure there is higher level; no descendent; then next dose level +1
      pmat2[pmat2[,tc] >= C & pmat2[,tc] <=D & pmat2[,dc]<=doses & pmat2[,"desc"]==0, dd]<-pmat2[pmat2[,tc] >= C & pmat2[,tc]<=D & pmat2[,dc]<=doses & pmat2[,"desc"]==0, dc]
      #process 'how many dlts at dj', when dlts between C and D; ensure the current dose level <= highest level because it is already increased; no descendent; 
      #then next dose level = current level because need another cohort at same lavel
      pmat2[pmat2[,tc] > D & pmat2[,dc]-1>=1 & pmat2[,"desc"]==0, dd]<-pmat2[pmat2[,tc] > D & pmat2[,dc]-1>=1 & pmat2[,"desc"]==0, dc]-1
      #process 'how many dlts at dj' in case that dlts >D; ensure there is lower level; no descendent; then next level -1
      pmat2[pmat2[,tc] > D & pmat2[,dc]<=doses & pmat2[,"desc"]==0, "desc"]<-1
      #set desc as 1
    }else{pmat2<-NULL}
    
    pmat3 <- pmat[rep(which(pmat[, "stop"] == 0 & pmat[, "desc"] == 1), each=B+1),]#repeat the row which is desendent for B+1 times each
    if(dim(pmat3)[1]>0){
      
      pmat3[,tc]<-0:B
      pmat3[,sc]<-B
      otherc.drug<- sapply(1:nrow(pmat3), function(j) which(pmat3[,grep("d ", names(pmat3))][j,1:which(colnames(pmat3[,grep("d ", names(pmat3))])==db)]==pmat3[j,dc]))
      for(j in 1:length(otherc.drug)){if(length(otherc.drug[[j]])==0){otherc.drug[[j]]<-NA}}
      otherc.drug<-unlist(otherc.drug)
      for(j in 1:length(otherc.drug)){
        #if(!is.na(otherc.drug[j])==TRUE){
        #   pmat3[pmat3[,tc] + pmat3[,as.character(paste("tox", otherc.drug[j]))] > E & pmat3[,dc] - 1>=1, dd] <- pmat3[pmat3[,tc] + pmat3[,as.character(paste("tox", otherc.drug[j]))] > E & pmat3[,dc] - 1>=1, dc] - 1}
        #}
        
        if(!is.na(otherc.drug[j])==TRUE){
          if(pmat3[j,tc]+pmat3[j, as.character(paste("tox", otherc.drug[j]))] > E & pmat3[j,dc] - 1 >= 1){pmat3[j,dd]<-pmat3[j, dc] - 1}
          #process 'are there more than E dlts at dj' in case that dlts>E; ensure there is lower level; then next level-1
        }
      }
      #browser()
    }else{pmat3<-NULL}
    
    pmat<-rbind(pmat1,pmat2,pmat3)
    
    #if(i!=3){
    #da <- as.character(paste("d", i - 3))
    #ta <- as.character(paste("tox", i - 3))
    #sa <- as.character(paste("ssize", i - 3))
    #pmat3<-pmat[pmat[, "stop"] == 0 & pmat[, "desc"] == 1,]
    #if(dim(pmat3)[1]!=0){
    #otherc.drug<- sapply(1:nrow(pmat3), function(j) which(pmat3[,grep("d ",names(pmat3))][j,1:max(pmat3[,da])]==pmat3[j,db]))#column number
    #for(j in 1:length(otherc.drug)){if(length(otherc.drug[[j]])==0){otherc.drug[[j]]<-NA}}
    #otherc.drug<-unlist(otherc.drug)
    #pmat4<-pmat5<-NULL
    #for(j in 1:length(otherc.drug)){
    #           ifelse(!is.na(otherc.drug[j])==TRUE, ifelse(pmat3[j,tc] + pmat3[j,as.character(paste("tox", otherc.drug[j]))] > E, pmat4<-rbind(pmat4,pmat3[j,]), pmat5<-rbind(pmat5,pmat3[j,])), pmat5<-rbind(pmat5,pmat3[j,]))
    #}
    #            #browser()
    #pmat4null<-as.numeric(is.null(pmat4))+as.numeric(is.null(dim(pmat4)[1]))
    #dimpmat4null<-length(as.numeric(dim(pmat4)[1]))
    #pmat4 <- pmat4[rep(which(pmat4[, "stop"] == 0 & pmat4[, "desc"] == 1), each=B+1),]
    #if(pmat4null==0 & dimpmat4null>0){
    #pmat4[, tc] <- 0:B
    #pmat4[, sc] <- B
    #otherc.drug<- sapply(1:nrow(pmat3), function(j) which(pmat3[,grep("d ",names(pmat3))][j,1:max(pmat3[,db])]==pmat3[j,dc]))#column number
    #for(j in 1:length(otherc.drug)){if(length(otherc.drug[[j]])==0){otherc.drug[[j]]<-NA}}
    #otherc.drug<-unlist(otherc.drug)
    #                for(j in 1:length(otherc.drug)){
    #if(!is.na(otherc.drug[j])==TRUE){pmat4[pmat4[,tc] + pmat4[,as.character(paste("tox", otherc.drug[j]))] > E & pmat4[,dc] - 1>=1, dd] <- pmat4[pmat4[,tc] + pmat4[,as.character(paste("tox", otherc.drug[j]))] > E & pmat4[,dc] - 1>=1, dc] - 1}
    #}
    #pmat3[pmat3[,tc] + pmat[]]<-pmat3[]
    #}else{pmat4<-NULL}
    #pmat<-rbind(pmat1,pmat2, pmat4, pmat5)
    #}}else{pmat<-rbind(pmat1,pmat2)}
    
    
    #pmat[pmat[, tc] == 0 & pmat[, "desc"] == 0 & pmat[, dc] + 1 <= doses, dd] <- pmat[pmat[, tc] == 0 & pmat[,
    #"desc"] == 0 & pmat[, dc] + 1 <= doses, dc] + 1
    #pmat[pmat[, tc] == 1 & pmat[, "desc"] == 0 & pmat[, tb] == 0, dd] <- pmat[pmat[, tc] == 1 & pmat[, "desc"] ==
    #0 & pmat[, tb] == 0, dc]
    #pmat[pmat[, tc] == 1 & pmat[, "desc"] == 0 & pmat[, tb] == 1 & pmat[, dc] - 1 >= 1, dd] <- pmat[pmat[, tc] ==
    #1 & pmat[, "desc"] == 0 & pmat[, tb] == 1, dc] - 1
    #pmat[pmat[, tc] == 1 & pmat[, "desc"] == 0 & pmat[, tb] == 1 & pmat[, dc] - 1 >= 1, "desc"] <- 1
    #pmat[pmat[, tc] > 1 & pmat[, dc] - 1 >= 1, dd] <- pmat[pmat[, tc] > 1 & pmat[, dc] - 1 >= 1, dc] - 1
    #pmat[pmat[, tc] > 1 & pmat[, dc] - 1 >= 1, "desc"] <- 1
    
    excluding.dd <- names(pmat)[grepl("d ", names(pmat)) & names(pmat) != dd] #all d except di
    cnt <- apply(pmat[!is.na(pmat[, dd]), dd] == pmat[!is.na(pmat[, dd]), excluding.dd], 1, sum, na.rm = T) #pick the row whose d3 is not a missing value, compare
    #the unmissing values to d3(dd), if they equal dd, add up.# whether a dose level has been tested twice
    pmat[!is.na(pmat[, dd]), dd][cnt > 1] <- NA #process 'how many patients given dj-1 when there is A+B patients'
    pmat[is.na(pmat[, dd]), "stop"] <- 1 #stop
    stopped.pmat <- pmat[pmat$stop == 1, -2]
    
    if(dim(stopped.pmat)[1]>0){
      path.mat<-rbind(path.mat,stopped.pmat) # Edit from published threep3 - keeps record of all paths generated
      #browser()
      dose.mat <- stopped.pmat[, grep("d ", names(stopped.pmat))[1:(i - 1)]]
      tox.mat <- stopped.pmat[, grep("tox", names(stopped.pmat))[1:(i - 1)]]
      ssize.mat <- stopped.pmat[, grep("ssize", names(stopped.pmat))[1:(i - 1)]]
      prob.new <- apply(matrix(dbinom(as.matrix(tox.mat), as.matrix(ssize.mat),
                                      truep[as.matrix(dose.mat)]), nrow = nrow(dose.mat)),1, prod, na.rm = T)
      prob <- c(prob, prob.new)
      ssize.new <- apply(ssize.mat,1,sum)#rep(3 * (i - 1), nrow(stopped.pmat))
      ssize <- c(ssize, ssize.new)
      last.drug <- stopped.pmat[, dc]
      previous.drug <- stopped.pmat[, db]
      otherc.drug <- sapply(1:nrow(stopped.pmat), function(j) which(stopped.pmat[,grep("d ", names(stopped.pmat))][j,1:which(colnames(stopped.pmat[,grep("d ", names(stopped.pmat))])==db)]==stopped.pmat[j,dc]))
      #otherc.drug <- sapply(1:nrow(stopped.pmat), function(j) which(stopped.pmat[,grep("d",names(stopped.pmat))][j,1:max(stopped.pmat[,db])]==stopped.pmat[j,dc]))#column number
      for(j in 1:length(otherc.drug)){if(length(otherc.drug[[j]])==0){otherc.drug[[j]]<-NA}}
      otherc.drug<-unlist(otherc.drug)
      last.tox <- stopped.pmat[, tc]
      previous.tox <- stopped.pmat[, tb]
      otherc.tox<-rep(NA,length(otherc.drug))
      for(j in 1:length(otherc.tox)){ifelse(is.na(otherc.drug[j])==TRUE,otherc.tox[j]<-NA,otherc.tox[j]<-stopped.pmat[j,grep("tox", names(stopped.pmat))][,otherc.drug[j]])}
      #otherc.drug<- sapply(1:nrow(stopped.pmat), function(j) which(stopped.pmat[,grep("tox",names(stopped.pmat))][j,1:max(stopped.pmat[,dc])-1]==stopped.pmat[j,dc]))
      #for(i in 1:length(otherc.drug)){if(length(otherc.drug[[i]])==0){otherc.drug[[i]]<-NA}}
      #otherc.drug<-unlist(otherc.drug)
      last.ssize <- stopped.pmat[, sc]
      previous.ssize <- stopped.pmat[, sb]
      mtd.new <- rep(NA, nrow(stopped.pmat))
      #browser()
      mtd.new[last.tox < C & last.drug!=previous.drug] <- last.drug[last.tox < C & last.drug!=previous.drug]
      mtd.new[last.tox + previous.tox > E & previous.drug == last.drug] <- last.drug[last.tox + previous.tox > E & previous.drug == last.drug] - 1
      mtd.new[last.tox + previous.tox <= E & previous.drug == last.drug] <- last.drug[last.tox + previous.tox <= E & previous.drug == last.drug]
      mtd.new[last.tox > D & previous.drug != last.drug] <- last.drug[last.tox > D & previous.drug != last.drug] - 1
      for(j in 1:length(otherc.tox)){
        if(!is.na(otherc.tox[j])==TRUE){mtd.new[last.tox + otherc.tox <= E & last.drug == dose.mat[,otherc.drug[j]]] <- last.drug[last.tox + otherc.tox <= E & last.drug == dose.mat[,otherc.drug[j]]]}
      }
      #mtd.new[last.tox  1 & previous.drug != last.drug] <- last.drug[last.tox == 1 & previous.drug != last.drug]
      mtd <- c(mtd, mtd.new)
      exp <- exp + sapply(1:doses, function(j) {
        sum((stopped.pmat[,grep("ssize", names(stopped.pmat))]) * (stopped.pmat[, grep("d", names(stopped.pmat))] == 1) * prob.new/ssize.new, na.rm = T) })
      dlt.no <- c(dlt.no, apply(stopped.pmat[, grep("tox",
                                                    names(stopped.pmat))], 1, sum, na.rm = T))
    }
  }
  path.mat<-cbind(path.mat,prob,ssize,mtd) 
  return(path.mat)
}
 

