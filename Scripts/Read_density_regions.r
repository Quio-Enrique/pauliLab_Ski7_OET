#######################
# Luis Enrique Cabrera Quio
# Script to calculate relative read distribution of 5' UTRs, 3' UTRs, and CDSs per transcript
#######################

library(ggplot2)

maindir<-"ReadDensityRegion/"
samp<- c("Oocytes", "Egg", "Embryos")
gt<-c("WT", 'Ski7')
reps<-c(paste0("R", 1:3))

for(s in samp){
  if(s == "Embryos" | s == "Oocytes"){ tps<-c(paste0("T", 1:4)) } else { tps<-c("ACT", "INA", "FER") }
    for(p in tps){
      for(g in gt){
        for(r in reps){
          five<-read.table(paste0(maindir, s, '/', g, '/', p, '/', r, '/', '<OUTht_5UTR>'), sep = '\t', stringsAsFactors = F)
          cds<-read.table(paste0(maindir, s, '/', g, '/', p, '/', r, '/', 'OUTht_CDS'), sep = '\t', stringsAsFactors = F)
          three<-read.table(paste0(maindir, s, '/', g, '/', p, '/', r, '/', 'OUTth_3UTR'), sep = '\t', stringsAsFactors = F)
          comptrx<-intersect(intersect(five$V1, cds$V1), three$V1)
          comptrx<-data.frame(Trx = comptrx, Leader = five$V2[match(comptrx, five$V1)], CDS = cds$V2[match(comptrx, cds$V1)], Trailer = three$V2[match(comptrx, three$V1)])
          comptrxtotal<-rowSums(comptrx)
          dir.create(paste0(maindir, s, '/', g, '/', p, '/', r), recursive = T)
          write.table(comptrxtotal, paste0(maindir, s, '/', g, '/', p, '/', r, '/Read_number_', g, '_', p, '_', r), sep = '\t', quote = F)
          tmpdir<-'DEGs/'
          up<-read.table(paste0(tmpdir, s, '/', p, '/Tables/Up_', s, '_DEGs_alph01_', p, '.txt'), sep = '\t', stringsAsFactors = F, header = T)[,1]
          down<-read.table(paste0(tmpdir, s, '/', p, '/Tables/Down_', s, '_DEGs_alph01_', p, '.txt'), sep = '\t', stringsAsFactors = F, header = T)[,1]
          unc<-read.table(paste0(tmpdir, s, '/', p, '/Tables/Unchanged_', s, '_DEGs_alph01_', p, '.txt'), sep = '\t', stringsAsFactors = F, header = T)[,1]
          updist<-avgreps[na.omit(match(paste0("gene:", up), rownames(avgreps))),]
          updist[is.na(updist)]<-0
          updist<-melt(updist)
          downdist<-avgreps[na.omit(match(paste0("gene:", down), rownames(avgreps))),]
          downdist[is.na(downdist)]<-0
          downdist<-melt(downdist)
          uncdist<-avgreps[na.omit(match(paste0("gene:", unc), rownames(avgreps))),]
          uncdist[is.na(uncdist)]<-0
          uncdist<-melt(uncdist)
          all<-rbind(updist, downdist, uncdist)
          forsum<-match(unique(all$Var1), comptrx))
          totread<-readnum[forsum,]
          names(totread)<-unique(all$Var1)
          write.table(totread, paste0(maindir, s, '/', g, '/', p, '/', r, '/Read_number_used_trxs_', g, '_', p, '_', r), sep = '\t', quote = F)
        }
      }
    }
  }
}

# Calculate library size from the transcripts used to use as a normalising factor for the plots

nreadall<-vector()
for(s in samp){
  if(s == "Embryos" | s == "Oocytes"){ tps<-c(paste0("T", 1:4)) } else { tps<-c("ACT", "INA", "FER") }
  for(p in tps){
    for(g in gt){
      for(r in reps){
        if(g == 'WT' & (p == 'INA' | p == "FER") & (r == "R3")){ next
        }else{
          readsnall<-read.table(paste0(maindir, s, '/', g, '/', p, '/', r, '/Read_number_used_trxs_', g, '_', p, '_', r), sep = '\t', stringsAsFactors = F, row.names = 1)
          nreadall<-c(nreadall, colSums(readsnall))
          names(nreadall)<-c(names(nreadall)[1:length(nreadall)-1], paste(s, p, g, r, sep = '_'))
        }
      }
    }
  }
  nfactall<-min(nreadall)/nreadall
  for(p in tps){
    for(g in gt){
      for(r in reps){
        five<-read.table(paste0(maindir, s, '/', g, '/', p, '/', r, '/', '<OUTht_5UTR>'), sep = '\t', stringsAsFactors = F)
        cds<-read.table(paste0(maindir, s, '/', g, '/', p, '/', r, '/', 'OUTht_CDS'), sep = '\t', stringsAsFactors = F)
        three<-read.table(paste0(maindir, s, '/', g, '/', p, '/', r, '/', 'OUTth_3UTR'), sep = '\t', stringsAsFactors = F)
        comptrx<-intersect(intersect(five$V1, cds$V1), three$V1)
        comptrx<-data.frame(Trx = comptrx, UTR5 = five$V2[match(comptrx, five$V1)], CDS = cds$V2[match(comptrx, cds$V1)], UTR3 = three$V2[match(comptrx, three$V1)])
        rownames(comptrx)<-comptrx[,1]
        comptrx<-comptrx[,-1]
        if(g == 'WT' & (p == 'INA' | p == "FER") & (r == "R3")){ nfall<-1
        } else { nfall<-nfactall[which(names(nfactall) == paste(s, p, g, r, sep = '_'))] }
        comptrxall<-comptrx*nfall
        fivelen<-read.table('UTR5_length.txt', sep = '\t', stringsAsFactors = F, header = T)
        cdslen<-read.table('CDS_length.txt', sep = '\t', stringsAsFactors = F, header = T)
        threelen<-read.table('UTR3_length.txt', sep = '\t', stringsAsFactors = F, header = T)
        tmp<-intersect(intersect(fivelen$Transcript, cdslen$Transcript), threelen$Transcript)
        lenlong<-data.frame(Trx = tmp, UTR5 = fivelen$Length[match(tmp, fivelen$Transcript)], CDS = cdslen$Length[match(tmp, cdslen$Transcript)], UTR3 = threelen$Length[match(tmp, threelen$Transcript)])
        rownames(lenlong)<-lenlon$Trx
        lenlong<-lenlong[,-1]
        densiall<-comptrxlongall/lenlong
        densiall[is.na(densiall)]<-0
        write.table(densiall, paste0(maindir, s, '/', g, '/', p, '/', r, '/Density_prop_nfact_all_', g, '_', p, '_', r), sep = '\t', quote = F)
      }
      r1<-as.matrix(read.table(paste0(maindir, s, '/', g, '/', p, '/R1/Density_prop_nfact_all_', g, '_', p, '_R1'), stringsAsFactors = F, sep = '\t', header = T))
      r2<-as.matrix(read.table(paste0(maindir, s, '/', g, '/', p, '/R2/Density_prop_nfact_all_', g, '_', p, '_R2'), stringsAsFactors = F, sep = '\t', header = T))
      r3<-as.matrix(read.table(paste0(maindir, s, '/', g, '/', p, '/R3/Density_prop_nfact_all_', g, '_', p, '_R3'), stringsAsFactors = F, sep = '\t', header = T))
      if(g == 'WT' & (p == 'INA' | p == "FER")){ replist<-list(r1, r2) } else { replist<-list(r1, r2, r3) }
      repbi<-do.call(cbind, replist)
      repbi<-array(repbi, dim = c(dim(replist[[1]]), length(replist)))
      avgreps<-apply(repbi, c(1, 2), mean, na.rm = TRUE)
      colnames(avgreps)<-c("Leader", "CDS", "Trailer")
      rownames(avgreps)<-rownames(r1)
      dir.create(paste0(maindir, s, '/', g, '/', p, '/Rs'), recursive = T)
      write.table(avgreps, paste0(maindir, s, '/', g, '/', p, '/Rs/Density_prop_nfact_all_average_Reps'), sep = '\t', quote = F)
      tmpdir<-'DEGs/'
      up<-read.table(paste0(tmpdir, s, p, '/Tables/Up_', s, '_DEGs_alph01_', p, '.txt'), sep = '\t', stringsAsFactors = F, header = T)[,1]
      down<-read.table(paste0(tmpdir, s, p, '/Tables/Down_', s, '_DEGs_alph01_', p, '.txt'), sep = '\t', stringsAsFactors = F, header = T)[,1]
      unc<-read.table(paste0(tmpdir, s, p, '/Tables/Unchanged_', s, '_DEGs_alph01_', p, '.txt'), sep = '\t', stringsAsFactors = F, header = T)[,1]
      updist<-avgreps[na.omit(match(paste0("gene:", up), rownames(avgreps))),]
      updist[is.na(updist)]<-0
      updist<-melt(updist)
      updist$DEG <- rep("Up")
      updist$GT <- rep(g)
      downdist<-avgreps[na.omit(match(paste0("gene:", down), rownames(avgreps))),]
      downdist[is.na(downdist)]<-0
      downdist<-melt(downdist)
      downdist$DEG <- rep("Down")
      downdist$GT <- rep(g)
      uncdist<-avgreps[na.omit(match(paste0("gene:", unc), rownames(avgreps))),]
      uncdist[is.na(uncdist)]<-0
      uncdist<-melt(uncdist)
      uncdist$DEG <- rep("Unch")
      uncdist$GT <- rep(g)
      all<-rbind(updist, downdist, uncdist)
      write.table(all, paste0(maindir, s, '/' g, '/', p, '/Rs/Density_prop_nfact_all_', g, '_avgReps_DEGs'), sep = '\t', quote = F, row.names = F)
    }
  }
}
					
# Calculate the ratio of ski7/WT of every gene body region

for(s in samp){
  if(s == "Embryos" | s == "Oocytes"){ tps<-c(paste0("T", 1:4)) } else { tps<-c("ACT", "INA", "FER") }
  tpsratio<-c()
  wilsall<-c()
  for(p in tps){
    wt<-read.table(paste0(maindir, s, '/WT/', p, '/Rs/Density_prop_nfact_all_WT_avgReps_DEGs'), sep = '\t', stringsAsFactors = F, header = T)
    ski<-read.table(paste0(maindir, s, '/Ski7/', p, '/Rs/Density_prop_nfact_all_Ski7_avgReps_DEGs'), sep = '\t', stringsAsFactors = F, header = T)
    ratio<-data.frame(Gene = wt$Var1, Region = wt$Var2, Ratio = ski$value/wt$value, DEG = wt$DEG)
    ratio<-ratio[-(which(is.na(ratio$Ratio))),]
    ratio<-ratio[-(which(is.infinite(ratio$Ratio))),]
    ratio$Region<-factor(ratio$Region, levels = c("UTR5", "CDS", "UTR3"))
    ratio$DEG<-factor(ratio$DEG, levels = c("Up", "Down"))
    meanlu<-mean(ratio$Ratio[which(ratio$Region == "UTR5" & ratio$DEG == "Up")])
    meancu<-mean(ratio$Ratio[which(ratio$Region == "CDS" & ratio$DEG == "Up")])
    meantu<-mean(ratio$Ratio[which(ratio$Region == "UTR3" & ratio$DEG == "Up")])
    meanld<-mean(ratio$Ratio[which(ratio$Region == "UTR5" & ratio$DEG == "Down")])
    meancd<-mean(ratio$Ratio[which(ratio$Region == "CDS" & ratio$DEG == "Down")])
    meantd<-mean(ratio$Ratio[which(ratio$Region == "UTR3" & ratio$DEG == "Down")])
    meanreg<-rbind(cbind(meanlu, meancu, meantu), cbind(meanld, meancd, meantd))
    tpsratio<-cbind(tpsratio, meanreg)
    wiltestcdsup<-wilcox.test(ratio$Ratio[which(ratio$Region == "UTR5" & ratio$DEG == "Up")], ratio$Ratio[which(ratio$Region == "CDS" & ratio$DEG == "Up")])$p.value
    wiltestcdsdo<-wilcox.test(ratio$Ratio[which(ratio$Region == "UTR5" & ratio$DEG == "Down")], ratio$Ratio[which(ratio$Region == "CDS" & ratio$DEG == "Down")])$p.value
    wiltesttraup<-wilcox.test(ratio$Ratio[which(ratio$Region == "UTR5" & ratio$DEG == "Up")], ratio$Ratio[which(ratio$Region == "UTR3" & ratio$DEG == "Up")])$p.value
    wiltesttrado<-wilcox.test(ratio$Ratio[which(ratio$Region == "UTR5" & ratio$DEG == "Down")], ratio$Ratio[which(ratio$Region == "UTR3" & ratio$DEG == "Down")])$p.value
    wils<-rbind(c(0, wiltestcdsup, wiltesttraup), c(0, wiltestcdsdo, wiltesttrado))
    wilsall<-cbind(wilsall, wils)
  }
  colnames(wilsall)<-colnames(tpsratio)<-paste0(c("5UTR_", "CDS_", "3UTR_"), rep(1:4, each = 3))
  rownames(tpsratio)<-rownames(tpsratio)<-c("Up", "Down")
  wilsall<-melt(wilsall)
  tpsratio<-melt(tpsratio)
  tpsratio$pval <- wilsall$value
  upsi<-ggplot(tpsratio[which(tpsratio$Var1 == "Up"),], aes(Var2, (as.numeric(value)))) + geom_bar() + geom_text(aes(label = scientific(pval, digits = 2)))
  dosi<-ggplot(tpsratio[which(tpsratio$Var1 == "Down"),], aes(Var2, (as.numeric(value)))) + geom_bar() + geom_text(aes(label = scientific(pval, digits = 2)))
  print(ggarrange(upsi, dosi, ncol = 1))
}
