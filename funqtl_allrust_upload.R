#load libraries
library(readr)
library(ggplot2)
library(data.table)
library(funqtl)
library(reshape2)

setwd("~/Desktop/pavi_newmap/Pathogen/all_rust")

#load cleaned data
load("~/Desktop/pavi_newmap/Rdata/csvs_clean2_june13.rda")


#creating cross objects
load("~/Desktop/pavi_newmap/Rdata/pv2017_cross_out.rda")
crosstest<-cross_out
#clean out old phenotypes
crosstest$pheno[,13:20]<-NULL

#create function to add phenos to the cross
pheno_add<-function(rustpheno, crosstest){
  #id column for cross
  cross.ids = crosstest$pheno$PLOT_GL 
  #id column for each site
  rust.ids = rustpheno$PLOT_GL 
  pheno2cross = rustpheno
  pheno2cross_order = pheno2cross[order(pheno2cross$PLOT_GL),]
  #rownames(pheno2cross)<-rust.ids
  pheno2crossmini = pheno2cross_order[pheno2cross_order$PLOT_GL %in% cross.ids,]
  phenomini<-as.data.frame(pheno2crossmini)
  #remove duplicates
  phenomini2<-phenomini[!duplicated(phenomini$PLOT_GL),]
  #spread to cross.ids
  phenofinal = phenomini2[match(cross.ids,phenomini2$PLOT_GL),]
  crosstest$pheno<-data.frame(crosstest$pheno, phenofinal)
  crosstest
}

#your data starts on column 15

crosslist_all<-lapply(csvs_clean2, function(x){pheno_add(x,crosstest)})

#test that it worked. all columns should have nonzero values
sapply(crosslist_all, function(x){sum(x$pheno[,15:17], na.rm=T)})

dates<-lapply(crosslist_all, FUN=function(x){names(x$pheno[,15:ncol(x$pheno)])})
times<-lapply(dates, FUN=function(x){substr(x, 6,8)})
timesnumeric<-lapply(dates, FUN=function(x){as.numeric(substr(x, 6,8))})

##funqtl
crosslist_keep<-crosslist_all

#position of desired phenos
phes <- lapply(crosslist_keep, FUN=function(x){seq(15, nphe(x), by=1)})
#calc genotype probabilities ~5m
genoprobs <- lapply(crosslist_keep, FUN=function(x){calc.genoprob(x, step=0)})

#scanone 
outs <- vector("list", length(phes))
for (i in 1:length(phes)){
  outs[[i]]<- scanone(genoprobs[[i]], pheno.col = phes[[i]], method="hk")}
#est qtl effects at each time for each locus
effs <- mapply(FUN=function(x,y){geteffects(x, pheno.cols=y)}, x=genoprobs, y=phes)

csv_names<-names(crosslist_all)
timesnumeric_keep<-timesnumeric
sitenames_yr_keep<-as.list(substr(csv_names, 1,9))
sitenames_keep<-substr(csv_names, 1,4)

#test plots
plot(outs[[1]])
plot(outs[[12]])
plot(outs[[21]])

#plot
pdf(file="effectsplots_rustall.pdf")
mapply(FUN=function(x,y,z,t){
  plotlod(output=x, effects=y, y=z, gap=15,
          ylab=paste("Julian Day",t), ylim=c(100:150))
}, x=outs, y=effs, z=timesnumeric_keep, t=sitenames_yr_keep)
dev.off()

##SLOD and MLOD
outsscanoneF <- lapply(genoprobs, FUN=function(x){scanoneF(x, pheno.cols = 15:nphe(x), method="hk")})
#save(outsscanoneF, file="outsscanoneF_july8.Rdata")
#check SLOD and MLOD
pdf(file="SLODMLOD_test.pdf")
mapply(FUN=function(x,
                    #y,
                    t){
  par(mfrow=c(2,1))
  plot(x, 
       ylim=c(0,8), 
       main=paste("The SLOD curve for",t),
       bandcol="gray90")
   plot(x, lodcolumn=2, 
       ylim=c(0,15),
       main=paste("The MLOD curve for",t), 
       bandcol="gray90")
}, x=outsscanoneF, 
t=sitenames_yr_keep)
dev.off()

#slow step to get permutation threshold
permouts <- lapply(crosslist_keep, FUN=function(x){
  scanoneF(x, pheno.cols=15:(nphe(x)),method = "hk", n.perm=1000)
})
#save(permouts, file="permouts_rustall.Rdata")

permoutsums<-lapply(permouts, FUN=function(x){summary(x)})
permoutsums

#save(permoutsums, file = "permoutsums_rustall_july8.Rdata")

#estimate allele effects
allele_fx <- vector("list", length(outsscanoneF))
for (i in 1:length(outsscanoneF)){
  allele_fx[[i]]<- summary(outsscanoneF[[i]], format="tabByCol", pvalues=TRUE, perm=permouts[[i]], ci.function="lodint")}
names(allele_fx)<-names(outsscanoneF)

pdf(file="SLODMLOD_all.pdf")
mapply(FUN=function(x,
                    y,
                    t){
  par(mfrow=c(2,1))
  plot(x, 
       ylim=c(0,4.5), 
       main=paste("The SLOD curve for",t),
       bandcol="gray90")
  abline(h=y[1,1], col="red", lty=3)
  plot(x, lodcolumn=2, 
       ylim=c(0,15),
       main=paste("The MLOD curve for",t), 
       bandcol="gray90")
  abline(h=y[1,2], col="red", lty=3)
}, x=outsscanoneF, 
y=permoutsums, 
t=sitenames_yr_keep)
dev.off()

#clean up outputs for plotting
outplotdf<-lapply(outsscanoneF, FUN=function(x){as.data.frame(x)})
outplot1<-Map(cbind, outplotdf, site=sitenames_keep, year=substr(sitenames_yr_keep,6,9))
plotmax<- vector("list", length(genoprobs))
for(i in 1:length(genoprobs))
  plotmax[i]<-as.data.frame(ifelse(outplot1[[i]]$slod>permoutsums[[i]][1,1],outplot1[[i]]$slod, permoutsums[[i]][1,1]))
outplot2<-Map(cbind, outplot1, siglod=plotmax)
allplots<-do.call(rbind, outplot2)
##fix the order
allplots$site = factor(allplots$site, levels=c('KBSM','LINC','MNHT','CLMB','OVTN','TMPL','PKLE','BFLG', 'KING'))

##combined qtl plot
library(tidyverse)

permnames<-vector("list", 21)
for(i in 1:21){permnames[[i]]<-cbind(as.data.frame(permouts[[i]]),names(permouts)[i])}
permnames_long<-do.call(rbind, permnames)
colnames(permnames_long)[3]<-"site_yr"
permnames_long$rank<-rank(permnames_long$slod)
#take the 5% * 21
crit<-permnames_long[which(permnames_long$rank==(.05*2100)),1]*21

allplots_site<-allplots %>%
  group_by(chr,pos) %>%
  summarise(slod=sum(slod),
            mlod=sum(mlod))

ggplot(allplots_site, aes(x=pos, y=slod))+
  geom_line()+
  geom_hline(aes(yintercept=crit), lty="dashed")+
  facet_grid(.~chr, scales="free", space="free")+
  theme_classic(base_size = 15)+
  labs(x="Position", y="SLOD")+
  scale_x_continuous(breaks=NULL)
#ggsave("combined_slod.png", width=10, height=5)

library(grid)
library(gridExtra)

#split ny N and S
norths<-c(1:6,10:13)
souths<-c(7:9, 14:21)
perm_north<-vector("list", 10)
perm_south<-vector("list", 11)
for(i in 1:10){
  perm_north[[i]]<-permouts[[norths[i]]]
  perm_north[[i]]<-cbind(as.data.frame(perm_north[[i]]),names(permouts)[norths[i]])
}
perm_north_long<-do.call(rbind, perm_north)
colnames(perm_north_long)[3]<-"site_yr"
perm_north_long$rank<-rank(perm_north_long$slod)
#find the critical value
crit_north<-permnames_long[which(perm_north_long$rank==(.95*length(perm_north_long$rank))),1]*length(perm_north)

for(i in 1:11){
  perm_south[[i]]<-permouts[[souths[i]]]
  perm_south[[i]]<-cbind(as.data.frame(perm_south[[i]]),names(permouts)[souths[i]])
}
perm_south_long<-do.call(rbind, perm_south)
colnames(perm_south_long)[3]<-"site_yr"
perm_south_long$rank<-rank(perm_south_long$slod)
#find the critical value
crit_south<-permnames_long[which(perm_south_long$rank==(.95*length(perm_south_long$rank))),1]*length(perm_south)

allplots2<-allplots
allplots2$ns<-plyr::mapvalues(allplots2$site, 
                        c("KBSM", "LINC", "MNHT", "CLMB", "OVTN", "TMPL" ,"PKLE" , "KING"),
                        c("N","N","N","N","S","S","S","S"))
#detach(package:plyr)
allplots_site2<-allplots2 %>%
  group_by(chr,pos, ns) %>%
  summarise(slod=sum(slod),
            mlod=sum(mlod))

crits<-data.frame(crit=c(crit_north, crit_south), ns=c("North", "South"))

ggplot(allplots_site2, aes(x=pos, y=slod))+
  geom_line()+
  #geom_hline(data=crits,aes(yintercept=crit), lty="dashed")+
  facet_grid(ns~chr, space="free", scales="free")+
  theme_classic(base_size = 15)+
  labs(x="Position", y="SLOD")+
  scale_x_continuous(breaks=NULL)

