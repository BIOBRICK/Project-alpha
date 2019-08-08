# TCGA mRNA-Seq data preprocess

In order to make TCGA Data into pesudotime series data, we use cancer "Clinical Stage" as time tag. Thus, we have to make the expression matrix clean together with clinical data: we use GDC Data transfer tool and R to do this:

```bash
gdc-client.exe download -m manifest.txt 
```

```R
# R script
#Merge clinical data
Your_work_path = 'Your_work_path'

setwd(Your_work_path)
library(dplyr)
library(hash)
library(GDCRNATools)
options(stringsAsFactors = F)


clinicaldataPath = 'Your_clinical_data_path'
clinicalDa <- gdcClinicalMerge(path = clinicaldataPath, key.info = T)
Patient_id<-row.names(clinicalDa)
rownames(clinicalDa)<-NULL
clinicalDa<-cbind(Patient_id,clinicalDa)

#输出之前最好确认一下第6列确实是stage，一般来说没问题，但是有时候也有可能第六列不是的，这时候要改成对应的列数。
clinical<-clinicalDa[,c(1,6)]
names(clinical)<-c('Id','stage')

write.csv(clinical,file = "clinical.csv")

#为了进行后续生存分析，这里另外提取一组生存信息用的表：
clinical_vital<-subset(clinicalDa,select = c("Patient_id","days_to_death","days_to_last_followup","vital_status","pathologic_stage"))
write.csv(clinical_vital,file = "clinical_vital.csv")

##################################################################################################
#Merge RNA data and clinical data.
clinicalFile="clinical.csv"
expFile="symbol.txt"

rt=read.table(clinicalFile,header=T,check.names=F,sep=",")

#transform clinical data into hash table
h = hash(keys = rt$Id, values = paste(rt$stage))
exp=read.table(expFile,header=T,check.names=F,sep="\t")

#eradicate replicated gene symbol, use max data to substitude.
exp<- exp %>%
  group_by(id) %>%
  summarise_all(max)

exprownames<-exp$id

exp<-exp[,-1]
row.names(exp)<-exprownames

write.table("sample\tStage",file="survivalInput.txt",sep="\t"
            ,quote=F,row.names=F,col.names=F)
for(i in names(exp)){
  j=unlist(strsplit(i,"\\-"))
  if(grepl("^0",j[4])){
    name4=paste(j[1],j[2],j[3],j[4],j[5],j[6],j[7],sep="-")
    name3=paste(j[1],j[2],j[3],sep="-")
    if(has.key(name3,h)){
      write.table(paste(name4,h[[name3]],sep="\t"),file="survivalInput.txt",sep="\t",
                  quote=F,append=TRUE,row.names=F,col.names=F)
    }
  }else if(grepl("^1",j[4])){
    name4=paste(j[1],j[2],j[3],j[4],j[5],j[6],j[7],sep="-")
    name3=paste(j[1],j[2],j[3],sep="-")
    if(has.key(name3,h)){
      write.table(paste(name4,"Normal",sep="\t"),file="survivalInput.txt",sep="\t",
                  quote=F,append=TRUE,row.names=F,col.names=F)
    }
  }
}

expt<-as.data.frame(t(exp))
exptrownames<-row.names(expt)
expt2<-cbind(exptrownames,expt)
stagedata<-read.table(file = "survivalInput.txt",header = T,sep = "\t")
exp_clinc<-merge(expt2,stagedata,by.x = "exptrownames",by.y = 'sample')
row.names(exp_clinc)<-exp_clinc$exptrownames
exp_clinc2<-exp_clinc[,-1]
stage<-exp_clinc2$Stage
final<-subset(exp_clinc2,select = -Stage)
final<-cbind(stage,final)

#Patient_id<-rownames(final)
#final<-cbind(Patient_id,final)
#rownames(final)<-NULL

final<-filter(final,final$stage!='NA')
final<-as.data.frame(t(final))
write.table(final,file = 'exp_with_stage.txt',sep = '\t',col.names = FALSE)

rm(list = ls())

```

Here is a very important poblem, if we want to open the final matrix table an have a look(no matter what reason). It's possible that we will use Excel-like software. This however, will cause a very sever problem !

Such action will make some of your symbol ID be transferred into "Date" format and will never be able to get them back !

So, I suggest that you use [escape_excel](https://github.com/pstew/escape_excel) to handle this problem, and to be honest, this is the **only** method......

Cite: *Welsh E A , Stewart P A , Kuenzi B M , et al. Escape Excel: A tool for preventing gene symbol and accession conversion errors[J]. PLOS ONE, 2017, 12(9):e0185207-.*

Here I have to say:

****ING EXCEL!!!

