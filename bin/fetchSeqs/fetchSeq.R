#install the package
#source("https://bioconductor.org/biocLite.R")
#biocLite("SRAdb")

library(SRAdb)
#get sradblite, this file update regularly, get online version
sqlfile <- 'SRAmetadb.sqlite'
if(!file.exists()) sqlfile <<- getSRAdbFile
dbcon <- dbConnect(SQLite(), sqlfile)
sra_tables <- dbListTables(dbcon)
sra_tables;
mydata <- read.table("input.txt")
studies <- names(table(mydata))

for (i in 1:length(studies)){
 #print(studies[i])
 projid <- studies[i]
 print(projid)
 dir.create(projid)
 #select * from sra, sample where study_alias='PRJNA264976' and sra.sample_accession=sample.sample_accession
 query <- paste("select * from sra, sample where study_alias='",projid , "'and sra.sample_accession=sample.sample_accession;", sep="")
 print(query)
 rs <- dbGetQuery(dbcon, query)
 #print(rs)
 print(colnames(rs))
 write.table(rs,file=paste(projid, "/table.", studies[i], ".txt", sep=""), sep="\t", row.names=FALSE)
 #sra_address <- listSRAfile(rs$run_accession, dbcon, fileType='sra')
 sra_info <- getSRAinfo(rs$run_accession,dbcon, sraType='sra')
 print(sra_info)
 #you have options to download fastq or sra file
 getSRAfile(rs$run_accession, dbcon, fileType='sra', destDir=projid)
}
