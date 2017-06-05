motifid=read.csv("ATtRACT_db.txt",stringsAsFactors=F,sep="\t",skip=0,header = TRUE)
motifid=cbind(motifid$Gene_name,motifid$Organism,motifid$Matrix_id)
motifid=motifid[!duplicated(motifid),]
xx=which(motifid[,2]=="Homo_sapiens")
Human=motifid[xx,]
RBP_gene=c()
Biofunction <- file("RBP_result4.txt", "r")
line=readLines(Biofunction,n=1)
i=0
while( length(line) != 0 ) {
  i <- i+1
  print(i)
  msig.go <- strsplit(line,"\t")
  tt=msig.go[[1]]
  ID=msig.go[[1]][1]
  xa=which(Human[,3]==ID)
  if(length(xa)>0){
    RT=cbind(Human[xa,1],tt[3],tt[4],tt[5],tt[8])
    RBP_gene=rbind(RBP_gene,RT)
  }
  line=readLines(Biofunction,n=1);
}
close(Biofunction)
colnames(RBP_gene)=c("RBP_gene","Entrez","start","end","p_value")
write.table(RBP_gene,"RBP_gene_binding.txt",sep = "\t",row.names = FALSE, col.names = T,quote = FALSE)