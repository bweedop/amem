library(ape)

csv.to.nexus<-function(csv){
  spreadsheet<-read.csv(csv, row.names = 1, header = TRUE)
  spreadsheet<-as.matrix(spreadsheet)
  colnames(spreadsheet)<-seq(1, ncol(spreadsheet))
  write.nexus.data(spreadsheet, "character_matrix.nex", format = 'dna')
  nexus.file <- "character_matrix.nex"
  for (f in nexus.file){
    x <- readLines(f)
    y <- gsub("DATATYPE=DNA", "DATATYPE=STANDARD", x)
    cat(y, file = f, sep = "\n")
  }
}