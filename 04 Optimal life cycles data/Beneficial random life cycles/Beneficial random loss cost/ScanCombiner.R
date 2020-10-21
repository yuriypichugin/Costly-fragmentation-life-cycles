setwd('/Volumes/OSX/Users/pichugin/Projects/Fragmentation with costs/Data/Beneficial random loss cost 2020/')

# FileList=list.files(pattern="^[D]") #here is regular expression, which I don't understand
FileList=list.files(path = "Results_01/") #here is regular expression, which I don't understand

FilePath = paste("Results_01/", FileList[1], sep ="")
Data2Write=read.csv(FilePath)
 
for(FileName in FileList[2:length(FileList)]){
  FilePath = paste("Results_01/", FileName, sep ="")
  Data=read.csv(FilePath)
  Data2Write = rbind(Data2Write, Data)
}

write.csv(Data2Write, file="Results_02/Combined_Loss_Scan.txt")
