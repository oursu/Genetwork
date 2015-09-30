
args=commandArgs(trailingOnly=TRUE)

f=args[1]
out=args[2]

print("starting")
data=read.table(f,header=TRUE)
print(data)

pdf(out)
par(las=3)
barplot(as.matrix(data),col='pink',cex.lab=0.5,cex.axis=0.5,cex.names=0.5)
dev.off()