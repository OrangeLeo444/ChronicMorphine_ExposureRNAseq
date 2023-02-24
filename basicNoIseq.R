library(NOISeq)
myDataMyFactors<-readData(data=cts_clean, factors=cData_filtered)
myDataMyFactors

myTMM <- tmm(assayData(myDataMyFactors)$exprs,long=1000,lc=0)

myfilt <- filtered.data(cts_clean, factor = myDataMyFactors$group, norm=FALSE, depth=NULL, method=1, cv.cutoff=100, cpm=1, p.adj="fdr")


myDataMyFactors<-readData(data=myfilt, factors=cData_filtered)
myTMM <- tmm(assayData(myDataMyFactors)$exprs,long=1000,lc=0)

myDataMyFactors<-readData(data=myTMM, factors=cData_filtered)

mynoiseqSim <- noiseq(myDataMyFactors, k = 0.5, norm = "n", factor = "group", pnr = 0.2, nss = 5, v = 0.02)
mynoiseqSim.deg <- degenes(mynoiseqSim, q=0.9, M=NULL)




#Plot expression with DEGs highlighted in red
DE.plot(mynoiseqSim, q = 0.9, graphic = "expr", log.scale = TRUE)
DE.plot(mynoiseqSim, q = 0.9, graphic = "MD", log.scale = TRUE)
mynoiseqSim.deg1 = degenes(mynoiseqSim, q = 0.9, M = "up")
mynoiseqSim.deg2 = degenes(mynoiseqSim, q = 0.9, M = "down")
