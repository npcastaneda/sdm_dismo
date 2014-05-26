dismo.mod <- function(sp,registros=reg.filt,predictors=MyExpl,MaxEnt=T,Bioclim=T,Mahal=F,part=3){
  library(dismo)
  #library(XML)
  library(raster)
  library(rgdal)
  print(date())
  
  cat(paste("Modelando",sp,"...",'\n'))
  #extrai as coordenadas de cada especie 
  coord <- registros[registros[,1]==sp,2:3]
  #tabela de valores
  #presvals <- extract(predictors, coord)
  backgr <- randomPoints(predictors, 500)
  
  group <- kfold(coord, part)
  
  for (i in unique(group)){
    cat(paste("Partição número",i,'\n'))
    pres_train <- coord[group != i, ]
    pres_test <- coord[group == i, ]
    
    #faz os modelos
    cat(paste("Fazendo os modelos...",sp,i,'\n'))
    if (Bioclim==T){
      cat(paste("Bioclim",'\n'))
      bc <- bioclim (predictors, pres_train, progress="text")
    }
    if (MaxEnt==T){ 
      cat(paste("Maxent",'\n'))
      mx <- maxent (predictors, pres_train, progress="text")
    }
    if (Mahal==T){ 
      cat(paste("Mahalanobis",'\n'))
      ma <- mahal  (predictors, pres_train, progress="text")
    }
    
    #avaliação
    cat(paste("Avaliando os modelos...",sp,i,'\n'))
    if(exists("bc")) {
      ebc <- evaluate(pres_test, backgr, bc, predictors)
      thresholdbc <- ebc@t[which.max(ebc@TPR + ebc@TNR)]
    } 
    if(exists("mx")){
      emx <- evaluate(pres_test,backgr,mx,predictors)
      thresholdmx <- emx@t[which.max(emx@TPR + emx@TNR)]
    }
    
    if(exists("ma")){
      ema <- evaluate(pres_test,backgr,ma,predictors)
      thresholdma <- ema@t[which.max(ema@TPR + ema@TNR)]
    }
    
    #prediz para a área toda
    cat(paste("Projetando os modelos...",sp,i,'\n'))
    if(exists("bc")) 
      contbc <- predict(predictors, bc, progress='text',filename=paste('./models/',sp,i,'bioclim_prediction.grd'),overwrite=T)
    if(exists("mx")) 
      contmx <- predict(mx, predictors,progress='text',filename=paste('./models/',sp,i,'maxent_prediction.grd'),overwrite=T) 
    if(exists("ma")) 
      contma <- predict(ma, predictors,progress='text',filename=paste('./models/',sp,i,'mahalanobis_prediction.grd'),overwrite=T) 
    
    #cria os modelos binarios
    cat(paste("Criando os modelos binários...",sp,i,'\n'))
    if(exists("contbc")) 
      binbc <- contbc > thresholdbc
    if(exists("contmx")) 
      binmx <- contmx > thresholdmx
    if(exists("contma")) 
      binma <- contma > thresholdma
    
    #cria o modelo final
    cat(paste("Criando os modelos cortados...",sp,i,'\n'))
    if(exists("contbc")) 
      bcfinal <- contbc * binbc
    if(exists("contmx"))
      mxfinal <- contmx * binmx
    if(exists("contma"))
      mafinal <- contma * binma
    
    #escreve continuos
    cat(paste("Escrevendo os rasters...",sp,i,'\n'))
    if(exists("contbc"))
      writeRaster(x=contbc,filename=paste0("./models/",sp,i,"bioclim_cont.grd"),progress="text",overwrite=T)
    if(exists("contmx"))
      writeRaster(x=contmx,filename=paste0("./models/",sp,i,"maxent_cont.grd"),progress="text",overwrite=T)
    if(exists("contma")) 
      writeRaster(x=contma,filename=paste0("./models/",sp,i,"mahal_cont.grd"),progress="text",overwrite=T)
    
    #escreve binarios
    if(exists("binbc"))
      writeRaster(x=binbc,filename=paste0("./models/",sp,i,"bioclim_binario.grd"),progress="text",overwrite=T)
    if(exists("binmx"))
      writeRaster(x=binmx,filename=paste0("./models/",sp,i,"maxent_binario.grd"),progress="text",overwrite=T)
    if(exists("binma")) 
      writeRaster(x=binma,filename=paste0("./models/",sp,i,"mahal_binario.grd"),progress="text",overwrite=T)
    
    #escreve finais
    if(exists("bcfinal"))
      writeRaster(x=bcfinal,filename=paste0("./models/",sp,i,"bioclim_cut.grd"),overwrite=T)
    if(exists("mxfinal"))
      writeRaster(x=mxfinal,filename=paste0("./models/",sp,i,"maxent_cut.grd"),overwrite=T)
    if(exists("mafinal"))    
      writeRaster(x=mafinal,filename=paste0("./models/",sp,i,"mahal_cut.grd"),overwrite=T)
    
    cat(paste("Gerando o arquivo de avaliação...",sp,i,'\n'))
    sink(file=paste0("./models/evaluate_",sp,i,".txt"),split=T)
    print(paste(sp,i))
    if(exists("ebc")){
      print("BIOCLIM")
      print(ebc)
    }
    if(exists("emx")){
      print("MAXENT")
      print(emx)
    }
    if(exists("ema")){
      print("MAHALANOBIS")
      ema
    }
    sink()   
    
  }
  
  
  #totalconsensus
  #cat(paste("Gerando os modelos de consenso absoluto...",sp,i,'\n'))  
  #if(exists(c("binbc","binmx"))){
  # totalconsensus <- binbc * binmx * binma
  #writeRaster(x=totalconsensus,filename=paste0("./models/totalconsensus",sp,".grd"),overwrite=T)
  #}
  cat("====FIN====",'\n')
  print(date())
}

