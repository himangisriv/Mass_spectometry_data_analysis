library(NPARC)
library(tidyverse)
library(readxl)
library(tidyverse)
library(ggplot2)

#upload te data 
my_data <- read_excel("2022.03.23-2.xlsx")
my_data <- na.omit(my_data) 
my_data<-my_data %>% select("protein id","37C (TMT-128C) m128.134436","40C (TMT-127N) m127.124761","43C (TMT-126) m126.127726","46C 127C (TMT-127C) m127.131081",
                            "49C (TMT-131) m131.13818","52C (TMT-130C) m130.141145","55C (TMT-130N) m130.134825","58C (TMT-129C) m129.13779","61C (TMT-129N) m129.131471",
                            "64C (TMT-128N) m128.128116")

my_data <-my_data %>% group_by(my_data$`protein id`)%>%
  summarise_at(vars("37C (TMT-128C) m128.134436","40C (TMT-127N) m127.124761","43C (TMT-126) m126.127726","46C 127C (TMT-127C) m127.131081",
                    "49C (TMT-131) m131.13818","52C (TMT-130C) m130.141145","55C (TMT-130N) m130.134825","58C (TMT-129C) m129.13779","61C (TMT-129N) m129.131471",
                    "64C (TMT-128N) m128.128116"), mean)

cols <- c("protein",seq(37,64,3))
colnames(my_data) <- cols
df<-my_data
df <- as.data.frame(t(apply(my_data[,-1], 1, function(x)(x/max(x)))))
df <- cbind(my_data[,1], df)
df<-na.omit(df)
mydata2<-df %>% gather(., "Temp", "Relative abundance", 2:11)
mydata2$Temp<-as.numeric(mydata2$Temp)
protein=unique(mydata2$protein)

#function to check the curve
EXAMINE_NULL<-function(data,protein_name)
{
  stk42<- filter(data, protein == protein_name)
  # Fit the data
  nullFit <- NPARC:::fitSingleSigmoid(x = stk42$Temp, y = stk42$`Relative abundance`)
  if(is(nullFit)!="try-error")
  {
    #retrieve the parameters
    Pl=as.numeric(strsplit(capture.output(nullFit)[5], " +")[[1]][2])
    a=as.numeric(strsplit(capture.output(nullFit)[5], " +")[[1]][3])
    b=as.numeric(strsplit(capture.output(nullFit)[5], " +")[[1]][4])
    RSS.p <- sum(residuals(nullFit)^2)
    TSS <- sum((stk42$`Relative abundance` - mean(stk42$`Relative abundance`))^2)
    #findind the R2 value for the data fitted 
    R2<-1 - (RSS.p/TSS)
    val<-c(protein_name,as.numeric(Pl),as.numeric(a),as.numeric(b),as.numeric(R2))
  }
  else
  {
    Pl<-NA
    a<-NA
    b<-NA
    RSS.p <- NA
    TSS <- NA
    R2<-NA
    print(paste0("convergence not possible at ", protein_name, ".") )
    val<-c(protein_name,Pl,a,b,R2)
  }
  
  return(val)
}

#function to plot the graph 
plot_ms_data<-function(value,data,protein_name)
{
  stk42<- filter(data, protein == protein_name)
  if(is.na(value[2])==FALSE & is.na(value[3])==FALSE & 
     is.na(value[4])==FALSE & is.na(value[5])==FALSE)
  {
    if(value[5]>=0.9)
    {
      Pl=as.numeric(value[2])
      a=as.numeric(value[3])
      b=as.numeric(value[4])
      #png(file = paste0("/Users/himangisrivastava/Desktop/spectometry_data/",paste(protein_name, '.png', sep = '')))
      stk4_plot_orig<-ggplot(data = stk42, mapping=aes(x = stk42$Temp, y = stk42$`Relative abundance`)) +
        geom_point()+
        theme_bw() +
        stat_function(fun = function(x) ((1 - Pl)/(1 + exp((b - a/x))) + Pl ))+
        ggtitle(protein,value[5])
      print(stk4_plot_orig)
      #dev.off()
    }
  }
  
}

#function to retrieve the molten temperature
retrieve_molten_temp<-function(data,protein_name)
  {
  stk42<- filter(data, protein == protein_name)
  nullFit <- NPARC:::fitSingleSigmoid(x = stk42$Temp, y = stk42$`Relative abundance`)
  if(is(nullFit)!="try-error")
  {
  nullPredictions <- broom::augment(nullFit)
  dat<-nullPredictions %>% filter(.fitted <= 0.60 & .fitted >=0.50)%>% select(x)
  }
  else{
    print("convergence error")
  }
  return(dat)
}
protein_name = protein[8]
data_f = data.frame(matrix(nrow=length(protein),ncol=5))
for (i in 1: length(protein))
{
  protein_name=protein[i]
  
  moten_temp<-retrieve_molten_tem(mydata2,protein_name)$x
  data_f[i, ] <- value
  print(plot_ms_data(value,mydata2,protein_name))
  
}
