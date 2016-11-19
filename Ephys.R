##Load necessary working packages and graphics packages
library("ggplot2")
library("reshape2")
library("gridExtra")
library("minpack.lm")
library("nls2")


############################################################################################
###################################   Functions   ##########################################
############################################################################################

################################### File Maintenance ######################################
uploadfiles<-function(path){
  # Loads all atf files in a directory into the workspace
  #
  # Args:
  #   path: path to directory containing atf files
  #
  # Returns:
  #   A list of data frames for each atf file in the directory
  setwd(path)
  data<-lapply(dir(),read.table,skip=10,header=TRUE) 
}

change_headers<-function(atf_data){
  # Changes default atf headers with simplified headers denoted T for time, V# for voltage 
  # and number of the trace, or C# for current and number of the trace
  #
  # Args:
  #   atf_data: a data frame in atf format
  #
  # Returns:
  #   The same data frame but with simplified headers
  names_col<-vector("list",length(atf_data))
  
  for(i in 1:length(atf_data)){
    if(i==1){
      names_col[[i]]<-'T'
    }else if(i%%2==0){
      names_col[[i]]<-paste('C',i-1,sep="")
    }else{
      names_col[[i]]<-paste('V',i-2,sep="")
    }
  }
  names_col
}

organize<-function(atf_data){
  # Separates current and voltage data with respect to time and returns in both long and 
  # wide format
  #
  # Args:
  #   atf_data: a data frame in atf format
  #
  # Returns:
  #   A list of 4 data frames [current,voltage,voltage(long format),current(long format)]
  
  colnames(atf_data)<-change_headers(atf_data)
  
  odd_indices<-seq(3,length(atf_data),2)
  even_indices<-seq(2,length(atf_data),2)
  voltage<-data.frame(atf_data[1],atf_data[,odd_indices])
  current<-data.frame(atf_data[1],atf_data[,even_indices])
  
  current_melt<-melt(current,id.vars="T")
  voltage_melt<-melt(voltage,id.vars="T")
  org_data<-list(voltage,current,voltage_melt,current_melt)
}

#################################### Normalization #######################################

readcap <- function()
{ 
  n <- readline(prompt="Enter the capacitance: ")
  return(as.numeric(n))
}

align_peak<-function(cell_data){
  zero_mark<-subset(cell_data,cell_data$Inorm==max(Inorm))$time
  aligned_cell_data<-cell_data
  aligned_cell_data$time<-aligned_cell_data$time-zero_mark
  return(aligned_cell_data)
}

normalize_current<-function(cell_data){
  peak<-max(cell_data$Inorm)
  normalized<-cell_data
  normalized$Inorm<-normalized$Inorm/peak
  return(normalized)
}

normalize_multiple<-function(cell_data_list){
  normpeaks<-by(cell_data_list,cell_data_list$cell,normalize)
  df<-do.call("rbind",normpeaks)
  return(df)
}

align_multiple<-function(cell_data_list){
  aligned_data<-by(cell_data_list,cell_data_list$cell,align_peak)
  df <- do.call("rbind", aligned_data)
  return(df)
}

cdensity<-function(atf_data){
  cap<-readcap()
  list_format<-organize(atf_data)
  current<-list_format[[4]]
  current$variable<-as.character(current$variable)
  act<-subset(current,T>0.5730 & T<0.8682)
  max_trace_num<-subset(act,value==max(value))$variable
  max_trace<-subset(act,variable==max_trace_num)
  itrace<-data.frame(max_trace$T,max_trace$value)
  colnames(itrace)<-c("t","i")
  itrace$t<-itrace$t*1000
  itrace$i<-itrace$i/cap
  
  testPlot<-ggplot(itrace, aes(x=t, y=i)) + 
    geom_point(colour="black",shape=1, size=1) +
    xlab('time (ms)')+
    ylab('Current Density (pA/pF)') +
    theme(legend.title=element_blank()) +
    theme(panel.margin=unit(2,"cm"),
          panel.background=element_blank(),
          axis.line.x=element_line(colour="black"),
          axis.line.y=element_line(colour="black"),
          axis.title=element_text(face="bold"),
          plot.margin=unit(c(1,1,1,1),"cm"))
  
  print(testPlot)
  return(itrace)
}

##################################### Protocols ########################################

rundown_io<-function(atf_data){
  list_format<-organize(atf_data)
    
  current_melted<-subset(list_format[[4]], T>0.01 & T<0.1)
  mean_i<-aggregate(current_melted$value,by=list(current_melted$variable),mean)

  rundown<-data.frame((seq(1,length(mean_i$x),1)-1)*5,mean_i$x)
  colnames(rundown)<-c("t","i")
    
  plt<-ggplot(rundown,aes(x=t,y=i))+
    geom_point(color="blue",shape=1,size=2) +
    ylab('I (pA)') +
    xlab("time (s)")+
    ylim(c(0,1.2*max(rundown$i)))+
    geom_segment(x=5*5,xend=12*5,y=1.1*max(rundown$i),yend=1.1*max(rundown$i), colour="orange",size=3)+
    geom_segment(x=19*5,xend=33*5,y=1.1*max(rundown$i),yend=1.1*max(rundown$i), colour="orange",size=3)+
    #geom_segment(x=58*5,xend=67*5,y=1.1*max(rundown$i),yend=1.1*max(rundown$i), colour="orange",size=3)+
    ggtitle("Rundown")+
    annotate("text",x=mean(c(5,12))*5,y=1.18*max(rundown$i),label=paste('Ani9'))+
    annotate("text",x=mean(c(19,33))*5,y=1.18*max(rundown$i),label=paste('Ani9'))+
    #annotate("text",x=mean(c(58,67))*5,y=1.18*max(rundown$i),label=paste('Ani9'))+
    theme(legend.title=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          axis.title=element_text(size=14),
          axis.line.x = element_line(colour = "black"),
          axis.line.y= element_line(colour="black"),
          plot.margin=unit(c(1,1,1,1),'cm'))
    
  print(plt)
  return(rundown)
}

rundown<-function(atf_data){
  # Extracts max current from each trace normalizes each max current to the overall max
  # The normalized max current is plotted as a function of time
  # Assumes regular intervals of 20 sec between traces when constructing plot
  # Leak current is taken as the current remaining after a 300 ms step down
  # to -80 mV after the 80 mV test step
  #
  # Args:
  #   atf_data: a data frame in atf format
  #
  # Returns:
  #   A plot of time vs normalized max current
  
  list_format<-organize(atf_data)
  
  current_melted<-subset(list_format[[4]], T>0.2740 & T<0.8659)
  leak<-subset(current_melted,T==0.8658)
  
  max_I<-aggregate(current_melted$value, by=list(current_melted$variable), max)
  max_I$x<-max_I$x+leak$value
  rundown<-data.frame((seq(1,length(max_I$x),1)-1)*5,max_I$x)
  rundown<-cbind(rundown,rep("c1",length(rundown[,1])))
  colnames(rundown)<-c("time","Inorm","cell")
  rundown<-normalize_current(rundown)
  
  plt<-ggplot(rundown,aes(x=time,y=Inorm))+
    geom_point(color="blue",shape=1,size=2) +
    #ylab(expression('I (pA)')) +
    ylab(expression('I/I'[max]))+
    #ylim(c(0,max(rundown$Inorm*1.1)))+
    scale_y_continuous(limits=c(0,1.2),breaks=seq(0,1,0.5))+
    xlab("time (s)")+
    #geom_segment(x=0,xend=33*20,y=1.1,yend=1.1, colour="purple",size=3)+
    #geom_segment(x=33*20,xend=49*20,y=1.1,yend=1.1, colour="orange",size=3)+
    #geom_segment(x=49*20,xend=max(rundown$time),y=1.1,yend=1.1, colour="purple",size=3)+
    ggtitle("Rundown")+
    #annotate("text",x=mean(c(0,33))*20,y=1.2,label=paste('control'))+
    #annotate("text",x=mean(c(33,49))*20,y=1.2,label=paste('dtt'))+
    #annotate("text",x=mean(c(49,max(rundown$time)/20))*20,y=1.2,label=paste('control'))+
    theme(legend.title=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          axis.title=element_text(size=14),
          axis.line.x = element_line(colour = "black"),
          axis.line.y= element_line(colour="black"),
          plot.margin=unit(c(1,1,1,1),'cm'))

  print(plt)
  return(rundown)
  
  
}

plot_rundown<-function(rundown_data){
  ggplot(rundown_data,aes(x=time,y=Inorm,colour=cell))+
    geom_line(size=1) +
    ylab(expression('I/I'[max]))+
    xlab("time (s)")+
    theme(legend.title=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          axis.title=element_text(size=14),
          axis.line.x = element_line(colour = "black"),
          axis.line.y= element_line(colour="black"),
          plot.margin=unit(c(1,1,1,1),'cm'))
}

tracePlot<-function(atf_data){
  
  list_format<-organize(atf_data)
  
  current<-list_format[[4]]
  voltage<-list_format[[3]]
  
  current<-subset(current,T<=0.3 & T>=0.075)
  voltage<- subset(voltage,T<=0.3 & T>=0.075)
  
  c_ymax=max(current$value)
  c_ymin=min(current$value)
  c_xmax=max(current$T)
  c_xmin=min(current$T)
  
  c_scale_barx<-round((c_xmax-c_xmin)/10,2)
  c_scale_bary<-round((c_ymax-c_ymin)/4,-3)
  
  v_ymax=max(voltage$value)
  v_ymin=min(voltage$value)
  v_xmax=max(voltage$T)
  v_xmin=min(voltage$T)
  
  v_scale_barx<-round((v_xmax-v_xmin)/10,2)
  v_scale_bary<-round((v_ymax-v_ymin)/4,-1)
  

  
  currentPlot<-ggplot(current, aes(x=T, y=value, group = variable)) + 
    geom_line() +
    geom_segment(aes(x=1.2*c_xmax-c_scale_barx,y=c_ymax-c_scale_bary,xend=1.2*c_xmax-c_scale_barx,yend=c_ymax))+
    geom_segment(aes(x=1.2*c_xmax-c_scale_barx,y=c_ymax-c_scale_bary,xend=1.2*c_xmax,yend=c_ymax-c_scale_bary))+
    annotate("text",x=1.18*c_xmax-1.4*c_scale_barx,y=c_ymax-0.5*c_scale_bary,label=paste(c_scale_bary,"pA"))+
    annotate("text",x=1.2*c_xmax-0.5*c_scale_barx,y=c_ymax-1.2*c_scale_bary,label=paste(c_scale_barx*1000,"ms"))+
    theme(legend.title=element_blank()) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text = element_blank(),
          axis.title=element_blank())

  voltagePlot<-ggplot(voltage, aes(x=T,y=value, group = variable)) +
    geom_line() +
    geom_segment(aes(x=1.2*v_xmax-v_scale_barx,y=v_ymax-v_scale_bary,xend=1.2*v_xmax-v_scale_barx,yend=v_ymax))+
    geom_segment(aes(x=1.2*v_xmax-v_scale_barx,y=v_ymax-v_scale_bary,xend=1.2*v_xmax,yend=v_ymax-v_scale_bary))+
    annotate("text",x=1.18*v_xmax-1.4*v_scale_barx,y=v_ymax-0.5*v_scale_bary,label=paste(v_scale_bary,"mV"))+
    annotate("text",x=1.2*v_xmax-0.5*v_scale_barx,y=v_ymax-1.4*v_scale_bary,label=paste(v_scale_barx*1000,"ms"))+
    theme(legend.title=element_blank()) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text = element_blank(),
          axis.title=element_blank())
  
  grid.arrange(currentPlot,voltagePlot, heights=c(0.7,0.3),ncol=1)
}

rampIV<-function(atf_data){
  list_format<-organize(atf_data)
  
  current<-list_format[[2]]
  current<-subset(current,T>=0.045 & T<=0.54)
  voltage<-list_format[[1]]
  voltage<-subset(voltage,T>=0.045 & T<=0.54)
  meani<-rowMeans(current)
  meanv<-rowMeans(voltage)
  ci95i<-apply(current,1,sd)/sqrt(length(current))
  ci95i<-1.96*ci95i
  ramp<-data.frame(meanv,meani,ci95i)
  
  plt<-ggplot(ramp,aes(x=meanv, y=meani))+
    geom_line(colour="blue",size=1)+
    ylab("pA")+
    xlab("mV")+
    ylim(c(-500,2000))+
    scale_x_continuous(limits=c(-80,80),breaks=seq(-80,80,10))+
    geom_hline(yintercept=0)+
    geom_vline(xintercept=0)+
    theme_bw()+
    theme(panel.margin=unit(2,"cm"),
          axis.title=element_text(face="bold"),
          plot.margin=unit(c(1,1,1,1),"cm"))
  print(plt)
  return(ramp)
}

IV<-function(atf_data){
  list_format<-organize(atf_data)
  
  current<-list_format[[4]]
  voltage<-list_format[[3]]
  
  ss_current<-subset(current, T==0.258)
  ss_voltage<- subset(voltage, T==0.258)
  
  IV_data<-data.frame(ss_voltage$value,ss_current$value)
  colnames(IV_data)<-c("voltage","current")
  
  IVplot<-ggplot(IV_data, aes(x=voltage, y=current)) + 
    geom_point(colour="blue",shape=1, size=3)+
    xlab("Voltage (mV)")+
    ylab("Current (pA)")+
    geom_hline(yintercept=0, colour="red", lty="dashed")+
    scale_x_continuous(limits=c(-90,90),breaks=seq(-80,80,20))+
    theme(panel.margin=unit(2,"cm"),
          panel.background=element_blank(),
          axis.line.x=element_line(colour="black"),
          axis.line.y=element_line(colour="black"),
          axis.title=element_text(face="bold"),
          plot.margin=unit(c(1,1,1,1),"cm"))
  
  print(IVplot)
}


####################################### Kinetics #########################################

act_kinetics<-function(atf_data){
  list_format<-organize(atf_data)
  current<-list_format[[4]]
  current$variable<-as.character(current$variable)
  act<-subset(current,T>0.2730 & T<0.5686)
  max_trace_num<-subset(act,value==max(value))$variable
  max_trace<-subset(act,variable==max_trace_num)
  itrace<-data.frame(max_trace$T,max_trace$value)
  colnames(itrace)<-c("t","i")
  itrace$t<-itrace$t*1000
  write.csv(itrace, paste(wdir,"/activation_kinetics.txt",sep=""),row.names=FALSE)
  
  st1<-expand.grid(a=seq(-10000,10000,1000),tau=seq(20,200,20),c=seq(-5000,5000,500))
  roughmodel<-nls2(i~a*exp(-(t)/tau)+c,data=itrace,start=st1,
                   algorithm="brute-force")
  guesses<-coef(roughmodel)
  fitmodel<- nlsLM(i ~ a*exp(-t/tau)+c,data=itrace,start=list(a=guesses[1],tau=guesses[2],c=guesses[3]))
  params<-coef(fitmodel)
  
  testPlot<-ggplot(itrace, aes(x=t, y=i)) + 
    geom_point(colour="black",shape=1, size=1) +
    geom_line(aes(x=t,y=predict(fitmodel)),linetype="dashed",colour="red",size=1)+
    annotate("text",x=1.1*min(itrace$t),y=0.9*max(itrace$i),label=paste('tau=',round(as.numeric(params[2]),3),"ms"))+
    xlab('time (ms)')+
    ylab('I (pA)') +
    theme(legend.title=element_blank()) +
    theme(panel.margin=unit(2,"cm"),
          panel.background=element_blank(),
          axis.line.x=element_line(colour="black"),
          axis.line.y=element_line(colour="black"),
          axis.title=element_text(face="bold"),
          plot.margin=unit(c(1,1,1,1),"cm"))
  
  print(testPlot)
  return(as.numeric(params[2]))
}

deact_kinetics<-function(atf_data){
  list_format<-organize(atf_data)
  current<-list_format[[4]]
  current$variable<-as.character(current$variable)
  act<-subset(current,T>0.5730 & T<0.8682)
  min_trace_num<-subset(act,value==min(value))$variable
  min_trace<-subset(act,variable==min_trace_num)
  itrace<-data.frame(min_trace$T,min_trace$value)
  colnames(itrace)<-c("t","i")
  itrace$t<-itrace$t*1000
  write.csv(itrace, paste(wdir,"/activation_kinetics.txt",sep=""),row.names=FALSE)
  
  st1<-expand.grid(a=seq(-10000,10000,1000),tau=seq(20,300,20),c=seq(-5000,5000,500))
  roughmodel<-nls2(i~a*exp(-(t)/tau)+c,data=itrace,start=st1,
                   algorithm="brute-force")
  guesses<-coef(roughmodel)
  fitmodel<- nlsLM(i ~ a*exp(-t/tau)+c,data=itrace,start=list(a=guesses[1],tau=guesses[2],c=guesses[3]),
                   control=nls.lm.control(maxiter=500))
  params<-coef(fitmodel)
  
  #testPlot<-ggplot(itrace, aes(x=t, y=i)) + 
  #geom_point(colour="black",shape=1, size=1) +
  #geom_line(aes(x=t,y=predict(fitmodel)),linetype="dashed",colour="red",size=1)+
  #annotate("text",x=1.1*min(itrace$t),y=0.9*max(itrace$i),label=paste('tau=',round(as.numeric(params[2]),3),"ms"))+
  #xlab('time (ms)')+
  #ylab('I (pA)') +
  #theme(legend.title=element_blank()) +
  #theme(panel.margin=unit(2,"cm"),
  #panel.background=element_blank(),
  #axis.line.x=element_line(colour="black"),
  #axis.line.y=element_line(colour="black"),
  #axis.title=element_text(face="bold"),
  #plot.margin=unit(c(1,1,1,1),"cm"))
  
  #print(testPlot)
  return(as.numeric(params[2]))
}

dtau_all<-function(atf_data){
  list_format<-organize(atf_data)
  current<-list_format[[4]]
  current$variable<-as.character(current$variable)
  split_current<-split(current,current$variable)
  taus<-lapply(split_current,trace_dtau)
  return(as.numeric(taus))
}

plot_dtaus<-function(atf_data){
  taus<-dtau_all(atf_data)
  tau_data<-data.frame(seq(0,length(taus)-1,1),taus)
  colnames(tau_data)<-c("t","tau")
  tau_data$t<-tau_data$t*20
  write.csv(tau_data, paste(wdir,"/deact tau data.txt",sep=""),row.names=FALSE)
  
  
  plt<-ggplot(tau_data,aes(t,tau))+
    geom_point(colour="blue",shape=1,size=3)+
    ylim(0,1.2*max(tau_data$tau))+
    geom_segment(x=0,xend=18*20,y=1.1*max(tau_data$tau),yend=1.1*max(tau_data$tau), colour="purple",size=3)+
    geom_segment(x=18*20,xend=35*20,y=1.1*max(tau_data$tau),yend=1.1*max(tau_data$tau), colour="orange",size=3)+
    geom_segment(x=35*20,xend=max(tau_data$t),y=1.1*max(tau_data$tau),yend=1.1*max(tau_data$tau), colour="purple",size=3)+
    xlab("Time (s)")+
    ylab("tau (ms)")+
    ggtitle("Deactivation Kinetics")+
    annotate("text",x=mean(c(0,18))*20,y=1.2*max(tau_data$tau),label=paste('control'))+
    annotate("text",x=mean(c(18,35))*20,y=1.2*max(tau_data$tau),label=paste('dtt'))+
    annotate("text",x=mean(c(35,max(tau_data$t)/20))*20,y=1.2*max(tau_data$tau),label=paste('control'))+
    theme(panel.margin=unit(2,"cm"),
          panel.background=element_blank(),
          axis.line.x=element_line(colour="black"),
          axis.line.y=element_line(colour="black"),
          axis.title=element_text(face="bold"),
          plot.margin=unit(c(1,1,1,1),"cm"))
  
  print(plt)
  
  return(tau_data)
}

plot_ataus<-function(atf_data){
  taus<-atau_all(atf_data)
  tau_data<-data.frame(seq(0,length(taus)-1,1),taus)
  colnames(tau_data)<-c("t","tau")
  tau_data$t<-tau_data$t*20
  write.csv(tau_data, paste(wdir,"/act tau data.txt",sep=""),row.names=FALSE)
  
  
  plt<-ggplot(tau_data,aes(t,tau))+
    geom_point(colour="blue",shape=1,size=3)+
    ylim(0,1.2*max(tau_data$tau))+
    geom_segment(x=0,xend=33*20,y=1.1*max(tau_data$tau),yend=1.1*max(tau_data$tau), colour="purple",size=3)+
    geom_segment(x=33*20,xend=49*20,y=1.1*max(tau_data$tau),yend=1.1*max(tau_data$tau), colour="orange",size=3)+
    geom_segment(x=49*20,xend=max(tau_data$t),y=1.1*max(tau_data$tau),yend=1.1*max(tau_data$tau), colour="purple",size=3)+
    xlab("Time (s)")+
    ylab("tau (ms)")+
    ggtitle("Activation Kinetics")+
    annotate("text",x=mean(c(0,33))*20,y=1.2*max(tau_data$tau),label=paste('control'))+
    annotate("text",x=mean(c(33,49))*20,y=1.2*max(tau_data$tau),label=paste('dtt'))+
    annotate("text",x=mean(c(49,max(tau_data$t)/20))*20,y=1.2*max(tau_data$tau),label=paste('control'))+
    theme(panel.margin=unit(2,"cm"),
          panel.background=element_blank(),
          axis.line.x=element_line(colour="black"),
          axis.line.y=element_line(colour="black"),
          axis.title=element_text(face="bold"),
          plot.margin=unit(c(1,1,1,1),"cm"))
  
  print(plt)
  
  return(tau_data)
}
  
trace_dtau<-function(trace){
  dact<-subset(trace,T>0.5730 & T<0.8682)
  colnames(dact)<-c("t","num","i")
  dact$t<-dact$t*1000
  
  st1<-expand.grid(a=seq(-10000,0,500),tau=seq(10,200,20),c=seq(-1000,1000,200))
  roughmodel<-nls2(i~a*exp(-(t)/tau)+c,data=dact,start=st1,
                   algorithm="brute-force")
  guesses<-coef(roughmodel)
  fitmodel<- nlsLM(i ~ a*exp(-t/tau)+c,data=dact,start=list(a=guesses[1],tau=guesses[2],c=guesses[3]),
                   control=nls.lm.control(maxiter=500))
  params<-coef(fitmodel)
  
  print(params[2])
  testPlot<-ggplot(dact, aes(x=t, y=i)) + 
    geom_point(colour="black",shape=1, size=1) +
    geom_line(aes(x=t,y=predict(fitmodel)),linetype="dashed",colour="red",size=1)+
    annotate("text",x=1.1*min(dact$t),y=0.9*max(dact$i),label=paste('tau=',round(as.numeric(params[2]),3),"ms"))+
    xlab('time (ms)')+
    ylab('I (pA)') +
    theme(legend.title=element_blank()) +
    theme(panel.margin=unit(2,"cm"),
          panel.background=element_blank(),
          axis.line.x=element_line(colour="black"),
          axis.line.y=element_line(colour="black"),
          axis.title=element_text(face="bold"),
          plot.margin=unit(c(1,1,1,1),"cm"))
  
  print(testPlot)
  return(as.numeric(params[2]))
}

atau_all<-function(atf_data){
  list_format<-organize(data)
  current<-list_format[[4]]
  current$variable<-as.character(current$variable)
  current<-subset(current,variable != 'C1')
  split_current<-split(current,current$variable)
  taus<-lapply(split_current,trace_atau)
  return(as.numeric(taus))
  
}

trace_atau<-function(trace){
  act<-subset(trace,T>0.2730 & T<0.5686)
  colnames(act)<-c("t","num","i")
  act$t<-act$t*1000
  
  st1<-expand.grid(a=seq(-3000,0,300),tau=seq(20,300,20),c=seq(0,10000,1000))
  roughmodel<-nls2(i~a*exp(-(t)/tau)+c,data=act,start=st1,
                   algorithm="brute-force")
  guesses<-coef(roughmodel)
  fitmodel<- nlsLM(i ~ a*exp(-t/tau)+c,data=act,start=list(a=guesses[1],tau=guesses[2],c=guesses[3]))
  params<-coef(fitmodel)
  
  print(params[2])
  testPlot<-ggplot(act, aes(x=t, y=i)) + 
    geom_point(colour="black",shape=1, size=1) +
    geom_line(aes(x=t,y=predict(fitmodel)),linetype="dashed",colour="red",size=1)+
    annotate("text",x=1.1*min(act$t),y=0.9*max(act$i),label=paste('tau=',round(as.numeric(params[2]),3),"ms"))+
    xlab('time (ms)')+
    ylab('I (pA)') +
    theme(legend.title=element_blank()) +
    theme(panel.margin=unit(2,"cm"),
          panel.background=element_blank(),
          axis.line.x=element_line(colour="black"),
          axis.line.y=element_line(colour="black"),
          axis.title=element_text(face="bold"),
          plot.margin=unit(c(1,1,1,1),"cm"))
  
  print(testPlot)
  return(as.numeric(params[2]))
}

plot_tau<-function(ctrl_tau,cond_tau){
  ctrl_dtau<-data.frame(ctrl_dtau,rep('ctrl',length(ctrl_dtau)))
  colnames(ctrl_dtau)<-c('tau','cond')
  dtt_dtau<-data.frame(dtt_dtau,rep('dtt',length(ddt_dtau)))
  colnames(dtt_dtau)<-c('tau','cond')
  dact_kin<-rbind(ctrl_dtau,dtt_dtau)
  
  dtt_mean<- mean(dact_kin$tau[dact_kin$cond=='dtt'])
  ctrl_mean<- mean(dact_kin$tau[dact_kin$cond=='ctrl'])
  dtt_sem<-sd(dact_kin$tau[dact_kin$cond=='dtt'])/sqrt(length(dact_kin$cond[dact_kin$cond=='dtt']))
  ctrl_sem<-sd(dact_kin$tau[dact_kin$cond=='ctrl'])/sqrt(length(dact_kin$cond[dact_kin$cond=='ctrl']))
  
  mann_whit<-wilcox.test(tau~cond,data=dact_kin)
  
  agg_plot<-ggplot(dact_kin,aes(x=tau,y=cond,group=cond,colour=cond))+
    geom_point(size=3,shape=1)+
    scale_colour_manual(breaks=c('ctrl','dtt'),values=c('red','blue'))+
    geom_point(aes(x=ctrl_mean,y=cond[1]),shape=18,size=6,colour="red")+
    ggtitle("Activation Kinetics")+
    geom_point(aes(x=dtt_mean,y=cond[9]),shape=18,size=6,colour="blue")+
    geom_segment(aes(x=dtt_mean-dtt_sem,xend=dtt_mean+dtt_sem,y=cond[9],yend=cond[9]),colour='blue',size=1)+
    geom_segment(aes(x=ctrl_mean-ctrl_sem,xend=ctrl_mean+ctrl_sem,y=cond[1],yend=cond[1]),colour='red',size=1)+
    ylab("Condition")+
    xlab("tau (ms)")+
    theme_bw()+
    theme(panel.margin=unit(2,"cm"),
          legend.position="none",
          axis.title=element_text(face="bold"),
          plot.margin=unit(c(1,1,1,1),"cm"))
  
  print(agg_plot)
  return(mann_whit)
}

peak90<-function(cell_data){
  peak<-max(cell_data$Inorm)
  peaktime<-subset(cell_data,cell_data$Inorm == peak)
  sub_peak<-subset(cell_data,Inorm<0.80*peak & time>=peaktime$time)
  return(min(sub_peak$time)-peaktime$time)
}



############################################################################################
############################################################################################
############################################################################################

wdir<-"C:/Users/mny3/Desktop/Yang Lab/Data/Patching/TMEM16B/11-7-16 TMEM16B-HEK I-O"
filename<-"16n07016.txt"

setwd(wdir)
data<-read.table(filename,skip=10,header=TRUE)

trace8<-cbind(trace8,rep("c8",length(trace8[1])))



ggplot(sum, aes(x=t, y=i,colour=cond)) + 
  geom_point(size=1) +
  scale_colour_manual(breaks=c('ctrl','dtt'),values=c('red','blue'))+
  xlab('time (ms)')+
  ylab('pA/pF') +
  theme(legend.title=element_blank()) +
  theme(panel.margin=unit(2,"cm"),
        panel.background=element_blank(),
        axis.line.x=element_line(colour="black"),
        axis.line.y=element_line(colour="black"),
        axis.title=element_text(face="bold"),
        plot.margin=unit(c(1,1,1,1),"cm"))

