PANEL_supplemental_figure1=function(plot_all=FALSE){

	library(xlsx)
	options(stringsAsFactors=FALSE)

	print("INFO|supplemental_figure1|init")
	
	if(plot_all){ 
		pdf("PANEL_Figure2NEW_sym_immun_variables.pdf",height=15,width=15) 
	}else{
		pdf("PANEL_Figure2NEW_sym_immun_variables.pdf",height=8,width=8) }
	if(plot_all){ 
		par(mfrow=c(3,4)) 
	}else{	
		par(mfrow=c(3,4)) 
	}
	# ------------------------------------------------------------------------------
	table_file=XLSX_FILE
	# ------------------------------------------------------------------------------
	data=read.xlsx(table_file,encoding = "UTF-8",3)
	df=as.data.frame(data); dd=read.xlsx(table_file,encoding = "UTF-8",4)

	#
	# Immune variables
	#	
	colnames(dd)=gsub("\\."," ",colnames(dd))

	# 
	#	Combining of the values from the free light chains
	#	for making an extra column nasal.FLC
	#
	nasal.FLC=c()
	for(x in 1:dim(dd)[1]){
		arr=c(dd[x,"FLC  U 03BB "],dd[x,"FLC  U 03BA "])
		arr=arr[!is.na(arr)]
		if(length(arr)>0){ nasal.FLC=c(nasal.FLC,sum(arr))
		}else{ nasal.FLC=c(nasal.FLC,NA) }
	}

	#
	# now we add the FLC column into the dataframe
	#
	dd=as.data.frame(dd); dd=cbind(dd,nasal.FLC)

	#
	# vector containing the patients
	#
	X=dd[,2]; X=unlist(strsplit(X,"_")); X=unique(X[seq(1,2*dim(dd)[1],2)])
	cn=colnames(dd)
	for(my_cn in 3:length(cn)){
		
		print(cn[my_cn])
		
		b_in_sym=c(); 
		b_in_immune=c(); 
		b_out_sym=c(); 
		b_out_immune=c(); 
		
		na_in_sym=c(); 
		na_in_immune=c(); 
		na_out_sym=c(); 
		na_out_immune=c()
		
		c11=list(); 
		c11$estimate=NA; 
		c11$p.value=NA;
		
		c21=list(); 
		c21$estimate=NA; 
		c21$p.value=NA;
		
		c31=list(); 
		c31$estimate=NA; 
		c31$p.value=NA;
		
		c41=list(); 
		c41$estimate=NA; 
		c41$p.value=NA;
		
		for(p in 1:length(X)){
			df_filt=subset(df,df[,1]==tolower(X[p]))
			my_sym_in=df_filt[,"Symptom.birch..avg..in.season."]
			my_sym_out=df_filt[,"Symptom.birch..avg..out.of.season."]
			
			u6=paste(X[p],"_6",sep=""); 
			u7=paste(X[p],"_7",sep=""); 
			u8=paste(X[p],"_8",sep="")
			
			ix6=which(dd[,2]==u6); 
			ix7=which(dd[,2]==u7); 
			ix8=which(dd[,2]==u8); 
			ix_b_in=(c(ix6,ix7,ix8))
			
			u1=paste(X[p],"_1",sep=""); 
			u2=paste(X[p],"_2",sep=""); 
			u13=paste(X[p],"_13",sep=""); 
			u14=paste(X[p],"_14",sep=""); 
			u15=paste(X[p],"_15",sep="")
			
			ix1=which(dd[,2]==u1); 
			ix2=which(dd[,2]==u2); 
			ix13=which(dd[,2]==u13); 
			ix14=which(dd[,2]==u14); 
			ix15=which(dd[,2]==u15); 
			ix_b_out=(c(ix1,ix2,ix13,ix14,ix15))
			
			mu1=(dd[ix_b_in,my_cn]); 
			mu1=(as.numeric(mu1)); 
			mu1=mu1[!is.na(mu1)]; 
			mu1=mean(mu1);
			
			if(length(grep("^b",X[p]))){ 
				b_in_immune=c(b_in_immune,mu1); 
				b_in_sym=c(b_in_sym,my_sym_in); 
			}else if(length(grep("^na",X[p]))){ 
				na_in_immune=c(na_in_immune,mu1); 
				na_in_sym=c(na_in_sym,my_sym_in); 
			}
			
			mu2=(dd[ix_b_out,my_cn]); 
			mu2=(as.numeric(mu2)); 
			mu2=mu2[!is.na(mu2)]; 
			mu2=mean(mu2)
			
			if(length(grep("^b",X[p]))){
				b_out_immune=c(b_out_immune,mu2); 
				b_out_sym=c(b_out_sym,my_sym_out) 
			}else if(length(grep("^na",X[p]))){ 
				na_out_immune=c(na_out_immune,mu2); 
				na_out_sym=c(na_out_sym,my_sym_out); 
			}
		}
		c_x_max=dd[,my_cn]; 
		c_x_max=c_x_max[!is.na(c_x_max)]; 
		c_x_max=max(c_x_max);
		
		c_y_max=c(df[,"Symptom.birch..avg..in.season."],df[,"Symptom.birch..avg..out.of.season."]); 
		c_y_max=c_y_max[!is.na(c_y_max)]; 
		c_y_max=max(c_y_max)
		
		if(plot_all){
			plot(na_in_immune,na_in_sym,col="blue",main=cn[my_cn],pch=20,xlab="abundance [in season]",ylab="symptom [in season]", xlim=c(0,c_x_max), ylim=c(0,c_y_max)) 
		}else{
			plot(na_in_immune,na_in_sym,col="blue",main=cn[my_cn],pch=20,xlab="abundance [in season]",ylab="symptom [in season]") 
		}
		
		sdf=data.frame(na_in_immune,na_in_sym); sdf=na.omit(sdf);
		if(sd(sdf[,1])>0){
			abline(lm(sdf[,2]~sdf[,1]),col="blue",main=cn[my_cn])
			c11=cor.test(sdf[,1],sdf[,2],method="pearson"); c12=cor.test(sdf[,1],sdf[,2],method="spearman")
			if(!(plot_all)){
				legend("topright",legend=substitute(paste(italic(r)[spearman], "=", r1, ", ", italic(p)[spearman], "=", p1),list(r1=round(c11$estimate,2),p1=round(c11$p.value,2))),bty="n")
				legend("bottomright",legend=substitute(paste(italic(r)[pearson], "=", r2, ", ", italic(p)[pearson], "=", p2),list(r2=round(c12$estimate,2),p2=round(c12$p.value,2))),bty="n")
			}
		}
		
		if(plot_all){ points(na_out_immune,na_out_sym,col="#0000ff55",main=cn[my_cn],pch=20,xlab="abundance [out of season]",ylab="symptom [out of season]") }
		else{ plot(na_out_immune,na_out_sym,col="#0000ff55",main=cn[my_cn],pch=20,xlab="abundance [out of season]",ylab="symptom [out of season]") }
		sdf=data.frame(na_out_immune,na_out_sym); sdf=na.omit(sdf);
		if(sd(sdf[,1])>0){
			abline(lm(sdf[,2]~sdf[,1]),col="blue",main=cn[my_cn])
			c21=cor.test(sdf[,1],sdf[,2],method="pearson"); c22=cor.test(sdf[,1],sdf[,2],method="spearman")
			if(!(plot_all)){
				legend("topright",legend=substitute(paste(italic(r)[spearman], "=", r1, ", ", italic(p)[spearman], "=", p1),list(r1=round(c21$estimate,2),p1=round(c21$p.value,2))),bty="n")
				legend("bottomright",legend=substitute(paste(italic(r)[pearson], "=", r2, ", ", italic(p)[pearson], "=", p2),list(r2=round(c22$estimate,2),p2=round(c22$p.value,2))),bty="n")
			}
		}
		
		if(plot_all){ points(b_in_immune,b_in_sym,col="red",main=cn[my_cn],pch=20,xlab="abundance [in season]",ylab="symptom [in season]") }
		else{ plot(b_in_immune,b_in_sym,col="red",main=cn[my_cn],pch=20,xlab="abundance [in season]",ylab="symptom [in season]") }
		sdf=data.frame(b_in_immune,b_in_sym); sdf=na.omit(sdf);
		if(sd(sdf[,1])>0){
			abline(lm(sdf[,2]~sdf[,1]),col="red",main=cn[my_cn])
			c31=cor.test(sdf[,1],sdf[,2],method="pearson"); c32=cor.test(sdf[,1],sdf[,2],method="spearman")
			if(!(plot_all)){
				legend("topright",legend=substitute(paste(italic(r)[spearman], "=", r1, ", ", italic(p)[spearman], "=", p1),list(r1=round(c31$estimate,2),p1=round(c31$p.value,2))),bty="n")
				legend("bottomright",legend=substitute(paste(italic(r)[pearson], "=", r2, ", ", italic(p)[pearson], "=", p2),list(r2=round(c32$estimate,2),p2=round(c32$p.value,2))),bty="n")
			}
		}
		
		if(plot_all){ points(b_out_immune,b_out_sym,col="#ff000055",main=cn[my_cn],pch=20,xlab="abundance [out of season]",ylab="symptom [out of season]") }
		else{ plot(b_out_immune,b_out_sym,col="#ff000055",main=cn[my_cn],pch=20,xlab="abundance [out of season]",ylab="symptom [out of season]") }
		sdf=data.frame(b_out_immune,b_out_sym); sdf=na.omit(sdf);
		if(sd(sdf[,1])>0){
			abline(lm(sdf[,2]~sdf[,1]),col="#ff000055",main=cn[my_cn])
			c41=cor.test(sdf[,1],sdf[,2],method="pearson"); c42=cor.test(sdf[,1],sdf[,2],method="spearman")
			if(!(plot_all)){
				legend("topright",legend=substitute(paste(italic(r)[spearman], "=", r1, ", ", italic(p)[spearman], "=", p1),list(r1=round(c41$estimate,2),p1=round(c41$p.value,2))),bty="n")
				legend("bottomright",legend=substitute(paste(italic(r)[pearson], "=", r2, ", ", italic(p)[pearson], "=", p2),list(r2=round(c42$estimate,2),p2=round(c42$p.value,2))),bty="n")
			}else{
				legend("topright",legend=c(paste("r1=",round(c11$estimate,2),"p1=",round(c11$p.value,2)),paste("r2=",round(c21$estimate,2),"p2=",round(c21$p.value,2)),paste("r3=",round(c31$estimate,2),"p3=",round(c31$p.value,2)),paste("r4=",round(c41$estimate,2),"p4=",round(c41$p.value,2))),fill=c("blue","#0000ff55","red","#ff000055"))
			}
		}
	}
	dev.off()
}
