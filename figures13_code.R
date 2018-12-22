library(xlsx)

options(stringsAsFactors=FALSE)

#XLSX_FILE="../../PANEL_STUDY_REPORT/Tables_25apr18.xlsx"
XLSX_FILE="20180731_Panelstudy_Tables.xlsx"

#
#	PANEL_figures3D_ZOOM
#
PANEL_figures3D_ZOOM=function(show_z_scale){
	options(stringsAsFactors=FALSE)
	svg(paste("PANEL_figure3D_ZOOM_",show_z_scale,".svg",sep=""),width=20,height=10)

	library(ggplot2);
	library(grid);
	library(xlsx);
	library(gridExtra);
	library(ggpubr)

	#grid.newpage()
	if(show_z_scale==TRUE){
		data=read.xlsx(XLSX_FILE,encoding = "UTF-8",5)
		for(x in 3:dim(data)[2]){ data[,x]=as.numeric(data[,x]) }
	}else{
		data=read.xlsx(XLSX_FILE,encoding = "UTF-8",4) 
	}

	print(colnames(data))

	my_cols=c("sIgA","sIgG4","sIgG4","IgA","FLC.Î»","FLC.Îº","IL.33","MIP.1Î²","Eotaxin.2","Eotaxin.2","Eotaxin.2","Eotaxin.2","IL.1Î²","MCP.1","MCP.1")
	my_index=c(8,6,11,2,8,8,9,8,1,5,8,15,8,8,13)

	pushViewport(viewport(layout = grid.layout(nrow = length(my_cols)/2, ncol = 2)))
	define_region <- function(row, col){ viewport(layout.pos.row = row, layout.pos.col = col) }
	my_plots=list()
	for(x in 1:length(my_cols)){
		to_sel=paste("_",my_index[x],sep="")
		my_pat=paste("b.*",to_sel,"$",sep="")
		ix1=(grep(my_pat,as.character(data[,2])))
		A=data[ix1,my_cols[x]]; 
		A=A[!is.na(A)];
		print(ix1)

		my_pat=paste("na.*",to_sel,"$",sep=""); 
		ix2=(grep(my_pat,as.character(data[,2])))
		B=data[ix2,my_cols[x]]; 
		B=B[!is.na(B)]; 

		B=as.numeric(B); 
		A=as.numeric(A)

		if(length(A)>0 & length(B)>0){
			col1=c(rep("allergic",length(A)),rep("non-allergic",length(B)))
			col2=c(A,B)
			col3=c(rep("allergic",length(A)),rep("non-allergic",length(B)))
			ddf<-data.frame(A=col1,B=col2,C=col3)
			ddf[,1]=factor(ddf[,1]);ddf[,1]=factor(ddf[,1],levels(ddf[,1])[c(2,1)])
			ddf[,3]=factor(ddf[,3]);ddf[,3]=factor(ddf[,3],levels(ddf[,3])[c(2,1)])
			my_p_val=wilcox.test(ddf[,2]~ddf[,1])$p.value
			if(show_z_scale==TRUE){
				my_p=ggplot(ddf, aes(x=A, y=B, fill=C)) +
				geom_boxplot(outlier.shape=NA)+geom_jitter(position=position_jitter(width=.1, height=0))+ggtitle(paste(my_cols[x],my_index[x],my_p_val))+xlab("")+stat_summary(fun.y=mean, geom="point", shape=20, size=10, color="red", fill="red")+ylab("abundance")
			}else{
				my_p=ggplot(ddf, aes(x=A, y=B, fill=C))+geom_boxplot(outlier.shape=NA)+geom_jitter(position=position_jitter(width=.1, height=0))+ggtitle(paste(my_cols[x],my_index[x],my_p_val))+stat_summary(fun.y=mean, geom="point", shape=20, size=10, color="red", fill="red")+xlab("")+ylab("abundance")
			}
			if(x %% 2 == 0){ my_col=1; my_row=(x+1)/2 }else{ my_col=2; my_row=(x+1)/2 }
			my_plots[[x]]=my_p
		}else{
			print(c(my_cols[x],"ignored"))
		}
	}
	print(length(my_cols))
	plot(ggarrange(plotlist=my_plots,ncol=length(my_cols)/3,nrow=3))
	dev.off()
}

#
#	PANEL_supplemental_figure1
#
PANEL_supplemental_figure1=function(plot_all=FALSE){
	library(xlsx)
	
	if(plot_all){ pdf("PANEL_Figure2NEW_sym_immun_variables.pdf",height=15,width=15) }
	else{	pdf("PANEL_Figure2NEW_sym_immun_variables.pdf",height=8,width=8) }
	if(plot_all){ par(mfrow=c(3,4)) }
	else{	par(mfrow=c(3,4)) }
	# ------------------------------------------------------------------------------
	table_file=XLSX_FILE
	# ------------------------------------------------------------------------------
	options(stringsAsFactors=FALSE)
	data=read.xlsx(table_file,encoding = "UTF-8",3)
	df=as.data.frame(data); dd=read.xlsx(table_file,encoding = "UTF-8",4)
	colnames(dd)=gsub("\\."," ",colnames(dd))
	print(colnames(dd))

	nasal.FLC=c()
	for(x in 1:dim(dd)[1]){
		arr=c(dd[x,"FLC Î»"],dd[x,"FLC Îº"])
		arr=arr[!is.na(arr)]
		if(length(arr)>0){ nasal.FLC=c(nasal.FLC,sum(arr))
		}else{ nasal.FLC=c(nasal.FLC,NA) }
	}
	dd=as.data.frame(dd); dd=cbind(dd,nasal.FLC)
	X=dd[,2]; X=unlist(strsplit(X,"_")); X=unique(X[seq(1,2*dim(dd)[1],2)])
	cn=colnames(dd)
	for(my_cn in 3:length(cn)){
		b_in_sym=c(); b_in_immune=c(); b_out_sym=c(); b_out_immune=c(); na_in_sym=c(); na_in_immune=c(); na_out_sym=c(); na_out_immune=c()
		c11=list(); c11$estimate=NA; c11$p.value=NA;
		c21=list(); c21$estimate=NA; c21$p.value=NA;
		c31=list(); c31$estimate=NA; c31$p.value=NA;
		c41=list(); c41$estimate=NA; c41$p.value=NA;
		for(p in 1:length(X)){
			df_filt=subset(df,df[,1]==tolower(X[p]))
			my_sym_in=df_filt[,"Symptom.birch..avg..in.season."]
			my_sym_out=df_filt[,"Symptom.birch..avg..out.of.season."]
			u6=paste(X[p],"_6",sep=""); u7=paste(X[p],"_7",sep=""); u8=paste(X[p],"_8",sep="")
			ix6=which(dd[,2]==u6); ix7=which(dd[,2]==u7); ix8=which(dd[,2]==u8); ix_b_in=(c(ix6,ix7,ix8))
			u1=paste(X[p],"_1",sep=""); u2=paste(X[p],"_2",sep=""); u13=paste(X[p],"_13",sep=""); u14=paste(X[p],"_14",sep=""); u15=paste(X[p],"_15",sep="")
			ix1=which(dd[,2]==u1); ix2=which(dd[,2]==u2); ix13=which(dd[,2]==u13); ix14=which(dd[,2]==u14); ix15=which(dd[,2]==u15); ix_b_out=(c(ix1,ix2,ix13,ix14,ix15))
			mu1=(dd[ix_b_in,my_cn]); mu1=(as.numeric(mu1)); mu1=mu1[!is.na(mu1)]; mu1=mean(mu1);
			if(length(grep("^b",X[p]))){ b_in_immune=c(b_in_immune,mu1); b_in_sym=c(b_in_sym,my_sym_in); }else if(length(grep("^na",X[p]))){ na_in_immune=c(na_in_immune,mu1); na_in_sym=c(na_in_sym,my_sym_in); }
			mu2=(dd[ix_b_out,my_cn]); mu2=(as.numeric(mu2)); mu2=mu2[!is.na(mu2)]; mu2=mean(mu2)
			if(length(grep("^b",X[p]))){ b_out_immune=c(b_out_immune,mu2); b_out_sym=c(b_out_sym,my_sym_out) }else if(length(grep("^na",X[p]))){ na_out_immune=c(na_out_immune,mu2); na_out_sym=c(na_out_sym,my_sym_out); }
		}
		c_x_max=dd[,my_cn]; c_x_max=c_x_max[!is.na(c_x_max)]; c_x_max=max(c_x_max);
		c_y_max=c(df[,"Symptom.birch..avg..in.season."],df[,"Symptom.birch..avg..out.of.season."]); c_y_max=c_y_max[!is.na(c_y_max)]; c_y_max=max(c_y_max)
		if(plot_all){ plot(na_in_immune,na_in_sym,col="blue",main=cn[my_cn],pch=20,xlab="abundance [in season]",ylab="symptom [in season]", xlim=c(0,c_x_max), ylim=c(0,c_y_max)) }
		else{ plot(na_in_immune,na_in_sym,col="blue",main=cn[my_cn],pch=20,xlab="abundance [in season]",ylab="symptom [in season]") }
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

#
#	PANEL_figure3B
#
PANEL_figure3B=function(req_samples=5){
	library(pheatmap); 
	library(xlsx);
	options(stringsAsFactors=FALSE)

	do_normalise=function(){
		data=read.xlsx(XLSX_FILE,encoding = "UTF-8",4); 
		data2=data
		write.table(data,"before.txt")
		normalize_data=function(data){
			pat=unique(data[,2])
			for(p in 1:length(pat)){ pat[p]=strsplit(pat[p],"_"," ")[[1]][1] }
			pat=unique(pat)
			for(s in 3:dim(data)[2]){
				for(p in 1:length(pat)){
					ix=grep(paste(pat[p],"_",sep=""),data[,2])
					my_val=data[ix,s]; 
					my_val=as.double(as.character(my_val)); 
					my_val=my_val[!is.na(my_val)]
					if(length(my_val)<=req_samples){ 
						data[ix,s]=NA 
					} else{ 
						data[ix,s]=(data[ix,s]-mean(my_val))/(sd(my_val))
					}
				}
			}
			return(data)
		}
		data=normalize_data(data); 
		write.table(data,"val_ztrans.txt"); 
		write.table(data2,"raw_vals.txt")
	}
	do_normalise()

	calc_log_change=function(c_f){
		data=read.csv(c_f,sep=" ",header=T)

		df_f_matx=data.frame(); 
		df_p_matx=data.frame()
		cn=data[,2]; norm_fac=3
		for(y in 3:dim(data)[2]){
			arr_f=c(); arr_p=c()
			for(x in 1:15){
				c1=cn[grep(paste("*_",x,"$",sep=""),cn)]
				cg1=c1[grep("b",c1)]; cg2=c1[grep("n",c1)]
				D1=subset(data,data[,2] %in% cg1); D2=subset(data,data[,2] %in% cg2)
				D1=D1[,colnames(data)[y]]; D2=D2[,colnames(data)[y]]
				D1=D1[!is.na(D1)]; D2=D2[!is.na(D2)]
				if(length(D1)>0 && length(D2)>0){
					cur_p=t.test(D1,D2)$p.value
					arr_p=c(arr_p,cur_p)
					if(mean(D2)==0){ d2_m=0.001 } else{ d2_m=mean(D2) }
					if(mean(D1)==0){ d1_m=0.001 } else{ d1_m=mean(D1) }
     		                 	cur_fac=log2(d1_m/d2_m)
					if(!is.na(cur_fac)){
						if(cur_fac>=norm_fac){ cur_fac=norm_fac }
						else if(cur_fac<=(-norm_fac)){ cur_fac=-norm_fac }
					}else{ cur_fac=0 }
     	      				arr_f=c(arr_f,cur_fac)
				}else{
					arr_f=c(arr_f,NA); 
					arr_p=c(arr_p,NA)
				}
			}
			df_f_matx=rbind(df_f_matx,arr_f); 
			df_p_matx=rbind(df_p_matx,arr_p)
		}	
		rownames(df_f_matx)=colnames(data)[3:dim(data)[2]]
		rownames(df_p_matx)=colnames(data)[3:dim(data)[2]]
		colnames(df_f_matx)=c("Nov 18-25","Dec 7-11","Jan 11-15","Feb 1-12","Mar 7-29","Apr 4-8","Apr 11-19","Apr 25-28","May 6-23","May 20 - Jun 2","Jun 27 - Jul 1","Jul 25-29","Aug 30 - Sep 2","Sep 26 - Oct 10","Oct 24-27")
		colnames(df_p_matx)=c("Nov 18-25","Dec 7-11","Jan 11-15","Feb 1-12","Mar 7-29","Apr 4-8","Apr 11-19","Apr 25-28","May 6-23","May 20 - Jun 2","Jun 27 - Jul 1","Jul 25-29","Aug 30 - Sep 2","Sep 26 - Oct 10","Oct 24-27")
		output=list()
		output$first=df_f_matx
		output$second=df_p_matx
		return(output)
	}

	do_convert=function(df_f_matx, df_p_matx){
		c_n2=c("serum.sIgE","serum.sIgG4","serum.sIgA")

		M_f=df_f_matx[c_n2,]; M_f=as.matrix(df_f_matx[c_n2,])
		M_f2=df_f_matx[!(rownames(df_f_matx)%in%(c_n2)),]; M_f2=as.matrix(df_f_matx[!(rownames(df_f_matx)%in%(c_n2)),])

		M_p=df_p_matx[c_n2,]; M_p=as.matrix(df_p_matx[c_n2,])
		M_p2=df_p_matx[!(rownames(df_p_matx)%in%(c_n2)),]; M_p2=as.matrix(df_p_matx[!(rownames(df_p_matx)%in%(c_n2)),])

		output=list()
		output$first=M_f
		output$second=M_f2
		output$third=M_p
		output$fourth=M_p2
		return(output)
	}

	output=calc_log_change("raw_vals.txt")
	df_f_matx=output[[1]]; df_p_matx=output[[2]]
	output=do_convert(df_f_matx,df_p_matx)
	M1_r<<-output[[1]]; M2_r<<-output[[2]]; M1p_r<<-output[[3]]; M2p_r<<-output[[4]]
	output=calc_log_change("val_ztrans.txt")
	df_f_matx=output[[1]]; df_p_matx=output[[2]]
	output=do_convert(df_f_matx,df_p_matx)
	M1_nr<<-output[[1]]; M2_nr<<-output[[2]]; M1p_nr<<-output[[3]]; M2p_nr<<-output[[4]]

}

#
#	PANEL_figure3
#
PANEL_figure3=function(opt){
	library("plotrix")
	library(vcd)
	options(stringsAsFactors=FALSE)

	PANEL_figure3B()

	rm(M2_r); rm(M2p_r); rm(M2_nr); rm(M2p_nr);
	rm(M1_r); rm(M1_pr); rm(M1_nr); rm(M1p_nr);

	if(opt=="nasal"){
		Mr=M2_r[dim(M2_r):1,]; Mpr=M2p_r[dim(M2p_r):1,]; Mnr=M2_nr[dim(M2_nr):1,]; Mpnr=M2p_nr[dim(M2p_nr):1,] 
	}
	else if(opt=="serum"){ 
		Mr=M1_r[dim(M1_r):1,]; Mpr=M1p_r[dim(M1p_r):1,]; Mnr=M1_nr[dim(M1_nr):1,]; Mpnr=M1p_nr[dim(M1p_nr):1,] 
	}

	svg("figure3.svg",width=15,height=15)
	plot(1:50,1:50,type="n")
	for(y in 1:ncol(Mr)){
		text(y,25,colnames(Mr)[y],adj=0,srt=44,cex=0.75)
	}
	norm_factor=3
	for(x in 1:nrow(Mr)){
		legend(16,x+0.75,cex=0.8,rownames(Mr)[x],bty="n")
		for(y in 1:ncol(Mr)){
			if(is.na(Mpnr[x,y])){ c_r=0.15 }
			else if(Mpnr[x,y]<0.001){ c_r=0.8 }
			else if(Mpnr[x,y]<0.01){ c_r=0.6 }
			else if(Mpnr[x,y]<=0.05){ c_r=0.4 }
			else{ c_r=0.2 }
			if(!is.na(Mr[x,y])){
				if(Mr[x,y]<0){ col="#029bf4" }
				else{ col="#ff0000" }
				rid=abs(as.integer(Mr[x,y]))
				intens=as.integer(20+235*min(abs(as.integer(Mr[x,y])),norm_factor)/norm_factor); intens=as.hexmode(intens)
				if(nchar(as.character(intens))==1){ intens=paste("0",intens,sep="") }
				my_col=paste(col,intens,sep="")
				draw.circle(y,x,c_r,col=my_col)
			}
		}
	}	
	offset_x=22
	offset_y=40
	text(offset_x+4,offset_y+0,"not significant (n.s.)",adj=0,cex=0.6)
	text(offset_x+4,offset_y+2,"P<=0.05",adj=0,cex=0.8)
	text(offset_x+4,offset_y+4,"P<=0.01",adj=0,cex=0.8)
	text(offset_x+4,offset_y+6,"P<=0.001",adj=0,cex=0.8)
	draw.circle(offset_x+2,offset_y+0,0.2)
	draw.circle(offset_x+2,offset_y+2,0.4)
	draw.circle(offset_x+2,offset_y+4,0.5)
	draw.circle(offset_x+2,offset_y+6,0.6)
	rect(offset_x+2,offset_y+8,offset_x+3,offset_y+9,col="red")
	rect(offset_x+2,offset_y+10,offset_x+3,offset_y+11,col="blue")
	text(offset_x+4,offset_y+8,"higher in allergic patients",adj=0)
	text(offset_x+4,offset_y+10,"higher in non-allergic patients",adj=0)
	dev.off()
}

#
#	PANEL_figure2NEW
#
PANEL_figure2NEW=function (my_opt) {
	svg("PANEL_figure2NEW.svg")

	# --------------------------------------------------------------
	cor_method="pearson"; 
	my_genome="birch"; 
	prick_type="whole_season"; 
	sel_opt="all"
        # --------------------------------------------------------------

	options(stringsAsFactors=FALSE)
	data=read.xlsx(XLSX_FILE,3)
	print(colnames(data))
	df=as.data.frame(data); par(mfrow=c(3,2)); my_opt=c("all")
	make_plot=function(label1, label2,my_col, label1_show, label2_show, y_max){
		print(c(label1,label2))
		if(mo=="all"){ ix=1:dim(df)[1]; ix1=grep("b",df[,1]); ix2=grep("na",df[,1]) }
		else if(mo=="allergic"){ ix=grep("b",df[,1]);ix1=grep("b",df[,1]) }
		else if(mo=="non-allergic"){ ix=grep("na",df[,1]);ix1=grep("na",df[,1]) }
		df2=data.frame(df[,label1],df[,label2])
		cl=cor.test(df2[ix1,1],df2[ix1,2],method=cor_method)
		b_sel=df2[ix1,]; b_sel=na.omit(b_sel);
		plot(b_sel[,1],b_sel[,2],main=as.expression(bquote(italic(r)[spearman]==.(round(cl$estimate,2)) ~ italic(p)==.(round(cl$p.value,2)))),pch=20,col=my_col,bg="red",cex=1.5,xlab=label1_show,ylab=label2_show,cex.main=1.5,cex.lab=1.5,ylim=c(0,y_max),xlim=c(0,20))
		na_sel=df2[ix2,]; na_sel=na.omit(na_sel); points(na_sel[,1],na_sel[,2],col="#029bf444",pch=20)
		abline(lm(df2[ix1,2]~df2[ix1,1]),col="red",lty=2)
	}
	if(my_genome=="birch"){
		for(mo in my_opt){
			if(sel_opt=="nose"){
				label2A="Symptom.birch..avg.nose."; label2B="Symptom.birch..avg..nose.out.of.season."; label2_show="symptom nose (avg)"
			}else{
				label2A="Symptom.birch..avg..in.season."; label2B="Symptom.birch..avg..out.of.season."; label2_show="symptom (avg)"
			}
			label1A="Pricktest..wheal..mm..in.season."; label1B="Pricktest..wheal..mm..out.of.season."
			if(mo=="all"){
				my_col=c()
				for(x in 1:dim(df)[1]){
					if(substr(df[x,1],1,1)=="b"){ my_col=c(my_col,"#ff0000dd")
					}else if(substr(df[x,1],1,1)=="n"){ my_col=c(my_col,"#029bf444") }
				}
			}
			else if(mo=="allergic"){ my_col="#FFA500FF" } else if(mo=="non-allergic"){ my_col="#FFA50033"; }
			make_plot(label1B,label2B,my_col,"wheal [mm]",label2_show,5)
			make_plot(label1A,label2A,my_col,"wheal [mm]",label2_show,10)
		}
	}
	if(my_genome=="birch"){
		for(mo in my_opt){
			if(sel_opt=="nose"){
				label2A="Symptom.birch..avg..nose."; label2B="Symptom.birch.nose.out.of.season."; label2_show="symptom nose (avg)"
			}else{
				label2A="Symptom.birch..avg..in.season."; label2B="Symptom.birch..avg..out.of.season."; label2_show="symptom (avg)"
			}
			label1A="Pricktest..flare..mm..in.season."; label1B="Pricktest..flare..mm..out.of.season."
			if(mo=="all"){
				my_col=c()
				for(x in 1:dim(df)[1]){
					if(substr(df[x,1],1,1)=="b"){ my_col=c(my_col,"#ff0000dd") }else if(substr(df[x,1],1,1)=="n"){ my_col=c(my_col,"#029bf444") }
				}
			}
			else if(mo=="allergic"){ my_col="#FFA500FF" } else if(mo=="non-allergic"){ my_col="#FFA50033" }
			make_plot(label1B,label2B,my_col,"flear [mm]",label2_show,5)
			make_plot(label1A,label2A,my_col,"flear [mm]",label2_show,10)
		}
	}
	if(my_genome=="birch" || my_genome=="grass" || my_genome=="hazel"){
		for(mo in my_opt){
			if(sel_opt=="nose"){
				label2A="Symptom..birch..avg.nose."; label2B="Symptom.birch..avg.nose.out.of.season."; label2_show="symptom nose (avg)"
			}else{
				label2A="Symptom.birch..avg..in.season."; label2B="Symptom.birch..avg..out.of.season."; label2_show="symptom (avg)"
			}
			label1="ImmunoCAP..birch..ISU.E..in.season."
			if(mo=="all"){
				my_col=c()
				for(x in 1:dim(df)[1]){ if(substr(df[x,1],1,1)=="b"){ my_col=c(my_col,"#ff0000dd") }else if(substr(df[x,1],1,1)=="n"){ my_col=c(my_col,"#029bf444") } }
			}
			else if(mo=="allergic"){ my_col="#FFA500FF" } else if(mo=="non-allergic"){ my_col="#FFA50033" }
			make_plot(label1,label2B,my_col,"ImmunoCAP [ISU-E]",label2_show,5)
			plot(1:10,1:10,type="n",xlab="",ylab="",xaxt="n",yaxt="n")
		}
	}
	dev.off()
}

#
#	PANEL_figure1
#
PANEL_figure1=function (my_opt) {
	options(stringsAsFactors=FALSE)
	library(png); library(RCurl); library(grid); library(xlsx)
	svg("figure1.svg",width=10,height=10)
	man=readPNG("man.png");
	woman=readPNG("woman.png")
	my_opt=c("hazel","alder","birch","grass")
	plot(1:1000,1:1000,type="n"); offset=0;
	
	for(mo in 1:length(my_opt)){
		FAC=3
		data=read.xlsx(XLSX_FILE,2)
		my_df=as.data.frame(data)
		my_df=subset(my_df,my_df[,6]==my_opt[mo])
		tmp=my_df[,3]
		symptome_col=as.double(tmp)
		symptome_col[as.double(tmp)>1]="#f41a02ff"
		symptome_col[as.double(tmp)>0.5 & as.double(tmp)<=1]="#f41a02ff"
		symptome_col[as.double(tmp)>0 & as.double(tmp)<=0.5]="#f41a0244"
            symptome_col[as.double(tmp)<=0]="#f41a0211"
		tmp=my_df[,2]
		immunoCAP_col=as.double(tmp)
		immunoCAP_col[as.double(tmp)>5]="#029bf4ff"
		immunoCAP_col[as.double(tmp)>3 & as.double(tmp)<=5]="#029bf4ff"
		immunoCAP_col[as.double(tmp)>2 & as.double(tmp)<=3]= "#029bf444"
		immunoCAP_col[as.double(tmp)<=2]="#029bf411"
		for(x in 1:dim(my_df)[1]){
			c_e=my_df[x,2];c_e=FAC*((c_e/5)/pi)^0.5; points(offset+50,1000-50*x,cex=c_e,pch=20,col=immunoCAP_col[x])
		}
		for(x in 1:dim(my_df)[1]){
			points(offset+75,1000-50*x,cex=FAC*(my_df[x,3]/pi)^0.5,pch=20,col=symptome_col[x])
            }
		FAC=1
		for(x in 1:dim(my_df)[1]){
			points(offset+100,1000-50*x,col='#31A72B',cex=FAC*(my_df[x,4]/pi)^0.5,pch=20)
            }
		if(mo==1){
			for(x in 1:dim(my_df)[1]){
	            	text(offset,1000-50*x,my_df[x,1],col="black")
			}
		}
		text(offset+100,100,my_opt[mo])
		offset=offset+150
		legend(700,500,lty=1,c("ImmunoCAP","Symptoms","Pollen amount"),col=c("#029bf4ff","#f41a02ff","#31A72B"),pch=20)
		if(mo==1){
			for(x in 1:dim(my_df)[1]){
				if(my_df[x,5]=="woman"){ rasterImage(woman, 600,1000-50*x-25, 620,1000-50*x+25,interpolate=FALSE)
				}else{ rasterImage(man, 600,1000-50*x-25, 620,1000-50*x+25,interpolate=FALSE)}
			}
		}
	}
	rect(25,620,150,1000); rect(150,620,300,1000); rect(300,620,450,1000); rect(450,620,600,1000)
	rect(25,200,150,620); rect(150,200,300,620); rect(300,200,450,620); rect(450,200,600,620)
	dev.off()
}

PANEL_figure1()
PANEL_figure2NEW()
PANEL_figure3("nasal")
PANEL_figure3("serum")
PANEL_figures3D_ZOOM(TRUE)
PANEL_supplemental_figure1()


