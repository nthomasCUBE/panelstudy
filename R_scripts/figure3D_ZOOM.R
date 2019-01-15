PANEL_figures3D_ZOOM=function(show_z_scale){
	options(stringsAsFactors=FALSE)
	svg(paste("PANEL_figure3D_ZOOM_",show_z_scale,".svg",sep=""),width=20,height=10)

	#grid.newpage()
	if(show_z_scale==TRUE){
		data=read.xlsx(XLSX_FILE,encoding = "UTF-8",5)
		for(x in 3:dim(data)[2]){ data[,x]=as.numeric(data[,x]) }
	}else{
		data=read.xlsx(XLSX_FILE,encoding = "UTF-8",4) 
	}

	print(colnames(data))

	my_cols=c("sIgA","sIgG4","sIgG4","IgA","FLC.λ","FLC.κ","IL.33","MIP.1β","Eotaxin.2","Eotaxin.2","Eotaxin.2","Eotaxin.2","IL.1β","MCP.1","MCP.1")
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
