PANEL_figure3B=function(req_samples=5){
	
#install.packages("pheatmap")
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