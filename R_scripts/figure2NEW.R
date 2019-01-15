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