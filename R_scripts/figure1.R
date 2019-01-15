PANEL_figure1=function (my_opt) {
	options(stringsAsFactors=FALSE)
	
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
