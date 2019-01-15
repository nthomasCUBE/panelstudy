PANEL_figure3=function(opt){
#install.packages("plotrix")
#install.packages("vcd")
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