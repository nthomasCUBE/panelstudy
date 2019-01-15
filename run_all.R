getwd()
#setwd("C:/Users/goekkame/Desktop/PAPER/R-Code/Folder_of_R_code")
options(stringsAsFactors=FALSE)
options(java.parameters = "-Xmx4000m")

library(bitops)
library(ggpubr)
library(grid);
library(gridExtra);
library(magrittr)
library(png);
library(RCurl);
library(ggplot2);
library(grid);
library(xlsx);

#XLSX_FILE="../../PANEL_STUDY_REPORT/Tables_25apr18.xlsx"
XLSX_FILE="20180731_Panelstudy_Tables.xlsx"

source("R_scripts/figure1.R")
source("R_scripts/figure2NEW.R")
source("R_scripts/figure3.R")
source("R_scripts/figure3B.R")
source("R_scripts/figure3D_ZOOM.R")
source("R_scripts/suppl_figure1.R")

run_all=function(){
	system("rm *svg")

	#PANEL_figure1()
	#PANEL_figure2NEW()
	#PANEL_figure3("nasal") 
	#PANEL_figure3("serum") 
	#PANEL_figures3D_ZOOM(TRUE)
	PANEL_supplemental_figure1(TRUE) 
}
run_all()

