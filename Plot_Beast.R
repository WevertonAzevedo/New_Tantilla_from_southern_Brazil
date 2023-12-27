##########################################################################################
# R script used to plot calibrated trees (BEAST).
# Wrote to generate figures in Azevedo et al. - Integrated evidence sheds light on the taxonomy of the widespread Tantilla melanocephala species
complex (Serpentes: Colubridae) and indicates the existence of a new species from southern South America
# Felipe G. Grazziotin - fgrazziotin@gmail.com
# &
# Weverton S. Azevedo - weverton.azevedo@hotmail.com
##########################################################################################
# This is a modular general script combining original functions and commands with several other, derived from blogs, forums, books and regular papers.
# Feel free to change, use and identify mistakes.
##########################################################################################

##########################################################################################
# Load libraries
##########################################################################################
library(ape)
library(phytools)
library(strap)
##########################################################################################

##########################################################################################
# Externally manipulate tree files before open them in R
##########################################################################################

# first: open the summary tree of Beast in Figtree and save as nexus.

# second: run the commands below in a pain text editor with grep.

 (\[&CAheight_95%_HPD={[\d\.,E-]+)},[a-z_=\d+\.,{}%-]+\]
 replace
 \1]


 \[&height[a-z_=\d+\.,{}%-]+\]
 replace


 \[&CAheight_95%_HPD={([\d\.E-]+),([\d\.E-]+)}\]
 replace
 \1_\2
##########################################################################################

##########################################################################################
# Plot complete calibrated tree
##########################################################################################

# read tree file
t1<-read.nexus("out_ed_clean.tre")

# get max 95% CI
max95<-gsub("\\d+.[0-9E-]+_","",t1$node.label)

# get min 95% CI
min95<-gsub("(\\d+.[0-9E-]+)_\\d+.[0-9E-]+","\\1",t1$node.label)

# set root time
t1$root.time<-max(node.depth.edgelength(t1)) 
lim<-t1$root.time

#######################################
# hard option to plot using plot.phylo
#######################################

#limit timescales chart to pruned tree root + some value to assure space for root CI bar 
redTimeScale<-timescales$ICS2015[timescales$ICS2015[,2]<=lim+15,]

# create vector with tips name
tips_for_tree<-gsub("_"," ",t1$tip.label)

# load table with tips colors
tips_col_table<-read.csv("tip_colors.csv")

# selects unique values of tips colors
unique_colors<-unique(tips_col_table[,2])

# create vector to select the colors of tips
tip_color<-rep("black",Ntip(t1))

# loop to set tips colors based on table information
for(i in seq(unique_colors)){
	tip_color[t1$tip.label%in%tips_col_table[,1][tips_col_table[,2]==unique_colors[i]]]<-unique_colors[i]
}

 # create vector to select branch colors
branch_color<-rep("black",Nedge(t1))

# loop to set branch colors based on table information
for(i in seq(unique_colors)){
	branch_color[which(t1$edge[,2]%in%c(getMRCA(t1, which(t1$tip.label%in%tips_col_table[tips_col_table[,2]==unique_colors[i],1])),getDescendants(t1,getMRCA(t1, which(t1$tip.label%in%tips_col_table[tips_col_table[,2]==unique_colors[i],1])))))]<-unique_colors[i]
}
#######################################


#######################################
# plot tree
#######################################
pdf("Best.pdf",family="ArialMT", width=8.25, height=11.75)

	# set layouts
	layout(matrix(c(1,3,4,2,5,6),3,2),widths=c(5,3),heights=c(11,0.8,0.2))

	# set margins
	par(oma=rep(0,4),mar=rep(0,4))

	# plot tree to get coordinates
	plot.phylo(t1, x.lim=c(0,lim), y.lim=c(1,Ntip(t1)+1), type="phy", no.margin=TRUE, edge.width=1, show.tip.label=F, plot=F)

	# get coordinates
	obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)

	# plot age bars
	for(i in 1:length(redTimeScale[,5][redTimeScale[,5]=="Epoch"])){
		if(i%%2==0){
			rect(max(obj$xx)-redTimeScale[,1][redTimeScale[,5]=="Epoch"][i],-3.5, max(obj$xx)-redTimeScale[,2][redTimeScale[,5]=="Epoch"][i], max(obj$yy), col=rgb(240,240,240, maxColorValue=255), border=NA)
		}else{
			rect(max(obj$xx)-redTimeScale[,1][redTimeScale[,5]=="Epoch"][i],-3.5, max(obj$xx)-redTimeScale[,2][redTimeScale[,5]=="Epoch"][i], max(obj$yy), col=rgb(255,255,255, maxColorValue=255), border=NA)
		}
		segments(max(obj$xx)-redTimeScale[,2][redTimeScale[,5]=="Epoch"][i],-3.5,max(obj$xx)-redTimeScale[,2][redTimeScale[,5]=="Epoch"][i],-1)
		if(redTimeScale[,2][redTimeScale[,5]=="Epoch"][i]<10){
			text(max(obj$xx)-redTimeScale[,2][redTimeScale[,5]=="Epoch"][i],-0.8,round(redTimeScale[,2][redTimeScale[,5]=="Epoch"][i],1),cex=0.6,srt=90,adj=c(0,0.5))
		}else{
		text(max(obj$xx)-redTimeScale[,2][redTimeScale[,5]=="Epoch"][i],-0.8,round(redTimeScale[,2][redTimeScale[,5]=="Epoch"][i],0),cex=0.6,srt=90,adj=c(0,0.5))
		}
	}

	# # set parameters
	par(new=T)

	# plot tree without tip labels
	plot.phylo(t1, x.lim=c(0,lim), y.lim=c(1,Ntip(t1)+1), type="phy", no.margin=TRUE, edge.width=1, edge.color=branch_color, show.tip.label=F)

	for(i in seq(min95)){
		rect(t1$root.time-as.numeric(min95[i]),obj$yy[i+Ntip(t1)]+.175,t1$root.time-as.numeric(max95[i]),obj$yy[i+length(t1$tip.label)]-.175, col=rgb(0.1,0.9,1,0.3), lwd=.25, border=NA)
	}

	par(new=F, mar=rep(0,4),xpd = NA)
	# plot to set layout coordinates
	plot.phylo(t1, x.lim=c(0,lim), y.lim=c(1,Ntip(t1)+1), type="phy", no.margin=TRUE, edge.width=1, show.tip.label=F, plot=F)
	
	# plot tip labels
	text(-8,seq(Ntip(t1)),tips_for_tree,cex=0.9,pos=4,font=3,col=tip_color)

	# set parameters
	par(new=F,mar=rep(0,4))
	# plot to set layout coordinates
	plot(0,xlim=obj$x.lim,ylim=obj$y.lim,axes=F,type="n")
	# plot geo chart
	for(i in 1:length(redTimeScale[,5][redTimeScale[,5]=="Period"])){
		rect(max(obj$xx)-redTimeScale[,1][redTimeScale[,5]=="Period"][i], 0, max(obj$xx)-redTimeScale[,2][redTimeScale[,5]=="Period"][i], max(obj$yy)/3, col=rgb(redTimeScale[,10][redTimeScale[,5]=="Period"][i], redTimeScale[,11][redTimeScale[,5]=="Period"][i], redTimeScale[,12][redTimeScale[,5]=="Period"][i], maxColorValue=255))
		if(redTimeScale[,8][redTimeScale[,5]=="Period"][i]!="Quat."){
			text(((max(obj$xx)-redTimeScale[,1][redTimeScale[,5]=="Period"][i])+(max(obj$xx)-redTimeScale[,2][redTimeScale[,5]=="Period"][i]))/2, (max(obj$yy)/3)/2,redTimeScale[,8][redTimeScale[,5]=="Period"][i],cex=1)
		}
	}	
	for(i in 1:length(redTimeScale[,5][redTimeScale[,5]=="Epoch"])){
		rect(max(obj$xx)-redTimeScale[,1][redTimeScale[,5]=="Epoch"][i], max(obj$yy)/3, max(obj$xx)-redTimeScale[,2][redTimeScale[,5]=="Epoch"][i], max(obj$yy)+5, col=rgb(redTimeScale[,10][redTimeScale[,5]=="Epoch"][i], redTimeScale[,11][redTimeScale[,5]=="Epoch"][i], redTimeScale[,12][redTimeScale[,5]=="Epoch"][i], maxColorValue=255))
		if(redTimeScale[,7][redTimeScale[,5]=="Epoch"][i]!="Quaternary"){
			text(((max(obj$xx)-redTimeScale[,1][redTimeScale[,5]=="Epoch"][i])+(max(obj$xx)-redTimeScale[,2][redTimeScale[,5]=="Epoch"][i]))/2, ((max(obj$yy)/3)+max(obj$yy))/2,redTimeScale[,8][redTimeScale[,5]=="Epoch"][i],cex=0.8,srt=90)
		}
	}

	# plot to set layout coordinates
	plot(0,xlim=obj$x.lim,ylim=obj$y.lim,axes=F,type="n")

	# set tick number
	Nticks<-round(max(obj$xx)/10,0)

	# plot ticks for the geo chart
	for(i in 0:Nticks){
		segments(max(obj$xx)-(i*10),max(obj$yy)/2,max(obj$xx)-(i*10),max(obj$yy)+5)
		text(max(obj$xx)-(i*10),(max(obj$yy)/3),i*10,srt=90,adj=c(0.6,0.5),cex=0.9)
	}
	
dev.off()
#######################################

##########################################################################################

##########################################################################################
# Plot reduced calibrated tree
##########################################################################################
# remove all objects from the environment 
rm(list = ls())

# read tree file
t1<-read.nexus("out_ed_clean.tre")

# get the MRCA for the ingroup

t1p<-extract.clade(t1,getMRCA(t1,c("Scolecophis_atrocinctus","Tantilla_supracincta")))
temp<-extract.clade(t1,getMRCA(t1,c("Scolecophis_atrocinctus","Tantilla_supracincta")))
# write tree to manually rotate clades

t1p<-rotate(t1p,getMRCA(t1p,c("Tantilla_gracilis","Tantilla_supracincta")))
t1p<-rotate(t1p,getMRCA(t1p,c("Tantilla_coronata","Tantilla_relicta")))
t1p<-rotate(t1p,getMRCA(t1p,c("Tantilla_melanocephala_BRA_RR_Boa_Vista_IBSP87381","Tantilla_melanocephala_Guyana1")))
t1p<-rotate(t1p,getMRCA(t1p,c("Tantilla_aff_melanocephala_BRA_RS_Rosario_do_Sul_IBSP90139","Tantilla_melanocephala_BRA_SP_Itu_IBSP86599")))

# get max 95% CI
max95<-gsub("\\d+.[0-9E-]+_","",t1p$node.label)

# get min 95% CI
min95<-gsub("(\\d+.[0-9E-]+)_\\d+.[0-9E-]+","\\1",t1p$node.label)

# set root time
t1p$root.time<-max(node.depth.edgelength(t1p)) 
lim<-t1p$root.time

node.dates<-round(t1p$root.time-node.depth.edgelength(t1p)[(Ntip(t1p)+1):(Ntip(t1p)+Nnode(t1p))],1)

#######################################
# hard option to plot using plot.phylo
#######################################

#limit timescales chart to pruned tree root + some value to assure space for root CI bar 
redTimeScale<-timescales$ICS2015[timescales$ICS2015[,2]<=lim+15,]

# create vector with tips name
tips_for_tree<-t1p$tip.label

tips_for_tree<-tips_for_tree[t1p$edge[,2][t1p$edge[,2]<=Ntip(t1p)]]

# load table with tips colors
tips_col_table<-read.csv("tip_colors.csv")

# selects unique values of tips colors
unique_colors<-unique(tips_col_table[,2])

# create vector to select the colors of tips
tip_color<-rep("black",Ntip(t1p))

# loop to set tips colors based on table information
for(i in seq(unique_colors)){
	tip_color[tips_for_tree%in%tips_col_table[,1][tips_col_table[,2]==unique_colors[i]]]<-unique_colors[i]
}

 # create vector to select branch colors
branch_color<-rep("black",Nedge(t1p))

# loop to set branch colors based on table information
for(i in seq(unique_colors)){
	branch_color[which(t1p$edge[,2]%in%c(getMRCA(t1p, which(t1p$tip.label%in%tips_col_table[tips_col_table[,2]==unique_colors[i],1])),getDescendants(t1p,getMRCA(t1p, which(t1p$tip.label%in%tips_col_table[tips_col_table[,2]==unique_colors[i],1])))))]<-unique_colors[i]
}

# change tip labels
tips_for_tree<-gsub("_"," ",tips_for_tree)
tips_for_tree<-gsub(" IBSP\\d+","",tips_for_tree)
tips_for_tree<-gsub(" UFMG\\d+","",tips_for_tree)
tips_for_tree<-gsub(" SMF97636","",tips_for_tree)
tips_for_tree<-gsub(" JAC29265","",tips_for_tree)
tips_for_tree<-gsub(" Sao "," São ",tips_for_tree)
tips_for_tree<-gsub(" Rosario"," Rosário",tips_for_tree)
tips_for_tree<-gsub(" Jeronimo"," Jerônimo",tips_for_tree)
tips_for_tree<-gsub(" Carapicuiba"," Carapicuíba",tips_for_tree)
tips_for_tree<-gsub(" aff"," aff.",tips_for_tree)

tips_for_tree<-gsub(" melanocephala BRA SP"," cf. melanocephala BRA SP",tips_for_tree)

tips_for_tree<-gsub("Guyana1","GUY Aishalton",tips_for_tree)
tips_for_tree<-gsub("Guyana2","GUY Dubulay",tips_for_tree)
tips_for_tree<-gsub("Kourou French Guiana","FGU Kourou",tips_for_tree)
tips_for_tree<-gsub("Trinidad","TTO Trinidad",tips_for_tree)
tips_for_tree<-gsub("Tobago","TTO Tobago",tips_for_tree)
tips_for_tree<-gsub("Macuro Venezuela","VEN Macuro",tips_for_tree)

tips_for_tree<-gsub("melanocephala$","cf. melanocephala BRA SP São Paulo",tips_for_tree)

tips_for_tree<-gsub("cf boipiranga"," boipiranga",tips_for_tree)


#######################################


#######################################
# plot tree
#######################################
pdf("Best_reduced.pdf",family="ArialMT", width=8.25, height=7.83)

	# set layouts
	layout(matrix(c(1,3,4,2,5,6),3,2),widths=c(5,2.2),heights=c(11,1.0,0.5))
	
	####################
	# set margins
	par(oma=rep(0,4),mar=rep(0,4))

	# plot tree to get coordinates
	plot.phylo(t1p, x.lim=c(0,lim), y.lim=c(1,Ntip(t1p)+1), type="phy", no.margin=TRUE, edge.width=1, show.tip.label=F, plot=F)
	
	# get coordinates
	obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)

	# plot age bars
	for(i in 1:length(redTimeScale[,5][redTimeScale[,5]=="Age"])){
		if(i%%2==0){
			rect(max(obj$xx)-redTimeScale[,1][redTimeScale[,5]=="Age"][i],-3.5, max(obj$xx)-redTimeScale[,2][redTimeScale[,5]=="Age"][i], max(obj$yy), col=rgb(240,240,240, maxColorValue=255), border=NA)
		}else{
			rect(max(obj$xx)-redTimeScale[,1][redTimeScale[,5]=="Age"][i],-3.5, max(obj$xx)-redTimeScale[,2][redTimeScale[,5]=="Age"][i], max(obj$yy), col=rgb(255,255,255, maxColorValue=255), border=NA)
		}
		segments(max(obj$xx)-redTimeScale[,2][redTimeScale[,5]=="Age"][i],0.05, max(obj$xx)-redTimeScale[,2][redTimeScale[,5]=="Age"][i],-1)
		if(redTimeScale[,2][redTimeScale[,5]=="Age"][i]<10){
			text(max(obj$xx)-redTimeScale[,2][redTimeScale[,5]=="Age"][i],0.15,round(redTimeScale[,2][redTimeScale[,5]=="Age"][i],1),cex=0.7,srt=90,adj=c(0,0.5))
		}else{
		text(max(obj$xx)-redTimeScale[,2][redTimeScale[,5]=="Age"][i],0.15,round(redTimeScale[,2][redTimeScale[,5]=="Age"][i],0),cex=0.7,srt=90,adj=c(0,0.5))
		}
	}

	####################
	# # set parameters
	par(new=T)
	
	# plot tree without tip labels
	plot.phylo(t1p, x.lim=c(0,lim), y.lim=c(1,Ntip(t1p)+1), type="phy", no.margin=TRUE, edge.width=1.4, edge.color=branch_color, show.tip.label=F)
	
	# plot CI bars
	for(i in seq(min95)){
		rect(t1p$root.time-as.numeric(min95[i]),obj$yy[i+Ntip(t1p)]+.1,t1p$root.time-as.numeric(max95[i]),obj$yy[i+length(t1p$tip.label)]-.1, col=rgb(0.1,0.9,1,0.8), lwd=.25, border=NA)
	}
	# plot node dates
	nodelabels(node.dates,frame="none", adj=c(1.2,1.4), cex=1)

	####################
	par(new=F, mar=rep(0,4),xpd = NA)

	# plot to set layout coordinates
	plot.phylo(t1p, x.lim=c(0,lim), y.lim=c(1,Ntip(t1p)+1), type="phy", no.margin=TRUE, edge.width=1, show.tip.label=F, plot=F)
	
	# plot tip labels
	text(-3.2,seq(Ntip(t1p)),tips_for_tree,cex=1,pos=4,font=3,col=tip_color)

	####################
	# set parameters
	par(new=F,mar=rep(0,4))
	
	# plot to set layout coordinates
	plot(0,xlim=obj$x.lim,ylim=obj$y.lim,axes=F,type="n")
	
	# plot geo chart
	for(i in 1:length(redTimeScale[,5][redTimeScale[,5]=="Epoch"])){
		if(i%%2==0){			
			rect(max(obj$xx)-redTimeScale[,1][redTimeScale[,5]=="Epoch"][i], 0, max(obj$xx)-redTimeScale[,2][redTimeScale[,5]=="Epoch"][i], max(obj$yy)/3, col=rgb(240,240,240, maxColorValue=255))
		}else{			
			rect(max(obj$xx)-redTimeScale[,1][redTimeScale[,5]=="Epoch"][i], 0, max(obj$xx)-redTimeScale[,2][redTimeScale[,5]=="Epoch"][i], max(obj$yy)/3, col=rgb(255,255,255, maxColorValue=255))
		}
		if(redTimeScale[,8][redTimeScale[,5]=="Epoch"][i]!="Quat."){
			text(((max(obj$xx)-redTimeScale[,1][redTimeScale[,5]=="Epoch"][i])+(max(obj$xx)-redTimeScale[,2][redTimeScale[,5]=="Epoch"][i]))/2, (max(obj$yy)/3)/2,redTimeScale[,8][redTimeScale[,5]=="Epoch"][i],cex=1)
		}
	}




	for(i in 1:length(redTimeScale[,5][redTimeScale[,5]=="Age"])){
		if(i%%2==0){			
			rect(max(obj$xx)-redTimeScale[,1][redTimeScale[,5]=="Age"][i], max(obj$yy)/3, max(obj$xx)-redTimeScale[,2][redTimeScale[,5]=="Age"][i], max(obj$yy)+5, col=rgb(240,240,240, maxColorValue=255))
		}else{			
			rect(max(obj$xx)-redTimeScale[,1][redTimeScale[,5]=="Age"][i], max(obj$yy)/3, max(obj$xx)-redTimeScale[,2][redTimeScale[,5]=="Age"][i], max(obj$yy)+5, col=rgb(255,255,255, maxColorValue=255))
		}
		text(((max(obj$xx)-redTimeScale[,1][redTimeScale[,5]=="Age"][i])+(max(obj$xx)-redTimeScale[,2][redTimeScale[,5]=="Age"][i]))/2, ((max(obj$yy)/3)+max(obj$yy))/2,redTimeScale[,8][redTimeScale[,5]=="Age"][i],cex=0.6,srt=90)
	}

	####################
	# plot to set layout coordinates
	plot(0,xlim=obj$x.lim,ylim=obj$y.lim,axes=F,type="n")

	# set tick number
	Nticks<-round(max(obj$xx)/10,0)
	
dev.off()

#######################################

##########################################################################################
