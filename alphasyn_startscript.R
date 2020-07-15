################################################################################
###Transcriptomics - alphasyn timeseries
################################################################################


###Set your work directory
setwd("~/Desktop/WUR Eleni/R/Gegenereerde_data/") ###set this to the folder you want to use in R
workwd <- getwd()
filename <- "alphasynRIL"

###Load pre-made functions
    #uses eQTL pipeline functions https://git.wur.nl/mark_sterken/eQTL_pipeline
    #     transcriptomics functions https://git.wur.nl/mark_sterken/Transcriptomics.workbench
git_dir <- "~/Desktop/WUR Eleni/R/R_scripts/"  #the folder in which you place the NEMA_functions folder
source(paste(git_dir,"/NEMA_functions/Loader.R",sep=""))



################################################################################
###Dependencies
################################################################################

install <- FALSE  ####first time you run it; set to true to install packages!
if(install){
    install.packages("tidyverse")
    install.packages("colorspace")
    install.packages("RColorBrewer")
    install.packages("BiocManager")
    BiocManager::install("limma")
    BiocManager::install("statmod")
    install.packages("gridExtra")
    install.packages("VennDiagram")
    install.packages("openxlsx")
    install.packages("rmarkdown")
}

###load
    library("colorspace")
    library("RColorBrewer")
    library(limma)
    library(gridExtra)
    library("VennDiagram")
    library(openxlsx)
    library("rmarkdown")
    library(tidyverse)


################################################################################
###Plotting theme, colours
################################################################################


###Set plotting theme
    presentation <- theme(axis.text.x = element_text(size=10, face="bold", color="black"),
                          axis.text.y = element_text(size=10, face="bold", color="black"),
                          axis.title.x = element_text(size=12, face="bold", color="black"),
                          axis.title.y = element_text(size=12, face="bold", color="black"),
                          strip.text.x = element_text(size=12, face="bold", color="black"),
                          strip.text.y = element_text(size=12, face="bold", color="black"),
                          plot.title = element_text(size=14, face="bold"),
                          strip.background = element_rect(fill= "grey80",color="black"),
                          panel.background = element_rect(fill = "white",color="black"),
                          panel.grid.major = element_line(colour = "grey80"),
                          panel.grid.minor = element_blank(),
                          legend.position = "right")


    blank_theme <- theme(plot.background = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.border = element_blank(),
                         panel.background = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.ticks = element_blank())



###Here you can set colours for plotting in theme using ggplot2
    #display.brewer.all()
    myColors <- c(brewer.pal(9,"Set1")[c(2,5)],brewer.pal(9,"Purples")[c(4,6,6)],"black","darkgrey","black",brewer.pal(12,"Paired")[c(1,7)])
    
    
    names(myColors) <- c("CB4856","N2","aSlight","aS","trans","cis","notsig","RIL","SCH4856","NL5901")
    
    colScale <- scale_colour_manual(name = "Treatment",values = myColors)
    fillScale <- scale_fill_manual(name = "Treatment",values = myColors)

################################################################################
###load all the data
################################################################################

    ###Load normalized data
        ###all quality is good
        load(file = "../Raw_data/aSRIL_datafiles/obj_list.data.Rdata")

    ###populations
        load(file = "../Raw_data/aSRIL_datafiles/obj_popmap.Rdata")
        load(file = "../Raw_data/aSRIL_datafiles/obj_popmrk.Rdata")
    
    ###eQTL mapping
        load(file="../Raw_data/aSRIL_datafiles/obj_aS.eQTL.table.Rdata")
        load(file="../Raw_data/aSRIL_datafiles/obj_peak.aS.QTL.Rdata")
        load(file="../Raw_data/aSRIL_datafiles/obj_aS.QTL.Rdata")
        
    ###load qPCR data
        load(file = "../Raw_data/aSRIL_datafiles/obj_qpcr_data.Rdata")
        load(file = "../Raw_data/aSRIL_datafiles/obj_qpcr.QTL.Rdata")
        
################################################################################
###plot ILs
################################################################################

    strain1 <- wur.pop.map[,53]; strain1[strain1==-1] <- 2
    strain2 <- wur.pop.map[,69]
    
    pdf("F1.pdf",width = 4,height=4)
        plot.genotype.strain(wur.pop.marker,strain.map1 = strain1,strain.map2 = strain2,both.chr = TRUE,col.lines=c("#1F78B4","grey","#FF7F00","darkgreen"))
    dev.off()
        
    pdf("P1.pdf",width = 8,height=4)   
        par(mfrow=c(1,2))
        plot.genotype.strain(wur.pop.marker,strain.map1 = strain1,both.chr = TRUE,col.lines=c("#1F78B4","grey","#FF7F00","darkgreen"))
        plot.genotype.strain(wur.pop.marker,strain.map1 = strain2,both.chr = TRUE,col.lines=c("#1F78B4","grey","#FF7F00","darkgreen"))
    dev.off()
          
    pdf("F2_after_selfing.pdf",width = 16,height=4)
        par(mfrow=c(1,4))
        plot.genotype.strain(wur.pop.marker,strain.map1 = rep(1,length(strain1)),both.chr = TRUE,col.lines=c("#1F78B4","grey","#FF7F00","darkgreen"))
        plot.genotype.strain(wur.pop.marker,strain.map1 = strain1*strain2,both.chr = TRUE,col.lines=c("#1F78B4","grey","#FF7F00","darkgreen"))
        plot.genotype.strain(wur.pop.marker,strain.map1 = strain1,strain.map2 = strain2,both.chr = TRUE,col.lines=c("#1F78B4","grey","#FF7F00","darkgreen"))
        plot.genotype.strain(wur.pop.marker,strain.map1 = strain2,strain.map2 = strain1,both.chr = TRUE,col.lines=c("#1F78B4","grey","#FF7F00","darkgreen"))
    dev.off()
                        

################################################################################
###Polymorphic genes in a region
################################################################################

    CB4856.DB[[1]]
    
    aS.eQTL.table[1,]
    
    dplyr::filter(Agilent.Ce.V2,chromosome=="V",!is.na(chromosome),
                                gene_bp_end > 15397385, gene_bp_start < 18222746) %>%
    dplyr::select(gene_sequence_name) %>%
    dplyr::filter(!duplicated(gene_sequence_name))
    

    genes <- dplyr::filter(Agilent.Ce.V2,chromosome=="V",!is.na(chromosome),
                                         gene_bp_end > 15397385, gene_bp_start < 18222746) %>%
             dplyr::select(gene_sequence_name) %>%
             dplyr::filter(!duplicated(gene_sequence_name)) %>%
             unlist() %>%
             as.character()

    filter(CB4856.DB[[2]],Sequence_name %in% genes)            
        
################################################################################
###plot QTL profile
################################################################################
        
    ###ignore the warnings; you can select any gene
        head(aS.eQTL.table)
    
    
        data.plot <- prep.ggplot.QTL.profile(peak.aS.QTL,aS.QTL,"AGIWUR10019")
        data.plot[[2]] <- mutate(data.plot[[2]], geno_strain=ifelse(genotype==-1,"CB4856","N2"))
        
        sf2a <- ggplot(data.plot$QTL_profile,aes(x=qtl_bp,y=qtl_significance,alpha=0.2)) +
                geom_line(size=1.5,colour=brewer.pal(9,"Set1")[3]) + facet_grid(~qtl_chromosome,scales="free",space="free_x") + presentation + theme(legend.position = "none") +
                geom_abline(intercept=3.1,slope=0,linetype=2,size=1) + labs(x="QTL position (Mb)",y="significance (-log10(p))") +
                scale_x_continuous(breaks=c(5,10,15,20)*10^6,labels=c(5,10,15,20)) + ylim(0,5.5)
        
        sf2b <- ggplot(data.plot[[2]],aes(x=geno_strain,y=trait_value)) +
                geom_jitter(height=0,width=0.25,aes(colour=geno_strain),alpha=0.2) + geom_boxplot(outlier.shape=NA,alpha=0.2,aes(fill=geno_strain)) +
                xlab("Genotype\nat marker") + ylab(data.plot[[2]]$trait[1]) + facet_grid(~Chromosome+Position) +
                presentation + colScale + fillScale + theme(legend.position = "none") +
                annotate("text",x=1.5,y=max(data.plot[[2]]$trait_value,na.rm=T),label=paste("italic(R)^{2}==",round(data.plot[[2]]$R_squared[1],digits=2),sep=""),parse=TRUE)
             
        annotation.grobA <- title.grob <- textGrob(label = "A",x = unit(0, "lines"),y = unit(0, "lines"),hjust = 0, vjust = 0,gp = gpar(fontsize = 20,fontface="bold"))
        annotation.grobB <- title.grob <- textGrob(label = "B",x = unit(0, "lines"),y = unit(0, "lines"),hjust = 0, vjust = 0,gp = gpar(fontsize = 20,fontface="bold"))

        sf2a <- arrangeGrob(sf2a,top=annotation.grobA)
        sf2b <- arrangeGrob(sf2b,top=annotation.grobB)

        grid.arrange(sf2a,sf2b,widths=c(2,1))
        
        pdf("Figure-eQTL_profile.pdf",width=8,height=4)        
            grid.arrange(sf2a,sf2b,widths=c(2,1))
        dev.off()  

################################################################################
###Select trans-bands
################################################################################
        
transtable <- dplyr::filter(aS.eQTL.table,trans_band != "none") %>%
        dplyr::filter(!duplicated(trans_band)) %>%
        dplyr::select(trans_band) %>%
        tidyr::separate(trans_band,into=c("chromosome","rest"),sep=":",remove = FALSE) %>%
        tidyr::separate(rest,into=c("loc_left","loc_right"),sep="-") %>%
        dplyr::mutate(chromosome=gsub("chr","",chromosome), loc_right=gsub("Mb","",loc_right)) %>%
        dplyr::mutate(loc_left=1e6*as.numeric(loc_left),loc_right=1e6*as.numeric(loc_right))
        
        output <- NULL       
        for(i in 1:nrow(transtable)){
            genes <- dplyr::filter(Agilent.Ce.V2,chromosome==transtable[i,2],!is.na(chromosome),
                                   gene_bp_end > transtable[i,3], gene_bp_start < transtable[i,4]) %>%
                dplyr::select(gene_sequence_name) %>%
                dplyr::filter(!duplicated(gene_sequence_name)) %>%
                unlist() %>%
                as.character()
            
            tmp <- filter(CB4856.DB[[2]],Sequence_name %in% genes)  
            
            output <- rbind(output,
                            cbind(transtable[i,],tmp)
            )
        }
        
        ###check if is in list
        filter(output,trans_band %in% c("chrV:14.5-17.0Mb","chrV:5.5-6.0Mb"))
        filter(output,Sequence_name %in% filter(aS.eQTL.table, qtl_type == "cis")$gene_public_name)
        
        ###check if is not in a list
        filter(output,!trans_band %in% c("chrV:14.5-17.0Mb","chrV:5.5-6.0Mb"))
        
        ###check if larger or smaller then
        filter(output, loc_left < 1000000)
        
        ###check if not equal to
        filter(output, chromosome != "V")
        
        
        
################################################################################
###check gene expression at trans-band
################################################################################
        
        transtable <- dplyr::filter(aS.eQTL.table,trans_band != "none") %>%
            dplyr::filter(!duplicated(trans_band)) %>%
            dplyr::select(trans_band) %>%
            tidyr::separate(trans_band,into=c("chromosome","rest"),sep=":",remove = FALSE) %>%
            tidyr::separate(rest,into=c("loc_left","loc_right"),sep="-") %>%
            dplyr::mutate(chromosome=gsub("chr","",chromosome), loc_right=gsub("Mb","",loc_right)) %>%
            dplyr::mutate(loc_left=1e6*as.numeric(loc_left),loc_right=1e6*as.numeric(loc_right))
        
        
        for(i in 1:nrow(transtable)){ 
            
            data.plot <- filter(aS.eQTL.table,trans_band == transtable[i,1])
            
            #ggplot(data.plot,aes(x=qtl_effect,y=qtl_significance,size=qtl_R2_sm)) +
            #geom_point()
            print(
                ggplot(data.plot,aes(x=qtl_effect,y=PL_alpha_effect)) +
                    geom_point() + geom_smooth(method="lm")
            )
        }
        
        
        CB4856.DB[[1]]
        
        aS.eQTL.table[1,]
        
        dplyr::filter(Agilent.Ce.V2,chromosome=="V",!is.na(chromosome),
                      gene_bp_end > 15397385, gene_bp_start < 18222746) %>%
            dplyr::select(gene_sequence_name) %>%
            dplyr::filter(!duplicated(gene_sequence_name))
        
        
        genes <- dplyr::filter(Agilent.Ce.V2,chromosome=="V",!is.na(chromosome),
                               gene_bp_end > 15397385, gene_bp_start < 18222746) %>%
            dplyr::select(gene_sequence_name) %>%
            dplyr::filter(!duplicated(gene_sequence_name)) %>%
            unlist() %>%
            as.character()
        
        filter(CB4856.DB[[2]],Sequence_name %in% genes)            
        
################################################################################
###check gene expression at trans-band
################################################################################
        
        filter(aS.eQTL.table,trait=="AGIWUR1000") %>%
            merge(list.data,by.x=1,by.y=2)        
        
        
        
##############################################################################
###Transbanden uit Wormnet 
###########################################################################        
       
###transband ChI
        
        WormNetChI <- read.delim("Trans_Band_ChI.txt", header = FALSE)        
        head(WormNetChI)
        filter(output, chromosome == "I", Sequence_name %in% WormNetChI$V1)

        outputChI <- filter(output, chromosome == "I", Sequence_name %in% WormNetChI$V1)        
      
        ###Filter op p-value 0.05
        WormNetChI_p_value <- dplyr::filter(WormNetChI,V5<0.20)
            
        outputChI_p_value <-  dplyr::filter(outputChI, Sequence_name %in% dplyr::filter(WormNetChI,V5<0.20)$V1)
        
        outputChI_p_value <-  merge(outputChI, WormNetChI_p_value,by.x=5,by.y=1)

        
###transband ChII        
        
        WormNetChII <- read.delim("Trans_Band_ChII", header = FALSE)        
        head(WormNetChII)
        filter(output, chromosome == "II", Sequence_name %in% WormNetChII$V1)
        
        outputChII <- filter(output, chromosome == "II", Sequence_name %in% WormNetChII$V1)        
        
        ###Filter op p-value 0.05
        WormNetChII_p_value <- dplyr::filter(WormNetChII,V5<0.20)
        
        outputChII_p_value <-  dplyr::filter(outputChII, Sequence_name %in% dplyr::filter(WormNetChII,V5<0.20)$V1)
        outputChII_p_value <-  merge(outputChII, WormNetChII_p_value,by.x=5,by.y=1)
        
        
###transband ChIV
        
        WormNetChIV <- read.delim("Trans_Band_ChIV", header = FALSE)        
        head(WormNetChIV)
        filter(output, chromosome == "IV", Sequence_name %in% WormNetChIV$V1)
        
        outputChIV <- filter(output, chromosome == "IV", Sequence_name %in% WormNetChIV$V1)        
        
        ###Filter op p-value 0.05
        WormNetChIV_p_value <- dplyr::filter(WormNetChIV,V5<0.20)
        
        outputChIV_p_value <-  dplyr::filter(outputChIV, Sequence_name %in% dplyr::filter(WormNetChIV,V5<0.20)$V1)
        outputChIV_p_value <-  merge(outputChIV, WormNetChIV_p_value,by.x=5,by.y=1)
        
        
###transband ChV11 deze was > 300
        
        ###Een lang lijst met 4 x random de eerste 300 in WormNet
        WormNetChV11 <- read.delim("Trans_Band_ChV11", header = FALSE)        
        head(WormNetChV11)
       
        ###Filter ook op transband ipv alleen chromosome V    
        outputChV11 <- filter(output, trans_band == "chrV:11.0-12.5Mb", Sequence_name %in% WormNetChV11$V1) 
    
        ###Filter op p-value 0.05
        WormNetChV11_p_value <- dplyr::filter(WormNetChV11,V5<0.20)
        
        outputChV11_p_value <-  dplyr::filter(outputChV11, Sequence_name %in% dplyr::filter(WormNetChV11,V5<0.20)$V1)
        outputChV11_p_value <-  merge(outputChV11, WormNetChV11_p_value,by.x=5,by.y=1) %>%
                                dplyr::filter(!duplicated(Sequence_name))
               
###transband ChV5.5
        
        WormNetChV5.5 <- read.delim("Trans_Band_ChV5.5", header = FALSE)        
        head(WormNetChV5.5)    
       
        ###Filter ook op transband ipv alleen chromosome V    
        outputChV5.5 <- filter(output, trans_band == "chrV:5.5-6.0Mb", Sequence_name %in% WormNetChV5.5$V1) 

        ###Filter op p-value 0.05
        WormNetChV5.5_p_value <- dplyr::filter(WormNetChV5.5,V5<0.20)
        
        outputChV5.5_p_value <-  dplyr::filter(outputChV5.5, Sequence_name %in% dplyr::filter(WormNetChV5.5,V5<0.20)$V1)
        outputChV5.5_p_value <-  merge(outputChV5.5, WormNetChV5.5_p_value,by.x=5,by.y=1)
        
        
###transband ChV14.5 moet nog want die is te groot
    
        WormNetChV14.5 <- read.delim("Trans_Band_ChV14.5", header = FALSE)        
        head(WormNetChV14.5)
      
        ###Filter ook op transband ipv alleen chromosome V    
        outputChV14.5 <- filter(output, trans_band == "chrV:14.5-17.0Mb", Sequence_name %in% WormNetChV14.5$V1) 
        
        ###Filter op p-value 0.20
        WormNetChV14.5_p_value <- dplyr::filter(WormNetChV14.5,V5<0.20)
        
        outputChV14.5_p_value <-  dplyr::filter(outputChV14.5, Sequence_name %in% dplyr::filter(WormNetChV14.5,V5<0.20)$V1)
        outputChV14.5_p_value <-  merge(outputChV14.5, WormNetChV14.5_p_value,by.x=5,by.y=1) %>%
                                  dplyr::filter(!duplicated(Sequence_name))
        
###transband ChX14
        
        WormNetChX14 <- read.delim("Trans_Band_ChX14", header = FALSE)        
        head(WormNetChX14)
        
        ###Filter ook op transband ipv alleen chromosome X    
        outputChX14 <- filter(output, trans_band == "chrX:14.0-14.5Mb", Sequence_name %in% WormNetChX14$V1) 

        ###Filter op p-value 0.20
        WormNetChX14_p_value <- dplyr::filter(WormNetChX14,V5<0.20)
        
        outputChX14_p_value <-  dplyr::filter(outputChX14, Sequence_name %in% dplyr::filter(WormNetChX14,V5<0.20)$V1)
        outputChX14_p_value <-  merge(outputChX14, WormNetChX14_p_value,by.x=5,by.y=1)
        
        
###transband ChX8.5
        
        WormNetChX8.5 <- read.delim("Trans_Band_ChX8.5", header = FALSE)        
        head(WormNetChX8.5)
        
        ###Filter ook op transband ipv alleen chromosome X    
        outputChX8.5 <- filter(output, trans_band == "chrX:8.5-9.0Mb", Sequence_name %in% WormNetChX8.5$V1) 
        
        ###Filter op p-value 0.20
        WormNetChX8.5_p_value <- dplyr::filter(WormNetChX8.5,V5<0.20)
        
        outputChX8.5_p_value <-  dplyr::filter(outputChX8.5, Sequence_name %in% dplyr::filter(WormNetChX8.5,V5<0.20)$V1)
        outputChX8.5_p_value <-  merge(outputChX8.5, WormNetChX8.5_p_value,by.x=5,by.y=1)

        
        
###transband ChX6.0
        
        WormNetChX6.0 <- read.delim("Trans_Band_ChX6.0", header = FALSE)        
        
        ###Filter ook op transband ipv alleen chromosome X    
        outputChX6.0 <- filter(output, trans_band == "chrX:6.0-6.5Mb", Sequence_name %in% WormNetChX6.0$V1) 
        
        ###Filter op p-value 0.20
        WormNetChX6.0_p_value <- dplyr::filter(WormNetChX6.0,V5<0.20)
        
        outputChX6.0_p_value <-  dplyr::filter(outputChX6.0, Sequence_name %in% dplyr::filter(WormNetChX6.0,V5<0.20)$V1)
        outputChX6.0_p_value <-  merge(outputChX6.0, WormNetChX6.0_p_value,by.x=5,by.y=1)%>%       
                                 dplyr::filter(!duplicated(Sequence_name))
        
        
#################Write excel bestand van je objects om een tabel te krijgen van alle outputCh.._p_value
        
        write.xlsx(outputChI_p_value, outputChII_p_value, outputChIV_p_value, outputChV11_p_value, outputChV5.5_p_value, outputChV14.5_p_value, outputChX14_p_value, outputChX8.5_p_value, outputChX6.0_p_value, file = "Table_filtered_trans_bands.xlsx", asTable = TRUE)  
        
        write.xlsx(outputChII_p_value, file = "Table_filtered_trans_bandI.xlsx", asTable = TRUE)  
        
        write.xlsx(outputChIV_p_value, file = "Table_filtered_trans_bandIV.xlsx", asTable = TRUE)  
        
        write.xlsx(outputChV11_p_value, file = "Table_filtered_trans_bandV11.xlsx", asTable = TRUE)  
        
        write.xlsx(outputChV5.5_p_value, file = "Table_filtered_trans_bandV5.5.xlsx", asTable = TRUE)  
        
        write.xlsx(outputChV14.5_p_value, file = "Table_filtered_trans_bandV14.5.xlsx", asTable = TRUE)  
        
        write.xlsx(outputChX14_p_value, file = "Table_filtered_trans_bandX14.xlsx", asTable = TRUE)  
        
        write.xlsx(outputChX8.5_p_value, file = "Table_filtered_trans_bandX8.5.xlsx", asTable = TRUE)  
        
        write.xlsx(outputChX6.0_p_value, file = "Table_filtered_trans_bandX6.0.xlsx", asTable = TRUE)  
        
        
        
        
        
        
        
############ Inladen WormExp excel sheet en genen koppelen
   
        ####Chromosome 1
        
             
        read.xlsx("WormExp trans-bands.xlsx", sheet = 1)
        
        wormexp1 <- read.xlsx("WormExp trans-bands.xlsx", sheet = 1)
        
        genes <- dplyr::select(Agilent.Ce.V2, gene_WBID, gene_sequence_name, gene_public_name)%>% 
                 dplyr::filter(!is.na(gene_WBID), gene_WBID!="", !duplicated(gene_WBID))
        
        grepl("tut-1", wormexp1[,2])
        
        match.nu <- function(x,y){which(grepl(x, y))}
        tapply(genes[,3], genes[,1], match.nu, y=wormexp1[,2])
        selectedgenesChI <- tapply(genes[,3], genes[,1], match.nu, y=wormexp1[,2])

        length(selectedgenesChI[[1]])        

        unlist(lapply(selectedgenesChI, length))
        
        unlist(lapply(selectedgenesChI, length))[unlist(lapply(selectedgenesChI, length))!=0]
        
        selectedgenesChI <- lapply(selectedgenesChI, paste, collapse = ";")
        
        selectedgenesChI <- data.frame(cbind(gene_WBID = names(selectedgenesChI), gene_match = unlist(selectedgenesChI)))
        
        head(selectedgenesChI)
    
        selectedgenesChI <- filter(selectedgenesChI, gene_match !="") %>%
                         merge(genesChI, by.x = 1, by.y = 1)
        
        
            
            #####genes polymorphisms and in trans-band region
       
            selectedgenesChIfiltered <- filter(output, trans_band == "chrI:12.5-13.0Mb", Sequence_name %in% selectedgenesChI$gene_sequence_name)%>%
                                        merge(output, selectedgenesChI, by.x = 5, by.y = 3)
                           
                           
        ####Chromosome 2
        
        read.xlsx("WormExp trans-bands.xlsx", sheet = 2)
            
        wormexp2 <- read.xlsx("WormExp trans-bands.xlsx", sheet = 2)
            
        ###grepl("tut-1", wormexp2[,2])
            
        match.nu <- function(x,y){which(grepl(x, y))}
        tapply(genes[,3], genes[,1], match.nu, y=wormexp2[,2])
        selectedgenesChII <- tapply(genes[,3], genes[,1], match.nu, y=wormexp2[,2])
            
        length(selectedgenesChII[[1]])        
            
        unlist(lapply(selectedgenesChII, length))
            
        unlist(lapply(selectedgenesChII, length))[unlist(lapply(selectedgenesChII, length))!=0]
            
        selectedgenesChII <- lapply(selectedgenesChII, paste, collapse = ";")
            
        selectedgenesChII <- data.frame(cbind(gene_WBID = names(selectedgenesChII), gene_match = unlist(selectedgenesChII)))
            
        head(selectedgenesChII)
            
        selectedgenesChII <- filter(selectedgenesChII, gene_match !="") %>%
        merge(genes, by.x = 1, by.y = 1)
            
            
            
                  #####genes polymorphisms and in trans-band region
            
                  selectedgenesChIIfiltered <- filter(output, trans_band == "chrII:6.0-6.5Mb", Sequence_name %in% selectedgenesChII$gene_sequence_name)%>%
                  merge(output, selectedgenesChII, by.x = 5, by.y = 3)
            
  
            
            ####Chromosome 4     
          
            
        read.xlsx("WormExp trans-bands.xlsx", sheet = 3)
        
        wormexp4 <- read.xlsx("WormExp trans-bands.xlsx", sheet = 3)
        
        ###grepl("tut-5", wormexp2[,2])
        
        
        match.nu <- function(x,y){which(grepl(x, y))}
        tapply(genes[,3], genes[,1], match.nu, y=wormexp4[,2])
        selectedgenesChIV <- tapply(genes[,3], genes[,1], match.nu, y=wormexp4[,2])
        
        length(selectedgenesChIV[[1]])        
        
        unlist(lapply(selectedgenesChIV, length))
        
        unlist(lapply(selectedgenesChIV, length))[unlist(lapply(selectedgenesChIV, length))!=0]
        
        selectedgenesChIV <- lapply(selectedgenesChIV, paste, collapse = ";")
        
        selectedgenesChIV <- data.frame(cbind(gene_WBID = names(selectedgenesChIV), gene_match = unlist(selectedgenesChIV)))
        
        head(selectedgenesChIV)
        
        selectedgenesChIV <- filter(selectedgenesChIV, gene_match !="") %>%
          merge(genes, by.x = 1, by.y = 1)
        
        #####genes polymorphisms and in trans-band region
        
        selectedgenesChIVfiltered <- filter(output, trans_band == "chrIV:14.5-15.0Mb", Sequence_name %in% selectedgenesChIV$gene_sequence_name)%>%
          merge(output, selectedgenesChIV, by.x = 5, by.y = 3)
        
        
        
    ####Chromosome V 11.0-12.5Mb    

        
        read.xlsx("WormExp trans-bands.xlsx", sheet = 4)
        
        wormexp5_11 <- read.xlsx("WormExp trans-bands.xlsx", sheet = 4)
        
        ###grepl("tut-5", wormexp5_11[,2])
        
        
        match.nu <- function(x,y){which(grepl(x, y))}
        tapply(genes[,3], genes[,1], match.nu, y=wormexp5_11[,2])
        selectedgenesChV11 <- tapply(genes[,3], genes[,1], match.nu, y=wormexp5_11[,2])
        
        length(selectedgenesChV11[[1]])        
        
        unlist(lapply(selectedgenesChV11, length))
        
        unlist(lapply(selectedgenesChV11, length))[unlist(lapply(selectedgenesChV11, length))!=0]
        
        selectedgenesChV11 <- lapply(selectedgenesChV11, paste, collapse = ";")
        
        selectedgenesChV11 <- data.frame(cbind(gene_WBID = names(selectedgenesChV11), gene_match = unlist(selectedgenesChV11)))
        
        head(selectedgenesChV11)
        
        selectedgenesChV11 <- filter(selectedgenesChV11, gene_match !="") %>%
          merge(genes, by.x = 1, by.y = 1)
        
        #####genes polymorphisms and in trans-band region
        
        selectedgenesChV11filtered <- filter(output, trans_band == "chrV:11.0-12.5Mb", Sequence_name %in% selectedgenesChV11$gene_sequence_name)%>%
          merge(output, selectedgenesChV11, by.x = 5, by.y = 3)
      
        
        
  ####Chromosome V 14.5-17.0Mb    
        
        
        read.xlsx("WormExp trans-bands.xlsx", sheet = 5)
        
        wormexp5_14 <- read.xlsx("WormExp trans-bands.xlsx", sheet = 5)
        
        ###grepl("tut-5", wormexp5_14[,2])
        
        
        match.nu <- function(x,y){which(grepl(x, y))}
        tapply(genes[,3], genes[,1], match.nu, y=wormexp5_14[,2])
        selectedgenesChV14 <- tapply(genes[,3], genes[,1], match.nu, y=wormexp5_14[,2])
        
        length(selectedgenesChV14[[1]])        
        
        unlist(lapply(selectedgenesChV14, length))
        
        unlist(lapply(selectedgenesChV14, length))[unlist(lapply(selectedgenesChV14, length))!=0]
        
        selectedgenesChV14 <- lapply(selectedgenesChV14, paste, collapse = ";")
        
        selectedgenesChV14 <- data.frame(cbind(gene_WBID = names(selectedgenesChV14), gene_match = unlist(selectedgenesChV14)))
        
        head(selectedgenesChV14)
        
        selectedgenesChV14 <- filter(selectedgenesChV14, gene_match !="") %>%
          merge(genes, by.x = 1, by.y = 1)
        
        #####genes polymorphisms and in trans-band region
        
        selectedgenesChV14filtered <- filter(output, trans_band == "chrV:14.5-17.0Mb", Sequence_name %in% selectedgenesChV14$gene_sequence_name)%>%
          merge(output, selectedgenesChV14, by.x = 5, by.y = 3)
        
        
  
        
  ####Chromosome V 5.5-6.0Mb    
        
        
        read.xlsx("WormExp trans-bands.xlsx", sheet = 6)
        
        wormexp5_5 <- read.xlsx("WormExp trans-bands.xlsx", sheet = 6)
        
        ###grepl("tut-5", wormexp5_5[,2])
        
        
        match.nu <- function(x,y){which(grepl(x, y))}
        tapply(genes[,3], genes[,1], match.nu, y=wormexp5_5[,2])
        selectedgenesChV5 <- tapply(genes[,3], genes[,1], match.nu, y=wormexp5_5[,2])
        
        length(selectedgenesChV5[[1]])        
        
        unlist(lapply(selectedgenesChV5, length))
        
        unlist(lapply(selectedgenesChV5, length))[unlist(lapply(selectedgenesChV5, length))!=0]
        
        selectedgenesChV5 <- lapply(selectedgenesChV5, paste, collapse = ";")
        
        selectedgenesChV5 <- data.frame(cbind(gene_WBID = names(selectedgenesChV5), gene_match = unlist(selectedgenesChV5)))
        
        head(selectedgenesChV5)
        
        selectedgenesChV5 <- filter(selectedgenesChV5, gene_match !="") %>%
          merge(genes, by.x = 1, by.y = 1)
        
        #####genes polymorphisms and in trans-band region
        
        selectedgenesChV5filtered <- filter(output, trans_band == "chrV:5.5-6.0Mb", Sequence_name %in% selectedgenesChV5$gene_sequence_name)%>%
          merge(output, selectedgenesChV5, by.x = 5, by.y = 3)
        
  
        
  ####Chromosome X 14.0-14.5Mb    
        
        
        read.xlsx("WormExp trans-bands.xlsx", sheet = 7)
        
        wormexpX_14 <- read.xlsx("WormExp trans-bands.xlsx", sheet = 7)
        
        grepl("tut-5", wormexpX_14[,2])
        
        
        match.nu <- function(x,y){which(grepl(x, y))}
        tapply(genes[,3], genes[,1], match.nu, y=wormexpX_14[,2])
        selectedgenesChX14 <- tapply(genes[,3], genes[,1], match.nu, y=wormexpX_14[,2])
        
        length(selectedgenesChX14[[1]])        
        
        unlist(lapply(selectedgenesChX14, length))
        
        unlist(lapply(selectedgenesChX14, length))[unlist(lapply(selectedgenesChX14, length))!=0]
        
        selectedgenesChX14 <- lapply(selectedgenesChX14, paste, collapse = ";")
        
        selectedgenesChX14 <- data.frame(cbind(gene_WBID = names(selectedgenesChX14), gene_match = unlist(selectedgenesChX14)))
        
        head(selectedgenesChX14)
        
        selectedgenesChX14 <- filter(selectedgenesChX14, gene_match !="") %>%
                              merge(genes, by.x = 1, by.y = 1)
        
        #####genes polymorphisms and in trans-band region
        
        selectedgenesChX14filtered <- filter(output, trans_band == "chrX:14.0-14.5Mb", Sequence_name %in% selectedgenesChX14$gene_sequence_name)%>%
          merge(output, selectedgenesChX14, by.x = 5, by.y = 3)

        
        
       
####Chromosome X 6.0-6.5Mb    
        
        
        read.xlsx("WormExp trans-bands.xlsx", sheet = 8)
        
        wormexpX_6 <- read.xlsx("WormExp trans-bands.xlsx", sheet = 8)
        
        grepl("tut-5", wormexpX_6[,2])
        
        
        match.nu <- function(x,y){which(grepl(x, y))}
        tapply(genes[,3], genes[,1], match.nu, y=wormexpX_6[,2])
        selectedgenesChX6 <- tapply(genes[,3], genes[,1], match.nu, y=wormexpX_6[,2])
        
        length(selectedgenesChX6[[1]])        
        
        unlist(lapply(selectedgenesChX6, length))
        
        unlist(lapply(selectedgenesChX6, length))[unlist(lapply(selectedgenesChX6, length))!=0]
        
        selectedgenesChX6 <- lapply(selectedgenesChX6, paste, collapse = ";")
        
        selectedgenesChX6 <- data.frame(cbind(gene_WBID = names(selectedgenesChX6), gene_match = unlist(selectedgenesChX6)))
        
        head(selectedgenesChX6)
        
        selectedgenesChX6 <- filter(selectedgenesChX6, gene_match !="") %>%
          merge(genes, by.x = 1, by.y = 1)
        
        #####genes polymorphisms and in trans-band region
        
        selectedgenesChX6filtered <- filter(output, trans_band == "chrX:6.0-6.5Mb", Sequence_name %in% selectedgenesChX6$gene_sequence_name)%>%
          merge(output, selectedgenesChX6, by.x = 5, by.y = 3) 
        
        
        
####Chromosome X 8.5-9.0Mb    
        
        
        read.xlsx("WormExp trans-bands.xlsx", sheet = 9)
        
        wormexpX_8 <- read.xlsx("WormExp trans-bands.xlsx", sheet = 9)
        
        grepl("tut-5", wormexpX_8[,2])
        
        
        match.nu <- function(x,y){which(grepl(x, y))}
        tapply(genes[,3], genes[,1], match.nu, y=wormexpX_8[,2])
        selectedgenesChX8 <- tapply(genes[,3], genes[,1], match.nu, y=wormexpX_8[,2])
        
        length(selectedgenesChX8[[1]])        
        
        unlist(lapply(selectedgenesChX8, length))
        
        unlist(lapply(selectedgenesChX8, length))[unlist(lapply(selectedgenesChX8, length))!=0]
        
        selectedgenesChX8 <- lapply(selectedgenesChX8, paste, collapse = ";")
        
        selectedgenesChX8 <- data.frame(cbind(gene_WBID = names(selectedgenesChX8), gene_match = unlist(selectedgenesChX8)))
        
        head(selectedgenesChX8)
        
        selectedgenesChX8 <- filter(selectedgenesChX8, gene_match !="") %>%
          merge(genes, by.x = 1, by.y = 1)
        
        #####genes polymorphisms and in trans-band region
        
        selectedgenesChX8filtered <- filter(output, trans_band == "chrX:8.5-9.0Mb", Sequence_name %in% selectedgenesChX8$gene_sequence_name)%>%
          merge(output, selectedgenesChX8, by.x = 5, by.y = 3)   
        