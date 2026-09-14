library(circlize)
library(stringr)
integrationBed <- read.delim("/Users/wenjingu/Library/CloudStorage/SynologyDrive-macbook13/Google_Drive_Backup/research/hpv fusion/project/newSampleCircos/integrationBed.txt", header=FALSE)
colnames(integrationBed) = c('chr','start','end','value1')
integrationBed$value1 = -integrationBed$value1



circos.initializeWithIdeogram()
bed = generateRandomBed(nr = 200)
col_fun = colorRamp2(breaks = c(-5,0), colors = c("red", "grey"))
circos.genomicTrack(integrationBed, bg.border = NA, track.height = 0.5,
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom  = 0,
                                         col = col_fun(value[[1]]), border = '#ffffff00' ,...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lwd = 1, lty = 2, col = "grey")
                      circos.lines(CELL_META$cell.xlim, c(-2, -2), lwd = 1, lty = 2, col = "grey")
                      circos.lines(CELL_META$cell.xlim, c(-4, -4), lwd = 1, lty = 2, col = "grey")
                      circos.lines(CELL_META$cell.xlim, c(-6, -6), lwd = 1, lty = 2, col = "grey")
                      circos.lines(CELL_META$cell.xlim, c(-8, -8), lwd = 1, lty = 2, col = "grey")
                      circos.lines(CELL_META$cell.xlim, c(-10, -10), lwd = 1, lty = 2, col = "grey")
                    })

set.seed(123)
bed1 = generateRandomBed(nr = 100)
bed1 = bed1[sample(nrow(bed1), 20), ]
bed2 = generateRandomBed(nr = 100)
bed2 = bed2[sample(nrow(bed2), 20), ]

circos.clear()

track <- read.delim("/Volumes/GoogleDrive/My Drive/research/hpv fusion/project/newSampleCircos/track.bed",header = FALSE,stringsAsFactors = FALSE)
colnames(track) = c('gene','start','end')
hpvSite = read.delim("/Volumes/GoogleDrive/My Drive/research/hpv fusion/project/newSampleCircos/hpvIntegration_part3.bed",header = FALSE,stringsAsFactors = FALSE)
geneSite = read.delim("/Volumes/GoogleDrive/My Drive/research/hpv fusion/project/newSampleCircos/genomeIntegration_part3.bed",header = FALSE,stringsAsFactors = FALSE)
colnames(hpvSite) = c('gene','start','end')
colnames(geneSite) = c('gene','start','end')
label = read.delim("/Volumes/GoogleDrive/My Drive/research/hpv fusion/project/newSampleCircos/label.bed",header = FALSE,stringsAsFactors = FALSE)
colnames(label) = c('chr','start','end','value1')

col1 =rand_color(13,transparency = 0.5)
col2 = rep("lightBlue",37-13)
col = c(col1,col2)
circos.genomicInitialize(track,labels.cex = 0.5,track.height = 0.1,plotType = NULL)
circos.track(ylim = c(0, 1), 
             bg.col = col, 
             bg.border = NA, track.height = 0.1, panel.fun = function(x, y) {
               chr = CELL_META$sector.index
               xlim = CELL_META$xlim
               ylim = CELL_META$ylim
               if(str_detect(chr,"chr") == FALSE){
                 label = paste("------",chr,sep = " ")
               }
               else{
                 label = chr
               }
               
               circos.text(mean(xlim),2,label , cex = 0.6, col = "Black",
                           facing = "clockwise", niceFacing = TRUE)
             })

hpvGenes = unique(track$gene[1:18])

linkColor = rep(0,nrow(hpvSite))
i = 1
for(each in hpvSite$gene){
  j = 1
  for(eachGene in hpvGenes){
    if(each == eachGene){
      linkColor[i] = col1[j]
    }
    j = j + 1
  }
  i = i + 1
}

circos.genomicLink(hpvSite,geneSite, col = linkColor, border = NA,lwd = 0.5)



