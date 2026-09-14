library(circlize)
library(stringr)

integrationBed <- read.delim("/Users/guwenjin/Documents/manuscript/figure_3/integrationBed_Sup_3_type2.txt", header=FALSE)
colnames(integrationBed) = c('chr','start','end','value1')
integrationBed$value1 = -integrationBed$value1



circos.initializeWithIdeogram()
bed = generateRandomBed(nr = 200)
col_fun = colorRamp2(breaks = c(-5,0), colors = c("red", "grey"))
circos.genomicTrack(integrationBed, bg.border = NA, track.height = 0.8,
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom  = 0,
                                         col = col_fun(value[[1]]), border = '#ffffff00' ,...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lwd = 1, lty = 2, col = "grey")
                      circos.lines(CELL_META$cell.xlim, c(-5, -5), lwd = 1, lty = 2, col = "grey")
                      circos.lines(CELL_META$cell.xlim, c(-10, -10), lwd = 1, lty = 2, col = "grey")
                      circos.lines(CELL_META$cell.xlim, c(-15, -15), lwd = 1, lty = 2, col = "grey")
                      circos.lines(CELL_META$cell.xlim, c(-20, -20), lwd = 1, lty = 2, col = "grey")
                      circos.lines(CELL_META$cell.xlim, c(-25, -25), lwd = 1, lty = 2, col = "grey")
                      
                    })

circos.clear()
