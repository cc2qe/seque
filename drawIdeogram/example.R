source('~/code/seque/drawIdeogram/src/drawIdeogram.R')

## Human chromosome example
# define the plot area
plot(c(0, getChromLength(10, genome='hg19')), c(-2,2), type='n')
# draw the ideogram for human chromosome 16
drawIdeogram(10, genome="hg19")


## Mouse chromosome example
# define the plot area
plot(c(0, getChromLength(4, genome='mm9')), c(-2,2), type='n')
# draw the ideogram for mouse chromosome 4
drawIdeogram(4, genome="mm9")