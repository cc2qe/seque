source('~/code/seque/drawIdeogram/src/drawIdeogram.R')

# define the plot area
plot(c(0, getChromLength(4, genome='mm9')), c(-2,2), type='n')

# draw the ideogram for mouse chromosome 4
drawIdeogram(4, genome="mm9")