library(MASS) # smoothing
library(RColorBrewer)
library(hexbin) # honeycomb plot

pal <- brewer.pal(11,'Spectral')
rf <- colorRampPalette(rev(pal[3:11]))
r <- rf(32)

cov <- read.table('coverages_2.tsv')

minor_variant_rel_cov <- cov$V1 / (cov$V1 + cov$V2)
total_pair_cov <- cov$V1 + cov$V2

high_cov_filt <- quantile(total_pair_cov, 0.99) > total_pair_cov
minor_variant_rel_cov <- minor_variant_rel_cov[high_cov_filt]
total_pair_cov <- total_pair_cov[high_cov_filt]

h1 <- hist(minor_variant_rel_cov, breaks = 100, plot = F)
h2 <- hist(total_pair_cov, breaks = 100, plot = F)

top <- max(h1$counts, h2$counts)
k <- kde2d(minor_variant_rel_cov, total_pair_cov, n=30)
k$z <- sqrt(k$z)

png('~/Desktop/Avag1_kmer_pair_coverage.png')
# margins 'c(bottom, left, top, right)'
par(mar=c(4.8,4.8,1,1))
layout(matrix(c(2,4,1,3), 2, 2, byrow=T), c(3,1), c(1,3))
#plot the image
# 2D HISTOGRAM
image(k, col = r,
      xlab = 'Normalized minor kmer coverage',
      ylab = 'Total coverage of the kmer pair', cex.lab = 1.4
)
# TODO image ticks by the coverage counts 2n, 3n, 4n, 5n, 6n, 7n, 8n...

n <- 92
# for(i in 2:6){
#       lines(c(0, 0.6), rep(i * n, 2), lwd = 1.4)
#       text(0.1, i * n, paste0(i,'x'), pos = 3)
# }

# EXPECTED COMPOSITIONS
text(1/2 - 0.01, 2 * n, 'AB', offset = 0, cex = 1.4)
text(1/3, 3 * n, 'AAB', offset = 0, cex = 1.4)
text(1/4, 4 * n, 'AAAB', offset = 0, cex = 1.4)
text(1/2 - 0.022, 4 * n, 'AABB', offset = 0, cex = 1.4)
text(2/5, 5 * n, 'AAABB', offset = 0, cex = 1.4)
text(1/5, 5 * n, 'AAAAB', offset = 0, cex = 1.4)
text(3/6 - 0.035, 6 * n, 'AAABBB', offset = 0, cex = 1.4)
text(2/6, 6 * n, 'AAAABB', offset = 0, cex = 1.4)
text(1/6, 6 * n, 'AAAAAB', offset = 0, cex = 1.4)

# minor_variant_rel_cov HISTOGRAM
par(mar=c(0,3.8,1,0))
barplot(h1$counts, axes=F, ylim=c(0, top), space=0, col = pal[2])
# total pair coverage HISTOGRAM
par(mar=c(3.8,0,0.5,1))
barplot(h2$counts, axes=F, xlim=c(0, top), space=0, col = pal[2], horiz=T)

# LEGEND
par(mar=c(0,0,2,1))
plot.new()
title('Coverage (^2 scale)')
for(i in 1:32){
      rect(0,(i - 0.01) / 33, 0.5, (i + 0.99) / 33, col = r[i])
}
kmer_max <- (length(total_pair_cov) * max((k$z)^2)) / sum((k$z)^2)
for(i in 0:6){
      text(0.75, i / 6, round((sqrt(kmer_max) * i)^2 / 6000) * 1000, offset = 0)
}

dev.off()

# alternative plotting
# h <- hexbin(df)
# # h@count <- sqrt(h@count)
# plot(h, colramp=rf)
# gplot.hexbin(h, colorcut=10, colramp=rf)
