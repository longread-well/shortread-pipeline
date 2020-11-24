library( wesanderson )

samples = c( "NA12878", "HV31" )
coverage = "30X"
technology = c( "coolmps", "illumina" )
k = 31

data = list()
for( sample in samples ) {
   data[[sample]] = list()
   for( coverage in c( "30X", "full" )) {
      print( coverage)
      tryCatch({
          data[[sample]][[sprintf( "whg_illumina_%s", coverage )]] = read.table( sprintf( "results/kmer/histograms/illumina/whg/whg_illumina_%s_%s.k=%d.bq=any.histogram.gz", sample, coverage, k), header = F, as.is=T )
      }, error = function(e){} )
      tryCatch({
          data[[sample]][[sprintf( "nygc_illumina_%s", coverage )]] = read.table( sprintf( "results/kmer/histograms/illumina/nygc/nygc_illumina_%s_150bp_%s.k=%d.bq=any.histogram.gz", sample, coverage, k), header = F, as.is=T )
      }, error = function(e){} )
      tryCatch({
         data[[sample]][[sprintf( "whg_coolmps_%s", coverage )]] = read.table( sprintf( "results/kmer/histograms/coolmps/whg/whg_coolmps_%s_%s.k=%d.bq=any.histogram.gz", sample, coverage, k), header = F, as.is=T )
      }, error = function(e){} )
   }
}

blank.plot <- function( xlim = c( 0, 1 ), ylim = c( 0, 1 ), xlab = '', ylab = '' ) {
   plot(
       0, 0,
       col = 'white', bty = 'n', xaxt = 'n', yaxt = 'n',
       xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab
   )
}

panel1 <- function( data, sample, whats = c( 'whg_illumina_30X', 'nygc_illumina_30X', 'whg_coolmps_30X' ), k = 31 ) {

    palette = wes_palette( "Zissou1" )
    colours = palette[c(1,4,5)]
    names(colours) = whats

    par( mar = c( 4.1, 4.1, 2.1, 1.1 ))
    blank.plot(
        xlim = c( -2, 35 ),
        ylim = c( 0, 0.09 ),
    #    ylim = c( 0, 5E9 ),
        xlab = sprintf( "k-mer coverage (k=%d)", k ),
        ylab = "proportion of read k-mers"
    )
    estimates = data.frame()
    for( i in 1:length( whats )) {
       what = whats[i]
       X = data[[sample]][[what]]
       denominator = sum(as.numeric(X[,2])*X[,1])
       
       
       points(
           X[,1],
           as.numeric(X[,2])*X[,1] / denominator,
           type = 'l',
           col = colours[i]
       )
       points(
          X[1,1],
          X[1,2]*X[1,1] / denominator,
          pch = 19,
          col = colours[i]
       )
       polygon(
          c( 1, X[1:5,1] ),
          c( 0, X[1:5,1] * X[1:5,2] ) / denominator,
          col = colours[i],
          border = NA
       )
       errors = data.frame(
          name = what,
          total = denominator,
          error = sum( X[1:5,1] * X[1:5,2] )
       )
       errors$rate = errors$error / errors$total
       estimates = rbind( estimates, errors )
       grid()
       text(
          0.5,
          X[1,2]*X[1,1] / denominator - 0.01,
          adj = 1,
          label = sprintf( "%.2f%%", errors$rate * 100 / k ),
          col = colours[i]
       )
       segments(
          x0 = 0.65,
          x1 = 0.85,
          y0 = X[1,2]*X[1,1] / denominator - 0.01,
          y1 = X[1,2]*X[1,1] / denominator - 0.01,
          col = colours[i]
       )
    }
    abline( v = 5, col = 'black', lty = 2 )
    axis(1)
    axis(2, at = c( 0, 0.02, 0.04, 0.06, 0.08, 0.1 ), labels = c( "0", "2%", "4%", "6%", "8%", "10%" ), las = 1 )

    legend.names = c(
       "whg_illumina_30X" = "Illumina Novaseq (WHG) 150bp @ 30X",
       "nygc_illumina_30X" = "Illumina Novaseq (NYGC) 151bp @ 30X",
       "whg_coolmps_30X" = "MGI CoolMPS (WHG) 100bp @ 30X"
    )

    legend(
       "top",
       legend = legend.names[ estimates$name ],
       col = colours[whats],
       lty = 1,
       pch = 19,
       bty = 'n'
    )
}

panel2 <- function( data, sample, whats = c( 'whg_illumina_full', 'whg_coolmps_full' ), k = 31 ) {

    palette = wes_palette( "Zissou1" )
    colours = palette[c(1,5)]
    names(colours) = whats

    par( mar = c( 4.1, 4.1, 2.1, 1.1 ))
    blank.plot(
        xlim = c( 0, 100 ),
        ylim = c( 0, 1.5E8 ),
    #    ylim = c( 0, 5E9 ),
        xlab = sprintf( "k-mer coverage (k=%d)", k ),
        ylab = "Number of distinct k-mers"
    )
    estimates = data.frame()
    for( i in 1:length( whats )) {
       what = whats[i]
       X = data[[sample]][[what]]
       denominator = sum(as.numeric(X[,2])*X[,1])
       points(
           X[,1],
           as.numeric(X[,2]),
           type = 'l',
           col = colours[i]
       )
    }
    grid()
    axis(1)
    axis(2, at = c( 0, 0.02, 0.04, 0.06, 0.08, 0.1 ), labels = c( "0", "2%", "4%", "6%", "8%", "10%" ), las = 1 )
}

pdf( file = "paper/figures/Figure_01_kmers.pdf", width = 8, height = 7 )
par( mfrow = c( 2, 1 ))
panel1( data, "NA12878", whats = c( 'whg_illumina_30X', 'nygc_illumina_30X', 'whg_coolmps_30X' ) )
panel2( data, "NA12878", whats = c( 'whg_illumina_full', 'whg_coolmps_full' ) )
dev.off()

