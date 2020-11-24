normal.model <- function( x, coeffs, which = NULL ) {
    mu = coeffs[1]
    variance = coeffs[2]
    alpha = coeffs[3:length(coeffs)]
    v = variance * 1:length(alpha)
    result = 0
    if( is.null( which )) {
        which = 1:length(alpha)
    }
    for( i in which ) {
        result = result + alpha[i] * dnorm( x, mean = i * mu, sd = sqrt( v[i] ))
    }
    return(result)
}

t.model <- function( x, coeffs, t.df = 3, which = NULL ) {
    mu = coeffs[1]
    variance = coeffs[2]
    df = t.df
    alpha = coeffs[3:length(coeffs)]
    v = variance * 1:length(alpha)
    result = 0
    if( is.null( which )) {
        which = 1:length(alpha)
    }
    for( i in which ) {
        result = result + alpha[i] * dt( (x - i*mu)/sqrt(v[i]), df = df )
    }
    return(result)
}

# From genomescope R script:
genomescope1.model <- function( x, coeffs, k = 31 ) {
    d = coeffs[1]
    r = coeffs[2]
    kmercov = coeffs[3]
    bias = coeffs[4]
    length = coeffs[5]
    result = (((2*(1-d)*(1-(1-r)^k)) + (2*d*(1-(1-r)^k)^2) + (2*d*((1-r)^k)*(1-(1-r)^k))) * dnbinom(x, size = kmercov   / bias, mu = kmercov)     * length +
                              (((1-d)*((1-r)^k)) + (d*(1-(1-r)^k)^2))                                      * dnbinom(x, size = kmercov*2 / bias, mu = kmercov * 2) * length + 
                              (2*d*((1-r)^k)*(1-(1-r)^k))                                                  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov * 3) * length + 
                              (d*(1-r)^(2*k))                                                              * dnbinom(x, size = kmercov*4 / bias, mu = kmercov * 4) * length)
    return( result )
}

plot.it <- function( X, model, coeffs, which = NULL, xlim = c( 0, nrow(X)), ylim = c( 0, 5E7) ) {
    plot( X[,1], X[,2], type = 'l', xlim = xlim, ylim = ylim )
    points(
        X[,1],
        model(X[,1], coeffs ),
        type =  'l',
        col = 'red', 
        lty = 1
    )
    if( !is.null( which )) {
        for( i in which ) {
            points(
                X[,1],
                model(X[,1], coeffs, i ),
                type =  'l',
                col = 'grey', 
                lty = 2
            )
        }
    }
}

mean.residual.squared.error <- function( fit ) {
    mean( summary( fit )$residuals^2 )
}

plot.fits  <-function( fits, filename ) {
    pdf( file = filename, width = 6, height = 10 )
    layout(
        matrix( c( 1, 1, 1:(length(fits)*2)+1), ncol = 2, byrow = T ),
        heights = c( 0.1, rep( 1, length(fits)))
    )
    par( mar = c( 0.1, 1, 0.2, 1 ))
    plot( 0, 0, col = 'white', xlab = '', ylab = '', bty = 'n', xaxt = 'n', yaxt = 'n', xlim = c( 0, 1), ylim = c( 0, 1 ) )
    text( 0.5, 0.5, "Parametric fits to kmer distribution", cex = 2, adj = c( 0.5, 0.5 ), xpd = NA )
    for( i in 1:length(fits)) {
        par( mar = c( 2, 2, 1, 1 ))
        plot.it( X, fits[[i]]$model, fits[[i]]$coeffs[,1], xlim = c( 0, mu * 4 ), ylim = c( 0, X[2*mu,2] * 1.1 ) )
        grid()
        legend(
            "topright",
            legend = c(
                sprintf( "%s", names(fits)[i] )
            ),
            cex = 2,
            bty = 'n'
        )
        legend(
            "topleft",
            c( "empirical", "fit" ),
            lty = 1,
            col = c( "black", "red" ),
            bty = 'n'
        )
        par( mar = c( 2, 2, 1, 1 ))
        plot.it( X, fits[[i]]$model, fits[[i]]$coeffs[,1], xlim = c( 0, mu * 2 ), ylim = c( 0, X[2*mu,2] * 0.25 ) )
    
        values = fits[[i]]$model( X[,1], fits[[i]]$coeffs[,1] )
        residual = X[,2] - values
    
        points( X[1:80,1], residual[1:80], type = 'l', lty = 3, col = 'gold2' )
        abline( v = min( which( residual[1:80] <= values[1:80] )), col = 'gold2', lty = 1 )

        # P(error|cov=x) = P(cov=x|error) P(error) / ( P(cov=x|error)P(error) + P(cov=x|!error) P(!error) )
        residual[ 60:length(residual)] = 0
        error_distribution =  residual / sum( residual )
        nonerror_distribution = values / sum( values )

        error.prior = 1E-3
        posteriors1 = (error_distribution * error.prior) / ( error_distribution * error.prior+ nonerror_distribution * (1-error.prior ))
        error.prior = 1E-4
        posteriors2 = (error_distribution * error.prior) / ( error_distribution * error.prior+ nonerror_distribution * (1-error.prior ))
        min( which( residual <= values ))

        legend(
            "topright",
            c(
                sprintf( "50/50 split point = %d", min( which( residual <= values )) ),
                sprintf( "...w/1 in 1000 prior = %d", min( which( posteriors1 <= 0.5 )) ),
                sprintf( "...w/1 in 10000 prior = %d", min( which( posteriors2 <= 0.5 )) ),
                sprintf( "mrse = %.2g", mean.residual.squared.error( fits[[i]] ))
            ),
            bty = 'n'
        )
        grid()
    }
    dev.off()
}


#################################
# Load data and fit distribution.

fit.models <- function( X ) {
    mu = (which.max( X[10:nrow(X),2] )+9)/2

    fits = list()
    fits[['normal']] = nls(
        count ~ normal.model( cov, coeffs ),
        start = list(
            coeffs = c(
                mu = mu,
                variance = mu,
                alpha = c( 1E8, 5E8, 1E5, 1E7, 1E5, 1E6 )
            )
        ),
        data = X[10:1000,],
        trace = TRUE,
        control = list(
            minFactor = 1E-12,
            maxiter = 10000
        )
    )
    fits[['normal']]$model = normal.model
    fits[['normal']]$coeffs = summary(fits[['normal']])$coeff

    for( df in c( 5, 10, 20, 30, 50, 100 )) {
        name = sprintf( "t.df=%d", df )
        fits[[name]] = nls(
            count ~ t.model( cov, coeffs, t.df = df ),
            start = list(
                coeffs = c(
                    mu = mu,
                    variance = mu,
                    alpha = c( 1E8, 5E8, 1E5, 1E7, 1E5, 1E6 )
                )
            ),
            data = X[10:1000,],
            trace = TRUE,
            control = list(
                minFactor = 1E-12,
                maxiter = 10000
            )
        )
        fits[[name]]$model = function( cov, coeffs ) { return( t.model( cov, coeffs, t.df = df ) ) }
        fits[[name]]$coeffs = summary(fits[[name]])$coeff
    }

    # This tends to fail:
    fits[['genomescope.1']] = nls(
        count ~ genomescope1.model( cov, coeffs ),
        start = list(
            coeffs = c(
                d = 0,
                r = 0,
                kmercov = (which.max( X[10:nrow(X),2] )+9)/2,
                bias = 0.5,
                length = 3.2E9
            )
        ),
        data = X[10:1000,],
        trace = TRUE,
        control = list(
            minFactor = 1E-12,
            maxiter = 100,
            warnOnly = TRUE
        )
    )
    fits[['genomescope.1']]$model = genomescope1.model
    fits[['genomescope.1']]$coeffs = summary(fits[['genomescope.1']])$coeff
    return( fits )
}

# X = read.table( "/path/to/histogram", head = F, as.is = T )
# colnames(X) = c( "cov", "count" )

fits = fit.models( X )

plot.fits( fits, "/tmp/kmer_distributions.pdf" )
