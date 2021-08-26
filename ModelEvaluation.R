###########################################################################
###   This script was developed by Dr. Jen Cruz to evaluate the model  ###
### described in Sauer et al. (in review ) Multiscale drivers of amphibian #
### occupancy in urban ponds. The model evaluation includes species-level #
##  deviance residuals and an overall Bayesian p value, as posterior  ###
##   model checks. The model is a single season, multispecies occupancy ##
###  model, fit under a Bayesian framework.                           ###
###########################################################################

########### clean workspace and load required packages ####################
###########################################################################

#clean workspace to improve efficiency:
rm(list = ls() ) 

####### load relevant packages ###
library( tidyverse ) #dataframe manipulations.
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( jagsUI ) #to run RJAGS

#######################    import relevant data   ##############

#load relevant workspace 
load( "OccResults.RData" )

################################################################################
#################### viewing model results ####################################
##################################################################################
#define model object
mr <- cm1

#Total number of iterations ran:
N <- mr$mcmc.info$n.samples
#I = number of ponds
#S = number of species
#labels for occupancy models
labs <- modcovs
#labels for detection model
detlabs <- c( "Floating vegetation", "Day of year") 

########################################################################
## Calculating deviance from likelihood of data conditional on model ###
########################################################################
# Define function to calculate Model deviances 
CalcDevs <- function(lik_yobs = mr$sims.list$lik_yobs,
                     lik_yhat = mr$sims.list$lik_yhat ){

  #assign temporary objects
  ll_yobs_is <- ll_yhat_is <- array( 0, dim = c(  N, I, S ) )
  ll_yobs_s <-  ll_yhat_s <-  matrix( nrow = N, ncol = S )
  ll_yobs <- ll_yhat <- Dev_obs <- Dev_hat <- rep( NA, N )
  Dev_obs_s <<- Dev_hat_s <<- matrix( nrow = N, ncol = S )
  
  for( n in 1:N ){ #loop over MCMC iterations
    for( s in 1:S ){ #loop over species
      for( i in 1:I ){ #loop over sites
        # sum log likelihoods
        ll_yobs_is[ n, i, s ] <- sum( log( lik_yobs[ n, i,  s,  ] ) )
        ll_yhat_is[ n, i, s ] <- sum( log( lik_yhat[ n, i, s,  ] ) ) 
        
      } #close I loop
      # sum log likelihoods
      ll_yobs_s[ n, s ] <- sum( ll_yobs_is[ n, , s ] )
      ll_yhat_s[ n, s ] <- sum( ll_yhat_is[ n, , s ] ) 
      # Calculate species-level deviances
      Dev_obs_s[ n, s ] <<- -2 * ( ll_yobs_s[ n, s ])
      Dev_hat_s[ n, s ] <<- -2 * ( ll_yhat_s[ n, s ] )
      
    } #close S loop

    #sum log likelihoods across detected species 
    ll_yobs[ n ] <- sum( ll_yobs_s[ n,  ] ) #for observed y
    ll_yhat[ n ] <- sum( ll_yhat_s[ n,  ] ) #for estimated y
    # Calculate overall deviances
    Dev_obs[ n ] <- -2 * ( ll_yobs[ n ])
    Dev_hat[ n ] <- -2 * ( ll_yhat[ n ] )
  } # close n loop
  
return( data.frame( Dev_obs, Dev_hat, ll_yobs, ll_yhat ) )
  
} # close function

# Run function to calculate deviances
ModDevs <- CalcDevs( lik_yobs = mr$sims.list$lik_yobs,
                     lik_yhat = mr$sims.list$lik_yhat )

#Plot the differences between observed and predicted deviance for each species:
colnames(Dev_obs_s ) <- spp[1:S]
gather( as.data.frame(Dev_obs_s - Dev_hat_s ), 
        key = Species, value = Deviance ) %>% 
  ggplot( ., aes( x = Species, y = Deviance ) ) +
  theme_classic() +
  geom_boxplot( ) + geom_hline( yintercept = 0, size = 2 )
#Note that Grey and Leopard frogs had the most uncertainty

#######################################################################
# Calculate Bayesian p-value as mean number of times that Deviance of #
# observed data was greater than the deviance of predicted model #####
#########################################################################
Baypvalue <- mean( ModDevs$Dev_obs > ModDevs$Dev_hat )
print( Baypvalue )

#plot Bayesian p value for the model
par( mfrow = c(1,1), cex.axis = 1.5, mar=c(3,4,2,2)+1)
plot( ModDevs$Dev_obs, ModDevs$Dev_hat, 
      main = paste0( "Bayesian P value = ", 
                    round( Baypvalue, 3 ) ), 
      xlab = "Deviance of observed data", 
      ylab = "Deviance of predicted data", tcl = 0.2, 
      cex.lab = 1.5, bty = "l"  ) 
abline( 0, 1, col = 'black', lwd = 3 )

#############       END OF SCRIPT         ###########################
