#############################################################################
## This script was developed by Dr. Jen Cruz for Sauer et al. (in review) #
# Multiscale drivers of amphibian occupancy in urban ponds.             ###
# The script fits a multispecies, single-season occupancy model in JAGS ###
#                                                                       ###
# Data: Detection of amphibian species at 96 ponds in Wisconsin         ##    
##   surveyed in two occasions over one breeding season.                 # 
#                                                                        #
#  Predictors include habitat cover, water chemistry and pond metrics.   #
# Model includes group category that divide species based on life history #
##############################################################################

########### clean workspace and load required packages ####################
###########################################################################

#####clean workspace to improve efficiency: ###
rm(list = ls() ) 

####### load relevant packages ###
library( tidyverse ) #dataframe manipulations.
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( jagsUI ) #to run RJAGS

########## end of package loading ###########
#######################    import relevant data   ##############
#start with workspace created in DataPrep.R, which has objects #
# required for analysis:
load( "DataWorkspace.RData" )

#########################################################################################
###########################################################################
####################### define MCMC settings ##############################
#thinning, burnin, chains
nt <- 5; nb <- 15000; nc <- 3 

##### end of MCMC parameters definition ##############

#### Define single season, multispecies occupancy model in JAGS    ###
sink( "cm1.txt" )
cat( "
     model{
     
      #priors
      #for detection model: 
      #define intercept as mean probs:
      int.det <- log( mean.int.det / ( 1 - mean.int.det ) ) 
      mean.int.det ~ dbeta( 4, 4 )
      
      #priors for detection coefficients:
      for( n in 1:2 ){
        #define as a slightly informative prior
        alpha.det[n] ~ dnorm( 0, 0.1 ) T(-7, 7 )
       }
      
      #random species intercept for detection
      for ( s in 1:S ){  #loop over species
        eps.det[ s ] ~ dnorm( 0, pres.det ) T(-7, 7) 
      } #s
      
     #associated variance of random intercepts:     
     pres.det <- 1/ ( sigma.det * sigma.det )
     #sigma prior specified as a student t half-normal:
     sigma.det ~ dt( 0, 2.5, 7 ) T( 0, )
     
      #for occupancy model:
      #define intercept mean prob:
      # with separate intercepts for each group
      for( c in 1:2 ){
        int.occ[c] <- log( mean.int.occ[c] / ( 1 - mean.int.occ[c] ) )
        mean.int.occ[c] ~ dbeta( 4, 4 )
      }
       
      #fixed coefs representing group effects
      for( n in 1:OP ) { #loop over all occupancy predictors
        for( c in 1:2 ){ #loop over groups         
          beta.spp [ c, n ] ~ dnorm( 0, 0.1 ) #slighly informative prior
       }#c
      }#n
      
      #ecological model for occupancy:
      for( s in 1:S ){ #loop over species
          for ( i in 1:I ){  #loop over sites
            #latent, estimated true occupancy 
            z[ i, s ] ~ dbern( psi[ i, s ] ) 
            #prob of occupancy at site i for species s
            logit( psi[ i, s ] ) <- int.occ[ gp[s] ] + 
                    #fixed, group effects:
                    inprod( beta.spp[ gp[s],1:OP ], XI[i,] ) 

          } #close i loop
        } #close s loop
      
      #observation model:
      for ( i in 1:I ){  #loop over sites
        for( j in 1:J ){  #loop over surveys
          for( s in 1:S ){ #loop over species
                   
            logit( p[ i, s, j ] ) <- int.det + eps.det[ s ] +
                  # floating vegetation cover %:
                  alpha.det[1] * FV.pSTD[i,j] +
                  #date of survey
                  alpha.det[2] * DateSTD[i,j] 
          
            #linking both model outputs to observations
             y_obs[ i, s, j ] ~ dbern( z[ i, s ] * p[ i, s, j ] ) 
             
            #estimated detections from model: 
             yhat[ i, s, j ] ~ dbern( z[ i, s ] * p[ i, s, j ] ) 
             
            } #close s loop
          } #close j loop
        } #close i loop
        
       for ( i in 1:I ){  #loop over sites
         for( j in 1:J ){  #loop over surveys
           for( s in 1:S ){ #loop over species

              # Bernoulli likelihood of observations:
              lik_yobs[ i, s, j ] <- ( ( psi[ i, s ] * p[ i, s, j ] )^y_obs[ i, s, j ] ) *
                            ( ( 1 - psi[ i, s ] * p[ i, s, j ] )^( 1 - y_obs[ i, s, j ] ) )

              # likelihood of estimated detections:
              lik_yhat[ i, s, j ] <- ( ( psi[ i, s ]* p[ i, s, j ] )^yhat[ i, s, j ] ) *
                            ( ( 1 - psi[ i, s ] * p[ i, s, j ] )^( 1 - yhat[ i, s, j ] ) )
          }#s
        } #j
      } #i

        
    } #model close
     
     ", fill = TRUE )

sink()
################ end of model specification  #####################################
modelname <- "cm1.txt"
#parameters monitored #only keep those relevant for model comparisons (with different variances)
params <- c(  'int.det' #intercept for detection
              ,'int.occ'  #intercept for occupancy
              , 'eps.det' #random species intercepts in detection
              , 'alpha.det'  #fixed coefficients for detection
              , 'beta.spp' #random slopes in occupancy
              , 'psi' #occupancy probability
              , 'p' #detection probability
              , 'sigma.eps', 'sigma.det'#,'beta.sig'#std devs for random intercepts
               , 'lik_yobs' #likelihood for each occupancy observation
               , 'lik_yhat' #likelihood for occupancy observations predicted by the model
               , 'yhat' #estimated occurrence from model
)

#Predictors excluded from original list:
exc <-c( "predB","year","nearest",
         "Grass", "Wetlands" )
modcovs <- Statvarnames2[ !( Statvarnames2 %in% exc ) ]
modcovs 
#removed rare species:
y_obs <- AO3D
dim(y_obs )
spp
y_obs <- y_obs[,1:5,]
dim(y_obs)
S <- dim(y_obs)[2]
#create gp vector that groups species into two groups based on life-history
gp <- c( 2,1,1,2,2)
#Initial values for the model coefficients
zst <- matrix(data = 1, nrow = I, ncol = S )
inits <- function(){ list( beta.occ = rnorm( length(modcovs) ), 
                           alpha.det =  rnorm( 2 ), 
                           z = zst) }
#Replace missing date and FV with mean
DateSTD[ is.na(DateSTD[,2] ),2 ] <- 0
FV.pSTD[ is.na(FV.pSTD[,2] ),2 ] <- 0

#define data to go in model:
str( win.data <- list( y_obs = y_obs, #detections for each species
                       I = I , J = J, S = S, #sites, surveys, species
                       OP = length(modcovs), #occupancy predictors
                       #Static predictors
                       XI = as.matrix( StatVarSTD2[ , modcovs] ),
                       #IXJ predictors:
                       FV.pSTD = FV.pSTD, DateSTD = DateSTD
                       ,gp = gp
) )                
#call JAGS and summarize posteriors:
cm1 <- autojags( win.data, inits = inits, params, modelname, #
                 n.chains = nc, n.thin = nt, n.burnin = 0,
                 iter.increment = 5000, max.iter = 300000, 
                 Rhat.limit = 1.01,
                 save.all.iter = FALSE, parallel = TRUE ) 

########## end  ######################
##################################################################
### save workspace ###
save.image( "OccResults.RData" )

#############################end of script #######################