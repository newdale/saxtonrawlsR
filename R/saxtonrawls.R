#'
#' Use particle size analysis data to determine
#' Soil Water Characteristics based on Saxton and Rawls (2006)
#' Available at: https://doi.org/10.2136/sssaj2005.0117
#'
#' @param sand numeric, sand (0.5 - 2 mm) content in decimal format (e.g. 0.25)
#' @param clay numeric, clay (<0.002 mm) content in decimal format (e.g. 0.25)
#' @param OM numeric, Organic matter content in percentage (e.g., 2.5)
#' @param DF numeric, density factor, default value 1
#' @param Rw numeric, Volume fraction of gravel in decimal format, g/cm3, default is 0.
#' @param Rv numeric, Weight fraction of gravel in decimal format, g/g, default is NULL.
#' @return matrix
#' @export
#'
#' @examples
#' #Determine soil water characteristics for a single observation
#' saxtonrawls(sand=0.45, clay=0.35, OM=4, DF=1, Rw=0.10)
#'
#'
#' #Determine soil water characteristics for multiple observations
#' dat<- data.frame(sand=c(20,40,80), clay= c(65,30,10), OM=c(2,2,2), DF=c(1,1,1), Rw=c(0.1, 0.1, 0.1))
#' mapply(saxtonrawls,sand=dat$sand, clay=dat$clay, OM=dat$OM, DF=dat$DF, Rw=dat$Rw)
#'
#'
saxtonrawls<- function(sand,clay,OM,DF=1,Rw=0,Rv=NULL){

  if(is.na(sand)){stop('Sand content is a required argument, please provide a value')}
  if(is.na(clay)){stop('Clay content is a required argument, please provide a value')}
  if(is.na(OM)){stop('Organic matter content is a required argument, please provide a value')}
  if(DF<0.9 | DF>1.3){stop('DF values must be between 0.9 and 1.3')}
  if(OM>8){warning('Equations were calibrated with OM<8. Calculations with OM>8 should be scrutinized.')}
  if(clay>60){warning('Equations were calibrated with clay<60. Calculations with clay>60 should be scrutinized.')}
  if(clay>1 | sand >1){stop('Sand and clay content should be in decimal format but exceeds 1, please check input data')}
  #if(!is.null(Rw) & Rw>1){stop('Gravel content should be in decimal format but exceeds 1, please check input data')}
  #if(!is.null(Rv) & Rv>1){stop('Gravel content should be in decimal format but exceeds 1, please check input data')}
  if(clay+sand >1){stop('Sand+Clay exceeds 1 but should be in decimal format, please check input data')}

  # if gravel content is provided both as weight and volume, write warning and default to weight
  if(!is.null(Rw) & !is.null(Rv)){warning('Gravel content was provided on weight and volume basis, using weight basis')}
  if(!is.null(Rw) & !is.null(Rv)){Rv<- NULL}

  if(!is.na(sand) & !is.na(clay) & !is.na(OM)){

    # set the variables
    sand<- sand
    clay<- clay
    OM<- OM
    DF<- DF

    #estimate moisture at 1500 kPa
    o1500t<- -0.024*sand + 0.487*clay + 0.006*OM + 0.005*(sand*OM) - 0.013*(clay*OM) + 0.068*(sand*clay) + 0.031
    o1500<- o1500t + (0.14*o1500t - 0.02)

    #estimate moisture at 33 kPa
    o33t<- -0.251*sand + 0.195*clay + 0.011*OM + 0.006*(sand*OM) - 0.027*(clay*OM) + 0.452*(sand*clay) + 0.299
    o33<- o33t + (1.283*o33t^2 - 0.374*o33t - 0.015)

    #estimate moisture 0-33 kPa
    os_33t<- 0.278*sand + 0.034*clay + 0.022*OM - 0.018*(sand*OM) - 0.027*(clay*OM) - 0.584*(sand*clay) + 0.078
    os_33<- os_33t + (0.636*os_33t - 0.107)

    #estimate air entry suction
    psi_et<- -21.67*sand - 27.93*clay - 81.97*os_33 + 71.12*(sand*os_33) + 8.29*(clay*os_33) + 14.05*(sand*clay) + 27.16
    psi_e<- psi_et + (0.02*psi_et^2 - 0.113*psi_et - 0.70)

    #estimate moisture at 0 kPa (saturation)
    os<- o33 + os_33 - 0.097*sand + 0.043

    #estimate the bulk density based on saturation/porosity
    Dnorm<- (1-os)*2.65

    #now we can adjust for density factor
    Dadj<- Dnorm*DF
    osdf<- 1 - (Dadj/2.65)
    o33df<- o33 - 0.2*(os-osdf)
    os_33df<- osdf - o33df

    #estimate plant available water
    PAW<- o33df - o1500

    #estimate lambda, the slope or logarithmic tension-moisture curve
    B<- (log(1500)-log(33))/(log(o33df)-log(o1500))
    A<- exp(log(33) + B*log(o33))
    lambda<- 1/B

    #estimate saturated hydraulic conductivity
    ksat<- 1930*(osdf - o33df)^(3-lambda) #removed the df

    #modification for gravel content
    alpha<- Dadj/2.65

    # gravel can be provided either as Rv or Rw
    # if provided as Rw, Rv is estimated below
    # if provided as Rv, values should be used directly
    if(is.null(Rv)){
      Rv<- (alpha*Rw)/(1-(Rw*(1-alpha)))
    }else{Rv<- Rv}

    if(is.null(Rw)){
      Rw<- -Rv/(-Rv+(Rv*alpha)-alpha)
    }else{Rw<-Rw}

    Db<- Dadj*(1-Rv)+(Rv*2.65)
    PAWb<- PAW*(1-Rv)
    ksatb_s<- (1-Rw)/(1-Rw*(1-3*alpha/2))
    ksatb<- ksat*ksatb_s

    out<- (c(Wilting_Point=round(o1500*100,3), Field_Capacity=round(o33df*100,3), Saturation=round(osdf*100,3),
             Plant_Avail_Water =round(PAWb,3), KSat= round(ksatb,3), MatricBulkDensity= round(Dadj,3)))

    return(out)

  } else {
    return(c(Wilting_Point=NA, Field_Capacity=NA, Saturation=NA,
             Plant_Avail_Water=NA, KSat= NA, MatricBulkDensity= NA))

  }
}
