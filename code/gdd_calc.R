#CALCULATING GROWING DEGREE DAYS
#FUNCTION: degreedays() v 0.3 has been tested for accuracy but has not been optimized for performance.
#By: Elise Larsen, 2013-11-05, SESYNC
#Updated: 2021-02-23 for Pheno mismatch project 
#Contact: eal109 [at] georgetown.edu

##CHANGES: Default LDT set to 10 (assumed temperatures given in degrees Celsius).

#Single sine wave approximation for growing degree days with lower & upper development thresholds.
#Based on Degree-Days: The Calculation and Use of Heat Units in Pest Management FROM Zalom et al. Uc Coop. Extension 1983 (2m-12/83-WC/FB), Baskerville & Emin 1969.


##BEGIN degreedays() FUNCTION

#Calculated degree days using single sine wave approximation
degreedays=function(tmin,         #minimum daily temperature
                tmax,         #maximum daily temperature 
                ldt = 10,     #lower developmental threshold, default = 10C
                udt = 33) {   #upper developmental threshold, default = 33C
  
  ## CHECK FOR APPROPRIATE PARAMETERS
    if(missing(tmin) | missing(tmax)) {
    warning("No calculation: Missing Temperature Parameter(s)", immediate. = TRUE)
    return(NA)
    stop
  }
  ## Check that temps are numeric
  if(!is.numeric(tmin) | !is.numeric(tmax)){
    warning("No calculation: Temperature(s) not numeric", immediate. = TRUE)
    return(NA)
    stop
  }
  ## Check that tmax > tmin
  if(tmax < tmin){
    warning("No calculation: Maximum Temperature Less Than Minimum Temperature", immediate. = TRUE)
    return(NA)
    stop
  }
  ##CHECK LDT, UDT PARAMETERS
  if(udt < ldt){
    warning("Upper Temperature Threshold Less Than Lower Temperature Threshold", immediate. = TRUE)
    return(NA)
    stop
  }

  ##BEGIN CALCULATIONS
  ##Calculation case 1:  LDT < UDT < Tmin < Tmax: d = UDT-LDT
  if(tmin>=udt) {return(udt-ldt)} #max gdd
  ##Calculation case 2:  Tmin < Tmax < LDT < UDT: d = 0
  if(tmax<ldt) {return(0)} else {
    ##tmax >=ldt: some gdd will be calculated
    if(tmax<udt) {  
      ##Calculation case 3: LDT < Tmin < Tmax < UDT: d = beta - LDT  
      if(tmin>=ldt) {return(((tmax + tmin) / 2) - ldt) }
      #Calculation case 4: Tmin < LDT < Tmax < UDT: 
      if(tmin<ldt) {
        alpha<-((tmax - tmin) / 2)
        beta<-((tmax + tmin) / 2)
        theta1<-(asin((ldt - beta) / alpha))
        return( (1/pi) *((beta - ldt) * (pi/2 - theta1) + alpha * (cos(theta1))) )
      } 
    } else {     #tmax>=udt
    #Calculation: case 5: LDT < Tmin < UDT < Tmax: 
    if(tmin>=ldt & tmin<udt) {
      alpha<-((tmax - tmin) / 2)
      beta<-((tmax + tmin) / 2)
      theta2<-(asin((udt - beta) / alpha))
      return( (1/pi) * ((beta - ldt) * (pi/2 + theta2) + (udt - ldt) * (pi/2 - theta2) - alpha * (cos(theta2))))
    } 
    #Calculation: case 6: Tmin < LDT < UDT < Tmax: 
    if(tmin<ldt) {
      alpha<-((tmax - tmin) / 2)
      beta<-((tmax + tmin) / 2)
      theta1<-(asin((ldt - beta) / alpha))
      theta2<-(asin((udt - beta) / alpha))
      return( (1/pi) * ((beta - ldt) * (theta2 - theta1) + alpha * (cos(theta1) - cos(theta2)) + (udt - ldt) * (pi/2 - theta2)) )
     }
    }
  }
}
#END FUNCTION

