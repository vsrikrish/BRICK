##==============================================================================
## Define the coupled BRICK model and its inputs/outputs
## The inputs/outputs will vary based on which components you use, but an overview
## is given below.
##
## Input:
##  parameters.in          input vector of model parameters
##  parnames.in            vector of parameter names
##  forcing.in             matrix of radiative forcing input
##  l.project              making projections or hindcasts?
##  slope.Ta2Tg.in         slope of Antarctic vs global mean temperature regression
##  intercept.Ta2Tg.in     intercept of Antarctic vs global mean temperature regression
##  lws.params             mean and stdev of lws trend (in m/year)
##  mod.time               time (in years) of the model simulation
##  obs.temp               mean global surface temperature anomalies, used if no climate module
##  obs.emis               emissions by gas, used with Hector for projections
##  hector.params          name of Hector .ini file (for use with the R package) and scenario name
##  ind.norm.data          indices within the model output for setting zero anomaly of calibration data fields
##  ind.norm.sl            indices within model output for setting zero sea level
##  tstep                  model tstep [years]
##  i0                     index of reference year, within mod.time. For initial conditions to sub-models.
##  l.aisfastdy            logical, whether or not to use AIS fast dynamics emulator
##
## Requires:
##  luse.brick, includes: luse.doeclim, luse.gsic, luse.te, luse.tee, luse.simple,
##			  luse.dais, luse.lws, and luse.XXX, where XXX
##                        may be replaced with your favorite model component
##
## Questions? Tony Wong <twong@psu.edu>
##==============================================================================
## Copyright 2016 Tony Wong, Alexander Bakker
## This file is part of BRICK (Building blocks for Relevant Ice and Climate
## Knowledge). BRICK is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## BRICK is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with BRICK.  If not, see <http://www.gnu.org/licenses/>.
##==============================================================================

# function to find
is.error <- function(x) inherits(x, "try-error")

brick_model = function(parameters.in,
                       parnames.in,
                       forcing.in,
                       l.project = FALSE,
                       slope.Ta2Tg.in = 1,
                       intercept.Ta2Tg.in = 0,
                       lws.params = NULL,
                       tstep = 1,
                       mod.time,
                       obs.temp = NULL,
                       obs.emis = NULL,
                       hector.params = NULL,
                       ind.norm.data = NULL,
                       ind.norm.sl = NULL,
                       luse.brick,
                       i0,
                       l.aisfastdy
){

  # Initialize the list of output (do NOT grow lists/arrays in R)
  # The +1 is to have total global mean sea level (relative to ind.relative) in
  # the output.
  brick.out = vector('list',sum(luse.brick)+1)
  slr.out = rep(0,length(mod.time))
  outcnt=1

  # Initialize temperature for coupling. Most models require normalization to
  # preindustrial (1850-1870) values, so that should be the default normalization
  # period as set in ind.norm.data. Exceptions are handled at their model calls.
  temp.preindustrial = rep(NA, length(mod.time))

  # Initialize ocean heat change for coupling
  deltaH.couple = rep(NA, length(mod.time))

  #=============================================================================
  # SNEASY - Simple Nonlinear EArth SYstem model (including DOECLIM and CCM)

  if (luse.brick[,"luse.sneasy"]) {

    ## Grab the SNEASY parameters
    S            =parameters.in[match("S"           ,parnames.in)]
    kappa.sneasy =parameters.in[match("kappa.sneasy",parnames.in)]
    alpha.sneasy =parameters.in[match("alpha.sneasy",parnames.in)]
    Q10          =parameters.in[match("Q10"         ,parnames.in)]
    beta.sneasy  =parameters.in[match("beta.sneasy" ,parnames.in)]
    eta          =parameters.in[match("eta"         ,parnames.in)]
    h.sneasy     =parameters.in[match("h.sneasy"    ,parnames.in)]
    T0           =parameters.in[match("T0"          ,parnames.in)]
    H0           =parameters.in[match("H0"          ,parnames.in)]
    CO20         =parameters.in[match("CO20"        ,parnames.in)]
    MOC0         =parameters.in[match("MOC0"        ,parnames.in)]

    sneasy.out = sneasy(S=S, kappa=kappa.sneasy, alpha=alpha.sneasy, Q10=Q10,
                        beta=beta.sneasy, eta=eta, hydsens=h.sneasy, init.CO2=CO20,
                        init.MOC=MOC0, tstep=tstep, mod.time=mod.time,
                        forcing.co2=forcing.in$co2, forcing.aero=forcing.in$aero,
                        forcing.other=forcing.in$other)

    ## Normalize temperature to match the observations
    itmp = ind.norm.data[match("temp",ind.norm.data[,1]),2]:ind.norm.data[match("temp",ind.norm.data[,1]),3]
    sneasy.out$temp = sneasy.out$temp - mean(sneasy.out$temp[itmp])

    #itmp = ind.norm.data[match("ocheat",ind.norm.data[,1]),2]:ind.norm.data[match("ocheat",ind.norm.data[,1]),3]
    #sneasy.out$ocheat = sneasy.out$ocheat - mean(sneasy.out$ocheat[itmp])

    temp.preindustrial = sneasy.out$temp + T0
    deltaH.couple = diff(sneasy.out$ocheat) * 10^22 #in Joules
    brick.out[[outcnt]] = sneasy.out; names(brick.out)[outcnt]="sneasy.out"; outcnt=outcnt+1;

  }

  #=============================================================================
  # DOECLIM - climate and ocean energy balance

    if (luse.brick[,"luse.doeclim"]) {

    ## Grab the DOECLIM parameters
    S            =parameters.in[match("S.doeclim"            ,parnames.in)]
    kappa.doeclim=parameters.in[match("kappa.doeclim",parnames.in)]
    alpha.doeclim=parameters.in[match("alpha.doeclim",parnames.in)]
    T0           =parameters.in[match("T0"           ,parnames.in)]
    H0           =parameters.in[match("H0"           ,parnames.in)]

    ## Set up the radiative forcing
    forcing.total = forcing_total(forcing=forcing.in,
                                  alpha.doeclim=alpha.doeclim,
                                  l.project=l.project,
                                  begyear=mod.time[1],
                                  endyear=mod.time[length(mod.time)])

    ## Run DOECLIM at these parameter values
    doeclim.out = doeclimF(S=S, kappa=kappa.doeclim, forcing.total=forcing.total, mod.time=mod.time)

    ## Normalize temperature and ocean heat to match the observations
    itmp = ind.norm.data[match("temp",ind.norm.data[,1]),2]:ind.norm.data[match("temp",ind.norm.data[,1]),3]
    doeclim.out$temp = doeclim.out$temp - mean(doeclim.out$temp[itmp])

    #itmp = ind.norm.data[match("ocheat",ind.norm.data[,1]),2]:ind.norm.data[match("ocheat",ind.norm.data[,1]),3]
    #doeclim.out$ocheat = doeclim.out$ocheat - mean(doeclim.out$ocheat[itmp])

    temp.preindustrial = doeclim.out$temp + T0
    deltaH.couple = diff(doeclim.out$ocheat) * 10^22 #in Joules
    brick.out[[outcnt]] = doeclim.out; names(brick.out)[outcnt]="doeclim.out"; outcnt=outcnt+1;

    }
  
  #=============================================================================
  # Hector - climate and ocean energy balance
  
  if (luse.brick[,"luse.hector"]) {
    require(hector)
    
    # convert annual ocean heat flux (W/m^2) to cumulative ocean heat content anomaly (10^22 J)
    flux.to.heat = function(heatflux.mixed, heatflux.interior)
    {
      flnd = 0.29 # area land fraction
      fso = 0.95 # ocean area fraction of interior
      secs.per.year = 31556926
      earth.area = 510065600 * 10^6
      ocean.area = (1-flnd)*earth.area
      powtoheat = ocean.area*secs.per.year / 10^22 # in 10^22 J/yr
      
      heat.mixed = cumsum(heatflux.mixed) * powtoheat
      heat.interior = fso * cumsum(heatflux.interior) * powtoheat
      ocean.heat = heat.mixed + heat.interior
      
      return(list(ocean.heat=ocean.heat, heat.mixed=heat.mixed, heat.interior=heat.interior))
    }
    
    ## Grab the Hector parameters
    S            =as.numeric(parameters.in[match("S.temperature",parnames.in)])
    kappa        =as.numeric(parameters.in[match("diff.temperature",parnames.in)])
    alpha        =as.numeric(parameters.in[match("alpha.temperature",parnames.in)])
    T0           =as.numeric(parameters.in[match("offset.Tgav_obs", parnames.in)])
    ## start the Hector core with some initialization file
    inifile <- system.file(file.path('input', hector.params$inifile), package='hector', mustWork=TRUE)
    hcore <- newcore(inifile, suppresslogging=TRUE, name=hector.params$scenario)
    
    # if emissions are passed in, set those values
    if (!is.null(obs.emis)) {
      setvar(hcore, obs.emis$year, EMISSIONS_BC(), obs.emis[, 'BC'], getunits(EMISSIONS_BC()))
      setvar(hcore, obs.emis$year, EMISSIONS_C2F6(), obs.emis[, 'C2F6'], getunits(EMISSIONS_C2F6()))
      setvar(hcore, obs.emis$year, EMISSIONS_CF4(), obs.emis[, 'CF4'], getunits(EMISSIONS_CF4()))
      setvar(hcore, obs.emis$year, EMISSIONS_CH4(), obs.emis[, 'CH4'], getunits(EMISSIONS_CH4()))
      setvar(hcore, obs.emis$year, EMISSIONS_CO(), obs.emis[, 'CO'], getunits(EMISSIONS_CO()))
      setvar(hcore, obs.emis$year, FFI_EMISSIONS(), obs.emis[, 'CO2'], getunits(FFI_EMISSIONS()))
      setvar(hcore, obs.emis$year, EMISSIONS_HFC125(), obs.emis[, 'HFC125'], getunits(EMISSIONS_HFC125()))
      setvar(hcore, obs.emis$year, EMISSIONS_HFC134A(), obs.emis[, 'HFC134a'], getunits(EMISSIONS_HFC134A()))
      setvar(hcore, obs.emis$year, EMISSIONS_HFC143A(), obs.emis[, 'HFC143a'], getunits(EMISSIONS_HFC143A()))
      setvar(hcore, obs.emis$year, EMISSIONS_HFC227EA(), obs.emis[, 'HFC227ea'], getunits(EMISSIONS_HFC227EA()))
      setvar(hcore, obs.emis$year, EMISSIONS_HFC23(), obs.emis[, 'HFC23'], getunits(EMISSIONS_HFC23()))
      setvar(hcore, obs.emis$year, EMISSIONS_HFC245FA(), obs.emis[, 'HFC245fa'], getunits(EMISSIONS_HFC245FA()))
      setvar(hcore, obs.emis$year, EMISSIONS_HFC32(), obs.emis[, 'HFC32'], getunits(EMISSIONS_HFC32()))
      setvar(hcore, obs.emis$year, EMISSIONS_N2O(), obs.emis[, 'N2O'], getunits(EMISSIONS_N2O()))
      setvar(hcore, obs.emis$year, EMISSIONS_NMVOC(), obs.emis[, 'NMVOC'], getunits(EMISSIONS_NMVOC()))
      setvar(hcore, obs.emis$year, EMISSIONS_NOX(), obs.emis[, 'NOx'], getunits(EMISSIONS_NOX()))
      setvar(hcore, obs.emis$year, EMISSIONS_OC(), obs.emis[, 'OC'], getunits(EMISSIONS_OC()))
      setvar(hcore, obs.emis$year, EMISSIONS_SF6(), obs.emis[, 'SF6'], getunits(EMISSIONS_SF6()))
      setvar(hcore, obs.emis$year, EMISSIONS_SO2(), obs.emis[, 'SO2'], getunits(EMISSIONS_SO2()))
      setvar(hcore, obs.emis$year, LUC_EMISSIONS(), obs.emis[, 'LUC'], getunits(LUC_EMISSIONS()))
    }
    
    # run Hector at the parameter values
    setvar(hcore, NA, ECS(), S, 'degC')
    setvar(hcore, NA, DIFFUSIVITY(), kappa, 'cm2/s')
    setvar(hcore, NA, AERO_SCALE(), alpha, '(unitless)')
    begyear <- mod.time[1]
    endyear <- mod.time[length(mod.time)]
    reset(hcore)
    # since weird parameter values might get passed to Hector, we will try to run the core and catch an error if it fails
    run_out <- try(run(hcore,  endyear), silent=TRUE)
    # if the run failed, return NA
    if (is.error(run_out)) {
      return(NA)
    }
    temp <- fetchvars(hcore, begyear:endyear, GLOBAL_TEMP())$value
    flux_mixed <- fetchvars(hcore, begyear:endyear, FLUX_MIXED())$value
    flux_interior <- fetchvars(hcore, begyear:endyear, FLUX_INTERIOR())$value
    shutdown(hcore)
    ocheat <- flux.to.heat(flux_mixed, flux_interior)$ocean.heat
    hector.out <- list(temp=temp, ocheat=ocheat)
      
    ## Normalize temperature and ocean heat to match the observations
    itmp = ind.norm.data[match("temp",ind.norm.data[,1]),2]:ind.norm.data[match("temp",ind.norm.data[,1]),3]
    hector.out$temp = hector.out$temp - mean(hector.out$temp[itmp])

    temp.preindustrial = hector.out$temp + T0
    deltaH.couple = diff(hector.out$ocheat) * 10^22 #in Joules
    brick.out[[outcnt]] = hector.out; names(brick.out)[outcnt]="hector.out"; outcnt=outcnt+1;
  }
  

  #=============================================================================
  # Establish coupling if there is no module to estimate global mean surface temperature
  # Normalization period for obs.temp should match those passed by SNEASY and DOECLIM (1850-1870)

  if (!luse.brick[,"luse.doeclim"] & !luse.brick[, "luse.hector"] & !luse.brick[,"luse.sneasy"]) {temp.preindustrial=obs.temp}

  #=============================================================================
  # GSIC-MAGICC - glaciers and small ice caps

  if (luse.brick[,"luse.gsic"]) {

    ## Grab the GSIC parameters
    beta0  = parameters.in[match("beta0_gsic.slr_brick"  ,parnames.in)]
    V0.gsic= parameters.in[match("V0_gsic.slr_brick",parnames.in)]
    n      = parameters.in[match("n_gsic.slr_brick"      ,parnames.in)]
    Gs0    = parameters.in[match("Gs0_gsic.slr_brick"    ,parnames.in)]

    ## Normalize temperature to match what the sub-model expects (the parameters
    ## may assume a particular time period associated with Tg=0, for example)
    # GSIC-MAGICC expects temp.preindustrial relative to late 1800s, which it already is
    # with ind.norm.data for temperature rel to 1850-70 (Wigley and Raper 2005)

    ## Run GSIC-MAGICC at these parameter values, using temperature output from DOECLIM
    gsic.out = gsic_magiccF(beta0=beta0, V0=V0.gsic, n=n, Gs0=Gs0 , Tg=temp.preindustrial, i0=i0$gsic)

    ## Subtract off normalization period model GSIC output as the zero point
    itmp = ind.norm.data[match("gsic",ind.norm.data[,1]),2]:ind.norm.data[match("gsic",ind.norm.data[,1]),3]

    gsic.out.norm = gsic.out - mean(gsic.out[itmp])

    brick.out[[outcnt]] = gsic.out.norm; names(brick.out)[outcnt]="gsic.out"; outcnt=outcnt+1;

    ## Add this contribution to the total sea level rise
    slr.out = slr.out + (gsic.out - mean(gsic.out[ind.norm.sl]))

  }

  #=============================================================================
  # TE - thermal expansion

  if (luse.brick[,"luse.te"]) {

    ## Grab the BRICK-TE parameters
    a.te     =parameters.in[match("a_te.slr_brick"     ,parnames.in)]
    b.te     =parameters.in[match("b_te.slr_brick"     ,parnames.in)]
    invtau.te=parameters.in[match("invtau_te.slr_brick",parnames.in)]
  #  TE0      =parameters.in[match("TE0"      ,parnames.in)]

    ## Normalize temperature to match what the sub-model expects (the parameters
    ## may assume a particular time period associated with Tg=0, for example)
    ## TE expects temp.preindustrial relative to late 1800s, which it already is
    ## with ind.norm.data for temperature rel to 1850-70

    ## Run BRICK-TE (thermosteric expansion) model, using temp output from DOECLIM
    ## i0$te=1
    te.out = brick_te_F(a=a.te , b=b.te, invtau=invtau.te, Tg=temp.preindustrial)

    ## Subtract off normalization period
    itmp = ind.norm.data[match("te",ind.norm.data[,1]),2]:ind.norm.data[match("te",ind.norm.data[,1]),3]
    te.out.norm = te.out - mean(te.out[itmp])

    brick.out[[outcnt]] = te.out.norm; names(brick.out)[outcnt]="te.out"; outcnt=outcnt+1;

    ## Add this contribution to the total sea level rise
    slr.out = slr.out + (te.out.norm - mean(te.out.norm[ind.norm.sl]))

  }

  #=============================================================================
  # TEE - explicit thermosteric expansion,
  #       optional replacement for TE, requires a deltaH time series from doeclim

  if (luse.brick[,"luse.tee"]) {

    ## Grab the BRICK-TEE parameters
    a.tee  =parameters.in[match("a_tee.slr_brick",parnames.in)]
  #  TE0    =parameters.in[match("TE0"  ,parnames.in)]

    ## Run BRICK-TEE (explicit thermosteric expansion) model, using OHC output from DOECLIM
    ## i0$te=1
    te.out = brick_tee_F(a=a.tee, deltaH=deltaH.couple)

    ## Subtract off normalization period
    itmp = ind.norm.data[match("te",ind.norm.data[,1]),2]:ind.norm.data[match("te",ind.norm.data[,1]),3]
    te.out.norm = te.out - mean(te.out[itmp])

    brick.out[[outcnt]] = te.out.norm; names(brick.out)[outcnt]="te.out"; outcnt=outcnt+1;

    ## Add this contribution to the total sea level rise
    slr.out = slr.out + (te.out.norm - mean(te.out.norm[ind.norm.sl]))
  }

  #=============================================================================
  # SIMPLE - Greenland ice sheet

  if (luse.brick[,"luse.simple"]) {

    ## Grab SIMPLE parameters
    a.simple    =parameters.in[match("a_simple.slr_brick"    ,parnames.in)]
    b.simple    =parameters.in[match("b_simple.slr_brick"    ,parnames.in)]
    alpha.simple=parameters.in[match("alpha_simple.slr_brick",parnames.in)]
    beta.simple =parameters.in[match("beta_simple.slr_brick" ,parnames.in)]
    V0          =parameters.in[match("V0_simple.slr_brick"          ,parnames.in)]

    ## Normalize temperature to match what the sub-model expects (the parameters
    ## may assume a particular time period associated with Tg=0, for example)
    # SIMPLE expects temp.simple relative to 1960-1990. i0$gis should match this.
    temp.simple = temp.preindustrial - mean(temp.preindustrial[i0$gis])

    ## Run SIMPLE (Greenland Ice Sheet model)
    simple.out = simpleF(a=a.simple, b=b.simple, alpha=alpha.simple,
                         beta=beta.simple, V0=V0, Tg=temp.simple, i0=i0$gis)

    ## Add this contribution to the total sea level rise
    slr.out = slr.out + (simple.out$sle.gis - mean(simple.out$sle.gis[ind.norm.sl]))

    ## Subtract off normalization period
    itmp = ind.norm.data[match("gis",ind.norm.data[,1]),2]:ind.norm.data[match("gis",ind.norm.data[,1]),3]
    simple.out$sle.gis = simple.out$sle.gis - mean(simple.out$sle.gis[itmp])

    brick.out[[outcnt]] = simple.out; names(brick.out)[outcnt]="simple.out"; outcnt=outcnt+1;

  }

  #=============================================================================
  # DAIS - Antarctic ice sheet

  if (luse.brick[,"luse.dais"]) {

    ## Grab DAIS parameters
    anto.a=parameters.in[match("anto.a.slr_brick",parnames.in)]
    anto.b=parameters.in[match("anto.b.slr_brick",parnames.in)]
    gamma =parameters.in[match("gamma.slr_brick" ,parnames.in)]
    alpha.dais =parameters.in[match("alpha.dais.slr_brick" ,parnames.in)]
    mu =parameters.in[match("mu.slr_brick" ,parnames.in)]
    nu =parameters.in[match("nu.slr_brick" ,parnames.in)]
    P0 =parameters.in[match("P0.slr_brick" ,parnames.in)]
    kappa.dais =parameters.in[match("kappa.dais.slr_brick" ,parnames.in)]
    f0 =parameters.in[match("f0.slr_brick" ,parnames.in)]
    h0 =parameters.in[match("h0.slr_brick" ,parnames.in)]
    c =parameters.in[match("c.slr_brick" ,parnames.in)]
    b0 =parameters.in[match("b0.slr_brick" ,parnames.in)]
    slope =parameters.in[match("slope.slr_brick" ,parnames.in)]
    if(l.aisfastdy) {
      Tcrit = parameters.in[match("Tcrit.slr_brick",parnames.in)]
      lambda = parameters.in[match("lambda.slr_brick",parnames.in)]
    } else {
      Tcrit = NULL
      lambda = NULL
    }

		## Calculate the sea level updated from the other model components'
		## contributions. From Shaffer (2014), SL should be relative to 1961-1990
		## mean. Implement fingerprinting of local sea-level sources on AIS?

    l.fprint=TRUE

    SL.couple = slr.out
    if(l.fprint) {
      dSL.gis = diff(simple.out$sle.gis)
      dSL.gsic= diff(gsic.out)
      dSL.te  = diff(te.out)
      for (i in 2:length(mod.time)) {
        SL.couple[i] = SL.couple[i-1] + tstep*(1.0*dSL.gis[i-1] +
                                               1.0*dSL.gsic[i-1]+
                                               1.0*dSL.te[i-1]   )
      }
    }

    ## Normalize
    itmp = ind.norm.data[match("sl",ind.norm.data[,1]),2]:ind.norm.data[match("sl",ind.norm.data[,1]),3]
    SL.couple = SL.couple - mean(SL.couple[itmp])
    dSL.couple = c(-999,diff(SL.couple))
    include_dSLais = 0    # in coupled model, feeding AIS dSL without AIS contribution

    ## Check to make sure output from other models was reasonable
    if(any(is.na(SL.couple))) {
      slr.out=rep(NA,length(mod.time))
      brick.out[[outcnt]] = slr.out; names(brick.out)[outcnt]="dais.out"; outcnt=outcnt+1;
    } else {
      dais.out = daisanto_fastdynF(anto.a=anto.a , anto.b=anto.b,
                                   gamma=gamma   , alpha=alpha.dais,
                                   mu=mu         , nu=nu        ,
                                   P0=P0         , kappa=kappa.dais,
                                   f0=f0         , h0=h0        ,
                                   c=c           , b0=b0        ,
                                   slope=slope   , l.aisfastdy=l.aisfastdy,
                                   Tcrit=Tcrit   , lambda=lambda,
                                   slope.Ta2Tg=slope.Ta2Tg.in, intercept.Ta2Tg=intercept.Ta2Tg.in,
                                   Tg=temp.preindustrial, SL=SL.couple, dSL=dSL.couple ,
                                   includes_dSLais = include_dSLais)
      ## Subtract off normalization period
      itmp = ind.norm.data[match("ais",ind.norm.data[,1]),2]:ind.norm.data[match("ais",ind.norm.data[,1]),3]
      dais.out$Vais = dais.out$Vais - mean(dais.out$Vais[itmp])
      brick.out[[outcnt]] = dais.out; names(brick.out)[outcnt]="dais.out"; outcnt=outcnt+1;
      ## Add this contribution to the total sea level rise
      slr.out = slr.out + (dais.out$Vais - mean(dais.out$Vais[ind.norm.sl]))
    }
  }

  #=============================================================================
  # LWS - land water storage contributions

  if (luse.brick[,"luse.lws"]) {

    ## Grab LWS parameters
    # (fed in from lws.params)

    lws.out <- brick_lws(lws.mean=lws.params$mean, lws.sd=lws.params$sd,
                         lws0=lws.params$lws0, Tg=temp.preindustrial)

    ## Subtract off normalization period
    lws.out <- lws.out - mean(lws.out[ind.norm.sl])

    ## Add this contribution to the total sea level rise
    slr.out = slr.out + lws.out

    brick.out[[outcnt]] = lws.out; names(brick.out)[outcnt]="lws.out"; outcnt=outcnt+1;

  }

  #=============================================================================
  # Total sea-level rise

  ## Add the SLR to the output
  brick.out[[sum(luse.brick)+1]] = slr.out - mean(slr.out[ind.norm.sl])
  names(brick.out)[sum(luse.brick)+1]="slr.out"

  ## Check to make sure all the output made it
  if(outcnt!=sum(luse.brick)+1) print('ERROR - missing model output!')

  #=============================================================================

  return(brick.out)
}

##==============================================================================
## End
##==============================================================================
