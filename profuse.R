############################################ profuseAllStar.R ############################################
profuseAllStarFound2Fit = function(image,
                           sigma = NULL,
                           locs = NULL,
                           segim = NULL,
                           magzero = 0,
                           resamp = NULL,
                           psf_dim = c(51,51),
                           star_con = 2,
                           star_con_fit = TRUE,
                           star_circ = TRUE,
                           rough = FALSE,
                           star_dom_mag = NULL,
                           Nstar = 4,
                           ...) {

  message('    Running ProFound')
  if(!requireNamespace("ProFound", quietly = TRUE)){stop('The ProFound package is required to run this function!')}

  image[image < quantile(image, 0.01, na.rm=TRUE)*10] = NA
  image[image > quantile(image, 0.99, na.rm=TRUE)*10] = NA

  if(is.null(segim)){
    mini_profound = ProFound::profoundProFound(
      image = image,
      sky = 0,
      redosky = FALSE,
      magzero = magzero,
      verbose = FALSE,
      boundstats = TRUE,
      ...
    )
    segim = mini_profound$segim
  }else{
    mini_profound = ProFound::profoundProFound(
      image = image,
      segim = segim,
      sky = 0,
      redosky = FALSE,
      magzero = magzero,
      verbose = FALSE,
      boundstats = TRUE,
      iters = 0,
      ...
    )
  }


  if(is.null(sigma)){
    gain = ProFound::profoundGainEst(image, objects = mini_profound$objects, sky = 0)
    sigma = ProFound::profoundMakeSigma(
      image = image,
      objects = mini_profound$objects,
      gain = gain,
      sky = 0,
      skyRMS = mini_profound$skyRMS,
      plot = FALSE
    )
  }

  if(!is.null(resamp)){
    if(resamp > 1){
      sigma = sigma*resamp
    }
  }

  if(is.null(locs)){
    if(is.null(star_dom_mag)){
      star_dom_mag = median(mini_profound$segstats$mag, na.rm=TRUE)
    }
    cut_R50 = median(mini_profound$segstats[mini_profound$segstats$mag < star_dom_mag,'R50'], na.rm=TRUE)
    locs = mini_profound$segstats[mini_profound$segstats$mag < star_dom_mag & mini_profound$segstats$R50 < cut_R50,c("xcen","ycen")]
  }
  locs = as.matrix(rbind(locs))

  segID_tar = unique(mini_profound$segim[locs])
  segID_tar[segID_tar > 0]

  if(length(segID_tar) == 0){
    message('All stars are detected as sky (lower skycut)!')
    return(NULL)
  }

  segID_tar = mini_profound$segstats[mini_profound$segstats$segID %in% segID_tar & mini_profound$segstats$Nobject==0 & mini_profound$segstats$Nborder==0 & mini_profound$segstats$Nmask==0,'segID']

  if(length(segID_tar) < Nstar){
    segID_tar = unique(mini_profound$segim[locs])
    segID_tar[segID_tar > 0]
    segID_tar = mini_profound$segstats[mini_profound$segstats$segID %in% segID_tar & mini_profound$segstats$Nobject==0,'segID']
  }

  if(length(segID_tar) > Nstar){
    segID_tar = segID_tar[1:Nstar]
  }else{
    Nstar = length(segID_tar)
  }

  if(Nstar == 0){
    message('All stars too nearby other objects!')
    return(NULL)
  }

  loc_tar = which(mini_profound$segstats$segID %in% segID_tar)

  xcen = mini_profound$segstats[loc_tar, 'xcen']
  ycen = mini_profound$segstats[loc_tar, 'ycen']

  region = matrix(mini_profound$segim %in% segID_tar, nrow=dim(image)[1], ncol=dim(image)[2])

  gridNy = ceiling(sqrt(Nstar))
  gridNx = ceiling(Nstar/gridNy)

  grid = expand.grid((1:gridNx - 0.5) * psf_dim[1], (1:gridNy - 0.5) * psf_dim[2])

  image_psf = matrix(0, nrow=psf_dim[1]*gridNx, ncol=psf_dim[2]*gridNy)
  sigma_psf = matrix(0, nrow=psf_dim[1]*gridNx, ncol=psf_dim[2]*gridNy)
  region_psf = matrix(0, nrow=psf_dim[1]*gridNx, ncol=psf_dim[2]*gridNy)

  for(i in 1:Nstar){
    subx = 1:psf_dim[1] + grid[i,1] - psf_dim[1]/2
    suby = 1:psf_dim[2] + grid[i,2] - psf_dim[2]/2
    image_psf[subx,suby] = magcutout(image, loc=c(xcen[i],ycen[i]), box=psf_dim)$image
    sigma_psf[subx,suby] = magcutout(sigma, loc=c(xcen[i],ycen[i]), box=psf_dim)$image
    region_psf[subx,suby] = magcutout(region, loc=c(xcen[i],ycen[i]), box=psf_dim)$image
  }

  region_psf[is.na(image_psf) | is.na(region_psf)] = 0

  if(star_circ){
    ang = 0
    axrat = 1
  }else{
    ang = median(mini_profound$segstats[loc_tar, 'ang'],na.rm=TRUE)
    axrat = median(mini_profound$segstats[loc_tar, 'axrat'],na.rm=TRUE)
  }

  size = median(mini_profound$segstats[loc_tar, 'R50'],na.rm=TRUE)

  modellist = list(
    moffat = list(
      xcen = grid[1:Nstar,1],
      ycen = grid[1:Nstar,2],
      mag = mini_profound$segstats[loc_tar, 'mag'],
      fwhm = rep(size, Nstar),
      con = rep(star_con, Nstar),
      ang = rep(ang, Nstar),
      axrat = rep(axrat, Nstar)
    )
  )

  tofit = list(
    moffat = list(
      xcen = rep(TRUE,Nstar),
      ycen = rep(TRUE,Nstar),
      mag = rep(TRUE,Nstar),
      fwhm = c(TRUE, rep(NA, Nstar-1)),
      con = c(star_con_fit, rep(NA, Nstar-1)),
      ang = c(!star_circ, rep(NA, Nstar-1)),
      axrat = c(!star_circ, rep(NA, Nstar-1))
    )
  )

  tolog = list(
    moffat = list(
      xcen = rep(FALSE, Nstar),
      ycen = rep(FALSE, Nstar),
      mag = rep(FALSE, Nstar),
      fwhm = rep(TRUE, Nstar), #fwhm is best fit in log space
      con = rep(TRUE, Nstar), #con is best fit in log space
      ang = rep(FALSE, Nstar),
      axrat = rep(TRUE, Nstar) #axrat is best fit in log space
    )
  )

  intervals = list(moffat = list(
    xcen = rep(list(c(-50, dim(image)[1] + 50)), Nstar),
    ycen = rep(list(c(-50, dim(image)[2] + 50)), Nstar),
    mag = rep(list(c(0, 40)), Nstar),
    fwhm = rep(list(c(0.5, size*4)), Nstar),
    con = rep(list(c(1, 10)), Nstar),
    ang = rep(list(c(-180, 360)), Nstar),
    axrat = rep(list(c(0.1, 1)), Nstar)
  ))

  Data = profitSetupData(
    image = image_psf,
    region = region_psf,
    sigma = sigma_psf,
    segim = segim,
    psf = NULL,
    modellist = modellist,
    tofit = tofit,
    tolog = tolog,
    intervals = intervals,
    magzero = magzero,
    algo.func = 'LD',
    verbose = FALSE,
    rough = rough
  )
  Data$Nmod = Nstar
  return(invisible(list(profound = mini_profound, Data = Data)))
}

profuseAllStarDoFit = function(image,
                       sigma = NULL,
                       locs = NULL,
                       magzero = 0,
                       psf_dim = c(51,51),
                       rough = FALSE,
                       plot = FALSE,
                       seed = 666,
                       ...) {

  timestart = proc.time()[3] # start timer
  call = match.call(expand.dots=TRUE)

  message('Running Found2Fit')
  found2fit = profuseAllStarFound2Fit(
    image = image,
    sigma = sigma,
    locs = locs,
    magzero = magzero,
    psf_dim = psf_dim,
    rough = rough,
    ...
  )
  Data = found2fit$Data

  if (plot) {
    profitLikeModel(parm = Data$init,
                    Data = Data,
                    makeplots = TRUE)
    legend('topright', legend = 'Start')
  }

  lowers = unlist(Data$intervals)[c(T, F)]
  lowers[unlist(Data$tolog) == T] = log10(lowers[unlist(Data$tolog) == T])
  lowers = lowers[which(unlist(Data$tofit))]
  uppers = unlist(Data$intervals)[c(F, T)]
  uppers[unlist(Data$tolog) == T] = log10(uppers[unlist(Data$tolog) == T])
  uppers = uppers[which(unlist(Data$tofit))]

  message('Running Highander')
  if(!requireNamespace("ProFound", quietly = TRUE)){stop('The Highander package is required to run this function!')}
  highfit = Highlander::Highlander(
    parm = Data$init,
    Data = Data,
    likefunc = profitLikeModel,
    seed = seed,
    lower = lowers,
    upper = uppers,
    applyintervals = FALSE,
    applyconstraints = FALSE
  )
  names(highfit$parm) = names(Data$init)

  if (plot) {
    profitLikeModel(highfit$parm, Data = Data, makeplots = TRUE)
    legend('topright', legend = 'After')
  }

  highfit$profound = found2fit$profound
  highfit$Data = Data
  highfit$initmodel = profitRemakeModellist(Data$init, Data = Data)
  highfit$finalmodel = profitRemakeModellist(highfit$parm, Data = Data)

  highfit$parm[highfit$parm < lowers] = lowers[highfit$parm < lowers]
  highfit$parm[highfit$parm > uppers] = uppers[highfit$parm > uppers]

  highfit$error = apply(highfit$LD_last$Posterior1,
                        MARGIN = 2,
                        FUN = 'sd')

  if(psf_dim[1] %% 2 == 0){psf_dim[1] = psf_dim[1] + 1}
  if(psf_dim[2] %% 2 == 0){psf_dim[2] = psf_dim[2] + 1}

  temp_modellist = list(
    moffat = list(
      xcen = psf_dim[1]/2,
      ycen = psf_dim[2]/2,
      mag = 0,
      fwhm = highfit$finalmodel$modellist$moffat$fwhm[1],
      con = highfit$finalmodel$modellist$moffat$con[1],
      ang = highfit$finalmodel$modellist$moffat$ang[1],
      axrat = highfit$finalmodel$modellist$moffat$axrat[1]
    )
  )

  highfit$psf = profitMakeModel(temp_modellist, dim=psf_dim)
  highfit$psf_modellist = temp_modellist

  psf_fluxcheck = sum(highfit$psf$z)

  if(psf_fluxcheck < 0.95){
    message('WARNING: psf output image contains less than 95% of the total model flux! Consider increasing the size of psf_dim.')
  }

  if(psf_fluxcheck > 0.999){
    message('WARNING: psf output image contains more than 99.9% of the total model flux! Consider decreasing the size of psf_dim.')
  }
  highfit$psf = highfit$psf$z / sum(highfit$psf$z)
  highfit$psf_fluxcheck

  highfit$time = (proc.time()[3]-timestart)/60
  highfit$date = date()
  highfit$call = call
  highfit$ProFit.version = packageVersion('ProFit')
  highfit$ProFound.version = packageVersion('ProFound')
  highfit$Highlander.version = packageVersion('Highlander')
  highfit$R.version = R.version

  return(highfit)
}

############################################ profuseFound2Fit.R ############################################

profuseFound2Fit = function(image,
                           sigma = NULL,
                           loc = NULL,
                           segim = NULL,
                           Ncomp = 2,
                           cutbox = dim(image),
                           psf = NULL,
                           magdiff = 2.5,
                           magzero = 0,
                           gain = NULL,
                           resamp = NULL,
                           loc_use = FALSE,
                           loc_fit = TRUE,
                           mag_fit = TRUE,
                           sing_nser = 2,
                           bulge_nser = 4,
                           disk_nser = 1,
                           sing_nser_fit = TRUE,
                           bulge_nser_fit = FALSE,
                           disk_nser_fit = FALSE,
                           bulge_circ = TRUE,
                           nser_upper=5.3,
                           star_rough = TRUE,
                           fit_rough = FALSE,
                           psf_dim = c(51, 51),
                           star_con = 2,
                           star_con_fit = TRUE,
                           star_circ = TRUE,
                           offset = NULL,
                           tightcrop = TRUE,
                           deblend_extra = TRUE,
                           fit_extra = FALSE,
                           pos_delta=10,
                           ...) {
  if(Ncomp == 0.5){psf = NULL}
  if((Ncomp >= 1 & is.null(psf)) | (Ncomp == 0 & is.null(psf))){
    message('Need PSF for Ncomp >= 1 or Ncomp == 0. Running AllStarDoFit!')
    psf = profuseAllStarDoFit(image = image,
                              resamp = resamp,
                              psf_dim = psf_dim,
                              star_con = star_con,
                              star_con_fit = star_con_fit,
                              star_circ = star_circ,
                              magzero = magzero,
                              rough = star_rough,
                              skycut = 2, #works well for stars
                              SBdilate = 2)$psf #works well for stars
  }

  image[image < quantile(image, 0.01, na.rm=TRUE)*10] = NA
  #image[image > quantile(image, 0.99, na.rm=TRUE)*10] = NA

  if(!is.null(loc)){
    cutim = magicaxis::magcutout(image, loc = loc, box = cutbox)
    loc_cut = cutim$loc
    cutim = cutim$image
  }else{
    loc_cut = dim(image) / 2
    cutim = image

  }

  if (!is.null(sigma) & !is.null(loc)) {
    cutsigma = magicaxis::magcutout(sigma, loc = loc, box = cutbox)$image
  } else{
    cutsigma = NULL
  }

  if(is.null(segim)){
    cutseg = NULL
  }else if(dim(segim)[1] == dim(cutim)[1] & dim(segim)[2] == dim(cutim)[2]){
    cutseg = segim
  }else if (dim(segim)[1] == dim(image)[1] & dim(segim)[2] == dim(image)[2] & !is.null(loc)) {
    cutseg = magicaxis::magcutout(segim, loc = loc, box = cutbox)$image
  }else{
    message('No input segim that matches the input image- will create one using ProFound!')
    cutseg = NULL
  }

  loc = loc_cut

  message('    Running ProFound')
  if(!requireNamespace("ProFound", quietly = TRUE)){stop('The ProFound package is required to run this function!')}

  if(is.null(cutseg)){
    mini_profound = ProFound::profoundProFound(
      image = cutim,
      sky = 0,
      redosky = FALSE,
      nearstats = TRUE,
      groupby = 'segim',
      magzero = magzero,
      verbose = FALSE,
      ...
    )
    cutseg = mini_profound$segim
  }else{
    mini_profound = ProFound::profoundProFound(
      image = cutim,
      segim = cutseg,
      sky = 0,
      redosky = FALSE,
      nearstats = TRUE,
      groupby = 'segim',
      magzero = magzero,
      verbose = FALSE,
      iters = 0,
      ...
    )
  }

  if (is.null(cutsigma)) {
    if(is.null(gain)){
      gain = ProFound::profoundGainEst(cutim, objects = mini_profound$objects, sky = 0)
    }
    cutsigma = ProFound::profoundMakeSigma(
      image = cutim,
      objects = mini_profound$objects,
      gain = gain,
      sky = 0,
      skyRMS = mini_profound$skyRMS,
      plot = FALSE
    )
  }

  if(!is.null(resamp)){
    if(resamp > 1){
      cutsigma = cutsigma*resamp
    }
  }

  if(deblend_extra & fit_extra==FALSE){
    cutim = ProFound::profoundFluxDeblend(mini_profound, image_reweight=TRUE)$image
  }

  segID_tar = mini_profound$segim[cutbox[1] / 2, cutbox[2] / 2]
  if (segID_tar == 0) {
    message('Target appears to be sky! Consider using different ProFound parameters.')
    return(NULL)
  }
  loc_tar = which(mini_profound$segstats$segID == segID_tar)
  magID_tar = mini_profound$segstats[mini_profound$segstats$segID == segID_tar, 'mag']

  if(fit_extra){
    segID_ext = unlist(mini_profound$near$nearID[mini_profound$near$segID == segID_tar])
    if (length(segID_ext) > 0) {
      loc_ext = match(segID_ext, mini_profound$segstats$segID)
      loc_ext = loc_ext[which(mini_profound$segstats[loc_ext, "mag"] < magID_tar + magdiff)]
      segID_ext = mini_profound$segstats[loc_ext, 'segID']
      N_ext = length(loc_ext)
    } else{
      N_ext = 0
    }
  }else{
    segID_ext = {}
    N_ext = 0
  }

  region = matrix(mini_profound$segim %in% c(segID_tar, segID_ext), nrow=cutbox[1], ncol=cutbox[2])
  regionlim = which(region, arr.ind=TRUE)

  if(tightcrop){
    xlo = min(regionlim[,1])
    xhi = max(regionlim[,1])
    ylo = min(regionlim[,2])
    yhi = max(regionlim[,2])

    cutim = cutim[xlo:xhi, ylo:yhi]
    region = region[xlo:xhi, ylo:yhi]
    cutsigma = cutsigma[xlo:xhi, ylo:yhi]

    if(loc_use){
      xcen = loc[1] - xlo + 1
      ycen = loc[2] - ylo + 1
      xcen_int = xcen + c(-pos_delta,pos_delta)
      ycen_int = ycen + c(-pos_delta,pos_delta)
    }else{
      xcen = mini_profound$segstats[loc_tar, 'xmax'] - xlo + 1
      ycen = mini_profound$segstats[loc_tar, 'ymax'] - ylo + 1
      xcen_int = xcen + c(-pos_delta,pos_delta)
      ycen_int = ycen + c(-pos_delta,pos_delta)
    }
  }else{
    xlo = 1L
    ylo = 1L
    xcen = mini_profound$segstats[loc_tar, 'xmax']
    ycen = mini_profound$segstats[loc_tar, 'ymax']
    xcen_int = xcen + c(-pos_delta,pos_delta)
    ycen_int = ycen + c(-pos_delta,pos_delta)
  }

  if (Ncomp == 0) {
    modellist = list(
      pointsource = list(
        xcen = mini_profound$segstats[loc_tar, 'xmax'],
        ycen = mini_profound$segstats[loc_tar, 'ymax'],
        mag = mini_profound$segstats[loc_tar, 'mag']
      )
    )
  }else if(Ncomp == 0.5) {
    if(star_circ){
      ang = 0
      axrat = 1
    }else{
      ang = mini_profound$segstats[loc_tar, 'ang']
      axrat = mini_profound$segstats[loc_tar, 'axrat']
    }

    modellist = list(
      moffat = list(
          xcen = xcen,
          ycen = ycen,
          mag = mini_profound$segstats[loc_tar, 'mag'],
          fwhm = mini_profound$segstats[loc_tar, 'R50'] * 2,
          con = star_con,
          ang = ang,
          axrat = axrat
      )
    )
  } else if (Ncomp == 1) {
    modellist = list(
      sersic = list(
        xcen = xcen,
        ycen = ycen,
        mag = mini_profound$segstats[loc_tar, 'mag'],
        re = mini_profound$segstats[loc_tar, 'R50'],
        nser = sing_nser,
        ang = mini_profound$segstats[loc_tar, 'ang'],
        axrat = mini_profound$segstats[loc_tar, 'axrat']
      )
    )
  } else if (Ncomp == 2) {
    modellist = list(
      sersic = list(
        xcen = rep(xcen, 2),
        ycen = rep(ycen, 2),
        mag = rep(mini_profound$segstats[loc_tar, 'mag'], 2) + 0.752575,
        re = mini_profound$segstats[loc_tar, 'R50'] * c(0.5, 1.5),
        nser = c(bulge_nser, disk_nser),
        ang = c(ifelse(bulge_circ, 0, mini_profound$segstats[loc_tar, 'ang']), mini_profound$segstats[loc_tar, 'ang']),
        axrat = c(1, mini_profound$segstats[loc_tar, 'axrat'])
      )
    )
  } else if (Ncomp == 1.5) {
    modellist = list(
      pointsource = list(
        xcen = xcen,
        ycen = ycen,
        mag = mini_profound$segstats[loc_tar, 'mag'] + 0.752575
      ),
      sersic = list(
        xcen = xcen,
        ycen = ycen,
        mag = mini_profound$segstats[loc_tar, 'mag'] + 0.752575,
        re = mini_profound$segstats[loc_tar, 'R50'],
        nser = disk_nser,
        ang = mini_profound$segstats[loc_tar, 'ang'],
        axrat = mini_profound$segstats[loc_tar, 'axrat']
      )
    )
  } else if (Ncomp == 3) {
    modellist = list(
      sersic = list(
        xcen = rep(xcen, 3),
        ycen = rep(ycen, 3),
        mag = rep(mini_profound$segstats[loc_tar, 'mag'], 3) + 0.4771213,
        re = mini_profound$segstats[loc_tar, 'R50'] * c(0.5, 2, 1),
        nser = c(bulge_nser, disk_nser, disk_nser),
        ang = c(ifelse(bulge_circ, 0, mini_profound$segstats[loc_tar, 'ang']), rep(mini_profound$segstats[loc_tar, 'ang'],2)),
        axrat = c(1, rep(mini_profound$segstats[loc_tar, 'axrat'],2))
      )
    )
  }

  if (Ncomp == 0) {
    tofit = list(
      pointsource = list(
        xcen = loc_fit,
        ycen = loc_fit,
        mag = mag_fit
      )
    )
    constraints = NULL
  }else if (Ncomp == 0.5) {
    tofit = list(
      moffat = list(
        xcen = loc_fit,
        ycen = loc_fit,
        mag = mag_fit,
        fwhm = TRUE,
        con = star_con_fit,
        ang = !star_circ,
        axrat = !star_circ
      )
    )
    constraints = NULL
  }else if (Ncomp == 1) {
    tofit = list(
      sersic = list(
        xcen = loc_fit,
        ycen = loc_fit,
        mag = mag_fit,
        re = TRUE,
        nser = sing_nser_fit,
        ang = TRUE,
        axrat = TRUE
      )
    )
    constraints = NULL
  } else if (Ncomp == 2) {
    tofit = list(sersic = list(
      xcen = c(loc_fit, NA), #The NA couples the components together
      ycen = c(loc_fit, NA), #The NA couples the components together
      mag = rep(mag_fit, 2),
      re = rep(TRUE, 2),
      nser = c(bulge_nser_fit, disk_nser_fit),
      ang = c(!bulge_circ, TRUE),
      axrat = c(!bulge_circ, TRUE)
    ))
    constraints = NULL
  } else if (Ncomp == 1.5) {
    tofit = list(
      pointsource = list(
        xcen = FALSE,
        ycen = FALSE,
        mag = mag_fit
      ),
      sersic = list(
        xcen = loc_fit,
        ycen = loc_fit,
        mag = mag_fit,
        re = TRUE,
        nser = disk_nser_fit,
        ang = TRUE,
        axrat = TRUE
      )
    )
    constraints = .constraints_func1_5
  } else if (Ncomp == 3) {
    tofit = list(sersic = list(
      xcen = c(loc_fit, NA, NA), #The NA couples the components together
      ycen = c(loc_fit, NA, NA), #The NA couples the components together
      mag = rep(mag_fit, 3),
      re = rep(TRUE, 3),
      nser = c(bulge_nser_fit, disk_nser_fit, NA),
      ang = c(!bulge_circ, TRUE, NA),
      axrat = c(!bulge_circ, TRUE, NA)
    ))
    constraints = NULL
  }

  if (Ncomp == 0) {
    tolog = list(
      pointsource = list(
        xcen = FALSE,
        ycen = FALSE,
        mag = FALSE
      )
    )
    constraints = NULL
  }else if (Ncomp == 0.5) {
    tolog = list(
      moffat = list(
        xcen = rep(FALSE, Ncomp),
        ycen = rep(FALSE, Ncomp),
        mag = rep(FALSE, Ncomp),
        fwhm = rep(TRUE, Ncomp),
        #fwhm is best fit in log space
        con = rep(TRUE, Ncomp),
        #con is best fit in log space
        ang = rep(FALSE, Ncomp),
        axrat = rep(TRUE, Ncomp) #axrat is best fit in log space
      )
    )
  } else if (Ncomp == 1 | Ncomp == 2) {
    tolog = list(
      sersic = list(
        xcen = rep(FALSE, Ncomp),
        ycen = rep(FALSE, Ncomp),
        mag = rep(FALSE, Ncomp),
        re = rep(TRUE, Ncomp),
        #re is best fit in log space
        nser = rep(TRUE, Ncomp),
        #nser is best fit in log space
        ang = rep(FALSE, Ncomp),
        axrat = rep(TRUE, Ncomp) #axrat is best fit in log space
      )
    )
  } else if (Ncomp == 1.5) {
    tolog = list(
      pointsource = list(
        xcen = FALSE,
        ycen = FALSE,
        mag = FALSE
      ),
      sersic = list(
        xcen = FALSE,
        ycen = FALSE,
        mag = FALSE,
        re = TRUE,
        #re is best fit in log space
        nser = TRUE,
        #nser is best fit in log space
        ang = FALSE,
        axrat = TRUE #axrat is best fit in log space
      )
    )
  } else if (Ncomp == 3) {
    tolog = list(
      sersic = list(
        xcen = rep(FALSE, 3),
        ycen = rep(FALSE, 3),
        mag = rep(FALSE, 3),
        re = rep(TRUE, 3),
        #re is best fit in log space
        nser = rep(TRUE, 3),
        #nser is best fit in log space
        ang = rep(FALSE, 3),
        axrat = rep(TRUE, 3) #axrat is best fit in log space
      )
    )
  }

  #maxsize = sqrt(dim(cutim)[1]^2 + dim(cutim)[2]^2)
  maxsize = mini_profound$segstats[loc_tar, 'R50'] * 4

  if (Ncomp == 0) {
    intervals = list(pointsource = list(
      xcen = list(xcen_int),
      ycen = list(ycen_int),
      mag = list(c(0, 40))
    ))
  } else if (Ncomp == 0.5) {
    intervals = list(moffat = list(
      xcen = list(xcen_int),
      ycen = list(ycen_int),
      mag = list(c(0, 40)),
      fwhm = list(c(0.5, maxsize)),
      con = list(c(1, 10)),
      ang = list(c(-180, 360)),
      axrat = list(c(0.5, 1))
    ))
  } else if (Ncomp == 1) {
    intervals = list(sersic = list(
      xcen = list(xcen_int),
      ycen = list(ycen_int),
      mag = list(c(0, 40)),
      re = list(c(1, maxsize)),
      nser = list(c(0.5, nser_upper)),
      ang = list(c(-180, 360)),
      axrat = list(c(0.01, 1))
    ))
  } else if (Ncomp == 2) {
    intervals = list(
      sersic = list(
        xcen = list(xcen_int, xcen_int),
        ycen = list(ycen_int, ycen_int),
        mag = list(c(0, 40), c(0, 40)),
        re = list(c(1, maxsize), c(1, maxsize)),
        nser = list(c(2, nser_upper), c(0.5, 2)),
        ang = list(c(-180, 360), c(-180, 360)),
        axrat = list(c(0.01, 1), c(0.01, 1))
      )
    )
  } else if (Ncomp == 1.5) {
    intervals = list(
      pointsource = list(
        xcen = list(xcen_int),
        ycen = list(ycen_int),
        mag = list(c(0, 40))
      ),
      sersic = list(
        xcen = list(xcen_int),
        ycen = list(ycen_int),
        mag = list(c(0, 40)),
        re = list(c(1, maxsize)),
        nser = list(c(0.5, nser_upper)),
        ang = list(c(-180, 360)),
        axrat = list(c(0.01, 1))
      )
    )
  } else if (Ncomp == 3) {
    intervals = list(
      sersic = list(
        xcen = list(xcen_int, xcen_int, xcen_int),
        ycen = list(ycen_int, ycen_int, ycen_int),
        mag = list(c(0, 40), c(0, 40), c(0, 40)),
        re = list(c(1, maxsize), c(1, maxsize), c(1, maxsize)),
        nser = list(c(2, nser_upper), c(0.5, 2), c(0.5, 2)),
        ang = list(c(-180, 360), c(-180, 360), c(-180, 360)),
        axrat = list(c(0.01, 1), c(0.01, 1), c(0.01, 1))
      )
    )
  }

  if (fit_extra & N_ext > 0) {
    modellist = c(modellist,
                  list(
                    sersic = list(
                      xcen = mini_profound$segstats[loc_ext, 'xmax'] - xlo + 1L,
                      ycen = mini_profound$segstats[loc_ext, 'ymax'] - ylo + 1L,
                      mag = mini_profound$segstats[loc_ext, 'mag'],
                      re = mini_profound$segstats[loc_ext, 'R50'],
                      nser = rep(2, N_ext),
                      ang = mini_profound$segstats[loc_ext, 'ang'],
                      axrat = mini_profound$segstats[loc_ext, 'axrat']
                    )
                  )
                )

    tofit = c(tofit,
              list(
                sersic = list(
                  xcen = rep(FALSE, N_ext),
                  ycen = rep(FALSE, N_ext),
                  mag = rep(TRUE, N_ext),
                  re = rep(TRUE, N_ext),
                  nser = rep(TRUE, N_ext),
                  ang = rep(FALSE, N_ext),
                  axrat = rep(TRUE, N_ext)
                )
              )
            )

    tolog = c(tolog,
              list(
                sersic = list(
                  xcen = rep(FALSE, N_ext),
                  ycen = rep(FALSE, N_ext),
                  mag = rep(FALSE, N_ext),
                  re = rep(TRUE, N_ext),
                  #re is best fit in log space
                  nser = rep(TRUE, N_ext),
                  #nser is best fit in log space
                  ang = rep(FALSE, N_ext),
                  axrat = rep(TRUE, N_ext) #axrat is best fit in log space
                )
              )
            )

    maxsize = max(mini_profound$segstats[loc_ext, 'R50']*4, na.rm=TRUE)

    intervals = c(intervals,
                  list(
                    sersic = list(
                      xcen = rep(list(xcen_int), N_ext),
                      ycen = rep(list(ycen_int), N_ext),
                      mag = rep(list(c(0, 40)), N_ext),
                      re = rep(list(c(1, maxsize)), N_ext),
                      nser = rep(list(c(0.5, nser_upper)), N_ext),
                      ang = rep(list(c(-180, 360)), N_ext),
                      axrat = rep(list(c(0.01, 1)), N_ext)
                    )
                  )
                )
  }

  Data = profitSetupData(
    image = cutim,
    region = region,
    sigma = cutsigma,
    segim = cutseg,
    psf = psf,
    modellist = modellist,
    tofit = tofit,
    tolog = tolog,
    intervals = intervals,
    constraints = constraints,
    magzero = magzero,
    algo.func = 'LD',
    verbose = FALSE,
    offset = offset,
    rough = fit_rough
  )
  Data$Nmod = Ncomp + N_ext
  return(invisible(list(profound = mini_profound, Data = Data)))
}

profuseDoFit = function(image,
                       sigma = NULL,
                       loc = NULL,
                       F2F = NULL,
                       Ncomp = 2,
                       cutbox = dim(image),
                       psf = NULL,
                       magdiff = 2.5,
                       magzero = 0,
                       psf_dim = c(51,51),
                       star_rough = TRUE,
                       fit_rough = FALSE,
                       plot = FALSE,
                       seed = 666,
                       optim_iters = 2,
                       Niters = c(1000,1000),
                       NfinalMCMC = Niters[2],
                       ...) {

  timestart = proc.time()[3] # start timer
  #call = match.call(expand.dots=TRUE)

  if(is.null(F2F)){
    message('Running Found2Fit')
    F2F = profuseFound2Fit(
      image = image,
      sigma = sigma,
      loc = loc,
      Ncomp = Ncomp,
      cutbox = cutbox,
      psf = psf,
      magdiff = magdiff,
      magzero = magzero,
      fit_rough = fit_rough,
      ...
    )
  }else{
    F2F = profuseRegenPSF_F2F(F2F)
  }

  Data = F2F$Data

  if (plot) {
    profitLikeModel(parm = Data$init,
                    Data = Data,
                    makeplots = TRUE)
    legend('topright', legend = 'Start')
  }

  lowers = unlist(Data$intervals)[c(T, F)]
  lowers[unlist(Data$tolog) == T] = log10(lowers[unlist(Data$tolog) == T])
  lowers = lowers[which(unlist(Data$tofit))]
  uppers = unlist(Data$intervals)[c(F, T)]
  uppers[unlist(Data$tolog) == T] = log10(uppers[unlist(Data$tolog) == T])
  uppers = uppers[which(unlist(Data$tofit))]

  message('Running Highlander')
  if(!requireNamespace("ProFound", quietly = TRUE)){stop('The Highander package is required to run this function!')}
  highfit = Highlander::Highlander(
    parm = Data$init,
    Data = Data,
    likefunc = profitLikeModel,
    seed = seed,
    lower = lowers,
    upper = uppers,
    applyintervals = FALSE,
    applyconstraints = FALSE,
    optim_iters = optim_iters,
    Niters = Niters,
    NfinalMCMC = NfinalMCMC,
    parm.names = Data$parm.names
  )
  names(highfit$parm) = names(Data$init)

  if (plot) {
    profitLikeModel(highfit$parm, Data = Data, makeplots = TRUE)
    legend('topright', legend = 'After')
  }

  highfit$profound = F2F$profound
  highfit$Data = Data
  highfit$initmodel = profitRemakeModellist(Data$init, Data = Data)
  highfit$finalmodel = profitRemakeModellist(highfit$parm, Data = Data)

  highfit$parm[highfit$parm < lowers] = lowers[highfit$parm < lowers]
  highfit$parm[highfit$parm > uppers] = uppers[highfit$parm > uppers]

  highfit$error = apply(highfit$LD_last$Posterior1,
                         MARGIN = 2,
                         FUN = 'sd')

  if(Ncomp == 0.5){
    if(psf_dim[1] %% 2 == 0){psf_dim[1] = psf_dim[1] + 1}
    if(psf_dim[2] %% 2 == 0){psf_dim[2] = psf_dim[2] + 1}
    temp_modellist = highfit$finalmodel$modellist[[1]]
    temp_modellist$moffat$mag = 0
    temp_modellist$moffat$xcen = psf_dim[1]/2
    temp_modellist$moffat$ycen = psf_dim[2]/2

    highfit$psf = profitMakeModel(temp_modellist, dim=psf_dim)
    highfit$psf_modellist = temp_modellist

    psf_fluxcheck = sum(highfit$psf$z)

    if(psf_fluxcheck < 0.95){
      message('WARNING: psf output image contains less than 95% of the total model flux! Consider increasing the size of psf_dim.')
    }

    if(psf_fluxcheck > 0.999){
      message('WARNING: psf output image contains more than 99.9% of the total model flux! Consider decreasing the size of psf_dim.')
    }
    highfit$psf = highfit$psf$z / sum(highfit$psf$z)
    highfit$psf_fluxcheck
  }

  highfit$time = (proc.time()[3]-timestart)/60
  highfit$date = date()
  #highfit$call = call
  highfit$ProFit.version = packageVersion('ProFit')
  highfit$ProFound.version = packageVersion('ProFound')
  highfit$Highlander.version = packageVersion('Highlander')
  highfit$R.version = R.version
  highfit$LD_last$Call = NULL
  highfit$LD_last$Model = NULL

  return(highfit)
}

.constraints_func1_5 = function(modellist=NULL) {
  modellist$pointsource$xcen = modellist$sersic$xcen
  modellist$pointsource$ycen = modellist$sersic$ycen
  return(modellist)
}

# .constraints_func3 = function(modellist=NULL) { #I don't think I actually need this, can use NA
#   modellist$sersic$nser[3] = modellist$sersic$nser[2]
#   modellist$sersic$ang[3] = modellist$sersic$ang[2]
#   modellist$sersic$axrat[3] = modellist$sersic$axrat[2]
#   return(modellist)
# }


############################################ profuseMultiBand.R ############################################

profuseMultiBandFound2Fit = function(image_list,
                                     segim_list = NULL,
                                     segim_global = NULL,
                                    sky_list = NULL,
                                    skyRMS_list = NULL,
                                    loc = NULL,
                                    parm_global = c("sersic.xcen1", "sersic.ycen1", "sersic.re1", "sersic.ang2", "sersic.axrat2"),
                                    Ncomp = 2,
                                    cutbox = dim(image),
                                    psf_list = NULL,
                                    magdiff = 2.5,
                                    magzero = NULL,
                                    gain = NULL,
                                    resamp = NULL,
                                    sing_nser = 2,
                                    bulge_nser = 4,
                                    disk_nser = 1,
                                    sing_nser_fit = TRUE,
                                    bulge_nser_fit = FALSE,
                                    disk_nser_fit = FALSE,
                                    bulge_circ =  TRUE,
                                    nser_upper=5.3,
                                    star_rough = TRUE,
                                    fit_rough = FALSE,
                                    psf_dim = c(51, 51),
                                    star_con = 2,
                                    star_con_fit = TRUE,
                                    star_circ = TRUE,
                                    tightcrop = TRUE,
                                    wave = NULL,
                                    smooth.parm = NULL,
                                    parm_ProSpect = NULL,
                                    data_ProSpect = NULL, #perhaps need a way to specify extra data going to bulge/disk. Naming or list?
                                    logged_ProSpect = NULL,
                                    intervals_ProSpect = NULL,
                                    ...){
  Nim = length(image_list)

  if(is.null(magzero)){
    magzero = rep(0, Nim)
  }

  for(i in 1:Nim){
    if(is.null(sky_list[i][[1]]) | is.null(skyRMS_list[i][[1]])){ #[i][[1]] looks silly, but it will return NULL when sky_list = NULL for any i (default). [[i]] will error in this case
      message("Image ",i,": running initial ProFound")
      profound = ProFound::profoundProFound(image = image_list[[i]],
                                            sky = sky_list[i][[1]],
                                            skyRMS = skyRMS_list[i][[1]],
                                            magzero = magzero[i],
                                            ...)
      if(is.null(sky_list[i])){
        image_list[[i]] = image_list[[i]] - profound$sky
      }else{
        image_list[[i]] = image_list[[i]] - sky_list[i]
      }
      if(is.null(skyRMS_list[i][[1]])){
        skyRMS_list[[i]] = profound$skyRMS
      }
    }

    if(is.null(psf_list[i][[1]])){
      message("Image ",i,": running AllStarDoFit")
      psf_list[[i]] = profuseAllStarDoFit(image = image_list[[i]],
                                          resamp = resamp[i][[1]],
                                         psf_dim = psf_dim,
                                         star_con = star_con,
                                         star_con_fit = star_con_fit,
                                         star_circ = star_circ,
                                         magzero = magzero[i],
                                         rough = star_rough,
                                         skycut = 2, #works well for stars
                                         SBdilate = 2)$psf #works well for stars
    }
  }

  message("Making image stack")

  multi_stack = ProFound::profoundMakeStack(
    image_list = image_list,
    skyRMS_list = skyRMS_list,
    magzero_in = magzero,
    magzero_out = 0
  )

  message("Running ProFound on stack")

  if(is.null(segim_global) & !is.null(segim_list)){
    segim_global = profuseSegimGlobal(segim_list)
  }

  multi_stack_pro = ProFound::profoundProFound(image=multi_stack$image,
                                               segim=segim_global,
                                               sky=0,
                                               skyRMS=multi_stack$skyRMS,
                                               redosky=FALSE,
                                               static_photom=TRUE,
                                               ...)

  message("Running Found2Fit on stack")

  F2Fstack = profuseFound2Fit(image = multi_stack$image,
                             sigma = multi_stack$skyRMS, #not quite a sigma map, but doesn't matter for the stack F2F
                             loc = loc,
                             segim = multi_stack_pro$segim,
                             Ncomp = Ncomp,
                             psf = matrix(1,1,1), #Doesn't matter what we pass in here
                             magzero = 0,
                             mag_fit = is.null(parm_ProSpect),
                             sing_nser = sing_nser,
                             bulge_nser = bulge_nser,
                             disk_nser = disk_nser,
                             sing_nser_fit = sing_nser_fit,
                             bulge_nser_fit = bulge_nser_fit,
                             disk_nser_fit = disk_nser_fit,
                             bulge_circ =  bulge_circ,
                             nser_upper = nser_upper,
                             tightcrop = FALSE,
                             fit_extra = FALSE
  )

  if(tightcrop){
    regionlim = which(F2Fstack$Data$region, arr.ind=TRUE)

    xlo = min(regionlim[,1])
    xhi = max(regionlim[,1])
    ylo = min(regionlim[,2])
    yhi = max(regionlim[,2])

    if(!is.null(F2Fstack$Data$modellist$sersic)){
      F2Fstack$Data$modellist$sersic$xcen = F2Fstack$Data$modellist$sersic$xcen - xlo + 1
      F2Fstack$Data$modellist$sersic$ycen = F2Fstack$Data$modellist$sersic$ycen - ylo + 1

      for(i in 1:length(F2Fstack$Data$intervals$sersic$xcen)){
        F2Fstack$Data$intervals$sersic$xcen[[i]] = F2Fstack$Data$intervals$sersic$xcen[[i]] - xlo + 1
        F2Fstack$Data$intervals$sersic$ycen[[i]] = F2Fstack$Data$intervals$sersic$ycen[[i]] - ylo + 1
      }
    }

    if(!is.null(F2Fstack$Data$modellist$moffat)){
      F2Fstack$Data$modellist$moffat$xcen = F2Fstack$Data$modellist$moffat$xcen - xlo + 1
      F2Fstack$Data$modellist$moffat$ycen = F2Fstack$Data$modellist$moffat$ycen - ylo + 1

      for(i in 1:length(F2Fstack$Data$intervals$moffat$xcen)){
        F2Fstack$Data$intervals$moffat$xcen[[i]] = F2Fstack$Data$intervals$moffat$xcen[[i]] - xlo + 1
        F2Fstack$Data$intervals$moffat$ycen[[i]] = F2Fstack$Data$intervals$moffat$ycen[[i]] - ylo + 1
      }
    }

    if(!is.null(F2Fstack$Data$modellist$pointsource)){
      F2Fstack$Data$modellist$pointsource$xcen = F2Fstack$Data$modellist$pointsource$xcen - xlo + 1
      F2Fstack$Data$modellist$pointsource$ycen = F2Fstack$Data$modellist$pointsource$ycen - ylo + 1

      for(i in 1:length(F2Fstack$Data$intervals$pointsource$xcen)){
        F2Fstack$Data$intervals$pointsource$xcen[[i]] = F2Fstack$Data$intervals$pointsource$xcen[[i]] - xlo + 1
        F2Fstack$Data$intervals$pointsource$ycen[[i]] = F2Fstack$Data$intervals$pointsource$ycen[[i]] - ylo + 1
      }
    }

    xcenloc = grep('xcen', F2Fstack$Data$parm.names)
    if(length(xcenloc) > 0){
      F2Fstack$Data$init[xcenloc] = F2Fstack$Data$init[xcenloc] - xlo + 1
    }

    ycenloc = grep('ycen', F2Fstack$Data$parm.names)
    if(length(ycenloc) > 0){
      F2Fstack$Data$init[ycenloc] = F2Fstack$Data$init[ycenloc] - ylo + 1
    }

  }else{
    xlo = 1L
    ylo = 1L
    xhi = dim(multi_stack$image)[1]
    yhi = dim(multi_stack$image)[2]
  }

  if(!is.null(parm_ProSpect)){
    F2Fstack$Data$tofit$sersic$mag[] = FALSE
  }

  mag_stack = ProFound::profoundFlux2Mag(flux=sum(multi_stack$image[F2Fstack$Data$region], na.rm=TRUE), magzero=0)

  MF2F = list()

  for(i in 1:Nim){
    if(is.null(gain[[i]])){
     gain_loc = ProFound::profoundGainEst(image_list[[i]], objects=multi_stack_pro$objects, sky = 0)
    }else{
      gain_loc = gain[[i]]
    }
    sigma = ProFound::profoundMakeSigma(
      image = image_list[[i]],
      objects = multi_stack_pro$objects,
      gain = gain_loc,
      sky = 0,
      skyRMS = skyRMS_list[[i]],
      plot = FALSE
    )

    if(!is.null(resamp[i][[1]])){
      if(resamp[i][[1]] > 1){
        sigma = sigma*resamp[i][[1]]
      }
    }

    if(is.null(segim_list[i][[1]])){
      segim_use = F2Fstack$Data$segim
    }else{
      segim_use = segim_list[i][[1]]
    }

    message("Image ",i,": running SetupData")

    region = (segim_use == which.max(tabulate(segim_use[F2Fstack$Data$region])))

    MF2F[[i]] = profitSetupData(
      image = image_list[[i]][xlo:xhi,ylo:yhi],
      region = region[xlo:xhi,ylo:yhi],
      sigma = sigma[xlo:xhi,ylo:yhi],
      segim = segim_use[xlo:xhi,ylo:yhi],
      psf = psf_list[[i]],
      modellist = F2Fstack$Data$modellist,
      tofit = F2Fstack$Data$tofit,
      tolog = F2Fstack$Data$tolog,
      intervals = F2Fstack$Data$intervals,
      constraints = F2Fstack$Data$constraints,
      magzero = magzero[i],
      algo.func = 'LD',
      verbose = FALSE,
      rough = fit_rough
    )
  }

  names(MF2F) = names(image_list)

  #Check below for ProFuse- how this all works could be complicated... see also profitLikeModel and profitRemakeModelList
  if(is.null(parm_global)){
    parm = F2Fstack$Data$init
    for(i in 1:Nim){
      MF2F[[i]]$parmuse = 1:length(parm)
    }
  }else{
    parm_init = F2Fstack$Data$init
    if(is.character(parm_global)){
      parm_global = match(parm_global, names(parm_init))
    }
    parm = parm_init[parm_global]
    Nparm = length(parm_init)
    parm_local = 1:Nparm
    parm_local = parm_local[-parm_global] #anything not global is local - check for ProFuse
    for(i in 1:Nim){
      if(length(parm_local) > 0){
        parm_temp = F2Fstack$Data$init[parm_local] #extract local - check for ProFuse
        names(parm_temp) = paste0(names(parm_temp),'_',i) #mod names for local by adding the band number - check for ProFuse
        parm = c(parm, parm_temp) #create parent parm object (with all unique parm_local vector appended) - check for ProFuse
        parmuse = 1:Nparm
        parmuse[parm_global] = 1:length(parm_global)
        parmuse[parm_local] = length(parm_global) + 1:length(parm_local) + (i-1)*length(parm_local) #define parmuse location in parent parm object - check for ProFuse
        MF2F[[i]]$parmuse = parmuse
      }else{
        MF2F[[i]]$parmuse = 1:Nparm
      }
    }
  }

#Ideas for ProFuse. Need to calculate the parmuse positions as per they will be after all the ProSpect related parameters are removed and the relevant per band magnitudes added and named. This probably means we need a parm_ProSpect object that exhaustively identifies all of these (and need to enforce them being at the end). Probably cannot make it any more flexible just to be safe. In principle then we just need remove these ProSpect related arguments and replace that part of the parm with the mag inside profitLikeModel, which we then pass into profitRemakeModelList to make the target model images.

  #Create mag offsets based on magzero points and average mag of the stack.
  for(i in 1:Nim){
    mag_image = ProFound::profoundFlux2Mag(flux=sum(image_list[[i]][F2Fstack$Data$region], na.rm=TRUE), magzero=magzero[i])
    mag_diff = mag_stack - mag_image
    sel = grep(paste0('.*mag.*\\_',i), names(parm))
    parm[sel] = parm[sel] - mag_diff
  }

  if(is.null(data_ProSpect$LumDist_Mpc)){
    data_ProSpect$LumDist_Mpc = cosdistLumDist(data_ProSpect$z, H0 = 67.8, OmegaM = 0.308)
  }

  if(is.null(data_ProSpect$magemax)){
    if(!is.null(data_ProSpect$agemax)){
      data_ProSpect$magemax = data_ProSpect$agemax/1e9
    }
  }

  if(is.null(data_ProSpect$Zagemax)){
    if(!is.null(data_ProSpect$agemax)){
      data_ProSpect$Zagemax = data_ProSpect$agemax/1e9
    }
  }

  MF2F$init = c(parm, unlist(parm_ProSpect))
  MF2F$parm.names = names(MF2F$init)
  MF2F$mon.names = F2Fstack$Data$mon.names
  MF2F$Nim = Nim #Number of images
  MF2F$Ncomp = ceiling(Ncomp) #Number of components 0.5 -> 1 and 1.5 -> 2
  MF2F$N = F2Fstack$Data$N #This is the number of fitting pixels (cannot rename)
  MF2F$wave = wave
  MF2F$smooth.parm = smooth.parm
  MF2F$parm_ProSpect = parm_ProSpect
  MF2F$data_ProSpect = data_ProSpect
  MF2F$logged_ProSpect = logged_ProSpect
  MF2F$intervals_ProSpect = intervals_ProSpect

  return(MF2F)
}

profuseMultiBandDoFit = function(image_list,
                                 segim_list = NULL,
                                 segim_global = NULL,
                                sky_list = NULL,
                                skyRMS_list = NULL,
                                loc = NULL,
                                MF2F = NULL,
                                parm_global = c("sersic.xcen1", "sersic.ycen1", "sersic.re1", "sersic.ang2", "sersic.axrat2"),
                                Ncomp = 2,
                                cutbox = dim(image),
                                psf_list = NULL,
                                magzero = NULL,
                                psf_dim = c(51,51),
                                star_rough = TRUE,
                                fit_rough = FALSE,
                                seed = 666,
                                optim_iters = 2,
                                Niters = c(1000,1000),
                                NfinalMCMC = Niters[2],
                                ...) {

  timestart = proc.time()[3] # start timer
  #call = match.call(expand.dots=TRUE)

  if(is.null(MF2F)){
    message('Running MultiBandFound2Fit')
    MF2F = profuseMultiBandFound2Fit(
      image_list = image_list,
      segim_list = segim_list,
      segim_global = segim_global,
      sky_list = sky_list,
      skyRMS_list = skyRMS_list,
      loc = loc,
      parm_global = parm_global,
      Ncomp = Ncomp,
      cutbox = cutbox,
      psf_list = psf_list,
      magzero = magzero,
      star_rough = star_rough,
      fit_rough = fit_rough,
      ...
    )
  }else{
    MF2F = profuseRegenPSF_MF2F(MF2F)
  }

  lower_profit = {}
  upper_profit = {}
  logged_profit = {}

  for(i in 1:length(MF2F[[1]]$intervals)){ #loop over profiles
    for(j in 1:length(MF2F[[1]]$intervals[[i]])){ #loop over parameters
      for(k in 1:length(MF2F[[1]]$intervals[[i]][[1]])){ #loop over components
        if(isTRUE(MF2F[[1]]$tofit[[i]][[j]][[k]])){
          lower_profit = c(lower_profit, MF2F[[1]]$intervals[[i]][[j]][[k]][1])
          upper_profit = c(upper_profit, MF2F[[1]]$intervals[[i]][[j]][[k]][2])
          logged_profit = c(logged_profit, MF2F[[1]]$tolog[[i]][[j]][[k]])
        }
      }
    }
  }

  lower_profit[logged_profit] = log10(lower_profit[logged_profit])
  upper_profit[logged_profit] = log10(upper_profit[logged_profit])

  lower = c(lower_profit, MF2F$intervals_ProSpect$lo)
  upper = c(upper_profit, MF2F$intervals_ProSpect$hi)

  message('Running Highander on multi-band data')
  if(!requireNamespace("ProFound", quietly = TRUE)){stop('The Highander package is required to run this function!')}
  highfit = Highlander::Highlander(
    parm = MF2F$init,
    Data = MF2F,
    likefunc = profitLikeModel,
    seed = seed,
    lower = lower,
    upper = upper,
    applyintervals = TRUE,
    applyconstraints = FALSE,
    optim_iters = optim_iters,
    Niters = Niters,
    NfinalMCMC = NfinalMCMC,
    parm.names = MF2F$parm.names
  )

  highfit$MF2F = MF2F
  highfit$error = apply(highfit$LD_last$Posterior1,
                        MARGIN = 2,
                        FUN = 'sd')

  if(!is.null(MF2F$smooth.parm) & !is.null(MF2F$wave)){
    namevec = names(MF2F$smooth.parm)
    highfit$parm_smooth = highfit$parm
    for(i in 1:length(MF2F$smooth.parm)){
      highfit$parm_smooth = .smooth_parm(parm=highfit$parm_smooth, MF2F$parm.names, extract=namevec[i], wave=MF2F$wave, func=MF2F$smooth.parm[[i]])
    }
  }else{
    highfit$parm_smooth = NULL
  }

  highfit$time = (proc.time()[3]-timestart)/60
  highfit$date = date()
  #highfit$call = call
  highfit$ProFit.version = packageVersion('ProFit')
  highfit$ProFound.version = packageVersion('ProFound')
  highfit$Highlander.version = packageVersion('Highlander')
  highfit$R.version = R.version
  highfit$LD_last$Call = NULL
  highfit$LD_last$Model = NULL
  highfit$RedChi2 = highfit$RedChi2/MF2F$Nim

  class(highfit) = 'profusemulti'

  return(highfit)
}

.smooth_parm = function(parm, parm.names, extract='mag1', wave, func=smooth.spline){
  parm_loc = grep(extract,parm.names)

  parm[parm_loc] = func(log(wave),parm[parm_loc])$y
  return(parm)
}

profuseSegimGlobal = function(segim_list){
  if(requireNamespace("imager", quietly = TRUE)){
    i=NULL
    Csegim_list = foreach(i=1:length(segim_list))%do%{
      temp=segim_list[[i]]
      temp[temp==0L]=NA
      imager::as.cimg(temp)
    }
    Csegim_list = imager::as.imlist(Csegim_list)
    segim_global = imager::parmed(Csegim_list, na.rm=TRUE)
    segim_global = as.matrix(segim_global)
    segim_global[is.na(segim_global)] = 0L
    return(segim_global)
  }else{
    stop('The imager package is needed for segim merging to work. Please install from CRAN.', call. = FALSE)
  }
}


############################################ profuseRegen.R ############################################

profuseRegenPSF_F2F = function(F2F){
  F2F$convopt = list(convolver = profitMakeConvolver('brute', dim(F2F$image), F2F$psf), openclenv=NULL)
  return(F2F)
}

profuseRegenPSF_MF2F = function(MF2F){
  for(i in 1:MF2F$Nim){
    MF2F[[i]]$convopt = list(convolver = profitMakeConvolver('brute', dim(MF2F[[i]]$image), MF2F[[i]]$psf), openclenv=NULL)
  }
  return(MF2F)
}
