# Selection of some functions from the magicaxis package, used for downloading SDSS images and cutting them

magWCSradec2xy=function (RA, Dec, header, CRVAL1 = 0, CRVAL2 = 0, CRPIX1 = 0, 
                         CRPIX2 = 0, CD1_1 = 1, CD1_2 = 0, CD2_1 = 0, CD2_2 = 1, CTYPE1 = "RA--TAN", 
                         CTYPE2 = "DEC--TAN", loc.diff = c(0, 0), coord.type = "deg", 
                         sep = ":") 
{
  if (length(dim(RA)) == 2) {
    Dec = RA[, 2]
    RA = RA[, 1]
  }
  if (coord.type == "sex") {
    RA = hms2deg(RA, sep = sep)
    Dec = dms2deg(Dec, sep = sep)
  }
  RA = as.numeric(RA)
  Dec = as.numeric(Dec)
  tempxy = radec2xy(RA, Dec, header = header, CRVAL1 = CRVAL1, 
                    CRVAL2 = CRVAL2, CRPIX1 = CRPIX1, CRPIX2 = CRPIX2, CD1_1 = CD1_1, 
                    CD1_2 = CD1_2, CD2_1 = CD2_1, CD2_2 = CD2_2, CTYPE1 = CTYPE1, 
                    CTYPE2 = CTYPE2)
  tempxy[, 1] = tempxy[, 1] - 0.5 - loc.diff[1]
  tempxy[, 2] = tempxy[, 2] - 0.5 - loc.diff[2]
  return = tempxy
}

magWCSxy2radec=function (x, y, header, CRVAL1 = 0, CRVAL2 = 0, CRPIX1 = 0, CRPIX2 = 0, 
                         CD1_1 = 1, CD1_2 = 0, CD2_1 = 0, CD2_2 = 1, CTYPE1 = "RA--TAN", 
                         CTYPE2 = "DEC--TAN", loc.diff = c(0, 0)) 
{
  if (length(dim(x)) == 2) {
    y = x[, 2]
    x = x[, 1]
  }
  x = as.numeric(x)
  y = as.numeric(y)
  tempradec = xy2radec(x + 0.5 + loc.diff[1], y + 0.5 + loc.diff[2], 
                       header = header, CRVAL1 = CRVAL1, CRVAL2 = CRVAL2, CRPIX1 = CRPIX1, 
                       CRPIX2 = CRPIX2, CD1_1 = CD1_1, CD1_2 = CD1_2, CD2_1 = CD2_1, 
                       CD2_2 = CD2_2, CTYPE1 = CTYPE1, CTYPE2 = CTYPE2)
  tempradec[, 1] = tempradec[, 1]
  tempradec[, 2] = tempradec[, 2]
  return = tempradec
}

magcutout=function (image, loc = dim(image)/2, box = c(100, 100), shiftloc = FALSE, 
                    paddim = TRUE, plot = FALSE, ...) 
{
  loc = as.numeric(loc)
  xcen = loc[1]
  ycen = loc[2]
  loc = ceiling(loc)
  if (length(box) == 1) {
    box = rep(box, 2)
  }
  xlo = ceiling(loc[1] - (box[1]/2 - 0.5))
  xhi = ceiling(loc[1] + (box[1]/2 - 0.5))
  ylo = ceiling(loc[2] - (box[2]/2 - 0.5))
  yhi = ceiling(loc[2] + (box[2]/2 - 0.5))
  loc.diff = c(x = xlo - 1, y = ylo - 1)
  expand = paddim && shiftloc
  diffxlo = xlo - 1
  if (diffxlo < 0) {
    xlo = 1
    if (expand) 
      xhi = xlo + (box[1] - 1)
  }
  diffxhi = xhi - dim(image)[1]
  if (diffxhi > 0) {
    xhi = dim(image)[1]
    if (expand) {
      xlo = xlo - diffxhi
      if (xlo < 1) 
        xlo = 1
    }
  }
  diffylo = ylo - 1
  if (diffylo < 0) {
    ylo = 1
    if (expand) 
      yhi = ylo + (box[2] - 1)
  }
  diffyhi = yhi - dim(image)[2]
  if (diffyhi > 0) {
    yhi = dim(image)[2]
    if (expand) {
      ylo = ylo - diffyhi
      if (ylo < 1) 
        ylo = 1
    }
  }
  if (!paddim && !shiftloc) {
    if (diffxlo < 0 && (-diffxlo > diffxhi)) 
      xhi = xhi - max(diffxhi, 0) + diffxlo
    if (diffxhi > 0 && (-diffxlo < diffxhi)) 
      xlo = xlo + diffxhi - min(diffxlo, 0)
    if (diffylo < 0 && (-diffylo > diffyhi)) 
      yhi = yhi - max(diffyhi, 0) + diffylo
    if (diffyhi > 0 && (-diffylo < diffyhi)) 
      ylo = ylo + diffyhi - min(diffylo, 0)
  }
  xsel = as.integer(xlo:xhi)
  ysel = as.integer(ylo:yhi)
  xsel = xsel[xsel > 0]
  ysel = ysel[ysel > 0]
  if (length(xsel) == 0 | length(ysel) == 0) {
    image = matrix(NA, box[1], box[2])
  }
  else {
    image = image[xsel, ysel]
    if (paddim && !shiftloc && any(c(diffxlo, -diffxhi, diffylo, 
                                     -diffyhi) < 0)) {
      padded = matrix(NA, box[1], box[2])
      padded[xsel - diffxlo, ysel - diffylo] = image
      image = padded
    }
  }
  if (paddim & shiftloc == FALSE) {
    loc = c(x = xcen - diffxlo, y = ycen - diffylo)
  }
  else {
    loc = c(x = xcen - xlo + 1, y = ycen - ylo + 1)
  }
  loc.orig = c(x = xcen, y = ycen)
  loc.diff = c(x = loc.orig[1] - loc[1], y = loc.orig[2] - 
                 loc[2])
  output = list(image = image, loc = loc, loc.orig = loc.orig, 
                loc.diff = loc.diff, xsel = xsel, ysel = ysel)
  if (plot) {
    if (all(is.na(image))) {
      image[] = 0
      magimage(image, ...)
    }
    else {
      magimage(image, ...)
    }
  }
  invisible(output)
}

magcutoutWCS=function(image, header, loc, box = c(100, 100), shiftloc = FALSE, 
                      paddim = TRUE, plot = FALSE, CRVAL1 = 0, CRVAL2 = 0, CRPIX1 = 0, 
                      CRPIX2 = 0, CD1_1 = 1, CD1_2 = 0, CD2_1 = 0, CD2_2 = 1, coord.type = "deg", 
                      sep = ":", loc.type = c("coord", "coord"), approx.map = FALSE, 
                       ...) 
{
  if (length(loc.type) == 1) {
    loc.type = rep(loc.type, 2)
  }
  if (length(box) == 1) {
    box = rep(box, 2)
  }
  if (!missing(image)) {
    if (any(names(image) == "imDat") & missing(header)) {
      imtype = "FITSio"
      header = image$hdr
      image = image$imDat
    }
    if (any(names(image) == "dat") & missing(header)) {
      imtype = "astro"
      header = image$hdr[[1]]
      header = data.frame(key = header[, 1], value = header[, 
                                                            2], stringsAsFactors = FALSE)
      image = image$dat[[1]]
    }
    if (any(names(image) == "image") & missing(header)) {
      header = image$header
      image = image$image
      if (is.matrix(header) | is.data.frame(header)) {
        imtype = "astro"
      }
      else {
        imtype = "FITSio"
      }
    }
    if (!missing(header)) {
      if (is.matrix(header) | is.data.frame(header)) {
        imtype = "astro"
      }
      else {
        imtype = "FITSio"
      }
    }
  }
  if (missing(loc)) {
    loc = magWCSxy2radec(dim(image)[1]/2, dim(image)[2]/2, 
                         header = header, CRVAL1 = CRVAL1, CRVAL2 = CRVAL2, 
                         CRPIX1 = CRPIX1, CRPIX2 = CRPIX2, CD1_1 = CD1_1, 
                         CD1_2 = CD1_2, CD2_1 = CD2_1, CD2_2 = CD2_2)[1, ]
    tempxy = cbind(dim(image)[1]/2, dim(image)[2]/2)
  }
  else {
    if (loc.type[1] == "coord") {
      if (coord.type == "sex") {
        loc[1] = hms2deg(loc[1], sep = sep)
        loc[2] = dms2deg(loc[2], sep = sep)
      }
      loc = as.numeric(loc)
      tempxy = magWCSradec2xy(loc[1], loc[2], header = header, 
                              CRVAL1 = CRVAL1, CRVAL2 = CRVAL2, CRPIX1 = CRPIX1, 
                              CRPIX2 = CRPIX2, CD1_1 = CD1_1, CD1_2 = CD1_2, 
                              CD2_1 = CD2_1, CD2_2 = CD2_2)
    }
    else if (loc.type[1] == "image") {
      tempxy = rbind(loc)
      loc = magWCSxy2radec(loc[1], loc[2], header = header, 
                           CRVAL1 = CRVAL1, CRVAL2 = CRVAL2, CRPIX1 = CRPIX1, 
                           CRPIX2 = CRPIX2, CD1_1 = CD1_1, CD1_2 = CD1_2, 
                           CD2_1 = CD2_1, CD2_2 = CD2_2)[1, ]
    }
  }
  xcen = tempxy[1, 1]
  ycen = tempxy[1, 2]
  if (loc.type[2] == "coord") {
    box = box/3600
    tempxy = magWCSradec2xy(loc[1] - box[1]/2/cos(loc[2] * 
                                                    pi/180), loc[2], header = header, CRVAL1 = CRVAL1, 
                            CRVAL2 = CRVAL2, CRPIX1 = CRPIX1, CRPIX2 = CRPIX2, 
                            CD1_1 = CD1_1, CD1_2 = CD1_2, CD2_1 = CD2_1, CD2_2 = CD2_2)
    xlo = xcen - sqrt((tempxy[1, 1] - xcen)^2 + (tempxy[1, 
                                                        2] - ycen)^2)
    tempxy = magWCSradec2xy(loc[1] + box[1]/2/cos(loc[2] * 
                                                    pi/180), loc[2], header = header, CRVAL1 = CRVAL1, 
                            CRVAL2 = CRVAL2, CRPIX1 = CRPIX1, CRPIX2 = CRPIX2, 
                            CD1_1 = CD1_1, CD1_2 = CD1_2, CD2_1 = CD2_1, CD2_2 = CD2_2)
    xhi = xcen + sqrt((tempxy[1, 1] - xcen)^2 + (tempxy[1, 
                                                        2] - ycen)^2)
    tempxy = magWCSradec2xy(loc[1], loc[2] - box[2]/2, header = header, 
                            CRVAL1 = CRVAL1, CRVAL2 = CRVAL2, CRPIX1 = CRPIX1, 
                            CRPIX2 = CRPIX2, CD1_1 = CD1_1, CD1_2 = CD1_2, CD2_1 = CD2_1, 
                            CD2_2 = CD2_2)
    ylo = ycen - sqrt((tempxy[1, 1] - xcen)^2 + (tempxy[1, 
                                                        2] - ycen)^2)
    tempxy = magWCSradec2xy(loc[1], loc[2] + box[2]/2, header = header, 
                            CRVAL1 = CRVAL1, CRVAL2 = CRVAL2, CRPIX1 = CRPIX1, 
                            CRPIX2 = CRPIX2, CD1_1 = CD1_1, CD1_2 = CD1_2, CD2_1 = CD2_1, 
                            CD2_2 = CD2_2)
    yhi = ycen + sqrt((tempxy[1, 1] - xcen)^2 + (tempxy[1, 
                                                        2] - ycen)^2)
    xtemp = sort(c(xlo, xhi))
    xlo = ceiling(xtemp[1])
    xhi = ceiling(xtemp[2])
    ytemp = sort(c(ylo, yhi))
    ylo = ceiling(ytemp[1])
    yhi = ceiling(ytemp[2])
    box = c(xhi - xlo + 1, yhi - ylo + 1)
  }
  else {
  }
  cutout = magcutout(image, loc = c(xcen, ycen), box = box, 
                     shiftloc = shiftloc, paddim = paddim, plot = FALSE)
  cut_image = cutout$image
  xlo = cutout$loc.diff[1] + 1
  xhi = xlo + dim(cut_image)[1] - 1
  ylo = cutout$loc.diff[2] + 1
  yhi = ylo + dim(cut_image)[2] - 1
  xcen.new = xcen - xlo + 1
  ycen.new = ycen - ylo + 1
  pixscale = getpixscale(header = header, CD1_1 = CD1_1, CD1_2 = CD1_2, 
                         CD2_1 = CD2_1, CD2_2 = CD2_2)
  loc.diff = c(xlo - 1, ylo - 1)
  cut_xlo = 1
  cut_xhi = dim(cut_image)[1]
  cut_ylo = 1
  cut_yhi = dim(cut_image)[2]
  usr.WCS = rbind(magWCSxy2radec(xlo - 1, ylo - 1, header = header, 
                                 CRVAL1 = CRVAL1, CRVAL2 = CRVAL2, CRPIX1 = CRPIX1, CRPIX2 = CRPIX2, 
                                 CD1_1 = CD1_1, CD1_2 = CD1_2, CD2_1 = CD2_1, CD2_2 = CD2_2), 
                  magWCSxy2radec(xlo - 1, yhi, header = header, CRVAL1 = CRVAL1, 
                                 CRVAL2 = CRVAL2, CRPIX1 = CRPIX1, CRPIX2 = CRPIX2, 
                                 CD1_1 = CD1_1, CD1_2 = CD1_2, CD2_1 = CD2_1, CD2_2 = CD2_2), 
                  magWCSxy2radec(xhi, ylo - 1, header = header, CRVAL1 = CRVAL1, 
                                 CRVAL2 = CRVAL2, CRPIX1 = CRPIX1, CRPIX2 = CRPIX2, 
                                 CD1_1 = CD1_1, CD1_2 = CD1_2, CD2_1 = CD2_1, CD2_2 = CD2_2), 
                  magWCSxy2radec(xhi, yhi, header = header, CRVAL1 = CRVAL1, 
                                 CRVAL2 = CRVAL2, CRPIX1 = CRPIX1, CRPIX2 = CRPIX2, 
                                 CD1_1 = CD1_1, CD1_2 = CD1_2, CD2_1 = CD2_1, CD2_2 = CD2_2))
  usr.WCS = cbind(x.cut = c(cut_xlo - 1, cut_xlo - 1, cut_xhi, 
                            cut_xhi), y.cut = c(cut_ylo - 1, cut_yhi, cut_ylo - 1, 
                                                cut_yhi), x.orig = c(xlo - 1, xlo - 1, xhi, xhi), y.orig = c(ylo - 
                                                                                                               1, yhi, ylo - 1, yhi), usr.WCS)
  if (approx.map) {
    approx.map.RA = approxfun(seq(usr.WCS[1, "RA"], usr.WCS[4, 
                                                            "RA"], len = 100), seq(usr.WCS[1, "x.cut"], usr.WCS[4, 
                                                                                                                "x.cut"], len = 100))
    approx.map.Dec = approxfun(seq(usr.WCS[1, "Dec"], usr.WCS[4, 
                                                              "Dec"], len = 100), seq(usr.WCS[1, "y.cut"], usr.WCS[4, 
                                                                                                                   "y.cut"], len = 100))
    approx.map = function(RA, Dec) {
      if (length(dim(RA)) == 2) {
        Dec = RA[, 2]
        RA = RA[, 1]
      }
      invisible(cbind(x = approx.map.RA(RA), y = approx.map.Dec(Dec)))
    }
  }
  else {
    approx.map = NULL
  }
  if (!missing(header)) {
    dimdiff = dim(cut_image) - dim(image)
    hdradd = list(CRPIX1 = -loc.diff[1], CRPIX2 = -loc.diff[2], 
                  NAXIS1 = dimdiff[1], NAXIS2 = dimdiff[2])
    if (imtype == "FITSio") {
      for (hdrname in names(hdradd)) {
        if (hdradd[[hdrname]] != 0) {
          hdrrow = which(header == hdrname) + 1
          header[hdrrow] = as.character(as.numeric(header[hdrrow]) + 
                                          hdradd[[hdrname]])
        }
      }
    }
    else if (imtype == "astro") {
      for (hdrname in names(hdradd)) {
        if (hdradd[[hdrname]] != 0) {
          hdrrow = which(header[, "key"] == hdrname)
          header[hdrrow, "value"] = as.character(as.numeric(header[hdrrow, 
                                                                   "value"]) + hdradd[[hdrname]])
        }
      }
    }
    else {
      header = NULL
    }
  }
  else {
    header = NULL
  }
  output = list(image = cut_image, loc = c(x = as.numeric(xcen.new), 
                                           y = as.numeric(ycen.new)), loc.orig = c(x = as.numeric(xcen), 
                                                                                   y = as.numeric(ycen)), loc.diff = c(as.numeric(loc.diff[1]), 
                                                                                                                       as.numeric(loc.diff[2])), xsel = xlo:xhi, ysel = ylo:yhi, 
                loc.WCS = loc, scale.WCS = pixscale, usr.WCS = usr.WCS, 
                approx.map = approx.map, header = header)
  if (plot) {
    if (all(is.na(cut_image))) {
      cut_image[] = 0
      magimageWCS(image = cut_image, header = header, CRVAL1 = CRVAL1, 
                  CRVAL2 = CRVAL2, CRPIX1 = CRPIX1, CRPIX2 = CRPIX2, 
                  CD1_1 = CD1_1, CD1_2 = CD1_2, CD2_1 = CD2_1, 
                  CD2_2 = CD2_2, ...)
      cut_image[] = NA
    }
    else {
      magimageWCS(image = cut_image, header = header, CRVAL1 = CRVAL1, 
                  CRVAL2 = CRVAL2, CRPIX1 = CRPIX1, CRPIX2 = CRPIX2, 
                  CD1_1 = CD1_1, CD1_2 = CD1_2, CD2_1 = CD2_1, 
                  CD2_2 = CD2_2, ...)
    }
  }
  invisible(output)
}

magimage=function (x, y, z, zlim, xlim, ylim, col = grey((0:1000)/1000), 
                   add = FALSE, useRaster = TRUE, asp = 1, magmap = TRUE, locut = 0.4, 
                   hicut = 0.995, flip = FALSE, range = c(0, 1), type = "quan", 
                   stretch = "asinh", stretchscale = "auto", bad = NA, clip = "", 
                   axes = TRUE, frame.plot = TRUE, sparse = "auto", qdiff = FALSE, 
                   ...) 
{
  dots = list(...)
  dotskeepimage = c("xaxs", "yaxs", "breaks", "oldstyle")
  if (length(dots) > 0) {
    dotsimage = dots[names(dots) %in% dotskeepimage]
  }
  else {
    dotsimage = {
    }
  }
  if (!missing(x)) {
    if (is.list(x)) {
      if ("y" %in% names(x)) {
        y = x$y
      }
      if ("z" %in% names(x)) {
        z = x$z
      }
      if ("x" %in% names(x)) {
        x = x$x
      }
    }
    if (is.matrix(x)) {
      z = x
      x = seq(0.5, dim(z)[1] - 0.5)
      if (missing(y)) {
        y = seq(0.5, dim(z)[2] - 0.5)
      }
    }
  }
  if (!missing(z)) {
    if (is.matrix(z)) {
      if (missing(x)) {
        x = seq(0.5, dim(z)[1] - 0.5)
      }
      if (missing(y)) {
        y = seq(0.5, dim(z)[2] - 0.5)
      }
    }
  }
  if (is.vector(z) & !missing(x)) {
    z = matrix(z, nrow = length(x))
    if (missing(y)) {
      y = seq(0.5, dim(z)[2] - 0.5)
    }
  }
  if (is.vector(z) & !missing(y)) {
    z = matrix(z, ncol = length(y))
    if (missing(x)) {
      x = seq(0.5, dim(z)[1] - 0.5)
    }
  }
  if (missing(xlim) & length(x) == dim(z)[1]) {
    xlim = range(x, na.rm = TRUE) + c(-0.5, 0.5) * diff(range(x, 
                                                              na.rm = TRUE))/(dim(z)[1] - 1)
  }
  if (missing(ylim) & length(y) == dim(z)[2]) {
    ylim = range(y, na.rm = TRUE) + c(-0.5, 0.5) * diff(range(y, 
                                                              na.rm = TRUE))/(dim(z)[2] - 1)
  }
  if (missing(xlim) & length(x) == (dim(z)[1] + 1)) {
    xlim = range(x, na.rm = TRUE)
  }
  if (missing(ylim) & length(y) == (dim(z)[2] + 1)) {
    ylim = range(y, na.rm = TRUE)
  }
  if (x[1] > x[length(x)]) {
    x = rev(x)
    xlim = rev(xlim)
  }
  if (y[1] > y[length(y)]) {
    y = rev(y)
    ylim = rev(ylim)
  }
  if (sparse == "auto") {
    sparse = ceiling(max(dim(z)/1000))
  }
  if (sparse > 1) {
    samplex = seq(sparse/2, length(x), by = sparse)
    sampley = seq(sparse/2, length(y), by = sparse)
    x = x[samplex]
    y = y[sampley]
    z = z[samplex, sampley]
  }
  if (qdiff) {
    col = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100))
    if (missing(hicut)) {
      maximg = max(abs(z), na.rm = TRUE)
      locut = -maximg
      hicut = maximg
    }
    else {
      locut = -hicut
    }
    type = "num"
    zlim = c(0, 1)
  }
  if (magmap) {
    if (type == "quan") {
      if (quantile(z, locut, na.rm = T) != quantile(z, 
                                                    hicut, na.rm = T)) {
        z = magmap(data = z, locut = locut, hicut = hicut, 
                   flip = flip, range = range, type = type, stretch = stretch, 
                   stretchscale = stretchscale, bad = bad, clip = clip)$map
      }
      else {
        print("Too many same valued pixels: turning off magmap scaling!")
      }
    }
    else {
      z = magmap(data = z, locut = locut, hicut = hicut, 
                 flip = flip, range = range, type = type, stretch = stretch, 
                 stretchscale = stretchscale, bad = bad, clip = clip)$map
    }
  }
  if (missing(zlim)) {
    zlim = range(z, na.rm = TRUE)
  }
  do.call("image", c(list(x = x, y = y, z = z, zlim = zlim, 
                          xlim = xlim, ylim = ylim, col = col, add = add, useRaster = useRaster, 
                          axes = FALSE, asp = asp, xlab = "", ylab = "", main = ""), 
                     dotsimage))
  if (add == FALSE) {
    if (axes) {
      magaxis(...)
      if (frame.plot) {
        box()
      }
    }
  }
  return = list(x = x, y = y, z = z)
}

magmap=function (data, locut = 0, hicut = 1, flip = FALSE, range = c(0, 
                                                                     2/3), type = "quan", stretch = "lin", stretchscale = 1, bad = NA, 
                 clip = "") 
{
  if (stretchscale == "auto") {
    good = is.na(data) == FALSE & is.nan(data) == FALSE & 
      is.infinite(data) == FALSE & is.null(data) == FALSE
    if (length(which(good)) == 0) {
      stop("There is no numeric data!")
    }
    absdata = abs(data[good] - median(data[good], na.rm = TRUE))
    stretchscale = 1/median(absdata, na.rm = TRUE)
    if (!is.finite(stretchscale)) {
      stretchscale = 1
    }
  }
  if (stretch == "log" | stretch == "sqrt") {
    good = is.na(data) == FALSE & is.nan(data) == FALSE & 
      is.infinite(data) == FALSE & is.null(data) == FALSE & 
      data > 0
    if (length(which(good)) == 0) {
      stop("There is no numeric data with a value greater than 0!")
    }
  }
  else {
    good = is.na(data) == FALSE & is.nan(data) == FALSE & 
      is.infinite(data) == FALSE & is.null(data) == FALSE
    if (length(which(good)) == 0) {
      stop("There is no numeric data!")
    }
  }
  if (type == "quan") {
    locut = quantile(data[good], locut)
    hicut = quantile(data[good], hicut)
  }
  else if (type == "num") {
    locut = locut
    hicut = hicut
  }
  else if (type == "sig") {
    locut = quantile(data[good], pnorm(locut))
    hicut = quantile(data[good], pnorm(hicut))
  }
  else if (type == "rank") {
    locut = 1
    hicut = length(data[good])
    data[good][order(data[good])] = locut:hicut
  }
  else {
    stop(type, "is not a valid type option!")
  }
  loreturn = locut
  hireturn = hicut
  if (stretch == "log" & locut <= 0) {
    stop("locut <= 0 and stretch='log'- this is not allowed!")
  }
  if (stretch == "log" & hicut <= 0) {
    stop("hicut <=0 and stretch='log'- this is not allowed!")
  }
  if (locut > hicut) {
    stop("locut>hicut is not allowed")
  }
  if (locut == hicut) {
    data[good] = (range[2] + range[1])/2
  }
  if (locut < hicut) {
    if (stretch == "lin") {
    }
    else if (stretch == "log") {
      locut = log10(locut)
      hicut = log10(hicut)
      data = suppressWarnings(log10(data))
    }
    else if (stretch == "atan") {
      locut = atan(locut * stretchscale)
      hicut = atan(hicut * stretchscale)
      data = atan(data * stretchscale)
    }
    else if (stretch == "asinh") {
      locut = asinh(locut * stretchscale)
      hicut = asinh(hicut * stretchscale)
      data = asinh(data * stretchscale)
    }
    else if (stretch == "sqrt") {
      locut = sqrt(locut)
      hicut = sqrt(hicut)
      data = suppressWarnings(sqrt(data))
    }
    else if (stretch == "cdf") {
      cdf = ecdf(data[good])
      locut = cdf(locut)
      hicut = cdf(hicut)
      data[good] = cdf(data[good])
    }
    else {
      stop(paste(stretch, "is not a valid stretch option!"))
    }
    losel = data < locut & good
    hisel = data > hicut & good
    data[losel] = locut
    data[hisel] = hicut
    data[good] = data[good] - locut
    data[good] = range[1] + (data[good] * (range[2] - range[1])/(hicut - 
                                                                   locut))
    if (flip) {
      data[good] = range[2] - data[good] + range[1]
    }
    if (clip == "NA") {
      data[losel] = NA
      data[hisel] = NA
    }
  }
  data[!good] = bad
  return(list(map = data, datalim = c(loreturn, hireturn), 
              maplim = range, loclip = length(which(data[good] == range[1]))/length(data[good]), 
              hiclip = length(which(data[good] == range[2]))/length(data[good])))
}

magaxis=function (side = 1:2, majorn = 5, minorn = "auto", tcl = 0.5, 
                  ratio = 0.5, labels = TRUE, unlog = "auto", mgp = c(2, 0.5, 
                                                                      0), mtline = 2, xlab = NULL, ylab = NULL, crunch = TRUE, 
                  logpretty = TRUE, prettybase = 10, powbase = 10, hersh = FALSE, 
                  family = "sans", frame.plot = FALSE, usepar = FALSE, grid = FALSE, 
                  grid.col = "grey", grid.lty = 1, grid.lwd = 1, lwd.axis = 1, 
                  lwd.ticks = lwd.axis, ...) 
{
  dots = list(...)
  dotskeepaxis = c("cex.axis", "col.axis", "font.axis", "xaxp", 
                   "yaxp", "tck", "las", "fg", "xpd", "xaxt", "yaxt", "col.ticks", 
                   "tick")
  dotskeepmtext = c("cex.lab", "col.lab", "font.lab")
  if (length(dots) > 0) {
    dotsaxis = dots[names(dots) %in% dotskeepaxis]
    dotsmtext = dots[names(dots) %in% dotskeepmtext]
  }
  else {
    dotsaxis = {
    }
    dotsmtext = {
    }
  }
  if (length(mtline) == 1) {
    mtline = rep(mtline, 2)
  }
  majornlist = majorn
  minornlist = minorn
  labelslist = labels
  unloglist = unlog
  crunchlist = crunch
  logprettylist = logpretty
  prettybaselist = prettybase
  powbaselist = powbase
  gridlist = grid
  if (length(majorn) == 1 & length(side) > 1) {
    majornlist = rep(majorn, length(side))
  }
  if (length(minorn) == 1 & length(side) > 1) {
    minornlist = rep(minorn, length(side))
  }
  if (length(labels) == 1 & length(side) > 1) {
    labelslist = rep(labels, length(side))
  }
  if (length(unlog) == 1 & length(side) > 1 & (unlog[1] == 
                                               T | unlog[1] == F | unlog[1] == "auto")) {
    unloglist = rep(unlog, length(side))
  }
  if (length(crunch) == 1 & length(side) > 1) {
    crunchlist = rep(crunch, length(side))
  }
  if (length(logpretty) == 1 & length(side) > 1) {
    logprettylist = rep(logpretty, length(side))
  }
  if (length(prettybase) == 1 & length(side) > 1) {
    prettybaselist = rep(prettybase, length(side))
  }
  if (length(powbase) == 1 & length(side) > 1) {
    powbaselist = rep(powbase, length(side))
  }
  if (length(grid) == 1 & length(side) > 1) {
    gridlist = rep(grid, length(side))
  }
  if (unlog[1] == "") {
    unloglist = rep(FALSE, length(side))
  }
  if (unlog[1] == "x") {
    unloglist = rep(FALSE, length(side))
    unloglist[side %in% c(1, 3)] = TRUE
  }
  if (unlog[1] == "y") {
    unloglist = rep(FALSE, length(side))
    unloglist[side %in% c(2, 4)] = TRUE
  }
  if (unlog[1] == "xy" | unlog[1] == "yx") {
    unloglist = rep(TRUE, length(side))
  }
  if (length(majornlist) != length(side)) {
    stop("Length of majorn vector mismatches number of axes!")
  }
  if (length(minornlist) != length(side)) {
    stop("Length of minorn vector mismatches number of axes!")
  }
  if (length(labelslist) != length(side)) {
    stop("Length of labels vector mismatches number of axes!")
  }
  if (length(unloglist) != length(side)) {
    stop("Length of unlog vector mismatches number of axes!")
  }
  if (length(crunchlist) != length(side)) {
    stop("Length of crunch vector mismatches number of axes!")
  }
  if (length(logprettylist) != length(side)) {
    stop("Length of logpretty vector mismatches number of axes!")
  }
  if (length(prettybaselist) != length(side)) {
    stop("Length of prettybase vector mismatches number of axes!")
  }
  if (length(powbaselist) != length(side)) {
    stop("Length of powbase vector mismatches number of axes!")
  }
  if (length(gridlist) != length(side)) {
    stop("Length of grid vector mismatches number of axes!")
  }
  currentfamily = par("family")
  if (hersh & family == "serif") {
    par(family = "HersheySerif")
  }
  if (hersh & family == "sans") {
    par(family = "HersheySans")
  }
  if (hersh == F & family == "serif") {
    par(family = "serif")
  }
  if (hersh == F & family == "sans") {
    par(family = "sans")
  }
  if (missing(lwd.axis)) {
    lwd.axis = par()$lwd
  }
  if (missing(lwd.ticks)) {
    lwd.ticks = par()$lwd
  }
  if (usepar) {
    if (missing(tcl)) {
      tcl = par()$tcl
    }
    if (missing(mgp)) {
      mgp = par()$mgp
    }
  }
  for (i in 1:length(side)) {
    currentside = side[i]
    majorn = majornlist[i]
    minorn = minornlist[i]
    labels = labelslist[i]
    unlog = unloglist[i]
    crunch = crunchlist[i]
    logpretty = logprettylist[i]
    prettybase = prettybaselist[i]
    powbase = powbaselist[i]
    grid = gridlist[i]
    lims = par("usr")
    if (currentside %in% c(1, 3)) {
      lims = lims[1:2]
      if (par("xlog")) {
        logged = T
      }
      else {
        logged = F
      }
    }
    else {
      lims = lims[3:4]
      if (par("ylog")) {
        logged = T
      }
      else {
        logged = F
      }
    }
    lims = sort(lims)
    if (unlog == "auto") {
      if (logged) {
        unlog = T
      }
      else {
        unlog = F
      }
    }
    if ((logged | unlog) & powbase == 10) {
      usemultloc = (10^lims[2])/(10^lims[1]) < 50
    }
    else {
      usemultloc = F
    }
    if (unlog) {
      sci.tick = maglab(10^lims, n = majorn, log = T, exptext = T, 
                        crunch = crunch, logpretty = logpretty, usemultloc = usemultloc, 
                        prettybase = prettybase, powbase = powbase, hersh = hersh)
      major.ticks = log(sci.tick$tickat, powbase)
      uselabels = sci.tick$exp
      labloc = log(sci.tick$labat, powbase)
      if (usemultloc == F) {
        if (minorn == "auto") {
          splitmin = (powbase^major.ticks[2])/(powbase^major.ticks[1])
        }
        else {
          splitmin = minorn + 1
        }
        if (splitmin > 10) {
          minors = seq(major.ticks[1], major.ticks[2]) - 
            major.ticks[1]
        }
        else {
          minors = log(seq(powbase^major.ticks[1], powbase^major.ticks[2], 
                           len = splitmin), powbase) - major.ticks[1]
        }
      }
    }
    if (logged & unlog == F) {
      sci.tick = maglab(10^lims, n = majorn, log = T, exptext = F, 
                        crunch = crunch, logpretty = logpretty, usemultloc = usemultloc, 
                        prettybase = prettybase, powbase = powbase, hersh = hersh)
      major.ticks = log(sci.tick$tickat, powbase)
      uselabels = sci.tick$exp
      labloc = log(sci.tick$labat, powbase)
      if (usemultloc == F) {
        if (minorn == "auto") {
          splitmin = (powbase^major.ticks[2])/(powbase^major.ticks[1])
        }
        else {
          splitmin = minorn + 1
        }
        if (splitmin > 10) {
          minors = seq(major.ticks[1], major.ticks[2]) - 
            major.ticks[1]
        }
        else {
          minors = log(seq(powbase^major.ticks[1], powbase^major.ticks[2], 
                           len = splitmin), powbase) - major.ticks[1]
        }
      }
    }
    if (logged == F & unlog == F) {
      sci.tick = maglab(lims, n = majorn, log = F, exptext = F, 
                        prettybase = prettybase, hersh = hersh)
      major.ticks = sci.tick$tickat
      uselabels = sci.tick$exp
      labloc = sci.tick$labat
      if (minorn == "auto") {
        splitmin = length(pretty(major.ticks[1:2]))
      }
      else {
        splitmin = minorn + 1
      }
      minors = seq(major.ticks[1], major.ticks[2], len = splitmin) - 
        major.ticks[1]
    }
    if (grid) {
      if (currentside == 1) {
        if (logged) {
          abline(v = powbase^labloc, col = grid.col, 
                 lty = grid.lty, lwd = grid.lty)
        }
        else {
          abline(v = labloc, col = grid.col, lty = grid.lty, 
                 lwd = grid.lty)
        }
      }
      if (currentside == 2) {
        if (logged) {
          abline(h = powbase^labloc, col = grid.col, 
                 lty = grid.lty, lwd = grid.lty)
        }
        else {
          abline(h = labloc, col = grid.col, lty = grid.lty, 
                 lwd = grid.lty)
        }
      }
    }
    if (logged) {
      do.call("axis", c(list(side = currentside, at = powbase^major.ticks, 
                             tcl = tcl, labels = FALSE, mgp = mgp, lwd = lwd.axis, 
                             lwd.ticks = lwd.ticks), dotsaxis))
    }
    else {
      do.call("axis", c(list(side = currentside, at = major.ticks, 
                             tcl = tcl, labels = FALSE, mgp = mgp, lwd = lwd.axis, 
                             lwd.ticks = lwd.ticks), dotsaxis))
    }
    if (labels) {
      if (logged) {
        do.call("axis", c(list(side = currentside, at = powbase^labloc, 
                               tick = F, labels = uselabels, mgp = mgp, lwd = lwd.axis, 
                               lwd.ticks = lwd.ticks), dotsaxis))
      }
      else {
        do.call("axis", c(list(side = currentside, at = labloc, 
                               tick = F, labels = uselabels, mgp = mgp, lwd = lwd.axis, 
                               lwd.ticks = lwd.ticks), dotsaxis))
      }
    }
    if (usemultloc == F & minorn > 1) {
      minors = minors[-c(1, length(minors))]
      minor.ticks = c(outer(minors, major.ticks, `+`))
      if (logged) {
        do.call("axis", c(list(side = currentside, at = powbase^minor.ticks, 
                               tcl = tcl * ratio, labels = FALSE, mgp = mgp, 
                               lwd = lwd.axis, lwd.ticks = lwd.ticks), dotsaxis))
      }
      else {
        do.call("axis", c(list(side = currentside, at = minor.ticks, 
                               tcl = tcl * ratio, labels = FALSE, mgp = mgp, 
                               lwd = lwd.axis, lwd.ticks = lwd.ticks), dotsaxis))
      }
    }
  }
  if (length(dotsmtext) > 0) {
    names(dotsmtext) = c("cex", "col", "font")[match(names(dotsmtext), 
                                                     dotskeepmtext)]
  }
  if (is.null(xlab) == FALSE) {
    do.call("mtext", c(list(text = xlab, side = ifelse(side[1] %in% 
                                                         c(1, 3), side[1], side[2]), line = mtline[1]), dotsmtext))
  }
  if (is.null(ylab) == FALSE) {
    do.call("mtext", c(list(text = ylab, side = ifelse(side[2] %in% 
                                                         c(2, 4), side[2], side[1]), line = mtline[2]), dotsmtext))
  }
  if (frame.plot) {
    box()
  }
  par(family = currentfamily)
}

maglab=function (lims, n, log = FALSE, exptext = TRUE, crunch = TRUE, 
                 logpretty = TRUE, usemultloc = FALSE, multloc = c(1, 2, 5), 
                 prettybase = 10, powbase = 10, hersh = FALSE, trim = FALSE) 
{
  if (usemultloc & log == F) {
    stop("If using multloc then log must be TRUE!")
  }
  lims = lims/(prettybase/10)
  if (log & usemultloc == F) {
    lims = log(lims, powbase)
  }
  if (usemultloc == F) {
    if (missing(n)) {
      labloc = pretty(lims)
    }
    else {
      labloc = pretty(lims, n)
    }
  }
  if (log) {
    if (usemultloc == F) {
      labloc = labloc + log10(prettybase/10)
      labloc = labloc[round(labloc - log(prettybase/10, 
                                         powbase), 10)%%1 == 0]
      if (min(labloc) > lims[1]) {
        labloc = c(min(labloc) - 1, labloc)
      }
      if (max(labloc) < lims[2]) {
        labloc = c(labloc, max(labloc) + 1)
      }
      labloc = round(labloc, 10)
      labloc = powbase^labloc
      tickloc = labloc
    }
    if (usemultloc) {
      labloc = {
      }
      for (i in 1:length(multloc)) {
        labloc = c(labloc, multloc[i] * powbase^seq(ceiling(log(lims[1], 
                                                                powbase)) - 1, floor(log(lims[2], powbase)) + 
                                                      1))
      }
      labloc = sort(labloc)
      tickloc = {
      }
      for (i in 1:9) {
        tickloc = c(tickloc, i * powbase^seq(ceiling(log(lims[1], 
                                                         powbase)) - 1, floor(log(lims[2], powbase)) + 
                                               1))
      }
      tickloc = sort(tickloc)
    }
    char = {
    }
    if (exptext) {
      for (i in 1:length(labloc)) {
        char = c(char, format(labloc[i]))
      }
    }
    if (!exptext) {
      for (i in 1:length(labloc)) {
        char = c(char, format(log(labloc[i], powbase)))
      }
    }
  }
  else {
    labloc = labloc * (prettybase/10)
    tickloc = labloc
    char = {
    }
    for (i in 1:length(labloc)) {
      char = c(char, format(labloc[i]))
    }
  }
  if (log & usemultloc == F) {
    lims = powbase^(lims)
  }
  if (trim) {
    char = char[labloc >= lims[1] & labloc <= lims[2]]
    labloc = labloc[labloc >= lims[1] & labloc <= lims[2]]
    tickloc = tickloc[tickloc >= lims[1] & tickloc <= lims[2]]
  }
  check = grep("e", char)
  if (length(check) > 0) {
    char = format(labloc, scientific = T)
    check = grep("0e+00", char, fixed = T)
    char[check] = "0"
    if (hersh) {
      check = grep("e+0", char, fixed = T)
      char[check] = sub("e+0", "e+", char[check], fixed = T)
      check = grep("e-0", char, fixed = T)
      char[check] = sub("e-0", "e-", char[check], fixed = T)
      check = grep("e+", char, fixed = T)
      char[check] = paste(sub("e+", "\\mu10\\sp", char[check], 
                              fixed = T), "\\ep", sep = "")
      check = grep("e-", char, fixed = T)
      char[check] = paste(sub("e-", "\\mu10\\sp-", char[check], 
                              fixed = T), "\\ep", sep = "")
    }
    else {
      check = grep("e+", char, fixed = T)
      char[check] = paste(sub("e+", "*x*10^{", char[check], 
                              fixed = T), "}", sep = "")
      check = grep("e-", char, fixed = T)
      char[check] = paste(sub("e-", "*x*10^{-", char[check], 
                              fixed = T), "}", sep = "")
    }
  }
  if (crunch) {
    check = grep("1*x*", char)
    if (length(check) > 0) {
      if (hersh) {
        char[check] = sub("1\\mu", "", char[check], fixed = T)
      }
      else {
        char[check] = sub("1*x*", "", char[check], fixed = T)
      }
    }
  }
  if (hersh) {
    exp = char
  }
  else {
    exp = parse(text = char)
  }
  return(list(tickat = tickloc, labat = labloc, exp = exp))
}

magimageRGB=function (x, y, R, G, B, saturation = 1, zlim, xlim, ylim, add = FALSE, 
                      useRaster = TRUE, asp = 1, magmap = TRUE, locut = 0.4, hicut = 0.995, 
                      flip = FALSE, range = c(0, 1), type = "quan", stretch = "asinh", 
                      stretchscale = "auto", bad = range[1], clip = "", axes = TRUE, 
                      frame.plot = TRUE, sparse = "auto", ...) 
{
  dots = list(...)
  dotskeepimage = c("xaxs", "yaxs", "breaks", "oldstyle")
  if (length(dots) > 0) {
    dotsimage = dots[names(dots) %in% dotskeepimage]
  }
  else {
    dotsimage = {
    }
  }
  if (!missing(x)) {
    if (is.list(x)) {
      if ("R" %in% names(x)) {
        R = x$R
      }
      if ("G" %in% names(x)) {
        G = x$G
      }
      if ("B" %in% names(x)) {
        B = x$B
      }
      if ("y" %in% names(x)) {
        y = x$y
      }
      else {
        y = seq(0.5, dim(R)[2] - 0.5)
      }
      if ("x" %in% names(x)) {
        x = x$x
      }
      else {
        x = seq(0.5, dim(R)[1] - 0.5)
      }
    }
    if (is.array(x)) {
      R = x[, , 1]
      G = x[, , 2]
      B = x[, , 3]
      x = seq(0.5, dim(R)[1] - 0.5)
      if (missing(y)) {
        y = seq(0.5, dim(R)[2] - 0.5)
      }
    }
  }
  if (!missing(R)) {
    if (is.matrix(R)) {
      if (missing(x)) {
        x = seq(0.5, dim(R)[1] - 0.5)
      }
      if (missing(y)) {
        y = seq(0.5, dim(R)[2] - 0.5)
      }
    }
  }
  if (missing(xlim) & length(x) == dim(R)[1]) {
    xlim = range(x, na.rm = TRUE) + c(-0.5, 0.5) * diff(range(x, 
                                                              na.rm = TRUE))/(dim(R)[1] - 1)
  }
  if (missing(ylim) & length(y) == dim(R)[2]) {
    ylim = range(y, na.rm = TRUE) + c(-0.5, 0.5) * diff(range(y, 
                                                              na.rm = TRUE))/(dim(R)[2] - 1)
  }
  if (missing(xlim) & length(x) == (dim(R)[1] + 1)) {
    xlim = range(x, na.rm = TRUE)
  }
  if (missing(ylim) & length(y) == (dim(R)[2] + 1)) {
    ylim = range(y, na.rm = TRUE)
  }
  if (is.numeric(x)) {
    if (x[1] > x[length(x)]) {
      x = rev(x)
      xlim = rev(xlim)
    }
  }
  if (is.numeric(y)) {
    if (y[1] > y[length(y)]) {
      y = rev(y)
      ylim = rev(ylim)
    }
  }
  if (sparse == "auto") {
    sparse = ceiling(max(dim(R)/1000))
  }
  if (sparse > 1) {
    samplex = seq(sparse/2, length(x), by = sparse)
    sampley = seq(sparse/2, length(y), by = sparse)
    x = x[samplex]
    y = y[sampley]
    R = R[samplex, sampley]
    G = G[samplex, sampley]
    B = B[samplex, sampley]
  }
  if (missing(zlim)) {
    zlim = c(0, length(R))
  }
  if (length(locut) < 3) {
    locut = rep(locut[1], 3)
  }
  if (length(hicut) < 3) {
    hicut = rep(hicut[1], 3)
  }
  if (length(stretchscale) < 3) {
    stretchscale = rep(stretchscale[1], 3)
  }
  if (magmap) {
    R = magmap(data = R, locut = locut[1], hicut = hicut[1], 
               flip = flip, range = range, type = type, stretch = stretch, 
               stretchscale = stretchscale[1], bad = bad, clip = clip)$map
    G = magmap(data = G, locut = locut[2], hicut = hicut[2], 
               flip = flip, range = range, type = type, stretch = stretch, 
               stretchscale = stretchscale[2], bad = bad, clip = clip)$map
    B = magmap(data = B, locut = locut[3], hicut = hicut[3], 
               flip = flip, range = range, type = type, stretch = stretch, 
               stretchscale = stretchscale[3], bad = bad, clip = clip)$map
  }
  if (saturation != 1) {
    Y = (R + G + B)/3
    Ra = Y + saturation * (2 * R - G - B)/3
    Ga = Y + saturation * (2 * G - R - B)/3
    Ba = Y + saturation * (2 * B - R - G)/3
    rm(Y)
    R = Ra
    rm(Ra)
    G = Ga
    rm(Ga)
    B = Ba
    rm(Ba)
  }
  do.call("image", c(list(x = x, y = y, z = matrix(1:length(R), 
                                                   dim(R)[1]), zlim = zlim, xlim = xlim, ylim = ylim, col = rgb(R, 
                                                                                                                G, B), add = add, useRaster = useRaster, axes = FALSE, 
                          asp = asp, xlab = "", ylab = "", main = ""), dotsimage))
  if (add == FALSE) {
    if (axes) {
      magaxis(...)
      if (frame.plot) {
        box()
      }
    }
  }
  return = list(x = x, y = y, R = R, G = G, B = B)
}

getpixscale=function (header, CD1_1 = 1, CD1_2 = 0, CD2_1 = 0, CD2_2 = 1) 
{
  if (!missing(header)) {
    if (is.data.frame(header) | is.matrix(header)) {
      locs = match(c("CD1_1", "CD1_2", "CD2_1", "CD2_2", 
                     "CDELT1", "CDELT2"), header[, 1])
      locs = locs[is.na(locs) == FALSE]
      headerWCS = data.frame(header[locs, 1], as.numeric(header[locs, 
                                                                2]), stringsAsFactors = FALSE)
      if ("CD1_1" %in% headerWCS[, 1]) {
        CD1_1 = headerWCS[headerWCS[, 1] == "CD1_1", 
                          2]
        if ("CD1_2" %in% headerWCS[, 1]) {
          CD1_2 = headerWCS[headerWCS[, 1] == "CD1_2", 
                            2]
        }
        else {
          message("Missing CD1_2")
        }
      }
      else {
        if ("CDELT1" %in% headerWCS[, 1]) {
          CD1_1 = headerWCS[headerWCS[, 1] == "CDELT1", 
                            2]
        }
        else {
          message("Missing CD1_1 and CDELT1")
        }
      }
      if ("CD2_2" %in% headerWCS[, 1]) {
        CD2_2 = headerWCS[headerWCS[, 1] == "CD2_2", 
                          2]
        if ("CD2_1" %in% headerWCS[, 1]) {
          CD2_1 = headerWCS[headerWCS[, 1] == "CD2_1", 
                            2]
        }
        else {
          message("Missing CD2_1")
        }
      }
      else {
        if ("CDELT2" %in% headerWCS[, 1]) {
          CD2_2 = headerWCS[headerWCS[, 1] == "CDELT2", 
                            2]
        }
        else {
          message("Missing CD2_2 and CDELT2")
        }
      }
    }
    else {
      if ("CD1_1" %in% header) {
        CD1_1 = as.numeric(header[which(header == "CD1_1") + 
                                    1])
        if ("CD1_2" %in% header) {
          CD1_2 = as.numeric(header[which(header == "CD1_2") + 
                                      1])
        }
        else {
          message("Missing CD1_2")
        }
      }
      else {
        if ("CDELT1" %in% header) {
          CD1_1 = as.numeric(header[which(header == "CDELT1") + 
                                      1])
        }
        else {
          message("Missing CD1_1 and CDELT1")
        }
      }
      if ("CD2_2" %in% header) {
        CD2_2 = as.numeric(header[which(header == "CD2_2") + 
                                    1])
        if ("CD2_1" %in% header) {
          CD2_1 = as.numeric(header[which(header == "CD2_1") + 
                                      1])
        }
        else {
          message("Missing CD2_1")
        }
      }
      else {
        if ("CDELT1" %in% header) {
          CD2_2 = as.numeric(header[which(header == "CDELT2") + 
                                      1])
        }
        else {
          message("Missing CD2_2 and CDELT2")
        }
      }
    }
  }
  return(3600 * (sqrt(CD1_1^2 + CD1_2^2) + sqrt(CD2_1^2 + CD2_2^2))/2)
}
radec2xy=function (RA, Dec, header, CRVAL1 = 0, CRVAL2 = 0, CRPIX1 = 0, 
                   CRPIX2 = 0, CD1_1 = 1, CD1_2 = 0, CD2_1 = 0, CD2_2 = 1, CTYPE1 = "RA--TAN", 
                   CTYPE2 = "DEC--TAN") 
{
  if (length(dim(RA)) == 2) {
    Dec = RA[, 2]
    RA = RA[, 1]
  }
  RA = as.numeric(RA)
  Dec = as.numeric(Dec)
  if (!missing(header)) {
    if (is.data.frame(header) | is.matrix(header)) {
      locs = match(c("CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2", 
                     "CRPIX1", "CRPIX2", "CD1_1", "CD1_2", "CD2_1", 
                     "CD2_2", "CDELT1", "CDELT2"), header[, 1])
      locs = locs[is.na(locs) == FALSE]
      headerWCS = data.frame(header[locs, 1], as.character(header[locs, 
                                                                  2]), stringsAsFactors = FALSE)
      if ("CTYPE1" %in% headerWCS[, 1]) {
        CTYPE1 = headerWCS[headerWCS[, 1] == "CTYPE1", 
                           2]
      }
      else {
        message("Missing CTYPE1")
      }
      if ("CTYPE2" %in% headerWCS[, 1]) {
        CTYPE2 = headerWCS[headerWCS[, 1] == "CTYPE2", 
                           2]
      }
      else {
        message("Missing CTYPE2")
      }
      if ("CRVAL1" %in% headerWCS[, 1]) {
        CRVAL1 = as.numeric(headerWCS[headerWCS[, 1] == 
                                        "CRVAL1", 2])
      }
      else {
        message("Missing CRVAL1")
      }
      if ("CRVAL2" %in% headerWCS[, 1]) {
        CRVAL2 = as.numeric(headerWCS[headerWCS[, 1] == 
                                        "CRVAL2", 2])
      }
      else {
        message("Missing CRVAL2")
      }
      if ("CRPIX1" %in% headerWCS[, 1]) {
        CRPIX1 = as.numeric(headerWCS[headerWCS[, 1] == 
                                        "CRPIX1", 2])
      }
      else {
        message("Missing CRPIX1")
      }
      if ("CRPIX2" %in% headerWCS[, 1]) {
        CRPIX2 = as.numeric(headerWCS[headerWCS[, 1] == 
                                        "CRPIX2", 2])
      }
      else {
        message("Missing CRPIX2")
      }
      if ("CD1_1" %in% headerWCS[, 1]) {
        CD1_1 = as.numeric(headerWCS[headerWCS[, 1] == 
                                       "CD1_1", 2])
        if ("CD1_2" %in% headerWCS[, 1]) {
          CD1_2 = as.numeric(headerWCS[headerWCS[, 1] == 
                                         "CD1_2", 2])
        }
        else {
          message("Missing CD1_2")
        }
      }
      else {
        if ("CDELT1" %in% headerWCS[, 1]) {
          CD1_1 = as.numeric(headerWCS[headerWCS[, 1] == 
                                         "CDELT1", 2])
        }
        else {
          message("Missing CD1_1 and CDELT1")
        }
      }
      if ("CD2_2" %in% headerWCS[, 1]) {
        CD2_2 = as.numeric(headerWCS[headerWCS[, 1] == 
                                       "CD2_2", 2])
        if ("CD2_1" %in% headerWCS[, 1]) {
          CD2_1 = as.numeric(headerWCS[headerWCS[, 1] == 
                                         "CD2_1", 2])
        }
        else {
          message("Missing CD2_1")
        }
      }
      else {
        if ("CDELT2" %in% headerWCS[, 1]) {
          CD2_2 = as.numeric(headerWCS[headerWCS[, 1] == 
                                         "CDELT2", 2])
        }
        else {
          message("Missing CD2_2 and CDELT2")
        }
      }
    }
    else {
      if ("CTYPE1" %in% header) {
        CTYPE1 = as.character(header[which(header == 
                                             "CTYPE1") + 1])
      }
      else {
        message("Missing CTYPE1")
      }
      if ("CTYPE2" %in% header) {
        CTYPE2 = as.character(header[which(header == 
                                             "CTYPE2") + 1])
      }
      else {
        message("Missing CTYPE2")
      }
      if ("CRVAL1" %in% header) {
        CRVAL1 = as.numeric(header[which(header == "CRVAL1") + 
                                     1])
      }
      else {
        message("Missing CRVAL1")
      }
      if ("CRVAL2" %in% header) {
        CRVAL2 = as.numeric(header[which(header == "CRVAL2") + 
                                     1])
      }
      else {
        message("Missing CRVAL2")
      }
      if ("CRPIX1" %in% header) {
        CRPIX1 = as.numeric(header[which(header == "CRPIX1") + 
                                     1])
      }
      else {
        message("Missing CRPIX1")
      }
      if ("CRPIX2" %in% header) {
        CRPIX2 = as.numeric(header[which(header == "CRPIX2") + 
                                     1])
      }
      else {
        message("Missing CRPIX2")
      }
      if ("CD1_1" %in% header) {
        CD1_1 = as.numeric(header[which(header == "CD1_1") + 
                                    1])
        if ("CD1_2" %in% header) {
          CD1_2 = as.numeric(header[which(header == "CD1_2") + 
                                      1])
        }
        else {
          message("Missing CD1_2")
        }
      }
      else {
        if ("CDELT1" %in% header) {
          CD1_1 = as.numeric(header[which(header == "CDELT1") + 
                                      1])
        }
        else {
          message("Missing CD1_1 and CDELT1")
        }
      }
      if ("CD2_2" %in% header) {
        CD2_2 = as.numeric(header[which(header == "CD2_2") + 
                                    1])
        if ("CD2_1" %in% header) {
          CD2_1 = as.numeric(header[which(header == "CD2_1") + 
                                      1])
        }
        else {
          message("Missing CD2_1")
        }
      }
      else {
        if ("CDELT1" %in% header) {
          CD2_2 = as.numeric(header[which(header == "CDELT2") + 
                                      1])
        }
        else {
          message("Missing CD2_2 and CDELT2")
        }
      }
    }
  }
  RA0 = CRVAL1
  Dec0 = CRVAL2
  x0 = CRPIX1
  y0 = CRPIX2
  x1 = CD1_1
  x2 = CD1_2
  y1 = CD2_1
  y2 = CD2_2
  RA0 = RA0 * (pi/180)
  Dec0 = Dec0 * (pi/180)
  RA = RA * (pi/180)
  Dec = Dec * (pi/180)
  scalemat = matrix(c(x1, x2, y1, y2), 2) * (pi/180)
  if (grepl("TAN", CTYPE1)) {
    cosc1 = sin(Dec0) * sin(Dec) + (cos(Dec0) * cos(Dec) * 
                                      cos(RA - RA0))
  }
  else if (grepl("SIN", CTYPE1) | grepl("NCP", CTYPE1)) {
    if (grepl("NCP", CTYPE1)) {
      message("Approximating deprecated CTYPE1 NCP with SIN!")
    }
    cosc1 = 1
  }
  else {
    stop("Projection system is not recognised. Must be either TAN, SIN or NCP!")
  }
  if (grepl("TAN", CTYPE2)) {
    cosc2 = sin(Dec0) * sin(Dec) + (cos(Dec0) * cos(Dec) * 
                                      cos(RA - RA0))
  }
  else if (grepl("SIN", CTYPE2) | grepl("NCP", CTYPE2)) {
    if (grepl("NCP", CTYPE2)) {
      message("Approximating deprecated CTYPE2 NCP with SIN!")
    }
    cosc2 = 1
  }
  else {
    stop("Projection system is not recognised. Must be either TAN, SIN or NCP!")
  }
  xxfunc = function(RA0, Dec0, RA, Dec) {
    (cos(Dec) * sin(RA - RA0))/cosc1
  }
  yyfunc = function(RA0, Dec0, RA, Dec) {
    ((cos(Dec0) * sin(Dec)) - (sin(Dec0) * cos(Dec) * cos(RA - 
                                                            RA0)))/cosc2
  }
  XX = xxfunc(RA0, Dec0, RA, Dec)
  YY = yyfunc(RA0, Dec0, RA, Dec)
  raw = cbind(XX, YY)
  output = raw %*% solve(scalemat)
  output[, 1] = output[, 1] + x0
  output[, 2] = output[, 2] + y0
  colnames(output) = c("x", "y")
  return(output)
}

xy2radec=function (x, y, header, CRVAL1 = 0, CRVAL2 = 0, CRPIX1 = 0, CRPIX2 = 0, 
                   CD1_1 = 1, CD1_2 = 0, CD2_1 = 0, CD2_2 = 1, CTYPE1 = "RA--TAN", 
                   CTYPE2 = "DEC--TAN") 
{
  if (length(dim(x)) == 2) {
    y = x[, 2]
    x = x[, 1]
  }
  x = as.numeric(x)
  y = as.numeric(y)
  if (!missing(header)) {
    if (is.data.frame(header) | is.matrix(header)) {
      locs = match(c("CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2", 
                     "CRPIX1", "CRPIX2", "CD1_1", "CD1_2", "CD2_1", 
                     "CD2_2", "CDELT1", "CDELT2"), header[, 1])
      locs = locs[is.na(locs) == FALSE]
      headerWCS = data.frame(header[locs, 1], as.character(header[locs, 
                                                                  2]), stringsAsFactors = FALSE)
      if ("CTYPE1" %in% headerWCS[, 1]) {
        CTYPE1 = headerWCS[headerWCS[, 1] == "CTYPE1", 
                           2]
      }
      else {
        message("Missing CTYPE1")
      }
      if ("CTYPE2" %in% headerWCS[, 1]) {
        CTYPE2 = headerWCS[headerWCS[, 1] == "CTYPE2", 
                           2]
      }
      else {
        message("Missing CTYPE2")
      }
      if ("CRVAL1" %in% headerWCS[, 1]) {
        CRVAL1 = as.numeric(headerWCS[headerWCS[, 1] == 
                                        "CRVAL1", 2])
      }
      else {
        message("Missing CRVAL1")
      }
      if ("CRVAL2" %in% headerWCS[, 1]) {
        CRVAL2 = as.numeric(headerWCS[headerWCS[, 1] == 
                                        "CRVAL2", 2])
      }
      else {
        message("Missing CRVAL2")
      }
      if ("CRPIX1" %in% headerWCS[, 1]) {
        CRPIX1 = as.numeric(headerWCS[headerWCS[, 1] == 
                                        "CRPIX1", 2])
      }
      else {
        message("Missing CRPIX1")
      }
      if ("CRPIX2" %in% headerWCS[, 1]) {
        CRPIX2 = as.numeric(headerWCS[headerWCS[, 1] == 
                                        "CRPIX2", 2])
      }
      else {
        message("Missing CRPIX2")
      }
      if ("CD1_1" %in% headerWCS[, 1]) {
        CD1_1 = as.numeric(headerWCS[headerWCS[, 1] == 
                                       "CD1_1", 2])
        if ("CD1_2" %in% headerWCS[, 1]) {
          CD1_2 = as.numeric(headerWCS[headerWCS[, 1] == 
                                         "CD1_2", 2])
        }
        else {
          message("Missing CD1_2")
        }
      }
      else {
        if ("CDELT1" %in% headerWCS[, 1]) {
          CD1_1 = as.numeric(headerWCS[headerWCS[, 1] == 
                                         "CDELT1", 2])
        }
        else {
          message("Missing CD1_1 and CDELT1")
        }
      }
      if ("CD2_2" %in% headerWCS[, 1]) {
        CD2_2 = as.numeric(headerWCS[headerWCS[, 1] == 
                                       "CD2_2", 2])
        if ("CD2_1" %in% headerWCS[, 1]) {
          CD2_1 = as.numeric(headerWCS[headerWCS[, 1] == 
                                         "CD2_1", 2])
        }
        else {
          message("Missing CD2_1")
        }
      }
      else {
        if ("CDELT2" %in% headerWCS[, 1]) {
          CD2_2 = as.numeric(headerWCS[headerWCS[, 1] == 
                                         "CDELT2", 2])
        }
        else {
          message("Missing CD2_2 and CDELT2")
        }
      }
    }
    else {
      if ("CTYPE1" %in% header) {
        CTYPE1 = as.character(header[which(header == 
                                             "CTYPE1") + 1])
      }
      else {
        message("Missing CTYPE1")
      }
      if ("CTYPE2" %in% header) {
        CTYPE2 = as.character(header[which(header == 
                                             "CTYPE2") + 1])
      }
      else {
        message("Missing CTYPE2")
      }
      if ("CRVAL1" %in% header) {
        CRVAL1 = as.numeric(header[which(header == "CRVAL1") + 
                                     1])
      }
      else {
        message("Missing CRVAL1")
      }
      if ("CRVAL2" %in% header) {
        CRVAL2 = as.numeric(header[which(header == "CRVAL2") + 
                                     1])
      }
      else {
        message("Missing CRVAL2")
      }
      if ("CRPIX1" %in% header) {
        CRPIX1 = as.numeric(header[which(header == "CRPIX1") + 
                                     1])
      }
      else {
        message("Missing CRPIX1")
      }
      if ("CRPIX2" %in% header) {
        CRPIX2 = as.numeric(header[which(header == "CRPIX2") + 
                                     1])
      }
      else {
        message("Missing CRPIX2")
      }
      if ("CD1_1" %in% header) {
        CD1_1 = as.numeric(header[which(header == "CD1_1") + 
                                    1])
        if ("CD1_2" %in% header) {
          CD1_2 = as.numeric(header[which(header == "CD1_2") + 
                                      1])
        }
        else {
          message("Missing CD1_2")
        }
      }
      else {
        if ("CDELT1" %in% header) {
          CD1_1 = as.numeric(header[which(header == "CDELT1") + 
                                      1])
        }
        else {
          message("Missing CD1_1 and CDELT1")
        }
      }
      if ("CD2_2" %in% header) {
        CD2_2 = as.numeric(header[which(header == "CD2_2") + 
                                    1])
        if ("CD2_1" %in% header) {
          CD2_1 = as.numeric(header[which(header == "CD2_1") + 
                                      1])
        }
        else {
          message("Missing CD2_1")
        }
      }
      else {
        if ("CDELT1" %in% header) {
          CD2_2 = as.numeric(header[which(header == "CDELT2") + 
                                      1])
        }
        else {
          message("Missing CD2_2 and CDELT2")
        }
      }
    }
  }
  RA0 = CRVAL1
  Dec0 = CRVAL2
  x0 = CRPIX1
  y0 = CRPIX2
  x1 = CD1_1
  x2 = CD1_2
  y1 = CD2_1
  y2 = CD2_2
  RA0 = RA0 * (pi/180)
  Dec0 = Dec0 * (pi/180)
  scalemat = matrix(c(x1, x2, y1, y2), 2) * (pi/180)
  xytran = cbind(x - x0, y - y0) %*% scalemat
  x = xytran[, 1]
  y = xytran[, 2]
  rad = sqrt(x^2 + y^2)
  if (grepl("TAN", CTYPE1)) {
    radproj1 = atan(rad)
  }
  else if (grepl("SIN", CTYPE1) | grepl("NCP", CTYPE1)) {
    if (grepl("NCP", CTYPE1)) {
      message("Approximating deprecated CTYPE1 NCP with SIN!")
    }
    radproj1 = asin(rad)
  }
  else {
    stop("Projection system is not recognised. Must be either TAN, SIN or NCP!")
  }
  if (grepl("TAN", CTYPE2)) {
    radproj2 = atan(rad)
  }
  else if (grepl("SIN", CTYPE2) | grepl("NCP", CTYPE2)) {
    if (grepl("NCP", CTYPE2)) {
      message("Approximating deprecated CTYPE2 NCP with SIN!")
    }
    radproj2 = asin(rad)
  }
  else {
    stop("Projection system is not recognised. Must be either TAN, SIN or NCP!")
  }
  rafunc = function(RA0, Dec0, x, y) {
    RA0 + atan2(x * sin(radproj1), rad * cos(Dec0) * cos(radproj1) - 
                  y * sin(Dec0) * sin(radproj1))
  }
  decfunc = function(Dec0, x, y) {
    asin(cos(radproj2) * sin(Dec0) + (y * sin(radproj2) * 
                                        cos(Dec0)/rad))
  }
  RA = rafunc(RA0, Dec0, x, y) * 180/pi%%360
  Dec = decfunc(Dec0, x, y) * 180/pi%%90
  Dec[which(is.nan(Dec))] = Dec0 * 180/pi
  output = cbind(as.numeric(RA), as.numeric(Dec))
  colnames(output) = c("RA", "Dec")
  return(output)
}




