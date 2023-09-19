# utility functions (miscellaneous but useful)

getMomentFromSlip = function(slips, rigidity=4*10^10, doTaper=FALSE, lambda=1, fault=csz, depths=getFaultCenters(fault)[,3], 
                             normalizeTaper=FALSE, dStar=28000) {
  # get fault total area
  areas = fault$length*fault$width
  totArea = sum(areas)
  
  if(doTaper)
    slips = slips*taper(depths, lambda, dStar=dStar, normalize=normalizeTaper)
  
  Mo = sum(slips*rigidity*areas)
  
  if(Mo <= 0) {
    warning("negative seismic moment observed")
    return(0)
  }
  
  (log10(Mo) - 9.05)/1.5
}

# Projects from lon/lat (EPSG:4326) to utm coordinates in km based on UTM 10 
# (EPSG:32610) or back if inverse==TRUE.
# x: a matrix of points, where each row is the coordinates of a single point
# inverse: if FALSE, converts from lon/lat to utm, else does the reverse
# units: currently not used
projCSZ = function(x, inverse=FALSE, units=c("m", "km")) {
  units = match.arg(units)
  
  # determine the "to" projection systems
  utmProj = st_crs("EPSG:32610")
  lonLatProj = st_crs("EPSG:4326")
  
  if(!inverse) {
    toProj = utmProj
    fromProj = lonLatProj
  } else {
    toProj = lonLatProj
    fromProj = utmProj
  }
  
  # convert to an sf object
  x = st_multipoint(x, dim="XY")
  if(inverse) {
    # make sure we convert from km back to m as the projection is expecting
    x = x*1000
  }
  x = st_sfc(x, crs=fromProj)
  
  # we already know it is utm zone 10
  if(FALSE) {
    # calculate UTM zone based on slab 2 geometry
    slab = loadSlab2()
    lonRange = range(slab$lon)
    latRange = range(slab$lat)
    midLon = mean(lonRange)
    
    utmZone = ceiling((midLon-180)/6) # UTM zone 10
  }
  
  # transform coordinates and convert back to a matrix
  out = st_transform(x, toProj)
  out = st_coordinates(out)[,1:2]
  
  # convert from m to km if necessary
  if(toProj == utmProj) {
    out = out/1000
  }
  
  out
}




