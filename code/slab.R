# script for loading and discretizing Slab 2.0 data

# loads the Slab 2.0 model output points and depths
loadSlab2 = function() {
  slabDir = "data/cas_surf_09.21_slab2_output/"
  
  # doesn't work/bad depth values due to NAs:
  # slab = read.csv(paste0(slabDir, "cas_slab2_nod_09.04.23.csv"))
  # slab
  
  require(ncdf4)
  ncdat <- nc_open("data/cas_surf_09.21_slab2_output/cas_slab2_dep_09.04.23.grd")
  depth = ncvar_get(ncdat, "z")
  lon = ncvar_get(ncdat, "x")
  lat = ncvar_get(ncdat, "y")
  
  allLats = rep(lat, each=length(lon))
  allLons = rep(lon, length(lat))
  allDepths = c(depth)
  
  ncdat <- nc_open("data/cas_surf_09.21_slab2_output/cas_slab2_dip_09.04.23.grd")
  dip = ncvar_get(ncdat, "z")
  allDips = c(dip)
  
  ncdat <- nc_open("data/cas_surf_09.21_slab2_output/cas_slab2_str_09.04.23.grd")
  strike = ncvar_get(ncdat, "z")
  allStrikes = c(strike)
  
  out = data.frame(list(lon=allLons, lat=allLats, depth=allDepths, dip=allDips, strike=allStrikes))
  
  extentI = is.finite(allDepths)
  
  out[extentI,]
}

# Discretizes the Slab 2.0 geometry into a triangulation mesh.
# NOTE: this function does not currently vary the resolution of the mesh as a 
#     function of the depth, since this will likely not substantially impact the 
#     subsidence along the coast.
# n: number of nodes to form triangulation from (each node is the corner of 
#     several triangles)
# max.edge: maximum triangle edge length in km (for making triangle regular)
# maxDepth: maximum depth of the fault geometry
# ...: other inputs passed to inla.mesh.2d
discretizeSlab2 = function(n=2000, max.n=-1, max.edge=c(15, 100), maxDepth=30, ...) {
  
  # first load in the Slab 2.0 geometry
  slab = loadSlab2()
  lonLat = cbind(slab$lon, slab$lat)
  depths = slab$depth
  
  # get projected coordinates
  xy = projCSZ(lonLat)
  
  if(FALSE) {
    plotWithColor(xy[,1], xy[,2], depths, pch=19, cex=.3, xlab="Easting (km)", 
                  ylab="Northing (km)")
  }
  
  # keep only points with appropriate depth
  goodI = abs(depths) <= maxDepth
  lonLat = lonLat[goodI,]
  xy = xy[goodI,]
  depths = depths[goodI]
  
  # compute convex hull of slab geometry
  # hullI = chull(xy)
  # xyHull = xy[hullI,]
  hullInt = inla.nonconvex.hull.basic(xy, resolution=350, convex=-.012)
  xyHull = hullInt$loc
  
  hullExt = inla.nonconvex.hull.basic(xy, resolution=150, convex=-.4)
  
  if(FALSE) {
    plotWithColor(xy[,1], xy[,2], depths, pch=19, cex=.3, xlab="Easting (km)", 
                  ylab="Northing (km)")
    polygon(xyHull[,1], xyHull[,2], border="green")
  }
  
  # construct mesh with INLA
  # mesh = inla.mesh.2d(n=n, loc.domain=xyHull, max.edge=max.edge, ...)
  # mesh = inla.mesh.2d(n=n, boundary=hullInt, max.edge=max.edge, ...)
  mesh = inla.mesh.2d(n=n, boundary=list(hullInt, hullExt), max.n=max.n, max.edge=max.edge, ...)
  
  if(FALSE) {
    # plot the mesh
    plot(mesh, asp=1)
    # points(xy[,1], xy[,2], col="red", pch=".")
  }
  
  faultGeom = getGeomFromMesh(mesh, extent=hullInt$loc, maxDepth=maxDepth)
  browser()
  faultGeom
}

getGeomFromMesh = function(mesh, extent, maxDepth=30) {
  corners = mesh$loc[,1:2] # corners of the triangles, i.e. vertices
  
  # t: triangles, v: vertices
  tv = mesh$graph$tv
  vt = mesh$graph$vt
  tt = mesh$graph$tt
  tti = mesh$graph$tti
  vv = mesh$graph$vv # sparse matrix, 1 if connected
  
  nV = nrow(corners)
  nT = nrow(tv)
  
  inds = vt
  
  # for each triangle, calculate its center (for calculating spatial 
  # covariances)
  centers = matrix(nrow=nT, ncol=2)
  triCorners = list()
  for(ti in 1:nT) {
    vInds = tv[ti,]
    thisCoords = corners[vInds,]
    
    centers[ti,] = colMeans(thisCoords)
    triCorners = c(triCorners, list(thisCoords))
  }
  
  # get depths at all points
  allDepths = getPointDepths(rbind(centers, corners), maxDepth=maxDepth)
  centerDepths = allDepths[1:nrow(centers)]
  allCornerDepths = allDepths[-(1:nrow(centers))]
  
  # add depths to coordinates
  centers = data.frame(list(lon=centers[,1], lat=centers[,2], depth=centerDepths))
  for(i in 1:nT) {
    vInds = tv[ti,]
    thisDepths = allCornerDepths[vInds]
    triCorners[[i]] = data.frame(lon=triCorners[[i]][,1], lat=triCorners[[i]][,2], 
                                 depth=thisDepths)
  }
  
  # determine which triangles are in the fault extent
  internalI = fields::in.poly(centers, extent)
  internalCenters = centers[internalI,]
  externalCenters = centers[!internalI,]
  internalTriCorners = list()
  externalTriCorners = list()
  for(i in 1:nrow(centers)) {
    thisTriCorners = triCorners[[i]]
    if(internalI[i]) {
      internalTriCorners = c(internalTriCorners, list(thisTriCorners))
    } else {
      externalTriCorners = c(externalTriCorners, list(thisTriCorners))
    }
  }
  
  list(centers=internalCenters, corners=internalTriCorners, 
       externalTriangulation=list(centers=externalCenters, corners=externalTriCorners))
}

# simple nearest neighbor algorithm
# assume pts is in utm10
getPointDepths = function(pts, maxDepth=Inf) {
  # first load in the Slab 2.0 geometry
  slab = loadSlab2()
  lonLat = cbind(slab$lon, slab$lat)
  depths = slab$depth
  
  # get projected coordinates
  xy = projCSZ(lonLat)
  
  # keep only points with appropriate depth
  goodI = abs(depths) <= (maxDepth + 1)
  xy = xy[goodI,]
  depths = depths[goodI]
  
  # system.time(out <- interp::interp(x=xy[,1], y=xy[,2], z=depths, 
  #                                   xo=pts[,1], yo=pts[,2], output="points", 
  #                                   method="akima"))
  # 
  # system.time(out <- interp::interp(x=xy[,1], y=xy[,2], z=depths, 
  #                                   xo=pts[,1], yo=pts[,2], output="points", 
  #                                   method="akima"))
  
  distMat = rdist(pts, xy)
  nearestInds = apply(distMat, 1, which.min)
  out = depths[nearestInds]
  
  out
}






