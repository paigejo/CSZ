# testing script

# get two triangular subfaults forming a rectangle for testing 
# the okada model on
getTestTriangularSubfaults = function(nDown=150, nStrike=200, strike=0, dip=14) {
  
  # get old fault geometry data (convert from km to m and to longitude format consistent with other data)
  faultGeom = read.csv("data/CSZe01.csv")
  faultGeom$longitude = faultGeom$longitude - 360
  kmCols = c(4, 6, 7)
  faultGeom[,kmCols] = faultGeom[,kmCols] * 10^3
  
  csz = divideFault(faultGeom, nDown=nDown, nStrike=nStrike)
  rectSubfault = csz[1,]
  rectSubfault$strike = strike
  rectSubfault$dip = dip
  
  # extract first rectangular subfault
  geom = calcGeom(rectSubfault)
  rectCorners = geom$corners
  
  # divide the rectangle into 2 triangles
  tri1 = rectCorners[1:3,]
  tri2 = rectCorners[c(1, 3, 4),]
  
  # construct testing points
  testLon = range(rectCorners[,1])
  testLat = range(rectCorners[,2])
  
  list(rect=rectSubfault, tri1=tri1, tri2=tri2, testLon=testLon, testLat=testLat, rectGeom=geom)
}

# obtain okada model predictions for the first test triangular subfault
testTri1 = function() {
  # get the subfault. It must contain the following:
  #     fix_orientation: whether the orientation is fixed. Set by 
  #                      calculate_geometry_triangles
  #     corners: 3x3 matrix where each row is (lon, lat, depth) of a single corner
  #     centers: TODO
  #     lon: longitude TODO
  #     lat: latitude TODO
  #     depth: depth TODO
  #     strike: strike TODO
  #     dip: dip TODO
  #     slip: average coseismic slip for the subfault
  #     rake: rake of the subfault
  tri1 = getTestTriangularSubfaults()$tri1
  
  triGeom = calculate_geometry_triangles(tri1)
  
  # set test points
  lonTest = c(-124.2, -124.1)
  latTest = c(40.15, 40.30)
  
  # run okada model
  dtopo = okadaSubfaultTri(triGeom, x=lonTest, y=latTest)
  
  # plot results
  xlim = range(triGeom$corners[,1])
  ylim = range(triGeom$corners[,2])
  plotWithColor(triGeom$corners[,1], triGeom$corners[,2], triGeom$corners[,3], main="depth", 
                xlim=xlim, ylim=ylim, pch=19)
  
  plotWithColor(c(dtopo$X), c(dtopo$Y), c(dtopo$dX), main="dx", 
                xlim=xlim, ylim=ylim, pch=19)
  points(triGeom$corners[,1], triGeom$corners[,2])
  
  plotWithColor(c(dtopo$X), c(dtopo$Y), c(dtopo$dY), main="dy", 
                xlim=xlim, ylim=ylim, pch=19)
  points(triGeom$corners[,1], triGeom$corners[,2])
  
  plotWithColor(c(dtopo$X), c(dtopo$Y), c(dtopo$dZ), main="dz", 
                xlim=xlim, ylim=ylim, pch=19)
  points(triGeom$corners[,1], triGeom$corners[,2])
  
  # do the results change after reordering the corners?
  tri12 = tri1
  tri12 = tri12[c(2, 3, 1),]
  triGeom2 = calculate_geometry_triangles(tri12)
  dtopo2 = okadaSubfaultTri(triGeom2, x=lonTest, y=latTest)
  
  browser()
  triGeom
  triGeom2
  
  dtopo
  dtopo2
}

testTrianglesVsRect = function(nDown=150, nStrike=200, slip=100) {
  out = getTestTriangularSubfaults(nDown=nDown, nStrike=nStrike)
  tri1 = out$tri1
  tri2 = out$tri2
  rect = out$rect
  lonTest = out$testLon
  latTest = out$testLat
  rectGeom = out$rectGeom
  
  # project corners into E/N and calculate geometries
  triGeom1 = calculate_geometry_triangles(cbind(projCSZ(tri1[,1:2], units="m"), tri1[,3]))
  triGeom2 = calculate_geometry_triangles(cbind(projCSZ(tri2[,1:2], units="m"), tri2[,3]))
  # triGeom2 = calculate_geometry_triangles(tri2)
  
  # convert x/y back to lon/lat
  triGeom1$corners = tri1
  triGeom2$corners = tri2
  lonLat1 = projCSZ(cbind(triGeom1$lon, triGeom1$lat), inverse=TRUE, units="m")
  lonLat2 = projCSZ(cbind(triGeom2$lon, triGeom2$lat), inverse=TRUE, units="m")
  triGeom1$lon = lonLat1[1]
  triGeom1$lat = lonLat1[2]
  triGeom2$lon = lonLat2[1]
  triGeom2$lat = lonLat2[2]
  
  triGeom1$slip = slip
  triGeom2$slip = slip
  
  triGeomFull = list(triGeom1, triGeom2)
  
  # set test points
  # lonTest = c(-124.2, -124.1)
  # latTest = c(40.15, 40.30)
  
  # run okada model for individual triangular subfaults and 
  # full triangulated fault
  dtopo1 = okadaSubfaultTri(triGeom1, x=lonTest, y=latTest)
  dtopo2 = okadaSubfaultTri(triGeom2, x=lonTest, y=latTest)
  dtopoAll = okadaTri(triGeomFull, x=lonTest, y=latTest, slip=slip)
  
  # run okada model for individual rectangular subfault
  dtopoRect = okadaSubfaultRect(rect, x=lonTest, y=latTest, slip=slip)
  
  browser()
  
  print("2 triangular Okada:")
  print(dtopoAll$dZ)
  print("rectangular Okada:")
  print(dtopoRect)
}

testFaultDepthCalc = function() {
  browser()
  triGeom1 = discretizeSlab2(method="kernel")
  triGeom2 = discretizeSlab2(method="NN")
  
  for(i in 1:length(triGeom1$corners)) {
    triGeom1$corners[[i]][,1:2] = projCSZ(as.matrix(triGeom1$corners[[i]][,1:2]), inverse=TRUE, units="km")
  }
  for(i in 1:length(triGeom2$corners)) {
    triGeom2$corners[[i]][,1:2] = projCSZ(as.matrix(triGeom2$corners[[i]][,1:2]), inverse=TRUE, units="km")
  }
  
  if(FALSE) {
    centerDepths1 = triGeom1$centers[,3]
    plotPolyDat(triGeom1$corners, centerDepths1, border=rgb(0,0,0,0))
    
    plotPolyDat(triGeom1$corners, centerDepths1, border=rgb(0,0,0,0), project=TRUE, myProjection=projCSZ)
    # plotPolyDat(triGeom1$corners, centerDepths1, border=rgb(0,0,0,0), project=TRUE, myProjection=projCSZ2)
    
    centerDepths2 = triGeom2$centers[,3]
    plotPolyDat(triGeom2$corners, centerDepths2, border=rgb(0,0,0,0))
  }
  
  browser()
  
  triGeomFull1 = getFullFaultGeom(triangulatedGeom=triGeom1)
  triGeomFull2 = getFullFaultGeom(triangulatedGeom=triGeom2)
  
  # plot fault
  if(FALSE) {
    centerDepths1 = sapply(triGeomFull1, function(x) {x$depth})
    strikes1 = sapply(triGeomFull1, function(x) {x$strike})
    dips1 = sapply(triGeomFull1, function(x) {x$dip})
    
    centerDepths2 = sapply(triGeomFull2, function(x) {x$depth})
    strikes2 = sapply(triGeomFull2, function(x) {x$strike})
    dips2 = sapply(triGeomFull2, function(x) {x$dip})
    
    plotPolyDat(triGeom1$corners, centerDepths1, border=rgb(0,0,0,0))
    plotPolyDat(triGeom1$corners, strikes1, border=rgb(0,0,0,0))
    plotPolyDat(triGeom1$corners, dips1, border=rgb(0,0,0,0))
    
    plotPolyDat(triGeom2$corners, centerDepths2, border=rgb(0,0,0,0))
    plotPolyDat(triGeom2$corners, strikes2, border=rgb(0,0,0,0))
    plotPolyDat(triGeom2$corners, dips2, border=rgb(0,0,0,0))
  }
}

testFullFault = function() {
  browser()
  triGeom = discretizeSlab2(method="kernel")
  
  # convert coordinates back to lon/lat as expected by getFullFaultGeom()
  for(i in 1:length(triGeom$corners)) {
    triGeom$corners[[i]][,1:2] = projCSZ(as.matrix(triGeom$corners[[i]][,1:2]), inverse=TRUE, units="km")
  }
  
  triGeomFull = getFullFaultGeom(triangulatedGeom=triGeom)
  
  # plot fault
  if(FALSE) {
    centerDepths = sapply(triGeomFull, function(x) {x$depth})
    strikes = sapply(triGeomFull, function(x) {x$strike})
    dips = sapply(triGeomFull, function(x) {x$dip})
    
    plotPolyDat(triGeom$corners, centerDepths, border=rgb(0,0,0,0))
    plotPolyDat(triGeom$corners, cornerDepth1, border=rgb(0,0,0,0))
    plotPolyDat(triGeom$corners, cornerDepth2, border=rgb(0,0,0,0))
    plotPolyDat(triGeom$corners, cornerDepth3, border=rgb(0,0,0,0))
    plotPolyDat(triGeom$corners, strikes, border=rgb(0,0,0,0))
    plotPolyDat(triGeom$corners, dips, border=rgb(0,0,0,0))
  }
  
  # get lon/lat range
  allExtCorners = do.call("rbind", triGeom$externalTriangulation$corners)
  lonLat = projCSZ(as.matrix(allExtCorners[,1:2]), inverse=TRUE)
  lonRange = range(lonLat[,1])
  latRange = range(lonLat[,2])
  
  # run Okada model on full fault
  lonTest = seq(lonRange[1], lonRange[2], l=100)
  latTest = seq(latRange[1], latRange[2], l=100)
  
  totTime = system.time(dtopoAll <- okadaTri(triGeomFull, x=lonTest, y=latTest, slip=10))[3]
  # totTime/60
  # elapsed 
  # 4.178517
  
  depths = sapply(triGeomFull, function(x) {x$depth})
  strikes = sapply(triGeomFull, function(x) {x$strike})
  dips = sapply(triGeomFull, function(x) {x$dip})
  plotPolyDat(triGeom$corners, depths, border=rgb(0,0,0,0), 
              xlab="Easting (km)", ylab="Northing (km)", addColorBar = TRUE, 
              leaveRoomForLegend = TRUE)
  
  plotPolyDat(triGeom$corners, strikes, border=rgb(0,0,0,0), 
              xlab="Easting (km)", ylab="Northing (km)", addColorBar = TRUE, 
              leaveRoomForLegend = TRUE)
  
  plotPolyDat(triGeom$corners, dips, border=rgb(0,0,0,0), 
              xlab="Easting (km)", ylab="Northing (km)", addColorBar = TRUE, 
              leaveRoomForLegend = TRUE)
  
  plotPolyDat(triGeom$corners, border=rgb(0, 0, 0, 0.5), lwd=.3, 
              xlab="Easting (km)", ylab="Northing (km)", addColorBar = FALSE, 
              leaveRoomForLegend = TRUE)
  
  ENTest = projCSZ(cbind(c(dtopoAll$X), c(dtopoAll$Y)))
  
  plotWithColor(c(dtopoAll$X), c(dtopoAll$Y), c(dtopoAll$dZ), 
                pch=19, cex=.1, asp=1)
  # plotPolyDat(triGeom$corners, border=rgb(0, 0, 0, 0.05), lwd=.3, new=FALSE)
  # map("usa", add=TRUE)
  world(add=TRUE)
  polygon(projCSZ(triGeom$extent, inverse=TRUE), border="purple")
  
  plotWithColor(ENTest[,1], ENTest[,2], c(dtopoAll$dZ), 
                xlab="Easting (km)", ylab="Northing (km)", 
                pch=19, cex=.1, asp=1)
  # plotPolyDat(triGeom$corners, border=rgb(0, 0, 0, 0.05), lwd=.3, new=FALSE)
  polygon(triGeom$extent, border="purple")
}


