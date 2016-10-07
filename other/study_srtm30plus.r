#----
# initialize environment

# set working directory
wd = 'E:/bbest/CetMap'
setwd(wd)

# load packages
require('maptools') # readShapePoly
require('raster') # rasterFromXYZ
#require('RCurl')

#----
# read in bounding box extent
eez = readShapePoly('data/study/study_bbox_gcs.shp', proj4string=CRS('+proj=longlat +datum=WGS84'))
#summary(eez)
b = bbox(eez) # x:-98,-60; y:23,46
apply(bbox(eez),1,mean)# x: -79; y: 34.5 # centroid
#b = c( floor(b['x','min']), ceiling(b['x','max']), floor(b['y','min']), ceiling(b['y','max']) )

#----
# formulate URLs for UCSD Topex site and setup for append to XYZ file on per degree cell basis within extent

# get integer cell extents, base URL and XYZ file headers
u0 = 'http://topex.ucsd.edu/cgi-bin/get_srtm30.cgi?submit=get+data'
csv = 'data/bathy/xyz_test.csv'
xs = seq(floor(b['x','min']), ceiling(b['x','max']), by=1)
ys = floor(b['y','min']):ceiling(b['y','max'])
d = data.frame(x=numeric(),y=numeric(),z=integer())
write.csv(d, csv, row.names=F)

# iterate over x and y to fetch and append
for (ix in 1:(length(xs)-1)){
  for (iy in 1:(length(ys)-1)){  
    if (!(xs[ix] < -82 & ys[iy] > 32)){ # NOTE: hacked exceptions here to skip land areas

      # construct request URL
      u = sprintf('%s&west=%g&east=%g&south=%g&north=%g',u0,xs[ix],xs[ix+1],ys[iy],ys[iy+1])
      print(u) # show progress

      # write as read so memory doesn't get used up
      write.table(read.delim(url(u), header=F, col.names=c('x','y','z')), 
                  file=csv, append=T, sep=',', row.names=F, col.names=F)
    }
  }
}

#---
# read in data as raster from XYZ
r = rasterFromXYZ(d, crs='+proj=longlat +datum=WGS84') # res=c(NA,NA), digits=5)


#----
# plot
# Error: cannot allocate vector of size 83.4 Mb
spplot(r)
plot(r, col=colorRampPalette(c("red", "white", "blue"))(255))
