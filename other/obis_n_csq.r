# todo:
#  * check/resolve points on the line?

# INITIALIZE VARIABLES

# packages to load (and install if needed)
packages = c('sp','rgdal','maps','fields','maptools','PBSmapping')

# url prefixes 
hits.pre = 'http://iobis-gis.marine.rutgers.edu/geoserver/wfs?REQUEST=GetFeature&ResultType=hits&VERSION=1.0.0&TYPENAME=obis:expandeddigircache&BBOX=' 
wfs.pre = 'http://iobis-gis.marine.rutgers.edu/geoserver/wfs?REQUEST=GetFeature&VERSION=1.0.0&TYPENAME=obis:expandeddigircache&BBOX='
wms.pre = 'http://iobis-gis.marine.rutgers.edu/geoserver/wms?REQUEST=GetMap&LAYERS=obis:bluemarble,obis:expandeddigircache&FORMAT=image/png&WIDTH=500&HEIGHT=500&SRS=EPSG:4326&BBOX='

## wfs hits vars
hits.tmp = tempfile()
hits.max = 10000 # (GeoServer default for WFS max record return is 10,000)

# csquares spec vars
resolutions = c(100,10,5,1,0.5,0.1,0.05,0.1,0.05,0.01) # in decimal degrees, output 1, and 4 to 16 characters, per c-squares spec

# load required packages, installing if need be
for (p in packages){
	if (!p %in% installed.packages()[,'Package']){
		print(sprintf('installing: %s',p))
		install.packages(p, repos=repos, dependencies=T); Sys.sleep(2)
    }
	require(p, character.only=T)
}

# UTILITY FUNCTIONS

bb2poly = function(bb, ID='bbox', proj=CRS('+proj=longlat +datum=WGS84'), data=NULL) {
  # convert a bounding box defined as c(lonmin, latmin, lonmax, latmax) to a SpatialPolygon object
  #   or optionally adding data to a data frame, hence SpatialPolygonDataFrame
  poly = Polygon(data.frame(rbind(c(bb[1], bb[2]),
                                c(bb[1], bb[4]),
                                c(bb[3], bb[4]),
                                c(bb[3], bb[2]),
                                c(bb[1], bb[2]))))
  polys = SpatialPolygons(list(Polygons(list(poly), ID=ID)), proj4string=proj)
  if (is.null(data)){
  	return(polys)
  } else {
  	return(SpatialPolygonsDataFrame(polys, data=data.frame(data,row.names=ID)))
  }
}

get.hits.by.bb = function(bb){
	# determine hits given:
	#  input: bounding box (bb) as string of form "x.min, y.min, x.max.x, y.max" per BBOX WMS/WFS spec
	#  output: integer of hits for query
	#  example: get.hits('0,0,10,10')
	hits.url = paste(hits.pre, paste(bb, collapse=','), sep='')
	response = try(download.file(hits.url, hits.tmp))
	if (class(response) == 'try-error'){
		hits = Inf  # hack for now, since presume download timing out b/c of large am't of data
	} else {
		hits = readLines(hits.tmp, warn=F)[2]
		hits = as.integer(strsplit(hits,'"')[[1]][2])
	}
	return(hits)
}

get.finer.csquares = function(csq, id){
	# update csquares SpatialPolygonDataFrame with next finer resolution of csquares

	# get easy reference to string and bounding box of current csquare
	idx = which(rownames(csq@data) == id)
	b = as.vector(bbox(csq@polygons[[idx]]))
	s = csq@polygons[[idx]]@ID
	s1 = substr(s,1,1)
	csq@data[id,'use'] = F # check current csquare as unused, to be replaced with finer resolution csquares
	
	# get next resolution (res) based on number of characters in csquares id string
	if (nchar(s) == 1){
		res = resolutions[2]
	} else {
		res = resolutions[nchar(s)/2+1] # eg 6 characters = 2nd resolution, or 5 degree
	}
	
	# get sequence of geographic x and y along distance from origin (0,0), 
	#   according to quadrant and resolution
	if (s1=='1'){          # NE quadrant
		gx = seq(b[1], b[3], res)
		gy = seq(b[2], b[4], res)
	} else if (s1=='3'){   # SE quadrant
		gx = seq(b[1], b[3], res)
		gy = seq(b[4], b[2],-res)
	} else if (s1=='5'){   # SW quadrant
		gx = seq(b[3], b[1],-res)
		gy = seq(b[4], b[2],-res)
	} else if (s1=='7'){   # NW quadrant
		gx = seq(b[3], b[1],-res)
		gy = seq(b[2], b[4], res)
	}
	nx = length(gx)-1
	ny = length(gy)-1
	# intermediate quadrant contribution per x and y location
	qx = c(rep(1,nx/2), rep(2,nx/2))
	qy = c(rep(0,nx/2), rep(2,nx/2))

	for (ix in 1:nx){  # iterate by row 
		for (iy in 1:ny){ # then by col
		
			# get bounding box
			rx = gx[c(ix,ix+1)]
			ry = gy[c(iy,iy+1)]
			bb = c(min(rx), min(ry), max(rx), max(ry))
			
			# define string
			q = qx[ix] + qy[iy] # cell quadrant
			if (res == 10){               # add 10-degree cell
				string = sprintf('%s%d%02d',s,iy-1,ix-1)
			} else if (ns %% 4 == 2) {	# add intermediate quadrant, if modulus 2
				string = sprintf('%s:%d',s,q)
			} else {					# add x,y cell identifiers for 3-digit cycle
				string = sprintf('%s%d%d',s,iy-1,ix-1)
			}
			
			# get hits
			print(paste('getting hits for csquare', string,'bbox', paste(bb, collapse=',')))
			t0 = Sys.time()
			hits   = get.hits.by.bb(bb)
			
			# append to csq polygons
			csq = rbind(csq, bb2poly(bb, ID=string, data=list(hits=hits,use=T,res=res,time=Sys.time()-t0)))
		}
	}
	return(csq)
}

# INITIALIZE CSQUARES as SpatialPolygonsDataFrame with global quadrants
d = list(hits=Inf,use=T,res=100,time=0)  # set inital hits to Infinity so goes to next finer level in while loop
csq = rbind(bb2poly(c(   0,  0,180,90), ID='1', data=d),
			bb2poly(c(   0,-90,180, 0), ID='3', data=d),
			bb2poly(c(-180,-90,  0, 0), ID='5', data=d),
			bb2poly(c(-180,  0,  0,90), ID='7', data=d))


# ITERATE into finer c-squares resolution where reaching hits.max limit
while (max(csq@data[csq@data$use==T,'hits']) >= hits.max){
	id = rownames(csq@data)[which(csq$hits >= hits.max & csq$use==T)[1]]
	csq = get.finer.csquares(csq,id)
}

# PLOT CSQUARES

# plot initial csquares quadrant world map
csq.m = csq[which(csq$res==100),]
map('world', col='grey', xlim=bbox(csq.m)[1,], ylim=bbox(csq.m)[2,])
plot(csq.m, add=T)
text(coordinates(csq.m), labels=rownames(csq.m@data), cex=5)

# plot csq map within region of given csq string
csq.id = '1'
bb = bbox(csq@polygons[[which(rownames(csq@data)==csq.id)]])
map('world', col='grey', xlim=bb[1,], ylim=bb[2,])
plot(csq[which(csq$use==T),], add=T)

# plot full csq map, coloring by density
csq$area = sapply(csq@polygons, function(i) i@area)
csq$dens = log(csq$hits / csq$area)
csq.m = csq[which(csq$use==T),]   # subset
nbrks = 7
brks = quantile(csq.m$dens, seq(0,1,1/(nbrks-1)), na.rm=T)
cols = tim.colors(nbrks)[findInterval(csq.m$dens, brks)] 
plotPolys(SpatialPolygons2PolySet(csq.m), col=cols)
map('world', fill=T, col='grey', add=T)
legend(c(-180,-150), c(-60,-90), fill=tim.colors(nbrks), legend=leglabs(format(brks,digits=2)), bty="o",bg='white',cex=1.7,title=expression(paste('obs density (n/',degree,')')))

# generate KML Google Earth overlay
d = '/Users/bbest/Desktop'
f = 'obis_n_csq_'
kml = sprintf('%s/%s.kml',d,f)
png = sprintf('%s/%s.png',d,f)
ge = GE_SpatialGrid(csq.m, maxPixels=5000)
png(file=png, width=ge$width, height=ge$height, bg="transparent")
par(mar=c(0,0,0,0), xaxs="i", yaxs="i")
plot(csq.m, col=cols[idx], xlim=ge$xlim, ylim=ge$ylim)
map('world', fill=T, col='grey', add=T)
dev.off()
kmlOverlay(ge, kml, sprintf('./%s',basename(png)))


# DOWNLOAD WFS DATA BY CSQUARES

# initalize vars
wfs.tmp = tempfile()
pts.df = data.frame(datelastmodified=character(), institutioncode=character(), collectioncode=character(), scientificname=character(), 
           catalognumber=character(), basisofrecord=character(), kingdom=character(), phylum=character(), class=character(), 
           ordername=character(), family=character(), genus=character(), species=character(), subspecies=character(), scientificnameauthor=character(), 
           collector=character(), latitude=numeric(), longitude=numeric(), slatitude=numeric(),
           recordlastcached=character(), res_name=character(), recordurl=character(), csquare=character(), locality=character(), 
           country=character(), fieldnumber=character(), individualcount=integer(), gid=integer(), yearcollected=integer(),
           monthcollected=integer(), daycollected=integer(), minimumdepth=numeric(), maximumdepth=numeric(), source=character(), 
           citation=character(), datecollected=character(), avgdep=numeric(), deperr=numeric(), identifiedby=character(), 
           julianday=integer(), typestatus=character(), coordinateprecision=numeric(), timeofday=numeric(), timezone=character(),
           stringsAsFactors=F)  # unique by institutioncode, collectioncode, catalognumber? remove duplicates
stats = data.frame(bbox=character(), hits=numeric(), time=character(), url.hits=character(), url.wfs=character(), url.wms=character())

idx = which(csq$use==T)
csq.polys = csq[idx,]@polygons
csq.data = csq[idx,]@data

for (i in 1:length(csq.polys)){
	p = csq.polys[[i]]
	bb = as.vector(bbox(p))
	dat = csq.data[i,]
	if (dat$hits > 0){
						
		# download WFS and read in as GML
		wfs.url  = paste(wfs.pre, bb, sep='')
		wms.url  = paste(wms.pre, bb, sep='')
		t0 = Sys.time()
		download.file(wfs.url, wfs.tmp)
		t = Sys.time()-t0
		pts.tmp = readOGR(wfs.tmp, 'expandeddigircache')
			
		# modify/add/delete fields so match
		for (f in names(pts.df)){
			if (f %in% names(pts.tmp)){
				# convert to same data type
				pts.tmp@data[,f] = as( pts.tmp@data[,f], class(pts.df[,f]) )
			} else {
				# add field with NAs
				pts.tmp@data[,f] = rep(NA, nrow(pts.tmp))
			}
		}
		# exclude fields
		pts.tmp@data = pts.tmp@data[, which(names(pts.tmp) %in% names(pts.df))]
			
		# append or set new			
		if (exists('pts')){
			pts = rbind(pts, pts.tmp)
		} else {
			pts = pts.tmp
		}
	}
		
	# track fetch stats
	stats = rbind(stats, data.frame(bbox=bb, hits=hits, time=t,  url.hits=hits.url, url.wfs=wfs.url, url.wms=wms.url))
	}
}

# SUMMARIZE AND PLOT WFS DATA BY CSQUARES

# text summary
summary(pts)
dim(pts)

# plot map
map('world', xlim=bf[c(1,3)], ylim=bf[c(2,4)])
points(pts)

# plot fetch stats
hist(stats$hits)
plot(stats$time, stats$hits)

# FUTURE...

# look at date field for observations through time

# write to pgsql db, either with RODBC or ogr2ogr

#  * better to paginate, and have general xml data srxr, not precast
#  * wrap as kepler tool.  could even do postgresql pl/r
# long term:
#   * Catalog of full dataset: GeoNetwork (OGC Catalog Service) vs MetaCat, GCMD
#   * 3-tiered mapping, prob as pgsql procedure with numeric sorted csqu indecies

