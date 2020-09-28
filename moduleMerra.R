#' @title        Merra2 grid
#' @description  Returns Merra2 grid. If OPeNDAP params used, returns the equivalent coordinates
#' @param        lat_odap Latitude as used in OPeNDAP [0:1:360]. For debugging use.
#' @param        lon_odap Longitude as used in OPeNDAP [0:1:575]. For debugging use.
#' @param        rangeOnly returns a list of vectors with all the latitude and longitude
#' @return       data.frame
#' @note         source: http://goldsmr4.sci.gsfc.nasa.gov/dods/M2C0NXASM.das
gridMerra2 <- function(lat_odap = NA, lon_odap = NA, rangeOnly = F){
  lon_res <- 0.625
  lat_res <- 0.5
  longitude <- seq(-180, 179.375, by = lon_res)
  latitude  <- seq(-90, 90, by = lat_res)
  if(rangeOnly) return(list(lat = latitude, lon = longitude))
  if(all(is.na(c(lat_odap, lon_odap)))){
    grid <- expand.grid(lat = latitude, lon = longitude, KEEP.OUT.ATTRS = FALSE)
  } else {
    grid <- data.frame(lat = latitude[lat_odap + 1],  # OPeNDAP vector indices start at 0
                       lon = longitude[lon_odap + 1]) # and in R start at 1
  }

  attr(grid, 'lon_res') <- lon_res
  attr(grid, 'lat_res') <- lat_res
  return(grid)
}

#' @title       Generate Merra Links
#'
#' @description This function generate Merra daily file links from any of the 8 different products listed below. It allows
#'              to generate a single node or a box, at different time intervals. See the arguments for more details.
#'
#' @param       Product Character string. Copy & paste only one of the following options (To see all the available variables click the product name):\itemize{
#'              \item \href{http://disc.sci.gsfc.nasa.gov/mdisc/data-holdings/merra/inst1_2d_int_Nx.shtml}{"inst1_2d_int_Nx"} MERRA IAU 2d Vertical integrals
#'              \item \href{http://disc.sci.gsfc.nasa.gov/mdisc/data-holdings/merra/tavg1_2d_flx_Nx.shtml}{"tavg1_2d_flx_Nx"} MERRA IAU 2d surface turbulent flux diagnostics
#'              \item \href{http://disc.sci.gsfc.nasa.gov/mdisc/data-holdings/merra/tavg1_2d_int_Nx.shtml}{"tavg1_2d_int_Nx"} MERRA IAU 2d Vertical integrals
#'              \item \href{http://disc.sci.gsfc.nasa.gov/mdisc/data-holdings/merra/tavg1_2d_lnd_Nx.shtml}{"tavg1_2d_lnd_Nx"} MERRA IAU 2d land surface diagnostics
#'              \item \href{http://disc.sci.gsfc.nasa.gov/mdisc/data-holdings/merra/tavg1_2d_ocn_nx.shtml}{"tavg1_2d_ocn_Nx"} MERRA IAU 2d ocean surface diagnostics
#'              \item \href{http://disc.sci.gsfc.nasa.gov/mdisc/data-holdings/merra/tavg1_2d_rad_Nx.shtml}{"tavg1_2d_rad_Nx"} MERRA IAU 2d surface and TOA radiation fluxes
#'              \item \href{http://disc.sci.gsfc.nasa.gov/mdisc/data-holdings/merra/tavg1_2d_slv_Nx.shtml}{"tavg1_2d_slv_Nx"} MERRA IAU 2d atmospheric single-level diagnostics
#'              \item \href{http://disc.sci.gsfc.nasa.gov/mdisc/data-holdings/merra/tavg1_2d_mld_nx.shtml}{"tavg1_2d_mld_Nx"} MERRA-Land 2d land surface diagnostics
#'              \item \bold{Note:} A single file with all the above information can be found \href{https://www.dropbox.com/s/da5d2loyr7keyqz/Subset_2012.Nasa.File Specification for MERRA Products.pdf?dl=0}{here}
#'              }
#' @param       Variables Character vector. Desired variables are case sensitive. \itemize{
#'              \item \code{c('U2M','V2M')} Will work only if \code{Product = "tavg1_2d_slv_Nx"}
#'              }
#' @param       Time Character vector which can take several options. See the following examples:\itemize{
#'              \item \code{"2014-04-21"} Will generate all the hourly values of that day
#'              \item \code{"2014-04"} Will generate all the hourly values of that month
#'              \item \code{"2014"} Will generate all the hourly values of that year
#'              \item \code{c("2013-06-15","2014-10-02")} Will generate all the hourly values for that custom range. Both elements must be YYYY-MM-DD.
#'              \item \code{"all"} Will generate all the hourly values from 1979-01-01 to the last day of the latest available month.
#'              \item \bold{Note:} To check the latest available date click \href{http://goldsmr2.sci.gsfc.nasa.gov/opendap/hyrax/MERRA/contents.html}{here} and navigate through the product.
#'              }
#' @param       Coordinates Numeric vector. Can be a single point or a box, which will download the closest node or the interior nodes rescpectively. See the following examples:\itemize{
#'              \item \code{c(42.394, -76.668)} (Latitude, Longitude) Will download the data of the closest node
#'              \item \code{c(42.394, -76.668, 40.515, -74.085)} \bold{Must} have this format: First pair is the \bold{top left} corner of the box. Second pair is the \bold{bottom right} corner of the box.
#'              }
#' @param       ExtraVars Default TRUE. In addition to the variables requested above, will download Merra's \code{'XDim','YDim','TIME'} (longitude, latitude and time respectively)
#' @param       MERRA2 latest version of Merra
#' @return      Character vector containing the urls.
#'
#' @references  A more comprehensive description of MERRA variables can be found \href{http://gmao.gsfc.nasa.gov/products/documents/GEOS-5_Filespec_Glossary.pdf}{here}.
#'
#' @examples \dontrun{
#' # Example:
#' urls <- generateMerraLinks(Product = 'tavg1_2d_slv_Nx', Variables = c('U2M','V2M'),
#'                            Time = c('2013-11-01','2013-12-31'),
#'                            Coordinates = c(42.394, -76.668, 40.515, -74.085),
#'                            ExtraVars = T, MERRA2 = T)
#' }
#'
#' @author      Christian Rivero \email{crivero@@resurety.com}
generateMerraLinks <- function(Product = NA, Variables = NA, Time = NA, Coordinates = NA, ExtraVars = TRUE, MERRA2 = FALSE){

  lookup <- data.frame(Products     = c('inst1_2d_int_Nx','tavg1_2d_flx_Nx','tavg1_2d_int_Nx','tavg1_2d_lnd_Nx','tavg1_2d_ocn_Nx','tavg1_2d_rad_Nx','tavg1_2d_slv_Nx'),
                       ShortName_M1 = c('MAI1NXINT',      'MAT1NXFLX',      'MAT1NXINT',      'MAT1NXLND',      'MAT1NXOCN',      'MAT1NXRAD',      'MAT1NXSLV'),
                       ShortName_M2 = c('M2I1NXINT',      'M2T1NXFLX',      'M2T1NXINT',      'M2T1NXLND',      'M2T1NXOCN',      'M2T1NXRAD',      'M2T1NXSLV'),
                       stringsAsFactors=F)
  if(MERRA2){
    rootUrl <- 'https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/hyrax/MERRA2'
    shName <- lookup$ShortName_M2[lookup$Products == Product]
    str <- '^MERRA2.*nc4.html$'
    str2 <- '.nc' # could also be .nc4 which has a slightly better compression than netCDF V3. To use netCDF V4, NCO tool should be compiled with --enable-netcdf4
    lon_str <- ',lon'
    lat_str <- ',lat'
    time_str <- ',time[0:1:23]'

  } else {
    rootUrl <- 'https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/hyrax/MERRA'
    shName <- lookup$ShortName_M1[lookup$Products == Product]
    str <- '^MERRA.*hdf.html$'
    str2 <- '.nc'
    lon_str <- ',XDim'
    lat_str <- ',YDim'
    time_str <- ',TIME[0:1:23]'
  }

  # 1. Generate sequence of dates
  if(length(Time)==2){
    seque <- seq(as.Date(Time[1]), as.Date(Time[2]), "day")
    seque <- gsub('-', '', seque) # removes dashes
  } else if(length(Time)==1){
    if(nchar(Time)==4){ # entire year
      seqFrom <- paste0(Time,'-01-01')
      seqTo   <- paste0(Time,'-12-31')
      seque <- seq(as.Date(seqFrom), as.Date(seqTo), "day")
      seque <- gsub('-', '', seque) # removes dashes
    } else if(nchar(Time)==7){ # entire month
      seqFrom <- paste0(Time,'-01')
      nmbd    <- numberOfDays(Time)
      seqTo   <- paste0(Time,'-',nmbd)
      seque <- seq(as.Date(seqFrom), as.Date(seqTo), "day")
      seque <- gsub('-', '', seque) # removes dashes
    } else if(nchar(Time)==10){ # single day
      seque <- gsub('-', '', as.Date(Time))
    } else if(Time == "all"){
      seqFrom <- '1979-01-01'
      seqTo <- lastWebMerraMonth()
      seqTo <- paste(seqTo, numberOfDays(seqTo), sep = "-")
      seque <- seq(as.Date(seqFrom), as.Date(seqTo), "day")
      seque <- gsub('-', '', seque) # removes dashes
    } else return(stop('Parameter "Time" has a bad format'))
  } else return (stop('More than 2 parameters is not allowed in parameter "Time"'))

  # 2. Construct the body (left part) of the url by using: Product, Date sequence and crawling OPeNDAP
  ## First: Access to the Product Section
  contentsUrl <- file.path(rootUrl, 'contents.html')
  rawXml <- RCurl::getURL(contentsUrl)
  rootDoc   <- XML::htmlParse(rawXml)
  rootLinks <- XML::xpathSApply(rootDoc, "//a/@href") # Selects all attributes that are named href
  XML::free(rootDoc) # errase it

  prodLink  <- rootLinks[grepl(paste0("^",shName), rootLinks)] # Find string that starts with the corresponding ShortName
  prodUrl   <- file.path(rootUrl, prodLink)
  ## Use date sequence to get the daily URLs. Crawl by YEAR/MONTH pair values.
  unqYM <- unique(substr(seque,1,6))
  prodUrl2 <- gsub('/contents.html','',prodUrl) # errasing the tail
  unqYMurl <- file.path(prodUrl2, substr(unqYM,1,4),substr(unqYM,5,6),'contents.html') # YEAR/MONTH URLs to be crawled
  ## Crawl and extract all daily urls in every month. This is done because the header of the daily urls can start with either MERRA100, MERRA200, MERRA300 or MERRA301 (3 consecutive months that varies depending on the product)
  bodylefturl <- sapply(unqYMurl, function(x){
    rawXml <- RCurl::getURL(x)
    rootDoc2   <- XML::htmlParse(rawXml)
    rootLinks2 <- XML::xpathSApply(rootDoc2, "//a/@href") # Selects all attributes that are named href
    XML::free(rootDoc2) # errase it
    tempLinks   <- rootLinks2[grepl(str, rootLinks2)]# Find daily files that must start with 'MERRA' and end with 'hdf.html'
    tempLinks   <- unique(tempLinks) # *.hdf.html is found twice
    fullLinks   <- file.path(gsub('/contents.html','',x), tempLinks)
    return(fullLinks)
  }, USE.NAMES = FALSE)
  # Checkpoint 1: write.csv(bodylefturl, '/bodylefturl.csv')

  # Crop the list only if Time is of lenght = 2, which means custom length
  #if(length(Time)==2){
    # TODO crop bodylefturl which contains full months of data
  #}
  # remove .html to paste the right part of the url
  bodylefturl2 <- gsub('.html', str2, unlist(bodylefturl))

  # 3. Construct the right side of the URL by using: Variables and Coordinates
  # Find Merra's OPeNDAP indices
  if(length(Coordinates)==2){ # Retrieve Closest Node
    ijmn <- latlonToMerraij(inpoint=Coordinates, MERRA2 = MERRA2)
    # Output strings
    #YDim <- paste0(ydim,':1:',ydim)
    #XDim <- paste0(xdim,':1:',xdim)
    #YDim <- paste0('[',YDim,']') # Adding brackets
    #XDim <- paste0('[',XDim,']')
  } else if(length(Coordinates)==4){ # Retrieve Interior Node
    topl_ij <- latlonToMerraij(inpoint=Coordinates[c(1,2)], MERRA2 = MERRA2)
    botr_ij <- latlonToMerraij(inpoint=Coordinates[c(3,4)], MERRA2 = MERRA2)
    topbot  <- rbind(topl_ij, botr_ij)
    i_ordrd <- topbot$imn[order(topbot$imn)]
    j_ordrd <- topbot$jmn[order(topbot$jmn)]

    if(FALSE){ # debug plot
      tmp <- expand.grid(j = seq(j_ordrd[1],j_ordrd[2]),
                         i = seq(i_ordrd[1],i_ordrd[2]))

      tmp2 <- apply(tmp, 1, function(x) gridMerra2(x[1],x[2]))

      ylim <- range(c(sapply(tmp2, '[[',1), Coordinates[c(1,3)]))
      xlim <- range(c(sapply(tmp2, '[[',2), Coordinates[c(2,4)]))

      plot(Coordinates[c(2,4)],Coordinates[c(1,3)],xlim=xlim, ylim=ylim, pch=3,col=2)
      sapply(tmp2,function(x) points(x[2],x[1]))
    }

    # Out strings
    YDim <- paste(j_ordrd, collapse=':1:')
    XDim <- paste(i_ordrd, collapse=':1:')
    YDim <- paste0('[',YDim,']') # Adding brackets
    XDim <- paste0('[',XDim,']')
  } else return(stop('Parameter "Coordinates" has a bad format'))

  # 4. Paste the Variables
  # The OPeNDAP syntax is: url?VAR1[TIME][YDim][XDim],VAR2[TIME][YDim][XDim],...
  # Variables XDim, YDim and TIME are special cases, they only take one argument
  rightside <- paste0(Variables,'[0:1:23]', YDim, XDim, collapse=',')
  if(ExtraVars){
    rightside <- paste0(rightside, lon_str, XDim, lat_str, YDim, time_str)
  }

  outUrls <- paste0(bodylefturl2, '?', rightside)

  return(outUrls)
}

#' @title       Convert Lat/lon to Merra indices i/j
#' @description Make the conversition by finding the closest node
#' @param       inpoint Input numeric vector. Eg c(42.394, -76.668)
#' @param       MERRA2 for MERRA2
#' @return      data.frame
latlonToMerraij <- function(inpoint=NA, MERRA2 = FALSE){
  if(MERRA2){
    gr <- gridMerra2(rangeOnly = T)
    lons <- gr$lon
    lats <- gr$lat
  } else {
    IMn <- 1:540 # Longitude indices
    JMn <- 1:361 # Latitude indices   [pg 7](http://gmao.gsfc.nasa.gov/products/documents/MERRA_File_Specification.pdf)

    lons <- -180 + (2/3)*(IMn - 1)  # [pg 7] idem
    lats <- -90  + (1/2)*(JMn - 1)
  }

  lowestLatidx <- findInterval(inpoint[1], lats) # Finds the index of the closest lat. Always gets the lowest lat.
  lowestLonidx <- findInterval(inpoint[2], lons)
  # The closest nodes is found by calculating the distance to each of the 16 neighbors
  rangeLatidx <- seq(lowestLatidx-1, lowestLatidx+2)
  rangeLonidx <- seq(lowestLonidx-1, lowestLonidx+2)

  rangeLats <- lats[(lowestLatidx-1):(lowestLatidx+2)] # 4 points
  rangeLons <- lons[(lowestLonidx-1):(lowestLonidx+2)] # 4 points
  # TODO if any value in rangelats/lons is equal to NA or numeric(0) STOP! because the requested coordinates are out of the planet
  allCombs <- expand.grid(lon = rangeLons, lat = rangeLats) # Find all combinations of the 16 points
  allCombs$distVect <- raster::pointDistance(p1=allCombs, p2=c(inpoint[2], inpoint[1]), lonlat = T)
  minIdx   <- which.min(allCombs$distVect)
  minLon <- allCombs$lon[minIdx]
  minLat <- allCombs$lat[minIdx]
  # Find minLon/Lat corresponding indeces in JMn, IMn
  imn <- which(lons==minLon) - 1 # -1 to make it OPeNDAP compatible
  jmn <- which(lats==minLat) - 1

  outdf <- data.frame(imn, jmn)
  return(outdf)
}

#' @title       Assign Record Dimension
#' @description Assign TIME as record dimension to Merra .nc files
#' @param       ncfileIn  Input .nc file
#' @param       ncfileOut Output .nc file
#' @param       recDim    default is "TIME"
#' @return      char string with output file name
assignRecordDim <- function(ncfileIn = NA, ncfileOut = NA, recDim = "TIME"){
  if(is.na(ncfileOut)) ncfileOut <- ncfileIn
  status <- system2(command = "ncks",
                    args    = paste("-O --mk_rec_dmn", recDim, ncfileIn, ncfileOut))
  attr(ncfileOut, "cmd_status") <- status
  return(ncfileOut)
}

#' @title       Stack NetCDF files
#' @description Stack multiple .nc files through tis record dimension to a single file
#' @param       ncfiles   char vector with .nc files to be stacked
#' @param       ncfileOut filename of final stacked .nc file
#' @return      char string with output file name
stackNcFiles <- function(ncfiles = NA, ncfileOut = NA){
  status <- system2(command = "ncrcat",
                    args    = paste(paste(ncfiles, collapse = " "), "-o", ncfileOut))
  attr(ncfileOut, "cmd_status") <- status
  return(ncfileOut)
}

#' @title       Last Web Merra Month
#' @description Find most recent Merra month for all their products
#' @return      char string with the latest month and attribute 'more' with detail info
#'              return NA means that MERRA is still updating all the products, which doesn't happen at the same time
lastWebMerraMonth <- function(){
  summdf <- data.frame(
    servers = c('https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/hyrax/MERRA2/M2T1NXSLV.5.12.4/catalog.xml',
                'https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/hyrax/MERRA2/M2T1NXFLX.5.12.4/catalog.xml'),
    LastMonth        = NA,
    filesLastMonth   = NA,
    stringsAsFactors = FALSE
  )

  for(idx in 1:nrow(summdf)){
    url        <- summdf$servers[idx]

    rawXml     <- RCurl::getURL(url)
    tableTop   <- XML::xmlRoot(XML::xmlTreeParse(rawXml))[['dataset']]
    lastYear   <- XML::xmlGetAttr(tableTop[length(tableTop)-2]$catalogRef, 'name')
    urlYear    <- gsub('catalog.xml', paste0(lastYear,'/catalog.xml'), url)

    rawXmlYear <- RCurl::getURL(urlYear)
    tableTop   <- XML::xmlRoot(XML::xmlTreeParse(rawXmlYear))[['dataset']]
    lastMonth  <- XML::xmlGetAttr(tableTop[length(tableTop)]$catalogRef, 'name')
    urlMonth   <- gsub(lastYear, paste0(lastYear,'/',lastMonth), urlYear)

    rawXmlMonth<- RCurl::getURL(urlMonth)
    monthTop   <- XML::xmlRoot(XML::xmlTreeParse(rawXmlMonth))[['dataset']]
    dailyFiles <- length(monthTop)/2 # it counts both the nc4 and nc4.xml files
    summdf$LastMonth[idx]      <- paste0(lastYear,'-',lastMonth)
    summdf$filesLastMonth[idx] <- dailyFiles
  }

  latestMonth <- unique(summdf$LastMonth)
  if(length(latestMonth) > 1) latestMonth <- NA # Merra products are not updated all the same day
  attr(latestMonth, "more") <- summdf

  return(latestMonth)
}

#' @title       Last Local Merra Month
#' @description Find the latest Merra available month in our server by picking 4 file randomly per folder
#' @param       files2sample Amount of files to check their last reading
#' @return      char string with the most recent month & meta 'more": data.frame containing: server folder | last month | files per folder
#'              returns NA when at least one Merra variable in our server has files with different time lengths.
lastLocalMerraMonth <- function(files2sample = 10){
  merra2Dir  <- buildData('MERRAV2')
  merra2Grid <- gridMerra2(rangeOnly = TRUE)
  production_vars  <- c('U2M','U10M','U50M','T2M')
  production_areas <- c('USbig', 'AUSbig'); names(production_areas) <- production_areas

  randLatLon <- lapply(production_areas, function(x){
    area_coords <- getMerraBox()[[x]]

    area_lats <- data.table::between(x = merra2Grid$lat, lower = area_coords[3], upper = area_coords[1])
    area_lats <- merra2Grid$lat[area_lats]

    area_lons <- data.table::between(x = merra2Grid$lon, lower = area_coords[2], upper = area_coords[4])
    area_lons <- merra2Grid$lon[area_lons]

    rand_lats <- sample(x = area_lats, size = files2sample)
    rand_lons <- sample(x = area_lons, size = files2sample)

    rand_coords <- paste0(rand_lats, '_', rand_lons, '.gz')
    rand_paths  <- expand.grid(x, production_vars, rand_coords)
    return(rand_paths)
  })

  files2check <- do.call(rbind, randLatLon)
  row.names(files2check) <- NULL
  files2check$fpath <- file.path(merra2Dir, files2check$Var2, files2check$Var3)

  files2check$LastMonth <- sapply(files2check$fpath, function(x){
    raw <- readRDS(x)
    last_month <- raw$Time[nrow(raw)]
    last_month <- substr(last_month, 1, 7)
    return(last_month)
  })

  files2check$TotalHours <- sapply(files2check$fpath, function(x){
    raw <- readRDS(x)
    return(nrow(raw))
  })

  latestMonth <- unique(files2check$LastMonth)
  if(length(latestMonth) > 1 || all(is.na(files2check$LastMonth))){
    cat('Following variables have different last timestamp: ', fill = T)
    print(files2check)
    latestMonth <- NA
  }

  totalHours <- unique(files2check$TotalHours)
  if(length(totalHours) > 1 || all(is.na(files2check$TotalHours))){
    cat('Following variables have more hours than expexted: ', fill = T)
    print(files2check)
    latestMonth <- NA
  }

  attr(latestMonth, "more") <- files2check
  return(latestMonth)
}

#' @title        Merra's update function
#' @description  Run the steps to update Merra: Find latest month, download files, assign record dimension, stack, drill and append
#' @param        web What month to update? Default NA will find the latest Merra month available
#' @return       list containing all the files were updated (drilled)
updateMerra <- function(web = NA){
  merra2Dir <- buildData('MERRAV2')
  if(is.na(web)){
    # Determine the next month to update
    web <- lastWebMerraMonth()
    if(is.na(web)) return("Merra is updating their products. Delaying update until they finish.")
    srv <- lastLocalMerraMonth()
    if(is.na(srv)) stop("At least one Merra variable in our server has files with different time lengths. Run lastLocalMerraMonth() for details.")
    srv <- substr(srv, 1, 7)
    if(web == srv) return("Merra in server has the latest data. No need to sync.")
    # Calculating day difference between months
    dif <- difftime(as.Date(paste0(web, "-15")),
                    as.Date(paste0(srv, "-15")),
                    units = 'days')
    if(dif < 28 | dif > 31) stop("Error. Last Merra month local and web are not 1 month apart. Check lastWebMerraMonth() and lastLocalMerraMonth() for details.")
  }

  thisMonth <- web

  # Define Merra variables to download
  slv_prod  <- "tavg1_2d_slv_Nx"
  slv_vars  <- c("U10M", "V10M", "U50M", "V50M", "T2M","U2M","V2M")
  flx_prod  <- "tavg1_2d_flx_Nx"
  flx_vars  <- c("PBLH")
  area      <- getMerraBox()$USbig
  # Generate links to download
  l1 <- generateMerraLinks(Product = slv_prod,  Variables   = slv_vars,
                           Time    = thisMonth, Coordinates = area, MERRA2 = TRUE)
  l2 <- generateMerraLinks(Product = flx_prod,  Variables   = flx_vars,
                           Time    = thisMonth, Coordinates = area, MERRA2 = TRUE)
  # Download files & generate filenames to write
  landingDir <- file.path(merra2Dir, "landing")
  if(!dir.exists(landingDir)) tmp <- dir.create(landingDir)

  # Issue R16-1218
  myopts <- RCurl::curlOptions(netrc          = TRUE, # uses netrc
                               netrc.file     = file.path(getwd(), "inst", "netrc"), # define netrc file path
                               cookiejar      = 'cookies.txt', # writes in cookie jar
                               cookie         = 'cookies.txt', # uses cookies
                               followlocation = TRUE) # follow locations. Without cookies curl hops infinitely


  basename_slv <- paste(paste(area, collapse = "_"), slv_prod, paste(slv_vars, collapse = "_"), thisMonth, sep = "_")
  fnames_slv   <- file.path(landingDir, paste0(paste(basename_slv,  seq(1, length(l1)), sep = "_"), ".nc"))
  d1           <- mapply(downloadSafe, l1, fnames_slv, MoreArgs = list(.opts = myopts))

  basename_flx <- paste(paste(area, collapse = "_"), flx_prod, paste(flx_vars, collapse = "_"), thisMonth, sep = "_")
  fnames_flx   <- file.path(landingDir, paste0(paste(basename_flx,  seq(1, length(l2)), sep = "_"), ".nc"))
  d2           <- mapply(downloadSafe, l2, fnames_flx, MoreArgs = list(.opts = myopts))
  # Check consistency of scale, offset factors
  is_consistent_d1 <- check_MERRAV2_ScaleOffset(ncfiles = d1, vars = slv_vars)
  if(!is_consistent_d1) {stop('SLV daily files have different set of scale and offset factors')}
  is_consistent_d2 <- check_MERRAV2_ScaleOffset(ncfiles = d2, vars = flx_vars)
  if(!is_consistent_d2) {stop('FLX daily files have different set of scale and offset factors')}
  # Assign Record dimension
  recdimSlv <- sapply(d1, assignRecordDim, recDim = "time")
  recdimFlx <- sapply(d2, assignRecordDim, recDim = "time")
  # Stack files
  stackSlv <- stackNcFiles(ncfiles = recdimSlv, ncfileOut = file.path(landingDir, paste0(basename_slv,".nc")))
  stackFlx <- stackNcFiles(ncfiles = recdimFlx, ncfileOut = file.path(landingDir, paste0(basename_flx,".nc")))
  # Drill & Append
  drillSlv <- drillNcVar(ncfile = stackSlv, append2existing = TRUE, CORES = 6)
  drillFlx <- drillNcVar(ncfile = stackFlx, append2existing = TRUE, CORES = 6)
  # Writing Power Shear Parameters
  #cat('Re calculating new alpha parameters at',date(), fill = TRUE)
  #alphas <- writeAlpha()
  msg1 <- paste('Success!', 'area', 'USbig', 'month', thisMonth, 'was updated:', length(unlist(drillSlv)), 'time series')
  msg2 <- paste('Success!', 'area', 'USbig', 'month', thisMonth, 'was updated:', length(unlist(drillFlx)), 'time series')

  msg3 <- updateMerraMonthly_generic('AUSbig', thisMonth)

  allmsg <- c(msg1, msg2, msg3)
  for(i in allmsg){message(i, appendLF = T)}

  return(TRUE)
}

#' @title        Merra's bounding box
#' @description  Return the geo regions that are used to download Merra
#'               Northwest point, Southeast point
#' @return       data.frame
getMerraBox <- function(){

  merraBox <- list(USbig  = c(72, -172, 15, -60), # For this area, roughly 60 MB per variable per month
                   AUStst = c(-35.313862,141.127742,-39.836308,146.882650), # DEV-2507
                   AUSbig = c(-10.5, 113, -44, 156) #DEV-4352
                   )

  return(merraBox)
}

#' @title        Experimental
#' @description  Merra's first run function
#' @param        tempIn var to control for loop
#' @return       logical TRUE if finish
firstRunMerra <- function(tempIn=1){

  merraDir <- NULL # Appease R CMD CHECK

  time <- "all"
  box  <- getMerraBox()$USbig

  vars <- list()
  vars[[1]] <- list(prod = "tavg1_2d_slv_Nx",
                    vars = c("U500", "V500", "TS"))
  vars[[2]] <- list(prod = "tavg1_2d_flx_Nx",
                    vars = c("PBLH"))

  links2down <- readRDS('/media/WD1/temp/links2down.rds')
  #links2down <- lapply(vars, function(x) generateMerraLinks(Product     = x$prod,
  #                                                          Variables   = x$vars,
  #                                                          Time        = time,
  #                                                          Coordinates = box))
  # saveRDS(links2down, '/media/WD1/temp/links2down.rds')
  # The first run should be done yearly beacuse of memory contraints
  for(j in tempIn:length(links2down)){
    urls  <- links2down[[j]]
    links <- data.frame(url  = urls,
                        year = sapply(strsplit(urls, "/"), "[[", 8),
                        stringsAsFactors = FALSE)
    years <- unique(links$year)
    append <- TRUE # only for first year
    for(y in 1997:2015){
      cat("Processing year", y, fill = T)
      # Subset links
      ylinks <- links$url[links$year == y]
      # Create .nc file names
      basename <- paste(paste(box, collapse = "_"), vars[[j]]$prod, paste(vars[[j]]$vars, collapse = "_"), y, sep = "_")
      fnames   <- file.path(merraDir, "landing", paste0(paste(basename,  seq(1, length(ylinks)), sep = "_"), ".nc"))
      # Downlad files
      d1 <- mapply(downloadSafe, ylinks, fnames)
      # Assign Record dimension
      recdim <- sapply(d1, assignRecordDim)
      # Stack in annual files
      stack <- stackNcFiles(ncfiles = recdim, ncfileOut = file.path(merraDir, "landing", paste0(basename,".nc")))
      # Errasing daily files
      tmp <- file.remove(d1)
      # Drill & Append
      drill <- drillNcVar(ncfile = stack, append2existing = append)
      append <- TRUE
    }
  }
  return(TRUE)
}

#' @title        Update Merra monthly generic
#' @description  Generic function to update Merra in a given area every month
#' @param        area_name From getMerraBox() Eg. "AUSbig"
#' @param        month Eg. "2018-09", could be also a year "2017"
#' @return       Drilled time series path
updateMerraMonthly_generic <- function(area_name, month){
  # INPUTS
  # area_name = "AUSbig"
  # month = '2018-09'
  #
  mb <- getMerraBox()
  area <- mb[[area_name]]
  if(is.null(area)) {stop('area_name ', area_name, ' cannot be found in getMerraBox()')}

  raw_dir <- buildData('MERRAV2', 'landing', area_name, 'raw')
  links_dir <- buildData('MERRAV2', 'landing', area_name, 'links')

  if(!dir.exists(raw_dir)) {dir.create(raw_dir, recursive = T)}
  if(!dir.exists(links_dir)) {dir.create(links_dir, recursive = T)}

  slv_prod  <- "tavg1_2d_slv_Nx"
  slv_vars  <- c("U10M", "V10M", "U50M", "V50M", "T2M","U2M","V2M")

  ml_all <- generateMerraLinks(Product     = slv_prod,
                               Variables   = slv_vars,
                               Time        = month,
                               Coordinates = area,
                               MERRA2      = TRUE)

  netrc_file <- system.file('netrc', package = 'REdata')
  myopts <- RCurl::curlOptions(netrc          = TRUE, # uses netrc
                               netrc.file     = netrc_file, # define netrc file path
                               cookiejar      = 'cookies.txt', # writes in cookie jar
                               cookie         = 'cookies.txt', # uses cookies
                               followlocation = TRUE) # follow locations. Without cookies curl hops infinitely

  basename_slv <- paste(paste0('area',area_name), slv_prod, paste(slv_vars, collapse = "_"), month, sep = "_")
  fnames_slv   <- file.path(raw_dir, paste0(paste(basename_slv,  seq(1, length(ml_all)), sep = "_"), ".nc"))
  d1 <- mapply(downloadSafe, ml_all, fnames_slv, MoreArgs = list(.opts = myopts))

  is_consistent_d1 <- check_MERRAV2_ScaleOffset(ncfiles = d1, vars = slv_vars)
  stopifnot(is_consistent_d1)

  recdimSlv <- sapply(d1, assignRecordDim, recDim = "time")

  stackSlv <- stackNcFiles(ncfiles   = recdimSlv,
                           ncfileOut = file.path(raw_dir, paste0(basename_slv,".nc")))

  drillSlv <- drillNcVar(ncfile = stackSlv, append2existing = TRUE, CORES = 6)
  msg <- paste('Success!', 'area', area_name, 'month', month, 'was updated:', length(unlist(drillSlv)), 'time series')
  return(msg)
}

#' @title        First Run Generic
#' @description  First Run to get an arbritary geographical area
#' @param        area_name From getMerraBox()
#' @return       Drilled time series path
firstRunMerra_generic <- function(area_name){
  #area <- getMerraBox()$AUS_SE
  #area <- getMerraBox()$AUS_QLD
  area <- getMerraBox()[[area_name]]
  thisMonth <- '1980-01_2018-08'

  rfl_idx <- landingDir <- NULL # Appease R CMD CHECK

  raw_dir <- buildData('MERRAV2', 'landing', area_name, 'raw')
  links_dir <- buildData('MERRAV2', 'landing', area_name, 'links')

  if(!dir.exists(raw_dir)) {dir.create(raw_dir, recursive = T)}
  if(!dir.exists(links_dir)) {dir.create(links_dir, recursive = T)}

  slv_prod  <- "tavg1_2d_slv_Nx"
  slv_vars  <- c("U10M", "V10M", "U50M", "V50M", "T2M","U2M","V2M")

  ml1 <- generateMerraLinks(Product     = slv_prod,
                            Variables   = slv_vars,
                            Time        = c('1980-01-01','1999-12-31'),
                            Coordinates = area,
                            MERRA2      = TRUE)
  ml2 <- generateMerraLinks(Product     = slv_prod,
                            Variables   = slv_vars,
                            Time        = c('2000-01-01','2010-12-31'),
                            Coordinates = area,
                            MERRA2      = TRUE)
  ml3 <- generateMerraLinks(Product     = slv_prod,
                            Variables   = slv_vars,
                            Time        = c('2011-01-01','2018-08-31'),
                            Coordinates = area,
                            MERRA2      = TRUE)

  ml_all <- c(ml1, ml2, ml3)
  #ml_all <- ml11
  saveRDS(ml_all, file.path(links_dir, 'merra_links.rds'))

  netrc_file <- system.file('netrc', package = 'REdata')
  myopts <- RCurl::curlOptions(netrc          = TRUE, # uses netrc
                               netrc.file     = netrc_file, # define netrc file path
                               cookiejar      = 'cookies.txt', # writes in cookie jar
                               cookie         = 'cookies.txt', # uses cookies
                               followlocation = TRUE) # follow locations. Without cookies curl hops infinitely


  basename_slv <- paste(paste0('area',area_name), slv_prod, paste(slv_vars, collapse = "_"), thisMonth, sep = "_")
  fnames_slv   <- file.path(raw_dir, paste0(paste(basename_slv,  seq(1, length(ml_all)), sep = "_"), ".nc"))
  saveRDS(fnames_slv, file.path(links_dir, 'fnames_slv.rds'))

  #NOTE: if download will take a lot of time see the script inst/Rscript/merra_pull.R which is useful for running via Rscript in the background
  d1           <- mapply(downloadSafe, ml_all, fnames_slv, MoreArgs = list(.opts = myopts))
  #saveRDS(fnames_slv,'/media/NAS1/MERRAV2/landing/AUS_SE/links/fnames_slv.rds')
  # ad-hoc check if all files were donwloaded
  #rfl <- list.files(landingDir, pattern = '*.nc')
  #rfl_idx <- sapply(strsplit(rfl, '_'), '[[', 16)
  #rfl_idx <- as.numeric(gsub('.nc', '', rfl_idx))
  #rfl <- data.frame(rfl, rfl_idx, stringsAsFactors = FALSE)
  #rfl <- rfl[order(rfl$rfl_idx),]
  #which(diff(rfl$rfl_idx) != 1)

  d1 <- list.files(raw_dir, pattern = paste0(thisMonth, '.*.nc'), full.names = TRUE)
  is_consistent_d1 <- check_MERRAV2_ScaleOffset(ncfiles = d1, vars = slv_vars)
  stopifnot(is_consistent_d1)

  recdimSlv <- sapply(d1, assignRecordDim, recDim = "time")

  # below didn't work: had to use command line
  # ls | sort -V | ncrcat -o areaAUS_SE_tavg1_2d_slv_Nx_U10M_V10M_U50M_V50M_T2M_U850_V850_1981-01_2017-08.nc
  # ls areaAUS_SE_tavg1_2d_slv_Nx_U10M_V10M_U50M_V50M_T2M_U850_V850_2017-09_2017-10* | sort -V | ncrcat -o areaAUS_SE_tavg1_2d_slv_Nx_U10M_V10M_U50M_V50M_T2M_U850_V850_2017-09_2017-10.nc
  # ls areaAUS_SE_tavg1_2d_slv_Nx_U10M_V10M_U50M_V50M_T2M_U850_V850_2017-11_2017-12* | sort -V | ncrcat -o areaAUS_SE_tavg1_2d_slv_Nx_U10M_V10M_U50M_V50M_T2M_U850_V850_2017-11_2017-12.nc
  # ls areaAUS_SE_tavg1_2d_slv_Nx_U10M_V10M_U50M_V50M_T2M_U850_V850_2018-03* | sort -V | ncrcat -o areaAUS_SE_tavg1_2d_slv_Nx_U10M_V10M_U50M_V50M_T2M_U850_V850_2018-03.nc

  stackSlv <- stackNcFiles(ncfiles   = rfl_idx$recdimSlv,
                           ncfileOut = file.path(landingDir,"areaAUS_SE_tavg1_2d_slv_Nx_U10M_V10M_U50M_V50M_T2M_U850_V850_1981-01_2017-08.nc"))
  #stackSlv <- file.path(landingDir,"areaAUS_SE_tavg1_2d_slv_Nx_U10M_V10M_U50M_V50M_T2M_U850_V850_2017-09_2017-10.nc")
  #stackSlv <- file.path(landingDir,"areaAUS_SE_tavg1_2d_slv_Nx_U10M_V10M_U50M_V50M_T2M_U850_V850_2017-11_2017-12.nc")

  #NOTE: For one-off drills use the following script: inst/Rscript/one_time_merra.R which is useful for running via Rscript in the background
  stackSlv <- file.path(landingDir,"areaAUS_SE_tavg1_2d_slv_Nx_U10M_V10M_U50M_V50M_T2M_U850_V850_2018-06.nc")

  drillSlv <- drillNcVar(ncfile = stackSlv, append2existing = TRUE)
  return(drillSlv) # added 100 nodes per variable
}

#' @title        First Run Brasil
#' @description  First Run to get Brasilian project Merra nodes
#'               More info: https://resurety.atlassian.net/browse/DEV-2539
#' @return       Drilled time series path
firstRunMerra_BRAtst <- function(){
  area <- getMerraBox()$BRAtst
  slv_prod  <- "tavg1_2d_slv_Nx"
  slv_vars  <- c("U10M", "V10M", "U50M", "V50M", "T2M","U850","V850")

  rfl_idx <- NULL # Appease R CMD CHECK

  ml <- generateMerraLinks(Product     = slv_prod,
                           Variables   = slv_vars,
                           Time        = c('1981-01-01','1999-12-31'),
                           Coordinates = area,
                           MERRA2      = TRUE)
  ml2 <- generateMerraLinks(Product     = slv_prod,
                            Variables   = slv_vars,
                            Time        = c('2000-01-01','2010-12-31'),
                            Coordinates = area,
                            MERRA2      = TRUE)
  ml3 <- generateMerraLinks(Product     = slv_prod,
                            Variables   = slv_vars,
                            Time        = c('2011-01-01','2017-07-31'),
                            Coordinates = area,
                            MERRA2      = TRUE)
  ml_all <- c(ml, ml2, ml3)
  #saveRDS(ml_all,'/media/NAS1/MERRAV2/tmp_Brasil_DEV2539/merra_links.rds')

  myopts <- RCurl::curlOptions(netrc          = TRUE, # uses netrc
                               netrc.file     = file.path(getwd(), "inst", "netrc"), # define netrc file path
                               cookiejar      = 'cookies.txt', # writes in cookie jar
                               cookie         = 'cookies.txt', # uses cookies
                               followlocation = TRUE) # follow locations. Without cookies curl hops infinitely

  thisMonth <- '1981-01_2017-07'
  landingDir <- '/media/NAS1/MERRAV2/tmp_Brasil_DEV2539/raw'

  basename_slv <- paste('areaBRAtst', slv_prod, paste(slv_vars, collapse = "_"), thisMonth, sep = "_")
  fnames_slv   <- file.path(landingDir, paste0(paste(basename_slv,  seq(1, length(ml_all)), sep = "_"), ".nc"))
  d1           <- mapply(downloadSafe, ml_all, fnames_slv, MoreArgs = list(.opts = myopts))

  # ad-hoc check if all files were donwloaded
  #rfl <- list.files(landingDir, pattern = '*.nc')
  #rfl_idx <- sapply(strsplit(rfl, '_'), '[[', 15)
  #rfl_idx <- as.numeric(gsub('.nc', '', rfl_idx))
  #rfl <- data.frame(rfl, rfl_idx, stringsAsFactors = FALSE)
  #rfl <- rfl[order(rfl$rfl_idx),]
  #which(diff(rfl$rfl_idx) != 1)

  d1 <- list.files(landingDir, pattern = '*.nc', full.names = TRUE)
  is_consistent_d1 <- check_MERRAV2_ScaleOffset(ncfiles = d1, vars = slv_vars)
  stopifnot(is_consistent_d1)

  recdimSlv <- sapply(d1, assignRecordDim, recDim = "time")

  # below didn't work: had to use command line
  # ls | sort -V | ncrcat -o areaBRAtst_tavg1_2d_slv_Nx_U10M_V10M_U50M_V50M_T2M_U850_V850_1981-01_2017-07.nc
  stackSlv <- stackNcFiles(ncfiles   = rfl_idx$recdimSlv,
                           ncfileOut = file.path(landingDir,"areaBRAtst_tavg1_2d_slv_Nx_U10M_V10M_U50M_V50M_T2M_U850_V850_1981-01_2017-07.nc"))
  stackSlv <- file.path(landingDir,"areaBRAtst_tavg1_2d_slv_Nx_U10M_V10M_U50M_V50M_T2M_U850_V850_1981-01_2017-07.nc")

  drillSlv <- drillNcVar(ncfile = stackSlv, append2existing = FALSE)
  return(drillSlv) # added 100 nodes per variable
}

#' @title        First Run Glascow
#' @description  First Run to get UK project Merra nodes
#'               More info: https://resurety.atlassian.net/browse/DEV-3848
#' @return       Drilled time series path
firstRunMerra_Glascow <- function(){

  rfl_idx <- NULL # Appease R CMD CHECK

  area <- getMerraBox()$Glascow
  slv_prod  <- "tavg1_2d_slv_Nx"
  slv_vars  <- c("U10M", "V10M", "U50M", "V50M", "T2M","U850","V850")
  #
  # ml <- generateMerraLinks(Product     = slv_prod,
  #                          Variables   = slv_vars,
  #                          Time        = c('1981-01-01','1999-12-31'),
  #                          Coordinates = area,
  #                          MERRA2      = TRUE)
  # ml2 <- generateMerraLinks(Product     = slv_prod,
  #                           Variables   = slv_vars,
  #                           Time        = c('2000-01-01','2010-12-31'),
  #                           Coordinates = area,
  #                           MERRA2      = TRUE)
  # ml3 <- generateMerraLinks(Product     = slv_prod,
  #                           Variables   = slv_vars,
  #                           Time        = c('2011-01-01','2018-06-30'),
  #                           Coordinates = area,
  #                           MERRA2      = TRUE)
  # ml_all <- c(ml, ml2, ml3)
  #saveRDS(ml_all,'/media/NAS1/MERRAV2/tmp_Glascow_DEV3848/merra_links.rds')
  ml_all <- readRDS('/media/NAS1/MERRAV2/tmp_Glascow_DEV3848/merra_links.rds')
  myopts <- RCurl::curlOptions(netrc          = TRUE, # uses netrc
                               netrc.file     = file.path(getwd(), "inst", "netrc"), # define netrc file path
                               cookiejar      = 'cookies.txt', # writes in cookie jar
                               cookie         = 'cookies.txt', # uses cookies
                               followlocation = TRUE) # follow locations. Without cookies curl hops infinitely

  thisMonth <- '1981-01_2018-06'
  landingDir <- '/media/NAS1/MERRAV2/tmp_Glascow_DEV3848/raw'

  basename_slv <- paste('areaGlascow', slv_prod, paste(slv_vars, collapse = "_"), thisMonth, sep = "_")
  fnames_slv   <- file.path(landingDir, paste0(paste(basename_slv,  seq(1, length(ml_all)), sep = "_"), ".nc"))
  d1           <- mapply(downloadSafe, ml_all[1412:length(ml_all)], fnames_slv[1412:length(fnames_slv)], MoreArgs = list(.opts = myopts))
  return('DONE')
  # ad-hoc check if all files were donwloaded
  #rfl <- list.files(landingDir, pattern = '*.nc')
  #rfl_idx <- sapply(strsplit(rfl, '_'), '[[', 15)
  #rfl_idx <- as.numeric(gsub('.nc', '', rfl_idx))
  #rfl <- data.frame(rfl, rfl_idx, stringsAsFactors = FALSE)
  #rfl <- rfl[order(rfl$rfl_idx),]
  #which(diff(rfl$rfl_idx) != 1)

  d1 <- list.files(landingDir, pattern = '*.nc', full.names = TRUE)
  is_consistent_d1 <- check_MERRAV2_ScaleOffset(ncfiles = d1, vars = slv_vars)
  stopifnot(is_consistent_d1)

  recdimSlv <- sapply(d1, assignRecordDim, recDim = "time")

  # below didn't work: had to use command line
  # ls | sort -V | ncrcat -o areaBRAtst_tavg1_2d_slv_Nx_U10M_V10M_U50M_V50M_T2M_U850_V850_1981-01_2017-07.nc
  stackSlv <- stackNcFiles(ncfiles   = rfl_idx$recdimSlv,
                           ncfileOut = file.path(landingDir,"areaBRAtst_tavg1_2d_slv_Nx_U10M_V10M_U50M_V50M_T2M_U850_V850_1981-01_2017-07.nc"))
  stackSlv <- file.path(landingDir,"areaBRAtst_tavg1_2d_slv_Nx_U10M_V10M_U50M_V50M_T2M_U850_V850_1981-01_2017-07.nc")

  drillSlv <- drillNcVar(ncfile = stackSlv, append2existing = FALSE)
  return(drillSlv) # added 100 nodes per variable
}

#' @title        Calculate Shear Alpha
#' @description  Power law alpha parameter using wind 10, 50 and monthly displacement height
#' @param        nodeLat MERRA2 latitude
#' @param        nodeLon MERRA2 longitude
#' @return       data.frame with time and alpha parameters
shearNode <- function(nodeLat, nodeLon){
  merra2Dir <- buildData('MERRAV2')
  # Raw Merra files
  dspfile <- file.path(merra2Dir,'monthly2015_DISPH',paste0(nodeLat,'_',nodeLon,'.gz'))
  u10file <- file.path(merra2Dir,'U10M',paste0(nodeLat,'_',nodeLon,'.gz'))
  u50file <- file.path(merra2Dir,'U50M',paste0(nodeLat,'_',nodeLon,'.gz'))
  # Load files
  disph <- readRDS(dspfile)
  u10 <- readRDS(u10file)
  u50 <- readRDS(u50file)
  # Add monthly displacement height to u10
  u10$Month <- as.numeric(substr(u10$Time,6,7))
  u10$disph <- disph$DISPH[u10$Month]
  # Calculating power law shear parameter (alpha)
  alpha <- log(u50$Speed/u10$Speed)/log(50/(10+u10$disph)) # 10 meters above displacement height
  outdf <- data.frame(Time = u10$Time, Alpha = alpha, stringsAsFactors = FALSE)
  return(outdf)
}

#' @title        Writes Alpha Shear
#' @description  Calculates shear using shearNode() and saves it in .../MERRAV2/Alpha10_50_Disph
#' @return       file.path of written nodes
writeAlpha <- function(){
  merra2Dir <- buildData('MERRAV2')
  alphaDir <- file.path(merra2Dir,'Alpha10_50_Disph')
  if(!dir.exists(alphaDir)) dir.create(alphaDir)
  allNodes <- list.files(file.path(merra2Dir, 'U2M'), pattern = "*.gz")
  ov <- sapply(allNodes, function(n){
    # Striping lat/lon values
    n <- gsub(".gz", "", n)
    n <- unlist(strsplit(n, "_"))
    lat <- as.numeric(n[1])
    lon <- as.numeric(n[2])
    outfile <- file.path(alphaDir,paste0(lat,'_',lon,'.gz'))
    # Calculate alpha
    alpha <- shearNode(nodeLat = lat, nodeLon = lon)
    # Extract desired columns: time and shear
    saveRDS(alpha, outfile)
    return(outfile)
  })
  return(ov)
}

#' @title        Checks MERRAV2 Scale & Offset
#' @description  Checks that all the variables `vars` in the `ncfiles` have consistent scale_factor
#'               and add_offset parameters.
#' @param        ncfiles Char vector with full file paths to nc files
#' @param        vars Char vector with all the variables, as expressed inside the nc files, to check
#' @return       TRUE if scale_factor and add_offset are consistent in all files, FALSE otherwise
#'               Attribute 'all_info' contains each file's parameters for visual inspection
#' @note         This function was created because while testing ERA5, it was found that different
#'               annual files have different set of scale and offset values.
#'               This is not the case for MERRAV2 based on their documentation:
#'               https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf
#'               Section 2.2 Variables
check_MERRAV2_ScaleOffset <- function(ncfiles, vars = c('U10M','V10M','U50M','V50M','T2M')){

  out <- lapply(vars, function(var){
    per_var <- lapply(ncfiles, function(y){
      con1 <- ncdf4::nc_open(y, readunlim = FALSE)
      scale_factor <- con1$var[[var]]$scaleFact
      add_offset <-  con1$var[[var]]$addOffset
      ncfile <- y
      ncdf4::nc_close(con1)
      return(data.frame(ncfile=ncfile,var=var,scale_factor=scale_factor, add_offset=add_offset, stringsAsFactors = F))
    })
    per_var <- do.call(rbind, per_var)
    return(per_var)
  })
  out <- do.call(rbind, out)

  same_scale <- length(unique(out$scale_factor)) == 1
  same_offset <- length(unique(out$add_offset)) == 1

  is_scale_and_offset_ok <- same_scale & same_offset
  attr(is_scale_and_offset_ok, 'all_info') <- out

  return(is_scale_and_offset_ok)
}
