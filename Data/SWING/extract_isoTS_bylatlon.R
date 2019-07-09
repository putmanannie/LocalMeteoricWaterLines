# extract_isoTS_byLatLon.R
# richfiorella 181116.

# using annie's list of lat/lons, extract timeseries of monthly
# data from 21 years of climate model data. write out one file
# per model simulation.

# clean old data
rm(list=ls())

# load required packages
library(raster)
library(readxl)

# open list of lat/lons
df <- read_excel("Lat_Lon_list.xlsx",sheet="Sheet1")

# load model data.
model <- brick("~/Dropbox/L4/nudged.gissE.1980-2001.nc",varname="d18Op")

# check model coordinates, and see if lons are listed as 0 to 360 or -180 to 180
if (model@extent@xmax > 185) {
  df$Lon[df$Lon < 0 ] <- df$Lon[df$Lon < 0] + 360
}

# create spatial datapoints coordinate based on data frame.
df_points <- SpatialPoints(data.frame(df$Lon,df$Lat),proj4string = model@crs)

# set up output dataframe.
out_df <- list()

for (i in 1:model@data@nlayers) {
  out_df[[i]] <- extract(model[[i]],df_points)
} 

# okay now bind along columns. should be 398 sites (rows) by 252 columns (months)
out_df_restruct <- as.data.frame(do.call(cbind,out_df))

colnames(out_df_restruct) <- seq(as.Date("1980-01-15"),as.Date("2000-12-15"),by="months")

# okay, now cbind this to df.
data_out <- cbind(df,out_df_restruct)

# save data and variable. 
write.csv(data_out,"nudged_gissE_1x1_d18O.csv")
