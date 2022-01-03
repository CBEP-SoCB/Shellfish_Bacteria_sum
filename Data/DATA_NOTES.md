# Data Notes For Shellfish Sanitation Data

## Sampling Locations (known as 'Stations' in our anlysis)
GIS data were assembled from publicly available DMR "P90" data, which included
latitudes and longitudes.  Those files do not contain raw observations, so we
only used them to generate location data. The online location of the "P90" files
has changed over the past few years, and older data is no longer available.

DMR expose related data 
[here](https://dmr-maine.opendata.arcgis.com/search?tags=shellfish).
We assembled a dataset based on the union of all locations from 2016, 2017 and 
2018 P90 data. We exported locations both as a CSV file and as a shapefile (in
the "GIS" folder).

After merging three years of data, statewide data was read into ArcGIS as a 
CSV table. Casco Bay Locations were selected by locating points within 500 m of 
the Casco Bay watershed layer, and resulting data was exported as a 
shapefile.

In principle such a simple geographic search criteria could miss sites
more than 1 KM from land, or pick up extra locations just outside the watershed.
We checked, and found only a single Station Location not in the Growing Areas 
in Casco Bay.  We remove that point from further analysis.

().

### `Sellfish p90 Locations.csv`

Column Name  | Contents                                                      
-------------|-----------------------------------------------
Station      | Alphanumeric Station code from DMR
Lat          | Latitude, WGS 1984, decimal degrees 
Long         | Longitude, WGS 1984, decimal degrees

The first two characters of the Station ID represent the DMR Growing Area.  We
calculated a Grow_Area Field for the shapefile based on that relationship.


### Shapefile `Casco_Bay_Shellfish_p90_Locations`
Column Name  | Contents                                                      
-------------|-----------------------------------------------
Station      | Alphanumeric Station code from DMR
Lat          | Latitude, WGS 1984, decimal degrees 
Long         | Longitude, WGS 1984, decimal degrees
Grow_Area    | Two letter code for DMR Growing Area

##  Near Impervious Cover Estimates
Impervious cover estimates (calculated only for Station locations) were
based on Maine IF&W one meter pixel impervious cover data, which is based
largely on data from 2007.  CBEP has a version of this impervious cover data for
the Casco Bay watershed towns in our GIS data archives. Analysis followed the
following steps. 

1. Town by town IC data in a Data Catalog were assembled into a large `tif` 
   file using the "Mosaic Raster Catalog"  item from the context menu from the
   ArcGIS table of contents.

2. We created a polygon that enclosed all of the Casco Bay Station locations and
   a 2000 meter buffer around them.  Because our version of the Impervious Cover
   layer is limited to Casco Bay Watershed Towns, we can not develop impervious
   cover statistics for nearby sites outside the watershed towns.

3. We used "Extract by Mask" to extract a smaller version of the impervious
   cover data for just our buffered sample region.  

4. We used "Aggregate" to reduce the resolution of the impervious cover raster
   to a 5 meter resolution, summing the total impervious cover within the
   5m x 5 m area, generating a raster with values from zero to 25. This
   speeds later processing, with a negligible reduction in precision.

5. We used "Focal Statistics" to generate rasters that show the cumulative area
   of impervious cover (in meters) within 100 m, 500 m, and 1000 m. 

6. We clipped the `cnty24p` data layer, to the mask polygon, merged all 
   polygons to a single multipolygon and added a dummy attribute with a value 
   of one.  We converted that to a raster, with a 5 meter pixel, and a value of 
   one everywhere there was land (Note that each pixel has value of one, but
   covers 5m x 5m = 25 meters square, so this needs to be taken into account
   later.  

7. We used "Focal Statistics" to generate rasters that show the cumulative sum
   (NOT area) of the land cover raster within 100 m, 500 m, and 1000 m.
   (to get true area, we still need to multiply values by 25).

8. We extracted the values of the three rasters produced in step 5 and three
   rasters produced in step 7 at each Station location. We used  'Extract 
   Multi Values to Points'. (variable names are imperv_[radius] and 
   land[_radius] respectively).  For IC, but not land cover, we replaced any 
   null values with zeros, to account for points that lie more than the specified 
   distance from impervious cover using the field calculator.

9. We calculated (two versions of) percent cover with the Field Calculator.   
   *   We divided the total impervious cover within the specified distance by the 
       area of the circle sampled under step (5) ($\pi \cdot r^2$).  
   *   We divided the total impervious cover within the specified distance by the 
       estimated land area within each circle, for a percent impervious per unit 
       land. (Land area is 25 times the extracted value from the raster).  
   *   Variable names are pct_[radius] and pct_l_[radius], respectively for percent
       based on total area and land area.  

10.  Impervious cover data was exported in a text file "station_imperviousness.csv".

### `station_imperviousness.csv`

Column Name  | Contents                                             
-------------|-----------------------------------------------  
OBJECTID     | Arbitrary ID assigned by ArcGIS  
Station      | Alphanumeric DMR Code  
Lat          | Latitude, WGS84, decimal degrees  
Long         | longitude, WGS84, decimal degrees  
Grow_Area    | Two letter abbreviation for DMR Growing Areas  
land_1000    | Land within 1000 meters of sample location  
land_500     | Land within 500 meters of sample location  
land_100     | Land within 100 meters of sample location  
imperv_1000  | Impervious area within 1000 meters of sample location  
imperv_500   | Impervious area within 500 meters of sample location  
imperv_100   | Impervious area within 100 meters of sample location  
pct_100      | Percent of 100 m circle that is impervious  
pct_500      | Percent of 500 m circle that is impervious  
pct_1000     | Percent of 1000 m circle that is impervious  
pct_l_100    | Percent of land area within 100 meters that is impervious  
pct_l_500    | Percent of land area within 500 meters that is impervious  
pct_l_1000   | Percent of land area within 1000 meters that is impervious  

# Data on Bacteria at Shellfish Sampling Locations
This is data reworked from the original DMR data.  it contains values for 
coliform abundance in individual water samples. DMR uses membrane filtration to 
count coliform colony forming units.

Column Name  | Contents                             | Units
-------------|--------------------------------------|-------------
SDate        | Sample date                          | yyyy-mm-dd
STime        | Sample time                          | HH:MM:SS Local time?
SDateTime    | Sample date and time                 | yyyy-mm-ddTHH:MM:SSZ
Station      | DMR Station Code                     | Alphanumeric
GROW_AREA    | DMR Growing area that contains the sample location        |  
Tide         | Tide stage: "L" = low;  "LF" = low flood; "F" = Flood; "HF" = high flood; "H" = high; "HE" = high ebb; "E" = ebb;   "LE" = low ebb |  
Class        | Harvesting area classification 'A' = Approved, 'CA'= Conditionally Approved, 'CR' = Conditionally Restricted, 'R' = Restricted, 'P' = Prohibited  |  
Temp         | Water temperature                    | Degrees C 
Sal          | Salinity                             | Original metadata said "pct" but values suggest PPT
ColiScore    | Coliform data as interpreted by DMR  | Coliform bacteria colonies per 100 ml
RawColi      | Raw coli data, including indicators of censoring | Same
YEAR         | Year of sample collection            | Four digit integer
LCFlag       | Flag for left censored values        | TRUE / FALSE
RCFlag       | Flag for right censored values       | TRUE/FALSE
ColiVal      | Coliform numbers again, including censoring limits for censored values. | Coliform bacteria colonies per 100 ml
 
# Weather Data
Daily weather data for the Portland Jetport was downloaded via a NOAA weather
data API using a small python script.  More details on the NOAA weather API and
on the python programs we used to access data are available at a companion
github archives on
[climate change](https://github.com/CBEP-SoCB-Details/CDO_Portland_Jetport.git).
Documentation on specific NOAA data sets that are available through their API is
available at https://www.ncdc.noaa.gov/cdo-web/datasets.  

The version of the data included here was derived from downloaded data by 
selecting data from the years 2015 through 2019, and selectively removing 
data columns of no interest here.

Here, we downloaded daily (GHCND) weather summaries via API v2. Information
on this API is available here: https://www.ncdc.noaa.gov/cdo-web/webservices/v2

Documentation on specific datasets is available at
https://www.ncdc.noaa.gov/cdo-web/

## `Portland_Jetportt_2015-2019.csv`
Column Name     | Contents                | Units                         
----------------|-------------------------|------
date   | Date of weather observations     |  mm/dd/yyyy
AWND   | Average daily wind speed         |meters per second
PRCP   | Pprecipitation, rainfall equivalents | mm
SNOW   | Snowfall                         | mm
SNWD   | Snow depth                       | mm
TAVG   | Nominal daily average temperature; average of TMIN and TMAX? | Celsius
TMAX   | Maximum daily temperature        | Celsius
TMIN   | Minimum daily temperature        | Celsius
WDF2   | W Direction of fastest 2-minute wind | degrees
WDF5   | Direction of fastest 5-second wind   | degrees
WSF2   | Fastest 2-minute wind speed      | meters per second
WSF5   | Fastest 5-second wind speed      | meters per second

## Units
Data is in SI units, except that NOAA provides some data in tenths of the
nominal units. This is not well documented through the API, but obvious in 
context. Temperatures are reported in tenths of degrees C, and precipitation in
tenths of a millimeter. For this analysis, we disregard trace rainfall.
