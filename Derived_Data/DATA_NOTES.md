# Data Notes For Shellfish Sanitation Data

# Initial GIS Data
## Sampling Locations (known as 'Stations" in our anlysis)
GIS data were assembled from the p90 data, which included lats and longs. 
Reviewing the p90 data files, makes it  clear the p90 data files from each year 
are not independent.  Therefore we only used them to generate GIS layers for
displaying results and assessing geospatial relationships in these data.

We merged three p90 data files (for 2016, 2017 and 2018) to generate a more 
complete list of sampling locations than could be assembled from any 
one year.

After merging three years of data, statewide data was read into ArcGIS as a 
CSV table. Casco Bay Locations were selected by locating points within 500 m of 
the Casco Bay watershed layer, and resulting data was exported as a 
shapefile (`Casco_Bay_Shellfish_p90_Locations`). 

Although in principle such simple geographic search criteria could miss sites
more than 1 KM from land, or pick up extra locations just outside the watershed, 
the process appears to have worked (based on visual inspection).  We found only a 
single Station Location not in the Growing Areas in Casco Bay.  We remove 
those points from further analysis, both in GIS, and in our R code.

The first two characters of the Station ID represent the DMR Growing Area.  We
calculated a GrowingArea Field based on that relationship.

## Growing Areas
Complementary growing area coverage was derived from data in the 
MaineDMR_Public_Health_-_Current_NSSP_Classifications shapefile.  Data was
selected by identifying polygons that intersect the Casco Bay watershed.
Most attributes were removed to simplify the coverage for our purposes.

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
       estimated land area within each circle, for a percent impervious per unit land.
       (Land area is 25 times the extracted value from the raster).  
   *   Variable names are pct_[radius] and pct_l_[radius], respectively for percent
       based on total area nad land area.  

10.  Impervious cover data was exported in a text file "station_imperviousness.csv".

