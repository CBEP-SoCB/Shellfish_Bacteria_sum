# Primary Data
Curtis C. Bohlen received an e-mail From Benjamin Wahle of DMR on November 27,
2019, containing the following file: "Casco Bay WQ 15 19.csv"  This file
contains the primary data on which our analyses are based.

Curtis C. Bohlen received another e-mail in Decemberof 2019, from Wahle, 
containing useful metadata:
"Geomean and P90 Explanation.docx"

#"P90" data
After accessing several years' worth of "P90" data from DMR's public data
repositories,  we realized the data were not complete enough for our analytic 
needs, and in particular, that the mapped p90 values released each year were
based on the most recent 30 observations at each site, which extend over
multiple years.  Thus the annual "p90" data are not sequentially independent.

We ended up using the "p90" data principally for it's geographic content.
Data was accessed as follows:

## "Public" 2018 "P90" data
2018 P90 Data were downloaded from the following URL:
https://dmr-maine.opendata.arcgis.com/datasets/mainedmr-public-health-2018-p90-scores/data
as a Shapefile and a CSV file.

## Older "P90" data
P90 datasets prior to 2018 could not be accessed from that site. Years 2017 and
2016 appear here: (https://dmr-maine.opendata.arcgis.com/search?q=p90) (or could
be accessed by editing the URL for the 2018 data in the obvious way)

2017 P90 Data URL: (https://dmr-maine.opendata.arcgis.com/datasets/mainedmr-public-health-2017-p90-scores)

2016 P90 Data URL: (https://dmr-maine.opendata.arcgis.com/datasets/mainedmr-public-health-2016-p90-scores)

But clicking on the 'Download" tabs from those pages did not work.

We requested the historical P90 data from DMR via e-mail.  Curtis C. Bohlen 
received an e-mail from Benjamin Wahle of DMR on November 26, 2019, containing 
the following two files:

"2016 p90 for CBEP.xlsx"
"2017 p90 for CBEP.xlsx"

# Weather Data
CBEP uses a custom Python program to download data from NOAA's online data
repositories.  Specifically, data were accessed through NOAA's National Centers
for Environmental Information.

Here, we have downloaded daily (GHCND) weather summaries via API v2. Information
on this API is available here: https://www.ncdc.noaa.gov/cdo-web/webservices/v2

Documentation on specific datasets is available at
https://www.ncdc.noaa.gov/cdo-web/datasets

Portland Jetport weather data was downloaded by Curtis C. Bohlen using a custom
python script, titled "noaaweatherdataGUI.py" on April 23, 2020, when analyzing 
Long Creek Watershed Management District data.  The version of the data included 
here was derived from that download by selecting data from the years 2015 through
2019.

## Units
Data is in SI units, except that NOAA provides some data in tenths of the
nominal units. This is not well documented through the API, but obvious in 
context. Temperatures are reported in tenths of degrees C, and precipitation in
tenths of a millimeter. For this analysis, we disregard trace rainfall.

# Growing Areas
The shapefile `MaineDMR_Public_Health_-_Current_NSSP_Classifications-shp` was 
downloaded from DMR's Open Data Portal here:
https://dmr-maine.opendata.arcgis.com/datasets/mainedmr-public-health-current-nssp-classifications
By Curtis C. Bohlen on February 12, 2021.
