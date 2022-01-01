# Shellfish_Bacteria
Analysis of several years of data on bacteria levels in shellfish, from Maine's 
Shellfish Sanitation Program.

# Technical Background
Maine's Bureau of Public Health and the Department of Marine Resources and
manage the state's Shellfish Sanitation Program to ensure that molluscan
shellfish like clams and oysters are safe for human consumption.  The program is
a wide ranging one, including environmental surveillance for potential health
hazards, certification and inspection of shellfish dealers, and more.

Our interest here is specifically on the level of certain bacteria 
observed in Maine's coastal waters.  Data is collected as part of Maine's 
Shellfish Growing Area Classification System, which determines whether or not
specific areas are open for harvest of shellfish.  Shellfish areas can be either
permanently or episodically closed, based on assessment of risk due to
contamination, especially by pathogens or by phytotoxins (e.g., "red tide").

Shellfish harvesting areas like clam flats may be closed for harvest after
conducting a complex risk assessment process.  Generally, no single
measurement of bacteria levels is enough to close a tidal flat for shellfish
harvests, but a pattern of elevated bacteria levels would be. 

Conversely, sites may be closed based other information besides results of field
sampling. For example, most waters near Portland are permanently closed because
of proximity to wastewater discharges and combined sewer overflow discharge
points. (Paradoxically, that means DMR collects little bacteria data near the
City, where we might expect elevated levels).  A tidal flat may also be closed
after heavy rains, if it is in a location where rainfall has been observed to 
lead to elevated levels of bacteria.

Harvesting areas can independently be closed due to phytotoxins.  Closures due
to "red tide" (a bloom of the toxic dinoflagellate, *Alexandrium fundyense*) 
happen every year in the eastern part of Casco Bay, especially in May and June. 
Closures due to other phytotoxins are less common, but have been increasing in 
frequency.  We do not consider phytotoxin related closures further here, as our
interest focuses on available bacteria data.

Maine DMR collects data on concentrations of *E. coli* bacteria in coastal
waters at shellfish growing areas throughout the state.  Generally, data is
collected from each station about six times a year, and regulatory decisions
about safety of local shellfish harvests are based on evaluating the most recent 
thirty samples (thus approximately five years of data). 

Standards for whether a Station is open for harvest, conditionally open, etc.
are based in part on whether local bacteria levels meet certain thresholds.
Because bacteria data is highly skewed (many samples detect little or no
bacteria, but a few detect extremely high levels), these threshold are based on
either a long-term geometric mean of observations, or the 90th percentile of
observations.

## Growing Area Classification Standards
Growing Area Classification | Activity Allowed              |	Geometric mean FC/100ml	| 90th Percentile (P90) FC/100ml
----------------------------|-------------------------------|-------------------------|-------------------------------
Approved	               | Harvesting allowed	                  | ≤ 14	              | ≤ 31
Conditionally Approved	 | Harvesting allowed except during specified conditions | ≤ 14 in open status	| ≤ 31 in open status
Restricted	             | Depuration harvesting or relay only	| ≤ 88 and >15	      | ≤ 163 and > 31
Conditionally Restricted | Depuration harvesting or relay allowed except during specified conditions	| ≤ 88 in open status	| ≤ 163 in open status
Prohibited	             | Aquaculture seed production only	    | >88	                |>163


## *E. coli*, Enterococci, and Health Risk
while the Maine Beaches Program measures "enterococci" bacteria, the
shellfish program measures "*E. coli*"" bacteria.  These two measures address
similar problems, from slightly different perspectives. The different methods 
can be confusing, but mostly reflect different scientific and regulatory 
histories of concern about health risks due to consumption of seafood versus
water contact recreation.  The tests have a similar purpose.

Both *E. coli* and enterococci tests quantify groups of bacteria (the  tests
measure some of the same bacteria) that are common in vertebrate digestive
systems, especially the digestive systems of warm blooded organisms. These
bacteria are less common in aquatic environments. Thus elevated levels of these
bacteria (as documented by either test) indicate a higher likelihood of recent
contamination of surface waters by fecal material. 

Many human pathogens are transmitted efficiently through fecal contamination.
The bacteria being tested for under each test are not the primary pathogen of
concern (although some bacteria detected by the tests can be pathogenic).  
Instead, test results provide an indicator of risk of presence of pathogens more 
generally. While that is a concern for swimmers, bivalve shellfish, as filter 
feeders, can concentrate some pathogens in their bodies, posing even higher
risk to anyone consuming them.

Unfortunately, because of different methods, the *E coli* and enterococci tests 
are not directly comparable.  However, the  two measures are usually highly 
correlated. 

=======

Analysis of several years of data on bacteria levels in shellfish waters, from
Maine's Shellfish Sanitation Program

# Introduction
This archive includes R scripts and related products analyzing levels of
bacteria in waters where shellfish are harvested in Casco Bay, Maine.  The data
was provided  by the Department of Marine Resources, which is responsible for
tracking potential hazards to human health via consumption of shellfish.

Marine bivalve mollusks like clams, mussels, and oysters are filter feeders.
They ingest plankton and bacterioplankton from their environment, and under
certain circumstances can concentrate toxins or pathogens from the marine
environment.  Bivalves with elevated levels of toxins or pathogens can pose a
serious risk to human health.

The Bureau of Public Health oversees the application of the National Shellfish
Sanitation Program (NSSP) in Maine. The Department of Marine Resources (DMR)
collects related data that helps ensure safety of consumption of shellfish from
Maine waters.  Data collected includes data related to potential human health
hazards from consuming shellfish, such as the presence of phytotoxins in
shellfish tissue, presence of potentially toxic phytoplankton in marine waters,
and presence of elevated levels of indicator bacteria in local wasters.

Elevated levels of indicator bacteria do not necessarily indicate an imminent
threat to human health, they are an indicator of likely contamination of coastal
waters with fecal material, especially from warm blooded organisms. While risks
of exposure to (human) pathogens are highest if pollution comes from human
sources, any source of such contamination can be problematic.

# Statement of Purpose
CBEP is committed to the ideal of open science.  Our State of the Bay data
archives ensure the science underlying the 2020 State of the Bay report is
documented and reproducible by others. The purpose of these archives is to
release raw data and data analysis code whenever possible to allow others to
review, critique, learn from, and build upon CBEP science.

# Archive Structure
CBEP 2020 State of the Bay data analysis repositories are divided into from two
to four sub-folders.  All archives contain at least an "Original_Data" and a
"Graphics" folder.  The other two folders are only included if strictly
necessary.

- Original Data.  Original data, with a "DATA_SOURCES.md" or "READ ME.txt" 
file that documents data sources.
**DATA IN THIS FOLDER IS AS ORIGINALLY PROVIDED OR ACCESSED.** 

- Derived Data.  Data derived from the original raw data.  Includes
documentation of data reorganization steps, either in the form of files (R
notebooks, Excel files, etc.) that embody data transformations, or via README.md
or DATA_NOTES.md files.

- Analysis.  Contains one or more R Notebooks proceeding through the data
analysis steps. This often includes both preliminary data analysis --
principally graphical, and detailed analysis where necessary.

- Graphics.  Contains R Notebooks stepping through development of graphics, and
also copies of resulting graphics, usually in \*.png and \*.pdf formats.  These
graphics may differ from graphics as they appear in final State of the Bay
graphical layouts.

# Summary of Data Sources
We accessed publicly available data on bacteria levels in Maine coastal waters
from the web site of Maine DMR.  These public data sets provide a synopsis of
recent observations of bacteria levels, usually from the most recent 30
observations collected from any particular location.  The summary data includes
median and 90th percentiles of those most recent observations.  These summary
values relate directly to Maine and NSSP safety standards.

Bacteria data is notoriously variable, with a handful of very high observations
scattered within a much larger number of low values.  Thus certain
statistical analyses of bacteria levels requires access to raw data, not these
sumamries.  In late 2019, We requested data directly from Maine DMR, and
received raw data  for the years 2015 through 2019 via e-mail.  All analyses are
based on these original observations.

