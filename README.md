
# Explanation
#The R Code, used to build Fig 8 of Martin et al.,2026. It provides an astronomical parametrization for analyzing how the Earth's orbit influenced daylight during the Pliocene epoch at the Iberian Margin.

#The script uses the Laskar (2004) astronomical solution to reconstruct Earth's orbital parameters (eccentricity, obliquity, and precession) between 4.4 and 5.4 Ma.  It targets the extreme versions of the Pliocene orbit—comparing periods of maximum versus minimum eccentricity—to see how they affected the "daylight budget" at 40°N latitude. Seasons are geometrically defined based  on True Solar Longitude (λ).  The parametrization integrates Keplerian Physics, calculating the duration of each season and the instantaneous day length, while  accounting for the fact that Earth moves faster near perihelion (closest to the sun) and slower near aphelion(farthest from the sun).


# OUTPUTS:   Two plots 
    # p_hours :  plot of Total Daylight Hours Difference
    # p_absolute :  plot of  Total Accumulated Daylight Hours
# ---------------------------------------------------------------------------

#Author: OchoaD
#To be Published in Martin-Garcia et al., (Fig. 8)

# Libraries required
library(dplyr)
library(ggplot2)
library(astrochron)
library(palinsol)
library(tibble)



