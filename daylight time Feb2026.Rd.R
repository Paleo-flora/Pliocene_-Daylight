
###############################################################################
# Pliocene Seasonal Daylight Lengths: Perihelion vs Aphelion
# La2004 solution | 40°N | 4.4–5.4 Ma
###############################################################################

# Explanation
#The R Code, used to build Fig XX, provides an astronomical parametrization for analyzing how the Earth's orbit influenced daylight during the Pliocene epoch at the Iberian Margin.

#The script uses the Laskar (2004) astronomical solution to reconstruct Earth's orbital parameters (eccentricity, obliquity, and precession) between 4.4 and 5.4 Ma.  It targets the extreme versions of the Pliocene orbit—comparing periods of maximum versus minimum eccentricity—to see how they affected the "daylight budget" at 40°N latitude. Seasons are geometrically defined based  on True Solar Longitude (λ).  The parametrization integrates Keplerian Physics, calculating the duration of each season and the instantaneous day length, while  accounting for the fact that Earth moves faster near perihelion (closest to the sun) and slower near aphelion(farthest from the sun).

# ---------------------------------------------------------------------------
# OUTPUTS:   Two plots 
    # p_hours :  plot of Total Daylight Hours Difference
    # p_absolute :  plot of  Total Accumulated Daylight Hours
# ---------------------------------------------------------------------------


library(dplyr)
library(ggplot2)
library(astrochron)
library(palinsol)
library(tibble)


# ---------------------------------------------------------------------------
# 1. DATA
# ---------------------------------------------------------------------------

la04 <- getLaskar(sol = "la04", verbose = FALSE)

# Pliocene interval (4.4–5.4 Ma)
pliocene <- la04 %>%
  filter(Time_ka >= 4400 & Time_ka <= 5400) %>%
  mutate(
    ecc  = ecc_LA04,
    eps  = tilt_LA04,
    # Correct longitude of perihelion (ϖ), radians
    varpi = (pi - prec_LA04) %% (2*pi)
  ) %>%
  select(Time_ka, ecc, eps, varpi)

# Eccentricity extrema
scenarios <- bind_rows(
  pliocene %>% slice_max(ecc, n = 1) %>% mutate(scenario = "Eccentricity Max"),
  pliocene %>% slice_min(ecc, n = 1) %>% mutate(scenario = "Eccentricity Min")
)

# ---------------------------------------------------------------------------
# 2. SEASON DEFINITIONS (TRUE SOLAR LONGITUDE)
# ---------------------------------------------------------------------------

season_defs <- tibble(
  season = c("Spring", "Summer", "Autumn", "Winter"),
  l_start = c(0, pi/2, pi, 3*pi/2),
  l_end   = c(pi/2, pi, 3*pi/2, 2*pi),
  l_mid   = c(pi/4, 3*pi/4, 5*pi/4, 7*pi/4)
)

# Identify perihelion vs aphelion side
orbital_side <- function(l_mid, varpi) {
  d <- abs((l_mid - varpi + pi) %% (2*pi) - pi)
  ifelse(d < pi/2, "Perihelion-side", "Aphelion-side")
}

# ---------------------------------------------------------------------------
# 3. DAYLIGHT INTEGRATION (KEPLERIAN)
# ---------------------------------------------------------------------------

integrate_daylight <- function(l_start, l_end, ecc, eps, varpi, lat_deg = 40) {
  
  ### Day length for every small step of solar longitude (solar angle) - usign Kepler to consider speed changes (peri vs aph)
  
  lat_rad <- lat_deg * pi / 180
  lambda <- seq(l_start, l_end, length.out = 1000)
  d_lambda <- lambda[2] - lambda[1]
  
  # A. Calculate Instantaneous Day Length (Geometric)
  # Declination: delta = asin(sin(eps) * sin(lambda))
  delta <- asin(sin(eps) * sin(lambda))
  arg <- -tan(lat_rad) * tan(delta)
  arg <- pmax(-1, pmin(1, arg)) 
  omega <- acos(arg)
   day_hours <- (24 / pi) * omega
  
  # B. Calculate Time Derivative (dt/dlambda) - Kepler's 2nd Law
  nu <- lambda - varpi # True anomaly
  year_len <- 365.2422 # Tropical year approx
    term1 <- (year_len / (2 * pi))
  term2 <- (1 - ecc^2)^1.5
  term3 <- (1 + ecc * cos(nu))^2
    dt_dlambda <- term1 * (term2 / term3)
  
  # C. Integrate: Sum(DayLength * dt)
  # Integral = Sum( DayHours(lambda) * (dt/dlambda) * d_lambda )
  total_hours <- sum(day_hours * dt_dlambda * d_lambda)
  
  return(total_hours)
}

# ---------------------------------------------------------------------------
# 4. SEASON CALCULATION
# ---------------------------------------------------------------------------

compute_seasons <- function(row) {
    orbit <- c(eps = row$eps, ecc = row$ecc, varpi = row$varpi)
  
  season_defs %>%
    rowwise() %>%
    mutate(
      # 1. Season Length (Days) -
      d_start = l2day(orbit, l_start), 
      d_end   = l2day(orbit, l_end),
      length_days = (d_end - d_start) %% 365.2422, #
      
      # 2. Total Daylight Hours 
      total_light_hours = integrate_daylight(l_start, l_end, row$ecc, row$eps, row$varpi, lat_deg=40),  #total "budget" of light for that season.
      
      # 3. Metadata
      orbital_position = orbital_side(l_mid, row$varpi),
      scenario = row$scenario
    ) %>%
    ungroup() %>%
    select(scenario, season, orbital_position, length_days, total_light_hours)
}

df_results <- scenarios %>%
  split(1:nrow(.)) %>%
  lapply(compute_seasons) %>%
  bind_rows()


compute_seasons_absolute <- function(row) { # Seasons Absolute 
  season_defs %>%
    rowwise() %>%
    mutate(
      total_light_hours = integrate_daylight(l_start, l_end, row$ecc, row$eps, row$varpi, lat_deg=40),
      scenario = row$scenario
    ) %>%
    ungroup() %>%
    select(scenario, season, total_light_hours)
}

df_absolute <- scenarios %>%
  split(1:nrow(.)) %>%
  lapply(compute_seasons_absolute) %>%
  bind_rows()

# Set factor order for logical plotting
df_absolute$season <- factor(df_absolute$season, levels = c("Spring", "Summer", "Autumn", "Winter"))


# ---------------------------------------------------------------------------
# 5. ORBTAL COMPARISON 
# ---------------------------------------------------------------------------

# Calculate Anomalies (Max Eccentricity - Min Eccentricity)
ref_vals <- df_results %>%
  filter(scenario == "Eccentricity Min") %>%
  select(season, ref_hours = total_light_hours, ref_days = length_days)

df_anomalies <- df_results %>%
  filter(scenario == "Eccentricity Max") %>%
  left_join(ref_vals, by = "season") %>%
  mutate(
    diff_hours = total_light_hours - ref_hours,
    diff_days = length_days - ref_days
  )

print(df_anomalies)

# Plot: Total Daylight Hours Difference
p_hours <- ggplot(df_anomalies, aes(x = season, y = diff_hours, fill = orbital_position)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("Perihelion-side" = "#d73027", "Aphelion-side" = "#4575b4")) +
  labs(
    title = "Difference in Total Daylight Hours (40°N)",
    subtitle = "Eccentricity Max vs Min (Combined effect of Day Length + Season Duration)",
    y = "Change in Total Daylight Hours",
    x = NULL
  ) +
  theme_minimal(base_size = 14)

print(p_hours)


# Plot: Total Accumulated Daylight Hours
p_absolute<- ggplot(df_absolute, aes(x = season, y = total_light_hours, fill = scenario)) +
  geom_col(width = 0.7, color = "black", alpha = 0.8) +
  geom_text(aes(label = round(total_light_hours, 0)), 
            vjust = -0.5, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = c("Eccentricity Min" = "#4575b4",      "Eccentricity Max" = "#d73027"  )) +
   facet_wrap(~scenario) + 
  
  labs(
    title = "Total Daylight Hours per Season: Pliocene (40°N)",
    subtitle = "Comparing Eccentricity Min (Blue) vs Eccentricity Max (Red)",
    y = "Total Accumulated Daylight Hours",
    x = NULL
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none", 
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold", size = 12)
  )

print(p_absolute)

# -------------- The end
