# This is a script for plotting UVB radiation over a 12 month calendar year in Polson, Montana
# The UVB radiation data was calculated using the Tropospheric Ultraviolet-Visible (TUV) Model v.5.3 (https://www.acom.ucar.edu/Models/TUV/Interactive_TUV/)
# I used the latitude, longitude, and elevation of Polson, MT on the 15th day of each month in 2019 for this, 
# based on methods described in: 
 # Bromage, S., Enkhmaa, D., Baatar, T., Garmaa, G., Bradwin, G., Yondonsambuu, B., Sengee, T., Jamts, E.,
 # Suldsuren, N., McElrath, T. F., Cantonwine, D. E., Hoover, R. N., Troisi, R., & Ganmaa, D. (2019). 
 # Comparison of seasonal serum 25-hydroxyvitamin D concentrations among pregnant women in Mongolia and Boston. 
 # The Journal of steroid biochemistry and molecular biology, 193, 105427. https://doi.org/10.1016/j.jsbmb.2019.105427 
    # polson latitude: 47.6932 degrees north
    # polson longitude: 114.1631 degrees west
    # polson elevation: 2,927 feet == 0.8921 km above sea level


# load necessary packages
source("../scripts/load_R_packages.R")

# load the UVB (wavelengths: 280 nm to 320 nm) radiation data
UVBradiationData <-
read.csv(
  file = "../polson_UVB_radiation_data/UVB_radiation_per_month.csv",
  header = TRUE
  )

UVBradiationPlot <-
# plot the UVB radiation for each month with a smooth trend line
UVBradiationData %>%
  ggplot(
    data = .,
    mapping = aes(
      x = Month_number,
      y = UVB_rad_W_perSqMeter
      )
    ) +
  geom_line() +
  geom_vline(xintercept = 6,linetype = "dashed",color = "red") +
  scale_x_continuous(
    name = "Month", 
    breaks = seq(1,12,1),
    labels = UVBradiationData$Month_abbrev
      ) +
  theme_classic() +
  ylab(label = bquote(bold('UVB irradiance'~(W/m^2)))) +
  theme(
    axis.title.x = element_text(
                                family = "Arial",
                                face = "bold",
                                size = 11,
                                margin = margin(t = 20,r = 0,b = 0,l = 0)
                                ),
    axis.text = element_text(
                             family = "Arial",
                             face = "bold",
                             size = 11
                             ),
    axis.title.y = element_text(
                                family = "Arial",
                                face = "bold",
                                size = 11,
                                margin = margin(t = 0, r = 20, b = 0, l = 0)
                                ),
    axis.text.y = element_text(
                              family = "Arial",
                              face = "bold",
                              size = 11
                              )
    ) 

# save the results to the file system as a .tiff
tiff(
  filename = "../polson_UVB_radiation_data/UVBradiationPlot.tiff",
  compression = "none",
  res = 300,
  units = "in",
  width = 7,
  height = 5
  )

UVBradiationPlot
dev.off()



