# This is a Pipeline to look for asteroids using the Blanco or VST telescopes
# This has been developed by Sebastián Hermosilla Canobra
# If you have any question write to: s.hermosilla@ing.uchile.cl
import fileOrganizer as file_org

# Directory to the observations
directory = '/home/raw/isoul/blanco/202309'

# Directory where the information of the observations is going to be saved
output_folder = 'obs_202309'

band = 'r'

# Detection parameter
fwhm= 7.66 # fit width half maximum in gaussian fit
threshold = 3.7 # threshold in units of standar deviation

# Find the useful files in the directory
file_org.find_science_observations(directory, output_folder)

# Determine the observations that were taken the same night and create trios
file_org.create_trios_of_observations(output_folder)

# Determine bad pixel from observations
#file_org.create_bad_pixels_mask(directory, output_folder)

# Detect the sources in every observation
file_org.detect_sources_in_observations(output_folder, directory, fwhm, threshold, band=band)