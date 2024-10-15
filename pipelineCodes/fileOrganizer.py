from astropy.io import fits
import detection_functions as det_f
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.time import Time
import astropy.units as u
import pandas as pd
import numpy as np
import os
from tqdm import tqdm


def process_fits_file(file_path):
    with fits.open(file_path) as hdulist:
        header = hdulist[0].header
        proctype_value = header.get('PROCTYPE', '').strip()
        prodtype_value = header.get('PRODTYPE', '').strip()

        return 'InstCal' in proctype_value and 'image' in prodtype_value

def get_coordinates(file_path):
    with fits.open(file_path) as hdulist:
        header = hdulist[0].header
        ra = header.get('RA', None)
        dec = header.get('DEC', None)
        c = SkyCoord(ra + " " + dec, unit=(u.hourangle, u.deg))
        return c.ra.degree, c.dec.degree

def get_mjd(file_path):
    with fits.open(file_path) as hdulist:
        header = hdulist[0].header
        mjd = header.get('MJD-OBS', None)
        return mjd

def extract_observation_info(fits_path, degree=False):
    # Open the FITS file
    with fits.open(fits_path) as hdul:
        # Extract header information
        header = hdul[0].header

        # Extract observation date
        observation_date = header['DATE-OBS']
        observation_time = Time(observation_date, format='isot', scale='utc')

        # Extract RA and DEC
        ra = header['RA']
        dec = header['DEC']

        # Convert RA and DEC to HH MM SS.d format
        sky_coord = SkyCoord(ra + ' ' + dec, unit=(u.hourangle, u.deg))
        ra_str = sky_coord.ra.to_string(u.hour, sep=' ', precision=1, pad=True)
        dec_str = sky_coord.dec.to_string(u.degree, sep=' ', precision=1, alwayssign=True, pad=True)
        time_string = observation_time.datetime.strftime('%Y %m %d')
        decimal_hours = (observation_time.datetime.hour + observation_time.datetime.minute/60 + observation_time.datetime.second/(60*60))/24
        time_string += f'{decimal_hours:5f}'.lstrip('0')

        if degree:
            return  observation_time, sky_coord.ra.deg, sky_coord.dec.deg

    return time_string, ra_str, dec_str

def search_hdu_by_coord(path_to_output_folder, coord, file_fits, night):
    path_to_detections = '../outputs/' + path_to_output_folder + '/detections_' + night + '/' + file_fits
    det = pd.read_csv(path_to_detections, header=0)

    coord_det = SkyCoord(det['RA'], det['Dec'], unit='deg')
    sep = coord.separation(coord_det)
    n_min = np.argmin(sep.arcsecond)

    id_min = det['id'].iloc[n_min]
    hdu_str = id_min.split('_')[0]
    return int(hdu_str[3:])

def find_science_observations(directory, output_folder):
    print('')
    print(80 * '-')
    print('')
    print('WELCOME TO THE ASTEROID PIPELINE.')
    print('')
    print(80 * '-')
    print('')
    print('In the output folder you will find all your observations organized in .txt files')
    print('')
    print('STEP 1:')
    print('Running the function: find_science_observations')
    output_directory = '../outputs/' + output_folder
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    output_filename = output_directory + '/science_observations_file.txt'
    with open(output_filename, 'w') as output_file:
        output_file.write("Filename\tRA (deg)\tDEC (deg)\tMJD\n")  # Header

        for filename in os.listdir(directory):
            if filename.endswith('.fits.fz'):
                file_path = os.path.join(directory, filename)

                if process_fits_file(file_path):
                    ra, dec = get_coordinates(file_path)
                    mjd = get_mjd(file_path)
                    output_file.write(f'{filename}\t{ra:.6f}\t{dec:.6f}\t{mjd}\n')
                    #print(f'File {filename} meets the conditions.')
    print(f'The science observations have been found and saved into {output_filename }')

def time_convert(time, time_format="mjd"):
    t = Time(time, format=time_format)
    return t


def group_observations(obs_data):
    # Initialize an empty list to store observation groups
    observation_groups = []

    while not obs_data.empty:
        # Take out the first observation from the catalog
        current_observation = obs_data.iloc[0]

        # Create a SkyCoord object for the current observation
        current_coord = SkyCoord(ra=current_observation['RA (deg)'] * u.deg,
                                 dec=current_observation['DEC (deg)'] * u.deg)

        # Create a group for the current observation
        current_group = [current_observation]

        # Remove the current observation from the catalog
        obs_data = obs_data.iloc[1:]

        # Look for observations within 4 hours with similar coordinates
        for _, obs in obs_data.iterrows():
            # Check if the observation is within 4 hours
            if (obs['datetime'] - current_observation['datetime']).value <= 0.6:
                # Create a SkyCoord object for the comparison observation
                comparison_coord = SkyCoord(ra=obs['RA (deg)'] * u.deg,
                                            dec=obs['DEC (deg)'] * u.deg)

                # Check if the coordinates are similar (within 0.1 degree)
                if current_coord.separation(comparison_coord).degree <= 0.1:
                    # Add the observation to the current group
                    current_group.append(obs)

                    # Remove the observation from the catalog
                    obs_data = obs_data.drop(_)

        # Add the current group to the list of observation groups
        observation_groups.append(current_group)

    return observation_groups


def create_trios_of_observations(folder_with_science_observations):
    print('')
    print('STEP 2:')
    print('Running the function: create_trios_of_observations')
    input_directory = '../outputs/' + folder_with_science_observations + '/science_observations_file.txt'
    obs_data = pd.read_csv(input_directory, delimiter='\t')

    # Convert MJD to astropy Time objects for easier manipulation
    obs_data['datetime'] = obs_data['MJD'].apply(time_convert)

    # Order observations by datetime
    obs_data = obs_data.sort_values(by='datetime')

    # Call the function with your observations
    observation_groups = group_observations(obs_data)

    output_file_path = '../outputs/' + folder_with_science_observations + '/trios_of_observations.txt'

    with open(output_file_path, 'w') as output_file:
        for i, group in enumerate(observation_groups, start=1):
            output_file.write(f"Group {i}:\n")
            for obs in group:
                output_file.write(f"{obs['Filename']} {obs['RA (deg)']:.4f} {obs['DEC (deg)']:.4f} {obs['datetime'].iso}\n")
            output_file.write("###GROUP###" + "\n")
    print(f"Observation trios written to {output_file_path}")


def read_dates_of_trios_files(trios_document):
    with open(trios_document, "r") as file:
        # Split the file into groups
        groups = (file.read().split("###GROUP###\n"))[:-1]
        # Process each group
        dates = []
        for group in groups:
            group_lines = group.strip().split('\n')
            dates.append(group_lines[1].split()[-2].replace("-", ""))
        return dates

def create_folder(date, folder):
    folder_name = folder + f"detections_{date}"
    os.makedirs(folder_name, exist_ok=True)
    return folder_name

def split_groups_of_trios(trios_document):
    with open(trios_document, "r") as file:
        # Split the file into groups
        groups = (file.read().split("###GROUP###\n"))[:-1]
        # Process each group
    return groups

def identify_groups(dates, date):
    bool = (np.array(dates) == date)
    return bool
def create_detection_file(folder_name, csv_filename, groups, groups_bool):
    print(f'Creating file: {csv_filename}')
    with open(os.path.join(folder_name, csv_filename), "w") as csv_file:
        csv_file.write("filea\tfileb\tfilec\n")
        for i, group in enumerate(groups):
            if groups_bool[i] :
                group_lines = group.strip().split('\n')
                if len(group_lines) == 4:
                    for observation_line in group_lines[1:]:
                        fits_filename = observation_line.split()[0][:-8]
                        csv_file.write(f"{fits_filename}.csv\t")
                    csv_file.write(f'\n')
                if len(group_lines) > 4:
                    for observation_line in group_lines[1:4]:
                        fits_filename = observation_line.split()[0][:-8]
                        csv_file.write(f"{fits_filename}.csv\t")
                    csv_file.write(f'\n')


def create_bad_pixels_mask(fits_path, output_folder, obs_n=None):
    print('')
    print('STEP 3:')
    print('Running the function: create_bad_pixels_mask')
    print('If this take too long, try adjusting the parameter obs_n')

    calibrated_obs_files = '../outputs/' + output_folder + '/science_observations_file.txt'
    obs_data = pd.read_csv(calibrated_obs_files, delimiter='\t')
    if obs_n != None:
        obs_data = obs_data.iloc[:obs_n]

    output_directory = '../outputs/' + output_folder + '/bad_pixels_masks'
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    output_file_path = output_directory

    # Open a file to save the number of CCDs that contain the observations
    fits_file = fits_path + '/' + obs_data['Filename'].iloc[0]
    hdulist = fits.open(fits_file)
    n_hdus = len(hdulist)               # Number of CCDs
    N_files = len(obs_data) # Number of observations

    #Loop to create the maske in every CCD
    # 1 is not entered as it is usually empty
    total_bad_pix = 0 # Calculate the amount of bad pixels in the observation
    total_good_pix = 0  # Calculate the amount of bad pixels in the observation
    for i in tqdm(range(1, n_hdus), desc='Progress detecting bad pixels'):
        acum_data = []

        for j in range(N_files):
            fits_file = fits_path + '/' + obs_data['Filename'].iloc[j]
            hdulist = fits.open(fits_file)
            data = hdulist[i].data / np.nanmedian(hdulist[i].data)
            acum_data.append(data)
            hdulist.close()
        acum_data = np.array(acum_data)
        median_ccd = np.nanmedian(acum_data, axis=0)
        mediana_pix = np.nanmedian(median_ccd)
        iqr = np.nanpercentile(median_ccd, 75) - np.nanpercentile(median_ccd, 25)
        dif = np.abs(median_ccd - mediana_pix)
        bool_ccd = dif > 3 * iqr
        total_bad_pix += np.sum(bool_ccd)
        total_good_pix += np.sum(~bool_ccd)
        median_ccd[bool_ccd] = 1.0
        median_ccd[~bool_ccd] = 0.0
        ccd_df = pd.DataFrame(median_ccd)
        hdu_file = output_file_path + '/hdu' + str(i) + '.csv'
        ccd_df.to_csv(hdu_file, index=False, header=False)
    print(f'A mask for every CCD has been saved in {output_file_path}')
    percentage = np.round(100*total_bad_pix/(total_bad_pix+total_good_pix), 2)
    print(f'A total of {total_bad_pix} bad pixels have been found in the instrument, which corresponds to {percentage} %')

def run_source_detection(group_lines, output_folder, folder_name, data_directory, fwhm, threshold, band='r'):
    for observation_line in group_lines[1:]:
        fits_filename = data_directory + '/' + observation_line.split()[0]
        output_filename = f"{observation_line.split()[0][:-8]}.csv"
        extrainfo_filename = "extra_info_" + output_filename
        print(f'-Detecting sources in: {os.path.join(folder_name, output_filename)}')

        # Run source_detection function and save the result in the folder
        det_f.process_fits_file(fits_filename, output_folder, fwhm, threshold, output_csv_filename=os.path.join(folder_name, output_filename), extra_info_csv_filename=os.path.join(folder_name, extrainfo_filename), band=band, plot_images=False)


def detect_sources_in_observations(folder_with_trios, data_directory, fwhm, threshold, band='r'):
    print("")
    print("STEP 4:")
    print("Building the necesary to run findPOTATOs")
    trios_document = '../outputs/' + folder_with_trios + '/trios_of_observations.txt'
    dates = read_dates_of_trios_files(trios_document)
    groups = split_groups_of_trios(trios_document)
    current_date = None
    for i, date in enumerate(dates):
        group = groups[i]
        if current_date != date:
            current_date = date
            folder_name = create_folder(date, '../outputs/' + folder_with_trios + '/')
            print('')
            print(f'New folder: {folder_name}')
            csv_filename = f"image_triplets_{date}.csv"
            groups_bool = identify_groups(dates, date)
            create_detection_file(folder_name, csv_filename, groups, groups_bool)

        group_lines = group.strip().split('\n')
        if len(group_lines) == 4:
            print('')
            print(f'Meet the conditions {group_lines[0]}')
            run_source_detection(group_lines, folder_with_trios, folder_name, data_directory, fwhm, threshold, band)

        if len(group_lines) > 4:
            print('')
            print(f'Meet the conditions {group_lines[0]}')
            run_source_detection(group_lines[:4], folder_with_trios, folder_name, data_directory, fwhm, threshold, band)