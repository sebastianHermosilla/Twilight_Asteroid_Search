from astropy.io import fits
import numpy as np
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from astropy.table import vstack, QTable
from photutils.aperture import CircularAperture
from astropy.wcs import WCS
from astropy.visualization import SqrtStretch
import matplotlib.pyplot as plt
import pandas as pd

def process_fits_file(fits_filename, output_folder, fwhm, threshold, output_csv_filename="detection_file.csv",
                      extra_info_csv_filename="extra_information.csv", band='r', plot_images=True):


    # Open the FITS file
    hdulist = fits.open(fits_filename)

    # Get the total number of HDUs
    total_hdus = len(hdulist)

    # Initialize empty lists to store sources from all HDUs
    all_sources = []
    all_extra_info = []

    # Iterate through all HDUs in the FITS file
    for hdu_index in range(total_hdus):
        # Get the HDU
        hdu = hdulist[hdu_index]

        # Skip non-image HDUs
        if hdu.data is None or len(hdu.data.shape) != 2:
            continue

        # Extract image data
        path_to_badPixels = '../outputs/' + output_folder + '/bad_pixels_masks/hdu' + str(hdu_index) + '.csv'
        df_badPixels = pd.read_csv(path_to_badPixels, header=None, index_col=False)
        # print(df_badPixels.values)
        badPixels_to_mask = df_badPixels.values == 1.0

        # Extract image data
        image = hdu.data

        # Calculate mean, median, and std for sigma clipping
        mean, median, std = sigma_clipped_stats(image, sigma=3.0)

        # Detect sources using DAOStarFinder
        daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold * std)
        mask = np.zeros(image.shape, dtype=bool)
        mask[:25, :] = True
        mask[-25:, :] = True
        mask[:, :25] = True
        mask[:, -25:] = True
        mask[image > 50000] = True
        mask[badPixels_to_mask] = True
        sources = daofind(image)

        # Create custom IDs based on HDU number, xcentroid, and ycentroid
        hdu_number_str = f'HDU{hdu_index}_'
        custom_ids = [f"{hdu_number_str}x{int(x)}_y{int(y)}" for x, y in
                      zip(sources['xcentroid'], sources['ycentroid'])]
        sources['id'] = custom_ids

        # Filter detections
        pix_value = sources['peak'] + sources['sky']
        flux = sources['flux']
        non_saturated = pix_value < 40000
        saturated = pix_value > 50000
        x = pix_value[non_saturated]
        y = flux[non_saturated]

        slope, intercept = np.polyfit(x, y, 1)
        line = pix_value * slope + intercept

        deviation = flux - line
        std = np.std(deviation)

        upper_rare = (deviation - 5 * std) > 0
        lower_rare = (deviation + 5 * std) < 0
        bool_rare = upper_rare | lower_rare

        filtered_sources = sources[~bool_rare & ~saturated]

        # Append sources and extra information to the lists
        all_sources.append(QTable(filtered_sources[['id', 'xcentroid', 'ycentroid', 'flux', 'mag']]))
        all_extra_info.append(QTable(filtered_sources[['id', 'peak', 'sharpness', 'roundness1', 'roundness2', "flux", "mag", 'npix', 'sky']]))

        # Plot the image with apertures if plot_images is True
        if plot_images:
            positions = np.transpose((filtered_sources['xcentroid'], filtered_sources['ycentroid']))
            apertures = CircularAperture(positions, r=50.0)
            fig, ax = plt.subplots(1, 1, figsize=(12, 12))
            norm = ImageNormalize(stretch=SqrtStretch())
            ax.imshow(image, cmap='Greys', origin='lower', norm=norm,
                      interpolation='nearest')
            apertures.plot(color='blue', lw=1.5, alpha=0.5)
            ax.set_title(f'HDU {hdu_index}')
            plt.show()

    # Concatenate sources and extra information from all HDUs into single tables
    all_sources = vstack(all_sources, join_type='exact')
    all_extra_info = vstack(all_extra_info, join_type='exact')

    # Extract source coordinates
    positions = np.transpose((all_sources['xcentroid'], all_sources['ycentroid']))

    # Create CircularApertures for visualization
    apertures = CircularAperture(positions, r=50.0)

    # Plot the final image with apertures if plot_images is True
    if plot_images:
        fig, ax = plt.subplots(1, 1, figsize=(12, 12))
        norm = ImageNormalize(stretch=SqrtStretch())
        ax.imshow(image, cmap='Greys', origin='lower', norm=norm,
                  interpolation='nearest')
        apertures.plot(color='blue', lw=1.5, alpha=0.5)
        ax.set_title('All HDUs Combined')
        plt.show()

    # Extract and add RA and Dec to the detection file
    ra_vec, dec_vec = [], []
    for i in range(len(all_sources)):
        x, y = all_sources["xcentroid"][i], all_sources["ycentroid"][i]
        hdu_str = all_sources['id'][i]
        hdu_ind = hdu_str.split('_')[0]
        hdu_num = int(float(hdu_ind[3:]))
        mywcs = WCS(hdulist[hdu_num].header)  # You may need to adjust this index if your WCS is in a different HDU
        ra, dec = mywcs.all_pix2world([[x, y]], 0)[0]
        ra_vec.append(ra)
        dec_vec.append(dec)

    all_sources["RA"] = ra_vec
    all_sources["Dec"] = dec_vec

    # Add MJD-OBS and observatory_code to the detection file
    mjd_time = hdulist[0].header["MJD-OBS"]
    all_sources["mjd"] = mjd_time * np.ones(len(ra_vec))
    all_sources["observatory_code"] = 807 * np.ones(len(ra_vec))
    all_sources["band"] =  np.full(len(ra_vec), band, dtype=str)
    all_sources["mag"] = all_sources["mag"].copy()

    # Close the FITS file
    hdulist.close()

    # Save the detection file and extra information as CSVs
    all_sources[['id', 'RA', 'Dec', 'mjd', 'observatory_code', 'xcentroid', 'ycentroid', 'mag', 'band']].write(output_csv_filename, format='csv', overwrite=True)
    all_extra_info.write(extra_info_csv_filename, format='csv', overwrite=True)

    return all_sources, all_extra_info