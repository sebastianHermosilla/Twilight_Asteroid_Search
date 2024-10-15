from astroquery.jplhorizons import Horizons
from astropy.coordinates import SkyCoord
import visualize_tracks as visual
import fileOrganizer as file_org
from astropy import units as u
from astropy.time import Time
from sbident import SBIdent
from fpdf import FPDF
from tqdm import tqdm
import pandas as pd
import numpy as np
import os

print(80*'-')
print('')
print('WELCOME to the post-findPOTATOEs pipeline')
print('This code is meant to identify what have you found in the observations')
print('')
print(80*'-')

def determine_the_nights(path_to_output_folder):
    path_to_trios = '../outputs/' + path_to_output_folder + '/trios_of_observations.txt'
    dates = file_org.read_dates_of_trios_files(path_to_trios)
    return np.unique(np.array(dates))

def correct_dec(string):
    short  = string.replace("'", ":")
    short_no_space = short.replace(" ", ":")
    return short_no_space[:-1]

def match_tracklets_with_observations(path_to_tracklets, path_to_output_folder, night, path_to_fits_files, observation_run_id='_4s'):
    print('Running function: match_tracklets_with_observations()')
    triplets_document = '../outputs/' + path_to_output_folder + '/detections_' + night + '/image_triplets_' + night + '.csv'
    image_triplet = pd.read_csv(triplets_document, sep='\t', index_col=False)
    ra_obs = []
    dec_obs = []
    date_obs = []
    for n, obs in image_triplet.iterrows():
        filea = obs.filea[:-4]
        fileb = obs.fileb[:-4]
        filec = obs.filec[:-4]
        filea_path = path_to_fits_files + filea + '.fits.fz'
        obs_date, ra, dec = file_org.extract_observation_info(filea_path)
        ra_obs.append(ra)
        dec_obs.append(dec)
        date_obs.append(obs_date)
    image_triplet['ra'] = ra_obs
    image_triplet['dec'] = dec_obs
    image_triplet['date'] = date_obs
    total_path_of_tracklets = path_to_tracklets +'tracklets_' + night + observation_run_id + '.txt'
    if os.path.isfile(total_path_of_tracklets):
        tracklet_ids = []
        coordinates_track = []
        with open(total_path_of_tracklets, "r") as file:
            for line in file:
                # Split each line into individual components
                components = line.split()

                # Store the components in respective variables
                tracklet_ids.append(night + '_' + components[0])
                coordinates_track.append(" ".join(components[4:10]))
        coord_det = SkyCoord(coordinates_track, unit=(u.hourangle, u.deg))
        str_coord = image_triplet['ra'] + ' ' + image_triplet['dec']
        coord_obs = SkyCoord(str_coord.values, unit=(u.hourangle, u.deg))
        file_of_track = []
        HDUs_tracks = []
        for t in tqdm(range(len(tracklet_ids)), desc='Progress matching tracklets with observations'):
            seps = coord_obs.separation(coord_det[t])
            minimo = np.min(seps.arcsecond)
            bool_minimo = minimo == seps.arcsecond
            if t % 3 == 0:
                file_fits = image_triplet['filea'].loc[bool_minimo].values[0]

            if t % 3 == 1:
                file_fits = image_triplet['fileb'].loc[bool_minimo].values[0]

            if t % 3 == 2:
                file_fits = image_triplet['filec'].loc[bool_minimo].values[0]

            file_of_track.append(file_fits[:-4])
            if t%3 == 0:
                n_hdu = file_org.search_hdu_by_coord(path_to_output_folder, coord_det[t], file_fits, night)
            HDUs_tracks.append(n_hdu)


        data_tracks = {'id': tracklet_ids, 'coord': coordinates_track, 'file': file_of_track, 'hdu':HDUs_tracks}
        df_track = pd.DataFrame(data_tracks)
        return df_track


def match_with_known_asteroids(tracks_df, path_to_fits_files, path_to_output_folder, night):
    print('Running function: match_with_known_asteroids()')
    triplets_document = '../outputs/' + path_to_output_folder + '/detections_' + night + '/image_triplets_' + night + '.csv'
    image_triplet = pd.read_csv(triplets_document, sep='\t', index_col=False)

    asteroids_df = pd.DataFrame()

    tracks_df['asteroid name'] = [None] * len(tracks_df)
    tracks_df['asteroid coord'] = [None] * len(tracks_df)
    tracks_df['asteroid V mag'] = [None] * len(tracks_df)
    for f, file in image_triplet.iterrows():
        filea_path = path_to_fits_files + file['filea'][:-4] + '.fits.fz'
        print(' ')
        print(f'-Looking asteroids in field: {file["filea"][:4] }')
        obs_date, ra, dec = file_org.extract_observation_info(filea_path, degree=True)
        center = SkyCoord(ra, dec, unit='deg')
        epoch = Time(obs_date)

        # Works for blanco
        search_radius = 70  # * u.arcmin
        observatory = str(807)  # observatory code for Cerro Tololo

        sbid = SBIdent(observatory, epoch, center, maglim=21, hwidth=search_radius / 60)

        if len(sbid.results) == 0:
            print('No Asteroids in this field')
            continue

        sbid = sbid.results.to_pandas()

        asteroids_names = sbid.iloc[:, 0]
        asteroids_ras = sbid.iloc[:, 1]
        asteroids_decs = sbid.iloc[:, 2].apply(correct_dec)
        asteroids_mag = sbid["Visual magnitude (V)"].values

        coords_string = asteroids_ras + ' ' + asteroids_decs
        coord_asteroids = SkyCoord(coords_string.values, unit=(u.hourangle, u.deg))
        asteroids_files = [file['filea'][:-4]]*len(asteroids_names)
        asteroids_det = [False] * len(asteroids_names)

        asteroids_data = {'Name': asteroids_names, 'coord': coords_string, 'mag': asteroids_mag, 'file': asteroids_files, 'detected':asteroids_det}
        asteroids_provitional_df = pd.DataFrame(asteroids_data)

        tracks_in_file = tracks_df.loc[tracks_df['file'] == file['filea'][:-4]]

        for t, track in tracks_in_file.iterrows():
            coord_track = SkyCoord(track['coord'], unit=(u.hourangle, u.deg))
            seps = coord_asteroids.separation(coord_track)
            n_min = np.argmin(seps.arcsecond)
            if seps.arcsecond[n_min] < 2:
                print(f' MATCH: {track["id"]} correspond to {asteroids_names[n_min]} of magnitude {asteroids_mag[n_min]}')
                tracks_df.loc[tracks_df['id'] == track['id'], 'asteroid name'] = asteroids_names[n_min]
                tracks_df.loc[tracks_df['id'] == track['id'], 'asteroid coord'] = coords_string[n_min]
                tracks_df.loc[tracks_df['id'] == track['id'], 'asteroid V mag'] = asteroids_mag[n_min]
                asteroids_provitional_df.loc[asteroids_provitional_df['Name'] == asteroids_names[n_min], 'detected'] = True

        asteroids_df = pd.concat([asteroids_df, asteroids_provitional_df])
    return tracks_df, asteroids_df

def complete_asteroids_info(asteroids_df, path_to_fits_files, path_to_output_folder, night):
    print('Running function: complete_asteroids_info()')
    coord_asteroids = SkyCoord(asteroids_df['coord'].values, unit=(u.hourangle, u.deg))

    hdu_number = []
    a = []
    e = []
    i = []
    omega = []
    w = []
    for n, ast in asteroids_df.iterrows():
        coord = coord_asteroids[n]
        n_hdu = file_org.search_hdu_by_coord(path_to_output_folder, coord,  ast['file'] + '.csv', night)
        hdu_number.append(n_hdu)

        name = ast['Name'].split()[0]
        obj = Horizons(id=name, location='500@10')
        ele = obj.elements()
        a.append(float(ele['a']))
        e.append(float(ele['e']))
        i.append(float(ele['incl']))
        omega.append(float(ele['Omega']))
        w.append(float(ele['w']))

    asteroids_df['hdu'] = hdu_number
    asteroids_df['eccentricity'] = e
    asteroids_df['a sm axis'] = a
    asteroids_df['inclination'] = i
    asteroids_df['long of ascending node'] = omega
    asteroids_df['arg of perifocus'] = w
    return asteroids_df


def create_images(df_tracks, path_to_fits_files, path_to_output_folder, night, img_size=300, radius=None):
    coord_det = SkyCoord(df_tracks['coord'].values, unit=(u.hourangle, u.deg))
    save_folder = '../outputs/' + path_to_output_folder + '/detections_' + night + '/images_tracklets'
    print(f'You can find your images in: {save_folder}')

    if not os.path.exists(save_folder):
        print('Creando el directorio')
        os.makedirs(save_folder)

    for t in range(len(df_tracks)//3):
        file1 = df_tracks['file'].iloc[t*3]
        fits_path1 = path_to_fits_files + file1 + '.fits.fz'

        saving_path = save_folder + '/' + df_tracks['id'].iloc[t * 3]
        n_hdu = df_tracks['hdu'].iloc[t * 3]

        coords_list = [coord_det[3*t], coord_det[3*t+1], coord_det[3*t+2]]

        visual.image_with_path(fits_path1, coords_list, n_hdu=n_hdu, save_path=saving_path, img_size=img_size, radius=radius)

def create_gifs_of_tracklets(df_tracks, path_to_fits_files, path_to_output_folder, night, img_size=300, radius=None):
    coord_det = SkyCoord(df_tracks['coord'].values, unit=(u.hourangle, u.deg))
    save_folder = '../outputs/' + path_to_output_folder + '/detections_' + night + '/gifs_tracklets'
    print(f'You can find your gifs in: {save_folder}')
    if not os.path.exists(save_folder):
        print('Creando el directorio')
        os.makedirs(save_folder)

    for t in range(len(df_tracks)//3):
        file1 = df_tracks['file'].iloc[t*3]
        file2 = df_tracks['file'].iloc[t*3+1]
        file3 = df_tracks['file'].iloc[t*3+2]
        fits_path1 = path_to_fits_files + file1 + '.fits.fz'
        fits_path2 = path_to_fits_files + file2 + '.fits.fz'
        fits_path3 = path_to_fits_files + file3 + '.fits.fz'
        images_list = [fits_path1, fits_path2, fits_path3]

        saving_path = save_folder + '/' + df_tracks['id'].iloc[t*3]
        n_hdu = df_tracks['hdu'].iloc[t*3]

        visual.animation_by_coord(images_list, coord_det[t*3], n_hdu=n_hdu, save_path=saving_path, img_size=img_size, radius=radius)

def pdf_creator(tracks_df, path_to_output_folder, night):

    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()

    # Title
    pdf.set_font('Arial', 'B', 18)
    pdf.cell(200, 10, txt=f"Statistical resume of observation {night}", ln=True, align='C')

    # Paragraph
    pdf.set_font('Arial', '', 12)
    pdf.ln(10)
    text_paragraph = ('The following table has a statical resume of all the detections per night.')
    pdf.multi_cell(0, 10, txt=text_paragraph)

    # Table
    pdf.ln(10)
    pdf.set_font('Arial', 'B', 12)
    pdf.cell(0, 10, f'Findings per file', ln=True)

    pdf.set_font('Arial', 'B', 10)
    pdf.cell(40, 10, 'File', 1)
    pdf.cell(40, 10, 'detected tracklets', 1)
    pdf.cell(40, 10, 'known asteroids in the field', 1)
    pdf.cell(40, 10, 'known asteroids detected', 1)
    pdf.ln()

    pdf.ln(10)
    #pdf.image("", x=10, y=None, w=100)

    pdf.output(path_to_output_folder + 'statics_of_{night}.pdf')




