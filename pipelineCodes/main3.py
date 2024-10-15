import identify_findings as find
import os
import pandas as pd

### Variables ###

# Directory to the observations
path_to_fits_files = '/home/raw/isoul/blanco/202309/'

# Directory to the output folder
path_to_output_folder = 'obs_202309'

# Size of the images and gifs
image_pixels = 100
gif_pixels = 200

### Code ###

#Determine the nights of the observations runs
nights = find.determine_the_nights(path_to_output_folder)

for night in nights:
    print(' ')
    print(f'##### Night: {night} #####')
    print(' ')
    path_to_track_features = '../outputs/' + path_to_output_folder + '/detections_' + night + '/track_features.csv'

    if not os.path.exists(path_to_track_features):
        print('No tracks on this night')
        continue

    tracks_df = pd.read_csv(path_to_track_features, header=0, index_col=0)
    find.create_gifs_of_tracklets(tracks_df, path_to_fits_files, path_to_output_folder, night, img_size=gif_pixels)
    find.create_images(tracks_df, path_to_fits_files, path_to_output_folder, night, img_size=image_pixels)
    find.pdf_creator(tracks_df, path_to_output_folder, night)

#find.create_images_of_tracklets()