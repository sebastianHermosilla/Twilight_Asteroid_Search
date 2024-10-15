import identify_findings as find
import pandas as pd

### Variables ###

# Directory to the observations
path_to_fits_files = '/home/raw/isoul/blanco/202309/'

# Directory to the output folder
path_to_output_folder = 'obs_202309'

# Directory to the tracklets output of findPOTATOs
path_to_tracklets = '/home/shermosilla/Desktop/asteroids_blanco/findPOTATOs/outputs/'

# ID of the findPOTATOs-run (this can be just set as '')
id_observations = '_37s'

#### Code ####

# Determine the nights of the observation run
nights = find.determine_the_nights(path_to_output_folder)

for night in nights[:3]:
    print(' ')
    print(f'##### Night {night} #####')
    print(' ')

    saving_folder = '../outputs/' + path_to_output_folder + '/detections_' + night + '/'

    tracks_night = find.match_tracklets_with_observations(path_to_tracklets, path_to_output_folder, night, path_to_fits_files, observation_run_id=id_observations)
    #tracks_night = pd.read_csv('../outputs/obs_202309/detections_'+night+'/track_features.csv')

    if type(tracks_night) == type(None):
        print(f'no data on {night}')
        print(' ')
        continue

    # Deetermine if the tracklets corresponds to a known asteroid
    tracks_night, asteroids_df = find.match_with_known_asteroids(tracks_night, path_to_fits_files, path_to_output_folder, night)
    if type(tracks_night) != type(None):
        # Completes the information for the asteroids detected
        asteroids_df = find.complete_asteroids_info(asteroids_df, path_to_fits_files, path_to_output_folder, night)

    tracks_night.to_csv(saving_folder + 'track_features.csv', header=True)
    asteroids_df.to_csv(saving_folder + 'asteroids_features.csv', header=True)
