from astropy.visualization import ZScaleInterval
from astropy.coordinates import SkyCoord
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import pandas as pd
import numpy as np

def image_limits(image, pixel_x, pixel_y, half_size_x, half_size_y):
    Ny, Nx = np.shape(image)
    if (pixel_x <= half_size_x) and (pixel_y >= half_size_y):
        from_x = (0,pixel_x + half_size_x)
        from_y = (pixel_y - half_size_y,pixel_y + half_size_y)
        sub_image = image[pixel_y - half_size_y:pixel_y + half_size_y,
                                0:pixel_x + half_size_x]
    if (pixel_y <= half_size_y) and (pixel_x >= half_size_x):
        from_x = (pixel_x - half_size_x, pixel_x + half_size_x)
        from_y = (0, pixel_y + half_size_y)
        sub_image = image[0:pixel_y + half_size_y,
                                pixel_x - half_size_x:pixel_x + half_size_x]
    if (pixel_y <= half_size_y) and (pixel_x <= half_size_x):
        from_x = (0, pixel_x + half_size_x)
        from_y = (0, pixel_y + half_size_y)
        sub_image = image[0:pixel_y + half_size_y,
                                0:pixel_x + half_size_x]
    if (Nx - pixel_x <= half_size_x) and (Ny - pixel_y >= half_size_y):
        from_x = (pixel_x - half_size_x,Nx)
        from_y = (pixel_y - half_size_y,pixel_y + half_size_y)
        sub_image = image[pixel_y - half_size_y:pixel_y + half_size_y,
                                pixel_x - half_size_x:Nx]
    if (Ny - pixel_y <= half_size_y) and (Nx - pixel_x >= half_size_x):
        from_x = (pixel_x - half_size_x, pixel_x + half_size_x)
        from_y = (pixel_y - half_size_y, Ny)
        sub_image = image[pixel_y - half_size_y:Ny,
                                pixel_x - half_size_x:pixel_x + half_size_x]
    if (Ny - pixel_y <= half_size_y) and (Nx - pixel_x <= half_size_x):
        from_x = (pixel_x - half_size_x, Nx)
        from_y = (pixel_y - half_size_y, Ny)
        sub_image = image[pixel_y - half_size_y:Ny,
                                pixel_x - half_size_x:Nx]
    if (Ny - pixel_y <= half_size_y) and (pixel_x <= half_size_x):
        from_x = (0, pixel_x + half_size_x)
        from_y = (pixel_y - half_size_y, Ny)
        sub_image = image[pixel_y - half_size_y:Ny,
                                0:pixel_x + half_size_x]
    if (pixel_y <= half_size_y) and (Nx - pixel_x <= half_size_x):
        from_x = (pixel_x - half_size_x, Nx)
        from_y = (0, pixel_y + half_size_y)
        sub_image = image[0:pixel_y + half_size_y,
                                pixel_x - half_size_x:Nx]
    bool_low = (pixel_y > half_size_y) and (pixel_x > half_size_x)
    bool_up = (Ny - pixel_y > half_size_y) and (Nx - pixel_x > half_size_x)
    if bool_low and bool_up:
        from_x = (pixel_x - half_size_x, pixel_x + half_size_x)
        from_y = (pixel_y - half_size_y, pixel_y + half_size_y)
        sub_image = image[pixel_y - half_size_y:pixel_y + half_size_y,
                    pixel_x - half_size_x:pixel_x + half_size_x]
    return sub_image, from_x, from_y

def find_hdu_with_coordinate(detections_file, coordinate):
    # Read detections data from CSV file
    detections = pd.read_csv(detections_file, sep=',', header=0)

    # Find the detection closest to the given coordinate
    min_distance = float('inf')
    closest_detection = None

    for index, row in detections.iterrows():
        ra = row['RA']
        dec = row['Dec']
        det_coord = SkyCoord(ra, dec, unit='deg')
        distance = det_coord.separation(coordinate)

        if distance < min_distance:
            min_distance = distance
            closest_detection = row

    # Check if a detection was found
    if closest_detection is not None:
        # Get HDU information from the detection's id
        hdu_id = closest_detection['id']
        hdu_info = hdu_id.split('_')[0]
        x_center = int(closest_detection['xcentroid'])
        y_center = int(closest_detection['ycentroid'])
        hdu_index = int(float(hdu_info[3:]))

        # Return HDU information
        return hdu_index, x_center, y_center, hdu_id
    else:
        print("No detection found for the given coordinate.")
        return None, None, None

def image_with_path(fits_path, coords_list, n_hdu=0, save_path='animation', img_size=200, radius=None):
    coord1 = coords_list[0]
    coord2 = coords_list[1]
    coord3 = coords_list[2]

    hdulist = fits.open(fits_path)

    wcs = WCS(hdulist[n_hdu].header)

    x1_, y1_ = wcs.world_to_pixel(coord1)
    x2_, y2_ = wcs.world_to_pixel(coord2)
    x3_, y3_ = wcs.world_to_pixel(coord3)

    x1 = int(float(x1_))
    x2 = int(float(x2_))
    x3 = int(float(x3_))
    y1 = int(float(y1_))
    y2 = int(float(y2_))
    y3 = int(float(y3_))

    half_size_x, half_size_y = img_size // 2, img_size // 2
    image = hdulist[n_hdu].data

    sub_image, from_x1, from_y1 = image_limits(image, x1, y1, half_size_x, half_size_y)

    if np.isnan(sub_image).all():
        print(f'all nan in image in pix: ({x1}, {y1})')
        return

    z = ZScaleInterval()
    z1, z2 = z.get_limits((sub_image))
    plt.clf()
    plt.imshow(sub_image, cmap="gray", origin='lower', interpolation='none', vmin=z1, vmax=z2)
    plt.plot(x1 - from_x1[0], y1 - from_y1[0], 'x', color='red', markersize=3)
    plt.plot(x2 - from_x1[0], y2 - from_y1[0], 'x', color='red', markersize=3)
    plt.plot(x3 - from_x1[0], y3 - from_y1[0], 'x', color='red', markersize=3)

    theta = np.linspace(0, 2*np.pi, 100)
    if radius == None:
        radius = img_size*0.2
    r_x = radius * np.cos(theta)
    r_y = radius * np.sin(theta)
    plt.plot(x1 - from_x1[0] + r_x, y1 - from_y1[0] + r_y, color='orange')

    plt.savefig(save_path)

def animation_by_coord(images_path, coord, n_hdu=0, save_path='animation', img_size=200, radius=None):
    image_path1 = images_path[0]
    image_path2 = images_path[1]
    image_path3 = images_path[2]

    hdulist1 = fits.open(image_path1)
    hdulist2 = fits.open(image_path2)
    hdulist3 = fits.open(image_path3)

    wcs1 = WCS(hdulist1[n_hdu].header)
    wcs2 = WCS(hdulist2[n_hdu].header)
    wcs3 = WCS(hdulist3[n_hdu].header)

    x1_, y1_ = wcs1.world_to_pixel(coord)
    x2_, y2_ = wcs2.world_to_pixel(coord)
    x3_, y3_ = wcs3.world_to_pixel(coord)

    x1 = int(float(x1_))
    x2 = int(float(x2_))
    x3 = int(float(x3_))
    y1 = int(float(y1_))
    y2 = int(float(y2_))
    y3 = int(float(y3_))

    half_size_x, half_size_y = img_size // 2, img_size // 2
    image1 = hdulist1[n_hdu].data
    image2 = hdulist2[n_hdu].data
    image3 = hdulist3[n_hdu].data

    sub_image1, from_x1, from_y1 = image_limits(image1, x1, y1, half_size_x, half_size_y)
    sub_image2, from_x2, from_y2 = image_limits(image2, x2, y2, half_size_x, half_size_y)
    sub_image3, from_x3, from_y3 = image_limits(image3, x3, y3, half_size_x, half_size_y)

    if np.isnan(sub_image1).all() or np.isnan(sub_image2).all() or np.isnan(sub_image2).all():
        print(f'all nan in image in pix: ({x1}, {y1})')
        return

    images = [sub_image1, sub_image2, sub_image3]
    from_x = [from_x1, from_x2, from_x3]
    from_y = [from_y1, from_y2, from_y3]
    pos_x = [x1, x2, x3]
    pos_y = [y1, y2, y3]

    fig, ax = plt.subplots()
    def animate(i):
        ax.clear()
        z = ZScaleInterval()
        z1, z2 = z.get_limits((images[i]))
        ax.imshow(images[i], cmap="gray", origin='lower', interpolation='none', vmin=z1, vmax=z2)
        #anotate_coord(ax, images[i], from_x[i], from_y[i], wcs_i)
        if radius != None:
            fix_x = -(pos_x[i] - from_x[0][0]) + int(half_size_x*1.4)
            fix_y = -(pos_y[i] - from_y[0][0]) + int(half_size_y*1.4)
            ax.plot(fix_x, fix_y, 'x', color='red')
            ax.plot()
            if i == 0:
                theta = np.linspace(0, 2 * np.pi, 50)
                center_x = x1 - from_x[i][0]
                center_y = y1 - from_y[i][0]
                x = radius * np.cos(theta)
                y = radius * np.sin(theta)
                ax.plot(x + center_x, y + center_y, color='red')

    ani = animation.FuncAnimation(fig, animate, interval=600, frames=3)

    ani_name =  save_path + '.gif'
    ani.save(filename=ani_name, writer='pillow')
    plt.close(fig)