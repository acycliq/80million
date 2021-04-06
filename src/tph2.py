"""
Code to generate the Tph2 tiles
"""
import os
import shutil
import json
from pathlib import Path
import numpy as np
import pandas as pd
from functools import partial
import datashader as ds
from datashader import transfer_functions as tf
from datashader.utils import export_image
import src.config as config
import io
import matplotlib.pyplot as plt
import pyvips
import logging

logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)



def map_size(z):
    '''
    return the image size for each zoom level. Assumes that each map tile is 256x256
    :param z:
    :return:
    '''
    return 256 * 2 ** z


def tile_maker(zoom_levels, out_dir, img_path, z_depth='onetile'):
    # img_path = os.path.join(dir_path, 'demo_data', 'background_boundaries.tif')

    dim = map_size(zoom_levels)
    # remove the dir if it exists
    # if os.path.exists(out_dir):
    #     shutil.rmtree(out_dir)

    # now make a fresh one
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    im = pyvips.Image.new_from_file(img_path, access='sequential')

    # The following two lines add an alpha component to rgb which allows for transparency.
    # Is this worth it? It adds quite a bit on the execution time, about x2 increase
    im = im.colourspace('srgb')
    im = im.addalpha()

    logger.info('Resizing image: %s' % img_path)
    factor = dim / max(im.width, im.height)
    im = im.resize(factor)
    logger.info('Done! Image is now %d by %d' % (im.width, im.height))
    pixel_dims = [im.width, im.height]

    # sanity check
    assert max(im.width, im.height) == dim, 'Something went wrong. Image isnt scaled up properly. ' \
                                            'It should be %d pixels in its longest side' % dim

    # im = im.gravity('south-west', dim, dim) # <---- Uncomment this if the origin is the bottomleft corner

    # now you can create a fresh one and populate it with tiles
    logger.info('Started doing the image tiles ')
    im.dzsave(out_dir, layout='google', suffix='.png', background=0, skip_blanks=0, depth=z_depth)
    logger.info('Done. Pyramid of tiles saved at: %s' % out_dir)

    return pixel_dims


def load_data():
    cfg = config.DEFAULT

    # read the spots from the raw files
    cached_file = os.path.join(config.ROOT, 'src', 'parquet', 'Tph2_data.parquet.gzip')
    if Path(cached_file).is_file():
        logger.info('Loading cached file from %s' % cached_file)
        Tph2_data = pd.read_parquet(cached_file)
    else:
        spots_path = cfg['detected_transcripts']
        logger.info('Reading raw data from %s' % spots_path)
        chunks = pd.read_csv(spots_path, chunksize=100000)
        data = pd.concat(chunks)
        Tph2_data = data[data.gene == 'Tph2'][['gene', 'global_x', 'global_y']]

        target_file = os.path.join(config.ROOT, 'src', 'parquet', 'Tph2_data.parquet.gzip')
        Tph2_data.to_parquet(target_file, compression='gzip')
        logger.info('File save at %s' % target_file)
    return Tph2_data


def transformation(bbox, img):
    """
    Converts micron coordinates to pixel coordinates
    """

    # Affine transformation: a set of coefficients a, b, c, d for transforming
    # a point of a form (x, y) into (a*x + b, c*y + d)
    a = img['width'] / (bbox['x1'] - bbox['x0'])
    b = -1 * img['width'] / (bbox['x1'] - bbox['x0']) * bbox['x0']
    c = img['height'] / (bbox['y1'] - bbox['y0'])
    d = -1 * img['height'] / (bbox['y1'] - bbox['y0']) * bbox['y0']

    tx = lambda x: a * x + b
    ty = lambda y: c * y + d
    return tx, ty


def proj_factor(z, img):
    """
    calcs projected coords of point d at zoom level z
    :param d:
    :return:
    """
    # calc the factor
    dim = map_size(z)
    factor = dim / max(img['width'], img['height'])
    # out = d * factor
    return factor


def master_tile(data, img, hex_code, z, dot_size=None):
    """
    makes the map at zoom level = z to be broken down into tiles
    :param zoom_level:
    :return:
    """
    dim = map_size(z)
    proj_data = proj_factor(z, img) * data
    proj_data.y = dim - proj_data.y

    scene = ds.Canvas(x_range=[0, dim], y_range=[0, dim], plot_width=dim, plot_height=dim)
    aggregation = scene.points(proj_data, 'x', 'y')
    image = tf.shade(aggregation, cmap=[hex_code], alpha=100)
    if dot_size:
        image = tf.spread(image, px=1, shape='circle', name="spread square")
    export_image(image, 'master_tile', background=None)
    return image


def manifest(cfg):
    with open(cfg['manifest']) as f:
        settings = json.load(f)

    # bounding box in microns
    bbox = {'x0': settings['bbox_microns'][0],
            'x1': settings['bbox_microns'][2],
            'y0': settings['bbox_microns'][1],
            'y1': settings['bbox_microns'][3]}

    # image width and height in pixels
    img_shape = {'width': settings['mosaic_width_pixels'],
                 'height': settings['mosaic_height_pixels']}
    return bbox, img_shape


def get_color(gene):
    hex = [d["color"] for d in glyphsConfig['settings'] if d['gene'] == gene][0]
    return hex


def tile_generator(gene, dot_color, z):
    cfg = config.DEFAULT
    bbox, img_shape = manifest(cfg)
    Tph2_data = load_data()
    tx, ty = transformation(bbox, img_shape)
    _x = Tph2_data.global_x.apply(tx).values
    _y = Tph2_data.global_y.apply(ty).values
    point_px = pd.DataFrame({'x': _x,
                             'y': _y})
    mt = master_tile(point_px, img_shape, dot_color, z)

    tile_maker(z, 'pyramid', 'master_tile.png', z_depth='one')


if __name__ == "__main__":
    with open(r"glyphConfig.json") as f:
        glyphsConfig = json.load(f)

    gene = 'Tph2'
    hex_code = get_color(gene)
    tile_generator(gene, "#00FF00", 1)


    print('ok')



