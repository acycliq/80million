"""
functions to create the geneData flatfiles
"""
import json
import numpy as np
import pandas as pd
import logging
import src.config as config

logger = logging.getLogger()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s"
)


def transformation(cgf):
    with open(cgf['manifest']) as f:
        settings = json.load(f)

    # bounding box in microns
    bbox = {'x0': settings['bbox_microns'][0],
            'x1': settings['bbox_microns'][2],
            'y0': settings['bbox_microns'][1],
            'y1': settings['bbox_microns'][3]}

    # image width and height in pixels
    img = {'width': settings['mosaic_width_pixels'],
           'height': settings['mosaic_height_pixels']}

    # Affine transformation: a set of coefficients a, b, c, d for transforming
    # a point of a form (x, y) into (a*x + b, c*y + d)
    a = img['width'] / (bbox['x1'] - bbox['x0'])
    b = -1 * img['width'] / (bbox['x1'] - bbox['x0']) * bbox['x0']
    c = img['height'] / (bbox['y1'] - bbox['y0'])
    d = -1 * img['height'] / (bbox['y1'] - bbox['y0']) * bbox['y0']

    tx = lambda x: a * x + b
    ty = lambda y: c * y + d
    return tx, ty


def main():
    cfg = config.DEFAULT
    # tx, ty = transformation(cfg)

    # read the spots from the raw files
    spots_path = cfg['detected_transcripts']
    chunks = pd.read_csv(spots_path, chunksize=100000)
    data = pd.concat(chunks)

    data_z3 = data[data.global_z == 3]
    [unique_gene_names, gene_id, counts_per_gene] = np.unique(data_z3.gene.values, return_inverse=True,
                                                              return_counts=True)
    data_z3 = data_z3[['gene', 'global_x', 'global_y']].rename(
        columns={'gene': "Gene", 'global_x': "x", 'global_y': "y"})

    data_z3['Gene_id'] = gene_id
    data_z3['neighbour'] = np.ones(len(gene_id)).astype(np.int)
    data_z3['neighbour_array'] = [[1] for i in range(len(gene_id))]
    data_z3['neighbour_prob'] = [[1.0] for i in range(len(gene_id))]

    logger.info('Data size: %d, %d' % (data_z3.shape[0], data_z3.shape[1]))
    data_z3 = data_z3.sort_values(by=['x', 'y'])
    splitter_mb(data_z3, 'geneData', 99)


if __name__=="__main__":
    main()
    logger.info('Done!')





