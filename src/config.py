import os
from pathlib import Path

ROOT = Path(os.path.dirname(os.path.realpath(__file__))).parent.absolute()
print(ROOT)

DATASTORE = os.path.join('E:\\', 'Neuroscience', 'data', 'Stefano', 'aws', 'rawData', 'vz-data-to vs6', 'analyzed_data', 'v1', 'VS6_MsBrain_A3_VS6library_V3_LH_02-07-21')

DEFAULT = {
    'cell_by_gene': os.path.join(DATASTORE, 'cell_by_gene.csv'),

    'cell_metadata': os.path.join(DATASTORE, 'cell_metadata.csv'),

    'detected_transcripts': os.path.join(DATASTORE, 'detected_transcripts.csv'),

    'manifest': os.path.join(DATASTORE, 'images', 'manifest.json'),
}

