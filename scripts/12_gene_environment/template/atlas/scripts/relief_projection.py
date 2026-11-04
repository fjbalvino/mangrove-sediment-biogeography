"""Reproject the local global relief image for Leaflet's EPSG:3857 overview.

Input is the existing Natural Earth plate-carree raster (180 N-to-S degrees,
360 W-to-E degrees). Resampling changes projection, not spatial resolution.
The output is embedded in each viewer and never fetched from a map server.
"""
from pathlib import Path
import numpy as np
from PIL import Image

MAX_LATITUDE = float(np.degrees(np.arctan(np.sinh(np.pi))))


def row_latitudes(height):
    """Latitude at the centre of every north-to-south Web Mercator pixel."""
    if height <= 0:
        raise ValueError('Positive output height required')
    y = (np.arange(height, dtype=float) + 0.5) / height
    return np.degrees(np.arctan(np.sinh(np.pi * (1 - 2*y))))


def reproject_rgb(rgb):
    """Bilinear latitude sampling; retain the original longitude resolution."""
    rgb = np.asarray(rgb)
    if rgb.ndim != 3 or rgb.shape[2] != 3 or rgb.shape[1] != 2*rgb.shape[0]:
        raise ValueError('Expected a global 2:1 RGB plate-carree raster')
    height, width, _ = rgb.shape
    source_y = np.clip((90-row_latitudes(width))/180*height - 0.5, 0, height-1)
    y0 = np.floor(source_y).astype(int)
    y1 = np.minimum(y0+1, height-1)
    fraction = (source_y-y0)[:, None, None]
    values = rgb[y0].astype(float)*(1-fraction) + rgb[y1].astype(float)*fraction
    return np.clip(np.rint(values), 0, 255).astype('uint8')


def write_web_mercator(source, target):
    image = np.asarray(Image.open(source).convert('RGB'))
    output = reproject_rgb(image)
    Image.fromarray(output).save(target, format='JPEG', quality=90, optimize=True)
    return {'source_projection':'EPSG:4326 plate carree',
            'display_projection':'EPSG:3857', 'bounds':[[-MAX_LATITUDE,-180],[MAX_LATITUDE,180]],
            'source_pixels':[image.shape[1],image.shape[0]],
            'display_pixels':[output.shape[1],output.shape[0]],
            'interpretation':'Embedded overview; no extra spatial detail is inferred by resampling'}
