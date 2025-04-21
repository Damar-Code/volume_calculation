import os
import numpy as np
import geopandas as gpd
from osgeo import gdal
from shapely.geometry import shape


def sample_dem(point, band, transform):
    x, y = point
    x_origin = transform[0]
    y_origin = transform[3]
    pixel_width = transform[1]
    pixel_height = transform[5]
    x_offset = int((x - x_origin) / pixel_width)
    y_offset = int((y - y_origin) / pixel_height)
    return band.ReadAsArray(x_offset, y_offset, 1, 1)[0, 0]


def get_min_height(polygon, band, transform):
    return min(sample_dem((x, y), band, transform) for x, y in polygon.exterior.coords)


def get_avg_height(polygon, band, transform):
    heights = [sample_dem((x, y), band, transform) for x, y in polygon.exterior.coords]
    return sum(heights) / len(heights) if heights else 0


def calculate_polygon_volume(polygon, dem_band, dem_transform, base_band, base_level, counting_method, step_size):
    volume = 0.0
    neg_volume = 0.0

    min_x, min_y, max_x, max_y = polygon.bounds
    x_coords = np.arange(min_x, max_x, step_size)
    y_coords = np.arange(min_y, max_y, step_size)

    for x in x_coords:
        for y in y_coords:
            point = (x, y)
            if polygon.contains(shape({'type': 'Point', 'coordinates': point})):
                height = sample_dem(point, dem_band, dem_transform)
                base_height = sample_dem(point, base_band, dem_transform) if base_band else base_level
                height_diff = abs(height - base_height)
                prism_vol = step_size * step_size * height_diff

                if counting_method == 'count_only_above' and height >= base_height:
                    volume += prism_vol
                elif counting_method == 'count_only_below' and height < base_height:
                    volume += prism_vol
                elif counting_method == 'count_above_and_below':
                    if height >= base_height:
                        volume += prism_vol
                    else:
                        neg_volume += prism_vol

    return volume, neg_volume


def run_volume_calculation(polygon_file, dem_file, output_path,
                           counting_method='count_only_above',
                           base_level_method='approximate_via_min',
                           base_level_value=None,
                           base_level_dem=None,
                           step_size=1.0,
                           output_column='volume_m3',
                           merge_result=False):
    """
    Run volume calculation for given polygon and DEM file.

    Args:
        polygon_file (str): Path to the input shapefile.
        dem_file (str): Path to the DEM raster file.
        output_path (str): Output shapefile path.
        counting_method (str): Method to calculate volume.
        base_level_method (str): How to determine base level.
        base_level_value (float): Manual base level if needed.
        base_level_dem (str): Optional secondary DEM file for base.
        step_size (float): Sampling resolution (in map units).
        output_column (str): Name of the column to store volume values.
        merge_result (bool): Placeholder if you want to expand batch processing.

    Returns:
        GeoDataFrame: Output polygon with volume data.
    """
    polygons = gpd.read_file(polygon_file)
    dem = gdal.Open(dem_file)
    dem_band = dem.GetRasterBand(1)
    dem_transform = dem.GetGeoTransform()

    base_band = None
    if base_level_method == 'use_dem_layer' and base_level_dem:
        base_dem = gdal.Open(base_level_dem)
        base_band = base_dem.GetRasterBand(1)

    results = {}

    for idx, row in polygons.iterrows():
        geom = row.geometry

        if base_level_method == 'approximate_via_min':
            base_level = get_min_height(geom, dem_band, dem_transform)
        elif base_level_method == 'approximate_via_avg':
            base_level = get_avg_height(geom, dem_band, dem_transform)
        elif base_level_method == 'use_dem_layer':
            base_level = None
        elif base_level_method == 'manual':
            base_level = base_level_value
        else:
            raise ValueError("Invalid base level method")

        volume, neg_volume = calculate_polygon_volume(
            geom, dem_band, dem_transform,
            base_band, base_level,
            counting_method, step_size
        )

        results[idx] = {
            'volume': volume,
            'neg_volume': neg_volume if counting_method == 'count_above_and_below' else None
        }

    polygons[output_column] = polygons.index.map(lambda idx: results[idx]['volume'])

    if counting_method == 'count_above_and_below':
        polygons[f'{output_column}_below'] = polygons.index.map(lambda idx: results[idx]['neg_volume'])

    polygons.to_file(output_path)
    print(f'✅ Volume calculation complete. Shapefile saved to: {output_path}')

    return polygons


# HOW TO RUN!!!
# result = run_volume_calculation(
#     polygon_file= 'G:/TEAM SHARE/Working Directory/Damar/project/wood disposal/vector/wood_disposal.shp',
#     dem_file= 'G:/TEAM SHARE/Working Directory/Damar/project/wood disposal/raster/dsm_05.tif',
#     output_path= 'G:/TEAM SHARE/Working Directory/Damar/project/wood disposal/vector/wood_disposal_volume.shp',
#     counting_method='count_only_above',
#     base_level_method='approximate_via_min',
#     step_size=1.0,
#     output_column='volume_m3'
# )

# print(result.head())
# print(f'Total volume: {result["volume_m3"].sum():.2f} m³')
