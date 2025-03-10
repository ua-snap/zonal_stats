import requests
import geopandas as gpd
import pandas as pd
import xarray as xr
import numpy as np
import rasterio
import rioxarray
from rasterio.features import rasterize


# config
# GS_BASE_URL = "https://gs.earthmaps.io/geoserver/"
GS_BASE_URL = "https://gs.mapventure.org/geoserver/"


# function to get the WFS polygon request url
def generate_wfs_huc8_url(workspace, huc8_id):
    wfs_url = (
        GS_BASE_URL
        + f"wfs?service=WFS&version=2.0.0&request=GetFeature&typeName={workspace}&outputFormat=application/json"
    )
    wfs_url += "&propertyName=(id,area_type,the_geom)"
    wfs_url += f"&filter=<Filter><PropertyIsEqualTo><PropertyName>id</PropertyName><Literal>{huc8_id}</Literal></PropertyIsEqualTo></Filter>"
    return wfs_url


# function to get the WFS polygon request url for HUC10 and HUC12 contained by the HUC8
def generate_wfs_huc10_huc12_within_bounds_url(workspace, huc8_bounds_str):
    wfs_url = (
        GS_BASE_URL
        + f"wfs?service=WFS&version=1.0.0&request=GetFeature&typeName={workspace}&propertyName=(id,area_type,the_geom)&outputFormat=application/json&cql_filter=BBOX(the_geom, {huc8_bounds_str}) AND area_type IN ('HUC10', 'HUC12')"
    )
    return wfs_url


# function to get all huc10 and huc12 polygons that are within a huc8, given the huc8 id
def get_huc10_and_huc12_polygons_within_huc8(huc8_id):

    all_gdfs = []

    # first get the huc8 polygon
    huc8_url = generate_wfs_huc8_url(
        "playground:all_areas",
        huc8_id,
    )

    with requests.get(huc8_url) as r:
        huc8_gdf = gpd.GeoDataFrame.from_features(r.json())
        huc8_gdf.set_geometry("geometry", inplace=True)
        huc8_gdf.set_crs("EPSG:4326", inplace=True)

    huc8_bounds = huc8_gdf.geometry[0].bounds
    huc8_bounds_str = (
        f"{huc8_bounds[0]}, {huc8_bounds[1]}, {huc8_bounds[2]}, {huc8_bounds[3]}"
    )

    # now get the huc10 and huc12 polygons that are within the huc8
    huc10_huc12_url = generate_wfs_huc10_huc12_within_bounds_url(
        "playground:all_areas", huc8_bounds_str
    )

    with requests.get(huc10_huc12_url) as r:
        huc10_huc12_gdf = gpd.GeoDataFrame.from_features(r.json()["features"])
        huc10_huc12_gdf.set_crs("EPSG:4326", inplace=True)

    # add all geodataframes to the list
    all_gdfs.append(huc8_gdf)
    all_gdfs.append(huc10_huc12_gdf)

    # find the union of all geodataframes and calculate the bounding box
    union_gdf = gpd.GeoDataFrame(pd.concat(all_gdfs, ignore_index=True))
    union_gdf = union_gdf.to_crs("EPSG:3338")
    bbox = union_gdf.total_bounds

    # for each polygon in the geodataframe, extract an individual geodataframe, the area in km2, and the area_type
    # populate a dictionary using the id of the polygon as the key
    polygons = {}
    for _, row in union_gdf.iterrows():
        polygons[row["id"]] = {
            "gdf": union_gdf[union_gdf["id"] == row["id"]],
            "area_km2": row["geometry"].area / 1e6,  # convert from m2 to km2
            "area_type": row["area_type"],
        }

    print(f"Found {len(polygons)} polygons within the HUC8 {huc8_id} bounding box.")
    # get mean area of all polygons with area_type HUC10 and print a message
    huc10_polygons = [
        poly for poly in polygons.values() if poly["area_type"] == "HUC10"
    ]
    if huc10_polygons:
        mean_area_km2 = np.mean([poly["area_km2"] for poly in huc10_polygons])
        print(
            f"Mean area of {len(huc10_polygons)} HUC10 polygons: {mean_area_km2:.2f} km²"
        )
    # get mean area of all polygons with area_type HUC12 and print a message
    huc12_polygons = [
        poly for poly in polygons.values() if poly["area_type"] == "HUC12"
    ]
    if huc12_polygons:
        mean_area_km2 = np.mean([poly["area_km2"] for poly in huc12_polygons])
        print(
            f"Mean area of {len(huc12_polygons)} HUC12 polygons: {mean_area_km2:.2f} km²"
        )

    return polygons, bbox


# function to create an empty xarray dataset with band name "data"
# dataset will have user-defined grid resolution (in meters),
# CRS, and bounding box
def create_empty_dataset(resolution, bbox, crs="EPSG:3338"):

    minx = bbox[0]
    miny = bbox[1]
    maxx = bbox[2]
    maxy = bbox[3]

    ds = xr.DataArray(
        np.zeros(
            (
                int(maxy - miny + resolution) // resolution,
                int(maxx - minx + resolution) // resolution,
            )
        ),
        dims=("y", "x"),
        coords={
            "y": np.arange(
                miny,
                maxy,
                resolution,
            ),
            "x": np.arange(
                minx,
                maxx,
                resolution,
            ),
        },
        attrs={"crs": crs},
    ).to_dataset(name="data", promote_attrs=True)

    return ds


# function to populate a dataset with random numbers to simulate climate data
# use a normal distribution and user-defined anomaly rate
def populate_dataset(ds, anomaly_rate=0.25):

    ds["data"] = xr.DataArray(
        np.random.normal(0, 1, ds.data.shape),
        dims=("y", "x"),
        coords={
            "y": ds.y,
            "x": ds.x,
        },
    )
    # add anomalies to the data
    ds["data"] = ds["data"].where(
        np.random.rand(*ds.data.shape) > anomaly_rate,
        10 * np.random.normal(0, 1, ds.data.shape),
    )

    # add some additional structure to the data so we can orient it visually:
    # make the top row the max value
    ds["data"][0, :] = 1
    # make the left column the min value
    ds["data"][:, 0] = 0

    return ds


# function to interpolate the dataset using a multiplier
def interpolate_dataset(ds, x_dim, y_dim, multiplier, method):
    """Interpolate the xarray dataset to a higher resolution for more accurate zonal stats.

    Args:
        ds (xarray.DataSet): xarray dataset returned from fetching a bbox from a coverage
        x_dim (str): name of the x dimension
        y_dim (str): name of the y dimension
        mulitplier (int): multiplier to increase the resolution by
        method (str): method to use for interpolation

    Returns:
        ds_new (xarray.DataSet): xarray dataset interpolated to higher resolution
    """
    X = x_dim
    Y = y_dim

    new_lon = np.linspace(ds[X][0].item(), ds[X][-1].item(), ds.sizes[X] * multiplier)
    new_lat = np.linspace(ds[Y][0].item(), ds[Y][-1].item(), ds.sizes[Y] * multiplier)

    ds_new = ds["data"].interp(method=method, coords={X: new_lon, Y: new_lat})

    return ds_new


# function to rasterize a polygon geodataframe at the resolution of a given dataset
def rasterize_polygon(ds, x_dim, y_dim, poly_gdf):

    rasterized_polygon_array = rasterize(
        [(poly_gdf.geometry.iloc[0].geoms[0], 1)],
        out_shape=(
            ds[y_dim].values.shape[0],
            ds[x_dim].values.shape[0],
        ),  # must be YX order for numpy array!
        transform=ds.rio.transform(
            recalc=True
        ),  # must recalc since we interpolated, otherwise the old stored transform is used and rasterized polygon is not aligned
        fill=0,
        all_touched=False,
    )

    return rasterized_polygon_array


# function to calculate zonal stats for a dataset and a rasterized polygon array
# dataset and array must have the same shape, and the polygon area must be value of 1 in the array
def calculate_zonal_stats(ds, poly_array):
    data_array = ds.values
    zonal_stats = {}

    # extract values from the dataset array using the rasterized polygon array
    values = data_array[poly_array == 1].tolist()

    # calculate zonal stats
    if values:
        # zonal_stats["min"] = np.nanmin(values)
        # zonal_stats["max"] = np.nanmax(values)
        zonal_stats["mean"] = np.nanmean(values)
        zonal_stats["count"] = len(values)
    else:
        # zonal_stats["min"] = np.nan
        # zonal_stats["max"] = np.nan
        zonal_stats["mean"] = np.nan
        zonal_stats["count"] = 0

    return zonal_stats
