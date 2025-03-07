# Zonal Statistics Experiment :test_tube:


### Background

Zonal statistics is the summarization of raster data cells that intersect a vector polygon. When SNAP runs zonal statistics using the `rasterstats` package, we have [two methods](https://pythonhosted.org/rasterstats/manual.html#rasterization-strategy) to choose from: 
- the `all_touched` method selects all grid cells that touch the polygon boundary
- the `centroid` method selects only the grid cells that have their centroid contained by the polygon boundary

SNAP has chosen to use the `centroid` method, as it is more conservative; the `all_touched` method tends to grab too much data from outside of the polygon boundary. However, problems arise when the polygon is small relative to the grid cells of the raster dataset - much of the polygon area can be unaccounted for if the pixel centroids do not fall within the polygon. And sometimes no pixel centroids fall within the polygon and the result is null, even though dataset values exist under the polygon shape.
  
  
| <small>*For example, overlaying this HUC-12 polygon on a 10km resolution dataset shows that only one pixel centroid falls within the polygon, so only one data value would be used to compute zonal statistics:*</small> |
|:--:| 
|![Alt text](img/image-2.png)|

One way to remedy this is to resample the input dataset to a higher resolution (i.e., a smaller grid cell size) in order to capture more data falling within the polygon boundary. At a high enough resolution, this method also has the effect of implementing a spatially weighted average where  pixels are effectively "split" by the polygon boundary. This is called resampling, and the scale factor is the value used to divide the grid axes when creating a grid of smaller cells. 


|<small>*For example, resampling by a scale factor of 10 divides each grid axis into 10 smaller parts. One grid cell now turns into 100 smaller grid cells. This would allow 5 unique data values to be used in zonal statistics, with their frequency being determined by the area of overlap with the polygon (i.e., spatial weighting). The resampled values are populated using the nearest neighbor method, which is both computationally inexpensive and conservative because it does not interpolate any new values:*</small> |
|:--:|
|![Alt text](img/image-3.png)|

### Purpose

>This document describes an experiment designed to find a relationship between polygon area, grid cell area, and the scale factor (i.e., resampling rate) in order to develop a widely applicable and defensible resampling method for zonal statistics. 

**Why not just pick one scale factor and apply it across all zonal statistics operations?** Our heterogenous collection of polygons and datasets makes a "one size fits all" solution to resampled zonal statistics unintuitive: we have too much variation to guess at an effective value to apply universally. If we want to use a single scale factor, we first need to justify that its the right one.

**Why not choose a scale factor based on features of the vector dataset (e.g., mean or median polygon size) and/or the raster dataset (e.g., resolution)?** Setting the resampling parameters based on vector or raster dataset metadata is appealing, but is unworkable for several reasons. For one, vector data polygon sizes within some datasets (e.g., parks and protected areas) vary wildly, so mean or median values are not very meaningful or representative of the individual features and we would likely not solve the issue. We can't assume that the polygon areas in any given dataset follow a normal distribution. For two, we would need to establish explicit rules for each vector + raster dataset combination, which is not exactly scaleable. Every time we add new vector or raster data, we would need to define a new set of parameters. And maybe most importantly, our methodology would be difficult to communicate and defend to the data users, because it would vary depending on the polygons + datasets used.

>What we need is a function that uses the area of the individual polygon feature (regardless of source) & the area of the grid cells in the input dataset to determine the scale factor that will be used to resample the input data before zonal statistics:
>
>**$f$(polygon area,  dataset grid cell area) = scale factor**


## Phase 1 Experimental Dataset

One HUC-8 polygon was chosen to be the spatial domain of the experiment, and all HUC-10 and HUC-12 level polygons within this HUC-8's bounding box were fetched to make up our set of polygons. 

Empty `xarray` datasets were created at 1km, 2km, 4km, and 10km resolutions. The datasets were populated using random numbers from 0 to 1 in a normal distribution, with 25% of those values being multiplied by 10 to create "anomolies" in the dataset. (This helps to simulate real data, and guards against the fact that zonal means will always tend towards 0.5 in a random dataset without anomolies.)

|<small>*Visual representation of our dataset of HUC-10 & HUC-12 polygons with simulated gridded data at 4 resolutions:*</small> |
|:--:|
|![Alt text](img/image-4.png)|

### The experiment

For each individual polygon + dataset combination, we resampled the dataset at each scale factor from 1 (original resolution) to 12 (high resolution) and computed a zonal mean for the polygon. 

We then plotted the values for 40 randomly chosen individual polygon + dataset combinations to visually determine the scale factor at which the mean seemed to "flatline", which would indicate the optimal scale factor for resampling. (In other words, this is the scale factor where any further increase in scale factor would not meaningfully alter the result.)

For each observation, we recorded the polygon area, grid cell area, and optimal scale factor in a table for analysis.

|<small>*Example of an individual observation, where the optimal scale factor for each combination was determined by visual inspection of the plots. Here, the optimal scale factors would be 1 (1km grid), 1 (2km grid), 4 (4km grid), and 8 (10km grid):*</small>|
|:--:|
|![Alt text](img/image-9.png)|


|<small>*Plotting the whole dataset together, we can see that the zonal means become most sensitive to scale factor when resolutions are coarse relative to polygon size:*</small>|
|:--:|
|![Alt text](img/image-5.png)|
|![Alt text](img/image-6.png)|
|![Alt text](img/image-7.png)|
|![Alt text](img/image-8.png)|


### Results

 The polygon area and grid cell area were computed to a simple ratio, and a polylogarithmic function was fit to the dataset relating the ratio to the optimal scale factor. 


![Alt text](img/image-10.png)

 This would be the line of best fit if we really needed the *optimal* scale factor, but in our case we are really seeking the *maximum* scale factor. We want to capture as many observations as possible *under* this line. That way, we are sure that we are resampling to a fine enough resolution in all cases.

 To accomplish this, we use a hyperbolic function to create a line that includes the vast majority of the data points:

![Alt text](img/image-12.png)

 Hyperbolic functions can have horizontal asymptotes, and in this case we use 1 (scale factor of 1 = original resolution) so as the polygon area increases relative to grid size, we tend towards keeping the original resolution. As the polygon area decreases relative to grid size, we increase the scale factor to a maximum of about 15. 

 Looking at compute times for the resampling operation, we see the longest times for the fine resolution (1km) dataset since there are more values to compute. However, our function as defined shows that we would rarely run into a case where we are going to resample a 1km resolution dataset. These look like very reasonable compute times (<6ms) for all other cases in the experiment.

![Alt text](img/image-15.png)


### Conclusion

By using the hyperbolic function below, we can determine an appropriate scale factor to use in resampling any input dataset to match the summarizing polygon (so long as we have the area of the polygon and the area of the dataset grid cells in the same units):

> **scale factor = 350 / ( (polygon area / grid cell area) + 24 ) + 1**

Scale factors will end up being between 1 and 15, and should probably be rounded up before use in interpolation. 

For a dataset with bounding box the size of a HUC-8, interpolation using the nearest neighbor method should not be computationally expensive, even if our function leads to "oversampling" beyond what is necessary. It's important to note that only subsets of the coverage (data fetched using a bounding box) will be interpolated - not the entire dataset! - but performance may decrease if we request zonal stats for very large polygon + fine resolution dataset combination. More testing will be needed to see if this resampling methodology breaks down under those extreme cases.

## Phase 2 Experimental Dataset

TBD
