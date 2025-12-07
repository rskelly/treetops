# Treetops Project Plan

## Problem
There is no tool available to locate tree tops and delineate tree crowns in the DSM of a forest canopy, derived from aerial LiDAR or photographic imagery (via SFM).

## Objective
- Develop a software application that can:
 Locate tree tops in the DSM of a forest canopy, derived from aerial LiDAR or photographic imagery (via SFM).
- Delineate a tree crown associated with each tree top, based on certain thresholds, such a crown depth and maximum extent.
- Aggregate tree tops which reside within single crowns (for trees with multiple tops).
- Produce geospatial data products, including rasters and vectors (points, polygons) which represent the location and extent of the derived features, and associated attributes.

## Requirements

- A software program with a graphical user interface which allows the user to specify input output files, select configurations, execute the process and view errors and status reports.

### Program Requirements

The user interface provides the ability to:
 - Enter all configuration parameters.
 - Select input and output files.
 - View status updates during processing.
 - View error messages.
 - Cancel processing

#### Algorithm

1. Gaussian Smoothing -- the input raster is optionally smoothed before subsequent processing.

2. Detection of "tree tops": pixels with the maximum value within a moveable region are identified as tree tops and saved for subsequent processing. Each tree top is labeled with a unique integer ID.

3. Detection of "tree crowns": a new integer raster -- the target raster -- is created with the same extents as the source raster. Beginning with each tree top, a modified flood-fill algorithm "fills" the region around each tree top, extending outwards until the pixel values reach a minimum threshold, a maximum radius or touch a neighbouring filled pixel.  

    A pixel takes the value of the associated tree top's ID and is assigned at the same position to a new integer raster with the same extent as the original raster. The fill is performed iteratively: one flood-fill step is performed at each tree top in sequence, before returning to the first tree top and performing the next fill step. This is repeated until no more pixels are available to be filled.

4. Tree crowns are merged [under which conditions] to accomodate instances where a single tree contains more than one top, such as in the case of "candelabra" cedars. For each moved polygon, the highest tree top is maintained, and the others are deleted. The merged crown taks the ID of the highest top.

5. The filled tree crown raster is vectorized: polygonal regions of pixels are converted to polygons, labeled with the associated tree top's ID and other relevant, such as the minium and maximum heights of the crown.

6. The tree crown rasters are optionally cleaned: holes in the crown polygons are removed. "Dangles" -- chains of single pixels extending from the boundary of the polygon -- are trimmed.

7. The smoothed raster is optionally saved as a new raster and the tree tops and crowns are saved to vector files.

#### Program Inputs

##### Input File 

- The program requires a single file input: a raster (grid) file which represents a forest canopy surface, in any format readable by the GDAL library. 
- The input raster's CRS and resolution will dictate the CRS and resolution of the outputs.
- Grid cell values will represent the surface height, i.e., the height of the forest canopy or ground.

##### Output Files

- When the input file is selected, output file paths are automatically filled by modifying the name of the original file to indicate the type and contents of the new files. By default, output files are saved in the same directory as the input file.
- The smoothed raster may be optionally saved as a raster with the same characteristics as the input raster.
- The tree tops and crowns are saved in vector form, as points and polygons, respectively. Vector files use the same CRS as the input file.
- The configuration is saved as a JSON file.

#### Configurations

- Configurations will be saved in a JSON file with a standard format, which can be re-loaded into any instance of the program to replicate the original process.
- Configuration files shall be versioned. Newer versions of the program will modify deprecated configurations to satisfy their requirements.

##### Input Raster

- The user will select which band from the input raster is to be processed. Defaults to 1.

##### Smoother Configuration

- The Gaussian smoother requires a standard deviation, mean and radius. Defaults to [].

##### Tree Top Configuration

- A radius is given to define the size of the window used to locate maxima. Defaults to 5 metres.

##### Tree Crown Configuration

- A radius is given to restrict the horizontal extent of a tree crown. Defaults to 10 metres.
- A percentage is given to restrict the height of the crown, measured from the tree top downwards, as a proportion of the overall height of the tree from ground to top. That is, a 100m tall tree with a crown height restriction of 0.75, will have a 75m crown extending from 25m above ground to 100m, if no other restrictions apply. Defaults to 75%.

##### Tree Crown Cleanup Configuration

- A boolean is given to dictate whether holes should be removed from crown polygons. Defaults to true.
- A boolean is given to dictate whether "dangles" should be removed from tree crown polygons. Defaults to true.

# Execution

The algorithm will be executed by a single button, and canceled with a single button.

If an error occurs, a pop-up dialog will display readable text explaining the error.

An embedded text field and status bar will indicate, in text, the operation currently being performed and, graphically, the overall progress of processing.

On completion, the program will indicate whether processing completed successfully or failed.

If processing failed, provide the user the opportunity to change any implicated inputs or configurations and restart from any previous successful stage in processing.

# Data Storage

## Treetops

- A struct representing a pixel, with a column and row field, each a 32-bit unsigned integer.
- A struct contining the cel row and column index as an unsigned, 32-bit integer, an unsigned 16-bit integer representing the ID, and a queue for storing instances of the pixel struct.
- A vector is initialized to store treetops.

## Input Raster

- A memory-mapped input segment is initialized to store pixels from the input raster. 
- The segment's datatype corresponds to that of the raster, i.e., Float32 -> float; Float64 -> double. 
- The length of the segment equals the total number of pixels in the raster band, i.e., number of columns x number of rows.

## Tree Crowns Raster

- Tree crowns are stored in a memory-mapped crown segment of equal size to the input raster's segment. The crown segment is typed as a 32-bit unsigned integer. All segment values are set to an initial value of 0.

# Algorithm

## Input Loading

1. The input file is read by GDAL using the appropriate raster driver. The following properties are extracted:
  - The coordinate reference system.
  - The width and height of the raster in pixels.
  - The affine transform which describes the geospatial bounds of the raster.
  - The data type of the raster's pixels.
2. The input raster is loaded into the input segment.
3. Using the raster's transform, a conversion factor from metres to pixels is calculated.

## Gaussian Smoothing

1. A circular Gaussian filter with the configured parameters is passed over the input segment, pixel-by-pixel. The filter's radius is converted from metres by the calculated conversion factor. The segment is modified in-place.

## Tree Top Detection

1. A circular window of the configured size is passed over the input segment. The size of the window in pixels is determined by applying the conversion factor to the radius in metres. Within each window the column and row coordinates of the pixel with the maximum value are determined. A Treetop struct is instantiated with the next available ID value, and the row and column index. The instance is added to the tree tops list. A pixel struct with the column and row indexes is added to the tree top's queue.

## Tree Crown Delineation

1. Initialize a counter.
2. For each tree top instance in the list:
    1. Increment the counter by the number of pixels in the top's queue.
    1. While the the queue's length is greater than zero:
        1. Remove the first pixel.
        2. Check that the pixel satisfies the following conditions:
    
            - It is within the maximum crown radius.
            - It is within the crown height threshold.  
   
           If these conditions are satisfied:
        
            1. Fill the pixel with the value of the top's ID.
            2. For each each of the pixel's 8 neighbouring pixels whose value is 0:
                1. Add a pixel instance with its coordinates to the top's queue.
3. If the counter's value is 0, exit.
4. Clear the counter.
5. Return to 2.

## Tree Crown Merging

?

## Tree Crown Cleanup

...