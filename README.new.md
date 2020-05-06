Treetops 

User Interface

Files Section

The Files section enables the user to select input and output files that will be used by each of the three processing phases, Smoothing, Tree Tops and Tree Crowns.

The settings file will load previous settings from a file so that a session can be resumed. When an existing settings file is selected, a pop-up will warn the user that the settings currently set in the UI will be overwritten by the saved values. Because the program automatically configures the settings file name when an input raster is selected, it may happen that the dynamically-created file already exists. When this occurs, the warning is raised. The user can choose to load the settings from the extant file, or start over with new settings.

If a non-existent settings file is entered, no warning is raised.

Whether a file is required depends on which of the three processing phases are selected. The Treetops Database and Crowns Raster, for example, are not required for smoothing, but if crowns are to be processed, the Original CHM, Smoothed CHM and Treetops Database are all required. If the required fields are not filled, the Run button will not be enabled.

Original CHM (Input)

This is the original canopy height model as a raster. 

Any raster format can be read, so long as there is an available GDAL driver (https://www.gdal.org/formats_list.html). Note that not all GDAL installations come with all drivers. The Treetops installer may not have the driver you need for a particular file, but this would be unusual as most current formats are represented.

The default band is set to 1, but the user may select another. Selecting a band that is not available will raise an error.

Smoothed CHM (Output)

This is the file that will result from the smoothing operation. This is a 32-bit floating-point raster. Allowed formats are IMG (Erdas Imagine) and GTiff (GeoTiff).

The file will be overwritten, and data written to band 1.

Treetops Database (Output)

This is the database that will contain tree top point geometries and other characteristics. Allowed formats are ESRI Shapefile and SQLite (Spatialite).

The maximum size of a Shapefile is about 2GB. If this format is selected and the output size exceeds the limit, the file will be saved as SQLite and the extension and driver selection will be updated to ".sqlite" and "SQLite", respectively. For extremely large or high-resolution datasets, SQLite is a better choice.

Crowns Raster (Output)

This is the raster that will contain delineated crowns derived from the treetops and the CHM. This is a 32-bit integer raster. Allowed formats are IMG (Erdas Imagine) and GTiff (GeoTiff).

Crowns Database (Output)

If given, the Crowns Raster will be polygonized and its attributes transferred to the database along with the crown geometries. Allowed formats are ESRI Shapefile and SQLite (Spatialite).

The maximum size of a Shapefile is about 2GB. If this format is selected and the output size exceeds the limit, the file will be saved as SQLite and the extension and driver selection will be updated to ".sqlite" and "SQLite", respectively. For extremely large or high-resolution datasets, SQLite is a better choice. If the Treetops Database exceeds the limit, the Crowns Database will also, by definition, so it will be converted automatically if the treetops phase fails.




 