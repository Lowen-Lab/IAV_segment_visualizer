# IAV segment visualizer

## What this script does
* This script generates scale genome diagrams for Influenza A virus, adding the location of features to scale.
* One of the required inputs is a file that has the locations and desired labels, colors, etc (if any) of genetic features that you want to include.
* It has multiple color schemes built in, all picked to be red/green color blindness friendly.

![alt text](https://github.com/Lowen-Lab/IAV_segment_visualizer/blob/main/example_IAV_diagram_figure.png)

# How to use this script
The interface for using this script is to open it in a text editor like sublime text and edit the options in the "USER DEFINED VARIABLES" section
* Once you've edited the appropriate variables in that section and made a *"features file"*, run the script from terminal using the following command:
```
python visualize_IAV_segments.py
```
* Check the output to make sure everything looks correct.
* The PDF should be easily edited in a software like inkscape or adobe illustrator to make minor changes for visual preferences

## Generating a *"features file"*
The feature file should have this format:

|#segment | feature_start | feature_stop | featureID | feature_color | strand|
| --------|---------------|--------------|-----------|---------------|-------|
|NA | 379 | 533 | amplicon | blue | 1|
|PA | 35 | 426 | deletion | red | -1|
|PB1 | 1435 | nan | nan | nan | None|

- The only two required columns in this feature file is the first two columns, like this:
|#segment | feature_start |
| --------|---------------|
|NA | 379 |
|PA | 35 |
|PB1 | 1435 |

In this case, the script will generate a genome diagram with grey vertical lines at all of the locations specified in column 2

#### Notes:
- If you don't want to include any values in your final figure, such as labels for a subset of the features included in a table, replace the featureID with "nan"
  - If you simply leave that column out, it will create errors while attempting read in any columns that follow the value you skip (e.g. you can't specify a color for a feature without also adding adding place-holder "nan" values for feature_stop and featureID)
- all lines beginning with a hash character (#) are ignored

## Pick a color scheme
The following color_schemes are built in: 
* rainbow
* greyscale"
* grey(2-tone)
* blue(2-tone)
* red(2-tone)
* green(2-tone)
* grey
* blue
* red
* green
* *other* 
   - if you use "other", you also need to supply either a color name from python or a color library, or the color hex value such as "#b7da9c" to the variables *"other_primary_color"* and *"other_secondary_color"* below. Otherwise these values will be ignored

# How to install the dependencies of this script
If you haven't already, [download and install miniconda](https://docs.conda.io/en/latest/miniconda.html)

Once you have conda installed, run these commands in the appropriate terminal to numpy, biopython, and reportlab
```
conda install numpy
conda install -c conda-forge biopython
conda install -c anaconda reportlab
```
