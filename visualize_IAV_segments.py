'''
####################################################################################################################
IAV segment visualizer v0.1

WHAT THIS SCRIPT DOES
This script generates scale genome diagrams for Influenza A virus, adding the location of features to scale. One of 
the inputs for this script is a file that has the locations and desired labels, colors, etc (if any) of the features
you want to include. It has multiple color schemes built in, all picked to be red/green color blindness friendly.

HOW TO USE THIS SCRIPT
Create a feature file as described below and direct the script to that file by editing the fields in the "USER DEFINED VARIABLES"
section. Edit any other relevant fields as well before running. 

HOW TO INSTALL THE DEPENDANCIES OF THIS SCRIPT USING ANACONDA:
enter these commands in your terminal:
conda install numpy
conda install -c conda-forge biopython
conda install -c anaconda reportlab

If you have any problems, please contact the author of this script: Dave VanInsberghe at dvanins@emory.edu

#################################################   UPDATE NOTES   #################################################

'''
############################################## USER DEFINED VARIABLES ##############################################
output_name = "labeled_IAV_genome.pdf"
segment_alignment = "center" #left, center, or Right
page_width = 12 #The page width of your PDF output (in cm)
page_height = 12 #(in cm)
feature_annot_angle = 30 #The angle you want the annotations of your features to appear relative to the gene diagram (e.g. 0, 45, 90, or -90)
add_segment_labels = False #If True, the segments will be labled
add_feature_labels = False #If True, the featureID of the table will be written into the diagram. If False, no labels will be added. If a feature_ID is listed as "nan", no label will be added even if add_feature_labels == True
overwrite_existing_files = True

feature_filename = "features_to_add.txt"
'''
An example of how the feature_file should look is below:

#segment	feature_start	feature_stop	featureID	feature_color	strand
NA			379				533				amplicon	gainsboro		1
PA			35				426				deletion	red				-1
PB1			1435			nan				nan			nan				None

- This feature file can be as simple as two columns (segment and feature start site), but any values in this table that
you don't want to include in the figure need to be replaced with "nan" for the script to ignore it properly while also 
including information from any column to the right of that value you don't want to include
- all lines beginning with a hash character (#) are ignored
'''

color_scheme = "rainbow"
'''Available color_scheme values: 
- "rainbow"
- "greyscale"
- "grey(2-tone)"
- "blue(2-tone)"
- "red(2-tone)"
- "green(2-tone)"
- "grey"
- "blue"
- "red"
- "green"

 - "other" - if you use "other", you also need to supply either a color name from python or a color library, or the color hex value such as "#b7da9c" below
'''
#These are ignored unless color_scheme = "other"
other_primary_color = 'blue' #color of primary gene products (examples: #a4dbe0, "red", "dodgerblue")
other_secondary_color = 'lightblue' #color of all alternative gene products

############################################  END OF USER INPUT SECTION ############################################
import os
import sys
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib.units import cm
from random import randint
import numpy as np

project_dir = "./"

segment_length_dict = {'PB2':2329,'PB1':2323,'PA':2200,'HA':1750,'NP':1546,'NA':1459,'M':1031,'NS':887}
segment_order_dict = {0:'PB2',1:'PB1',2:'PA',3:'HA',4:'NP',5:'NA',6:'M',7:'NS'}
segment_CDS_location_dict = {'PB2':{'PB2':[(25,2304)]},
	'PB1':{'PB1':[(25,2298)],'PB1-F2':[(119,358)]},
	'PA':{'PA':[(25,2175)],'PA-X':[(25,594),(596,724)]},
	'HA':{'HA':[(25,1725)]},
	'NP':{'NP':[(25,1521)]},
	'NA':{'NA':[(25,1434)]},
	'M':{'M1':[(25,783)],'M2':[(25,50),(739,1006)]},
	'NS':{'NS1':[(25,684)],'NS2':[(25,54),(527,862)]}}
CDS_strand_dict = {'PB2':1, 'PB1':1, 'PB1-F2':-1, 'PA':1, 'PA-X':-1, 'HA':1, 'NP':1, 'NA':1, 'M1':1, 'M2':-1, 'NS1':1, 'NS2':-1}
max_stop = segment_length_dict['PB2']
primary_gene_products = ['PB2','PB1','PA','HA','NP','NA','M1','NS1']
alternative_gene_products = ['PB1-F2','PA-X','M2','NS2']

grey_2tone_color_dict = {'PB2':'#bbbdc0','PB1':'#bbbdc0','PB1-F2':'#e6e7e8','PA':'#bbbdc0','PA-X':'#e6e7e8','HA':'#bbbdc0','NP':'#bbbdc0','NA':'#bbbdc0','M1':'#bbbdc0','M2':'#e6e7e8','NS1':'#bbbdc0','NS2':'#e6e7e8'}
blue_2tone_color_dict = {'PB2':'#6aadc8','PB1':'#6aadc8','PB1-F2':'#a4dbe0','PA':'#6aadc8','PA-X':'#a4dbe0','HA':'#6aadc8','NP':'#6aadc8','NA':'#6aadc8','M1':'#6aadc8','M2':'#a4dbe0','NS1':'#6aadc8','NS2':'#a4dbe0'}
red_2tone_color_dict = {'PB2':'#c85147','PB1':'#c85147','PB1-F2':'#da8f7e','PA':'#c85147','PA-X':'#da8f7e','HA':'#c85147','NP':'#c85147','NA':'#c85147','M1':'#c85147','M2':'#da8f7e','NS1':'#c85147','NS2':'#da8f7e'}
green_2tone_color_dict = {'PB2':'#88b17d','PB1':'#88b17d','PB1-F2':'#b7da9c','PA':'#88b17d','PA-X':'#b7da9c','HA':'#88b17d','NP':'#88b17d','NA':'#88b17d','M1':'#88b17d','M2':'#b7da9c','NS1':'#88b17d','NS2':'#b7da9c'}

grey_color_dict = {'PB2':'#bbbdc0','PB1':'#bbbdc0','PB1-F2':'#bbbdc0','PA':'#bbbdc0','PA-X':'#bbbdc0','HA':'#bbbdc0','NP':'#bbbdc0','NA':'#bbbdc0','M1':'#bbbdc0','M2':'#bbbdc0','NS1':'#bbbdc0','NS2':'#bbbdc0'}
blue_color_dict = {'PB2':'#6aadc8','PB1':'#6aadc8','PB1-F2':'#6aadc8','PA':'#6aadc8','PA-X':'#6aadc8','HA':'#6aadc8','NP':'#6aadc8','NA':'#6aadc8','M1':'#6aadc8','M2':'#6aadc8','NS1':'#6aadc8','NS2':'#6aadc8'}
red_color_dict = {'PB2':'#c85147','PB1':'#c85147','PB1-F2':'#c85147','PA':'#c85147','PA-X':'#c85147','HA':'#c85147','NP':'#c85147','NA':'#c85147','M1':'#c85147','M2':'#c85147','NS1':'#c85147','NS2':'#c85147'}
green_color_dict = {'PB2':'#88b17d','PB1':'#88b17d','PB1-F2':'#88b17d','PA':'#88b17d','PA-X':'#88b17d','HA':'#88b17d','NP':'#88b17d','NA':'#88b17d','M1':'#88b17d','M2':'#88b17d','NS1':'#88b17d','NS2':'#88b17d'}

rainbow_color_dict = {'PB2':'#be1e2d','PB1':'#ef4036','PB1-F2':'#f37b60','PA':'#f7931d','PA-X':'#fbaf3f','HA':'#f1ec6e','NP':'#cfe0af','NA':'#6fae69','M1':'#3d63ab','M2':'#6688c0','NS1':'#a378a0','NS2':'#ccbcd5'}
greyscale_color_dict = {'PB2':'#58585b','PB1':'#656669','PB1-F2':'#727376','PA':'#7e8082','PA-X':'#8a8c8f','HA':'#96989b','NP':'#a3a5a8','NA':'#b0b2b4','M1':'#bdbfc1','M2':'#cbcdce','NS1':'#d8dadb','NS2':'#e6e7e8'}

other_color_dict = {}
if color_scheme == "other":
	for segment_name in segment_CDS_location_dict:
		for CDS_product in segment_CDS_location_dict[segment_name]:
			if CDS_product in primary_gene_products:
				other_color_dict[CDS_product] = other_primary_color
			else:
				other_color_dict[CDS_product] = other_secondary_color

color_options = {"rainbow":rainbow_color_dict,"grey(2-tone)":grey_2tone_color_dict,"blue(2-tone)":blue_2tone_color_dict,"red(2-tone)":red_2tone_color_dict,"green(2-tone)":green_2tone_color_dict,"grey":grey_color_dict,"blue":blue_color_dict,"red":red_color_dict,"green":green_color_dict,"greyscale":greyscale_color_dict,"other":other_color_dict}
CDS_color_dict = color_options[color_scheme]

offset_dict = {}
for segment_name in segment_length_dict:
	if max_stop-segment_length_dict[segment_name] >0 or segment_alignment == "left":
		if segment_alignment == "center":
			offset = int((max_stop-segment_length_dict[segment_name])/2)
		elif segment_alignment == "right":
			offset = int((max_stop-segment_length_dict[segment_name]))
	else:
		offset = 0
	offset_dict[segment_name] = offset

####################################################################################################################
#Load in feature file
feature_dictionary = {}
infile = open(project_dir+feature_filename,"r")
for line in infile:
	line = line.strip()
	if len(line)>0:
		if line[0]!="#":
			line = line.split("\t")
			if len(line) >=2:
				segment_name = line[0]
				feature_start = line[1]
				try:
					feature_start = int(feature_start)
				except:
					sys.exit('"feature_start" variable is not integer: "'+feature_start+'"')
			if len(line)>=3:
				feature_stop = line[2]
				if feature_stop == "nan":
					feature_stop = feature_start
				else:
					try:
						feature_stop = int(feature_stop)
					except:
						sys.exit('"feature_stop" variable is not integer: "'+feature_stop+'"')
					if feature_stop < feature_start:
						sys.exit('"feature_stop" value is smaller than "feature_start". Exiting.')
			else:
				feature_stop = feature_start

			if len(line)>=4:
				feature_ID = line[3]
			else:
				feature_ID = "nan"

			if len(line)>=5:
				feature_color = line[4]
				if feature_color == "nan":
					feature_color = 'grey'
			else:
				feature_color = 'grey'

			if len(line)>=6:
				strand_val = line[5]
				if strand_val == "None":
					strand_val = None
				else:
					try:
						strand_val = int(strand_val)
					except:
						strand_val = None
			else:
				strand_val = None

			tup = (feature_ID,feature_start,feature_stop,feature_color,strand_val)
			try:
				feature_dictionary[segment_name].append(tup)
			except:
				feature_dictionary[segment_name] = []
				feature_dictionary[segment_name].append(tup)
infile.close()


gdd = GenomeDiagram.Diagram('Genome Diagram')
for num in range(0,8):
	segment_name = segment_order_dict[num]
	print(segment_name)
	segment_offset = offset_dict[segment_name]
	segment_start = 0+segment_offset
	segment_stop = segment_length_dict[segment_name]+segment_offset
	gdt_features = gdd.new_track(track_level=1, name=segment_name,start=segment_start,end=segment_stop, greytrack=False,scale_largetick_interval=1000,scale_largeticks=0.25,scale_smalltick_interval=100,scale_smallticks=0.12,scale_largetick_labels=0, scale_smalltick_labels=0)
	for CDS_product in segment_CDS_location_dict[segment_name]:
		gds_features = gdt_features.new_set()
		for loc_tup in segment_CDS_location_dict[segment_name][CDS_product]:
			CDS_start = loc_tup[0]+segment_offset
			CDS_stop = loc_tup[1]+segment_offset
			CDS_strand = CDS_strand_dict[CDS_product]
			label_angle_phase = 0
			if CDS_strand == -1:
				label_angle_phase = 180
			feature = SeqFeature(FeatureLocation(CDS_start, CDS_stop), strand=CDS_strand)
			gds_features.add_feature(feature, sigil="BOX",color=CDS_color_dict[CDS_product],label=add_segment_labels, label_angle=(0-label_angle_phase),label_position="start",label_size=10,name=CDS_product,border=None)
	try:
		feature_list = feature_dictionary[segment_name]
	except:
		feature_list = []
	for feature_tup in feature_list:
		feature_ID = feature_tup[0]
		feature_start = feature_tup[1]+segment_offset
		feature_stop = feature_tup[2]+segment_offset
		feature_color = feature_tup[3]
		feature_strand = feature_tup[4]
		label_angle_phase = 0
		if feature_strand == -1:
			label_angle_phase = 180
		feature = SeqFeature(FeatureLocation(feature_start, feature_stop), strand=feature_strand)
		gds_features.add_feature(feature, sigil="BOX",color=feature_color,label=add_feature_labels, label_angle=(feature_annot_angle-label_angle_phase),name=feature_ID,border=None)
gdd.draw(format='linear', pagesize=(page_width*cm,page_height*cm), fragments=1, start=0, end=max_stop,track_size=0.8,xl=0.15, xr=0.15, yt=0, yb=0, tracklines=False)
write_output_file = False
if os.path.isfile(project_dir+output_name) == True and overwrite_existing_files == False:
	user_input = input('The specified file already exists, overwrite? (yes or no)\n')
	if user_input == "yes" or user_input == "Yes" or user_input == "y" or user_input == "YES" or user_input == "Y":
		write_output_file = True
else:
	write_output_file = True
if write_output_file == True:
	gdd.write(project_dir+output_name, "pdf")
