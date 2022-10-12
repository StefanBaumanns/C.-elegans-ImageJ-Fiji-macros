// Macro for automatic measurment of ROS and MMP in the head regions of C. elegans nematodes using red light emitting fluorescent probes (e.g. Mitotracker Red or TMRE).
// RFP and TRANS Images are needed (e. g. at 10x magnification).
// The images of the channels must be sorted in individual folders and the mean background of the fluorescence channels has to be determined beforehand.
// Plugins: Biovoxxel, Morpholibj, Bio-Formats

#@ File (label = "TRANS input directory", style = "directory") input
#@ File (label = "RFP input directory", style = "directory") rfpInput
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".tif") suffix
#@ Integer (label="Mean background value RFP", min=0, max=1000) MeanBackgroundValueRFP
#@ String (label = "Group name", value = "") GroupName

// opening of RFP images as image stacks
run("Image Sequence...", "open=["+rfpInput+"] sort");
run("Subtract...", "value=MeanBackgroundValueRFP stack");

processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
	// saves results for all images in a single file
	run("Close All");
	selectWindow("Results");
	run("Summarize");
	saveAs("Results", output + "/" + GroupName + "_Results.tsv"); 
	selectWindow("Results");
	run("Close");
	selectWindow("ROI Manager");
	run("Close");
	}


function processFile(input, output, file) {

    // creation of a ROI based on TRANS images
	run("Bio-Formats", "open=[" + input + "/" + file +"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT"); // opening of TRANS images
	originalImage = getTitle();
	rename(replace(originalImage, "TRANS/", ""));
	renamedImage = getTitle();
	run("Duplicate...", "title=binary"); // duplicate original image and work on the copy
	
   	waitForUser ("Make a rough polygon selection around the head without inclusion of dark particles in the background. Draw a straight line at the terminal bulb in a right angle to the body. Press OK afterwards.");
	run("Select None");
	run("Gaussian Weighted Median", "radius=2");
	run("Morphological Filters", "operation=Erosion element=Disk radius=2");
	run("Enhance Contrast...", "saturated=0.01 normalize");
	run("Threshold...");
	setAutoThreshold("Huang");
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("Close");
	run("Fill Holes");
	run("EDM Binary Operations", "iterations=1 operation=erode");
	//run("EDM Binary Operations", "iterations=1 operation=open");
	run("Restore Selection");
	run("Clear Outside");
	run("Select None");
	run("EDM Binary Operations", "iterations=100 operation=close");
	run("Fill Holes");


	run("Analyze Particles...", "size=1000-100000 add"); // size exclusion to filter objects smaller or larger than the head
	selectWindow("binary");
	run("Close");
	selectWindow("binary-Erosion");
	run("Close");
	selectWindow(renamedImage);
	selectWindow(renamedImage);
	roiManager("Show All without labels");
	
	if (getBoolean("Is the segmentation result OK? If the extracted ROI is missing, close all open windows and restart the macro")) 	{
		roiManager("Deselect");
		roiManager("Save", output+ "/" + file + "_ROI.zip"); // saves Rois zip file

		// analysing RFP
		roiManager("Deselect");
		selectWindow("RFP");
		run("Set Measurements...", "area mean min integrated display redirect=None decimal=6");
		roiManager("Measure");
    }

   	else {
		roiManager("Deselect");
		roiManager("Delete");
  		waitForUser ("Make a precise polygon selection of the head. Draw a straight line at the terminal bulb in a right angle to the body. Press OK afterwards."); 
		roiManager("Add");
    	
    	roiManager("Deselect");
		roiManager("Save", output+ "/" + file + "_ROI.zip"); // saves Rois zip file

		// analysing RFP
		roiManager("Deselect");
		selectWindow("RFP");
		run("Set Measurements...", "area mean min integrated display redirect=None decimal=6");
		roiManager("Measure");
    } 
    
	selectWindow("RFP");
	run("Next Slice [>]");
	selectWindow(renamedImage);
	run("Close");	
	roiManager("Delete"); // clear ROI Manager for next image	
}






