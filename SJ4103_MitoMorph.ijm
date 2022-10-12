// Macro for automatic measurment of mitochondrial morphology using C. elegans SJ4103, expressing GFP in mitochondria of body wall muscle cells (for alternative tagging of mitchondria the macro may be adjusted)
// The input images need a high spatial resolution and a good SNR (e.g. through a high NA objective, deconvolution, matching refractive indices) and should be scaled for correct measurement of mitochondrial size
// Plugins: adjustable watershed, Biovoxxel, Morpholibj, Bio-Formats


#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".tif") suffix
#@ String (label = "Group name", value = "") GroupName


// Creation of a summary

measurement = newArray("Total Cell Area [square µm]",
		           	   "Mitochondrial Coverage [percent]",
		               "Mitochondria [count]",
		               "Mitochondria Per Area [count/square µm]",
		               "Mean Size [square µm]",
		               "Mean Aspect Ratio [major- /minor axis]",
		               "Mean Circularity [4pi*area/perimeter^2]",
		               "Mean Form Factor [perimeter^2/4pi*area]");
Table.create("Summary");
Table.setColumn(" ", measurement);


processFolder(input);

// Function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
		processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
		processFile(input, output, list[i]);
}

}


function processFile(input, output, file) {
	run("Bio-Formats", "open=[" + input + "/" + file +"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	rename(file);
	originalName = getTitle(); // Get original image name

	// Creation of a Z-max projection
	run("Z Project...", "projection=[Max Intensity]");
	selectWindow(originalName);
	run("Close");

	


	// Pre-processing
	run("Enhance Contrast...", "saturated=0.1 normalize");
	setOption("ScaleConversions", true);
	run("8-bit");
	run("Bilateral Filter", "spatial=3 range=20"); // The range may be adjusted to a lower value for images with poor contrast
	run("Bilateral Filter", "spatial=3 range=20");
	run("Subtract Background...", "rolling=15"); // The radius (rolling) may be adjusted based on the size of the largest mitochondria
	run("Unsharp Mask...", "radius=2 mask=0.60");
	run("Subtract Background...", "rolling=15"); 
	
	rename(originalName);
	
	run("Duplicate...", "title=Binary");
	waitForUser ("Make a polygon selection of a single muscle cell without slicing through mitochondria. Regions near the edge of the image or unsharp regions should be excluded. Press OK afterwards.");

	run("Enhance Contrast...", "saturated=1.0 normalize");
	run("Select None");

	// Extraction of objects
	run("Auto Local Threshold", "method=Phansalkar radius=35 parameter_1=0 parameter_2=0 white");


	run("Restore Selection");
	run("Clear Outside");

	// Post-processing
	run("Fill Holes");
	run("Adjustable Watershed", "tolerance=1.5"); // The tolerance may be adjusted to a higher value or the command may be excluded to avoid oversegmentation of large elongated mitochondria

	
	// Measurement of cell size based on closing of the extracted objects
	run("Select None");
	run("Duplicate...", " ");
	run("EDM Binary Operations", "iterations=100 operation=close"); 
	run("Set Measurements...", "area perimeter shape area_fraction display redirect=None decimal=6");
	run("Analyze Particles...", "size=100-10000 display add clear"); 
	TotalCellArea = getResult("Area", nResults()-1);

	// Measurment of mitochondrial coverage
	selectWindow("Binary");
	run("Clear Results");
	roiManager("Measure");
	MitochondrialCoverage = getResult("%Area", nResults()-1);

	roiManager("Deselect"); 
	roiManager("Save", output+ "/" + originalName + "_ROI_MuscleCell.zip"); // Saves the ROI of the muscle cell
	roiManager("Delete");

	// Measurement of objects
	run("Analyze Particles...", "size=4-Infinity pixel add clear display"); // The size exclusion can be removed
	roiManager("Measure");

	

n = nResults;
	Circ = newArray(n);
	Perim = newArray(n);
	for (i=0; i<n; i++) {
      Circ[i] = getResult("Circ.", i);
      Perim[i]= getResult("Perim.", i);
    }
	for (i=0; i<n; i++) {
      setResult("FF", i, 1/Circ[i]);
      setResult("AWFF", i, (Perim[i] * Perim[i]) / (4 *3.14159265358979));
    }
    updateResults();
	run("Summarize");


	// Read out values from results table
	selectWindow("Results"); rows = nResults;
	ObjectCounts = newArray(rows);
	Objects = (parseFloat(CountForMe(ObjectCounts)));
	ObjectsPerArea = Objects / TotalCellArea;
	MeanArea = getResult("Area", nResults() - 4);
	MeanAR = getResult("AR", nResults() - 4);
	MeanCircularity = getResult("Circ.", nResults() - 4);
	Label = getResult("Label", nResults() - 5);
	MeanFormFactor = getResult("FF", nResults () -4);
	
			     
	Numbers = newArray( TotalCellArea,
                        MitochondrialCoverage,
				        Objects,
			 	        ObjectsPerArea,
				        MeanArea,
				        MeanAR,
				        MeanCircularity,
				        MeanFormFactor);

				   


// Control of the result through the user
	selectWindow(originalName);
	roiManager("Show All without labels");
	waitForUser ("Check the segmentation result and press OK afterwards. If the result is good, also press OK in the next pop-up.");
if (getBoolean("Is the segmentation result OK?")) 	{
	selectWindow(originalName);
	roiManager("Show None");
	save(output  + "/" + originalName + "_Z-max"); // save Z-max projection
	selectWindow("Binary");
	roiManager("Show None");
	save(output  + "/" + originalName + "_Binary" ); // Save current binary image
	roiManager("Deselect");
	roiManager("Save", output+ "/" + originalName + "_ROI_Mitochondria.zip"); // Saves the ROIs of mitochondria

	label = originalName;
	selectWindow("Summary");
	Table.setColumn(label, Numbers);
 	
}
else {
	File.delete(output+ "/" + originalName + "_ROI_MuscleCell.zip");
}
	// clean up
	roiManager("Deselect");
	roiManager("Delete"); 
	selectWindow("Results");
	run("Close");
	selectWindow("ROI Manager");
	run("Close");

	selectWindow(originalName);
	close("\\Others");
	

while (getBoolean("Do you want to analyse the muscle cell again or another muscle cell in the same picture?")) 	{

	// Second analysis
	run("Duplicate...", " ");
	selectWindow(originalName);
	run("Close");
    originalName = getTitle();

	run("Duplicate...", "title=Binary");
	// Pre-processing


	waitForUser ("Make a polygon selection of a single muscle cell without slicing through mitochondria. Regions near the edge of the image or unsharp regions should be excluded. Press OK afterwards.");
	run("Enhance Contrast...", "saturated=1.0 normalize");
	run("Select None");

	// Extraction of objects
	run("Auto Local Threshold", "method=Phansalkar radius=35 parameter_1=0 parameter_2=0 white");

	run("Restore Selection");
	run("Clear Outside");

	// Post-processing
	run("Fill Holes");
	//run("Adjustable Watershed", "tolerance=1.5"); // The tolerance may be adjusted to a higher value to avoid oversegmentation of large elongated mitochondria

	
	// Measurement of cell size based on closing of the extracted objects
	run("Select None");
	run("Duplicate...", " ");
	run("EDM Binary Operations", "iterations=100 operation=close"); 
	run("Set Measurements...", "area perimeter shape area_fraction display redirect=None decimal=6");
	run("Analyze Particles...", "size=100-10000 display add clear"); 
	TotalCellArea = getResult("Area", nResults()-1);

	// Measurment of mitochondrial coverage
	selectWindow("Binary");
	run("Clear Results");
	roiManager("Measure");
	MitochondrialCoverage = getResult("%Area", nResults()-1);

	roiManager("Deselect"); 
	roiManager("Save", output+ "/" + originalName + "_ROI_MuscleCell.zip"); // Saves the ROI of the muscle cell
	roiManager("Delete");

	// Measurement of objects
	run("Analyze Particles...", "size=4-Infinity pixel add clear display"); // The size exclusion can be removed
	roiManager("Measure");

	


	n = nResults;
	Circ = newArray(n);
	Perim = newArray(n);
	for (i=0; i<n; i++) {
      Circ[i] = getResult("Circ.", i);
      Perim[i]= getResult("Perim.", i);
    }
	for (i=0; i<n; i++) {
      setResult("FF", i, 1/Circ[i]);
    }
    updateResults();
	run("Summarize");


	// Read out values from results table
	selectWindow("Results"); rows = nResults;
	ObjectCounts = newArray(rows);
	Objects = (parseFloat(CountForMe(ObjectCounts)));
	ObjectsPerArea = Objects / TotalCellArea;
	MeanArea = getResult("Area", nResults() - 4);
	MeanAR = getResult("AR", nResults() - 4);
	MeanCircularity = getResult("Circ.", nResults() - 4);
	Label = getResult("Label", nResults() - 5);
	MeanFormFactor = getResult("FF", nResults () -4);
	

	Numbers = newArray(	TotalCellArea,
                        MitochondrialCoverage,
				        Objects,
			 	        ObjectsPerArea,
				        MeanArea,
				        MeanAR,
				        MeanCircularity,
				        MeanFormFactor);
       

				   


// Control of the result through the user
	selectWindow(originalName);
	roiManager("Show All without labels");
	waitForUser ("Check the segmentation result and press OK afterwards. If the result is good, also press OK in the next pop-up.");
if (getBoolean("Is the segmentation result OK?")) 	{
	selectWindow(originalName);
	roiManager("Show None");
	save(output  + "/" + originalName + "_Z-max"); // save Z-max projection
	selectWindow("Binary");
	roiManager("Show None");
	save(output  + "/" + originalName + "_Binary"); // Save current binary image

	roiManager("Deselect");
	roiManager("Save", output+ "/" + originalName + "_ROI_Mitochondria.zip"); // Saves the ROIs of mitochondria

	label = originalName;
	selectWindow("Summary");
	Table.setColumn(label, Numbers);

}
else {
	File.delete(output+ "/" + originalName + "_ROI_MuscleCell.zip");
}

	// clean up
	roiManager("Deselect");
	roiManager("Delete"); 
	selectWindow("Results");
	run("Close");
	selectWindow("ROI Manager");
	run("Close");

	selectWindow(originalName);
	close("\\Others");

}
 run("Close All");

}

// Saves the summary

selectWindow("Summary");
dir = getDirectory("Choose a directory for the summary.");
saveAs("Results", dir + "Summary_" + GroupName + ".tsv"); 
run("Close");

if (isOpen("Log")) {
 selectWindow("Log");
 run("Close");
}

// Maths
function CountForMe(data) {
	entries = data.length;
	total = 0.0;
	for (i=0; i<entries-4; i++) 
		{
			total = total + 1;
		}
		return(total);
}
