/////////////////////////////////////////////////////////////////////////////
//////////////////// Lesion segmentation with WEKA //////////////////////////
//
// author: Niccolo` Bassetti - niccolo.bassetti@protonmail.com
//
// principal investigator: Nina Fatouros - nina.fatouros@wur.nl

//////////////////////////////////
// PREPARATION
//////////////////////////////////
//
// Input pictures to be analysed need to be placed in folder "Input_pic" 
// An empty folder called "Results" is needed to store the results.
//
// This pipeline will generate these results
// 1) "your_pic_lesion.jpg", a .jpeg copy of your image to be analysed
// 2) "your_pic_output.jpg", a .jpeg version of your image with a yellow border indicating the area measured
// 3) "your_image.RoiSet.zip", a .zip that can be loaded on Roi Manager to re-analysed the are measured.
//     This area can be adjusted (remove certain particles) and remeasured. 
// 4) "Results.txt", an overview of all the measurements done after each run of this script. 
//    NOTE: this file needs to be copied elsewehere, otherwise will be overwritten (you will loose your results!). 



//////////////////////////////////
// !!! IMPORTANT !!!
//////////////////////////////////
//
// 1) Empty "Results" folder before running this pipeline.
// 
// 2) Increase memory available for Fiji before starting
//   Change maimum/parallel values to allow Fiji to speed up image classificatio (CPU intensive process).
//   Click Edit > Options > Memory & Threads to select more memory


// TO START the pipeline, click Run (Ctrl + R)


//////////////////////////////////
// ANALYSIS
/////////////////////////////////

// This is the input directory
input = getDirectory("Select Input folder");

// This is the output directory
output = getDirectory("Select Results folder");

// Pipeline is repeated for every image of "Input" folder
list = getFileList(input);
for (i=0; i < list.length; i++) {		
	open(input+list[i]);

	// Resize image to speed up analysis, Step necessary only if it takes long to use WEKA 
	// CHANGE width & length based on your preference
	// BE AWARE smaller size = reduced quality!
	run("Size...", "width=1296 height=972 constrain average interpolation=Bilinear"); 

	// Set scale measurement
	setTool("rectangle");	// select area around scale bar
	waitForUser( "Pause","select area with scale bar, press OK when done");
	run("To Selection");	// zoom on selected area around scale bar 
	run("Select None");
	setTool("line");
	waitForUser( "Pause","Select scale bar of known distance. Enter known distance inside next window. Press OK to continue");
	run("Set Scale..."); // Enter the known distance of the selection
	run("Original Scale");	// Zoom out back to original scale
	
	// Rename and duplicate image.
	Name=getTitle();
	imgName=File.nameWithoutExtension;
	run("Duplicate...", "title=&imgName");

	// Reduce the image content to speed up the calculation time
	// Select area with lesion (leave out background with agar and leaf damage from borer)
	setTool("freehand");
	waitForUser("Pause","Select zone with lesion. If no lesion is visible, don`t select anything. Press OK when done."); // select area with HR
	if (selectionType() < 0) { 	// This returns a black image if no lesion is visible
		run("Select All");
		run("Clear", "slice");
		run("Select None");
		run("8-bit");
		setAutoThreshold("Default");
		setThreshold(85, 85);
		run("Analyze Particles...", " size=20-Infinity pixel summarize add slice"); //Measure thresholded area

	// Save pics results
		close(); // close copy
		path2 = output+File.nameWithoutExtension;
		saveAs("JPEG", path2 + "_lesion.jpeg");
		saveAs("JPEG", path2 + "_output.jpeg");
		close(); // close original
		roiManager("reset"); 		// clean ROI manager
	} else {						// If area with HR-like cell death is selected, this part prepares the selection for WEKA segmentation
		run("Make Inverse");
		run("Clear", "slice");
		run("Select None");
	
	// Image preprocessing. Can enhance segmentation results. Not used here. 
		//run("Enhance Contrast...", "saturated=0.3");
		//run("Remove Outliers...", "radius=2 threshold=50 which=Bright");
			
	// Calling WEKA with appropriate settings
		run("Trainable Weka Segmentation");
		selectWindow(Name);		// change line based on setup data analysis
		setLocation(2000, 200);		// this relocates the WEKA window for comfort. It can be skipped or change value to relocate window within your screen

	//selectWindow("Trainable Weka Segmentation v3.2.34"); 
		setTool("oval");		
		wait(1000);
/// This part selects the classifiers (called by WEKA "Feature") used to train the model. 
/// Change classifiers based on what works best on your pictures.
/// Write "false" to deselect or "true" to select a specific classifier. 
		// Deselect basic features
		waitForUser("Select appropriate features in Settings. This part can be automated in the pipeline (see script). Press OK to continue");
		call("trainableSegmentation.Weka_Segmentation.setFeature", "Gaussian_blur=false");
		call("trainableSegmentation.Weka_Segmentation.setFeature", "Sobel_filter=false");
		call("trainableSegmentation.Weka_Segmentation.setFeature", "Hessian=false");
		call("trainableSegmentation.Weka_Segmentation.setFeature", "Difference_of_gaussians=false");
		call("trainableSegmentation.Weka_Segmentation.setFeature", "Membrane_projections=true");
		//	Select additional features
		call("trainableSegmentation.Weka_Segmentation.setFeature", "Variance=true");
		call("trainableSegmentation.Weka_Segmentation.setFeature", "Mean=true");
		call("trainableSegmentation.Weka_Segmentation.setFeature", "Minimum=true");
		call("trainableSegmentation.Weka_Segmentation.setFeature", "Maximum=true");
		call("trainableSegmentation.Weka_Segmentation.setFeature", "Median=true");
		//call("trainableSegmentation.Weka_Segmentation.setFeature", "Bilateral=false");
		//call("trainableSegmentation.Weka_Segmentation.setClassBalance", "false");  
		call("trainableSegmentation.Weka_Segmentation.changeClassName", "0", "cell death");
		call("trainableSegmentation.Weka_Segmentation.changeClassName", "1", "healthy tissue");
		wait(1000);
		waitForUser( "Pauze", "Select representative cell death spots. Add selected spots to the list. Press OK when training is OK.");
		waitForUser( "Pauze", "Select representative healthy leaf tissue spots. Add selected spots to the list. Press OK when training is OK.");
		waitForUser( "Pauze", "Train classifier. Use overlay button to check classification.If unsufficient, add more training spots and repeat training. Advisable not to train more than 3-4 times to avoid overfitting. Press OK when training is OK.");
		//select real HR, no outliers
		selectWindow("Trainable Weka Segmentation v3.2.34"); // prevents error of weka was not selected
		setTool("freehand");
		waitForUser("Pause","Remove outliers and select real HR. Press OK when done."); // select area with HR
		
		// measure HR on classified image
		call("trainableSegmentation.Weka_Segmentation.getResult");
		selectWindow("Classified image");
		imgName=File.nameWithoutExtension;
		run("Rename...", "title=&imgName");
		run("Restore Selection");			
		run("8-bit"); 	//sets classified image to 8 bit for measuring
		setAutoThreshold("Default");
		setThreshold(85, 85);
		run("Analyze Particles...", " size=20-Infinity pixel summarize add slice"); //Measure thresholded area
		close(); // close classified image
		selectWindow("Trainable Weka Segmentation v3.2.34");
		run("Close"); // close WEKA & copy picture

		// Save pics results
		path2 = output+File.nameWithoutExtension;
		selectWindow(Name);
		//saveAs("PNG", path2 + "_lesion.png");	ignore saving original image to save space
		roiManager("Show All without labels");
		run("Flatten");
		saveAs("PNG", path2 + "_output.png");
		close(); // close flattened copy
		if (roiManager("count")!=0) {
			roiManager("Save", path2 + "_RoiSet.zip")
			roiManager("reset"); // clean ROI manager
			close(); // close original
		}else{close(); // close original 	
		}	//If function prevents error if no HR 
		run("Close All"); // close original avoid ask for saving modifications // TO CHECK IF IT WORKS
	}
}

selectWindow("Summary");
saveAs("Text", output + "Results.txt"); //Save measurement results as txt file
		
	
