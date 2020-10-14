/*  PiQSARS : Ratiometric probes monitoring 
 *  Written by Elise Lévy, Université Paris-Saclay, UVSQ, INRAE, VIM, 78350, Jouy-en-Josas, France.    
 *  Use and distribution of this macro is free for academic purposes only, not for commercial use.
 *  
 *  For more information see publication: 
 *  Lévy, E., Jaffrézic, F., Laloë, D., Rezaei, H., Huang, M.-E., Béringue, V., Martin, D., and Vernis, L. (2020). 
 *  PiQSARS: A pipeline for quantitative and statistical analyses of ratiometric fluorescent biosensors. 
 *  MethodsX 7, 101034.
 *  
 *  Contact address: elise.levy@agroparistech.fr, davy.martin@inrae.fr, laurence.vernis@cnrs.fr                                               
 *  
 *  
 *  Note 1 : The acquisition must be a hyperstack with 2 channels corresponding to both acquisitions.
 *  Note 2 : Before running the macro, it is recommended that you process your negative control image in the same way as the image to analyse (background subtraction, filtering ...).
 *           You will thus be able to use this image as a negative control.
 *  Note 3 : Prefer using only non-accentuated letters, numerals and underscores for images names and paths.
*/


// --- Default parameters of some dialogboxes

DefaultRatio="C2/C1 (e.g. HyPer)";		// alternatively : "C2/C1 (e.g. HyPer)" or "C1/C2 (e.g. roGFP)"
DefaultBackgroundSubtraction="yes";		// alternatively : "yes" or "no"
DefaultBinning="yes";					// alternatively : "yes" or "no"
DefaultMed="no";						// alternatively : "yes" or "no"
DefaultGauss="no";						// alternatively : "yes" or "no"
DefaultMedianRadius=2;					// alternatively : any other number
DefaultGaussianRadius=4;				// alternatively : any other number
DefaultBinaryImage="Sum";				// alternatively : "Sum", "C1" or "C2"
DefaultIncludeEdges="yes";				// alternatively : "yes" or "no"

// --- Size and citcularity parameters for cells segmentation

MinCellSize=100;						// alternatively : any other number (your cells' minimum size in pixels² or µm² depending if your images are calibrated or not)
MaxCellSize="Infinity";					// alternatively : "Infinity" or any number (your cells' maximum size in pixels² or µm² depending if your images are calibrated or not)
MinCellCircularity=0;					// alternatively : any other number between 0 and 1 (the highest is the number you choose, the more circular the filtered particles will be)
MaxCellCircularity=1;					// alternatively : any other number between 0 and 1 (the lowest is the number you choose, the less circular the filtered particles will be)
IncludeHoles="";						// alternatively : "include " (if you want to include holes) (notice the space at the end of the string).

// --- Name of the report containing all the user's choices

ReportName="AnalysisReport"				// alternatively : any character string

// --- Images opening and processing - ratio calculation - user's choice saving

run("Misc...", "divide=NaN run"); // The result of a division by zero will be NaN
EmptyWorkspace();
waitForUser("Chose the folder in which the control data will be saved");
directory=getDirectory("Results folder");
res=OpenSplitDisplay();
	nbSlices=res.length-1;
	times=Array.deleteIndex(res,nbSlices);
	ImageName=res[nbSlices];
waitForUser("Observe your images and decide if you want \nto apply an additive binning to increase their dynamics. \nCaution: it will decrease their resolution. \nAfterwards, click on OK");
res = Dialogbox1(DefaultRatio,DefaultBackgroundSubtraction,DefaultBinning,DefaultMed,DefaultGauss);
	CalculatedRatio=res[0];
	SubtractBG=res[1];
	AdditiveBinning=res[2];
	Med=res[3];
	Gauss=res[4];
Binning("C1",AdditiveBinning,directory);
Binning("C2",AdditiveBinning,directory);
BackgroundSubtraction(SubtractBG,directory,nbSlices);
Radii=FilterRadius(Med,Gauss,DefaultMedianRadius,DefaultGaussianRadius);
Filtering("C1",Med,Gauss,Radii,directory);
Filtering("C2",Med,Gauss,Radii,directory);
nb_ROI=LetTheUserThink();
RatioCalculation(CalculatedRatio,directory);
UsersChoiceSaving(ImageName,CalculatedRatio,SubtractBG,AdditiveBinning,Med,Gauss,Radii,nb_ROI,MinCellSize,MaxCellSize,MinCellCircularity,MaxCellCircularity,IncludeHoles,directory,ReportName);

// --- Cells segmentation and intensities measurments - user's choice completing

SelectDuplicate(nb_ROI,directory);
i_C1=newArray(0); i_C2=newArray(0); ratio=newArray(ratio); 
for(i=0;i<nb_ROI;i++){
	ThresholdChoiceC1=Segmentation("C1",directory,i);
	ThresholdChoiceC2=Segmentation("C2",directory,i);
	ParmsMask=DialogboxMask(DefaultBinaryImage,DefaultIncludeEdges);
		//Mask_choice=ParmsMask[0];
		//Include_edges=ParmsMask[1];
	UsersChoiceCompleting(i,ThresholdChoiceC1[0],ThresholdChoiceC1[1],ThresholdChoiceC2[0],ThresholdChoiceC2[1],ParmsMask[0],ParmsMask[1],directory,ReportName);
	MaskGeneration(ParmsMask[0],ParmsMask[1],i,MinCellSize,MaxCellSize,MinCellCircularity,MaxCellCircularity,IncludeHoles,directory);
	IntensitiesC1=IntensityMeasurment("C1",i,nbSlices); 
		i_C1=Array.concat(i_C1,IntensitiesC1);
	IntensitiesC2=IntensityMeasurment("C2",i,nbSlices); 
		i_C2=Array.concat(i_C2,IntensitiesC2);
	IntensitiesRatio=IntensityMeasurment("ratio",i,nbSlices); 
		ratio=Array.concat(ratio,IntensitiesRatio);
	CloseImagesAndClean(i);
}
ResultsSaving(nb_ROI,times,i_C1,i_C2,ratio);
EmptyWorkspace();


// --------------------------------------------------------

function EmptyWorkspace(){
	
	run("Close All");
	roiManager("reset");
  	Table.deleteRows(0, nResults-1);	
}

// ---------------------------------------------------------

function OpenSplitDisplay(){
	waitForUser("Open the raw image you want to analyze");
	open(); 
	imageName=getTitle();
	run("Split Channels");
	selectWindow("C1-"+imageName); 
	rename("C1.tif");
	nbSlices=nSlices;
	run("HiLo"); // HiLo LUT --> easier detection of saturated cells
	run("Enhance Contrast", "saturated=0.35");	
	selectWindow("C2-"+imageName); 
	rename("C2.tif");
	run("HiLo"); // HiLo LUT --> easier detection of saturated cells
	run("Enhance Contrast", "saturated=0.35");	
	if(nbSlices==1){
		time_vect=newArray(1);
	} else {
		run("Plot Z-axis Profile");
		Plot.showValues(); 
		vectHeaders=split(String.getResultsHeadings);
		header=vectHeaders[0]; // header of the 1st column of the 'Results' table (usually "X" or "X0")
		time_vect=newArray(nbSlices+1);
		for(j=0;j<nbSlices;j++){
	    	time_vect[j]=getResult(header,j)/60;
		}
		Table.deleteRows(0, nResults-1); 
		time_vect[nbSlices]=imageName;
	}
	return time_vect;
}

// ----------------------------------------------------------

function Dialogbox1(DefaultRatio,DefaultBackgroundSubtraction,DefaultBinning,DefaultMed,DefaultGauss){

	Dialog.create("Information about the image treatment");
	Dialog.addRadioButtonGroup("Which ratio do you want to calculate?", newArray("C2/C1 (e.g. HyPer)", "C1/C2 (e.g. roGFP)") , 1, 2, DefaultRatio);
	Dialog.addRadioButtonGroup("Do you want to use a 2*2 additive binning\nto increase the signal intensity? \nCaution: it will decrease the resolution.",newArray("yes","no"),1,2,DefaultBinning);
	Dialog.addRadioButtonGroup("Do you want to subtract \nthe background mean value from the images?", newArray("yes", "no"), 1, 2, DefaultBackgroundSubtraction);
	Dialog.addRadioButtonGroup("Do you want to apply a median filter \nto reduce the salt & pepper noise?", newArray("yes", "no"), 1, 2, DefaultMed);
	Dialog.addRadioButtonGroup("Do you want to apply a gaussian filter \nto smooth the signal and make the segmentation easier?", newArray("yes", "no"), 1, 2, DefaultGauss);
	Dialog.show; 

	return newArray(Dialog.getRadioButton,Dialog.getRadioButton,Dialog.getRadioButton,Dialog.getRadioButton,Dialog.getRadioButton); // Gets the answers to the 5 questions in the right order
}


// -----------------------------------------------------------

function Binning(ImageName,AdditiveBinning,res_dir){
	
	if(AdditiveBinning=="yes"){
		selectWindow(ImageName+".tif"); 
		run("Bin...", "x=2 y=2 z=1 bin=Sum");
		saveAs("Tiff",res_dir+ImageName+"_binning.tif"); 
		rename(ImageName+".tif");
	}
}

// -----------------------------------------------------------

function BackgroundSubtraction(SubtractBG,res_dir,nb_slices){

	if(SubtractBG=="yes"){
		Dialog.create("Acquisition choice"); items = newArray("C1", "C2"); 
		Dialog.addRadioButtonGroup("On which channel do you want to \ndraw a background-containing rectangle ?", items, 1, 2, "C1"); 
		Dialog.show;
		Acquisition_bg=Dialog.getRadioButton;
		if(Acquisition_bg=="C1"){
			selectWindow("C1.tif"); 
			} else{
			selectWindow("C2.tif"); 	
		}

		//  The user draws a rectangle on the chosen acquisition, which is reported on the other one.
		setTool("rectangle");
		waitForUser ("Background definition", "Draw a no-cell-containing rectangle to define the background \non the chosen acquisition.");
		if(Acquisition_bg=="C1"){
			selectWindow("C2.tif");
			run("Restore Selection");
		}else{
			selectWindow("C1.tif");
			run("Restore Selection");
		}
		roiManager("Add");
		roiManager("Save",directory+"background_rectangle.zip");
		run("Set Measurements...", "mean redirect=None decimal=3");
		
		// Background measurment and subtraction on the C1 acquisition
		selectWindow("C1.tif");
		for (i=0;i<nb_slices;i++){
			setSlice(i+1); 
			run("Measure");
		}
		saveAs ("Results", res_dir+"C1_background.csv");
		selectWindow("C1.tif");
		run("Select All");
		for (i=0;i<nb_slices;i++){
			meanbg=getResult('Mean',i); 
			setSlice(i+1); 
			run("Subtract...", "value="+meanbg);
		}
		selectWindow("C1.tif"); 
		saveAs("Tiff",res_dir+"C1_background_subtraction.tif"); 
		rename("C1.tif");
		Table.deleteRows(0, nResults-1);
		
		// Background measurment and subtraction on the C2 acquisition
		selectWindow("C2.tif");
		for (i=0;i<nb_slices;i++){
			setSlice(i+1);
			run("Measure");
		}
		saveAs ("Results", directory+"C2_background.csv");
		selectWindow("C2.tif");run("Select All");
		for (i=0;i<nb_slices;i++){
			meanbg=getResult('Mean',i);
			setSlice(i+1);
			run("Subtract...", "value="+meanbg);
		}
		selectWindow("C2.tif"); 
		saveAs("Tiff",res_dir+"C2_background_subtraction.tif"); 
		rename("C2.tif");
		Table.deleteRows(0, nResults-1);
	}
	roiManager("reset");

}

// -----------------------------------------------------------

function FilterRadius(Med,Gauss,DefaultMed,DefaultGauss){
	
	if(Med=="no"){
		radmed=newArray(0);
	}
	if(Gauss=="no"){
		radgauss=newArray(0);
	}
		
	if(Med=="yes"||Gauss=="yes"){
		Dialog.create("Filter(s) Radius(es)");
		if(Med=="yes"){
			Dialog.addNumber("Median filter radius", DefaultMed);
		}
		if(Gauss=="yes"){
			Dialog.addNumber("Gaussian filter radius", DefaultGauss);
		}
		Dialog.show();
		if(Med=="yes"){
			radmed=Dialog.getNumber();
			if(Gauss=="yes"){
				radgauss=Dialog.getNumber();
			}
		}else{
			if(Gauss=="yes"){
				radgauss=Dialog.getNumber();
			}
		}
	}

	return Array.concat(radmed,radgauss);
}

// -----------------------------------------------------------

function Filtering(ImageName,Med,Gauss,radiuses,res_dir){
	
	if(Med=="yes"){
		selectWindow(ImageName+".tif");
		run("Median...", "radius="+radiuses[0]+" stack"); 
		saveAs("Tiff",res_dir+ImageName+"_med.tif"); 
		rename (ImageName+".tif");
		if(Gauss=="yes"){
			selectWindow(ImageName+".tif"); 
			run("Duplicate...", "duplicate");
			run("Gaussian Blur...", "sigma="+radiuses[1]+" stack"); 
			saveAs ("Tiff", res_dir+ImageName+"_gauss.tif"); 
			rename(ImageName+"segm.tif");
		}else{
			selectWindow(ImageName+".tif");
			run("Duplicate...", "duplicate");
			rename(ImageName+"segm.tif");
		}
	} else {
		if(Gauss=="yes"){
			selectWindow(ImageName+".tif"); 
			run("Duplicate...", "duplicate");
			run("Gaussian Blur...", "sigma="+radiuses[0]+" stack"); 
			saveAs ("Tiff", res_dir+ImageName+"_gauss.tif"); 
			rename(ImageName+"segm.tif");
		}else{
			selectWindow(ImageName+".tif");
			run("Duplicate...", "duplicate");
			rename(ImageName+"segm.tif");
		}
	}
}

// ----------------------------------------------------------

function LetTheUserThink(){
	
	selectWindow("C2segm.tif"); 
	run("Brightness/Contrast...");
	waitForUser ("nb cells", "Please count the cells you wish to analyze. \nThen, click on OK");
	Dialog.create("Nb ROI"); 
	Dialog.addNumber("Number of cells to be analyzed", 1); Dialog.show();
	
	return Dialog.getNumber();
}

// ------------------------------------------------------------
function RatioCalculation(probe,res_dir){

	if(probe=="C2/C1 (e.g. HyPer)"){
		imageCalculator("Divide create 32-bit stack", "C2.tif","C1.tif"); 
		selectWindow("Result of C2.tif");
		saveAs ("Tiff", res_dir+"ratio.tif");
	}else{
		imageCalculator("Divide create 32-bit stack", "C1.tif","C2.tif"); 
		selectWindow("Result of C1.tif");
		saveAs ("Tiff", res_dir+"ratio.tif");
	}
}

// ------------------------------------------------------------

function SelectDuplicate(nb_ROI,res_dir){

	// --- Rectangles selection
	
	Dialog.create("Choice of the image on which to define the rectangles");                     
	items = newArray("C1", "C2");
	Dialog.addRadioButtonGroup("On which channel do you want to \ndefine the rectangles around the cells?", items, 1, 2, "C1");
	Dialog.show;
	
	if(Dialog.getRadioButton=="C1"){
		
		for(i=0;i<nb_ROI;i++){ // number of the analysed cell = i+1
			selectWindow("C1segm.tif"); 
			setTool("rectangle");
			waitForUser ("Rectangle "+i+1+" selection", "Draw a rectangle containing the cell "+i+1+ " on every slice of the stack \nif possible so that the cell does not touch the edge of the rectangle. \nThen, click on OK.");
			roiManager("Add"); 
			roiManager("Show All with labels"); 
		}
	} else {
		for(i=0;i<nb_ROI;i++){ // n° de la cellule analysée = i+1
			selectWindow("C2segm.tif"); setTool("rectangle");
			waitForUser ("Rectangle "+i+1+" selection", "Draw a rectangle containing the cell "+i+1+ " on every slice of the stack \nif possible so that the cell does not touch the edge of the rectangle. \nThen, click on OK");
			roiManager("Add");
			roiManager("Show All with labels");
		}
	}
	
	// --- Rectangles duplication
	
	for(i=0;i<nb_ROI;i++){
		
		// --- ROI duplication on the C1 acquisition (smoothed or not)
		
		selectWindow ("C1segm.tif");
		roiManager ("Select",i);
		run ("Duplicate...", "duplicate"); 
		rename(i+1+"_ROI_C1_segm.tif"); 
		
		// --- ROI duplication on the C2 acquisition (smoothed or not)
		
		selectWindow ("C2segm.tif"); 
		roiManager ("Select",i);
		run ("Duplicate...", "duplicate"); 
		rename(i+1+"_ROI_C2_segm.tif"); 
		
		// --- ROI duplication on the C1 acquisition (denoised or not)
		
		selectWindow ("C1.tif"); 
		roiManager ("Select",i); 
		run ("Duplicate...", "duplicate"); 
		rename(i+1+"_ROI_C1.tif");
		
		// --- ROI duplication on the C2 acquisition (denoised or not)
		
		selectWindow ("C2.tif"); 
		roiManager ("Select",i); 
		run ("Duplicate...", "duplicate"); 
		rename(i+1+"_ROI_C2.tif");
		
		// --- ROI duplication on the ratio image
		
		selectWindow ("ratio.tif"); 
		roiManager ("Select",i); 
		run ("Duplicate...", "duplicate");
		rename(i+1+"_ROI_ratio.tif");
	}
	
	// --- Saving the rectangles' coordinates
	roiManager("Save", res_dir+"Analyzed cells.zip"); 
}

// -------------------------------------------------------------

function Segmentation (imagename,res_dir,i){

	selectWindow(i+1+"_ROI_"+imagename+"_segm.tif"); 
	run("Enhance Contrast", "saturated=0.35");	
	run("Threshold..."); 
	getThreshold(lower,upper);
	waitForUser("Cell "+i+1+" segmentation", "Adjust the threshold of the cell nr "+i+1+ "\non the "+imagename+" acquisition and click on OK"); 
	run("Convert to Mask", " background=Dark black"); 
	saveAs("Tiff",res_dir+i+1+"_"+imagename+"_binary.tif");

	return newArray(lower,upper);
}

// -------------------------------------------------------------

function DialogboxMask(DefaultBinaryChoice,DefaultIncludeEdges){

	Dialog.create("Choice of the final binary image");                   
	Dialog.addRadioButtonGroup("From which binary image \ndo you want to create the mask?", newArray("Sum", "C1","C2"), 1,3, DefaultBinaryChoice);
	Dialog.addRadioButtonGroup("Do you want to analyse edges-touching particles?", newArray("no", "yes"), 1, 2, DefaultIncludeEdges);
	Dialog.show; 
	MaskChoice=Dialog.getRadioButton();
	Include_edges=Dialog.getRadioButton();
	
	return newArray(MaskChoice,Include_edges);
}

// -------------------------------------------------------------

function MaskGeneration(maskchoice, include_edges,i,MinSize,MaxSize,MinCirc,MaxCirc,IncludeHoles,res_dir){

	// --- Generating the final binary image depending on the user's choice
	
	if(maskchoice=="C1"){
		selectWindow(i+1+"_C1_binary.tif") ;
		run ("Duplicate...", "duplicate");
		rename(i+1+"_binary.tif");
	}
	if(maskchoice=="C2"){
		selectWindow(i+1+"_C2_binary.tif") ;
		run ("Duplicate...", "duplicate");
		rename(i+1+"_binary.tif");
	}
	if(maskchoice=="Sum"){
		imageCalculator("Add create stack", toString(i+1)+"_C1_binary.tif",toString(i+1)+"_C2_binary.tif");
		selectWindow("Result of "+i+1+"_C1_binary.tif"); 
		rename(i+1+"_binary.tif");
	}
		
	selectWindow (i+1+"_C1_binary.tif"); 
	close(); 
	selectWindow (i+1+"_C2_binary.tif"); 
	close();

	// --- Filtering the objects > 100 µm² touching the edges or not to generate the final mask

	selectWindow(i+1+"_binary.tif");
	if(include_edges=="no"){
		run("Analyze Particles...", "size="+MinSize+"-"+MaxSize+" circularity="+MinCirc+"-"+MaxCirc+" show=Masks exclude clear "+IncludeHoles+"stack");
			} else {
		run("Analyze Particles...", "size="+MinSize+"-"+MaxSize+" circularity="+MinCirc+"-"+MaxCirc+" show=Masks clear "+IncludeHoles+"stack");
	}
	selectWindow("Mask of "+i+1+"_binary.tif");
	saveAs("Tiff",res_dir+i+1+"_mask.tif"); 
	
	// Division of the mask image by 255 so that the pixels are worth 1 if inside the segmented cell and 0 if outside.
	
	selectWindow(i+1+"_mask.tif"); 
	run("Divide...", "value=255 stack");
}

// -------------------------------------------------------------

function IntensityMeasurment(imagename,i,nb_slices){
	
	// --- Division of the rectangle by the mask so that the pixels outside the cell are replaced by "NaN".

	imageCalculator("Divide create 32-bit stack", toString(i+1)+"_ROI_"+imagename+".tif",toString(i+1)+"_mask.tif");

	// --- Measuring the mean intensity on each slice

	run("Set Measurements...", "mean redirect=None decimal=3");
	selectWindow("Result of "+i+1+"_ROI_"+imagename+".tif");
	for(j=0;j<nb_slices;j++){
		setSlice(j+1);
		run("Measure");
	}

	// --- Storing the results in an arra and returning it
	
	vect_intensities=newArray(nResults);
	for(j=0;j<nResults;j++){
		vect_intensities[j]=getResult("Mean", j);
	}
	Table.deleteRows(0, nResults-1);
	
	return vect_intensities;	
}

// -------------------------------------------------------------
function CloseImagesAndClean(i){
	selectWindow (i+1+"_ROI_C1.tif"); close(); 
	selectWindow ("Result of "+i+1+"_ROI_C1.tif"); close(); 
	selectWindow ("Result of "+i+1+"_ROI_C2.tif"); close() ;
	selectWindow ("Result of "+i+1+"_ROI_ratio.tif"); close();
	selectWindow (i+1+"_binary.tif"); close(); 
	selectWindow (i+1+"_mask.tif"); close(); 
	Table.deleteRows(0, nResults-1);
	roiManager("reset");
}

// -------------------------------------------------------------
function ResultsSaving(nb_ROI,times,i_C1,i_C2,ratio){
	
	// --- Array of the stacked images numbers

	ArrayCellNumber=newArray(0);
	for(i=0;i<nb_ROI;i++){
		CellNumber=newArray(lengthOf(times));
		for(j=0;j<lengthOf(times);j++){
			CellNumber[j]=i+1;
		}
		ArrayCellNumber=Array.concat(ArrayCellNumber,CellNumber);
	}
	
	// --- Array of the stacked times
	
	ArrayTime=newArray(0);
	for(i=0;i<nb_ROI;i++){
		ArrayTime=Array.concat(ArrayTime,times);
	}

	// --- Table filling with every measured value
	
	Table.create("Summary");
	Table.setColumn("cell_number", ArrayCellNumber);
	Table.setColumn("time",ArrayTime);
	Table.setColumn("i_C1",i_C1);
	Table.setColumn("i_C2",i_C2);
	Table.setColumn("ratio",ratio);
	waitForUser("Choose the file where you want to save the summary matrix"); 
	saveAs("Results");
}


// -------------------------------------------------------------
function UsersChoiceSaving(imagename,ratio,BG,binning,med,gauss,radii,nbcells,mincellsize,maxcellsize,mincirc,maxcirc,includeholes,dir,filename){
	Report=File.open(dir+filename+".txt");
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);	
	print(Report,"Ratio analysis of the "+imagename+" image");
	print(Report,"Performed on the "+year+"."+month+"."+dayOfMonth+" at "+hour+"h "+minute+"min");
	print(Report,"Calculated ratio = "+ratio);
	print(Report,"");
	print(Report,"Image preprocessing:");
	if(BG=="yes") { print(Report,"	Background subtraction performed");
	} else { print(Report,"	No background subtraction performed");
	}
	if(binning=="yes") { print(Report,"	2*2 additive binning performed");
	} else { print(Report,"	No additive binning performed");	
	}
	if(med=="yes"){ print(Report,"	"+radii[0]+"-pixel-radius median filtering");
	} else { print(Report,"	No median filtering");
	}
	if(gauss=="yes"){ print(Report,"	"+radii[1]+"-pixel-radius gaussian filtering");
	} else { print(Report,"	No gaussian filtering");	
	}
	print(Report, nbcells+" Analyzed cells");
	print(Report,"Analyzed particles:");
	print(Report,"	Minimum cells size = "+mincellsize+" µm²");
	print(Report,"	Maximum cell size = "+maxcellsize+" µm²");
	print(Report,"	Minimal circularity = "+mincirc);
	print(Report,"	Maximal circularity = "+maxcirc);
	if(includeholes==""){ print(Report,"	No holes inclusion");
	} else { print(Report,"	Holes inclusion");
	}
	print(Report,"Parameters for individual cell segmentation");
	File.close(Report);
}

// -------------------------------------------------------------

function UsersChoiceCompleting(cellnum,minthresholdC1,maxthresholdC1,minthresholdC2,maxthresholdC2,maskchoice,IncludeEdges,dir,filename){
	File.append("	Cell "+cellnum+1+" segmentation", dir+filename+".txt");
	File.append("		Chosen binary image: "+maskchoice, dir+filename+".txt");
	if(IncludeEdges=="yes"){
		File.append("		Edge-touching particles included", dir+filename+".txt"); 
	}else{
		File.append("		Edge-touching particles not included", dir+filename+".txt");
	}
	File.append("		Threshold values:", dir+filename+".txt");
	if(maskchoice=="C1"){
		File.append("			C1: "+minthresholdC1+" - "+maxthresholdC1, dir+filename+".txt");
	}
	if(maskchoice=="C2"){
		File.append("			C2: "+minthresholdC2+" - "+maxthresholdC2, dir+filename+".txt");
	}
	if(maskchoice=="Sum"){
		File.append("			C1: "+minthresholdC1+" - "+maxthresholdC1, dir+filename+".txt");
		File.append("			C2: "+minthresholdC2+" - "+maxthresholdC2,dir+filename+".txt");
	}
	
}