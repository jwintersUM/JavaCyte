/*
JavaCyte: AUTOMATIC ANALYSIS OF STRUCTURAL REMODELING PROPERTIES IN CARDIAC TISSUE SAMPLES
version: 1.0
created by J.Winters, S.Verheule, M. Von Braunmuhl
Maastricht University, department of physiology

Correspondence to:
s.verheule@maastrichtuniversity.nl
-------------------------------------------
*/

setBatchMode(true); 

// ====== Generate a dialog box and retrieve users choices ==============================================================================================================================
// Lets the user specify which analysis needs to be performed, in which magnification the images are made etc.

// Image properties
Dialog.create("Settings");
Dialog.addMessage('Settings:');
Dialog.addChoice("Input",																									newArray("Single file", "Folder"),"Folder" );
Dialog.addString("File suffix: ", ".bmp", 5);
Dialog.addNumber("Scale in px/um", 5.3364); //This allows the user to determine the scaling. It is important that this is tested for each microscope separately, once.

//Analysis properties
Dialog.addMessage('Analysis:');
//Dialog.addCheckbox("cdidiameter", false);
Dialog.addCheckbox("Total ECM content", false); 
Dialog.addCheckbox("Cell to cell distance", false);
Dialog.addCheckbox("Cardiomyocyte dissociation index", false);
Dialog.addCheckbox("Cell diameter", false);
Dialog.addCheckbox("Capillary count and size", false);
Dialog.addCheckbox("Fibroblast density", false);
Dialog.addCheckbox("Myocyte count", false);

Dialog.addHelp('url');
Dialog.show();

// retrieve values 
// The preferences entered by the user in the dialoge box 'Settings' are retrieved and stored in variables.
map = Dialog.getChoice();
suffix = Dialog.getString();
scale = Dialog.getNumber();
scaleratio = scale/5.3364;

//cdidiam = Dialog.getCheckbox();
fibr = Dialog.getCheckbox();
dis = Dialog.getCheckbox();
myocyteIsolationIndex = Dialog.getCheckbox();
diam = Dialog.getCheckbox();
capil = Dialog.getCheckbox();
fbCount = Dialog.getCheckbox();
myocount = Dialog.getCheckbox();

//=====================================================================================================================
// If the user wants to run a batch of images and there are still images open, a message is displayed to close all images.
// if the user only wants to analyse 1 image, but has multiple (or no) images open, a message is displayed stating that the user should open only 1 image.
if (map == 'Folder' && nImages != 0){
	exit("STOP! \nclose all images!");
}
if (map != 'Folder' && nImages != 1){
	exit("STOP! \nopen 1 image!");
}

// You have to select at least one test to run.
if (dis == false && diam == false && myocyteIsolationIndex == false && capil == false && fibr == false  && fbCount == false && myocount == false){
	exit("STOP! \nPlease choose at least 1 option");
}

output = getDirectory("Output directory for result table and control images"); // the user selects the output directory for exported images

// Creating list of images in directory
if (map == 'Folder'){
	if (fibr == true || dis == true || myocyteIsolationIndex ==true || diam == true || capil ==true || myocount ==true){
		input = getDirectory("input directory"); 
	}
	
	if (fbCount == true){
		viminput = getDirectory("input directory for vimentin images");
		vimgsiinput = getDirectory("input directory for gsi images");
	}
}

//============================= Part 2: initiate processing of input files============================================================================================================================================================================
if (map=='Folder'){
	if (fibr == true || dis == true || myocyteIsolationIndex ==true || diam == true || capil == true || myocount==true){processFolder(input);}
	if (fbCount == true){processFolderVim(viminput);}
}else{file = getTitle();}

// a forloop is created. This forloop is excecuted on the files in the input folder.
//The loop opens up a file with correct suffix, calls the appropriate analysis function to analyse the file, than creates a result file and control picture with the same name. 

function processFolder(input){
	list = getFileList(input);
	total = list.length;
    	for (l = 0; l < list.length; l++) {
        	if(endsWith(list[l], suffix) == true){
        		open(input+list[l]);
        		file = getTitle();
				print("Start processing of image " + l+1 + "out of " + total);
				
				createWI(output,file);
				if (dis == true || myocyteIsolationIndex == true || diam == true || myocount == true) {
					noise = noiseSetting(output,file);
					myocyteRecognition(output,file);
					if (dis == true || myocyteIsolationIndex == true) {celltocelldistance(output, file);}
					if (diam == true ){diameter(output, file);}	
					if (myocount == true) {myocytecount(output,file);}
				}

				if (fibr==true){fibrosis(output,file);}
			
			close("recimage");
			close("WI");
			close(file);		
		}
	}
	print("all files analyzed, results are saved");
}
		

// fibroblasts
function processFolderVim(viminput) {
    vimlist = getFileList(viminput);
	vimgsilist = getFileList(vimgsiinput);
    for (l = 0; l < vimlist.length; l++) {
        if(endsWith(vimlist[l], suffix3) == true){
        	open(viminput+vimlist[l]);
        	vimfile = getTitle();
			open(vimgsiinput + vimgsilist[l]);
			vimgsifile = getTitle();
            fibroblast(output, vimfile);
			close(vimfile);
			close(vimgsifile);
           }
    }
}


// =======================================Noise setting=======================================================================

function noiseSetting(output,file){
//===setting ideal prominence is required for optimal recognition of myocytes. The program determines the ideal setting based on picture specifications

selectWindow(file);
getStatistics(area, mean, min, max, std);
mean = mean;
stdev = std;

gb = 7;

selectWindow(file);
run("Duplicate...", "title=origins");
selectWindow("WI");
run("Duplicate...", "title=Black");
selectWindow("Black");

run("Create Selection");
selectWindow("origins");
run("Restore Selection");

run("Set Measurements...", "mean standard min median redirect=None decimal=3");
run("Measure");

meanECM = getResult("Mean", 0);
SDECM = getResult("StdDev", 0);
MinECM = getResult("Min", 0);
MaxECM = getResult("Max",0);
MedianECM = getResult("Median", 0);

run("Clear Results");
selectWindow("origins");
run("Make Inverse");
run("Measure");
meanCell = getResult("Mean",0);
SDCell = getResult("StdDev", 0);
MinCell = getResult("Min", 0);
MaxCell = getResult("Max",0);
MedianCell = getResult("Median", 0);
run("Clear Results");

close("Black");
close("origins");
noise = 2.46291+(0.17304*stdev)+(0.39708*mean)-(0.10713*meanECM)-(0.07971*SDECM)-(0.33095*meanCell)+(0.45177*SDCell);
print(noise);

return noise;
}




// ====================================== part 3A: create a working image===========================================
function createWI(output, file) {

print("creating a working image for " + file);	
//run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");

//parameters	:
getStatistics(area, mean, min, max, std);
mean1 = mean;
stdev1 = std;
phansalkarradius  = 0.8* mean1 + stdev1; 
minsize = 1000*scaleratio;
maxsize = 27000*scaleratio;


if (dis == true || myocyteIsolationIndex ==true || diam == true || myocount==true){
		
		run("Duplicate...", "title=WI");
		run("Clear Results");
		run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=1.50 mask=*None* fast_(less_accurate)");
		run("Subtract Background...", "rolling=750");
		run("Gaussian Blur...", "sigma=3");
		run("8-bit");
		run("Auto Local Threshold", "method=Phansalkar radius=phansalkarradius parameter_1=0 parameter_2=0 white");
		run("Median...", "radius=2");
		run("Invert");
		run("Fill Holes");
		run("Options...", "iterations=5 count=1 black do=Open");
		run("Clear Results");
		
		run("Set Measurements...", "area centroid redirect=None decimal=3");
		run("Analyze Particles...", "size=0-minsize circularity=0.00-1.00 exclude");

		Spx = newArray(nResults);
    		Spy = newArray(nResults);
		for (i = 0 ;i <nResults; i++){
			Spx[i] = getResult("X", i);
   			Spy[i] = getResult("Y", i);
   		}	
	
		run("Select None");

		for (i=0; i<nResults; i++){
			doWand(Spx[i],Spy[i]);
			run("Enlarge...", "enlarge=2");
			run("Colors...", "foreground=black background=white selection=yellow");
			run("Fill", "slice");
		}
		run("Select None");
		run("Clear Results");
		
		selectWindow("WI");
		
		print("Creating WI complete");
		run("Input/Output...", "jpeg=85");
	 	saveAs("Jpeg", output + file +"wi");
	 	rename("WI");
		
		run("Select None");
		selectWindow(file);
}

if (fibr==true){
	run("Duplicate...", "title=WIfibrosis");
	run("8-bit");
	run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=2 mask=*None* fast_(less_accurate)");
	run("Gaussian Blur...", "sigma=2");
	run("Auto Local Threshold", "method=Phansalkar radius=phansalkarradius parameter_1=0 parameter_2=0 white");
	run("Median...", "radius=2");
	run("Options...", "iterations=2 count=1 black do=Close");
	saveAs("jpeg", output + file + "overallfibrosis");
	selectWindow("WIfibrosis");
	WIfibrosis = getTitle();
	
}
		
if (capil == true){
		
	}

if (fbCount == true){
		
	}


}


// ====================================== part 3B: create myocyte recognition image===========================================
function myocyteRecognition(output,file) {
print("creating image for recognition of myocytes of " + file);	
//run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");

selectWindow(file);
run("Duplicate...", "title=originsri");
selectWindow("WI");
run("Duplicate...", "title=Blackri");
selectWindow("Blackri");

run("Create Selection");
selectWindow("originsri");
run("Restore Selection");

run("Set Measurements...", "mean standard min median redirect=None decimal=3");
run("Measure");

meanECM = getResult("Mean", 0);
SDECM = getResult("StdDev", 0);
threshold = meanECM - SDECM;

gb = 7;
close("Blackri");
close("originsri");

selectWindow(file);

//Step 2: image preparation
	run("Duplicate...", "title=recimage");
	run("Duplicate...", "title=recimage2");
	run("Enhance Local Contrast (CLAHE)","blocksize=127 histogram=256 maximum=1.5 mask=*None* fast_(less_accurate)");
	run("Subtract Background...","rolling=500");
	
//step 4: find local maxima on recimage2
	run("Clear Results");
	selectWindow("recimage");
	run("Duplicate...", "title=recognize");	
	run("Gaussian Blur...", "sigma=gb");
	run("Find Maxima...", "prominence=noise exclude light output=[Point Selection]");
	getSelectionCoordinates(x, y);
	

// Step 5: check if the local maxima are actually in a cell, and not just a dark spot within the matrix
// this can be done by checking the size of the particle and by checking the mean pixel value of a circle around the maximum and its stdev
// change the pixel value of local maxima that are not cells, so they are not recognized as a separate segment anymore

selectWindow("recimage2");
run("Duplicate...", "title=measure");
run("Gaussian Blur...", "sigma=gb");

for (i=0; i<x.length;i++){
	selectWindow("measure");
	run("Clear Results");
	setTool("oval");
	b = x[i]-5;
	c = y[i]-5;
	makeOval(b,c,10, 10);
	getStatistics(area, mean, min, max, std);
	meangray = mean;
	if (meangray >threshold){
		selectWindow("recimage");
		setTool("oval");
		d = x[i]-30;
		e = y[i]-30;
		makeOval(d,e,60, 60);
		setForegroundColor(93,0,14);   
		run("Fill", "slice");
	}
	selectWindow("WI");
	val = getPixel(x[i],y[i]);
	if ( val <25 ){
		selectWindow("recimage");
		setTool("oval");
		d = x[i]-30;
		e = y[i]-30;
		makeOval(d,e,60, 60);
		setForegroundColor(93,0,14);   
		run("Fill", "slice");
		
	}
}

print("Creating rI complete");
		selectWindow("recimage");
		run("Select None");
		run("Input/Output...", "jpeg=85");
	 	saveAs("Jpeg", output + file +"ri");
	 	rename("recimage");
	 	
		
selectWindow("recimage2");
close();
		
selectWindow("measure");
close();

selectWindow("recognize");
close();
selectWindow(file);



}



//=======overal fibrosis================================================================================================================================================================
function fibrosis(output, file) {

	print("Processing: overall fibrosis of " + file);	
	run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");

//Calculate percentage of overall fibrosis

	selectWindow("WIfibrosis");
	
	run("Clear Results");
	nBins = 256;
	row = 0;
	getHistogram(values, counts, nBins);
	for (i=0; i<nBins; i++) {
	    setResult("Value", row, values[i]);
	    setResult("Count", row, counts[i]);
		row++;
	}
	updateResults();

	blackpx = getResult("Count",0);
	whitepx = getResult("Count",255);
	totalpx = blackpx + whitepx;
	percentage = (whitepx/totalpx)*100;
	
	run("Clear Results");
	setResult("% Overall fibrosis", 0, percentage);
	setResult("Name", 0, file);;
	
updateResults;	

saveAs("Results", output + file +"overallfibrosis.csv");
run("Clear Results");
close("WIfibrosis");


}



// ========= Intermyocyte distance and cardiomyocyte dissociation index  ========================================================================================================================
function celltocelldistance(output, file) {

print("Processing: cell to cell distance of " + file);	
//run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");

extremectc = 14;

//  initialize variables
var CxR = newArray();
var CyR = newArray();
var NxR = newArray();
var NyR = newArray();
var WidthR = newArray();
var TotalwidthR = newArray();
var PeaksR = newArray();
var cmisindex = newArray();
var cdi50xR = newArray();
var cdi50yR = newArray();
var cdi50R = newArray();
gb=7;
// getting recimage and WI ready
		selectWindow("recimage");
		run("Duplicate...", "title=recimageBlurred");
		run("Gaussian Blur...", "sigma=gb");
		//run("Gaussian Blur...", "sigma=2");

		selectWindow("WI");
		run("Invert");
		run("Duplicate...", "title=CTCWI");


// Step 6: segment particles to get all cells
	
	selectWindow("recimage");
	run("Duplicate...", "title=CTC_ctr");
	run("Duplicate...", "title=CTC");
	selectWindow("CTC");
	run("Gaussian Blur...", "sigma=gb");
	//run("Gaussian Blur...", "sigma=2");
	run("Find Maxima...", "prominence=noise exclude light output=[Segmented Particles]");
	rename("Basis");
	selectWindow("CTC");
	close();	
	selectWindow("Basis");


//Phase A: Appointing coordinates to myocytes	
	run("Set Measurements...", "area centroid redirect=None decimal=3");
	run("Analyze Particles...", "size=minsize-maxsize circularity=0.00-0.95 exclude display");
	Cx = newArray(nResults);
    Cy = newArray(nResults);
	for (i = 0 ;i <nResults; i++){
		Cx[i] = getResult("X", i);
   		Cy[i] = getResult("Y", i);
   	}	
  
	
//Phase B: find neighbors for every cell (center cell); C= center cell; N= Neighbor cell
	//create a duplication (="wand") of the skeleton image for the purpose to identify all neighbor cells of the centrum cell in the loop
	// identify neighbors by using the 'wand' function of imageJ to select the entire centrum cell. Then the selection is slightly enlarged and filled with white color. 
	// Therefore the boundaries of the centrum cell diminish. When using the 'wand' function again, the area of the centrum cell together with the neighborcells is selected. 
	// By restoring this selection on the "basis" imagej we can identify the centroid coordinates of the cells within this selection. 

	run("Clear Results");
   	
	for (i=0; i<Cx.length; i++){
		var iso = 0;
		var noniso = 0;
		var remain = 0;
		var totalsepta = 0;
		var isoindex = 0;
		var NumberOfPeaks = 50;
		//var ncount = 0;
		
		//Select one cell
		selectWindow("Basis");
		doWand(Cx[i],Cy[i]);
		selectWindow("recimageBlurred");
		run("Restore Selection");
		run("Find Maxima...", "prominence=noise exclude light output=[Point Selection]");
		getSelectionCoordinates(xc,yc);	
		xc = Array.trim(xc, 1);
		yc = Array.trim(yc, 1);
		

		run("Clear Results");
		selectWindow("recimageBlurred");
		run("Select None");


		// select neighbours of the cell
		selectWindow("Basis");
		run("Select None");
		run("Duplicate...", "title=wand");
		selectWindow("wand");
		doWand(Cx[i],Cy[i]);
		run("Enlarge...", "enlarge=2");
		run("Colors...", "foreground=white background=black selection=yellow");
		run("Fill", "slice");
		doWand(Cx[i],Cy[i]);
		selectWindow("Basis");
		run("Restore Selection");
		run("Set Measurements...", "centroid redirect=None decimal=3");
		run("Analyze Particles...", "size=minsize-maxsize circularity=0.00-0.95 display");
		close("wand");
	//=============================================	
		

		//Resulttable -> Array (Nx and Ny)
		//Centroid of all cells in the restored selection		
		Nx = newArray(nResults);
    	Ny = newArray(nResults);
    	
		for (j = 0 ;j <nResults(); j++){
			Nx[j] = getResult("X", j);
    		Ny[j] = getResult("Y", j);
		}

		//print("Number of neighbors:    " + nResults);
		
		run("Clear Results");
		run("Select None");

		
// Phase C: Loop through all the neighbor cells to make a line selection between the centroid coordinate of the neighbor cell of the loop and the centrum cell of the underlying loop.
				
		for (j=0; j<Nx.length;j++){

			if (Cx[i] != Nx[j]){
				
			selectWindow("Basis");
			doWand(Nx[j],Ny[j]);
			selectWindow("recimageBlurred");
			run("Restore Selection");
			run("Find Maxima...", "prominence=noise exclude light output=[Point Selection]");
			getSelectionCoordinates(xn,yn);
			xn = Array.trim(xn, 1);
			yn = Array.trim(yn, 1);
			
			run("Clear Results");
										
			run("Select None");

			selectWindow("WI");
				
				for(t=0; t<xc.length; t++){
				makeLine(xc[t], yc[t], xn[t], yn[t]);}
				original_plot = getProfile(); //saving plot into array
				tolerance = 25;
				original_peaks= Array.findMaxima(original_plot, tolerance);
				if (original_peaks.length ==0) {NumberOfPeaks =0;}
				for(m=1; m <=original_peaks.length ; m++) {NumberOfPeaks = m;}
			
			//measuring distances
			// peaks that end with y=255 will disturb the loop to measure width. In addition, this measurement is likely not accurate. It should therefore be ignored.
			
			if (NumberOfPeaks != 0){
			
				width = 0;
				begin = 0;
				end = original_plot.length -1;
				if (original_plot[begin]<250 && original_plot[end] <250){
					k = 0;
					beginpoint = 0;
					while(k < original_plot.length){
						if (original_plot[k] > 250) {
							peak_begin = k;
							if (beginpoint == 0){beginpoint=k;}	
							while (original_plot[k] > 250 && k < original_plot.length) {k=k+1;}
							peak_end = k-1;
							distance = peak_end-peak_begin;
							distance = distance/scale;
							width = width + distance;
							totalwidth = peak_end-beginpoint;
							totalwidth = totalwidth/scale;
						}else{k=k+1;}
					}
				} else{break;}
			} else {width = 0;}

			// Adding results to result arrays and draw the measurement onto the CTC and WI image
			if (width >= 1 && width < 50) {
				setForegroundColor(0, 255, 0);	
				setLineWidth(8);
				// edge ctr picture
				selectWindow("CTC_ctr");
				for(l=0; l<xc.length; l++){
				drawLine(xc[l], yc[l], xn[l], yn[l]);}
				
								
				// wi ctr
				selectWindow("CTCWI");
				setForegroundColor(0, 255, 0);	
				setLineWidth(8);
				for(l=0; l<xc.length; l++){
				drawLine(xc[l], yc[l], xn[l], yn[l]);}
				
				CxR = Array.concat(CxR,Cx[i]);
				CyR = Array.concat(CyR,Cy[i]);
				NxR = Array.concat(NxR,Nx[j]);
				NyR = Array.concat(NyR,Ny[j]);
				WidthR = Array.concat(WidthR,width);
				TotalwidthR = Array.concat(TotalwidthR,totalwidth);
				PeaksR = Array.concat(PeaksR, NumberOfPeaks);
				
			}

			if (myocyteIsolationIndex == true){
				if (width > extremectc){
					if(width < 50){iso = iso +1;} else {remain=remain+1;}}
				if (width <= extremectc){
					if (width >= 1){noniso = noniso + 1;}}
				if (width < 1){remain=remain+1;}
			}
				
			run("Clear Results");
		
			}
	} // end for every neighbor
	
	if(myocyteIsolationIndex == true) {
			totalsepta = iso + noniso;
			
			if (totalsepta != 0){
				isoindex = iso/totalsepta;
				cmisindex = Array.concat(cmisindex,isoindex);
			} 
			else {
				isoindex=999;
				cmisindex = Array.concat(cmisindex,isoindex);
			}
			ncount = Array.concat(ncount,totalsepta);
			if (isoindex>0.24){
				if (isoindex<5){
					cdi50xR = Array.concat(cdi50xR,Cx[i]);
					cdi50yR = Array.concat(cdi50yR,Cy[i]);
					cdi50R = Array.concat(cdi50R, isoindex);
				}				
			} 	
	}

}// end for each cell

// Phase D:  writing coordinates and distances in result table and save results as csv
//only if not same combination of coordinates

for (l = 0; l<WidthR.length; l++){
		setResult("Cx", l, CxR[l]);
		setResult("Cy", l, CyR[l]);
		setResult("Nx", l, NxR[l]);
		setResult("Ny", l, NyR[l]);
		setResult("Width", l, WidthR[l]);
		setResult("TotalWidth", l, TotalwidthR[l]);
		setResult("Peaks", l, PeaksR[l]);
} 

updateResults;
run("Input/Output...", "jpeg=85 save_column save_row");
saveAs("Results", output + file +"distance.csv");
run("Clear Results");

if (myocyteIsolationIndex == true){
for (i=0; i<Cx.length; i++){
	setResult("Cx", i, Cx[i]);
	setResult("Cy", i, Cy[i]);				
	setResult("myocyteIsolationIndex", i, cmisindex[i]);
	setResult("neighborcount", i, ncount[i]);
	}

	updateResults;
	run("Input/Output...", "jpeg=85 save_column save_row");
	saveAs("Results", output + file +"isoIndex.csv");
	run("Clear Results");
	
}


//Phase E: creating the controlpicture
	selectWindow("Basis");
	run("Create Selection");
	run("Colors...", "foreground=black background=white selection=gray");
	selectWindow("CTC_ctr");
	run("Restore Selection");
	run("Colors...", "foreground=gray background=white selection=gray");
	run("Draw", "slice");
	run("Input/Output...", "jpeg=85");
	saveAs("Jpeg", output + file +"edgectr");
	run("Close");

	selectWindow("CTCWI");
	run("Input/Output...", "jpeg=85");
	saveAs("Jpeg", output + file +"wictr");
	run("Close");
	close("Basis");

if (myocyteIsolationIndex == true){
	selectWindow("recimage");
	run("Duplicate...", "title=CDI_ctr");
	selectWindow("CDI_ctr");
	run("Gaussian Blur...", "sigma=gb");
	run("Find Maxima...", "prominence=noise exclude light output=[Segmented Particles]");
	rename("cdi_ct");
	close("CDI_ctr");
	selectWindow("cdi_ct");

	for (i=0; i<cdi50xR.length; i++){
		run("Select None");
		selectWindow("cdi_ct");
		run("RGB Color");
		doWand(cdi50xR[i],cdi50yR[i]);
				
		if (cdi50R[i]<0.40){setForegroundColor(255,255,0);}
		if (cdi50R[i]<0.55 && cdi50R[i]>0.39){setForegroundColor(255,201,0);} 
		if (cdi50R[i]<0.70 && cdi50R[i]>0.54){setForegroundColor(255,153,0);}
		if (cdi50R[i]<0.85 && cdi50R[i]>0.65){setForegroundColor(255,105,0);}
		if (cdi50R[i]>0.85){setForegroundColor(255,0,0);}
		
		run("Fill","slice");
		
	}

	run("Input/Output...", "jpeg=85");
	 	saveAs("Jpeg", output + file +"ctr");
	 	run("Close");
}

if (isOpen("Results")) {
   selectWindow("Results");
   run("Close");
} 

run("Select None");

if(isOpen("CTC")) {
	selectWindow("CTC");
	run("Close");
}

close("recimageBlurred");

} 	


//==============================================Cell Diameter=========================================
function diameter(output, file){
print("Processing: Cell diameter of " + file);	
run("Set Scale...", "distance=scale known=1 pixel=1 unit=um");

selectWindow("recimage");
//selectWindow(file);
//run("Duplicate
gb = 7;
run("Duplicate...", "title=recdiam");
run("Gaussian Blur...", "sigma=gb");
run("Find Maxima...", "prominence=noise exclude light output=[Segmented Particles]");
rename("Basisdiam");
run("Create Selection");
//run("Enlarge...", "enlarge=4");
selectWindow("WI");
run("Duplicate...", "title=diam");
run("Restore Selection");
setForegroundColor(0, 0, 0);
run("Fill", "slice");
run("Select None");
setOption("BlackBackground", true);
run("Make Binary");
//run("Invert");
run("Set Measurements...", "area centroid mean shape feret's redirect=None decimal=3");
run("Set Scale...", "distance=scale known=1 pixel=1 unit=um");
run("Analyze Particles...", "size=minsize-maxsize  circularity=0.00-0.95 show=[Overlay Masks] display exclude");
selectWindow("diam");
run("Restore Selection");
setForegroundColor(0, 0, 0);
run("Fill", "slice");
run("Select None");
saveAs("jpeg", output + file + "overlaymask");
close();
selectWindow("recdiam");
close();	
selectWindow("Basisdiam");
close();

run("Select None");
	
DCx = newArray(nResults);
DCy = newArray(nResults);
area = newArray(nResults);
circ = newArray(nResults);
mf = newArray(nResults);
for (i = 0 ;i <nResults(); i++){
	DCx[i] = getResult("X", i);
	DCy[i] = getResult("Y", i);
   	area[i] = getResult("Area", i);
	circ[i] = getResult("Circ.",i);
  	mf[i] = getResult("MinFeret", i);
  		
	}
	run("Clear Results");
		
	for (g = 0; g<DCx.length; g++){
		setResult("Cx", g, DCx[g]);
		setResult("Cy", g, DCy[g]);
		setResult("Area", g, area[g]);
		setResult("Circ", g, circ[g]);
		setResult("MinFeret", g, mf[g]);
		
	}
	saveAs("Results", output + file +"diameter.csv");
	run("Clear Results");


if (isOpen("Results")){
  selectWindow("Results");
  close();
}
  
}


//=========process Myocyte count===============================================================================================
function myocytecount(output, file) {

print("Processing: Counting myocytes of " + file);	
//run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");

gb = 7;
selectWindow("recimage");
run("Duplicate...", "title=count");
run("Gaussian Blur...", "sigma=gb");
run("Find Maxima...", "prominence=noise exclude light output=[Segmented Particles]");
rename("Basis");
selectWindow("count");
close();	
run("Clear Results");
selectWindow("Basis");


//Phase A: Appointing coordinates to myocytes	
	run("Set Measurements...", "area centroid redirect=None decimal=3");
	run("Analyze Particles...", "size=minsize-maxsize circularity=0.00-0.95 display exclude");
	count = nResults;
	CCx = newArray(nResults);
    CCy = newArray(nResults);
	mf = newArray(nResults);
	for (i = 0 ;i <nResults(); i++){
		CCx[i] = getResult("X", i);
   		CCy[i] = getResult("Y", i);
   		mf[i] = getResult("MinFeret", i);
  	}
	run("Clear Results");
	setResult("Count", 0, count);
	saveAs("Results", output + file +"myocytecount.csv");

	run("Clear Results");
	
	for (i = 0; i<CCx.length; i++){
		setResult("Cx", i, CCx[i]);
		setResult("Cy", i, CCy[i]);
		setResult("MinFeret", i, mf[i]);
	}

	saveAs("Results", output + file +"myocytecountCentroids.csv");
	
	selectWindow("Basis");
	close();	

}

// ======= process image to calculate amount and size of capillaries ========================================================
function capillaries(output, file){

print("Processing: Capillaries of " + file);	
run("Set Scale...", "distance=scale known=1 pixel=1 unit=um");
sf=scale/5.336;

//======phase A: image preparation===========================================
selectWindow(file);
run("Duplicate...", "title=gsi"); 
run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=2.50 mask=*None* fast_(less_accurate)");
run("Gaussian Blur...", "sigma=5");

min=newArray(3);
max=newArray(3);
filter=newArray(3);
a=getTitle();
run("HSB Stack");
run("Convert Stack to Images");
selectWindow("Hue");
rename("0");
selectWindow("Saturation");
rename("1");
selectWindow("Brightness");
rename("2");
getStatistics(area, mean, min[2], max[2], std);
meanB = mean;
stdevB = std;
B = meanB + (1.25*stdevB);     // The value of B will be used as a brightness value to threshold the capillary image

min[0]=34;
max[0]=129;
filter[0]="pass";
min[1]=0;
max[1]=255;
filter[1]="pass";
min[2]=B;
max[2]=255;
filter[2]="pass";
for (i=0;i<3;i++){
  selectWindow(""+i);
  setThreshold(min[i], max[i]);
  run("Convert to Mask");
  if (filter[i]=="stop")  run("Invert");
}
imageCalculator("AND create", "0","1");
imageCalculator("AND create", "Result of 0","2");
for (i=0;i<3;i++){
  selectWindow(""+i);
  close();
}
selectWindow("Result of 0");
close();
selectWindow("Result of Result of 0");
rename("capillaries");
setOption("BlackBackground", true);
run("Make Binary");
run("Median...", "radius=10");
run("Options...", "iterations=5 count=1 black do=Close");
run("Fill Holes");
if (coimage == true) {
	run("Duplicate...", "title=shapecirc");
}



// =======================Phase B: analysis of capillary count and shape=============================================
selectWindow("capillaries");
run("Set Measurements...", "centroid area shape feret's redirect=None decimal=3");
run("Analyze Particles...", "size=2-Infinity show=Overlay display");
if(nResults > 10){
	run("Input/Output...", "jpeg=85 save_column save_row");
	saveAs("Results", output + file +"capillaries.csv");
	run("Clear Results");
	}

if (coimage == true) {
	selectWindow("capillaries");
	run("Overlay Options...", "stroke=none width=0 fill=green set apply show");
	saveAs("jpeg", output + file + "sizeFilter");
}

run("Close");

// ============Phase C: creating the control image==========================
/*if (coimage == true){
	selectWindow("shapecirc");
	run("Set Measurements...", "centroid area shape feret's redirect=None decimal=3");
	run("Analyze Particles...", "size=4-Infinity circularity=0.01-1.00 show=Overlay display");
	run("Overlay Options...", "stroke=none width=0 fill=magenta set apply show");
	saveAs("jpeg", output + file + "sizecircFilterGB2B55");
	run("Distribution...", "parameter=MinFeret automatic");
	saveAs("jpeg", output + file + "MinFeretDistribution");
	run("Close");
	selectWindow("shapecirc");
	run("Close");
	run("Clear Results");
	}
*/
run("Clear Results");
close(file);
}

//========= Fibroblast analysis ==================================================================================
function fibroblast(output, vimfile) {

print("Processing: Fibroblast count of " + vimfile);	
run("Set Scale...", "distance=scale known=1 pixel=1 unit=um");

C=0;
D=0;
meanD = 0;
stdevD = 0;
mean = 0;
std =0;
meanC = 0;
stdevC = 0;


// ==========================Phase A: color tresholding vimfile=============================================================

selectWindow(vimfile);
//run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=1.50 mask=*None* fast_(less_accurate)");
run("Duplicate...", "title=vim"); 
run("Color Threshold...");
min=newArray(3);
max=newArray(3);
filter=newArray(3);
a=getTitle();
run("HSB Stack");
run("Convert Stack to Images");
selectWindow("Hue");
rename("0");
selectWindow("Saturation");
rename("1");
selectWindow("Brightness");
rename("2");
getStatistics(area, mean, min[2], max[2], std);
meanC = mean;
stdevC = std;
C = meanC + (1.40*stdevC);
print("The blue threshold value for " + vimfile + "is..." + C);
print("The mean value for " + vimfile + "is..." + meanC);
print("The std for " + vimfile + "is..." +stdevC);
min[0]=125;
max[0]=205;
filter[0]="pass";
min[1]=0;
max[1]=255;
filter[1]="pass";
min[2]=C;
max[2]=255;
filter[2]="pass";
for (i=0;i<3;i++){
  selectWindow(""+i);
  setThreshold(min[i], max[i]);
  run("Convert to Mask");
  if (filter[i]=="stop")  run("Invert");
}
imageCalculator("AND create", "0","1");
imageCalculator("AND create", "Result of 0","2");
for (i=0;i<3;i++){
  selectWindow(""+i);
  close();
}
selectWindow("Result of 0");
close();
selectWindow("Result of Result of 0");
rename("blue");
setOption("BlackBackground", true);
run("Make Binary");
run("Options...", "iterations=5 count=1 black do=Close");
run("Fill Holes");


//==================Phase B: tresholding vimgsifile=================================================
selectWindow(vimgsifile);
run("Duplicate...", "title=black");
run("Duplicate...", "title=vimgsi"); 
min=newArray(3);
max=newArray(3);
filter=newArray(3);
a=getTitle();
run("HSB Stack");
run("Convert Stack to Images");
selectWindow("Hue");
rename("0");
selectWindow("Saturation");
rename("1");
selectWindow("Brightness");
rename("2");
getStatistics(area, std, min[2], max[2], std);
meanD = mean;
stdevD = std;
D = meanD + (1.25*stdevD);
print("The green threshold value for " + vimfile + "is..." + D);
rename("2");
min[0]=34;
max[0]=129;
filter[0]="pass";
min[1]=0;
max[1]=255;
filter[1]="pass";
min[2]=D;
max[2]=255;
filter[2]="pass";
for (i=0;i<3;i++){
  selectWindow(""+i);
  setThreshold(min[i], max[i]);
  run("Convert to Mask");
  if (filter[i]=="stop")  run("Invert");
}
imageCalculator("AND create", "0","1");
imageCalculator("AND create", "Result of 0","2");
for (i=0;i<3;i++){
  selectWindow(""+i);
  close();
}
selectWindow("Result of 0");
close();
selectWindow("Result of Result of 0");
rename("green");
setOption("BlackBackground", true);
run("Make Binary");run("Options...", "iterations=5 count=1 black do=Close");
run("Fill Holes");

//=========== phase C: creating black image for stack======================================
selectWindow("black");
run("Color Threshold...");
min=newArray(3);
max=newArray(3);
filter=newArray(3);
a=getTitle();
run("HSB Stack");
run("Convert Stack to Images");
selectWindow("Hue");
rename("0");
selectWindow("Saturation");
rename("1");
selectWindow("Brightness");
rename("2");
min[0]=0;
max[0]=255;
filter[0]="stop";
min[1]=0;
max[1]=255;
filter[1]="pass";
min[2]=86;
max[2]=255;
filter[2]="pass";
for (i=0;i<3;i++){
  selectWindow(""+i);
  setThreshold(min[i], max[i]);
  run("Convert to Mask");
  if (filter[i]=="stop") {run("Invert");}
}
imageCalculator("AND create", "0","1");
imageCalculator("AND create", "Result of 0","2");
for (i=0;i<3;i++){
  selectWindow(""+i);
  close();
}
selectWindow("Result of 0");
close();
selectWindow("Result of Result of 0");
rename(a);
setOption("BlackBackground", true);
run("Make Binary");
run("Invert");

//====================Phase C: putting all images together in a stack=======================================
//selectWindow(vimfile);
//close();
//selectWindow(vimgsifile);
//close();
run("Concatenate...", " title=Stack keep open image1=black image2=green image3=blue");
run("Stack to RGB");
selectWindow("Stack");
close();
selectWindow("black");
close();
selectWindow("Stack (RGB)");
saveAs("jpeg", output + vimfile +"ctrstack");
close();

//===================Phase D: Selecting only vimentin positive items for analysis
fibroblastCountR = newArray(nResults);
imageCalculator("Subtract create", "blue","green");
selectWindow("Result of blue");
run("Dilate");
run("Median...", "radius=5");
run("Options...", "iterations=6 count=1 black do=Close");
run("Fill Holes");

run("Set Measurements...", "area centroid redirect=None decimal=3");
minsizefibroblast = (37.42*scaleratio)*(37.42*scaleratio);
run("Analyze Particles...", "size=minsizefibroblast-Infinity show=Overlay");
//run("Input/Output...", "jpeg=85 save_column save_row");
saveAs("Results", output + vimfile +"fibroblastCount.csv");
run("Overlay Options...", "stroke=none width=0 fill=magenta set apply show");
run("Labels...", "color=white font=10 show draw");
saveAs("jpeg", output + vimfile +"ctrfibroblastCount");
close();



selectWindow("green");
run("Close");
selectWindow("blue");
run("Close");
run("Clear Results");
}

setBatchMode(false);
waitForUser("Results are saved");