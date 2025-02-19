// This macro runs Bio-formats to convert ND2 files to tif and perform operations e.g: Gaussian blur
// This macro prepares data for SuperSegger and Morphometrics analysis (MorphoSegger)
// Author: Andres Florez     April 27, 2020, Harvard University


//Select input, and create output directory
inDir = getDirectory("Choose input Directory");
outDir = inDir+File.separator+"Analysis"+File.separator;
if (File.exists(outDir))
      exit("Directory already exists");
else
File.makeDirectory(outDir);
filenames = getFileList(inDir); 

//Batch Mode
setBatchMode(true);

//Set parameters for image files
timePrefix="_t";
extension=".nd2";

//Select the frames to analyze (ajdust this for each experiment)
timeStart=1;
timeEnd=89;

pattern = ".*"; // for selecting all the files in the folder
// different examples of regex-based selection 
//pattern = "01-03.*Pos [3-9].*";
//pattern = ".*Pos [7-9].*";
//pattern = "01-02.*";

count = 0;
for (i = 0; i < filenames.length; i++) {
	currFile = inDir+filenames[i];
	if(endsWith(currFile, extension) && matches(filenames[i], pattern)) { // process matching regex, 
		count++;
        run("Bio-Formats Importer", "open=currFile autoscale color_mode=Default rois_import=[ROI manager] specify_range split_channels view=Hyperstack stack_order=XYCZT t_begin=&timeStart t_end=&timeEnd t_step=1");
        fluor = getImageID(); //these steps are necessary to select the image window properly
        selectImage(fluor);
        fluorName= getTitle();
        titleFluor = replace(fluorName, extension, "xy"); // Adjust the renaming according to the initial ND2 file name
        rename(titleFluor);
        run("Gaussian Blur...", "sigma=1 stack");
        run("Image Sequence... ", "format=TIFF name=["+titleFluor+timePrefix+"] save=&outDir");
        //saveAs("Tiff", outDir+titleFluor);
        close();
        phase = getImageID();
        selectImage(phase);
        phaseName= getTitle();
        titlePhase = replace(phaseName, extension, "xy");
        rename(titlePhase);
        run("Image Sequence... ", "format=TIFF name=["+titlePhase+timePrefix+"] save=&outDir");
        close();
        
        
	}
}
print("Number of files processed: "+count);

