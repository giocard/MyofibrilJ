// ImageJ macro to analyse myofibril cross sections
//
// Max Planck Institute of Biochemistry - Giovanni Cardone and Maria Spletter
//
/* Usage
 * - analysis can be run on entire image or on a selected ROI
 * - analysis is run on all the images opened
*/
/* Results
 * Parameter computed: 
 * 	 diameter
 *   density
 *   cluster distribution
 * Diameter:
 *   Initial estimate of diameter is obtained by autocorrelation, and this value is
 *   used to calibrate the filter parameters for the final estimate.
 *   The diameter is eventually estimated by averaging all the corss sections and 
 *   calculating from the radial profile of the average the width at 26% of the maximum height
 * Density:
 *   The density is estimated in two different ways:
 *   - from the nearest neighbor distance d (see below) as 1/d^2 (optimal for regularly distributed points)
 *   - as the ratio between points and area (optimal if all points are detected)
 * Cluster distribution
 *   The cluster distribution is estimated by means of Nearest Neighbor Distance.
 *   The values reported are :
 *   - nearest neighbor distance (average and standard deviation)
 *   - nearest neighbor index (NNI) and p-value (NNI<=0.5 -> clustered; 0.5<NNI<1.5 -> randomly distributed; NNI>=1.5 -> regularly distributed) 
 */
/* Acknowledgements
 *   procedure for calculating autocorrelation is based on script from Michael Schmid (version: 2008-May-14)
 */

requires("1.42p");
setBatchMode(true);

ResultsTable = "Fibril cross sections";
totImages = nImages;

if (totImages==0) {
	exit("No images open!");
}

initializeResultsTable(ResultsTable);

for (in=1; in<=totImages; in++){

    selectImage(in);
    title=getTitle();
    // process only images loaded from file
    // to avoid analysis of results from previous run
    dirname = getInfo("image.directory");
    if (dirname==""){
        continue;
    }

    doSelection = false;
    selType = selectionType();
    if (selType !=-1) doSelection = true;
    
    Stack.getDimensions(width, height, NChannels, NSlices, NFrames);
    Stack.getPosition(currentChannel, currentSlice, currentFrame);
    inputID=getImageID();

    outTitle = title;
    if (NChannels>1)  outTitle = outTitle+"Ch"+currentChannel;
    if (NSlices>1)  outTitle = outTitle+"Sl"+currentSlice;
    if (doSelection == true) outTitle = outTitle+"-crop";

    print("\n"+"--- "+outTitle+" ---");


    // Step 1: initial estimate of diameter by autocorrelation
    // initialize parameters
    if (doSelection==true) {
        getSelectionBounds(x, y, width, height);
    } else {
        width=getWidth;
        height=getHeight;
    }
    minAutoCorrSize = minOf(128,maxOf(width, height));

    getPixelSize(unit, pixelWidth, pixelHeight);
    if(pixelWidth!=pixelHeight) {
        exit("Error: can not handle different pixel size along x and y ("+pixelWidth+","+pixelHeight+")");
    }

    // Extract multiple subregions of decreasing size from the image and calculate the autocorrelation
    // and estimate the diameter as the lowest of the first minima in the radial profile from each autocorrelation
    tmpWidth = width;
    tmpHeight = height;
    initialDiameter = maxOf(width,height);
    for (k=1; k<4 && maxOf(tmpHeight,tmpWidth)>=minAutoCorrSize; k++) {
        for (r=1; r<=k*k; r++) {
            //maxRadius may be modified, should not be larger than 0.3*minOf(width, height);
            maxRadius=0.3*minOf(tmpWidth, tmpHeight);
            kernel = maxOf(maxOf(tmpWidth, tmpHeight)/10,20);
            minFFTsize=1.3*maxOf(tmpWidth, tmpHeight);
            size=4;
            while(size<minFFTsize) size*=2;
            
            hannID = HanningWindow(tmpWidth,tmpHeight);
            selectImage(hannID);
            run("Canvas Size...", "width="+ size+" height="+ size+" position=Center zero");
            makeRectangle(floor((size-tmpWidth)/2), floor((size-tmpHeight)/2), tmpWidth, tmpHeight);
            run("Make Inverse");
            run("Set...", "value=0");
            run("Select None");

            selectImage(inputID);
            Stack.setPosition(currentChannel, currentSlice, currentFrame);
            //make a copy and prefilter image
            tempTitle="temp-"+k+"-"+random();
            run("Duplicate...", "title="+tempTitle);
            tempID=getImageID();
            run("Subtract Background...", "rolling="+kernel);
            // constrast enhancing tends to suppress the minimum in some cases
            //run("Enhance Local Contrast (CLAHE)", "blocksize="+kernel+" histogram=256 maximum=3 mask=*None* fast_(less_accurate)");
            run("Gaussian Blur...", "sigma=1");
            run("32-bit");
            // crop a random subregion in the image
            makeRectangle(floor(random()*(width-tmpWidth)), floor(random()*(height-tmpHeight)), tmpWidth, tmpHeight);
            run("Crop");
            getRawStatistics(nPixels, mean);
            run("Canvas Size...", "width="+ size+" height="+ size+" position=Center zero");
            makeRectangle(floor((size-tmpWidth)/2), floor((size-tmpHeight)/2), tmpWidth, tmpHeight);
            run("Make Inverse");
            run("Set...", "value="+mean);
            run("Select None");
            // apodize using Hann window
            run("Subtract...", "value="+mean);
            imageCalculator("Multiply 32-bit", tempTitle,"Hanning Window");
            run("Add...", "value="+mean);
            getRawStatistics(nPixels, mean);

            //make autocorrelation of particle image
            run("FD Math...", "image1=["+tempTitle+"] operation=Correlate image2=["+tempTitle+"] result=AutoCorrelation do");
            ACimgID=getImageID();
            run("Subtract...", "value="+(nPixels*mean*mean));
            selectImage(tempID);
            close();

            // make autocorrelation reference to correct finite image size effects
            // ** not sure it is needed in this situation **
            newImage("frame", "8-bit White", tmpWidth, tmpHeight, 1);
            run("Set...", "value=255");
            tempID=getImageID();
            tempTitle="tempref-"+k+"-"+random();
            rename(tempTitle);
            run("Canvas Size...", "width="+ size+" height="+ size+" position=Center zero");
            run("FD Math...", "image1=["+tempTitle+"] operation=Correlate image2=["+tempTitle+"] result=AutoCorrReference do");
            refID=getImageID();
            imageCalculator("Divide", ACimgID, refID);
            selectImage(refID);
            close();
            selectImage(tempID);
            close();
            
            selectImage("Hanning Window");
            close();

            //prepare normalized power spectrum for radial averaging
            selectImage(ACimgID);
            circleSize=2*floor(maxRadius)+1;
            norm = getPixel(size/2, size/2);
            run("Divide...", "value="+norm);
            norm = getPixel(size/2, size/2);
            plotRadialProfile(size/2,size/2,floor(maxRadius)-1);
            Plot.getValues(x, y);
            close();
            selectImage(ACimgID);
            close();

            yFinal = newArray(y.length+1);
            yFinal[0] = 1.;
            for (i=0; i<y.length; i++)
                yFinal[i+1]+=y[i];

            xFinal = newArray(y.length+1);
            xFinal[0] = 0;
            for (i=0; i<x.length; i++)
                xFinal[i+1] = x[i];

            //Plot.create("Autocorrelation of "+title, "Distance (pixels)", "Normalized Autocorrelation", xFinal, yFinal);
            //setBatchMode("show");

            // determine size as location of minimum
            lastY = 1;
            minY = 0;
            minX = -1;
            done = false;
            for (i=1; i<yFinal.length && !done; i++) {
                if (yFinal[i] > 1 || yFinal[i] < lastY) {lastY = yFinal[i];}
                else {
                    minY = yFinal[i-1];
                    minX = xFinal[i-1];
                    done=true;
                }
            }

			// keep the minimum value as the initial estimate for the diameter
            if (done==true && minX<initialDiameter){
                initialDiameter = minX;
            }
        }
        tmpHeight/=2;
        tmpWidth/=2;
    }

    // check if pre-estimate of diameter (in pixels) is found
    if (initialDiameter == maxOf(width,height)) {
        if (doSelection) {
            showMessage("ERROR: could not determine an approximate size of the sections! Please select a different subregion.");
        } else {
            showMessage("ERROR: could not determine an approximate size of the sections! Please select a subregion while avoiding empty space.");
        }
        exit();
    }


    // initialize parameters for final estimate of diameter and other parameters
    sigmaGauss = maxOf(1,initialDiameter/8);
    minKernel = 3*initialDiameter;
    boxsize = floor((4*initialDiameter)/2)*2+1;
    area = width * height * pixelHeight * pixelHeight;

    kernelSize = maxOf(minKernel, minOf(width,height)/10);
    contrastSize = minOf(width,height)/2;
    histoSize = minOf(256,width*height/2);
    tempTitle = title+"_processed";
    selectImage(inputID);
    run("Duplicate...", "title="+tempTitle);
    copyImgID = getImageID();
    run("Subtract Background...", "rolling="+kernelSize+" slice");
    run("Gaussian Blur...", "sigma="+sigmaGauss);
    run("Enhance Local Contrast (CLAHE)", "blocksize="+contrastSize+" histogram="+histoSize+" maximum=3 mask=*None* fast_(less_accurate)");
    getRawStatistics(nPixels, mean, min, max, std);
    // find position of all the cross-sections
    threshold = (max-mean)/10;
    run("Find Maxima...", "noise="+threshold+" output=[List]");

    selectImage(copyImgID);
    close();

    // read maxima coordinates and eliminate those touching the edges
    nPoints = 0;
    x = newArray();
    y = newArray();
    for(c=0; c<nResults; c++) {
        xtmp = getResult("X", c);
        ytmp = getResult("Y", c);
        xstart = xtmp-floor(initialDiameter/2);
        ystart = ytmp-floor(initialDiameter/2);
        if (xstart>0 && ystart>0 && (xstart+initialDiameter)<width && (ystart+initialDiameter)<height) {
            x = Array.concat(x,xtmp);
            y = Array.concat(y,ytmp);
            nPoints++;
        }
    }

    // calculate Nearest Neighbor distance for each point/cross-section
    count = 0;
    NN = newArray(nPoints);
    NNdist = newArray((nPoints*nPoints)-nPoints);
    NNsum = 0;

    maxDist = maxOf(width,height);
    for(i1=0; i1<nPoints; i1++) {
       NN[i1] = maxDist;
    }
    for(i1=0; i1<nPoints; i1++) {
        for(i2=0; i2<nPoints; i2++) {
            if(i1!=i2) {
                NNdist[count] = sqrt(pow((x[i2]-x[i1]), 2) + pow((y[i2]-y[i1]), 2));
                if (NNdist[count]<NN[i1]) {
                    NN[i1] = NNdist[count];
                }
                count++;
            }
        }
        NNsum = NNsum + NN[i1];
    }

    NNmean = NNsum / nPoints;

    //determine the sum of the differences
    NNvariance = 0;

    for(i=0; i<nPoints; i++) {
        NNvariance = NNvariance + pow((NN[i]-NNmean), 2);
    }

    //determine variance and standard deviaiton
    NNvariance = NNvariance/(nPoints-1);
    NNstdDev = sqrt(NNvariance);

    NNmean = NNmean * pixelHeight;
    NNstdDev = NNstdDev * pixelHeight;
    densityNN = 1./NNmean/NNmean;
    NNrandom = 0.5 * sqrt(area/nPoints);
    NNI = NNmean/NNrandom;
    NNIse = 0.26136 * sqrt(area/nPoints/nPoints);
    NNIzScore = (NNmean - NNrandom) / NNIse;

    print("\nNumber of spots: "+nPoints);
    print("\n== Nearest neighbor distance analysis ==");
    print("   Average nearest neighbor distance: "+NNmean+" "+unit);
    print("   Standard Deviation: "+NNstdDev+" "+unit);
    print("   Density: "+densityNN+" spots / "+unit+"^2");
    print("   Nearest Neighbor Index: "+NNI);
    if (NNI<=0.5) {
        print("  ->  Spots are presumably clustered");
    } else if (NNI>0.5 && NNI < 1.5) {
        print("  ->  Spots are presumably randomly distributed");
    } else if (NNI >= 1.5) {
        print("  ->  Spots are presumably dispersed, regularly distributed");
    }
    if (NNIzScore<-2.58 || NNIzScore > 2.58) {
        pVal = "< 0.01";
    } else if (NNIzScore<-1.96 || NNIzScore > 1.96) {
        pVal = "< 0.05";
    } else if (NNIzScore<-1.65 || NNIzScore > 1.65) {
        pVal = "< 0.1";
    } else {
        pVal = "> 0.1";
    }
    print("   z-score: "+NNIzScore+" (p-value "+pVal+")");

    nBins = maxOf(15,floor(nPoints/30));
    makeHistogram(title+" - nearest distance histogram", pixelHeight, unit, nBins, NN);
    setBatchMode("show");

    // average spots and estimate diameter
    print("\n== Diameter estimate ==");
    selectImage(inputID);
    run("Duplicate...", "title=refAvg");
    refAvgImgID = getImageID();
    run("Gaussian Blur...", "sigma="+maxOf(1,sigmaGauss/2));

    caaResult = cropAndAverage(refAvgImgID,x,y,boxsize);
    avgImgID = caaResult[0];
    ncount = caaResult[1];
    print("   Number of spots averaged: "+ncount);
    selectImage(avgImgID);
    norm = getPixel(floor(boxsize/2), floor(boxsize/2));
    run("Divide...", "value="+norm);
    circleSize=boxsize;
    circleCentre = floor(boxsize/2);
    circleRadius = floor(boxsize/2);
    plotRadialProfile(circleCentre,circleCentre,circleRadius);
    Plot.getValues(xp, yp);
    close();
    selectImage(refAvgImgID);
    close();
    selectImage(avgImgID);
    close();

    yFinal = newArray(yp.length+1);
    yFinal[0]+=1.;
    for (i=0; i<yp.length; i++)
      yFinal[i+1]+=yp[i];

    xFinal = newArray(yp.length+1);
    xFinal[0] = 0;
    for (i=0; i<xp.length; i++)
       xFinal[i+1] = xp[i];

    //Plot.create("Average spot from "+title, "Distance (pixels)", "Intensity", xFinal, yFinal);

    lastY = 1;
    minY = 0;
    minX = -1;
    done = false;
    for (i=0; i<yFinal.length && !done; i++) {
        if (yFinal[i] > 1 || yFinal[i] <= lastY) lastY = yFinal[i];
        else {
            minY = yFinal[i-1];
            minX = xFinal[i-1];
            done=true;
        }
    }
    if (done!=true) minY = yFinal[yFinal.length-1];
    //targetFWHM = 0.5*(1.+minY);
    // diameter is estimated as the full width at 26% of the maximum heigth
    targetFWHM = minY+0.26*(1.-minY);
    done = false;
    for (i=0; i<yFinal.length && !done; i++) {
        if (yFinal[i] < targetFWHM) {
            xInterp = interpolate(xFinal[i-1],xFinal[i],yFinal[i-1],yFinal[i],targetFWHM);
            xFWHM = xInterp*2;
            done = true;
        }
    }

    finalDiameter=xFWHM;

    selectImage(inputID);
    run("Duplicate...", "title="+title+"_spots");
    imgMarkedID=getImageID();
    plotCircles(imgMarkedID,x,y,finalDiameter);
    run("Enhance Contrast", "saturated=0.35");
    setBatchMode("show");

    initialDiameter = initialDiameter*pixelHeight;
    finalDiameter = finalDiameter*pixelHeight;
    print("   Diameter: "+finalDiameter+" "+unit+"  [initial estimate by autocorrelation: "+initialDiameter+" "+unit+" ]");

    addMeasure(ResultsTable, outTitle, nPoints, area, finalDiameter, NNmean, NNstdDev, densityNN, NNI, pVal, unit);

    selectWindow("Results");
    run("Close");
    selectWindow("Log");
}
setBatchMode(false);
if(totImages>1) showMessage("Done!");


/*****************  FUNCTIONS ********************/

function initializeResultsTable(ResultsTable) {
/* 
 * create an empty table if not already open
 * 
 */

  if (!isOpen(ResultsTable)){
    run("Table...", "name=["+ResultsTable+"]");
    print("["+ResultsTable+"]", "\\Headings: Image \t N \t Area \t Diameter \t NearestDistanceMean \t NearestDistanceStdDev \t DensityByNearestDistance \t DensityByArea \t NNIndex \t p-val \t unit");
  }

}


function addMeasure(ResultsTable, title, nPoints, area, diameter, NNmean, NNstdDev, densityNN, NNI, pVal, unit) {
/* 
 * add a row to existing table with results from last image analysed
 * 
 */
	areaDensity = d2s(nPoints/area,3);
    print("["+ResultsTable+"]", title+"\t"+nPoints+"\t"+area+"\t"+diameter+"\t"+NNmean+"\t"+NNstdDev+"\t"+densityNN+"\t"+areaDensity+"\t"+NNI+"\t"+pVal+"\t"+unit+"\n");
	
}


function plotCircles(imgID,x,y,diameter) {
/*
 * Draw circles at given coordinates
 * 
 */

    roiManager("Reset");
    selectImage(imgID);
    radius=diameter/2;;
    for (i=0; i<x.length; i++) {
        makeOval(x[i]-radius, y[i]-radius, 2*radius, 2*radius);
        roiManager("Add");
    }
    roiManager("Show All without labels");
}


function interpolate(x1,x2,y1,y2,y0){
/*
 * Estimate x0 at given y0 by linear interpolation
 * 
 */

   return x1+(x2-x1)/(y2-y1)*(y0-y1);
}


function cropAndAverage(imgID,xcoords,ycoords,boxsize) {
/*
 * crop boxes around detected spots and average all of them
 * 
 */
    selectImage(imgID);
    run("Duplicate...", "title=tmpImage"+imgID);
    tmpImgID = getImageID();
    width=getWidth;
    height=getHeight;

    if (isOpen("averageImg") ){
        selectImage("averageImg");
        close();
    }
    if (isOpen("partImg") ){
        selectImage("partImg");
        close();
    }

    istart = 0;
    done = false;
    for ( i = 0; i< xcoords.length && done==false; i++){
        xstart = xcoords[i]-floor(boxsize/2);
        ystart = ycoords[i]-floor(boxsize/2);
        if (xstart>0 && ystart>0 && (xstart+boxsize)<width && (ystart+boxsize)<height) {
            makeRectangle(xstart,ystart,boxsize,boxsize);
            run("Duplicate...", "title=averageImg");
            run("32-bit");
            istart = i;
            done = true;
        }
    }

    ncount = 1;
    for ( i = istart+1; i< xcoords.length; i++){
        xstart = xcoords[i]-floor(boxsize/2);
        ystart = ycoords[i]-floor(boxsize/2);
        if (xstart>0 && ystart>0 && (xstart+boxsize)<width && (ystart+boxsize)<height) {
            selectImage(tmpImgID);
            makeRectangle(xstart,ystart,boxsize,boxsize);
            run("Duplicate...", "title=partImg");
            run("32-bit");
            imageCalculator("Add", "averageImg","partImg");
            selectImage("partImg");
            close();
            ncount++;
        }
    }
    selectImage(tmpImgID);
    close();
    selectImage("averageImg");
    run("Divide...", "value="+ncount);
    if (ncount < 50) print("WARNING: low number of spots ("+ncount+")! The estimate of the diameter could not be accurate");
    selectImage("averageImg");
    imgID = getImageID();
    results = newArray(imgID,ncount);
    return results;
}


function plotRadialProfile(X0,Y0,mR) {
/*
 * Plot radial profile within given radius
 * 
 * Adapted from Radial Profile Plot plugin (https://imagej.nih.gov/ij/plugins/radial-profile.html)
 * 
 */

	nBins = floor(3*mR/4);

	xmin=floor(X0-mR)+1;
	xmax=floor(X0+mR);
	ymin=floor(Y0-mR)+1;
	ymax=floor(Y0+mR);
	AccumX = newArray(nBins);
	AccumY = newArray(nBins);

	for (i=xmin; i<xmax; i++) {
		for (j=ymin; j<ymax; j++) {
			R = sqrt((i-X0)*(i-X0)+(j-Y0)*(j-Y0));
			thisBin = floor((R/mR)*nBins);
			if (thisBin==0) thisBin=1;
			thisBin=thisBin-1;
			if (thisBin>nBins-1) thisBin=nBins-1;
			AccumX[thisBin]=AccumX[thisBin]+1;
			AccumY[thisBin]=AccumY[thisBin]+getPixel(i,j);
		}
	}
	for (i=0; i<nBins;i++) {
		AccumY[i] = AccumY[i] / AccumX[i];
		AccumX[i] = mR*((i+1)/nBins);
	}
	Plot.create("Radial Profile Plot", "Radius [pixels]", "Normalized Integrated Intensity",  AccumX, AccumY);
	Plot.show();

}


function makeHistogram(title, pixel, unit, nBins, ArrayName){
/* 
 * plot histogram of NN distances from array
 * 
 */
    histogramBinArrayMin=newArray(nBins);
    histogramBinArrayMax=newArray(nBins);
    histogramBinArray=newArray(nBins);
    histogramPlotXArray=newArray(nBins);
    histogramPlotFreqArray=newArray(nBins);
    Array.getStatistics(ArrayName, ArrayMin, ArrayMax, ArrayMean, ArraystdDev);
    binSize=(ArrayMax-ArrayMin)/(nBins+1);

    //set the bin values
    for(i=0;i<nBins;i++){
        histogramBinArrayMin[i]=ArrayMin+(i*binSize);
        histogramBinArrayMax[i]=ArrayMin+((i+1)*binSize);
    }

    //determine the number plot profile values in each bin
    for(i=0;i<lengthOf(ArrayName);i++){
        chkE=ArrayName[i];
        for(j=0;j<nBins;j++){
            if(chkE>=histogramBinArrayMin[j] && chkE<histogramBinArrayMax[j]){
                chkBinCount=histogramBinArray[j];
                histogramBinArray[j]=chkBinCount+1;
            }
            if(chkE==ArrayMax){
                chkBinCount=histogramBinArray[nBins-1];
            }
        }
    }

    //generates the values for plotting
    for(i=0;i<nBins;i++){
        histogramPlotXArray[i]=(histogramBinArrayMin[i]+histogramBinArrayMax[i])/2;//for x axis center in the middle of the bin
        histogramPlotFreqArray[i]=histogramBinArray[i]/(lengthOf(ArrayName)); //for y axis is number in each bin divided by the total number of plot values
    }


    for(i=0;i<histogramPlotXArray.length;i++){
        histogramPlotXArray[i] = histogramPlotXArray[i] * pixel;
    }

    //make a plot of the histogram
    PlotAsColumns(title, "distance ("+unit+")", "Normalized count", histogramPlotXArray, histogramPlotFreqArray,binSize*pixel);
    Plot.addText(d2s(ArrayMean*pixel,2)+" "+fromCharCode(177)+" "+d2s(ArraystdDev*pixel,2)+" "+unit,0.05,0.1);
    Plot.show();
}


function PlotAsColumns(PlotName,xAxisText,yAxisText,xArray,yArray,columnWidth){
/* 
 * subroutine of makeHistogram to draw bars
 * 
 */

    Array.getStatistics(xArray,xmin,xmax,xmean,xstdDev);
    Array.getStatistics(yArray,ymin,ymax,ymean,ystdDev);

    Plot.create(PlotName, xAxisText, yAxisText);
    Plot.setLimits(xmin-columnWidth/2,xmax+columnWidth/2,0,ymax);
    for(i=0;i<xArray.length;i++){
        x1=xArray[i]-columnWidth/2;
        x2=xArray[i]+columnWidth/2;
        y1=yArray[i];
        y2=y1;
        Plot.drawLine(x1,y1,x2,y2);
        Plot.drawLine(x1,y1,x1,0);
        Plot.drawLine(x2,y2,x2,0);
    }
}


function HanningWindow(width,height) {
/*
 *  generate a 2D Hanning apodizing window
 *  
 */
    newImage("Hanning Window", "32-bit Black", width, height, 1);
    imgID = getImageID();
    for (v=0; v<height; v++) {
        for (u=0; u<width; u++) {
            ru=2*u/width-1;
            rv = 2*v/height-1;
            ruv=pow(ru,2)+pow(rv,2);
            ruv = sqrt(ruv);
            wuv =  0.5*(cos(PI*ruv)+1);
            if (ruv >= 0 && ruv < 1)
                setPixel(u, v, 1*wuv);
            else
                setPixel(u, v, 0);
        }
    }
    return imgID;
}