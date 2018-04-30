// ImageJ macro to analyse myofibril longitudinal sections
//
// Max Planck Institute of Biochemistry - Giovanni Cardone and Maria Spletter
//
/* Usage
 - analysis can be run on entire image or on a selected ROI (not recommended when analyzing thickness)
 - analysis is run on all the images opened
*/
/* Constraints
 *   the fibrils must be oriented horizontally
*/
/* Results
 * Parameters estimated:
 *   sarcomere length
 *   fibril thickness
 * Both measures are estimated in different ways
 * Sarcomere length/repeat:
 *   - A global estimate is obtained by measuring the periodicity in the Fourier Transform (FT) from the entire image.
 *   - A localized estimate is obtained similarly but after dividing the image in subregions and averaging their FT
 * Fibril thickness:
 *   - A global estimate is obtained by analysing the profile of the autocorrelation function calculated along the
 *     direction perpendicular to the fibrils
 *   - A localized estimated is obtained by looking for profiles perpendicular to fibrils and averaging their widths
*/

setBatchMode(true);

ResultsTable = "Fibril longitudinal sections";
// parameters on full size image
fulloriented = false;   // not implemented
fullapodized = false;
// parameters on subregions of image
oriented = false;       // not implemented
apodized = true;
pad = 2;
totImages = nImages;

if (totImages==0) {
	exit("No images open!");
}

initializeResultsTable(ResultsTable);

for (in=1; in<=totImages; in++){

    selectImage(in);
    imgID = getImageID();
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


    Stack.getDimensions(width, height, Nchannels, Nslices, Nframes);
    Stack.getPosition(currentChannel, currentSlice, currentFrame);
    if (doSelection==true) {
        getSelectionBounds(x, y, width, height);
    }
    Stack.getPosition(c,s,f);
    getPixelSize(unit, pixelWidth, pixelHeight);

    imgName = title;
    if (Nchannels>1)  imgName = imgName+"Ch"+currentChannel;
    if (Nslices>1)  imgName = imgName+"Sl"+currentSlice;
    if (doSelection==true) imgName +="-crop";
    print(" ");
    print("--- "+imgName+" ---");

    // work on duplicated image to restrict analysis to selection, if any
    tempTitle="temp-"+random();
    run("Duplicate...", "title="+tempTitle);
    tempID=getImageID();
    if (fulloriented == true) {
        rotAngle = orientImage(tempID);
        print("Initial rotation: "+rotAngle);
    }
    // calculate FFT from full image
    selectImage(tempID);
    if (fullapodized == true) {
        tempInputFFTID = apodizeHanning(tempID);
        selectImage(tempInputFFTID);
    }
    run("FFT");
    //run("32-bit");
    rename(imgName+"-FFT");
    fftID=getImageID();
    if (fullapodized == true) {
    	// remove intermediate image
        selectImage(tempInputFFTID);
        close();
    }
    globalRepeat = findRepeat(fftID);
    if (globalRepeat < 0){
        exit("Error: can not determine the repeat!");
    }
    selectImage(fftID);
    fftSize=getWidth();
    drawRepeatSpots(fftID,globalRepeat);
    physicalGlobalRepeat = (fftSize*pixelWidth)/globalRepeat;
    setBatchMode("show");

    print("\n== Repeat ==");

    print("   Estimate on entire image/selection: "+physicalGlobalRepeat+" ("+unit+")");

    // repeat value used for thickess estimation
    referenceRepeat = physicalGlobalRepeat/pixelWidth;

    // set parameters for boxwise analysis
    boxX = physicalGlobalRepeat/pixelWidth*8;
    boxY = floor(boxX/2)+1;

    fftResults = periodogramAverage(tempID,imgName,boxX,boxY,oriented,apodized);
    fftAvgID = fftResults[0];
    NfftAvg = fftResults[1];

    if (fftAvgID<0){
        selectImage(fftAvgID);
        setBatchMode("show");

        refinedRepeat = findRepeat(fftAvgID);
        if (refinedRepeat < 0){
            exit("Error: can not determine the repeat!");
        }

        selectImage(fftAvgID);
        fftSize=getWidth;
        drawRepeatSpots(fftAvgID,refinedRepeat);
        physicalRefinedRepeat = (fftSize*pixelWidth)/refinedRepeat;
        setBatchMode("show");

        print("   Localized estimate: "+physicalRefinedRepeat+" ("+unit+")");
        print("      Subregions averaged: "+NfftAvg);

    } else {
        refinedRepeat = -1;
        physicalRefinedRepeat = -1;
    }

    results = estimateThickness(tempID, imgName, referenceRepeat);
    thicknessGlobalACF = results[0];
    thicknessRealSpace = results[1];
    thicknessRealSpaceStd = results[2];
    Ntrs = results[3];

    selectImage(tempID);
    close();

    addMeasure(ResultsTable, imgName, physicalGlobalRepeat, physicalRefinedRepeat, NfftAvg, thicknessGlobalACF, thicknessRealSpace, thicknessRealSpaceStd, Ntrs, unit);
}

setBatchMode(false);
if(totImages>1) showMessage("Done!");


/*****************  FUNCTIONS ********************/

function periodogramAverage(imgID,refTitle,boxX,boxY,oriented,apodized) {
/*
 * 	generate mean FFT from image subregions by periodogram averaging 
 *  
 */
    // generate 2D Hanning window
    boxHX = boxX*pad;
    boxHY = boxY*pad;
    newImage("Hanning Window", "32-bit Black", boxHX, boxHY, 1);
    for (v=0; v<boxHY; v++) {
        for (u=0; u<boxHX; u++) {
            ru=2*u/boxHX-1;
            rv = 2*v/boxHY-1;
            ruv=pow(ru,2)+pow(rv,2);
            ruv = sqrt(ruv);
            wuv =  0.5*(cos(PI*ruv)+1);
            if (ruv >= 0 && ruv < 1)
                setPixel(u, v, 1*wuv);
            else
                setPixel(u, v, 0);
        }
    }

    selectImage(imgID);
    run("Select None");
    run("Duplicate...", "title=tmpImage"+imgID);
    tmpImgID = getImageID();
    width=getWidth;
    height=getHeight;

    if (isOpen("averageFFTImg") ){
        selectImage("averageFFTImg");
        close();
    }
    if (isOpen("boxFFTImg") ){
        selectImage("boxFFTImg");
        close();
    }
    if (isOpen("boxImg") ){
        selectImage("boxImg");
        close();
    }

    // check if multiple subregions can be extracted with current settings
    multiBoxes = ((width-boxX)>0) && ((height-boxY)>0);

    if (!multiBoxes) {
        print("Selection is too small. Skipping localized estimate of the repeat");
        selectImage(tmpImgID);
        close();
        selectImage("Hanning Window");
        close();
        return newArray(1,0);
    }
    // generate an empty image, used to accumulate the FFT
    makeRectangle(0,0,boxX,boxY);
    run("Duplicate...", "title=boxImg");
    if (pad>1) run("Canvas Size...", "width="+pad*boxX+" height="+pad*boxY+" position=Center zero");
    boxImgID = getImageID();
    run("FFT");
    run("32-bit");
    rename("averageFFTImg");
    run("Set...", "value=0");
    selectImage("boxImg");
    close();
    // generate average FFT image
    ncount = 0;
    weight = 0.0;
    for (x = 0; x< width-boxX; x+=boxX/2){
        for (y = 0; y< height-boxY; y+=boxY/2){
            selectImage(tmpImgID);
            makeRectangle(x,y,boxX,boxY);
            run("Duplicate...", "title=boxImg");
            boxID = getImageID();
            run("32-bit");
            getRawStatistics(nPixels, mean);
            if (oriented==true){
                rotAngle = orientImage(boxID);
            }
            if (pad>1) run("Canvas Size...", "width="+pad*boxX+" height="+pad*boxY+" position=Center zero");
            // smooth edges before FFT
            if (apodized == true) imageCalculator("Multiply 32-bit", "boxImg", "Hanning Window");
            run("FFT");
            run("32-bit");
            rename("boxFFTImg");
            run("Multiply...", "value="+mean);
            // accumulate FFTs
            imageCalculator("Add", "averageFFTImg","boxFFTImg");
            selectImage("boxImg");
            close();
            selectImage("boxFFTImg");
            close();
            ncount++;
            weight+=mean;
        }
    }
    selectImage(tmpImgID);
    close();
    if (ncount == 0) {
        print("Could not extract subregions for estimating the repeat");
        selectImage("Hanning Window");
        close();
        return newArray(1,0);
    }
    // normalize accumulated FFT to obtain average
    selectImage("averageFFTImg");
    run("Divide...", "value="+weight);
    selectImage("averageFFTImg");
    rename(refTitle+"-FFTavg");
    run("Enhance Contrast", "saturated=0.35");

    avgImgID = getImageID();

    selectImage("Hanning Window");
    close();

    return newArray(avgImgID,ncount);
}


function findRepeat(fftID){
/*
 * Determine periodicity in FTT image by analyzing average line profile along
 * the X axis and finding the repeat best in agreement with the profile
 * 
*/

    selectImage(fftID);
    fftSize=getWidth;

    // Line profile is obtained as average from multiple lines
    // Set initial values for crop region and optimize the width
    boxWidth = fftSize/4;
    boxStartX = fftSize/2;
    maxBoxHeight = floor(fftSize/2);
    boxHeight = findOptimumProfileThickess(fftID,boxStartX,boxWidth,maxBoxHeight);
    boxStartY = fftSize/2 - boxHeight/2;
    selectImage(fftID);
    makeRectangle(boxStartX, boxStartY, boxWidth, boxHeight);
    profile = getProfile();

    // find periodicity by testing all plausible values
    minRepeat = 3;
    maxRepeat = floor((profile.length-1)/1.5);
    bestMetric = 0;
    repeat = -1;
    // The metric is calculated, for each repeat value, as the sum of the
    // intensity difference between each expected peak position and
    // its closest expected minimum position (the one with lower intensity)
    // Negative differences are excluded from the sum
    for (r=minRepeat; r<=maxRepeat; r++){
        metric = 0.;
        halfRepeat = floor(r/2);
        imax = floor((profile.length-1)/(1.5*r));
        for (i=1; i<=imax; i++) {
            if ((profile[i*r]-maxOf(profile[i*r-halfRepeat],profile[i*r+halfRepeat]))>0){
                metric += (profile[i*r]-minOf(profile[i*r-halfRepeat],profile[i*r+halfRepeat]));
            }
        }
        if (metric>bestMetric) {
            bestMetric = metric;
            repeat = r;
        }
    }
    return repeat;
}


function orientImage(imgID) {

    //function not yet implemented

    return(0);
}


function apodizeHanning(imgID){
/*
 * apply Hanning window function to image to dampen values close to the edges
 * 
*/

	selectImage(imgID);
    title = getTitle();
    boxX = getWidth;
    boxY = getHeight;
    // generate 2D Hanning function
    newImage("Hanning Window", "32-bit Black", boxX, boxY, 1);
    for (v=0; v<boxY; v++) {
        for (u=0; u<boxX; u++) {
            ru=2*u/boxX-1;
            rv = 2*v/boxY-1;
            ruv=pow(ru,2)+pow(rv,2);
            ruv = sqrt(ruv);
            wuv =  0.5*(cos(PI*ruv)+1);
            if (ruv >= 0 && ruv < 1)
                setPixel(u, v, 1*wuv);
            else
                setPixel(u, v, 0);
        }
    }
    // apply apodization
    imageCalculator("Multiply 32-bit", title, "Hanning Window");
    newID = getImageID();
    selectImage("Hanning Window");
    close();
    return(newID);
}


function findOptimumProfileThickess(imgID,centerCoord,boxWidth,maxBoxHeight) {
/*
 * ideal thickness is defined as the one with least maxima (to reduce ambiguity detection)
 * 
*/


    selectImage(imgID);
    boxStartX = centerCoord;
    minMaxima = boxWidth;
    optHalfHeight = 1;
    for (y=1; y<=floor(maxBoxHeight/2); y++){
        boxStartY = centerCoord -y;
        boxHeight = 2*y+1;
        makeRectangle(boxStartX, boxStartY, boxWidth, boxHeight);
        profile = getProfile();
        nMaxima = 0;
        for (p=1; p<profile.length-1; p++){
            if( (profile[p-1] <= profile[p] && profile[p] >= profile[p+1])){
                nMaxima++;
            }
        }
        if (nMaxima < minMaxima){
            minMaxima = nMaxima;
            optHalfHeight = y;
        }
    }
    return 2*optHalfHeight+1;

}


function drawRepeatSpots(fftID,repeat){
/*
 * draw spots over the expected line peaks used to quantify the repeat
 * 
*/
    radius = repeat / 3;
    selectImage(fftID);
    run("8-bit");
    fftsize = getWidth;
    // crop the image to a significant size
    boxHeight = floor(fftsize/8)*2;
    boxWidth = fftsize/2;
    boxStartX = floor(0.5*(fftsize - boxWidth))+1;
    boxStartY = floor(0.5*(fftsize - boxHeight))+1;
    makeRectangle(boxStartX, boxStartY, boxWidth, boxHeight);
    run("Crop");
    centerX = boxWidth/2-radius-1;
    centerY = boxHeight/2-radius-1;
    Nrepeats = floor(boxWidth/2/repeat);
    for (i=0; i<=Nrepeats; i++) {
        cx = centerX + i*repeat;
        Overlay.drawEllipse(cx, centerY, 2*radius+1, 2*radius+1);
        cx = centerX - i*repeat;
        Overlay.drawEllipse(cx, centerY, 2*radius+1, 2*radius+1);
    }
    Overlay.show
}


function estimateThickness(imgID, refTitle, period){
/*
 * Estimate thickness of fibrils both in Fourier space (AutoCorrelation)
 * and in real space
 * 
 */
 
    selectImage(imgID);
    title = getTitle();
    width = getWidth;
    height = getHeight;
    refSize = maxOf(width,height);

    getPixelSize(unit, pixelWidth, pixelHeight);

	// Method 1:
	// Estimate thickness by analyzing the autocorrelation line profile along Y
	// for 4 subregions around the centre (from full size to 1/8 size)
	// and selecting the smallest estimate
	// The thickness is determined as the firt minimum in the autocorrelation function
    minAutoCorrSize = minOf(64,maxOf(width, height));
    tmpWidth = width;
    tmpHeight = height;
    optimalThickness = maxOf(width,height);
    for (k=1; k<4 && maxOf(tmpHeight,tmpWidth)>=minAutoCorrSize; k++) {

        selectImage(imgID);
        minFFTsize=maxOf(tmpWidth, tmpHeight);
        size=4;
        while(size<minFFTsize) size*=2;

        tempTitle="temp-"+random();
        run("Duplicate...", "title="+tempTitle);
        tempID=getImageID();
        getRawStatistics(nPixels, mean);
        if (!(size==width && size==height)) {
            run("Canvas Size...", "width="+ size+" height="+ size+" position=Center zero");
            makeRectangle(floor((size-tmpWidth)/2), floor((size-tmpHeight)/2), tmpWidth, tmpHeight);
            run("Make Inverse");
            run("Set...", "value="+mean);
            run("Select None");
        }
        //compute autocorrelation of particle image
        run("FD Math...", "image1=["+tempTitle+"] operation=Correlate image2=["+tempTitle+"] result=AutoCorrelation do");
        ACimgID=getImageID();
        run("Subtract...", "value="+(nPixels*mean*mean));
        selectImage(tempID);
        close();
        // plot profile along Y and determine the first minimum
        selectImage(ACimgID);
        makeLine(size/2, size/2, size/2, 1);
        profile = getProfile();

        selectImage(ACimgID);
        close();
        minX= findFirstMinimum(profile);

        if (minX>0 && minX<optimalThickness){
            optimalThickness = minX;
            bestprofile = Array.copy(profile);
        }

        tmpHeight/=2;
        tmpWidth/=2;
    }

    if (optimalThickness == maxOf(width,height)) {
        exit("Error: can not determine the thickness!");
    }

    print("\n== Thickness ==");
    physicalOptimalThickness = optimalThickness*pixelWidth;
    print("   Estimate by autocorrelation: "+physicalOptimalThickness+" ("+unit+")");

    xCoord = newArray(bestprofile.length);
    for (k=0; k<xCoord.length; k++) {
        xCoord[k] = k*pixelWidth;
    }
    Plot.create(refTitle+" - thickness by autocorrelation profile", unit, "a.u.", xCoord, bestprofile);
    Plot.setColor("red");
    if (optimalThickness>0) {
        Plot.drawLine(xCoord[optimalThickness], bestprofile[0], xCoord[optimalThickness], bestprofile[optimalThickness]);
        Plot.setColor("blue");
    }
    Plot.show();
    setBatchMode("show");


	// Method 2: estimate thickness by dividing the image in optimal subregions
	// and analyzing the intensity profile after projecting along the longitudinal
	// direction. Thickness of single fibrils is estimated as the half distance 
	// between two intensity minima
    boxW = minOf(width,1*period);
    boxH = minOf(maxOf(5*optimalThickness,height/2),height);
    nBoxesX = maxOf(1,floor(width/boxW));
    nBoxesY = maxOf(1,floor(height/boxH));

    thickArray = newArray();

    for (i=0; i<nBoxesX; i++){
        for (j=0; j<nBoxesY; j++){
            selectImage(imgID);
            makeRectangle(i*boxW, j*boxH, boxW, boxH);
            run("Duplicate...", "title=boxedImg");
            boxID=getImageID();
            getRawStatistics(nPixels, mean);
            run("Rotate 90 Degrees Right");
            run("Gaussian Blur...", "sigma=1");
            run("Select All");
            profile = getProfile();
            maxPos = newArray();
            minPos = newArray();
            // collect all maxima and minima
            for (k=1; k<profile.length-1; k++){
                if( (profile[k-1] <= profile[k] && profile[k] >= profile[k+1])){
                    maxPos = Array.concat(maxPos,k);
                }
                if( (profile[k-1] >= profile[k] && profile[k] <= profile[k+1])){
                    minPos = Array.concat(minPos,k);
                }
            }

            // Scan maxima and if peak is above the average intensity in the original subregion
            // then add half distance between previous and next minimum (only if both below average)
            // to the thickess estimates
            for (mx=0; mx<maxPos.length; mx++){
                if (profile[maxPos[mx]] > mean){
                    done = false;
                    for (mn=1; mn<minPos.length && !done; mn++){
                        if (minPos[mn]>maxPos[mx]){
                            done = true;
                            if (profile[minPos[mn]]<mean && profile[minPos[mn-1]]<mean){
                                 thickArray = Array.concat(thickArray,(minPos[mn]-minPos[mn-1])*0.5);
                            }
                        }
                    }
                }
            }
            selectImage(boxID);
            close();
        }
    }

    // Determine the thickness from the average of all the measures
    if (thickArray.length>0) {
	    Array.getStatistics(thickArray, tmin, tmax, tmean, tstdDev);
    	realspaceThickness = tmean;
    	realspaceThicknessError = tstdDev;

    	physicalRealspaceThickness = realspaceThickness*pixelWidth;
    	physicalRealspaceThicknessError = realspaceThicknessError*pixelWidth;
    } else {
    	physicalRealspaceThickness = -1.;
    	physicalRealspaceThicknessError = -1.;
    }
    print("   Estimate by averaging: "+physicalRealspaceThickness+" +/- "+physicalRealspaceThicknessError+" ("+unit+")");
    print("      Number of elements averaged: "+thickArray.length);

    results = newArray(physicalOptimalThickness, physicalRealspaceThickness, physicalRealspaceThicknessError, thickArray.length);
    return results;

}


function findFirstMinimum(profile) {
/*
 * Determine position of first minimum in profile curve
 * 
 */
    lastY = profile[0];
    minY = 0;
    minX = -1;
    done = false;
    for (i=1; i<profile.length && !done; i++) {
        if (profile[i] > profile[0] || profile[i] < lastY) {lastY = profile[i];}
        else {
            minY = profile[i-1];
            minX = i-1;
            done=true;
        }
    }

    return minX;
}


function initializeResultsTable(ResultsTable) {
/* 
 * create an empty table if not already open
 * 
 */

  if (!isOpen(ResultsTable)){
	  run("Table...", "name=["+ResultsTable+"]");
	print("["+ResultsTable+"]", "\\Headings: Image \t repeatByFullFFT \t repeatByLocalFFT \t Ncrops \t thicknessByAutoCorrelation \t thicknessByProfileMean \t thicknessByProfileStdDev \t Nprofiles \t unit");
  }

}


function addMeasure(ResultsTable, title, repeatGlobalFFT, repeatLocalFFT, Nrl, thicknessGlobalACF, thicknessRealSpace, thicknessRealSpaceStd, Ntrs, unit) {
/*
 * add measure(s) to the results table
 */

  print("["+ResultsTable+"]", title+"\t"+repeatGlobalFFT+"\t"+repeatLocalFFT+"\t"+Nrl+"\t"+thicknessGlobalACF+"\t"+thicknessRealSpace+"\t"+thicknessRealSpaceStd+"\t"+Ntrs+"\t"+unit+"\n");

}
