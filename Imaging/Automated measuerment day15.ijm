// Set the folder path where the images are stored
folderPath = "V:/AG Thiel/Mitarbeiter/Alexander/LAMA5/CRISPR/USC (LAMA5 KO)/KOs/LAMA5 KO C Nr.2 (FACS) - Erfolgreich/SpheroidMessung/Tag 15/KO3/";

// Get the list of all files in the folder
list = getFileList(folderPath);

// Loop through all files in the folder
for (i = 0; i < list.length; i++) {
    file = list[i];

    // Check if the file is a .tiff image (you can add other image formats if needed)
    if (endsWith(file, ".tiff")) {
        // Open the image
        open(folderPath + file);
        selectWindow(file);

        // Set the scale (distance=1 µm, known=0.37 for pixel size)
        run("Set Scale...", "distance=1 known=0.37 unit=µm");

        // Adjust contrast (set minimum to 90 and maximum to 100)
        setMinAndMax(30, 70);

        // Set a fixed threshold range (100-255)
        setThreshold(55, 255);

        setOption("BlackBackground", true);
        run("Convert to Mask");
        run("Despeckle");

        // Analyze particles and show the outlines in a separate window
        run("Analyze Particles...", "size=5000-Infinity circularity=0.1-1.00 show=Outlines display add");

        // Close the current image after processing
        close();
    }
}
