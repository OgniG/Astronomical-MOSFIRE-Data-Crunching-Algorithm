Version 1 8/13/16


1. Disclaimer: This program only works for MOSFIRE data from the Keck Archives although it may be easily edited to work for other instruments
2. Requirements:
   1. The following python packages installed: __future__, sys, decimal, os, os.path, itertools, re, csv, numpy, urllib2, gzip, astropy, datetime, dateutil.parser, MOSFIRE, glob, aplpy, warnings, tabulate, and matplotlib
   2. The program needs to be in the same directory as the FITS files
   3. An internet connection
   4. A python environment (or command line)
   5. Available storage space for tables to be saved in
   6. Current version of MOSFIRE Data Reduction Pipeline which can be found: https://keck-datareductionpipelines.github.io/MosfireDRP/ .
1. INSTRUCTIONS
   1. Download MASTER.py
   2. Download appropriate FITS files (if necessary)
   3. Double-click MASTER.py
   4. Follow on-screen instructions
1. PROGRAM DESIGN
   1. Feature 1: PI (Principal Investigator) Info: User inputs url link of ASCII table from Keck Archive and program outputs list of PIs, number of PIs, and list of projects under a selected PI.
   2. Feature 2: Project Info: User inputs url link of ASCII table from Keck Archive and program outputs list of projects, number of projects, number and listing of masks used in a selected project, and number of nights assigned to a selected project.
   3. Feature 3: Mask(Target) Info: User inputs url link of ASCII table from Keck Archive and program outputs the total exposure time of a selected target. 
   4. Feature 4: Summary Table Generation:
      1. Generation of a summary table without filter info: user inputs url link of ASCII table from Keck Archive and program outputs a summary table of all entries in ASCII table with info about PI, project, mask, and exposure time.
      2. Generation of a summary table with filter info based on FITS headers downloaded in directory: program reads all FITS headers the user has downloaded (whether compressed or uncompressed) and outputs a summary table with info about PI, project, mask, filter, and exposure time; ideal for users who do not have access to all FITS headers and who want to only look at images taken under a specific PI, project, or mask. 
   1. Feature 5: Object Search: 
      1. User must have relevant FITS files downloaded
      2. Categorizes objects from the files by mask and presents the information in tables which may be saved as .csv files (can be opened by Microsoft Excel).
      3. Allows user to search for specific objects, Right Ascensions, Declinations, file names, mask names, and dates and results are presented in tables which may be saved in the .csv format
   1. Feature 6: Plot of mask footprints:
      1. User must have: FITS image of field and desired FITS headers downloaded.
      2. Optional: 3-color image of field with exact coordinates of FITS image of field, FITS file of field contours
      3. Program reads FITS headers and takes coordinates of each mask and stores coordinates in ‘footcoord.txt’ file, then generates a plot of mask footprints in field. 
      4. Plot has basic aesthetics with grayscale/colorscale/3-color image as background, grid, and coordinate labels. To add more artistic features to plot, user can look to the documentation for the APLpy package. 
   1. Features 1, 2, 3 need only the link to ASCII table from Keck Archive. Feature 4 can generate summary table with either a link of ASCII table (for summary table without filter info) or downloaded FITS headers (for summary table with filter info). Feature 5 requires downloaded FITS headers. Feature 6 requires downloaded FITS headers and FITS image of field.
1. AREAS FOR FUTURE DEVELOPMENT
   1. Plot mask footprints with rotation angle information
   2. Object Catalogue Matching
   3. Add a GUI
   4. Work for instruments other than MOSFIRE (Tip: search for “edit the following line for switching instruments” in the script)
   5. Make the program more “flexible” to accommodate for typos, etc