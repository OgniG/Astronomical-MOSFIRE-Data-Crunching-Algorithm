#!/usr/bin/env python

from __future__ import division

print ("\n" + "Initializing" + "\n")

import sys
import decimal
import os
import os.path
import itertools as IT
import re
import csv
import numpy as np
import urllib2
import gzip
from astropy.table import Table, Column
from astropy.io import ascii
from astropy.extern.six.moves.urllib import request
from itertools import izip
from itertools import islice
from os import listdir
from os.path import isfile, join
from astropy.io import fits
from astropy.time import Time
from datetime import datetime
from dateutil.parser import parse
from astropy.coordinates import ICRS
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
import MOSFIRE
from MOSFIRE import IO, Options
import glob as glob
import aplpy
import matplotlib
import warnings
from tabulate import tabulate

def main_option():
    print "Press 1 for information about PIs: "
    print "Press 2 for information about projects: "
    print "Press 3 for information about specific mask/target: "
    print "Press 4 to generate a summary table: "
    print "Press 5 for Object Information or to Search: "
    print "Press 6 for a plot of mask footprints: "
    print "Press 0 to Quit: "
    x = raw_input("->")
    return x 

def pi_info():   
        pianswer = raw_input('Press 1 for a list of PIs\nPress 2 for the number of PIs\nPress 5 to quit:\n->')
        if pianswer == '1':
            print uniquepilist
            listorquit = raw_input('Press 1 for a list of the projects under a certain PI\nPress 5 to quit:\n->')
            if listorquit == "1":
                pichoice = raw_input('Enter the name of the PI whose projects you want listed: ')
                if pichoice in pilist:
                    piprojlist = []
                    for i, j in enumerate(pilist):
                        if j == pichoice:
                        	#edit the following line for switching instruments
                            pidatalist = info['progtitl']
                            newlist = pidatalist[i]
                            piprojlist.append(newlist)
                    uniquelist = list(set(piprojlist))
                    print uniquelist
                    numorquit = raw_input('Press 1 for the number of projects under this PI\nPress 5 to quit:\n->')
                    if numorquit == "1":
                        print len(uniquelist)
                    elif numorquit == "5":
                        sys.exit()
                else:
                    print('That PI is not part of this database. Please try again')
            if listorquit == "5":
                sys.exit()
        elif pianswer == '2':
            print len(uniquepilist)
        elif pianswer == "5":
            sys.exit()

def project_info():
        projanswer = raw_input('Press 1 for a list of projects\nPress 2 for the number of projects\nPress 5 to quit:\n->')
        if projanswer == "1":
            print uniqueprojlist
            maskanswer = raw_input('Press 1 for the masks used in a project\nPress 2 for the number of nights assigned to a project\nPress 5 to quit\n->')
            if maskanswer == "1":
                projname = raw_input('Enter the name of the project of which you want the masks used listed: ')
                if projname in projlist:
                    masklist = []
                    for i, j in enumerate(projlist):
                        if j == projname:
                        	#edit the following line for switching instruments
                            plist = info['targname']
                            newlist = plist[i]
                            masklist.append(newlist)
                    uniquelist = list(set(masklist))
                    print uniquelist
                    numorquit = raw_input('For the number of masks used, Press 1\nTo quit, Press 5:\n->')
                    if numorquit == "1":
                        print len(uniquelist)
                        sys.exit()
                    elif numorquit == "5":
                        sys.exit()
                else:
                    print('That project is not part of this database. Please try again')
            elif maskanswer == "2":
                projname = raw_input('Enter the name of the project of which you want the to know the number of nights assigned: ')
                if projname in projlist:
                    masklist = []
                    for i, j in enumerate(projlist):
                        if j == projname:
                        	#edit the following line for switching instruments
                            plist = info['date_obs']
                            newlist = plist[7]
                            masklist.append(newlist)
                    uniquelist = list(set(masklist))
                    print (len(uniquelist))
                else:
                    print('That project is not part of this database. Please try again')
            elif maskanswer == "5":
                sys.exit()
        elif projanswer == "2":
            print len(uniqueprojlist)
        elif projanswer == "5":
            sys.exit()

def mask_info():
                targetname = raw_input('Enter the name of the specific mask(target name) of which you want to know the total exposure time: ')
                with open('data.txt') as f: #use data from txt file
                    lines = f.readlines()
                    totaltimes = []
                    for line in lines:
                      if targetname in line:
                      	#edit the following line for switching instruments
                        totaltimes.append(line.split()[9]) #append value from elapsedtime/exposuretime column (index 11)
                f.close()
                totaltimesfloat = [] #convert items in totaltimes list from strings to floats
                for item in totaltimes:
                   totaltimesfloat.append(float(item))
                print sum(totaltimesfloat), "seconds"
                sys.exit()


def summary_table():
    answer = raw_input('Press 2 to generate summary table WITHOUT filter info\nPress 3 to generate summary table based on files downloaded WITH filter info\n->')
    if answer == '2':
       create_table2()
    elif answer == '3':
       create_table3()

def create_table2(): #generates complete table based on ASCII without filter info
    print ("Initializing")
    Tfull = []
    Efull = []
    for i in range(0, len(pilist)):
       pi = pilist[i]
       project_name = projlist[i]
       mask = masklist[i]
       exposure_time = extimelist[i]
       y = [pi, project_name, mask]
       Tfull.append(y)
       Efull.append(exposure_time)
    T = []
    E = []
    for j in range(0, len(Tfull)):
       state = False
       for i in range(0, len(T)):
          if Tfull[j] == T[i]:
              E[i] += Efull[j]
              state = True
              break
       if state == False:
          T.append(Tfull[j])
          E.append(Efull[j])
    for h in range(0, len(T)):
       T[h].append(E[h])
    #edit the following line for switching instruments
    print tabulate(T, headers=["PI", "Project Name", "Mask/Target Name", "Exposure Time (seconds)"])
    
def create_table3(): #generates table based on files available
	print ("Initializing")
	xdata_rows = []
	data_rows = []
	filedata = []
	totalexps = []
	dupes = []
	reallist = []
	tot_exps = []
	test_lists = []
	filenames = [f for f in listdir(os.getcwd()) if isfile(join(os.getcwd(), f))]
	for x in filenames:
		f = open(x)
		file_name = str(f.name)
		posofx = filenames.index(x) + 1
		lenlist = len(filenames)
		print ("Reading " + str(posofx) + " of " + str(lenlist) + " files: " + str(file_name))
		if posofx == lenlist:
			print("\n")				 
		if ".fits" in x or ".fits.gz" in x:
		   hdulist = fits.open(file_name)
		   try:
		   	#edit the following line for switching instruments
			   filepi = hdulist[0].header['PROGPI']
		   except:
			   KeyError
		   try:
		   	#edit the following line for switching instruments
			   fileprojname = hdulist[0].header['PROGTL1']
		   except:
			   KeyError
		   try:
		   	#edit the following line for switching instruments
			   filemask = hdulist[0].header['TARGNAME']
		   except:
			   KeyError
		   try:
		   	#edit the following line for switching instruments
			   filefilter = hdulist[0].header['MF2NAME']
		   except:
			   KeyError
		   try:
		   	#edit the following line for switching instruments
			   filestartexp = hdulist[0].header['TIME-OBS']
		   except:
			   KeyError
		   try:
		   	#edit the following line for switching instruments
			   fileendexp = hdulist[0].header['TIME-END']
		   except:
			   KeyError
		   try:
		   	#edit the following line for switching instruments
			   filedate = hdulist[0].header['DATE-OBS']
		   except:
			   KeyError
		   startexp = [[x.strip() for x in filestartexp.split(':')]]
		   endexp = [[x.strip() for x in fileendexp.split(':')]]
		   firstelement = int(endexp[0][0]) - int(startexp[0][0])
		   secondelement = int(endexp[0][1]) - int(startexp[0][1])
		   thirdelement = float(endexp[0][2]) - float(startexp[0][2])
		   if secondelement < 0:
			   firstelement = firstelement - 1
			   secondelement = secondelement + 60
		   if thirdelement < 0:
			   secondelement = secondelement - 1
			   thirdelement = thirdelement + 60
		   secondelement = round(secondelement, 3)
		   thirdelement = round(thirdelement, 3)
		   firstelement = round(firstelement, 3)
		   totalexp = [int(firstelement), int(secondelement), thirdelement]
		   totalexps.extend([totalexp])
		   fileinfo = [filedate, filepi, fileprojname, filemask, filefilter]
		   xdata_rows.append(fileinfo)
	for x in xdata_rows:
		start_at = -1
		locs = []
		while True:
			try:
				loc = xdata_rows.index(x,start_at+1)
			except ValueError:
				break
			else:
				locs.append(loc)
				start_at = loc
		dupes.append(locs)
	uniq_dupes = list(set(map(tuple, dupes)))
	udupelist = [list(elem) for elem in uniq_dupes]
	for i in xrange(0, len(udupelist)):
		exec("price%d = %s" % (i + 1, repr(udupelist[i]))) in globals(), locals()
	for w in range(0, len(udupelist)):
		addingnums = udupelist[w]
		testlist = []
		for q in addingnums:
			numslist = totalexps[q]
			testlist.append(numslist)
		reallist.append(testlist)
	for i in xrange(0, len(reallist)):
		exec("price%d = %s" % (i + 1, repr(reallist[i]))) in globals(), locals()
	for i in range(0, len(reallist)):
		firstelement = sum(i[0] for i in reallist[i])
		secondelement = sum(i[1] for i in reallist[i])
		thirdelement = sum(i[2] for i in reallist[i])
		if thirdelement > 60:
			elem = thirdelement//60
			thirdelement = thirdelement - (elem * 60)
			secondelement = secondelement + elem
		if secondelement > 60:
			elem = secondelement//60
			secondelement = secondelement - (elem * 60)
			firstelement = firstelement + elem
		secondelement = round(secondelement, 3)
		thirdelement = round(thirdelement, 3)
		firstelement = round(firstelement, 3)
		thelist = [[int(firstelement), int(secondelement), thirdelement]]
		tot_exps.extend(thelist)
	for i in xrange(0, len(tot_exps)):
		exec("price%d = %s" % (i + 1, repr(tot_exps[i]))) in globals(), locals()
	u_list = [list(x) for x in set(tuple(x) for x in xdata_rows)]
	for p in xrange(0, len(tot_exps)):
		if tot_exps[p][0] < 10:
			tot_exps[p][0] = "0" + str(tot_exps[p][0])
		if tot_exps[p][1] < 10:
			tot_exps[p][1] = "0" + str(tot_exps[p][1])
		exposure_lists = [str(tot_exps[p][0]) + ":" + str(tot_exps[p][1]) + ":" + str(tot_exps[p][2])]
		test_lists.extend(exposure_lists)
	for i in xrange(0, len(u_list)):
		exec("price%d = %s" % (i + 1, repr(u_list[i]))) in globals(), locals()
	for i in xrange(0, len(test_lists)):
		exec("price%d = %s" % (i + 1, repr(test_lists[i]))) in globals(), locals()
	for d in range(0, len(test_lists)):
		u_list[d].extend([test_lists[d]])
	def my_parse(lis):
		try: 
			return parse(lis[0])
		except ValueError:
			return datetime(1, 1, 1)
	n_list = sorted(u_list, key=my_parse)
	#edit the following line for switching instruments
	summtab = Table(rows=n_list, names=('Date', 'PI', 'ProjectName', 'Mask/TargetName', 'Filter', 'ExposureTime(sec)'))
	print summtab
	newans = raw_input("Do you want to save these tables? Press 1 for yes. Press 2 for no: ")
        if newans == "1":
            with open('table.txt', 'w') as file:
                file.writelines('\t'.join(i) + '\n' for i in n_list)
            with open('table.txt','r') as infile, open('table.csv', 'wb') as outfile:
                in_txt = csv.reader(infile, delimiter = '\t')
                out_csv = csv.writer(open('table.csv', 'wb'))
                out_csv.writerows(in_txt)
                print ("Saved" + "\n")
	
def grab_dimensions():
   print ("\n" + "Initializing" + "\n")
   ralist = []
   declist = []
   filenames = [f for f in listdir(os.getcwd()) if isfile(join(os.getcwd(), f))]
   for x in filenames:
	  f = open(x)
	  file_name = str(f.name)
	  if ".fits" in x or ".fits.gz" in x:
		 hdulist = fits.open(file_name)
		 try:
			 ralist.append(hdulist[0].header['CRVAL1'])
		 except:
			 KeyError
		 try:
			 declist.append(hdulist[0].header['CRVAL2'])
		 except:
			 KeyError
   if os.path.exists('footcoord.txt'):
       os.remove('footcoord.txt')
   p = open('footcoord.txt', 'w')
   for i in range(0, len(ralist)):
      p.write(str(ralist[i]))
      p.write(" ")
      p.write(str(declist[i]))
      p.write("\n")
   print "Mask coordinates saved in 'footcoord.txt' file. Now plotting..."
   
def create_plot():
   if os.path.exists('mycandelsplot.png'):
	  os.remove('mycandelsplot.png')
   fitscanvas = raw_input("What is the name of the FITS file to which you want your data plotted?\n ->")
   color = raw_input("Press 1 to have grayscale image.\nPress 2 to have colorscale image.\nPress 3 to have a 3-color image.\n ->")
   if color == "1":
	  fig = aplpy.FITSFigure(fitscanvas)
	  data = np.loadtxt('footcoord.txt')
	  ra, dec = data[:, 0], data[:, 1]
	  fig.show_rectangles(ra, dec, 0.05, 0.1, edgecolor='blue', facecolor='none')
	  fig.show_grayscale()
	  fig.tick_labels.set_font(size='small')
	  fig.add_grid()
	  contour = raw_input("Press 1 to add contours to plot.\nPress 2 to proceed to save.\n ->")
	  if contour == "1":
		contourfile = raw_input("What is the name of FITS file with contour data?\n ->")
		fig.show_contour(contourfile, colors='white')
		fig.save('mycandelsplot.png')
		print "Plot is saved as 'mycandelsplot.png' "
	  elif contour == "2":
		fig.save('mycandelsplot.png')
		print "Plot is saved as 'mycandelsplot.png' "
   elif color == "2":
	  fig = aplpy.FITSFigure(fitscanvas)
	  data = np.loadtxt('footcoord.txt')
	  ra, dec = data[:, 0], data[:, 1]
	  fig.show_rectangles(ra, dec, 0.05, 0.1, edgecolor='blue', facecolor='none')
	  fig.show_colorscale()
	  fig.tick_labels.set_font(size='small')
	  fig.add_grid()
	  contour = raw_input("Press 1 to add contours to plot.\nPress 2 to proceed to save.\n ->")
	  if contour == "1":
		contourfile = raw_input("What is the name of FITS file with contour data?\n ->")
		fig.show_contour(contourfile, colors='white')
		fig.save('mycandelsplot.png')
		print "Plot is saved as 'mycandelsplot.png' "
	  elif contour == "2":
		fig.save('mycandelsplot.png')
		print "Plot is saved as 'mycandelsplot.png' "
   elif color == "3":
	  pic = raw_input("What is the name of PNG file containing the 3-color image?\n ->")
	  fig = aplpy.FITSFigure(fitscanvas)
	  data = np.loadtxt('footcoord.txt')
	  ra, dec = data[:, 0], data[:, 1]
	  fig.show_rectangles(ra, dec, 0.05, 0.1, edgecolor='blue', facecolor='none')
	  fig.show_rgb(pic)
	  fig.tick_labels.set_font(size='small')
	  fig.add_grid()
	  contour = raw_input("Press 1 to add contours to plot.\nPress 2 to proceed to save.\n ->")
	  if contour == "1":
		contourfile = raw_input("What is the name of FITS file with contour data?\n ->")
		fig.show_contour(contourfile, colors='white')
		fig.save('mycandelsplot.png')
		print "Plot is saved as 'mycandelsplot.png' "
	  elif contour == "2":
		fig.save('mycandelsplot.png')
		print "Plot is saved as 'mycandelsplot.png' "
   sys.exit()


if __name__== "__main__":
    url = raw_input('URL of ASCII table off Keck Archive\n->')
    open('data.txt', 'wb').write(request.urlopen(url).read())
    info = ascii.read("data.txt")
    #edit the following line for switching instruments
    pilist = info['progpi']
    uniquepilist = list(set(pilist))
    #edit the following line for switching instruments
    projlist = info['progtitl']
    uniqueprojlist = list(set(projlist))
    #edit the following line for switching instruments
    masklist = info['targname']
    #edit the following line for switching instruments
    extimelist = info['elaptime']
    #edit the following line for switching instruments
    filenamelist = info['koaid']
    answer = main_option()  
    while answer != "0":        
        if answer == "1":
            pi_info()
        elif answer == "2":
            project_info()
        elif answer == "3":
            mask_info()  
        elif answer == "4":
            summary_table()
        elif answer == "5":
            path = os.getcwd()
            path = path + "/"
            filenames = [f for f in listdir(path) if isfile(join(path, f))]
            allfiles = sorted(glob.glob(path+'*fits.gz')) + sorted(glob.glob(path+'*fits'))
            inflist = []
            RA = []
            DEC = []
            h = 1
            for file in allfiles:
                f = open(file, 'r')
                y = f.name
                x = y[-25:]
                posofx = h
                fitsname = str(x)
                lenlist = len(filenames)
                print ("Reading " + str(posofx) + " of " + str(lenlist) + " files: " + str(x))
                try:
                    hdr, dat, bs = IO.readmosfits(file,Options.flat)
                    #edit the following line for switching instruments
                    for i, obj in enumerate(bs.targs.field("Target_Name")):
                    	#edit the following line for switching instruments
                        hrs = float(bs.targs.field("RA_Hours")[i])
                        #edit the following line for switching instruments
                        mins = float(bs.targs.field("RA_Minutes")[i])
                        #edit the following line for switching instruments
                        secs = float(bs.targs.field("RA_Seconds")[i])
                        time = (((((secs/60) + mins)/60) + hrs)*15)
                        time = str(round(time, 8))
                        #edit the following line for switching instruments
                        newhrs = float(bs.targs.field("Dec_Degrees")[i])
                        #edit the following line for switching instruments
                        newmins = float(bs.targs.field("Dec_Minutes")[i])
                        #edit the following line for switching instruments
                        newsecs = float(bs.targs.field("Dec_Seconds")[i])
                        temphrs = -newhrs
                        if newhrs > 0:
                            othertime = ((((newsecs/60) + newmins)/60) + newhrs)
                        elif newhrs < 0:
                            othertime = ((((newsecs/60) + newmins)/60) + temphrs)
                        othertime = str(round(othertime, 8))
                        RA.extend(time)
                        DEC.extend(othertime)
                        #edit the following line for switching instruments
                        mask = hdr['maskname']
                        head, sep, tail = mask.partition(' ')
                        #edit the following line for switching instruments
                        if head == "MIRA" or head == "align":
                            continue
                        else:
                            x = [head, fitsname, obj, time, othertime]
                            inflist.extend([x])
                except:
                    Exception
                h = h + 1
            ##saved here just in case anyone in the future figures out how to make catalog matching work
            ##stuck at matching portion at "MARKED LINE" below
            ##catfiles = sorted(glob.glob(path+'*cat'))
            ##catlist = list(catfiles)
            ##if len(catlist) > 0:
            ##    catans = raw_input("Catalog files have been found in the local directory. Should the catalog be matched with the current data set?" + "\n" + "Press 1 for yes. Press 2 for no: ")
            ##    if catans == "1": 
            ##        bans = raw_input("Choose a catalog file: " + str(catlist) + ": ")
            ##        cat=ascii.read(bans,data_start=1)     
            ##        with open(bans,'r') as infile, open('CatData.csv', 'wb') as outfile:
            ##            in_txt = csv.reader(infile, delimiter = ' ')
            ##            out_csv = csv.writer(open('CatData.csv', 'wb'))
            ##            out_csv.writerows(in_txt)
            ##        with open('CatData.csv') as f:
            ##            catdata = [[],[],[],[]] 
            ##            reader = csv.reader(f)
            ##            for row in reader:
            ##                for col in range(2):
            ##                    catdata[col].append(row[col])
            ##        if os.path.exists('CatData.csv'):
            ##            os.remove('CatData.csv')      
            ##        RAtemp = catdata[1]
            ##        RAtot = [float(i) for i in RAtemp]  
            ##        DECtemp = catdata[2]
            ##        DECtot = [float(i) for i in DECtemp]
            ##MARKED LINE
            ##        fullcat=SkyCoord(RAtot*u.degree,DECtot*u.degree,frame='icrs')
            ##        obscat=SkyCoord(RA*u.degree,DEC*u.degree,frame='icrs')
            ##        index, dist2d, dist3d=obscat.match_to_catalog_sky(fullcat)
            ##        print cat['ID'][(index)]
                
            def grouper(sequence):
                group, members = [], set()
                for item in inflist:
                    if group and members.isdisjoint(item):
                        yield group
                        group, members = [], set()
                    group.append(item)
                    members.update(item)
                yield group

            rightlist = []
            for group in grouper(inflist):
                rightlist.extend([group])
            for i in xrange(0, len(rightlist)):
                exec("price%d = %s" % (i + 1, repr(rightlist[i])));
            for e in rightlist:
                title = e[0][0]
                #edit the following line for switching instruments
                objtab = Table(rows=e, names=('Mask', 'File', 'Object', 'Right Acension', 'Declination'))
                print "\n"
                print objtab
                print "\n"
            newans = raw_input("Do you want to save these tables? Press 1 for yes. Press 2 for no: ")
            if newans == "1":
                for k in rightlist:
                    title = k[0][0]
                    with open(title + '_Mask_Objects.txt', 'w') as file:
                        file.writelines('\t'.join(i) + '\n' for i in k)
                    with open(title + '_Mask_Objects.txt','r') as infile, open(title + '_Mask_Objects.csv', 'wb') as outfile:
                        in_txt = csv.reader(infile, delimiter = '\t')
                        out_csv = csv.writer(open(title +  '_Mask_Objects.csv', 'wb'))
                        out_csv.writerows(in_txt)

                    if os.path.exists(title + '_Mask_Objects.txt'):
                        os.remove(title + '_Mask_Objects.txt') 
                print ("\n" + "Tables Saved" + "\n")
                theans = "9"
                while theans != "2":
                    theans = raw_input("Do you want to search for a specific item? Press 1 for Yes. Press 2 for No: ")
                    if theans == "1":
                        print "\n"
                        searchobj = raw_input("Search: ")
                        print "\n"
                        matches = []
                        for i in xrange(0, len(rightlist)):
                            match = [x for x in rightlist[i] if searchobj in x]
                            matches.extend(match)
                        if len(matches)>0:
                        	#edit the following line for switching instruments
                            matchtable = Table(rows=matches, names=('Mask', 'File', 'Object', 'Right Acension', 'Declination'))
                            print "\n"
                            print matchtable
                            print "\n"
                            newans = raw_input("Do you want to save these tables? Press 1 for Yes. Press 2 for No: ")
                            if newans == "1":
                                with open(searchobj + '_Query_Results_Table.txt', 'w') as file:
                                    file.writelines('\t'.join(i) + '\n' for i in matches)
                                with open(searchobj + '_Query_Results_Table.txt','r') as infile, open(searchobj + '_Query_Results_Table.csv', 'wb') as outfile:
                                    in_txt = csv.reader(infile, delimiter = '\t')
                                    out_csv = csv.writer(open(searchobj + '_Query_Results_Table.csv', 'wb'))
                                    out_csv.writerows(in_txt)
                                if os.path.exists(searchobj + '_Query_Results_Table.txt'):
                                    os.remove(searchobj + '_Query_Results_Table.txt')
                                print "\n"
                                print ("Files Saved")
                            elif newans == "2":
                                sys.exit()
                        elif len(matches) == 0:
                            print "\n"
                            print("No matches found")
                    elif theans == "2":
                        sys.exit()
            else:
                theans = "9"
                while theans != "2":
                    theans = raw_input("Do you want to search for a specific item? Press 1 for Yes. Press 2 for No: ")
                    if theans == "1":
                        print "\n"
                        searchobj = raw_input("Search: ")
                        print "\n"
                        matches = []
                        for i in xrange(0, len(rightlist)):
                            match = [x for x in rightlist[i] if searchobj in x]
                            matches.extend(match)
                        if len(matches)>0:
                        	#edit the following line for switching instruments
                            matchtable = Table(rows=matches, names=('Mask', 'File', 'Object', 'Right Acension', 'Declination'))
                            print "\n"
                            print matchtable
                            print "\n"
                            newans = raw_input("Do you want to save this table? Press 1 for Yes. Press 2 for No: ")
                            if newans == "1":
                                with open(searchobj + '_Query_Results_Table.txt', 'w') as file:
                                    file.writelines('\t'.join(i) + '\n' for i in matches)
                                with open(searchobj + '_Query_Results_Table.txt','r') as infile, open(searchobj + '_Query_Results_Table.csv', 'wb') as outfile:
                                    in_txt = csv.reader(infile, delimiter = '\t')
                                    out_csv = csv.writer(open(searchobj + '_Query_Results_Table.csv', 'wb'))
                                    out_csv.writerows(in_txt)
                                print "\n"
                                print ("Files Saved") 
                            elif newans == "2":
                                sys.exit()     
                        elif len(matches) == 0:
                            print "\n"
                            print("No matches found")
                    elif theans == "2":
                        sys.exit()
        elif answer == "6":
            grab_dimensions()
            create_plot()
        elif answer == "0":
            sys.exit()
