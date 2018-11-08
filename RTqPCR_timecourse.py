#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 20:01:00 2017

@author: Rory
"""

"""
Program for automated analysis and graphing of 96-well qPCR timecourse data.


Takes in qPCR data from 96-well plate reader. Assumes time course times are
0,15,30,45,60 and 75 minutes. Takes in sample names and data file, creates
new data file with modified and normalised data points as per usual de Bruin
lab protocol. Saves dat file as excel file, saves a pdf of graphs representing
data, and can email data and graphs to a specified email address.

"""


### DEFINITIONS:

import os
from os import listdir
import sys
import pandas as pd
import numpy as np
from os.path import isfile, join
import datetime
import xlsxwriter
import matplotlib.pyplot as plt

def create_excel_conversion():
    alpha = "abcdefghijklmnopqrstuvwxyz"
    h = 97
    w = 26
    exmat = [[0 for x in range(w)] for y in range(h)]
    for i in range(h):
        for j in range(w):
            exmat[i][j] = alpha[j].upper() + str(i+1)
    return exmat

def ec(datpos):
    h = 97
    w = 26
    uppos = datpos.upper()
    for i in range(h):
        for j in range(w):
            if exmat[i][j] == uppos:
                return(((i-1),(j-1)))

def choose_datafiles():
    os.chdir("/Users/Rory/Desktop")
    mypath = os.getcwd()
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    excels = []
    for file in onlyfiles:
        for c in range(len(file)-3):
            test = file[c]+file[c+1]+file[c+2]+file[c+3]
            if test == ".xls":
                excels.append(file)
    print("Excel sheets found in working directory:\n")
    for i in range(len(excels)):
        if i % 2 == 0:
            print(str(i+1) + "." + excels[i] +"   ", end='')
        else:
            print(str(i+1) + "." + excels[i])        
    print(" ")
    filechoice = "0"
    options = list(range(1,(len(excels)+1)))
    while int(filechoice) not in options:
        filechoice = input("Choose data sheet: ")
    filename = excels[int(filechoice)-1]
    print("Chosen: " + filename +".")
    excels.remove(filename)
    print("Remaining excel sheets: ")
    for i in range(len(excels)):
        if i % 2 == 0:
            print(str(i+1) + "." + excels[i] +"   ", end='')
        else:
            print(str(i+1) + "." + excels[i])
    loadingchoice = "0"
    options = str(list(range(1,(len(excels)+1))))
    while loadingchoice not in options:
        loadingchoice = input("Choose data for loading control: ")
    loadingname = excels[int(loadingchoice)-1]
    print("Chosen: " + loadingname +".")
    return filename, loadingname
    
def loadfile(filename):
    file = pd.ExcelFile(filename)
    dat = file.parse('0')
    print("successfully loaded " + filename)
    for i in range(12):
        dat['column' + str(i+16)] = pd.Series(np.array([None]*96))
    return dat

def detect_gene_name(filename):
    skip = ["!","£","$","%","&", "(",")","'",
            '"', "/",".",",","?","[","]","\n"]
    for c in skip:
        filename = filename.replace(c,"")
    numbers = "0,1,2,3,4,5,6,7,8,9"
    num = numbers.split(",")
    name = filename.lower()
    text = name.split(" ")
    for word in text:
        if len(word) >= 4:
            if (word[0:3].isalpha()) and (word[3] in num):
                return(word)
    return(None)
    
def enter_gene_name(filename):
    gene = detect_gene_name(filename).upper()
    if gene != None:
        q = ""
        while q != "y" and q != "n":
            q = input("Gene name detected: "+gene+". Is this correct? (y/n) ")
        if q == "y":
            return gene
        if q == "n":
            gene = None
    if gene == None:
        gene = input("Enter name of target gene here: ").upper()
    return gene

def enter_cell_strains():
    a = ""
    while a != "y" and a != "n":
        a = input("""Are these the strains you are using: 
            WT, Swi4Δ, Mbp1Δ, SVS1pr-CDC21, CDC21pr-SVS1? (y/n) """)
    if a == "y":
        strainlist = ["WT", "Swi4Δ", "Mbp1Δ", "SVS1pr-CDC21", "CDC21pr-SVS1"]
    if a == "n":
        strains = input("Enter names of strains, separated by spaces (max 5):\n")
        strainlist = strains.split(" ")
        if len(strainlist) > 5:
            sys.exit("Error. Please enter no more than five strains.")    
    return(strainlist)

def create_output_file(gene, dataframe):
    now = datetime.datetime.now()
    timestamp = now.strftime("%Y%m%d %H%M%S")
    title = gene + " qPCR Timecourse Data (" + timestamp +").xlsx"
    outfile = pd.ExcelWriter(title, engine = 'xlsxwriter')
    dataframe.to_excel(outfile)
    outfile.save()
    return title
    
### END OF DEFINITIONS.
    
### START OF PROGRAM:
exmat = create_excel_conversion() 
# creates a matrix for excel/A1 and Python/[0][0] conversions.
# use ec() to convert A1 to [0][0]
filename, controlname = choose_datafiles()
# allows you to select which of the excel sheets in pwd 
# you want to use for data and for loading control
dat = loadfile(filename)
lcd = loadfile(controlname)
# loads the data into a Pandas DataFrame.
gene = enter_gene_name(filename)
# gathers, through title parsing or user input, the name of the test gene.
strainlist = enter_cell_strains()
strno = len(strainlist)
if strno == 5:
    cycles = 4
    fifth = 1
elif strno < 5:
    cycles = strno
    fifth = 0
 # gathers name of 5 different cell strains

# Formatting:
setindex = np.arange(len(dat.index))
dat = dat.set_index(setindex)
lcd = lcd.set_index(setindex)
headers = ["Mean","SD","ACT1","Normalised to ACT1", "Normalised to WT0",
               "Normalised to WT30", "Normalised to 0 min"]
times = ["0","15","30","45","60","75"]
ind = dat.index
col = dat.columns
for i in range(len(strainlist)):
    index = ind[i*8]
    column = col[17]
    dat.loc[index,column] = strainlist[i]
    for j in range(len(headers)):
        index = ind[(i*8)+1]
        column = col[17+j]
        dat.loc[index,column] = headers[j]
    for k in range(len(times)):
        index = ind[(i*8)+2+k]
        column = col[16]
        dat.loc[index,column] = times[k]

# Calculations:

# first, create the loading control array:
lclist = []
lcsds = []
datacol = col[6]
editcount = 0
lcedits = []
for i in range(cycles):
    for j in range(len(times)):
        datastore = np.array([0.0,0.0,0.0])
        newdatastore = []
        diffstore = np.array([0.0,0.0,0.0])
        dataind = ind[(12*j)+(3*i)]
        for k in range(3):
            datastore[k] = lcd.loc[(dataind+k),datacol]
        datasd = np.std(datastore, ddof=1)
        if datasd > 0.3:
            bc = abs(datastore[2] - datastore[1])
            diffstore[0] = bc
            ac = abs(datastore[2] - datastore[0])
            diffstore[1] = ac
            ab = abs(datastore[1] - datastore[0])
            diffstore[2] = ab
            m = np.where(diffstore == np.min(diffstore))
            n = int(m[0])
            for val in range(len(datastore)):
                if datastore[val] != datastore[n]:
                    newdatastore.append(datastore[val])
            newdatastore = np.array(newdatastore)
            mean = np.mean(newdatastore)
            sd = np.std(newdatastore, ddof=1)
            editcount = editcount + 1
            editindex = j + 6*i
            lcedits.append(editindex)
        else:
            mean = np.mean(datastore)
            sd = np.std(datastore, ddof=1)
        lclist.append(mean)
        lcsds.append(sd)
if fifth == 1:       
    indices = [72,84,75,87,78,90]
    for t in range(len(times)):
        datastore = np.array([0.0,0.0,0.0])
        newdatastore = []
        diffstore = np.array([0.0,0.0,0.0])
        dataind = indices[t]
        for k in range(3):
                datastore[k] = lcd.loc[(dataind+k),datacol]
        datasd = np.std(datastore, ddof=1)
        if datasd > 0.3:
            bc = abs(datastore[2] - datastore[1])
            diffstore[0] = bc
            ac = abs(datastore[2] - datastore[0])
            diffstore[1] = ac
            ab = abs(datastore[1] - datastore[0])
            diffstore[2] = ab
            m = np.where(diffstore == np.min(diffstore))
            n = int(m[0])
            for val in range(len(datastore)):
                if datastore[val] != datastore[n]:
                    newdatastore.append(datastore[val])
            newdatastore = np.array(newdatastore)
            mean = np.mean(newdatastore)
            sd = np.std(newdatastore, ddof=1)
            editcount = editcount + 1
            editindex = t + 4*6
            lcedits.append(editindex)
        else:
            mean = np.mean(datastore)
            sd = np.std(datastore, ddof=1)
        lclist.append(mean)
        lcsds.append(sd)

print("NB: datapoints removed from loading control data: " + str(editcount))
    
# then load the LC data into dat dataframe:
newlclist = np.array(lclist)
for i in range(len(lclist)):
    column = col[19]
    index = 2+i+((i-i%6)/6)*2
    dat.loc[index,column] = np.around(newlclist[i], decimals=2)
    
for i in range(len(lcedits)):
    column = col[24]
    index = 2+lcedits[i]+((lcedits[i]-lcedits[i]%6)/6)*2
    dat.loc[index,column] = "NB: Loading control edited"

# then create array for means and sds from actual data:
datlist = []
datsds = []
datacol = col[6]
editcount = 0
datedits = []
for i in range(cycles):
    for j in range(len(times)):
        datastore = np.array([0.0,0.0,0.0])
        newdatastore = []
        diffstore = np.array([0.0,0.0,0.0])
        dataind = ind[(12*j)+(3*i)]
        for k in range(3):
            datastore[k] = dat.loc[(dataind+k),datacol]
        datasd = np.std(datastore, ddof=1)
        if datasd > 0.3:
            bc = abs(datastore[2] - datastore[1])
            diffstore[0] = bc
            ac = abs(datastore[2] - datastore[0])
            diffstore[1] = ac
            ab = abs(datastore[1] - datastore[0])
            diffstore[2] = ab
            m = np.where(diffstore == np.min(diffstore))
            n = int(m[0])
            for val in range(len(datastore)):
                if datastore[val] != datastore[n]:
                    newdatastore.append(datastore[val])
            newdatastore = np.array(newdatastore)
            mean = np.mean(newdatastore)
            sd = np.std(newdatastore, ddof=1)
            editcount = editcount + 1
            editindex = j + 6*i
            datedits.append(editindex)
        else:
            mean = np.mean(datastore)
            sd = np.std(datastore, ddof=1)
        datlist.append(mean)
        datsds.append(sd)
if fifth == 1:       
    indices = [72,84,75,87,78,90]
    for t in range(len(times)):
        datastore = np.array([0.0,0.0,0.0])
        newdatastore = []
        diffstore = np.array([0.0,0.0,0.0])
        dataind = indices[t]
        for k in range(3):
                datastore[k] = dat.loc[(dataind+k),datacol]
        datasd = np.std(datastore, ddof=1)
        if datasd > 0.3:
            bc = abs(datastore[2] - datastore[1])
            diffstore[0] = bc
            ac = abs(datastore[2] - datastore[0])
            diffstore[1] = ac
            ab = abs(datastore[1] - datastore[0])
            diffstore[2] = ab
            m = np.where(diffstore == np.min(diffstore))
            n = int(m[0])
            for val in range(len(datastore)):
                if datastore[val] != datastore[n]:
                    newdatastore.append(datastore[val])
            newdatastore = np.array(newdatastore)
            mean = np.mean(newdatastore)
            sd = np.std(newdatastore, ddof=1)
            editcount = editcount + 1
            editindex = t + 4*6
            datedits.append(editindex)
        else:
            mean = np.mean(datastore)
            sd = np.std(datastore, ddof=1)
        datlist.append(mean)
        datsds.append(sd)

print("NB: datapoints removed from target gene data: " + str(editcount))

# then load means and sds into dataframe:
newdatlist = np.array(datlist)
newdatsds = np.array(datsds)
for i in range(len(datlist)):
    column = col[17]
    index = 2+i+((i-i%6)/6)*2
    dat.loc[index,column] = np.around(newdatlist[i], decimals=2)
for i in range(len(datsds)):
    column = col[18]
    index = 2+i+((i-i%6)/6)*2
    dat.loc[index,column] = np.around(newdatsds[i], decimals=2) 

for i in range(len(datedits)):
    column = col[25]
    index = 2+datedits[i]+((datedits[i]-datedits[i]%6)/6)*2
    dat.loc[index,column] = "NB: Data edited"

# then perform requisite calculations and import into dataframe:

# first, normalise to loading control:
ratlist = []
for i in range(len(datlist)):
    column = col[20]
    index = 2+i+((i-i%6)/6)*2
    tar = datlist[i]
    ref = lclist[i]
    rat = (2**-tar)/(2**-ref)
    dat.loc[index,column] = rat
    ratlist.append(rat)
    
# next, normalise to WT0:
wtzlist = []
for i in range(len(datlist)):
    column = col[21]
    index = 2+i+((i-i%6)/6)*2
    num = ratlist[i]
    den = ratlist[0]
    wtz = num/den
    dat.loc[index,column] = wtz
    wtzlist.append(wtz)

# then, normalise to WT30:
wttlist = []
for i in range(len(datlist)):
    column = col[22]
    index = 2+i+((i-i%6)/6)*2
    top = ratlist[i]
    bot = ratlist[2]
    wtt = top/bot
    dat.loc[index,column] = wtt
    wttlist.append(wtt)   

# finally, normalise to t0:
tznlist = []
for i in range(len(datlist)):
    column = col[23]
    index = 2+i+((i-i%6)/6)*2
    up = ratlist[i]
    down = ratlist[i - (i%6)]
    tzn = up/down
    dat.loc[index,column] = tzn
    tznlist.append(tzn)   

# then create and export the chosen time-course graph:
now = datetime.datetime.now()
timestamp = now.strftime("%Y%m%d %H%M%S")
g = ""
while g != "y" and g != "n":
    g = input("Create Graph? (y/n) ")    
if g == "y":
    options = ['1','2','3']
    c = ""
    print("Graph options:")
    print("1. Expression normalised to WT0")
    print("2. Expression normalised to WT30")
    print("3. Expression normalised to t0")
    while c not in options:
        c = input("Please choose: ")
    if c == '1':
        x = [0,15,30,45,60,75]
        fir = wtzlist[:6]
        if len(strainlist) >= 2:
            sec = wtzlist[6:12]
        if len(strainlist) >= 3:
            thi = wtzlist[12:18]
        if len(strainlist) >= 4:
            fou = wtzlist[18:24]
        if len(strainlist) == 5:
            fif = wtzlist[24:]
        f = plt.figure()
        plt.plot(x,fir,color='blue')
        if len(strainlist) >= 2:
            plt.plot(x,sec,color='red')
        if len(strainlist) >= 3:
            plt.plot(x,thi,color='green')
        if len(strainlist) >= 4:
            plt.plot(x,fou,color='yellow')
        if len(strainlist) == 5:
            plt.plot(x,fif,color='orange')
        plt.legend(strainlist, loc='upper left')
        plt.ylabel(str(gene)+" expression relative to WT0")
        plt.xlabel("time (minutes)")
        plt.title(str(gene)+" Expression Levels")
        plt.xlim(0, 75)
        plt.xticks(x)
        plt.grid(axis='y', linestyle='dotted')
        plt.grid(axis='x', linestyle='dotted')
        plt.show()
        f.savefig(str(gene)+" Expression Levels (" +timestamp+").pdf",
                  bbox_inches='tight')
    if c == '2':
        x = [0,15,30,45,60,75]
        fir = wttlist[:6]
        if len(strainlist) >= 2:
            sec = wttlist[6:12]
        if len(strainlist) >= 3:
            thi = wttlist[12:18]
        if len(strainlist) >= 4:
            fou = wttlist[18:24]
        if len(strainlist) == 5:
            fif = wttlist[24:]
        f = plt.figure()
        plt.plot(x,fir,color='blue')
        if len(strainlist) >= 2:
            plt.plot(x,sec,color='red')
        if len(strainlist) >= 3:
            plt.plot(x,thi,color='green')
        if len(strainlist) >= 4:
            plt.plot(x,fou,color='yellow')
        if len(strainlist) == 5:
            plt.plot(x,fif,color='orange')
        plt.legend(strainlist, loc='upper left')
        plt.ylabel(str(gene)+" expression relative to WT30")
        plt.xlabel("time (minutes)")
        plt.title(str(gene)+" Expression Levels")
        plt.xlim(0, 75)
        plt.xticks(x)
        plt.grid(axis='y', linestyle='dotted')
        plt.grid(axis='x', linestyle='dotted')
        plt.show()
        f.savefig(str(gene)+" Expression Levels (" +timestamp+").pdf",
                  bbox_inches='tight')
    if c == '3':
        x = [0,15,30,45,60,75]
        fir = tznlist[:6]
        if len(strainlist) >= 2:
            sec = tznlist[6:12]
        if len(strainlist) >= 3:
            thi = tznlist[12:18]
        if len(strainlist) >= 4:
            fou = tznlist[18:24]
        if len(strainlist) == 5:
            fif = tznlist[24:]
        f = plt.figure()
        plt.plot(x,fir,color='blue')
        if len(strainlist) >= 2:
            plt.plot(x,sec,color='red')
        if len(strainlist) >= 3:
            plt.plot(x,thi,color='green')
        if len(strainlist) >= 4:
            plt.plot(x,fou,color='yellow')
        if len(strainlist) == 5:
            plt.plot(x,fif,color='orange')
        plt.legend(strainlist, loc='upper left')
        plt.ylabel(str(gene)+" expression relative to t0")
        plt.xlabel("time (minutes)")
        plt.title(str(gene)+" Expression Levels")
        plt.xlim(0, 75)
        plt.xticks(x)
        plt.grid(axis='y', linestyle='dotted')
        plt.grid(axis='x', linestyle='dotted')
        plt.show()
        f.savefig(str(gene)+" Expression Levels (" +timestamp+").pdf",
                  bbox_inches='tight')


# then export the data sheet 
title = create_output_file(gene, dat)

# then email it to someone: [credit goes to https://gist.github.com/rdempsey]
e = ""
while e != 'y' and e != 'n':
    e = input("Would you like me to email the graph and data to you? (y/n) ")
if e == 'y':
    message = """Hello,
    
    Attached are the qPCR data and graph from your recent analysis 
    using Python.
    
    Rory's Program.
    
    
    
    
    """
    graphtitle = str(gene)+" Expression Levels (" + timestamp+").pdf"
    graphpath = os.path.abspath(graphtitle)
    datapath = os.path.abspath(title)                           
    attachpaths = [graphpath,datapath]
    import os
    import smtplib
    from email import encoders
    from email.mime.base import MIMEBase
    from email.mime.multipart import MIMEMultipart
    from email.mime.text import MIMEText
    
    COMMASPACE = ', '
    
    sender = "roryspythonprograms@gmail.com"
    gmail_password = "computersdontneedpasswords"
    address = str(input("Enter email address: "))
    recipients = address
    
    # Create the enclosing (outer) message
    outer = MIMEMultipart()
    outer['Subject'] = 'Python qPCR Analysis Data and Graph'
    outer['To'] = COMMASPACE.join(recipients)
    outer['From'] = sender
    outer.preamble = 'You will not see this in a MIME-aware mail reader.\n'
    text = message
    outer.attach(MIMEText(text, 'plain'))
    
    # List of attachments
    attachments = attachpaths
    
    # Add the attachments to the message
    for file in attachments:
        try:
            with open(file, 'rb') as fp:
                msg = MIMEBase('application', "octet-stream")
                msg.set_payload(fp.read())
            encoders.encode_base64(msg)
            msg.add_header('Content-Disposition', 'attachment', filename=os.path.basename(file))
            outer.attach(msg)
        except:
            print("Unable to open one of the attachments. Error: ", sys.exc_info()[0])
            raise
    
    composed = outer.as_string()
    
    # Send the email
    try:
        with smtplib.SMTP('smtp.gmail.com', 587) as s:
            s.ehlo()
            s.starttls()
            s.ehlo()
            s.login(sender, gmail_password)
            s.sendmail(sender, recipients, composed)
            s.close()
        print("Email sent!")
    except:
        print("Unable to send the email. Error: ", sys.exc_info()[0])
        raise
if e == 'n':
    print("Okay, goodbye")





























