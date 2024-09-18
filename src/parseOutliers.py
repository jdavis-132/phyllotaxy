import re
import pandas as pd
import numpy as np

# Save path to files as a string
path='/Users/jensinadavis/Library/CloudStorage/GoogleDrive-davisjensina@gmail.com/My Drive/Schnable-Lab/in-silico/Data'

# Open file with JS_IDs matched to PI#s
JS_PI_FILE = open(path + '/JSname_PIname.csv', 'r')

# Iterate through JS_PI_FILE and save JS_ID as key and PI# as value in a dictionary
JSDict = {}

for line in JS_PI_FILE:
    # Strip whitespace characters at end of line
    line = line.strip()
    # Split into list by commas
    line = line.split(',')
    # Save to dictionary
    JSDict[line[1]] = line[2]

# Save list of files to process
fileList = ['outliers.csv']

outData = pd.DataFrame(np.nan, index=[0], columns=['plant_name', 'JS_ID', 'plant_num', 'PI_num',
'img_date', 'reason'])

# Iterate through list of files
for i in range(0, len(fileList), 1):
    # Save file name as a string
    fileName = fileList[i]
    
    # Open file for reading
    currFile = open(path + '/' +fileName, 'r')
    
    # Iterate through lines of file
    for line in currFile:
        # Strip whitespace characters at end of line
        line = line.strip()
        # Split line by commas
        line = line.split(',')
        
        # If the line doesn't contain 'Schnable' in the plant name field this is a header row, so ignore it
        if not re.search('Schnable', line[0]):
            continue
        # Initialize empty temp list and fill with info for one row of outData
        tempList = ['plant_name', 'JS_ID', 'plant_num','PI_num','img_date', 'reason']
        
        # Get plant name ID to allow back matching to original data file if needed
        plantName = line[0]
        
        # Split plantName to get the JS ID
        currJS_IDList = plantName.split('-')
        currJS_ID = ''
        plant_num = 1
            
        i = 0
        while not (re.search('JS', currJS_IDList[i])) and not (re.search('js', currJS_IDList[i])):
            i = i + 1
            
        currJS_ID = currJS_IDList[i]
        plant_num = currJS_IDList[i - 1]
        
        # Make all letters uppercase
        currJS_ID = currJS_ID.upper()
        
        # Use JS_ID to get the matching PI# from dictionary
        currPINum = JSDict[currJS_ID]
        
        # Split plantName to get imgDate
        currImgDate = plantName.split('_')
        currImgDate = currImgDate[3]
        
        # Add metadata for this skeleton to the tempList
        tempList[0] = plantName
        tempList[1] = currJS_ID
        tempList[2] = plant_num
        tempList[3] = currPINum
        tempList[4] = currImgDate
        tempList[5] = line[1]
        
        
        # Append tempList to outData
        outData.loc[len(outData)] = tempList

# Write outData to a csv file
outData.to_csv(path + '/parsedOutliers.csv', index=False)
        