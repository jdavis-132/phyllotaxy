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
fileList = ['angles_three_days_4cm.csv', 'angles_three_days_6cm.csv', 'angles_three_days_8cm.csv', 'angles_three_days_12cm.csv']

outData = pd.DataFrame(np.nan, index=[0], columns=['plant_name', 'JS_ID', 'plant_num', 'PI_num',
'img_date', 'voxel_len', 'phi_0', 'phi_1', 'phi_2', 'phi_3', 'phi_4', 'phi_5', 'phi_6',
'phi_7', 'phi_8', 'phi_9', 'phi_10', 'phi_11', 'phi_12', 'phi_13', 'phi_14'])

# Iterate through list of files
for i in range(0, len(fileList), 1):
    # Save file name as a string
    fileName = fileList[i]
    # Save voxel length as a string
    voxelLen = fileName.split('_')
    voxelLen = voxelLen[3].split('.')
    voxelLen = voxelLen[0]
    
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
        # If the array has less than 10 columns it doesn't have useful data (i.e. <=1 phi value)
        elif len(line) < 10:
            continue
        # If the accuracy is less than 0.70 or topology skeleton is false, ignore this line
        elif float(line[1]) < 0.7 or not line[2]:
            continue
        
        # Initialize empty temp list and fill with info for one row of outData
        tempList = ['plant_name', 'JS_ID', 'plant_num','PI_num','img_date', 'voxel_len',
        'phi_0','phi_1', 'phi_2', 'phi_3', 'phi_4', 'phi_5', 'phi_6','phi_7', 'phi_8',
        'phi_9', 'phi_10', 'phi_11', 'phi_12', 'phi_13', 'phi_14']
        
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
        tempList[5] = voxelLen
        
        k = 6
        
        for j in range(6, len(line), 3):
            tempList[k] = line[j]
            k = k + 1
            
            if k == len(tempList):
                break
        
        # Append tempList to outData
        outData.loc[len(outData)] = tempList

# Write outData to a csv file
outData.to_csv(path + '/processedData.csv', index=False)
        