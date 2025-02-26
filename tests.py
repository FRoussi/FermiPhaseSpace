import subprocess
import json

systFolder = 'systemsFolder'
exName = 'hc2'
filePath = systFolder + '/' + exName + '.json'

print('### FilePath : ' + filePath + ' ###')

# Create the file :
print("\n### makeSiteSystem ###")
subprocess.run(['./makeSitesSystem',exName, 'honeyCombs', '2'])


# Test the HF :
print('\n### fps HF ###')
subprocess.run(['./fps', filePath])

# Modify the json file
with open(filePath, 'r') as jsonFile:
    jsondt = json.load(jsonFile)


jsondt['system']['IC']['type'] = 'tbGroundState'
jsondt['parameters']['dt'] = 0.05
jsondt['parameters']['n_traj'] = 100

jsondt['method']='ED'

with open(filePath, 'w') as jsonFile:
    json.dump(jsondt, jsonFile, indent=4)


# Test ED
print('\n### fps ED ###')
subprocess.run(['./fps', filePath])


# Test fTWA
print('\n### fps fTWA ###')
with open(filePath, 'r') as jsonFile:
    jsondt = json.load(jsonFile)

jsondt['method']='fTWA'

with open(filePath, 'w') as jsonFile:
    json.dump(jsondt, jsonFile, indent=4)

subprocess.run(['./fps', filePath])


# Test GPSR
print('\n### fps GPSR ###')

with open(filePath, 'r') as jsonFile:
    jsondt = json.load(jsonFile)

jsondt['method']='GPSR'

with open(filePath, 'w') as jsonFile:
    json.dump(jsondt, jsonFile, indent=4)

subprocess.run(['./fps', filePath])
