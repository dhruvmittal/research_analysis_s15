import os
import averageData

VARIABLES = ['Bare_coupling', 'fugacity']
DATA_FILENAME = 'logfile0.out'

def listSubDir(path):
    result = []
    for dirname, dirnames, filenames in os.walk(path):
        for subdirname in dirnames:
            result.append(os.path.join(dirname, subdirname))
    return result

# Returns tuple: (last_parameter_in_path_value, average_data_for_parameter_value)
def getDataFromSubDirs(list_paths,error=False):
    all_subdir_data = []
    for path in list_paths:
        last_parameter_value = float(path.split('/')[-1].split('_')[-1])
        if error:
            data = averageData.getError(os.path.join(path, DATA_FILENAME))
        else:
            data = averageData.get(os.path.join(path, DATA_FILENAME))
        all_subdir_data.append((last_parameter_value, data))
    sorted_all_subdir_data = sorted(all_subdir_data, key=lambda data:data[0])
    return sorted_all_subdir_data
    

#print getDataFromSubDirs(listSubDir('../exec_coupling_fugacity/Bare_coupling_0.01/'))
