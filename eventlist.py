import os
import numpy as np

# Function to list files in a specific folder with a given extension
def list_files(folder_path, extension='SAC'):
    file_list = []
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith(f'.{extension}'):
                file_list.append(os.path.join(root, file))
    return file_list

# Function to count the number of traces for each frequency
def count_traces_by_frequency(lht_folder):
    trace_count_by_frequency = {}
    lht_files = list_files(lht_folder)
    for file_path in lht_files:
        # Extract frequency information from the file path
        frequency_info = os.path.basename(file_path).split('.')[3]
        trace_count_by_frequency[frequency_info] = trace_count_by_frequency.get(frequency_info, 0) + 1
    return trace_count_by_frequency

# Specify the base path for events up to "2007-12-31"
base_path = "./2017-2018_processed"
# Specify the target frequency
target_frequency = "0.03"
# write the list to a file
fw = open('filelist' + '.f' + target_frequency,'w')

total_events = 0  # Variable to store the total number of events
flist = []
evnumb = []
# Loop through all event folders and print the traces and counts for each frequency
for root, dirs, files in os.walk(base_path):
    for dir_name in dirs:
        lht_folder = os.path.join(root, dir_name, 'LHZ', f'f-{target_frequency}')
        
        # Check if the folder with the target frequency exists
        if os.path.exists(lht_folder):
            total_events += 1  # Increment the total number of events
            # Get the list of files for the LHT folder
            lht_files = list_files(lht_folder)
            # Get the count of traces for each frequency
            trace_count_by_frequency = count_traces_by_frequency(lht_folder)
            
            # Print the total number of traces and event number in the specified format
            print(f"{len(lht_files)}  {total_events}")
            evnumb.append(len(lht_files))
            for file in lht_files:
                flist.append(file)
                print(file)
# write the total event number to the output file
fw.writelines(str(len(evnumb)) + '\n')

# write the traces by event to the output file
ttrace = 0
for iev in range(total_events):
    fw.writelines(str(evnumb[iev]) + ' ' + str((iev+1)) + '\n')
    
# order the trace file names for each event
# the "41:44" in 'flist[ttrace -1][41:44]' is for the numbers in the file name(e.g. 001, 011)
# change "41:44" accordingly for your files

    trnumb = []
    trindx = []
    eflist = []
    for itrace in range(evnumb[iev]):
        ttrace += 1
        #print(flist[ttrace -1][36:39])
        trnumb.append(int(flist[ttrace-1][35:37]))
        eflist.append(flist[ttrace-1])
    trindx = np.argsort(trnumb)
    for j in range (evnumb[iev]):
        fw.writelines(eflist[trindx[j]] + '\n')
    
fw.close()
