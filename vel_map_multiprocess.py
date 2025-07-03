################################################################################
# IMPORT MODULES
#-------------------------------------------------------------------------------
import datetime
START = datetime.datetime.now()

import os.path, warnings

import numpy as np

from astropy.table import Table
import astropy.units as u

from multiprocessing import Process, Queue, Value
from queue import Empty

from ctypes import c_long


import sys
sys.path.insert(1, ) # CHANGE PATH AS NECESSARY
from mapSmoothness_functions import how_smooth

#warnings.simplefilter('ignore', np.RankWarning)
#warnings.simplefilter('ignore', RuntimeWarning)
################################################################################


def process_1_galaxy(job_queue, i, 
                     return_queue, 
                     num_masked_gal, 
                     num_not_smooth, 
                     num_missing_photo,
                     VEL_MAP_FOLDER, 
                     IMAGE_DIR, 
                     IMAGE_FORMAT, 
                     DRP_index, 
                     map_smoothness_max, 
                     DRP_table, 
                     vel_function, 
                     NSA_index, 
                     NSA_table):
    '''
    Main body of for-loop for processing one galaxy.
    '''

    ############################################################################
    # Open file to which we can write all the print statements.
    #---------------------------------------------------------------------------
    outfile = open('Process_' + str(i) + '_output.txt', 'wt')
    sys.stdout = outfile
    sys.stderr = outfile
    ############################################################################

    while True:
        try: 
            gal_ID = job_queue.get(timeout=1.0)
        except Empty:
        
            print('Queue is empty!', flush=True)
        
            outfile.close()
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__
            
            print('Worker', i, 'redirected stdout and stderr.', flush=True)
            
            return_queue.close()
            
            print('Worker', i, 'closed the return queue.', flush=True)
            
            return_queue.join_thread()
            
            print('Worker', i, 'joined the return queue.', flush=True)
            
            job_queue.close()
            
            print('Worker', i, 'closed the job queue.', flush=True)
            
            job_queue.join_thread()
            
            print('Worker', i, 'returned successfully', datetime.datetime.now(), flush=True)
            return
    
        start = datetime.datetime.now()

        ############################################################################
        # Fitting process here
        #---------------------------------------------------------------------------
        
        try:

            # main function call here

        except:
            print(gal_ID, 'CRASHED! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<', 
                            flush=True)
                        
            raise
        
        fit_time = datetime.datetime.now() - start

                

        print(gal_ID, "velocity map fit", fit_time, flush=True)
        ############################################################################

        output_tuple = (best, fit, params, and, uncertainties, idx)
        return_queue.put(output_tuple)



################################################################################
################################################################################
################################################################################

job_queue = Queue()
return_queue = Queue()

num_tasks = len(FILE_IDS)

# Load jobs into queue
for i,gal_ID in enumerate(FILE_IDS):
        
    job_queue.put(gal_ID)

    #if i > 10:
    #    num_tasks = 12
    #    break


print('Starting processes', datetime.datetime.now(), flush=True)

processes = []

for i in range(12):

    p = Process(target=process_1_galaxy, args=(job_queue, i, 
                                               return_queue, 
                                               num_masked_gal, 
                                               num_not_smooth, 
                                               num_missing_photo,
                                               VEL_MAP_FOLDER, 
                                               IMAGE_DIR, 
                                               IMAGE_FORMAT, 
                                               DRP_index, 
                                               map_smoothness_max, 
                                               DRP_table, 
                                               'tail', 
                                               NSA_index, 
                                               NSA_table))
    
    p.start()

    processes.append(p)

print('Populating output table', datetime.datetime.now(), flush=True)

################################################################################
# Iterate through the populated return queue to fill in the table
#-------------------------------------------------------------------------------
num_processed = 0

print(num_tasks)

while num_processed < num_tasks:

    try:
        return_tuple = return_queue.get(timeout=1.0)
    except:
        continue

    ############################################################################
    # Write the best-fit values and calculated parameters to a text file in 
    # ascii format.
    #---------------------------------------------------------------------------

    if PARAMETER is not None:

        # put PARAMETER in table

    num_processed += 1

    if num_processed % 5 == 0:
        DRP_table.write(TABLE_NAME.FITS, 
                format='fits', overwrite=True)
        print('Table written ', num_processed, flush=True)
    
    print(num_processed, ': ', idx)


print('Finished populating output table', datetime.datetime.now(), flush=True)
# Go through all the processes and join them back to the parent.
for p in processes:
    p.join(None)

DRP_table.write(TABLE_NAME.FITS, 
                format='fits', overwrite=True)
print('Table written', flush=True)

FINISH = datetime.datetime.now()
print("Runtime:", FINISH - START, flush=True)