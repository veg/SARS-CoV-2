import sys
import json
import argparse
import operator

arguments = argparse.ArgumentParser(description='Summarize selection analysis results.')

arguments.add_argument('-f', '--file',   help = 'File to process', required = True, type = str)
import_settings = arguments.parse_args()
hpd = 0.95

with open (import_settings.file, "r") as fh:
    fubar = json.load (fh)
    grid = fubar ['grid']
    for site, p in fubar['posterior']["0"].items():
        posteriors = [(i,v) for i,v in enumerate (p[0])]
        hdpc = sorted(posteriors, key=operator.itemgetter(1), reverse = True)
        psum = 0.
        i = 0
                
        alpha = sum ([k[0] * posteriors[i][1] for i, k in enumerate (grid)])
        beta  = sum ([k[1] * posteriors[i][1] for i, k in enumerate (grid)])
        ppp = sum ([posteriors[i][1] for i, k in enumerate (grid) if k[0]<k[1]])
        ppn = sum ([posteriors[i][1] for i, k in enumerate (grid) if k[0]>k[1]])
        non0 = sum ([posteriors[i][1] for i, k in enumerate (grid) if k[0] > 0])
        hpd0 = hpd * non0
        mean_omega  = sum ([k[1]/k[0] * posteriors[i][1] for i, k in enumerate (grid) if k[0]>0]) / non0
        
        omega = []
        
        while (psum <= hpd0):
            idx = hdpc[i][0]
            if grid[idx][0] > 0:
                psum += hdpc[i][1]
                omega.append (grid[idx][1]/grid[idx][0])
            i+=1
            
        #print (alpha, beta, mean_omega, min(omega), max (omega), ppn, ppp)
        if ppn > 0.95 or ppp > 0.95:
            print (mean_omega, min(omega), max (omega), ppn, ppp)
        
        #sys.exit (0)
        
    