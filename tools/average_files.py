#!/usr/bin/python
import os
import os.path


#**********************************************************************
def isNumber(s):
    try:
        i = float(s)
    except ValueError:
        return False;
    else:
        return True;

#**********************************************************************

dirs = [ '../platelet/16x12_Ht0.20_Ca0.25/run0/tmp',
	 '../platelet/16x12_Ht0.20_Ca0.25/run1/tmp',
	 '../platelet/16x12_Ht0.20_Ca0.25/run2/tmp',
	 '../platelet/16x12_Ht0.20_Ca0.25/run3/tmp',
	 '../platelet/16x12_Ht0.20_Ca0.25/run4/tmp',
	 '../platelet/16x12_Ht0.20_Ca0.25/run5/tmp' ];
out_dir = '../platelet/16x12_Ht0.20_Ca0.25/mean';

#fn_base = 'cell_Ht.000000_to_100000.dat';
#fn_base = 'probe_vel.dat';
#fn_base = 'tracer_vel.dat';
#fn_base = 'tracer_corr.dat';
fn_base = 'tracer_dz2_mean.dat';

nfiles = 0
nlines = 0
is_vals = [ ];
vals = [ ];

for dir in dirs:
    fn = dir + '/' + fn_base;

    if (not (os.path.isfile(fn))):
        continue
    else:
        print fn;
        file = open(fn, 'r');

    lines = file.readlines();

    if (nfiles == 0):
        # First run
	nlines = len(lines);
	is_vals = [None]*nlines;
	vals = [None]*nlines;

	for i in range(nlines):
	    line = lines[i];
	    vals_tmp = line.split();

	    if (len(line) == 0):
		is_vals[i] = False;
	    else:
		is_vals[i] = isNumber(vals_tmp[0]);

	    if (not(is_vals[i])): 
	        vals[i] = line.rstrip();
	    else:
		for j in range(len(vals_tmp)):
		    vals_tmp[j] = float(vals_tmp[j]);
		vals[i] = vals_tmp;
    else:
	for i in range(nlines):
	    if (not(is_vals[i])):
	        continue;

	    line = lines[i];
	    vals_tmp = line.split();

	    for j in range(len(vals_tmp)):
		vals_tmp[j] = float(vals_tmp[j]);
		vals[i][j] += float(vals_tmp[j]);

    file.close();
    nfiles += 1;

# Averaging
for i in range(nlines):
    if (is_vals[i]):
        for j in range(len(vals[i])):
	    vals[i][j] /= max(nfiles,1);

if ( not(os.path.isdir(out_dir)) ):
    os.mkdir(out_dir);
fn = out_dir + '/' + fn_base;
file = open(fn, 'w');

for i in range(nlines):
    line = '';

    if (not(is_vals[i])):
        line = vals[i];
    else:
	for j in range(len(vals[i])):
	    line += ' %12.5E' %( vals[i][j] );
    file.write(line + '\n');
file.close();

print "Number of files =", nfiles;
print "Averaged data is written to", fn;
