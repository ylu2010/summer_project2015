# print_selected_parameters.py

a2 = isel_rdisk
a1 = isel_mzcold
ip = numpy.nonzero(numpy.in1d(a2, a1))[0]

p_ip = p[:,ip]

fout = rundir+ '/p_mzcold_rdisk.csv'
header = ', '.join(param_names[:ndim])

nline = np.shape(p_ip)[1]
for i in range(nline):
    print p_ip[:, i]
    
np.savetxt(fout, np.transpose(p_ip), fmt='%g', delimiter=',', header=header)