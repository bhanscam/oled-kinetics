import numpy as np

# possible initial rates: {1e-6,1e-3,1,1e3,1e6,1e9}
m = 1e-9
mags = []
while m < 1e9:
	m *= 1000
	mags.append(m)

print("possible initial rates: ",mags)

# TODO: add a for loop for every rate in model
paramscan = []
for i1 in mags:
	for i2 in mags:
		for i3 in mags:
			for i4 in mags:
				for i5 in mags:
					paramscan.append([i1,i2,i3,i4,i5])
print("number of sets of rates to try: ",len(paramscan))
np.savetxt("scan_params.txt",paramscan)
