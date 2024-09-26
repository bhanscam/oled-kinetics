import numpy as np

pl = []
with open('exciplex_al_pl.txt', 'r') as file:
	for line in file:
		line = line.strip()
		words = line.split()
		pl.append(float(words[1]))

noise = np.random.normal(0,0.1,len(pl))

pl_noise = pl + noise

np.savetxt('exciplex_al_pl_noise.txt',pl_noise)
