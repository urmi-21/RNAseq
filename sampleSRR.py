import sys
import random
```
```
smap={}
with open(sys.argv[1],'r') as f:
	data=f.read().splitlines()
	

for l in data:
	thisSrp=l.split('\t')[0]
	thisSrr=l.split('\t')[1]
	if thisSrp in smap:
		smap[thisSrp].append(thisSrr)
	else:
		newL=[thisSrr]
		smap[thisSrp]=newL

#print smap
random.seed(9001)
#for each srp select a random srr
for srp in smap:
	print(random.choice(smap[srp]))
