import sys

sample = ''
sequences = 0
mapped = 0
unmapped = 0
mq0 = 0
for line in sys.stdin :
    s = line.split()

    if line.startswith('# The command line was:') :
        sample = s[-1].split('/')[-1].split('.')[0]

    if line.startswith('SN\traw total sequences:') :
        sequences = int(s[-1])
    
    if line.startswith('SN\treads mapped:') :
        mapped = int(s[-1])

    if line.startswith('SN\treads unmapped:') :
        unmapped = int(s[-1])

    if line.startswith('SN\treads MQ0:') :
        mq0 = int(s[3])

assert sequences == mapped + unmapped
assert mq0 <= mapped

useful = mapped - mq0
rate = float(0)
if sequences :
    rate = float(useful) * 100.0 / sequences

print(sample, rate, sep = ',')
