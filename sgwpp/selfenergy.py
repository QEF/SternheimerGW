import sys
import re

def plot_gw(plot_file, plot_data): 
  sigma = plot_data
#top/bottom/efstate:=highest/lowest/valence band top
  top     = 14
  efstate = 10
  bottom  = 2
  plot  = open(plot_file,'w')

  #print >> plot, "unset key"

  print >> plot, "plot '%s' u 1:%d w lp pt 7 lt -1" %(sigma, bottom)

  for i in range(bottom+1, efstate):
      print >> plot, "replot '%s' u 1:%d w lp pt 7 lt -1" %(sigma, i)

  print >> plot, "replot '%s' u 1:%d w lp pt 7 lt 3" %(sigma, efstate)

  for i in range(efstate+1, top+1):
      print >> plot, "replot '%s' u 1:%d w lp pt 7 lt 1" %(sigma, i)

  return

f = open(sys.argv[1]).read()

res   = open('resigma.dat', 'w')
ims   = open('imsigma.dat', 'w')
aspec = open('aspec.dat', 'w')

sigma_regex = re.compile(r'GW qp renorm.*?\n\n(.*?)\Z', re.M | re.S)
sigmare_regex   = re.compile(r'REsigma\n(.*?)IMsigma', re.M | re.S)
sigmaim_regex   = re.compile(r'IMsigma\n(.*?)ASpec', re.M | re.S)
sigmaspec_regex = re.compile(r'ASpec\n(.*?)\n\s{0,}\n', re.M | re.S)

#blocks_text = sigma_regex.findall(f)

block = sigmare_regex.findall(f)
#for line in block[0].split('\n'):
#print block[0]
for line in block[0].split('\n'):
  print >>res, line 

block = sigmaim_regex.findall(f)
for line in block[0].split('\n'):
  print >>ims, line 

block = sigmaspec_regex.findall(f)
for line in block[0].split('\n'):
  print >>aspec, line 

datafiles = ['resigma.dat', 'imsigma.dat', 'aspec.dat'] 
plotfiles = ['plotre.gnu', 'plotim.gnu', 'plotaspec.gnu'] 


for a, b in zip(plotfiles,datafiles):
    plot_gw(a,b)


