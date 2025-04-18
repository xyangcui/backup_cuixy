

ddir = "/home/sunming/data5/cuixy/Subpre_NPJ/data"
DIR = ddir+"/eof_"+var+".nc"
fout = addfile(DIR,"r")
eof  = fout->x
 
printVarSummary(eof) 
ntim = 4840
prinfo = True

pcvar = eof@pcvar
print(sum(pcvar(0:51)))

sig_pcv = eofunc_north(eof@pcvar, ntim, prinfo)