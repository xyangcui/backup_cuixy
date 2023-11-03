;------------------------------------------------------------
load "/home/sunming/data5/cuixy/script/ncl_func/Filter.ncl"

INDIR  = getenv("DATADIR")
OUTDIR = getenv("OUTDIR")
var    = getenv("var")

latmin = tofloat(getenv("latmin"))
latmax = tofloat(getenv("latmax"))
lonmin = tofloat(getenv("lonmin"))
lonmax = tofloat(getenv("lonmax"))

f = addfile(INDIR,"r")

u = f->$var$(:,0,{latmin:latmax},{lonmin:lonmax})

uanom = iso_filter(u,201,100,20,1)

system("rm -f "+OUTDIR)

f1    = addfile(OUTDIR,"c")
f1->u = uanom