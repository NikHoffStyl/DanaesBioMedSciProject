import os
import math

def mkdirIfNotExists(outDirectory_, opts=''):
    if not os.path.exists(outDirectory_):
        os.system('mkdir ' + opts + ' ' + outDirectory_)

def convertRatioToFloat(x):
    a,b = x.split("/")
    c = float(int(a)/int(b))
    return c

def getLoge(x):
    return math.log(x)

def countGenes(x):
    return len(x.split(";"))

def rmGOkey(x):
    return x.split("(GO")[0]
