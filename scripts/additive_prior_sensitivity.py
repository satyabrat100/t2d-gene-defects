import numpy as np
from scipy.stats import binom
rng=np.random.default_rng(20261018); n=200000
B=np.clip(rng.normal(.4,.1,n),0,1); S=np.clip(rng.normal(.4,.1,n),0,1); H=np.clip(rng.normal(.65,.1,n),0.3,1)
sdB=np.log(4.43/3.78)/1.645
def run(sB,sS,sH,mS=3.0,mH=2.0):
    pB=np.exp(rng.normal(np.log(3.78),sB,n)); pS=np.exp(rng.normal(np.log(mS),sS,n)); pH=np.exp(rng.normal(np.log(mH),sH,n))
    dX=pB*(1-B)/1.59; dY=pS*(1-S)/1.59; dZ=pH*(H-.3)/1.59
    return dX,dY,dZ,pB,pS,pH
dX,dY,dZ,pB,pS,pH=run(sdB,.2,.25)
print("medians",np.median(dX),np.median(dY),np.median(dZ))
print("P X>Y",(dX>dY).mean(),"P Y>Z",(dY>dZ).mean(),"P X>Y>Z",((dX>dY)&(dY>dZ)).mean())
# rank correlation of (dX-dY) with priors
from scipy.stats import spearmanr
for nm,v in [("phiB",pB),("phiS",pS),("B0",B),("S0",S)]: print(nm,spearmanr(v,dX-dY)[0])
dX,dY,dZ,*_=run(sdB*2,.4,.5); print("wide priors P X>Y",(dX>dY).mean(),"P Y>Z",(dY>dZ).mean())
# equal-prior check: phiS mean = phiB
dX,dY,dZ,*_=run(sdB,.2,.25,mS=3.78); print("phiS=phiB P X>Y",(dX>dY).mean())
# phiH needed for Z>=Y at median: phiH*0.35 = 3.0*0.6
print("phiH needed",3.0*.6/.35)
print("binom P(<=3 of 11)",binom.cdf(3,11,.5),"two-sided",2*binom.cdf(3,11,.5))
