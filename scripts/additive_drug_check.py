import numpy as np
rng=np.random.default_rng(20261018); n=200000
B=np.clip(rng.normal(.4,.1,n),0,1); S=np.clip(rng.normal(.4,.1,n),0,1); H=np.clip(rng.normal(.65,.1,n),0.3,1)
pB=np.exp(rng.normal(np.log(3.78),np.log(4.43/3.78)/1.645,n)); pS=np.exp(rng.normal(np.log(3.0),.2,n)); pH=np.exp(rng.normal(np.log(2.0),.25,n))
q=lambda a:np.round(np.percentile(a,[50,5,95]),2)
# TZD: insulin-stimulated disposal +54% (Inzucchi 1998) -> S*1.54
tzd=pS*(np.minimum(S*1.54,1)-S)/1.59
# metformin: EGP -19%, disposal +13% (Inzucchi 1998)
met=(pH*0.19*H + pS*(np.minimum(S*1.13,1)-S))/1.59
print("TZD",q(tzd),"trial ~1.0-1.6 (pioglitazone vs placebo)")
print("MET",q(met),"trial 1.12 (0.92-1.32)")
print("ratio TZD",1.37/np.median(tzd),"MET",1.12/np.median(met))
# needed weights
print("phiS needed for 1.37:",1.37*1.59/(0.4*.54),"phiH needed for 1.12 given phiS 3:",(1.12*1.59-3*0.4*.13)/(0.19*.65))
print("baseline excess glucose mmol/L:",np.median(pB*(0.9-B)+pS*(0.9-S)+pH*(H-.3)))
