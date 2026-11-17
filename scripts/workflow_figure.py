import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
fig,ax=plt.subplots(figsize=(10,5.6),dpi=300); ax.set_xlim(0,100); ax.set_ylim(0,56); ax.axis('off')
def box(x,y,w,h,t,c,fs=8.5,b=False):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.4,rounding_size=1.2",fc=c,ec='#444',lw=0.8))
    ax.text(x+w/2,y+h/2,t,ha='center',va='center',fontsize=fs,weight='bold' if b else 'normal',wrap=True)
def arr(x1,y1,x2,y2): ax.annotate('',(x2,y2),(x1,y1),arrowprops=dict(arrowstyle='->',lw=1.1,color='#333'))
# column headers
for x,t in [(2,'Defect X\nbeta-cell insulin release\n(islets)'),(35,'Defect Y\ninsulin sensitivity\n(skeletal muscle)'),(68,'Defect Z\nhepatic glucose output\n(liver)')]:
    box(x,47,30,7,t,'#dbe8f5',8.5,True)
box(2,36,30,7,'4 GEO cohorts (n = 172)\nStouffer meta-analysis, FDR','#f2f2f2')
box(35,36,30,7,'2 GEO cohorts (n = 142)\nStouffer meta-analysis, FDR','#f2f2f2')
box(68,36,30,7,'3 GEO cohorts (n = 56)\nStouffer meta-analysis, FDR','#f2f2f2')
for x in (17,50,83): arr(x,46.4,x,43.6)
box(2,24,30,7,'Cause test (RQ2): islet eQTL MR\nInsPIRE (420 donors) → DIAMANTE T2D','#fde9d9')
arr(17,35.4,17,31.6)
box(35,24,63,7,'RQ1 result: 1,256 islet genes (SLC2A2, IAPP, ABCC8 ↓); muscle 2 genes;\nliver none.  RQ2 result: 332 MR genes; expression and MR agree for 3 of 11','#fff8dc')
arr(32.6,27.5,34.4,27.5)
box(2,9.5,96,9,'Monte Carlo disease model: one state per defect (B, S, H) → glucose → HbA1c (ADAG)\nCalibration: sulfonylurea (−1.40% vs −1.51%).  External check: thiazolidinedione and\nmetformin under-predicted 3–4-fold.  Prior-sensitivity analysis','#e2f0d9')
arr(17,23.4,17,18.6); arr(66,23.4,66,18.6)
box(2,0.5,96,6,'RQ3: beta-cell restoration −1.47% HbA1c; Y and Z effects likely underestimated (additive model)','#d9ead3',9,True)
arr(50,9.4,50,7.1)
plt.savefig('figures/fig1_workflow.png',bbox_inches='tight',facecolor='white')
