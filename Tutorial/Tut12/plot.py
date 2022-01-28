import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


fig,ax=plt.subplots(1,2,figsize=(6,3))

FCS=nx.read_gpickle('fcs.pickle')
for i in range(len(FCS)):
    if (FCS[i,1:]<0).all():
        break
print(i)
FCS=FCS[:i,:]


FCCS=nx.read_gpickle('fccs.pickle')
for i in range(len(FCCS)):
    if (FCCS[i,1:]<0).all():
        break
FCCS=FCCS[:i,:]

ax[0].semilogx(FCS[:,0],FCS[:,1],'k',label='670 nm')
ax[0].semilogx(FCS[:,0],FCS[:,2],'r',label='518 nm')
ax[0].set_xlabel(r'$\tau$ (ns)')
ax[0].set_ylabel(r'$G_{avg}^{MS}(n_O,\lambda,\tau)$ (ns)')
ax[0].set_xlim(0.2,1000)
ax[0].set_ylim(0,max(FCS[1,:])*1.2)
ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)

ax[1].semilogx(FCCS[:,0],FCCS[:,1],'b',label='670,518 nm')
ax[1].set_xlabel(r'$\tau$ (ns)')
ax[1].set_ylabel(r'$G_{avg}^{MS}(n_O,\lambda,\tau)$ (ns)')
ax[1].set_xlim(0.2,1000)
ax[1].set_ylim(0,max(FCCS[1,:])*1.2)
ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)
fig.tight_layout(rect=[0,0,1,0.9])
fig.legend(ncol=3,edgecolor='k',framealpha=1.0)
plt.show()

