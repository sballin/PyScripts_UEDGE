import datetime
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib
from scipy.optimize import curve_fit
from scipy.special import erfc
plt.rcParams['font.size'] = 12
# Set figsize to the default in multiple ways for consistency across environments
# e.g. figsize is sensitive to whether the grid preview figure is shown or not
# figsize = (6.4, 4.8)
# plt.rcParams['figure.figsize'] = figsize


def plotvar(var, title="UEDGE data", log=False):
    patches = []

    for iy in np.arange(0,com.ny+2):
        for ix in np.arange(0,com.nx+2):
            rcol=com.rm[ix,iy,[1,2,4,3]]
            zcol=com.zm[ix,iy,[1,2,4,3]]
            rcol.shape=(4,1)
            zcol.shape=(4,1)
            polygon = Polygon(np.column_stack((rcol,zcol)), True)
            patches.append(polygon)

    vals=np.zeros((com.nx+2)*(com.ny+2))

    for iy in np.arange(0,com.ny+2):
        for ix in np.arange(0,com.nx+2):
            k=ix+(com.nx+2)*iy
            vals[k] = var[ix,iy]
    
    norm = matplotlib.colors.LogNorm() if log else None
    p = PatchCollection(patches, norm=norm)
    p.set_array(np.array(vals))

    # fig = plt.figure(figsize=figsize)
    ax = plt.gca()

    ax.add_collection(p)
    ax.set_facecolor('lightgray');
    ax.autoscale_view()
    plt.colorbar(p)
    plt.grid(True, which='both', color='#dddddd'); plt.gca().set_axisbelow(True);
    plt.axis('equal')
    
    plt.title(title)
    plt.xlabel(r'$R$ [m]')
    plt.ylabel(r'$Z$ [m]')
    
    plt.tight_layout()

rlabel = r'$R_{omp}-R_{sep}$ [m]'

with PdfPages('plots.pdf') as pdf:
    # Input table
    page = plt.figure(figsize=(8.5, 3.5))
    page.clf()
    txt = 'Simulation settings\n'
    #
    txt += '\nCore ne, ni: '
    if bbb.isnicore[0] == 0:
        txt += 'set flux to curcore/sy locally in ix'
    elif bbb.isnicore[0] == 1:
        txt += 'fixed uniform core density n = %.1e /m^3' % bbb.ncore[0]
    elif bbb.isnicore[0] == 2:
        txt += 'set flux & ni over range'
    elif bbb.isnicore[0] == 3:
        txt += 'set icur=curcore-recycc*fngy, const ni'
    elif bbb.isnicore[0] == 4:
        txt += 'use impur. source terms (impur only)'
    elif bbb.isnicore[0] == 5:
        txt += 'set d(ni)/dy=-ni/lynicore at midp & ni constant poloidally'
    #
    txt += '\nCore neutral density: '
    if bbb.isnicore[1] == 0:
        txt += 'set flux to curcore/sy locally in ix'
    elif bbb.isnicore[1] == 1:
        txt += 'fixed uniform core density ng = %.1e /m^3' % bbb.ncore[1]
    elif bbb.isnicore[1] == 2:
        txt += 'set flux & ng over range'
    elif bbb.isnicore[1] == 3:
        txt += 'set icur=curcore-recycc*fngy, const ng'
    elif bbb.isnicore[1] == 4:
        txt += 'use impur. source terms (impur only)'
    elif bbb.isnicore[1] == 5:
        txt += 'set d(ng)/dy=-ng/lynicore at midp & ng constant poloidally'
    #
    txt += '\nCore temp/power: '
    if bbb.iflcore == 0:
        txt += 'fixed Te = %.1e eV, Ti = %.1e eV' % (bbb.tcoree, bbb.tcorei)
    elif bbb.iflcore == 1:
        txt += 'fixed Pe = %.1e W, Pi = %.1e W' % (bbb.pcoree, bbb.pcorei)
    #
    txt += '\nCore ion parallel velocity (up): '
    if bbb.isupcore[0] == 0:
        txt += 'up=upcore at core boundary'
    elif bbb.isupcore[0] == 1:
        txt += 'd(up)/dy=0 at core boundary'
    elif bbb.isupcore[0] == 2:
        txt += 'd^2(up)/dy^2 = 0'
    elif bbb.isupcore[0] == 3:
        txt += 'fmiy = 0'
    elif bbb.isupcore[0] == 4:
        txt += 'tor. ang mom flux = lzflux & n*up/R=const'
    elif bbb.isupcore[0] == 5:
        txt += 'ave tor vel = utorave & n*up/R=const'
    #
    txt += '\nWall Te: '
    if bbb.istewc == 0:
        txt += 'set by zero energy flux'
    elif bbb.istewc == 1:
        if np.all(bbb.tewallo):
            txt += 'fixed Te = %.1e eV' % bbb.tewallo[0]
        else: 
            txt += 'fixed to 1D spatially varying profile (bbb.tewallo)'
    elif bbb.istewc == 2:
        txt += 'extrapolated'
    elif bbb.istewc == 3:
        txt += 'set Te scale length to lyte'
    elif bbb.istewc == 4:
        txt += 'set feey = bceew*fniy*te'
    #
    txt += '\nWall Ti: '
    if bbb.istiwc == 0:
        txt += 'set by zero energy flux'
    elif bbb.istiwc == 1:
        if np.all(bbb.tiwallo):
            txt += 'fixed Ti = %.1e eV' % bbb.tiwallo[0]
        else:
            txt += 'fixed to 1D spatially varying profile (bbb.tiwallo)'
    elif bbb.istiwc == 2:
        txt += 'extrapolated'
    #
    txt += '\nWall ne, ni: '
    if bbb.isnwcono[0] == 0:
        txt += 'old case; if ifluxni=0, dn/dy=0; if ifluxni=1, fniy=0'
    elif bbb.isnwcono[0] == 1:
        if np.all(bbb.nwallo):
            txt += 'fixed n = %.1e /m^3' % bbb.nwallo[0]
        else:
            txt += 'fixed to 1D spatially varying profile (bbb.nwallo)'
    elif bbb.isnwcono[0] == 2:
        txt += 'extrapolated'
    elif bbb.isnwcono[0] == 3:
        txt += 'approx grad-length lyni, but limited by nwomin'
    #
    txt += '\nWall ng: '
    if bbb.isnwcono[1] == 0:
        txt += 'old case; if ifluxni=0, dn/dy=0; if ifluxni=1, fniy=0'
    elif bbb.isnwcono[1] == 1:
        if np.all(bbb.nwallo):
            txt += 'fixed ng = %.1e /m^3' % bbb.nwallo[0]
        else:
            txt += 'fixed to 1D spatially varying profile (bbb.nwallo)'
    elif bbb.isnwcono[1] == 2:
        txt += 'extrapolated'
    elif bbb.isnwcono[1] == 3:
        txt += 'approx grad-length lyni, but limited by nwomin'
    #
    txt += '\nPlates H recycling coefficient: %.3e' % bbb.recycp[0]
    #
    txt += '\nImpurity model: '
    if bbb.isimpon == 0:
        txt += 'no impurities'
    elif bbb.isimpon == 2:
        txt += 'fixed-fraction model'
    elif bbb.isimpon == 3:
        txt += 'average-impurity-ion model (disabled)'
    elif bbb.isimpon == 4:
        txt += 'INEL multi-charge-state model (disabled)'
    elif bbb.isimpon == 5:
        txt += "Hirshman's reduced-ion model"
    elif bbb.isimpon == 6:
        txt += 'force-balance model or nusp_imp > 0; see also isofric for full-Z drag term'
    elif bbb.isimpon == 7:
        txt += 'simultaneous fixed-fraction and multi-charge-state (isimpon=6) models'
    #
    txt += '\nImpurity: %s' % (globals['imp'] if 'imp' in globals().keys() else 'none')
    #
    txt += '\nImpurity concentration: '
    if np.all(bbb.afracs):
        txt += 'fixed to uniform %.3e' % bbb.afracs[0,0]
    else:
        txt += 'fixed to 2D spatially varying profile (bbb.afracs)'
    #
    page.text(0.1, 0.9, txt, transform=page.transFigure, size=12, horizontalalignment='left', verticalalignment='top')
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
    
    # Input graphs
    # D OMP
    plt.figure(figsize=(8.5, 8.5))
    plt.subplot(221)
    plt.plot(com.yyc, bbb.dif_use[:,:,0][bbb.ixmp], c='C0', label=r'$D_{omp}$')
    plt.plot(com.yyc, bbb.dif_use[:,:,0][bbb.ixmp-1], c='C1', label=r'$D_{imp}$')
    plt.legend()
    plt.xlabel(rlabel)
    plt.title(r"$D$")
    plt.yscale('log')
    plt.grid(True, which='both', color='#dddddd'); plt.gca().set_axisbelow(True);
    # Chi OMP
    plt.subplot(222)
    plt.plot(com.yyc, bbb.kye_use[bbb.ixmp], c='C0', lw=1, label=r'$\chi_{e,omp}$')
    plt.plot(com.yyc, bbb.kyi_use[bbb.ixmp], c='C0', lw=2, label=r'$\chi_{i,omp}$')
    plt.plot(com.yyc, bbb.kye_use[bbb.ixmp-1], c='C1', lw=1, label=r'$\chi_{e,imp}$')
    plt.plot(com.yyc, bbb.kyi_use[bbb.ixmp-1], c='C1', lw=2, label=r'$\chi_{i,imp}$')
    plt.ylabel(r'$\chi$ [m$^2$/s]')
    plt.xlabel(rlabel)
    plt.yscale('log')
    plt.legend()
    plt.grid(True, which='both', color='#dddddd'); plt.gca().set_axisbelow(True);
    # D 2D image
    plt.subplot(223)
    plotvar(bbb.dif_use[:,:,0],  title=r"$D$ [m$^2$/s]", log=True)
    # Chi 2D image
    plt.subplot(224)
    plotvar(bbb.kye_use,  title=r"$\chi$ [m$^2$/s]", log=True)
    # Finish
    plt.tight_layout()
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
    
    plt.figure(figsize=(8.5, 8.5))
    plt.subplot(221)
    plotvar(bbb.te/bbb.ev,  title=r"$T_e$ [eV]", log=True)
    plt.subplot(222)
    plotvar(bbb.ti/bbb.ev,  title=r"$T_i$ [eV]", log=True)
    plt.subplot(223)
    plotvar(bbb.ni[:,:,0],  title=r"$n_{i,e}$ [m$^{-3}$]", log=True)
    plt.subplot(224)
    plotvar(bbb.ni[:,:,1],  title=r"$n_g$ [m$^{-3}$]", log=True)
    plt.tight_layout()
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
    
    plt.figure(figsize=(8.5, 15))
    #
    lineNi = {'c': 'C0',     'ls': '-',  'lw': 2, 'label': r'$n_{i,e}$'}
    lineNg = {'c': 'C0',     'ls': '--', 'lw': 2, 'label': r'$n_g$'}
    lineTi = {'c': 'tomato', 'ls': '-',  'lw': 2, 'label': r'$T_i$'}
    lineTe = {'c': 'tomato', 'ls': '-',  'lw': 1, 'label': r'$T_e$'}
    #
    plt.subplot(421)
    plt.title('Inner midplane')
    plt.plot(com.yyc, bbb.ni[:,:,0][bbb.ixmp-1], **lineNi)
    plt.plot(com.yyc, bbb.ni[:,:,1][bbb.ixmp-1], **lineNg)
    plt.ylabel(r'$n$ [m$^{-3}$]')
    plt.xlabel(rlabel)
    plt.grid(True, which='both', color='#dddddd'); plt.gca().set_axisbelow(True);
    plt.legend()
    #
    plt.subplot(422)
    plt.title('Outer midplane')
    plt.plot(com.yyc, bbb.ni[:,:,0][bbb.ixmp], **lineNi)
    plt.plot(com.yyc, bbb.ni[:,:,1][bbb.ixmp], **lineNg)
    plt.ylabel(r'$n$ [m$^{-3}$]')
    plt.xlabel(rlabel)
    plt.grid(True, which='both', color='#dddddd'); plt.gca().set_axisbelow(True);
    plt.legend()
    #
    plt.subplot(423)
    plt.title('Inner midplane')
    plt.plot(com.yyc, bbb.te[bbb.ixmp-1]/bbb.ev, **lineTe)
    plt.plot(com.yyc, bbb.ti[bbb.ixmp-1]/bbb.ev, **lineTi)
    plt.ylabel(r'$T$ [eV]')
    plt.xlabel(rlabel)
    plt.grid(True, which='both', color='#dddddd'); plt.gca().set_axisbelow(True);
    plt.legend()
    #
    plt.subplot(424)
    plt.title('Outer midplane')
    plt.plot(com.yyc, bbb.te[bbb.ixmp]/bbb.ev, **lineTe)
    plt.plot(com.yyc, bbb.ti[bbb.ixmp]/bbb.ev, **lineTi)
    plt.ylabel(r'$T$ [eV]')
    plt.xlabel(rlabel)
    plt.grid(True, which='both', color='#dddddd'); plt.gca().set_axisbelow(True);
    plt.legend()
    #
    plt.subplot(425)
    plt.title('Inner plate')
    plt.plot(com.yyc, bbb.ni[:,:,0][1], **lineNi)
    plt.plot(com.yyc, bbb.ni[:,:,1][1], **lineNg)
    plt.ylabel(r'$n$ [m$^{-3}$]')
    plt.xlabel(rlabel)
    plt.grid(True, which='both', color='#dddddd'); plt.gca().set_axisbelow(True);
    plt.legend()
    #
    plt.subplot(426)
    plt.title('Outer plate')
    plt.plot(com.yyc, bbb.ni[:,:,0][com.nx], **lineNi)
    plt.plot(com.yyc, bbb.ni[:,:,1][com.nx], **lineNg)
    plt.ylabel(r'$n$ [m$^{-3}$]')
    plt.xlabel(rlabel)
    plt.grid(True, which='both', color='#dddddd'); plt.gca().set_axisbelow(True);
    plt.legend()
    #
    plt.subplot(427)
    plt.title('Inner plate')
    plt.plot(com.yyc, bbb.te[1]/bbb.ev, **lineTe)
    plt.plot(com.yyc, bbb.ti[1]/bbb.ev, **lineTi)
    plt.ylabel(r'$T$ [eV]')
    plt.xlabel(rlabel)
    plt.grid(True, which='both', color='#dddddd'); plt.gca().set_axisbelow(True);
    plt.legend()
    # plt.yscale('log')
    #
    plt.subplot(428)
    plt.title('Outer plate')
    plt.plot(com.yyc, bbb.te[com.nx]/bbb.ev, **lineTe)
    plt.plot(com.yyc, bbb.ti[com.nx]/bbb.ev, **lineTi)
    plt.ylabel(r'$T$ [eV]')
    plt.xlabel(rlabel)
    plt.grid(True, which='both', color='#dddddd'); plt.gca().set_axisbelow(True);
    plt.legend()
    # Finish
    plt.tight_layout()
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()
    
    # qpar
    bbb.fetx = bbb.feex+bbb.feix
    bpol_local = 0.5*(com.bpol[:,:,2] + com.bpol[:,:,4])
    bphi_local = 0.5*(com.bphi[:,:,2] + com.bphi[:,:,4])
    btot_local = np.sqrt(bpol_local**2+bphi_local**2)
    rrf = bpol_local/btot_local
    iy = com.iysptrx+1
    xq = com.yyc[iy:]
    ixo = com.ixpt2[0]
    ixi = com.ixpt1[0]
    #-radial profile of qpar below entrance to the outer leg
    qparo = bbb.fetx[ixo,iy:]/com.sx[ixo,iy:]/rrf[ixo,iy:]
    intqo = np.sum(bbb.fetx[ixo+1,iy:]) # integral along first set of edges that *enclose* the outer divertor
    #-radial profile of qpar below entrance to the inner leg
    qpari = -bbb.fetx[ixi,iy:]/com.sx[ixi,iy:]/rrf[ixi,iy:]
    intqi = np.sum(-bbb.fetx[ixi-1,iy:])
    # lamda_q fits
    expfun = lambda x, A, lamda_q_inv: A*np.exp(-x*lamda_q_inv) # needs to be in this form for curve_fit to work
    omax = np.argmax(qparo) # only fit stuff to right of max
    qofit, _ = curve_fit(expfun, xq[omax:], qparo[omax:])
    lqo = 1000./qofit[1] # lamda_q in mm
    imax = np.argmax(qpari) # only fit stuff to right of max
    qifit, _ = curve_fit(expfun, xq[imax:], qpari[imax:])
    lqi  = 1000./qifit[1] # lamda_q in mm
    #
    plt.figure(figsize=(8.5, 8.5))
    plt.subplot(211)
    plt.title(r'$q_\parallel$ at divertor entrance ($P_{inner}:P_{outer}$ = 1:%.1f)' % (intqo/intqi))
    plt.plot(xq*1000, qparo, c='C0', label=r'Outer divertor ($P_{xpt}$ = %.1e W)' % intqo)
    plt.plot(xq*1000, qpari, c='C1', label=r'Inner divertor ($P_{xpt}$ = %.1e W)' % intqi)
    # plt.yscale('log')
    ylim = plt.gca().get_ylim()
    plt.plot(xq[omax:]*1000, expfun(xq, *qofit)[omax:], c='C0', ls=':', 
             label='Outer exp. fit ($\lambda_q$ = %.3f mm)' % lqo)
    plt.plot(xq[imax:]*1000, expfun(xq, *qifit)[imax:], c='C1', ls=':', 
             label='Inner exp. fit ($\lambda_q$ = %.3f mm)' % lqi)
    plt.ylim(ylim)
    max_lamda_q = max(1./qofit[1], 1./qifit[1])
    plt.xlim([-max_lamda_q*1000, 10*max_lamda_q*1000])
    plt.xlabel( r'$R_{omp}-R_{sep}$ [mm]')
    plt.ylabel(r'$q_\parallel$ [W/m$^2$]')
    plt.grid(True, which='both', color='#dddddd'); plt.gca().set_axisbelow(True);
    plt.legend()
    
    # qsurf
    bbb.fetx = bbb.feex+bbb.feix
    #-radial profile of qpar below entrance to the outer leg
    qsurfo = bbb.fetx[com.nx,:]/com.sx[com.nx,:]
    intqo = np.sum(bbb.fetx[com.nx,:])
    #-radial profile of qpar below entrance to the inner leg
    qsurfi = -bbb.fetx[0,:]/com.sx[0,:]
    intqi = np.sum(-bbb.fetx[0,:])
    # lamda_q fits
    def qEich(rho, q0, S, lqi, qbg, rho_0):
        rho = rho - rho_0
        # lqi is inverse lamda_q
        return q0/2*np.exp((S*lqi/2)**2-rho*lqi)*erfc(S*lqi/2-rho/S)+qbg
    oguess = (np.max(qsurfo)-np.min(qsurfo), 1./qofit[1], qofit[1], np.min(qsurfo), 0)
    qsofit, _ = curve_fit(qEich, com.yyc, qsurfo, p0=oguess)
    lqeo, So = 1000./qsofit[2], qsofit[1]*1000 # lamda_q and S in mm
    iguess = (np.max(qsurfi)-np.min(qsurfi), 1./qifit[1], qifit[1], np.min(qsurfi), 0)
    qsifit, _ = curve_fit(qEich, com.yyc, qsurfi, p0=iguess)
    lqei, Si = 1000./qsifit[2], qsifit[1]*1000 # lamda_q and S in mm 
    #
    plt.subplot(212)
    plt.title(r'$q_{surf}$ ($P_{inner}:P_{outer}$ = 1:%.1f)' % (intqo/intqi))
    plt.plot(com.yyc*1000, qsurfo, c='C0', label=r'Outer divertor ($P_{surf}$ = %.1e W)' % intqo)
    plt.plot(com.yyc*1000, qEich(com.yyc, *qsofit), c='C0', ls=':',
             label=r'Outer Eich fit ($\lambda_q$ = %.3f mm, $S$ = %.3f mm)' % (lqeo, So))
    plt.plot(com.yyc*1000, qsurfi, c='C1', label=r'Inner divertor ($P_{surf}$ = %.1e W)' % intqi)
    plt.plot(com.yyc*1000, qEich(com.yyc, *qsifit), c='C1', ls=':',
             label=r'Inner Eich fit ($\lambda_q$ = %.3f mm, $S$ = %.3f mm)' % (lqei, Si))
    plt.xlabel(r'$R_{omp}-R_{sep}$ [mm]')
    plt.ylabel(r'$q_{surf}$ [W/m$^2$]')
    plt.grid(True, which='both', color='#dddddd'); plt.gca().set_axisbelow(True);
    plt.yscale('log')
    plt.legend()
    plt.tight_layout()
    pdf.savefig()  # saves the current figure into a pdf page
    plt.close()

    d = pdf.infodict()
    d['Title'] = 'UEDGE plots'
    d['CreationDate'] = datetime.datetime.today()
