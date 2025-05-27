from scipy.interpolate import griddata
import numpy as np

def contourPlotRZ(fbmr, fig, r_grid, z_grid, f_in, title=''):

    fbmfile  = fbmr.fileName.split('/')[-1]
    runid    = fbmfile[:8]
    title += ', Run %s, t =%6.3f s' %(runid, fbmr.time)
    fig.clf()
    ax = fig.add_subplot(111, aspect='equal')
    fig.text(0.5, 0.95, title, ha='center')
    if 'gc_d' in globals():
        for gc in gc_d.values():
            ax.plot(m2cm*gc.r, m2cm*gc.z, 'b-')
    for irho in range(fbmr.r_surf.shape[0]):
        ax.plot(fbmr.r_surf[irho, :], fbmr.z_surf[irho, :], 'r-', linewidth=0.5)
    for jbar, myr in enumerate(fbmr.rbar):
        ax.plot(myr, fbmr.zbar[jbar], 'r-')
    f_grid = griddata((fbmr.r2d, fbmr.z2d), f_in, (r_grid, z_grid), method='cubic')
    plot2d = ax.contourf(r_grid, z_grid, f_grid, levels=50, cmap='viridis')
    fig.colorbar(plot2d, aspect=20)
    ax.set_xlabel('R [cm]')
    ax.set_ylabel('Z [cm]')
# Draw limiter, set plot boundaries
    xext =  (fbmr.rlim_pts.min() + fbmr.rlim_pts.max())*0.025
    yext = (-fbmr.ylim_pts.min() + fbmr.ylim_pts.max())*0.025
    xpol_lim=np.array ([fbmr.rlim_pts.min() - xext, fbmr.rlim_pts.max() + xext])
    ypol_lim=np.array ([fbmr.ylim_pts.min() - yext, fbmr.ylim_pts.max() + yext])
    xpol_lim=np.round(xpol_lim).astype(int)
    ypol_lim=np.round(ypol_lim).astype(int)
    rlin=fbmr.rlim_pts
    zlin=fbmr.ylim_pts
    ax.set_xlim(xpol_lim)
    ax.set_ylim(ypol_lim)
    ax.plot(rlin, zlin, 'g-', linewidth=2.5)
    canvas = fig.canvas
    canvas.draw()

    return ax
