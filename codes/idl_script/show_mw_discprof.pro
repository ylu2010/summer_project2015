plotps=0
plotpng=0

tag_name='Model-I'
tag_name='Model-EJ'
tag_name='Model-PR'
;tag_name='Model-II'
;tag_name='Model-III'
dir='/Users/luyu/project/disc/results/'+tag_name+'/m12_smooth_new2/'

fpng='disc_'+tag_name+'.png'

;file=dir+'disc_m3_f.2_t1_p2.dat'
file=dir+'disc_z0.dat'
readcol, file, comment='#', i, r, sigcold, sigstar, sigmole, sigatom, rout, mrhalo, mrstar, mrcold
sigstar_z0=sigstar
sigcold_z0=sigcold
sigatom_z0=sigatom
rhalf_star=interpol(rout, mrstar, 0.5*mrstar[n_elements(mrstar)-1])
rhalf_cold=interpol(rout, mrcold, 0.5*mrcold[n_elements(mrcold)-1])
print,rhalf_star, rhalf_cold, 0.015*300.
;;;;;;;;;;;;;;

file=dir+'disc_z1.dat'
readcol, file, comment='#', i, r, sigcold, sigstar, sigmole, sigatom, rout, mrhalo, mrstar, mrcold
sigstar_z1=sigstar
sigcold_z1=sigcold
sigatom_z1=sigatom
rhalf_star=interpol(rout, mrstar, 0.5*mrstar[n_elements(mrstar)-1])
rhalf_cold=interpol(rout, mrcold, 0.5*mrcold[n_elements(mrcold)-1])
print,rhalf_star, rhalf_cold, 0.015*300.

;;;;;;;;;;;;;;

file=dir+'disc_z2.dat'
readcol, file, comment='#', i, r, sigcold, sigstar, sigmole, sigatom, rout, mrhalo, mrstar, mrcold
sigstar_z2=sigstar
sigcold_z2=sigcold
sigatom_z2=sigatom
rhalf_star=interpol(rout, mrstar, 0.5*mrstar[n_elements(mrstar)-1])
rhalf_cold=interpol(rout, mrcold, 0.5*mrcold[n_elements(mrcold)-1])
print,rhalf_star, rhalf_cold, 0.015*300.



;;;;;;;;;;;;;;
;;;;;;;;;;;;;;
@set_plot_env
window,2, xsize=1100,ysize=400,retain=2     ;plot on Xwindow
!x.margin=[9.0,1.0]
!p.charsize=2.5

multiplot, /initialize
multiplot, [3,1]

xrange=[0.0, 27]
yrange=[1e-0, 1e4]
plot, [0], [0], xtitle='R (kpc)', ytitle=TexToIDL('\Sigma')+' (M'+sunsymbol()+'/pc!u2!n)', xrange=xrange, yrange=yrange, /ylog
oplot, r, sigstar_z0, color=fsc_color('red')
oplot, r, sigcold_z0, color=fsc_color('blue'), linestyle=5
oplot, r, sigatom_z0, color=fsc_color('green'), linestyle=2
xyouts, 3.0+xrange[0], 0.3*yrange[1], 'z=0', charsize=2.5
xyouts, 14+xrange[0], 0.3*yrange[1], tag_name, charsize=2.5
;legend, ['star', 'cold gas', 'atomic gas'], $
;color=[fsc_color('red'), fsc_color('blue'), fsc_color('green')], $
;linestyle=[0,5,2],/right, box=0
multiplot

xrange=[0.01,27]
plot, [0], [0], xtitle='R (kpc)', xrange=xrange, yrange=yrange, /ylog
oplot, r, sigstar_z1, color=fsc_color('red')
oplot, r, sigcold_z1, color=fsc_color('blue'), linestyle=5
oplot, r, sigatom_z1, color=fsc_color('green'), linestyle=2
xyouts, 3.0+xrange[0], 0.3*yrange[1], 'z=1', charsize=2.5
multiplot

plot, [0], [0], xtitle='R (kpc)', xrange=xrange, yrange=yrange, /ylog
oplot, r, sigstar_z2, color=fsc_color('red')
oplot, r, sigcold_z2, color=fsc_color('blue'), linestyle=5
oplot, r, sigatom_z2, color=fsc_color('green'), linestyle=2
xyouts, 3.0+xrange[0], 0.3*yrange[1], 'z=2', charsize=2.5
;xyouts, 14+xrange[0], 0.3*yrange[1], tag_name, charsize=2.5
legend, ['star', 'cold gas', 'atomic gas'], $
color=[fsc_color('red'), fsc_color('blue'), fsc_color('green')], $
linestyle=[0,5,2],/right, box=0
multiplot,/reset
multiplot, /default

@reset_plot_env
end
