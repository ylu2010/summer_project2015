plotps=1
plotpng=0

lgmvir=11.0
mass_name='m'+string(lgmvir, format='(f4.1)')
mass_name='m'+string(lgmvir, format='(i02)')
nhalo=10

fpng='metal_hist_'+mass_name+'.png'
fps='metal_hist_'+mass_name+'.eps'

yield_const=0.03/(1.-0.43)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
tag_name='Model-EJ'
dir='/Users/luyu/project/disc/metal/'+tag_name+'/'+mass_name+'_smooth/'

file=dir+'sample_z0.0.dat'
readcol, file, comment='#', mhalo, rhalo, vhalo, mstar, rstar, mcold
ngal=n_elements(mhalo)
mhalo0=mhalo

file_hist=dir+'hist.dat'
readcol, file_hist, comment='#', z, thubble, mhalo, mhot, mcold, mstar, meject, $
        matom, mmole, mion, $
        rhalo, rdisc, rhalfcold, rhalfstar, rcooling, concen, $
        mratehalo, mratecooling, mratestar, mroutflow,$
        vhalo, shalo, tcool, $
		mzhot, mzcold, mzstar, mzeject
ntime=n_elements(z)/ngal
a=reform(mstar,ntime,ngal)
meanmstar=mean(a, dimension=2)
medianmstar=median(a, dimension=2)

a=reform(mzhot,ntime,ngal)
meanmzhot=mean(a, dimension=2)
medianmzhot=median(a, dimension=2)

a=reform(mzcold,ntime,ngal)
meanmzcold=mean(a, dimension=2)
medianmzcold=median(a, dimension=2)

a=reform(mzstar,ntime,ngal)
meanmzstar=mean(a, dimension=2)
medianmzstar=median(a, dimension=2)

a=reform(mzeject,ntime,ngal)
meanmzeject=mean(a, dimension=2)
medianmzeject=median(a, dimension=2)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
mstar_ej = meanmstar
mzhot_ej = meanmzhot
mzcold_ej = meanmzcold
mzstar_ej = meanmzstar
mzeject_ej = meanmzeject
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
tag_name='Model-RI'
dir='/Users/luyu/project/disc/metal/'+tag_name+'/'+mass_name+'_smooth/'

file=dir+'sample_z0.0.dat'
readcol, file, comment='#', mhalo, rhalo, vhalo, mstar, rstar, mcold
ngal=n_elements(mhalo)
mhalo0=mhalo

file_hist=dir+'hist.dat'
readcol, file_hist, comment='#', z, thubble, mhalo, mhot, mcold, mstar, meject, $
        matom, mmole, mion, $
        rhalo, rdisc, rhalfcold, rhalfstar, rcooling, concen, $
        mratehalo, mratecooling, mratestar, mroutflow,$
        vhalo, shalo, tcool, $
        mzhot, mzcold, mzstar, mzeject
ntime=n_elements(z)/ngal
a=reform(mstar,ntime,ngal)
meanmstar=mean(a, dimension=2)
medianmstar=median(a, dimension=2)

a=reform(mzhot,ntime,ngal)
meanmzhot=mean(a, dimension=2)
medianmzhot=median(a, dimension=2)

a=reform(mzcold,ntime,ngal)
meanmzcold=mean(a, dimension=2)
medianmzcold=median(a, dimension=2)

a=reform(mzstar,ntime,ngal)
meanmzstar=mean(a, dimension=2)
medianmzstar=median(a, dimension=2)

a=reform(mzeject,ntime,ngal)
meanmzeject=mean(a, dimension=2)
medianmzeject=median(a, dimension=2)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
mstar_ri = meanmstar
mzhot_ri = meanmzhot
mzcold_ri = meanmzcold
mzstar_ri = meanmzstar
mzeject_ri = meanmzeject
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
tag_name='Model-PR'
dir='/Users/luyu/project/disc/metal/'+tag_name+'/'+mass_name+'_smooth/'

file=dir+'sample_z0.0.dat'
readcol, file, comment='#', mhalo, rhalo, vhalo, mstar, rstar, mcold
ngal=n_elements(mhalo)
mhalo0=mhalo

file_hist=dir+'hist.dat'
readcol, file_hist, comment='#', z, thubble, mhalo, mhot, mcold, mstar, meject, $
        matom, mmole, mion, $
        rhalo, rdisc, rhalfcold, rhalfstar, rcooling, concen, $
        mratehalo, mratecooling, mratestar, mroutflow,$
        vhalo, shalo, tcool, $
        mzhot, mzcold, mzstar, mzeject
ntime=n_elements(z)/ngal
a=reform(mstar,ntime,ngal)
meanmstar=mean(a, dimension=2)
medianmstar=median(a, dimension=2)

a=reform(mzhot,ntime,ngal)
meanmzhot=mean(a, dimension=2)
medianmzhot=median(a, dimension=2)

a=reform(mzcold,ntime,ngal)
meanmzcold=mean(a, dimension=2)
medianmzcold=median(a, dimension=2)

a=reform(mzstar,ntime,ngal)
meanmzstar=mean(a, dimension=2)
medianmzstar=median(a, dimension=2)

a=reform(mzeject,ntime,ngal)
meanmzeject=mean(a, dimension=2)
medianmzeject=median(a, dimension=2)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
mstar_pr = meanmstar
mzhot_pr = meanmzhot
mzcold_pr = meanmzcold
mzstar_pr = meanmzstar
mzeject_pr = meanmzeject
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
thubble=hubble_time_z(z)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; observations
;;;;;;;;;;;;;;
;;
;;;;;;;;;;;;;;
t0=thubble[n_elements(z)-1]
t=thubble
t_lookback=t0-t

color=0
;;;;;;;;;;;;;;
@set_plot_env
if plotps then device,filename=fps,xsize=10,ysize=4.3,/inches,xoffset=1.0,yoffset=1,/color 
;window,2, xsize=900,ysize=350,retain=2     ;plot on Xwindow
!p.multi=[0,3,1]
!x.margin=[7.0,3.0]
!y.margin=[3.5, 2.5]
!p.charsize=2.0
!x.style=8+1

zp=[0.0, 0.5, 1.0, 2.0, 4.0]
tp=hubble_time_z(0.0)-hubble_time_z(zp)
string_tage=string(tp, format='(f5.2)')
string_z=string(zp, format='(f3.1)')
string_z=['0.0', '0.5', '1', '2', '4']

multiplot, /initialize
multiplot, [3,1]

yrange=[1.e-3,50.0]*10.^(lgmvir)*0.17*0.01*yield_const
xrange=[0,13.6]
plot, [0], [0], xtitle='', ytitle='M!iz!n (M!i'+sunsymbol()+'!n)', $
xrange=xrange, yrange=yrange, /ylog ;, /xlog
axis, xaxis=1, xtitle='', xticks=n_elements(zp),xtickv=tp, xtickn=string_z,charsize=1.8
oplot, t_lookback, mstar_ej*yield_const, linestyle=2
oplot, t_lookback, mzstar_ej, color=fsc_color('red')
oplot, t_lookback, mzcold_ej, color=fsc_color('blue'), linestyle=5
oplot, t_lookback, exp(alog(smooth(mzhot_ej,1))), color=fsc_color('orange'), linestyle=4
oplot, t_lookback, exp(alog(smooth(mzeject_ej,1))), color=fsc_color('purple'), linestyle=3

xyouts, 0.15, 0.77, 'Model-EJ', /norm
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
multiplot

np=n_elements(t_lookback)
nf=10.0
ns=np/nf
is=[fix(findgen(ns)*nf), np-2]

plot, [0], [0], xtitle='!6 look back time (Gyr)', $
xrange=xrange, yrange=yrange, /ylog ;, /xlog
axis, xaxis=1, xtitle='z', xticks=n_elements(zp[1:*]),xtickv=tp[1:*], xtickn=string_z[1:*],charsize=1.8
oplot, t_lookback[is], mstar_ri[is]*yield_const, linestyle=2
oplot, t_lookback[is], mzstar_ri[is], color=fsc_color('red')
oplot, t_lookback[is], mzcold_ri[is], color=fsc_color('blue'), linestyle=5
oplot, t_lookback[is], exp(alog(smooth(mzhot_ri[is],1))), color=fsc_color('orange'), linestyle=4
oplot, t_lookback[is], exp(alog(smooth(mzeject_ri[is],1))), color=fsc_color('purple'), linestyle=3

xyouts, 0.44, 0.77, 'Model-RI', /norm
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
multiplot

plot, [0], [0], xtitle='', $
xrange=xrange, yrange=yrange, /ylog ;, /xlog
axis, xaxis=1, xtitle='', xticks=n_elements(zp[1:*]),xtickv=tp[1:*], xtickn=string_z[1:*],charsize=1.8
np=n_elements(t_lookback)
nf=60.0
ns=np/nf
is=[fix(findgen(ns)*nf), np-2]

oplot, t_lookback, mstar_pr*yield_const, linestyle=2
oplot, t_lookback, mzstar_pr, color=fsc_color('red')
oplot, t_lookback, mzcold_pr, color=fsc_color('blue'), linestyle=5
oplot, t_lookback[is], exp(alog(smooth(mzhot_pr[is],1))), color=fsc_color('orange'), linestyle=4
oplot, t_lookback, exp(alog(smooth(mzeject_pr,1))), color=fsc_color('purple'), linestyle=3

xyouts, 0.7, 0.77, 'Model-PR', /norm
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if lgmvir ge 12.0 then legend, ['star', 'cold', 'hot', 'eject'], linestyle=[0,5,4,3], $
color=fsc_color(['red','blue','orange','purple']),box=0,bottom=1, charsize=1.5
xyouts, 0.13, 0.30, 'log M!ivir!n/M!i'+sunsymbol()+'!n='+ string(alog10(median(mhalo0)),format='(f4.1)'), /norm, charsize=1.8

multiplot,/reset
multiplot, /default

@reset_plot_env
end

