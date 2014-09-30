d=mrdfits('camb_28841594_scalcls.fits',1)

ncmb=n_elements(d.temperature)
xcmb=findgen(ncmb)+1

; this is roughly equal to the number of CMB modes but can also be arbitrary
n_modes=ncmb  ; for now use the same number of modes, can be arbitrary.
yamp=sqrt(double(d[0:n_modes-1].temperature))  ; assume file contains T^2 in K
yamp[0:9]=yamp[0:9]*0.0 ; decrease lowest modes since these cause non-gaussian structure

samp=10.0 ; sample spectra at x10 the number of modes

n=n_modes*samp
xax=findgen(n)/float(samp)

for jj=0,100 do begin

xfld=replicate(0.0,n)

for i=1L,n_modes-1 do begin
   phase=randomu(seed)*2.0*!PI
   ytmp=yamp[i]*sin(xax*!PI*i/(float(ncmb))+phase)
   xfld=xfld+ytmp
;   plot,xax,ytmp
;   tmp=''
;   read,tmp
endfor

ft=fft(xfld)

xaxnew=xax[0:n/2]*2

xax_new=xcmb*2
yplt=((abs(ft)^2)*xax_new*(xax_new+1)*1E12)/!PI

; note that the first 1000 correspond to the first 2000 multipoles

plot,xcmb*2,yplt,psym=3,xrange=[0,2000],yrange=[0,8000]
oplot,xcmb,yamp^2*xcmb*(xcmb+1)*1E12/(2.0*!PI)
if jj eq 0 then ytotpow = yplt else  ytotpow=ytotpow+yplt
oplot,xcmb*2,(ytotpow/(float(jj)+1)),psym=1

tmp=''
read,tmp

endfor

 
end
