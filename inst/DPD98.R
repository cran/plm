# @---------------------------- one-step estimates -----------------------------@
  
vc=s*inxmx;                                               @ covariance matrix @
  stderr=sqrt(diag(vc));                                      @ standard errors @
  trat=b./stderr;                                                    @ t-ratios @
    pvt=2*cdfnc(abs(trat));                                            @ p-values @
    
    if m0125;
  
  if sc2;                                            @ serial correlation tests @
    v2v=y2y+b'x2x*b-y2x*b-b'x2y;
  v2x=y2x-b'x2x;
  vhat2=(yky2+b'xkx2*b-2*ykx2*b-2*v2x*inxmx*xz*inzkz*(zky2-zkx2*b)+
    v2x*inxmx*v2x')*s;
  m2=v2v./sqrt(vhat2);
  pvm2=2*cdfnc(abs(m2));
endif;

if sc1;
  v1v=y1y+b'x1x*b-y1x*b-b'x1y;
  v1x=y1x-b'x1x;
  vhat1=(yky1+b'xkx1*b-2*ykx1*b-2*v1x*inxmx*xz*inzkz*(zky1-zkx1*b)+
         v1x*inxmx*v1x')*s;
  m1=v1v./sqrt(vhat1);
  pvm1=2*cdfnc(abs(m1));
  endif;

#  @------------------ robust statistics for one-step estimates -----------------@
    
    xax2=xzinzkz*zkz2*xzinzkz';

vc2=inxmx*xax2*inxmx;                                     @ covariance matrix @
stderr2=sqrt(diag(vc2));                                    @ standard errors @
trat2=b./stderr2;                                                  @ t-values @
pvt2=2*cdfnc(abs(trat2));                                          @ p-values @

if sc2;                                            @ serial correlation tests @
  if m34;
    v2v=y2y+b'x2x*b-y2x*b-b'x2y;
    v2x=y2x-b'x2x;
    endif;
    vhat2w=vomv2w-(2.*v2x*inxmx*xzinzkz*zomv2w)+(v2x*vc2*v2x');
  m2w=v2v./sqrt(vhat2w);
  pvm2w=2*cdfnc(abs(m2w));
endif;

if sc1;
  if m34;
    v1v=y1y+b'x1x*b-y1x*b-b'x1y;
    v1x=y1x-b'x1x;
                                                 endif;
                                                 vhat1w=vomv1w-(2.*v1x*inxmx*xzinzkz*zomv1w)+(v1x*vc2*v1x');
  m1w=v1v./sqrt(vhat1w);
  pvm1w=2*cdfnc(abs(m1w));
endif;

# @---------------------------- Two-step estimates -----------------------------@

if (iiv==1 and ijus==0);

  if imike==0;
    inzkz2=invpd(zkz2);
  else;
    inzkz2=pinv(zkz2);
  endif;
  clear zkz2;

  bxmx=xz*inzkz2*zx; bxmy=xz*inzkz2*zy;

  bvc=invpd(bxmx);                                        @ covariance matrix @

  bb=bvc*bxmy;                                           @ coefficient vector @
  bstderr=sqrt(diag(bvc));                                  @ standard errors @
  btrat=bb./bstderr;                                               @ t-values @
  bpvt=2*cdfnc(abs(btrat));                                        @ p-values @

  if sc2;                                          @ serial correlation tests @
    bv2v=y2y+bb'x2x*bb-y2x*bb-bb'x2y;
    bvhat2w=vomv2w-(2.*v2x*bvc*xz*inzkz2*zomv2w)+(v2x*bvc*v2x');
    bm2w=bv2v./sqrt(bvhat2w);
    pvbm2w=2*cdfnc(abs(bm2w));
  endif;
                                                 
  if sc1;
    bv1v=y1y+bb'x1x*bb-y1x*bb-bb'x1y;
    bvhat1w=vomv1w-(2.*v1x*bvc*xz*inzkz2*zomv1w)+(v1x*bvc*v1x');
    bm1w=bv1v./sqrt(bvhat1w);
    pvbm1w=2*cdfnc(abs(bm1w));
  endif;
