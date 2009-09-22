c $Header:$
cmf 27-06-03 find coulomb tail of ionic pseudopotential (epstail)
cmf 04-09-02 added kli exchange functional
c***********************************************************************
c normconserving pseudopotentials
c gradient-corrections
c nonlinear core-valence XC
c self-interaction correction (definitely experimental!)
c
c input
c im ............ channel no. for monitoring output
c irl ........... relativity mode of all-electron calculation
c                 1 = scalar-relativistic, 2 = non-relativistic
c iexc .......... exchange-correlation mode
c zin ........... full nuclear charge
c zc ............ number of core electrons
c zv ............ number of valence electrons
c lmax .......... highest angular momentum quantum number to generate
c                 pseudopotentials for 
c                 lmax < 4
c t_pp_mode() ... generate pseudowavefunction 
c                 .false. = using the eigenstate
c                 .true.  = from input reference energy (Hamann procedure)
c i_pp_type() ... pseudopotential scheme
c                 1 = Hamann, 2 = Troullier-Martins 
c np() .......... principal quantum number of reference state
c e() ........... eigenvalues
c f() ........... occupancy
c rc() .......... pseudoization cutoff radii
c v() ........... effective potential
c mmax .......... maximum used index radial grid 
c r() ........... logarithmic radial grid
c rnlc .......... partial core cutoff radius
c dc() .......... core charge density
c dcp() ......... dto. 1st derivative
c dcpp() ........ dto. 2nd derivative
c tsic .......... flag for self-interaction corrected pseudopotentials
c vorb() ........ orbitwise sic corrections to the effective potential
c
c output
c dc() .......... model core charge density (if any)
c dcp() ......... dto.
c dcpp() ........ dto.
c ups() ......... pseudo wavefunctions
c
c files
c fort.[iu] ..... monitoring data
c fort.25 ....... pseudo valence density (integrated over solid angle)
c fort.27 ....... model core density (integrated over solid angle)
c fort.28 ....... monitoring unscreening 
c fort.4[0,1,...] ionic pseudopotentials & wavefunctions for l=0,1,...
c fort.4[5,6,...] screened pseudopotentials 
c
c Martin Fuchs, FHI der MPG, Berlin, 09-1996
c***********************************************************************

      subroutine ncpp(iu,irl,iexc,zin,zc,zv,lmax,t_pp_mode,i_pp_type
     1,  np,e,f,rc,v
     1,  mmax,r
     1,  rnlc,dc,dcp,dcpp,rho,rhop,rhopp
     1,  tsic,vorb,ups
     1,  nval,nae,lae,eae,fae,svkli)

      implicit none
      include  'parameter.h'
      include  'default.h'
   
      real*8    epstail
      parameter (epstail=1.d-3)

      character spp*20,symz*8,symxc*26
      logical   t_pp_mode(ms)
      integer   iu,i,iexc,j,l,l1,lmax,mch,mchf,mmax,mode,nin,ninmx,nrc
      integer   ir,irl,iequidensity,io
      integer   np(ms),i_pp_type(ms),ninu(ms),modmap(7)
      real*8    al,als,amesh,ekin,ex,ec,eeel,eeig,eexc,etest,gam,gpr
      real*8    umch,uld,etot,evps,zin,zc,zv,cval,tc,t1,t2,fmom,dmelm
      real*8    rnlc,wm,upmch,uldf,sgrad,pi4,rcut_global,raetail
      real*8    f(ms),e(ms),rc(ms)
      real*8    r(mx),rho(mx),rhop(mx),rhopp(mx),u(mx),up(mx),upp(mx)
      real*8    v(mx),vin(mx),ves(mx),vxc(mx),dc(mx),dcp(mx),dcpp(mx)
      real*8    drh(mx,ms),drhp(mx,ms),drhpp(mx,ms),ups(mx,ms)
      real*8    vps(mx,ms),wa(mx),wb(mx),uae(mx,ms),vscr(mx,ms)
      real*8    vunscr(mx),vl,vh,bl,bh

      integer   nval,nae(ms),lae(ms)
      real*8    eae(ms),fae(ms)

      logical   tsic,tdo(10)
      integer   icut(10)
      real*8    rshift,rcut,e_sic,epot_sic,epot_sic_ref,emax,fwidth
      real*8    vesi(mx),vxc_spin(mx,2),dorb(mx,2),dnil(mx,2)
      real*8    esic(5),vorb(mx,5),vsic(mx,5),wsic(mx,5),fcut(mx)
ckli
      logical   tkli
      integer   iukli,opt_kli_shift,opt_kli_unscreen
      integer   index(ms),mp(ms),norb_kli,lp(ms)
      integer   lp_(ms),mp_(ms),np_(ms)
      real*8    fp(ms),fp_(ms),dens(mx),svm(ms),vsl(mx),rho_ae(mx)
      real*8    vcp(mx),vdummy(mx),svkli(ms),svkli_(ms),utmp(mx,ms)
      real*8    vckli(mx)
      real*8    range,rwidth,eshift
      parameter (iukli=80)
cspecial for l == 3
      real*8    alpha,beta,vir1

csic exchange-correlation energy pp
      real*8    des
      common/xcpp/des(mx)

      external  fmom,dmelm,spp
      data (modmap(i),i=1,7)/2,2,4,4,5,7,7/

      al=0.1d0*log(r(11)/r(1))
      als=al*al
      amesh=exp(al)
      tkli=.false.
      if(iexc.eq.12 .or. iexc.eq.13 .or. 
     1   iexc.eq.18 .or. iexc.eq.19) then
        tkli=.true.
        opt_kli_unscreen=2
        if((iexc.eq.12 .or. iexc.eq.13) .and. rnlc.ne.0.d0) then
          write(ie,'(3a,i3,a)') '%ncpp: ERROR - nonlinear core-valence',
     1     ' exchange-correlation (input rnlc > 0) not implemented',
     1     ' for iexc =',iexc,'. Stopping.'
          stop
        endif
      endif

      if(tsic) stop 'ncpp - sic not implemented'
      do i=1,10
        tdo(i)=.true.
        icut(i)=0
      enddo
      wm=1.d0
      do i=1,mmax
        wa(i)=0.d0
        wb(i)=1.d0
        dnil(i,1)=0.d0
        dnil(i,2)=0.d0
        dorb(i,2)=0.d0
        vxc_spin(i,2)=0.d0
      enddo

c save full potential
      call dcpv(mx,1,mmax,v,vin)

c find ionic cutoff radius: 
c unscreen full potential 
      do ir=1,mmax
        rho_ae(ir)=rho(ir)
      enddo
      do i=mmax,1,-1
        if(rho(i) .gt. dc(i)) iequidensity=i
      enddo
      call vestat(mmax,zv,eeel,r,rho,ves,.false.)
      if(rnlc .gt. 0.0) then
        do i=mmax,1,-1
          rho(i)=rho(i)+dc(i)
          rhop(i)=rhop(i)+dcp(i)
          rhopp(i)=rhopp(i)+dcpp(i)
        enddo
      endif
      call vexcor(iexc,mmax,r,rho,rhop,rhopp,vxc,tc,ex,ec,.false.)
c
c     open(11,file='tmp.unscr',status='unknown')
c     do i=1,mmax
c      vunscr(i)=vin(i)-vxc(i)-ves(i)
c      write(11,*) r(i),(1.d0-0.*exp(-r(i)**4))*vunscr(i)
c     enddo
c     close(11)

      do i=mmax,1,-1
        tc=abs(1.d0+(vin(i)-vxc(i)-ves(i))*r(i)/(zin-zc))
        if(tc .gt. epstail) goto 100
      enddo
  100 rcut_global=r(i)
      raetail=rcut_global
c switch off the global cutoff
      rcut_global=r(mmax)+1.0

c null eigenvalue energy accumulator
      ninmx=0
      eeig=0.0d0
      evps=0.0d0

c null charge accumulator & derivatives
      cval=0.d0
      do i=1,mmax
        rhopp(i)=0.d0
        rhop(i)=0.d0
        rho(i)=0.d0
      enddo

      write(iu,630)
  630 format(/'    === pseudo atom ===',
     1       //'  l  type  rcore       rmatch          eigenvalue(eV)',
     1       '      norm test   slope test',
     1      /'                                 all-electron     pseudo',
     1       '     1 =         1 =')

c loop to construct pseudopotential for each state
c     do j=1,mmax
c      rho_ae(j)=0.d0
c     enddo

      do l1=1,lmax+1
        l=l1-1
        etest=e(l1)
        nrc=log(rc(l1)/r(1))/al+1
        if(.not. t_pp_mode(l1)) then
          mode=irl
        else if( t_pp_mode(l1) .and. f(l1) .eq. 0.d0) then
          mode=3*irl
          if(i_pp_type(l1) .eq. 1) mch=log(2.5*rc(l1)/r(1))/al
          if(i_pp_type(l1) .eq. 2) mch=nrc+1
        endif

c full potential solution
        call dftseq(mode,zin,mmax,r,np(l1),l,wm,v,wa,wb,
     1    ninu(l1),mch,uldf,e(l1),uae(1,l1),up,upp)

c       if(f(l1) .gt. 0.0) then
c         do j=1,ninu(l1)
c           rho_ae(j)=rho_ae(j)+f(l1)*(uae(j)/r(j))**2
c         enddo
c       endif

        ninmx=max(ninmx,ninu(l1))

        if(i_pp_type(l1) .eq. 1) mchf=mch
        if(i_pp_type(l1) .eq. 2) mchf=nrc+1
        uldf=up(mchf)/uae(mchf,l1)
        upmch=up(mchf)
        umch=uae(mchf,l1)

        mode=modmap(mode)

d       write(ie,*) '& ncpp -- mode,t_pp_mode ',mode,t_pp_mode(l1)
d       write(ie,*) '& ncpp -- ae wfct        ',l,etest,e(l1)
d       write(ie,*) '& ncpp -- mchf,umch,upmch',mchf,umch,upmch
d       write(ie,*) '& ncpp -- i_pp_type      ',i_pp_type(l1)

c call potential subroutines
        call dcpv(mx,1,mmax,uae(1,l1),u)
        if(i_pp_type(l1) .eq. 1) then
          call hamann(iu,0.d0,epspp,l,rc(l1),e(l1),mode,mchf,mmax,r,
     1      u,up,upp,v,vps(1,l1))
        else if(i_pp_type(l1) .eq. 2) then
          call tromar(l,rc(l1),nrc,e(l1),mmax,r,u,up,upp,v,vps(1,l1))
        else
          stop 'ncpp - pseudopotential scheme not implemented'
        endif

        if(.false. .and. l == 3 ) then
          write(ie,*) 'ncpp - warning: special l = 3 construction'
          write(iu,*) 'ncpp - warning: special l = 3 construction'
          vir1=(v(nrc-2)+8*(v(nrc+1)-v(nrc-1))-v(nrc+2))/12.0d0
          vir1=vir1/(al*r(nrc))
          beta = -vir1/(2.d0*r(nrc)*v(nrc))
          alpha = v(nrc)*exp(beta*r(nrc)*r(nrc))
          write(ie,*) 'ncpp - alpha,beta',alpha,beta
          do i=1,nrc
            vps(i,l1) = alpha*exp(-beta*r(i)*r(i))
          enddo
          do i=nrc+1,mmax
            vps(i,l1) = v(i)
          enddo
          do i = nrc-10,nrc+10
            print *, r(i),vps(i,l1),v(i)
          enddo
        endif

        call dcpv(mx,1,mmax,u,ups(1,l1))

c accumulate charge, derivatives & eigenvalues
        cval=cval+f(l1)
        if(f(l1) .gt. 0.0) then
          do j=1,ninu(l1)
            drh(j,l1)=(u(j)/r(j))**2
            rho(j)=rho(j)+f(l1)*drh(j,l1)
          enddo
          do j=1,ninu(l1)
            drhp(j,l1)=2.d0*(up(j)/al-u(j))*u(j)/(r(j)**3)
            drhpp(j,l1)=2.d0*((up(j)/al-u(j))*(up(j)/al
     1               -3.d0*u(j))+(upp(j)-al*up(j))*u(j)/als )/r(j)**4
            rhop(j)=rhop(j)+f(l1)*drhp(j,l1)
            rhopp(j)=rhopp(j)+f(l1)*drhpp(j,l1)
          enddo
        endif

        mch=mchf
c @@@@@@@@@@@@@@
c
c This is to use the s pseudo also in the p channel
c
        do j=1,mmax
          vps(j,2) = vps(j,1)
        enddo
c @@@@@@@@@@@@@@
        call dftseq(mode,0.d0,mmax,r,l+1,l,wm,vps(1,l1),wa,wb,
     1    nin,mch,uldf,etest,u,up,upp)
        call dcpv(mx,1,mmax,u,ups(1,l1))
        eeig=eeig+f(l1)*e(l1)
        evps=evps+f(l1)*dmelm(mmax,al,r,vps(1,l1),drh(1,l1))
c
c diagnostic output
        gam=abs(umch/u(mchf))
        gpr=abs(upmch/up(mchf))
        rc(l1)=r(nrc)

        write(iu,'(i3,2x,1a1,(2x,6f12.7))') 
     1   l,spp(i_pp_type(l1)),rc(l1),r(mchf),e(l1)*ry2,etest*ry2,gam,gpr

c excess occupied states
        do i=1,nval
          
          if(l.eq.lae(i) .and. np(l1).ne.nae(i)) then

            if(tkli) then
              write(ie,'(a,a)') '%ncpp: ERROR - excess occupied states',
     1           ' not implemented for KLI exchange option. Stopping.'
              stop
            endif

            cval=cval+fae(i)
            etest=eae(i)
            call dftseq(2,0.d0,mmax,r,l+2,l,wm,vps(1,l1),wa,wb,
     1        nin,mch,uldf,etest,u,up,upp)
            gam=etest
c for consistent unscreening: augment density w/ full potential wavefunction
            call dftseq(irl,0.d0,mmax,r,nae(i),l,wm,v,wa,wb,
     1        nin,mch,uldf,gam,u,up,upp)

            do j=1,nin
              rhopp(j)=rhopp(j)+fae(i)*2.d0*((up(j)/al-u(j))*(up(j)/al
     1               -3.d0*u(j))+(upp(j)-al*up(j))*u(j)/als )/r(j)**4
              rhop(j)=rhop(j)+fae(i)*2.d0*(up(j)/al-u(j))*u(j)/(r(j)**3)
              u(j)=(u(j)/r(j))**2
              rho(j)=rho(j)+fae(i)*u(j)
            enddo

            eeig=eeig+fae(i)*etest
            evps=evps+fae(i)*dmelm(mmax,al,r,vps(1,l1),u)

            write(iu,'(i3,2x,1a,26x,2f12.7,3x,16a)') 
     1      l,spp(i_pp_type(l1)),eae(i)*ry2,etest*ry2,'(excited state!)'

          endif

        enddo

      enddo

c pseudopotentials done

c test: 1 = Integral(rho)/charge
      tc=fmom(0,mmax,al,1.d0,r,rho)

c test: 1 = -Integral(r*rho')/(3*charge)
      t1=-fmom(1,mmax,al,cval,r,rhop)/3.d0

c test: 1 = Integral(r**2 * rho'')/(12*charge)
      t2=fmom(2,mmax,al,cval,r,rhopp)/12.d0

      if(abs(cval-zv) .gt. 1e-12) 
     1  write(ie,*) '& ncpp -- warning: valency mismatch?',zv,cval

c potential energy (before unscreening)
      ekin=eeig-evps 

c unscreen pseudopotentials
      call outcore(25,mmax,1,r,rho,rhop,rhopp)
      if(cval .gt. 0.0) then
        
        call vestat(mmax,cval,eeel,r,rho,ves,.true.)
        if(rnlc .gt. 0.0) then

c nonlinear core correction
c can subject core + pseudo valence density the partial core construction TOGETHER and
c obtain very smoothly varying (GGA) pseudopotentials (restricting the shape of the
c total pseudo density)
c this option is taken out as it may lead to different pseudo core densities for different
c valence configurations and hence to transferability problems when PPs from different
c configurations are combined 
c         do i=1,mmax
c           dc(i)=rho(i)+dc(i)
c           dcp(i)=rhop(i)+dcp(i)
c           dcpp(i)=rhopp(i)+dcpp(i)
c         enddo

          do i=mmax,1,-1
            if(dc(i)-rho(i)*0.0 .gt. rho(i)) goto 81
          enddo
   81     iequidensity=i+1
          if(rnlc .gt. r(1)) call dnlcc7(mmax,rnlc,r,rho,dc,dcp,dcpp,24)
c
c         open(11,file='tmp.core')
c         do i=1,mmax
c           read(11,*,end=88) r(mx),dc(i),dcp(i),dcpp(i)
c         enddo
c         close(11)
c
c         write(23,*) 'warning - nlcv xc uses external core density'
c  88     continue
c

          do i=1,mmax
c output pseudo core density
c when restoring the core + pseudo valence density partial core construction version
c           write(27,'(e20.14,3(1x,e20.14))') r(i),dc(i)-rho(i),dcp(i)-rhop(i),dcpp(i)-rhopp(i)
            write(27,'(e20.14,3(1x,e20.14))') r(i),dc(i),dcp(i),dcpp(i)
            dc(i)=dc(i)+rho(i)
            dcp(i)=dcp(i)+rhop(i)
            dcpp(i)=dcpp(i)+rhopp(i)
          enddo
          close(27)
        else
          do i=1,mmax
            dc(i)=rho(i)
            dcp(i)=rhop(i)
            dcpp(i)=rhopp(i)
          enddo
        endif
        if(.not. tkli) then
          call vexcor(iexc,mmax,r,dc,dcp,dcpp,vxc,eexc,ex,ec,.true.)
c         call ecp(zin,mmax,r,vcp)
c         open(11,file='vkli_ps.dat')
c         do i=1,mmax
c           write(11,*) r(i),vcp(i)
c         enddo
c         close(11)
        else
c case of kli pseudopotential
          norb_kli=0
          do i=1,lmax+1
            if( f(i).gt.epstail ) norb_kli=norb_kli+1
            fp(i)=0.5*f(i)
            lp(i)=i-1
            mp(i)=0
          enddo

          if(tkli) write(iu,'(/a,1x,i2)') 
     1     ' %ncpp: kli - number of occupied orbitals = ',norb_kli  

c normally states are order by l1, for vklix they must be
c ordered w.r.t. to increasing eigenvalues, this is done here
          if(tkli) then
c create ordered index distinguishing cases, should be made automatic
c s<p<d
            index(1)=1
            index(2)=2
            index(3)=3
            index(4)=4
            index(5)=5
c use all-electron wfct in unscreening
            opt_kli_unscreen=2
c d<s<p
            if(e(3) .lt. e(1)) then
              index(1)=3
              index(2)=1
              index(3)=2
              index(4)=4
              index(5)=5
c use pseudo wfct in unscreening
              opt_kli_unscreen=3
            endif


            write(iukli,'(a)') '--- state ordering for vklix ---'
            do l1=1,lmax+1
             write(iukli,*) l1,e(l1),index(l1),e(index(l1))
             np_(l1)=np(index(l1))
             lp_(l1)=lp(index(l1))
             mp_(l1)=mp(index(l1))
             fp_(l1)=fp(index(l1))
             svkli_(l1)=svkli(index(l1))
            enddo

          endif

          if(opt_kli_unscreen == 1) then
c using xkli in sic style (not the kli pseudopotential!)
           write(iu,'(1x,a)') '%kli: unscreen in LDA'
           call vexcor(8,mmax,r,dc,dcp,dcpp,vxc,eexc,ex,ec,.true.)

          else if(tkli) then
           if(opt_kli_unscreen == 2) then
c using all-electron valence wavefunctions for unscreening
c the full ae density is used for damping contributions from
c the core region
c           use ae kli shifts:
            opt_kli_shift=2
            write(iu,'(1x,a,a,i2)') 
     1        '%ncpp: kli - unscreen all-electron valence,',
     1        ' opt_kli_shift =',opt_kli_shift
            do l1=1,lmax+1
              do ir=1,mmax 
                utmp(ir,l1)=uae(ir,index(l1))
              enddo
            enddo
            do ir=1,mmax
              dens(ir)=0.5d0*rho_ae(ir)
            enddo
            call vklix(norb_kli,np_,lp_,mp_,fp_,mmax,r,utmp,dens,
     1                vxc,vsl,svkli_,svm,opt_kli_shift,
     1                .true.,.false.,.true.)
c using pseudo valence wavefunctions for unscreening
           else if(opt_kli_unscreen == 3) then 
c           use pseudo kli shifts:
            opt_kli_shift=1
            write(iu,'(1x,a,a,i2)') 
     1        '%ncpp: kli - unscreen pseudo valence', 
     1        ' opt_kli_shift =',opt_kli_shift
            do l1=1,lmax+1
              do ir=1,mmax 
                utmp(ir,l1)=ups(ir,index(l1))
              enddo
            enddo
            do ir=1,mmax
              dens(ir)=0.5d0*rho(ir)
            enddo
            call vklix(norb_kli,np_,lp_,mp_,fp_,mmax,r,utmp,dens,
     1                vxc,vsl,svkli,svm,1,
     1                .true.,.false.,.true.)
            do l1=1,lmax+1
              svkli_(index(l1))=svkli_(l1)
            enddo

           else
            stop '%ncpp - invalid opt_kli_unscreen.'
           endif

c kli + lda or gga, if applicable, vckli is zeroed by vexcor
           call vexcor(iexc,mmax,r,dc,dcp,dcpp,vckli,eexc,ex,ec,.true.)
d          write(ie,*) '%ncpp: kli+lda',iexc,vckli(1),vckli(100)
           do ir=1,mmax
             vxc(ir)=vxc(ir)+vckli(ir)
           enddo
           
           write(iukli,'(a)') ' --- kli pseudo shifts (Ha) --- '
           write(iukli,'(a,i2)') 'opt_kli_unscreen =',opt_kli_unscreen
           do i=1,norb_kli
            write(iukli,'(a2,3i3,1x,e14.8)')
     1        '%',i,np_(i),lp_(i),svkli_(i)
           enddo
           write(iukli,*)

          endif
          

          if(iexc .eq. 13) then
c adding the core polarization potential
            write(iu,'(/,1x,a)') '%ncpp: add cp potential.'
            call ecp(zin,mmax,r,vcp)
            do ir=1,mmax
              vxc(ir) = vxc(ir) - vcp(ir)
            enddo
          endif

          do ir=1,mmax
            vcp(ir) = 0.d0
          enddo
          if(opt_kli_unscreen == 1) then
            write(6,*) 'Enter shifting range (0 < range < 20) ...'
            read(5,*) range
            if(range .gt. 0.d0 .and. range .lt. 20.d0) then
             write(6,*) '%ncpp: apply global shifting range.'
             write(iu,'(a,1x,e14.6)') '%ncpp: apply shifting range=',
     1         range
             eshift=1.d0/range
             rwidth=range/50.d0
             do i=1,mmax
             vcp(i) = eshift/(exp((r(i)-range)/rwidth)+1.d0)
     &         +1.d0/(exp(-(r(i)-range)/rwidth)+1.d0)/r(i)
             enddo
             do ir=1,mmax
               vxc(ir) = vxc(ir) - vcp(ir)            
             enddo
            endif
          endif
          open(11,file='vkli_ps.dat')
          write(11,'(a)') '#ps r(i) vxc(i) vslater(i) vcp(i)'
          do i=1,mmax
            write(11,*) r(i),vxc(i),vsl(i),vcp(i)
          enddo
          close(11)
        endif

        call dadv(mx,1,mmax,ves,vxc,vin)
        do l1=1,lmax+1
          do i=1,mmax
            vscr(i,l1)=vps(i,l1)
            vps(i,l1)=vps(i,l1)-vin(i)
          enddo
c enforce coulomb like tail beyond rcut_global
          do i=mmax,1,-1
            if(r(i) .lt. rcut_global) goto 110
            vps(i,l1)=vunscr(i)+exp(-2*(r(i)-rcut_global)**2)
     1              *(vps(i,l1)-vunscr(i))
            if(i .ne. mmax) 
     1        write(ie,*) '& ncpp - warning: global tail cutoff applied'
          enddo
  110     continue  
        enddo

      endif

c special:
c     call overlap(0,mx,r,vunscr)
c     call overlap(1,mx,r,vps(1,2))
c     call overlap(2,mx,r,vps(1,3))

      if(tkli) then
       do ir=1,mmax
         dens(ir)=0.5d0*rho(ir)
       enddo
       call vklix(norb_kli,np,lp,mp,fp,mmax,r,ups,dens,
     1              vdummy,vdummy,svkli,svm,2,
     1              .false.,.true.,.false.)
       ex=0.d0
       do i=1,norb_kli
        ex=ex+fp(i)*svm(i)
       enddo
      endif

c external potential energy (i.e. unscreened pseudopotentials)
      evps=evps-dmelm(mmax,al,r,vin,rho)
      etot=ekin+.5*eeel+evps+ex+ec

c pseudoatom output
      write(iu,*)
      if(rcut_global .lt. r(mmax))
     1  write(iu,608)'--- global tail cutoff beyond',rcut_global
      if(rnlc .gt. 0.) then
        write(iu,609)'--- nonlinear core-valence XC ---'
        write(iu,608)'core density cutoff radius',rnlc
        write(iu,608)'c-v equidensity radius',r(iequidensity)
        write(iu,608)'integrated model core density',
     1    fmom(0,mmax,al,1.d0,r,dc)-tc
        write(iu,608)'true core charge',zc
      endif
      if(rnlc .eq. 0.) then
        write(iu,609)'--- linearized core-valence XC ---'
        write(iu,608)'c-v equidensity radius',r(iequidensity)
      endif

      write(iu,631)
      write(iu,608) 'total energy',etot
      write(iu,608) 'kinetic energy',ekin
      write(iu,608) 'potential energy',evps
      write(iu,608) 'hartree energy',0.5*eeel
      write(iu,608) 'xc energy',ex+ec
      write(iu,608) ' c energy',ec
      write(iu,608) 'integrated valence density',tc
c
  631 format(/34x,'(Hartree a.u.)')
  608 format(a30,2x,f16.5)
  609 format(a36)
  632 format(a20,f18.6)
  633 format(a20,f18.6,2x,f18.6)
  634 format(a30,f16.3,5x,a,1p,e8.1,a)
      write(iu,608)'... 1st derivative test 1 =',t1
      write(iu,608)'... 2nd derivative test 1 =',t2
      write(iu,634)'coulomb tail beyond radius', raetail, 
     1  '(tolerance ',epstail,')'  
c
  635 format(a20,2x,i1,2x,a20)

c output pseudopotentials & pseudowavefunctions
      call outpot(40,mx,mmax,lmax+1,zin-zc,rc,r,ups,vps,vscr)
      call out39(mx,lmax,ninu,np,rc,r,ups,uae)
c for plot x y range
      vl=1.d20
      vh=-1.d20
      do l1=1,lmax+1
        call dextv(mx,1,mmax,vps(1,l1),bl,bh)
        vl=min(bl,vl)
        vh=max(bh,vh)
      enddo
      i=max(-100,int(vl-int(vl/5)*5)-1+int(vl/5)*5)
      j=min(50,int(vh-int(vh/5)*5)+1+int(vh/5)*5)
      l1=max(1,(j-i)/4)
      write(iu,637) 'y range plot',i,j,l1
  637 format(a30,6x,3i4)

c output for monitoring unscreening in case of GGAs
      pi4=16*atan(1.0)
      do i=1,mmax,max(mmax/100,1)
        sgrad=0.16*abs(dcp(i)/pi4)/(max(dc(i),1d-12)/pi4)**1.33
        write(28,'(5(e12.6,1x))') r(i),dc(i),sgrad,ves(i),vxc(i)
      enddo
      close(28)
      do i=1,mmax
        write(95,*) r(i),(f(3)-1.d0)*drh(i,3)
     1,(f(3)-1.d0)*drhp(i,3),(f(3)-1.d0)*drhpp(i,3)
      enddo

      return
      end
c
