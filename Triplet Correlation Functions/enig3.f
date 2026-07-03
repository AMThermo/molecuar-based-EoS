C *************************************************************************************
      program enig3
        implicit none
        character(len=6) :: string
        integer :: force, run, ncell, ncycle, nequi, nblock, i, j, n
        integer :: nconf, nacct, nove, ix, iy, iz
        double precision :: pi, temp, rho, lam, lama, lamr, lbox, hlbox
        double precision :: cut, maxds, xa(4), ya(4), za(4), sig
        double precision :: dx, dy, dz, r2, C, upot, sig2, invt, cut2
        double precision :: hlamr, hlama, lam2, invr2
        double precision, dimension(:), allocatable :: rx0, ry0, rz0
        common /isimvar/ force, run, n, ncycle, nequi, nblock, nconf
        common /denergy/ lam, lamr, lama, C, upot
        common /dsimvar/ temp, rho, cut, maxds, sig, lbox, hlbox
        common /factors/ sig2, invt, cut2, lam2, hlamr, hlama
        common /constants/ pi
        pi=dacos(-1.d0)
        open(unit=1, file="mc.in.txt", status="old")
        read(1,*) string, force
        read(1,*) string, ncell
        read(1,*) string, temp
        read(1,*) string, rho
        read(1,*) string, lam
        read(1,*) string, lamr, lama
        read(1,*) string, cut
        read(1,*) string, string, string, ncycle, nequi, nblock
        read(1,*) string, maxds
        read(1,*)
        read(1,*) (xa(i), ya(i), za(i), i=1,4)
        close(1)
        n=4*ncell*ncell*ncell
        allocate(rx0(n), ry0(n), rz0(n))
        nconf=(ncycle-nequi)/nblock
        lbox=(dble(n)/rho)**(1.0d0/3.0d0)
        hlbox=0.5d0*lbox
        sig=1.0d0/lbox
        maxds=maxds*sig
        invt=1.0d0/temp
        sig2=sig*sig
        if (force.eq.2) then
            lam=lam*sig
            lam2=lam*lam
        else if (force.eq.3) then
            cut=cut*sig
            if (cut.ge.0.5d0) then
                write(*,*) "INVALID CUT-OFF: HIGHER THAN HALF THE BOX!"
                stop
            end if
            cut2=cut*cut
            hlamr=lamr*0.5d0
            hlama=lama*0.5d0
            C=lamr/(lamr-lama)*(lamr/lama)**(lama/(lamr-lama))
        end if
        ! ------------------- INIT. FCC LATTICE ------------------
        i=1
        do j=1,4
            do ix=1,ncell
                do iy=1,ncell
                    do iz=1,ncell
                        rx0(i)=((ix-1)+xa(j))/ncell
                        ry0(i)=((iy-1)+ya(j))/ncell
                        rz0(i)=((iz-1)+za(j))/ncell
                        i=i+1
                    end do
                end do
            end do
        end do
        ! ------------------- INIT. TOTAL POT. ENERGY ------------------
        nove=0
        upot=0.0d0
        do i=1,n-1
            do j=i+1,n
                dx=rx0(j)-rx0(i)
                dy=ry0(j)-ry0(i)
                dz=rz0(j)-rz0(i)
                ! apply the MIC:
                dx=dx-dnint(dx)
                dy=dy-dnint(dy)
                dz=dz-dnint(dz)
                r2=dx*dx+dy*dy+dz*dz
                if (force.eq.1) then
                    if (r2.lt.sig2) then
                        write(*,*) "INVALID INIT.: OVLP. PARTICLES!"
                        stop
                    else
                        upot=upot+0.0d0
                    end if
                else if (force.eq.2) then
                    if (r2.lt.cut2) then
                        write(*,*) "INVALID INIT.: OVLP. PARTICLES!"
                        stop
                    else if (r2.lt.lam2) then
                        upot=upot-1.0d0
                    end if
                else if (force.eq.3) then
                    if (r2.lt.cut2) then
                        invr2=sig2/r2
                        upot=upot+C*invr2**hlamr
                        upot=upot-C*invr2**hlama
                    end if
                end if
            end do
        end do
        open(unit=2, file="mc.inc.txt")
        write(2,*) force
        write(2,*) n
        write(2,*) temp
        write(2,*) rho
        write(2,*) lam
        write(2,*) lamr
        write(2,*) lama
        write(2,*) C
        write(2,*) cut
        write(2,*) ncycle
        write(2,*) nequi
        write(2,*) nblock
        write(2,*) maxds
        write(2,*) lbox
        write(2,*) sig
        close(2)
        open(unit=7, file="upot.dat")
        open(unit=8, file="ratio.dat")
        open(unit=31, file="conrdf_1.dat")
        open(unit=32, file="conrdf_2.dat")
        open(unit=33, file="conrdf_3.dat")
        open(unit=34, file="conrdf_4.dat")
        open(unit=35, file="conrdf_5.dat")
        open(unit=41, file="contcf_1.dat")
        open(unit=42, file="contcf_2.dat")
        open(unit=43, file="contcf_3.dat")
        open(unit=44, file="contcf_4.dat")
        open(unit=45, file="contcf_5.dat")
        call metropolis(rx0, ry0, rz0)
        close(7)
        close(31)
        close(32)
        close(33)
        close(34)
        close(35)
        close(41)
        close(42)
        close(43)
        close(44)
        close(45)
        call statistics
        stop
      end

      function ran2(idum)
        implicit double precision(a-h,o-z)
        parameter(im1=2147483563, im2=2147483399, am=1.0d00/im1)
        parameter(imm1=im1-1, ia1=40014, ia2=40692, iq1=53668)
        parameter(iq2=52774, ir1=12211, ir2=3791, ntab=32)
        parameter(ndiv=1+imm1/ntab, eps=1.2d-14, rnmx=1.0d00-eps)
        dimension iv(ntab)
        save iv, iy, idum2
        data idum2/123456789/, iv/ntab*0/, iy/0/
        if (idum.le.0) then
            idum=max(-idum,1)
            idum2=idum
            do j=ntab+8,1,-1
                k=idum/iq1
                idum=ia1*(idum-k*iq1)-k*ir1
                if (idum.lt.0) idum=idum+im1
                if (j.le.ntab) iv(j)=idum
            end do
            iy=iv(1)
        end if
        k=idum/iq1
        idum=ia1*(idum-k*iq1)-k*ir1
        if (idum.lt.0) idum=idum+im1
        k=idum2/iq2
        idum2=ia2*(idum2-k*iq2)-k*ir2
        if (idum2.lt.0) idum2=idum2+im2
        j=1+iy/ndiv
        iy=iv(j)-idum2
        iv(j)=idum
        if (iy.lt.1) iy=iy+imm1
        ran2=min(am*iy, rnmx)
        return
      end

      function dva(r, s, t)
        implicit none
        double precision :: r, s, t, dva, pi
        common /constants/ pi
        dva=r**4-9.0d00*r**2*(s**2+t**2)+16.0d00*r*(s**3+t**3)
        dva=dva-9.0d00*s**4+18.0d0*s**2*t**2-9.0d0*t**4
        dva=pi*r**2*dva/72.0d00
        return
      end

      function dvb(r, s, t, ds)
        implicit none
        double precision :: r, s, t, ds, dvb, pi
        common /constants/ pi
        dvb=-9.0d00*ds**4-36.0d00*ds**3*s-54.0d00*(s*ds)**2
        dvb=dvb+18.0d00*ds**2*t**2-36.0d00*ds*s**3+36.0d00*ds*s*t**2
        dvb=dvb+r**4-9.0d00*r**2*(ds**2+2.0d00*ds*s+s**2+t**2)
        dvb=dvb+16.0d00*r*(ds**3+3.0d00*s*ds*(s+ds)+s**3+t**3)
        dvb=dvb-9.0d00*s**4+18.0d00*s**2*t**2-9.0d00*t**4
        dvb=pi*r**2/72.d00*dvb
        return
      end

      function dvc(r, s, t, dt)
        implicit none
        double precision :: r, s, t, dt, dvc, pi
        common /constants/ pi
        dvc=-9.0d00*dt**4-36.0d00*dt**3*t-54.0d00*dt**2*t**2
        dvc=dvc+18.0d00*dt**2*s**2-36.0d00*dt*t**3+36.0d00*dt*t*s**2
        dvc=dvc+r**4-9.0d00*r**2*(dt**2+2.0d00*dt*t+s**2+t**2)
        dvc=dvc+16.0d00*r*(dt**3+3.0d00*t*dt*(t+dt)+s**3+t**3)
        dvc=dvc-9.0d00*s**4+18.0d00*s**2*t**2-9.0d00*t**4
        dvc=pi*r**2/72.d00*dvc
        return
      end

      function dvd(r, s, t, ds, dt)
        implicit none
        double precision :: r, s, t, ds, dt, dvd, pi
        common /constants/ pi
        dvd=-9.0d00*ds**4-36.0d00*ds**3*s+18.0d00*ds**2*dt**2
        dvd=dvd+36.0d00*ds**2*dt*t-54.0d00*ds**2*s**2+18.d00*ds**2*t**2
        dvd=dvd+36.0d00*ds*dt**2*s+72.0d00*ds*dt*s*t-36.0d00*ds*s**3
        dvd=dvd+36.0d00*ds*s*t**2-9.0d00*dt**4-36.0d00*dt**3*t
        dvd=dvd+18.0d00*dt**2*s**2-54.0d00*(t*dt)**2+36.0d00*dt*s**2*t
        dvd=dvd-36.0d00*dt*t**3+r**4-9.0d00*(r*ds)**2-18.0d00*r**2*s*ds
        dvd=dvd-9.0d00*r**2*(dt**2+2.0d00*t*dt+s**2+t**2)
        dvd=dvd+16.0d00*r*(ds**3+3.0d00*ds**2*s+3.0d00*ds*s**2+dt**3)
        dvd=dvd+16.0d00*r*(3.0d00*dt**2*t+3.0d00*dt*t**2+s**3+t**3)
        dvd=dvd-9.0d00*s**4+18.0d00*s**2*t**2-9.0d00*t**4
        dvd=pi*r**2/72.d00*dvd
        return
      end

      subroutine metropolis(rx, ry, rz)
        implicit none
        integer :: run, n, ncycle, nequi, nblock, i, j, k, nove, iseed
        integer :: nmove, nacct, nconf, conf, force
        double precision :: press, temp, rho, lam, cut, ratio, C, upotn
        double precision :: lamr, lama, maxds, lbox, sig, sig2, hlbox
        double precision :: upot, vpot0, vpot, xnew, ynew, znew, ran2
        double precision :: dx, dy, dz, r2, rx, ry, rz, dvpot, invt
        double precision :: boltz, invr2, hlamr, hlama, lam2, cut2
        common /isimvar/ force, run, n, ncycle, nequi, nblock, nconf
        common /denergy/ lam, lamr, lama, C, upot
        common /dsimvar/ temp, rho, cut, maxds, sig, lbox, hlbox
        common /factors/ sig2, invt, cut2, lam2, hlamr, hlama
        dimension rx(n), ry(n), rz(n)
        iseed=-123456789
        nmove=0
        nacct=0
        do run=1,ncycle
            do k=1,n
                nmove=nmove+1
                i=int(ran2(iseed)*n)+1
                ! ---------------- OLD POT. ENERGY ---------------
                vpot0=0.0d0
                do 100 j=1,n
                    if (j.eq.i) go to 100
                    dx=rx(j)-rx(i)
                    dy=ry(j)-ry(i)
                    dz=rz(j)-rz(i)
                    ! apply the MIC:
                    dx=dx-dnint(dx)
                    dy=dy-dnint(dy)
                    dz=dz-dnint(dz)
                    r2=dx*dx+dy*dy+dz*dz
                    if (force.eq.2) then
                        if (r2.lt.lam2) vpot0=vpot0-1.0d00
                    else if (force.eq.3) then
                        if (r2.lt.cut2) then
                            invr2=sig2/r2
                            vpot0=vpot0+C*invr2**hlamr
                            vpot0=vpot0-C*invr2**hlama
                        end if
                    end if
100             continue
                ! displace particle i:
                xnew=rx(i)+maxds*(ran2(iseed)-0.5d00)
                ynew=ry(i)+maxds*(ran2(iseed)-0.5d00)
                znew=rz(i)+maxds*(ran2(iseed)-0.5d00)
                ! apply PBC:
                xnew=xnew-dble(floor(xnew))
                ynew=ynew-dble(floor(ynew))
                znew=znew-dble(floor(znew))
                ! ---------------- NEW POT. ENERGY ---------------
                nove=0
                vpot=0.0d0
                do 101 j=1,n
                    if (j.eq.i) go to 101
                    dx=rx(j)-xnew
                    dy=ry(j)-ynew
                    dz=rz(j)-znew
                    ! apply the MIC:
                    dx=dx-dnint(dx)
                    dy=dy-dnint(dy)
                    dz=dz-dnint(dz)
                    r2=dx*dx+dy*dy+dz*dz
                    if (force.eq.1) then
                        if (r2.lt.sig2) nove=nove+1
                    else if (force.eq.2) then
                        if (r2.lt.sig2) then
                            nove=nove+1
                        else if (r2.lt.lam2) then
                            vpot=vpot-1.0d00
                        end if
                    else if (force.eq.3) then
                        if (r2.lt.cut2) then
                            invr2=sig2/r2
                            vpot=vpot+C*invr2**hlamr
                            vpot=vpot-C*invr2**hlama
                        end if
                    end if
101             continue
                ! -------- METROPOLIS' ACCEPTANCE CRITERION --------
                if (force.eq.1) then
                    if (nove.eq.0) then
                        nacct=nacct+1
                        rx(i)=xnew
                        ry(i)=ynew
                        rz(i)=znew
                    end if
                else if (force.eq.2) then
                    if (nove.eq.0) then
                        dvpot=vpot-vpot0
                        boltz=dexp(-dvpot*invt)
                        if (boltz.gt.ran2(iseed)) then
                            nacct=nacct+1
                            rx(i)=xnew
                            ry(i)=ynew
                            rz(i)=znew
                            upot=upot+dvpot
                        end if
                    end if
                else if (force.eq.3) then
                    dvpot=vpot-vpot0
                    boltz=dexp(-dvpot*invt)
                    if (boltz.gt.ran2(iseed)) then
                        nacct=nacct+1
                        rx(i)=xnew
                        ry(i)=ynew
                        rz(i)=znew
                        upot=upot+dvpot
                    end if
                end if
            end do
            upotn=upot/dble(n)
            ratio=dble(nacct)/dble(nmove)*100.0d00
            write(7,*) run, upotn
            write(8,*) run, ratio
            do j=1,nconf
                conf=nequi+j*nblock
                if (run.eq.conf) then
                    call pdf(rx, ry, rz)
                end if
            end do
            write(*,*) "RUN:", run
        end do
        return
      end

      subroutine pdf(rx, ry, rz)
        implicit none
        character(len=9), dimension(2,5) :: filename
        integer :: force, run, n, ncycle, nequi, nblock
        integer :: nmove, binh, nbin(5), conf, i, j, l, nconf, w
        integer :: binr, bins, bint, flag(n,n), list(n,n)
        double precision :: press, temp, rho, cut, bin0(5), g3s(5,6)
        double precision :: maxds, lbox, sig, k, g2, ang(5), g2s(5)
        double precision :: dx, dy, dz, d(5), hlbox, pi, rr, nig2, maxr
        double precision :: nig3(6), g3(6), rm, r, s, t, ig2
        double precision :: rx(n), ry(n), rz(n), invfac(5), inm, ig3(6)
        double precision :: dv3a, dv3b, dv3c, dv3d, dva, dvb, dvc, dvd
        double precision :: rrdr, ss, tt, dr, ds, dt, fac(5), invd(5)
        double precision :: dxij, dyij, dzij, dxjk, dyjk, dzjk
        double precision :: dxik, dyik, dzik
        double precision, dimension(:,:), allocatable :: cumhp, dV2
        double precision, dimension(:,:,:), allocatable :: cumht, dV3
        double precision, dimension(:,:), allocatable :: hp
        double precision, dimension(:,:,:), allocatable :: ht
        common /isimvar/ force, run, n, ncycle, nequi, nblock, nconf
        common /dsimvar/ temp, rho, cut, maxds, sig, lbox, hlbox
        common /constants/ pi
        logical initialized
        save initialized
        save d, invd, invfac, filename, nbin, bin0, dV2, dV3, hp, ht
        save cumhp, cumht
        data initialized /.false./
        if (.not. initialized) then
            d(1)=0.010d00
            d(2)=0.020d00
            d(3)=0.040d00
            d(4)=0.050d00
            d(5)=0.0625d00
            ang(1)=60.0d0/180.0d0*pi
            ang(2)=90.0d0/180.0d0*pi
            ang(3)=109.5d0/180.0d0*pi
            ang(4)=120.0d0/180.0d0*pi
            ang(5)=180.0d0/180.0d0*pi
            filename(1,1)="rdf_1.dat"
            filename(1,2)="rdf_2.dat"
            filename(1,3)="rdf_3.dat"
            filename(1,4)="rdf_4.dat"
            filename(1,5)="rdf_5.dat"
            filename(2,1)="tcf_1.dat"
            filename(2,2)="tcf_2.dat"
            filename(2,3)="tcf_3.dat"
            filename(2,4)="tcf_4.dat"
            filename(2,5)="tcf_5.dat"
            do i=1,5
                fac(i)=dsqrt(2.0d0*(1.0d0-dcos(ang(i))))
                invfac(i)=1.0d0/fac(i)
            end do
            do l=1,5
                invd(l)=1.0d0/d(l)
                bin0(l)=int(invd(l))
                nbin(l)=int(3.0d0*invd(l))
            end do
            allocate(hp(5,0:maxval(nbin)-1))
            allocate(ht(5,0:maxval(nbin)-1,6))
            allocate(cumhp(5,0:maxval(nbin)-1))
            allocate(cumht(5,0:maxval(nbin)-1,6))
            allocate(dV2(5,0:maxval(nbin)-1))
            allocate(dV3(5,0:maxval(nbin)-1,6))
            do l=1,5
                ! ------------------ VOLUME ELEMENTS ----------------
                do k=0,nbin(l)-1
                    ! ------------------- PAIR --------------------
                    dV2(l,k)=((k+1)**3-k**3)*d(l)**3
                    dV2(l,k)=4.0d0/3.0d0*pi*dV2(l,k)
                    ! -------------- RIGID TRIPLET ----------------
                    do j=1,5
                        ds=d(l)
                        dt=d(l)
                        dr=d(l)*fac(j)
                        ss=k*ds
                        tt=k*dt
                        rr=k*dr
                        rrdr=rr+dr
                        ! compute volume elements:
                        if (rr.lt.ss+tt) then
                            dv3a=dva(rrdr, ss, tt)
                            dv3a=dv3a-dva(rr, ss, tt)
                        else
                            dv3a=0.0d0
                        end if
                        if (rrdr.lt.ss+tt+ds) then
                            dv3b=dvb(rrdr, ss, tt, ds)
                            dv3b=dv3b-dvb(rr, ss, tt, ds)
                        else
                            dv3b=dvb(ss+tt+ds, ss, tt, ds)
                            dv3b=dv3b-dvb(rr, ss, tt, ds)
                        end if
                        if (rrdr.lt.ss+tt+dt) then
                            dv3c=dvc(rrdr, ss, tt, dt)
                            dv3c=dv3c-dvc(rr, ss, tt, dt)
                        else
                            dv3c=dvc(ss+tt+dt, ss, tt, dt)
                            dv3c=dv3c-dvc(rr, ss, tt, dt)
                        end if
                        dv3d=dvd(rrdr, ss, tt, ds, dt)
                        dv3d=dv3d-dvd(rr, ss, tt, ds, dt)
                        dV3(l,k,j)=dv3a+dv3d-dv3b-dv3c
                        dV3(l,k,j)=4.0d0*pi*dV3(l,k,j)
                    end do
                    ! --------------- FLEX. TRIPLET -----------------
                    dv3a=(k+1)**4-k**4
                    dv3a=dv3a*((k+1)**2-k**2)
                    dv3b=(k+1)**3-k**3
                    dv3b=8.0d0/9.0d0*dv3b**2
                    dv3c=(k+1)**2-k**2
                    dv3c=dv3c**2*k**2
                    dV3(l,k,6)=dv3a+dv3b-dv3c
                    dV3(l,k,6)=pi*pi*d(l)**6*dV3(l,k,6)
                end do
            end do
            cumhp=0
            cumht=0
            initialized=.true.
        end if
        hp=0
        ht=0
        maxr=3.0d0
        ! ------------------- DISTANCE HISTOGRAM -------------------
        ! PAIRS:
        do i=1,n-1
            do j=i+1,n
                dx=rx(j)-rx(i)
                dy=ry(j)-ry(i)
                dz=rz(j)-rz(i)
                ! apply the MIC:
                dx=dx-dnint(dx)
                dy=dy-dnint(dy)
                dz=dz-dnint(dz)
                r=dsqrt(dx*dx+dy*dy+dz*dz)
                r=r*lbox
                if (r.lt.maxr) then
                    flag(i,j)=1
                    flag(j,i)=1
                    do l=1,5
                        binh=int(r*invd(l))
                        hp(l,binh)=hp(l,binh)+2
                        cumhp(l,binh)=cumhp(l,binh)+2
                    end do
                end if
            end do
        end do
        ! TRIPLETS:
        do i=1,n
            do 101 j=1,n-1
                if (flag(i,j).eq.0) go to 101
                do 102 k=j+1,n
                    if (flag(i,k).eq.0) go to 102
                    dxij=rx(j)-rx(i)
                    dyij=ry(j)-ry(i)
                    dzij=rz(j)-rz(i)
                    ! apply the MIC:
                    dxij=dxij-dnint(dxij)
                    dyij=dyij-dnint(dyij)
                    dzij=dzij-dnint(dzij)
                    r=dsqrt(dxij*dxij+dyij*dyij+dzij*dzij)
                    r=r*lbox
                    dxjk=rx(k)-rx(j)
                    dyjk=ry(k)-ry(j)
                    dzjk=rz(k)-rz(j)
                    ! apply the MIC:
                    dxjk=dxjk-dnint(dxjk)
                    dyjk=dyjk-dnint(dyjk)
                    dzjk=dzjk-dnint(dzjk)
                    s=dsqrt(dxjk*dxjk+dyjk*dyjk+dzjk*dzjk)
                    s=s*lbox
                    dxik=rx(k)-rx(i)
                    dyik=ry(k)-ry(i)
                    dzik=rz(k)-rz(i)
                    ! apply the MIC:
                    dxik=dxik-dnint(dxik)
                    dyik=dyik-dnint(dyik)
                    dzik=dzik-dnint(dzik)
                    t=dsqrt(dxik*dxik+dyik*dyik+dzik*dzik)
                    t=t*lbox
                    do l=1,5
                        binr=int(r*invd(l))
                        bint=int(t*invd(l))
                        if (binr.eq.bint) then
                            ht(l,binr,6)=ht(l,binr,6)+2
                            cumht(l,binr,6)=cumht(l,binr,6)+2
                            do w=1,5
                                bins=int(s*invd(l)*invfac(w))
                                if (bins.eq.bint.and.binr.eq.bins) then
                                    ht(l,binr,w)=ht(l,binr,w)+2
                                    cumht(l,binr,w)=cumht(l,binr,w)+2
                                end if
                            end do
                        end if
                    end do
102             continue
101         continue
        end do
        conf=(run-nequi)/nblock
        do l=1,5
            open(unit=l+10, file=filename(1,l))
            open(unit=l+20, file=filename(2,l))
            do i=0,nbin(l)-1
                rm=(dble(i)+0.5d0)*d(l)
                ! ----------------- RDF ------------------
                nig2=dV2(l,i)*rho
                ig2=dble(hp(l,i))/dble(n)/nig2
                g2=dble(cumhp(l,i))/dble(n)/nig2/dble(conf)
                ! ----------------- TCF ------------------
                nig3=dV3(l,i,:)*rho*rho
                ig3=dble(ht(l,i,:))/dble(n)/nig3
                g3=dble(cumht(l,i,:))/dble(n)/nig3/dble(conf)
                write(l+10,*) rm, g2
                write(l+20,*) rm, g3
                if (i.eq.bin0(l)) then
                    write(l+30,*) ig2 ! RDF at contact
                    write(l+40,*) ig3 ! TCF at contact
                end if
            end do
            close(l+10)
            close(l+20)
        end do
        return
      end

      subroutine statistics
        implicit none
        character(len=9), dimension(2,5) :: filename
        integer :: force, run, n, ncycle, nequi, nblock, i, l, nbin(5)
        integer :: bin0(5), nconf
        double precision :: temp, rho, cut, maxds, sig, lbox, hlbox
        double precision :: d(5), z, g2, g3(6), g2s(5), g3s(5,6)
        double precision :: Gam2LA(5,6), Gam2SA(5,6), lam2(5,6)
        common /isimvar/ force, run, n, ncycle, nequi, nblock, nconf
        common /dsimvar/ temp, rho, cut, maxds, sig, lbox, hlbox
        d(1)=0.01d0
        d(2)=0.02d0
        d(3)=0.04d0
        d(4)=0.05d0
        d(5)=0.0625d0
        filename(1,1)="rdf_1.dat"
        filename(1,2)="rdf_2.dat"
        filename(1,3)="rdf_3.dat"
        filename(1,4)="rdf_4.dat"
        filename(1,5)="rdf_5.dat"
        filename(2,1)="tcf_1.dat"
        filename(2,2)="tcf_2.dat"
        filename(2,3)="tcf_3.dat"
        filename(2,4)="tcf_4.dat"
        filename(2,5)="tcf_5.dat"
        open(unit=51, file="con_rdf.dat")
        open(unit=52, file="con_tcf.dat")
        open(unit=53, file="con_G2la.dat")
        open(unit=54, file="con_G2sa.dat")
        open(unit=55, file="con_lam2.dat")
        do l=1,5
            bin0(l)=int(1.0d0/d(l))
            nbin(l)=int(3.0d0/d(l))
            open(unit=l+10, file=filename(1,l))
            open(unit=l+20, file=filename(2,l))
            do i=0,nbin(l)-1
                read(l+10,*) z, g2
                read(l+20,*) z, g3
                if (i.eq.bin0(l)) then
                    g2s(l)=g2
                    g3s(l,:)=g3
                end if
            end do
            close(l+10)
            close(l+20)
            lam2(l,:)=g3s(l,:)/g2s(l)/g2s(l)-1.0d0 ! lambda in TPT2
            Gam2LA(l,:)=lam2(l,:)+1.0d0 ! linear approximation (LA)
            Gam2SA(l,:)=g3s(l,:)/g2s(l)/g2s(l)/g2s(l) ! superposition approximation (SA)
            write(51,*) d(l), g2s(l)
            write(52,*) d(l), g3s(l,:)
            write(53,*) d(l), Gam2LA(l,:)
            write(54,*) d(l), Gam2SA(l,:)
            write(55,*) d(l), lam2(l,:)
        end do
        close(51)
        close(52)
        close(53)
        close(54)
        close(55)
        return
      end
C *************************************************************************************
