c X+L few-level autoionizing system with b-b, b-f, f-f transitions
c Electron spectra for partial wave I & II
      MODULE Global
      IMPLICIT NONE
      REAL,PARAMETER::
     & cFT=0.39894,small=1.0E-6,large=1.0E6,clight=137.0,
     & aeV=27.211384,afs=0.024188843,In0=3.51E16,lam=45.563353,
     & amm=5.2917721E-8,acm3=6.7483346E24,torr=9.6564746E18
      COMPLEX,PARAMETER::i=(0.0,1.0)
      CHARACTER(LEN=48)::title
      INTEGER::
     & it,kt,Nkt,Nt,iE,jE,NE,m,Nm,n,Nn,ibk,Nbk,iEwx
      REAL::
     & pi,Emin,Emax,dE, ! dimension
     & EL0,tL,wL0,CEPL,t0, ! L
     & EX0,tX,wX0,CEPX, ! X
     & jg2,jg3,Dg1,D12,D13,Ip ! atomic
      INTEGER,ALLOCATABLE::nh(:)
      REAL,ALLOCATABLE::
     & Em(:),Vm(:),Km(:),En(:),Vn(:),Kn(:),E(:),jgn(:),
     & Dgm(:),Dmn(:,:),Dm2(:),D1n(:),
     & t(:),dt(:),b(:),EL(:),EL2(:),
     & tmE(:,:),cth1(:),sth1(:),tnE(:,:),cth2(:),sth2(:)
      COMPLEX,ALLOCATABLE::
     & FX(:),cm(:,:),cn(:,:),c1(:,:),c2(:,:),c3(:,:),cE1(:,:),cE2(:,:),
     & Hgm(:),Hmn(:,:),H1n(:),Hm2(:)
      CONTAINS

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE Input
      REAL::InX,lamL,InL
      REAL,ALLOCATABLE::bfs(:),Gm(:),Gn(:)
      NAMELIST/ATM/Dg1,D12,D13,Ip
      NAMELIST/XUV/wX0,tX,InX,CEPX
      NAMELIST/LSR/lamL,tL,InL,CEPL,t0
      NAMELIST/DMS/Nbk,Emin,Emax,NE,Nkt
      OPEN(11,FILE='xle.cfg',STATUS='OLD')
      READ(11,'(A48)') title
      READ(11,NML=ATM)
      READ(11,*) Nm
      ALLOCATE(Em(Nm),Gm(Nm),Km(Nm),Vm(Nm),Dgm(Nm),Dm2(Nm),
     & Hgm(Nm),Hm2(Nm))
      DO m=1,Nm
        READ(11,*) Em(m),Gm(m),Dgm(m),Dm2(m)
      END DO
      READ(11,*) Nn
      ALLOCATE(En(Nn),Gn(Nn),Kn(Nn),Vn(Nn),D1n(Nn),H1n(Nn),jgn(Nn))
      DO n=1,Nn
        READ(11,*) En(n),Gn(n),D1n(n)
      END DO
      ALLOCATE(Dmn(Nm,Nn),Hmn(Nm,Nn))
      DO m=1,Nm
        READ(11,*) (Dmn(m,n),n=1,Nn)
      END DO
      READ(11,NML=XUV)
      READ(11,NML=LSR)
      READ(11,NML=DMS)
      ALLOCATE(b(Nbk+1),bfs(Nbk+1),nh(Nbk)) ! b[fs] = block boundaries, nh = block number
      READ(11,*) (bfs(ibk),nh(ibk),ibk=1,Nbk),bfs(Nbk+1)
      CLOSE(11)
c atomic
      Ip=Ip/aeV
      Em=Em/aeV
      Km=Gm/2.0/aeV
      Vm=SQRT(Km/pi)
      En=En/aeV
      Kn=Gn/2.0/aeV
      Vn=SQRT(Kn/pi)
c field
      wX0=wX0/aeV
      wL0=lam/lamL
      tX=tX/afs/1.144
      tL=tL/afs/1.144
      t0=t0/afs
      EX0=SQRT(InX/In0)
      EL0=SQRT(InL/In0)
c dimension
      b=bfs/afs
      Emin=Emin/aeV
      Emax=Emax/aeV
c composite atomic parameters
      jg2=pi*Dg1*D12
      jg3=pi*Dg1*D13
      Hgm=Dgm-i*pi*Vm*Dg1
      H1n=D1n-i*pi*D12*Vn
      Hm2=Dm2-i*pi*Vm*D12 
      jgn=pi*Dg1*D1n
      DO m=1,Nm
        DO n=1,Nn
          Hmn(m,n)=Dmn(m,n)-i*pi*(Vm(m)*D1n(n)+Dm2(m)*Vn(n))
     &            -pi*pi*Vm(m)*D12*Vn(n)
        END DO
      END DO
 5    FORMAT(20F14.6)
 8    FORMAT(20A14)
      WRITE(150,8) 'Em(eV)','Gm(eV)'
      DO m=1,Nm
        WRITE(150,5) Em(m)*aeV,Gm(m)*aeV
      END DO
      WRITE(150,8) 'En(eV)','Gn(eV)'
      DO n=1,Nn
        WRITE(150,5) En(n)*aeV,Gn(n)*aeV
      END DO
      WRITE(150,8) 'qm(eV)','qmn(eV)'
      DO m=1,Nm
        WRITE(150,5) Dgm(m)/pi/Vm(m)/Dg1,
     &               ((Dmn(m,n)-pi*pi*Vm(m)*D12*Vn(n))/
     &                pi/(Vm(m)*D1n(n)+Dm2(m)*Vn(n)),n=1,Nn)
      END DO
      RETURN
      END SUBROUTINE Input

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c initial values
      SUBROUTINE Dimen
      INTEGER::ih,mh
      REAL::hh,xx,theta
c time
      Nt=SUM(nh)
      IF (Nkt==0) THEN
        kt=1
      ELSE
        kt=INT(Nt/Nkt)
      END IF
      Nt=Nt+1
      ALLOCATE(t(Nt),dt(Nt),FX(Nt),EL(Nt),EL2(Nt),
     & cm(Nm,Nt),cn(Nn,Nt))
      mh=0
      DO ibk=1,Nbk
        hh=(b(ibk+1)-b(ibk))/REAL(nh(ibk))
        DO ih=1,nh(ibk)
          t(mh+ih)=b(ibk)+hh*REAL(ih-1)
          dt(mh+ih)=hh
        END DO
        mh=mh+nh(ibk)
      END DO
      t(Nt)=t(Nt-1)+dt(Nt-1)
c continuum E
      dE=(Emax-Emin)/REAL(NE)
      NE=NE+1
      ALLOCATE(E(NE),c1(NE,Nt),c2(NE,Nt),c3(NE,Nt),
     & cE1(NE,Nt),cE2(NE,Nt))
      DO iE=1,NE
        E(iE)=Emin+dE*REAL(iE-1)
      END DO
      iEwx=INT((wX0-Emin)/dE)+1
c eigenstate phase theta
      ALLOCATE(tmE(Nm,NE),cth1(NE),sth1(NE),
     & tnE(Nn,NE),cth2(NE),sth2(NE))
      DO iE=1,NE
        xx=0.0
        DO m=1,Nm
          tmE(m,iE)=-Km(m)/(E(iE)-Em(m))
          xx=xx+tmE(m,iE)
        END DO
        theta=ATAN(xx)
        cth1(iE)=COS(theta)
        sth1(iE)=SIN(theta)
        xx=0.0
        DO n=1,Nn
          tnE(n,iE)=-Kn(n)/(E(iE)-En(n))
          xx=xx+tnE(n,iE)
        END DO
        theta=ATAN(xx)
        cth2(iE)=COS(theta)
        sth2(iE)=SIN(theta)
      END DO
      RETURN
      END SUBROUTINE Dimen

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE FieldX
      COMPLEX::a,ex1,ex2
      DO it=1,nt
        FX(it)=ES2(tX,EX0/2.0,t(it))*EXP(i*CEPX)
      END DO
      RETURN
      END SUBROUTINE FieldX

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE FieldL ! incident fields
      REAL::ts,ra
      DO it=1,nt
        ts=t(it)-t0
        EL(it)=ES2(tL,EL0,ts)*COS(wL0*ts+CEPL)
        EL2(it)=EL(it)*EL(it)
      END DO
      CALL FieldLw
      RETURN
      END SUBROUTINE FieldL

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE FieldLw
      COMPLEX::cw,ex1,ex2
      OPEN(26,FILE='Fw.dat',STATUS='UNKNOWN')
      DO iE=1,NE
        cw=(0.0,0.0)
        DO it=1,Nt-1
          ex1=EXP(-i*E(iE)*t(it))
          ex2=EXP(-i*E(iE)*t(it+1))
          cw=cw+dt(it)*(ex1*EL(it)+ex2*EL(it+1))/2.0
        END DO
        cw=cFT*cw
        WRITE(26,10) E(iE)*aeV,ABS(cw)**2.0/aeV
      END DO
      CLOSE(26)
 10   FORMAT(20ES14.6)
      RETURN
      END SUBROUTINE FieldLw

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c General sin^2 type pulse shape of
c (tau=FWHM/1.144, amplitude, time relative to peak)
      FUNCTION ES2(tau,amp,tt)
      REAL::ES2,tau,amp,tt
      IF (ABS(tt)<=pi/2.0*tau) THEN
        ES2=amp*COS(tt/tau)**2
      ELSE
        ES2=0.0
      END IF
      RETURN
      END FUNCTION ES2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE System
      CALL Dimen
      CALL FieldX
      CALL FieldL
      RETURN
      END SUBROUTINE System

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE Micro
      COMPLEX::x0,x1,x2,x3,x4,x5,xx,y1,y2,y3,y4,yy,der1,der2,der3
c bound states
      cm=(0.0,0.0)
      cn=(0.0,0.0)
      DO it=1,Nt-1
c cm
        DO m=1,Nm
          x0=i*(wX0-Em(m))-Km(m)
          x1=i*CONJG(FX(it))*Hgm(m)
          x2=(0.0,0.0)
          DO n=1,Nn
            x2=x2+Hmn(m,n)*cn(n,it)
          END DO
          xx=x1+i*EL(it)*x2
          y1=x0*cm(m,it)+xx
          y2=x0*(cm(m,it)+0.5*dt(it)*y1)+xx
          y3=x0*(cm(m,it)+0.5*dt(it)*y2)+xx
          y4=x0*(cm(m,it)+dt(it)*y3)+xx
          yy=(y1+2.0*y2+2.0*y3+y4)/6.0
          cm(m,it+1)=cm(m,it)+dt(it)*yy
        END DO
c cn
        DO n=1,Nn
          x0=i*(wX0-En(n))-Kn(n)
          x1=-CONJG(FX(it))*EL(it)*jgn(n)
          x2=(0.0,0.0)
          DO m=1,Nm
            x2=x2+Hmn(m,n)*cm(m,it)
          END DO
          xx=x1+i*EL(it)*x2
          y1=x0*cn(n,it)+xx
          y2=x0*(cn(n,it)+0.5*dt(it)*y1)+xx
          y3=x0*(cn(n,it)+0.5*dt(it)*y2)+xx
          y4=x0*(cn(n,it)+dt(it)*y3)+xx
          yy=(y1+2.0*y2+2.0*y3+y4)/6.0
          cn(n,it+1)=cn(n,it)+dt(it)*yy
        END DO
      END DO
c continuum states
      c1=(0.0,0.0)
      c2=(0.0,0.0)
      c3=(0.0,0.0)
      der1=(0.0,0.0)
      der2=(0.0,0.0)
      der3=(0.0,0.0)
      DO it=1,Nt-1
        DO iE=1,NE
          x0=i*(wX0-E(iE))
c c1
          x1=i*CONJG(FX(it))*Dg1
          x2=(0.0,0.0)
          DO m=1,Nm
            x2=x2+Vm(m)*cm(m,it)
          END DO
          x3=(0.0,0.0)
          DO n=1,Nn
            x3=x3+H1n(n)*cn(n,it)
          END DO
          x3=x3+pi*(D12*der2+D13*der3)
          xx=x1-i*x2+i*EL(it)*x3
          y1=x0*c1(iE,it)+xx
          y2=x0*(c1(iE,it)+0.5*dt(it)*y1)+xx
          y3=x0*(c1(iE,it)+0.5*dt(it)*y2)+xx
          y4=x0*(c1(iE,it)+dt(it)*y3)+xx
          yy=(y1+2.0*y2+2.0*y3+y4)/6.0
          c1(iE,it+1)=c1(iE,it)+dt(it)*yy
c c2
          x1=-EL(it)*CONJG(FX(it))*jg2
          x2=(0.0,0.0)
          DO n=1,Nn
            x2=x2+Vn(n)*cn(n,it)
          END DO
          x3=(0.0,0.0)
          DO m=1,Nm
            x3=x3+Hm2(m)*cm(m,it)
          END DO
          x3=x3+pi*D12*der1
          xx=x1-i*x2+i*EL(it)*x3
          y1=x0*c2(iE,it)+xx
          y2=x0*(c2(iE,it)+0.5*dt(it)*y1)+xx
          y3=x0*(c2(iE,it)+0.5*dt(it)*y2)+xx
          y4=x0*(c2(iE,it)+dt(it)*y3)+xx
          yy=(y1+2.0*y2+2.0*y3+y4)/6.0
          c2(iE,it+1)=c2(iE,it)+dt(it)*yy
c c3
          x1=-EL(it)*CONJG(FX(it))*jg3
          x3=pi*D13*der1
          xx=x1+i*EL(it)*x3
          y1=x0*c3(iE,it)+xx
          y2=x0*(c3(iE,it)+0.5*dt(it)*y1)+xx
          y3=x0*(c3(iE,it)+0.5*dt(it)*y2)+xx
          y4=x0*(c3(iE,it)+dt(it)*y3)+xx
          yy=(y1+2.0*y2+2.0*y3+y4)/6.0
          c3(iE,it+1)=c3(iE,it)+dt(it)*yy
        END DO
        der1=(c1(iEwx,it+1)-c1(iEwx,it))/dt(it)
        der2=(c2(iEwx,it+1)-c2(iEwx,it))/dt(it)
        der2=(c3(iEwx,it+1)-c3(iEwx,it))/dt(it)
      END DO
c eigenstates
      cE1=(0.0,0.0)
      cE2=(0.0,0.0)
      DO it=1,Nt
        DO iE=1,NE
          xx=(0.0,0.0)
          DO m=1,Nm
            xx=xx+tmE(m,iE)/pi/Vm(m)*cm(m,it)
          END DO
          cE1(iE,it)=cth1(iE)*xx-(cth1(iE)-i*sth1(iE))*c1(iE,it)
          DO m=1,Nm
            IF (E(iE)==Em(m)) THEN
              cE1(iE,it)=-cm(m,it)/pi/Vm(m)-i*c1(iE,it)
            END IF
          END DO
          xx=(0.0,0.0)
          DO n=1,Nn
            xx=xx+tnE(n,iE)/pi/Vn(n)*cn(n,it)
          END DO
          cE2(iE,it)=cth2(iE)*xx-(cth2(iE)-i*sth2(iE))*c2(iE,it)
          DO n=1,Nn
            IF (E(iE)==En(n)) THEN
              cE2(iE,it)=cn(n,it)/pi/Vn(n)+i*c2(iE,it)
            END IF
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE Micro

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE Output
      REAL,ALLOCATABLE::
     & pm(:,:),pn(:,:),p1(:,:),p2(:,:),p3(:,:),pE1(:,:),pE2(:,:),
     & pma(:),pna(:),p1a(:),p2a(:),p3a(:),pE1a(:),pE2a(:)
      ALLOCATE(pm(Nm,Nt),pn(Nn,Nt),p1(NE,Nt),p2(NE,Nt),p3(NE,Nt),
     & pma(Nt),pna(Nt),p1a(Nt),p2a(Nt),p3a(Nt),pE1(NE,Nt),pE2(NE,Nt),
     & pE1a(Nt),pE2a(Nt))
      pm=ABS(cm)**2
      pn=ABS(cn)**2
      p1=ABS(c1)**2
      p2=ABS(c2)**2
      p3=ABS(c3)**2
      pE1=ABS(cE1)**2
      pE2=ABS(cE2)**2
      pma=SUM(pm,1)
      pna=SUM(pn,1)
      p1a=SUM(p1,1)*dE
      p2a=SUM(p2,1)*dE
      p3a=SUM(p3,1)*dE
      pE1a=SUM(pE1,1)*dE
      pE2a=SUM(pE2,1)*dE
      DEALLOCATE(cm,cn,c1,c2,c3)
      OPEN(22,FILE='F.dat',STATUS='UNKNOWN')
      OPEN(30,FILE='pall.dat',STATUS='UNKNOWN')
      OPEN(32,FILE='pbound.dat',STATUS='UNKNOWN')
      OPEN(36,FILE='pcont.dat',STATUS='UNKNOWN')
      OPEN(38,FILE='peigen.dat',STATUS='UNKNOWN')
      WRITE(30,8) 't(fs)','m','n','1','2','3','E1','E2'
      DO it=1,Nt,kt
        WRITE(22,10) t(it)*afs,2.0*ABS(FX(it)),EL(it)
        WRITE(30,10) t(it)*afs,pma(it),pna(it),
     &               p1a(it),p2a(it),p3a(it),pE1a(it),pE2a(it)
        WRITE(32,10) t(it)*afs,(pm(m,it),m=1,Nm),(pn(n,it),n=1,Nn)
        DO iE=1,NE
          WRITE(36,10) t(it)*afs,E(iE)-Ip,p1(iE,it),p2(iE,it),p3(iE,it)
          WRITE(38,10) t(it)*afs,E(iE)-Ip,pE1(iE,it),pE2(iE,it)
        END DO
        WRITE(36,*)
        WRITE(38,*)
      END DO
      CLOSE(22)
      CLOSE(30)
      CLOSE(32)
      CLOSE(34)
      CLOSE(36)
      CLOSE(38)
      OPEN(46,FILE='spec.dat',STATUS='UNKNOWN')
      DO iE=1,NE
        WRITE(46,10) (E(iE)-Ip)*aeV,pE1(iE,nt),pE2(iE,nt),p3(iE,nt)
      END DO
      CLOSE(46)
 8    FORMAT(20A14)
 10   FORMAT(20ES14.6)
      RETURN
      END SUBROUTINE Output


      END MODULE Global

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      PROGRAM Main
      USE Global
      IMPLICIT NONE
      pi=ACOS(-1.0)
      CALL Input
      CALL System
      CALL Micro
      CALL Output
      STOP
      END PROGRAM Main

