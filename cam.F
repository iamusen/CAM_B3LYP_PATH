      MODULE cam 
        USE prec 
       !INTEGER, PARAMETER :: q     = SELECTED_REAL_KIND(10)
        REAL(q), SAVE      :: MU    = 0._q
        REAL(q), SAVE      :: CAMFA = 0._q
        REAL(q), SAVE      :: CAMFB = 0._q
      CONTAINS
      SUBROUTINE NSPCAM(R,D,F,V1,VD1,CU,CA,CB,CLYP,CVWN)
        IMPLICIT REAL(q) (A-Z)
        REAL(q), INTENT(IN)  :: R,D,CU,CA,CB,CLYP,CVWN
        REAL(q), INTENT(OUT) :: F,V1,VD1
       !Local variables
        REAL(q) R1,D1
        REAL(q) FCAM,FVWN,FLYP
        REAL(q) VCAM1,VCAM2,DCAM1,DCAM2
        REAL(q) VVWN1,VLYP1,DLYP1
        REAL(q) VVWN2,VLYP2,DLYP2
        R1    = R/2._q
        D1    = D/2._q
        CAMFA = CA
        CAMFB = CB
        MU    = CU
        CALL NSPB88(R1,R1,D1,D1,FCAM,VCAM1,VCAM2,DCAM1,DCAM2)
        CALL   LYP(R1,R1,D1,D1,FLYP,VLYP1,VLYP2,DLYP1,DLYP2,CLYP)
        CALL   VWN(R1,R1,FVWN,VVWN1,VVWN2,CVWN)
        F     = FCAM  + FLYP  + FVWN
        V1    = VCAM1 + VLYP1 + VVWN1
        VD1   = DCAM1 + DLYP1 
        RETURN
      END SUBROUTINE NSPCAM
      SUBROUTINE SPNCAM(R1,R2,D1,D2,F,V1,V2,VD1,VD2,CU,CA,CB,CLYP,CVWN)
        IMPLICIT REAL(q) (A-Z)
        REAL(q), INTENT(IN)  :: R1,R2,D1,D2,CU,CA,CB,CLYP,CVWN
        REAL(q), INTENT(OUT) :: F,V1,V2,VD1,VD2
       !Local variables
        REAL(q) FCAM,FVWN,FLYP
        REAL(q) VCAM1,VCAM2,DCAM1,DCAM2
        REAL(q) VVWN1,VLYP1,DLYP1
        REAL(q) VVWN2,VLYP2,DLYP2
        CAMFA = CA
        CAMFB = CB
        MU    = CU
        CALL   VWN(R1,R2,FVWN,VVWN1,VVWN2,CVWN)
        CALL   LYP(R1,R2,D1,D2,FLYP,VLYP1,VLYP2,DLYP1,DLYP2,CLYP)
        CALL SPNB88(R1,R2,D1,D2,FCAM,VCAM1,VCAM2,DCAM1,DCAM2)
        F     = FCAM  + FLYP  + FVWN
        V1    = VCAM1 + VLYP1 + VVWN1
        VD1   = DCAM1 + DLYP1
        V2    = VCAM2 + VLYP2 + VVWN2
        VD2   = DCAM2 + DLYP2 
        RETURN
      END SUBROUTINE SPNCAM
      SUBROUTINE SPNB88(R1,R2,D1,D2,F,V1,V2,VV1,VV2)
        IMPLICIT REAL(q) (A-Z)
        REAL(q), INTENT(IN)  :: R1,R2,D1,D2
        REAL(q), INTENT(OUT) :: F,V1,V2,VV1,VV2
        REAL(q) F1,F2
        F1  = 0._q
        F2  = 0._q
        V1  = 0._q
        V2  = 0._q
        VV1 = 0._q
        VV2 = 0._q
        CALL CAM_B88(R1,D1,F1,V1,VV1)
        CALL CAM_B88(R2,D2,F2,V2,VV2)
        F   = F1 + F2
        RETURN
      END SUBROUTINE SPNB88
      SUBROUTINE NSPB88(R1,R2,D1,D2,F,V1,V2,VV1,VV2)
        IMPLICIT REAL(q) (A-Z)
        REAL(q), INTENT(IN)  :: R1,R2,D1,D2
        REAL(q), INTENT(OUT) :: F,V1,V2,VV1,VV2
        REAL(q) F1
        V2  = 0._q
        VV2 = 0._q
        CALL CAM_B88(R1,D1,F1,V1,VV1)
        F   = F1 + F1
        RETURN
      END SUBROUTINE NSPB88
      SUBROUTINE CAM_B88(R,DR,F,V,VV)
        IMPLICIT REAL(q) (A-Z)
        REAL(q), INTENT(IN)  :: R,DR
        REAL(q), INTENT(OUT) :: F,V,VV
        REAL(q) Ep,Ef,DKN,DKG
        REAL(q) KG,K,H,HR,HD,SQ,DAR
        REAL(q) A,B,C,AA,AAA,AAAA,DEN,ARC
        REAL(q) R13,R43,X,XX,PIS,DENH,DEND
        REAL(q) DSQ,DKX,DXN,DXG,TRM
        REAL(q), PARAMETER :: F1o2 = 1._q/2._q 
        REAL(q), PARAMETER :: F1o3 = 1._q/3._q 
        REAL(q), PARAMETER :: F1o4 = 1._q/4._q 
        REAL(q), PARAMETER :: F3o4 = 3._q/4._q 
        REAL(q), PARAMETER :: F1o6 = 1._q/6._q 
        REAL(q), PARAMETER :: F3o2 = 3._q/2._q
        REAL(q), PARAMETER :: F2o3 = 2._q/3._q 
        REAL(q), PARAMETER :: F4o3 = 4._q/3._q 
        REAL(q), PARAMETER :: F8o3 = 8._q/3._q 
        REAL(q), PARAMETER :: BETA = 0.0042_q
        REAL(q), PARAMETER :: BET2 = 0.0042_q*0.0042_q
        REAL(q), PARAMETER :: B6   = 0.0252_q
        REAL(q), PARAMETER :: PI   = 4._q*ATAN(1._q) 
        REAL(q), PARAMETER :: CX   = F3o2*((F3o4/PI)**F1o3)
       !REAL(q), PARAMETER :: MU   = 0.33_q 
       !REAL(q), PARAMETER :: CAMFA= 0.19_q
       !REAL(q), PARAMETER :: CAMFB= 0.46_q
        PIS = SQRT(PI)
        R13 = R**F1o3
        R23 = R13*R13
        R43 = R13*R
        R63 = R23*R43
        R53 = R13*R43
        R73 = R*R43
        ADR = ABS(DR)
        X   = ADR/R43 
        XX  = X*X
        XXX = XX*X
        SQ  = SQRT(XX+1._q)
        ARC = LOG(X+SQ)
        DEN = 1._q + B6*X*ARC 
        KG  = 2._q*BETA*XX/DEN
        K   = KG+2._q*CX
        F   = -F1o2*K*R43
        DKX = KG*(2._q/X-(B6*(X/SQ+ARC)/DEN))
        DXR = -F4o3*X/R
        DXG = 1._q/R43
        DKR = DKX*DXR
        DKG = DKX*DXG
        V   =-F1o2*(F4o3*K*R13+R43*DKR) 
        VV  =-F1o2*R43*DKG
        VV  = SIGN(1._q,DR)*VV
        SQK = SQRT(K)
        A   = SQK*MU*F1o6/PIS/R13
        AA  = A*A
        AAA = AA*A
        AAAA= AA*AA
        Ep  = Exp(-F1o4/AA)
        Ef  = Erf(F1o2/A)
        Efp = PIS*Ef
        B   = Ep - 1._q
        C   = 2._q*AA*B + F1o2
        H   = 1._q - CAMFA - CAMFB*F8o3*A*(Efp+2._q*A*(B-C))
        IF (H.LT.0.35_q.OR.Ef.LT.0.1_q) THEN
          H = 0.35_q
          F = F *H
          V = V *H
          VV= VV*H
          RETURN
        ELSE
          F   = H*F
          V   = V*H 
          VV  = VV*H
          Efa = -1._q*Ep/AA
          DBA = F1o2*Ep/AAA
          DCA = 4._q*A*B+2._q*AA*DBA
          DHA = -F8o3*((Efp +2._q*A*(B-C)) + &
                     A*(Efa+(2._q*(B-C)+2._q*A*(DBA-DCA))))
          DAR = MU/(6._q*PIS)*(F1o2*DKR/SQK/R13-F1o3*SQK/R43)
          DAG = MU/(6._q*PIS*R13)*(F1o2*DKG/SQK)
          DHR = DHA*DAR*CAMFB
          DHG = DHA*DAG*CAMFB
          V   = V - DHR*K*R43*F1o2
          VV  = VV- DHG*K*R43*F1o2
        ENDIF
        RETURN 
      END SUBROUTINE CAM_B88
      SUBROUTINE VWN(r1,r2,energy,va,vb,cvwn)
        IMPLICIT NONE
        REAL(q),  INTENT(in)  :: r1, r2, cvwn
        REAL(q),  INTENT(out) :: energy, va, vb 
        REAL(q),  PARAMETER   :: facx=0.787623317899743D00 ! (3/(4*pi))**(1/6)
        REAL(q),  PARAMETER   :: fach=1.709920934161366D00 ! 4/(9*(2**(1/3)-1))
        REAL(q),  PARAMETER   :: prmqp=6.15199081975908D00 ! Qp=sqrt(4*cp-bp*bp)
        REAL(q),  PARAMETER   :: prmqf=4.73092690956011D00 ! Qf=sqrt(4*cf-bf*bf)
        REAL(q),  PARAMETER   :: prmqa=7.12310891781812D00 ! Qa=sqrt(4*ca-ba*ba)
        REAL(q),  PARAMETER   :: prmqpinv=0.1625490071910D00 ! 1/Qp
        REAL(q),  PARAMETER   :: prmqfinv=0.211375068589462D0 ! 1/Qf
        REAL(q),  PARAMETER   :: prmqainv=0.14038813831676D0 ! 1/Qa
        REAL(q),  PARAMETER   :: prmxx0pinv=7.9650086660586D-2 ! 1/Xx0p=1/(x0p*x0p+bp*x0p+cp)
        REAL(q),  PARAMETER   :: prmxx0finv=6.3016782913201D-2 ! 1/Xx0f=1/(x0f*x0f+bf*x0f+cf)
        REAL(q),  PARAMETER   :: prmxx0ainv=7.6928162709702D-2 ! 1/Xx0a=1/(x0a*x0a+ba*x0a+ca)
        REAL(q),  PARAMETER   :: prmap=3.10907D-02, prmbp=3.72744D0, prmcp=1.29352D1
        REAL(q),  PARAMETER   :: prmaf=1.554535D-02, prmbf=7.06042D0, prmcf=1.80578D1
        REAL(q),  PARAMETER   :: prmaa=-1.688686D-02, prmba=1.13107D0, prmca=1.30045D1
        REAL(q),  PARAMETER   :: prmx0p=-0.10498D+00, prmx0f=-0.325D0, prmx0a=-0.0047584D0
        REAL(q),  PARAMETER   :: prmp2=9.690227711544374D-4  ! -ap*bp*x0p/Xx0p
        REAL(q),  PARAMETER   :: prmf2=2.247867095542611D-3  ! -af*bf*x0f/Xx0f
        REAL(q),  PARAMETER   :: prma2=-6.991730719309994D-6 ! -aa*ba*x0a/Xx0a
        REAL(q),  PARAMETER   :: prmp3=3.878329487811301D-2  ! ap*2*bp/Qp*(cp-x0p*x0p)/Xx0p
        REAL(q),  PARAMETER   :: prmf3=5.249139316978093D-2  ! af*2*bf/Qf*(cf-x0f*x0f)/Xx0f
        REAL(q),  PARAMETER   :: prma3=-5.365090596990275D-3 ! aa*2*ba/Qa*(ca-x0a*x0a)/Xx0a
        REAL(q),  PARAMETER   :: one=1.0D+00, two=2.D0, onethird=0.3333333333333333D0
        REAL(q),  PARAMETER   :: onesixth=0.1666666666666666D0, four=4.D0
        REAL(q),  PARAMETER   :: fourthird=1.333333333333333D0, nineeighth=1.125D0
        REAL(q),  PARAMETER   :: onehalf=1.5D0
        REAL(q) :: rhoa, rhob, rhot, x, xxp, xxf, xxa, xxpinv, xxfinv, xxainv, tmp1p, tmp2p, tmp3p
        REAL(q) :: tmp1f, tmp2f, tmp3f, tmp1a, tmp2a, tmp3a, epsp, epsf, epsa, xbq2p, xbq2f, xbq2a
        REAL(q) :: depsp, depsf, depsa, xbp, xbf, xba, bracketp, bracketf, bracketa, hx, dhx
        REAL(q) :: zeta, zeta2, zeta3, zeta4, gzeta, dgzeta, vwnpot, dvwnpot1, dvwnpot2
        rhoa= MAX(1.E-10_q,r1)
        rhob= MAX(1.E-10_q,r2)
        rhot= rhoa+rhob
        x= facx*(rhot**(-onesixth))
        xxp= x*x+prmbp*x+prmcp
        xxf= x*x+prmbf*x+prmcf
        xxa= x*x+prmba*x+prmca
        xxpinv=one/xxp
        xxfinv=one/xxf
        xxainv=one/xxa
        xbp= two*x+prmbp
        xbf= two*x+prmbf
        xba= two*x+prmba
        tmp1p= log(x*x*xxpinv)
        tmp1f= log(x*x*xxfinv)
        tmp1a= log(x*x*xxainv)
        tmp2p= log((x-prmx0p)*(x-prmx0p)*xxpinv)
        tmp2f= log((x-prmx0f)*(x-prmx0f)*xxfinv)
        tmp2a= log((x-prmx0a)*(x-prmx0a)*xxainv)
        tmp3p= atan(prmqp/xbp)
        tmp3f= atan(prmqf/xbf)
        tmp3a= atan(prmqa/xba)
        epsp= prmap*tmp1p+prmp2*tmp2p+prmp3*tmp3p
        epsf= prmaf*tmp1f+prmf2*tmp2f+prmf3*tmp3f
        epsa= prmaa*tmp1a+prma2*tmp2a+prma3*tmp3a
        xbq2p= one/(xbp*xbp+prmqp*prmqp)
        xbq2f= one/(xbf*xbf+prmqf*prmqf)
        xbq2a= one/(xba*xba+prmqa*prmqa)
        bracketp= two/(x-prmx0p)-xbp*xxpinv-four*(two*prmx0p+prmbp)*xbq2p
        bracketf= two/(x-prmx0f)-xbf*xxfinv-four*(two*prmx0f+prmbf)*xbq2f
        bracketa= two/(x-prmx0a)-xba*xxainv-four*(two*prmx0a+prmba)*xbq2a
        depsp=prmap*(two/x-xbp*xxpinv-four*prmbp*xbq2p-prmbp*prmx0p*prmxx0pinv*bracketp)
        depsf=prmaf*(two/x-xbf*xxfinv-four*prmbf*xbq2f-prmbf*prmx0f*prmxx0finv*bracketf)
        depsa=prmaa*(two/x-xba*xxainv-four*prmba*xbq2a-prmba*prmx0a*prmxx0ainv*bracketa)
        hx= fach*(epsf-epsp)/epsa-one
        dhx= fach*(depsf-depsp-(epsf-epsp)*depsa/epsa)/epsa
        zeta=(rhoa-rhob)/rhot
        zeta2= zeta *zeta
        zeta3= zeta2*zeta
        zeta4= zeta2*zeta2
        gzeta= nineeighth*((one+zeta)**fourthird+(one-zeta)**fourthird-two)
        dgzeta= onehalf*((one+zeta)**onethird-(one-zeta)**onethird)
        vwnpot= epsp+epsa*gzeta*(one+hx*zeta4)
        energy= cvwn*rhot*vwnpot
        dvwnpot1=-x*onesixth*(depsp+depsa*gzeta*(one+hx*zeta4)+epsa*gzeta*dhx*zeta4)
        dvwnpot2= epsa*(dgzeta*(one+hx*zeta4)+four*gzeta*hx*zeta3)
        va= cvwn*(vwnpot+dvwnpot1+dvwnpot2*(one-zeta))
        vb= cvwn*(vwnpot+dvwnpot1-dvwnpot2*(one+zeta))
        return
      END SUBROUTINE VWN 
      SUBROUTINE LYP(r1,r2,grhoa,grhob,energy,va,vb,vda,vdb,clyp)
        IMPLICIT NONE
        REAL(q),  INTENT(in)  :: r1, r2, grhoa, grhob, clyp
        REAL(q),  INTENT(out) :: va,vda,vb,vdb,energy
        REAL(q),  PARAMETER   :: face=3.646239897876478D+01 ! 2**(11/3)*3/10*(3*pi*pi)**(2/3)
        REAL(q),  PARAMETER   :: prma=0.04918D+00, prmb=0.132D+00, prmc=0.2533D+00, prmd=0.349D+00
        REAL(q),  PARAMETER   :: half=0.5D+00, one=1.0D+00, two=2.0D+00, three=3.0D+00, four=4.0D+00
        REAL(q),  PARAMETER   :: five=5.0D+00, seven=7.0D+00, eleven=11.0D+00, p47=47.0D+00
        REAL(q),  PARAMETER   :: onethird=0.3333333333333333D+00, oneninth=0.1111111111111111D+00
        REAL(q),  PARAMETER   :: fourthird=1.333333333333333D+00, eightthird=2.6666666666666666D+00
        REAL(q),  PARAMETER   :: p11third=3.666666666666666D+00
        REAL(q) :: rhoa,rhob,rhot, rhoab, rhot13, rhotinv, rhot13inv, rhot43inv, rhot53inv, rhoa83, rhob83
        REAL(q) :: rhoatinv, rhobtinv, denom, omega, delta, delta11a, delta11b, omegap, deltap
        REAL(q) :: gradaa, gradab, gradbb, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
        REAL(q) :: tmp10, tmp11,rhoa13,rhob13,tmp12
        REAL(q) :: dlypaa, dlypab, dlypbb, dlypaaa, dlypaab, dlypabb, dlypbaa, dlypbab, dlypbbb
        rhoa       = MAX(r1,0.1E-10_q)
        rhob       = MAX(r2,0.1E-10_q)
        rhot       = rhoa+rhob
        rhoab      = rhoa*rhob
        rhot13     = rhot**onethird
        rhoa13     = rhoa**onethird
        rhob13     = rhob**onethird
        rhotinv    = one/rhot
        rhot13inv  = rhot13*rhot13*rhotinv
        rhot43inv  = rhotinv*rhot13inv
        rhot53inv  = rhotinv*rhot13inv*rhot13inv
        rhoa83     = rhoa*rhoa*rhoa13*rhoa13
        rhob83     = rhob*rhob*rhob13*rhob13
        rhoatinv   = rhoa*rhotinv
        rhobtinv   = rhob*rhotinv
        denom      = one/(one+prmd*rhot13inv)
        omega      = exp(-prmc*rhot13inv)*denom*rhot53inv*rhot53inv*rhot13inv*prma*prmb
        delta      = prmc*rhot13inv+prmd*rhot13inv*denom
        delta11a   = (delta-eleven)*rhoatinv
        delta11b   = (delta-eleven)*rhobtinv
        gradaa     = grhoa*grhoa   
        gradab     = grhoa*grhob   
        gradbb     = grhob*grhob   
        tmp1       = oneninth*rhoab
        tmp2       = four*prma*denom*rhob*rhoatinv
        dlypaa     =-omega*(tmp1*(one-three*delta-delta11a)-rhob*rhob)
        dlypab     =-omega*(tmp1*(p47-seven*delta)-fourthird*rhot*rhot)
        dlypbb     =-omega*(tmp1*(one-three*delta-delta11b)-rhoa*rhoa)
        tmp3       = dlypaa*gradaa+dlypab*gradab+dlypbb*gradbb
        energy     = clyp*(-tmp2-face*omega*rhoab*(rhoa83+rhob83)+tmp3)
        omegap     =-onethird*rhot43inv*(eleven*rhot13-prmc-prmd*denom)
        deltap     = onethird*(prmd*prmd*rhot53inv*denom*denom-rhotinv*delta)*oneninth*rhoab
        tmp4       = oneninth*rhob*(one-three*delta-delta11a)-(three+rhoatinv)*deltap
        tmp5       = oneninth*rhob*(one-three*delta-delta11b)-(three+rhobtinv)*deltap
        dlypaaa    = omegap*dlypaa-omega*(tmp4-oneninth*rhoab*delta11b*rhotinv)
        dlypaab    = omegap*dlypab-omega*(oneninth*rhob*(p47-seven*delta)-seven*deltap-eightthird*rhot)
        dlypabb    = omegap*dlypbb-omega*(tmp5+oneninth*rhoab*delta11b*rhotinv-two*rhoa)
        tmp6       =-four*prma*denom*rhobtinv*(onethird*prmd*rhot43inv*denom*rhoa+one-rhoatinv)
        tmp7       =-face*omega*rhob*(omegap*rhoa*(rhoa83+rhob83)+(p11third*rhoa83+rhob83))
        tmp8       = dlypaaa*gradaa+dlypaab*gradab+dlypabb*gradbb
        tmp9       = clyp*(dlypaa+dlypab+dlypbb)
        tmp4       = oneninth*rhoa*(one-three*delta-delta11a)-(three+rhoatinv)*deltap
        tmp5       = oneninth*rhoa*(one-three*delta-delta11b)-(three+rhobtinv)*deltap
        dlypbaa    = omegap*dlypaa-omega*(tmp4+oneninth*rhoab*delta11a*rhotinv-two*rhob)
        dlypbab    = omegap*dlypab-omega*(oneninth*rhoa*(p47-seven*delta)-seven*deltap-eightthird*rhot)
        dlypbbb    = omegap*dlypbb-omega*(tmp5-oneninth*rhoab*delta11a*rhotinv)
        tmp10      =-four*prma*denom*rhoatinv*(onethird*prmd*rhot43inv*denom*rhob+one-rhobtinv)
        tmp11      =-face*omega*rhoa*(omegap*rhob*(rhoa83+rhob83)+(rhoa83+p11third*rhob83))
        tmp12      = dlypbaa*gradaa+dlypbab*gradab+dlypbbb*gradbb
        va         = clyp*(tmp6+tmp7+tmp8)
        vda        = clyp*(two*dlypaa*grhoa+dlypab*grhob)
        vb         = clyp*(tmp10+tmp11+tmp12)
        vdb        = clyp*(two*dlypbb*grhob+dlypab*grhoa)
        return
      END SUBROUTINE LYP
      END MODULE cam 
