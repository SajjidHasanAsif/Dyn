C +--------------------------------------------------------------------+
C |                        VUMAT BEGINS HERE                           |
C +--------------------------------------------------------------------+
      subroutine vumat( 
! Read only -
     $     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     $     stepTime, totalTime, dt, cmname, coordMp, charLength,
     $     props, density, strainInc, relSpinInc,
     $     tempOld, stretchOld, defgradOld, fieldOld,
     $     stressOld, stateOld, enerInternOld, enerInelasOld,
     $     tempNew, stretchNew, defgradNew, fieldNew,
! Write only -
     $     stressNew, stateNew, enerInternNew, enerInelasNew )

      include 'vaba_param.inc'

      dimension props(nprops), density(nblock),
     $coordMp(nblock,*),
     $charLength(*), strainInc(nblock,ndir+nshr),
     $relSpinInc(*), tempOld(*),
     $stretchOld(*), defgradOld(*),
     $fieldOld(*), stressOld(nblock,ndir+nshr),
     $stateOld(nblock,nstatev), enerInternOld(nblock),
     $enerInelasOld(nblock), tempNew(*),
     $stretchNew(*), defgradNew(*), fieldNew(*),
     $stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     $enerInternNew(nblock), enerInelasNew(nblock)

      character*80 cmname
C      write(*,'(A," HAS BEEN CALLED.")') CMNAME

      IF(CMNAME(1:5).EQ."X2DWC") THEN

         CALL MaterialWC( 
     $ nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     $ stepTime, totalTime, dt, cmname, coordMp, charLength,
     $ props, density, strainInc, relSpinInc,
     $ tempOld, stretchOld, defgradOld, fieldOld,
     $ stressOld, stateOld, enerInternOld, enerInelasOld,
     $ tempNew, stretchNew, defgradNew, fieldNew,
     $ stressNew, stateNew, enerInternNew, enerInelasNew )

C      WRITE(*,'(A," HAS BEEN CALLED.")') CMNAME

      END IF
      RETURN
      END 
c! +-----------------------------------------------------------------+
c! |                   SUBROUTINE MATERIALWC                         |
c! +-----------------------------------------------------------------+
      subroutine MaterialWC(
     $     nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     $     stepTime, totalTime, dt, cmname, coordMp, charLength,
     $     props, density, strainInc, relSpinInc,
     $     tempOld, stretchOld, defgradOld, fieldOld,
     $     stressOld, stateOld, enerInternOld, enerInelasOld,
     $     tempNew, stretchNew, defgradNew, fieldNew,
     $     stressNew, stateNew, enerInternNew, enerInelasNew )

      include 'vaba_param.inc'

c  ! Explanations of the arguments of subroutine VUMAT
c  !
c  ! Variables to be defined
c  !
c  ! stressNew(:) : stresses at the end of the increment,                   (output)
c  ! stateNew(:)  : state variables at the end of the increment,            (output)
c  !
c  ! Variables that can be updated
c  !
c  ! enerInternNew(:): Internal energy per unit mass at each mat. point at
c  !                the end of the increment                                (output)
c  ! enerInelasNew(:): Dissipated inelastic energy per unit mass at each
c  !                mat. point at the end of the increment                  (output)
c  !
c  ! Variables for information (or input only)
c  !
c  ! nblock       : no. of mat. pts to be processed in this call to vumat    (input)
c  ! ndir         : no. of direct components in a symmetric tensor           (input)
c  ! nshr         : no. of indirect components in a symmetric tensor         (input)
c  ! nstatev      : no. of user-defined state vars for this material         (input)
c  ! nfieldv      : no. of user-defined external field vars                  (input)
c  ! nporps       : no. of user-defined mat. properties                      (input)
c  ! lanneal      : =1 indicates that vumat is called during annealing       (input)   
c  ! stepTime     : time since the beginning of step                         (input)
c  ! totalTime    : total time since the beginning of first step             (input)
c  ! dt           : time increment size                                      (input)
c  ! cmname       : user specified material name, left justified.   
c  ! coordMp(:,:) : mat. pt. coords.                                         (input)
c  ! charLength(:): characteristic element length                            (input)
c  ! props(:)     : user-supplied material props                             (input)
c  ! density(:)   : current density at mat. pts. in the midstep conf.        (input)
c  ! strainInc(:,:) : strain increment tensor at each material point         (input)
c  ! relSpinInc(:,:): incremental relative rotation vector at each mat. pt.  (input)
c  ! tempOld(:,:)   : temperatures at the mat. pts. at the beginning of step (input)
c  ! stretchOld(:,:): stretch tensor U, at each mat. pt. at the beg. of step (input)
c  ! defgradOld(:,:): def. grad. tensor F, at each mat.pt. at beg. of step   (input)
c  ! fieldOld(:,:)  : values of the user defined field vars at each mat. pt.
c  !                at the beginning of the step                             (input)
c  ! stressOld(:,:) : stress tensor at each mat. pt. at the beg. of the step (input)
c  ! stateOld(:,:)  : state variables at each mat. pt. at the beg. of step   (input)
c  ! enerInternOld(:) : internal energy per unit mass at each mat. pt. at  
c  !                    the beg. of the incr.                                (input)
c  ! enerInelasOld(:) : dissipated inelastic energy per unit mass at each  
c  !                    mat. pt. at the beg. of the incr.                    (input)
c  ! tempNew(:)   : temperatures at each mat. pt. at the end of increment    (input)
c  ! stretchNew(:,:): stretch tensor, U, at each mat. pt. at the end of incr.(input)
c  ! defgradNew(:,:):  of how to provide some of the material specific
c  ! fieldNew(:,:)  : values of user-defined field vars at each mat.pt. at the
c  !                end of increment
c  !
c  ! OTHER USEFUL INFORMATION 
c  ! (1) shear strains are used instead of shear angles as opposed to UMAT.  
c  ! (2) strain order is \eps_1, \eps_2, \eps_3, \eps_{12}, \eps_{23}, \eps_{13}.
c  ! (3) when the local system is defined over the elements for which this routine
c  !     is called at every Gauss point, the rotation matrix is not needed, because
c  !     stresses and strains in this subroutine are defined with respect to that
c  !     local coordinate system.
c  !
c  ! All arrays dimensioned by (*) are not used in this algorithm
      dimension props(nprops), density(nblock),
     $     coordMp(nblock,*),
     $     charLength(*), strainInc(nblock,ndir+nshr),
     $     relSpinInc(*), tempOld(*),
     $     stretchOld(*), defgradOld(nblock,ndir+nshr+nshr),
     $     fieldOld(*), stressOld(nblock,ndir+nshr),
     $     stateOld(nblock,nstatev), enerInternOld(nblock),
     $     enerInelasOld(nblock), tempNew(*),
     $     stretchNew(*), defgradNew(nblock,ndir+nshr+nshr),
     $     fieldNew(*),
     $     stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     $     enerInternNew(nblock), enerInelasNew(nblock)

      common /matrix_prop/Em0, xnum0, Gm0, eta, xeta
      common /fiber_prop/Eaxf0, Elatf0, G12f0, G23f0, xnu12f0, xnu23f0
      common /matparam/xR1, xR2, xq, xpar, epsn1_f, epsn0_fc, xH, xn
      common /damage_prop/epsc_m, eps0_m, epsnc_f, epsn0_f, epstc_f, 
     $epst0_f, epsac_f, epsa0_f, epsbc_f, epsb0_f, epscc_f, epsc0_f,
     $epsm_f  
      common /common_vars/ PI, alpha, vf_my, vf_y, vf_m, vf
      common /debonding/xa, xp
      common /globaltolocal/ RRg2l1, RRg2l2, RRl2g1, RRl2g2, 
     $xInvRRg2l1, xInvRRl2g1, RG, theta_g2l, theta_l2g
      common /projection_vars/ xfyNijNkl, xfyMijMkl, xfyLijLkl,
     $xfyNijMkl, xfyNijLkl, xfyMijLkl,
     $xfyAijAkl, xfyBijBkl, xfyCijCkl,
     $xwyNijNkl, xwyMijMkl, xwyLijLkl,
     $xwyNijMkl, xwyNijLkl, xwyMijLkl,
     $xwyAijAkl, xwyBijBkl, xwyCijCkl,
     $xfyNij, xfyMij, xfyLij, xfyAij, xfyBij, xfyCij,
     $xwyNij, xwyMij, xwyLij, xwyAij, xwyBij, xwyCij
      common /changed_model_prop/ eps0, epsk, epsf, sigk, Ett, Ftt
      common /changed_model_prop/ eps0c, epskc, epsfc, sigkc, Ecc, Fcc
      common /shape_parameter/ sp_m, sp_t, sp_c

      integer np
      parameter (np=6)

      dimension RRg2l1(6,6), RRl2g1(6,6)
      dimension RRg2l2(6,6), RRl2g2(6,6)
      dimension xInvRRg2l1(6,6), xInvRRl2g1(6,6)
      dimension RG(3,3), eps33(3,3), deps33(3,3)
      dimension eps33_local(3,3), deps33_local(3,3)

      dimension xfyNij(6,6), xfyMij(6,6), xfyLij(6,6)
      dimension xfyAij(6,6), xfyBij(6,6), xfyCij(6,6)

      dimension xwyNij(6,6), xwyMij(6,6), xwyLij(6,6)
      dimension xwyAij(6,6), xwyBij(6,6), xwyCij(6,6)

      dimension xfyNijNkl(6,6,6), xfyMijMkl(6,6,6), xfyLijLkl(6,6,6)
      dimension xfyNijMkl(6,6,6), xfyNijLkl(6,6,6), xfyMijLkl(6,6,6)
      dimension xfyAijAkl(6,6,6), xfyBijBkl(6,6,6), xfyCijCkl(6,6,6)

      dimension xwyNijNkl(6,6,6), xwyMijMkl(6,6,6), xwyLijLkl(6,6,6)
      dimension xwyNijMkl(6,6,6), xwyNijLkl(6,6,6), xwyMijLkl(6,6,6)
      dimension xwyAijAkl(6,6,6), xwyBijBkl(6,6,6), xwyCijCkl(6,6,6)

      character*80 cmname
      dimension sig(6), eps(6)
      dimension eps_local(6), deps_local(6)
      dimension sig_old(6), eps_old(6)
      dimension deps(6), dsig(6)
      dimension stiff_lamina(6,6), compg(6,6)
      dimension stiff_lamina_local(6,6), stiff_lamina_global(6,6)
      dimension stiff_lamina_local_dam(6,6)
      dimension stiff_lamina_rot1(6,6), stiff_lamina_rot2(6,6)
      dimension comp_rot(6,6), comp_fyp(6,6)
      dimension stiff_m(6,6), stiff_fyp(6,6), stiff_wyp(6,6)
      dimension epsn_fy(6), depsn_fy(6), epsn_wy(6),depsn_wy(6)
      dimension wdtfy0(6),  wdtfy(6)
      dimension wdtfyi(6), wdtwyi(6)
      dimension wdtwy0(6), wdtwy(6) 
      dimension wdcfy0(6), wdcwy0(6), wdcfy(6), wdcwy(6)
      dimension wdcfyi(6), wdcwyi(6)
      dimension epsfy_mxt(6), epsfy_mxc(6)
      dimension epswy_mxt(6), epswy_mxc(6)


      integer jp, i, ind, IC, iflag

      PI=3.1415926535897932

      NTENS = 6

c--------------------------------------------------------------------------
c  Set material properties and the microplane system
c--------------------------------------------------------------------------
      theta_g2l = props(1)*PI/180
      theta_l2g = -props(1)*PI/180
      
      if (totalTime<=dt) then
         call inputparams()
         call setsystem()
      end if

      do i = 1, nblock
         strainInc(i,3)=-0.21863*(strainInc(i,1)+strainInc(i,2))
         do IC=1,4
            eps_old(IC)=stateOld(i,IC+38)
         end do
         do IC=1,4
            eps(IC)=eps_old(IC)+strainInc(i,IC)
            stateNew(i,IC+38)=eps(IC)
         end do
         eps(5)=0.0d0
         eps(6)=0.0d0
         eps_old(5)=0.0d0
         eps_old(6)=0.0d0

         DO IC = 1, NTENS
            IF (IC > 3) THEN    ! SHEAR STRAINS ARE SHEAR ANGLES IN MICROPLANE MODEL
               eps(IC)=2.0*eps(IC)
               eps_old(IC)=2.0*eps_old(IC)
            END IF
         END DO

         DO IC = 1, NTENS
            IF (IC > 3) THEN    ! SHEAR STRAINS ARE SHEAR ANGLES IN MICROPLANE MODEL
               deps(IC)=2.0*strainInc(i,IC)
            ELSE
               deps(IC)=strainInc(i,IC)
            END IF
         END DO
         deps(5)=0.d0
         deps(6)=0.d0

c  Rotate strain tensor to material coordinate system (1 = fill, 2 = warp)
C -- Convert eps and deps to 3x3 matrix
      call voigt_to_tensor(eps, eps33)
      call voigt_to_tensor(deps, deps33)

C -- Rotate eps and deps to local coordinates, abt z axis
c      eps33_local = matmul(transpose(RG),(matmul(eps33,RG)))
c      deps33_local = matmul(transpose(RG),(matmul(deps33,RG)))
        eps33_local = eps33
        deps33_local = deps33
    

C -- Convert eps33_local and deps33_local to voigt form
      call tensor_to_voigt(eps33_local, eps_local)
      call tensor_to_voigt(deps33_local, deps_local)
      

c     ! +---------------------------------------------------------------------------+
c     ! |                              PART I: MATRIX                               |
c     ! +---------------------------------------------------------------------------+
        eps_eff0  = stateOld(i,1)
        wdm0      = stateOld(i,2)
        wdm_eq0   = stateOld(i,45)                                                        ! Matrix Damage Variable SDV 45
        if(totalTime<=dt) then
           wdm0    = 0.0d0
           wdm_eq0 = 0.0d0
        end if
        xsum = 0.0d0
		xsum1 = 0.0d0
        do IC= 1,6
           xsum=xsum+(eps_local(IC)*eps_local(IC))
           xsum1=xsum1+(deps_local(IC)*deps_local(IC))
        end do
        eps_eff   = sqrt(2.d0*xsum/3.d0)                                                    ! Effective strian
        deps_eff  = sqrt(2.d0*xsum1/3.d0)
        eps_max       = max(eps_eff, stateOld(i,70))
        stateNew(i,70)= eps_max
        deps_effi     = abs(eps_eff - stateOld(i,70))                                                     ! Increment in effective strain
!        deps_eff      = max(deps_effi,stateOld(i,71))
        !stateNew(i,71)= deps_eff
        t1a           = 0.d0
        t2a           = 0.d0
        t3            = 0.d0
        eps_eff_nl    = 0.0d0
		!sr_m = deps_eff/dt
        sr_m   = (stateOld(i,75)+(deps_eff/dt))/2.0d0                                 ! matrix strain rate
        sr_m  = min(sr_m, 6.5d2)
        stateNew(i,75) = sr_m
        
C Viscoelastic formulation
        if(totalTime <= dt) then
            E0    = Em0
			sr_m  = 0.0d0
			eps_max=eps_eff
!			xHD  =  3.0d1
			wdm2 = 0.0d0  
!            sigma0= eps_eff*Em0
        else 
            E0    = stateold(i,73)
!            sigma0= stateOld(i,72)
        end if
!        dsigma        =Em0*(deps_eff-((sigma0*dt)/eta))
!        sigma_i       = stateold(i,72) + dsigma
!        StateNew(i,72)= sigma_i
!        if (dsigma==0 .or. deps_eff==0) then
!            E_vis = stateold(i,74)
!        else 
!            E_vis  =(abs(dsigma)/abs(deps_eff))
!        end if 
         xeta = 1.21D0
		 if(eps_max<1.0d-8) then
		   E_vis = Em0
          else 
		   E_vis = Em0 !+ (xeta*sr_m/eps_max) ! Kelvin Voigt model
		 end if
        if(totalTime < 3.d0*dt) then
          write(*,'(a,f18.5)') "E_vis = ", E_vis
          write(*,'(a,f18.5)') "eta = ", xeta
          write(*,'(a,f18.5)') "sr_m = ", sr_m
          write(*,'(a,f18.5)') "eps_eff = ", eps_eff		  
        end if
         stateNew(i,74) = E_vis
! Matrix rate dependent formulation
      sp_m1   = 2.0D3
	  sp_m2   = 5.0D1
	  sp_m3   = 3.0D0
	  sp_m4   = 5.0D1
      eps0_mD = eps0_m!*((1.d0+(abs(sr_m)/sp_m1)))!,0.95*epsc_mD)    ! Critical strain at damage initiation
      epsc_mD = epsc_m!/((1.d0+(abs(sr_m)/sp_m2)))     ! Critical strain at failure
      epsf_mD = epsm_f!*(1.d0+(abs(sr_m)/sp_m3))
      xH = 8.0D1 !116.67D0
	  xHD = xH!*(1.d0+(abs(sr_m)/sp_m4))
!        if(totalTime <= 3*dt) then
           EbyH = E_vis/xHD
           stateNew(i,80)=EbyH
!        else
!		   stateNew(i,80)=stateOld(i,80)
!        end if
      xsl     = stateNew(i,80)
      stateNew(i,76)=eps0_mD
      stateNew(i,77)=epsc_mD
      stateNew(i,78)=epsf_mD
	  !stateNew(i,80)=EbyH

C Damage variable for matrix part
      if(eps_max < eps0_mD) then
          wdmi         = 0.0d0
		  wdm2         = 0.0d0
		  wdm2 = max(wdm2, stateOld(i,79))
          stateNew(i,71)= 1.d0 
      else if (eps_max >= eps0_mD .and. eps_max< epsc_mD) then
          eps_eff_nl  = max((eps_max-eps0_mD),0.d0)
          t1a  = eps0_mD/eps_max
          t2a  = ((eps_eff_nl)**xn)/(xsl*eps_max)
          wdmi = max(1.d0-(t1a+t2a),0.d0)
          wdmi = min(wdmi, 0.9999d0)
          stateNew(i,71)= 2.d0 
		  wdm2 = 0.d0
 		  wdm2 = max(wdm2, stateOld(i,79))
!      else if (eps_max >= epsc_mD) then
       else if (eps_max >= epsc_mD .and. eps_max <= epsf_mD) then
!          wdmi  = min(stateOld(i,45) + ((eps_max-epsc_mD)*
!     $            (stateOld(i,45)-0.9999d0))/(epsc_mD-epsf_mD),0.9999d0)       
          wdm2 = min((eps_max - epsc_mD)/(epsf_mD-epsc_mD), 0.9999d0)!
		  wdm2 = max(wdm2, stateOld(i,79))
           wdmi = stateOld(i,45)
          stateNew(i,71)= 3.d0 
       else if (eps_max > = epsm_fD) then 
           wdmi = 0.9999d0        
 		  wdm2 = max(wdm2, stateOld(i,79))
           wdmi = stateOld(i,45)
           stateNew(i,71)= 4.d0 
      end if
C Prevention of damage reversibility in Matrix
      wdm      = max(wdmi,stateOld(i,45))
      stateNew(i,45)  = wdm
      stateNew(i,79)  = wdm2

C Calculate the elastic stress tensor
        call kelast_isotropic(wdm, wdm2, stiff_m, E_vis, eee)
!        call kelast_isotropic(wdm, stiff_m)
        stateNew(i,73) = eee
        
c     ! +---------------------------------------------------------------------------+
c     ! |              PART II: FILL YARN PLATE                                     |
c     ! +---------------------------------------------------------------------------+
         iplate=1

        do jp=1,np
           wdtfy0(jp)     = stateOld(i,2+jp)                                                ! Fill yarn damage variable under tension SDV- 3 to 8
           wdcfy0(jp)     = stateOld(i,8+jp)                                                ! Fill yarn damage variable under compression SDV- 9 to 14
        end do     
C Damage varible at the begining of time increment        
        if(totalTime<=dt) then
          wdtfy0         = 0.d0                                                                ! Damage variable in n direction is zero at the begining of time
          wdcfy0         = 0.d0                                                                ! Damage variable in t direction is zero at the begining of time
          wdtfy          = 0.d0
		wdcfy          = 0.d0
   		wdtfyi         = 0.d0
		wdcfyi         = 0.d0
        end if
C Initial strain and strain increment at the begining of the simulation in n,m,l direction-
      do jp=1,6
          epsn_fy(jp)    = 0.d0
          depsn_fy(jp)   = 0.d0
C Defining strain; and strain strain increment in the normal direction
          do IC=1,6
             epsn_fy(jp) = epsn_fy(jp)+(eps_local(IC)*xfyNij(IC,jp))
             depsn_fy(jp)= depsn_fy(jp)+(deps_local(IC)*xfyNij(IC,jp))
          end do
C Determination of max tensile strain
          epsn_fypos = max(epsn_fy(jp),0.d0)
          epsfy_mxt(jp) = max(epsn_fypos, stateOld(i,45+jp))
          stateNew(i,45+jp) = epsfy_mxt(jp)
      end do
C Calculation of damage variable for strain softening
      do jp=1,np
C Fill yarn in tension
!         if(epsn_fy(jp)>=0)then												        ! Fill yarn in tension 
C Increment in effective strain
              dep		= max(epsn_fy(jp)-epsn0_f,0.d0)                 
              deps_pp	= max(depsn_fy(jp),0.d0)
              sr		= deps_pp/dt                                        ! strain rate in fill yarn
C Rate Dependent formulation                                                           
			epskD	= epsk!*(1+(abs(sr)/sp_t))
			epsfD	= epsf!*(1+(abs(sr)/sp_t))
			sigkD	= sigk!*(1+(abs(sr)/sp_t))
			eps0D	= eps0!min(eps0*(1+(abs(sr)/sp_t)),0.9999*epsfD)
C Calculation of damage varible for tension in yarn plate
             if ((epsfy_mxt(jp)<eps0D))then
                 wdtfyi(jp)= 0d0
          else if((epsfy_mxt(jp)>=eps0D).and.(epsfy_mxt(jp)<epskD))then                 ! strain goes beyond initial critical strain but smaller than knee point
                 wdtfyi(jp)= max(((epskD-sigkD/Ett)/(epskD-eps0D))
     $						*(1-(eps0D/epsn_fy(jp))),0.d0)
             else if (epsfy_mxt(jp)>=epskD.and.epsfy_mxt(jp)<espfD)then                ! strian goes beyond knee point
                 wdtfyi(jp)= max(1-((sigkD/Ett)/(epsfD-epskD))
     $						*((epsfD/epsn_fy(jp))-1),0.d0)
             else if (epsfy_mxt(jp)>=epsfD) then                                      ! strain goes beyond final damage threshold
                  wdtfyi(jp) = 0.9999d0
             end if
!             wdcfyi(jp) = 0d0
!          end if
C Condition to prevent riversibility of damage
              wdtfy(jp) = max(wdtfyi(jp),stateOld(i,2+jp))
              stateNew(i,2+jp) = wdtfy(jp)
c Fill yarn in compression (modify)
!          if(epsn_fy(jp)<0) then
C Determination of max compressive strain
              epsn_fyneg  = min(epsn_fy(jp), 0.d0)
              epsfy_mxc(jp)       = min(epsn_fyneg, stateOld(i,51+jp))
              stateNew(i,51+jp)   = epsfy_mxc(jp)
              xabse               = abs(epsfy_mxc(jp))
          
C Damage Varibale for fill yarn under compression       
              if((xabse <= epskc).and.(xabse >=eps0c)) then
                  wdcfyi(jp)= max(0.d0,((1.d0-(eps0c/xabse))
     $                        *((epskc-(sigkc/Ecc))/(epskc-eps0c))))
              else if (xabse > epskc) then
                  wdcfyi(jp)= max(0.d0,(1.d0-(sigkc/Ecc)/xabse))  
              else
                  wdcfyi(jp)=0d0 
              end if
!              wdtfyi(jp) = 0d0
!          end if
C Condition to prevent reversibility of damage variable
              !wdcfy(jp) = 0d0!min(max(wdcfyi(jp),stateOld(i,8+jp)),0.9999d0)    ! T&C overlap 
              wdcfy(jp) = 0
              stateNew(i,8+jp) = wdcfy(jp)                                               ! no compression
        end do

        call kelast_yarnplate(wdm, wdm2, wdtfy, wdcfy, iplate, stiff_fyp)
c     ! +---------------------------------------------------------------------------+
c     ! |              PART III: WARP YARN PLATE                                    |
c     ! +---------------------------------------------------------------------------+
        iplate=2

        do jp=1,np
           wdtwy0(jp)     = stateOld(i,14+jp)                                         ! warp yarn under tension SDV 15-20
           wdcwy0(jp)     = stateOld(i,20+jp)                                         ! warp yarn under compression SDV 21-27
        end do
C Damage varible at the begining of time
        if(totalTime<=dt) then
          wdtwy0         = 0.d0
          wdcwy0         = 0.d0
          wdtwyi         = 0.d0
		wdcwyi         = 0.d0
   		wdtwy          = 0.d0
		wdcwy          = 0.d0
        end if
C Definining initial strain and strain increment in n,m,l direction for warp yarn- 
        do jp=1,6
           epsn_wy(jp) =0.d0
           depsn_wy(jp)=0.d0
C Calculation of strain and strain increment in triad direction-
           do IC=1,6
              epsn_wy(jp)=epsn_wy(jp)+(eps_local(IC)*xwyNij(IC,jp))
              depsn_wy(jp)=depsn_wy(jp)+(deps_local(IC)*xwyNij(IC,jp))
           end do
C Determination of maximum strain under tension in warp yarn
          epsn_wypos = max(epsn_wy(jp), 0.d0)
		  epswy_mxt(jp)       = max(epsn_wypos, stateOld(i,57+jp))
          stateNew(i,57+jp)   = epswy_mxt(jp)
        end do

        do jp=1,np
C Warp yarn in tension
!          if(epsn_wy(jp)>=0)then
C Strain increment
			dep             = max(epsn_wy(jp)-epsn0_f,0.d0)
			deps_pp	        = max(depsn_wy(jp),0.d0)
              sr		        = deps_pp/dt                                ! strian rate 
C Rate Dependent formulation                                                       
			epskD	= epsk!*(1+(abs(sr)/sp_t))
			epsfD	= epsf!*(1+(abs(sr)/sp_t))
			sigkD	= sigk!*(1+(abs(sr)/sp_t))
			eps0D	= eps0!min(eps0*(1+(abs(sr)/sp_t)),0.9999*epsfD)
C Calculation of damage varibele for warp warn under tension
             if((epswy_mxt(jp)<eps0D)) then
                  wdtwyi(jp) = 0d0 
		    else if((epswy_mxt(jp)>=eps0D).and.(epswy_mxt(jp)<epskD)) then
                  wdtwyi(jp) = max (((epskD-sigkD/Ett)/(epskD-eps0D))
     $						    *(1-(eps0D/epsn_wy(jp))),0.d0)
		    else if (epswy_mxt(jp) >= epskD.and.epswy_mxt(jp)<epsfD) then
			    wdtwyi(jp) = max (1-((sigkD/Ett)/(epsfD-epskD))
     $						    *((epsfD/epsn_wy(jp))-1),0.d0)
              else if (epswy_mxt(jp) >= epsfD) then
                  wdtwyi(jp) = 0.9999d0
              end if 
!              wdcwyi(jp) = 0d0
!          end if
C Prevention of damage reversibility
              wdtwy(jp)= max(wdtwyi(jp),stateOld(i,14+jp))
              stateNew(i,14+jp)   = wdtwy(jp)  
c Warp yarn in compression
!          if(epsn_wy(jp)<0)then
C Determination of maximum strain under compression in warp yarn              
              epsn_wyneg = min(epsn_wy(jp), 0.d0)
              epswy_mxc(jp)       = min(epsn_wyneg, stateOld(i,63+jp))
              stateNew(i,63+jp)   = epswy_mxc(jp)
              xabse               = abs(epswy_mxc(jp))
C Damage Varibale for warp yarn under compression         
              if((xabse <= epskc).and.(xabse >= eps0c))then
                  wdcwyi(jp)=max(0.d0,((1.d0-(eps0c/xabse))
     $                            *((epskc-(sigkc/Ecc))/(epskc-eps0c))))
              else if (xabse > epskc)then
                  wdcwyi(jp)=max(0.d0,(1.d0-(sigkc/Ecc)/xabse)) 
              else
                  wdcwyi(jp)=0d0 
              end if
!                 wdtwyi(jp)= 0d0
!           end if 
C Prention of damage reversibility        
              !wdcwy(jp)= 0d0 !min(max(wdcwyi(jp),stateOld(i,20+jp)),0.9999d0)
              wdcwy(jp) = 0.d0
              stateNew(i,20+jp) = wdcwy(jp)                            ! no compression damage
        end do

        call kelast_yarnplate(wdm, wdm2, wdtwy, wdcwy, iplate, stiff_wyp)

c     !     -------------------------------------------------
c     !     Lamina Stiffness tensor
c     !     -------------------------------------------------
         stiff_lamina_local =(1-vf)*stiff_m+(vf/2)*stiff_fyp+
     $(vf/2)*stiff_wyp

C   --Debonding damage--
      wddb0=stateOld(i,27)
      gam=abs(2.d0*eps33_local(1,2))
      dgam=abs(2.d0*deps33_local(1,2))
      if(totalTime<=dt) then
           wddb0=0.d0
           gam=0.d0
          dgam=0.d0
      end if
      dwddb=(0.0d0)*max((dgam*xp*xa*(gam**(xp-1.d0))*
     $exp(-xa*(gam**xp))),0.d0)
      wddb=min(wddb0+dwddb,0.999)
      do IC=1,6
         do JC=1,6
            stiff_lamina_local_dam(IC,JC)=stiff_lamina_local(IC,JC)
         end do
      end do
      stiff_lamina_local_dam(4,4)=(1.d0-wddb)*stiff_lamina_local(4,4)
      stiff_lamina_local_dam(5,5)=(1.d0-wddb)*stiff_lamina_local(5,5)
      stiff_lamina_local_dam(6,6)=(1.d0-wddb)*stiff_lamina_local(6,6)

      stateNew(i,27)=wddb
c   --End debonding damage --

C -- Back rotate stiff_lamina to global coordinates
c      stiff_lamina = matmul(xInvRRl2g1,
c     $(matmul(stiff_lamina_local,RRl2g2)))
       stiff_lamina = stiff_lamina_local

c---------write stiffness tensors---------------------------------------------
      call M66INV(stiff_lamina, compg)
      if(i<0)then
      if(totalTime<=dt) then
      write(*,'(a)') "In MaterialWC"
      write(*,'(a)') "stiff_m = "
      do IC=1,6
         write(*,'(6f18.6)') (stiff_m(IC,j), j=1,6)
      end do
      write(*,'(a)') "stiff_fyp = "
      do IC=1,6
         write(*,'(6f18.6)') (stiff_fyp(IC,j), j=1,6)
      end do
      write(*,'(a)') "stiff_wyp = "
      do IC=1,6
         write(*,'(6f18.6)') (stiff_wyp(IC,j), j=1,6)
      end do
      write(*,'(a)') "stiff_lamina_local = "
      do IC=1,6
         write(*,'(6f18.6)') (stiff_lamina_local(IC,j), j=1,6)
      end do
      write(*,'(a)') "stiff_lamina_local_dam = "
      do IC=1,6
         write(*,'(6f18.6)') (stiff_lamina_local_dam(IC,j), j=1,6)
      end do
      write(*,'(a)') "stiff_lamina  = "
      do IC=1,6
         write(*,'(6f18.6)') (stiff_lamina(IC,j), j=1,6)
      end do
      write(*,'(a)') "lamina compliance  = "
      do IC=1,6
         write(*,'(6g15.4)') (compg(IC,j), j=1,6)
      end do
      write(*,'(a,f18.5)') "E1 = ", 1.0d0/compg(1,1)
      write(*,'(a,f18.5)') "E2 = ", 1.0d0/compg(2,2)
      write(*,'(a,f18.5)') "E3 = ", 1.0d0/compg(3,3)
      write(*,'(a,f18.4)') "G12 = ", 1.0d0/compg(4,4)
      write(*,'(a,f18.4)') "G23 = ", 1.0d0/compg(5,5)
      write(*,'(a,f18.4)') "G31 = ", 1.0d0/compg(6,6)
      write(*,'(a,f18.3)') "nu21 = ", -1.0d0*compg(1,2)/compg(2,2)
      write(*,'(a,f18.3)') "nu12 = ", -1.0d0*compg(1,2)/compg(1,1)
      write(*,'(a,f18.3)') "nu31 = ", -1.0d0*compg(1,3)/compg(3,3)
      write(*,'(a,f18.3)') "nu13 = ", -1.0d0*compg(1,3)/compg(1,1)
      write(*,'(a,f18.3)') "nu32 = ", -1.0d0*compg(2,3)/compg(3,3)
      write(*,'(a,f18.3)') "nu23 = ", -1.0d0*compg(2,3)/compg(2,2)
      end if
      end if
      stateNew(i,28)=1.0d0/compg(1,1)
      stateNew(i,29)=1.0d0/compg(2,2)
      stateNew(i,30)=1.0d0/compg(3,3)
      stateNew(i,31)=1.0d0/compg(4,4)
      stateNew(i,32)=1.0d0/compg(5,5)
      stateNew(i,33)=1.0d0/compg(6,6)
      stateNew(i,34)=-1.0d0*compg(1,2)/compg(2,2)
      stateNew(i,35)=-1.0d0*compg(1,3)/compg(3,3)
      stateNew(i,36)=-1.0d0*compg(3,2)/compg(2,2)
c---------------------------------------------------------------------

      sig = matmul(stiff_lamina, eps)

      DO IC = 1, NTENS
         stressNew(i,IC) = sig(IC)
      END DO
      
      DO IC = 1,NTENS
         dsig(IC)=stressNew(i,IC)-stressOld(i,IC)
      END DO

c     ! Report the new values of the state variables
      stateNew(i,1)=eps_eff
      stateNew(i,2)=wdm
      !do jp=1,np
       ! stateNew(i,2+jp)=wdtfy(jp)
        !stateNew(i,8+jp)=wdcfy(jp)

!        stateNew(i,14+jp)=wdtwy(jp)
!        stateNew(i,20+jp)=wdcwy(jp)
      !end do
      stateNew(i,37)=1
c      wsum1=0.d0
c      wsum2=0.d0
c      wsum3=0.d0
c      wsum4=0.d0
c      do jp=1,6
c         wsum1=wsum1+stateNew(i,2+jp)
c         wsum2=wsum2+stateNew(i,8+jp)
c         wsum3=wsum3+stateNew(i,14+jp)
c         wsum4=wsum4+stateNew(i,20+jp)
c      end do
c      wsum=wsum1+wsum3
c      if(CMNAME(1:6).eq."X2DWC1")then
C        if(stateNew(i,2)>0.5) then
c          if(wsum1>3.0) then
c            stateNew(i,37)=0
c          end if
c          if(wsum3>3.0) then
c             stateNew(i,37)=0
c          end if
C        end if
c      end if
c      if(stateNew(i,2)>0.99) then
c        if(wsum1>5.9) then
c           stateNew(i,37)=0
c        end if
c        if(wsum3>5.9) then
c           stateNew(i,37)=0
c        end if
c      end if
!c      if(stateNew(i,1)>0.45d0) then
!          stateNew(i,37)=0
!c      end if

c     ! Update the specific internal energy -
         stressPower = (( stressOld(i,1)+stressNew(i,1) )*
     $strainInc(i,1) + ( stressOld(i,2)+stressNew(i,2) )*
     $strainInc(i,2) + 2.0*( stressOld(i,4)+stressNew(i,4) )*
     $strainInc(i,4) )/2.0
c     !
      enerInternNew(i) = enerInternOld(i) + stressPower / density(i)
      if(totalTime<dt) then
           enerInternNew(i)=0.d0
      end if

c     !
c     ! Update the dissipated inelastic specific energy -
      smean = (stressNew(i,1) + stressNew(i,2) +
     $     stressNew(i,3))/3.0
      equivStress = sqrt( 3.0/2.0 * ( (stressNew(i,1)-smean)**2 +
     $     (stressNew(i,2)-smean)**2 + (stressNew(i,3)-smean)**2 +
     $     two * stressNew(i,4)**2 ) )

C      enerInelasNew(i) = enerInelasOld(i)+(fractureWorkInc*4d0)
      enerInelasNew(i) = 0.d0
      stateNew(i,38)=deps_eff
      stateNew(i,43)=dwdm
      stateNew(i,44)=wdm_eq
      end do
      end 
C +--------------------------------------------------------------------+
C |                        SUBROUTINE INPUTPARAMS                      |
C +--------------------------------------------------------------------+
      subroutine inputparams() 
      include 'vaba_param.inc'
      common /matrix_prop/Em0, xnum0, Gm0, eta, xeta
      common /fiber_prop/Eaxf0, Elatf0, G12f0, G23f0, xnu12f0, xnu23f0
      common /matparam/xR1, xR2, xq, xpar, epsn1_f, epsn0_fc, xH, xn
      common /damage_prop/epsc_m, eps0_m, epsnc_f, epsn0_f, epstc_f,
     $epst0_f, epsac_f, epsa0_f, epsbc_f, epsb0_f, epscc_f, epsc0_f,  
     $epsm_f
      common /debonding/xa, xp
      common /common_vars/ PI, alpha, vf_my, vf_y, vf_m, vf
      common /changed_model_prop/ eps0, epsk, epsf, sigk, Ett, Ftt
      common /changed_model_prop/ eps0c, epskc, epsfc, sigkc, Ecc, Fcc
      common /shape_parameter/ sp_m, sp_t, sp_c
      
      integer np
      parameter (np=6)

      PI=3.1415926535897932
C Young Modulus and Poisson ratio of matrix in the 2DWC.
      Em0  = 3500.0 !4400
      xnum0 = 0.35
      Gm0=Em0/(2*(1+xnum0))

C Max. undulation angle of yarns
      alpha = 0.06831  !tan-1(b/a)
      Eaxf0  = 190000 !221000.0    !MPa
      Elatf0 = 40000.0      !MPa
      G12f0   = 24000.0     !MPa
      G23f0   = 14300.0     !MPa
      xnu12f0 = 0.26
      xnu23f0 = 0.49
      vf_my = 0.3            !Matrix volume fraction WITHIN YARN
      vf_y  = 0.7            !Fiber volume fraction WITHIN YARN
      vf    = 0.7661 !0.8571 !Yarn volume fraction WITHIN RUC
      vf_m  = 1.d0-vf        !Matrix volume fraction WITHIN RUC

      eps0_m = 6.5D-3          !changes the initial slope in the shear curve, can be used to achieve higher stress
      epsc_m = 1.55D-1!2.8d-2  !changes the nonlinear part also affect the initial slope
      epsm_f = 4.8D-1  


      epsn0_f = 1.09d-02
!      epsnc_f = 2.0d-02 

!      epst0_f = 5.0d-02
!      epstc_f = 1.0d-01

!      epsa0_f = 1.5d10
!      epsac_f = 2.0d0*epsa0_f

!      epsb0_f = 1.81d10
!      epsbc_f = 2.0d0*epsb0_f

!      epsc0_f = epsa0_f
!      epscc_f = epsac_f
C used for debonding damage
      xa=1.0d4
      xp=4.0d0

!      xR1 = 8.2d-03
!      xR2 = 1.6d5
!      xq = 4.25d0
!      xpar = 3.5d4/5.63d4! slope of initial plateau/E0; 3.5d4/5.63d4
!      epsn1_f = 1.37d0*epsn0_f
!      epsn0_fc = 1.03d0*epsn0_f

      xH = 3d1
      xn = 1/1.75d0

!parameters for changed model
      ! under tension
      Ftt	 = 866D0
      eps0 = 1.45D-2
      epsk = 0.01D0!20.d0*2.d0/(600.0d0*2.d0)
      epsf = 0.16D0
      sigk = 300.0D0!(((78.0d0-40.d0)/2.0d0)*2.0d0)/(epsf-epsk)
      Ett  = Ftt/eps0         ! Elastic modulus      

      ! under compression
      Fcc   = 201D0 
      eps0c = 0.011!(1.62e-02)
      epskc = 0.03!(20.d0*2.d0/(600.0d0*2.d0))
      epsfc = 0.148d0
      sigkc = 330!((((78.0d0-40.d0)/2.0d0)*2.0d0)/(epsfc-epskc))
      Ecc   = 53300.d0!
      
      ! yarn shape parameter for rate dependent model
      sp_t = 1.5D3 !2676!23.2D2                   ! under tension calibrated to 14.2D2
      sp_c = 7.45D2                        ! under compression calibrated to 7.45D2 [not been used]

      ! Matrix shape parameter
      sp_m   = 1.5D01
      eta    = 1.21D-2
      end 
      
C +--------------------------------------------------------------------+
C |                        SUBROUTINE SETSYSTEM                        |
C +--------------------------------------------------------------------+
      subroutine setsystem() 
      include 'vaba_param.inc'
      common /matrix_prop/Em0, xnum0, Gm0, eta, xeta
      common /fiber_prop/Eaxf0, Elatf0, G12f0, G23f0, xnu12f0, xnu23f0
      common /matparam/xR1, xR2, xq, xpar, epsn1_f, epsn0_fc, xH, xn
      common /common_vars/ PI, alpha, vf_my, vf_y, vf_m, vf
      common /globaltolocal/ RRg2l1, RRg2l2, RRl2g1, RRl2g2, 
     $xInvRRg2l1, xInvRRl2g1, RG, theta_g2l, theta_l2g
      common /projection_vars/ xfyNijNkl, xfyMijMkl, xfyLijLkl,
     $xfyNijMkl, xfyNijLkl, xfyMijLkl,
     $xfyAijAkl, xfyBijBkl, xfyCijCkl,
     $xwyNijNkl, xwyMijMkl, xwyLijLkl,
     $xwyNijMkl, xwyNijLkl, xwyMijLkl,
     $xwyAijAkl, xwyBijBkl, xwyCijCkl,
     $xfyNij, xfyMij, xfyLij, xfyAij, xfyBij, xfyCij,
     $xwyNij, xwyMij, xwyLij, xwyAij, xwyBij, xwyCij

      dimension xfyNij(6,6), xfyMij(6,6), xfyLij(6,6)
      dimension xfyAij(6,6), xfyBij(6,6), xfyCij(6,6)

      dimension xwyNij(6,6), xwyMij(6,6), xwyLij(6,6)
      dimension xwyAij(6,6), xwyBij(6,6), xwyCij(6,6)

      dimension xfyNijNkl(6,6,6), xfyMijMkl(6,6,6), xfyLijLkl(6,6,6)
      dimension xfyNijMkl(6,6,6), xfyNijLkl(6,6,6), xfyMijLkl(6,6,6)
      dimension xfyAijAkl(6,6,6), xfyBijBkl(6,6,6), xfyCijCkl(6,6,6)

      dimension xwyNijNkl(6,6,6), xwyMijMkl(6,6,6), xwyLijLkl(6,6,6)
      dimension xwyNijMkl(6,6,6), xwyNijLkl(6,6,6), xwyMijLkl(6,6,6)
      dimension xwyAijAkl(6,6,6), xwyBijBkl(6,6,6), xwyCijCkl(6,6,6)

      integer jp, ij(2,6), i, j, k
      integer np
      parameter (np=6)
      dimension te(3,36), te_loc(3,36)
      dimension xn1(3), xm1(3), xl1(3)
      dimension xn2(3), xm2(3), xl2(3)
      dimension xa1(3), xb1(3), xc1(3)
      dimension xa2(3), xb2(3), xc2(3)
      dimension RG(3,3)
 
      dimension RRg2l1(6,6), RRl2g1(6,6)
      dimension RRg2l2(6,6), RRl2g2(6,6)
      dimension xInvRRg2l1(6,6), xInvRRl2g1(6,6)

      PI=3.1415926535897932
      
      nvhi = np+1

      call getRotMat(theta_g2l, RRg2l1, RRg2l2, xInvRRg2l1)
      call getRotMat(theta_l2g, RRl2g1, RRl2g2, xInvRRl2g1)
      RG(1,1)=cos(theta_g2l)
      RG(1,2)=-sin(theta_g2l)
      RG(1,3)=0
      RG(2,1)=sin(theta_g2l)
      RG(2,2)=cos(theta_g2l)
      RG(2,3)=0
      RG(3,1)=0
      RG(3,2)=0
      RG(3,3)=1

c  ! +-                    -+   +-                    -+
c  ! | sig_11 sig_12 sig_13 | _ | sig_1  sig_4  sig_6  |
c  ! |        sig_22 sig_23 | = |        sig_2  sig_5  | according to array [ij] given below:
c  ! |               sig_33 |   |               sig_3  |
c  ! +-                    -+   +-                    -+                                                  
c
c  ! ABAQUS EXPLICIT ORDERING:
c  ! version f90:      
c      ij=reshape((/1,1,2,2,3,3,1,2,2,3,1,3/),(/2,6/)) ! 2 rows and 6 columns; filled in columnwise
      ij(1,1)=1
      ij(2,1)=1
      ij(1,2)=2
      ij(2,2)=2
      ij(1,3)=3
      ij(2,3)=3
      ij(1,4)=1
      ij(2,4)=2
      ij(1,5)=2
      ij(2,5)=3
      ij(1,6)=1
      ij(2,6)=3
c  ! +------------------------------------------------------------------------------+
c  ! |                     The 6 n vectors for the fill yarns                       |
c  ! +------------------------------------------------------------------------------+
      te_loc(1,1)= cos(alpha)
      te_loc(2,1)= 0.0
      te_loc(3,1)= sin(alpha) 

      te_loc(1,2)= 1.0
      te_loc(2,2)= 0.0
      te_loc(3,2)= 0.0

      te_loc(1,3)= cos(alpha)
      te_loc(2,3)= 0.0
      te_loc(3,3)= -sin(alpha)

      te_loc(1,4)= cos (alpha)
      te_loc(2,4)= 0.0
      te_loc(3,4)= -sin(alpha)

      te_loc(1,5)= 1.0
      te_loc(2,5)= 0.0
      te_loc(3,5)= 0.0

      te_loc(1,6)= cos(alpha)
      te_loc(2,6)= 0.0
      te_loc(3,6)= sin(alpha)
c  ! +------------------------------------------------------------------------------+
c  ! |                     The 6 m vectors for the fill yarns                       |
c  ! +------------------------------------------------------------------------------+
      te_loc(1,7)= sin(alpha)
      te_loc(2,7)= 0.0
      te_loc(3,7)= -cos(alpha) 

      te_loc(1,8)= 0.0
      te_loc(2,8)= 0.0
      te_loc(3,8)= -1.0

      te_loc(1,9)= -sin(alpha)
      te_loc(2,9)= 0.0
      te_loc(3,9)= -cos(alpha)

      te_loc(1,10)= sin(alpha)
      te_loc(2,10)= 0.0
      te_loc(3,10)= cos(alpha)

      te_loc(1,11)= 0.0
      te_loc(2,11)= 0.0
      te_loc(3,11)= 1.0

      te_loc(1,12)= -sin(alpha)
      te_loc(2,12)= 0.0
      te_loc(3,12)= cos(alpha)
c  ! +------------------------------------------------------------------------------+
c  ! |                     The 6 l vectors for the fill yarns                       |
c  ! +------------------------------------------------------------------------------+
      te_loc(1,13)= 0.0
      te_loc(2,13)= 1.0
      te_loc(3,13)= 0.0 

      te_loc(1,14)= 0.0
      te_loc(2,14)= 1.0
      te_loc(3,14)= 0.0

      te_loc(1,15)= 0.0
      te_loc(2,15)= 1.0
      te_loc(3,15)= 0.0

      te_loc(1,16)= 0.0
      te_loc(2,16)= -1.0
      te_loc(3,16)= 0.0

      te_loc(1,17)= 0.0
      te_loc(2,17)= -1.0
      te_loc(3,17)= 0.0

      te_loc(1,18)= 0.0
      te_loc(2,18)= -1.0
      te_loc(3,18)= 0.0
c  ! +------------------------------------------------------------------------------+
c  ! |                     The 6 n, m, l vectors for the warp yarns                 |
c  ! |                     rotate fill yarn by 90 deg about Z axis
c  ! +------------------------------------------------------------------------------+

      do k=1,18
         te_loc(1,k+18)=-te_loc(2,k)
         te_loc(2,k+18)=te_loc(1,k)
         te_loc(3,k+18)=te_loc(3,k)
      end do
            
c  ! +------------------------------------------------------------------------------+

      do k=1,36
          te(1,k)=te_loc(1,k)
          te(2,k)=te_loc(2,k)
          te(3,k)=te_loc(3,k)
      end do          
c  ! +------------------------------------------------------------------------------+
c  ! |                   Generate N_{ij}=n_i n_j  using Voight rule                 |
c  ! +------------------------------------------------------------------------------+
      do k=1,6
         do jp=1,np
            xfyNij(k,jp) = 0.0
            xfyMij(k,jp) = 0.0
            xfyLij(k,jp) = 0.0
            xfyAij(k,jp) = 0.0
            xfyBij(k,jp) = 0.0
            xfyCij(k,jp) = 0.0

            xwyNij(k,jp) = 0.0
            xwyMij(k,jp) = 0.0
            xwyLij(k,jp) = 0.0
            xwyAij(k,jp) = 0.0
            xwyBij(k,jp) = 0.0
            xwyCij(k,jp) = 0.0
         end do
      end do
      do jp=1, np
         xn1(1) = te(1,jp) 
         xn1(2) = te(2,jp) 
         xn1(3) = te(3,jp) 
         xm1(1) = te(1,jp+6)
         xm1(2) = te(2,jp+6)
         xm1(3) = te(3,jp+6)
         xl1(1) = te(1,jp+12)
         xl1(2) = te(2,jp+12)
         xl1(3) = te(3,jp+12)

         xn2(1) = te(1,jp+18)
         xn2(2) = te(2,jp+18)
         xn2(3) = te(3,jp+18)
         xm2(1) = te(1,jp+24)
         xm2(2) = te(2,jp+24)
         xm2(3) = te(3,jp+24)
         xl2(1) = te(1,jp+30)
         xl2(2) = te(2,jp+30)
         xl2(3) = te(3,jp+30)
         do k=1,6 
            i=ij(1,k) 
            j=ij(2,k) 
            xfyNij(k,jp) = xn1(i)*xn1(j) 
            xfyMij(k,jp) = xm1(i)*xm1(j)
            xfyLij(k,jp) = xl1(i)*xl1(j)
            xwyNij(k,jp) = xn2(i)*xn2(j)
            xwyMij(k,jp) = xm2(i)*xm2(j)
            xwyLij(k,jp) = xl2(i)*xl2(j)
            if(k>3)then
             xfyNij(k,jp)=xfyNij(k,jp)*2.d0
             xfyMij(k,jp)=xfyMij(k,jp)*2.d0
             xfyLij(k,jp)=xfyLij(k,jp)*2.d0
             xwyNij(k,jp)=xwyNij(k,jp)*2.d0
             xwyMij(k,jp)=xwyMij(k,jp)*2.d0
             xwyLij(k,jp)=xwyLij(k,jp)*2.d0
            end if

            xfyAij(k,jp) = 0.5d0*(xn1(i)*xl1(j) + xl1(i)*xn1(j))
            xfyBij(k,jp) = 0.5d0*(xm1(i)*xl1(j) + xl1(i)*xm1(j)) 
            xfyCij(k,jp) = 0.5d0*(xn1(i)*xm1(j) + xm1(i)*xn1(j))
            xwyAij(k,jp) = 0.5d0*(xn2(i)*xl2(j) + xl2(i)*xn2(j))
            xwyBij(k,jp) = 0.5d0*(xm2(i)*xl2(j) + xl2(i)*xm2(j))
            xwyCij(k,jp) = 0.5d0*(xn2(i)*xm2(j) + xm2(i)*xn2(j))
            if(k>3)then
               xfyAij(k,jp)=xfyAij(k,jp)*2.d0
               xfyBij(k,jp)=xfyBij(k,jp)*2.d0
               xfyCij(k,jp)=xfyCij(k,jp)*2.d0
               xwyAij(k,jp)=xwyAij(k,jp)*2.d0
               xwyBij(k,jp)=xwyBij(k,jp)*2.d0
               xwyCij(k,jp)=xwyCij(k,jp)*2.d0
            end if
         end do
      end do
c  ! +------------------------------------------------------------------------------+
c  ! |           Generate N_{ij} N_{kl}=n_i n_j n_k n_l  using Voigt rule           |
c  ! +------------------------------------------------------------------------------+
      do jp=1, np
         do i=1,6
            do j=1,6
               xfyNijNkl(i,j,jp) = xfyNij(i,jp)*xfyNij(j,jp)
               xfyMijMkl(i,j,jp) = xfyMij(i,jp)*xfyMij(j,jp)
               xfyLijLkl(i,j,jp) = xfyLij(i,jp)*xfyLij(j,jp)

               xfyNijMkl(i,j,jp) = 0.5d0*(xfyNij(i,jp)*xfyMij(j,jp) +
     $xfyMij(i,jp)*xfyNij(j,jp))
               xfyNijLkl(i,j,jp) = 0.5d0*(xfyNij(i,jp)*xfyLij(j,jp) +
     $xfyLij(i,jp)*xfyNij(j,jp))
               xfyMijLkl(i,j,jp) = 0.5d0*(xfyMij(i,jp)*xfyLij(j,jp) +
     $xfyLij(i,jp)*xfyMij(j,jp))

               xfyAijAkl(i,j,jp) = xfyAij(i,jp)*xfyAij(j,jp)
               xfyBijBkl(i,j,jp) = xfyBij(i,jp)*xfyBij(j,jp)
               xfyCijCkl(i,j,jp) = xfyCij(i,jp)*xfyCij(j,jp)

               xwyNijNkl(i,j,jp) = xwyNij(i,jp)*xwyNij(j,jp)
               xwyMijMkl(i,j,jp) = xwyMij(i,jp)*xwyMij(j,jp)
               xwyLijLkl(i,j,jp) = xwyLij(i,jp)*xwyLij(j,jp)

               xwyNijMkl(i,j,jp) = 0.5d0*(xwyNij(i,jp)*xwyMij(j,jp) + 
     $xwyMij(i,jp)*xwyNij(j,jp))
               xwyNijLkl(i,j,jp) = 0.5d0*(xwyNij(i,jp)*xwyLij(j,jp) +
     $xwyLij(i,jp)*xwyNij(j,jp))
               xwyMijLkl(i,j,jp) = 0.5d0*(xwyMij(i,jp)*xwyLij(j,jp) +
     $xwyLij(i,jp)*xwyMij(j,jp))

               xwyAijAkl(i,j,jp) = xwyAij(i,jp)*xwyAij(j,jp)
               xwyBijBkl(i,j,jp) = xwyBij(i,jp)*xwyBij(j,jp)
               xwyCijCkl(i,j,jp) = xwyCij(i,jp)*xwyCij(j,jp)
            end do
         end do
      end do

c  ! +-----------------------------------------------+
c  ! | xfyNijNkl is N_{ij}N_{kl} for the fill yarn   |
c  ! | xfyMijMkl is M_{ij}M_{kl} for the fill yarn   |
c  ! | xfyLijLkl is L_{ij}L_{kl} for the fill yarn   |
c  ! | xfyNijNkl is N_{ij}N_{kl} for the warp yarn   |
c  ! | xwyMijMkl is M_{ij}M_{kl} for the warp yarn   |
c  ! | xwyLijLkl is L_{ij}L_{kl} for the warp yarn   |
c  ! +-----------------------------------------------+

      end 
C +--------------------------------------------------------------------+
C |                        SUBROUTINE KELAST_ISOTROPIC                 |
C +--------------------------------------------------------------------+
      subroutine kelast_isotropic(wdm, wdm2, e, E_vis, eee)
   
      include 'vaba_param.inc'
      common /matrix_prop/Em0, xnum0, Gm0, eta, xeta
      dimension e(6,6)
      
      eee=(1.d0-wdm)*(1.d0-wdm2)*E_vis
!      eee=(1.d0-wdm)*Em0

      xnu=xnum0
      e = 0.0d0
C Fourth order elastic stiffness tensor
      e(1,1)=eee*(1.0D0-xnu)/((1.0D0+xnu)*(1.0D0-2.0D0*xnu))
      e(1,2)=eee*xnu/((1.0D0+xnu)*(1.0D0-2.0D0*xnu))
      e(1,3)=e(1,2)
      e(2,1)=e(1,2)
      e(2,2)=e(1,1)
      e(2,3)=e(1,2)
      e(3,1)=e(1,2)
      e(3,2)=e(1,2)
      e(3,3)=e(1,1)
      e(4,4)=eee/(2.0D0+2.0D0*xnu)
      e(5,5)=e(4,4)
      e(6,6)=e(4,4)

      return
      end 
C +--------------------------------------------------------------------+
C |                SUBROUTINE KELAST_YARNPLATE                         |
C +--------------------------------------------------------------------+
      subroutine kelast_yarnplate(wdm, wdm2, wdty, wdcy, iplate, stiff_yp)

      include 'vaba_param.inc'
      common /matrix_prop/Em0, xnum0, Gm0, eta, xeta
      common /fiber_prop/Eaxf0, Elatf0, G12f0, G23f0, xnu12f0, xnu23f0
      common /matparam/xR1, xR2, xq, xpar, epsn1_f, epsn0_fc, xH, xn
      common /common_vars/ PI, alpha, vf_my, vf_y, vf_m, vf
      common /projection_vars/ xfyNijNkl, xfyMijMkl, xfyLijLkl,
     $xfyNijMkl, xfyNijLkl, xfyMijLkl,
     $xfyAijAkl, xfyBijBkl, xfyCijCkl,
     $xwyNijNkl, xwyMijMkl, xwyLijLkl,
     $xwyNijMkl, xwyNijLkl, xwyMijLkl,
     $xwyAijAkl, xwyBijBkl, xwyCijCkl,
     $xfyNij, xfyMij, xfyLij, xfyAij, xfyBij, xfyCij,
     $xwyNij, xwyMij, xwyLij, xwyAij, xwyBij, xwyCij

      integer :: i, j, k, np, iplate
      parameter (np=6)
      dimension xfyNij(6,6), xfyMij(6,6), xfyLij(6,6)
      dimension xfyAij(6,6), xfyBij(6,6), xfyCij(6,6)

      dimension xwyNij(6,6), xwyMij(6,6), xwyLij(6,6)
      dimension xwyAij(6,6), xwyBij(6,6), xwyCij(6,6)

      dimension xfyNijNkl(6,6,6), xfyMijMkl(6,6,6), xfyLijLkl(6,6,6)
      dimension xfyNijMkl(6,6,6), xfyNijLkl(6,6,6), xfyMijLkl(6,6,6)
      dimension xfyAijAkl(6,6,6), xfyBijBkl(6,6,6), xfyCijCkl(6,6,6)

      dimension xwyNijNkl(6,6,6), xwyMijMkl(6,6,6), xwyLijLkl(6,6,6)
      dimension xwyNijMkl(6,6,6), xwyNijLkl(6,6,6), xwyMijLkl(6,6,6)
      dimension xwyAijAkl(6,6,6), xwyBijBkl(6,6,6), xwyCijCkl(6,6,6)

      dimension xNijNkl(6,6), xMijMkl(6,6), xLijLkl(6,6)
      dimension xNijMkl(6,6), xNijLkl(6,6), xMijLkl(6,6)
      dimension xAijAkl(6,6), xBijBkl(6,6), xCijCkl(6,6)

      dimension stiff_yp(6,6), stiff_f(6,6)
      dimension stiff_1(6,6), stiff_2(6,6), stiff_3(6,6)
      dimension stiff_4(6,6), stiff_5(6,6), stiff_6(6,6)
      dimension stiff_7(6,6), stiff_8(6,6), stiff_9(6,6)
      dimension stiff_sum(6,6)
      dimension AA(3,3)
      dimension wdny(np), wdmy(np), wdly(np)
      dimension wday(np), wdby(np), wdcy(np)
      dimension wdmly(np), wdacy(np)
      dimension wdty(np)
      dimension wdy(np)
      dimension xmpwt(np)

      PI=3.1415926535897932
      stiff_f=0
      xmpwt(1)=0.16666d0
      xmpwt(2)=0.16666d0
      xmpwt(3)=0.16666d0
      xmpwt(4)=0.16666d0
      xmpwt(5)=0.16666d0
      xmpwt(6)=0.16666d0
      do jp=1,np
	    wdcy(jp)=0.d0
        Eaxf=(1.d0-wdty(jp))*(1.d0-wdcy(jp))*Eaxf0
        Elatf=(1.d0-wdty(jp))*(1.d0-wdcy(jp))*Elatf0
C        wdacy(jp)=max(wdy(jp),wdcy(jp))
        G12f=(1.d0-wdty(jp))*(1.d0-wdcy(jp))*G12f0
        G23f=(1.d0-wdty(jp))*(1.d0-wdcy(jp))*G23f0

c        Eaxf=(1.d0-wdy(jp))*Eaxf0
c        Elatf=(1.d0-wdy(jp))*Elatf0
c        G12f=(1.d0-wdy(jp))*G12f0
c        G23f=(1.d0-wdy(jp))*G23f0

        Em=(1.d0-wdm)*(1.d0-wdm2)*Em0
        Gm=Em/(2*(1+xnum0))

        Eaxy=(vf_y*Eaxf)+(vf_my*Em);                              ! Eqn 17
c        Elaty = (vf_y*Elatf)+(vf_my*Em);
        Elaty=1/((vf_y/Elatf)+(vf_my/Em));                        ! Eqn 18
c        G12y = 1/((vf_y/G12f)+(vf_my/Gm));
        G12y = Gm * ((G12f+Gm)+vf_y*(G12f-Gm))/((G12f+Gm)-
     $vf_y*(G12f-Gm))                                             ! Eqn 19
        G23y = 1/((vf_y/G23f)+(vf_my/Gm));                        ! Eqn 20

        xnu12y =(vf_y*xnu12f0)+(vf_my*xnum0);                     ! Eqn 21
        xnu23y =(vf_y*xnu23f0)+(vf_my*xnum0);
        xnu13y =(vf_y*xnu12f0)+(vf_my*xnum0);

        G13y = G12y;
        xnu21y = xnu12y*Elaty/Eaxy;
        xnu31y = xnu13y*Elaty/Eaxy;
        xnu32y = xnu23y;

c     Define Cij
        AA(1,1) = 1.d0
        AA(1,2) = -1*xnu21y
        AA(1,3) = -1*xnu31y
        AA(2,1) =-1*xnu12y 
        AA(2,2) = 1.d0
        AA(2,3) = -1*xnu32y
        AA(3,1) = -1*xnu13y 
        AA(3,2) = -1*xnu23y
        AA(3,3) = 1

        DETAA = AA(1,1)*(AA(2,2)*AA(3,3)-AA(2,3)*AA(3,2))-
     $(AA(1,2)*(AA(2,1)*AA(3,3)-AA(3,1)*AA(2,3)))+
     $(AA(1,3)*(AA(2,1)*AA(3,2)-AA(3,1)*AA(2,2)))                 ! Eqn 26

        DELTA=DETAA/(Eaxy*Elaty*Elaty);                           ! Eqn 26
        
C Stiffness matrix in eqn 22,23
        C11 = (1-(xnu23y*xnu32y))/(Elaty*Elaty*DELTA);
        C22 = (1-(xnu13y*xnu31y))/(Eaxy*Elaty*DELTA);
        C33 = (1-(xnu12y*xnu21y))/(Eaxy*Elaty*DELTA);

        C12 = (xnu12y + (xnu13y*xnu32y))/(DELTA*Eaxy*Elaty);
        C13 = (xnu23y + (xnu21y*xnu13y))/(DELTA*Eaxy*Elaty);
        C23 = (xnu31y + (xnu21y*xnu32y))/(DELTA*Elaty*Elaty);
  
        C44 = G12y;
        C55 = G23y;
        C66 = G13y;

C   Define coeff of potential energy Eqn(27-34)

        D = C11 + C22 + C33 + 2*C12 + 2*C13 + 2*C23;          ! Eqn 26
        A1 = (1/(2*D)) *((C11+C12+C13)**2+
     $(C11*(C22+2*C23+C33)-(C12+C13)**2));
        A2 = (1/(2*D)) *((C12+C22+C23)**2+
     $(C22*(C11+2*C13+C33)-(C12+C23)**2));
        A3 = (1/(2*D)) *((C13+C23+C33)**2+
     $(C33*(C11+2*C12+C22)-(C13+C23)**2));

        A4 = (1/(2*D))*(2*(C11+C12+C13)*(C12+C22+C23)+
     $2*(C12**2-(C11+C13)*(C22+C23)+C12*(C13+C23+C33)));
        A5 = (1/(2*D))*(2*(C12+C22+C23)*(C13+C23+C33)+
     $2*(C23*(C11+C23)+C13*(-C22+C23)-C22*C33-C12*(C13-C23+C33)));
        A6 = (1/(2*D))*(2*(C13+C23+C33)*(C11+C12+C13)+
     $2*(C13*(C13+C22+C23)-C11*(C23+C33)-C12*(-C13+C23+C33)));

        A7 = (1/(2*D))*(D*C44);
        A8 = (1/(2*D))*(D*C55);
        A9 = (1/(2*D))*(D*C66);
  
        xNijNkl=0.d0
        xMijMkl=0.d0
        xLijLkl=0.d0
        xNijMkl=0.d0
        xNijLkl=0.d0
        xMijLkl=0.d0
        xAijAkl=0.d0
        xBijBkl=0.d0
        xCijCkl=0.d0
        stiff_1=0.d0
        stiff_2=0.d0
        stiff_3=0.d0
        stiff_4=0.d0
        stiff_5=0.d0
        stiff_6=0.d0
        stiff_7=0.d0
        stiff_8=0.d0
        stiff_9=0.d0
        stiff_sum=0.d0
        do k=1,6
          do j=1,6
             if(iplate==1) then
               xNijNkl(k,j)=xfyNijNkl(k,j,jp)
               xMijMkl(k,j)=xfyMijMkl(k,j,jp)
               xLijLkl(k,j)=xfyLijLkl(k,j,jp)
               xNijMkl(k,j)=xfyNijMkl(k,j,jp)
               xNijLkl(k,j)=xfyNijLkl(k,j,jp)
               xMijLkl(k,j)=xfyMijLkl(k,j,jp)
               xAijAkl(k,j)=xfyAijAkl(k,j,jp)
               xBijBkl(k,j)=xfyBijBkl(k,j,jp)
               xCijCkl(k,j)=xfyCijCkl(k,j,jp)
             else
               xNijNkl(k,j)=xwyNijNkl(k,j,jp)
               xMijMkl(k,j)=xwyMijMkl(k,j,jp)
               xLijLkl(k,j)=xwyLijLkl(k,j,jp)
               xNijMkl(k,j)=xwyNijMkl(k,j,jp)
               xNijLkl(k,j)=xwyNijLkl(k,j,jp)
               xMijLkl(k,j)=xwyMijLkl(k,j,jp)
               xAijAkl(k,j)=xwyAijAkl(k,j,jp)
               xBijBkl(k,j)=xwyBijBkl(k,j,jp)
               xCijCkl(k,j)=xwyCijCkl(k,j,jp)
             end if 
             stiff_1(k,j) = 2 * A1 * (xNijNkl(k,j))
             stiff_2(k,j) = 2 * A2 * (xMijMkl(k,j))
             stiff_3(k,j) = 2 * A3 * (xLijLkl(k,j))
             stiff_4(k,j) = A4 * xNijLkl(k,j)
             stiff_5(k,j) = A5 * xMijLkl(k,j)
             stiff_6(k,j) = A6 * xNijMkl(k,j)
             stiff_7(k,j) = 2 * A7 * (xAijAkl(k,j))
             stiff_8(k,j) = 2 * A8 * (xBijBkl(k,j))
             stiff_9(k,j) = 2 * A9 * (xCijCkl(k,j))

             stiff_sum(k,j) = xmpwt(jp)*(stiff_1(k,j)+stiff_2(k,j)+
     $stiff_3(k,j)+stiff_4(k,j)+stiff_5(k,j)+stiff_6(k,j)+
     $stiff_7(k,j)+stiff_8(k,j)+stiff_9(k,j))

            stiff_f(k,j)=stiff_f(k,j)+stiff_sum(k,j)
           end do
        end do
      end do

      stiff_yp=stiff_f
      end 
C +--------------------------------------------------------------------+
C |                  SUBROUTINE getRotMat                             |
C +--------------------------------------------------------------------+
      subroutine getRotMat(theta1, RR1, RR2, xInvRR1)
      include 'vaba_param.inc'
      dimension RR1(6,6), RR2(6,6), xInvRR1(6,6)

      xc=cos(theta1)
      xs=sin(theta1)

c     RR1 = Rotation matrix for fourth order tensor in Voigt form(6x6)

         do k=1,6
            do j=1,6
               RR1(k,j)=0.d0
c               xInvRR1(k,j)=0.d0
            end do
         end do
         RR1(1,1)=xc*xc
         RR1(1,2)=xs*xs
         RR1(1,4)=2*xc*xs
         RR1(2,1)=RR1(1,2)
         RR1(2,2)=RR1(1,1)
         RR1(2,4)=-RR1(1,4)
         RR1(3,3)=1.d0
         RR1(6,6)=xc
         RR1(6,5)=-xs
         RR1(5,6)=xs
         RR1(5,5)=xc
         RR1(4,1)=-xc*xs
         RR1(4,2)=xc*xs
         RR1(4,4)=(xc*xc)-(xs*xs)

         call M66INV(RR1, xInvRR1)

         do k=1,6
            do j=1,6
               RR2(k,j)=RR1(k,j)
            end do
         end do
         RR2(1,4)=xc*xs
         RR2(2,4)=-xc*xs
         RR2(4,1)=-2*xc*xs
         RR2(4,2)=2*xc*xs

c  ! Now, to rotate the stiffness matrix do: 
c  ! 1. Generate the regular rotation matrix that corresponds to the desired rotation;
c  ! 2. call this subroutine: call get_RotMat_for_2ndOrderTensor_in_Voight_form(R,RR,xInvRR)
c  ! 3. calculate stiff_rotated = matmul(RR,matmul(stiff,xInvRR))

      end
C +--------------------------------------------------------------------+
C |                  SUBROUTINE voigt_to_tensor                        |
C +--------------------------------------------------------------------+
      subroutine voigt_to_tensor(ABC, ABC33)
      include 'vaba_param.inc'
      dimension ABC(6)
      dimension ABC33(3,3)

      ABC33(1,1)=ABC(1)
      ABC33(2,2)=ABC(2)
      ABC33(3,3)=ABC(3)
      ABC33(1,2)=0.5d0*ABC(4)
      ABC33(2,1)=ABC33(1,2)
      ABC33(2,3)=0.5d0*ABC(5)
      ABC33(3,2)=ABC33(2,3)
      ABC33(1,3)=0.5d0*ABC(6)
      ABC33(3,1)=ABC33(1,3)

      end
C +--------------------------------------------------------------------+
C |                  SUBROUTINE tensor_to_voigt                        |
C +--------------------------------------------------------------------+
      subroutine tensor_to_voigt(ABC33, ABC)
      include 'vaba_param.inc'
      dimension ABC(6)
      dimension ABC33(3,3)

      ABC(1)=ABC33(1,1)
      ABC(2)=ABC33(2,2)
      ABC(3)=ABC33(3,3)
      ABC(4)=2*ABC33(1,2)
      ABC(5)=2*ABC33(2,3)
      ABC(6)=2*ABC33(1,3)

      end
C +--------------------------------------------------------------------+
C |                  SUBROUTINE M66INV                                 |
C +--------------------------------------------------------------------+
      SUBROUTINE M66INV (A, AINV)
!***********************************************************************************************************************************
! M66INV - Compute the inverse of a 6x6 matrix.
! A       = input 6x6 matrix to be inverted
! AINV  = output 6x6 inverse of matrix A
!***********************************************************************************************************************************
      include 'vaba_param.inc'

      DIMENSION A(6,6), AINV(6,6) 
      DIMENSION COFACTOR(6,6)
C      LOGICAL, INTENT(OUT) :: OK_FLAG
c      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10

      A11=A(1,1); A12=A(1,2); A13=A(1,3)
      A14=A(1,4); A15=A(1,5); A16=A(1,6)
      A21=A(2,1); A22=A(2,2); A23=A(2,3) 
      A24=A(2,4); A25=A(2,5); A26=A(2,6)
      A31=A(3,1); A32=A(3,2); A33=A(3,3)
      A34=A(3,4); A35=A(3,5); A36=A(3,6)
      A41=A(4,1); A42=A(4,2); A43=A(4,3)
      A44=A(4,4); A45=A(4,5); A46=A(4,6)
      A51=A(5,1); A52=A(5,2); A53=A(5,3)
      A54=A(5,4); A55=A(5,5); A56=A(5,6)
      A61=A(6,1); A62=A(6,2); A63=A(6,3)
      A64=A(6,4); A65=A(6,5); A66=A(6,6)

      DET = -(A16*A25*A34*A43*A52-A15*A26*A34*A43*A52-A16*A24*A35*A43*
     $A52+A14*A26*A35*A43*A52+A15*A24*A36*A43*A52-A14*A25*A36*A43*A52-
     $A16*A25*A33*A44*A52+A15*A26*A33*A44*A52+A16*A23*A35*A44*A52-A13*
     $A26*A35*A44*A52-A15*A23*A36*A44*A52+A13*A25*A36*A44*A52+A16*A24*
     $A33*A45*A52-A14*A26*A33*A45*A52-A16*A23*A34*A45*A52+A13*A26*A34*
     $A45*A52+A14*A23*A36*A45*A52-A13*A24*A36*A45*A52-A15*A24*A33*A46*
     $A52+A14*A25*A33*A46*A52+A15*A23*A34*A46*A52-A13*A25*A34*A46*A52-
     $A14*A23*A35*A46*A52+A13*A24*A35*A46*A52-A16*A25*A34*A42*A53+A15*
     $A26*A34*A42*A53+A16*A24*A35*A42*A53-A14*A26*A35*A42*A53-A15*A24*
     $A36*A42*A53+A14*A25*A36*A42*A53+A16*A25*A32*A44*A53-A15*A26*A32*
     $A44*A53-A16*A22*A35*A44*A53+A12*A26*A35*A44*A53+A15*A22*A36*A44*
     $A53-A12*A25*A36*A44*A53-A16*A24*A32*A45*A53+A14*A26*A32*A45*A53+
     $A16*A22*A34*A45*A53-A12*A26*A34*A45*A53-A14*A22*A36*A45*A53+A12*
     $A24*A36*A45*A53+A15*A24*A32*A46*A53-A14*A25*A32*A46*A53-A15*A22*
     $A34*A46*A53+A12*A25*A34*A46*A53+A14*A22*A35*A46*A53-A12*A24*A35*
     $A46*A53+A16*A25*A33*A42*A54-A15*A26*A33*A42*A54-A16*A23*A35*A42*
     $A54+A13*A26*A35*A42*A54+A15*A23*A36*A42*A54-A13*A25*A36*A42*A54-
     $A16*A25*A32*A43*A54+A15*A26*A32*A43*A54+A16*A22*A35*A43*A54-A12*
     $A26*A35*A43*A54-A15*A22*A36*A43*A54+A12*A25*A36*A43*A54+A16*A23*
     $A32*A45*A54-A13*A26*A32*A45*A54-A16*A22*A33*A45*A54+A12*A26*A33*
     $A45*A54+A13*A22*A36*A45*A54-A12*A23*A36*A45*A54-A15*A23*A32*A46*
     $A54+A13*A25*A32*A46*A54+A15*A22*A33*A46*A54-A12*A25*A33*A46*A54-
     $A13*A22*A35*A46*A54+A12*A23*A35*A46*A54-A16*A24*A33*A42*A55+A14*
     $A26*A33*A42*A55+A16*A23*A34*A42*A55-A13*A26*A34*A42*A55-A14*A23*
     $A36*A42*A55+A13*A24*A36*A42*A55+A16*A24*A32*A43*A55-A14*A26*A32*
     $A43*A55-A16*A22*A34*A43*A55+A12*A26*A34*A43*A55+A14*A22*A36*A43*
     $A55-A12*A24*A36*A43*A55-A16*A23*A32*A44*A55+A13*A26*A32*A44*A55+
     $A16*A22*A33*A44*A55-A12*A26*A33*A44*A55-A13*A22*A36*A44*A55+A12*
     $A23*A36*A44*A55+A14*A23*A32*A46*A55-A13*A24*A32*A46*A55-A14*A22*
     $A33*A46*A55+A12*A24*A33*A46*A55+A13*A22*A34*A46*A55-A12*A23*A34*
     $A46*A55+A15*A24*A33*A42*A56-A14*A25*A33*A42*A56-A15*A23*A34*A42*
     $A56+A13*A25*A34*A42*A56+A14*A23*A35*A42*A56-A13*A24*A35*A42*A56-
     $A15*A24*A32*A43*A56+A14*A25*A32*A43*A56+A15*A22*A34*A43*A56-A12*
     $A25*A34*A43*A56-A14*A22*A35*A43*A56+A12*A24*A35*A43*A56+A15*A23*
     $A32*A44*A56-A13*A25*A32*A44*A56-A15*A22*A33*A44*A56+A12*A25*A33*
     $A44*A56+A13*A22*A35*A44*A56-A12*A23*A35*A44*A56-A14*A23*A32*A45*
     $A56+A13*A24*A32*A45*A56+A14*A22*A33*A45*A56-A12*A24*A33*A45*A56-
     $A13*A22*A34*A45*A56+A12*A23*A34*A45*A56)*A61+(A16*A25*A34*A43*A51-
     $A15*A26*A34*A43*A51-A16*A24*A35*A43*A51+A14*A26*A35*A43*A51+A15*
     $A24*A36*A43*A51-A14*A25*A36*A43*A51-A16*A25*A33*A44*A51+A15*A26*
     $A33*A44*A51+A16*A23*A35*A44*A51-A13*A26*A35*A44*A51-A15*A23*A36*
     $A44*A51+A13*A25*A36*A44*A51+A16*A24*A33*A45*A51-A14*A26*A33*A45*
     $A51-A16*A23*A34*A45*A51+A13*A26*A34*A45*A51+A14*A23*A36*A45*A51-
     $A13*A24*A36*A45*A51-A15*A24*A33*A46*A51+A14*A25*A33*A46*A51+A15*
     $A23*A34*A46*A51-A13*A25*A34*A46*A51-A14*A23*A35*A46*A51+A13*A24*
     $A35*A46*A51-A16*A25*A34*A41*A53+A15*A26*A34*A41*A53+A16*A24*A35*
     $A41*A53-A14*A26*A35*A41*A53-A15*A24*A36*A41*A53+A14*A25*A36*A41*
     $A53+A16*A25*A31*A44*A53-A15*A26*A31*A44*A53-A16*A21*A35*A44*A53+
     $A11*A26*A35*A44*A53+A15*A21*A36*A44*A53-A11*A25*A36*A44*A53-A16*
     $A24*A31*A45*A53+A14*A26*A31*A45*A53+A16*A21*A34*A45*A53-A11*A26*
     $A34*A45*A53-A14*A21*A36*A45*A53+A11*A24*A36*A45*A53+A15*A24*A31*
     $A46*A53-A14*A25*A31*A46*A53-A15*A21*A34*A46*A53+A11*A25*A34*A46*
     $A53+A14*A21*A35*A46*A53-A11*A24*A35*A46*A53+A16*A25*A33*A41*A54-
     $A15*A26*A33*A41*A54-A16*A23*A35*A41*A54+A13*A26*A35*A41*A54+A15*
     $A23*A36*A41*A54-A13*A25*A36*A41*A54-A16*A25*A31*A43*A54+A15*A26*
     $A31*A43*A54+A16*A21*A35*A43*A54-A11*A26*A35*A43*A54-A15*A21*A36*
     $A43*A54+A11*A25*A36*A43*A54+A16*A23*A31*A45*A54-A13*A26*A31*A45*
     $A54-A16*A21*A33*A45*A54+A11*A26*A33*A45*A54+A13*A21*A36*A45*A54-
     $A11*A23*A36*A45*A54-A15*A23*A31*A46*A54+A13*A25*A31*A46*A54+A15*
     $A21*A33*A46*A54-A11*A25*A33*A46*A54-A13*A21*A35*A46*A54+A11*A23*
     $A35*A46*A54-A16*A24*A33*A41*A55+A14*A26*A33*A41*A55+A16*A23*A34*
     $A41*A55-A13*A26*A34*A41*A55-A14*A23*A36*A41*A55+A13*A24*A36*A41*
     $A55+A16*A24*A31*A43*A55-A14*A26*A31*A43*A55-A16*A21*A34*A43*A55+
     $A11*A26*A34*A43*A55+A14*A21*A36*A43*A55-A11*A24*A36*A43*A55-A16*
     $A23*A31*A44*A55+A13*A26*A31*A44*A55+A16*A21*A33*A44*A55-A11*A26*
     $A33*A44*A55-A13*A21*A36*A44*A55+A11*A23*A36*A44*A55+A14*A23*A31*
     $A46*A55-A13*A24*A31*A46*A55-A14*A21*A33*A46*A55+A11*A24*A33*A46*
     $A55+A13*A21*A34*A46*A55-A11*A23*A34*A46*A55+A15*A24*A33*A41*A56-
     $A14*A25*A33*A41*A56-A15*A23*A34*A41*A56+A13*A25*A34*A41*A56+A14*
     $A23*A35*A41*A56-A13*A24*A35*A41*A56-A15*A24*A31*A43*A56+A14*A25*
     $A31*A43*A56+A15*A21*A34*A43*A56-A11*A25*A34*A43*A56-A14*A21*A35*
     $A43*A56+A11*A24*A35*A43*A56+A15*A23*A31*A44*A56-A13*A25*A31*A44*
     $A56-A15*A21*A33*A44*A56+A11*A25*A33*A44*A56+A13*A21*A35*A44*A56-
     $A11*A23*A35*A44*A56-A14*A23*A31*A45*A56+A13*A24*A31*A45*A56+A14*
     $A21*A33*A45*A56-A11*A24*A33*A45*A56-A13*A21*A34*A45*A56+A11*A23*
     $A34*A45*A56)*A62-(A16*A25*A34*A42*A51-A15*A26*A34*A42*A51-A16*A24*
     $A35*A42*A51+A14*A26*A35*A42*A51+A15*A24*A36*A42*A51-A14*A25*A36*
     $A42*A51-A16*A25*A32*A44*A51+A15*A26*A32*A44*A51+A16*A22*A35*A44*
     $A51-A12*A26*A35*A44*A51-A15*A22*A36*A44*A51+A12*A25*A36*A44*A51+
     $A16*A24*A32*A45*A51-A14*A26*A32*A45*A51-A16*A22*A34*A45*A51+A12*
     $A26*A34*A45*A51+A14*A22*A36*A45*A51-A12*A24*A36*A45*A51-A15*A24*
     $A32*A46*A51+A14*A25*A32*A46*A51+A15*A22*A34*A46*A51-A12*A25*A34*
     $A46*A51-A14*A22*A35*A46*A51+A12*A24*A35*A46*A51-A16*A25*A34*A41*
     $A52+A15*A26*A34*A41*A52+A16*A24*A35*A41*A52-A14*A26*A35*A41*A52-
     $A15*A24*A36*A41*A52+A14*A25*A36*A41*A52+A16*A25*A31*A44*A52-A15
     $*A26*A31*A44*A52-A16*A21*A35*A44*A52+A11*A26*A35*A44*A52+A15*A21*
     $A36*A44*A52-A11*A25*A36*A44*A52-A16*A24*A31*A45*A52+A14*A26*A31*
     $A45*A52+A16*A21*A34*A45*A52-A11*A26*A34*A45*A52-A14*A21*A36*A45*
     $A52+A11*A24*A36*A45*A52+A15*A24*A31*A46*A52-A14*A25*A31*A46*A52-
     $A15*A21*A34*A46*A52+A11*A25*A34*A46*A52+A14*A21*A35*A46*A52-A11*
     $A24*A35*A46*A52+A16*A25*A32*A41*A54-A15*A26*A32*A41*A54-A16*A22*
     $A35*A41*A54+A12*A26*A35*A41*A54+A15*A22*A36*A41*A54-A12*A25*A36*
     $A41*A54-A16*A25*A31*A42*A54+A15*A26*A31*A42*A54+A16*A21*A35*A42*
     $A54-A11*A26*A35*A42*A54-A15*A21*A36*A42*A54+A11*A25*A36*A42*A54+
     $A16*A22*A31*A45*A54-A12*A26*A31*A45*A54-A16*A21*A32*A45*A54+A11*
     $A26*A32*A45*A54+A12*A21*A36*A45*A54-A11*A22*A36*A45*A54-A15*A22*
     $A31*A46*A54+A12*A25*A31*A46*A54+A15*A21*A32*A46*A54-A11*A25*A32*
     $A46*A54-A12*A21*A35*A46*A54+A11*A22*A35*A46*A54-A16*A24*A32*A41*
     $A55+A14*A26*A32*A41*A55+A16*A22*A34*A41*A55-A12*A26*A34*A41*A55-
     $A14*A22*A36*A41*A55+A12*A24*A36*A41*A55+A16*A24*A31*A42*A55-A14*
     $A26*A31*A42*A55-A16*A21*A34*A42*A55+A11*A26*A34*A42*A55+A14*A21*
     $A36*A42*A55-A11*A24*A36*A42*A55-A16*A22*A31*A44*A55+A12*A26*A31*
     $A44*A55+A16*A21*A32*A44*A55-A11*A26*A32*A44*A55-A12*A21*A36*A44*
     $A55+A11*A22*A36*A44*A55+A14*A22*A31*A46*A55-A12*A24*A31*A46*A55-
     $A14*A21*A32*A46*A55+A11*A24*A32*A46*A55+A12*A21*A34*A46*A55-A11*
     $A22*A34*A46*A55+A15*A24*A32*A41*A56-A14*A25*A32*A41*A56-A15*A22*
     $A34*A41*A56+A12*A25*A34*A41*A56+A14*A22*A35*A41*A56-A12*A24*A35*
     $A41*A56-A15*A24*A31*A42*A56+A14*A25*A31*A42*A56+A15*A21*A34*A42*
     $A56-A11*A25*A34*A42*A56-A14*A21*A35*A42*A56+A11*A24*A35*A42*A56+
     $A15*A22*A31*A44*A56-A12*A25*A31*A44*A56-A15*A21*A32*A44*A56+A11*
     $A25*A32*A44*A56+A12*A21*A35*A44*A56-A11*A22*A35*A44*A56-A14*A22*
     $A31*A45*A56+A12*A24*A31*A45*A56+A14*A21*A32*A45*A56-A11*A24*A32*
     $A45*A56-A12*A21*A34*A45*A56+A11*A22*A34*A45*A56)*A63+(A16*A25*A33*
     $A42*A51-A15*A26*A33*A42*A51-A16*A23*A35*A42*A51+A13*A26*A35*A42*
     $A51+A15*A23*A36*A42*A51-A13*A25*A36*A42*A51-A16*A25*A32*A43*A51+
     $A15*A26*A32*A43*A51+A16*A22*A35*A43*A51-A12*A26*A35*A43*A51-A15*
     $A22*A36*A43*A51+A12*A25*A36*A43*A51+A16*A23*A32*A45*A51-A13*A26*
     $A32*A45*A51-A16*A22*A33*A45*A51+A12*A26*A33*A45*A51+A13*A22*A36*
     $A45*A51-A12*A23*A36*A45*A51-A15*A23*A32*A46*A51+A13*A25*A32*A46*
     $A51+A15*A22*A33*A46*A51-A12*A25*A33*A46*A51-A13*A22*A35*A46*A51+
     $A12*A23*A35*A46*A51-A16*A25*A33*A41*A52+A15*A26*A33*A41*A52+A16*
     $A23*A35*A41*A52-A13*A26*A35*A41*A52-A15*A23*A36*A41*A52+A13*A25*
     $A36*A41*A52+A16*A25*A31*A43*A52-A15*A26*A31*A43*A52-A16*A21*A35*
     $A43*A52+A11*A26*A35*A43*A52+A15*A21*A36*A43*A52-A11*A25*A36*A43*
     $A52-A16*A23*A31*A45*A52+A13*A26*A31*A45*A52+A16*A21*A33*A45*A52-
     $A11*A26*A33*A45*A52-A13*A21*A36*A45*A52+A11*A23*A36*A45*A52+A15*
     $A23*A31*A46*A52-A13*A25*A31*A46*A52-A15*A21*A33*A46*A52+A11*A25*
     $A33*A46*A52+A13*A21*A35*A46*A52-A11*A23*A35*A46*A52+A16*A25*A32*
     $A41*A53-A15*A26*A32*A41*A53-A16*A22*A35*A41*A53+A12*A26*A35*A41*
     $A53+A15*A22*A36*A41*A53-A12*A25*A36*A41*A53-A16*A25*A31*A42*A53+
     $A15*A26*A31*A42*A53+A16*A21*A35*A42*A53-A11*A26*A35*A42*A53-A15*
     $A21*A36*A42*A53+A11*A25*A36*A42*A53+A16*A22*A31*A45*A53-A12*A26*
     $A31*A45*A53-A16*A21*A32*A45*A53+A11*A26*A32*A45*A53+A12*A21*A36*
     $A45*A53-A11*A22*A36*A45*A53-A15*A22*A31*A46*A53+A12*A25*A31*A46*
     $A53+A15*A21*A32*A46*A53-A11*A25*A32*A46*A53-A12*A21*A35*A46*A53+
     $A11*A22*A35*A46*A53-A16*A23*A32*A41*A55+A13*A26*A32*A41*A55+A16*
     $A22*A33*A41*A55-A12*A26*A33*A41*A55-A13*A22*A36*A41*A55+A12*A23*
     $A36*A41*A55+A16*A23*A31*A42*A55-A13*A26*A31*A42*A55-A16*A21*A33*
     $A42*A55+A11*A26*A33*A42*A55+A13*A21*A36*A42*A55-A11*A23*A36*A42*
     $A55-A16*A22*A31*A43*A55+A12*A26*A31*A43*A55+A16*A21*A32*A43*A55-
     $A11*A26*A32*A43*A55-A12*A21*A36*A43*A55+A11*A22*A36*A43*A55+A13*
     $A22*A31*A46*A55-A12*A23*A31*A46*A55-A13*A21*A32*A46*A55+A11*A23*
     $A32*A46*A55+A12*A21*A33*A46*A55-A11*A22*A33*A46*A55+A15*A23*A32*
     $A41*A56-A13*A25*A32*A41*A56-A15*A22*A33*A41*A56+A12*A25*A33*A41*
     $A56+A13*A22*A35*A41*A56-A12*A23*A35*A41*A56-A15*A23*A31*A42*A56+
     $A13*A25*A31*A42*A56+A15*A21*A33*A42*A56-A11*A25*A33*A42*A56-A13*
     $A21*A35*A42*A56+A11*A23*A35*A42*A56+A15*A22*A31*A43*A56-A12*A25*
     $A31*A43*A56-A15*A21*A32*A43*A56+A11*A25*A32*A43*A56+A12*A21*A35*
     $A43*A56-A11*A22*A35*A43*A56-A13*A22*A31*A45*A56+A12*A23*A31*A45*
     $A56+A13*A21*A32*A45*A56-A11*A23*A32*A45*A56-A12*A21*A33*A45*A56+
     $A11*A22*A33*A45*A56)*A64-(A16*A24*A33*A42*A51-A14*A26*A33*A42*
     $A51-A16*A23*A34*A42*A51+A13*A26*A34*A42*A51+A14*A23*A36*A42*A51-
     $A13*A24*A36*A42*A51-A16*A24*A32*A43*A51+A14*A26*A32*A43*A51+A16*
     $A22*A34*A43*A51-A12*A26*A34*A43*A51-A14*A22*A36*A43*A51+A12*A24*
     $A36*A43*A51+A16*A23*A32*A44*A51-A13*A26*A32*A44*A51-A16*A22*A33*
     $A44*A51+A12*A26*A33*A44*A51+A13*A22*A36*A44*A51-A12*A23*A36*A44*
     $A51-A14*A23*A32*A46*A51+A13*A24*A32*A46*A51+A14*A22*A33*A46*A51-
     $A12*A24*A33*A46*A51-A13*A22*A34*A46*A51+A12*A23*A34*A46*A51-A16*
     $A24*A33*A41*A52+A14*A26*A33*A41*A52+A16*A23*A34*A41*A52-A13*A26*
     $A34*A41*A52-A14*A23*A36*A41*A52+A13*A24*A36*A41*A52+A16*A24*A31*
     $A43*A52-A14*A26*A31*A43*A52-A16*A21*A34*A43*A52+A11*A26*A34*A43*
     $A52+A14*A21*A36*A43*A52-A11*A24*A36*A43*A52-A16*A23*A31*A44*A52+
     $A13*A26*A31*A44*A52+A16*A21*A33*A44*A52-A11*A26*A33*A44*A52-A13*
     $A21*A36*A44*A52+A11*A23*A36*A44*A52+A14*A23*A31*A46*A52-A13*A24*
     $A31*A46*A52-A14*A21*A33*A46*A52+A11*A24*A33*A46*A52+A13*A21*A34*
     $A46*A52-A11*A23*A34*A46*A52+A16*A24*A32*A41*A53-A14*A26*A32*A41*
     $A53-A16*A22*A34*A41*A53+A12*A26*A34*A41*A53+A14*A22*A36*A41*A53-
     $A12*A24*A36*A41*A53-A16*A24*A31*A42*A53+A14*A26*A31*A42*A53+A16*
     $A21*A34*A42*A53-A11*A26*A34*A42*A53-A14*A21*A36*A42*A53+A11*A24*
     $A36*A42*A53+A16*A22*A31*A44*A53-A12*A26*A31*A44*A53-A16*A21*A32*
     $A44*A53+A11*A26*A32*A44*A53+A12*A21*A36*A44*A53-A11*A22*A36*A44*
     $A53-A14*A22*A31*A46*A53+A12*A24*A31*A46*A53+A14*A21*A32*A46*A53-
     $A11*A24*A32*A46*A53-A12*A21*A34*A46*A53+A11*A22*A34*A46*A53-A16*
     $A23*A32*A41*A54+A13*A26*A32*A41*A54+A16*A22*A33*A41*A54-A12*A26*
     $A33*A41*A54-A13*A22*A36*A41*A54+A12*A23*A36*A41*A54+A16*A23*A31*
     $A42*A54-A13*A26*A31*A42*A54-A16*A21*A33*A42*A54+A11*A26*A33*A42*
     $A54+A13*A21*A36*A42*A54-A11*A23*A36*A42*A54-A16*A22*A31*A43*A54+
     $A12*A26*A31*A43*A54+A16*A21*A32*A43*A54-A11*A26*A32*A43*A54-A12*
     $A21*A36*A43*A54+A11*A22*A36*A43*A54+A13*A22*A31*A46*A54-A12*A23*
     $A31*A46*A54-A13*A21*A32*A46*A54+A11*A23*A32*A46*A54+A12*A21*A33*
     $A46*A54-A11*A22*A33*A46*A54+A14*A23*A32*A41*A56-A13*A24*A32*A41*
     $A56-A14*A22*A33*A41*A56+A12*A24*A33*A41*A56+A13*A22*A34*A41*A56-
     $A12*A23*A34*A41*A56-A14*A23*A31*A42*A56+A13*A24*A31*A42*A56+A14*
     $A21*A33*A42*A56-A11*A24*A33*A42*A56-A13*A21*A34*A42*A56+A11*A23*
     $A34*A42*A56+A14*A22*A31*A43*A56-A12*A24*A31*A43*A56-A14*A21*A32*
     $A43*A56+A11*A24*A32*A43*A56+A12*A21*A34*A43*A56-A11*A22*A34*A43*
     $A56-A13*A22*A31*A44*A56+A12*A23*A31*A44*A56+A13*A21*A32*A44*A56-
     $A11*A23*A32*A44*A56-A12*A21*A33*A44*A56+A11*A22*A33*A44*A56)*A65+
     $(A15*A24*A33*A42*A51-A14*A25*A33*A42*A51-A15*A23*A34*A42*A51+A13*
     $A25*A34*A42*A51+A14*A23*A35*A42*A51-A13*A24*A35*A42*A51-A15*A24*
     $A32*A43*A51+A14*A25*A32*A43*A51+A15*A22*A34*A43*A51-A12*A25*A34*
     $A43*A51-A14*A22*A35*A43*A51+A12*A24*A35*A43*A51+A15*A23*A32*A44*
     $A51-A13*A25*A32*A44*A51-A15*A22*A33*A44*A51+A12*A25*A33*A44*A51+
     $A13*A22*A35*A44*A51-A12*A23*A35*A44*A51-A14*A23*A32*A45*A51+A13*
     $A24*A32*A45*A51+A14*A22*A33*A45*A51-A12*A24*A33*A45*A51-A13*A22*
     $A34*A45*A51+A12*A23*A34*A45*A51-A15*A24*A33*A41*A52+A14*A25*A33*
     $A41*A52+A15*A23*A34*A41*A52-A13*A25*A34*A41*A52-A14*A23*A35*A41*
     $A52+A13*A24*A35*A41*A52+A15*A24*A31*A43*A52-A14*A25*A31*A43*A52-
     $A15*A21*A34*A43*A52+A11*A25*A34*A43*A52+A14*A21*A35*A43*A52-A11*
     $A24*A35*A43*A52-A15*A23*A31*A44*A52+A13*A25*A31*A44*A52+A15*A21*
     $A33*A44*A52-A11*A25*A33*A44*A52-A13*A21*A35*A44*A52+A11*A23*A35*
     $A44*A52+A14*A23*A31*A45*A52-A13*A24*A31*A45*A52-A14*A21*A33*A45*
     $A52+A11*A24*A33*A45*A52+A13*A21*A34*A45*A52-A11*A23*A34*A45*A52+
     $A15*A24*A32*A41*A53-A14*A25*A32*A41*A53-A15*A22*A34*A41*A53+A12*
     $A25*A34*A41*A53+A14*A22*A35*A41*A53-A12*A24*A35*A41*A53-A15*A24*
     $A31*A42*A53+A14*A25*A31*A42*A53+A15*A21*A34*A42*A53-A11*A25*A34*
     $A42*A53-A14*A21*A35*A42*A53+A11*A24*A35*A42*A53+A15*A22*A31*A44*
     $A53-A12*A25*A31*A44*A53-A15*A21*A32*A44*A53+A11*A25*A32*A44*A53+
     $A12*A21*A35*A44*A53-A11*A22*A35*A44*A53-A14*A22*A31*A45*A53+A12*
     $A24*A31*A45*A53+A14*A21*A32*A45*A53-A11*A24*A32*A45*A53-A12*A21*
     $A34*A45*A53+A11*A22*A34*A45*A53-A15*A23*A32*A41*A54+A13*A25*A32*
     $A41*A54+A15*A22*A33*A41*A54-A12*A25*A33*A41*A54-A13*A22*A35*A41*
     $A54+A12*A23*A35*A41*A54+A15*A23*A31*A42*A54-A13*A25*A31*A42*A54-
     $A15*A21*A33*A42*A54+A11*A25*A33*A42*A54+A13*A21*A35*A42*A54-A11*
     $A23*A35*A42*A54-A15*A22*A31*A43*A54+A12*A25*A31*A43*A54+A15*A21*
     $A32*A43*A54-A11*A25*A32*A43*A54-A12*A21*A35*A43*A54+A11*A22*A35*
     $A43*A54+A13*A22*A31*A45*A54-A12*A23*A31*A45*A54-A13*A21*A32*A45*
     $A54+A11*A23*A32*A45*A54+A12*A21*A33*A45*A54-A11*A22*A33*A45*A54+
     $A14*A23*A32*A41*A55-A13*A24*A32*A41*A55-A14*A22*A33*A41*A55+A12*
     $A24*A33*A41*A55+A13*A22*A34*A41*A55-A12*A23*A34*A41*A55-A14*A23*
     $A31*A42*A55+A13*A24*A31*A42*A55+A14*A21*A33*A42*A55-A11*A24*A33*
     $A42*A55-A13*A21*A34*A42*A55+A11*A23*A34*A42*A55+A14*A22*A31*A43*
     $A55-A12*A24*A31*A43*A55-A14*A21*A32*A43*A55+A11*A24*A32*A43*A55+
     $A12*A21*A34*A43*A55-A11*A22*A34*A43*A55-A13*A22*A31*A44*A55+A12*
     $A23*A31*A44*A55+A13*A21*A32*A44*A55-A11*A23*A32*A44*A55-A12*A21*
     $A33*A44*A55+A11*A22*A33*A44*A55)*A66

c      write(*,*)"Determinant = "
c      write(*,*)DET
      IF (ABS(DET).LE.1.0d-10) THEN
        AINV = 0.0D0
      END IF

      COFACTOR(1,1) = A26*A35*A44*A53*A62-A25*A36*A44*A53*A62-A26*A34*
     $A45*A53*A62+A24*A36*A45*A53*A62+A25*A34*A46*A53*A62-A24*A35*A46*
     $A53*A62-A26*A35*A43*A54*A62+A25*A36*A43*A54*A62+A26*A33*A45*A54*
     $A62-A23*A36*A45*A54*A62-A25*A33*A46*A54*A62+A23*A35*A46*A54*A62+
     $A26*A34*A43*A55*A62-A24*A36*A43*A55*A62-A26*A33*A44*A55*A62+A23*
     $A36*A44*A55*A62+A24*A33*A46*A55*A62-A23*A34*A46*A55*A62-A25*A34*
     $A43*A56*A62+A24*A35*A43*A56*A62+A25*A33*A44*A56*A62-A23*A35*A44*
     $A56*A62-A24*A33*A45*A56*A62+A23*A34*A45*A56*A62-A26*A35*A44*A52*
     $A63+A25*A36*A44*A52*A63+A26*A34*A45*A52*A63-A24*A36*A45*A52*A63-
     $A25*A34*A46*A52*A63+A24*A35*A46*A52*A63+A26*A35*A42*A54*A63-A25*
     $A36*A42*A54*A63-A26*A32*A45*A54*A63+A22*A36*A45*A54*A63+A25*A32*
     $A46*A54*A63-A22*A35*A46*A54*A63-A26*A34*A42*A55*A63+A24*A36*A42*
     $A55*A63+A26*A32*A44*A55*A63-A22*A36*A44*A55*A63-A24*A32*A46*A55*
     $A63+A22*A34*A46*A55*A63+A25*A34*A42*A56*A63-A24*A35*A42*A56*A63-
     $A25*A32*A44*A56*A63+A22*A35*A44*A56*A63+A24*A32*A45*A56*A63-A22*
     $A34*A45*A56*A63+A26*A35*A43*A52*A64-A25*A36*A43*A52*A64-A26*A33*
     $A45*A52*A64+A23*A36*A45*A52*A64+A25*A33*A46*A52*A64-A23*A35*A46*
     $A52*A64-A26*A35*A42*A53*A64+A25*A36*A42*A53*A64+A26*A32*A45*A53*
     $A64-A22*A36*A45*A53*A64-A25*A32*A46*A53*A64+A22*A35*A46*A53*A64+
     $A26*A33*A42*A55*A64-A23*A36*A42*A55*A64-A26*A32*A43*A55*A64+A22*
     $A36*A43*A55*A64+A23*A32*A46*A55*A64-A22*A33*A46*A55*A64-A25*A33*
     $A42*A56*A64+A23*A35*A42*A56*A64+A25*A32*A43*A56*A64-A22*A35*A43*
     $A56*A64-A23*A32*A45*A56*A64+A22*A33*A45*A56*A64-A26*A34*A43*A52*
     $A65+A24*A36*A43*A52*A65+A26*A33*A44*A52*A65-A23*A36*A44*A52*A65-
     $A24*A33*A46*A52*A65+A23*A34*A46*A52*A65+A26*A34*A42*A53*A65-A24*
     $A36*A42*A53*A65-A26*A32*A44*A53*A65+A22*A36*A44*A53*A65+A24*A32*
     $A46*A53*A65-A22*A34*A46*A53*A65-A26*A33*A42*A54*A65+A23*A36*A42*
     $A54*A65+A26*A32*A43*A54*A65-A22*A36*A43*A54*A65-A23*A32*A46*A54*
     $A65+A22*A33*A46*A54*A65+A24*A33*A42*A56*A65-A23*A34*A42*A56*A65-
     $A24*A32*A43*A56*A65+A22*A34*A43*A56*A65+A23*A32*A44*A56*A65-A22*
     $A33*A44*A56*A65+A25*A34*A43*A52*A66-A24*A35*A43*A52*A66-A25*A33*
     $A44*A52*A66+A23*A35*A44*A52*A66+A24*A33*A45*A52*A66-A23*A34*A45*
     $A52*A66-A25*A34*A42*A53*A66+A24*A35*A42*A53*A66+A25*A32*A44*A53*
     $A66-A22*A35*A44*A53*A66-A24*A32*A45*A53*A66+A22*A34*A45*A53*A66+
     $A25*A33*A42*A54*A66-A23*A35*A42*A54*A66-A25*A32*A43*A54*A66+A22*
     $A35*A43*A54*A66+A23*A32*A45*A54*A66-A22*A33*A45*A54*A66-A24*A33*
     $A42*A55*A66+A23*A34*A42*A55*A66+A24*A32*A43*A55*A66-A22*A34*A43*
     $A55*A66-A23*A32*A44*A55*A66+A22*A33*A44*A55*A66

      COFACTOR(2,1) = -A16*A35*A44*A53*A62+A15*A36*A44*A53*A62+A16*A34*
     $A45*A53*A62-A14*A36*A45*A53*A62-A15*A34*A46*A53*A62+A14*A35*A46*
     $A53*A62+A16*A35*A43*A54*A62-A15*A36*A43*A54*A62-A16*A33*A45*A54*
     $A62+A13*A36*A45*A54*A62+A15*A33*A46*A54*A62-A13*A35*A46*A54*A62-
     $A16*A34*A43*A55*
     $A62+A14*A36*A43*A55*A62+A16*A33*A44*A55*A62-A13*A36*A44*A55*A62-
     $A14*A33*
     $A46*A55*A62+A13*A34*A46*A55*A62+A15*A34*A43*A56*A62-A14*A35*A43*
     $A56*A62-A15*A33*A44*A56*A62+A13*A35*A44*A56*A62+A14*A33*A45*A56*
     $A62-A13*A34*
     $A45*A56*A62+A16*A35*A44*A52*A63-A15*A36*A44*A52*A63-A16*A34*A45*
     $A52*A63+A14*A36*A45*A52*A63+A15*A34*A46*A52*A63-A14*A35*A46*A52*
     $A63-A16*A35*
     $A42*A54*A63+A15*A36*A42*A54*A63+A16*A32*A45*A54*A63-A12*A36*A45*
     $A54*A63-A15*A32*A46*A54*A63+A12*A35*A46*A54*A63+A16*A34*A42*A55*
     $A63-A14*A36*A42*A55*A63-A16*A32*A44*A55*A63+A12*A36*A44*A55*A63+
     $A14*A32*A46*A55*
     $A63-A12*A34*A46*A55*A63-A15*A34*A42*A56*A63+A14*A35*A42*A56*A63+
     $A15*A32*
     $A44*A56*A63-A12*A35*A44*A56*A63-A14*A32*A45*A56*A63+A12*A34*A45*
     $A56*A63-A16*A35*A43*A52*A64+A15*A36*A43*A52*A64+A16*A33*A45*A52*
     $A64-A13*A36*A45*A52*A64-A15*A33*A46*A52*A64+A13*A35*A46*A52*A64+
     $A16*A35*A42*A53*A64-A15*A36*A42*A53*A64-A16*A32*A45*A53*A64+A12*
     $A36*A45*A53*A64+A15*A32*A46*A53*A64-A12*A35*A46*A53*A64-A16*A33*
     $A42*A55*A64+A13*A36*A42*A55*A64+A16*A32*A43*A55*A64-A12*A36*A43*
     $A55*A64-A13*A32*A46*A55*A64+A12*A33*A46*A55*A64+A15*A33*A42*A56*
     $A64-A13*A35*A42*A56*A64-A15*A32*A43*A56*A64+A12*A35*A43*A56*A64+
     $A13*A32*A45*A56*A64-A12*A33*A45*A56*A64+A16*A34*A43*A52*A65-A14*
     $A36*A43*A52*A65-A16*A33*A44*A52*A65+A13*A36*A44*A52*A65+A14*A33*
     $A46*A52*A65-A13*A34*A46*A52*A65-A16*A34*A42*A53*A65+A14*A36*
     $A42*A53*A65+A16*A32*A44*A53*A65-A12*A36*A44*A53*A65-A14*A32*A46*
     $A53*A65+A12*A34*A46*A53*A65+A16*A33*A42*A54*A65-A13*A36*A42*A54*
     $A65-A16*A32*A43*A54*A65+A12*A36*A43*A54*A65+A13*A32*A46*A54*A65-
     $A12*A33*A46*A54*A65-A14*A33*A42*A56*A65+A13*A34*A42*A56*A65+A14*
     $A32*A43*A56*A65-A12*A34*A43*A56*A65-A13*A32*A44*A56*A65+A12*A33*
     $A44*A56*A65-A15*A34*A43*A52*A66+A14*A35*A43*A52*A66+A15*A33*A44*
     $A52*A66-A13*A35*A44*A52*A66-A14*A33*A45*A52*A66+A13*A34*A45*A52*
     $A66+A15*A34*A42*A53*A66-A14*A35*A42*A53*A66-A15*A32*A44*A53*A66+
     $A12*A35*A44*A53*A66+A14*A32*A45*A53*A66-A12*A34*A45*A53*A66-A15*
     $A33*A42*A54*A66+A13*A35*A42*A54*A66+A15*A32*A43*A54*A66-A12*A35*
     $A43*A54*A66-A13*A32*A45*A54*A66+A12*A33*A45*A54*A66+A14*A33*A42*
     $A55*A66-A13*A34*A42*A55*A66-A14*A32*A43*A55*A66+A12*A34*A43*A55*
     $A66+A13*A32*A44*A55*A66-A12*A33*A44*A55*A66

      COFACTOR(3,1) = A16*A25*A44*A53*A62-A15*A26*
     $A44*A53*A62-A16*A24*A45*A53*A62+A14*A26*A45*A53*A62+A15*A24*A46*
     $A53*A62-A14*A25*A46*A53*A62-A16*A25*A43*A54*A62+A15*A26*A43*A54*
     $A62+A16*A23*A45*A54*A62-A13*A26*A45*A54*A62-A15*A23*A46*A54*A62+
     $A13*A25*A46*A54*A62+A16*A24*A43*A55*A62-A14*A26*A43*A55*A62-A16*
     $A23*A44*A55*A62+A13*A26*A44*A55*A62+A14*A23*A46*A55*A62-A13*A24*
     $A46*A55*A62-A15*A24*A43*A56*A62+A14*A25*A43*A56*A62+A15*A23*A44*
     $A56*A62-A13*A25*A44*A56*A62-A14*A23*A45*A56*A62+A13*A24*A45*A56*
     $A62-A16*A25*A44*A52*A63+A15*A26*A44*A52*A63+A16*A24*A45*A52*A63-
     $A14*A26*A45*A52*A63-A15*A24*A46*A52*A63+A14*A25*A46*A52*A63+A16*
     $A25*A42*A54*A63-A15*A26*A42*A54*A63-A16*A22*A45*A54*A63+A12*A26*
     $A45*A54*A63+A15*A22*A46*A54*A63-A12*A25*A46*A54*A63-A16*A24*A42*
     $A55*A63+A14*A26*A42*A55*A63+A16*A22*A44*A55*A63-A12*A26*A44*A55*
     $A63-A14*A22*A46*A55*A63+A12*A24*A46*A55*A63+A15*A24*A42*A56*A63-
     $A14*A25*A42*A56*A63-A15*A22*A44*A56*A63+A12*A25*A44*A56*A63+A14*
     $A22*A45*A56*A63-A12*A24*A45*A56*A63+A16*A25*A43*A52*A64-A15*A26*
     $A43*A52*A64-A16*A23*A45*A52*A64+A13*A26*A45*A52*A64+A15*A23*A46*
     $A52*A64-A13*A25*A46*A52*A64-A16*A25*A42*A53*A64+A15*A26*A42*A53*
     $A64+A16*A22*A45*A53*A64-A12*A26*A45*A53*A64-A15*A22*A46*A53*A64+
     $A12*A25*A46*A53*A64+A16*A23*A42*A55*A64-A13*A26*A42*A55*A64-A16*
     $A22*A43*A55*A64+A12*A26*A43*A55*A64+A13*A22*A46*A55*A64-A12*A23*
     $A46*A55*A64-A15*A23*A42*A56*A64+A13*A25*A42*A56*A64+A15*A22*A43*
     $A56*A64-A12*A25*A43*A56*A64-A13*A22*A45*A56*A64+A12*A23*A45*A56*
     $A64-A16*A24*A43*A52*A65+A14*A26*A43*A52*A65+A16*A23*A44*A52*A65-
     $A13*A26*A44*A52*A65-A14*A23*A46*A52*A65+A13*A24*A46*A52*A65+A16*
     $A24*A42*A53*A65-A14*A26*A42*A53*A65-A16*A22*A44*A53*A65+A12*A26*
     $A44*A53*A65+A14*A22*A46*A53*A65-A12*A24*A46*A53*A65-A16*A23*A42*
     $A54*A65+A13*A26*A42*A54*A65+A16*A22*A43*A54*A65-A12*A26*A43*A54*
     $A65-A13*A22*A46*A54*A65+A12*A23*A46*A54*A65+A14*A23*A42*A56*A65-
     $A13*A24*A42*A56*A65-A14*A22*A43*A56*A65+A12*A24*A43*A56*A65+A13*
     $A22*A44*A56*A65-A12*A23*A44*A56*A65+A15*A24*A43*A52*A66-A14*A25*
     $A43*A52*A66-A15*A23*A44*A52*A66+A13*A25*A44*A52*A66+A14*A23*A45*
     $A52*A66-A13*A24*A45*A52*A66-A15*A24*A42*A53*A66+A14*A25*A42*A53*
     $A66+A15*A22*A44*A53*A66-A12*A25*A44*A53*A66-A14*A22*A45*A53*A66+
     $A12*A24*A45*A53*A66+A15*A23*A42*A54*A66-A13*A25*A42*A54*A66-A15*
     $A22*A43*A54*A66+A12*A25*A43*A54*A66+A13*A22*A45*A54*A66-A12*A23*
     $A45*A54*A66-A14*A23*A42*A55*A66+A13*A24*A42*A55*A66+A14*A22*A43*
     $A55*A66-A12*A24*A43*A55*A66-A13*A22*A44*A55*A66+A12*A23*A44*A55*
     $A66

      COFACTOR(4,1) = -A16*A25*
     $A34*A53*A62+A15*A26*A34*A53*A62+A16*A24*A35*A53*A62-A14*A26*A35*
     $A53*A62-A15*A24*A36*A53*A62+A14*A25*A36*A53*A62+A16*A25*A33*A54*
     $A62-A15*A26*A33*A54*A62-A16*A23*A35*A54*A62+A13*A26*A35*A54*A62+
     $A15*A23*A36*A54*A62-A13*A25*A36*A54*A62-A16*A24*A33*A55*A62+A14*
     $A26*A33*A55*A62+A16*A23*A34*A55*A62-A13*A26*A34*A55*A62-A14*A23*
     $A36*A55*A62+A13*A24*A36*A55*A62+A15*A24*A33*A56*A62-A14*A25*A33*
     $A56*A62-A15*A23*A34*A56*A62+A13*A25*A34*A56*A62+A14*A23*A35*A56*
     $A62-A13*A24*A35*A56*A62+A16*A25*A34*A52*A63-A15*A26*A34*A52*A63-
     $A16*A24*A35*A52*A63+A14*A26*A35*A52*A63+A15*A24*A36*A52*A63-A14*
     $A25*A36*A52*A63-A16*A25*A32*A54*A63+A15*A26*A32*A54*A63+A16*A22*
     $A35*A54*A63-A12*A26*A35*A54*A63-A15*A22*A36*A54*A63+A12*A25*A36*
     $A54*A63+A16*A24*A32*A55*A63-A14*A26*A32*A55*A63-A16*A22*A34*A55*
     $A63+A12*A26*A34*A55*A63+A14*A22*A36*A55*A63-A12*A24*A36*A55*A63-
     $A15*A24*A32*A56*A63+A14*A25*A32*A56*A63+A15*A22*A34*A56*A63-A12*
     $A25*A34*A56*A63-A14*A22*A35*A56*A63+A12*A24*A35*A56*A63-A16*A25*
     $A33*A52*A64+A15*A26*A33*A52*A64+A16*A23*A35*A52*A64-A13*A26*A35*
     $A52*A64-A15*A23*A36*A52*A64+A13*A25*A36*A52*A64+A16*A25*A32*A53*
     $A64-A15*A26*A32*A53*A64-A16*A22*A35*A53*A64+A12*A26*A35*A53*A64+
     $A15*A22*A36*A53*A64-A12*A25*A36*A53*A64-A16*A23*A32*A55*A64+A13*
     $A26*A32*A55*A64+A16*A22*A33*A55*A64-A12*A26*A33*A55*A64-A13*A22*
     $A36*A55*A64+A12*A23*A36*A55*A64+A15*A23*A32*A56*A64-A13*A25*A32*
     $A56*A64-A15*A22*A33*A56*A64+A12*A25*A33*A56*A64+A13*A22*A35*A56*
     $A64-A12*A23*A35*A56*A64+A16*A24*A33*A52*A65-A14*A26*A33*A52*A65-
     $A16*A23*A34*A52*A65+A13*A26*A34*A52*A65+A14*A23*A36*A52*A65-A13*
     $A24*A36*A52*A65-A16*A24*A32*A53*A65+A14*A26*A32*A53*A65+A16*A22*
     $A34*A53*A65-A12*A26*A34*A53*A65-A14*A22*A36*A53*A65+A12*A24*A36*
     $A53*A65+A16*A23*A32*A54*A65-A13*A26*A32*A54*A65-A16*A22*A33*A54*
     $A65+A12*A26*A33*A54*A65+A13*A22*A36*A54*A65-A12*A23*A36*A54*A65-
     $A14*A23*A32*A56*A65+A13*A24*A32*A56*A65+A14*A22*A33*A56*A65-A12*
     $A24*A33*A56*A65-A13*A22*A34*A56*A65+A12*A23*A34*A56*A65-A15*A24*
     $A33*A52*A66+A14*A25*A33*A52*A66+A15*A23*A34*A52*A66-A13*A25*A34*
     $A52*A66-A14*A23*A35*A52*A66+A13*A24*A35*A52*A66+A15*A24*A32*A53*
     $A66-A14*A25*A32*A53*A66-A15*A22*A34*A53*A66+A12*A25*A34*A53*A66+
     $A14*A22*A35*A53*A66-A12*A24*A35*A53*A66-A15*A23*A32*A54*A66+A13*
     $A25*A32*A54*A66+A15*A22*A33*A54*A66-A12*A25*A33*A54*A66-A13*A22*
     $A35*A54*A66+A12*A23*A35*A54*A66+A14*A23*A32*A55*A66-A13*A24*A32*
     $A55*A66-A14*A22*A33*A55*A66+A12*A24*A33*A55*A66+A13*A22*A34*A55*
     $A66-A12*A23*A34*A55*A66

      COFACTOR(5,1) = A16*A25*A34*A43*A62-A15*A26*A34*A43*A62-A16*A24*
     $A35*A43*
     $A62+A14*A26*A35*A43*A62+A15*A24*A36*A43*A62-A14*A25*A36*A43*A62-
     $A16*A25*
     $A33*A44*A62+A15*A26*A33*A44*A62+A16*A23*A35*A44*A62-A13*A26*A35*
     $A44*A62-A15*A23*A36*A44*A62+A13*A25*A36*A44*A62+A16*A24*A33*A45*
     $A62-A14*A26*A33*A45*A62-A16*A23*A34*A45*A62+A13*A26*A34*A45*A62+
     $A14*A23*A36*A45*A62-A13*A24*A36*A45*A62-A15*A24*A33*A46*A62+A14*
     $A25*A33*A46*A62+A15*A23*A34*A46*A62-A13*A25*A34*A46*A62-A14*A23*
     $A35*A46*A62+A13*A24*A35*A46*A62-A16*A25*A34*A42*A63+A15*A26*A34*
     $A42*A63+A16*A24*A35*A42*A63-A14*A26*A35*A42*A63-A15*A24*A36*A42*
     $A63+A14*A25*A36*A42*A63+A16*A25*A32*A44*A63-A15*A26*A32*A44*A63-
     $A16*A22*A35*A44*A63+A12*A26*A35*A44*A63+A15*A22*A36*A44*A63-A12*
     $A25*A36*A44*A63-A16*A24*A32*A45*A63+A14*A26*A32*A45*A63+A16*A22*
     $A34*A45*A63-A12*A26*A34*A45*A63-A14*A22*A36*A45*A63+A12*A24*A36*
     $A45*A63+A15*A24*A32*A46*A63-A14*A25*A32*A46*A63-A15*A22*A34*A46*
     $A63+A12*A25*A34*A46*A63+A14*A22*A35*A46*A63-A12*A24*A35*A46*A63+
     $A16*A25*A33*A42*A64-A15*A26*A33*A42*A64-A16*A23*A35*A42*A64+A13*
     $A26*A35*A42*A64+A15*A23*A36*A42*A64-A13*A25*A36*A42*A64-A16*A25*
     $A32*A43*A64+A15*A26*A32*A43*A64+A16*A22*A35*A43*A64-A12*A26*A35*
     $A43*A64-A15*A22*A36*A43*A64+A12*A25*A36*A43*A64+A16*A23*A32*A45*
     $A64-A13*A26*A32*A45*A64-A16*A22*A33*A45*A64+A12*A26*A33*A45*A64+
     $A13*A22*A36*A45*A64-A12*A23*A36*A45*A64-A15*A23*A32*A46*A64+A13*
     $A25*A32*A46*A64+A15*A22*A33*A46*A64-A12*A25*A33*A46*A64-A13*A22*
     $A35*A46*A64+A12*A23*A35*A46*A64-A16*A24*A33*A42*A65+A14*A26*A33*
     $A42*A65+A16*A23*A34*A42*A65-A13*A26*A34*A42*A65-A14*A23*A36*A42*
     $A65+A13*A24*A36*A42*A65+A16*A24*A32*A43*A65-A14*A26*A32*A43*A65-
     $A16*A22*A34*A43*A65+A12*A26*A34*A43*A65+A14*A22*A36*A43*A65-A12*
     $A24*A36*A43*A65-A16*A23*A32*A44*A65+A13*A26*A32*A44*A65+A16*A22*
     $A33*A44*A65-A12*A26*A33*A44*A65-A13*A22*A36*A44*A65+A12*A23*A36*
     $A44*A65+A14*A23*A32*A46*A65-A13*A24*A32*A46*A65-A14*A22*A33*A46*
     $A65+A12*A24*A33*A46*A65+A13*A22*A34*A46*A65-A12*A23*A34*A46*A65+
     $A15*A24*A33*A42*A66-A14*A25*A33*A42*A66-A15*A23*A34*A42*A66+A13*
     $A25*A34*A42*A66+A14*A23*A35*A42*A66-A13*A24*A35*A42*A66-A15*A24*
     $A32*A43*A66+A14*A25*A32*A43*A66+A15*A22*A34*A43*A66-A12*A25*A34*
     $A43*A66-A14*A22*A35*A43*A66+A12*A24*A35*A43*A66+A15*A23*A32*A44*
     $A66-A13*A25*A32*A44*A66-A15*A22*A33*A44*A66+A12*A25*A33*A44*A66+
     $A13*A22*A35*A44*A66-A12*A23*A35*A44*A66-A14*A23*A32*A45*A66+A13*
     $A24*A32*A45*A66+A14*A22*A33*A45*A66-A12*A24*A33*A45*A66-A13*A22*
     $A34*A45*A66+A12*A23*A34*A45*A66

      COFACTOR(6,1) = -A16*A25*A34*A43*A52+A15*A26*A34*A43*
     $A52+A16*A24*A35*A43*A52-A14*A26*A35*A43*A52-A15*A24*A36*A43*A52+
     $A14*A25*A36*A43*A52+A16*A25*A33*A44*A52-A15*A26*A33*A44*A52-A16*
     $A23*A35*A44*A52+A13*A26*A35*A44*A52+A15*A23*A36*A44*A52-A13*A25*
     $A36*A44*A52-A16*A24*A33*A45*A52+A14*A26*A33*A45*A52+A16*A23*A34*
     $A45*A52-A13*A26*A34*A45*A52-A14*A23*A36*A45*A52+A13*A24*A36*A45*
     $A52+A15*A24*A33*A46*A52-A14*A25*A33*A46*A52-A15*A23*A34*A46*A52+
     $A13*A25*A34*A46*A52+A14*A23*A35*A46*A52-A13*A24*A35*A46*A52+A16*
     $A25*A34*A42*A53-A15*A26*A34*A42*A53-A16*A24*A35*A42*A53+A14*A26*
     $A35*A42*A53+A15*A24*A36*A42*A53-A14*A25*A36*A42*A53-A16*A25*A32*
     $A44*A53+A15*A26*A32*A44*A53+A16*A22*A35*A44*A53-A12*A26*A35*A44*
     $A53-A15*A22*A36*A44*A53+A12*A25*A36*A44*A53+A16*A24*A32*A45*A53-
     $A14*A26*A32*A45*A53-A16*A22*A34*A45*A53+A12*A26*A34*A45*A53+A14*
     $A22*A36*A45*A53-A12*A24*A36*A45*A53-A15*A24*A32*A46*A53+A14*A25*
     $A32*A46*A53+A15*A22*A34*A46*A53-A12*A25*A34*A46*A53-A14*A22*A35*
     $A46*A53+A12*A24*A35*A46*A53-A16*A25*A33*A42*A54+A15*A26*A33*A42*
     $A54+A16*A23*A35*A42*A54-A13*A26*A35*A42*A54-A15*A23*A36*A42*A54+
     $A13*A25*A36*A42*A54+A16*A25*A32*A43*A54-A15*A26*A32*A43*A54-A16*
     $A22*A35*A43*A54+A12*A26*A35*A43*A54+A15*A22*A36*A43*A54-A12*A25*
     $A36*A43*A54-A16*A23*A32*A45*A54+A13*A26*A32*A45*A54+A16*A22*A33*
     $A45*A54-A12*A26*A33*A45*A54-A13*A22*A36*A45*A54+A12*A23*A36*A45*
     $A54+A15*A23*A32*A46*A54-A13*A25*A32*A46*A54-A15*A22*A33*A46*A54+
     $A12*A25*A33*A46*A54+A13*A22*A35*A46*A54-A12*A23*A35*A46*A54+A16*
     $A24*A33*A42*A55-A14*A26*A33*A42*A55-A16*A23*A34*A42*A55+A13*A26*
     $A34*A42*A55+A14*A23*A36*A42*A55-A13*A24*A36*A42*A55-A16*A24*A32*
     $A43*A55+A14*A26*A32*A43*A55+A16*A22*A34*A43*A55-A12*A26*A34*A43*
     $A55-A14*A22*A36*A43*A55+A12*A24*A36*A43*A55+A16*A23*A32*A44*A55-
     $A13*A26*A32*A44*A55-A16*A22*A33*A44*A55+A12*A26*A33*A44*A55+A13*
     $A22*A36*A44*A55-A12*A23*A36*A44*A55-A14*A23*A32*A46*A55+A13*A24*
     $A32*A46*A55+A14*A22*A33*A46*A55-A12*A24*A33*A46*A55-A13*A22*A34*
     $A46*A55+A12*A23*A34*A46*A55-A15*A24*A33*A42*A56+A14*A25*A33*A42*
     $A56+A15*A23*A34*A42*A56-A13*A25*A34*A42*A56-A14*A23*A35*A42*A56+
     $A13*A24*A35*A42*A56+A15*A24*A32*A43*A56-A14*A25*A32*A43*A56-A15*
     $A22*A34*A43*A56+A12*A25*A34*A43*A56+A14*A22*A35*A43*A56-A12*A24*
     $A35*A43*A56-A15*A23*A32*A44*A56+A13*A25*A32*A44*A56+A15*A22*A33*
     $A44*A56-A12*A25*A33*A44*A56-A13*A22*A35*A44*A56+A12*A23*A35*A44*
     $A56+A14*A23*A32*A45*A56-A13*A24*A32*A45*A56-A14*A22*A33*A45*A56+
     $A12*A24*A33*A45*A56+A13*A22*A34*A45*A56-A12*A23*A34*A45*A56

      COFACTOR(1,2) = -A26*A35*A44*A53*
     $A61+A25*A36*A44*A53*A61+A26*A34*A45*A53*A61-A24*A36*A45*A53*A61-
     $A25*A34*A46*A53*A61+A24*A35*A46*A53*A61+A26*A35*A43*A54*A61-A25*
     $A36*A43*A54*A61-A26*A33*A45*A54*A61+A23*A36*A45*A54*A61+A25*A33*
     $A46*A54*A61-A23*A35*A46*A54*A61-A26*A34*A43*A55*A61+A24*A36*A43*
     $A55*A61+A26*A33*A44*A55*A61-A23*A36*A44*A55*A61-A24*A33*A46*A55*
     $A61+A23*A34*A46*A55*A61+A25*A34*A43*A56*A61-A24*A35*A43*A56*A61-
     $A25*A33*A44*A56*A61+A23*A35*A44*A56*A61+A24*A33*A45*A56*A61-A23*
     $A34*A45*A56*A61+A26*A35*A44*A51*A63-A25*A36*A44*A51*A63-A26*A34*
     $A45*A51*A63+A24*A36*A45*A51*A63+A25*A34*A46*A51*A63-A24*A35*A46*
     $A51*A63-A26*A35*A41*A54*A63+A25*A36*A41*A54*A63+A26*A31*A45*A54*
     $A63-A21*A36*A45*A54*A63-A25*A31*A46*A54*A63+A21*A35*A46*A54*A63+
     $A26*A34*A41*A55*A63-A24*A36*A41*A55*A63-A26*A31*A44*A55*A63+A21*
     $A36*A44*A55*A63+A24*A31*A46*A55*A63-A21*A34*A46*A55*A63-A25*A34*
     $A41*A56*A63+A24*A35*A41*A56*A63+A25*A31*A44*A56*A63-A21*A35*A44*
     $A56*A63-A24*A31*A45*A56*A63+A21*A34*A45*A56*A63-A26*A35*A43*A51*
     $A64+A25*A36*A43*A51*A64+A26*A33*A45*A51*A64-A23*A36*A45*A51*A64-
     $A25*A33*A46*A51*A64+A23*A35*A46*A51*A64+A26*A35*A41*A53*A64-A25*
     $A36*A41*A53*A64-A26*A31*A45*A53*A64+A21*A36*A45*A53*A64+A25*A31*
     $A46*A53*A64-A21*A35*A46*A53*A64-A26*A33*A41*A55*A64+A23*A36*A41*
     $A55*A64+A26*A31*A43*A55*A64-A21*A36*A43*A55*A64-A23*A31*A46*A55*
     $A64+A21*A33*A46*A55*A64+A25*A33*A41*A56*A64-A23*A35*A41*A56*A64-
     $A25*A31*A43*A56*A64+A21*A35*A43*A56*A64+A23*A31*A45*A56*A64-A21*
     $A33*A45*A56*A64+A26*A34*A43*A51*A65-A24*A36*A43*A51*A65-A26*A33*
     $A44*A51*A65+A23*A36*A44*A51*A65+A24*A33*A46*A51*A65-A23*A34*A46*
     $A51*A65-A26*A34*A41*A53*A65+A24*A36*A41*A53*A65+A26*A31*A44*A53*
     $A65-A21*A36*A44*A53*A65-A24*A31*A46*A53*A65+A21*A34*A46*A53*A65+
     $A26*A33*A41*A54*A65-A23*A36*A41*A54*A65-A26*A31*A43*A54*A65+A21*
     $A36*A43*A54*A65+A23*A31*A46*A54*A65-A21*A33*A46*A54*A65-A24*A33*
     $A41*A56*A65+A23*A34*A41*A56*A65+A24*A31*A43*A56*A65-A21*A34*A43*
     $A56*A65-A23*A31*A44*A56*A65+A21*A33*A44*A56*A65-A25*A34*A43*A51*
     $A66+A24*A35*A43*A51*A66+A25*A33*A44*A51*A66-A23*A35*A44*A51*A66-
     $A24*A33*A45*A51*A66+A23*A34*A45*A51*A66+A25*A34*A41*A53*A66-A24*
     $A35*A41*A53*A66-A25*A31*A44*A53*A66+A21*A35*A44*A53*A66+A24*A31*
     $A45*A53*A66-A21*A34*A45*A53*A66-A25*A33*A41*A54*A66+A23*A35*A41*
     $A54*A66+A25*A31*A43*A54*A66-A21*A35*A43*A54*A66-A23*A31*A45*A54*
     $A66+A21*A33*A45*A54*A66+A24*A33*A41*A55*A66-A23*A34*A41*A55*A66-
     $A24*A31*A43*A55*A66+A21*A34*A43*A55*A66+A23*A31*A44*A55*A66-A21*
     $A33*A44*A55*A66

      COFACTOR(2,2) = A16*A35*A44*A53*A61-A15*A36*A44*A53*A61-A16*A34*
     $A45*A53*A61+A14*A36*A45*A53*A61+A15*A34*A46*A53*A61-A14*A35*A46*
     $A53*A61-A16*A35*A43*A54*A61+A15*A36*A43*A54*A61+A16*A33*A45*A54*
     $A61-A13*A36*A45*A54*A61-A15*A33*A46*A54*A61+A13*A35*A46*A54*A61+
     $A16*A34*A43*A55*A61-A14*A36*A43*A55*A61-A16*A33*A44*A55*A61+A13*
     $A36*A44*A55*A61+A14*A33*A46*A55*A61-A13*A34*A46*A55*A61-A15*A34*
     $A43*A56*A61+A14*A35*A43*A56*A61+A15*A33*A44*A56*A61-A13*A35*A44*
     $A56*A61-A14*A33*A45*A56*A61+A13*A34*A45*A56*A61-A16*A35*A44*A51*
     $A63+A15*A36*A44*A51*A63+A16*A34*A45*A51*A63-A14*A36*A45*A51*A63-
     $A15*A34*A46*A51*A63+A14*A35*A46*A51*A63+A16*A35*A41*A54*A63-A15*
     $A36*A41*A54*A63-A16*A31*A45*A54*A63+A11*A36*A45*A54*A63+A15*A31*
     $A46*A54*A63-A11*A35*A46*A54*A63-A16*A34*A41*A55*A63+A14*A36*A41*
     $A55*A63+A16*A31*A44*A55*A63-A11*A36*A44*A55*A63-A14*A31*A46*A55*
     $A63+A11*A34*A46*A55*A63+A15*A34*A41*A56*A63-A14*A35*A41*A56*A63-
     $A15*A31*A44*A56*A63+A11*A35*A44*A56*A63+A14*A31*A45*A56*A63-A11*
     $A34*A45*A56*A63+A16*A35*A43*A51*A64-A15*A36*A43*A51*A64-A16*A33*
     $A45*A51*A64+A13*A36*A45*A51*A64+A15*A33*A46*A51*A64-A13*A35*A46*
     $A51*A64-A16*A35*A41*A53*A64+A15*A36*A41*A53*A64+A16*A31*A45*A53*
     $A64-A11*A36*A45*A53*A64-A15*A31*A46*A53*A64+A11*A35*A46*A53*A64+
     $A16*A33*A41*A55*A64-A13*A36*A41*A55*A64-A16*A31*A43*A55*A64+A11*
     $A36*A43*A55*A64+A13*A31*A46*A55*A64-A11*A33*A46*A55*A64-A15*A33*
     $A41*A56*A64+A13*A35*A41*A56*A64+A15*A31*A43*A56*A64-A11*A35*A43*
     $A56*A64-A13*A31*A45*A56*A64+A11*A33*A45*A56*A64-A16*A34*A43*A51*
     $A65+A14*A36*A43*A51*A65+A16*A33*A44*A51*A65-A13*A36*A44*A51*A65-
     $A14*A33*A46*A51*A65+A13*A34*A46*A51*A65+A16*A34*A41*A53*A65-A14*
     $A36*A41*A53*A65-A16*A31*A44*A53*A65+A11*A36*A44*A53*A65+A14*A31*
     $A46*A53*A65-A11*A34*A46*A53*A65-A16*A33*A41*A54*A65+A13*A36*A41*
     $A54*A65+A16*A31*A43*A54*A65-A11*A36*A43*A54*A65-A13*A31*A46*A54*
     $A65+A11*A33*A46*A54*A65+A14*A33*A41*A56*A65-A13*A34*A41*A56*A65-
     $A14*A31*A43*A56*A65+A11*A34*A43*A56*A65+A13*A31*A44*A56*A65-A11*
     $A33*A44*A56*A65+A15*A34*A43*A51*A66-A14*A35*A43*A51*A66-A15*A33*
     $A44*A51*A66+A13*A35*A44*A51*A66+A14*A33*A45*A51*A66-A13*A34*A45*
     $A51*A66-A15*A34*A41*A53*A66+A14*A35*A41*A53*A66+A15*A31*A44*A53*
     $A66-A11*A35*A44*A53*A66-A14*A31*A45*A53*A66+A11*A34*A45*A53*A66+
     $A15*A33*A41*A54*A66-A13*A35*A41*A54*A66-A15*A31*A43*A54*A66+A11*
     $A35*A43*A54*A66+A13*A31*A45*A54*A66-A11*A33*A45*A54*A66-A14*A33*
     $A41*A55*A66+A13*A34*A41*A55*A66+A14*A31*A43*A55*A66-A11*A34*A43*
     $A55*A66-A13*A31*A44*A55*A66+A11*A33*A44*A55*A66

      COFACTOR(3,2) = -A16*A25*A44*A53*A61+A15*A26*A44*A53*A61+A16*A24*
     $A45*A53*A61-A14*A26*A45*A53*A61-A15*A24*A46*A53*A61+A14*A25*A46*
     $A53*A61+A16*A25*A43*A54*A61-A15*A26*A43*A54*A61-A16*A23*A45*A54*
     $A61+A13*A26*A45*A54*A61+A15*A23*A46*A54*A61-A13*A25*A46*A54*A61-
     $A16*A24*A43*A55*A61+A14*A26*A43*A55*A61+A16*A23*A44*A55*A61-A13*
     $A26*A44*A55*A61-A14*A23*A46*A55*A61+A13*A24*A46*A55*A61+A15*A24*
     $A43*A56*A61-A14*A25*A43*A56*A61-A15*A23*A44*A56*A61+A13*A25*A44*
     $A56*A61+A14*A23*A45*A56*A61-A13*A24*A45*A56*A61+A16*A25*A44*A51*
     $A63-A15*A26*A44*A51*A63-A16*A24*A45*A51*A63+A14*A26*A45*A51*A63+
     $A15*A24*A46*A51*A63-A14*A25*A46*A51*A63-A16*A25*A41*A54*A63+A15*
     $A26*A41*A54*A63+A16*A21*A45*A54*A63-A11*A26*A45*A54*A63-A15*A21*
     $A46*A54*A63+A11*A25*A46*A54*A63+A16*A24*A41*A55*A63-A14*A26*A41*
     $A55*A63-A16*A21*A44*A55*A63+A11*A26*A44*A55*A63+A14*A21*A46*A55*
     $A63-A11*A24*A46*A55*A63-A15*A24*A41*A56*A63+A14*A25*A41*A56*A63+
     $A15*A21*A44*A56*A63-A11*A25*A44*A56*A63-A14*A21*A45*A56*A63+A11*
     $A24*A45*A56*A63-A16*A25*A43*A51*A64+A15*A26*A43*A51*A64+A16*A23*
     $A45*A51*A64-A13*A26*A45*A51*A64-A15*A23*A46*A51*A64+A13*A25*A46*
     $A51*A64+A16*A25*A41*A53*A64-A15*A26*A41*A53*A64-A16*A21*A45*A53*
     $A64+A11*A26*A45*A53*A64+A15*A21*A46*A53*A64-A11*A25*A46*A53*A64-
     $A16*A23*A41*A55*A64+A13*A26*A41*A55*A64+A16*A21*A43*A55*A64-A11*
     $A26*A43*A55*A64-A13*A21*A46*A55*A64+A11*A23*A46*A55*A64+A15*A23*
     $A41*A56*A64-A13*A25*A41*A56*A64-A15*A21*A43*A56*A64+A11*A25*A43*
     $A56*A64+A13*A21*A45*A56*A64-A11*A23*A45*A56*A64+A16*A24*A43*A51*
     $A65-A14*A26*A43*A51*A65-A16*A23*A44*A51*A65+A13*A26*A44*A51*A65+
     $A14*A23*A46*A51*A65-A13*A24*A46*A51*A65-A16*A24*A41*A53*A65+A14*
     $A26*A41*A53*A65+A16*A21*A44*A53*A65-A11*A26*A44*A53*A65-A14*A21*
     $A46*A53*A65+A11*A24*A46*A53*A65+A16*A23*A41*A54*A65-A13*A26*A41*
     $A54*A65-A16*A21*A43*A54*A65+A11*A26*A43*A54*A65+A13*A21*A46*A54*
     $A65-A11*A23*A46*A54*A65-A14*A23*A41*A56*A65+A13*A24*A41*A56*A65+
     $A14*A21*A43*A56*A65-A11*A24*A43*A56*A65-A13*A21*A44*A56*A65+A11*
     $A23*A44*A56*A65-A15*A24*A43*A51*A66+A14*A25*A43*A51*A66+A15*A23*
     $A44*A51*A66-A13*A25*A44*A51*A66-A14*A23*A45*A51*A66+A13*A24*A45*
     $A51*A66+A15*A24*A41*A53*A66-A14*A25*A41*A53*A66-A15*A21*A44*A53*
     $A66+A11*A25*A44*A53*A66+A14*A21*A45*A53*A66-A11*A24*A45*A53*A66-
     $A15*A23*A41*A54*A66+A13*A25*A41*A54*A66+A15*A21*A43*A54*A66-A11*
     $A25*A43*A54*A66-A13*A21*A45*A54*A66+A11*A23*A45*A54*A66+A14*A23*
     $A41*A55*A66-A13*A24*A41*A55*A66-A14*A21*A43*A55*A66+A11*A24*A43*
     $A55*A66+A13*A21*A44*A55*A66-A11*A23*A44*A55*A66

      COFACTOR(4,2) = A16*A25*A34*A53*A61-A15*A26*
     $A34*A53*A61-A16*A24*A35*A53*A61+A14*A26*A35*A53*A61+A15*A24*A36*
     $A53*A61-A14*A25*A36*A53*A61-A16*A25*A33*A54*A61+A15*A26*A33*A54*
     $A61+A16*A23*A35*A54*A61-A13*A26*A35*A54*A61-A15*A23*A36*A54*A61+
     $A13*A25*A36*A54*A61+A16*A24*A33*A55*A61-A14*A26*A33*A55*A61-A16*
     $A23*A34*A55*A61+A13*A26*A34*A55*A61+A14*A23*A36*A55*A61-A13*A24*
     $A36*A55*A61-A15*A24*A33*A56*A61+A14*A25*A33*A56*A61+A15*A23*A34*
     $A56*A61-A13*A25*A34*A56*A61-A14*A23*A35*A56*A61+A13*A24*A35*A56*
     $A61-A16*A25*A34*A51*A63+A15*A26*A34*A51*A63+A16*A24*A35*A51*A63-
     $A14*A26*A35*A51*A63-A15*A24*A36*A51*A63+A14*A25*A36*A51*A63+A16*
     $A25*A31*A54*A63-A15*A26*A31*A54*A63-A16*A21*A35*A54*A63+A11*A26*
     $A35*A54*A63+A15*A21*A36*A54*A63-A11*A25*A36*A54*A63-A16*A24*A31*
     $A55*A63+A14*A26*A31*A55*A63+A16*A21*A34*A55*A63-A11*A26*A34*A55*
     $A63-A14*A21*A36*A55*A63+A11*A24*A36*A55*A63+A15*A24*A31*A56*A63-
     $A14*A25*A31*A56*A63-A15*A21*A34*A56*A63+A11*A25*A34*A56*A63+A14*
     $A21*A35*A56*A63-A11*A24*A35*A56*A63+A16*A25*A33*A51*A64-A15*A26*
     $A33*A51*A64-A16*A23*A35*A51*A64+A13*A26*A35*A51*A64+A15*A23*A36*
     $A51*A64-A13*A25*A36*A51*A64-A16*A25*A31*A53*A64+A15*A26*A31*A53*
     $A64+A16*A21*A35*A53*A64-A11*A26*A35*A53*A64-A15*A21*A36*A53*A64+
     $A11*A25*A36*A53*A64+A16*A23*A31*A55*A64-A13*A26*A31*A55*A64-A16*
     $A21*A33*A55*A64+A11*A26*A33*A55*A64+A13*A21*A36*A55*A64-A11*A23*
     $A36*A55*A64-A15*A23*A31*A56*A64+A13*A25*A31*A56*A64+A15*A21*A33*
     $A56*A64-A11*A25*A33*A56*A64-A13*A21*A35*A56*A64+A11*A23*A35*A56*
     $A64-A16*A24*A33*A51*A65+A14*A26*A33*A51*A65+A16*A23*A34*A51*A65-
     $A13*A26*A34*A51*A65-A14*A23*A36*A51*A65+A13*A24*A36*A51*A65+A16*
     $A24*A31*A53*A65-A14*A26*A31*A53*A65-A16*A21*A34*A53*A65+A11*A26*
     $A34*A53*A65+A14*A21*A36*A53*A65-A11*A24*A36*A53*A65-A16*A23*A31*
     $A54*A65+A13*A26*A31*A54*A65+A16*A21*A33*A54*A65-A11*A26*A33*A54*
     $A65-A13*A21*A36*A54*A65+A11*A23*A36*A54*A65+A14*A23*A31*A56*A65-
     $A13*A24*A31*A56*A65-A14*A21*A33*A56*A65+A11*A24*A33*A56*A65+A13*
     $A21*A34*A56*A65-A11*A23*A34*A56*A65+A15*A24*A33*A51*A66-A14*A25*
     $A33*A51*A66-A15*A23*A34*A51*A66+A13*A25*A34*A51*A66+A14*A23*A35*
     $A51*A66-A13*A24*A35*A51*A66-A15*A24*A31*A53*A66+A14*A25*A31*A53*
     $A66+A15*A21*A34*A53*A66-A11*A25*A34*A53*A66-A14*A21*A35*A53*A66+
     $A11*A24*A35*A53*A66+A15*A23*A31*A54*A66-A13*A25*A31*A54*A66-A15*
     $A21*A33*A54*A66+A11*A25*A33*A54*A66+A13*A21*A35*A54*A66-A11*A23*
     $A35*A54*A66-A14*A23*A31*A55*A66+A13*A24*A31*A55*A66+A14*A21*A33*
     $A55*A66-A11*A24*A33*A55*A66-A13*A21*A34*A55*A66+A11*A23*A34*A55*
     $A66

      COFACTOR(5,2) = -A16*A25*
     $A34*A43*A61+A15*A26*A34*A43*A61+A16*A24*A35*A43*A61-A14*A26*A35*
     $A43*A61-A15*A24*A36*A43*A61+A14*A25*A36*A43*A61+A16*A25*A33*A44*
     $A61-A15*A26*A33*A44*A61-A16*A23*A35*A44*A61+A13*A26*A35*A44*A61+
     $A15*A23*A36*A44*A61-A13*A25*A36*A44*A61-A16*A24*A33*A45*A61+A14*
     $A26*A33*A45*A61+A16*A23*A34*A45*A61-A13*A26*A34*A45*A61-A14*A23*
     $A36*A45*A61+A13*A24*A36*A45*A61+A15*A24*A33*A46*A61-A14*A25*A33*
     $A46*A61-A15*A23*A34*A46*A61+A13*A25*A34*A46*A61+A14*A23*A35*A46*
     $A61-A13*A24*A35*A46*A61+A16*A25*A34*A41*A63-A15*A26*A34*A41*A63-
     $A16*A24*A35*A41*A63+A14*A26*A35*A41*A63+A15*A24*A36*A41*A63-A14*
     $A25*A36*A41*A63-A16*A25*A31*A44*A63+A15*A26*A31*A44*A63+A16*A21*
     $A35*A44*A63-A11*A26*A35*A44*A63-A15*A21*A36*A44*A63+A11*A25*A36*
     $A44*A63+A16*A24*A31*A45*A63-A14*A26*A31*A45*A63-A16*A21*A34*A45*
     $A63+A11*A26*A34*A45*A63+A14*A21*A36*A45*A63-A11*A24*A36*A45*A63-
     $A15*A24*A31*A46*A63+A14*A25*A31*A46*A63+A15*A21*A34*A46*A63-A11*
     $A25*A34*A46*A63-A14*A21*A35*A46*A63+A11*A24*A35*A46*A63-A16*A25*
     $A33*A41*A64+A15*A26*A33*A41*A64+A16*A23*A35*A41*A64-A13*A26*A35*
     $A41*A64-A15*A23*A36*A41*A64+A13*A25*A36*A41*A64+A16*A25*A31*A43*
     $A64-A15*A26*A31*A43*A64-A16*A21*A35*A43*A64+A11*A26*A35*A43*A64+
     $A15*A21*A36*A43*A64-A11*A25*A36*A43*A64-A16*A23*A31*A45*A64+A13*
     $A26*A31*A45*A64+A16*A21*A33*A45*A64-A11*A26*A33*A45*A64-A13*A21*
     $A36*A45*A64+A11*A23*A36*A45*A64+A15*A23*A31*A46*A64-A13*A25*A31*
     $A46*A64-A15*A21*A33*A46*A64+A11*A25*A33*A46*A64+A13*A21*A35*A46*
     $A64-A11*A23*A35*A46*A64+A16*A24*A33*A41*A65-A14*A26*A33*A41*A65-
     $A16*A23*A34*A41*A65+A13*A26*A34*A41*A65+A14*A23*A36*A41*A65-A13*
     $A24*A36*A41*A65-A16*A24*A31*A43*A65+A14*A26*A31*A43*A65+A16*A21*
     $A34*A43*A65-A11*A26*A34*A43*A65-A14*A21*A36*A43*A65+A11*A24*A36*
     $A43*A65+A16*A23*A31*A44*A65-A13*A26*A31*A44*A65-A16*A21*A33*A44*
     $A65+A11*A26*A33*A44*A65+A13*A21*A36*A44*A65-A11*A23*A36*A44*A65-
     $A14*A23*A31*A46*A65+A13*A24*A31*A46*A65+A14*A21*A33*A46*A65-A11*
     $A24*A33*A46*A65-A13*A21*A34*A46*A65+A11*A23*A34*A46*A65-A15*A24*
     $A33*A41*A66+A14*A25*A33*A41*A66+A15*A23*A34*A41*A66-A13*A25*A34*
     $A41*A66-A14*A23*A35*A41*A66+A13*A24*A35*A41*A66+A15*A24*A31*A43*
     $A66-A14*A25*A31*A43*A66-A15*A21*A34*A43*A66+A11*A25*A34*A43*A66+
     $A14*A21*A35*A43*A66-A11*A24*A35*A43*A66-A15*A23*A31*A44*A66+A13*
     $A25*A31*A44*A66+A15*A21*A33*A44*A66-A11*A25*A33*A44*A66-A13*A21*
     $A35*A44*A66+A11*A23*A35*A44*A66+A14*A23*A31*A45*A66-A13*A24*A31*
     $A45*A66-A14*A21*A33*A45*A66+A11*A24*A33*A45*A66+A13*A21*A34*A45*
     $A66-A11*A23*A34*A45*A66

      COFACTOR(6,2) = A16*A25*A34*A43*A51-A15*A26*A34*A43*A51-A16*A24*
     $A35*A43*A51+A14*A26*A35*A43*A51+A15*A24*A36*A43*A51-A14*A25*A36*
     $A43*A51-A16*A25*A33*A44*A51+A15*A26*A33*A44*A51+A16*A23*A35*A44*
     $A51-A13*A26*A35*A44*A51-A15*A23*A36*A44*A51+A13*A25*A36*A44*A51+
     $A16*A24*A33*A45*A51-A14*A26*A33*A45*A51-A16*A23*A34*A45*A51+A13*
     $A26*A34*A45*A51+A14*A23*A36*A45*A51-A13*A24*A36*A45*A51-A15*A24*
     $A33*A46*A51+A14*A25*A33*A46*A51+A15*A23*A34*A46*A51-A13*A25*A34*
     $A46*A51-A14*A23*A35*A46*A51+A13*A24*A35*A46*A51-A16*A25*A34*A41*
     $A53+A15*A26*A34*A41*A53+A16*A24*A35*A41*A53-A14*A26*A35*A41*A53-
     $A15*A24*A36*A41*A53+A14*A25*A36*A41*A53+A16*A25*A31*A44*A53-A15*
     $A26*A31*A44*A53-A16*A21*A35*A44*A53+A11*A26*A35*A44*A53+A15*A21*
     $A36*A44*A53-A11*A25*A36*A44*A53-A16*A24*A31*A45*A53+A14*A26*A31*
     $A45*A53+A16*A21*A34*A45*A53-A11*A26*A34*A45*A53-A14*A21*A36*A45*
     $A53+A11*A24*A36*A45*A53+A15*A24*A31*A46*A53-A14*A25*A31*A46*A53-
     $A15*A21*A34*A46*A53+A11*A25*A34*A46*A53+A14*A21*A35*A46*A53-A11*
     $A24*A35*A46*A53+A16*A25*A33*A41*A54-A15*A26*A33*A41*A54-A16*A23*
     $A35*A41*A54+A13*A26*A35*A41*A54+A15*A23*A36*A41*A54-A13*A25*A36*
     $A41*A54-A16*A25*A31*A43*A54+A15*A26*A31*A43*A54+A16*A21*A35*A43*
     $A54-A11*A26*A35*A43*A54-A15*A21*A36*A43*A54+A11*A25*A36*A43*A54+
     $A16*A23*A31*A45*A54-A13*A26*A31*A45*A54-A16*A21*A33*A45*A54+A11*
     $A26*A33*A45*A54+A13*A21*A36*A45*A54-A11*A23*A36*A45*A54-A15*A23*
     $A31*A46*A54+A13*A25*A31*A46*A54+A15*A21*A33*A46*A54-A11*A25*A33*
     $A46*A54-A13*A21*A35*A46*A54+A11*A23*A35*A46*A54-A16*A24*A33*A41*
     $A55+A14*A26*A33*A41*A55+A16*A23*A34*A41*A55-A13*A26*A34*A41*A55-
     $A14*A23*A36*A41*A55+A13*A24*A36*A41*A55+A16*A24*A31*A43*A55-A14*
     $A26*A31*A43*A55-A16*A21*A34*A43*A55+A11*A26*A34*A43*A55+A14*A21*
     $A36*A43*A55-A11*A24*A36*A43*A55-A16*A23*A31*A44*A55+A13*A26*A31*
     $A44*A55+A16*A21*A33*A44*A55-A11*A26*A33*A44*A55-A13*A21*A36*A44*
     $A55+A11*A23*A36*A44*A55+A14*A23*A31*A46*A55-A13*A24*A31*A46*A55-
     $A14*A21*A33*A46*A55+A11*A24*A33*A46*A55+A13*A21*A34*A46*A55-A11*
     $A23*A34*A46*A55+A15*A24*A33*A41*A56-A14*A25*A33*A41*A56-A15*A23*
     $A34*A41*A56+A13*A25*A34*A41*A56+A14*A23*A35*A41*A56-A13*A24*A35*
     $A41*A56-A15*A24*A31*A43*A56+A14*A25*A31*A43*A56+A15*A21*A34*A43*
     $A56-A11*A25*A34*A43*A56-A14*A21*A35*A43*A56+A11*A24*A35*A43*A56+
     $A15*A23*A31*A44*A56-A13*A25*A31*A44*A56-A15*A21*A33*A44*A56+A11*
     $A25*A33*A44*A56+A13*A21*A35*A44*A56-A11*A23*A35*A44*A56-A14*A23*
     $A31*A45*A56+A13*A24*A31*A45*A56+A14*A21*A33*A45*A56-A11*A24*A33*
     $A45*A56-A13*A21*A34*A45*A56+A11*A23*A34*A45*A56

      COFACTOR(1,3) = A26*A35*A44*A52*A61-A25*A36*A44*A52*
     $A61-A26*A34*A45*A52*A61+A24*A36*A45*A52*A61+A25*A34*A46*A52*A61-
     $A24*A35*A46*A52*A61-A26*A35*A42*A54*A61+A25*A36*A42*A54*A61+A26*
     $A32*A45*A54*A61-A22*A36*A45*A54*A61-A25*A32*A46*A54*A61+A22*A35*
     $A46*A54*A61+A26*A34*A42*A55*A61-A24*A36*A42*A55*A61-A26*A32*A44*
     $A55*A61+A22*A36*A44*A55*A61+A24*A32*A46*A55*A61-A22*A34*A46*A55*
     $A61-A25*A34*A42*A56*A61+A24*A35*A42*A56*A61+A25*A32*A44*A56*A61-
     $A22*A35*A44*A56*A61-A24*A32*A45*A56*A61+A22*A34*A45*A56*A61-A26*
     $A35*A44*A51*A62+A25*A36*A44*A51*A62+A26*A34*A45*A51*A62-A24*A36*
     $A45*A51*A62-A25*A34*A46*A51*A62+A24*A35*A46*A51*A62+A26*A35*A41*
     $A54*A62-A25*A36*A41*A54*A62-A26*A31*A45*A54*A62+A21*A36*A45*A54*
     $A62+A25*A31*A46*A54*A62-A21*A35*A46*A54*A62-A26*A34*A41*A55*A62+
     $A24*A36*A41*A55*A62+A26*A31*A44*A55*A62-A21*A36*A44*A55*A62-A24*
     $A31*A46*A55*A62+A21*A34*A46*A55*A62+A25*A34*A41*A56*A62-A24*A35*
     $A41*A56*A62-A25*A31*A44*A56*A62+A21*A35*A44*A56*A62+A24*A31*A45*
     $A56*A62-A21*A34*A45*A56*A62+A26*A35*A42*A51*A64-A25*A36*A42*A51*
     $A64-A26*A32*A45*A51*A64+A22*A36*A45*A51*A64+A25*A32*A46*A51*A64-
     $A22*A35*A46*A51*A64-A26*A35*A41*A52*A64+A25*A36*A41*A52*A64+A26*
     $A31*A45*A52*A64-A21*A36*A45*A52*A64-A25*A31*A46*A52*A64+A21*A35*
     $A46*A52*A64+A26*A32*A41*A55*A64-A22*A36*A41*A55*A64-A26*A31*A42*
     $A55*A64+A21*A36*A42*A55*A64+A22*A31*A46*A55*A64-A21*A32*A46*A55*
     $A64-A25*A32*A41*A56*A64+A22*A35*A41*A56*A64+A25*A31*A42*A56*A64-
     $A21*A35*A42*A56*A64-A22*A31*A45*A56*A64+A21*A32*A45*A56*A64-A26*
     $A34*A42*A51*A65+A24*A36*A42*A51*A65+A26*A32*A44*A51*A65-A22*A36*
     $A44*A51*A65-A24*A32*A46*A51*A65+A22*A34*A46*A51*A65+A26*A34*A41*
     $A52*A65-A24*A36*A41*A52*A65-A26*A31*A44*A52*A65+A21*A36*A44*A52*
     $A65+A24*A31*A46*A52*A65-A21*A34*A46*A52*A65-A26*A32*A41*A54*A65+
     $A22*A36*A41*A54*A65+A26*A31*A42*A54*A65-A21*A36*A42*A54*A65-A22*
     $A31*A46*A54*A65+A21*A32*A46*A54*A65+A24*A32*A41*A56*A65-A22*A34*
     $A41*A56*A65-A24*A31*A42*A56*A65+A21*A34*A42*A56*A65+A22*A31*A44*
     $A56*A65-A21*A32*A44*A56*A65+A25*A34*A42*A51*A66-A24*A35*A42*A51*
     $A66-A25*A32*A44*A51*A66+A22*A35*A44*A51*A66+A24*A32*A45*A51*A66-
     $A22*A34*A45*A51*A66-A25*A34*A41*A52*A66+A24*A35*A41*A52*A66+A25*
     $A31*A44*A52*A66-A21*A35*A44*A52*A66-A24*A31*A45*A52*A66+A21*A34*
     $A45*A52*A66+A25*A32*A41*A54*A66-A22*A35*A41*A54*A66-A25*A31*A42*
     $A54*A66+A21*A35*A42*A54*A66+A22*A31*A45*A54*A66-A21*A32*A45*A54*
     $A66-A24*A32*A41*A55*A66+A22*A34*A41*A55*A66+A24*A31*A42*A55*A66-
     $A21*A34*A42*A55*A66-A22*A31*A44*A55*A66+A21*A32*A44*A55*A66

      COFACTOR(2,3) = -A16*A35*A44*A52*
     $A61+A15*A36*A44*A52*A61+A16*A34*A45*A52*A61-A14*A36*A45*A52*A61-
     $A15*A34*A46*A52*A61+A14*A35*A46*A52*A61+A16*A35*A42*A54*A61-A15*
     $A36*A42*A54*A61-A16*A32*A45*A54*A61+A12*A36*A45*A54*A61+A15*A32*
     $A46*A54*A61-A12*A35*A46*A54*A61-A16*A34*A42*A55*A61+A14*A36*A42*
     $A55*A61+A16*A32*A44*A55*A61-A12*A36*A44*A55*A61-A14*A32*A46*A55*
     $A61+A12*A34*A46*A55*A61+A15*A34*A42*A56*A61-A14*A35*A42*A56*A61-
     $A15*A32*A44*A56*A61+A12*A35*A44*A56*A61+A14*A32*A45*A56*A61-A12*
     $A34*A45*A56*A61+A16*A35*A44*A51*A62-A15*A36*A44*A51*A62-A16*A34*
     $A45*A51*A62+A14*A36*A45*A51*A62+A15*A34*A46*A51*A62-A14*A35*A46*
     $A51*A62-A16*A35*A41*A54*A62+A15*A36*A41*A54*A62+A16*A31*A45*A54*
     $A62-A11*A36*A45*A54*A62-A15*A31*A46*A54*A62+A11*A35*A46*A54*A62+
     $A16*A34*A41*A55*A62-A14*A36*A41*A55*A62-A16*A31*A44*A55*A62+A11*
     $A36*A44*A55*A62+A14*A31*A46*A55*A62-A11*A34*A46*A55*A62-A15*A34*
     $A41*A56*A62+A14*A35*A41*A56*A62+A15*A31*A44*A56*A62-A11*A35*A44*
     $A56*A62-A14*A31*A45*A56*A62+A11*A34*A45*A56*A62-A16*A35*A42*A51*
     $A64+A15*A36*A42*A51*A64+A16*A32*A45*A51*A64-A12*A36*A45*A51*A64-
     $A15*A32*A46*A51*A64+A12*A35*A46*A51*A64+A16*A35*A41*A52*A64-A15*
     $A36*A41*A52*A64-A16*A31*A45*A52*A64+A11*A36*A45*A52*A64+A15*A31*
     $A46*A52*A64-A11*A35*A46*A52*A64-A16*A32*A41*A55*A64+A12*A36*A41*
     $A55*A64+A16*A31*A42*A55*A64-A11*A36*A42*A55*A64-A12*A31*A46*A55*
     $A64+A11*A32*A46*A55*A64+A15*A32*A41*A56*A64-A12*A35*A41*A56*A64-
     $A15*A31*A42*A56*A64+A11*A35*A42*A56*A64+A12*A31*A45*A56*A64-A11*
     $A32*A45*A56*A64+A16*A34*A42*A51*A65-A14*A36*A42*A51*A65-A16*A32*
     $A44*A51*A65+A12*A36*A44*A51*A65+A14*A32*A46*A51*A65-A12*A34*A46*
     $A51*A65-A16*A34*A41*A52*A65+A14*A36*A41*A52*A65+A16*A31*A44*A52*
     $A65-A11*A36*A44*A52*A65-A14*A31*A46*A52*A65+A11*A34*A46*A52*A65+
     $A16*A32*A41*A54*A65-A12*A36*A41*A54*A65-A16*A31*A42*A54*A65+A11*
     $A36*A42*A54*A65+A12*A31*A46*A54*A65-A11*A32*A46*A54*A65-A14*A32*
     $A41*A56*A65+A12*A34*A41*A56*A65+A14*A31*A42*A56*A65-A11*A34*A42*
     $A56*A65-A12*A31*A44*A56*A65+A11*A32*A44*A56*A65-A15*A34*A42*A51*
     $A66+A14*A35*A42*A51*A66+A15*A32*A44*A51*A66-A12*A35*A44*A51*A66-
     $A14*A32*A45*A51*A66+A12*A34*A45*A51*A66+A15*A34*A41*A52*A66-A14*
     $A35*A41*A52*A66-A15*A31*A44*A52*A66+A11*A35*A44*A52*A66+A14*A31*
     $A45*A52*A66-A11*A34*A45*A52*A66-A15*A32*A41*A54*A66+A12*A35*A41*
     $A54*A66+A15*A31*A42*A54*A66-A11*A35*A42*A54*A66-A12*A31*A45*A54*
     $A66+A11*A32*A45*A54*A66+A14*A32*A41*A55*A66-A12*A34*A41*A55*A66-
     $A14*A31*A42*A55*A66+A11*A34*A42*A55*A66+A12*A31*A44*A55*A66-A11*
     $A32*A44*A55*A66

      COFACTOR(3,3) = A16*A25*A44*A52*A61-A15*A26*A44*A52*A61-A16*A24*
     $A45*A52*A61+A14*A26*A45*A52*A61+A15*A24*A46*A52*A61-A14*A25*A46*
     $A52*A61-A16*A25*A42*A54*A61+A15*A26*A42*A54*A61+A16*A22*A45*A54*
     $A61-A12*A26*A45*A54*A61-A15*A22*A46*A54*A61+A12*A25*A46*A54*A61+
     $A16*A24*A42*A55*A61-A14*A26*A42*A55*A61-A16*A22*A44*A55*A61+A12*
     $A26*A44*A55*A61+A14*A22*A46*A55*A61-A12*A24*A46*A55*A61-A15*A24*
     $A42*A56*A61+A14*A25*A42*A56*A61+A15*A22*A44*A56*A61-A12*A25*A44*
     $A56*A61-A14*A22*A45*A56*A61+A12*A24*A45*A56*A61-A16*A25*A44*A51*
     $A62+A15*A26*A44*A51*A62+A16*A24*A45*A51*A62-A14*A26*A45*A51*A62-
     $A15*A24*A46*A51*A62+A14*A25*A46*A51*A62+A16*A25*A41*A54*A62-A15*
     $A26*A41*A54*A62-A16*A21*A45*A54*A62+A11*A26*A45*A54*A62+A15*A21*
     $A46*A54*A62-A11*A25*A46*A54*A62-A16*A24*A41*A55*A62+A14*A26*A41*
     $A55*A62+A16*A21*A44*A55*A62-A11*A26*A44*A55*A62-A14*A21*A46*A55*
     $A62+A11*A24*A46*A55*A62+A15*A24*A41*A56*A62-A14*A25*A41*A56*A62-
     $A15*A21*A44*A56*A62+A11*A25*A44*A56*A62+A14*A21*A45*A56*A62-A11*
     $A24*A45*A56*A62+A16*A25*A42*A51*A64-A15*A26*A42*A51*A64-A16*A22*
     $A45*A51*A64+A12*A26*A45*A51*A64+A15*A22*A46*A51*A64-A12*A25*A46*
     $A51*A64-A16*A25*A41*A52*A64+A15*A26*A41*A52*A64+A16*A21*A45*A52*
     $A64-A11*A26*A45*A52*A64-A15*A21*A46*A52*A64+A11*A25*A46*A52*A64+
     $A16*A22*A41*A55*A64-A12*A26*A41*A55*A64-A16*A21*A42*A55*A64+A11*
     $A26*A42*A55*A64+A12*A21*A46*A55*A64-A11*A22*A46*A55*A64-A15*A22*
     $A41*A56*A64+A12*A25*A41*A56*A64+A15*A21*A42*A56*A64-A11*A25*A42*
     $A56*A64-A12*A21*A45*A56*A64+A11*A22*A45*A56*A64-A16*A24*A42*A51*
     $A65+A14*A26*A42*A51*A65+A16*A22*A44*A51*A65-A12*A26*A44*A51*A65-
     $A14*A22*A46*A51*A65+A12*A24*A46*A51*A65+A16*A24*A41*A52*A65-A14*
     $A26*A41*A52*A65-A16*A21*A44*A52*A65+A11*A26*A44*A52*A65+A14*A21*
     $A46*A52*A65-A11*A24*A46*A52*A65-A16*A22*A41*A54*A65+A12*A26*A41*
     $A54*A65+A16*A21*A42*A54*A65-A11*A26*A42*A54*A65-A12*A21*A46*A54*
     $A65+A11*A22*A46*A54*A65+A14*A22*A41*A56*A65-A12*A24*A41*A56*A65-
     $A14*A21*A42*A56*A65+A11*A24*A42*A56*A65+A12*A21*A44*A56*A65-A11*
     $A22*A44*A56*A65+A15*A24*A42*A51*A66-A14*A25*A42*A51*A66-A15*A22*
     $A44*A51*A66+A12*A25*A44*A51*A66+A14*A22*A45*A51*A66-A12*A24*A45*
     $A51*A66-A15*A24*A41*A52*A66+A14*A25*A41*A52*A66+A15*A21*A44*A52*
     $A66-A11*A25*A44*A52*A66-A14*A21*A45*A52*A66+A11*A24*A45*A52*A66+
     $A15*A22*A41*A54*A66-A12*A25*A41*A54*A66-A15*A21*A42*A54*A66+A11*
     $A25*A42*A54*A66+A12*A21*A45*A54*A66-A11*A22*A45*A54*A66-A14*A22*
     $A41*A55*A66+A12*A24*A41*A55*A66+A14*A21*A42*A55*A66-A11*A24*A42*
     $A55*A66-A12*A21*A44*A55*A66+A11*A22*A44*A55*A66

      COFACTOR(4,3) = -A16*A25*A34*A52*A61+A15*A26*A34*A52*A61+A16*A24*
     $A35*A52*A61-A14*A26*A35*A52*A61-A15*A24*A36*A52*A61+A14*A25*A36*
     $A52*A61+A16*A25*A32*A54*A61-A15*A26*A32*A54*A61-A16*A22*A35*A54*
     $A61+A12*A26*A35*A54*A61+A15*A22*A36*A54*A61-A12*A25*A36*A54*A61-
     $A16*A24*A32*A55*A61+A14*A26*A32*A55*A61+A16*A22*A34*A55*A61-A12*
     $A26*A34*A55*A61-A14*A22*A36*A55*A61+A12*A24*A36*A55*A61+A15*A24*
     $A32*A56*A61-A14*A25*A32*A56*A61-A15*A22*A34*A56*A61+A12*A25*A34*
     $A56*A61+A14*A22*A35*A56*A61-A12*A24*A35*A56*A61+A16*A25*A34*A51*
     $A62-A15*A26*A34*A51*A62-A16*A24*A35*A51*A62+A14*A26*A35*A51*A62+
     $A15*A24*A36*A51*A62-A14*A25*A36*A51*A62-A16*A25*A31*A54*A62+A15*
     $A26*A31*A54*A62+A16*A21*A35*A54*A62-A11*A26*A35*A54*A62-A15*A21*
     $A36*A54*A62+A11*A25*A36*A54*A62+A16*A24*A31*A55*A62-A14*A26*A31*
     $A55*A62-A16*A21*A34*A55*A62+A11*A26*A34*A55*A62+A14*A21*A36*A55*
     $A62-A11*A24*A36*A55*A62-A15*A24*A31*A56*A62+A14*A25*A31*A56*A62+
     $A15*A21*A34*A56*A62-A11*A25*A34*A56*A62-A14*A21*A35*A56*A62+A11*
     $A24*A35*A56*A62-A16*A25*A32*A51*A64+A15*A26*A32*A51*A64+A16*A22*
     $A35*A51*A64-A12*A26*A35*A51*A64-A15*A22*A36*A51*A64+A12*A25*A36*
     $A51*A64+A16*A25*A31*A52*A64-A15*A26*A31*A52*A64-A16*A21*A35*A52*
     $A64+A11*A26*A35*A52*A64+A15*A21*A36*A52*A64-A11*A25*A36*A52*A64-
     $A16*A22*A31*A55*A64+A12*A26*A31*A55*A64+A16*A21*A32*A55*A64-A11*
     $A26*A32*A55*A64-A12*A21*A36*A55*A64+A11*A22*A36*A55*A64+A15*A22*
     $A31*A56*A64-A12*A25*A31*A56*A64-A15*A21*A32*A56*A64+A11*A25*A32*
     $A56*A64+A12*A21*A35*A56*A64-A11*A22*A35*A56*A64+A16*A24*A32*A51*
     $A65-A14*A26*A32*A51*A65-A16*A22*A34*A51*A65+A12*A26*A34*A51*A65+
     $A14*A22*A36*A51*A65-A12*A24*A36*A51*A65-A16*A24*A31*A52*A65+A14*
     $A26*A31*A52*A65+A16*A21*A34*A52*A65-A11*A26*A34*A52*A65-A14*A21*
     $A36*A52*A65+A11*A24*A36*A52*A65+A16*A22*A31*A54*A65-A12*A26*A31*
     $A54*A65-A16*A21*A32*A54*A65+A11*A26*A32*A54*A65+A12*A21*A36*A54*
     $A65-A11*A22*A36*A54*A65-A14*A22*A31*A56*A65+A12*A24*A31*A56*A65+
     $A14*A21*A32*A56*A65-A11*A24*A32*A56*A65-A12*A21*A34*A56*A65+A11*
     $A22*A34*A56*A65-A15*A24*A32*A51*A66+A14*A25*A32*A51*A66+A15*A22*
     $A34*A51*A66-A12*A25*A34*A51*A66-A14*A22*A35*A51*A66+A12*A24*A35*
     $A51*A66+A15*A24*A31*A52*A66-A14*A25*A31*A52*A66-A15*A21*A34*A52*
     $A66+A11*A25*A34*A52*A66+A14*A21*A35*A52*A66-A11*A24*A35*A52*A66-
     $A15*A22*A31*A54*A66+A12*A25*A31*A54*A66+A15*A21*A32*A54*A66-A11*
     $A25*A32*A54*A66-A12*A21*A35*A54*A66+A11*A22*A35*A54*A66+A14*A22*
     $A31*A55*A66-A12*A24*A31*A55*A66-A14*A21*A32*A55*A66+A11*A24*A32*
     $A55*A66+A12*A21*A34*A55*A66-A11*A22*A34*A55*A66

      COFACTOR(5,3) = A16*A25*A34*A42*A61-A15*A26*
     $A34*A42*A61-A16*A24*A35*A42*A61+A14*A26*A35*A42*A61+A15*A24*A36*
     $A42*A61-A14*A25*A36*A42*A61-A16*A25*A32*A44*A61+A15*A26*A32*A44*
     $A61+A16*A22*A35*A44*A61-A12*A26*A35*A44*A61-A15*A22*A36*A44*A61+
     $A12*A25*A36*A44*A61+A16*A24*A32*A45*A61-A14*A26*A32*A45*A61-A16*
     $A22*A34*A45*A61+A12*A26*A34*A45*A61+A14*A22*A36*A45*A61-A12*A24*
     $A36*A45*A61-A15*A24*A32*A46*A61+A14*A25*A32*A46*A61+A15*A22*A34*
     $A46*A61-A12*A25*A34*A46*A61-A14*A22*A35*A46*A61+A12*A24*A35*A46*
     $A61-A16*A25*A34*A41*A62+A15*A26*A34*A41*A62+A16*A24*A35*A41*A62-
     $A14*A26*A35*A41*A62-A15*A24*A36*A41*A62+A14*A25*A36*A41*A62+A16*
     $A25*A31*A44*A62-A15*A26*A31*A44*A62-A16*A21*A35*A44*A62+A11*A26*
     $A35*A44*A62+A15*A21*A36*A44*A62-A11*A25*A36*A44*A62-A16*A24*A31*
     $A45*A62+A14*A26*A31*A45*A62+A16*A21*A34*A45*A62-A11*A26*A34*A45*
     $A62-A14*A21*A36*A45*A62+A11*A24*A36*A45*A62+A15*A24*A31*A46*A62-
     $A14*A25*A31*A46*A62-A15*A21*A34*A46*A62+A11*A25*A34*A46*A62+A14*
     $A21*A35*A46*A62-A11*A24*A35*A46*A62+A16*A25*A32*A41*A64-A15*A26*
     $A32*A41*A64-A16*A22*A35*A41*A64+A12*A26*A35*A41*A64+A15*A22*A36*
     $A41*A64-A12*A25*A36*A41*A64-A16*A25*A31*A42*A64+A15*A26*A31*A42*
     $A64+A16*A21*A35*A42*A64-A11*A26*A35*A42*A64-A15*A21*A36*A42*A64+
     $A11*A25*A36*A42*A64+A16*A22*A31*A45*A64-A12*A26*A31*A45*A64-A16*
     $A21*A32*A45*A64+A11*A26*A32*A45*A64+A12*A21*A36*A45*A64-A11*A22*
     $A36*A45*A64-A15*A22*A31*A46*A64+A12*A25*A31*A46*A64+A15*A21*A32*
     $A46*A64-A11*A25*A32*A46*A64-A12*A21*A35*A46*A64+A11*A22*A35*A46*
     $A64-A16*A24*A32*A41*A65+A14*A26*A32*A41*A65+A16*A22*A34*A41*A65-
     $A12*A26*A34*A41*A65-A14*A22*A36*A41*A65+A12*A24*A36*A41*A65+A16*
     $A24*A31*A42*A65-A14*A26*A31*A42*A65-A16*A21*A34*A42*A65+A11*A26*
     $A34*A42*A65+A14*A21*A36*A42*A65-A11*A24*A36*A42*A65-A16*A22*A31*
     $A44*A65+A12*A26*A31*A44*A65+A16*A21*A32*A44*A65-A11*A26*A32*A44*
     $A65-A12*A21*A36*A44*A65+A11*A22*A36*A44*A65+A14*A22*A31*A46*A65-
     $A12*A24*A31*A46*A65-A14*A21*A32*A46*A65+A11*A24*A32*A46*A65+A12*
     $A21*A34*A46*A65-A11*A22*A34*A46*A65+A15*A24*A32*A41*A66-A14*A25*
     $A32*A41*A66-A15*A22*A34*A41*A66+A12*A25*A34*A41*A66+A14*A22*A35*
     $A41*A66-A12*A24*A35*A41*A66-A15*A24*A31*A42*A66+A14*A25*A31*A42*
     $A66+A15*A21*A34*A42*A66-A11*A25*A34*A42*A66-A14*A21*A35*A42*A66+
     $A11*A24*A35*A42*A66+A15*A22*A31*A44*A66-A12*A25*A31*A44*A66-A15*
     $A21*A32*A44*A66+A11*A25*A32*A44*A66+A12*A21*A35*A44*A66-A11*A22*
     $A35*A44*A66-A14*A22*A31*A45*A66+A12*A24*A31*A45*A66+A14*A21*A32*
     $A45*A66-A11*A24*A32*A45*A66-A12*A21*A34*A45*A66+A11*A22*A34*A45*
     $A66

      COFACTOR(6,3) = -A16*A25*
     $A34*A42*A51+A15*A26*A34*A42*A51+A16*A24*A35*A42*A51-A14*A26*A35*
     $A42*A51-A15*A24*A36*A42*A51+A14*A25*A36*A42*A51+A16*A25*A32*A44*
     $A51-A15*A26*A32*A44*A51-A16*A22*A35*A44*A51+A12*A26*A35*A44*A51+
     $A15*A22*A36*A44*A51-A12*A25*A36*A44*A51-A16*A24*A32*A45*A51+A14*
     $A26*A32*A45*A51+A16*A22*A34*A45*A51-A12*A26*A34*A45*A51-A14*A22*
     $A36*A45*A51+A12*A24*A36*A45*A51+A15*A24*A32*A46*A51-A14*A25*A32*
     $A46*A51-A15*A22*A34*A46*A51+A12*A25*A34*A46*A51+A14*A22*A35*A46*
     $A51-A12*A24*A35*A46*A51+A16*A25*A34*A41*A52-A15*A26*A34*A41*A52-
     $A16*A24*A35*A41*A52+A14*A26*A35*A41*A52+A15*A24*A36*A41*A52-A14*
     $A25*A36*A41*A52-A16*A25*A31*A44*A52+A15*A26*A31*A44*A52+A16*A21*
     $A35*A44*A52-A11*A26*A35*A44*A52-A15*A21*A36*A44*A52+A11*A25*A36*
     $A44*A52+A16*A24*A31*A45*A52-A14*A26*A31*A45*A52-A16*A21*A34*A45*
     $A52+A11*A26*A34*A45*A52+A14*A21*A36*A45*A52-A11*A24*A36*A45*A52-
     $A15*A24*A31*A46*A52+A14*A25*A31*A46*A52+A15*A21*A34*A46*A52-A11*
     $A25*A34*A46*A52-A14*A21*A35*A46*A52+A11*A24*A35*A46*A52-A16*A25*
     $A32*A41*A54+A15*A26*A32*A41*A54+A16*A22*A35*A41*A54-A12*A26*A35*
     $A41*A54-A15*A22*A36*A41*A54+A12*A25*A36*A41*A54+A16*A25*A31*A42*
     $A54-A15*A26*A31*A42*A54-A16*A21*A35*A42*A54+A11*A26*A35*A42*A54+
     $A15*A21*A36*A42*A54-A11*A25*A36*A42*A54-A16*A22*A31*A45*A54+A12*
     $A26*A31*A45*A54+A16*A21*A32*A45*A54-A11*A26*A32*A45*A54-A12*A21*
     $A36*A45*A54+A11*A22*A36*A45*A54+A15*A22*A31*A46*A54-A12*A25*A31*
     $A46*A54-A15*A21*A32*A46*A54+A11*A25*A32*A46*A54+A12*A21*A35*A46*
     $A54-A11*A22*A35*A46*A54+A16*A24*A32*A41*A55-A14*A26*A32*A41*A55-
     $A16*A22*A34*A41*A55+A12*A26*A34*A41*A55+A14*A22*A36*A41*A55-A12*
     $A24*A36*A41*A55-A16*A24*A31*A42*A55+A14*A26*A31*A42*A55+A16*A21*
     $A34*A42*A55-A11*A26*A34*A42*A55-A14*A21*A36*A42*A55+A11*A24*A36*
     $A42*A55+A16*A22*A31*A44*A55-A12*A26*A31*A44*A55-A16*A21*A32*A44*
     $A55+A11*A26*A32*A44*A55+A12*A21*A36*A44*A55-A11*A22*A36*A44*A55-
     $A14*A22*A31*A46*A55+A12*A24*A31*A46*A55+A14*A21*A32*A46*A55-A11*
     $A24*A32*A46*A55-A12*A21*A34*A46*A55+A11*A22*A34*A46*A55-A15*A24*
     $A32*A41*A56+A14*A25*A32*A41*A56+A15*A22*A34*A41*A56-A12*A25*A34*
     $A41*A56-A14*A22*A35*A41*A56+A12*A24*A35*A41*A56+A15*A24*A31*A42*
     $A56-A14*A25*A31*A42*A56-A15*A21*A34*A42*A56+A11*A25*A34*A42*A56+
     $A14*A21*A35*A42*A56-A11*A24*A35*A42*A56-A15*A22*A31*A44*A56+A12*
     $A25*A31*A44*A56+A15*A21*A32*A44*A56-A11*A25*A32*A44*A56-A12*A21*
     $A35*A44*A56+A11*A22*A35*A44*A56+A14*A22*A31*A45*A56-A12*A24*A31*
     $A45*A56-A14*A21*A32*A45*A56+A11*A24*A32*A45*A56+A12*A21*A34*A45*
     $A56-A11*A22*A34*A45*A56

      COFACTOR(1,4) = -A26*A35*A43*A52*A61+A25*A36*A43*A52*A61+A26*A33*
     $A45*A52*A61-A23*A36*A45*A52*A61-A25*A33*A46*A52*A61+A23*A35*A46*
     $A52*A61+A26*A35*A42*A53*A61-A25*A36*A42*A53*A61-A26*A32*A45*A53*
     $A61+A22*A36*A45*A53*A61+A25*A32*A46*A53*A61-A22*A35*A46*A53*A61-
     $A26*A33*A42*A55*A61+A23*A36*A42*A55*A61+A26*A32*A43*A55*A61-A22*
     $A36*A43*A55*A61-A23*A32*A46*A55*A61+A22*A33*A46*A55*A61+A25*A33*
     $A42*A56*A61-A23*A35*A42*A56*A61-A25*A32*A43*A56*A61+A22*A35*A43*
     $A56*A61+A23*A32*A45*A56*A61-A22*A33*A45*A56*A61+A26*A35*A43*A51*
     $A62-A25*A36*A43*A51*A62-A26*A33*A45*A51*A62+A23*A36*A45*A51*A62+
     $A25*A33*A46*A51*A62-A23*A35*A46*A51*A62-A26*A35*A41*A53*A62+A25*
     $A36*A41*A53*A62+A26*A31*A45*A53*A62-A21*A36*A45*A53*A62-A25*A31*
     $A46*A53*A62+A21*A35*A46*A53*A62+A26*A33*A41*A55*A62-A23*A36*A41*
     $A55*A62-A26*A31*A43*A55*A62+A21*A36*A43*A55*A62+A23*A31*A46*A55*
     $A62-A21*A33*A46*A55*A62-A25*A33*A41*A56*A62+A23*A35*A41*A56*A62+
     $A25*A31*A43*A56*A62-A21*A35*A43*A56*A62-A23*A31*A45*A56*A62+A21*
     $A33*A45*A56*A62-A26*A35*A42*A51*A63+A25*A36*A42*A51*A63+A26*A32*
     $A45*A51*A63-A22*A36*A45*A51*A63-A25*A32*A46*A51*A63+A22*A35*A46*
     $A51*A63+A26*A35*A41*A52*A63-A25*A36*A41*A52*A63-A26*A31*A45*A52*
     $A63+A21*A36*A45*A52*A63+A25*A31*A46*A52*A63-A21*A35*A46*A52*A63-
     $A26*A32*A41*A55*A63+A22*A36*A41*A55*A63+A26*A31*A42*A55*A63-A21*
     $A36*A42*A55*A63-A22*A31*A46*A55*A63+A21*A32*A46*A55*A63+A25*A32*
     $A41*A56*A63-A22*A35*A41*A56*A63-A25*A31*A42*A56*A63+A21*A35*A42*
     $A56*A63+A22*A31*A45*A56*A63-A21*A32*A45*A56*A63+A26*A33*A42*A51*
     $A65-A23*A36*A42*A51*A65-A26*A32*A43*A51*A65+A22*A36*A43*A51*A65+
     $A23*A32*A46*A51*A65-A22*A33*A46*A51*A65-A26*A33*A41*A52*A65+A23*
     $A36*A41*A52*A65+A26*A31*A43*A52*A65-A21*A36*A43*A52*A65-A23*A31*
     $A46*A52*A65+A21*A33*A46*A52*A65+A26*A32*A41*A53*A65-A22*A36*A41*
     $A53*A65-A26*A31*A42*A53*A65+A21*A36*A42*A53*A65+A22*A31*A46*A53*
     $A65-A21*A32*A46*A53*A65-A23*A32*A41*A56*A65+A22*A33*A41*A56*A65+
     $A23*A31*A42*A56*A65-A21*A33*A42*A56*A65-A22*A31*A43*A56*A65+A21*
     $A32*A43*A56*A65-A25*A33*A42*A51*A66+A23*A35*A42*A51*A66+A25*A32*
     $A43*A51*A66-A22*A35*A43*A51*A66-A23*A32*A45*A51*A66+A22*A33*A45*
     $A51*A66+A25*A33*A41*A52*A66-A23*A35*A41*A52*A66-A25*A31*A43*A52*
     $A66+A21*A35*A43*A52*A66+A23*A31*A45*A52*A66-A21*A33*A45*A52*A66-
     $A25*A32*A41*A53*A66+A22*A35*A41*A53*A66+A25*A31*A42*A53*A66-A21*
     $A35*A42*A53*A66-A22*A31*A45*A53*A66+A21*A32*A45*A53*A66+A23*A32*
     $A41*A55*A66-A22*A33*A41*A55*A66-A23*A31*A42*A55*A66+A21*A33*A42*
     $A55*A66+A22*A31*A43*A55*A66-A21*A32*A43*A55*A66

      COFACTOR(2,4) = A16*A35*A43*A52*A61-A15*A36*A43*A52*
     $A61-A16*A33*A45*A52*A61+A13*A36*A45*A52*A61+A15*A33*A46*A52*A61-
     $A13*A35*A46*A52*A61-A16*A35*A42*A53*A61+A15*A36*A42*A53*A61+A16*
     $A32*A45*A53*A61-A12*A36*A45*A53*A61-A15*A32*A46*A53*A61+A12*A35*
     $A46*A53*A61+A16*A33*A42*A55*A61-A13*A36*A42*A55*A61-A16*A32*A43*
     $A55*A61+A12*A36*A43*A55*A61+A13*A32*A46*A55*A61-A12*A33*A46*A55*
     $A61-A15*A33*A42*A56*A61+A13*A35*A42*A56*A61+A15*A32*A43*A56*A61-
     $A12*A35*A43*A56*A61-A13*A32*A45*A56*A61+A12*A33*A45*A56*A61-A16*
     $A35*A43*A51*A62+A15*A36*A43*A51*A62+A16*A33*A45*A51*A62-A13*A36*
     $A45*A51*A62-A15*A33*A46*A51*A62+A13*A35*A46*A51*A62+A16*A35*A41*
     $A53*A62-A15*A36*A41*A53*A62-A16*A31*A45*A53*A62+A11*A36*A45*A53*
     $A62+A15*A31*A46*A53*A62-A11*A35*A46*A53*A62-A16*A33*A41*A55*A62+
     $A13*A36*A41*A55*A62+A16*A31*A43*A55*A62-A11*A36*A43*A55*A62-A13*
     $A31*A46*A55*A62+A11*A33*A46*A55*A62+A15*A33*A41*A56*A62-A13*A35*
     $A41*A56*A62-A15*A31*A43*A56*A62+A11*A35*A43*A56*A62+A13*A31*A45*
     $A56*A62-A11*A33*A45*A56*A62+A16*A35*A42*A51*A63-A15*A36*A42*A51*
     $A63-A16*A32*A45*A51*A63+A12*A36*A45*A51*A63+A15*A32*A46*A51*A63-
     $A12*A35*A46*A51*A63-A16*A35*A41*A52*A63+A15*A36*A41*A52*A63+A16*
     $A31*A45*A52*A63-A11*A36*A45*A52*A63-A15*A31*A46*A52*A63+A11*A35*
     $A46*A52*A63+A16*A32*A41*A55*A63-A12*A36*A41*A55*A63-A16*A31*A42*
     $A55*A63+A11*A36*A42*A55*A63+A12*A31*A46*A55*A63-A11*A32*A46*A55*
     $A63-A15*A32*A41*A56*A63+A12*A35*A41*A56*A63+A15*A31*A42*A56*A63-
     $A11*A35*A42*A56*A63-A12*A31*A45*A56*A63+A11*A32*A45*A56*A63-A16*
     $A33*A42*A51*A65+A13*A36*A42*A51*A65+A16*A32*A43*A51*A65-A12*A36*
     $A43*A51*A65-A13*A32*A46*A51*A65+A12*A33*A46*A51*A65+A16*A33*A41*
     $A52*A65-A13*A36*A41*A52*A65-A16*A31*A43*A52*A65+A11*A36*A43*A52*
     $A65+A13*A31*A46*A52*A65-A11*A33*A46*A52*A65-A16*A32*A41*A53*A65+
     $A12*A36*A41*A53*A65+A16*A31*A42*A53*A65-A11*A36*A42*A53*A65-A12*
     $A31*A46*A53*A65+A11*A32*A46*A53*A65+A13*A32*A41*A56*A65-A12*A33*
     $A41*A56*A65-A13*A31*A42*A56*A65+A11*A33*A42*A56*A65+A12*A31*A43*
     $A56*A65-A11*A32*A43*A56*A65+A15*A33*A42*A51*A66-A13*A35*A42*A51*
     $A66-A15*A32*A43*A51*A66+A12*A35*A43*A51*A66+A13*A32*A45*A51*A66-
     $A12*A33*A45*A51*A66-A15*A33*A41*A52*A66+A13*A35*A41*A52*A66+A15*
     $A31*A43*A52*A66-A11*A35*A43*A52*A66-A13*A31*A45*A52*A66+A11*A33*
     $A45*A52*A66+A15*A32*A41*A53*A66-A12*A35*A41*A53*A66-A15*A31*A42*
     $A53*A66+A11*A35*A42*A53*A66+A12*A31*A45*A53*A66-A11*A32*A45*A53*
     $A66-A13*A32*A41*A55*A66+A12*A33*A41*A55*A66+A13*A31*A42*A55*A66-
     $A11*A33*A42*A55*A66-A12*A31*A43*A55*A66+A11*A32*A43*A55*A66

      COFACTOR(3,4) = -A16*A25*A43*A52*
     $A61+A15*A26*A43*A52*A61+A16*A23*A45*A52*A61-A13*A26*A45*A52*A61-
     $A15*A23*A46*A52*A61+A13*A25*A46*A52*A61+A16*A25*A42*A53*A61-A15*
     $A26*A42*A53*A61-A16*A22*A45*A53*A61+A12*A26*A45*A53*A61+A15*A22*
     $A46*A53*A61-A12*A25*A46*A53*A61-A16*A23*A42*A55*A61+A13*A26*A42*
     $A55*A61+A16*A22*A43*A55*A61-A12*A26*A43*A55*A61-A13*A22*A46*A55*
     $A61+A12*A23*A46*A55*A61+A15*A23*A42*A56*A61-A13*A25*A42*A56*A61-
     $A15*A22*A43*A56*A61+A12*A25*A43*A56*A61+A13*A22*A45*A56*A61-A12*
     $A23*A45*A56*A61+A16*A25*A43*A51*A62-A15*A26*A43*A51*A62-A16*A23*
     $A45*A51*A62+A13*A26*A45*A51*A62+A15*A23*A46*A51*A62-A13*A25*A46*
     $A51*A62-A16*A25*A41*A53*A62+A15*A26*A41*A53*A62+A16*A21*A45*A53*
     $A62-A11*A26*A45*A53*A62-A15*A21*A46*A53*A62+A11*A25*A46*A53*A62+
     $A16*A23*A41*A55*A62-A13*A26*A41*A55*A62-A16*A21*A43*A55*A62+A11*
     $A26*A43*A55*A62+A13*A21*A46*A55*A62-A11*A23*A46*A55*A62-A15*A23*
     $A41*A56*A62+A13*A25*A41*A56*A62+A15*A21*A43*A56*A62-A11*A25*A43*
     $A56*A62-A13*A21*A45*A56*A62+A11*A23*A45*A56*A62-A16*A25*A42*A51*
     $A63+A15*A26*A42*A51*A63+A16*A22*A45*A51*A63-A12*A26*A45*A51*A63-
     $A15*A22*A46*A51*A63+A12*A25*A46*A51*A63+A16*A25*A41*A52*A63-A15*
     $A26*A41*A52*A63-A16*A21*A45*A52*A63+A11*A26*A45*A52*A63+A15*A21*
     $A46*A52*A63-A11*A25*A46*A52*A63-A16*A22*A41*A55*A63+A12*A26*A41*
     $A55*A63+A16*A21*A42*A55*A63-A11*A26*A42*A55*A63-A12*A21*A46*A55*
     $A63+A11*A22*A46*A55*A63+A15*A22*A41*A56*A63-A12*A25*A41*A56*A63-
     $A15*A21*A42*A56*A63+A11*A25*A42*A56*A63+A12*A21*A45*A56*A63-A11*
     $A22*A45*A56*A63+A16*A23*A42*A51*A65-A13*A26*A42*A51*A65-A16*A22*
     $A43*A51*A65+A12*A26*A43*A51*A65+A13*A22*A46*A51*A65-A12*A23*A46*
     $A51*A65-A16*A23*A41*A52*A65+A13*A26*A41*A52*A65+A16*A21*A43*A52*
     $A65-A11*A26*A43*A52*A65-A13*A21*A46*A52*A65+A11*A23*A46*A52*A65+
     $A16*A22*A41*A53*A65-A12*A26*A41*A53*A65-A16*A21*A42*A53*A65+A11*
     $A26*A42*A53*A65+A12*A21*A46*A53*A65-A11*A22*A46*A53*A65-A13*A22*
     $A41*A56*A65+A12*A23*A41*A56*A65+A13*A21*A42*A56*A65-A11*A23*A42*
     $A56*A65-A12*A21*A43*A56*A65+A11*A22*A43*A56*A65-A15*A23*A42*A51*
     $A66+A13*A25*A42*A51*A66+A15*A22*A43*A51*A66-A12*A25*A43*A51*A66-
     $A13*A22*A45*A51*A66+A12*A23*A45*A51*A66+A15*A23*A41*A52*A66-A13*
     $A25*A41*A52*A66-A15*A21*A43*A52*A66+A11*A25*A43*A52*A66+A13*A21*
     $A45*A52*A66-A11*A23*A45*A52*A66-A15*A22*A41*A53*A66+A12*A25*A41*
     $A53*A66+A15*A21*A42*A53*A66-A11*A25*A42*A53*A66-A12*A21*A45*A53*
     $A66+A11*A22*A45*A53*A66+A13*A22*A41*A55*A66-A12*A23*A41*A55*A66-
     $A13*A21*A42*A55*A66+A11*A23*A42*A55*A66+A12*A21*A43*A55*A66-A11*
     $A22*A43*A55*A66

      COFACTOR(4,4) = A16*A25*A33*A52*A61-A15*A26*A33*A52*A61-A16*A23*
     $A35*A52*A61+A13*A26*A35*A52*A61+A15*A23*A36*A52*A61-A13*A25*A36*
     $A52*A61-A16*A25*A32*A53*A61+A15*A26*A32*A53*A61+A16*A22*A35*A53*
     $A61-A12*A26*A35*A53*A61-A15*A22*A36*A53*A61+A12*A25*A36*A53*A61+
     $A16*A23*A32*A55*A61-A13*A26*A32*A55*A61-A16*A22*A33*A55*A61+A12*
     $A26*A33*A55*A61+A13*A22*A36*A55*A61-A12*A23*A36*A55*A61-A15*A23*
     $A32*A56*A61+A13*A25*A32*A56*A61+A15*A22*A33*A56*A61-A12*A25*A33*
     $A56*A61-A13*A22*A35*A56*A61+A12*A23*A35*A56*A61-A16*A25*A33*A51*
     $A62+A15*A26*A33*A51*A62+A16*A23*A35*A51*A62-A13*A26*A35*A51*A62-
     $A15*A23*A36*A51*A62+A13*A25*A36*A51*A62+A16*A25*A31*A53*A62-A15*
     $A26*A31*A53*A62-A16*A21*A35*A53*A62+A11*A26*A35*A53*A62+A15*A21*
     $A36*A53*A62-A11*A25*A36*A53*A62-A16*A23*A31*A55*A62+A13*A26*A31*
     $A55*A62+A16*A21*A33*A55*A62-A11*A26*A33*A55*A62-A13*A21*A36*A55*
     $A62+A11*A23*A36*A55*A62+A15*A23*A31*A56*A62-A13*A25*A31*A56*A62-
     $A15*A21*A33*A56*A62+A11*A25*A33*A56*A62+A13*A21*A35*A56*A62-A11*
     $A23*A35*A56*A62+A16*A25*A32*A51*A63-A15*A26*A32*A51*A63-A16*A22*
     $A35*A51*A63+A12*A26*A35*A51*A63+A15*A22*A36*A51*A63-A12*A25*A36*
     $A51*A63-A16*A25*A31*A52*A63+A15*A26*A31*A52*A63+A16*A21*A35*A52*
     $A63-A11*A26*A35*A52*A63-A15*A21*A36*A52*A63+A11*A25*A36*A52*A63+
     $A16*A22*A31*A55*A63-A12*A26*A31*A55*A63-A16*A21*A32*A55*A63+A11*
     $A26*A32*A55*A63+A12*A21*A36*A55*A63-A11*A22*A36*A55*A63-A15*A22*
     $A31*A56*A63+A12*A25*A31*A56*A63+A15*A21*A32*A56*A63-A11*A25*A32*
     $A56*A63-A12*A21*A35*A56*A63+A11*A22*A35*A56*A63-A16*A23*A32*A51*
     $A65+A13*A26*A32*A51*A65+A16*A22*A33*A51*A65-A12*A26*A33*A51*A65-
     $A13*A22*A36*A51*A65+A12*A23*A36*A51*A65+A16*A23*A31*A52*A65-A13*
     $A26*A31*A52*A65-A16*A21*A33*A52*A65+A11*A26*A33*A52*A65+A13*A21*
     $A36*A52*A65-A11*A23*A36*A52*A65-A16*A22*A31*A53*A65+A12*A26*A31*
     $A53*A65+A16*A21*A32*A53*A65-A11*A26*A32*A53*A65-A12*A21*A36*A53*
     $A65+A11*A22*A36*A53*A65+A13*A22*A31*A56*A65-A12*A23*A31*A56*A65-
     $A13*A21*A32*A56*A65+A11*A23*A32*A56*A65+A12*A21*A33*A56*A65-A11*
     $A22*A33*A56*A65+A15*A23*A32*A51*A66-A13*A25*A32*A51*A66-A15*A22*
     $A33*A51*A66+A12*A25*A33*A51*A66+A13*A22*A35*A51*A66-A12*A23*A35*
     $A51*A66-A15*A23*A31*A52*A66+A13*A25*A31*A52*A66+A15*A21*A33*A52*
     $A66-A11*A25*A33*A52*A66-A13*A21*A35*A52*A66+A11*A23*A35*A52*A66+
     $A15*A22*A31*A53*A66-A12*A25*A31*A53*A66-A15*A21*A32*A53*A66+A11*
     $A25*A32*A53*A66+A12*A21*A35*A53*A66-A11*A22*A35*A53*A66-A13*A22*
     $A31*A55*A66+A12*A23*A31*A55*A66+A13*A21*A32*A55*A66-A11*A23*A32*
     $A55*A66-A12*A21*A33*A55*A66+A11*A22*A33*A55*A66

      COFACTOR(5,4) = -A16*A25*A33*A42*A61+A15*A26*A33*A42*A61+A16*A23*
     $A35*A42*A61-A13*A26*A35*A42*A61-A15*A23*A36*A42*A61+A13*A25*A36*
     $A42*A61+A16*A25*A32*A43*A61-A15*A26*A32*A43*A61-A16*A22*A35*A43*
     $A61+A12*A26*A35*A43*A61+A15*A22*A36*A43*A61-A12*A25*A36*A43*A61-
     $A16*A23*A32*A45*A61+A13*A26*A32*A45*A61+A16*A22*A33*A45*A61-A12*
     $A26*A33*A45*A61-A13*A22*A36*A45*A61+A12*A23*A36*A45*A61+A15*A23*
     $A32*A46*A61-A13*A25*A32*A46*A61-A15*A22*A33*A46*A61+A12*A25*A33*
     $A46*A61+A13*A22*A35*A46*A61-A12*A23*A35*A46*A61+A16*A25*A33*A41*
     $A62-A15*A26*A33*A41*A62-A16*A23*A35*A41*A62+A13*A26*A35*A41*A62+
     $A15*A23*A36*A41*A62-A13*A25*A36*A41*A62-A16*A25*A31*A43*A62+A15*
     $A26*A31*A43*A62+A16*A21*A35*A43*A62-A11*A26*A35*A43*A62-A15*A21*
     $A36*A43*A62+A11*A25*A36*A43*A62+A16*A23*A31*A45*A62-A13*A26*A31*
     $A45*A62-A16*A21*A33*A45*A62+A11*A26*A33*A45*A62+A13*A21*A36*A45*
     $A62-A11*A23*A36*A45*A62-A15*A23*A31*A46*A62+A13*A25*A31*A46*A62+
     $A15*A21*A33*A46*A62-A11*A25*A33*A46*A62-A13*A21*A35*A46*A62+A11*
     $A23*A35*A46*A62-A16*A25*A32*A41*A63+A15*A26*A32*A41*A63+A16*A22*
     $A35*A41*A63-A12*A26*A35*A41*A63-A15*A22*A36*A41*A63+A12*A25*A36*
     $A41*A63+A16*A25*A31*A42*A63-A15*A26*A31*A42*A63-A16*A21*A35*A42*
     $A63+A11*A26*A35*A42*A63+A15*A21*A36*A42*A63-A11*A25*A36*A42*A63-
     $A16*A22*A31*A45*A63+A12*A26*A31*A45*A63+A16*A21*A32*A45*A63-A11*
     $A26*A32*A45*A63-A12*A21*A36*A45*A63+A11*A22*A36*A45*A63+A15*A22*
     $A31*A46*A63-A12*A25*A31*A46*A63-A15*A21*A32*A46*A63+A11*A25*A32*
     $A46*A63+A12*A21*A35*A46*A63-A11*A22*A35*A46*A63+A16*A23*A32*A41*
     $A65-A13*A26*A32*A41*A65-A16*A22*A33*A41*A65+A12*A26*A33*A41*A65+
     $A13*A22*A36*A41*A65-A12*A23*A36*A41*A65-A16*A23*A31*A42*A65+A13*
     $A26*A31*A42*A65+A16*A21*A33*A42*A65-A11*A26*A33*A42*A65-A13*A21*
     $A36*A42*A65+A11*A23*A36*A42*A65+A16*A22*A31*A43*A65-A12*A26*A31*
     $A43*A65-A16*A21*A32*A43*A65+A11*A26*A32*A43*A65+A12*A21*A36*A43*
     $A65-A11*A22*A36*A43*A65-A13*A22*A31*A46*A65+A12*A23*A31*A46*A65+
     $A13*A21*A32*A46*A65-A11*A23*A32*A46*A65-A12*A21*A33*A46*A65+A11*
     $A22*A33*A46*A65-A15*A23*A32*A41*A66+A13*A25*A32*A41*A66+A15*A22*
     $A33*A41*A66-A12*A25*A33*A41*A66-A13*A22*A35*A41*A66+A12*A23*A35*
     $A41*A66+A15*A23*A31*A42*A66-A13*A25*A31*A42*A66-A15*A21*A33*A42*
     $A66+A11*A25*A33*A42*A66+A13*A21*A35*A42*A66-A11*A23*A35*A42*A66-
     $A15*A22*A31*A43*A66+A12*A25*A31*A43*A66+A15*A21*A32*A43*A66-A11*
     $A25*A32*A43*A66-A12*A21*A35*A43*A66+A11*A22*A35*A43*A66+A13*A22*
     $A31*A45*A66-A12*A23*A31*A45*A66-A13*A21*A32*A45*A66+A11*A23*A32*
     $A45*A66+A12*A21*A33*A45*A66-A11*A22*A33*A45*A66

      COFACTOR(6,4) = A16*A25*A33*A42*A51-A15*A26*
     $A33*A42*A51-A16*A23*A35*A42*A51+A13*A26*A35*A42*A51+A15*A23*A36*
     $A42*A51-A13*A25*A36*A42*A51-A16*A25*A32*A43*A51+A15*A26*A32*A43*
     $A51+A16*A22*A35*A43*A51-A12*A26*A35*A43*A51-A15*A22*A36*A43*A51+
     $A12*A25*A36*A43*A51+A16*A23*A32*A45*A51-A13*A26*A32*A45*A51-A16*
     $A22*A33*A45*A51+A12*A26*A33*A45*A51+A13*A22*A36*A45*A51-A12*A23*
     $A36*A45*A51-A15*A23*A32*A46*A51+A13*A25*A32*A46*A51+A15*A22*A33*
     $A46*A51-A12*A25*A33*A46*A51-A13*A22*A35*A46*A51+A12*A23*A35*A46*
     $A51-A16*A25*A33*A41*A52+A15*A26*A33*A41*A52+A16*A23*A35*A41*A52-
     $A13*A26*A35*A41*A52-A15*A23*A36*A41*A52+A13*A25*A36*A41*A52+A16*
     $A25*A31*A43*A52-A15*A26*A31*A43*A52-A16*A21*A35*A43*A52+A11*A26*
     $A35*A43*A52+A15*A21*A36*A43*A52-A11*A25*A36*A43*A52-A16*A23*A31*
     $A45*A52+A13*A26*A31*A45*A52+A16*A21*A33*A45*A52-A11*A26*A33*A45*
     $A52-A13*A21*A36*A45*A52+A11*A23*A36*A45*A52+A15*A23*A31*A46*A52-
     $A13*A25*A31*A46*A52-A15*A21*A33*A46*A52+A11*A25*A33*A46*A52+A13*
     $A21*A35*A46*A52-A11*A23*A35*A46*A52+A16*A25*A32*A41*A53-A15*A26*
     $A32*A41*A53-A16*A22*A35*A41*A53+A12*A26*A35*A41*A53+A15*A22*A36*
     $A41*A53-A12*A25*A36*A41*A53-A16*A25*A31*A42*A53+A15*A26*A31*A42*
     $A53+A16*A21*A35*A42*A53-A11*A26*A35*A42*A53-A15*A21*A36*A42*A53+
     $A11*A25*A36*A42*A53+A16*A22*A31*A45*A53-A12*A26*A31*A45*A53-A16*
     $A21*A32*A45*A53+A11*A26*A32*A45*A53+A12*A21*A36*A45*A53-A11*A22*
     $A36*A45*A53-A15*A22*A31*A46*A53+A12*A25*A31*A46*A53+A15*A21*A32*
     $A46*A53-A11*A25*A32*A46*A53-A12*A21*A35*A46*A53+A11*A22*A35*A46*
     $A53-A16*A23*A32*A41*A55+A13*A26*A32*A41*A55+A16*A22*A33*A41*A55-
     $A12*A26*A33*A41*A55-A13*A22*A36*A41*A55+A12*A23*A36*A41*A55+A16*
     $A23*A31*A42*A55-A13*A26*A31*A42*A55-A16*A21*A33*A42*A55+A11*A26*
     $A33*A42*A55+A13*A21*A36*A42*A55-A11*A23*A36*A42*A55-A16*A22*A31*
     $A43*A55+A12*A26*A31*A43*A55+A16*A21*A32*A43*A55-A11*A26*A32*A43*
     $A55-A12*A21*A36*A43*A55+A11*A22*A36*A43*A55+A13*A22*A31*A46*A55-
     $A12*A23*A31*A46*A55-A13*A21*A32*A46*A55+A11*A23*A32*A46*A55+A12*
     $A21*A33*A46*A55-A11*A22*A33*A46*A55+A15*A23*A32*A41*A56-A13*A25*
     $A32*A41*A56-A15*A22*A33*A41*A56+A12*A25*A33*A41*A56+A13*A22*A35*
     $A41*A56-A12*A23*A35*A41*A56-A15*A23*A31*A42*A56+A13*A25*A31*A42*
     $A56+A15*A21*A33*A42*A56-A11*A25*A33*A42*A56-A13*A21*A35*A42*A56+
     $A11*A23*A35*A42*A56+A15*A22*A31*A43*A56-A12*A25*A31*A43*A56-A15*
     $A21*A32*A43*A56+A11*A25*A32*A43*A56+A12*A21*A35*A43*A56-A11*A22*
     $A35*A43*A56-A13*A22*A31*A45*A56+A12*A23*A31*A45*A56+A13*A21*A32*
     $A45*A56-A11*A23*A32*A45*A56-A12*A21*A33*A45*A56+A11*A22*A33*A45*
     $A56

      COFACTOR(1,5) = A26*A34*A43*A52*A61-A24*A36*A43*A52*A61-A26*A33*
     $A44*A52*A61+A23*A36*A44*A52*A61+A24*A33*A46*A52*A61-A23*A34*A46*
     $A52*A61-A26*A34*A42*A53*A61+A24*A36*A42*A53*A61+A26*A32*A44*A53*
     $A61-A22*A36*A44*A53*A61-A24*A32*A46*A53*A61+A22*A34*A46*A53*A61+
     $A26*A33*A42*A54*A61-A23*A36*A42*A54*A61-A26*A32*A43*A54*A61+A22*
     $A36*A43*A54*A61+A23*A32*A46*A54*A61-A22*A33*A46*A54*A61-A24*A33*
     $A42*A56*A61+A23*A34*A42*A56*A61+A24*A32*A43*A56*A61-A22*A34*A43*
     $A56*A61-A23*A32*A44*A56*A61+A22*A33*A44*A56*A61-A26*A34*A43*A51*
     $A62+A24*A36*A43*A51*A62+A26*A33*A44*A51*A62-A23*A36*A44*A51*A62-
     $A24*A33*A46*A51*A62+A23*A34*A46*A51*A62+A26*A34*A41*A53*A62-A24*
     $A36*A41*A53*A62-A26*A31*A44*A53*A62+A21*A36*A44*A53*A62+A24*A31*
     $A46*A53*A62-A21*A34*A46*A53*A62-A26*A33*A41*A54*A62+A23*A36*A41*
     $A54*A62+A26*A31*A43*A54*A62-A21*A36*A43*A54*A62-A23*A31*A46*A54*
     $A62+A21*A33*A46*A54*A62+A24*A33*A41*A56*A62-A23*A34*A41*A56*A62-
     $A24*A31*A43*A56*A62+A21*A34*A43*A56*A62+A23*A31*A44*A56*A62-A21*
     $A33*A44*A56*A62+A26*A34*A42*A51*A63-A24*A36*A42*A51*A63-A26*A32*
     $A44*A51*A63+A22*A36*A44*A51*A63+A24*A32*A46*A51*A63-A22*A34*A46*
     $A51*A63-A26*A34*A41*A52*A63+A24*A36*A41*A52*A63+A26*A31*A44*A52*
     $A63-A21*A36*A44*A52*A63-A24*A31*A46*A52*A63+A21*A34*A46*A52*A63+
     $A26*A32*A41*A54*A63-A22*A36*A41*A54*A63-A26*A31*A42*A54*A63+A21*
     $A36*A42*A54*A63+A22*A31*A46*A54*A63-A21*A32*A46*A54*A63-A24*A32*
     $A41*A56*A63+A22*A34*A41*A56*A63+A24*A31*A42*A56*A63-A21*A34*A42*
     $A56*A63-A22*A31*A44*A56*A63+A21*A32*A44*A56*A63-A26*A33*A42*A51*
     $A64+A23*A36*A42*A51*A64+A26*A32*A43*A51*A64-A22*A36*A43*A51*A64-
     $A23*A32*A46*A51*A64+A22*A33*A46*A51*A64+A26*A33*A41*A52*A64-A23*
     $A36*A41*A52*A64-A26*A31*A43*A52*A64+A21*A36*A43*A52*A64+A23*A31*
     $A46*A52*A64-A21*A33*A46*A52*A64-A26*A32*A41*A53*A64+A22*A36*A41*
     $A53*A64+A26*A31*A42*A53*A64-A21*A36*A42*A53*A64-A22*A31*A46*A53*
     $A64+A21*A32*A46*A53*A64+A23*A32*A41*A56*A64-A22*A33*A41*A56*A64-
     $A23*A31*A42*A56*A64+A21*A33*A42*A56*A64+A22*A31*A43*A56*A64-A21*
     $A32*A43*A56*A64+A24*A33*A42*A51*A66-A23*A34*A42*A51*A66-A24*A32*
     $A43*A51*A66+A22*A34*A43*A51*A66+A23*A32*A44*A51*A66-A22*A33*A44*
     $A51*A66-A24*A33*A41*A52*A66+A23*A34*A41*A52*A66+A24*A31*A43*A52*
     $A66-A21*A34*A43*A52*A66-A23*A31*A44*A52*A66+A21*A33*A44*A52*A66+
     $A24*A32*A41*A53*A66-A22*A34*A41*A53*A66-A24*A31*A42*A53*A66+A21*
     $A34*A42*A53*A66+A22*A31*A44*A53*A66-A21*A32*A44*A53*A66-A23*A32*
     $A41*A54*A66+A22*A33*A41*A54*A66+A23*A31*A42*A54*A66-A21*A33*A42*
     $A54*A66-A22*A31*A43*A54*A66+A21*A32*A43*A54*A66

      COFACTOR(2,5) = -A16*A34*A43*A52*A61+A14*A36*A43*A52*A61+A16*A33*
     $A44*A52*A61-A13*A36*A44*A52*A61-A14*A33*A46*A52*A61+A13*A34*A46*
     $A52*A61+A16*A34*A42*A53*A61-A14*A36*A42*A53*A61-A16*A32*A44*A53*
     $A61+A12*A36*A44*A53*A61+A14*A32*A46*A53*A61-A12*A34*A46*A53*A61-
     $A16*A33*A42*A54*A61+A13*A36*A42*A54*A61+A16*A32*A43*A54*A61-A12*
     $A36*A43*A54*A61-A13*A32*A46*A54*A61+A12*A33*A46*A54*A61+A14*A33*
     $A42*A56*A61-A13*A34*A42*A56*A61-A14*A32*A43*A56*A61+A12*A34*A43*
     $A56*A61+A13*A32*A44*A56*A61-A12*A33*A44*A56*A61+A16*A34*A43*A51*
     $A62-A14*A36*A43*A51*A62-A16*A33*A44*A51*A62+A13*A36*A44*A51*A62+
     $A14*A33*A46*A51*A62-A13*A34*A46*A51*A62-A16*A34*A41*A53*A62+A14*
     $A36*A41*A53*A62+A16*A31*A44*A53*A62-A11*A36*A44*A53*A62-A14*A31*
     $A46*A53*A62+A11*A34*A46*A53*A62+A16*A33*A41*A54*A62-A13*A36*A41*
     $A54*A62-A16*A31*A43*A54*A62+A11*A36*A43*A54*A62+A13*A31*A46*A54*
     $A62-A11*A33*A46*A54*A62-A14*A33*A41*A56*A62+A13*A34*A41*A56*A62+
     $A14*A31*A43*A56*A62-A11*A34*A43*A56*A62-A13*A31*A44*A56*A62+A11*
     $A33*A44*A56*A62-A16*A34*A42*A51*A63+A14*A36*A42*A51*A63+A16*A32*
     $A44*A51*A63-A12*A36*A44*A51*A63-A14*A32*A46*A51*A63+A12*A34*A46*
     $A51*A63+A16*A34*A41*A52*A63-A14*A36*A41*A52*A63-A16*A31*A44*A52*
     $A63+A11*A36*A44*A52*A63+A14*A31*A46*A52*A63-A11*A34*A46*A52*A63-
     $A16*A32*A41*A54*A63+A12*A36*A41*A54*A63+A16*A31*A42*A54*A63-A11*
     $A36*A42*A54*A63-A12*A31*A46*A54*A63+A11*A32*A46*A54*A63+A14*A32*
     $A41*A56*A63-A12*A34*A41*A56*A63-A14*A31*A42*A56*A63+A11*A34*A42*
     $A56*A63+A12*A31*A44*A56*A63-A11*A32*A44*A56*A63+A16*A33*A42*A51*
     $A64-A13*A36*A42*A51*A64-A16*A32*A43*A51*A64+A12*A36*A43*A51*A64+
     $A13*A32*A46*A51*A64-A12*A33*A46*A51*A64-A16*A33*A41*A52*A64+A13*
     $A36*A41*A52*A64+A16*A31*A43*A52*A64-A11*A36*A43*A52*A64-A13*A31*
     $A46*A52*A64+A11*A33*A46*A52*A64+A16*A32*A41*A53*A64-A12*A36*A41*
     $A53*A64-A16*A31*A42*A53*A64+A11*A36*A42*A53*A64+A12*A31*A46*A53*
     $A64-A11*A32*A46*A53*A64-A13*A32*A41*A56*A64+A12*A33*A41*A56*A64+
     $A13*A31*A42*A56*A64-A11*A33*A42*A56*A64-A12*A31*A43*A56*A64+A11*
     $A32*A43*A56*A64-A14*A33*A42*A51*A66+A13*A34*A42*A51*A66+A14*A32*
     $A43*A51*A66-A12*A34*A43*A51*A66-A13*A32*A44*A51*A66+A12*A33*A44*
     $A51*A66+A14*A33*A41*A52*A66-A13*A34*A41*A52*A66-A14*A31*A43*A52*
     $A66+A11*A34*A43*A52*A66+A13*A31*A44*A52*A66-A11*A33*A44*A52*A66-
     $A14*A32*A41*A53*A66+A12*A34*A41*A53*A66+A14*A31*A42*A53*A66-A11*
     $A34*A42*A53*A66-A12*A31*A44*A53*A66+A11*A32*A44*A53*A66+A13*A32*
     $A41*A54*A66-A12*A33*A41*A54*A66-A13*A31*A42*A54*A66+A11*A33*A42*
     $A54*A66+A12*A31*A43*A54*A66-A11*A32*A43*A54*A66

      COFACTOR(3,5) = A16*A24*A43*A52*A61-A14*A26*A43*A52*
     $A61-A16*A23*A44*A52*A61+A13*A26*A44*A52*A61+A14*A23*A46*A52*A61-
     $A13*A24*A46*A52*A61-A16*A24*A42*A53*A61+A14*A26*A42*A53*A61+A16*
     $A22*A44*A53*A61-A12*A26*A44*A53*A61-A14*A22*A46*A53*A61+A12*A24*
     $A46*A53*A61+A16*A23*A42*A54*A61-A13*A26*A42*A54*A61-A16*A22*A43*
     $A54*A61+A12*A26*A43*A54*A61+A13*A22*A46*A54*A61-A12*A23*A46*A54*
     $A61-A14*A23*A42*A56*A61+A13*A24*A42*A56*A61+A14*A22*A43*A56*A61-
     $A12*A24*A43*A56*A61-A13*A22*A44*A56*A61+A12*A23*A44*A56*A61-A16*
     $A24*A43*A51*A62+A14*A26*A43*A51*A62+A16*A23*A44*A51*A62-A13*A26*
     $A44*A51*A62-A14*A23*A46*A51*A62+A13*A24*A46*A51*A62+A16*A24*A41*
     $A53*A62-A14*A26*A41*A53*A62-A16*A21*A44*A53*A62+A11*A26*A44*A53*
     $A62+A14*A21*A46*A53*A62-A11*A24*A46*A53*A62-A16*A23*A41*A54*A62+
     $A13*A26*A41*A54*A62+A16*A21*A43*A54*A62-A11*A26*A43*A54*A62-A13*
     $A21*A46*A54*A62+A11*A23*A46*A54*A62+A14*A23*A41*A56*A62-A13*A24*
     $A41*A56*A62-A14*A21*A43*A56*A62+A11*A24*A43*A56*A62+A13*A21*A44*
     $A56*A62-A11*A23*A44*A56*A62+A16*A24*A42*A51*A63-A14*A26*A42*A51*
     $A63-A16*A22*A44*A51*A63+A12*A26*A44*A51*A63+A14*A22*A46*A51*A63-
     $A12*A24*A46*A51*A63-A16*A24*A41*A52*A63+A14*A26*A41*A52*A63+A16*
     $A21*A44*A52*A63-A11*A26*A44*A52*A63-A14*A21*A46*A52*A63+A11*A24*
     $A46*A52*A63+A16*A22*A41*A54*A63-A12*A26*A41*A54*A63-A16*A21*A42*
     $A54*A63+A11*A26*A42*A54*A63+A12*A21*A46*A54*A63-A11*A22*A46*A54*
     $A63-A14*A22*A41*A56*A63+A12*A24*A41*A56*A63+A14*A21*A42*A56*A63-
     $A11*A24*A42*A56*A63-A12*A21*A44*A56*A63+A11*A22*A44*A56*A63-A16*
     $A23*A42*A51*A64+A13*A26*A42*A51*A64+A16*A22*A43*A51*A64-A12*A26*
     $A43*A51*A64-A13*A22*A46*A51*A64+A12*A23*A46*A51*A64+A16*A23*A41*
     $A52*A64-A13*A26*A41*A52*A64-A16*A21*A43*A52*A64+A11*A26*A43*A52*
     $A64+A13*A21*A46*A52*A64-A11*A23*A46*A52*A64-A16*A22*A41*A53*A64+
     $A12*A26*A41*A53*A64+A16*A21*A42*A53*A64-A11*A26*A42*A53*A64-A12*
     $A21*A46*A53*A64+A11*A22*A46*A53*A64+A13*A22*A41*A56*A64-A12*A23*
     $A41*A56*A64-A13*A21*A42*A56*A64+A11*A23*A42*A56*A64+A12*A21*A43*
     $A56*A64-A11*A22*A43*A56*A64+A14*A23*A42*A51*A66-A13*A24*A42*A51*
     $A66-A14*A22*A43*A51*A66+A12*A24*A43*A51*A66+A13*A22*A44*A51*A66-
     $A12*A23*A44*A51*A66-A14*A23*A41*A52*A66+A13*A24*A41*A52*A66+A14*
     $A21*A43*A52*A66-A11*A24*A43*A52*A66-A13*A21*A44*A52*A66+A11*A23*
     $A44*A52*A66+A14*A22*A41*A53*A66-A12*A24*A41*A53*A66-A14*A21*A42*
     $A53*A66+A11*A24*A42*A53*A66+A12*A21*A44*A53*A66-A11*A22*A44*A53*
     $A66-A13*A22*A41*A54*A66+A12*A23*A41*A54*A66+A13*A21*A42*A54*A66-
     $A11*A23*A42*A54*A66-A12*A21*A43*A54*A66+A11*A22*A43*A54*A66

      COFACTOR(4,5) = -A16*A24*A33*A52*
     $A61+A14*A26*A33*A52*A61+A16*A23*A34*A52*A61-A13*A26*A34*A52*A61-
     $A14*A23*A36*A52*A61+A13*A24*A36*A52*A61+A16*A24*A32*A53*A61-A14*
     $A26*A32*A53*A61-A16*A22*A34*A53*A61+A12*A26*A34*A53*A61+A14*A22*
     $A36*A53*A61-A12*A24*A36*A53*A61-A16*A23*A32*A54*A61+A13*A26*A32*
     $A54*A61+A16*A22*A33*A54*A61-A12*A26*A33*A54*A61-A13*A22*A36*A54*
     $A61+A12*A23*A36*A54*A61+A14*A23*A32*A56*A61-A13*A24*A32*A56*A61-
     $A14*A22*A33*A56*A61+A12*A24*A33*A56*A61+A13*A22*A34*A56*A61-A12*
     $A23*A34*A56*A61+A16*A24*A33*A51*A62-A14*A26*A33*A51*A62-A16*A23*
     $A34*A51*A62+A13*A26*A34*A51*A62+A14*A23*A36*A51*A62-A13*A24*A36*
     $A51*A62-A16*A24*A31*A53*A62+A14*A26*A31*A53*A62+A16*A21*A34*A53*
     $A62-A11*A26*A34*A53*A62-A14*A21*A36*A53*A62+A11*A24*A36*A53*A62+
     $A16*A23*A31*A54*A62-A13*A26*A31*A54*A62-A16*A21*A33*A54*A62+A11*
     $A26*A33*A54*A62+A13*A21*A36*A54*A62-A11*A23*A36*A54*A62-A14*A23*
     $A31*A56*A62+A13*A24*A31*A56*A62+A14*A21*A33*A56*A62-A11*A24*A33*
     $A56*A62-A13*A21*A34*A56*A62+A11*A23*A34*A56*A62-A16*A24*A32*A51*
     $A63+A14*A26*A32*A51*A63+A16*A22*A34*A51*A63-A12*A26*A34*A51*A63-
     $A14*A22*A36*A51*A63+A12*A24*A36*A51*A63+A16*A24*A31*A52*A63-A14*
     $A26*A31*A52*A63-A16*A21*A34*A52*A63+A11*A26*A34*A52*A63+A14*A21*
     $A36*A52*A63-A11*A24*A36*A52*A63-A16*A22*A31*A54*A63+A12*A26*A31*
     $A54*A63+A16*A21*A32*A54*A63-A11*A26*A32*A54*A63-A12*A21*A36*A54*
     $A63+A11*A22*A36*A54*A63+A14*A22*A31*A56*A63-A12*A24*A31*A56*A63-
     $A14*A21*A32*A56*A63+A11*A24*A32*A56*A63+A12*A21*A34*A56*A63-A11*
     $A22*A34*A56*A63+A16*A23*A32*A51*A64-A13*A26*A32*A51*A64-A16*A22*
     $A33*A51*A64+A12*A26*A33*A51*A64+A13*A22*A36*A51*A64-A12*A23*A36*
     $A51*A64-A16*A23*A31*A52*A64+A13*A26*A31*A52*A64+A16*A21*A33*A52*
     $A64-A11*A26*A33*A52*A64-A13*A21*A36*A52*A64+A11*A23*A36*A52*A64+
     $A16*A22*A31*A53*A64-A12*A26*A31*A53*A64-A16*A21*A32*A53*A64+A11*
     $A26*A32*A53*A64+A12*A21*A36*A53*A64-A11*A22*A36*A53*A64-A13*A22*
     $A31*A56*A64+A12*A23*A31*A56*A64+A13*A21*A32*A56*A64-A11*A23*A32*
     $A56*A64-A12*A21*A33*A56*A64+A11*A22*A33*A56*A64-A14*A23*A32*A51*
     $A66+A13*A24*A32*A51*A66+A14*A22*A33*A51*A66-A12*A24*A33*A51*A66-
     $A13*A22*A34*A51*A66+A12*A23*A34*A51*A66+A14*A23*A31*A52*A66-A13*
     $A24*A31*A52*A66-A14*A21*A33*A52*A66+A11*A24*A33*A52*A66+A13*A21*
     $A34*A52*A66-A11*A23*A34*A52*A66-A14*A22*A31*A53*A66+A12*A24*A31*
     $A53*A66+A14*A21*A32*A53*A66-A11*A24*A32*A53*A66-A12*A21*A34*A53*
     $A66+A11*A22*A34*A53*A66+A13*A22*A31*A54*A66-A12*A23*A31*A54*A66-
     $A13*A21*A32*A54*A66+A11*A23*A32*A54*A66+A12*A21*A33*A54*A66-A11*
     $A22*A33*A54*A66

      COFACTOR(5,5) = A16*A24*A33*A42*A61-A14*A26*A33*A42*A61-A16*A23*
     $A34*A42*A61+A13*A26*A34*A42*A61+A14*A23*A36*A42*A61-A13*A24*A36*
     $A42*A61-A16*A24*A32*A43*A61+A14*A26*A32*A43*A61+A16*A22*A34*A43*
     $A61-A12*A26*A34*A43*A61-A14*A22*A36*A43*A61+A12*A24*A36*A43*A61+
     $A16*A23*A32*A44*A61-A13*A26*A32*A44*A61-A16*A22*A33*A44*A61+A12*
     $A26*A33*A44*A61+A13*A22*A36*A44*A61-A12*A23*A36*A44*A61-A14*A23*
     $A32*A46*A61+A13*A24*A32*A46*A61+A14*A22*A33*A46*A61-A12*A24*A33*
     $A46*A61-A13*A22*A34*A46*A61+A12*A23*A34*A46*A61-A16*A24*A33*A41*
     $A62+A14*A26*A33*A41*A62+A16*A23*A34*A41*A62-A13*A26*A34*A41*A62-
     $A14*A23*A36*A41*A62+A13*A24*A36*A41*A62+A16*A24*A31*A43*A62-A14*
     $A26*A31*A43*A62-A16*A21*A34*A43*A62+A11*A26*A34*A43*A62+A14*A21*
     $A36*A43*A62-A11*A24*A36*A43*A62-A16*A23*A31*A44*A62+A13*A26*A31*
     $A44*A62+A16*A21*A33*A44*A62-A11*A26*A33*A44*A62-A13*A21*A36*A44*
     $A62+A11*A23*A36*A44*A62+A14*A23*A31*A46*A62-A13*A24*A31*A46*A62-
     $A14*A21*A33*A46*A62+A11*A24*A33*A46*A62+A13*A21*A34*A46*A62-A11*
     $A23*A34*A46*A62+A16*A24*A32*A41*A63-A14*A26*A32*A41*A63-A16*A22*
     $A34*A41*A63+A12*A26*A34*A41*A63+A14*A22*A36*A41*A63-A12*A24*A36*
     $A41*A63-A16*A24*A31*A42*A63+A14*A26*A31*A42*A63+A16*A21*A34*A42*
     $A63-A11*A26*A34*A42*A63-A14*A21*A36*A42*A63+A11*A24*A36*A42*A63+
     $A16*A22*A31*A44*A63-A12*A26*A31*A44*A63-A16*A21*A32*A44*A63+A11*
     $A26*A32*A44*A63+A12*A21*A36*A44*A63-A11*A22*A36*A44*A63-A14*A22*
     $A31*A46*A63+A12*A24*A31*A46*A63+A14*A21*A32*A46*A63-A11*A24*A32*
     $A46*A63-A12*A21*A34*A46*A63+A11*A22*A34*A46*A63-A16*A23*A32*A41*
     $A64+A13*A26*A32*A41*A64+A16*A22*A33*A41*A64-A12*A26*A33*A41*A64-
     $A13*A22*A36*A41*A64+A12*A23*A36*A41*A64+A16*A23*A31*A42*A64-A13*
     $A26*A31*A42*A64-A16*A21*A33*A42*A64+A11*A26*A33*A42*A64+A13*A21*
     $A36*A42*A64-A11*A23*A36*A42*A64-A16*A22*A31*A43*A64+A12*A26*A31*
     $A43*A64+A16*A21*A32*A43*A64-A11*A26*A32*A43*A64-A12*A21*A36*A43*
     $A64+A11*A22*A36*A43*A64+A13*A22*A31*A46*A64-A12*A23*A31*A46*A64-
     $A13*A21*A32*A46*A64+A11*A23*A32*A46*A64+A12*A21*A33*A46*A64-A11*
     $A22*A33*A46*A64+A14*A23*A32*A41*A66-A13*A24*A32*A41*A66-A14*A22*
     $A33*A41*A66+A12*A24*A33*A41*A66+A13*A22*A34*A41*A66-A12*A23*A34*
     $A41*A66-A14*A23*A31*A42*A66+A13*A24*A31*A42*A66+A14*A21*A33*A42*
     $A66-A11*A24*A33*A42*A66-A13*A21*A34*A42*A66+A11*A23*A34*A42*A66+
     $A14*A22*A31*A43*A66-A12*A24*A31*A43*A66-A14*A21*A32*A43*A66+A11*
     $A24*A32*A43*A66+A12*A21*A34*A43*A66-A11*A22*A34*A43*A66-A13*A22*
     $A31*A44*A66+A12*A23*A31*A44*A66+A13*A21*A32*A44*A66-A11*A23*A32*
     $A44*A66-A12*A21*A33*A44*A66+A11*A22*A33*A44*A66

      COFACTOR(6,5) = -A16*A24*A33*A42*A51+A14*A26*A33*A42*A51+A16*A23*
     $A34*A42*A51-A13*A26*A34*A42*A51-A14*A23*A36*A42*A51+A13*A24*A36*
     $A42*A51+A16*A24*A32*A43*A51-A14*A26*A32*A43*A51-A16*A22*A34*A43*
     $A51+A12*A26*A34*A43*A51+A14*A22*A36*A43*A51-A12*A24*A36*A43*A51-
     $A16*A23*A32*A44*A51+A13*A26*A32*A44*A51+A16*A22*A33*A44*A51-A12*
     $A26*A33*A44*A51-A13*A22*A36*A44*A51+A12*A23*A36*A44*A51+A14*A23*
     $A32*A46*A51-A13*A24*A32*A46*A51-A14*A22*A33*A46*A51+A12*A24*A33*
     $A46*A51+A13*A22*A34*A46*A51-A12*A23*A34*A46*A51+A16*A24*A33*A41*
     $A52-A14*A26*A33*A41*A52-A16*A23*A34*A41*A52+A13*A26*A34*A41*A52+
     $A14*A23*A36*A41*A52-A13*A24*A36*A41*A52-A16*A24*A31*A43*A52+A14*
     $A26*A31*A43*A52+A16*A21*A34*A43*A52-A11*A26*A34*A43*A52-A14*A21*
     $A36*A43*A52+A11*A24*A36*A43*A52+A16*A23*A31*A44*A52-A13*A26*A31*
     $A44*A52-A16*A21*A33*A44*A52+A11*A26*A33*A44*A52+A13*A21*A36*A44*
     $A52-A11*A23*A36*A44*A52-A14*A23*A31*A46*A52+A13*A24*A31*A46*A52+
     $A14*A21*A33*A46*A52-A11*A24*A33*A46*A52-A13*A21*A34*A46*A52+A11*
     $A23*A34*A46*A52-A16*A24*A32*A41*A53+A14*A26*A32*A41*A53+A16*A22*
     $A34*A41*A53-A12*A26*A34*A41*A53-A14*A22*A36*A41*A53+A12*A24*A36*
     $A41*A53+A16*A24*A31*A42*A53-A14*A26*A31*A42*A53-A16*A21*A34*A42*
     $A53+A11*A26*A34*A42*A53+A14*A21*A36*A42*A53-A11*A24*A36*A42*A53-
     $A16*A22*A31*A44*A53+A12*A26*A31*A44*A53+A16*A21*A32*A44*A53-A11*
     $A26*A32*A44*A53-A12*A21*A36*A44*A53+A11*A22*A36*A44*A53+A14*A22*
     $A31*A46*A53-A12*A24*A31*A46*A53-A14*A21*A32*A46*A53+A11*A24*A32*
     $A46*A53+A12*A21*A34*A46*A53-A11*A22*A34*A46*A53+A16*A23*A32*A41*
     $A54-A13*A26*A32*A41*A54-A16*A22*A33*A41*A54+A12*A26*A33*A41*A54+
     $A13*A22*A36*A41*A54-A12*A23*A36*A41*A54-A16*A23*A31*A42*A54+A13*
     $A26*A31*A42*A54+A16*A21*A33*A42*A54-A11*A26*A33*A42*A54-A13*A21*
     $A36*A42*A54+A11*A23*A36*A42*A54+A16*A22*A31*A43*A54-A12*A26*A31*
     $A43*A54-A16*A21*A32*A43*A54+A11*A26*A32*A43*A54+A12*A21*A36*A43*
     $A54-A11*A22*A36*A43*A54-A13*A22*A31*A46*A54+A12*A23*A31*A46*A54+
     $A13*A21*A32*A46*A54-A11*A23*A32*A46*A54-A12*A21*A33*A46*A54+A11*
     $A22*A33*A46*A54-A14*A23*A32*A41*A56+A13*A24*A32*A41*A56+A14*A22*
     $A33*A41*A56-A12*A24*A33*A41*A56-A13*A22*A34*A41*A56+A12*A23*A34*
     $A41*A56+A14*A23*A31*A42*A56-A13*A24*A31*A42*A56-A14*A21*A33*A42*
     $A56+A11*A24*A33*A42*A56+A13*A21*A34*A42*A56-A11*A23*A34*A42*A56-
     $A14*A22*A31*A43*A56+A12*A24*A31*A43*A56+A14*A21*A32*A43*A56-A11*
     $A24*A32*A43*A56-A12*A21*A34*A43*A56+A11*A22*A34*A43*A56+A13*A22*
     $A31*A44*A56-A12*A23*A31*A44*A56-A13*A21*A32*A44*A56+A11*A23*A32*
     $A44*A56+A12*A21*A33*A44*A56-A11*A22*A33*A44*A56

      COFACTOR(1,6) = -A25*A34*A43*A52*A61+A24*
     $A35*A43*A52*A61+A25*A33*A44*A52*A61-A23*A35*A44*A52*A61-A24*A33*
     $A45*A52*A61+A23*A34*A45*A52*A61+A25*A34*A42*A53*A61-A24*A35*A42*
     $A53*A61-A25*A32*A44*A53*A61+A22*A35*A44*A53*A61+A24*A32*A45*A53*
     $A61-A22*A34*A45*A53*A61-A25*A33*A42*A54*A61+A23*A35*A42*A54*A61+
     $A25*A32*A43*A54*A61-A22*A35*A43*A54*A61-A23*A32*A45*A54*A61+A22*
     $A33*A45*A54*A61+A24*A33*A42*A55*A61-A23*A34*A42*A55*A61-A24*A32*
     $A43*A55*A61+A22*A34*A43*A55*A61+A23*A32*A44*A55*A61-A22*A33*A44*
     $A55*A61+A25*A34*A43*A51*A62-A24*A35*A43*A51*A62-A25*A33*A44*A51*
     $A62+A23*A35*A44*A51*A62+A24*A33*A45*A51*A62-A23*A34*A45*A51*A62-
     $A25*A34*A41*A53*A62+A24*A35*A41*A53*A62+A25*A31*A44*A53*A62-A21*
     $A35*A44*A53*A62-A24*A31*A45*A53*A62+A21*A34*A45*A53*A62+A25*A33*
     $A41*A54*A62-A23*A35*A41*A54*A62-A25*A31*A43*A54*A62+A21*A35*A43*
     $A54*A62+A23*A31*A45*A54*A62-A21*A33*A45*A54*A62-A24*A33*A41*A55*
     $A62+A23*A34*A41*A55*A62+A24*A31*A43*A55*A62-A21*A34*A43*A55*A62-
     $A23*A31*A44*A55*A62+A21*A33*A44*A55*A62-A25*A34*A42*A51*A63+A24*
     $A35*A42*A51*A63+A25*A32*A44*A51*A63-A22*A35*A44*A51*A63-A24*A32*
     $A45*A51*A63+A22*A34*A45*A51*A63+A25*A34*A41*A52*A63-A24*A35*A41*
     $A52*A63-A25*A31*A44*A52*A63+A21*A35*A44*A52*A63+A24*A31*A45*A52*
     $A63-A21*A34*A45*A52*A63-A25*A32*A41*A54*A63+A22*A35*A41*A54*A63+
     $A25*A31*A42*A54*A63-A21*A35*A42*A54*A63-A22*A31*A45*A54*A63+A21*
     $A32*A45*A54*A63+A24*A32*A41*A55*A63-A22*A34*A41*A55*A63-A24*A31*
     $A42*A55*A63+A21*A34*A42*A55*A63+A22*A31*A44*A55*A63-A21*A32*A44*
     $A55*A63+A25*A33*A42*A51*A64-A23*A35*A42*A51*A64-A25*A32*A43*A51*
     $A64+A22*A35*A43*A51*A64+A23*A32*A45*A51*A64-A22*A33*A45*A51*A64-
     $A25*A33*A41*A52*A64+A23*A35*A41*A52*A64+A25*A31*A43*A52*A64-A21*
     $A35*A43*A52*A64-A23*A31*A45*A52*A64+A21*A33*A45*A52*A64+A25*A32*
     $A41*A53*A64-A22*A35*A41*A53*A64-A25*A31*A42*A53*A64+A21*A35*A42*
     $A53*A64+A22*A31*A45*A53*A64-A21*A32*A45*A53*A64-A23*A32*A41*A55*
     $A64+A22*A33*A41*A55*A64+A23*A31*A42*A55*A64-A21*A33*A42*A55*A64-
     $A22*A31*A43*A55*A64+A21*A32*A43*A55*A64-A24*A33*A42*A51*A65+A23*
     $A34*A42*A51*A65+A24*A32*A43*A51*A65-A22*A34*A43*A51*A65-A23*A32*
     $A44*A51*A65+A22*A33*A44*A51*A65+A24*A33*A41*A52*A65-A23*A34*A41*
     $A52*A65-A24*A31*A43*A52*A65+A21*A34*A43*A52*A65+A23*A31*A44*A52*
     $A65-A21*A33*A44*A52*A65-A24*A32*A41*A53*A65+A22*A34*A41*A53*A65+
     $A24*A31*A42*A53*A65-A21*A34*A42*A53*A65-A22*A31*A44*A53*A65+A21*
     $A32*A44*A53*A65+A23*A32*A41*A54*A65-A22*A33*A41*A54*A65-A23*A31*
     $A42*A54*A65+A21*A33*A42*A54*A65+A22*A31*A43*A54*A65-A21*A32*A43*
     $A54*A65

      COFACTOR(2,6) = A15*A34*
     $A43*A52*A61-A14*A35*A43*A52*A61-A15*A33*A44*A52*A61+A13*A35*A44*
     $A52*A61+A14*A33*A45*A52*A61-A13*A34*A45*A52*A61-A15*A34*A42*A53*
     $A61+A14*A35*A42*A53*A61+A15*A32*A44*A53*A61-A12*A35*A44*A53*A61-
     $A14*A32*A45*A53*A61+A12*A34*A45*A53*A61+A15*A33*A42*A54*A61-A13*
     $A35*A42*A54*A61-A15*A32*A43*A54*A61+A12*A35*A43*A54*A61+A13*A32*
     $A45*A54*A61-A12*A33*A45*A54*A61-A14*A33*A42*A55*A61+A13*A34*A42*
     $A55*A61+A14*A32*A43*A55*A61-A12*A34*A43*A55*A61-A13*A32*A44*A55*
     $A61+A12*A33*A44*A55*A61-A15*A34*A43*A51*A62+A14*A35*A43*A51*A62+
     $A15*A33*A44*A51*A62-A13*A35*A44*A51*A62-A14*A33*A45*A51*A62+A13*
     $A34*A45*A51*A62+A15*A34*A41*A53*A62-A14*A35*A41*A53*A62-A15*A31*
     $A44*A53*A62+A11*A35*A44*A53*A62+A14*A31*A45*A53*A62-A11*A34*A45*
     $A53*A62-A15*A33*A41*A54*A62+A13*A35*A41*A54*A62+A15*A31*A43*A54*
     $A62-A11*A35*A43*A54*A62-A13*A31*A45*A54*A62+A11*A33*A45*A54*A62+
     $A14*A33*A41*A55*A62-A13*A34*A41*A55*A62-A14*A31*A43*A55*A62+A11*
     $A34*A43*A55*A62+A13*A31*A44*A55*A62-A11*A33*A44*A55*A62+A15*A34*
     $A42*A51*A63-A14*A35*A42*A51*A63-A15*A32*A44*A51*A63+A12*A35*A44*
     $A51*A63+A14*A32*A45*A51*A63-A12*A34*A45*A51*A63-A15*A34*A41*A52*
     $A63+A14*A35*A41*A52*A63+A15*A31*A44*A52*A63-A11*A35*A44*A52*A63-
     $A14*A31*A45*A52*A63+A11*A34*A45*A52*A63+A15*A32*A41*A54*A63-A12*
     $A35*A41*A54*A63-A15*A31*A42*A54*A63+A11*A35*A42*A54*A63+A12*A31*
     $A45*A54*A63-A11*A32*A45*A54*A63-A14*A32*A41*A55*A63+A12*A34*A41*
     $A55*A63+A14*A31*A42*A55*A63-A11*A34*A42*A55*A63-A12*A31*A44*A55*
     $A63+A11*A32*A44*A55*A63-A15*A33*A42*A51*A64+A13*A35*A42*A51*A64+
     $A15*A32*A43*A51*A64-A12*A35*A43*A51*A64-A13*A32*A45*A51*A64+A12*
     $A33*A45*A51*A64+A15*A33*A41*A52*A64-A13*A35*A41*A52*A64-A15*A31*
     $A43*A52*A64+A11*A35*A43*A52*A64+A13*A31*A45*A52*A64-A11*A33*A45*
     $A52*A64-A15*A32*A41*A53*A64+A12*A35*A41*A53*A64+A15*A31*A42*A53*
     $A64-A11*A35*A42*A53*A64-A12*A31*A45*A53*A64+A11*A32*A45*A53*A64+
     $A13*A32*A41*A55*A64-A12*A33*A41*A55*A64-A13*A31*A42*A55*A64+A11*
     $A33*A42*A55*A64+A12*A31*A43*A55*A64-A11*A32*A43*A55*A64+A14*A33*
     $A42*A51*A65-A13*A34*A42*A51*A65-A14*A32*A43*A51*A65+A12*A34*A43*
     $A51*A65+A13*A32*A44*A51*A65-A12*A33*A44*A51*A65-A14*A33*A41*A52*
     $A65+A13*A34*A41*A52*A65+A14*A31*A43*A52*A65-A11*A34*A43*A52*A65-
     $A13*A31*A44*A52*A65+A11*A33*A44*A52*A65+A14*A32*A41*A53*A65-A12*
     $A34*A41*A53*A65-A14*A31*A42*A53*A65+A11*A34*A42*A53*A65+A12*A31*
     $A44*A53*A65-A11*A32*A44*A53*A65-A13*A32*A41*A54*A65+A12*A33*A41*
     $A54*A65+A13*A31*A42*A54*A65-A11*A33*A42*A54*A65-A12*A31*A43*A54*
     $A65+A11*A32*A43*A54*A65

      COFACTOR(3,6) = -A15*A24*A43*A52*A61+A14*A25*A43*A52*A61+A15*A23*
     $A44*A52*A61-A13*A25*A44*A52*A61-A14*A23*A45*A52*A61+A13*A24*A45*
     $A52*A61+A15*A24*A42*A53*A61-A14*A25*A42*A53*A61-A15*A22*A44*A53*
     $A61+A12*A25*A44*A53*A61+A14*A22*A45*A53*A61-A12*A24*A45*A53*A61-
     $A15*A23*A42*A54*A61+A13*A25*A42*A54*A61+A15*A22*A43*A54*A61-A12*
     $A25*A43*A54*A61-A13*A22*A45*A54*A61+A12*A23*A45*A54*A61+A14*A23*
     $A42*A55*A61-A13*A24*A42*A55*A61-A14*A22*A43*A55*A61+A12*A24*A43*
     $A55*A61+A13*A22*A44*A55*A61-A12*A23*A44*A55*A61+A15*A24*A43*A51*
     $A62-A14*A25*A43*A51*A62-A15*A23*A44*A51*A62+A13*A25*A44*A51*A62+
     $A14*A23*A45*A51*A62-A13*A24*A45*A51*A62-A15*A24*A41*A53*A62+A14*
     $A25*A41*A53*A62+A15*A21*A44*A53*A62-A11*A25*A44*A53*A62-A14*A21*
     $A45*A53*A62+A11*A24*A45*A53*A62+A15*A23*A41*A54*A62-A13*A25*A41*
     $A54*A62-A15*A21*A43*A54*A62+A11*A25*A43*A54*A62+A13*A21*A45*A54*
     $A62-A11*A23*A45*A54*A62-A14*A23*A41*A55*A62+A13*A24*A41*A55*A62+
     $A14*A21*A43*A55*A62-A11*A24*A43*A55*A62-A13*A21*A44*A55*A62+A11*
     $A23*A44*A55*A62-A15*A24*A42*A51*A63+A14*A25*A42*A51*A63+A15*A22*
     $A44*A51*A63-A12*A25*A44*A51*A63-A14*A22*A45*A51*A63+A12*A24*A45*
     $A51*A63+A15*A24*A41*A52*A63-A14*A25*A41*A52*A63-A15*A21*A44*A52*
     $A63+A11*A25*A44*A52*A63+A14*A21*A45*A52*A63-A11*A24*A45*A52*A63-
     $A15*A22*A41*A54*A63+A12*A25*A41*A54*A63+A15*A21*A42*A54*A63-A11*
     $A25*A42*A54*A63-A12*A21*A45*A54*A63+A11*A22*A45*A54*A63+A14*A22*
     $A41*A55*A63-A12*A24*A41*A55*A63-A14*A21*A42*A55*A63+A11*A24*A42*
     $A55*A63+A12*A21*A44*A55*A63-A11*A22*A44*A55*A63+A15*A23*A42*A51*
     $A64-A13*A25*A42*A51*A64-A15*A22*A43*A51*A64+A12*A25*A43*A51*A64+
     $A13*A22*A45*A51*A64-A12*A23*A45*A51*A64-A15*A23*A41*A52*A64+A13*
     $A25*A41*A52*A64+A15*A21*A43*A52*A64-A11*A25*A43*A52*A64-A13*A21*
     $A45*A52*A64+A11*A23*A45*A52*A64+A15*A22*A41*A53*A64-A12*A25*A41*
     $A53*A64-A15*A21*A42*A53*A64+A11*A25*A42*A53*A64+A12*A21*A45*A53*
     $A64-A11*A22*A45*A53*A64-A13*A22*A41*A55*A64+A12*A23*A41*A55*A64+
     $A13*A21*A42*A55*A64-A11*A23*A42*A55*A64-A12*A21*A43*A55*A64+A11*
     $A22*A43*A55*A64-A14*A23*A42*A51*A65+A13*A24*A42*A51*A65+A14*A22*
     $A43*A51*A65-A12*A24*A43*A51*A65-A13*A22*A44*A51*A65+A12*A23*A44*
     $A51*A65+A14*A23*A41*A52*A65-A13*A24*A41*A52*A65-A14*A21*A43*A52*
     $A65+A11*A24*A43*A52*A65+A13*A21*A44*A52*A65-A11*A23*A44*A52*A65-
     $A14*A22*A41*A53*A65+A12*A24*A41*A53*A65+A14*A21*A42*A53*A65-A11*
     $A24*A42*A53*A65-A12*A21*A44*A53*A65+A11*A22*A44*A53*A65+A13*A22*
     $A41*A54*A65-A12*A23*A41*A54*A65-A13*A21*A42*A54*A65+A11*A23*A42*
     $A54*A65+A12*A21*A43*A54*A65-A11*A22*A43*A54*A65

      COFACTOR(4,6) = A15*A24*A33*A52*A61-A14*A25*A33*A52*
     $A61-A15*A23*A34*A52*A61+A13*A25*A34*A52*A61+A14*A23*A35*A52*A61-
     $A13*A24*A35*A52*A61-A15*A24*A32*A53*A61+A14*A25*A32*A53*A61+A15*
     $A22*A34*A53*A61-A12*A25*A34*A53*A61-A14*A22*A35*A53*A61+A12*A24*
     $A35*A53*A61+A15*A23*A32*A54*A61-A13*A25*A32*A54*A61-A15*A22*A33*
     $A54*A61+A12*A25*A33*A54*A61+A13*A22*A35*A54*A61-A12*A23*A35*A54*
     $A61-A14*A23*A32*A55*A61+A13*A24*A32*A55*A61+A14*A22*A33*A55*A61-
     $A12*A24*A33*A55*A61-A13*A22*A34*A55*A61+A12*A23*A34*A55*A61-A15*
     $A24*A33*A51*A62+A14*A25*A33*A51*A62+A15*A23*A34*A51*A62-A13*A25*
     $A34*A51*A62-A14*A23*A35*A51*A62+A13*A24*A35*A51*A62+A15*A24*A31*
     $A53*A62-A14*A25*A31*A53*A62-A15*A21*A34*A53*A62+A11*A25*A34*A53*
     $A62+A14*A21*A35*A53*A62-A11*A24*A35*A53*A62-A15*A23*A31*A54*A62+
     $A13*A25*A31*A54*A62+A15*A21*A33*A54*A62-A11*A25*A33*A54*A62-A13*
     $A21*A35*A54*A62+A11*A23*A35*A54*A62+A14*A23*A31*A55*A62-A13*A24*
     $A31*A55*A62-A14*A21*A33*A55*A62+A11*A24*A33*A55*A62+A13*A21*A34*
     $A55*A62-A11*A23*A34*A55*A62+A15*A24*A32*A51*A63-A14*A25*A32*A51*
     $A63-A15*A22*A34*A51*A63+A12*A25*A34*A51*A63+A14*A22*A35*A51*A63-
     $A12*A24*A35*A51*A63-A15*A24*A31*A52*A63+A14*A25*A31*A52*A63+A15*
     $A21*A34*A52*A63-A11*A25*A34*A52*A63-A14*A21*A35*A52*A63+A11*A24*
     $A35*A52*A63+A15*A22*A31*A54*A63-A12*A25*A31*A54*A63-A15*A21*A32*
     $A54*A63+A11*A25*A32*A54*A63+A12*A21*A35*A54*A63-A11*A22*A35*A54*
     $A63-A14*A22*A31*A55*A63+A12*A24*A31*A55*A63+A14*A21*A32*A55*A63-
     $A11*A24*A32*A55*A63-A12*A21*A34*A55*A63+A11*A22*A34*A55*A63-A15*
     $A23*A32*A51*A64+A13*A25*A32*A51*A64+A15*A22*A33*A51*A64-A12*A25*
     $A33*A51*A64-A13*A22*A35*A51*A64+A12*A23*A35*A51*A64+A15*A23*A31*
     $A52*A64-A13*A25*A31*A52*A64-A15*A21*A33*A52*A64+A11*A25*A33*A52*
     $A64+A13*A21*A35*A52*A64-A11*A23*A35*A52*A64-A15*A22*A31*A53*A64+
     $A12*A25*A31*A53*A64+A15*A21*A32*A53*A64-A11*A25*A32*A53*A64-A12*
     $A21*A35*A53*A64+A11*A22*A35*A53*A64+A13*A22*A31*A55*A64-A12*A23*
     $A31*A55*A64-A13*A21*A32*A55*A64+A11*A23*A32*A55*A64+A12*A21*A33*
     $A55*A64-A11*A22*A33*A55*A64+A14*A23*A32*A51*A65-A13*A24*A32*A51*
     $A65-A14*A22*A33*A51*A65+A12*A24*A33*A51*A65+A13*A22*A34*A51*A65-
     $A12*A23*A34*A51*A65-A14*A23*A31*A52*A65+A13*A24*A31*A52*A65+A14*
     $A21*A33*A52*A65-A11*A24*A33*A52*A65-A13*A21*A34*A52*A65+A11*A23*
     $A34*A52*A65+A14*A22*A31*A53*A65-A12*A24*A31*A53*A65-A14*A21*A32*
     $A53*A65+A11*A24*A32*A53*A65+A12*A21*A34*A53*A65-A11*A22*A34*A53*
     $A65-A13*A22*A31*A54*A65+A12*A23*A31*A54*A65+A13*A21*A32*A54*A65-
     $A11*A23*A32*A54*A65-A12*A21*A33*A54*A65+A11*A22*A33*A54*A65

      COFACTOR(5,6) = -A15*A24*A33*A42*
     $A61+A14*A25*A33*A42*A61+A15*A23*A34*A42*A61-A13*A25*A34*A42*A61-
     $A14*A23*A35*A42*A61+A13*A24*A35*A42*A61+A15*A24*A32*A43*A61-A14*
     $A25*A32*A43*A61-A15*A22*A34*A43*A61+A12*A25*A34*A43*A61+A14*A22*
     $A35*A43*A61-A12*A24*A35*A43*A61-A15*A23*A32*A44*A61+A13*A25*A32*
     $A44*A61+A15*A22*A33*A44*A61-A12*A25*A33*A44*A61-A13*A22*A35*A44*
     $A61+A12*A23*A35*A44*A61+A14*A23*A32*A45*A61-A13*A24*A32*A45*A61-
     $A14*A22*A33*A45*A61+A12*A24*A33*A45*A61+A13*A22*A34*A45*A61-A12*
     $A23*A34*A45*A61+A15*A24*A33*A41*A62-A14*A25*A33*A41*A62-A15*A23*
     $A34*A41*A62+A13*A25*A34*A41*A62+A14*A23*A35*A41*A62-A13*A24*A35*
     $A41*A62-A15*A24*A31*A43*A62+A14*A25*A31*A43*A62+A15*A21*A34*A43*
     $A62-A11*A25*A34*A43*A62-A14*A21*A35*A43*A62+A11*A24*A35*A43*A62+
     $A15*A23*A31*A44*A62-A13*A25*A31*A44*A62-A15*A21*A33*A44*A62+A11*
     $A25*A33*A44*A62+A13*A21*A35*A44*A62-A11*A23*A35*A44*A62-A14*A23*
     $A31*A45*A62+A13*A24*A31*A45*A62+A14*A21*A33*A45*A62-A11*A24*A33*
     $A45*A62-A13*A21*A34*A45*A62+A11*A23*A34*A45*A62-A15*A24*A32*A41*
     $A63+A14*A25*A32*A41*A63+A15*A22*A34*A41*A63-A12*A25*A34*A41*A63-
     $A14*A22*A35*A41*A63+A12*A24*A35*A41*A63+A15*A24*A31*A42*A63-A14*
     $A25*A31*A42*A63-A15*A21*A34*A42*A63+A11*A25*A34*A42*A63+A14*A21*
     $A35*A42*A63-A11*A24*A35*A42*A63-A15*A22*A31*A44*A63+A12*A25*A31*
     $A44*A63+A15*A21*A32*A44*A63-A11*A25*A32*A44*A63-A12*A21*A35*A44*
     $A63+A11*A22*A35*A44*A63+A14*A22*A31*A45*A63-A12*A24*A31*A45*A63-
     $A14*A21*A32*A45*A63+A11*A24*A32*A45*A63+A12*A21*A34*A45*A63-A11*
     $A22*A34*A45*A63+A15*A23*A32*A41*A64-A13*A25*A32*A41*A64-A15*A22*
     $A33*A41*A64+A12*A25*A33*A41*A64+A13*A22*A35*A41*A64-A12*A23*A35*
     $A41*A64-A15*A23*A31*A42*A64+A13*A25*A31*A42*A64+A15*A21*A33*A42*
     $A64-A11*A25*A33*A42*A64-A13*A21*A35*A42*A64+A11*A23*A35*A42*A64+
     $A15*A22*A31*A43*A64-A12*A25*A31*A43*A64-A15*A21*A32*A43*A64+A11*
     $A25*A32*A43*A64+A12*A21*A35*A43*A64-A11*A22*A35*A43*A64-A13*A22*
     $A31*A45*A64+A12*A23*A31*A45*A64+A13*A21*A32*A45*A64-A11*A23*A32*
     $A45*A64-A12*A21*A33*A45*A64+A11*A22*A33*A45*A64-A14*A23*A32*A41*
     $A65+A13*A24*A32*A41*A65+A14*A22*A33*A41*A65-A12*A24*A33*A41*A65-
     $A13*A22*A34*A41*A65+A12*A23*A34*A41*A65+A14*A23*A31*A42*A65-A13*
     $A24*A31*A42*A65-A14*A21*A33*A42*A65+A11*A24*A33*A42*A65+A13*A21*
     $A34*A42*A65-A11*A23*A34*A42*A65-A14*A22*A31*A43*A65+A12*A24*A31*
     $A43*A65+A14*A21*A32*A43*A65-A11*A24*A32*A43*A65-A12*A21*A34*A43*
     $A65+A11*A22*A34*A43*A65+A13*A22*A31*A44*A65-A12*A23*A31*A44*A65-
     $A13*A21*A32*A44*A65+A11*A23*A32*A44*A65+A12*A21*A33*A44*A65-A11*
     $A22*A33*A44*A65

      COFACTOR(6,6) = A15*A24*A33*A42*A51-A14*A25*A33*A42*A51-A15*A23*
     $A34*A42*A51+A13*A25*A34*A42*A51+A14*A23*A35*A42*A51-A13*A24*A35*
     $A42*A51-A15*A24*A32*A43*A51+A14*A25*A32*A43*A51+A15*A22*A34*A43*
     $A51-A12*A25*A34*A43*A51-A14*A22*A35*A43*A51+A12*A24*A35*A43*A51+
     $A15*A23*A32*A44*A51-A13*A25*A32*A44*A51-A15*A22*A33*A44*A51+A12*
     $A25*A33*A44*A51+A13*A22*A35*A44*A51-A12*A23*A35*A44*A51-A14*A23*
     $A32*A45*A51+A13*A24*A32*A45*A51+A14*A22*A33*A45*A51-A12*A24*A33*
     $A45*A51-A13*A22*A34*A45*A51+A12*A23*A34*A45*A51-A15*A24*A33*A41*
     $A52+A14*A25*A33*A41*A52+A15*A23*A34*A41*A52-A13*A25*A34*A41*A52-
     $A14*A23*A35*A41*A52+A13*A24*A35*A41*A52+A15*A24*A31*A43*A52-A14*
     $A25*A31*A43*A52-A15*A21*A34*A43*A52+A11*A25*A34*A43*A52+A14*A21*
     $A35*A43*A52-A11*A24*A35*A43*A52-A15*A23*A31*A44*A52+A13*A25*A31*
     $A44*A52+A15*A21*A33*A44*A52-A11*A25*A33*A44*A52-A13*A21*A35*A44*
     $A52+A11*A23*A35*A44*A52+A14*A23*A31*A45*A52-A13*A24*A31*A45*A52-
     $A14*A21*A33*A45*A52+A11*A24*A33*A45*A52+A13*A21*A34*A45*A52-A11*
     $A23*A34*A45*A52+A15*A24*A32*A41*A53-A14*A25*A32*A41*A53-A15*A22*
     $A34*A41*A53+A12*A25*A34*A41*A53+A14*A22*A35*A41*A53-A12*A24*A35*
     $A41*A53-A15*A24*A31*A42*A53+A14*A25*A31*A42*A53+A15*A21*A34*A42*
     $A53-A11*A25*A34*A42*A53-A14*A21*A35*A42*A53+A11*A24*A35*A42*A53+
     $A15*A22*A31*A44*A53-A12*A25*A31*A44*A53-A15*A21*A32*A44*A53+A11*
     $A25*A32*A44*A53+A12*A21*A35*A44*A53-A11*A22*A35*A44*A53-A14*A22*
     $A31*A45*A53+A12*A24*A31*A45*A53+A14*A21*A32*A45*A53-A11*A24*A32*
     $A45*A53-A12*A21*A34*A45*A53+A11*A22*A34*A45*A53-A15*A23*A32*A41*
     $A54+A13*A25*A32*A41*A54+A15*A22*A33*A41*A54-A12*A25*A33*A41*A54-
     $A13*A22*A35*A41*A54+A12*A23*A35*A41*A54+A15*A23*A31*A42*A54-A13*
     $A25*A31*A42*A54-A15*A21*A33*A42*A54+A11*A25*A33*A42*A54+A13*A21*
     $A35*A42*A54-A11*A23*A35*A42*A54-A15*A22*A31*A43*A54+A12*A25*A31*
     $A43*A54+A15*A21*A32*A43*A54-A11*A25*A32*A43*A54-A12*A21*A35*A43*
     $A54+A11*A22*A35*A43*A54+A13*A22*A31*A45*A54-A12*A23*A31*A45*A54-
     $A13*A21*A32*A45*A54+A11*A23*A32*A45*A54+A12*A21*A33*A45*A54-A11*
     $A22*A33*A45*A54+A14*A23*A32*A41*A55-A13*A24*A32*A41*A55-A14*A22*
     $A33*A41*A55+A12*A24*A33*A41*A55+A13*A22*A34*A41*A55-A12*A23*A34*
     $A41*A55-A14*A23*A31*A42*A55+A13*A24*A31*A42*A55+A14*A21*A33*A42*
     $A55-A11*A24*A33*A42*A55-A13*A21*A34*A42*A55+A11*A23*A34*A42*A55+
     $A14*A22*A31*A43*A55-A12*A24*A31*A43*A55-A14*A21*A32*A43*A55+A11*
     $A24*A32*A43*A55+A12*A21*A34*A43*A55-A11*A22*A34*A43*A55-A13*A22*
     $A31*A44*A55+A12*A23*A31*A44*A55+A13*A21*A32*A44*A55-A11*A23*A32*
     $A44*A55-A12*A21*A33*A44*A55+A11*A22*A33*A44*A55

      AINV = TRANSPOSE(COFACTOR) / DET
      END
C-------------------------------------------------------------------------------------------
