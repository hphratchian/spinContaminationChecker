      module hphSpinFun_mod
!
!     This module supports the program scfEnergyTerms.
!
!     -H. P. Hratchian, 2020.
!
!
!     USE Connections
!
      use mqc_general
      use mqc_molecule
      use mqc_gaussian
      use mqc_algebra2
      use mqc_algebra
      use iso_fortran_env
!
!     Variable Declarations
!
      implicit none
      integer,parameter::IOut=6
!
!
!     Module Procedures
!
      CONTAINS
!
!
      subroutine formFock(nBasis,density,ERIs,coulomb)
!
!     This subroutine forms a Coulomb matrix from a density matrix and ERIs. The
!     input and output arrays are MQC variables.
!
!
      implicit none
      integer::nBasis
      type(MQC_Variable),intent(in)::density,ERIs
      type(MQC_Variable),intent(out)::coulomb
      integer::iMu,iNu,iLambda,iSigma
      real,dimension(nBasis,nBasis)::tempCoulomb
!
      call ERIs%print(IOut,' In formCoulomb: ERIs=',blankAtTop=.True.)
!
!     Do the work...
!
      tempCoulomb = float(0)
      do iSigma    = 1,nBasis
        do iLambda = 1,nBasis
          do iNu   = 1,nBasis
            do iMu = 1,nBasis
              write(IOut,'(/,1x,I3,I3,I3,I3,2x,F10.3,2x,F10.3,2x,I3,2x,I3)')  &
                iMu,iNu,iLambda,iSigma,  &
                float(ERIs%getVal([iMu,iNu,iLambda,iSigma])),  &
                float(density%getVal([iLambda,iSigma])),  &
                MQC_Variable_getArrayPosition(ERIs,[iMu,iNu,iLambda,iSigma])
              tempCoulomb(iMu,iNu) = tempCoulomb(iMu,iNu) +  &
                float(ERIs%getVal([iMu,iNu,iLambda,iSigma])) *  &
                float(density%getVal([iLambda,iSigma]))
            endDo
          endDo
        endDo
      endDo
      coulomb = tempCoulomb
!
      return
      end subroutine formFock
!
!
      subroutine formCoulomb(nBasis,density,ERIs,coulomb,initialize)
!
!     This subroutine forms a Coulomb matrix from a density matrix and ERIs. The
!     input and output arrays are MQC variables.
!
!
      implicit none
      integer::nBasis
      type(MQC_Variable),intent(in)::density,ERIs
      type(MQC_Variable),intent(out)::coulomb
      logical,optional::initialize
      integer::iMu,iNu,iLambda,iSigma
      real,dimension(nBasis,nBasis)::tempCoulomb
      logical::init
!
      call ERIs%print(IOut,' In formCoulomb: ERIs=',blankAtTop=.True.)
!
!     Do initial set-up work.
!
      init = .true.
      if(Present(initialize)) init = initialize
!
!     Work through the integral loops to build Coulomb contributions.
!
      if(init) tempCoulomb = float(0)
      do iSigma    = 1,nBasis
        do iLambda = 1,nBasis
          do iNu   = 1,nBasis
            do iMu = 1,nBasis
              write(IOut,'(/,1x,I3,I3,I3,I3,2x,F10.3,2x,F10.3,2x,I3,2x,I3)')  &
                iMu,iNu,iLambda,iSigma,  &
                float(ERIs%getVal([iMu,iNu,iLambda,iSigma])),  &
                float(density%getVal([iLambda,iSigma])),  &
                MQC_Variable_getArrayPosition(ERIs,[iMu,iNu,iLambda,iSigma])
              tempCoulomb(iMu,iNu) = tempCoulomb(iMu,iNu) +  &
                float(ERIs%getVal([iMu,iNu,iLambda,iSigma])) *  &
                float(density%getVal([iLambda,iSigma]))
            endDo
          endDo
        endDo
      endDo
      coulomb = tempCoulomb
!
      return
      end subroutine formCoulomb
!
!
      subroutine formExchange(nBasis,density,ERIs,exchange,initialize)
!
!     This subroutine forms an Exchange matrix from a density matrix and ERIs. The
!     input and output arrays are MQC variables.
!
!
      implicit none
      integer::nBasis
      type(MQC_Variable),intent(in)::density,ERIs
      type(MQC_Variable),intent(out)::exchange
      logical,optional::initialize
      integer::iMu,iNu,iLambda,iSigma
      real,dimension(nBasis,nBasis)::tempExchange
      logical::init
!
!     Do initial set-up work.
!
      init = .true.
      if(Present(initialize)) init = initialize
!
!     Work through the integral loops to build Exchange contributions.
!
      if(init) tempExchange = float(0)
      do iMu = 1,nBasis
        do iNu = 1,nBasis
          do iLambda = 1,nBasis
            do iSigma = 1,nBasis
              write(IOut,'(1x,I3,I3,I3,I3,2x,F10.3,2x,F10.3)')  &
                iMu,iNu,iLambda,iSigma,  &
                float(ERIs%getVal([iMu,iSigma,iLambda,iNu])),  &
                float(density%getVal([iLambda,iSigma]))
              tempExchange(iMu,iNu) = tempExchange(iMu,iNu) -  &
                float(ERIs%getVal([iMu,iSigma,iLambda,iNu])) *  &
                float(density%getVal([iLambda,iSigma]))
            endDo
          endDo
        endDo
      endDo
      exchange = tempExchange
!
      return
      end subroutine formExchange

!
!=====================================================================
!
!     PROCEDURE S2_MAT_ELEM    
! 
      Function S2_Mat_Elem(IOut,IPrint,NBasis,Alpha_String_1,Beta_String_1, &
      Alpha_String_2,Beta_String_2,MO_Overlap)
!
!     This function returns the CI S**2 matrix elements used for computing S**2
!     values of CI vectors for a given alpha and beta string combination.
!     The MO overlap matrix is required
!
!     Variable Declarations...
!
      Implicit None
      Integer(kind=int64)::IOut,IPrint,NBasis,IPos,JPos,IDiff,Det_Diff,NAlpha,NBeta, &
        IOcc,JOcc,KOcc,LOcc,Mat_Sign,Alpha_Diff_Cnt,Beta_Diff_Cnt,NBit_Ints, &
        I,J,II,JJ
      real(kind=real64)::S2_Mat_Elem,Zero=0.0d0,Quarter=0.25d0,ABTerm,One=1.0d0
      Integer(kind=int64),Dimension(4)::Orbs,Spin,Det
      Integer(kind=int64),Dimension(:)::Alpha_String_1,Alpha_String_2,Beta_String_1, &
        Beta_String_2
      Integer(kind=int64),Dimension(:),Allocatable::Alpha_Diff,Beta_Diff
      Real(kind=real64),Dimension(:,:),Allocatable::MO_Overlap
!
!      Write(IOut,*) 'alpha 1:'
!      Write(IOut,'(B64)') Alpha_String_1
!      Write(IOut,*) 'alpha 2:'
!      Write(IOut,'(B64)') Alpha_String_2
!      Write(IOut,*) 'alpha XOR results in slater condon'
!      Write(IOut,'(B64)') IEOR(Alpha_String_1,Alpha_String_2)   
!      Write(IOut,*) 'NDifa=', PopCnt(IEOR(Alpha_String_1,Alpha_String_2))
!      Write(IOut,*)
!      Write(IOut,*) 'beta 1:'
!      Write(IOut,'(B64)') Beta_String_1
!      Write(IOut,*) 'beta 2:'
!      Write(IOut,'(B64)') Beta_String_2
!      Write(IOut,*) 'beta XOR results in slater condon'
!      Write(IOut,'(B64)') IEOR(Beta_String_1,Beta_String_2)   
!      Write(IOut,*) 'NDifb=', PopCnt(IEOR(Beta_String_1,Beta_String_2))
!      Write(IOut,*)
!
      NBit_Ints = (NBasis/Bit_Size(0))+1 
      Allocate(Alpha_Diff(NBit_Ints),Beta_Diff(NBit_Ints))
      Det_Diff = 0
      Alpha_Diff_Cnt = 0
      Beta_Diff_Cnt = 0
      NAlpha = 0
      NBeta = 0
      Do I = 1,NBit_Ints
        Alpha_Diff(I) = IEOR(Alpha_String_1(I),Alpha_String_2(I))
!        Write(IOut,*) 'Alpha Diff',I,':'
!        Write(IOut,'(B64)') Alpha_Diff(I)
!        Write(IOut,*) '-------------'
        Alpha_Diff_Cnt = Alpha_Diff_Cnt + PopCnt(Alpha_Diff(I)) 
        NAlpha = NAlpha + PopCnt(Alpha_String_1(I))
        Beta_Diff(I) = IEOR(Beta_String_1(I),Beta_String_2(I))
!        Write(IOut,*) 'Beta Diff',I,':'
!        Write(IOut,'(B64)') Beta_Diff(I)
!        Write(IOut,*) '-------------'
        Beta_Diff_Cnt = Beta_Diff_Cnt + PopCnt(Beta_Diff(I))
        NBeta = NBeta + PopCnt(Beta_String_1(I))
      EndDo
!      Write(IOut,*)'Alpha_Diff_Cnt:',Alpha_Diff_Cnt,'Beta_Diff_Cnt:',Beta_Diff_Cnt
      Det_Diff = Alpha_Diff_Cnt/2 + Beta_Diff_Cnt/2

      If(Mod(Alpha_Diff_Cnt,2).ne.0.or.Mod(Beta_Diff_Cnt,2).ne.0) then
        Write(IOut,*) "ERROR: S2_Mat_Elem has been handed spin non-conserving &
        determinants"
        Call Exit()
      EndIf
      Select Case (Det_Diff)
!
        Case(3:)
          S2_Mat_Elem = Zero 
          Return
!
        Case(2)
          IDiff = 1
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(BTest(Alpha_Diff(I),J).eq..True.) then
              Orbs(IDiff) = IPos+1
              Spin(IDiff) = 0
              If(BTest(Alpha_String_1(I),J).eq..True.) then
                Det(IDiff) = 1
              Else 
                Det(IDiff) = 2
              EndIf
              IDiff = IDiff + 1
            EndIf
            If(BTest(Beta_Diff(I),J).eq..True.) then
              Orbs(IDiff) = IPos+1
              Spin(IDiff) = 1
              If(BTest(Beta_String_1(I),J).eq..True.) then
                Det(IDiff) = 1
              Else 
                Det(IDiff) = 2
              EndIf
              IDiff = IDiff + 1
            EndIf
          EndDo
!          Call Print_Vector(IOut,Orbs,'Orbs')
!          Call Print_Vector(IOut,Spin,'Spin')
!          Call Print_Vector(IOut,Det,'Det')
!
          IOcc = 0
          Do IPos = Orbs(1)-1, 0, -1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(Spin(1).eq.0) then
              If(Det(1).eq.1) then
                If(BTest(Alpha_String_1(I),J).eq..True.) IOcc = IOcc + 1
              ElseIf(Det(1).eq.2) then
                If(BTest(Alpha_String_2(I),J).eq..True.) IOcc = IOcc + 1
              EndIf
            ElseIf(Spin(1).eq.1) then
              If(Det(1).eq.1) then
                If(BTest(Beta_String_1(I),J).eq..True.) IOcc = IOcc + 1
              ElseIf(Det(1).eq.2) then
                If(BTest(Beta_String_2(I),J).eq..True.) IOcc = IOcc + 1
              EndIf
            EndIf
          EndDo
!          Write(IOut,*) 'IOcc:',IOcc
          JOcc = 0
          Do IPos = Orbs(2)-1, 0, -1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(Spin(2).eq.0) then
              If(Det(2).eq.1) then
                If(BTest(Alpha_String_1(I),J).eq..True.) JOcc = JOcc + 1
              ElseIf(Det(2).eq.2) then
                If(BTest(Alpha_String_2(I),J).eq..True.) JOcc = JOcc + 1
              EndIf
            ElseIf(Spin(2).eq.1) then
              If(Det(2).eq.1) then
                If(BTest(Beta_String_1(I),J).eq..True.) JOcc = JOcc + 1
              ElseIf(Det(2).eq.2) then
                If(BTest(Beta_String_2(I),J).eq..True.) JOcc = JOcc + 1
              EndIf
            EndIf
          EndDo
!          Write(IOut,*) 'JOcc:',JOcc
          KOcc = 0
          Do IPos = Orbs(3)-1, 0, -1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(Spin(3).eq.0) then
              If(Det(3).eq.1) then
                If(BTest(Alpha_String_1(I),J).eq..True.) KOcc = KOcc + 1
              ElseIf(Det(3).eq.2) then
                If(BTest(Alpha_String_2(I),J).eq..True.) KOcc = KOcc + 1
              EndIf
            ElseIf(Spin(3).eq.1) then
              If(Det(3).eq.1) then
                If(BTest(Beta_String_1(I),J).eq..True.) KOcc = KOcc + 1
              ElseIf(Det(3).eq.2) then
                If(BTest(Beta_String_2(I),J).eq..True.) KOcc = KOcc + 1
              EndIf
            EndIf
          EndDo
!          Write(IOut,*) 'KOcc:',KOcc
          LOcc = 0
          Do IPos = Orbs(4)-1, 0, -1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(Spin(4).eq.0) then
              If(Det(4).eq.1) then
                If(BTest(Alpha_String_1(I),J).eq..True.) LOcc = LOcc + 1
              ElseIf(Det(4).eq.2) then
                If(BTest(Alpha_String_2(I),J).eq..True.) LOcc = LOcc + 1
              EndIf
            ElseIf(Spin(4).eq.1) then
              If(Det(4).eq.1) then
                If(BTest(Beta_String_1(I),J).eq..True.) LOcc = LOcc + 1
              ElseIf(Det(4).eq.2) then
                If(BTest(Beta_String_2(I),J).eq..True.) LOcc = LOcc + 1
              EndIf
            EndIf
          EndDo
!          Write(IOut,*) 'LOcc:',LOcc
!          Mat_Sign = -1
!          Mat_Sign = (-1)**(IOcc+JOcc+KOcc+LOcc-3)
!          Write(IOut,*) 'Permutations:',(IOcc+JOcc+KOcc+LOcc-3)
          Mat_Sign = (-1)**(2*(NAlpha+NBeta)+1-IOcc-JOcc-KOcc-LOcc)
!          Write(IOut,*) 'Permutations:',(2*(NAlpha+NBeta)+1-IOcc-JOcc-KOcc-LOcc)
!          Write(IOut,*) 'Mat_Sign:',Mat_Sign
!
          If(Det(1).eq.Det(2).and.Det(3).eq.Det(4)) then
            If(Spin(1).eq.Spin(3).and.Spin(2).eq.Spin(4)) then
              If(Spin(1).eq.0.and.Spin(2).eq.1) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(1),Orbs(4)+NBasis)*MO_Overlap(Orbs(2)+NBasis,Orbs(3))
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.0) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(2),Orbs(3)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(4))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Elem = Zero
              EndIf
            ElseIf(Spin(1).eq.Spin(4).and.Spin(2).eq.Spin(3)) then
              If(Spin(1).eq.0.and.Spin(2).eq.1) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(1),Orbs(3)+NBasis)*MO_Overlap(Orbs(2)+NBasis,Orbs(4))
              ElseIf(Spin(1).eq.1.and.Spin(2).eq.0) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(2),Orbs(4)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(3))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Elem = Zero
              EndIf
            Else
!             This suggests that there are unbalanced spins between determinants 
              S2_Mat_Elem = Zero
            EndIf
          ElseIf(Det(1).eq.Det(3).and.Det(2).eq.Det(4)) then
            If(Spin(1).eq.Spin(2).and.Spin(3).eq.Spin(4)) then
              If(Spin(1).eq.0.and.Spin(3).eq.1) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(1),Orbs(4)+NBasis)*MO_Overlap(Orbs(3)+NBasis,Orbs(2))
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.0) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(3),Orbs(2)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(4))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Elem = Zero
              EndIf
            ElseIf(Spin(1).eq.Spin(4).and.Spin(2).eq.Spin(3)) then
              If(Spin(1).eq.0.and.Spin(3).eq.1) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(1),Orbs(2)+NBasis)*MO_Overlap(Orbs(3)+NBasis,Orbs(4))
              ElseIf(Spin(1).eq.1.and.Spin(3).eq.0) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(3),Orbs(4)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(2))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Elem = Zero
              EndIf
            Else
!             This suggests that there are unbalanced spins between determinants 
              S2_Mat_Elem = Zero
            EndIf
          ElseIf(Det(1).eq.Det(4).and.Det(2).eq.Det(3)) then
            If(Spin(1).eq.Spin(2).and.Spin(3).eq.Spin(4)) then
              If(Spin(1).eq.0.and.Spin(4).eq.1) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(1),Orbs(3)+NBasis)*MO_Overlap(Orbs(4)+NBasis,Orbs(2))
              ElseIf(Spin(1).eq.1.and.Spin(4).eq.0) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(4),Orbs(2)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(3))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Elem = Zero
              EndIf
            ElseIf(Spin(1).eq.Spin(3).and.Spin(2).eq.Spin(4)) then
              If(Spin(1).eq.0.and.Spin(4).eq.1) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(1),Orbs(2)+NBasis)*MO_Overlap(Orbs(4)+NBasis,Orbs(3))
              ElseIf(Spin(1).eq.1.and.Spin(4).eq.0) then
                S2_Mat_Elem = Mat_Sign*MO_Overlap(Orbs(4),Orbs(3)+NBasis)*MO_Overlap(Orbs(1)+NBasis,Orbs(2))
              Else
!             Setting anything that isn't alpha-beta --> alpha-beta excitations to zero
                S2_Mat_Elem = Zero
              EndIf
            Else
!             This suggests that there are unbalanced spins between determinants 
              S2_Mat_Elem = Zero
            EndIf
          EndIf
!          Write(IOut,*) 'S2_Mat_Elem:',S2_Mat_Elem

          Return

        Case(1)
          IDiff = 1
!          Allocate(Orbs(2),Spin(2))
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            If(BTest(Alpha_Diff(I),J).eq..True.) then
              Orbs(IDiff) = IPos+1
              Spin(IDiff) = 0
              If(BTest(Alpha_String_1(I),J).eq..True.) then
                Det(IDiff) = 1
              Else 
                Det(IDiff) = 2
              EndIf
              IDiff = IDiff + 1
            EndIf
            If(BTest(Beta_Diff(I),J).eq..True.) then
              Orbs(IDiff) = IPos+1
              Spin(IDiff) = 1
              If(BTest(Beta_String_1(I),J).eq..True.) then
                Det(IDiff) = 1
              Else 
                Det(IDiff) = 2
              EndIf
              IDiff = IDiff + 1
            EndIf
          EndDo
!          Write(IOut,*)'Orb 1:',Orbs(1),' Orb 2:',Orbs(2)
!          Write(IOut,*)'Spin 1:',Spin(1),' Spin 2:',Spin(2)
!          Write(IOut,*)'Det 1:',Det(1),' Det 2:',Det(2)
!
          S2_Mat_Elem = Zero 
          If(Spin(1).ne.Spin(2)) then
            S2_Mat_Elem = Zero
!
          ElseIf(Spin(1).eq.0) then
!
            IOcc = 0
            Do IPos = Orbs(1)-1, 0, -1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(Det(1).eq.1) then
                If(BTest(Alpha_String_1(I),J).eq..True.) IOcc = IOcc + 1
              ElseIf(Det(1).eq.2) then
                If(BTest(Alpha_String_2(I),J).eq..True.) IOcc = IOcc + 1
              EndIf
            EndDo
!            Write(IOut,*) 'IOcc:',IOcc
            JOcc = 0
            Do IPos = Orbs(2)-1, 0, -1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(Det(2).eq.1) then
                If(BTest(Alpha_String_1(I),J).eq..True.) JOcc = JOcc + 1
              ElseIf(Det(2).eq.2) then
                If(BTest(Alpha_String_2(I),J).eq..True.) JOcc = JOcc + 1
              EndIf
            EndDo
!            Write(IOut,*) 'JOcc:',JOcc
!
            Do IPos = 0, NBasis-1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(BTest(Beta_String_1(I),J).eq..True.) then
                S2_Mat_Elem = S2_Mat_Elem + MO_Overlap(IPos+1+NBasis,Orbs(1)) * MO_Overlap(Orbs(2),IPos+1+NBasis)
              EndIf
            EndDo
!            S2_Mat_Elem = - S2_Mat_Elem
!            Write(IOut,*) 'Permutations:',(2*NAlpha+1-IOcc-JOcc)
!            Write(IOut,*) 'Mat_Sign:',(-1)**(2*NAlpha+1-IOcc-JOcc)
            S2_Mat_Elem = (-1)**(2*NAlpha+1-IOcc-JOcc) * S2_Mat_Elem
!            S2_Mat_Elem = (-1)**(IOcc+JOcc-1) * S2_Mat_Elem
!
          ElseIf(Spin(1).eq.1) then

            IOcc = 0
            Do IPos = Orbs(1)-1, 0, -1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(Det(1).eq.1) then
                If(BTest(Beta_String_1(I),J).eq..True.) IOcc = IOcc + 1
              ElseIf(Det(1).eq.2) then
                If(BTest(Beta_String_2(I),J).eq..True.) IOcc = IOcc + 1
              EndIf
            EndDo
!            Write(IOut,*) 'IOcc:',IOcc
            JOcc = 0
            Do IPos = Orbs(2)-1, 0, -1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(Det(2).eq.1) then
                If(BTest(Beta_String_1(I),J).eq..True.) JOcc = JOcc + 1
              ElseIf(Det(2).eq.2) then
                If(BTest(Beta_String_2(I),J).eq..True.) JOcc = JOcc + 1
              EndIf
            EndDo
!            Write(IOut,*) 'JOcc:',JOcc
!
            Do IPos = 0, NBasis-1
              I = NBit_Ints - IPos/Bit_Size(0)
              J = Mod(IPos,Bit_Size(0)) 
              If(BTest(Alpha_String_1(I),J).eq..True.) then
                S2_Mat_Elem = S2_Mat_Elem + MO_Overlap(IPos+1,Orbs(1)+NBasis) * MO_Overlap(Orbs(2)+NBasis,IPos+1)
              EndIf
            EndDo
!            S2_Mat_Elem = - S2_Mat_Elem
!            Write(IOut,*) 'Permutations:',(2*NBeta+1-IOcc-JOcc)
!            Write(IOut,*) 'Mat_Sign:',(-1)**(2*NBeta+1-IOcc-JOcc)
            S2_Mat_Elem = (-1)**(2*NBeta+1-IOcc-JOcc) * S2_Mat_Elem
!            S2_Mat_Elem = (-1)**(IOcc+JOcc-1) * S2_Mat_Elem

          EndIf

!          Write(IOut,*) 'S2_Mat_Elem:',S2_Mat_Elem

          Return
!
        Case(0)
          ABTerm = Zero
          Do IPos = 0, NBasis-1
            I = NBit_Ints - IPos/Bit_Size(0)
            J = Mod(IPos,Bit_Size(0)) 
            Do JPos = 0, NBasis-1
              II = NBit_Ints - JPos/Bit_Size(0)
              JJ = Mod(JPos,Bit_Size(0)) 
              If((BTest(Alpha_String_2(I),J).eq..True.).and.(BTest(Beta_String_2(II),JJ).eq..True.)) then
                ABTerm = ABTerm + MO_Overlap(IPos+1,JPos+1+NBasis)*MO_Overlap(JPos+1+NBasis,IPos+1) 
              EndIf
            EndDo
          EndDo
          S2_Mat_Elem = Quarter*((NAlpha-NBeta)**2+2*(NAlpha+NBeta)) - ABTerm
!          Write(IOut,*) 'S2_Mat_Elem:',S2_Mat_Elem
!
          Return
!
      End Select
!
      End Function S2_Mat_Elem  

      
      
      
!
!
      end module hphSpinFun_mod
