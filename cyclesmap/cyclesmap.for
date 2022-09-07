* Calcualtes number of cycles in networks

      implicit double precision(x)
      character*90 filename,fileout,fileoutk
      parameter  (nodosmax=130000,Nedgesmax=2000000)

      dimension nconnectivity(1:Nedgesmax)
      dimension npunteroini(1:nodosmax)
      dimension npunterofin(1:nodosmax)
      dimension ndegree(1:nodosmax)
      dimension nedge(1:Nedgesmax,1:2)
      dimension nvec(1:nodosmax)
      dimension nvec1(1:nodosmax)
      dimension nvec2(1:nodosmax)

      data npunteroini/nodosmax*0/
      data npunterofin/nodosmax*0/
      data ndegree/nodosmax*0/

*Lee red de fichero y carga en vector de dos sentidos
*****************************************************
      CALL get_command_argument(1, filename)
      CALL get_command_argument(2, fileout)

      !fileoutk='test_edgelist_0_degrees.dat'
    
      open(1,file=filename,status='unknown')
      NODOS=0
      nlist=0
      do while(.true.)
        read(1,*,END=10) i,j
        ndegree(i)=ndegree(i)+1
        ndegree(j)=ndegree(j)+1
        nlist=nlist+1
        nedge(nlist,1)=i
        nedge(nlist,2)=j
        if(i.gt.NODOS)then
          NODOS=i
        endif
        if(j.gt.NODOS)then
          NODOS=j
        endif
      enddo
10    continue
      close(1)

      !      open(1,file=fileoutk,status='unknown')
            indexaux=1
            ndtot=0
            nodosreal=0
            xk=0.D0
            do i=1,NODOS
             if(ndegree(i).gt.0)then
               npunteroini(i)=indexaux
               npunterofin(i)=indexaux-1
               indexaux=indexaux+ndegree(i)
               ndtot=ndtot+ndegree(i)
               nodosreal=nodosreal+1
               xk=xk+ndegree(i)
        !       write(1,*)ndegree(i)
             endif
            enddo
       !     close(1)
            xk=xk/dble(nodosreal)

            do l=1,nlist
             i=nedge(l,1)
             j=nedge(l,2)
             npunterofin(i)=npunterofin(i)+1
             nconnectivity(npunterofin(i))=j
             npunterofin(j)=npunterofin(j)+1
             nconnectivity(npunterofin(j))=i
            enddo


*********clustering
*            write(6,*)'Calculating node clustering ...'
            nds1=0
            xckiter=0.D0
            do i=1,NODOS
            if(ndegree(i).gt.1)then
              nds1=nds1+1
              nclust=0
              do j=npunteroini(i),npunterofin(i)-1
               do k2=j+1,npunterofin(i)
                 do l=npunteroini(nconnectivity(j)),
     +              npunterofin(nconnectivity(j))
                    if(nconnectivity(l).eq.nconnectivity(k2))then
                      nclust=nclust+1
                    endif
                  enddo
                enddo
              enddo
              xckiter=xckiter
     +        +2*dble(nclust)/(dble(ndegree(i))*(dble(ndegree(i))-1.))
            endif
            enddo
            xckiter=xckiter/dble(nds1)

*********cycles of edge
            nl1=0
            squaresP=0
            pentagonsP = 0

            xceiter=0.D0
            xsqiter=0.D0
            xpiter=0.D0
            ntt=0
            nts=0
            ntp=0
            do i=1,nlist   !do edges
              nodo1=nedge(i,1)
              nodo2=nedge(i,2)
              mledge=0
              msq=0
              mpent=0

*Calculating edge clustering ...
              if((ndegree(nodo1).gt.1).and.(ndegree(nodo2).gt.1))then
              nl1=nl1+1
              nvec=0
              do j=npunteroini(nodo1),npunterofin(nodo1)
               nodo3=nconnectivity(j)
               do k2=npunteroini(nodo2),npunterofin(nodo2)
                 nodo4=nconnectivity(k2)
                 if(nodo3.eq.nodo4)then
                   mledge=mledge+1
                   nvec(mledge)=nodo3
                 endif
               enddo
              enddo
              ntt=ntt+mledge
              xceiter=xceiter+dble(mledge)/
     +         (DMIN1(dble(ndegree(nodo1)),dble(ndegree(nodo2)))-1.0D0)

*Calculating edge squares ...
              if((ndegree(nodo1).ge.(mledge+1)).and.
     +           (ndegree(nodo2).ge.(mledge+1)))then !if more than multiplicity
              
              do j=npunteroini(nodo1),npunterofin(nodo1)
               nodo3=nconnectivity(j)
               if(nodo3.ne.nodo2)then

                 nsn=0
                 do k2=1,mledge
                   if(nodo3.eq.nvec(k2))then
                     nsn=1
                     goto 40
                   endif
                 enddo
40               continue

                 if(nsn.eq.0)then
                   do k2=npunteroini(nodo2),npunterofin(nodo2)
                     nodo4=nconnectivity(k2)
                     if(nodo4.ne.nodo1)then

                       nsn2=0
                       do k3=1,mledge
                         if(nodo4.eq.nvec(k3))then
                           nsn2=1
                           goto 50
                         endif
                       enddo
50                     continue

                       if(nsn2.eq.0)then
                         do k4=npunteroini(nodo3),npunterofin(nodo3)
                           nodo5=nconnectivity(k4)
                           if(nodo5.eq.nodo4)then
                             msq=msq+1
                             goto 60
                           endif
                         enddo
60                       continue
                       endif

                     endif
                   enddo

                 endif
               endif
              enddo

              nts=nts+msq
              if((msq.gt.0))then
               nkc1=ndegree(nodo1)-mledge-1
               nkc2=ndegree(nodo2)-mledge-1
               if((nkc1.gt.0).and.(nkc2.gt.0))then
                   squaresP=squaresP+1
                   xledge=dble(msq)/(dble(nkc1)*dble(nkc2))
               else
                   xledge=0.D0
               endif
              else
                xledge=0.D0
              endif
              xsqiter=xsqiter+xledge

              endif !endif more than multiplicity

*Calculating edge pentagones ...
              if((ndegree(nodo1).ge.(mledge+1)).and.
     +           (ndegree(nodo2).ge.(mledge+1)))then !Si los dos nodos tienen vecinos que no sean comunes ni ellos

               ndsv1=0
               do j2=npunteroini(nodo1),npunterofin(nodo1) !iteramos por los vecinos de nodo1 (nodo3)
               nodo3=nconnectivity(j2)
               if(ndegree(nodo3).gt.1)then !si el grado de nodo 3 es mayor que 1 y no es el nodo2
               if(nodo3.ne.nodo2)then
                 nsn=0
                 do k3=1,mledge
                   if(nvec(k3).eq.nodo3)then !verificamos que no sea un vecino comun de nodo1 y nodo2
                     nsn=1
                     goto 70
                   endif
                 enddo
70               continue
                 nnei=0
                 if(nsn.eq.0)then
                   do k3=npunteroini(nodo3),npunterofin(nodo3) !si el nodo3 tiene un vecino
                   if(nconnectivity(k3).ne.nodo1)then ! que no es el nodo 1
                     do k4=npunteroini(nodo2),npunterofin(nodo2)
                     if(nconnectivity(k4).ne.nodo1)then
                     if(nconnectivity(k3).eq.nconnectivity(k4))then ! y que está conectado con el nodo2
                       nnei=1 !marcamos el nodo3 como malo
                       goto 71
                     endif
                     endif
                     enddo
                   endif
                   enddo
71                 continue
                 endif
                 if((nsn.eq.0).and.(nnei.eq.0))then !si el nodo3 no esta marcado como malo, o metemos en ndsv1
                   ndsv1=ndsv1+1
                   nvec1(ndsv1)=nodo3
                   !write(6,*)'vector 1 ',nodo1,nodo2,ndsv1,nodo3
                 endif
               endif
               endif
               enddo

               ndsv2=0
               do j2=npunteroini(nodo2),npunterofin(nodo2)
               nodo4=nconnectivity(j2)
               if(ndegree(nodo4).gt.1)then
               if(nodo4.ne.nodo1)then
                 nsn=0
                 do k3=1,mledge
                   if(nvec(k3).eq.nodo4)then
                     nsn=1
                     goto 80
                   endif
                 enddo
80               continue
                 nnei=0
                 if(nsn.eq.0)then
                   do k3=npunteroini(nodo4),npunterofin(nodo4)
                   if(nconnectivity(k3).ne.nodo2)then
                     do k4=npunteroini(nodo1),npunterofin(nodo1)
                     if(nconnectivity(k4).ne.nodo2)then
                     if(nconnectivity(k3).eq.nconnectivity(k4))then
                       nnei=1
                       goto 81
                     endif
                     endif
                     enddo
                   endif
                   enddo
81                 continue
                 endif
                 if((nsn.eq.0).and.(nnei.eq.0))then
                   ndsv2=ndsv2+1
                   nvec2(ndsv2)=nodo4
                   !write(6,*)'vector 2 ',nodo1,nodo2,ndsv2,nodo4
                 endif
               endif
               endif
               enddo

               do k3=1,ndsv1
               do k4=1,ndsv2 !sumamos posibles pentagonos
                 pentagonsP = pentagonsP + 
     +           MIN(ndegree(nvec1(k3))-1, ndegree(nvec1(k4))-1) 
                 do k5=npunteroini(nvec1(k3)),npunterofin(nvec1(k3))
                 do k6=npunteroini(nvec2(k4)),npunterofin(nvec2(k4))
                   if(nconnectivity(k5).eq.nconnectivity(k6))then
                     nin=0
                     do k7=1,mledge
                       if(nvec(k7).eq.nconnectivity(k5))then
                         nin=1
                         goto 90
                       endif
                     enddo
90                   continue
                     
                     if(nin.eq.0)then
                       mpent=mpent+1
                       goto 100
                     endif
                    endif
                  enddo
                  enddo
100               continue
                enddo
                enddo
                !write(6,*)ndsv1,ndsv2,mpent

               ntp=ntp+mpent
               if(mpent.gt.0)then
                 xpent=dble(mpent)/(dble(ndsv1)*dble(ndsv2))
               else
                 xpent=0.0D0
               endif
               xpiter=xpiter+xpent

              endif !endif more than m in p

              endif !endif degrees larger than 1

            enddo   !enddo edges

            xceiter=xceiter/dble(nl1)		
            xsqiterMA=xsqiter/dble(nl1)
            write(6,*)"pentagonsP"
            write(6,*)pentagonsP
            write(6,*)"nl1"
            write(6,*)nl1
            write(6,*)"xpiter"
            write(6,*)xpiter
			xsqiterP=xsqiter/dble(squaresP)
            xpiterMA=xpiter/dble(nl1)
            xpiterP=xpiter/dble(pentagonsP)
            open(1,file=fileout,status='unknown')
            write(1,*)xceiter,",",xsqiterMA,",",xsqiterP
     +           ,",",xpiterMA,",",xpiterP,",",NODOS,",",nlist
            close(1)
      end
