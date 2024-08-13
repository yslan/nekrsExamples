c-----------------------------------------------------------------------
      subroutine my_cbc_chk(s3)
c     This print types of current BC, counting as well
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer ncbc,lcbc,i,j,e,f,inid,mtype,idummy,ncbc_tmp,ifld,nel
     $      , nfldt
      parameter (lcbc=10) ! max # unique CBC
      integer cbc_count(lcbc),cbc_count_tmp(lcbc)
      character*3 cbc_list(lcbc),cbc_list_tmp(lcbc),cbc3,s3,tag3
      character*2 s2tmp
      logical cbc_in_list

      nfldt = nfield
      if (ifmhd) nfldt = ifldmhd   ! mhd always uses nscal+3, with or without temperature

      do ifld=0,nfldt

        if (ifld.eq.0) then
          tag3='PR '
          nel = nelv
        elseif (ifld.eq.1) then
          tag3='VEL'
          nel = nelv
        elseif (ifmhd.AND.ifld.eq.ifldmhd) then
          tag3='MHD'
          nel = nelv
        elseif (ifld.gt.1) then
          write(s2tmp,'(I2.2)') ifld-2
          tag3='S'//s2tmp
          nel = nelt
        endif

        call izero(cbc_count,lcbc)

        ncbc = 1                 ! # unique CBC
        cbc_list(1) = 'E  ' ! internal BC

        do e=1,nel
        do f=1,2*ldim
          cbc3 = cbc(f,e,ifld)
          cbc_in_list = .false.

          if (cbc3.eq.'E  '.or.cbc3.eq.'   ') then
            cbc_count(1) = cbc_count(1) + 1
            cbc_in_list = .true.
          else
            do i=2,ncbc ! go throught the registered CBC
              if (cbc3.eq.cbc_list(i)) then
                cbc_count(i) = cbc_count(i) + 1
                cbc_in_list = .true.
              endif
            enddo
          endif
          if (.not.cbc_in_list) then
            ncbc = ncbc + 1
            if (ncbc.gt.lcbc)
     $        call exitti('BCchk: Increase lcbc$',ncbc)

            cbc_list(ncbc) = cbc3
            cbc_count(ncbc) = 1
          endif
        enddo
        enddo

        call nekgsync()

        ! All reduce to nid=0
        if (nid.eq.0) then
          do inid=1,np-1
            mtype = inid
            call csend(mtype,idummy,4,inid,0) ! handshake

            call crecv(mtype,ncbc_tmp,4)
            call crecv(mtype,cbc_list_tmp(1),3*ncbc_tmp,0,0)
            call crecv(mtype,cbc_count_tmp(1),4*ncbc_tmp,0,0)

            cbc_count(1) = cbc_count(1)+cbc_count_tmp(1)
            do j=2,ncbc_tmp
              cbc_in_list = .false.
              do i=2,ncbc
                if (cbc_list(i).eq.cbc_list_tmp(j)) then
                  cbc_count(i) = cbc_count(i) + cbc_count_tmp(j)
                  cbc_in_list = .true.
                endif
              enddo
              if (.not.cbc_in_list) then
                ncbc = ncbc + 1
                if (ncbc.gt.lcbc)
     $            call exitti('BCchk: Increase lcbc$',ncbc)

                cbc_list(ncbc) = cbc_list_tmp(j)
                cbc_count(ncbc) = cbc_count_tmp(j)
              endif
            enddo
          enddo
        else
          mtype = nid
          call crecv(mtype,idummy,4) ! ! handshake

          call csend(mtype,ncbc,4,0,0)
          call csend(mtype,cbc_list,3*ncbc,0,0)
          call csend(mtype,cbc_count,4*ncbc,0,0)
        endif

        call nekgsync()

        ! print
        if (nid.eq.0) then
          do i=1,ncbc
            write(6,41) s3,ifld,tag3,cbc_list(i),cbc_count(i)
          enddo
        endif

      enddo ! ifld

   41 format(a3,' BC: 'i3,' ',a3,a5,i10)

      return
      end
c-----------------------------------------------------------------------
      subroutine my_bdid_chk(s3)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      
      integer nbid,lcbc,bid,i,j,e,f,inid,mtype,idummy,nbid_tmp,nel,ifld
      parameter (lcbc=10) ! max # unique CBC
      integer bid_count(lcbc),bid_count_tmp(lcbc)
      integer bid_list(lcbc),bid_list_tmp(lcbc),cbc3
      logical bid_in_list
      character*3 s3
      character*5 tag5
      
      ! pass 1: bdid
      do ifld=1,2
        
        nel = nelv
        tag5 = "bdidv"
        if (ifld.eq.2) then
          nel=nelt
          tag5 = "bdidt"
        endif

        call izero(bid_count,lcbc)
        
        nbid = 1         ! # unique CBC
        bid_list(1) = 0 ! internal BC
        
        do e=1,nel
        do f=1,2*ldim

          bid = boundaryID(f,e)
          if (ifld.eq.2) bid = boundaryIDt(f,e)  

          bid_in_list = .false.          
          if (bid.eq.0) then
            bid_count(1) = bid_count(1) + 1
            bid_in_list = .true.
          else 
            do i=2,nbid ! go throught the registered CBC
              if (bid.eq.bid_list(i)) then 
                bid_count(i) = bid_count(i) + 1
                bid_in_list = .true.
              endif
            enddo
          endif
          if (.not.bid_in_list) then
            nbid = nbid + 1
            if (nbid.gt.lcbc)
     $        call exitti('BCchk: Increase lcbc$',nbid)

            bid_list (nbid) = bid
            bid_count(nbid) = 1
          endif
        enddo
        enddo

        call nekgsync()

        ! All reduce to nid=0
        if (nid.eq.0) then
          do inid=1,np-1
            mtype = inid
            call csend(mtype,idummy,isize,inid,0) ! handshake
            
            call crecv(mtype,nbid_tmp,isize)
            call crecv(mtype,bid_list_tmp(1), isize*nbid_tmp)
            call crecv(mtype,bid_count_tmp(1),isize*nbid_tmp)
            
            bid_count(1) = bid_count(1)+bid_count_tmp(1)
            do j=2,nbid_tmp 
              bid_in_list = .false.
              do i=2,nbid
                if (bid_list(i).eq.bid_list_tmp(j)) then
                  bid_count(i) = bid_count(i) + bid_count_tmp(j)
                  bid_in_list = .true.
                endif
              enddo
              if (.not.bid_in_list) then
                nbid = nbid + 1
                if (nbid.gt.lcbc)
     $            call exitti('BCchk: Increase lcbc$',nbid)
                
                bid_list(nbid) = bid_list_tmp(j)
                bid_count(nbid) = bid_count_tmp(j)
              endif
            enddo
          enddo
        else
          mtype = nid
          call crecv(mtype,idummy,isize) ! ! handshake
          
          call csend(mtype,nbid,isize,0,0)
          call csend(mtype,bid_list,isize*nbid,0,0)
          call csend(mtype,bid_count,isize*nbid,0,0)
        endif
   
        call nekgsync()
   
        ! print
        if (nid.eq.0) then
          do i=1,nbid
            write(6,41) s3,tag5,bid_list(i),bid_count(i)
          enddo
        endif

      enddo !ifld

   41 format(a3,' BCid: ',a5,i5,i10)
      return
      end
c---------------------------------------------------------------------
      subroutine my_dump_bc(s3)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer e, f, bid
      real tmp1(lx1*ly1*lz1*lelt)
     $   , tmp2(lx1*ly1*lz1*lelt)
     $   , tmp3(lx1*ly1*lz1*lelt)
     $   , ttt (lx1*ly1*lz1,lelt,3)

      character*3 cb3, s3

      call rzero(tmp1,lx1*ly1*lz1*nelv)
      call rzero(tmp2,lx1*ly1*lz1*nelv)
      call rzero(ttt(1,1,1),lx1*ly1*lz1*nelt)
      call rzero(ttt(1,1,2),lx1*ly1*lz1*nelt)
      call rzero(ttt(1,1,3),lx1*ly1*lz1*nelt)

      do e=1,nelv
      do f=1,2*ldim
        cb3 = cbc(f,e,1)
        if (cb3.eq.'W  ') then
          call facev (tmp1,e,f,1.0,lx1,ly1,lz1)
        elseif (cb3.eq.'v  '.OR.cb3.eq.'V  ') then
          call facev (tmp1,e,f,2.0,lx1,ly1,lz1)
        elseif (cb3.eq.'o  '.OR.cb3.eq.'O  ') then
          call facev (tmp1,e,f,3.0,lx1,ly1,lz1)
        elseif (cb3.eq.'SYM') then
          call facev (tmp1,e,f,4.0,lx1,ly1,lz1)
        endif

        bid = boundaryID(f,e)
        if (bid.gt.0) call facev (tmp2,e,f,real(bid),lx1,ly1,lz1)
      enddo
      enddo

      do e=1,nelt
      do f=1,2*ldim
        cb3 = cbc(f,e,2)
        if (cb3.eq.'t  ') then
          call facev (ttt(1,1,1),e,f,1.0,lx1,ly1,lz1)
        elseif (cb3.eq.'I  '.OR.cb3.eq.'O  ') then
          call facev (ttt(1,1,1),e,f,2.0,lx1,ly1,lz1)
        elseif (cb3.eq.'f  ') then
          call facev (ttt(1,1,1),e,f,3.0,lx1,ly1,lz1)
        endif

        bid = boundaryIDt(f,e)
        if (bid.gt.0) call facev (ttt(1,1,2),e,f,real(bid),lx1,ly1,lz1)
      enddo

        call cfill(ttt(1,e,3),real(lglel(e)),lx1*ly1*lz1)
      enddo

c      call outpost(tmp1,tmp2,vz,pr,ttt,s3)
      call outpost2(tmp1,tmp2,vz,pr,ttt,3,s3)

      return
      end
c---------------------------------------------------------------------
      subroutine chk_pebble_rad(s3,id_peb,Rs)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      integer e,f,kx1,kx2,ky1,ky2,kz1,kz2,ix,iy,iz,ia,id_peb
      real cx,cy,cz,xd,yd,zd,r,rmin,rmax,glmin,glmax,Rs
      character*3 s3

      Rs = 0.5

      rmin = 1.E10
      rmax = -1.E10
      do e=1,nelv
      do f=1,2*ldim
         if (boundaryID(f,e).eq.id_peb) then
            ! TODO: sph center
            cx = 0.0 ! single pebble
            cy = 0.0
            cz = 0.0

            CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,lx1,ly1,lz1,f)
            do iz=lz1,lz2
            do iy=ly1,ly2
            do ix=lx1,lx2
               xd = xm1(ix,iy,iz,e) - cx
               yd = ym1(ix,iy,iz,e) - cy
               zd = zm1(ix,iy,iz,e) - cz
               r = xd**2 + yd**2 + zd**2
               if (r.gt.0) r = sqrt(r)
               rmin = min(rmin,r)
               rmax = max(rmax,r)
            enddo
            enddo
            enddo
         endif
      enddo 
      enddo 
      rmin = glmin(rmin,1)
      rmax = glmax(rmax,1)
      Rs = (rmax+rmin)/2.0

      if (nio.eq.0) then
         write(*,*)'rad (surface)',s3,id_peb,cx,cy,cz,Rs,abs(rmax-rmin)
      endif

      if (nelv.eq.nelt) return
      rmax = -1.E10
      do e=nelv+1,nelt
      do ia=1,lx1*ly1*lz1
        ! TODO: sph center
        cx = 0.0 ! single pebble
        cy = 0.0
        cz = 0.0
        xd = xm1(ia,1,1,e) - cx
        yd = ym1(ia,1,1,e) - cy
        zd = zm1(ia,1,1,e) - cz
        r = xd**2 + yd**2 + zd**2
        if (r.gt.0) r = sqrt(r)
        rmax = max(rmax,r)
      enddo
      enddo
      rmax = glmax(rmax,1)

      if (nio.eq.0) then
         write(*,*)'rad (solid)  ',s3,id_peb,cx,cy,cz,rmax,abs(Rs-rmax)
      endif

      return
      end
c---------------------------------------------------------------------
