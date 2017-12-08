! rannw.f - normal random number generator - version 1.02 of 28 Oct 2001
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      C
!  Copyright (C) 2001 R. P. Brent.                                     C
!                                                                      C
!  This program is free software; you can redistribute it and/or       C
!  modify it under the terms of the GNU General Public License,        C
!  version 2, June 1991, as published by the Free Software Foundation. C
!  For details see http://www.gnu.org/copyleft/gpl.html .              C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
        subroutine rannw (am, sd, ix, da, n, dwork, nwork)
!       implicit none
        integer ix, n, nwork
        double precision am, sd, da(*), dwork(*)
!
! Returns n normally distributed random numbers in da(1) .. da(n),
! using Wallace's method as modified by Brent.
!
! For a description of the method and other references, see:
!
!   R. P. Brent, "Random number generation and simulation on vector and 
!   parallel computers", Proc. Fourth International Euro-Par Conference,
!   Lecture Notes in Computer Science, Vol. 1470, Springer-Verlag, 
!   Berlin, 1998, 1-20.
!
! Time per normally distributed number on the Fujitsu VP2200/10
! is approximately (6.8*f + 3.2) nsec (i.e. (1.8*f + 1.0) cycles)
! if n and nwork are large (e.g. n = 100000, nwork = 660000).
! Here f is a small positive integer which can be as small as 1
! but for statistical reasons we prefer f = 3 or 4.
! To obtain the fastest possible generator, at the cost of some
! possible loss of quality, set f to 1 in the data statement below.
!
! With f = 1 the numbers should still be normally distributed, but there 
! may be some dependence between consecutive blocks of n numbers. This
! could be acceptable (depending on the application) if the blocks are 
! large.
!
! rannw calls rannwb (see below) and also requires
! a uniform random number generator (currently ranut).
!
! am is the mean of the normal distribution.
!
! sd is the standard deviation (should be positive).
!
! ix is a nonzero seed provided by the caller.
! On the first call, ix should be positive.  It is returned as zero
! and should remain zero for subsequent calls.
!
! da(*) is a double-precision array of dimension at least n.
! rannw fills da(1), da(2), ... , da(n) with pseudo-random number
! normally distributed with mean am, standard deviation sd.
! da need not be preserved between calls (i.e. da(1) .. da(n) can
! be modified by the caller).
!
! n is the number of random numbers to return.
! The dimension of da must be at least n.
! For efficiency n should be moderately large.
!
! n, am and sd may be different on each call to rannw.
! If n is not positive, no random numbers are returned.
!
! dwork(*) is a double-precision array of dimension at least nwork
! used as a work area and to save state information.
! dwork(1) .. dwork(nwork) must be preserved between calls to rannw.
!
! nwork must be at least 1350.
!
! For efficiency and the best statistical properties
! nwork should be large (preferably in the range 50000 to 1000000,
! depending on memory limitations).
!
! R. P. Brent, last revised 20011028
!
        integer f, j, jlim, jlow, flag, lbp1, lbp2, lbpf, nloc, nd, usiz
        double precision r, t, pi2, ssq, scale
        integer nbuf, nz2, uind, x1ind, y1ind, x2ind, y2ind, z2ind
        double precision checkt, errort
!
! f is an important "quality factor" (or "fudge factor").
!
! Only a fraction 1/f of numbers generated are returned to the user,
! the others are discarded. This improves the performance on statistical
! tests but has an obvious cost in time per number returned.
! f should be in the range 1 to 8 (recommended value is 2, 3, or 4).
!
        data f /4/
!
! To avoid "drift" due to rounding errors, a sum of squares is
! recomputed occasionally.   checkt determines how often this is done.
! The smaller checkt, the less often sum of squares is recomputed
! (in fact it is recomputed about once in 1/checkt calls to rannwb).
!
        data checkt /1.0D-3/
!
! If relative error exceeds errort regard as overwriting error
!
        data errort /1.0D-6/
!
! nz2 should be at least 58 and no larger than nbuf
!
        data nz2 /63/
!
!        write(6,*) ' in rannw, am,sd',  am, sd
        if (ix .ne. 0) then
!         print *, ' Using quality factor f = ', f, ' in rannw'
!
! First call or restart so determine buffer size nbuf
! which must be a power of two
!


          nbuf = 64
10        nbuf = 2*nbuf
          if (nbuf .le. (nwork-(nz2+7))/10) go to 10
! Check for nwork too small
          if (nbuf .lt. 256) then
            print *, ' Error 1 in rannw - nwork = ', nwork, &
      &   ' should be at least 1350'
            stop
            endif
          x1ind = 7
          y1ind = x1ind + nbuf
          x2ind = y1ind + nbuf
          y2ind = x2ind + nbuf
          z2ind = y2ind + nbuf
          uind  = z2ind + nz2
          usiz  = nwork - uind
! Call ranut to get 2*nbuf uniform numbers in dwork(x2ind+1..)
! dwork(uind+1..uind+usiz) used for ranut workspace (preserve in
! case a restart is required).
! ix is used as the seed for ranut.
! A different uniform random number generator could be used here.
          call ranut (ix, dwork(x2ind+1), 2*nbuf, dwork(uind+1), usiz)
          pi2 = 8.0d0*datan(1.0d0)
! Generate first 2*nbuf normal numbers by (slow) Box-Muller method.
! Also compute sum of squares
          ssq = 0.0d0
! Following loop vectorises
          do 20 j = 1, nbuf
          t = -2.0d0*dlog(1.0d0 - dwork(x2ind+j))
          ssq = ssq + t
          r = dsqrt(t)
          dwork(x1ind+j) = r*dsin(pi2*dwork(y2ind+j))
          dwork(y1ind+j) = r*dcos(pi2*dwork(y2ind+j))
20        continue
! flag indicates whether current normal numbers in (x1,y1) or (x2,y2)
! For convenience of indexing, flag = 0 or 2*nbuf
          flag = 0
! Locations lbp1+1 to lbp2 of (x1,y1) or (x2,y2) are still to be used
! (x1(1) and x2(1) are reserved for use by chi-square correction).
          lbp1 = 1
          lbp2 = 2*nbuf
! Save essential variables in first 4 words of dwork for later calls
          dwork(1) = dble(nbuf)
          dwork(2) = dble(lbp1)
          dwork(3) = dble(lbp2)
          dwork(4) = dble(flag)
! For later convenience save sqrt(sum of squares), A and B,
! where A**2 = 1/(1 + sqrt(1 - 1/nu)), B**2 = sqrt(nu(nu - 1)),
! nu = lbp2 = 2*nbuf is the number of degrees of freedom
! and (A*normal + B)**2 is a good approximation to chi-squared
! (it has the correct mean nu and variance 2*nu).
          dwork(5) = dsqrt(ssq)
          dwork(6) = dsqrt(1.0d0/(1.0d0 + &
     &                 dsqrt(1.0d0 - 1.0d0/dble(lbp2))))
          dwork(7) = dsqrt(dsqrt(dble(lbp2)*dble(lbp2-1)))
          endif
! Usual entry - restore variables saved in dwork
        nbuf = int(dwork(1))
        lbp1 = int(dwork(2))
        lbp2 = int(dwork(3))
        flag = int(dwork(4))
! Check consistency of dwork(1) to dwork(7)
        if ((nbuf .lt. 256) .or. (lbp2 .ne. 2*nbuf) .or. &
     &      (lbp1 .le. 0)   .or. (lbp1 .ge. lbp2) .or.   & 
     &      ((flag .ne. 0) .and. (flag .ne. lbp2)) .or.  &
     &      (dwork(5) .le. 0d0) .or.                     &
     &      (dwork(6) .gt. 0.7073d0) .or.                &
     &      (dwork(6) .lt. 0.7071d0) .or.                &
     &      (dwork(7) .lt. 22.61d0)) then
          print *, ' Error 2 in rannw - dwork overwritten', &
     &             ' internal consistency check failed'
          stop
          endif
! Reconstruct indices for "virtual" arrays x1, y1, x2, y2, z2
        x1ind = 7
        y1ind = x1ind + nbuf
        x2ind = y1ind + nbuf
        y2ind = x2ind + nbuf
        z2ind = y2ind + nbuf
        uind  = z2ind + nz2
        usiz  = nwork - uind
! nloc is local copy of n - reduced down to zero.
! nd is index into da array - increased as nloc is decreased.
        nloc = n
        nd = 0
! Only return a fraction 1/f of numbers generated by rannwb
        lbpf = lbp2/f
! Loop if necessary to fill da(1)..da(n)
30      if (nloc .le. 0) then
! Save variables which may have changed and return
          dwork(2) = dble(lbp1)
          dwork(4) = dble(flag)
          return
          endif
! Else need to copy more to user buffer (adjusting mean and sd)
        jlim = min(lbpf - lbp1, nloc)
        jlow = lbp1 + x1ind + flag
!
! Following loop takes one cycle per random number on VP2200/10 -
! most of this overhead could be avoided in the case am = 0
! by working directly in the array da rather than in dwork.
! The loop vectorises.
!
          do 40 j = 1, jlim
40        da(j+nd) = am + sd*dwork(j+jlow)
!
        lbp1 = lbp1 + jlim
        nd = nd + jlim
        nloc = nloc - jlim
        if (lbp1 .ge. lbpf) then
! Call ranut to get nz2 uniform numbers in dwork(z2ind+1..).
! dwork(uind+1..uind+usiz) is used for ranut workspace.
          call ranut (ix, dwork(z2ind+1), nz2, dwork(uind+1), usiz)
! Use (A*normal + B)**2 as approximation to chi-squared,
! where A = dwork(6), B = dwork(7) precomputed (see above).
! Indexing by +- flag swaps (x1,y1) <-> (x2,y2) if flag nonzero.
          ssq = dwork(6)*dwork(x1ind+flag+1) + dwork(7)
          if (ssq .le. 0d0) ssq = dwork(5)
! Compute scale factor to give correct sum of squares
          scale = ssq/dwork(5)
! Call rannwb to generate 2*nbuf more normal numbers
          call rannwb (nbuf, scale, &
     &         dwork(x1ind+flag+1), dwork(y1ind+flag+1), &
     &         dwork(x2ind+1-flag), dwork(y2ind+1-flag), &
     &         dwork(z2ind+1), nz2)
! Adjust sqrt(sum of squares).
          dwork(5) = ssq
! For a double-precision version this recomputation is hardly
! necessary, but is a check on overwriting of dwork so do it
! occasionally, depending on a pseudo-random uniform
          if (dwork(z2ind+1) .lt. checkt) then
            ssq = 0d0
! Following loop vectorises
            do 50 j = 1, nbuf
            ssq = ssq + dwork(j+x2ind-flag)**2 + &
     &                  dwork(j+y2ind-flag)**2
50          continue
            ssq = dsqrt(ssq)
! Error return if error seems too large to be accounted for
! by rounding errors so probably caused by overwriting buffer.
            if (dabs(ssq - dwork(5)) .ge. errort*ssq) then
              print *, ' Error 3 in rannw - buffer probably overwritten'
              print *, ' Relative error in norm ', (ssq-dwork(5))/ssq, &
     &                 ' is too large'
              stop
              endif
! Reset dwork(5) to accurate value.
            dwork(5) = ssq
            endif
! Reverse flag (alternates between 0 and lbp2 = 2*nbuf)
          flag = lbp2 - flag
          lbp1 = 1
          endif
        go to 30
        end
!
        subroutine rannwb (n, scale, x1, y1, x2, y2, z2, nz2)
!       implicit none
        integer n, nz2
        double precision scale, x1(*), y1(*), x2(*), y2(*), z2(*)
!
! This routine is called by rannw.
! Assumes 2*n normal numbers in x1(1..n), y1(1..n)
! and nz2 uniform numbers in z2(1..nz2),
! returns 2*n normal numbers in x2(1..n), y2(1..n)
! scaled by given scale factor.
! nz2 should be at least 3*max(alpha+beta) + 9, i.e. 57
!
! In call by rannw, x1(1) and x2(1) are reserved for determination
! of scale factor by chi-squared approximation, so should not be
! used independently.
!
! R. Brent, last revised 20011028
!
        integer j, j0, j1, k1, k2, k3
        integer alpha, beta, gamma, delta
        double precision t1, t2, a00, a01, t, tr
        double precision c0, c1
        c0 = 1.0d0/(2.0d0 + dsqrt(3.0d0))
        c1 = 1.0d0/dsqrt(3.0d0) - c0
        k3 = nz2
! Simple random choice of alpha, beta, gamma, delta.
! alpha, beta should be odd and greater than 1,
! gamma and delta should be nonnegative.
        k3 = k3 - 4
        if (z2(k3+4) .gt. 0.5d0) then
          alpha = 3
        else
          alpha = 9
          endif
        k3 = k3 - 1
        if (z2(k3+3) .gt. 0.5d0) then
          beta = 5
        else
          beta = 7
          endif
        gamma = int(dble(n)*z2(k3+2))
        delta = int(dble(n)*z2(k3+1))
        j0 = 0
! while (j0 .lt. n) ...
! Compute j1, k1, k2 so no "mod n" operations on indices in loop
! (so conceptual loop is split into about alpha + beta smaller loops)
10      k1 = gamma + 1 - n*((alpha*j0 + gamma)/n)
        k2 = delta + 1 - n*((beta*j0  + delta)/n)
        j1 = min((n - k1)/alpha, (n - k2)/beta)
        j1 = min(n-1, j1)
! Use uniform to compute pseudo-random scaled 2x2 orthogonal matrix.
! Matrix is scale*(orthogonal matrix defined by t = tan(theta/2)),
! where theta is in (pi/6, pi/3) so sin(theta), cos(theta) .ge. 0.5
! Following should not occur unless nz2 is too small.
        if (k3 .lt. 4) k3 = nz2
        k3 = k3 - 3
! Force t in [tan(pi/12), tan(pi/6)] to ensure mixing.
        t  = c0 + c1*z2(k3+3)
        tr = scale/(1.0d0 + t*t)
! If t = tan(theta/2) then a00 = scale*cos(theta)
        a00 = (1.0d0 - t*t)*tr
! Randomise sign of cos
        if (z2(k3+2) .lt. 0.5d0) a00 = -a00
! If t = tan(theta/2) then a01 = scale*sin(theta)
        a01 = 2.0d0*t*tr
! Randomise sign of sin
        if (z2(k3+1) .lt. 0.5d0) a01 = -a01
!
! This is the inner loop for normal random number generation
! (apart from linear transformation and copying to user buffer).
! The non-unit stride access is not too slow because we choose
! alpha and beta to be odd. The loop vectorises.
!
        do 20 j = j0, j1
          t1 = x1(alpha*j + k1)
          t2 = y1(beta*j  + k2)
          x2(j+1) = a00*t1 + a01*t2
          y2(j+1) = a00*t2 - a01*t1
20        continue
!
        j0 = j1 + 1
        if (j0 .lt. n) go to 10
        return
        end
!
! End of rannw.f
! =====================================================================
! Add ranut and associated routines start here.
!
