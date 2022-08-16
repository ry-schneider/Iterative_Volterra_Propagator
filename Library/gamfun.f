c $Header: /usr/home/bis/mesa/library/bismath/RCS/gamfun.f,v 1.1 2002/08/06 13:15:50 bis Exp bis $
*deck @(#)gamfun.f	1.1 9/9/91
      function  gamfun(z)
c  this is a procedure that evaluates gamma(z) for
c     0 lt z le 3 to 16 significant figures
c    it is based on a chebyshev-type polynomial
c   approximation given in h. werner and r. collinge, math. comput.
c    15 (1961), pp. 195-97.
c   approximations to the gamma function, accurate up to 18 significant
c   digits, may be found in the paper quoted above
c
c
c
        implicit real*8 (a-h,o-z)
        dimension  a(18)
c
        a(1)=1.0d0
        a(2)=.4227843350984678d0
        a(3)=.4118403304263672d0
        a(4)=.0815769192502609d0
        a(5)=.0742490106800904d0
        a(6)=-.0002669810333484d0
        a(7)=.0111540360240344d0
        a(8)=-.0028525821446197d0
        a(9)=.0021036287024598d0
        a(10)=-.0009184843690991d0
        a(11)=.0004874227944768d0
        a(12)=-.0002347204018919d0
        a(13)=.0001115339519666d0
        a(14)=-.0000478747983834d0
        a(15)=.0000175102727179d0
        a(16)=-.0000049203750904d0
        a(17)=.0000009199156407d0
        a(18)=-.0000000839940496d0
c
c
c
        if(z.le.1.0d0  ) go to 10
        if(z.le.2.0d0  ) go to 20
        t=z-2.0d0
        go to 30
10    t=z
        go to 30
20    t=z-1.0d0
30    p=a(18)
        do 40 k1=1,17
        k=18-k1
        p=t*p+a(k)
40    continue
c
        if(z.gt.2.0d0  ) go to 50
        if(z.gt.1.0d0  ) go to 60
        gamfun=p/(z*(z+1.0d0  ))
        return
60    gamfun=p/z
        return
50    gamfun=p
        return
        end
