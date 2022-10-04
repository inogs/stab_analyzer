import numpy as np
import sys
from scipy.integrate import odeint


class LYAP(object):
    '''Compute the largest Lyapunov exponent of a timeseries'''
    def __init__(self, data):
        self.data = data
        self.datcnt = len(data)

    def basgen(self, tau, ndim, ires, maxbox):
        """
        returns ndim, ires, tau, datcnt, boxcnt, datmax, datmin,
        boxlen, datptr[:boxcnt], nxtbox[:boxcnt, :ndim]
        , where[:boxcnt,:ndim], nxtdat[:datcnt], data
        in a dictionary
        """
        delay = np.array([0,tau,(ndim-1)*tau],dtype=int)

        nxtbox = np.zeros((maxbox,ndim),dtype=int)
        where = np.zeros((maxbox,ndim),dtype=int)
        datptr =  np.full(maxbox,-1,dtype=int)
        nxtdat = np.zeros(self.datcnt,dtype=int)

        datmin = min(self.data)
        datmax = max(self.data)

        datmin = datmin - 0.01*(datmax - datmin)
        datmax = datmax + 0.01*(datmax - datmin)
        boxlen = (datmax - datmin)/ires 
    
        boxcnt = 1

        for i in range(int(self.datcnt-(ndim-1)*tau)):
            target = np.floor((self.data[i+delay]-datmin)/boxlen).astype(int)
            runner = 1
            chaser = 0

            j = 0
            while j < ndim:
                tmp = where[int(runner),j]-target[j]
                if tmp < 0:
                    chaser = runner
                    runner = nxtbox[int(runner),j]
                    if runner != 0:
                        continue
                if tmp != 0 :
                    boxcnt += 1

                    if boxcnt == maxbox:
                        print('Grid overflow, increase number of box count')
                        sys.exit()

                    for k in range(ndim):
                        where[boxcnt,k] = where[int(chaser),k]
                    where[boxcnt,j] = target[j]
                    nxtbox[int(chaser),j] = boxcnt
                    nxtbox[boxcnt,j] = runner
                    runner = boxcnt
                j += 1
            nxtdat[i] = datptr[int(runner)]
            datptr[int(runner)] = i
            
        used = 0
        for i in range(boxcnt+1):
            if datptr[i] != -1:
                used += 1
        print('Created: ', boxcnt)
        print('Used: ', used)
        newdict = {'ndim':int(ndim), 'ires':int(ires), 'tau':int(tau), 'datcnt':self.datcnt, 'boxcnt':int(boxcnt), 'datmax':datmax, 'datmin':datmin, 'boxlen':boxlen, 'datptr':datptr[1:boxcnt+1], 'nxtbox':nxtbox[1:boxcnt+1, :ndim]-1, 'where':where[1:boxcnt+1,:ndim]-1, 'nxtdat':nxtdat[0:self.datcnt], 'data':self.data}
    
    
        return newdict


    def search(iflag,ndim,ires,datmin,
            boxlen,nxtbox,where,datptr,nxtdat,data,delay,oldpnt,newpnt,
            datuse,dismin,dismax,thmax,evolve):
        """
        searches for the most viable point for fet
        return bstpnt, bstdis, thbest
        """
        target = np.zeros(ndim)
        oldcrd = np.zeros(ndim)
        zewcrd = np.zeros(ndim)
    
        oldcrd = data[int(oldpnt)+delay]
        zewcrd = data[int(newpnt)+delay]
        igcrds = np.floor((oldcrd - datmin) / boxlen).astype(int) -1 #added -1
        oldist = np.sqrt(np.sum(np.power(oldcrd - zewcrd,2)))
        irange = int(np.round(dismin/boxlen))
        if irange == 0 :
            irange = 1
    
        thbest = thmax
        bstdis = dismax
        bstpnt = -1
    
        goto30 = 1
        while goto30 == 1:
            goto30 = 0
            for icnt in range(int(((2*irange+1)**ndim))):
                goto140 = 0
                icounter = icnt
                for i in range(ndim):
                    ipower = np.power(2*irange+1,ndim-(i+1))
                    ioff = int(np.floor(icounter/ipower))
                    icounter = icounter - ioff*ipower
                    target[i] = igcrds[i] - irange + ioff
    
                    if target[i] < -1:
                        goto140 = 1
                        break
                    if target[i] > ires:
                        goto140 = 1
                        break
    
                if goto140 ==1:
                    continue

                if irange != 1:
                    iskip = 1
                    for i in range(ndim):
                        if abs(int(np.round(target[i] - igcrds[i]))) == irange:
                            iskip = 0
                    if iskip == 1:
                        continue
                runner = 0
                for i in range(ndim):
                    goto80 = 0
                    goto70 = 1
                    while goto70 == 1:
                        goto70 = 0
                        if where[int(runner),i] == target[i]:
                            goto80 = 1
                            break
                        runner = nxtbox[int(runner),i]
                        if runner !=-1 :
                            goto70 = 1

                    if goto80 == 1:
                        continue
                    goto140 = 1
                    break
                if goto140 == 1:
                    continue
    
                if runner == -1:
                    continue
                runner = datptr[int(runner)]
                if runner == -1:
                    continue
                goto90 = 1
                while goto90 == 1:
                    goto90 = 0
                    while True:
                        if abs(int(np.round(runner-oldpnt))) < evolve:
                            break
                        if abs(int(np.round(runner - datuse))) < (2*evolve):
                            break

                        bstcrd = data[int(runner)+delay]
                        abc1 = oldcrd - bstcrd
                        abc2 = oldcrd - zewcrd
                        tdist = np.sum(abc1*abc1)
                        tdist = np.sqrt(tdist)
                        dot = np.sum(abc1*abc2)
    
                        if tdist < dismin:
                            break
                        if tdist >= bstdis:
                            break
                        if tdist == 0:
                            break
                        goto120 = 0
                        if iflag == 0 :
                            goto120 = 1
                        if goto120 == 0:
                            ctheta = min(abs(dot/(tdist*oldist)),1)
                            theta = 57.3*np.arccos(ctheta)
                            if theta >= thbest:
                                break
                            thbest = theta
                        bstdis = tdist
                        bstpnt = runner
                        break
                    runner = nxtdat[int(runner)]
                    if runner != -1:
                        goto90 = 1
    
            irange += 1
            if irange <= (0.5 + int(np.round((dismax/boxlen)))):
                goto30 = 1
                continue
            return bstpnt, bstdis, thbest

    def fet(db, dt, evolve, dismin, dismax, thmax):
    
        out = []
        
        ndim = db['ndim']
        ires = db['ires']
        tau = db['tau']
        datcnt = db['datcnt']
        datmin = db['datmin']
        boxlen = db['boxlen']
        
        datptr = db['datptr']
        nxtbox = db['nxtbox']
        where = db['where']
        nxtdat = db['nxtdat']
        data = db['data']
    
        delay = np.array([0,tau,(ndim-1)*tau])
        datuse = datcnt - (ndim-1)*tau - evolve
    
        its = 0
        SUM = 0
        savmax = dismax
    
        oldpnt = 0     #1 in matlab original
        newpnt = 0     #1 in matlab original
        goto50 = 1
        while goto50 == 1:
            goto50 = 0
            bstpnt, bstdis, thbest = LYAP.search(0, ndim, ires, datmin, boxlen, nxtbox, where, \
                    datptr, nxtdat, data, delay, oldpnt, newpnt, datuse, dismin, dismax, \
                    thmax, evolve)
            
            while bstpnt == -1 :
                dismax = dismax * 2
                bstpnt, bstdis, thbest = LYAP.search(0, ndim, ires, datmin, boxlen, nxtbox, where, \
                    datptr, nxtdat, data, delay, oldpnt, newpnt, datuse, dismin, dismax, \
                    thmax, evolve)
    
            dismax = savmax 
            newpnt = bstpnt
            disold = bstdis
            iang = -1
    
            goto60 = 1
            while goto60 == 1:
                goto60 = 0
    
                oldpnt += evolve
                newpnt += evolve
    
                if oldpnt >= datuse:
                    print('Lyapunov exponent: ', zlyap)
                    return out, SUM, zlyap
    
                if newpnt >= datuse:
                    oldpnt = oldpnt - evolve
                    goto50 = 1
                    break
    
                p1 = data[int(oldpnt) + delay]
                p2 = data[int(newpnt) + delay]
                disnew = np.sqrt(np.sum(np.power(p2-p1,2)))
    
                its = its + 1
    
                SUM = SUM + np.log(disnew/disold)
                zlyap =  SUM/(its*evolve*dt*np.log(2))    # base 2 Lyapunov exponent 
#               print('***********')
#               print('z_lyap:  ',zlyap)
    
                out = [out, its*evolve, disold, disnew, zlyap, (oldpnt-evolve), (newpnt-evolve)]
    
                #if iang == -1:
                #   fprintf(fileID, out[end,0:3])
                #else:
                #   fprintf(fileID, [out[end,0:3i],iang])
                
                if disnew <= dismax:
                    disold = disnew
                    iang = -1
                    goto60 = 1
                    continue
                
                bstpnt, bstdis, thbest = LYAP.search(1, ndim, ires, datmin, boxlen, nxtbox, where, \
                    datptr, nxtdat, data, delay, oldpnt, newpnt, datuse, dismin, dismax, \
                    thmax, evolve)
                
                if bstpnt != -1: 
                    newpnt = bstpnt
                    disold = bstdis
                    iang = np.floor(thbest)
                    goto60 = 1
                    continue
                else:
                    goto50 = 1
                    break


    def lyap_e(self,tau=10,ndim=3,ires=10,maxbox=6000,
            dt=0.01,evolve=20,dismin=0.001,dismax=0.3,thmax=30):
    
        db = LYAP.basgen(self,tau,ndim,ires,maxbox)
    
        l = LYAP.fet(db,dt,evolve,dismin,dismax,thmax)
    
        return l[-1]  #returns the Lyapunov exponent in base 2
    



