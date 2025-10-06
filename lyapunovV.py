import numpy as np
import sys
import warnings

# Make numpy warnings raise exceptions
np.seterr(all="raise")


class LYAP(object):
    '''Compute the largest Lyapunov exponent of a timeseries'''
    def __init__(self, data):
        self.data = data
        self.datcnt = len(data)
        if not np.all(np.isfinite(self.data)):
            print("[DEBUG] Warning: input data contains NaN or inf values!")

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

        datmin = np.nanmin(self.data)
        datmax = np.nanmax(self.data)

        datmin = datmin - 0.01*(datmax - datmin)
        datmax = datmax + 0.01*(datmax - datmin)
        boxlen = (datmax - datmin)/ires
        if not np.isfinite(boxlen):
            print(f"[DEBUG] boxlen is invalid: {boxlen}")
    
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
        newdict = {'ndim':int(ndim), 'ires':int(ires), 'tau':int(tau), 'datcnt':self.datcnt, 'boxcnt':int(boxcnt), 'datmax':datmax, 'datmin':datmin, 'boxlen':boxlen, 'datptr':datptr[1:boxcnt+1], 'nxtbox':nxtbox[1:boxcnt+1, :ndim]-1, 'where':where[1:boxcnt+1,:ndim]-1, 'nxtdat':nxtdat[0:self.datcnt], 'data':self.data}
    
    
        return newdict
    

    def searchV(iflag,ndim,ires,datmin,
            boxlen,nxtbox,where,datptr,nxtdat,data,delay,oldpnt,newpnt,
            datuse,dismin,dismax,thmax,evolve):
        def my_func(run,bool_stop,thbest, bstdis, bstpnt):
                select = True
                bstcrd = data[run+delay]
                abc1 = oldcrd - bstcrd
                abc2 = oldcrd - zewcrd
                tdist = np.sum(abc1*abc1,axis=0)
                tdist = np.sqrt(tdist)
                dot = np.sum(abc1*abc2,axis=0)
                if np.abs(np.round(run-oldpnt).astype(int)) < evolve:
                    select = False
                if np.abs(np.round(run-datuse).astype(int)) < (2*evolve):
                    select = False
                if tdist < dismin:
                    select = False
                if tdist >= bstdis:
                    select = False
                if tdist == 0:
                    select = False
                if bool_stop:
                    select = False

                if iflag != 0:
                    ctheta = min(np.abs(dot/(tdist*oldist)),1)
                    theta = 57.3*np.arccos(ctheta)
                    if theta >= thbest:
                        select = False
                    if select:
                        thbest = theta
                if select:
                    bstdis = tdist
                    bstpnt = run
                return thbest, bstdis, bstpnt

        SMAX=int(np.round((dismax/boxlen)))+1

        oldcrd = data[int(oldpnt)+delay]
        zewcrd = data[int(newpnt)+delay]


        igcrds = np.floor((oldcrd - datmin) / boxlen).astype(int) -1 #added -1
        oldist = np.sqrt(np.sum(np.power(oldcrd - zewcrd,2)))
        Irange = int(np.round(dismin/boxlen))

        pos = 0
        count=0

        if Irange == 0 :
            Irange = 1

        thbest1 = thmax
        bstdis1 = dismax
        bstpnt1 = -1

        indim=range(ndim)
        ipower=np.zeros(ndim)

        #do at least one cycle
        if SMAX <= Irange:
          SMAX = Irange+1
        if SMAX > 250:
            sys.exit()
#       SMAX = min(SMAX,150)
        for irange in range(Irange,SMAX):
            DMAX=int((2*irange+1)**ndim)
            Identity=np.arange(DMAX,dtype=int)
            target = np.zeros((ndim,DMAX),dtype=int)
            BOOLtarget = np.logical_not( np.zeros(DMAX,dtype=bool) )
            BOOLtarget_helper = np.logical_not( np.zeros(DMAX,dtype=bool) )
            ioff=np.zeros((ndim,DMAX),dtype=int)
            icounter=np.zeros((ndim,DMAX),dtype=int)

            for i in range(ndim):
                ipower[i] = int(np.power(2*irange+1,ndim-(i+1)))
                icounter[i,:]=np.arange(DMAX,dtype=int)

            ioff[0,:] = np.floor(icounter[0,:]/ipower[0]).astype(int)
            icounter[0,:]= icounter[0,:] - ioff[0,:]*ipower[0]
            target[0,:]=igcrds[0] - irange + ioff[0,:]


            for i in range(1,ndim):

                ioff[i,:] = np.floor(icounter[i-1,:]/ipower[i]).astype(int)
                icounter[i,:]= icounter[i-1,:] - ioff[i,:]*ipower[i]
                target[i,:]=igcrds[i] - irange + ioff[i,:]
            for i in range(ndim):
                BOOLtarget[target[i,:] <      - 1] = False # all dimensions NDIM should be nan when one entry is nan
                BOOLtarget[(target[i,:]) > ires - 2] = False

            if irange != 1:
                for i in range(ndim):
                    BOOLtarget_helper[(np.abs((np.round(target[i,:] - igcrds[i])).astype(int))) == irange] = False
                BOOLtarget = np.logical_and(BOOLtarget,np.logical_not(BOOLtarget_helper))

            AAAA=len(where)
            runner1D=np.zeros(DMAX,dtype=int)
            runnerFLAG=np.logical_not(np.zeros(DMAX,dtype=bool)) # initialized to True
            runnerFLAGgood=np.zeros((ndim,DMAX),dtype=bool) # initialized to True
            runnerFLAG = np.logical_and(BOOLtarget , runnerFLAG) # discard wrong icnt
            runner2D=np.zeros((AAAA,DMAX),dtype=int)

            for i in range(ndim):
                loc=np.zeros(DMAX,dtype=int)
                if i != 0:
                    runnerFLAG[runnerFLAGgood[i-1,:]] = True

                for jj in range(0,AAAA):
                    if jj ==0:
                        runner2D[jj,:]=runner1D
                    else:
                        runner2D[jj,:]=nxtbox[runner2D[jj-1,:],i]# -- > check

                    if np.all(runner2D[jj,:] == -1):
                        break
                    loc_bool=np.logical_and(where[runner2D[jj,:],i] == target[i,:] , runner2D[jj,:] != -1)
                    loc[np.logical_and(loc_bool,  runnerFLAG )] = jj
                    runnerFLAGgood[i,:] = np.logical_or(runnerFLAGgood[i,:],np.logical_and(loc_bool,BOOLtarget))
                    BOOLtarget[np.logical_and(runner2D[jj,:] == -1,np.logical_not(runnerFLAGgood[i,:]))] = False
                    runnerFLAG = np.logical_and(np.logical_not(loc_bool) , runnerFLAG)
                runnerFLAG=np.zeros(DMAX,dtype=bool)
                runner1D=runner2D[loc,Identity]


            runner1D[np.logical_not(runnerFLAGgood[-1,:])] = -1
            runner1D[runner1D != -1] = datptr[runner1D[runner1D != -1]]


            if runner1D.all() == -1:
                continue

#restrict to no -1 values
            restrict_bool = runner1D != -1
            runner1D = runner1D[restrict_bool]
            BOOLtarget = BOOLtarget[restrict_bool]
            DMAX = len(runner1D)


            DATCNT = len(nxtdat)
            loc_bool=np.full(DMAX,False,dtype=bool)
            not_loc_bool=np.full(DMAX,False,dtype=bool)
            BOOLtarget = np.logical_and(runner1D != -1 , BOOLtarget) # discard wrong icnt


            for icnt in range(DMAX):
                runner2D=np.zeros(DATCNT,dtype=int)
                loc=np.zeros(DMAX,dtype=int)
                jj=0
                runner2D[jj]=runner1D[icnt]
                loc_bool[icnt]=np.logical_and(BOOLtarget[icnt],loc_bool[icnt])
                #update bstdis bstpnt theta at jj = 0
                thbest, bstdis, bstpnt = thbest1, bstdis1, bstpnt1
                thbest1, bstdis1, bstpnt1 = my_func(runner2D[jj],loc_bool[icnt],thbest, bstdis, bstpnt)

####     start with jj > 0
                for jj in range(1,DATCNT):
                    not_loc_bool[icnt] = np.logical_not(loc_bool[icnt]) # elements still possible to update
                    runner2D[jj]=nxtdat[runner2D[jj-1]]
                    if loc_bool[icnt]:
                        runner2D[jj]=-1
                    if runner2D[jj] == -1:
                        break
                    thbest, bstdis, bstpnt = thbest1, bstdis1, bstpnt1
                    thbest1, bstdis1, bstpnt1 = my_func(runner2D[jj],loc_bool[icnt],thbest, bstdis, bstpnt)

                    loc_bool[icnt]=np.logical_or(runner2D[jj] == -1,loc_bool[icnt])
        return  bstpnt1, bstdis1, thbest1




    def search(iflag,ndim,ires,datmin,
            boxlen,nxtbox,where,datptr,nxtdat,data,delay,oldpnt,newpnt,
            datuse,dismin,dismax,thmax,evolve):
        """
        searches for the most viable point for fet
        return bstpnt, bstdis, thbest
        """
        target = np.zeros(ndim,dtype=int)
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
            DMAX = int(((2*irange+1)**ndim))
            runner1D = np.full(DMAX, 0,dtype=int)
            target2D = np.full((ndim,DMAX), 0,dtype=int)
            BOOL = np.full(DMAX,-1,dtype=int)
            for icnt in range(int(((2*irange+1)**ndim))):
                goto140 = 0
                icounter = icnt
                for i in range(ndim):
                    ipower = int(np.power(2*irange+1,ndim-(i+1)))
                    ioff = int(np.floor(icounter/ipower))
                    icounter = icounter - ioff*ipower
                    target[i] = igcrds[i] - irange + ioff
                    target2D[i,icnt] = target[i]
                    if target[i] < -1:
                        goto140 = 1
                        break
                    if target[i] > ires-2:
                        goto140 = 1
                        break
    
                if goto140 ==1:
                    continue

                if irange != 1:
                    iskip = 1
                    for i in range(ndim):
                        if abs(int(np.round(target[i] - igcrds[i]))) == irange:
                            iskip = 0
                            BOOL[icnt] = 0
                    if iskip == 1:
                        continue
                runner = 0
                for i in range(ndim):
                    goto80 = 0
                    goto70 = 1
                    jj=0
                    while goto70 == 1:
                        goto70 = 0
                        if where[int(runner),i] == target[i]:
                            runner1D[icnt] = runner
                            goto80 = 1
                            break
                        runner = nxtbox[int(runner),i]
                        if int(runner) !=-1 :
                            goto70 = 1
                        jj+=1

                    if goto80 == 1:
                        continue
                    goto140 = 1
                    break
                if goto140 == 1:
                    continue
    
                if int(runner) == -1:
                    continue
                runner = datptr[int(runner)]
                if int(runner) == -1:
                    continue
                goto90 = 1
                jjj=0
                while goto90 == 1:
                    goto90 = 0
                    iii=0
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
                    jjj +=1
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
            bstpnt, bstdis, thbest = LYAP.searchV(0, ndim, ires, datmin, boxlen, nxtbox, where, \
                    datptr, nxtdat, data, delay, oldpnt, newpnt, datuse, dismin, dismax, \
                    thmax, evolve)
            while bstpnt == -1 :
                dismax = dismax * 2
                bstpnt, bstdis, thbest = LYAP.searchV(0, ndim, ires, datmin, boxlen, nxtbox, where, \
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

                #if p1 or p2 are NaN, skip: goto50 = 1 break
                if np.any(np.isnan(p1)) or np.any(np.isnan(p2)):
                    oldpnt = oldpnt - evolve
                    goto50 = 1
                    break

                disnew = np.sqrt(np.sum(np.power(p2-p1,2)))
                if not np.isfinite(disnew):
                    print(f"[DEBUG] disnew became NaN/inf at its={its}, oldpnt={oldpnt}, newpnt={newpnt}")

    
                its = its + 1

                try:
                    SUM = SUM + np.log(disnew/disold)
                    if not np.isfinite(SUM):
                        print(f"[DEBUG] SUM became NaN/inf at its={its}, disnew={disnew}, disold={disold}")
                    zlyap =  SUM/(its*evolve*dt*np.log(2))    # base 2 Lyapunov exponent 
                    if not np.isfinite(zlyap):
                        print(f"[DEBUG] zlyap became NaN/inf at its={its}, SUM={SUM}, disnew={disnew}, disold={disold}")
                except (FloatingPointError, OverflowError, ZeroDivisionError, ValueError):
                    print('Numerical error encountered, returning NaN')
                    return [], 0, np.nan
    
                out = [out, its*evolve, disold, disnew, zlyap, (oldpnt-evolve), (newpnt-evolve)]
    
                
                if disnew <= dismax:
                    disold = disnew
                    iang = -1
                    goto60 = 1
                    continue
                

                bstpnt, bstdis, thbest = LYAP.searchV(1, ndim, ires, datmin, boxlen, nxtbox, where, \
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
    

    #new implementation of Palladin et al 1995 
    def search_delta(self, embedded, oldpnt, delta0, Delta):
        x_old = embedded[oldpnt]
        diffs = embedded - x_old
        dists = np.linalg.norm(diffs, axis=1)

        # exclude self
        dists[oldpnt] = -np.inf

        # mask invalid neighbors
        mask = (dists >= delta0) & (dists <= Delta)
        dists[~mask] = -np.inf   
        
        if np.all(~mask):
            return None, None, None

        best_pnt = np.argmax(dists)
        return best_pnt, dists[best_pnt], best_pnt

    def fet_temporal(self, db, dt, delta0=1e-5, Delta=0.3):
        """
        Implements Paladin et al., 1995 algorithm with vectorized neighbor search
        """
        out = []
        tau = db['tau']
        ndim = db['ndim']
        datcnt = db['datcnt']
        datuse = datcnt - (ndim - 1) * tau

        #Step 1: precompute embedded phase space
        embedded = np.array([
            db['data'][i + np.arange(ndim) * tau]
            for i in range(datuse)
        ])

        evolve = 1
        oldpnt = 0
        SUM = 0
        count = 0

        while oldpnt < datuse:
            #Step 2: vectorized search
            newpnt, dist, new_time = self.search_delta(embedded, oldpnt, delta0, Delta)

            if newpnt is None:
                oldpnt += evolve
                continue

            delta_time = new_time - oldpnt
            SUM += delta_time
            count += 1

            out.append([oldpnt, newpnt, count, delta_time])
            oldpnt = newpnt  # rescale trajectory

        # compute mean Lyapunov exponent
        if count > 0:
            delta_time_bar = SUM / count
            lyap_exp = (1 / delta_time_bar) * np.log(Delta / delta0) / dt
        else:
            lyap_exp = np.nan

        return out, lyap_exp

    def lyap_e_paladin(self, tau=10, ndim=3, ires=10, maxbox=6000,
                       dt=0.01, delta0=1e-5, Delta=0.3):
        """
        Wrapper function to compute Lyapunov exponent using Paladin et al., 1995 method
        """
        db = self.basgen(tau, ndim, ires, maxbox)
        out, lyap_exp = self.fet_temporal(db, dt, delta0, Delta)
        return out, lyap_exp