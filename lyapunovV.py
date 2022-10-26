import numpy as np
import sys
from scipy.integrate import odeint
import ipdb


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
                print('count8',flush=True)
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
#       print('Created: ', boxcnt,flush=True)
#       print('Used: ', used,flush=True)
#       print('nxtbox',nxtbox[1:boxcnt+1, :ndim]-1)
#       print('where',where[1:boxcnt+1,:ndim]-1)
#       print('datptr',datptr[1:boxcnt+1])
        newdict = {'ndim':int(ndim), 'ires':int(ires), 'tau':int(tau), 'datcnt':self.datcnt, 'boxcnt':int(boxcnt), 'datmax':datmax, 'datmin':datmin, 'boxlen':boxlen, 'datptr':datptr[1:boxcnt+1], 'nxtbox':nxtbox[1:boxcnt+1, :ndim]-1, 'where':where[1:boxcnt+1,:ndim]-1, 'nxtdat':nxtdat[0:self.datcnt], 'data':self.data}
    
    
        return newdict
    

    def searchV(debug,iflag,ndim,ires,datmin,
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
#                   if debug and( icnt == 2 ):
#                     ipdb.set_trace()
                    ctheta = min(np.abs(dot/(tdist*oldist)),1)
                    theta = 57.3*np.arccos(ctheta)
                    print('dot0',dot,'oldist',oldist,'oldpnt',oldpnt,'newpnt',newpnt,'irange',irange,'runner',run,'icnt',icnt)
                    if theta >= thbest:
                        select = False
                    if select:
                        print('dot',dot,'oldist',oldist,'oldpnt',oldpnt,'newpnt',newpnt,'irange',irange,'runner',run,'icnt',icnt)
                        print('dot1',dot,'bstcrd',bstcrd,'run',run,'delay',delay)
                        print('theta',theta,'thbest',thbest,'irange',irange, 'icnt',icnt)
                        thbest = theta
                if select:
                    bstdis = tdist
                    bstpnt = run
                return thbest, bstdis, bstpnt

        SMAX=int(np.round((dismax/boxlen)))+1
#       ipdb.set_trace()

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

        for irange in range(Irange,SMAX):
#       for irange in range(1,2):
            #       print('irange',irange)
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
#           np.savetxt('target.txt',np.reshape(target,(DMAX,ndim)).astype(int), fmt='%s', delimiter=",")
#           np.savetxt('where.txt',where, fmt='%s', delimiter=",")
#           np.savetxt('nxtbox.txt',nxtbox, fmt='%s', delimiter=",")
            for i in range(ndim):
                BOOLtarget[target[i,:] <      - 1] = False # all dimensions NDIM should be nan when one entry is nan
                BOOLtarget[(target[i,:]) > ires - 2] = False
#           if debug and (irange==6):
#             ipdb.set_trace()

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
#                   BOOLtarget[runner2D[jj,:] == -1] = False
                    loc[np.logical_and(loc_bool,  runnerFLAG )] = jj
                    runnerFLAGgood[i,:] = np.logical_or(runnerFLAGgood[i,:],np.logical_and(loc_bool,BOOLtarget))
                    BOOLtarget[np.logical_and(runner2D[jj,:] == -1,np.logical_not(runnerFLAGgood[i,:]))] = False
                    runnerFLAG = np.logical_and(np.logical_not(loc_bool) , runnerFLAG)
                runnerFLAG=np.zeros(DMAX,dtype=bool)
                runner1D=runner2D[loc,Identity]
#               runner1D[runnerFLAGgood]=runner2D[loc[runnerFLAGgood],Identity[runnerFLAGgood]]


            runner1D[np.logical_not(runnerFLAGgood[-1,:])] = -1
            runner1D[runner1D != -1] = datptr[runner1D[runner1D != -1]]

            kkk=0
            for iv,v in enumerate(runner1D):
                    if v!=-1:
                        print('v ',v,iv,irange,len(runner1D),flush=True)
                        kkk+=1
#           if debug and (irange == 6):
#              ipdb.set_trace()
#           print('runner2',runner1D,len(runner1D),kkk)

#           if irange == 8 :
#               print(irange,SMAX,flush=True)
#               sys.exit()

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

            print("runner1DP",runner1D)

            for icnt in range(DMAX):
#               if (( icnt == 740 ) and debug) and (irange == 6):
#                  ipdb.set_trace()
                runner2D=np.zeros(DATCNT,dtype=int)
                loc=np.zeros(DMAX,dtype=int)
                jj=0
                runner2D[jj]=runner1D[icnt]
                loc_bool[icnt]=np.logical_and(BOOLtarget[icnt],loc_bool[icnt])
                #update bstdis bstpnt theta at jj = 0
                thbest, bstdis, bstpnt = thbest1, bstdis1, bstpnt1
                print('Prunner2Dup',runner2D[jj],'jj',jj,'icnt',icnt)
#               if debug and (icnt==2):
#                   ipdb.set_trace()
                thbest1, bstdis1, bstpnt1 = my_func(runner2D[jj],loc_bool[icnt],thbest, bstdis, bstpnt)

####     start with jj > 0
                for jj in range(1,DATCNT):
                    not_loc_bool[icnt] = np.logical_not(loc_bool[icnt]) # elements still possible to update
                    print('Prunner2D',runner2D[jj-1],'Prunner2Dup',nxtdat[runner2D[jj-1]],'jj',jj,'icnt',icnt)
                    runner2D[jj]=nxtdat[runner2D[jj-1]]
                    if loc_bool[icnt]:
                        runner2D[jj]=-1
                    if runner2D[jj] == -1:
                        break
                    thbest, bstdis, bstpnt = thbest1, bstdis1, bstpnt1
#                   if debug and (icnt == 2):
#                     ipdb.set_trace()
                    thbest1, bstdis1, bstpnt1 = my_func(runner2D[jj],loc_bool[icnt],thbest, bstdis, bstpnt)
                    if (thbest1>5.871) and (thbest1 < 5.873):
                         print('thetacheck',irange,runner2D[jj],jj,icnt)

                    loc_bool[icnt]=np.logical_or(runner2D[jj] == -1,loc_bool[icnt])
            print(thbest1, bstdis1, bstpnt1,'irange',irange)
            print('********')
#           if irange == 8:
#               sys.exit()
#       print(runner1D[runner1D != -1],irange,len(runner1D))
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
            print('count4',flush=True)
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

                print('igcrds',igcrds,target,irange,abs(int(np.round(target[i] - igcrds[i]))))
                if irange != 1:
                    iskip = 1
                    for i in range(ndim):
                        if abs(int(np.round(target[i] - igcrds[i]))) == irange:
                            iskip = 0
                            BOOL[icnt] = 0
                    if iskip == 1:
                        continue
                print('BOOL2', BOOL)
                runner = 0
                for i in range(ndim):
                    goto80 = 0
                    goto70 = 1
                    jj=0
                    while goto70 == 1:
                        print('count5',flush=True)
                        goto70 = 0
                        if where[int(runner),i] == target[i]:
                            runner1D[icnt] = runner
                            if irange ==8:
                                print('runnertrue', 'ndim', i, 'icnt',icnt,'runner',runner,'jj',jj,
                                        'target',target,'where',where[int(runner),:])
                                if icnt == 289:
                                    print('target289',target)
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
    
#               runner1D[icnt] = runner
#               print('runner1D',runner1D,flush=True)
#               print('runner0D',runner,int(((2*irange+1)**ndim)),flush=True)
                print('datptr',datptr,flush=True)
                print('target2D',target2D,icnt,irange,int(((2*irange+1)**ndim)),flush=True)
                if int(runner) == -1:
                    continue
#               print('runner0D',runner,int(((2*irange+1)**ndim)),flush=True)
                runner = datptr[int(runner)]
#               print('runnerlast',runner,icnt,irange,int(((2*irange+1)**ndim)),flush=True)
                if int(runner) == -1:
                    continue
                print('runnerlast',runner,icnt,irange,int(((2*irange+1)**ndim)),flush=True)
                goto90 = 1
                jjj=0
                while goto90 == 1:
                    print('count6',flush=True)
                    goto90 = 0
                    iii=0
                    while True:
                        print('count7',flush=True)
                        if abs(int(np.round(runner-oldpnt))) < evolve:
                            if icnt ==24:
                               print('Stop0', runner, nxtdat[int(runner)],bstpnt)
                            break
                        if abs(int(np.round(runner - datuse))) < (2*evolve):
                            if icnt ==24:
                               print('Stop1', runner, nxtdat[int(runner)],bstpnt)
                            break

                        bstcrd = data[int(runner)+delay]
                        abc1 = oldcrd - bstcrd
                        abc2 = oldcrd - zewcrd
                        tdist = np.sum(abc1*abc1)
                        tdist = np.sqrt(tdist)
                        dot = np.sum(abc1*abc2)
                        print('tdist',tdist,'jj',jjj,'irange',irange,'icnt',icnt, 'runner', runner)
                        if tdist < dismin:
                            if icnt ==24:
                               print('Stop2', runner, nxtdat[int(runner)],bstpnt)
                            break
                        if tdist >= bstdis:
                            if icnt ==24:
                               print('Stop3', runner, nxtdat[int(runner)],bstpnt)
                            break
                        if tdist == 0:
                            if icnt ==24:
                               print('Stop4', runner, nxtdat[int(runner)],bstpnt)
                            break
                        goto120 = 0
                        if iflag == 0 :
                            goto120 = 1
                        if goto120 == 0:
                            ctheta = min(abs(dot/(tdist*oldist)),1)
                            theta = 57.3*np.arccos(ctheta)
                            print('dot0',dot,'oldist',oldist,'oldpnt',oldpnt,'newpnt',newpnt,'irange',irange,'runner',runner,'icnt',icnt)
                            if theta < thbest:
                                print('theta',theta,'thbest',thbest)
                            if theta >= thbest:
                                break
                            if (theta>13.804) and (theta < 13.805):
                                print('thetacheck',irange,runner)
                            print('dot',dot,'oldist',oldist,'oldpnt',oldpnt,'newpnt',newpnt,'irange',irange,'runner',runner,'icnt',icnt)
#                           print('dot1',dot,'abc1',abc1,'abc2',abc2,'oldcrd',oldcrd,'zewcrd',zewcrd)
#                           print('dot1',dot,'abc1',abc1,'abc2',abc2,'oldcrd',oldcrd,'zewcrd',zewcrd,'bstcrd',bstcrd)
                            print('dot1',dot,'bstcrd',bstcrd,'runner',runner,'delay',delay)
                            thbest = theta
                        bstdis = tdist
                        bstpnt = runner
                        break
                    print('00',bstdis,bstpnt,'jj',jjj)
                    if icnt ==24:
                        print('runner24_0', runner, nxtdat[int(runner)],bstpnt)
                    print('Prunner2D',runner,'Prunner2Dup',nxtdat[runner],'jj',jjj,'icnt',icnt)
                    runner = nxtdat[int(runner)]
                    jjj +=1
                    if runner != -1:
                        goto90 = 1
            
            
            #########
                print('bstpnt',bstpnt,'bstdis', bstdis,'thbest', thbest,'irange',irange,'icnt',icnt,'runner',runner)
#           if irange == 1:
#               print('bstpnt',bstpnt,'bstdis', bstdis,'thbest', thbest)
#               sys.exit() 
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
        mycount = 0
        debug = False
        while goto50 == 1:
            goto50 = 0
            print('count1',flush=True)
#P          np.save('datmin.npy', datmin)
#P          np.save('boxlen.npy', boxlen)
#P          np.save('nxtbox.npy', nxtbox)
#P          np.save('where.npy', where)
#P          np.save('datptr.npy', datptr)
#P          np.save('nxtdat.npy', nxtdat)
#P          np.savez_compressed('data.npz', data=data.data, mask=data.mask)
#P          np.save('delay.npy', delay)
#P          np.save('oldpnt.npy', oldpnt)
#P          np.save('newpnt.npy', newpnt)
#P          np.save('datuse.npy', datuse)
#P          np.save('dismin.npy', dismin)
#P          np.save('thmax.npy', thmax)
#P          np.save('evolve.npy', evolve)
#P          np.save('dismax.npy', dismax)
            if mycount == 6:
               debug = True
#           ipdb.set_trace()
            bstpnt, bstdis, thbest = LYAP.searchV(debug,0, ndim, ires, datmin, boxlen, nxtbox, where, \
                    datptr, nxtdat, data, delay, oldpnt, newpnt, datuse, dismin, dismax, \
                    thmax, evolve)
            print('mycount',mycount, bstpnt, bstdis, thbest,0)
#           ipdb.set_trace()
            mycount += 1
            while bstpnt == -1 :
                print('count2',flush=True)
                dismax = dismax * 2
                if mycount ==6 :
                    debug=True
#               if debug:
#               ipdb.set_trace()
                bstpnt, bstdis, thbest = LYAP.searchV(debug,0, ndim, ires, datmin, boxlen, nxtbox, where, \
                    datptr, nxtdat, data, delay, oldpnt, newpnt, datuse, dismin, dismax, \
                    thmax, evolve)
                print('mycount',mycount, bstpnt, bstdis, thbest,1)
#               if debug:
#               ipdb.set_trace()
                mycount +=1 
            print('fet','bstpnt',bstpnt,'bstdis',bstdis,'thbest',thbest)
            dismax = savmax 
            newpnt = bstpnt
            disold = bstdis
            iang = -1
    
            goto60 = 1
            while goto60 == 1:
                print('count3',flush=True)
                goto60 = 0
    
                oldpnt += evolve
                newpnt += evolve
    
                if oldpnt >= datuse:
                    print('Lyapunov exponent: ', zlyap)
                    return out, SUM, zlyap
    
                if newpnt >= datuse:
                    oldpnt = oldpnt - evolve
                    goto50 = 1
                    print('goto50=1_0')
                    break
    
                p1 = data[int(oldpnt) + delay]
                p2 = data[int(newpnt) + delay]
                disnew = np.sqrt(np.sum(np.power(p2-p1,2)))
    
                its = its + 1
    
                SUM = SUM + np.log(disnew/disold)
                zlyap =  SUM/(its*evolve*dt*np.log(2))    # base 2 Lyapunov exponent 
                print('***********',flush=True)
                print('z_lyap:  ',zlyap,flush=True)
    
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
                
                np.save('oldpnt.npy', oldpnt)
                np.save('newpnt.npy', newpnt)
                print('input', 'oldpnt', oldpnt, 'newpnt',newpnt)

                if mycount ==6 :
                     debug=True
#               ipdb.set_trace()
                bstpnt, bstdis, thbest = LYAP.searchV(debug,1, ndim, ires, datmin, boxlen, nxtbox, where, \
                    datptr, nxtdat, data, delay, oldpnt, newpnt, datuse, dismin, dismax, \
                    thmax, evolve)
#               bstpnt, bstdis, thbest = LYAP.search(1, ndim, ires, datmin, boxlen, nxtbox, where, \
#                   datptr, nxtdat, data, delay, oldpnt, newpnt, datuse, dismin, dismax, \
#                   thmax, evolve)
                print('mycount',mycount, bstpnt, bstdis, thbest,2)
#               ipdb.set_trace()
                mycount +=1
                
                
                if bstpnt != -1: 
                    newpnt = bstpnt
                    disold = bstdis
                    iang = np.floor(thbest)
                    goto60 = 1
                    continue
                else:
                    print('goto50=1_1')
                    goto50 = 1
                    break


    def lyap_e(self,tau=10,ndim=3,ires=10,maxbox=6000,
            dt=0.01,evolve=20,dismin=0.001,dismax=0.3,thmax=30):
    
        db = LYAP.basgen(self,tau,ndim,ires,maxbox)
    
        l = LYAP.fet(db,dt,evolve,dismin,dismax,thmax)
    
        return l[-1]  #returns the Lyapunov exponent in base 2
    



