# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 15:28:52 2024

@author: SSKJCFD004
"""

import matplotlib.pyplot as plt
import numpy as np
import os 
import copy

def cal_elastance(tt,Tvcp,Tacp,Tarp,tac,tar,tao,T,c_type:str):
    elas = 0.1
    t = tt%T
    if c_type == 'ra' or c_type == 'la':
        if t<=tar+Tarp-T:
            elas = 0.5* (1+ np.cos( np.pi*(t+T-tar)/Tarp ))
        elif t>tar+Tarp-T and t<=tac:
            elas = 0
        elif t>tac and t<=tac+Tacp:
            elas = 0.5* (1- np.cos( np.pi*(t-tac)/Tacp ))
        else:
            elas = 0.5* (1+ np.cos( np.pi*(t-tar)/Tarp))
    elif (c_type == 'rv' or c_type == 'lv'):
        if t<=Tvcp:
            elas = 0.5*(1-np.cos( np.pi*t/Tvcp )) + np.exp(-(T-Tvcp)/tao) * np.cos(0.5*np.pi*t/Tvcp)
        else:
            elas = np.exp(-(t-Tvcp)/tao)
    else:
        elas = 0.001
        
    return elas

    
    C_matj2 = np.mat(np.zeros((L.shape[0],1)))
    C_matj3 = np.mat(np.zeros((L.shape[0],1)))
    
    for i in range(0,L.shape[0]):
        if i in valve:
            C_matj2[i,0] = (vj[i,0]-vm[i,0])*cham_e[j,valve.index(i)]
        else:
            C_matj2[i,0] = (vj[i,0]-vm[i,0])/C[i]
        if cc_in[i+7]<0:
            C_matj3[i,0] = q6j-qj[cc_out[i+7],0]
        else:
            C_matj3[i,0] = qj[cc_in[i+7],0]-qj[cc_out[i+7],0]
            
    return C_matj2,C_matj3

if __name__ =="__main__":
    os.chdir(r'F:\linzhihong_doc\erke\python')
    if(1):
        T = 0.8
        n_cycle = 25
        dt=0.001
        t = np.arange(0,n_cycle*T,dt)
        # simulating LPN coupled with CFD, CFD is simplified as a compartment with r6,c6,l6
        ##0-6: ra, rv, pulmony artery, la,lv, sa,sv
        ##p6,q6:aorta
        P = np.mat(np.zeros((len(t),7))) #result storage array for LPN pressure mmHg
        Q = np.mat(np.zeros((len(t),7))) #result storage array for LPN flow ml/s
        V = np.mat(np.zeros((len(t),7))) ##volume of chambers
        p6 = np.zeros(len(t))  #result storage array for CFD
        q6 = np.zeros(len(t)) #result storage array for CFD
        
        #parameter are set randomly
        L = np.zeros(7) #parameter of indutance
        R = np.zeros(7) # parameter of resistance
        C = np.zeros(7) #parameter of compliance of elastance
        B = np.zeros(7) #parameter of flow coefficient of valve
        alpha = np.zeros(7) #parameter of viscoresistance
        r6 = 0.143 #parameter of RLC for CFD
        c6 = 1.13
        l6 = 0.015
        
        if (1):
            L[0],B[0],R[0],alpha[0] = 5e-4, 1e-5,1e-3,5e-4
            L[1],B[1],R[1],alpha[1] = 5e-4, 0.001,1.5e-3,5e-4
            L[2],R[2],C[2]  =         1e-3, 0.0761,9
            L[3],B[3],R[3],alpha[3] = 5e-4, 1e-5,1e-3,5e-4
            L[4],B[4],R[4],alpha[4] = 5e-4, 1.5e-5,1.5e-3,5e-4
            L[5],R[5],C[5]  =         5e-3, 0.779,0.82
            L[6],R[6],C[6]  =         2e-3, 0.21,93.2
            l6,r6,c6 =                0.015,0.143,1.13
        
        valve = [0,1,3,4] #Q node for valve
        visco = [0,1,3,4] #P node for valve
    
    if(1):
        HR = 60/T
        tao = (30.2*np.exp(-HR/81.2)+31.4)/1000
        RR = 60/HR
        QT = -0.33*RR**2+0.69*RR+0.029
        Tvcp = 0.714*QT
        Tacp = 0.6*Tvcp
        Tarp = Tacp
        tac = T-Tacp-0.05
        tar = tac+Tacp
        cham_e = np.zeros((len(t),4))
        for i in range(0,t.shape[0]):
            cham_e[i,2] = 0.25*cal_elastance(t[i], Tvcp, Tacp, Tarp, tac, tar, tao, T, 'la')+0.25
            cham_e[i,3] = 2.87*cal_elastance(t[i], Tvcp, Tacp, Tarp, tac, tar, tao, T, 'lv')+0.056
            cham_e[i,0] = 0.13*cal_elastance(t[i], Tvcp, Tacp, Tarp, tac, tar, tao, T, 'ra')+0.13
            cham_e[i,1] = 0.48*cal_elastance(t[i], Tvcp, Tacp, Tarp, tac, tar, tao, T, 'rv')+0.05
        
        # plt.figure()
        # plt.plot(t,cham_e[:,2],'k')
        # plt.plot(t,cham_e[:,0],'r')
        # plt.plot(t,lea,'b')
        # plt.plot(t,lev,'y')
        # plt.show()
        
        
        pn_in =  [0,1,2,3,4,  5,6]
        pn_out = [1,2,3,4,-10,6,0]
        qn_in =  [6,0,1,2,3,-10,5]
        qn_out = [0,1,2,3,4,5  ,6]
        cc_in = [0,1,2,3,4  ,5,6,6,0,1,2,3,-10,5] ## pin node for Qnode:0-6; Qin for Pnode: 7-13; node=-10 mean value calculate from CFD
        cc_out =[1,2,3,4,-10,6,0,0,1,2,3,4,5  ,6] ## pout node for Qnode:0-6;Qout for Pnode:7-13; node=-10 mean value calculate from CFD
    
    #initialization
    if (1): ##note!!!: initial volume must be set for ventricle or atrium, otherwise no flow would be generated 
        Q[0,:] = 0  ## initial flow  0ml/s
        V[0,0],V[0,1],V[0,2],V[0,3] = 38.46,75.82,192.932,68.42
        V[0,4],V[0,5],V[0,6]        = 117.4,93.32,41.98
        q6[0] = 0
        p6[0] = 19.64/c6
        for i in range(0,len(L)):
            if i in valve:
                P[0,i] = V[0,i]*cham_e[0,valve.index(i)]
            else:
                P[0,i] = V[0,i]/C[i]
        eps,ct,flag = 1,0,1
    
    buffer = np.zeros(L.shape[0])
    for j in range(1,len(t)):
        eps,ct,flag = 1,0,1
        q6m = q6[j-1]
        p6m = p6[j-1]
        pjm = P[j-1,:].T
        qjm = Q[j-1,:].T
        vm = V[j-1,:].T
        pj = copy.deepcopy(pjm)
        qj = copy.deepcopy(qjm)
        vj = copy.deepcopy(vm)
        p6j = copy.deepcopy(p6m)
        q6j = copy.deepcopy(q6m)
        pj2 = copy.deepcopy(pjm)
        qj2 = copy.deepcopy(qjm)
        vj2 = copy.deepcopy(vm)
        
        for i in range(0,L.shape[0]):
            if cc_out[i]<0:
                dp = pjm[cc_in[i],0]-p6m
                buffer[i] = (dp-R[i]*qjm[i,0])/L[i]
            else:
                dp = pjm[cc_in[i],0]-pjm[cc_out[i],0]
                buffer[i] = (dp-R[i]*qjm[i,0])/L[i]
            if i in valve:
                if qjm[i,0]>=0:
                    qj[i,0] = qjm[i,0] + buffer[i]*dt
                else:
                    qj[i,0] = (dp*dt+L[i]*qjm[i,0]) / (L[i]+(R[i]+100000)*dt)
            else:
                qj[i,0] = qjm[i,0] + buffer[i]*dt
        for i in range(0,L.shape[0]):
            if cc_in[i+7]<0:
                dq = (q6m-qj[cc_out[i+7],0] + q6m-qjm[cc_out[i+7],0])/2
            else:
                dq = (qj[cc_in[i+7],0]-qj[cc_out[i+7],0] + qjm[cc_in[i+7],0]-qjm[cc_out[i+7],0])/2
            vj[i,0] = vm[i,0] + dt*dq
            
        for i in range(0,L.shape[0]):
            if i in valve:
                pj[i,0] = vj[i,0]*cham_e[j,valve.index(i)]
            else:
                pj[i,0] = (vj[i,0]-vm[i,0])/C[i] + pjm[i,0]
        q5 = qj2[4,0]
        p6j = (q5-q6m)*dt/c6 + p6m
        q6j = ((p6j-pj[5,0])*dt+l6*q6m)/(r6*dt+l6)
                
        while(flag):
 
            for i in range(0,L.shape[0]):
                if cc_out[i]<0:
                    dp = pj[cc_in[i],0]-p6j
                else:
                    dp = pj[cc_in[i],0]-pj[cc_out[i],0]
                if i in valve:
                    ## note:the flow sign(>0 or <0) determines whether valve close or open
                    ## in some thesis, the pressure difference determines valve opening, because the valve is in series with only resistor
                    ## in this model, valve is in series with resistor and inductor, therefore the flow sign should be used
                    ## considering cases of valve regurgitation, the regurgitation resistor is set to be 1e6
                    ## for regurgitation, the flow should be solved in implicit way, with resistor in denominator; otherwise, solution divergence
                    if qj[i,0]>=0:
                        qj2[i,0] = qjm[i,0] + (buffer[i] + (dp-R[i]*qj[i,0])/L[i])*0.5*dt
                    else:
                        
                        qj2[i,0] = (dp*dt+L[i]*qjm[i,0]) / (L[i]+(R[i]+100000)*dt)
                else:
                    qj2[i,0] = qjm[i,0] + (buffer[i] + (dp-R[i]*qj[i,0])/L[i])*0.5*dt
            
            ## another way to block reverse flow
            # for i in range(0,L.shape[0]):
            #     if cc_out[i]<0:
            #         dp = pj[cc_in[i],0]-p6j
            #     else:
            #         dp = pj[cc_in[i],0]-pj[cc_out[i],0]
            #     if i in valve:
            #         if qj[i,0]>=0:
            #             qj2[i,0] = (dp*dt+L[i]*qjm[i,0]) / (L[i]+R[i]*dt)
            #         else:
            #             qj2[i,0] = 0
            #     else:
            #         qj2[i,0] = (dp*dt+L[i]*qjm[i,0]) / (L[i]+R[i]*dt)
            # # ##have to set zero flow again 
            # for i in range(0,L.shape[0]):
            #     if i in valve:
            #         if qj2[i,0]<0:
            #             qj2[i,0] = 0
            
            for i in range(0,L.shape[0]):
                if cc_in[i+7]<0:
                    dq = (q6j-qj2[cc_out[i+7],0] + q6m-qjm[cc_out[i+7],0])/2
                else:
                    dq = (qj2[cc_in[i+7],0]-qj2[cc_out[i+7],0] + qjm[cc_in[i+7],0]-qjm[cc_out[i+7],0])/2
                vj2[i,0] = vm[i,0] + dt*dq
            
            
            for i in range(0,L.shape[0]):
                if i in valve:
                    pj2[i,0] = vj2[i,0]*cham_e[j,valve.index(i)]
                else:
                    pj2[i,0] = (vj2[i,0]-vm[i,0])/C[i] + pjm[i,0]
            
            q5 = qj2[4,0]
            ## p6 and q6 are pressure and flow for the compartment representing CFD process
            ## p6 and q6 are solved at the last in a weak copuling way.
            p6j2 = (q5-q6j)*dt/c6 + p6m
            q6j2 = ((p6j2-pj2[5,0])*dt+l6*q6m)/(r6*dt+l6)
            
            eps = abs(p6j-p6j2)
            ct = ct+1
            if ct>30 or eps<0.0001:
                flag=0
            else:
                pj = copy.deepcopy(pj2)
                qj = copy.deepcopy(qj2)
                vj = copy.deepcopy(vj2)
                p6j = 0.5*(p6j+p6j2)
                q6j = 0.5*(q6j+q6j2)
                
        Q[j,:] =  qj2.T
        P[j,:] =  pj2.T
        V[j,:] = vj2.T
        p6[j] = p6j2
        q6[j] = q6j2


    
    plt.figure()
    # plt.plot(t,p6,'k')
    plt.plot(t[:],Q[:,4],'k')
    # plt.plot(t[14000:],p6[14000:],'k')
    # plt.plot(t[14000:],P[14000:,4],'r')
    # plt.plot(t,lea,'b')
    # plt.plot(t,lev,'y')
    plt.show()






































