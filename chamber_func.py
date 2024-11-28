# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 09:34:57 2024

@author: SSKJCFD004
"""

# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import os 
import copy

def cal_elastance(tt,c_type:str):
    elas = 0.1
    t = tt%Tc
    Ts1_a,Ts2_a = 0.92*Tc,0.08*Tc
    Ts1_v,Ts2_v = 0.375*Tc,0.5625*Tc;
    
    if c_type == 'ra' or c_type == 'la':
        if t<Ts1_a:
            elas = 0
        elif t>=Ts1_a and t<Ts1_a+Ts2_a:
            elas = 1-np.cos((t-Ts1_a)/Ts2_a * 2*np.pi)
        else:
            elas = 0
    elif (c_type == 'rv'):
        if t<Ts1_v:
            elas = (1-np.cos(t/Ts1_v* np.pi))
        elif t>=Ts1_v and t<Ts2_v:
            elas = (1+ np.cos((t-Ts1_v)/(Ts2_v-Ts1_v) * np.pi))
        else:
            elas = 0.
    elif (c_type == 'lv'):
        if t<Ts1_v:
            elas =  (1-np.cos(t/Ts1_v* np.pi))
        elif t>=Ts1_v and t<Ts2_v:
            elas =  (1+ np.cos((t-Ts1_v)/(Ts2_v-Ts1_v) * np.pi))
        else:
            elas = 0.
    else:
        elas = 0.001
        
    return elas

def cal_ela2(tt):
    t = tt%Tc
    A_,B_,C_ = 1.0,80.0,0.27
    elas = A_*np.exp(-B_*(t-C_)**2)
    return elas
    
def cal_Pes(v_,vd_,Ees_):
    return Ees_*(v_-vd_)

def cal_Ped(v_,v0_,A_,lam_):
    ped = A_*(np.exp(lam_*(v_-v0_))-1)
    if ped>250:
        ped = 250
    return ped
    
def cal_pump_r(pre):
    if pre<1:
        r = -3.5*(pre-1.0)
    else:
        r=0
    return r

def pump_HQ(flow_,rpm_):
    if flow_<=0.0:
        f_=0.0
    elif flow_>2000:
        f_ = 2000
    else:
        f_ = flow_
    a = 7.3e-6*rpm_**2
    b = -1.2e-5 * rpm_*f_
    c = -1.8e-3*f_*f_
    return a+b+c

def cal_buf_p2(pm_,qm_,p6m_,qpm_,p6lm_,rpm_):
    buffer_ = np.zeros(len(pn_in))
    bufferpp_ = 0
    for i in range(0,len(pn_in)):
        if L[i]>0:
            if pn_out[i]<0 :
                dp = pm_[pn_in[i]]-p6m_
                buffer_[i] = (dp-R[i]*qm_[i])/L[i]
            else:
                dp = pm_[pn_in[i]]-pm_[pn_out[i]]
                buffer_[i] = (dp-R[i]*qm_[i])/L[i]
    if rpm_>0:
        dpp_ = pm_[lv]-p6lm_ + pump_HQ(qpm_,rpm_)
        vpr_ = cal_pump_r(pm_[lv])
        bufferpp_ = (dpp_ - (PR+vpr_)*qpm_)/PL
    else:
        bufferpp_ = 0.0
    # print(buffer_)
    return buffer_,bufferpp_

def lpn(buffer_,bufferpp_,
        pj2_,pm_,qj2_,qm_,vj2_,vm_,
        qpj2_,qpm_,p6j2_,p6m_,q6j2_,q6m_,q6lj2_,q6lm_,p6lj2_,p6lm_,
        rpm_,j_):
    ## solve flow and pressure for one time step
    pj_ = copy.deepcopy(pj2_)
    qj_ = copy.deepcopy(qj2_)
    vj_ = copy.deepcopy(vj2_)
    p6j_ = copy.deepcopy(p6j2_)
    qpj_ = copy.deepcopy(qpj2_)
    q6j_ = copy.deepcopy(q6j2_)
    q6lj_ = copy.deepcopy(q6lj2_)
    p6lj_ = copy.deepcopy(p6lj2_)
    for i in range(0,len(pn_in)):
        if L[i]>0:
            if pn_out[i]<0:
                dp = pj_[pn_in[i]]-p6j_
            else:
                dp = pj_[pn_in[i]]-pj_[pn_out[i]]
        else:
            if pn_out[i]<0 :
                dp = pj_[pn_in[i]]-p6j_
            else:
                dp = pj_[pn_in[i]]-pj_[pn_out[i]]
        if L[i] >0:
            if i in valve:
                if qj_[i]>0:
                    # qj_[i] = qm_[i] + (buffer_[i] + (dp-R[i]*qj_[i])/L[i])*0.5*dt
                    qj_[i] = (dp*dt+L[i]*qm_[i]) / (L[i]+(R[i])*dt)
                else:
                    qj_[i] = (dp*dt+L[i]*qm_[i]) / (L[i]+(R[i]+1000000000.)*dt)
            else:
                # qj_[i] = qm_[i] + (buffer_[i] + (dp-R[i]*qj_[i])/L[i])*0.5*dt
                qj_[i] = (dp*dt+L[i]*qm_[i]) / (L[i]+(R[i])*dt)
        else:
            if i in valve:
                if qj_[i]>0:
                    qj_[i] = dp/R[i]
                else:
                    qj_[i] = dp/(R[i]+100000000)
            else:
                qj_[i] = dp/R[i]
    
    if rpm_>0:
        dpp = pj_[lv]-p6j_ + pump_HQ(qpj_,3000)
        vpr = cal_pump_r(pj_[lv])
        if qpj_>=0:
            qpj_ = qpm_ + (bufferpp_ + (dpp - (PR+vpr)*qpj_)/PL)*0.5*dt
        else:
            qpj_ = 0
    else:
        qpj_=0
    
    q6lj_ =  ((p6j_-p6lj_)*dt+l6*q6lm_) / (l6+r6*dt)
    q6j_[7] = ((p6lj_-pj_[7])*dt+l6_u*q6m_[7]) / (l6_u+r6_u*dt)
    q6j_[8] = ((p6lj_-pj_[8])*dt+l6_d*q6m_[8]) / (l6_d+r6_d*dt)
    # delta1 = abs(q6lj_+qpj_-q6j_[7]-q6j_[8]) ## q6l is only from the ventricle
    delta1 = abs(q6lj_-q6j_[7]-q6j_[8]) ## q6l is the sum of Qventricle and Qpump
    
    p6lj2 = p6lj_*1.01+0.01
    q6lj2 =  ((p6j_-p6lj2)*dt+l6*q6lm_) / (l6+r6*dt)
    q6j7 = ((p6lj2-pj_[7])*dt+l6_u*q6m_[7]) / (l6_u+r6_u*dt)
    q6j8 = ((p6lj2-pj_[8])*dt+l6_d*q6m_[8]) / (l6_d+r6_d*dt)
    # delta = abs(q6lj2+qpj_-q6j7-q6j8)
    delta = abs(q6lj2-q6j7-q6j8)
    k = (delta-delta1)/(p6lj2-p6lj_)
    flag_i,ci=1,0
    
    while(flag_i):
        if k==0:
            p6lj_ = p6lj_- k
        else:
            p6lj_ = p6lj_- delta1/k
        q6lj_ =  ((p6j_-p6lj_)*dt+l6*q6lm_) / (l6+r6*dt)
        q6j_[7] = ((p6lj_-pj_[7])*dt+l6_u*q6m_[7]) / (l6_u+r6_u*dt)
        q6j_[8] = ((p6lj_-pj_[8])*dt+l6_d*q6m_[8]) / (l6_d+r6_d*dt)
        # delta1 = abs(q6lj_+qpj_-q6j_[7]-q6j_[8])
        delta1 = abs(q6lj_-q6j_[7]-q6j_[8])
        
        p6lj2 = p6lj_*1.01+0.01
        q6lj2 =  ((p6j_-p6lj2)*dt+l6*q6lm_) / (l6+r6*dt)
        q6j7 = ((p6lj2-pj_[7])*dt+l6_u*q6m_[7]) / (l6_u+r6_u*dt)
        q6j8 = ((p6lj2-pj_[8])*dt+l6_d*q6m_[8]) / (l6_d+r6_d*dt)
        # delta = abs(q6lj2+qpj_-q6j7-q6j8)
        delta = abs(q6lj2-q6j7-q6j8)
        k = (delta-delta1)/(p6lj2-p6lj_)
        ci = ci+1
        # print(delta1,delta,k,percent,percent2)
        if abs(delta1)<0.001 or ci>20:
            flag_i = 0
        # print(k,delta1,j_,ci)

    for i in range(0,len(pn_in)):
        if para_p[i]==1:
            dq = (qj_[qn_in[i]]-qj_[qn_out[i]]-qpj_ + qm_[qn_in[i]]-qm_[qn_out[i]]-qpm_)/2
        elif para[i]==1:
            index = list(np.where(para==-para[i])[0])
            qin=0.0
            qinm = 0.0
            for k in index:
                qin += qj_[k]
                qinm += qm_[k]
            dq = (qin-qj_[qn_out[i]] + qinm-qm_[qn_out[i]])/2.0
        elif qn_in[i]>=0:
            dq = (qj_[qn_in[i]]-qj_[qn_out[i]] + qm_[qn_in[i]]-qm_[qn_out[i]])/2
        elif qn_in[i]==-10:
            dq = (q6j_[i]-qj_[qn_out[i]] + q6m_[i]-qm_[qn_out[i]])/2
        
        vj_[i] = vm_[i] + dt*dq
        
    for i in range(0,len(pn_in)):
        if i in ven:
            pj_[i] = vj_[i]*cham_e[j_,valve.index(i)]
        else:
            pj_[i] = (vj_[i]-vm_[i])/C[i] + pm_[i]
    
    # p6j_ = (qj_[lv] -q6lj_ + qm_[lv] -q6lm_)*0.5*dt/c6 + p6m_ ## q6l is only from the ventricle
    p6j_ = (qj_[lv]+qpj_ -q6lj_ + qm_[lv]+qpm_ -q6lm_)*0.5*dt/c6 + p6m_ ## q6l is the sum of Qventricle and Qpump
    
    return p6j_,qpj_,pj_,qj_,vj_,q6j_,q6lj_,p6lj_

def solve_lpn(Tc_,rpm_,n_tc_,n_branch_):
    ##solve the flow and pressure of a single cardiac cycle
    #%%time  parameter
    dt= Tc_/n_tc_
    n_cycle = 0
    
    #%%define PQV array
    P_ = np.zeros((len(time),n_branch_)) #result storage array for LPN pressure mmHg
    Q_ = np.zeros((len(time),n_branch_)) #result storage array for LPN flow ml/s
    V_ = np.zeros((len(time),n_branch_)) ##volume of chambers
    p6_ = np.zeros(len(time))
    q6_ = np.zeros((len(time),n_branch_))
    q6l_ = np.zeros(len(time))
    p6l_ = np.zeros(len(time))
    qp_ = np.zeros(len(time))
    buffer = np.zeros(n_branch_)

    #%%initial arrays
    Q_[0,:] = 0  ## initial flow  0ml/s
    V_[0,0],V_[0,1],V_[0,2],V_[0,3] = 90,236.25,0,0
    V_[0,4],V_[0,5],V_[0,6],V_[0,7] = 0,212,360,0
    
    v_rest = V_tot-sum(V_[0,0:8])
    
    V_[0,11],V_[0,12] = v_rest/2,v_rest/2
 
    for i in range(0,len(pn_in)):
        if i in ven:
            P_[0,i] = V_[0,i]*cham_e[0,valve.index(i)]
        else:
            P_[0,i] = V_[0,i]/C[i]
    
    
    #%% solve PQV for the first time
    eps_c = 100
    for j in range(1,n_tc+1):
        q6m = q6_[j-1,:]
        p6m = p6_[j-1]
        pm = P_[j-1,:]
        qm = Q_[j-1,:]
        vm = V_[j-1,:]
        qpm = qp_[j-1]
        p6lm = p6l_[j-1]
        q6lm = q6l_[j-1]
        pj = copy.deepcopy(pm)
        qj = copy.deepcopy(qm)
        vj = copy.deepcopy(vm)
        p6j = copy.deepcopy(p6m)
        q6j = copy.deepcopy(q6m)
        qpj = copy.deepcopy(qpm)
        p6lj = copy.deepcopy(p6lm)
        q6lj = copy.deepcopy(q6lm)
        
        buffer,bufferpp=cal_buf_p2(pm,qm,p6m,qpm,p6lm,rpm_)
        ##solve the pressure and flow of current timestep for the first time
        p6j,qpj,pj,qj,vj,q6j,q6lj,p6lj = lpn(buffer,bufferpp,pj,pm,qj,qm,
                                                vj,vm,qpj,qpm,p6j,p6m,q6j,
                                                q6m,q6lj,q6lm,p6lj,p6lm,rpm,j)
        ##solve the pressure and flow of current timestep for the second time
        p6j,qpj,pj,qj,vj,q6j,q6lj,p6lj = lpn(buffer,bufferpp,pj,pm,qj,qm,
                                                vj,vm,qpj,qpm,p6j,p6m,q6j,
                                                q6m,q6lj,q6lm,p6lj,p6lm,rpm,j)
        Q_[j,:] =  qj[:]
        P_[j,:] =  pj[:]
        V_[j,:] =  vj[:]
        qp_[j] = qpj
        p6_[j] = p6j
        q6_[j,:] = q6j
        q6l_[j] = q6lj
        p6l_[j] = p6lj
    
    #%% solve PQV till eps_c is small 
    while(eps_c>0.0005 and n_cycle<40):
        temp = p6_[1:n_tc+1]
        avgp1 = sum(temp)*dt/Tc
        temp = q6l_[1:n_tc+1]
        avgp1 = avgp1 + sum(temp)*dt/Tc
        
        Q_[0,:] =  Q_[n_tc,:]
        P_[0,:] =  P_[n_tc,:]
        V_[0,:] =  V_[n_tc,:]
        qp_[0] = qp_[n_tc]
        p6_[0] = p6_[n_tc]
        q6_[0,:] = q6_[n_tc,:]
        q6l_[0] = q6l_[n_tc]
        p6l_[0] = p6l_[n_tc]
        
        
        for j in range(1,n_tc+1):
            q6m = q6_[j-1,:]
            p6m = p6_[j-1]
            pm = P_[j-1,:]
            qm = Q_[j-1,:]
            vm = V_[j-1,:]
            qpm = qp_[j-1]
            p6lm = p6l_[j-1]
            q6lm = q6l_[j-1]
            pj = copy.deepcopy(pm)
            qj = copy.deepcopy(qm)
            vj = copy.deepcopy(vm)
            p6j = copy.deepcopy(p6m)
            q6j = copy.deepcopy(q6m)
            qpj = copy.deepcopy(qpm)
            p6lj = copy.deepcopy(p6lm)
            q6lj = copy.deepcopy(q6lm)
            
            buffer,bufferpp=cal_buf_p2(pm,qm,p6m,qpm,p6lm,rpm_)
            
            p6j,qpj,pj,qj,vj,q6j,q6lj,p6lj = lpn(buffer,bufferpp,pj,pm,qj,qm,
                                                    vj,vm,qpj,qpm,p6j,p6m,q6j,
                                                    q6m,q6lj,q6lm,p6lj,p6lm,rpm,j)
            
            p6j,qpj,pj,qj,vj,q6j,q6lj,p6lj = lpn(buffer,bufferpp,pj,pm,qj,qm,
                                                    vj,vm,qpj,qpm,p6j,p6m,q6j,
                                                    q6m,q6lj,q6lm,p6lj,p6lm,rpm,j)
            Q_[j,:] =  qj[:]
            P_[j,:] =  pj[:]
            V_[j,:] =  vj[:]
            qp_[j] = qpj
            p6_[j] = p6j
            q6_[j,:] = q6j
            q6l_[j] = q6lj
            p6l_[j] = p6lj
            
        
        temp = p6_[1:n_tc+1]
        avgp2 = sum(temp)*dt/Tc
        temp = q6l_[1:n_tc+1]
        avgp2 = avgp2 + sum(temp)*dt/Tc
        eps_c = abs(avgp1-avgp2)
        n_cycle += 1
        
        if (eps_c>10000 and n_cycle>=9):
            P_= None
            break
    print(eps_c,n_cycle)
    
    return P_,Q_,V_,qp_,p6_,q6_,p6l_,q6l_

def check_timestep(dt_input,n_branch_):
    c_=0.1
    flag_ = 1
    dt_min = np.zeros(n_branch_+1)+10.1
    dt_ = dt_input
    count_dt = 0
    count_max_ = []
    for i in range(0,n_branch_+1):
        flag_=1
        dt_ = dt_input
        if i==0 or i==5:
            c_ = 1/(Echam[0,0]+Echam[0,1]) if i==0 else 1/(Echam[0,2]+Echam[0,3])
        elif i==1 or i==6:
            c_ = 1/(Echam[1,0]) if i==1 else 1/(Echam[1,1])
        elif i==n_branch_:
            c_ = c6
        else:
            c_ = C[i]
        
        if i==n_branch_:
            li,ri,lo,ro = L[6],R[6],l6,r6
        elif para[i]==1:
            index = list(np.where(para==-para[i])[0])
            qin,qout = index[0],i
            li,ri,lo,ro = L[qin],R[qin],L[qout],R[qout]
        elif qn_in[i]>=0:
            qin,qout = qn_in[i],qn_out[i]
            li,ri,lo,ro = L[qin],R[qin],L[qout],R[qout]
        elif qn_in[i]==-10:
            qin,qout=-i,qn_out[i]
            if i==7:
                li,ri,lo,ro = l6_u,r6_u,L[qout],R[qout]
            elif i==8:
                li,ri,lo,ro = l6_d,r6_d,L[qout],R[qout]
        
        
        count_dt = 0
        while(flag_==1):
            if (li>0 and ri>0) and (lo>0 and ro>0):
                A = np.zeros((3,3))
                B = np.zeros((3,3))
                A[0,0],A[1,1],A[2,2] = li,lo,c_
                B[0,0],B[0,2] = -ri,-1.0
                B[1,1],B[1,2] = -ro,1.0
                B[2,0],B[2,1] = 1.0,-1.0
            elif (li>0 and ri>0) and (lo<0 and ro>0):
                A = np.zeros((2,2))
                B = np.zeros((2,2))
                A[0,0],A[1,1] = li,c_
                B[0,0],B[0,1] = -ri,-1.0
                B[1,0],B[1,1] = 1.0,-1.0/ro
            elif (li<0 and ri>0) and (lo>0 and ro>0):
                A = np.zeros((2,2))
                B = np.zeros((2,2))
                A[0,0],A[1,1] = lo,c_
                B[0,0],B[0,1] = -ro,1.0
                B[1,0],B[1,1] = -1.0,-1.0/ri
            elif (li<0 and ri>0) and (lo<0 and ro>0):
                A = np.zeros((1,1))
                B = np.zeros((1,1))
                A[0,0] = c_
                B[0,0] = -(1.0/ri + 1.0/ro)
            D=np.linalg.inv(A)*B*dt_
            I = np.diag(np.ones(A.shape[0]))
            G = 0.5*D*D+D+I
            e_value,e_vector = np.linalg.eig(G)
            count_dt +=1
            if (max(abs(e_value))<=1 or count_dt>5):
                flag_=0
            else:
                dt_ = dt_/2
            
        print(e_value,count_dt,i)
        dt_min[i] = dt_
        count_max_.append(count_dt)
    dt_ = min(dt_min)
    return dt_,count_max_
    
def cal_target():
    Vlv = V[0:n_tc,lv]
    Vrv = V[0:n_tc,rv]
    Pao = p6[:n_tc+1]
    target_e_ = {}
    target_ = {}
    
    edvi = max(Vlv)/1.41
    target_e_['EDVI'] = abs((EDVI-edvi)/(EDVI+edvi)*2)
    target_['EDVI'] = edvi
    
    rvef = (max(Vrv)-min(Vrv))/max(Vrv)
    target_e_['RVEF'] = abs((rvef-RVEF)/(rvef+RVEF)*2)
    target_['RVEF'] = rvef
    
    svi = (max(Vlv)-min(Vlv))/1.41
    target_e_['SVI'] = abs((svi-SVI)/(svi+SVI)*2)
    target_['SVI'] = svi
    
    dp = min(Pao)
    target_e_['DP'] = abs((dp-DP)/(dp+DP)*2)
    target_['DP'] = dp
    
    sp = max(Pao)
    target_e_['SP'] = abs((sp-SP)/(sp+SP)*2)
    target_['SP'] = sp
    
    co = (max(Vlv)-min(Vlv)) * (60/Tc) /1000
    target_e_['CO'] = abs((co - CO)/(co + CO)*2)
    target_['CO'] = co
    
    wf1 = (Pao-dp)/(sp-dp)
    t1 = np.arange(0,1+0.5/n_tc,1/n_tc)
    t2 = target_wf[:,0]
    wf2 = np.zeros_like(t2)
    for i in range(0,t2.shape[0]):
        dis_m=10
        for j in range(0,t1.shape[0]):
            dis = abs(t2[i]-t1[j])
            if dis<dis_m:
                dis_m=dis
                j_m=j
        if j_m == t1.shape[0]-1 or t2[i]<=t1[j_m]:
            wf2[i] = (wf1[j_m]-wf1[j_m-1])/(t1[j_m]-t1[j_m-1]) * (t2[i]-t1[j_m-1]) + wf1[j_m-1]
        elif j_m==0 or t2[i,0]>t1[j_m]:
            wf2[i] = (wf1[j_m]-wf1[j_m+1])/(t1[j_m]-t1[j_m+1]) * (t2[i]-t1[j_m+1]) + wf1[j_m+1]
    
    x = list(target_wf[0:-1,1])
    y = list(wf2[0:-1])
    x = x+x
    y = y+y
    mx,my = np.mean(x),np.mean(y)
    x = np.array(x)
    y = np.array(y)
    #NCC normalized cross correlation
    corr = np.correlate(x-mx,y-my, mode='full')
    norm_f = np.sqrt(np.sum((x-mx)**2)*np.sum((y-my)**2))
    ncc = corr/norm_f 

    lag = np.argmax(ncc) - (x.shape[0]-1)
    max_cor = ncc[lag+x.shape[0]-1]
    # print(lag,max_cor)

    y3 = []
    y2=list(y)
    y3=y3+y2[-(lag+0):]
    y3+=y2[0:-(lag+0)]
    err = 0
    sum1=0
    sum2 = 0
    for i in range(0,len(y3)):
        sum1 += x[i]
        sum2 += y3[i]
        err += abs(x[i]-y3[i])
    # target_e_['err'] = err/len(y3)
    target_e_['err'] = err/abs(sum1)
    target_e_['y3'] = y3
    target_e_['max_cor'] = max_cor
    
    return target_,target_e_

if __name__ =="__main__":
    os.chdir(r'F:\linzhihong_doc\erke\python')
    #%% define basic time related parameter
    if(1):
        Tc =0.8
        rpm= 3000
        n_tc = 800 #number of timestep in on cycle
        dt= Tc/n_tc #initial timestep or maximum timestep
        n_cycle = 0
        n_branch=13
    #%% target data: aortic pressure waveform, peak pressure SP, diastolic pressure DP, LEDVI, CVP, RVEF,SVI,CO,EDVI
    wf = np.load('wave_form.npz')
    target_wf = wf['wf_n1']
    wf.close()
    SP = 90
    DP = 55
    CO = 2.9
    EDVI = 140
    RVEF = 0.13
    SVI = 19
    
    #%%parameter fixed and not optimized
    ## LR for 3Dmodel,fixed
    l6,r6 = 5.5E-4,0.012
    l6_u,r6_u = 0.01558,0.1079
    l6_d,r6_d = 0.01042,0.07213
    ## LR for pump,fixed
    PL = 0.0127*2
    PR = 0.0677*2
    #%% PQ connectivity fixed
    pn_in =  [0,1,2,3,4,5,  6,7,  8, 9,10,11,12] ##pin node for Qnode
    pn_out = [1,2,3,4,5,6,-10,9, 10,11,12, 0, 0]  ##pout node for Qnode
    para   = [1,0,0,0,0,0,  0,0,  0, 0, 0,-1, -1]  ## 1 mean Qin node is sum of parallel branch
    para_p = [0,0,0,0,0,0,  1,0,  0, 0, 0,0, 0]  ## 1 mean Qin node is sum of parallel branch
    qn_in =  [-100,0,1,2,3,4,5, -10,-10,7, 8, 9,10] ##Qin node for Pnode
    qn_out = [0 ,1,2,3,4,5,6,   7,  8,9,10,11,12] ##Qout node for Pnode
    para = np.array(para)
    valve = [0,1,5,6] #Q node for valve
    rv,lv = 1,6
    ven = [0,1,5,6]
    
    #%% initial set chamber parameter, need a function to define Echam
    ##optimized Echam
    Echam = np.zeros((2,4))
    Echam[0,0],Echam[0,1],Echam[0,2],Echam[0,3]=0.05,0.15,0.1,0.15 ## elatance for atriun
    Echam[1,0],Echam[1,1],Echam[1,2],Echam[1,3] = 1.2,0.1,0.525,0.1 ## elatance for ventricle
    
    es_la,ed_la,es_ra,ed_ra = Echam[0,0],Echam[0,1],Echam[0,2],Echam[0,3]
    es_lv,ed_lv,es_rv,ed_rv = Echam[1,0],Echam[1,1],Echam[1,2],Echam[1,3]
    
    #%% initial LPN parameter, need a function to define LRC_para
    # 0-12:0ra, 1rv, 2pas, 3pat, 4pvn, 5la, 6lv, 7sat_up, 8sat_down, 9sar_up, 10sar_down, 11svc_up, 12svc_down
    # c6:Cao; l6_u,r6_u: RL_ao_up; l6_d,r6_d: RL_ao_down
    L = np.zeros(n_branch) #parameter of indutance
    R = np.zeros(n_branch) # parameter of resistance
    C = np.zeros(n_branch) #parameter of compliance of elastance
    L[0],R[0],C[0] = -1, 0.003,-1
    L[1],R[1],C[1] = 1e-5, 0.003,-1
    L[2],R[2],C[2]  =  1.56e-5, 0.002,0.18
    L[3],R[3],C[3] = 0.0017, 0.135,4.75
    L[4],R[4],C[4] = -1, 0.015,20.5
    L[5],R[5],C[5]  = -1, 0.003,-1
    L[6],R[6],C[6]  = 1e-5, 0.003,-1
    L[7],R[7],C[7]  = 0.0034, 0.15,0.0875
    L[8],R[8],C[8]  = 0.0034, 0.075,0.6125
    L[9],R[9],C[9]  = -1, 2.325,0.0625
    L[10],R[10],C[10]  = -1, 1.55,0.4375
    L[11],R[11],C[11]  = -1, 0.0375,8.2
    L[12],R[12],C[12]  = -1, 0.025,12.3
    c6 = 0.3 
    V_tot = 900
    # optimized parameter in LRC_para
    LRC_para = {}
    LRC_para['L'] = L
    LRC_para['R'] = R
    LRC_para['C'] = C
    LRC_para['C2'] = [c6]
    
    #%%
    dt2,count_max = check_timestep(dt,n_branch)
    # dt = dt2
    n_tc = int(np.ceil(Tc/dt))
    time = np.arange(0,3*Tc,dt) 
    cham_e = np.zeros((len(time),4))
    for i in range(0,len(time)):
        cham_e[i,2] = es_la*cal_elastance(time[i],'la')+ed_la
        cham_e[i,3] = es_lv*cal_elastance(time[i],'lv')+ed_lv
        cham_e[i,0] = es_ra*cal_elastance(time[i],'ra')+ed_ra
        cham_e[i,1] = es_rv*cal_elastance(time[i],'rv')+ed_rv
    
    P,Q,V,qp,p6,q6,p6l,q6l = solve_lpn(Tc,rpm,n_tc,n_branch)
    if P is None:
        print("program not converged, this group of parameter should be abandoned")
    
    target,target_e = cal_target()
    
    delta = 0.3*target_e['EDVI'] + 0.3*target_e['RVEF'] + 0.3*target_e['SVI'] + target_e['SP'] + target_e['DP'] + target_e['CO'] + target_e['err'] + 0.3*(1-target_e['max_cor'])
    
    # target=np.zeros(3)
    
    plt.figure()
    # plt.plot(time,cham_e[:,1]/2,'k')
    # plt.plot(time,cham_e[:,3]/2,'r')
    # plt.plot(time,cham_e[:,0]/2,'k')
    # plt.plot(time,cham_e[:,2]/2,'r')
    plt.plot(time,Q[:,6],'k')
    plt.plot(time,q6l,'r')
    plt.plot(time,qp,'b')
    plt.show()
    
    # np.savez('interm_ven.npz',Q = Q, P = P, V=V, 
    #           p6=p6, q6=q6, qp = qp,p6l=p6l,q6l=q6l,
    #           vspt=vspt,p_peri = p_peri,p3w = p3w)
    np.savez('interm_ven.npz',Q = Q, P = P, V=V, 
              p6=p6, q6=q6, qp = qp,p6l=p6l,q6l=q6l)
   





















