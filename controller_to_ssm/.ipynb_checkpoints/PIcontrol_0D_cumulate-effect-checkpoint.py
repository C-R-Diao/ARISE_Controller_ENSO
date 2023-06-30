
'''
changed based on:
/glade/campaign/cesm/collections/ARISE-SAI-1.5/b.e21.BW.f09_g17.SSP245-TSMLT-GAUSS-DEFAULT.001/controller/PIcontrol.py

Chenrui Diao; 2022-11-01

**********
fix T1 and T2
create T0 with fixed-freq "ENSO" variation based on CESM2-ARISE results (see main.ipynb)
create "online" controller based on the linear relationship between injection amount and SSP gmgst (looking at the ens-mean from ARISE and SSP)

**********
disable the buffer year

'''
import numpy as np
exec(open('./IO_module.py').read())


#### USER-SPECIFIED CONTROL PARAMETERS ####
#refvals=[288.13,0.76,-5.98] reference values for CESM2(WACCM6) (Tilmes et al. 2020) 2020-2039
refvals=[288.64,0.8767,-5.89] # updated to be average over years 2010-2029
kivals=[0.0183,0.0753,0.3120]
kpvals=[0.0183,0.0753,0.3120]

firstyear= 2035
baseyear = 2030
x_ramp = 5.0 # defines a range of years over which the feedback is ramped up

#### USER SPECIFIED CALCULATIONS ####
logfilename='ControlLog_'+runname+'.txt'

logheader=['Timestamp','dT0','sum(dT0)','dT1','sum(dT1)','dT2','sum(dT2)','L0','L1N','L1S','L2','30S(Tg)','15S(Tg)','15N(Tg)','30N(Tg)','Total(Tg)']

firsttime=0
if os.path.exists(maindir+'/'+logfilename)==False:
    firsttime=1
else:
    loglines=readlog(maindir+'/'+logfilename)

T1 = refvals[1]
T2 = refvals[2]

de=np.array([T0-refvals[0],T1-refvals[1],T2-refvals[2]]) # error terms

if firsttime==1:
    timestamp=firstyear
    sumde=de
    sumdt2=de[2]
else:
    timestamp=int(loglines[-1][0])+1
    sumdt0=float(loglines[-1][2])+(T0-refvals[0])
    sumdt1=float(loglines[-1][4])+(T1-refvals[1])
    sumdt2=float(loglines[-1][6])+(T2-refvals[2])
    sumde=np.array([sumdt0,sumdt1,sumdt2])

#print('year: ', timestamp, '; T0: ', T0)

dt=timestamp-baseyear
dt2=timestamp-firstyear
ramp_up = 1.0
if (dt2<x_ramp):
    ramp_up = dt2 / x_ramp

# updated based on feedback simulation
l0hat=0.00347*dt
l1hat=-0.000*dt
l2hat=0.00*dt

# feedback
l0kp1=(kpvals[0]*de[0]+kivals[0]*sumde[0])*ramp_up
l1kp1=(kpvals[1]*de[1]+kivals[1]*sumde[1]-0.5*l0kp1)*ramp_up
l2kp1=(kpvals[2]*de[2]+kivals[2]*sumde[2]-l0kp1)*ramp_up

# all of the feeds
l0step4=l0kp1+l0hat
l1step4=l1kp1+l1hat
l2step4=l2kp1+l2hat

l0=max(l0step4,0)
l1n=min(max(l1step4,0),l0)
l1s=min(max(-l1step4,0),l0)
l2=min(max(l2step4,0),l0-l1s-l1n)
ell=np.array([[l0],[l1n],[l1s],[l2]])

# preventing integrator wind-up
if (l2==(l0-l1s-l1n)):
    sumdt2=sumdt2-(T2-refvals[2])
    sumde[2]=sumdt2

M=np.array([[0,30,30,0],[0,0,45,20],[20,45,0,0],[40,0,0,40]])
F=np.array([[1,1,1,1],[0,1,0,0],[0,0,1,0],[0,0,0,1]])

q=np.dot(np.dot(np.transpose(M),np.linalg.inv(F)),ell)
total_s = sum(q)

for k in range(len(q)):
    q[k]=max(q[k],0)

newline=[str(timestamp),str(de[0]),str(sumde[0]),str(de[1]),str(sumde[1]),str(de[2]),str(sumde[2]),str(l0),str(l1n),str(l1s),str(l2),str(q[0])[1:-1],str(q[1])[1:-1],str(q[2])[1:-1],str(q[3])[1:-1], str(total_s)[1:-1]]
if firsttime==1:
    linestowrite=[logheader,newline]
else:
    linestowrite=[]
    for k in range(len(loglines)):
        linestowrite.append(loglines[k])
    linestowrite.append(newline)

writelog(maindir+'/'+logfilename,linestowrite)

