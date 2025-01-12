import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#import

Y=0.92
S=0.12

Ep=0  #mm/h 
Es=0 #mm/h 
Ep=0.43  #mm/h 
Es=0.43 #mm/h 

p=0.69  #splash pinning
B=0.9 #Retention
q=0 #0 is broad-leaved tree
Ω=0.9 #cluster

LAI=2.23*Ω
SAI=0.154
SAIe=SAI #SAIe
FVC=0.7 
m=0.02
pt=SAIe/(SAIe+LAI)#stem interception ratio

D_d=1  #1 is not considering stem flow

#excel，1 is yes
write=1

G=0.69
if q>0:
    G=0.5


I = [[10], [0.25]]  # intensity(mm/h)，time(h)

start_time = 0.0
end_time = sum(I[1])-0.000001
# t_values
total_duration = end_time-start_time  #duration
num_intervals = int(30)  # 等分的份数sum(I[1])*11
# intervals
interval_duration = total_duration / num_intervals
t_values = [float(interval_duration) * i for i in range(num_intervals + 1)]


#get I
def get_intensity_and_duration(t):
    current_time = 0
    for intensity, duration in zip(I[0], I[1]):
        # print(intensity)
        next_time = current_time + duration
        if current_time <= t <= next_time:
            return intensity
        current_time = next_time
    return None, None  # 如果 t 超出范围，则返回 None 或其他适当的值


#calcylation


#Leaf interception
def modelleaf(t_values, I):
    results1 = []
    water_storage = 0  # 初始化冠层蓄水量为0
    w1=0
    E1=0
    if q>0:
        K=FVC * (1 - (1 - q) ** (LAI*Ω / FVC)) * (1 - pt)
    else:
        K = FVC * (1 - (1 - B * p) ** (LAI*G / FVC)) * (1 - pt)
        # K = 0.1

    for t in t_values:
        current_I = get_intensity_and_duration(t)  # 取出降雨强度列表对应元素作为当前强度
        if t==0:
            t0=0
            t1=0
            evaporation_origin=0
        else:
            if water_storage >= (current_I * K * Y) / (current_I * K + Ep):
                t0=1000
            else:
                t0=-Y / (current_I * K + Ep) * np.log(1 - (water_storage * (current_I*K+Ep)) / (current_I * K * Y ))
            # print("t0",t0)
            t0=round(t0, 3)
            t1=t0+t_values[2]-t_values[1]
            # print("t1",t1)

        
        w1 = (current_I * K * Y) / (current_I * K + Ep) * (1 - np.exp(-((current_I * K + Ep) / Y) * t1))
        
        E0 = ((current_I * K * Ep) / (current_I * K + Ep)) * (t0 + (Y / (current_I * K + Ep)) * (np.exp(-((current_I * K + Ep) / Y) * t0) - 1))
        E2 = ((current_I * K * Ep) / (current_I * K + Ep)) * (t1 + (Y / (current_I * K + Ep)) * (np.exp(-((current_I * K + Ep) / Y) * t1) - 1))
        E1 = E2-E0+evaporation_origin

        result1 = w1 + E1
        results1.append(result1)

        water_storage = w1  # 更新冠层蓄水量
        # print("water_storage",water_storage)
        evaporation_origin=E1 #更新蒸发量
        # print("evaporation_origin",evaporation_origin)

    return results1

#stem interception
def modelstem(t_values, I):
    results2 = []
    Ppts2=[]
    water_storage = 0  # 初始化冠层蓄水量为0
    w1=0
    E1=0
    P0=0
    K = pt*FVC*p
    for t in t_values:
        current_I = get_intensity_and_duration(t)  # 取出降雨强度列表对应元素作为当前强度
        if t==0:
            t0=0
            t1=0
            evaporation_origin=0
        else:
            if water_storage >= (current_I * K * S) / (current_I * K + Es):
                t0=1000
            else:
                t0=-S / (current_I * K + Es) * np.log(1 - (water_storage * (current_I*K+Es)) / (current_I * K * S ))
            # print("t0",t0)
            t0=round(t0, 3)
            t1=t0+t_values[2]-t_values[1]
            # print("t1",t1)

        
        w1 = (current_I * K * S) / (current_I * K + Es) * (1 - np.exp(-((current_I * K + Es) / S) * t1))
        
        E0 = ((current_I * K * Es) / (current_I * K + Es)) * (t0 + (S / (current_I * K + Es)) * (np.exp(-((current_I * K + Es) / S) * t0) - 1))
        E2 = ((current_I * K * Es) / (current_I * K + Es)) * (t1 + (S / (current_I * K + Es)) * (np.exp(-((current_I * K + Es) / S) * t1) - 1))
        E1 = E2-E0+evaporation_origin

        result2 = w1 + E1
        results2.append(result2)

        #降雨量
        if t==0:
            P=0
        else:
            P=current_I*(t_values[2]-t_values[1])+P0
        
        #茎流量
        Ppt2=(P*K-result2)*(1-D_d)
        Ppts2.append(Ppt2)

        water_storage = w1  # 更新冠层蓄水量
        # print("water_storage",water_storage)
        evaporation_origin=E1 #更新蒸发量
        # print("evaporation_origin",evaporation_origin)
        P0=P #更新P0

    return results2,Ppts2


wl=modelleaf(t_values, I)
ws, Ppt=modelstem(t_values, I)


#splash evaporation
def splashevaporation(t_values, I):
    Ef_org=0
    Efs=[]
    i=0
    # print("Kef",Kef)
    for t in t_values:
        current_I = get_intensity_and_duration(t)  # 取出降雨强度列表对应元素作为当前强度
        if q>0:
            Kef=FVC*m*(1-q)*(1-pt)*((1-(1-q)**(LAI*G/FVC)))
            
        else:
            Kef=FVC*m*(1-pt)*(1-p)*((1-(1-B*p)**(LAI*G/FVC))/(B*p))
            
        i=i+1

        if t==0:
            Ef=0
        else:
            Ef=(Kef)*current_I*(t_values[2]-t_values[1])+Ef_org
        Ef_org=Ef
        Efs.append(Ef)
    return Efs

Ef=splashevaporation(t_values, I)


w = np.add(wl, ws)
w1 = np.add(w, Ef)


# DataFrame
df = pd.DataFrame({
    't_values/h': t_values,
    'wl/mm': wl,
    'ws/mm': ws,
    'w/mm':w,
    'w1/mm': w1,
    'Ppt/mm':Ppt,
    'Ef/mm':Ef
})

# output Excel
if write==1:
    df.to_excel('w_model_output.xlsx', index=False)