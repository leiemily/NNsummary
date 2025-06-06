import deepxde as dde
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from scipy.io import savemat

def gen_traindata():
    data = scipy.io.loadmat('dataset/frac2-10tpb-both1-free.mat')
    
    s, f = data["aS_save2"], data["aF_save2"] 
    dsdt, dfdt = data["adsdt_save2"], data["adFdt_save2"]  

    S = np.reshape(s, (-1, 1)) 
    F = np.reshape(f, (-1, 1))
    dS = np.reshape(dsdt, (-1, 1)) 
    dF = np.reshape(dfdt, (-1, 1))
    return np.hstack((F, S, dF, dS)) 

def ex_func(s):
    data2 = scipy.io.loadmat('dataset/frac2-10tpb-both1-free.mat')
    
    ds2dt2, df2dt2 = data2["ads2dt2_save2"], data2["adF2dt2_save2"] 

    #reshape   
    dS2dt2 = np.reshape(ds2dt2, (-1, 1))
    dF2dt2 = np.reshape(df2dt2, (-1, 1))
    return np.hstack((dS2dt2, dF2dt2))


def G34nn(x, y, ex):
    ds2dt2 = ex[:, 0:1]
    dF2dt2 = ex[:, 1:2]
    
    G3_output = y[:, 0:1]
    G4_output = y[:, 1:2]
    
    eq1 = dF2dt2 - (G3_output + G4_output * ds2dt2)
    
    return eq1


# define domain
geom = dde.geometry.geometry_nd.Hypercube(xmin=[0, 0, 0, 0], xmax=[1, 1, 1, 1])

# Get the train data
observe_x = gen_traindata()


data = dde.data.PDE(
    geom, 
    G34nn, 
    [],
    num_domain=0, 
    num_boundary=0,
    anchors=observe_x,
    auxiliary_var_function=ex_func,
)

layer_size = [4] + [20]*4 + [2]
activation = "tanh"
initializer = "Glorot uniform"
net = dde.maps.FNN(layer_size, activation, initializer)


model = dde.Model(data, net)
model.compile("adam", lr=1e-3, loss="MSE")


def pred1_f(internal_period, total_epoch, mat_file_path, llsi, test_trajectory_numbers):
    data = scipy.io.loadmat(mat_file_path)

    
    f, df = data["F_save2"], data["dFdt_save2"]
    s, ds, ds2 = data["S_save2"], data["dsdt_save2"], data["ds2dt2_save2"]
    
    l2, lls2 = data["Length2"].T, data["Length2lls"].T  
    dt, tscale = data["dt00"], data["t_scale"]
    dt = dt/tscale    
    num_f = data["Numberof_F"]
    num_len = int(num_f)
    
    
    S = np.reshape(s, (-1, 1)) 
    F = np.reshape(f, (-1, 1))
    dF = np.reshape(df, (-1, 1))
    dS = np.reshape(ds, (-1, 1))
    dS2 = np.reshape(ds2, (-1, 1))
    
    
    """test process"""
    len_all = len(s.T)
    
    f_pre = np.zeros((len_all,1))
    df_pre = np.zeros((len_all,1))
    df2_pre = np.zeros((len_all,1))
    g3_pre = np.zeros((len_all,1))
    g4_pre = np.zeros((len_all,1))
    model_number = int(total_epoch/internal_period)
    test_errors = np.zeros((model_number + 1, 1))
    
    for n in range(model_number + 1):
        model_number_n = n * internal_period
        model.restore("my2dfreeModel/2dfree_model.ckpt-" + str(model_number_n), verbose=1)

        for j in range(num_len):
            for i in range(llsi-1):  
                f_pre[l2[j]+lls2[i,j]] = F[l2[j]+lls2[i,j]]
                f_pre[l2[j]+lls2[i,j]+1] = F[l2[j]+lls2[i,j]+1]
                for k in range(int(l2[j]+lls2[i,j]+1),int(l2[j]+lls2[i+1,j]-1)):
                    f1 = np.array([f_pre[k]])
                    s1 = np.array([S[k]])
                
                    df_pre[k] = (f_pre[k]-f_pre[k-1])/dt
                    df1 = np.array([df_pre[k]])
                    ds1 = np.array([dS[k]])
                
                    X_pre = np.hstack((f1,s1,df1,ds1))
                    g_pre = model.predict(X_pre)
                    
                    g3_pre[k] = g_pre[0,0]
                    g4_pre[k] = g_pre[0,1]
                    df2_pre[k] = g3_pre[k] + dS2[k] * g4_pre[k]
                    
                    f_pre[k+1] = df2_pre[k] * dt**2 + 2*f_pre[k] - f_pre[k-1]
                    
        """"""        
        for j in range(num_len):
            for i in range(llsi-1): #start
                k = int(l2[j]+lls2[i,j])
                f1 = np.array([f_pre[k]])
                s1 = np.array([S[k]])
            
                df_pre[k] = (f_pre[k+1]-f_pre[k])/dt
                df1 = np.array([df_pre[k]])
                ds1 = np.array([dS[k]])
            
                X_pre = np.hstack((f1,s1,df1,ds1))
                g_pre = model.predict(X_pre)
                
                g3_pre[k] = g_pre[0,0]
                g4_pre[k] = g_pre[0,1]
                df2_pre[k] = g3_pre[k] + dS2[k] * g4_pre[k]
                
                
                k = int(l2[j]+lls2[i+1,j]-1)#end
                f1 = np.array([f_pre[k]])
                s1 = np.array([S[k]])
            
                df_pre[k] = (f_pre[k]-f_pre[k-1])/dt
                df1 = np.array([df_pre[k]])
                ds1 = np.array([dS[k]])
            
                X_pre = np.hstack((f1,s1,df1,ds1))
                g_pre = model.predict(X_pre)
            
                g3_pre[k] = g_pre[0,0]
                g4_pre[k] = g_pre[0,1]
                df2_pre[k] = g3_pre[k] + dS2[k] * g4_pre[k]
                
        """对每个saved NN，计算测试集的MSE"""
        test_start_j = num_len - test_trajectory_numbers
        f_true1 = F[int(l2[test_start_j]) : int(l2[num_len])]
        f_pre1 = f_pre[int(l2[test_start_j]) : int(l2[num_len])]
        
        
        mse = np.mean((f_true1 - f_pre1) ** 2)
        test_errors[n] = mse
    

    
    """train process 用的是最后状态的NN"""
    X_train = np.hstack((F, S, dF, dS))
    Y_train = model.predict(X_train)
    Y_all = np.hstack((df2_pre,df_pre,f_pre,g3_pre,g4_pre,Y_train))
    
    return test_errors, Y_all




def pred2_f(internal_period, total_epoch, mat_file_path, test_trajectory_numbers):
    """load data"""
    data = scipy.io.loadmat(mat_file_path)
    
    f, df = data["F_save"], data["dFdt_save"]
    s, ds, ds2 = data["S_save"], data["dsdt_save"], data["ds2dt2_save"]
    
    dt, tscale = data["dt00"], data["t_scale"]
    dt = dt/tscale
    num_f = data["Numberof_F"]
    
    j1 = 0
    j2 = int(num_f)
    
    f = f[j1:j2,:]
    df = df[j1:j2,:]
    s = s[j1:j2,:]
    ds = ds[j1:j2,:]
    ds2 = ds2[j1:j2,:]
    
    num_len = len(s)
    t_len = len(s.T)
    
    S = np.reshape(s, (-1, 1)) 
    F = np.reshape(f, (-1, 1))
    dF = np.reshape(df, (-1, 1))
    dS = np.reshape(ds, (-1, 1))
    dS2 = np.reshape(ds2, (-1, 1))
    
    
    f_pre = np.zeros((t_len*num_len,1))
    df_pre = np.zeros((t_len*num_len,1))
    df2_pre = np.zeros((t_len*num_len,1))
    g3_pre = np.zeros((t_len*num_len,1))
    g4_pre = np.zeros((t_len*num_len,1))
    
    model_number = int(total_epoch/internal_period)
    test_errors = np.zeros((model_number + 1, 1))
    
    for n in range(model_number + 1):
        model_number_n = n * internal_period
        model.restore("my2dfreeModel/2dfree_model.ckpt-" + str(model_number_n), verbose=1)

        for j in range(num_len):
            f_pre[j*t_len] = F[j*t_len]
            f_pre[j*t_len+1] = F[j*t_len+1]
            
            for i in range(j*t_len+1,(j+1)*t_len-1):            
                f1 = np.array([f_pre[i]])
                s1 = np.array([S[i]])
                
                df_pre[i] = (f_pre[i]-f_pre[i-1])/dt
                df1 = np.array([df_pre[i]])
                ds1 = np.array([dS[i]])
                
                X_pre = np.hstack((f1,s1,df1,ds1))
                g_pre = model.predict(X_pre)
                
                g3_pre[i] = g_pre[0,0]
                g4_pre[i] = g_pre[0,1]
                df2_pre[i] = g3_pre[i] + dS2[i] * g4_pre[i]
                
                f_pre[i+1] = df2_pre[i] * dt**2 + 2*f_pre[i] - f_pre[i-1]
                
                
        for j in range(num_len):           
            k = int(j*t_len)#start
            f1 = np.array([f_pre[k]])
            s1 = np.array([S[k]])
            df_pre[k] = (f_pre[k+1]-f_pre[k])/dt
            df1 = np.array([df_pre[k]])
            ds1 = np.array([dS[k]])
        
            X_pre = np.hstack((f1,s1,df1,ds1))
            g_pre = model.predict(X_pre)
            
            g3_pre[k] = g_pre[0,0]
            g4_pre[k] = g_pre[0,1]
            df2_pre[k] = g3_pre[k] + dS2[k] * g4_pre[k]
            
            
            k = int((j+1)*t_len - 1)#end
            f1 = np.array([f_pre[k]])
            s1 = np.array([S[k]])
            df_pre[k] = (f_pre[k]-f_pre[k-1])/dt
            df1 = np.array([df_pre[k]])
            ds1 = np.array([dS[k]])
        
            X_pre = np.hstack((f1,s1,df1,ds1))
            g_pre = model.predict(X_pre)
            
            g3_pre[k] = g_pre[0,0]
            g4_pre[k] = g_pre[0,1]
            df2_pre[k] = g3_pre[k] + dS2[k] * g4_pre[k]
            
        """对每个saved NN，计算测试集的MSE"""    
        test_start_j = num_len - test_trajectory_numbers
        f_true1 = F[test_start_j*t_len : num_len*t_len]
        f_pre1 = f_pre[test_start_j*t_len : num_len*t_len]
        
        mse = np.mean((f_true1 - f_pre1) ** 2)
        test_errors[n] = mse
        
    """train process 用的是最后状态的NN"""
    X_train = np.hstack((F, S, dF, dS))
    Y_train = model.predict(X_train)
    Y_all = np.hstack((df2_pre,df_pre,f_pre,g3_pre,g4_pre,Y_train))
 
    return test_errors, Y_all
    
"""model train"""
internal_period = 100
total_epoch = 1500

checkpointer = dde.callbacks.ModelCheckpoint(
     "my2dfreeModel/2dfree_model.ckpt", verbose=1, save_better_only=False,period=internal_period,
 )

losshistory, train_state = model.train(epochs=0, model_save_path = "my2dfreeModel/2dfree_model.ckpt")

losshistory, train_state = model.train(epochs=total_epoch,
                                       display_every=internal_period,
                                       callbacks=[checkpointer],
                                       model_restore_path = "my2dfreeModel/2dfree_model.ckpt-0")


dde.saveplot(losshistory, train_state, issave=True, isplot=False)


"""model test"""
model_number = int(total_epoch/internal_period) + 1# 假设你想要model_number个元素
start_value = 0  # 起始值
ep_number = [start_value + internal_period * n for n in range(model_number)]

"""model test"""
"""long-term prediction用的是最后状态的NN"""
mat_file_path1 = 'dataset/frac2-10tpb-step1.mat'
llsi = 11
test_trajectory_number1 = 4
test_error1, Y_all1 = pred1_f(internal_period, total_epoch, mat_file_path1,llsi,test_trajectory_number1)

mat_file_path2 = 'dataset/frac2-10tpb-mix1.mat'
test_trajectory_number2 = 2
test_error2, Y_all2 = pred2_f(internal_period, total_epoch, mat_file_path2,test_trajectory_number2)

mat_file_path3 = 'dataset/frac2-10tpb-expmix1.mat'
test_trajectory_number3 = 4
test_error3, Y_all3 = pred2_f(internal_period, total_epoch, mat_file_path3,test_trajectory_number3)


array_names0 = ['ep_number','test_error1', 'Y_all1', 'test_error2', 'Y_all2','test_error3', 'Y_all3']
mat_dict0 = {name: array for name, array in zip(array_names0, [ep_number,test_error1, Y_all1, test_error2, Y_all2, 
                                                             test_error3, Y_all3])}
savemat('frac2-10tps-save2d-free.mat', mat_dict0)
savemat('frac2-10tpt-save2d-free.mat', mat_dict0)


############训练好的NN预测s平移的测试集############
s_up_path1 = 'dataset/frac2-12tps-expmix1.mat'
s_up_error1, s_up_Y_all1 = pred2_f(internal_period, total_epoch, s_up_path1, test_trajectory_number3)

s_up_path2 = 'dataset/frac2-14tps-expmix1.mat'
s_up_error2, s_up_Y_all2 = pred2_f(internal_period, total_epoch, s_up_path2, test_trajectory_number3)

s_up_path3 = 'dataset/frac2-16tps-expmix1.mat'
s_up_error3, s_up_Y_all3 = pred2_f(internal_period, total_epoch, s_up_path3, test_trajectory_number3)

s_up_path4 = 'dataset/frac2-18tps-expmix1.mat'
s_up_error4, s_up_Y_all4 = pred2_f(internal_period, total_epoch, s_up_path4, test_trajectory_number3)

s_up_path5 = 'dataset/frac2-20tps-expmix1.mat'
s_up_error5, s_up_Y_all5 = pred2_f(internal_period, total_epoch, s_up_path5, test_trajectory_number3)

s_up_path6 = 'dataset/frac2-22tps-expmix1.mat'
s_up_error6, s_up_Y_all6 = pred2_f(internal_period, total_epoch, s_up_path6, test_trajectory_number3)

s_up_path7 = 'dataset/frac2-24tps-expmix1.mat'
s_up_error7, s_up_Y_all7 = pred2_f(internal_period, total_epoch, s_up_path7, test_trajectory_number3)

array_names1 = ['ep_number','test_error3', 'Y_all3']
s_up_mat_dict1 = {name: array for name, array in zip(array_names1, [ep_number,s_up_error1, s_up_Y_all1])}
s_up_mat_dict2 = {name: array for name, array in zip(array_names1, [ep_number,s_up_error2, s_up_Y_all2])}
s_up_mat_dict3 = {name: array for name, array in zip(array_names1, [ep_number,s_up_error3, s_up_Y_all3])}
s_up_mat_dict4 = {name: array for name, array in zip(array_names1, [ep_number,s_up_error4, s_up_Y_all4])}
s_up_mat_dict5 = {name: array for name, array in zip(array_names1, [ep_number,s_up_error5, s_up_Y_all5])}
s_up_mat_dict6 = {name: array for name, array in zip(array_names1, [ep_number,s_up_error6, s_up_Y_all6])}
s_up_mat_dict7 = {name: array for name, array in zip(array_names1, [ep_number,s_up_error7, s_up_Y_all7])}

savemat('frac2-12tps-save2d-free.mat', s_up_mat_dict1)
savemat('frac2-14tps-save2d-free.mat', s_up_mat_dict2)
savemat('frac2-16tps-save2d-free.mat', s_up_mat_dict3)
savemat('frac2-18tps-save2d-free.mat', s_up_mat_dict4)
savemat('frac2-20tps-save2d-free.mat', s_up_mat_dict5)
savemat('frac2-22tps-save2d-free.mat', s_up_mat_dict6)
savemat('frac2-24tps-save2d-free.mat', s_up_mat_dict7)

############训练好的NN预测s压缩的测试集############
s_left_path1 = 'dataset/frac2-12tpt-expmix1.mat'
s_left_error1, s_left_Y_all1 = pred2_f(internal_period, total_epoch, s_left_path1, test_trajectory_number3)

s_left_path2 = 'dataset/frac2-14tpt-expmix1.mat'
s_left_error2, s_left_Y_all2 = pred2_f(internal_period, total_epoch, s_left_path2, test_trajectory_number3)

s_left_path3 = 'dataset/frac2-16tpt-expmix1.mat'
s_left_error3, s_left_Y_all3 = pred2_f(internal_period, total_epoch, s_left_path3, test_trajectory_number3)

s_left_path4 = 'dataset/frac2-18tpt-expmix1.mat'
s_left_error4, s_left_Y_all4 = pred2_f(internal_period, total_epoch, s_left_path4, test_trajectory_number3)

s_left_path5 = 'dataset/frac2-20tpt-expmix1.mat'
s_left_error5, s_left_Y_all5 = pred2_f(internal_period, total_epoch, s_left_path5, test_trajectory_number3)


#array_names1 = ['ep_number','test_error3', 'Y_all3']
s_left_mat_dict1 = {name: array for name, array in zip(array_names1, [ep_number,s_left_error1, s_left_Y_all1])}
s_left_mat_dict2 = {name: array for name, array in zip(array_names1, [ep_number,s_left_error2, s_left_Y_all2])}
s_left_mat_dict3 = {name: array for name, array in zip(array_names1, [ep_number,s_left_error3, s_left_Y_all3])}
s_left_mat_dict4 = {name: array for name, array in zip(array_names1, [ep_number,s_left_error4, s_left_Y_all4])}
s_left_mat_dict5 = {name: array for name, array in zip(array_names1, [ep_number,s_left_error5, s_left_Y_all5])}


savemat('frac2-12tpt-save2d-free.mat', s_left_mat_dict1)
savemat('frac2-14tpt-save2d-free.mat', s_left_mat_dict2)
savemat('frac2-16tpt-save2d-free.mat', s_left_mat_dict3)
savemat('frac2-18tpt-save2d-free.mat', s_left_mat_dict4)
savemat('frac2-20tpt-save2d-free.mat', s_left_mat_dict5)


s_right_path1 = 'dataset/frac2-22tpt-expmix1.mat'
s_right_error1, s_right_Y_all1 = pred2_f(internal_period, total_epoch, s_right_path1, test_trajectory_number3)

s_right_path2 = 'dataset/frac2-24tpt-expmix1.mat'
s_right_error2, s_right_Y_all2 = pred2_f(internal_period, total_epoch, s_right_path2, test_trajectory_number3)

s_right_path3 = 'dataset/frac2-26tpt-expmix1.mat'
s_right_error3, s_right_Y_all3 = pred2_f(internal_period, total_epoch, s_right_path3, test_trajectory_number3)

s_right_path4 = 'dataset/frac2-28tpt-expmix1.mat'
s_right_error4, s_right_Y_all4 = pred2_f(internal_period, total_epoch, s_right_path4, test_trajectory_number3)

s_right_path5 = 'dataset/frac2-30tpt-expmix1.mat'
s_right_error5, s_right_Y_all5 = pred2_f(internal_period, total_epoch, s_right_path5, test_trajectory_number3)


#array_names1 = ['ep_number','test_error3', 'Y_all3']
s_right_mat_dict1 = {name: array for name, array in zip(array_names1, [ep_number,s_right_error1, s_right_Y_all1])}
s_right_mat_dict2 = {name: array for name, array in zip(array_names1, [ep_number,s_right_error2, s_right_Y_all2])}
s_right_mat_dict3 = {name: array for name, array in zip(array_names1, [ep_number,s_right_error3, s_right_Y_all3])}
s_right_mat_dict4 = {name: array for name, array in zip(array_names1, [ep_number,s_right_error4, s_right_Y_all4])}
s_right_mat_dict5 = {name: array for name, array in zip(array_names1, [ep_number,s_right_error5, s_right_Y_all5])}


savemat('frac2-22tpt-save2d-free.mat', s_right_mat_dict1)
savemat('frac2-24tpt-save2d-free.mat', s_right_mat_dict2)
savemat('frac2-26tpt-save2d-free.mat', s_right_mat_dict3)
savemat('frac2-28tpt-save2d-free.mat', s_right_mat_dict4)
savemat('frac2-30tpt-save2d-free.mat', s_right_mat_dict5)



"""draw figures"""
# Plot MSE err
plt.figure()
plt.subplot(2, 2, 1)
plt.plot(ep_number, np.log10(test_error1))
plt.xlabel("epoch")
plt.ylabel("step MSE err")
plt.show()

plt.subplot(2, 2, 2)
plt.plot(ep_number, np.log10(test_error2))
plt.xlabel("epoch")
plt.ylabel("mix MSE err")
plt.show()

plt.subplot(2, 2, 3)
plt.plot(ep_number, np.log10(test_error3))
plt.xlabel("epoch")
plt.ylabel("expmix MSE err")
plt.show()