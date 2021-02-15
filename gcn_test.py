import torch
import numpy as np
from torchvision.datasets import mnist
from torch import nn
from torch.autograd import Variable
import matplotlib.pyplot as plt
import torch.nn.functional as F
from torch.utils.data import DataLoader
from torch.utils.data import TensorDataset,DataLoader
from tqdm import tqdm
from scipy.stats import pearsonr
from model import Net,Net2
from sklearn.model_selection import train_test_split
from config import Config
def load_data(configs,k=True):
    print('loading data....')
    x_test = np.load(configs.test_x_file)
    y_test = np.load(configs.test_y_file)
    configs.input_dim=x_test.shape[1]
    configs.output_dim=y_test.shape[1]

    try:
        print('using knowledge from',configs.corelation_matrix)
        coef_m = np.load(configs.corelation_matrix)
        coef_list=[]
        for c in configs.corelation_list:
            coef_list.append(np.load(c))
        coef_list=np.array(coef_list)
        print(coef_list.shape)
    except:
        print('random init')
        coef_m=np.random.random((y_test.shape[1],y_test.shape[1]))
    print(x_test.shape,y_test.shape)
    print(x_test.shape,y_test.shape)

    x_test = torch.from_numpy(x_test)
    y_test = torch.from_numpy(y_test)
    datas_test = TensorDataset(x_test, y_test)
    data_loader_test = DataLoader(datas_test, batch_size=configs.batch_size, shuffle=False)
    if k:
        return data_loader_test, coef_m, configs,coef_list
    else:
        return data_loader_test,coef_m,configs

def test(test_loader,coef_m,coef_list,configs):
    label_shape=configs.output_dim
    if configs.net_model=='k':
        net = Net2(configs.input_dim, configs.output_dim,configs.device,coef_list)
    else:
        net = Net(configs.input_dim, configs.output_dim, configs.device, coef_m)
    net.load_state_dict(torch.load(configs.val_model))
    net.to(configs.device)
    print(net)
    criterion = nn.MSELoss()
    net.eval()
    output=[[] for _ in range(label_shape)]
    targets=[[] for _ in range(label_shape)]
    val_loss=0
    for data,target in test_loader:
        data, target = Variable(data).float(), Variable(target).float()
        # data.cuda()
        data = data.to(configs.device)
        target = target.to(configs.device)
        out = net(data)
        loss = criterion(out, target)
        val_loss+=loss.item()
        for i in range(label_shape):
            output[i].append(out[:,i])
            targets[i].append(target[:,i])
    val_loss/=len(test_loader)
    pearson_list=[]
    for i in range(label_shape):
        output[i]=torch.cat([x for x in output[i]]).cpu().detach().numpy().tolist()
        targets[i] = torch.cat([x for x in targets[i]]).cpu().detach().numpy().tolist()
        p=pearsonr(output[i],targets[i])[0]
        if p>=-1 and p<=1:
            pearson_list.append(p)
        else:
            pearson_list.append(0)
            print('0 pers at',i,max(output[i]),min(output[i]))

    pear_sum=np.mean(np.array(pearson_list))
    print('test_loss', val_loss, '\tpearsonr', pear_sum)



    with open(configs.logs,'a') as f:
        f.write('test\t')
        f.write(configs.net_model+'\t')
        f.write(str(configs.test)+'\t')
        f.write(configs.name+'\t')
        f.write(configs.result_save+'\n')
        f.write(str(val_loss)+'\t'+str(pear_sum)+'\n')

    print(configs.net_model,configs.test,configs.name,configs.result_save)