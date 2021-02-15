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
    print('using valid is',configs.test)
    x_train=np.load(configs.train_x_file)
    y_train=np.load(configs.train_y_file)
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
        coef_m=np.random.random((y_train.shape[1],y_train.shape[1]))
    print(x_test.shape,y_test.shape)
    print(x_test.shape,y_test.shape)

    x = torch.from_numpy(x_train)
    y = torch.from_numpy(y_train)
    x_test = torch.from_numpy(x_test)
    y_test = torch.from_numpy(y_test)
    datas_test = TensorDataset(x_test, y_test)
    data_loader_test = DataLoader(datas_test, batch_size=configs.batch_size, shuffle=False)
    datas = TensorDataset(x, y)
    data_loader = DataLoader(datas, batch_size=configs.batch_size, shuffle=True)
    if k:
        return data_loader, data_loader_test, coef_m, configs,coef_list
    else:
        return data_loader,data_loader_test,coef_m,configs


def train(train_loader,test_loader,coef_m,coef_list,configs):
    min_val_loss = 1000
    min_train_loss = 1000
    min_epoch = 0
    label_shape=configs.output_dim
    if configs.net_model=='k':
        net = Net2(configs.input_dim, configs.output_dim,configs.device,coef_list)
    else:
        net = Net(configs.input_dim, configs.output_dim, configs.device, coef_m)
    net.to(configs.device)
    print(net)
    criterion = nn.MSELoss()
    # 随机梯度下降
    optimizer = torch.optim.Adam(net.parameters(), lr=10e-6)
    pear_sum_max=0
    loss_epoch=0
    for e in range(configs.epoch):
        print(e,end='\t')
        train_loss = 0
        val_loss=0
        net.train()
        optimizer.zero_grad()
        for data,target in tqdm(train_loader,total=len(train_loader)):
            data, target = Variable(data).float(), Variable(target).float()
            # data.cuda()
            data=data.to(configs.device)
            target=target.to(configs.device)
            out=net(data)

            loss=criterion(out,target)
            loss.backward()
            optimizer.step()
            train_loss+=loss.item()
        train_loss/=len(train_loader)
        print('train_loss',train_loss,end='\t')

        # if configs.test:
        net.eval()
        output=[[] for _ in range(label_shape)]
        targets=[[] for _ in range(label_shape)]
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
        print('val_loss', val_loss, '\tpearsonr', pear_sum)


        print(train_loss,min_train_loss)
        if (configs.test and (val_loss<min_val_loss)) or(configs.test==False and (train_loss<min_train_loss)):
            if (configs.test and (val_loss<min_val_loss)):
                print('validation is smaller')
            min_val_loss=val_loss
            min_train_loss=train_loss
            loss_epoch=e
            torch.save(net.state_dict(),configs.val_model)
            min_epoch=e
            patient_count=0
            test_loss=0
            net.eval()
            output = [[] for _ in range(label_shape)]
            targets = [[] for _ in range(label_shape)]
            for data, target in test_loader:
                data, target = Variable(data).float(), Variable(target).float()
                data = data.to(configs.device)
                target = target.to(configs.device)
                out = net(data)
                loss = criterion(out, target)
                test_loss+= loss.item()
                for i in range(label_shape):
                    output[i].append(out[:, i])
                    targets[i].append(target[:, i])
            test_pearson_list = []
            for i in range(label_shape):
                output[i] = torch.cat([x for x in output[i]]).cpu().detach().numpy().tolist()
                targets[i] = torch.cat([x for x in targets[i]]).cpu().detach().numpy().tolist()
                p = pearsonr(output[i], targets[i])[0]
                if p >= -1 and p <= 1:
                    test_pearson_list.append(p)
                else:
                    test_pearson_list.append(0)
            test_pear_sum = np.mean(np.array(test_pearson_list))
            test_loss/=len(test_loader)
            print('\ntest',test_pear_sum,test_loss,e)
            print(test_pearson_list)
            output = np.array(output)
            print(output.shape)
            path = configs.result_save
            np.save(path, output)
            print('saved! at ' + path)
        print(pear_sum_max,min_epoch)
    if configs.test:
        print(test_pear_sum,test_loss)
    print(min_val_loss,loss_epoch)
    with open(configs.logs,'a') as f:
        f.write(configs.net_model+'\t')
        f.write(str(configs.test)+'\t')
        f.write(configs.name+'\t')
        f.write(configs.result_save+'\n')
        if configs.test:
            f.write(str(test_pear_sum)+'\t'+str(test_loss)+'\n')
        else:
            f.write(str(pear_sum_max)+'\t'+str(min_val_loss)+'\n')
    print(configs.net_model,configs.test,configs.name,configs.result_save)

if __name__ == '__main__':
    configs=Config()
    configs.update()
    data_loader,data_loader_test,coef_m,configs,coef_list=load_data(configs)
    train(data_loader,data_loader_test,coef_m,coef_list,configs)


