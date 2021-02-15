import torch.nn as nn
import torch.nn.functional as F
import torch
from torch.autograd import Variable

class Net(nn.Module):
    def __init__(self,input_dim,output_dim,device,coef_m=None):
        super(Net, self).__init__()
        self.output_dim=output_dim
        self.fc_input=nn.Linear(input_dim,1024)
        self.fc1=nn.Linear(1024,128)
        self.dropout=nn.Dropout(0.3)
        self.output_list=[]
        self.device=device
        self.output_list2=[]

        for i in range(output_dim):
            self.output_list.append(nn.Linear(128,64))
            self.output_list2.append(nn.Linear(64,1))
        self.output_list=nn.ModuleList(self.output_list)
        self.output_list2=nn.ModuleList(self.output_list2)
        # self.fc2=nn.Linear(output_dim,output_dim)

        self.trans_gcn=nn.Linear(64,64)
        self.nodes_trans_nn = nn.Linear(output_dim, output_dim, bias=False)
        self.nodes_trans=Variable(torch.from_numpy(coef_m).contiguous().float(),requires_grad=True).to(device)
        coef_w=torch.from_numpy(coef_m*0.25).contiguous().float()
        self.nodes_trans_nn.weight=torch.nn.Parameter(coef_w)

        self.output_mapping=nn.Linear(64,1)
        self.prelu=nn.PReLU()
    def forward(self,x):
        x.to(self.device)
        x=F.relu(self.fc_input(x))
        x=F.relu(self.fc1(x))
        hidden_list=[]
        output_list=[]
        for i in range(self.output_dim):
            hidden_list.append(F.relu(self.output_list[i](x)))
            hidden_list[i] = hidden_list[i].unsqueeze(1)
        hiddens = torch.cat([h for h in hidden_list], 1)
        nodes_hidden = hiddens.permute(0, 2, 1)
        nodes_trans = self.nodes_trans_nn(nodes_hidden)
        nodes_trans = nodes_trans.permute(0, 2, 1)
        # print(nodes_trans.shape)
        trans = F.sigmoid(self.trans_gcn(nodes_trans))

        for i in range(self.output_dim):
            output_list.append(F.elu(self.output_list2[i](trans[:,i:i+1,:].squeeze(1))))
        outputs=self.prelu(self.output_mapping(trans)).squeeze(-1)
        return outputs

class Net2(nn.Module):
    def __init__(self,input_dim,output_dim,device,coef_m=None):
        super(Net2, self).__init__()
        self.output_dim = output_dim
        self.device = device
        self.fc_input = nn.Linear(input_dim, 1024)
        self.fc1 = nn.Linear(1024, 128)
        self.dropout = nn.Dropout(0.3)
        self.output_list = []
        self.device = device
        self.output_list2 = []
        self.coef_m = Variable(torch.from_numpy(coef_m).contiguous().float(), requires_grad=True).to(self.device)
        self.coef_m = self.coef_m.permute(1, 2, 0)

        self.output_list = nn.ModuleList(self.output_list)
        self.output_list2 = nn.ModuleList(self.output_list2)
        # self.fc2=nn.Linear(output_dim,output_dim)

        self.feature_nums = self.coef_m.shape[2]
        self.coef_nn_list = []
        self.coef_nn_list_att = []
        self.feature_hidden = 32
        for i in range(self.feature_nums):
            self.coef_nn_list.append(nn.Linear(1, self.feature_hidden))
            self.coef_nn_list_att.append(nn.Linear(self.feature_hidden, self.feature_hidden, bias=False))
        self.coef_nn_list = nn.ModuleList(self.coef_nn_list)
        self.coef_nn_list_att = nn.ModuleList(self.coef_nn_list_att)

        self.output_list_att = []
        for i in range(self.output_dim):
            self.output_list_att.append(nn.Linear(self.feature_hidden, self.feature_hidden, bias=False))
        self.output_list_att = nn.ModuleList(self.output_list_att)

        self.trans_gcn = nn.Linear(64 + self.feature_hidden * self.output_dim,
                                   64 + self.feature_hidden * self.output_dim,bias=False)
        self.nodes_trans_nn = nn.Linear(output_dim, output_dim,bias=False)

        self.output_mapping = nn.Linear(64 + self.feature_hidden * self.output_dim, 1)
        self.prelu = nn.PReLU()

        # self.coef_feature=nn.Linear(self.feature_nums*self.feature_hidden,self.feature_nums*self.feature_hidden)

        # self.coef_feature=nn.Linear(self.feature_nums*self.feature_hidden*self.output_dim,self.feature_nums*self.feature_hidden)
        self.nodes_feature = nn.Linear(64 + self.feature_hidden * self.output_dim, 64)

        for i in range(output_dim):
            self.output_list.append(nn.Linear(128, 64))
            self.output_list2.append(nn.Linear(64 + self.feature_hidden * self.output_dim, 1))

    def forward(self, x):
        batch_size = x.shape[0]
        x.to(self.device)
        x = F.relu(self.fc_input(x))
        x = F.relu(self.fc1(x))
        hidden_list = []
        output_list = []
        for i in range(self.output_dim):
            hidden_list.append(F.relu(self.output_list[i](x)))
            hidden_list[i] = hidden_list[i].unsqueeze(1)
        coef_feature_list = []
        for i in range(self.feature_nums):
            coef_feature_list.append(F.elu(self.coef_nn_list[i](self.coef_m[:, :, i].unsqueeze(-1))))
        coef_feature = torch.cat([h for h in coef_feature_list], 2)
        sum = torch.sum(torch.exp(F.elu(coef_feature)), 2)
        sum = sum.unsqueeze(2)

        hiddens = torch.cat([h for h in hidden_list], 1)

        coef_feature_list_att = []
        coef_feature_s = torch.zeros((self.output_dim, self.output_dim, self.feature_hidden)).contiguous().float().to(
            self.device)
        for idx, f in enumerate(coef_feature_list):
            att = torch.exp(F.elu(f))
            p = F.elu(att * self.coef_nn_list_att[idx](f) / sum)
            coef_feature_list_att.append(p)
            coef_feature_s = coef_feature_s + p
        coef_feature_s /= self.feature_nums
        coef_feature_s = F.elu(coef_feature_s)

        # coef_feature = torch.cat([h for h in coef_feature_list_att], 2)

        nodes_feature_s = torch.zeros((self.output_dim, self.feature_hidden)).contiguous().float().to(self.device)
        coef_sum = torch.sum(torch.exp(F.elu(coef_feature_s)), 1)
        coef_feature_att = []
        for i in range(self.output_dim):
            nodes_feat = coef_feature_s[:, i, :]
            p = torch.exp(F.elu(nodes_feat)) / coef_sum
            att = p * self.output_list_att[i](nodes_feat)
            coef_feature_att.append(att)
        coef_feature_att = torch.cat([h for h in coef_feature_att], 1)
        coef_feature_att = coef_feature_att.unsqueeze(0)
        coef_feature_att = coef_feature_att.repeat(batch_size, 1, 1)
        hiddens = torch.cat([hiddens, coef_feature_att], 2)


        # graph neural network
        nodes_hidden = hiddens.permute(0, 2, 1)
        nodes_trans = self.nodes_trans_nn(nodes_hidden)
        nodes_trans = nodes_trans.permute(0, 2, 1)
        trans = F.sigmoid(self.trans_gcn(nodes_trans))

        for i in range(self.output_dim):
            output_list.append(self.prelu(self.output_list2[i](trans[:, i:i + 1, :].squeeze(1))))
        outputs = self.prelu(self.output_mapping(trans)).squeeze(-1)
        return outputs