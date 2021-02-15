'''
cite_pbmc_16000
./my_data/train_cite_cbmc_split
./my_data/train_cite_cbmc_y

cite init
trainx2.npy
valy2.npy

cite_dev_proteins
./my_data/train_cbmc_x.npy
./my_data/dev_cbmc_x.npy

cbmc
(4194, 20501)
(4194, 13)
(1797, 20501)
(1797, 13)

pbmc
(3897, 17114)
(3897, 10)
(1670, 17114)
(1670, 10)

'''
import os
class Config():
    def __init__(self):
        self.name='cbmc'
        self.net_model='k'
        self.type_list = 'combine_database_text_experiment_coe'
        # type_list=['coe']

        self.device='cuda:0'
        self.train_x_file=''#'./my_data/train_'+self.name+'_x_raw.npy'
        self.train_y_file=''#'./my_data/train_'+self.name+'_y.npy'
        self.val_x_file=''#'./my_data/dev_'+self.name+'_x_raw.npy'
        self.val_y_file=''#'./my_data/dev_'+self.name+'_y.npy'
        self.test_x_file=''#'./my_data/dev_'+self.name+'_x_raw.npy'
        self.test_y_file=''#'./my_data/dev_'+self.name+'_y.npy'
        self.corelation_list=[]

        self.corelation_matrix=''# './prior/'+self.name+'.npy'
        self.test=True

        self.model_save_file='./model_saver_v2/raw/'
        self.model_save_root=self.model_save_file+'dive_'+self.name+'/'
        if not os.path.exists(self.model_save_root):
            os.makedirs(self.model_save_root)
        self.val_model=self.model_save_root+'val_net.pt' # used
        self.peason_model=self.model_save_root+'net.pt'
        self.logs=self.model_save_root+'logs.txt'
        self.result_save=''

        self.input_dim=0
        self.output_dim=0
        self.epoch=450
        self.batch_size=32

    def update(self):
        self.type_list_t=self.type_list.split('_')
        self.train_x_file='./my_data/train_'+self.name+'_x.npy'
        self.train_y_file='./my_data/train_'+self.name+'_y.npy'
        self.val_x_file='./my_data/dev_'+self.name+'_x.npy'
        self.val_y_file='./my_data/dev_'+self.name+'_y.npy'
        self.test_x_file='./my_data/dev_'+self.name+'_x.npy'
        self.test_y_file='./my_data/dev_'+self.name+'_y.npy'
        self.corelation_matrix='./prior/'+self.name+'.npy'



        for t in self.type_list_t:
            self.corelation_list.append('./prior/'+t+'_'+self.name+'.npy')
        self.model_save_root = self.model_save_file + 'dive_' + self.name + '/'
        if not os.path.exists(self.model_save_root):
            os.makedirs(self.model_save_root)
        if self.net_model=='k':
            self.result_save=self.model_save_root+'cite_t_'
            for t in self.type_list_t:
                self.result_save+=t
            self.result_save+='.npy'
        else:
            self.result_save = self.model_save_root + 'cite_' + self.name + '.npy'

        self.val_model=self.model_save_root+self.type_list+'.pt' # used
        self.peason_model=self.model_save_root+self.type_list+'net_t.pt'
        self.logs=self.model_save_root+'logs.txt'



