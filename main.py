from argparse import ArgumentParser
from config import Config
from gcn_train import load_data,train
from gcn_test import load_data as test_data_load
from gcn_test import test
configs=Config()
parser=ArgumentParser()
parser.add_argument('--mode',type=str,default='train')
parser.add_argument('--name',type=str,default=configs.name)
parser.add_argument('--train_x',type=str,default=configs.train_x_file)
parser.add_argument('--train_y',type=str,default=configs.train_y_file)
parser.add_argument('--valid_x',type=str,default=configs.val_x_file)
parser.add_argument('--valid_y',type=str,default=configs.val_y_file)
parser.add_argument('--test_x',type=str,default=configs.test_x_file)
parser.add_argument('--test_y',type=str,default=configs.test_y_file)
parser.add_argument('--cor_m',type=str,default=configs.corelation_matrix)
parser.add_argument('--model_save',type=str,default=configs.val_model)
parser.add_argument('--model_save_root',type=str,default=configs.model_save_file)
parser.add_argument('--result_save',type=str,default=configs.result_save)
parser.add_argument('--use_valid',type=bool,default=configs.test)
parser.add_argument('--type_list',type=str,default=configs.type_list)
parser.add_argument('--net_model',type=str,default=configs.net_model)
parser.add_argument('--device',type=str,default=configs.device)

parser.add_argument('--epoch',type=int,default=450)
parser.add_argument('--batch_size',type=int,default=32)

args=parser.parse_args()
print(args)

mode=args.mode

if mode=='train':
    configs.train_x_file=args.train_x
    configs.train_y_file=args.train_y
    configs.val_x_file=args.valid_x
    configs.val_y_file=args.valid_y
    configs.test_x_file=args.valid_x
    configs.test_y_file=args.valid_y
    configs.corelation_matrix=args.cor_m
    configs.val_model=args.model_save
    configs.result_save=args.result_save
    configs.test=args.use_valid
    configs.epoch=args.epoch
    configs.batch_size=args.batch_size
    configs.type_list=args.type_list
    configs.net_model=args.net_model
    configs.model_save_file=args.model_save_root
    configs.name=args.name
    configs.device=args.device
else:
    configs.test_x_file=args.valid_x
    configs.test_y_file=args.valid_y
    configs.corelation_matrix=args.cor_m
    configs.val_model=args.model_save
    configs.result_save=args.result_save
    configs.test=args.use_valid
    configs.epoch=args.epoch
    configs.batch_size=args.batch_size
    configs.type_list=args.type_list
    configs.net_model=args.net_model
    configs.model_save_file=args.model_save_root
    configs.name=args.name
    configs.device=args.device

configs.update()

if __name__ == '__main__':
    if mode=='train':
        data_loader,data_loader_test,coef_m,configs,coef_list=load_data(configs)
        train(data_loader,data_loader_test,coef_m,coef_list,configs)
    else:
        print('test mode!')
        data_loader_test,coef_m,configs,coef_list=test_data_load(configs)
        test(data_loader_test,coef_m,coef_list,configs)
