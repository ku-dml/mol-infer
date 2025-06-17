import os, json, sys, argparse, time, datetime, warnings
import torch

# warnings.simplefilter('ignore')

RESULT_filename = "./result_file.txt"

from src import GNN, chemicalgraph, dataset, network, constants

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_sdf", help='input sdf file')
    parser.add_argument("input_value", help='input value file')
    # parser.add_argument("output_pre", help='input output directory for results of preliminary experiments')
    # parser.add_argument("output_eval", help='input output file for results of evaluation experiments')
    parser.add_argument('-e', '--exfv', nargs='+', type=str)
    parser.add_argument('-c', '--clsf', nargs='+', type=str)
    # parser.add_argument('-s', '--standardization', nargs='*', type=str)   
    parser.add_argument('-m', '--metric', nargs='+', type=str)  

    parser.add_argument('-val', '--val', nargs='+', type=str) 

    parser.add_argument('-o', '--output_prefix', nargs='+', type=str)

    parser.add_argument('-r', '--resume', nargs='+', type=str)

    args = parser.parse_args()

    return(args) 

def count_parameters(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)

def main(args):
    sdf_filename = args.input_sdf
    value_filename = args.input_value
    extra_fv_filename = None
    extra_x_length = None

    time_start = time.time()

    if args.exfv is not None:
        extra_fv_filename = args.exfv[0]

    prop_name = "QM9"

    # print(prop_name, args.val)

    # read sdf file
    if extra_fv_filename is not None:
        cg_list = chemicalgraph.construct_from_SDF_file(sdf_filename, extra_fv_used=True, check_qm9=True)
        extra_x_length = cg_list.read_extra_fv_file(extra_fv_filename)
    else:
        cg_list = chemicalgraph.construct_from_SDF_file(sdf_filename, use_rdkit=False, check_qm9=True)

    HAR2EV = 27.2113825435
    KCALMOL2EV = 0.04336414

    conversion = {
        "mu": 1.0, "alpha": 1.0, "homo": HAR2EV,
        "lumo": HAR2EV, "gap": HAR2EV, "r2": 1.0, "zpve": HAR2EV, "u0": HAR2EV,
        "u298": HAR2EV, "h298": HAR2EV, "g298": HAR2EV, "cv": 1.0, "u0_atom": KCALMOL2EV,
        "u298_atom": KCALMOL2EV, "h298_atom": KCALMOL2EV, "g298_atom": KCALMOL2EV
    }

    cg_list.read_value_file_csv(value_filename, args.val, conversion)

    

    # cg_list_without_H = chemicalgraph.create_H_suppresed_chemical_graph(cg_list)

    cg_list.get_fringe_tree()

    cg_dataset = dataset.ChemicalGraph_Datasets(root='./data/')
    cg_dataset.construct_from_ChemicalGraph_list(cg_list, FT_use=True)

    ft_dataset = dataset.FringeTree_Datasets(root='./data')
    ft_dataset.construct_from_FringeTree_list(cg_list.fringe_tree_list)

    num_cg_atom_features = cg_dataset.data.x.shape[1]
    num_cg_edge_features = cg_dataset.data.edge_attr.shape[1]
    num_ft_atom_features = ft_dataset.data.x.shape[1]
    num_ft_edge_features = ft_dataset.data.edge_attr.shape[1]

    print(f"cg_len: {len(cg_dataset)}, num_atom_features: {num_cg_atom_features}, num_edge_features: {num_cg_edge_features}")
    print(f"ft_len: {len(ft_dataset)}, num_atom_features: {num_ft_atom_features}, num_edge_features: {num_ft_edge_features}")
    print(f"Learning start, device={constants.get_default_device()}")

    # SEARCH_PARAMS = {
    #     # 'model': ['GNN_int_no_edge_type'],
    #     'model': ['GNN_2L_no_edge_type'],
    #     'num_atom_features': [num_cg_atom_features],
    #     # 'model': ['GNN_2L'],
    #     'hidden_size': [32],
    #     'num_layers': [1],
    #     # 'ANN_channels': [(16, 1), (32, 1), (64, 1), (16, 16, 1), (32, 16, 1), (64, 32, 1), (64, 32, 16, 1)],
    #     'ANN_channels': [(32, 1)],
    #     'weight_decay': [1e-4],
    #     'metric': ['r2'] # r2, aucroc, bacc
    # }

    # SEARCH_PARAMS = {
    #     'model': ['GNN_2L_no_edge_type'],
    #     'num_atom_features': [num_cg_atom_features],
    #     'hidden_size': [16, 32],
    #     'num_layers': [1,2,3],
    #     'ANN_channels': [(32, 7), (32,16,7)],
    #     'weight_decay': [1e-4],
    #     'metric': ['r2'] # r2, aucroc, bacc
    # }

    if args.metric is not None:
        SEARCH_PARAMS['metric'] = args.metric

    if args.clsf is not None:
        SEARCH_PARAMS['clsf'] = [n_class]

    if args.resume is not None:
        resume_filename = args.resume[0]
    else:
        resume_filename = None

    # train_size = int(len(cg_dataset) * 0.8)
    # val_size = int(len(cg_dataset) * 0.1)
    # test_size = len(cg_dataset) - train_size - val_size

    train_size = 110000
    val_size = 10000
    test_size = len(cg_dataset) - train_size - val_size

    generator = torch.Generator().manual_seed(42)

    train_dataset, val_dataset, test_dataset = torch.utils.data.random_split(cg_dataset, [train_size, val_size, test_size], generator=generator)

    # best_params, best_train_score, best_test_score, best_train_score_mae, best_test_score_mae, best_train_score_all, best_test_score_all, best_train_score_mae_all, best_test_score_mae_all = GNN.hyper_param_tuning_trainval(train_dataset, val_dataset, ft_dataset, extra_x_length, search_params=SEARCH_PARAMS)

    # print(f"best_params: {json.dumps(best_params)}")
    # print(f"best_train_score: r2:{best_train_score}, mae:{best_train_score_mae}")
    # print(f"best_test_score: r2:{best_test_score}, mae:{best_test_score_mae}")
    # print(f"best_train_score: r2_all:{best_train_score_all}")
    # print(f"best_test_score: r2_all:{best_test_score_all}")
    # print(f"best_train_score: mae_all:{best_train_score_mae_all}")
    # print(f"best_test_score: mae_all:{best_test_score_mae_all}")


    output_prefix = args.output_prefix[0]

    best_params = {
            'model': 'GNN_2L_no_edge_type',
            'num_atom_features': num_cg_atom_features,
            'hidden_size': 16,
            'p_size': 32,
            'num_layers': 3,
            'ANN_channels': (16,len(args.val)),
            'weight_decay': 0,
            'metric': 'r2' # r2, aucroc, bacc
        }

    best_val_loss, test_loss_mae, test_loss_r2 = \
        GNN.eval_GNN(train_dataset, val_dataset, test_dataset, ft_dataset, extra_x_length, best_params, num_epoch=300, output_prefix=output_prefix, resume_filename=resume_filename)

    # for name, param in model.named_parameters():
    #     if param.requires_grad:
    #         print(name, param.data)     


if __name__ == "__main__":
    main(get_args())





