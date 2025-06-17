import os, sys, datetime

import torch
import torch.nn.functional as F
from torch.nn.utils import clip_grad_norm_
from torch_geometric.loader import DataLoader
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score, roc_auc_score, balanced_accuracy_score, mean_absolute_error
from sklearn.utils import class_weight
from scipy.special import softmax
from warmup_scheduler import GradualWarmupScheduler

import numpy as np
import random

from src import chemicalgraph, dataset, network, constants

import itertools, json, time

RANDOM_SEED = 42
RANDOM_SEED_BASE = 10000

def set_seed(seed):
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)

def count_parameters(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)

def test(model, loader, ema, device, ft_dataset_used, ft_dataset_batch):
    mae = 0
    ema.assign(model)
    for data in loader:
        data = data.to(device)
        if ft_dataset_used:
            output = model(data, ft_dataset_batch)
        else:
            output = model(data)
        mae += (output - data.y).abs().sum().item()
        y_true = data.y.cpu().detach().numpy()
        y_pred = output.cpu().detach().numpy()
    ema.resume(model)
    model_mae = mae / len(loader.dataset)
    model_r2 = r2_score(y_true, y_pred)

    return model_mae, model_r2, y_true, y_pred

def learn_GNN_selfmade_trainvaltest(
    train_dataset, val_dataset, test_dataset, ft_dataset, params, 
    extra_x_length=None, batch_size=32, num_epoch=5000,
    random_seed=RANDOM_SEED, learning_rate=1e-4, verbose=True, save=False, 
    log_prefix=f"./log/{datetime.datetime.today()}", model_name_save="GNN"
):
    set_seed(random_seed)

    model_name = params['model'] if 'model' in params.keys() else 'GNN_selfmade'
    hidden_size = params['hidden_size']
    num_layers = params['num_layers']
    ANN_channels = params['ANN_channels']
    weight_decay = params['weight_decay']

    num_atom_features = params['num_atom_features']

    # metric_name = params['metric']

    # if metric_name == 'r2':
    metric = r2_score

    mae_train_summary = []
    mae_val_summary = []
    timestamp_summary = []
    best_val_loss = None

    # train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    train_loader = DataLoader(train_dataset, batch_size=batch_size)
    val_loader = DataLoader(val_dataset, batch_size=len(val_dataset))
    test_loader = DataLoader(test_dataset, batch_size=len(test_dataset))

    # train_all_loader = DataLoader(train_dataset, batch_size=len(train_dataset))

    ft_dataloader = DataLoader(ft_dataset, batch_size=len(ft_dataset))
    for ft_dataset_batch in ft_dataloader:
        break

    ft_dataset_used = False

    if model_name == "GNN_int":
        model = network.GNN_int(
            in_channels=num_atom_features,
            hidden_channels=hidden_size,
            num_layers=num_layers,
            p_channels=hidden_size,
            ANN_channels=ANN_channels,
            num_edge_types=3,
            extra_x_length=extra_x_length
        )
    elif model_name == "GNN_int_no_edge_type":
        model = network.GNN_int_no_edge_type(
            in_channels=num_atom_features,
            hidden_channels=hidden_size,
            num_layers=num_layers,
            p_channels=hidden_size,
            ANN_channels=ANN_channels,
            extra_x_length=extra_x_length
        )
    elif model_name == "GNN_2L":
        model = network.GNN_2L(
            in_channels_FT=ft_dataset.num_node_features,
            hidden_channels_FT=64,
            num_layers_FT=6,
            p_channels_FT=128,
            ANN_channels_FT=(128, 64, 32, 16, 8),
            in_channels_IN=num_atom_features,
            hidden_channels_IN=hidden_size,
            num_layers_IN=num_layers,
            p_channels_IN=hidden_size,
            ANN_channels_IN=ANN_channels,
            num_edge_types=3,
            extra_x_length=extra_x_length,
            dropout=0.2
        )
        ft_dataset_used = True
    elif model_name == "GNN_2L_no_edge_type":
        model = network.GNN_2L_no_edge_type(
            in_channels_FT=ft_dataset.num_node_features,
            hidden_channels_FT=256,
            num_layers_FT=10,
            p_channels_FT=256,
            ANN_channels_FT=(256, 256, 128, 64, 32, 16, 8),
            in_channels_IN=num_atom_features,
            hidden_channels_IN=hidden_size,
            num_layers_IN=num_layers,
            p_channels_IN=hidden_size,
            ANN_channels_IN=ANN_channels,
            extra_x_length=extra_x_length,
            dropout=0.2
        )
        ft_dataset_used = True
    # model = network.GraphSAGE(
    #     in_channels=train_dataset.num_node_features,
    #     hidden_channels=hidden_size,
    #     num_layers=num_layers,
    #     out_channels=1
    # )

    # model = network.GraphATConv(
    #     in_channels=train_dataset.num_node_features,
    #     hidden_channels=hidden_size,
    #     num_layers=num_layers,
    #     out_channels=1
    # )

    device = torch.device(constants.get_default_device())

    model.to(device)
    ft_dataset_batch.to(device)
    print("Number of model parameters: ", count_parameters(model))

    ema = network.EMA(model, decay=0.999)

    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, weight_decay=weight_decay, amsgrad=False)
    scheduler = torch.optim.lr_scheduler.ExponentialLR(optimizer, gamma=0.9961697)
    scheduler_warmup = GradualWarmupScheduler(optimizer, multiplier=1.0, total_epoch=1, after_scheduler=scheduler)

    if save is True:
        output_model_name = f"./checkpoint_model_output_{model_name_save}.pth"
        print(f"The model will be saved as {output_model_name}")
    else:
        output_model_name = None
    early_stopping = network.EarlyStopping(patience=500, verbose=False, path=output_model_name)

    print(f"Epoch, loss_train, score_train, score_test, time")

    time_epoch_start = time.time()

    for epoch in range(num_epoch):
        running_loss = 0
        model.train()

        for index, batch in enumerate(train_loader):
            batch = batch.to(device)
            optimizer.zero_grad()

            if ft_dataset_used:
                out = model(batch, ft_dataset_batch)
            else:
                out = model(batch)

            loss = F.l1_loss(out, batch.y)
            running_loss += loss.item() * len(batch)
            loss.backward()
            clip_grad_norm_(model.parameters(), max_norm=1000, norm_type=2)
            optimizer.step()

            curr_epoch = epoch + float(index) / (len(train_dataset) / batch_size)
            scheduler_warmup.step(curr_epoch)

            ema(model)

        running_loss_mae = running_loss / len(train_dataset)
        val_loss_mae, val_loss_r2, val_true, val_pred = test(model, val_loader, ema, device, ft_dataset_used, ft_dataset_batch)
        # train_loss_mae, train_loss_r2 = test(model, train_all_loader, ema, device, ft_dataset_used, ft_dataset_batch)
        mae_train_summary.append(running_loss_mae)
        mae_val_summary.append(val_loss_mae)
        timestamp_summary.append(epoch + 1)

        if best_val_loss is None or val_loss_mae <= best_val_loss:
            test_loss_mae, test_loss_r2, test_true, test_pred = test(model, test_loader, ema, device, ft_dataset_used, ft_dataset_batch)
            best_val_loss = val_loss_mae
            torch.save(model.state_dict(), output_model_name)

            print(f'Epoch: {epoch+1:04d}, Train MAE: {running_loss_mae:.7f}, Val MAE: {val_loss_mae:.7f}, Val R2: {val_loss_r2:.7f}, Test MAE: {test_loss_mae:.7f}, Val R2: {test_loss_r2:.7f}, Time: {time.time() - time_epoch_start:.7f}.')

    # print(f"{metric_name}_train: {r2_train_score_tmp}\n{metric_name}_test: {r2_test_score_tmp}\n")

    with open(f"{log_prefix}_timestamp.txt", 'a') as f:
        f.write(f"Epoch, mae_train, mae_test\n")
        for _epoch, _mae_train, _mae_test in zip(timestamp_summary, mae_train_summary, mae_val_summary):
            f.write(f"{_epoch}, {_mae_train}, {_mae_test}\n")

        f.write("\n")
        f.write(f"Val:\n")
        f.write(f"val_true, val_pred\n")
        for (v_t, v_p) in zip(list(val_true), list(val_pred)):
            f.write(f"{v_t}, {v_p}\n")
        f.write("\n")
        f.write(f"Test:\n")
        f.write(f"test_true, test_pred\n")
        for (t_t, t_p) in zip(list(test_true), list(test_pred)):
            f.write(f"{t_t}, {t_p}\n")


    print(f"Log file saved in {log_prefix}_timestamp.txt!!!")


    print('Best Validation MAE:', best_val_loss)
    print('Testing MAE:', test_loss_mae)
    print('Testing R2:', test_loss_r2)

    return best_val_loss, test_loss_mae, test_loss_r2

# def hyper_param_tuning_trainval(train_dataset, val_dataset, ft_dataset, extra_x_length, search_params):
#     best_train_score = None
#     best_test_score = None
#     best_train_score_mae = None
#     best_test_score_mae = None
#     best_train_score_all = None
#     best_test_score_all = None
#     best_train_score_mae_all = None
#     best_test_score_mae_all = None
#     best_params = None

#     keys = list(search_params)
#     for values in itertools.product(*map(search_params.get, keys)):
#         params = dict(zip(keys, values))

#         print(f"params: {json.dumps(params)}")

#         time_start_org = time.time()

#         # for split_seed in range(1, Times_pre+1):
#         r2_train_score, r2_test_score, mae_train_score, mae_test_score, _, r2_train_score_all, r2_test_score_all, mae_train_score_all, mae_test_score_all = learn_GNN_selfmade_trainval(train_dataset, val_dataset, ft_dataset, params, 
#             extra_x_length=extra_x_length, random_seed=RANDOM_SEED_BASE)

#         print(f"train score: r2: {r2_train_score}\ttest score: r2: {r2_test_score}")
#         print(f"train score: r2: {r2_train_score_all}\ntest score: r2: {r2_test_score_all}")
#         print(f"train score: mae: {mae_train_score}\ttest score: r2: {mae_test_score}")
#         print(f"train score: mae: {mae_train_score_all}\ntest score: r2: {mae_test_score_all}")
#         print(f"time: {time.time() - time_start_org}\n")

#         if best_test_score is None or mae_test_score > best_test_score_mae:
#             best_train_score = r2_train_score
#             best_test_score = r2_test_score
#             best_train_score_mae = mae_train_score
#             best_test_score_mae = mae_test_score
#             best_train_score_all = r2_train_score_all
#             best_test_score_all = r2_test_score_all
#             best_train_score_mae_all = mae_train_score_all
#             best_test_score_mae_all = mae_test_score_all

#             best_params = params.copy()

#     return best_params, best_train_score, best_test_score, best_train_score_mae, best_test_score_mae, best_train_score_all, best_test_score_all, best_train_score_mae_all, best_test_score_mae_all

def eval_GNN(train_dataset, val_dataset, test_dataset, ft_dataset, extra_x_length, best_params, num_epoch=5000, model_name_save="GNN"):
    # train_val_dataset = torch.utils.data.ConcatDataset([train_dataset, val_dataset])

    print(f"\nEVALUATION:")
    print(f"params: {json.dumps(best_params)}")

    time_start_org = time.time()

    best_val_loss, test_loss_mae, test_loss_r2 = learn_GNN_selfmade_trainvaltest(train_dataset, val_dataset, test_dataset, ft_dataset, best_params, 
        extra_x_length=extra_x_length, random_seed=RANDOM_SEED, save=True, num_epoch=num_epoch, model_name_save=model_name_save)

    print(f"time: {time.time() - time_start_org}\n")

    return model, best_val_loss, test_loss_mae, test_loss_r2

