import torch
import torch_sparse
from torch_geometric.nn.conv import MessagePassing
from torch_geometric.nn.models import GAT, GraphSAGE
import torch_geometric
import torch.nn.functional as F
# from torch_geometric.nn.dense.linear import Linear
from torch.nn import Linear, Softmax, Parameter

import numpy as np

import os, copy

from src import constants

# os.environ['PYTORCH_ENABLE_MPS_FALLBACK'] = '1'

# set manual seed for torch in intializing weights
RANDOM_SEED_WEIGHT = 101
torch.manual_seed(RANDOM_SEED_WEIGHT)

class Score:
    # A class to record the learning score
    # Here the score is assumed to be like: higher -> better
    def __init__(self, **kwargs):
        self.score_dict = copy.deepcopy(kwargs)

    def __gt__(self, other):
        for k, v in self.score_dict.items():
            if v < other.score_dict[k]:
                return False
        return True

    def __str__(self):
        return str(self.score_dict)

    def __repr__(self):
        return repr(self.score_dict)

    def update(self, other):
        for k in self.score_dict.keys():
            self.score_dict[k] = max(self.score_dict[k], other.score_dict[k])

class EarlyStopping:
    def __init__(self, patience=100, verbose=True, path="checkpoint_model.pth"):
        self.patience = patience
        self.verbose = verbose
        self.counter = 0
        self.best_score = None
        self.early_stop = False
        self.val_loss_min = np.Inf
        self.path = path
        self.to_output = False
        # self.log_filename = log_filename

    def __call__(self, val_loss, model, epoch_index):
        score = val_loss
        if self.best_score is None:
            self.best_score = score
            self.checkpoint(val_loss, model, epoch_index)
            self.to_output = True
        elif self.best_score > score:
            self.counter += 1
            if self.counter >= self.patience:
                self.early_stop = True
            self.to_output = False
        else:
            # self.best_score = score
            self.best_score.update(score)
            self.checkpoint(val_loss, model, epoch_index)
            self.counter = 0
            self.to_output = True

    def checkpoint(self, val_loss, model, epoch_index):
        if self.verbose:
            print(f'Epoch {epoch_index}: Validation loss decreased ({self.val_loss_min:.6f} --> {val_loss:.6f}).  Saving model ...')

        if self.path is not None:
            torch.save(model.state_dict(), self.path)
        self.val_loss_min = val_loss



# Implementation for computing exponential moving averages (EMA) of model parameters
# copy from https://github.com/XieResearchGroup/Physics-aware-Multiplex-GNN/blob/main/utils/ema.py

class EMA:
    def __init__(self, model, decay):
        self.decay = decay
        self.shadow = {}
        self.original = {}

        for name, param in model.named_parameters():
            if param.requires_grad:
                self.shadow[name] = param.data.clone()

    def __call__(self, model, num_updates=99999):
        decay = min(self.decay, (1.0 + num_updates) / (10.0 + num_updates))
        for name, param in model.named_parameters():
            if param.requires_grad:
                assert name in self.shadow
                new_average = \
                    (1.0 - decay) * param.data + decay * self.shadow[name]
                self.shadow[name] = new_average.clone()

    def assign(self, model):
        for name, param in model.named_parameters():
            if param.requires_grad:
                assert name in self.shadow
                self.original[name] = param.data.clone()
                param.data = self.shadow[name]

    def resume(self, model):
        for name, param in model.named_parameters():
            if param.requires_grad:
                assert name in self.shadow
                param.data = self.original[name]

# class self_conv(MessagePassing):
#     ## selfmade convolution layer

#     def __init__(
#         self, in_channels, out_channels, 
#         num_edge_types=1,
#         aggr='add', normalize=False, bias=True, device=constants.get_default_device(),
#         **kwargs,
#     ):
#         self.in_channels = in_channels
#         self.out_channels = out_channels
#         self.normalize = normalize

#         self.num_edge_types = num_edge_types

#         super().__init__(aggr, **kwargs)
        
#         self.lin_l = {i: Linear(in_channels, out_channels, bias=bias, device=device) for i in range(0, num_edge_types + 1)}

#         self.reset_parameters()

#     def reset_parameters(self):
#         self.aggr_module.reset_parameters()
#         for lin in self.lin_l.values():
#             lin.reset_parameters()

#     def forward(self, x, edge_index, edge_attr, size=None):
#         out = None
#         for i in range(1, self.num_edge_types + 1):
#             out_i = self.propagate(edge_index[:, edge_attr[:, i - 1]], x=x, size=size)
#             out_i = self.lin_l[i](out_i)

#             if out is None:
#                 out = out_i
#             else:
#                 out += out_i        

#         out = out + self.lin_l[0](x)
        
#         if self.normalize:
#             out = F.normalize(out, p=2., dim=-1)

#         return out

#     def message(self, x_j):
#         return x_j

#     # def message_and_aggregate(self, adj_t, x):
#     #     adj_t = adj_t.set_value(None, layout=None)
#     #     return torch_sparse.matmul(adj_t, x, reduce=self.aggr)

#     def __repr__(self):
#         return (f'{self.__class__.__name__}({self.in_channels}, '
#                 f'{self.out_channels}, aggr={self.aggr})')

class self_conv_no_edge_type(MessagePassing):
    ## selfmade convolution layer

    def __init__(
        self, in_channels, out_channels, 
        aggr='add', normalize=False, bias=False, device=constants.get_default_device(),
        **kwargs,
    ):
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.normalize = normalize

        # self.num_edge_types = 1

        super().__init__(aggr, **kwargs)
        
        # self.lin_l = {i: Linear(in_channels, out_channels, bias=bias, device=device) for i in range(0, num_edge_types + 1)}
        self.lin_l = Linear(in_channels, out_channels, bias=bias, device=device)
        self.lin_l_self = Linear(in_channels, out_channels, bias=bias, device=device)
        self.bias = Parameter(torch.empty(out_channels))

        self.reset_parameters()

    def reset_parameters(self):
        self.aggr_module.reset_parameters()
        # for lin in self.lin_l.values():
        #     lin.reset_parameters()
        self.lin_l.reset_parameters()
        self.lin_l_self.reset_parameters()
        self.bias.data.zero_()

    def forward(self, x, edge_index, size=None):
        out = self.propagate(edge_index, x=x, size=size)
        out = self.lin_l(out)       

        out = out + self.lin_l_self(x)

        out = out + self.bias
        
        if self.normalize:
            out = F.normalize(out, p=2., dim=-1)

        return out

    def message(self, x_j):
        return x_j

    # def message_and_aggregate(self, adj_t, x):
    #     adj_t = adj_t.set_value(None, layout=None)
    #     return torch_sparse.matmul(adj_t, x, reduce=self.aggr)

    def __repr__(self):
        return (f'{self.__class__.__name__}({self.in_channels}, '
                f'{self.out_channels}, aggr={self.aggr})')

# class GNN_int(torch.nn.Module):
#     def __init__(self, in_channels, hidden_channels, num_layers, p_channels, ANN_channels, num_edge_types, extra_x_length=None, dropout=None, device=constants.get_default_device()):
#         super().__init__()

#         # to do, maybe reduce the size of fv in hidden layer
#         # GNN layer
#         self.convs = torch.nn.ModuleList()
#         for i in range(num_layers):
#             if i == 0:
#                 self.convs.append(self_conv(in_channels, hidden_channels, num_edge_types, device=device))
#             else:
#                 self.convs.append(self_conv(hidden_channels, hidden_channels, num_edge_types, device=device))

#         # last layer of conv
#         self.convs_to_p = torch.nn.Linear(hidden_channels, p_channels)

#         # normal ANN
#         self.Linears = torch.nn.ModuleList()
#         if extra_x_length is None:
#             self.Linears.append(torch.nn.Linear(p_channels, ANN_channels[0]))
#         else:
#             self.Linears.append(torch.nn.Linear(p_channels + extra_x_length, ANN_channels[0]))
#         self.extra_x_length = extra_x_length

#         for i in range(len(ANN_channels) - 1):
#             self.Linears.append(torch.nn.Linear(ANN_channels[i], ANN_channels[i + 1]))

#         self.dropout = dropout

#         self.activation1 = torch.nn.ReLU()
#         self.activation2 = torch.nn.LeakyReLU(0.1)

#     def forward(self, data):
#         x, edge_index, edge_attr, batch = data.x, data.edge_index, data.edge_attr, data.batch

#         # x = self.embedding_layer(x)

#         for i in range(len(self.convs)):
#             x = self.convs[i](x, edge_index, edge_attr)
#             x = self.activation2(x)
#             if self.dropout is not None:
#                 x = F.dropout(x, p=self.dropout, training=self.training)

#         x = self.convs_to_p(x)
#         x = torch_geometric.nn.global_add_pool(x, batch, size=max(batch.tolist()) + 1)

#         if self.extra_x_length is not None:
#             extra_x = data.extra_x
#             x = torch.cat((extra_x, x), 1)

#         for i in range(len(self.Linears) - 1):
#             x = self.Linears[i](x)
#             x = self.activation1(x)

#         x = self.Linears[-1](x)

#         return x

# class GNN_int_no_edge_type(torch.nn.Module):
#     def __init__(self, in_channels, hidden_channels, num_layers, p_channels, ANN_channels, extra_x_length=None, dropout=None, device=constants.get_default_device()):
#         super().__init__()

#         # to do, maybe reduce the size of fv in hidden layer
#         # GNN layer
#         self.convs = torch.nn.ModuleList()
#         for i in range(num_layers):
#             if i == 0:
#                 self.convs.append(self_conv_no_edge_type(in_channels, hidden_channels, device=device))
#             else:
#                 self.convs.append(self_conv_no_edge_type(hidden_channels, hidden_channels, device=device))

#         # last layer of conv
#         self.convs_to_p = torch.nn.Linear(hidden_channels, p_channels)

#         # normal ANN
#         self.Linears = torch.nn.ModuleList()
#         if extra_x_length is None:
#             self.Linears.append(torch.nn.Linear(p_channels, ANN_channels[0]))
#         else:
#             self.Linears.append(torch.nn.Linear(p_channels + extra_x_length, ANN_channels[0]))
#         self.extra_x_length = extra_x_length

#         for i in range(len(ANN_channels) - 1):
#             self.Linears.append(torch.nn.Linear(ANN_channels[i], ANN_channels[i + 1]))

#         self.dropout = dropout

#         self.activation1 = torch.nn.ReLU()
#         self.activation2 = torch.nn.LeakyReLU(0.1)

#     def forward(self, data):
#         x, edge_index, batch = data.x, data.edge_index, data.batch

#         # x = self.embedding_layer(x)

#         for i in range(len(self.convs)):
#             x = self.convs[i](x, edge_index)
#             x = self.activation2(x)
#             if self.dropout is not None:
#                 x = F.dropout(x, p=self.dropout, training=self.training)

#         x = self.convs_to_p(x)
#         x = torch_geometric.nn.global_add_pool(x, batch, size=max(batch.tolist()) + 1)

#         if self.extra_x_length is not None:
#             extra_x = data.extra_x
#             x = torch.cat((extra_x, x), 1)

#         for i in range(len(self.Linears) - 1):
#             x = self.Linears[i](x)
#             x = self.activation1(x)

#         x = self.Linears[-1](x)

#         return x

# class GNN_FT(torch.nn.Module):
#     def __init__(self, in_channels, hidden_channels, num_layers, p_channels, ANN_channels, num_edge_types, dropout=None, device=constants.get_default_device()):
#         super().__init__()

#         self.embedding_layer = torch.nn.Linear(in_channels, hidden_channels)

#         # to do, maybe reduce the size of fv in hidden layer
#         # GNN layer
#         self.convs = torch.nn.ModuleList()
#         for i in range(num_layers):
#             # if i == 0:
#             #     self.convs.append(self_conv(in_channels, hidden_channels, num_edge_types, normalize=True, device=device))
#             # else:
#             self.convs.append(self_conv(hidden_channels, hidden_channels, num_edge_types, normalize=True, device=device))

#         # last layer of conv
#         self.convs_to_p = torch.nn.Linear(hidden_channels, p_channels)

#         # normal ANN
#         self.Linears = torch.nn.ModuleList()
#         self.Linears.append(torch.nn.Linear(p_channels, ANN_channels[0]))
        
#         for i in range(len(ANN_channels) - 1):
#             self.Linears.append(torch.nn.Linear(ANN_channels[i], ANN_channels[i + 1]))

#         self.dropout = dropout

#         self.activation1 = torch.nn.ReLU()
#         self.activation2 = torch.nn.LeakyReLU(0.1)

#         self.softmax = Softmax(dim=1)

#     def forward(self, data):
#         x, edge_index, edge_attr, batch = data.x, data.edge_index, data.edge_attr, data.batch

#         x = self.embedding_layer(x)

#         for i in range(len(self.convs)):
#             x = self.convs[i](x, edge_index, edge_attr)
#             x = self.activation2(x)
#             if self.dropout is not None:
#                 x = F.dropout(x, p=self.dropout, training=self.training)

#         x = self.convs_to_p(x)
#         x = self.activation2(x)
#         x = torch_geometric.nn.global_mean_pool(x, batch, size=max(batch.tolist()) + 1)

#         for i in range(len(self.Linears) - 1):
#             x = self.Linears[i](x)
#             x = self.activation1(x)

#         x = self.Linears[-1](x)

#         x = self.softmax(x)

#         return x

class GNN_FT_no_edge_type(torch.nn.Module):
    def __init__(self, in_channels, hidden_channels, num_layers, p_channels, ANN_channels, dropout=None, device=constants.get_default_device()):
        super().__init__()

        self.embedding_layer = torch.nn.Linear(in_channels, hidden_channels)

        # to do, maybe reduce the size of fv in hidden layer
        # GNN layer
        self.convs = torch.nn.ModuleList()
        for i in range(num_layers):
            # if i == 0:
            #     self.convs.append(self_conv_no_edge_type(in_channels, hidden_channels, normalize=True, device=device))
            # else:
            self.convs.append(self_conv_no_edge_type(hidden_channels, hidden_channels, normalize=True, device=device))

        # last layer of conv
        self.convs_to_p = torch.nn.Linear(hidden_channels, p_channels)

        # normal ANN
        self.Linears = torch.nn.ModuleList()
        self.Linears.append(torch.nn.Linear(p_channels, ANN_channels[0]))
        
        for i in range(len(ANN_channels) - 1):
            self.Linears.append(torch.nn.Linear(ANN_channels[i], ANN_channels[i + 1]))

        self.dropout = dropout

        self.activation1 = torch.nn.ReLU()
        self.activation2 = torch.nn.LeakyReLU(0.1)

        self.softmax = Softmax(dim=1)

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch

        x = self.embedding_layer(x)

        for i in range(len(self.convs)):
            x = self.convs[i](x, edge_index)
            x = self.activation2(x)
            if self.dropout is not None:
                x = F.dropout(x, p=self.dropout, training=self.training)

        x = self.convs_to_p(x)
        x = self.activation2(x)
        x = torch_geometric.nn.global_mean_pool(x, batch, size=max(batch.tolist()) + 1)

        for i in range(len(self.Linears) - 1):
            x = self.Linears[i](x)
            x = self.activation1(x)

        x = self.Linears[-1](x)

        x = self.softmax(x)

        return x

# class GNN_2L_int(torch.nn.Module):
#     def __init__(self, 
#         in_channels_atom, out_channels_FT, hidden_channels, num_layers, p_channels, ANN_channels,  
#         num_edge_types, extra_x_length=None, dropout=None, device=constants.get_default_device()
#     ):
#         super().__init__()

#         # GNN layer
#         self.convs = torch.nn.ModuleList()
#         for i in range(num_layers):
#             if i == 0:
#                 self.convs.append(self_conv(in_channels_atom + out_channels_FT, hidden_channels, num_edge_types, device=device))
#             else:
#                 self.convs.append(self_conv(hidden_channels, hidden_channels, num_edge_types, device=device))

#         # last layer of conv
#         self.convs_to_p = torch.nn.Linear(hidden_channels, p_channels)

#         # normal ANN
#         self.Linears = torch.nn.ModuleList()
#         if extra_x_length is None:
#             self.Linears.append(torch.nn.Linear(p_channels, ANN_channels[0]))
#         else:
#             self.Linears.append(torch.nn.Linear(p_channels + extra_x_length, ANN_channels[0]))
#         self.extra_x_length = extra_x_length

#         for i in range(len(ANN_channels) - 1):
#             self.Linears.append(torch.nn.Linear(ANN_channels[i], ANN_channels[i + 1]))

#         self.dropout = dropout

#         self.activation1 = torch.nn.ReLU()
#         self.activation2 = torch.nn.LeakyReLU(0.1)

#     def forward(self, data, FT_output):
#         x, x_FT_index, edge_index, edge_attr, batch = \
#             data.x, data.pos, data.edge_index, data.edge_attr, data.batch

#         x = torch.cat((x, FT_output[x_FT_index]), 1)

#         for i in range(len(self.convs)):
#             x = self.convs[i](x, edge_index, edge_attr)
#             x = self.activation2(x)
#             if self.dropout is not None:
#                 x = F.dropout(x, p=self.dropout, training=self.training)

#         x = self.convs_to_p(x)
#         x = self.activation2(x)
#         x = torch_geometric.nn.global_add_pool(x, batch, size=max(batch.tolist()) + 1)

#         if self.extra_x_length is not None:
#             extra_x = data.extra_x
#             x = torch.cat((extra_x, x), 1)

#         for i in range(len(self.Linears) - 1):
#             x = self.Linears[i](x)
#             x = self.activation1(x)

#         x = self.Linears[-1](x)

#         return x

class GNN_2L_int_no_edge_type(torch.nn.Module):
    def __init__(self, 
        in_channels_atom, out_channels_FT, hidden_channels, num_layers, p_channels, ANN_channels,  
        extra_x_length=None, dropout=None, clsf=False, device=constants.get_default_device()
    ):
        super().__init__()

        # GNN layer
        self.convs = torch.nn.ModuleList()
        for i in range(num_layers):
            if i == 0:
                self.convs.append(self_conv_no_edge_type(in_channels_atom + out_channels_FT, hidden_channels, device=device))
            else:
                self.convs.append(self_conv_no_edge_type(hidden_channels, hidden_channels, device=device))

        # last layer of conv
        self.convs_to_p = torch.nn.Linear(hidden_channels, p_channels, bias=False)

        # normal ANN
        self.Linears = torch.nn.ModuleList()
        if extra_x_length is None:
            self.Linears.append(torch.nn.Linear(p_channels, ANN_channels[0]))
        else:
            self.Linears.append(torch.nn.Linear(p_channels + extra_x_length, ANN_channels[0]))
        self.extra_x_length = extra_x_length

        for i in range(len(ANN_channels) - 1):
            self.Linears.append(torch.nn.Linear(ANN_channels[i], ANN_channels[i + 1]))

        self.dropout = dropout

        self.activation1 = torch.nn.ReLU()
        self.activation2 = torch.nn.LeakyReLU(0.1)

    def forward(self, data, FT_output):
        x, x_FT_index, edge_index, batch = \
            data.x, data.pos, data.edge_index, data.batch

        x = torch.cat((x, FT_output[x_FT_index, :]), 1)

        for i in range(len(self.convs)):
            x = self.convs[i](x, edge_index)
            x = self.activation2(x)
            if self.dropout is not None:
                x = F.dropout(x, p=self.dropout, training=self.training)

        x = self.convs_to_p(x)
        x = torch_geometric.nn.global_add_pool(x, batch, size=max(batch.tolist()) + 1)
        x = self.activation2(x)

        if self.extra_x_length is not None:
            extra_x = data.extra_x
            x = torch.cat((extra_x, x), 1)

        for i in range(len(self.Linears) - 1):
            x = self.Linears[i](x)
            x = self.activation1(x)

        x = self.Linears[-1](x)

        return x

# class GNN_2L(torch.nn.Module):
#     def __init__(self, 
#         in_channels_FT, hidden_channels_FT, num_layers_FT, p_channels_FT, ANN_channels_FT,
#         in_channels_IN, hidden_channels_IN, num_layers_IN, p_channels_IN, ANN_channels_IN,  
#         num_edge_types, extra_x_length=None, dropout=None, device=constants.get_default_device()
#     ):
#         super().__init__()

#         self.GNN_2L_int = GNN_2L_int(
#             in_channels_IN, ANN_channels_FT[-1], hidden_channels_IN, num_layers_IN, p_channels_IN, ANN_channels_IN,  
#             num_edge_types, extra_x_length=extra_x_length, dropout=dropout, device=device
#         )

#         self.GNN_FT = GNN_FT(
#             in_channels_FT, hidden_channels_FT, num_layers_FT, p_channels_FT, ANN_channels_FT, num_edge_types, 
#             dropout=dropout, device=device
#         )

#     def forward(self, data, FT_data):
#         FT_output = self.GNN_FT.forward(FT_data)
#         x = self.GNN_2L_int(data, FT_output)

#         return x

class GNN_2L_no_edge_type(torch.nn.Module):
    def __init__(self, 
        in_channels_FT, hidden_channels_FT, num_layers_FT, p_channels_FT, ANN_channels_FT,
        in_channels_IN, hidden_channels_IN, num_layers_IN, p_channels_IN, ANN_channels_IN,  
        extra_x_length=None, dropout=None, device=constants.get_default_device()
    ):
        super().__init__()

        self.GNN_2L_int = GNN_2L_int_no_edge_type(
            in_channels_IN, ANN_channels_FT[-1], hidden_channels_IN, num_layers_IN, p_channels_IN, ANN_channels_IN,  
            extra_x_length=extra_x_length, dropout=dropout, device=device
        )

        self.GNN_FT = GNN_FT_no_edge_type(
            in_channels_FT, hidden_channels_FT, num_layers_FT, p_channels_FT, ANN_channels_FT,
            dropout=dropout, device=device
        )

    def forward(self, data, FT_data):
        FT_output = self.GNN_FT.forward(FT_data)
        x = self.GNN_2L_int(data, FT_output)

        return x

class GNN_2L_no_edge_type_GraphSAGE(torch.nn.Module):
    def __init__(self, 
        in_channels_FT, hidden_channels_FT, num_layers_FT, p_channels_FT, ANN_channels_FT,
        in_channels_IN, hidden_channels_IN, num_layers_IN, p_channels_IN, ANN_channels_IN,  
        extra_x_length=None, dropout=None, device=constants.get_default_device()
    ):
        super().__init__()

        self.GNN_2L_int = GNN_2L_int_no_edge_type(
            in_channels_IN, ANN_channels_FT[-1], hidden_channels_IN, num_layers_IN, p_channels_IN, ANN_channels_IN,  
            extra_x_length=extra_x_length, dropout=dropout, device=device
        )

        # self.GNN_FT = GNN_FT_no_edge_type(
        #     in_channels_FT, hidden_channels_FT, num_layers_FT, p_channels_FT, ANN_channels_FT,
        #     dropout=dropout, device=device
        # )

        # self.GNN_FT = GAT(in_channels=in_channels_FT, hidden_channels=hidden_channels_FT, 
        #     num_layers=num_layers_FT, out_channels=ANN_channels_FT[-1], dropout=dropout)

        self.GNN_FT = GraphSAGE(in_channels=in_channels_FT, hidden_channels=hidden_channels_FT, 
            num_layers=num_layers_FT, out_channels=ANN_channels_FT[-1], dropout=dropout)


    def forward(self, data, FT_data):
        FT_output = self.GNN_FT.forward(x=FT_data.x, edge_index=FT_data.edge_index)
        x = self.GNN_2L_int(data, FT_output)

        return x



# class GraphSAGE(torch.nn.Module):
#     def __init__(self, in_channels, hidden_channels, num_layers, out_channels, dropout=None):
#         super().__init__()

#         self.embedding_layer = torch.nn.Linear(in_channels, hidden_channels)

#         self.convs = torch.nn.ModuleList()
#         for i in range(num_layers):
#             self.convs.append(SAGEConv(hidden_channels, hidden_channels))

#         self.linear1 = torch.nn.Linear(hidden_channels, hidden_channels)
#         self.linear2 = torch.nn.Linear(hidden_channels, out_channels)

#         self.dropout = dropout

#         self.activation = torch.nn.LeakyReLU(0.1)

#     def forward(self, data):
#         x, edge_index, batch = data.x, data.edge_index, data.batch

#         x = self.embedding_layer(x)
#         for i in range(len(self.convs) - 1):
#             x = self.convs[i](x, edge_index)
#             x = self.activation(x)
#             if self.dropout is not None:
#                 x = F.dropout(x, p=self.dropout, training=self.training)

#         x = self.convs[-1](x, edge_index)

#         x = global_add_pool(x, batch)

#         x = self.linear1(x)
#         x = self.activation(x)
#         x = self.linear2(x)

#         return x

# class GraphATConv(torch.nn.Module):
#     def __init__(self, in_channels, hidden_channels, num_layers, out_channels, dropout=None):
#         super().__init__()

#         self.embedding_layer = torch.nn.Linear(in_channels, hidden_channels)

#         self.convs = torch.nn.ModuleList()
#         for i in range(num_layers):
#             self.convs.append(GATConv(hidden_channels, hidden_channels))

#         self.linear1 = torch.nn.Linear(hidden_channels, hidden_channels)
#         self.linear2 = torch.nn.Linear(hidden_channels, out_channels)

#         self.dropout = dropout

#         self.activation = torch.nn.LeakyReLU(0.1)

#     def forward(self, data):
#         x, edge_index, edge_attr, batch = data.x, data.edge_index, data.edge_attr, data.batch

#         x = self.embedding_layer(x)
#         for i in range(len(self.convs) - 1):
#             x = self.convs[i](x, edge_index, edge_attr)
#             x = self.activation(x)
#             if self.dropout is not None:
#                 x = F.dropout(x, p=self.dropout, training=self.training)

#         x = self.convs[-1](x, edge_index, edge_attr)

#         x = global_add_pool(x, batch)

#         x = self.linear1(x)
#         x = self.activation(x)
#         x = self.linear2(x)

#         return x

