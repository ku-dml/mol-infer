import torch, torchvision
import numpy as np

from torch_geometric.data import Data, InMemoryDataset

from src import constants

# .x, .edge_index, .edge_attr, .y

class ChemicalGraph_Datasets(InMemoryDataset):
    def __init__(self, root, transform=None):
        super().__init__(root, transform)

        self.transform = transform
        self.data_list = list()
        # self.data = None


    def __len__(self):
        return self.data_num

    def __getitem__(self, idx):
        return self.data_list[idx]

    @property
    def processed_file_names(self):
        return 'data.pt'

    def get(self, idx):
        return self.data_list[idx]

    def construct_from_ChemicalGraph_list(self, cg_list, clsf=False, FT_use=True, device=constants.get_default_device()):
        # u_idx_dict = dict()
        # u_idx = 0
        extra_x_used = False

        if cg_list.extra_fv_dict is not None:
            extra_x_used = True

        for i, cg in enumerate(cg_list):
            u_idx_dict = dict()
            u_idx = 0

            x = list()
            x_FT = list()
            x_FT_index = list()
            edge_index = [list(), list()]
            edge_attr = list()
            y = list()
            idx = list()

            for u in cg.interior_indices():
                index_u = u_idx
                u_idx_dict[u] = index_u
                u_idx += 1

            for u in cg.interior_indices():  
                idx_u = u_idx_dict[u]

                atom_u = cg[u]

                x.append(atom_u.get_feature())
                x_FT.append(atom_u.get_feature_FT(cg_list.TS_list))
                # x_FT_index.append([atom_u.fringe_tree_TS])
                x_FT_index.append(atom_u.fringe_tree_TS)

                ### ##
                # maybe no need to use this part
                # if FT_use:
                    # x[-1] = np.append(x[-1], x_FT[-1])
                ###

                for v in atom_u.adj:
                    if v not in cg.interior_indices():
                        continue

                    idx_v = u_idx_dict[v]
                    edge_index[0].append(idx_u)
                    edge_index[1].append(idx_v)
                    edge_attr_tmp = [0, 0, 0]
                    edge_attr_tmp[atom_u.beta[v] - 1] = 1
                    edge_attr.append(edge_attr_tmp)

            y.append(cg_list.get_value(cg.CID))
            
            x = torch.tensor(np.array(x), dtype=torch.float, device=device) 
            x_FT = torch.tensor(np.array(x_FT), dtype=torch.float, device=device)
            pos = torch.tensor(np.array(x_FT_index), dtype=torch.long, device=device)
            edge_index = torch.tensor(np.array(edge_index), dtype=torch.long, device=device)
            # edge_attr = torch.BoolTensor(np.array(edge_attr))
            edge_attr = torch.tensor(np.array(edge_attr), dtype=torch.bool, device=device)
            y = torch.tensor(np.array(y), dtype=torch.float, device=device) 

            if extra_x_used:
                extra_x = [cg_list.extra_fv_dict[cg.CID]]
                extra_x = torch.tensor(np.array(extra_x), dtype=torch.float, device=device)

            if extra_x_used:
                data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr, y=y, extra_x=extra_x, name=cg.CID, idx=i)
                data.x_FT = x_FT
                data.pos = pos
            else:
                data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr, y=y, name=cg.CID, idx=i)
                data.x_FT = x_FT
                data.pos = pos

            self.data_list.append(data)

        self.data_num = len(self.data_list)

        self.data, self.slices = self.collate(self.data_list)

class FringeTree_Datasets(InMemoryDataset):
    def __init__(self, root, transform=None):
        super().__init__(root, transform)

        self.transform = transform
        self.data_list = list()
        # self.data = None

    def __len__(self):
        return self.data_num

    def __getitem__(self, idx):
        return self.data_list[idx]

    @property
    def processed_file_names(self):
        return 'data.pt'

    def get(self, idx):
        return self.data_list[idx]

    def construct_from_FringeTree_list(self, ft_list, device=constants.get_default_device()):
        for i, rt in enumerate(ft_list):
            u_idx_dict = dict()
            u_idx = 0

            x = list()
            edge_index = [list(), list()]
            edge_attr = list()
            idx = list()

            for u in rt.vertex_list.keys():
                index_u = u_idx
                u_idx_dict[u] = index_u
                u_idx += 1

            for u in rt.vertex_list.keys():  
                idx_u = u_idx_dict[u]

                atom_u = rt[u]

                x.append(atom_u.get_feature_for_FT(rt.root))

                for v in atom_u.adj:
                    idx_v = u_idx_dict[v]
                    edge_index[0].append(idx_u)
                    edge_index[1].append(idx_v)
                    edge_attr_tmp = [0, 0, 0]
                    edge_attr_tmp[atom_u.beta[v] - 1] = 1
                    edge_attr.append(edge_attr_tmp)

            x = torch.tensor(np.array(x), dtype=torch.float, device=device) 
            edge_index = torch.tensor(np.array(edge_index), dtype=torch.long, device=device)
            # edge_attr = torch.BoolTensor(np.array(edge_attr))
            edge_attr = torch.tensor(np.array(edge_attr), dtype=torch.bool, device=device)
            data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr, name=f"FT_{i+1}", idx=i)
            self.data_list.append(data)

        self.data_num = len(self.data_list)

        self.data, self.slices = self.collate(self.data_list)






