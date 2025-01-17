import pandas as pd
import numpy as np
import sys, warnings, argparse

warnings.simplefilter('ignore')

ZERO_TOL = 0.00000001

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv_ch1", help='input csv file for ch1') # _desc, not _desc_norm
    parser.add_argument("input_csv_ch2", help='input csv file for ch2') # _desc, not _desc_norm
    parser.add_argument("output_prefix", help='output prefix')
    parser.add_argument('-aT', '--addT', nargs='*', type=int)  ## use T as a feature, -aT (k1) (k2) ...: use T^k1, T^k2, ...
    parser.add_argument('-T', '--T', nargs='*', type=float)   ## specified T of selected compounds
                                                            ## -T: select all (with T information)
                                                            ## -T 0: select the compounds with the most common T
                                                            ## -T (t): select the compounds with the specified T
                                                            ## -T (t1) (t2): select the compounds between t1 and t2
    # parser.add_argument('-sp1', '--sp1', nargs='*', type=float)   ## special feature 1, need -T
    # parser.add_argument('-sp2', '--sp2', nargs='*', type=float)   ## special feature 2, need -T

    args = parser.parse_args()

    return(args) 

def most_common(lst):
    lst_tmp = set(lst)
    lst_tmp.discard(0.0)
    return max(lst_tmp, key=lst.count)

class compound_desc_combiner():
    def __init__(self, x1, x2):
        self.x1 = x1
        self.x2 = x2
        self.x = None
        # self.y = None
        self.T = None
        self.T_drop = None

    def init_x(self):
        CID_list = [CID[:-2] for CID in self.x1.index]
        self.T = [float(CID.split('_')[-1]) for CID in CID_list]
        self.x = pd.DataFrame(index=pd.Index(CID_list, name='CID'))

        for fv in list(self.x1.columns):
            fv_name = fv + "_ch1"
            self.x[fv_name] = self.x1[fv].values
        for fv in list(self.x2.columns):
            fv_name = fv + "_ch2"
            self.x[fv_name] = self.x2[fv].values

    def add_T(self, k):
        T_list_tmp = [pow(T, k) if T != 0 else 0 for T in self.T]
        self.x[f"T^{k}"] = T_list_tmp

    def add_sp1(self):
        T_list_tmp = [T if T != 0 else 1.0 for T in self.T]
        T_list_tmp = np.array(T_list_tmp)
        for fv in self.x.columns:
            fv_name = fv + "_T^-1"
            self.x[fv_name] = self.x[fv].values / T_list_tmp

    def add_sp2(self):
        T_list_tmp = [T if T != 0 else 1.0 for T in self.T]
        T_list_tmp_2 = [pow(T, 2) if T != 0 else 1.0 for T in self.T]
        T_list_tmp = np.array(T_list_tmp)
        T_list_tmp_2 = np.array(T_list_tmp_2)
        for fv in self.x.columns:
            fv_name = fv + "_T-1"
            self.x[fv_name] = self.x[fv].values / T_list_tmp
            fv_name = fv + "_T-2"
            self.x[fv_name] = self.x[fv].values / T_list_tmp_2

    def select_from_T(self, T_opt):
        #### here only make a mask list of which compound to be selected

        if T_opt is None:
            # select all, no matter T information exists or not
            self.T_drop = [0] * len(self.x.index)
        elif len(T_opt) == 0:
            # select all with T information
            self.T_drop = [0 if T != 0 else 1 for T in self.T]
        elif len(T_opt) == 1:
            if T_opt[0] == 0:
                # select the most common T
                T_common = most_common(self.T)
                self.T_drop = [0 if T == T_common else 1 for T in self.T]
                print(f"The most common T is: {T_common}")
            else:
                # given T
                self.T_drop = [0 if T == T_opt[0] else 1 for T in self.T]
        elif len(T_opt) == 2:
            # between T1 and T2
            self.T_drop = [0 if T >= T_opt[0] and T <= T_opt[1] else 1 for T in self.T]
        else:
            # no definition
            self.T_drop = [0] * len(self.x.index)

    def normalize(self):
        min_x = np.amin(self.x, axis=0)
        max_x = np.amax(self.x, axis=0)
        self.x_norm = (self.x - min_x) / (max_x - min_x)

    def remove_var0(self):
        x_tmp = self.x.values
        K = x_tmp.shape[1]
        zero_column = [self.x.columns[i] for i in range(K) if abs(np.var(x_tmp[:, i])) < ZERO_TOL and 'T^' not in self.x.columns[i]]
        self.x = self.x.drop(columns=zero_column)

    def output(self, output_prefix):
        #### do the selection of compounds here
        compound_to_drop = [ind for (ind, T_mask) in zip(self.x.index, self.T_drop) if T_mask == 1]
        self.x = self.x.drop(index=compound_to_drop)
        print(f"#selected compounds: {self.x.values.shape[0]}")

        # output_descname = output_prefix + "_desc.csv"
        # self.x.to_csv(output_descname, sep=',')
        # print(f"#feature: {self.x.values.shape[1]}")
        
        self.remove_var0()
        output_descname = output_prefix + "_desc.csv"   
        self.x.to_csv(output_descname, sep=',')   
        # print(f"#var>0 feature: {self.x.values.shape[1]}")  
        
        self.normalize()
        output_descname = output_prefix + "_desc_norm.csv"
        self.x_norm.to_csv(output_descname, sep=',')   

def main(args):

    ch1_descname = args.input_csv_ch1
    ch2_descname = args.input_csv_ch2

    csv_data_1 = pd.read_csv(ch1_descname, index_col=0)
    csv_data_2 = pd.read_csv(ch2_descname, index_col=0)

    combiner = compound_desc_combiner(csv_data_1, csv_data_2)
    combiner.init_x()

    if args.addT is not None:
        for k in args.addT:
            combiner.add_T(int(k))

    # if args.sp1 is not None:
    #     if args.T is None:
    #         print("-sp1 needs -T !!!")
    #         sys.exit()
    #     else:
    #         combiner.add_sp1()

    # if args.sp2 is not None:
    #     if args.T is None:
    #         print("-sp2 needs -T !!!")
    #         sys.exit()
    #     else:
    #         combiner.add_sp2()

    combiner.select_from_T(args.T)

    combiner.output(args.output_prefix)


if __name__ == "__main__":
    main(get_args())
