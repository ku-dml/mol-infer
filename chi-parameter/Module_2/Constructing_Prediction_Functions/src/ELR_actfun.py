import math
import numpy as np
import numbers

def act_fun_0(x):
    return x

def act_fun_1(x):

    if x <= 1/3:
        y = 2 * x
    else:
        y = 0.5 * (x - 1/3) + 2/3

    return y

def act_fun_2(x):

    if x <= 2/3:
        y = 1.25 * x
    else:
        y = 0.5 * (x - 2/3) + 5/6

    return y

def act_fun_3(x):

    if x <= 1/3:
        y = 0.5 * x
    else:
        y = 1.25 * (x - 1/3) + 1/6

    return y

def act_fun_4(x):

    if x <= 2/3:
        y = 0.5 * x
    else:
        y = 2 * (x - 2/3) + 1/3

    return y

def act_fun_x2(x):
    return x ** 2

def act_fun_x2_alt(x):
    return 1 - (1 - x) ** 2

def act_fun_x3(x):
    return x ** 3

def act_fun_x3_alt(x):
    return 1 - (1 - x) ** 3

def act_fun_sp1(x):
    if x == 0.5:
        y = x
    else:
        y = ((x - 0.5) * 2) ** 3 / 2 + 0.5
    
    return y

def act_fun_sp2(x):
    if x == 0.5:
        y = x
    elif x > 0.5:
        y = ((x - 0.5) * 2) ** (1./3.) / 2 + 0.5
    else:
        y = -(((0.5 - x) * 2) ** (1./3.)) / 2 + 0.5
    
    return y

def act_fun_half_main(x, t):
    if x <= t:
        y = 0.5 / t * x
    else:
        y = 0.5 / (1 - t) * (x - 1) + 1

    return y

def act_fun_gen_main(x, t1, t2):


    if t1 == 0:
        y = x
    elif t1 == 1:
        y = x
    elif x <= t1:
        y = t2 / t1 * x
    else:
        y = (1 - t2) / (1 - t1) * (x - 1) + 1

    # print(x, t1, t2, y)

    return y

def act_fun_half_p4(x):

    return act_fun_half_main(x, 0.4)

def act_fun_half_p6(x):

    return act_fun_half_main(x, 0.6)

def act_fun_half_p25(x):

    return act_fun_half_main(x, 0.25)

def act_fun_half_p75(x):

    return act_fun_half_main(x, 0.75)

def act_fun_half_q1_main(x, t):

    if x <= t:
        y = -0.5 / (t ** 2) * ((x - t) ** 2) + 0.5
    else:
        y = 0.5 / ((1 - t) ** 2) * ((x - t) ** 2) + 0.5

    return y

def act_fun_half_q1_p4(x):

    return act_fun_half_q1_main(x, 0.4)

def act_fun_half_q1_p6(x):

    return act_fun_half_q1_main(x, 0.6)

def act_fun_half_q1_p25(x):

    return act_fun_half_q1_main(x, 0.25)

def act_fun_half_q1_p75(x):

    return act_fun_half_q1_main(x, 0.75)

def act_fun_half_q2_main(x, t):

    if x <= t:
        y = 0.5 / (t ** 2) * (x ** 2)
    else:
        y = -0.5 / ((1 - t) ** 2) * ((1 - x) ** 2) + 1

    return y

def act_fun_half_q2_p4(x):

    return act_fun_half_q2_main(x, 0.4)

def act_fun_half_q2_p6(x):

    return act_fun_half_q2_main(x, 0.6)

def act_fun_half_q2_p25(x):

    return act_fun_half_q2_main(x, 0.25)

def act_fun_half_q2_p75(x):

    return act_fun_half_q2_main(x, 0.75)


def act_fun_e_main(x, t):

    return math.exp(x * t) / math.exp(t)

def act_fun_e_1(x):

    return act_fun_e_main(x, 1)

def act_fun_e_2(x):

    return act_fun_e_main(x, 2)

def act_fun_e_3(x):

    return act_fun_e_main(x, 3)

def act_fun_e_4(x):

    return act_fun_e_main(x, 4)

def act_fun_e_5(x):

    return act_fun_e_main(x, 5)

def act_fun_e_6(x):

    return act_fun_e_main(x, 6)

def act_fun_e_7(x):

    return act_fun_e_main(x, 7)

def act_fun_e_8(x):

    return act_fun_e_main(x, 8)

def act_fun_e_9(x):

    return act_fun_e_main(x, 9)

def act_fun_e_10(x):

    return act_fun_e_main(x, 10)

def act_fun_sin_x(x, t):

    return math.sin(2 * math.pi * t * x)

def act_fun_cos_x(x, t):

    return math.cos(2 * math.pi * t * x)

def act_fun_f_list(n):

    ans = [act_fun_0]

    for i in range(1, n + 1):
        ans.append(lambda x: act_fun_sin_x(x, i))

    for i in range(1, n + 1):
        ans.append(lambda x: act_fun_cos_x(x, i))

    return ans

def act_fun_thun_x(x, p, t):
    ans = 0
    if x <= t:
        ans = x * p
    elif x >= 1- t:
        ans = 1 - p * (1 - x)
    else:
        ans = ((1 - 2 * p * t) * x - (1 - p) * t) / (1 - 2 * t)

    return ans

def act_fun_q_main(_list, x):
    ans = []

    for _l in _list:
        if _l == 0:
            ans.append(act_fun_0)
        elif isinstance(_l, numbers.Number):
            q = np.percentile(x, _l)
            # if q == 0 or q == 1:
            #     ans.append(act_fun_0)
            # else:
            ans.append(lambda x, qq=q, _ll=_l: act_fun_gen_main(x, qq, _ll))
        else:
            ans.append(_l)

    return ans

# def act_fun_f_diff(n, x):

#     ans = [1]

#     for i in range(1, n + 1):
#         ans.append(lambda x: 2 * math.pi * i * act_fun_cos_x(x, i))

#     for i in range(1, n + 1):
#         ans.append(lambda x: -2 * math.pi * i * act_fun_sin_x(x, i))

#     return ans

def act_fun_all(x, list_act_fun, coef_act_fun):
    return sum([f(x) * e for f, e in zip(list_act_fun, coef_act_fun)])