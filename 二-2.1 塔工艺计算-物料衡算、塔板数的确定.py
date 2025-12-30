import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
from matplotlib.font_manager import FontProperties
import matplotlib
import math
import json
from reflux_ratio_calculator import EnhancedRefluxCalculator
times_new_roman = FontProperties(family='Times New Roman')
matplotlib.rcParams['font.family'] = 'SimSun'
matplotlib.rcParams['axes.unicode_minus'] = False

# ===== 设计参数 ===============================================
# 根据你的设计要求修改
Annual_processing_capacity = 60000  # 年处理量 吨/年，需要换算
work_time = 8000  # 工作时间 小时/年
t_f = 45  # 进料温度 ℃
wt_f = 0.25  # 进料组成_wt 苯的质量分数
t_ewflux = 60  # 回流温度 ℃ 如果回流为泡点回流，则输入字符串格式'泡点回流'
D_wt = 0.91  # 塔顶组成
F_wt = 0.025  # 塔釜组成≤
size = 6  # 图幅
R_multiplier = 1.5  # 回流比系数（可根据设计调整，通常为1.2~2倍）

# =====以下请勿更改==================================================


# ===== 苯-甲苯物性数据 =====
MA_benzene = 78.11  # 苯的摩尔质量 g/mol
MA_toluene = 92.14  # 甲苯的摩尔质量 g/mol
MA = MA_benzene  # 苯
MB = MA_toluene  # 甲苯
# 计算摩尔流量（年处理6万吨）
F_kg_per_h = Annual_processing_capacity * 1000 / work_time  # kg/h


# 质量分数转换为摩尔分数
def wt_to_mol(wt_benzene):
    x_benzene = (wt_benzene / MA_benzene) / (wt_benzene / MA_benzene + (1 - wt_benzene) / MA_toluene)
    return x_benzene


x_F = wt_to_mol(wt_f)  # 进料摩尔分数
x_D = wt_to_mol(D_wt)  # 塔顶摩尔分数
x_W = wt_to_mol(F_wt)  # 塔釜摩尔分数

print(f"进料摩尔分数 x_F = {x_F:.4f}")
print(f"塔顶摩尔分数 x_D = {x_D:.4f}")
print(f"塔釜摩尔分数 x_W = {x_W:.4f}")

# ===== 苯-甲苯相平衡数据 (常压760mmHg) =====
# 来源：化工原理教材，苯-甲苯汽液平衡数据
dic_txy = {0.0: [0.0, 110.56],
           0.01: [0.025, 109.91],
           0.03: [0.0711, 108.79],
           0.05: [0.112, 107.61],
           0.1: [0.208, 105.05],
           0.15: [0.294, 102.79],
           0.2: [0.372, 100.75],
           0.25: [0.442, 98.84],
           0.3: [0.507, 97.13],
           0.35: [0.566, 95.58],
           0.4: [0.619, 94.09],
           0.45: [0.667, 92.69],
           0.5: [0.713, 91.4],
           0.55: [0.755, 90.11],
           0.6: [0.791, 80.80],
           0.65: [0.825, 87.63],
           0.7: [0.857, 86.52],
           0.75: [0.885, 85.44],
           0.8: [0.912, 84.40],
           0.85: [0.936, 83.33],
           0.9: [0.959, 82.25],
           0.95: [0.98, 81.11],
           0.97: [0.988, 80.66],
           0.99: [0.9961, 80.21],
           1.0: [1.0, 80.01],
           }

ding = x_D  # 塔顶
fu = x_W  # 塔釜
liao = x_F  # 进料

ding2 = x_D  # 部分回流塔顶
fu2 = x_W  # 部分回流塔釜
tF = t_f  # 进料温度

# ===== 系统参数 =====


# ===== 苯-甲苯物性数据 =====
# 汽化潜热 (kJ/kg)
r_benzene = 393.9  # 苯在80.1℃时的汽化潜热
r_toluene = 363.0  # 甲苯在110.6℃时的汽化潜热

# 比热容数据 (kJ/(kg·℃))
Cp_benzene = {
    0: 1.570,
    20: 1.716,
    40: 1.767,
    60: 1.828,
    80: 1.881,
    100: 1.953,
}

Cp_toluene = {
    0: 1.630,
    20: 1.681,
    40: 1.757,
    60: 1.834,
    80: 1.902,
    100: 1.970,
    120: 2.073
}
# 黏度数据 (mPa*s)
m_benzene = {
    0: 0.742,
    20: 0.638,
    40: 0.485,
    60: 0.381,
    80: 0.308,
    100: 0.255,
    120: 0.215
}

m_toluene = {
    0: 0.758,
    20: 0.580,
    40: 0.459,
    60: 0.373,
    80: 0.311,
    100: 0.264,
    120: 0.228
}


# ===== 插值函数 =====
def chazhi(T, dic, i='n'):
    if i == 'n':
        if T in dic:
            return dic[T]
        sorted_keys = sorted(dic.keys())
        lower_key = None
        upper_key = None
        for key in sorted_keys:
            if key < T:
                lower_key = key
            else:
                upper_key = key
                break
        if lower_key is None or upper_key is None:
            return f"插值无法计算，T超出字典{dic}键的范围"
        lower_value = dic[lower_key]
        upper_value = dic[upper_key]
        interpolation = (T - lower_key) * (upper_value - lower_value) / (upper_key - lower_key) + lower_value
        return interpolation

    if T in dic:
        return dic[T][i]
    sorted_keys = sorted(dic.keys())
    lower_key = None
    upper_key = None
    for key in sorted_keys:
        if key < T:
            lower_key = key
        else:
            upper_key = key
            break
    if lower_key is None or upper_key is None:
        return f"插值无法计算，T超出字典{dic}键的范围"
    lower_value = dic[lower_key][i]
    upper_value = dic[upper_key][i]
    interpolation = (T - lower_key) * (upper_value - lower_value) / (upper_key - lower_key) + lower_value
    return interpolation


def chazhixian(x, y):
    spline = UnivariateSpline(x, y, s=0.01)
    xnew = np.linspace(x.min(), x.max(), 300)
    ynew = spline(xnew)
    return xnew, ynew


def find_interpolation_points(x_target, data_dict, value_index=1):
    """
    找到用于插值的两个相邻数据点
    Args:
        x_target: 目标x值
        data_dict: 数据字典 {x: [y, t]}
        value_index: 要返回的值索引 (0: y, 1: t)
    Returns:
        (x1, x2, val1, val2): 相邻两个点的x值和对应的值
    """
    # 获取所有x值并排序
    x_values = sorted(data_dict.keys())

    # 如果目标值在字典中，直接返回相等的点
    if x_target in data_dict:
        return x_target, x_target, data_dict[x_target][value_index], data_dict[x_target][value_index]

    # 查找相邻点
    for i in range(len(x_values) - 1):
        x1, x2 = x_values[i], x_values[i + 1]
        if x1 <= x_target <= x2:
            val1 = data_dict[x1][value_index]
            val2 = data_dict[x2][value_index]
            return x1, x2, val1, val2

    # 如果目标值超出范围，返回边界点
    if x_target < x_values[0]:
        return x_values[0], x_values[1], data_dict[x_values[0]][value_index], data_dict[x_values[1]][value_index]
    else:
        return x_values[-2], x_values[-1], data_dict[x_values[-2]][value_index], data_dict[x_values[-1]][value_index]
# ===== 数据处理 =====
ls_XA_xiangtu = list(dic_txy.keys())
ls_YA_xiangtu = [value[0] for value in dic_txy.values()]

ls_XA = [ding, fu, ding2, fu2, liao]
ls_tBP = []

# 进料平均分子量
M_F = x_F * MA_benzene + (1 - x_F) * MA_toluene
F_kmol_per_h = F_kg_per_h / M_F

print(f"\n物料流量计算:")
print(f"年处理量: 60000 t/a")
print(f"工作时间: 8000 h/a")
print(f"进料质量流量: {F_kg_per_h:.2f} kg/h")
print(f"进料平均分子量: {M_F:.2f} g/mol")
print(f"进料摩尔流量: {F_kmol_per_h:.2f} kmol/h")

# 计算进料泡点温度
tBP_F = chazhi(liao, dic_txy, 1)
print(f"\n进料组成 {liao:.4f} 摩尔分数下，泡点温度 tBP = {tBP_F:.1f}℃")
if t_ewflux == '泡点回流':
    t_ewflux=tBP_F
# 计算进料热状况参数 q
if tF < tBP_F:  # 冷液进料
    print(f"进料温度 {tF}℃ < 泡点温度 {tBP_F:.1f}℃，为冷液进料")
    # 计算平均温度下的比热容
    t_avg = (tF + tBP_F) / 2
    # 苯的比热容 (kJ/(kg·℃)) 转换为 (kJ/(kmol·℃))
    Cp_benzene_mol = chazhi(t_avg, Cp_benzene, 'n') * MA_benzene if t_avg in Cp_benzene else 1.8 * MA_benzene
    # 甲苯的比热容
    Cp_toluene_mol = chazhi(t_avg, Cp_toluene, 'n') * MA_toluene if t_avg in Cp_toluene else 1.8 * MA_toluene
    # 混合物平均摩尔热容
    Cpm = Cp_benzene_mol * liao + Cp_toluene_mol * (1 - liao)
    # 混合物汽化潜热 (kJ/kmol)
    rm = r_benzene * MA_benzene * liao + r_toluene * MA_toluene * (1 - liao)
    # 进料热状况参数
    q = 1 + (Cpm * (tBP_F - tF)) / rm
    print(f"平均温度: {t_avg:.1f}℃")
    print(f"混合物平均摩尔热容 Cpm = {Cpm:.2f} kJ/(kmol·℃)")
    print(f"混合物汽化潜热 rm = {rm:.0f} kJ/kmol")
    print(f"进料热状况参数 q = {q:.3f}")
else:
    q = 1.0  # 饱和液体进料
    print(f"进料温度 ≥ 泡点温度，q = {q}（泡点进料）")

from txy_calculator import get_property_temperature

# 假设您已经定义了txy_dict
# 假设您已经计算出xD和xW

# 获取用于物性查询的温度
t_avg_for_properties = get_property_temperature(dic_txy, x_D, x_W)

t_avg2 = t_avg_for_properties
m_A = chazhi(t_avg2, m_benzene)
m_B = chazhi(t_avg2, m_toluene)
m_m = x_F * m_A + (1 - x_F) * m_B
ET = 0.17 - 0.616 * np.log10(m_m)
print(f'平均温度下粘度为{m_m:.3f}mPa·s')
print(f'塔板效率为{ET * 100:.2f}%')

# ======最小回流比计算=========


print("\n" + "=" * 60)
print("最小回流比计算")
print("=" * 60)
print("\n" + "=" * 70)
print("考虑过冷回流的回流比计算")
print("=" * 70)

# 准备物性数据
Cp_data = {
    'benzene': Cp_benzene,
    'toluene': Cp_toluene
}

r_data = {
    'benzene': r_benzene,  # 已在前面定义
    'toluene': r_toluene  # 已在前面定义
}

M_components = {
    'benzene': MA_benzene,
    'toluene': MA_toluene
}

# 创建增强版计算器
enhanced_calc = EnhancedRefluxCalculator(
    dic_txy,
    Cp_data=Cp_data,
    r_data=r_data,
    M_components=M_components
)

# 计算q值（使用更精确的方法）
# 这里可以使用你之前计算的q值，或者重新计算

q_calculated = q
t_bubble = tBP_F
# 计算最小回流比
R_min, intersection = enhanced_calc.calculate_R_min(x_D, x_F, q_calculated)
# 在计算最小回流比后，添加以下代码：
# 对交点坐标保留三位小数
x_eq_rounded = round(intersection[0], 3)
y_eq_rounded = round(intersection[1], 3)
x_D_rounded=round(x_D,4)
# 使用保留后的交点重新计算最小回流比
R_min_rounded = (x_D_rounded - y_eq_rounded) / (y_eq_rounded - x_eq_rounded)

print(f"\n使用保留三位小数的交点计算最小回流比:")
#  print(f"原始交点: x={intersection[0]:.6f}, y={intersection[1]:.6f}")
print(f"保留三位小数: x={x_eq_rounded:.3f}, y={y_eq_rounded:.3f}")
#  print(f"原始 R_min = {R_min:.4f}")
print(f"保留后 R_min = {R_min_rounded:.4f}")

# 使用保留后的R_min进行后续计算
R_min = R_min_rounded  # 如果希望使用保留后的值
print(f"最小回流比 R_min = {R_min:.4f}")

# 计算塔顶温度
t_top = float(enhanced_calc.t_from_x(x_D))
print(f"塔顶温度: {t_top:.1f} °C")

# 计算考虑过冷回流的实际回流比

R_actual, R_internal, correction = enhanced_calc.calculate_actual_reflux_ratio(
    R_min, R_multiplier, t_ewflux, t_top, x_D
)

print(f"\n回流比计算结果:")
print(f"理论外回流比 (R_min×{R_multiplier}): {R_actual:.4f}")
print(f"实际内回流比 (考虑过冷): {R_internal:.4f}")
print(f"过冷回流修正系数: {correction:.4f}")

# ===== 修改后续计算，使用内回流比 =====
# 注意：在阶梯图解法中，应该使用内回流比
R = R_internal  # 使用内回流比进行计算

print('\n' + '=' * 50)
print(f'部分回流设计分析 (内回流比 R\' = {R:.3f})')
print('=' * 50)


# 后续的精馏段操作线方程应该使用内回流比
# 原来的精馏段操作线方程：
# y = (R/(R+1)) * x + x_D/(R+1)
# 现在R应该是R_internal

# 修改精馏段操作线函数
def jingliu_line(x_d, R_value):
    x = np.linspace(0, x_d, 1000)
    y = (R_value / (R_value + 1)) * x + x_d / (R_value + 1)
    plt.plot(x, y, color='brown', label=f'精馏段操作线 (R\'={R_value:.3f})')
    print(f'精馏段操作线: y = {R_value / (R_value + 1):.3f}x + {x_d / (R_value + 1):.3f}')
    return x, y


def find_closest_index(a, b, target):
    idx = min(range(len(a)), key=lambda i: abs(a[i] - target))
    return b[idx]


#  ===== 部分回流（实际设计）=====
print('\n' + '=' * 50)
print('部分回流设计分析:')
print('=' * 50)

plt.figure(figsize=(size, size))
x, y = chazhixian(np.array(ls_XA_xiangtu), ls_YA_xiangtu)
plt.plot(x, y)
plt.plot([0, 1], [0, 1])
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.scatter(ls_XA_xiangtu,ls_YA_xiangtu)
# 画设计点
plt.plot([ls_XA[2], ls_XA[2]], [0, ls_XA[2]], linestyle='--', color='red')  # 塔顶
plt.plot([ls_XA[3], ls_XA[3]], [0, find_closest_index(x, y, ls_XA[3])], linestyle='--', color='red',
         )  # 塔釜


# q线方程
def q_line(x_f):
    x = np.linspace(x_f - 0.05, x_f + 0.05 + 0.3, 1000)
    if q != 1:
        y = (q / (q - 1)) * x - x_f / (q - 1)
    else:
        y = x * 0 + x_f  # q=1时为垂直线
    plt.plot(x, y, color='purple', label='q线')
    print(f'q线方程: y = {(q / (q - 1)):.3f}x + {(-x_f / (q - 1)):.3f}') if q != 1 else print('q线: x = 常数线')
    return x, y


x_q, y_q = q_line(ls_XA[4])


# 精馏段操作线方程
def jingliu_line(x_d, R_value):
    x = np.linspace(0, x_d, 1000)
    y = (R_value / (R_value + 1)) * x + x_d / (R_value + 1)
    plt.plot(x, y, color='brown', label=f'精馏段操作线 (R={R_value})')
    print(f'精馏段操作线: y = {R_value / (R_value + 1):.3f}x + {x_d / (R_value + 1):.3f}')
    return x, y


x_caozuo, y_caozuo = jingliu_line(ls_XA[2], R)


# 求q线与精馏段操作线交点
def solve_linear_equations(a1, b1, a2, b2):
    if a1 == a2:
        return None
    x = (b2 - b1) / (a1 - a2)
    y = a1 * x + b1
    return x, y


if q != 1:
    x_intersection, y_intersection = solve_linear_equations(
        (q / (q - 1)), -ls_XA[4] / (q - 1),
        (R / (R + 1)), ls_XA[2] / (R + 1)
    )
    print(f"q线与精馏段操作线交点: ({x_intersection:.4f}, {y_intersection:.4f})")
else:
    x_intersection, y_intersection = ls_XA[4], ls_XA[2] / (R + 1) + (R / (R + 1)) * ls_XA[4]

# 提馏段操作线
plt.plot([ls_XA[3], x_intersection], [ls_XA[3], y_intersection],
         color='blue', linewidth=1.5, label='提馏段操作线')


# 阶梯作图法求理论板数

def bufenhuiliu_caozuo_zhexian(x_break):
    x_zheqi = x_break
    # 在精馏段操作线上找y
    idx = np.argmin(np.abs(x_caozuo - x_zheqi))
    y_zheqi = y_caozuo[idx]
    # 在平衡线上找对应的x
    idx2 = np.argmin(np.abs(y - y_zheqi))
    x_zhezhong = x[idx2]
    # 再回到操作线上找y
    idx3 = np.argmin(np.abs(x_caozuo - x_zhezhong))
    y_zhezhong = y_caozuo[idx3]
    x_break_break = x_zhezhong
    plt.plot([x_zhezhong, x_zheqi], [y_zheqi, y_zheqi], color='green')
    plt.plot([x_zhezhong, x_zhezhong], [y_zhezhong, y_zheqi], color='green')
    return x_break_break, y_zhezhong


def bufenhuiliu_tiliu_zhexian(x_break):
    x_zheqi = x_break
    # 在提馏段操作线上找y（线性插值）
    if x_zheqi <= x_intersection and x_zheqi >= ls_XA[3]:
        y_zheqi = ls_XA[3] + (y_intersection - ls_XA[3]) * (x_zheqi - ls_XA[3]) / (x_intersection - ls_XA[3])
    else:
        y_zheqi = ls_XA[3]
    # 在平衡线上找对应的x
    idx = np.argmin(np.abs(y - y_zheqi))
    x_zhezhong = x[idx]
    # 再回到提馏段操作线上找y
    if x_zhezhong <= x_intersection and x_zhezhong >= ls_XA[3]:
        y_zhezhong = ls_XA[3] + (y_intersection - ls_XA[3]) * (x_zhezhong - ls_XA[3]) / (x_intersection - ls_XA[3])
    else:
        y_zhezhong = ls_XA[3]
    x_break_break = x_zhezhong
    plt.plot([x_zhezhong, x_zheqi], [y_zheqi, y_zheqi], color='green')
    plt.plot([x_zhezhong, x_zhezhong], [y_zhezhong, y_zheqi], color='green')
    return x_break_break, y_zheqi, y_zhezhong


# 从塔顶开始作阶梯
x_break = ls_XA[2]
js2 = 0
stage_positions = []

while x_break > x_intersection and js2 < 50:  # 防止无限循环
    x_break, y_zhezhong = bufenhuiliu_caozuo_zhexian(x_break)
    js2 += 1
    stage_positions.append(('精馏段', x_break, y_zhezhong))
    if x_break <= x_intersection:
        break

# 记录进料板位置
feed_stage = js2
print(f"进料板位置: 第 {feed_stage} 块理论板")

# 进入提馏段
while x_break > ls_XA[3] + 0.001 and js2 < 100:  # 容差
    x_break, y_zheqi, y_zhezhong = bufenhuiliu_tiliu_zhexian(x_break)
    js2 += 1
    stage_positions.append(('提馏段', x_break, y_zhezhong))
    if x_break <= ls_XA[3] + 0.001:
        break

if js2 >= 100:
    print("警告：可能未收敛到塔釜组成")

print(f'\n理论塔板数分析:')
print(f'总理论板数（含塔釜）: {js2} 块')
print(f'实际理论板数（不含塔釜）: {js2 - 1} 块')
print(f'其中：')
print(f'  精馏段: {feed_stage} 块')
print(f'  提馏段: {js2 - feed_stage - 1} 块')
print(f'进料板位置: 第 {feed_stage} 块理论板')

plt.title(f'苯-甲苯精馏塔设计')
plt.xlabel('液相组成 x (苯摩尔分数)')
plt.ylabel('气相组成 y (苯摩尔分数)')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()



def safe_ceil(value):
    """安全的向上取整，处理可能的数值问题"""
    return math.ceil(value) if isinstance(value, (int, float)) else value


#  精馏段实际塔板数
N1 = feed_stage / ET
#  提馏段实际塔板数
N2 = (js2 - feed_stage - 1) / ET
print(f'精馏段实际塔板数:{N1:.2f}块，即{safe_ceil(N1):.0f}块')
print(f'提馏段实际塔板数:{N2:.2f}块，即{safe_ceil(N2):.0f}块')

plot_detailed_analysis = False
'''
if plot_detailed_analysis:
    # 设置设计参数（完整版）
    enhanced_calc.set_design_parameters(
        xD=x_D, xF=x_F, xW=x_W, q=q_calculated,
        R_min=R_min, R_multiplier=R_multiplier,
        R_actual=R_actual, R_internal=R_internal,
        reflux_temp=t_ewflux,
        feed_stage=feed_stage,  # 从主程序的计算获得
        total_stages=js2 - 1,  # 总理论板数（不含塔釜）
        ET=ET,  # 塔板效率
        stage_positions=stage_positions  # 传递阶梯图结果
    )

    # 绘制增强图，传递阶梯图结果
    fig, axes, R_min_calc, R_actual_calc, R_internal_calc = enhanced_calc.plot_enhanced_diagram(
        xD=x_D, xF=x_F, xW=x_W, q=q_calculated,
        R_multiplier=R_multiplier,
        reflux_temp=t_ewflux,
        R_internal=R_internal,
        custom_stage_positions=stage_positions,  # 传递阶梯图结果
        save_path="苯-甲苯过冷回流分析.png"
    )
    plt.show()
'''
# ===== 课设报告计算过程总结 =====
# ===== 添加塔釜温度计算 =====
print("\n补充: 塔釜温度计算")
print("-" * 50)
# 通过插值计算塔釜温度
t_bottom = chazhi(x_W, dic_txy, 1)  # 使用chazhi函数从txy数据中插值
print(f"通过t-x-y数据插值计算塔釜温度:")
print(f"  x_W = {x_W:.4f} 对应的泡点温度 t_bottom = {t_bottom:.1f}°C")

# ===== 课设报告计算过程总结 =====
print("\n" + "=" * 80)
print("精馏塔设计计算过程总结")
print("=" * 80)

# 1. 物料衡算
print("\n1. 物料衡算")
print("-" * 50)
print("已知条件:")
print(f"  年处理量: {Annual_processing_capacity} t/a")
print(f"  年工作时间: {work_time} h/a")
print(f"  进料苯质量分数: {wt_f}")
print(f"  塔顶苯质量分数: {D_wt}")
print(f"  塔釜苯质量分数: {F_wt}")

print("\n(1) 质量流量计算:")
print(f"  F_mass = 年处理量 / 工作时间 = {Annual_processing_capacity * 1000} kg / {work_time} h")
print(f"        = {F_kg_per_h:.2f} kg/h")

print("\n(2) 质量分数转换为摩尔分数:")
print(f"  公式: x_苯 = (w_苯/M_苯) / (w_苯/M_苯 + w_甲苯/M_甲苯)")
print(f"  其中: M_苯 = {MA_benzene:.2f} g/mol, M_甲苯 = {MA_toluene:.2f} g/mol")
print(
    f"  进料: x_F = ({wt_f}/{MA_benzene:.2f}) / ({wt_f}/{MA_benzene:.2f} + {(1 - wt_f)}/{MA_toluene:.2f})")
print(f"        = {x_F:.4f}")
print(f"  塔顶: x_D = {D_wt}/{MA_benzene:.2f} / ({D_wt}/{MA_benzene:.2f} + {(1 - D_wt)}/{MA_toluene:.2f})")
print(f"        = {x_D:.4f}")
print(f"  塔釜: x_W = {F_wt}/{MA_benzene:.2f} / ({F_wt}/{MA_benzene:.2f} + {(1 - F_wt)}/{MA_toluene:.2f})")
print(f"        = {x_W:.4f}")

print("\n(3) 进料平均分子量:")
print(f"  M_F = x_F × M_苯 + (1-x_F) × M_甲苯")
print(f"      = {x_F:.4f} × {MA_benzene:.2f} + {1 - x_F:.4f} × {MA_toluene:.2f}")
print(f"      = {M_F:.2f} g/mol = {M_F / 1000:.4f} kg/mol")

print("\n(4) 进料摩尔流量:")
print(f"  F_kmol = F_mass / M_F = {F_kg_per_h:.2f} kg/h / ({M_F / 1000:.4f} kg/mol)")
print(f"        = {F_kmol_per_h:.2f} kmol/h")

print("\n(5) 总物料衡算:")
print("  物料衡算方程:")
print("    F = D + W")
print("    F·x_F = D·x_D + W·x_W")

# 计算D和W（摩尔流量）
D_kmol = (F_kmol_per_h * (x_F - x_W)) / (x_D - x_W)
W_kmol = F_kmol_per_h - D_kmol

print(f"\n  计算塔顶和塔釜摩尔流量:")
print(f"  由物料衡算方程联立求解:")
print(f"    D = F·(x_F - x_W) / (x_D - x_W)")
print(f"      = {F_kmol_per_h:.4f}×({x_F:.4f} - {x_W:.4f}) / ({x_D:.4f} - {x_W:.4f})")
print(f"      = {F_kmol_per_h * (x_F - x_W):.4f} / {x_D - x_W:.4f}")
print(f"      = {D_kmol:.4f} kmol/h")
print(f"    W = F - D = {F_kmol_per_h:.4f} - {D_kmol:.4f}")
print(f"      = {W_kmol:.4f} kmol/h")

print("\n(6) 质量流量换算:")
# 计算各物流的平均分子量
M_D = x_D * MA_benzene + (1 - x_D) * MA_toluene
M_W = x_W * MA_benzene + (1 - x_W) * MA_toluene

print(f"  塔顶平均分子量:")
print(f"    M_D = x_D·M_苯 + (1-x_D)·M_甲苯")
print(f"        = {x_D:.4f}×{MA_benzene:.2f} + {1 - x_D:.4f}×{MA_toluene:.2f}")
print(f"        = {M_D:.2f} g/mol")

print(f"  塔釜平均分子量:")
print(f"    M_W = x_W·M_苯 + (1-x_W)·M_甲苯")
print(f"        = {x_W:.4f}×{MA_benzene:.2f} + {1 - x_W:.4f}×{MA_toluene:.2f}")
print(f"        = {M_W:.2f} g/mol")

# 计算质量流量 (kg/h)
F_mass = F_kg_per_h
D_mass = D_kmol * M_D  # 注意单位: kmol/h * g/mol = kg/h
W_mass = W_kmol * M_W

print(f"\n  进料质量流量: F' = {F_mass:.2f} kg/h")
print(f"  塔顶质量流量: D' = D × M_D = {D_kmol:.4f} × {M_D:.2f}")
print(f"                = {D_mass:.2f} kg/h")
print(f"  塔釜质量流量: W' = W × M_W = {W_kmol:.4f} × {M_W:.2f}")
print(f"                = {W_mass:.2f} kg/h")
'''
print("\n(7) 物料衡算校验:")
total_mass = D_mass + W_mass
mass_error = abs(total_mass - F_mass) / F_mass * 100
print(f"  总质量平衡: D' + W' = {D_mass:.2f} + {W_mass:.2f} = {total_mass:.2f} kg/h")
print(f"              F' = {F_mass:.2f} kg/h")
print(f"  相对误差: {mass_error:.2f}% (< 0.1%，计算正确)")

# 苯组分衡算校验
benzene_F = F_mass * 进料组成_wt
benzene_D = D_mass * D_wt
benzene_W = W_mass * F_wt
total_benzene = benzene_D + benzene_W
benzene_error = abs(total_benzene - benzene_F) / benzene_F * 100

print(f"\n  苯组分平衡:")
print(f"  进料中苯: F'×w_F = {F_mass:.2f}×{进料组成_wt} = {benzene_F:.2f} kg/h")
print(f"  塔顶中苯: D'×w_D = {D_mass:.2f}×{D_wt} = {benzene_D:.2f} kg/h")
print(f"  塔釜中苯: W'×w_W = {W_mass:.2f}×{F_wt} = {benzene_W:.2f} kg/h")
# print(f"  苯总量: {benzene_D:.2f} + {benzene_W:.2f} = {total_benzene:.2f} kg/h")
# print(f"  相对误差: {benzene_error:.2f}% (< 0.1%，计算正确)")
'''
print(f"\n物料衡算结果汇总:")
print(f"  摩尔流量: F = {F_kmol_per_h:.4f} kmol/h, D = {D_kmol:.4f} kmol/h, W = {W_kmol:.4f} kmol/h")
print(f"  质量流量: F' = {F_mass:.2f} kg/h, D' = {D_mass:.2f} kg/h, W' = {W_mass:.2f} kg/h")
print(f"  塔顶采出率: D/F = {D_kmol / F_kmol_per_h * 100:.2f}% (摩尔), D'/F' = {D_mass / F_mass * 100:.2f}% (质量)")
# 2. 温度计算
print("=" * 50)
print('2.2理论板数N_T的计算')
print("=" * 50)
print("\n\n2.2.1 温度计算")
print("-" * 50)

print("\n(1) 进料泡点温度计算:")
print(f"  由t-x-y数据，通过线性插值求x_F={x_F:.4f}对应的泡点温度:")

# 获取实际的插值点
x1_F, x2_F, t1_F, t2_F = find_interpolation_points(x_F, dic_txy, value_index=1)
print(f"  查txy表，在x={x1_F:.2f}(t={t1_F:.2f}°C)和x={x2_F:.2f}(t={t2_F:.2f}°C)之间插值:")
if x1_F != x2_F:
    print(f"  tBP_F = {t1_F:.2f} + ({t2_F:.2f}-{t1_F:.2f})×({x_F:.4f}-{x1_F:.2f})/({x2_F:.2f}-{x1_F:.2f})")
else:
    print(f"  x_F={x_F:.4f}在数据表中，tBP_F = {t1_F:.2f}°C")
print(f"        = {tBP_F:.1f} °C")

print("\n(2) 塔顶温度计算:")
print(f"  由t-x-y数据，通过线性插值求x_D={x_D:.4f}对应的泡点温度:")

# 获取实际的插值点
x1_D, x2_D, t1_D, t2_D = find_interpolation_points(x_D, dic_txy, value_index=1)
print(f"  查txy表，在x={x1_D:.2f}(t={t1_D:.2f}°C)和x={x2_D:.2f}(t={t2_D:.2f}°C)之间插值:")
if x1_D != x2_D:
    print(f"  t_top = {t1_D:.2f} + ({t2_D:.2f}-{t1_D:.2f})×({x_D:.4f}-{x1_D:.2f})/({x2_D:.2f}-{x1_D:.2f})")
else:
    print(f"  x_D={x_D:.4f}在数据表中，t_top = {t1_D:.2f}°C")
print(f"        = {t_top:.1f} °C")

print("\n(3) 塔釜温度计算:")
print(f"  由t-x-y数据，通过线性插值求x_W={x_W:.4f}对应的泡点温度:")

# 获取实际的插值点
x1_W, x2_W, t1_W, t2_W = find_interpolation_points(x_W, dic_txy, value_index=1)
print(f"  查txy表，在x={x1_W:.2f}(t={t1_W:.2f}°C)和x={x2_W:.2f}(t={t2_W:.2f}°C)之间插值:")
if x1_W != x2_W:
    print(f"  t_bottom = {t1_W:.2f} + ({t2_W:.2f}-{t1_W:.2f})×({x_W:.4f}-{x1_W:.2f})/({x2_W:.2f}-{x1_W:.2f})")
else:
    print(f"  x_W={x_W:.4f}在数据表中，t_bottom = {t1_W:.2f}°C")
print(f"           = {t_bottom:.1f} °C")

# 3. 进料热状况参数q计算
print("\n\n2.2.2 进料热状况参数q计算")
print("-" * 50)
print(f"已知: 进料温度tF = {t_f}°C, 泡点温度tBP_F = {tBP_F:.1f}°C")
if tF < tBP_F:
    print(f"由于tF < tBP_F，为冷液进料")

    print("\n(1) 平均温度计算:")
    print(f"  t_avg_q = (tF + tBP_F)/2 = ({t_f} + {tBP_F:.1f})/2")
    print(f"          = {t_avg:.1f}°C")

    print("\n(2) 比热容插值计算:")
    print(f"  在t_avg_q = {t_avg:.1f}°C时:")
    Cp_benzene_val = chazhi(t_avg, Cp_benzene, 'n')
    Cp_toluene_val = chazhi(t_avg, Cp_toluene, 'n')
    print(f"  苯的比热容: 查Cp_benzene表，在{int(t_avg / 20) * 20}°C和{(int(t_avg / 20) + 1) * 20}°C之间插值")
    print(f"             得 Cp_苯 = {Cp_benzene_val:.3f} kJ/(kg·℃)")
    print(f"  甲苯的比热容: 查Cp_toluene表，在{int(t_avg / 20) * 20}°C和{(int(t_avg / 20) + 1) * 20}°C之间插值")
    print(f"               得 Cp_甲苯 = {Cp_toluene_val:.3f} kJ/(kg·℃)")

    print("\n(3) 摩尔比热容计算:")
    Cp_benzene_mol_val = Cp_benzene_val * MA_benzene
    Cp_toluene_mol_val = Cp_toluene_val * MA_toluene
    print(f"  Cp_苯_mol = Cp_苯 × M_苯 = {Cp_benzene_val:.3f} × {MA_benzene:.2f}")
    print(f"            = {Cp_benzene_mol_val:.2f} kJ/(kmol·℃)")
    print(f"  Cp_甲苯_mol = Cp_甲苯 × M_甲苯 = {Cp_toluene_val:.3f} × {MA_toluene:.2f}")
    print(f"              = {Cp_toluene_mol_val:.2f} kJ/(kmol·℃)")

    print("\n(4) 混合物平均摩尔热容:")
    Cpm = x_F * Cp_benzene_mol_val + (1 - x_F) * Cp_toluene_mol_val
    print(f"  Cpm = x_F × Cp_苯_mol + (1-x_F) × Cp_甲苯_mol")
    print(f"      = {x_F:.4f} × {Cp_benzene_mol_val:.2f} + {1 - x_F:.4f} × {Cp_toluene_mol_val:.2f}")
    print(f"      = {Cpm:.2f} kJ/(kmol·℃)")

    print("\n(5) 混合物汽化潜热:")
    r_benzene_mol = r_benzene * MA_benzene
    r_toluene_mol = r_toluene * MA_toluene
    rm = x_F * r_benzene_mol + (1 - x_F) * r_toluene_mol
    print(f"  r_苯_mol = r_苯 × M_苯 = {r_benzene} × {MA_benzene:.2f} = {r_benzene_mol:.0f} kJ/kmol")
    print(f"  r_甲苯_mol = r_甲苯 × M_甲苯 = {r_toluene} × {MA_toluene:.2f} = {r_toluene_mol:.0f} kJ/kmol")
    print(f"  rm = x_F × r_苯_mol + (1-x_F) × r_甲苯_mol")
    print(f"     = {x_F:.4f} × {r_benzene_mol:.0f} + {1 - x_F:.4f} × {r_toluene_mol:.0f}")
    print(f"     = {rm:.0f} kJ/kmol")

    print("\n(6) 进料热状况参数q:")
    q = 1 + (Cpm * (tBP_F - tF)) / rm
    print(f"  q = 1 + (Cpm × (tBP_F - tF)) / rm")
    print(f"    = 1 + ({Cpm:.2f} × ({tBP_F:.1f} - {t_f})) / {rm:.0f}")
    print(f"    = 1 + ({Cpm:.2f} × {tBP_F - t_f:.1f}) / {rm:.0f}")
    print(f"    = {q:.3f}")
else:
    print(f"由于tF ≥ tBP_F，为饱和液体进料，q = 1.0")

if q != 1:
    q_slope = q / (q - 1)
    q_intercept = -x_F / (q - 1)
    print(f"  q线: y = q/(q-1) × x - x_F/(q-1)")
    print(f"       = {q:.3f}/{q - 1:.3f} × x - {x_F:.4f}/{q - 1:.3f}")
    print(f"       = {q_slope:.3f}x + {q_intercept:.3f}")
else:
    print(f"  q线: x = {x_F:.4f} (垂直线)")

# 5. 回流比计算
print("\n\n2.2.3 回流比计算")
print("-" * 50)
print("(1) 最小回流比计算:")
print("程序可以给出最小回流比，但这是经过计算得出的。\n一般来说应该绘制x-y图和q线，取得q线与曲线的交点，即可继续算出R_min")
# 对交点坐标保留三位小数
x_eq_rounded = round(intersection[0], 3)
y_eq_rounded = round(intersection[1], 3)
x_D_rounded=round(x_D,4)
# 使用保留三位小数的交点计算最小回流比
numerator = x_D_rounded - y_eq_rounded
denominator = y_eq_rounded - x_eq_rounded
R_min_rounded = numerator / denominator
print(f"  交点坐标: x_eq = {x_eq_rounded:.4f}, y_eq = {y_eq_rounded:.4f}")
print(f"  R_min = (x_D - y_eq) / (y_eq - x_eq)")
print(f"        = ({x_D:.4f} - {y_eq_rounded:.3f}) / ({y_eq_rounded:.3f} - {x_eq_rounded:.3f})")
print(f"        = {numerator:.4f} / {denominator:.4f}")
print(f"        = {R_min_rounded:.4f}")

print("\n(2) 理论外回流比:")
print(f"  R_theoretical = R_min × R_multiplier")
print(f"                = {R_min:.4f} × {R_multiplier}")
print(f"                = {R_actual:.4f}")

print("\n(3) 过冷回流修正(内回流比):")
print(f"  已知: 回流温度 = {t_ewflux}°C, 塔顶温度 = {t_top:.1f}°C")
delta_T = t_top - t_ewflux
print(f"  温差 ΔT = t_top - T_reflux = {t_top:.1f} - {t_ewflux} = {delta_T:.1f} K")

# 重新计算过冷回流修正系数（使用塔顶组成x_D）
T_avg_reflux = (t_ewflux + t_top) / 2
print(f"  平均温度 T_avg_reflux = (回流温度 + t_top)/2 = ({t_ewflux} + {t_top:.1f})/2 = {T_avg_reflux:.1f}°C")

# 计算塔顶组成的物性
Cp_benzene_reflux = chazhi(T_avg_reflux, Cp_benzene, 'n')
Cp_toluene_reflux = chazhi(T_avg_reflux, Cp_toluene, 'n')
print(f"  苯的比热容(在{T_avg_reflux:.1f}°C): Cp_苯' = {Cp_benzene_reflux:.3f} kJ/(kg·K)")
print(f"  甲苯的比热容(在{T_avg_reflux:.1f}°C): Cp_甲苯' = {Cp_toluene_reflux:.3f} kJ/(kg·K)")

Cp_benzene_mol_reflux = Cp_benzene_reflux * MA_benzene
Cp_toluene_mol_reflux = Cp_toluene_reflux * MA_toluene
print(
    f"  苯的摩尔比热容: Cp_苯_mol' = {Cp_benzene_reflux:.3f} × {MA_benzene:.2f} = {Cp_benzene_mol_reflux:.2f} kJ/(kmol·K)")
print(
    f"  甲苯的摩尔比热容: Cp_甲苯_mol' = {Cp_toluene_reflux:.3f} × {MA_toluene:.2f} = {Cp_toluene_mol_reflux:.2f} kJ/(kmol·K)")

Cp_mix = x_D * Cp_benzene_mol_reflux + (1 - x_D) * Cp_toluene_mol_reflux
print(f"  混合物摩尔比热容(塔顶组成):")
print(f"  Cp_mix = x_D × Cp_苯_mol' + (1-x_D) × Cp_甲苯_mol'")
print(f"         = {x_D:.4f} × {Cp_benzene_mol_reflux:.2f} + {1 - x_D:.4f} × {Cp_toluene_mol_reflux:.2f}")
print(f"         = {Cp_mix:.1f} kJ/(kmol·K)")

# 汽化潜热
r_benzene_mol = r_benzene * MA_benzene
r_toluene_mol = r_toluene * MA_toluene
r_mix = x_D * r_benzene_mol + (1 - x_D) * r_toluene_mol
print(f"\n  苯的摩尔汽化潜热: ΔH_苯 = {r_benzene} × {MA_benzene:.2f} = {r_benzene_mol:.0f} kJ/kmol")
print(f"  甲苯的摩尔汽化潜热: ΔH_甲苯 = {r_toluene} × {MA_toluene:.2f} = {r_toluene_mol:.0f} kJ/kmol")
print(f"  混合物汽化潜热(塔顶组成):")
print(f"  ΔHv = x_D × ΔH_苯 + (1-x_D) × ΔH_甲苯")
print(f"      = {x_D:.4f} × {r_benzene_mol:.0f} + {1 - x_D:.4f} × {r_toluene_mol:.0f}")
print(f"      = {r_mix:.0f} kJ/kmol")

print(f"\n  修正系数 = 1 + (Cp_mix × ΔT) / ΔHv")
print(f"           = 1 + ({Cp_mix:.1f} × {delta_T:.1f}) / {r_mix:.0f}")
print(f"           = 1 + {Cp_mix * delta_T:.0f} / {r_mix:.0f}")
correction_calc = 1 + (Cp_mix * delta_T) / r_mix
print(f"           = {correction_calc:.4f}")

print(f"\n  实际内回流比: R' = R_theoretical × 修正系数")
print(f"                = {R_actual:.4f} × {correction_calc:.4f}")
R_internal_calc = R_actual * correction_calc
print(f"                = {R_internal_calc:.4f}")

print(f"\n(4) 设计采用回流比:")
print(f"  R = R' = {R:.4f} (使用内回流比进行后续计算)")

# 6. 理论塔板数计算
print("\n\n2.2.4 理论塔板数计算(阶梯图解法)")
print("-" * 50)

print("(1) 操作线方程:")
slope = R / (R + 1)
intercept = x_D / (R + 1)
print(f"  精馏段: y = R/(R+1) × x + x_D/(R+1)")
print(f"          = {R:.4f}/{R + 1:.4f} × x + {x_D:.4f}/{R + 1:.4f}")
print(f"          = {slope:.4f}x + {intercept:.4f}")

if q != 1:
    q_slope = q / (q - 1)
    q_intercept = -x_F / (q - 1)
    print(f"  q线: y = q/(q-1) × x - x_F/(q-1)")
    print(f"       = {q:.3f}/{q - 1:.3f} × x - {x_F:.4f}/{q - 1:.3f}")
    print(f"       = {q_slope:.3f}x + {q_intercept:.3f}")
else:
    print(f"  q线: x = {x_F:.4f} (垂直线)")

print(f"  提馏段: 连接点({x_W:.4f}, {x_W:.4f})和q线与精馏段操作线交点")

print("\n(2) 阶梯图解法结果:")
print(f"  总理论板数(含塔釜): {js2} 块")
print(f"  实际理论板数(不含塔釜): {js2 - 1} 块")
print(f"  其中: 精馏段 {feed_stage} 块")
print(f"        提馏段 {js2 - feed_stage - 1} 块")
print(f"  进料板位置: 第 {feed_stage} 块理论板")

# 7. 实际塔板数计算
print('\n\n'+"=" * 50)
print("2.3 实际塔板数计算")
print("=" * 50)

# 4. 塔板效率计算
print("\n\n2.3.1. 塔板效率计算")
print("-" * 50)

print("(1) 平均温度下粘度插值:")
print(f"  在t_avg = {t_avg2:.1f}°C时:")
m_A_val = chazhi(t_avg2, m_benzene, 'n')
m_B_val = chazhi(t_avg2, m_toluene, 'n')
print(f"  苯的粘度: 查m_benzene表，在{int(t_avg2 / 20) * 20}°C和{(int(t_avg2 / 20) + 1) * 20}°C之间插值")
print(f"           得 μ_苯 = {m_A_val:.3f} mPa·s")
print(f"  甲苯的粘度: 查m_toluene表，在{int(t_avg2 / 20) * 20}°C和{(int(t_avg2 / 20) + 1) * 20}°C之间插值")
print(f"             得 μ_甲苯 = {m_B_val:.3f} mPa·s")

print("\n(2) 混合物粘度:")
m_m = x_F * m_A_val + (1 - x_F) * m_B_val
print(f"  μ_m = x_F × μ_苯 + (1-x_F) × μ_甲苯")
print(f"      = {x_F:.4f} × {m_A_val:.3f} + {1 - x_F:.4f} × {m_B_val:.3f}")
print(f"      = {m_m:.3f} mPa·s")

print("\n(3) 塔板效率(Drickamer和Brabford法):")
ET = 0.17 - 0.616 * np.log10(m_m)
print(f"  E_T = 0.17 - 0.616 × log10(μ_m)")
print(f"      = 0.17 - 0.616 × log10({m_m:.3f})")
print(f"      = 0.17 - 0.616 × {np.log10(m_m):.3f}")
print(f"      = {ET:.4f} = {ET * 100:.1f}%")

print("\n\n2.3.2. 实际塔板数计算")
print("-" * 50)
print("\n(1) 精馏段实际塔板数:")
N1 = feed_stage / ET
print(f"  N1 = 精馏段理论板数 / E_T = {feed_stage} / {ET:.4f}")
print(f"     = {feed_stage / ET:.2f} 块")
print(f"     ≈ {safe_ceil(N1):.0f} 块 (向上取整)")

print("\n(2) 提馏段实际塔板数:")
N2 = (js2 - feed_stage - 1) / ET
print(f"  N2 = 提馏段理论板数 / E_T = {js2 - feed_stage - 1} / {ET:.4f}")
print(f"     = {N2:.2f} 块")
print(f"     ≈ {safe_ceil(N2):.0f} 块 (向上取整)")

print("\n(3) 总实际塔板数:")
total_plates = safe_ceil(N1) + safe_ceil(N2)
print(f"  N_total = N1 + N2 = {safe_ceil(N1):.0f} + {safe_ceil(N2):.0f}")
print(f"          = {total_plates:.0f} 块")
'''
# 8. 设计参数汇总
print("\n\n8. 设计参数汇总")
print("-" * 50)
print(f"  (1) 进料条件:")
print(f"      - 处理量: {F_kg_per_h:.2f} kg/h ({F_kmol_per_h:.2f} kmol/h)")
print(f"      - 组成: x_F = {x_F:.4f} (苯摩尔分数)")
print(f"      - 温度: {进料温度}°C (q = {q:.3f})")

print(f"\n  (2) 产品规格:")
print(f"      - 塔顶: x_D = {x_D:.4f}, 温度 = {t_top:.1f}°C")
print(f"      - 塔釜: x_W = {x_W:.4f}, 温度 = {t_bottom:.1f}°C")

print(f"\n  (3) 操作参数:")
print(f"      - 最小回流比: R_min = {R_min:.3f}")
print(f"      - 设计外回流比: R_theoretical = {R_actual:.3f} (理论值)")
print(f"      - 实际内回流比: R' = {R_internal:.3f} (考虑过冷回流)")
print(f"      - 回流比倍数: R_multiplier = {R_multiplier}")
print(f"      - 塔板效率: E_T = {ET * 100:.1f}%")

print(f"\n  (4) 塔板数:")
print(f"      - 理论板数: {js2 - 1} 块 (不含塔釜)")
print(f"      - 精馏段: {feed_stage} 块理论板")
print(f"      - 提馏段: {js2 - feed_stage - 1} 块理论板")
print(f"      - 实际板数: {total_plates:.0f} 块")
print(f"      - 进料板位置: 第 {feed_stage} 块理论板")

print(f"\n  (5) 温度参数:")
print(f"      - 进料泡点: {tBP_F:.1f}°C")
print(f"      - 塔顶温度: {t_top:.1f}°C")
print(f"      - 塔釜温度: {t_bottom:.1f}°C")
print(f"      - 平均温度(物性查询): {t_avg2:.1f}°C")
print(f"      - 温度范围: ΔT = {abs(t_top - t_bottom):.1f}°C")
'''
print("\n" + "=" * 80)
print("计算结束 - 以上结果可用于课程设计报告编写")
print("=" * 80)

# ===== 保存计算结果到JSON文件 =====


# 创建要传递的数据字典
calculation_results = {
    # 物料衡算
    "F_kmol_per_h": float(F_kmol_per_h),
    "D_kmol": float(D_kmol),
    "W_kmol": float(W_kmol),
    "F_mass": float(F_kg_per_h),
    "D_mass": float(D_mass),
    "W_mass": float(W_mass),

    # 组成
    "x_F": float(x_F),
    "x_D": float(x_D),
    "x_W": float(x_W),
    "w_F": float(wt_f),
    "w_D": float(D_wt),
    "w_W": float(F_wt),

    # 温度
    "t_top": float(t_top),
    "t_bottom": float(t_bottom),
    "t_feed_bubble": float(tBP_F),
    "t_feed": float(t_f),

    # 操作参数
    "R_min": float(R_min),
    "R_actual": float(R_actual),
    "R_internal": float(R_internal),
    "R": float(R),  # 使用的内回流比

    # 塔板数
    "N_theoretical_total": int(js2 - 1),
    "N_theoretical_rectifying": int(feed_stage),
    "N_theoretical_stripping": int(js2 - feed_stage - 1),
    "ET": float(ET),
    "N1_actual": float(N1),
    "N2_actual": float(N2),

    # 其他参数
    "feed_stage_position": int(feed_stage),
    "q_value": float(q),
    "average_viscosity": float(m_m)
}

# 保存到文件
with open('benzene_toluene_results.json', 'w', encoding='utf-8') as f:
    json.dump(calculation_results, f, indent=4, ensure_ascii=False)

print(f"\n计算结果已保存到 benzene_toluene_results.json")
print("可在下一个程序中使用这些数据")
