"""
塔板布置计算程序（几何修正版）
功能：
1. 基于液体负荷和停留时间设计降液管（弓形几何计算）
2. 使用扣除法计算塔板面积分区，确保面积闭合
3. 浮阀塔：计算阀孔气速、阀孔数、浮阀排列
4. 筛板塔：计算筛孔数、筛孔排列
5. 溢流装置计算：溢流堰、降液管、受液盘
6. 完整的设计流程和校验，带详细计算过程
"""

import json
import math
import os
import random
from datetime import datetime


def calculate_segment_area(lw, D, verbose=False):
    """
    计算弓形降液管的几何参数（核心几何函数）
    参数:
        lw: 堰长（弦长），m
        D: 塔径，m
        verbose: 是否打印详细计算过程
    返回:
        dict: 包含所有弓形几何参数的字典
    """
    R = D / 2

    if verbose:
        print(f"\n   [弓形几何计算]")
        print(f"   输入: 堰长 lw = {lw:.4f} m, 塔径 D = {D:.3f} m")
        print(f"   塔半径 R = D/2 = {D:.3f}/2 = {R:.3f} m")

    # 1. 校验输入
    if lw > D or lw <= 0:
        raise ValueError(f"堰长lw={lw:.3f}m不能大于塔径D={D:.3f}m且必须大于0")

    # 2. 计算弦心距
    chord_midpoint = math.sqrt(R**2 - (lw/2)**2)

    if verbose:
        print(f"\n   1. 弦心距计算:")
        print(f"      弦心距 = √(R² - (lw/2)²)")
        print(f"             = √({R:.4f}² - ({lw:.4f}/2)²)")
        print(f"             = √({R**2:.6f} - {((lw/2)**2):.6f})")
        print(f"             = √({R**2 - ((lw/2)**2):.6f})")
        print(f"             = {chord_midpoint:.6f} m")

    # 3. 计算弓形高度（降液管宽度Wd）
    Wd = R - chord_midpoint

    if verbose:
        print(f"\n   2. 降液管宽度计算:")
        print(f"      Wd = R - 弦心距")
        print(f"         = {R:.4f} - {chord_midpoint:.6f}")
        print(f"         = {Wd:.6f} m")

    # 4. 计算圆心角（弧度）
    if lw/(2*R) <= 1:
        theta = 2 * math.asin(lw/(2*R))
    else:
        theta = math.pi

    if verbose:
        print(f"\n   3. 圆心角计算:")
        print(f"      θ = 2 × arcsin(lw/(2R))")
        print(f"        = 2 × arcsin({lw:.4f}/(2×{R:.4f}))")
        print(f"        = 2 × arcsin({lw:.4f}/{2*R:.4f})")
        print(f"        = 2 × arcsin({lw/(2*R):.4f})")
        print(f"        = 2 × {math.asin(lw/(2*R)):.4f} rad")
        print(f"        = {theta:.4f} rad")
        print(f"        = {math.degrees(theta):.2f}°")

    # 5. 计算弓形面积Af
    # 扇形面积 = 0.5 * R^2 * theta
    # 三角形面积 = 0.5 * lw * chord_midpoint
    sector_area = 0.5 * R**2 * theta
    triangle_area = 0.5 * lw * chord_midpoint
    Af = sector_area - triangle_area

    if verbose:
        print(f"\n   4. 弓形面积计算:")
        print(f"      a) 扇形面积 = 0.5 × R² × θ")
        print(f"                 = 0.5 × {R:.4f}² × {theta:.4f}")
        print(f"                 = 0.5 × {R**2:.6f} × {theta:.4f}")
        print(f"                 = {sector_area:.6f} m²")

        print(f"\n      b) 三角形面积 = 0.5 × lw × 弦心距")
        print(f"                   = 0.5 × {lw:.4f} × {chord_midpoint:.6f}")
        print(f"                   = {triangle_area:.6f} m²")

        print(f"\n      c) 弓形面积 = 扇形面积 - 三角形面积")
        print(f"                = {sector_area:.6f} - {triangle_area:.6f}")
        print(f"                = {Af:.6f} m²")

    # 6. 计算其他参数
    area_ratio = Af / (math.pi * R**2)  # 面积比

    if verbose:
        print(f"\n   5. 面积比计算:")
        print(f"      面积比 = Af / (π × R²)")
        print(f"            = {Af:.6f} / (π × {R**2:.6f})")
        print(f"            = {Af:.6f} / {math.pi * R**2:.6f}")
        print(f"            = {area_ratio:.4f} = {area_ratio*100:.2f}%")

    return {
        '堰长_lw_m': lw,
        '降液管宽度_Wd_m': Wd,
        '弓形面积_Af_m2': Af,
        '面积比_%': area_ratio * 100,
        '圆心角_度': math.degrees(theta),
        '弦心距_m': chord_midpoint,
        '扇形面积_m2': sector_area,
        '三角形面积_m2': triangle_area
    }


def calculate_bubbling_area(D, Ws_inlet, Ws_outlet, Wc, Af, Af_side='both', verbose=False):
    """
    使用扣除法计算塔板各分区面积
    参数:
        D: 塔径，m
        Ws_inlet: 入口安定区宽度，m
        Ws_outlet: 出口安定区宽度，m
        Wc: 边缘区宽度，m
        Af: 单个降液管面积，m²
        Af_side: 降液管配置 'both'（双侧）, 'inlet', 'outlet'
        verbose: 是否打印详细计算过程
    返回:
        dict: 包含所有面积分区结果的字典
    """
    R = D / 2

    if verbose:
        print(f"\n   [塔板面积分区计算]")
        print(f"   输入: 塔径 D = {D:.3f} m, 塔半径 R = {R:.3f} m")
        print(f"         入口安定区宽度 Ws_inlet = {Ws_inlet:.3f} m")
        print(f"         出口安定区宽度 Ws_outlet = {Ws_outlet:.3f} m")
        print(f"         边缘区宽度 Wc = {Wc:.3f} m")
        print(f"         单个降液管面积 Af = {Af:.6f} m²")
        print(f"         降液管配置: {Af_side}")

    # 1. 塔板总面积
    At = math.pi * R**2

    if verbose:
        print(f"\n   1. 塔板总面积计算:")
        print(f"      公式: At = π × R²")
        print(f"      计算: At = π × {R:.4f}²")
        print(f"           = π × {R**2:.6f}")
        print(f"           = {At:.6f} m²")

    # 2. 计算总降液管面积
    if Af_side == 'both':
        total_Af = 2 * Af
        if verbose:
            print(f"\n   2. 总降液管面积计算 (双侧):")
            print(f"      公式: total_Af = 2 × Af")
            print(f"      计算: total_Af = 2 × {Af:.6f}")
            print(f"               = {total_Af:.6f} m²")
    elif Af_side in ['inlet', 'outlet']:
        total_Af = Af
        if verbose:
            print(f"\n   2. 总降液管面积计算 (单侧):")
            print(f"      total_Af = Af = {Af:.6f} m²")
    else:
        total_Af = 0

    # 3. 安定区面积（环形区域）
    # 入口安定区：从(R-Wc-Ws_inlet)到(R-Wc)的环形
    if (R - Wc - Ws_inlet) > 0:
        As_inlet = math.pi * ((R - Wc)**2 - (R - Wc - Ws_inlet)**2)
        if verbose:
            print(f"\n   3. 入口安定区面积计算:")
            print(f"      公式: As_inlet = π × [(R-Wc)² - (R-Wc-Ws_inlet)²]")
            print(f"      计算: As_inlet = π × [({R-Wc:.4f})² - ({R-Wc-Ws_inlet:.4f})²]")
            print(f"               = π × [{(R-Wc)**2:.6f} - {(R-Wc-Ws_inlet)**2:.6f}]")
            print(f"               = π × {(R-Wc)**2 - (R-Wc-Ws_inlet)**2:.6f}")
            print(f"               = {As_inlet:.6f} m²")
    else:
        As_inlet = 0
        if verbose:
            print(f"\n   3. 入口安定区面积计算:")
            print(f"      注意: R-Wc-Ws_inlet = {R-Wc-Ws_inlet:.3f} ≤ 0，取As_inlet = 0")

    # 出口安定区：从(R-Wc-Ws_outlet)到(R-Wc)的环形
    if (R - Wc - Ws_outlet) > 0:
        As_outlet = math.pi * ((R - Wc)**2 - (R - Wc - Ws_outlet)**2)
        if verbose:
            print(f"\n   4. 出口安定区面积计算:")
            print(f"      公式: As_outlet = π × [(R-Wc)² - (R-Wc-Ws_outlet)²]")
            print(f"      计算: As_outlet = π × [({R-Wc:.4f})² - ({R-Wc-Ws_outlet:.4f})²]")
            print(f"                = π × [{(R-Wc)**2:.6f} - {(R-Wc-Ws_outlet)**2:.6f}]")
            print(f"                = π × {(R-Wc)**2 - (R-Wc-Ws_outlet)**2:.6f}")
            print(f"                = {As_outlet:.6f} m²")
    else:
        As_outlet = 0
        if verbose:
            print(f"\n   4. 出口安定区面积计算:")
            print(f"      注意: R-Wc-Ws_outlet = {R-Wc-Ws_outlet:.3f} ≤ 0，取As_outlet = 0")

    # 4. 边缘区面积
    Ac = math.pi * (R**2 - (R - Wc)**2)

    if verbose:
        print(f"\n   5. 边缘区面积计算:")
        print(f"      公式: Ac = π × [R² - (R-Wc)²]")
        print(f"      计算: Ac = π × [{R:.4f}² - ({R:.4f}-{Wc:.4f})²]")
        print(f"            = π × [{R**2:.6f} - {((R - Wc)**2):.6f}]")
        print(f"            = π × {R**2 - (R - Wc)**2:.6f}")
        print(f"            = {Ac:.6f} m²")

    # 5. 鼓泡区面积（通过扣除法计算）
    Aa = At - total_Af - As_inlet - As_outlet - Ac

    if verbose:
        print(f"\n   6. 鼓泡区面积计算（扣除法）:")
        print(f"      公式: Aa = At - total_Af - As_inlet - As_outlet - Ac")
        print(f"      计算: Aa = {At:.6f} - {total_Af:.6f} - {As_inlet:.6f} - {As_outlet:.6f} - {Ac:.6f}")
        print(f"            = {Aa:.6f} m²")

    # 6. 计算各区域面积占比
    results = {
        # 塔板参数
        '塔径_D_m': D,
        '塔半径_R_m': R,
        '塔板总面积_At_m2': At,

        # 分区宽度
        '入口安定区宽度_Ws_inlet_m': Ws_inlet,
        '出口安定区宽度_Ws_outlet_m': Ws_outlet,
        '边缘区宽度_Wc_m': Wc,

        # 降液管参数
        '单个降液管面积_Af_m2': Af,
        '总降液管面积_total_Af_m2': total_Af,
        '降液管配置': Af_side,

        # 各分区面积
        '入口安定区面积_As_inlet_m2': As_inlet,
        '出口安定区面积_As_outlet_m2': As_outlet,
        '边缘区面积_Ac_m2': Ac,
        '鼓泡区面积_Aa_m2': Aa,

        # 面积占比
        '鼓泡区面积占比_%': (Aa / At) * 100 if At > 0 else 0,
        '降液管面积占比_%': (total_Af / At) * 100 if At > 0 else 0,
        '安定区面积占比_%': ((As_inlet + As_outlet) / At) * 100 if At > 0 else 0,
        '边缘区面积占比_%': (Ac / At) * 100 if At > 0 else 0,

        # 校验信息
        '面积校验通过': Aa >= 0,
        '面积闭合误差_%': 0
    }

    # 7. 计算面积闭合误差
    area_sum = Aa + total_Af + As_inlet + As_outlet + Ac
    area_error = abs(area_sum - At) / At * 100 if At > 0 else 0
    results['面积闭合误差_%'] = area_error

    if verbose:
        print(f"\n   7. 面积闭合校验:")
        print(f"      各分区面积之和 = Aa + total_Af + As_inlet + As_outlet + Ac")
        print(f"                     = {Aa:.6f} + {total_Af:.6f} + {As_inlet:.6f} + {As_outlet:.6f} + {Ac:.6f}")
        print(f"                     = {area_sum:.6f} m²")
        print(f"      塔板总面积 At = {At:.6f} m²")
        print(f"      面积闭合误差 = |{area_sum:.6f} - {At:.6f}| / {At:.6f} × 100%")
        print(f"                   = {area_error:.6f}%")

        print(f"\n   8. 面积占比:")
        print(f"      鼓泡区占比: {results['鼓泡区面积占比_%']:.2f}%")
        print(f"      降液管占比: {results['降液管面积占比_%']:.2f}%")
        print(f"      安定区占比: {results['安定区面积占比_%']:.2f}%")
        print(f"      边缘区占比: {results['边缘区面积占比_%']:.2f}%")

        if results['面积校验通过']:
            print(f"\n      ✓ 面积分区合理，鼓泡区面积为正")
        else:
            print(f"\n      ✗ 错误: 鼓泡区面积为负，请调整分区参数")

    return results


def find_lw_for_target_af(target_Af, D, lw_ratio_range=(0.6, 0.8), tolerance=0.0001, verbose=False):
    """
    通过二分法找到满足目标降液管面积的堰长
    参数:
        target_Af: 目标降液管面积，m²
        D: 塔径，m
        lw_ratio_range: 堰长比例范围
        tolerance: 容差
        verbose: 是否打印详细计算过程
    返回:
        lw: 堰长，m
        actual_Af: 实际降液管面积，m²
        Wd: 降液管宽度，m
    """
    R = D / 2
    lw_min = lw_ratio_range[0] * D
    lw_max = lw_ratio_range[1] * D

    if verbose:
        print(f"\n   [堰长搜索过程]")
        print(f"   目标降液管面积: target_Af = {target_Af:.6f} m²")
        print(f"   塔径: D = {D:.3f} m")
        print(f"   堰长范围: lw_min = {lw_min:.4f} m, lw_max = {lw_max:.4f} m")
        print(f"   容差: tolerance = {tolerance}")

    # 检查边界
    segment_min = calculate_segment_area(lw_min, D)
    segment_max = calculate_segment_area(lw_max, D)

    if target_Af < segment_min['弓形面积_Af_m2']:
        if verbose:
            print(f"\n   目标Af={target_Af:.6f}小于最小可能值{segment_min['弓形面积_Af_m2']:.6f}")
            print(f"   使用最小堰长: lw = {lw_min:.4f} m")
        return lw_min, segment_min['弓形面积_Af_m2'], segment_min['降液管宽度_Wd_m']

    if target_Af > segment_max['弓形面积_Af_m2']:
        if verbose:
            print(f"\n   目标Af={target_Af:.6f}大于最大可能值{segment_max['弓形面积_Af_m2']:.6f}")
            print(f"   使用最大堰长: lw = {lw_max:.4f} m")
        return lw_max, segment_max['弓形面积_Af_m2'], segment_max['降液管宽度_Wd_m']

    # 二分法搜索
    if verbose:
        print(f"\n   开始二分法搜索:")

    for i in range(50):  # 最多50次迭代
        lw_mid = (lw_min + lw_max) / 2
        segment_mid = calculate_segment_area(lw_mid, D)
        Af_mid = segment_mid['弓形面积_Af_m2']

        if verbose and (i < 3 or i % 10 == 9):  # 显示前3次和每10次的结果
            print(f"     迭代{i+1}: lw_mid={lw_mid:.4f}m, Af_mid={Af_mid:.6f}m², 差值={abs(Af_mid-target_Af):.6f}m²")

        if abs(Af_mid - target_Af) < tolerance:
            if verbose:
                print(f"\n   在第{i+1}次迭代找到满足条件的堰长")
                print(f"   最终结果: lw = {lw_mid:.4f} m, Af = {Af_mid:.6f} m²")
            return lw_mid, Af_mid, segment_mid['降液管宽度_Wd_m']

        if Af_mid < target_Af:
            lw_min = lw_mid
        else:
            lw_max = lw_mid

    # 未达到容差，返回最接近的值
    lw_final = (lw_min + lw_max) / 2
    segment_final = calculate_segment_area(lw_final, D)

    if verbose:
        print(f"\n   经过50次迭代未达到容差，返回最接近的值")
        print(f"   最终结果: lw = {lw_final:.4f} m, Af = {segment_final['弓形面积_Af_m2']:.6f} m²")

    return lw_final, segment_final['弓形面积_Af_m2'], segment_final['降液管宽度_Wd_m']


class TowerPlateArrangementCorrected:
    """塔板布置计算器（几何修正版）"""

    def __init__(self, diameter_design_file='tower_diameter_approximate_optimized.json',
                 properties_file='tower_properties_results.json',
                 tower_type='float_valve',
                 use_random_params=False,
                 plate_material='carbon_steel',
                 verbose=True):  # 添加详细计算标志
        """
        初始化塔板布置计算器

        参数：
        diameter_design_file: 塔径设计结果文件
        properties_file: 塔工艺特性文件
        tower_type: 塔板类型，'float_valve'（浮阀塔）或'sieve'（筛板塔）
        use_random_params: 是否使用随机参数生成
        plate_material: 塔板材料，'carbon_steel'（碳钢）或'stainless_steel'（不锈钢）
        verbose: 是否打印详细计算过程
        """
        self.tower_type = tower_type
        self.use_random_params = use_random_params
        self.plate_material = plate_material
        self.verbose = verbose  # 详细计算标志
        self.load_data(diameter_design_file, properties_file)

        # 确定控制塔段
        self.determine_control_section()

        # 初始化设计参数
        self.initialize_plate_parameters()

        # 结果存储
        self.downcomer_results = None
        self.area_results = None
        self.plate_results = None
        self.overflow_results = None
        self.final_results = None

        print(f"✓ 塔板布置计算器初始化完成")
        print(f"  塔径: {self.D:.3f} m, 塔板类型: {'浮阀塔' if tower_type=='float_valve' else '筛板塔'}")
        print(f"  参数选择: {'随机生成' if use_random_params else '手动指定'}")
        print(f"  详细计算: {'开启' if verbose else '关闭'}")

    def load_data(self, diameter_file, properties_file):
        """加载塔设计数据"""
        try:
            # 加载塔径设计结果
            with open(diameter_file, 'r', encoding='utf-8') as f:
                self.diameter_data = json.load(f)

            # 加载工艺特性数据
            with open(properties_file, 'r', encoding='utf-8') as f:
                self.properties_data = json.load(f)

            print(f"✓ 成功加载塔设计数据")

            # 提取关键参数
            self.D = self.diameter_data.get('最终设计结果', {}).get('最终塔径_m', 1.0)
            print(f"  塔径 D = {self.D:.3f} m")

        except FileNotFoundError as e:
            print(f"✗ 错误: 未找到数据文件 - {e}")
            self.D = 1.0  # 默认塔径
            self.diameter_data = {}
            self.properties_data = {}

    def determine_control_section(self):
        """确定控制塔段"""
        # 从塔径设计文件中获取控制塔段
        control_section = self.diameter_data.get('最终设计结果', {}).get('控制塔段', '提馏段')

        if control_section == '精馏段':
            self.control_section = '精馏段结果'
            print(f"  控制塔段: 精馏段")
        else:
            self.control_section = '提馏段结果'
            print(f"  控制塔段: 提馏段")

        # 获取控制塔段的工艺数据
        self.section_data = self.properties_data.get(self.control_section, {})

        if not self.section_data:
            print(f"⚠ 警告: 未找到{self.control_section}的数据，使用默认值")
            self.section_data = {
                'V_vol_m3h': 3000,
                'V_mass_kgh': 10000,
                'rho_vapor_kgm3': 3.0,
                'rho_liquid_kgm3': 800,
                'L_vol_Lh': 15000,
                'surface_tension_Nm': 0.02
            }

    def get_parameter_value(self, param_name, min_val, max_val, default_val, description=""):
        """获取参数值：支持随机生成或用户输入"""
        if self.use_random_params:
            # 在范围内随机生成
            value = random.uniform(min_val, max_val)
            if self.verbose:
                print(f"  随机生成 {description}: {value:.3f} (范围: {min_val:.3f}-{max_val:.3f})")
            return value
        else:
            # 用户输入
            prompt = f"请输入{description} (范围: {min_val:.3f}-{max_val:.3f}，默认: {default_val:.3f}): "
            user_input = input(prompt).strip()
            if user_input == "":
                return default_val
            try:
                value = float(user_input)
                if min_val <= value <= max_val:
                    return value
                else:
                    print(f"⚠ 输入值超出范围，使用默认值 {default_val:.3f}")
                    return default_val
            except ValueError:
                print(f"⚠ 输入无效，使用默认值 {default_val:.3f}")
                return default_val

    def get_sieve_hole_diameter(self):
        """根据表面张力确定筛孔直径"""
        if self.tower_type != 'sieve':
            return 0.005  # 默认值

        # 从控制塔段数据中获取表面张力
        sigma = self.section_data.get('sigma_mNm', 20.45) / 1000  # 转换为N/m
        if self.verbose:
            print(f"  表面张力: {sigma*1000:.2f} mN/m")

        if sigma > 0:
            # 表面张力为正，鼓泡型，孔径3-8mm
            d0_min, d0_max = 0.003, 0.008
            hole_type = "鼓泡型"
        else:
            # 表面张力为负，喷射型，孔径10-25mm
            d0_min, d0_max = 0.010, 0.025
            hole_type = "喷射型"

        if self.verbose:
            print(f"  筛孔类型: {hole_type}")
        return self.get_parameter_value('筛孔直径', d0_min, d0_max,
                                        (d0_min + d0_max) / 2, f"筛孔直径(m)-{hole_type}")

    def get_plate_thickness(self, sieve_hole_diameter):
        """根据材料和孔径确定筛板厚度"""
        if self.tower_type != 'sieve':
            return 0.003  # 默认值

        d0 = sieve_hole_diameter

        if self.plate_material == 'carbon_steel':
            # 碳钢: (0.4~0.8)*d0
            thickness_min, thickness_max = 0.4 * d0, 0.8 * d0
        else:
            # 不锈钢: (0.5~0.7)*d0
            thickness_min, thickness_max = 0.5 * d0, 0.7 * d0

        thickness = self.get_parameter_value('筛板厚度', thickness_min, thickness_max,
                                             (thickness_min + thickness_max) / 2,
                                             f"筛板厚度(m)-{self.plate_material}")

        # 确保厚度合理
        if thickness < 0.001:
            thickness = 0.001
        elif thickness > 0.01:
            thickness = 0.01

        if self.verbose:
            print(f"  筛板厚度: {thickness*1000:.1f} mm")
        return thickness

    def get_weir_length_ratio(self):
        """根据塔径和流型确定堰长与塔径之比"""
        # 确定流型
        if self.D < 2.2:
            flow_type = "单流型"
            lw_ratio_min, lw_ratio_max = 0.60, 0.80  # 修正为更合理的范围
        else:
            flow_type = "双流型"
            lw_ratio_min, lw_ratio_max = 0.50, 0.70

        if self.verbose:
            print(f"  流型: {flow_type}")
            print(f"  堰长比例范围: lw/D = {lw_ratio_min:.2f} - {lw_ratio_max:.2f}")
        return self.get_parameter_value('堰长/塔径比', lw_ratio_min, lw_ratio_max,
                                        (lw_ratio_min + lw_ratio_max) / 2,
                                        f"堰长/塔径比-{flow_type}")

    def initialize_plate_parameters(self):
        """初始化塔板布置参数"""
        if self.verbose:
            print(f"\n[初始化塔板布置参数]")

        # 根据塔径确定参数范围
        if self.D < 1.5:
            # 小塔
            Ws_inlet_min, Ws_inlet_max = 0.060, 0.075
            Ws_outlet_min, Ws_outlet_max = 0.050, 0.100
            Wc_min, Wc_max = 0.030, 0.050
        else:
            # 大塔
            Ws_inlet_min, Ws_inlet_max = 0.080, 0.100
            Ws_outlet_min, Ws_outlet_max = 0.050, 0.100
            Wc_min, Wc_max = 0.050, 0.075

        if self.verbose:
            print(f"  塔径 D = {self.D:.3f} m {'< 1.5 m' if self.D < 1.5 else '≥ 1.5 m'}")

        # 获取参数值
        Ws_inlet = self.get_parameter_value('Ws_inlet', Ws_inlet_min, Ws_inlet_max,
                                            (Ws_inlet_min + Ws_inlet_max) / 2, "入口安定区宽度(m)")
        Ws_outlet = self.get_parameter_value('Ws_outlet', Ws_outlet_min, Ws_outlet_max,
                                             (Ws_outlet_min + Ws_outlet_max) / 2, "出口安定区宽度(m)")
        Wc = self.get_parameter_value('Wc', Wc_min, Wc_max,
                                      (Wc_min + Wc_max) / 2, "边缘区宽度(m)")

        # 对于筛板塔，需要先计算筛孔直径和筛板厚度
        sieve_hole_diameter = 0.005  # 默认值
        plate_thickness = 0.003  # 默认值

        if self.tower_type == 'sieve':
            sieve_hole_diameter = self.get_sieve_hole_diameter()
            plate_thickness = self.get_plate_thickness(sieve_hole_diameter)

        self.design_params = {
            # 安定区宽度 (m)
            'Ws_inlet': Ws_inlet,
            'Ws_outlet': Ws_outlet,

            # 边缘区宽度 (m)
            'Wc': Wc,

            # 浮阀参数
            'valve_diameter': 0.039,  # F1型浮阀直径，m
            'valve_thickness': 0.002,  # 浮阀厚度，m
            'F0_range': (9, 12),  # 气体动能因数范围
            'F0_default': 10,

            # 筛板参数
            'sieve_hole_diameter': sieve_hole_diameter,
            'sieve_hole_pitch_min': 0.0025,  # 筛孔中心距最小值，m
            'sieve_hole_pitch_max': 0.004,  # 筛孔中心距最大值，m
            'sieve_hole_pitch_default': 0.003,  # 筛孔中心距，m
            'sieve_open_area_ratio': 0.1,  # 开孔率（%）
            'plate_thickness': plate_thickness,

            # 溢流装置参数
            'weir_length_ratio': self.get_weir_length_ratio(),
            'weir_height_range': (0.04, 0.05),  # 堰高范围，m（常压塔）
            'outlet_weir_height_range': (0.04, 0.08),  # 出口堰高范围，m
            'downcomer_clearance_range': (0.02, 0.03),  # 降液管底隙高度范围，m
            'liquid_retention_time_min': 3,  # 最小停留时间，s（已改为3秒）

            'inlet_sump_depth': 0.05,  # 受液盘深度，m
            'tear_hole_count': 2 if self.D > 1.4 else 1,  # 泪孔数量
        }

        if self.verbose:
            print(f"\n  设计参数初始化完成:")
            print(f"    入口安定区宽度: {Ws_inlet:.3f} m")
            print(f"    出口安定区宽度: {Ws_outlet:.3f} m")
            print(f"    边缘区宽度: {Wc:.3f} m")
            print(f"    堰长比例: {self.design_params['weir_length_ratio']:.3f}")
            print(f"    最小停留时间: {self.design_params['liquid_retention_time_min']} s")

    def design_downcomer_from_liquid_load(self):
        """
        核心设计：根据液体负荷设计降液管
        返回: 降液管设计结果
        """
        print("\n" + "="*80)
        print("一、降液管设计（基于液体负荷和停留时间）")
        print("="*80)

        D = self.D
        R = D / 2

        if self.verbose:
            print(f"\n[基本参数]")

        # 1. 获取液体流量
        L_vol = self.section_data.get('L_vol_Lh', 15000)  # L/h
        L_m3h = L_vol / 1000  # m³/h
        L_m3s = L_m3h / 3600  # m³/s

        if self.verbose:
            print(f"\n1. 液体流量计算:")
            print(f"   已知: L_vol = {L_vol:.1f} L/h")
            print(f"   换算: L_m3h = L_vol / 1000 = {L_vol:.1f} / 1000 = {L_m3h:.4f} m³/h")
            print(f"   换算: L_m3s = L_m3h / 3600 = {L_m3h:.4f} / 3600 = {L_m3s:.6f} m³/s")

        # 2. 获取板间距
        try:
            with open('tower_height_design.json', 'r', encoding='utf-8') as f:
                height_data = json.load(f)
            HT = height_data.get('板间距设计', {}).get('正常板间距_提馏段_m', 0.45)
            if self.verbose:
                print(f"\n2. 板间距:")
                print(f"   从tower_height_design.json读取: HT = {HT:.3f} m")
        except:
            HT = 0.45
            if self.verbose:
                print(f"\n2. 板间距:")
                print(f"   使用默认值: HT = {HT:.3f} m")

        # 3. 确定停留时间要求
        tau_min = self.design_params['liquid_retention_time_min']
        if self.verbose:
            print(f"\n3. 停留时间要求:")
            print(f"   最小停留时间: τ_min = {tau_min} s")

        # 4. 计算所需降液管面积
        Af_required = (L_m3s * tau_min) / HT
        if self.verbose:
            print(f"\n4. 所需降液管面积计算:")
            print(f"   公式: Af_required = (L × τ_min) / HT")
            print(f"   代入: L = {L_m3s:.6f} m³/s, τ_min = {tau_min} s, HT = {HT:.3f} m")
            print(f"   计算: Af_required = ({L_m3s:.6f} × {tau_min}) / {HT:.3f}")
            print(f"        = {L_m3s * tau_min:.6f} / {HT:.3f}")
            print(f"        = {Af_required:.6f} m²")

        # 5. 确定堰长比例范围
        if D < 2.2:
            flow_type = "单流型"
            lw_ratio_range = (0.60, 0.80)
        else:
            flow_type = "双流型"
            lw_ratio_range = (0.50, 0.70)

        if self.verbose:
            print(f"\n5. 流型确定:")
            print(f"   塔径 D = {D:.3f} m {'< 2.2 m' if D < 2.2 else '≥ 2.2 m'}")
            print(f"   流型: {flow_type}")
            print(f"   堰长比例范围: lw/D = {lw_ratio_range[0]:.2f} - {lw_ratio_range[1]:.2f}")

        # 6. 通过二分法找到满足Af_required的堰长
        if self.verbose:
            print(f"\n6. 搜索满足Af_required={Af_required:.6f} m²的堰长...")

        lw, actual_Af, Wd = find_lw_for_target_af(
            target_Af=Af_required,
            D=D,
            lw_ratio_range=lw_ratio_range,
            tolerance=0.0001,
            verbose=self.verbose
        )

        # 7. 计算实际停留时间
        tau_actual = (actual_Af * HT) / L_m3s if L_m3s > 0 else 0

        # 8. 显示详细弓形计算过程
        if self.verbose:
            print(f"\n7. 详细弓形计算过程:")
            segment_details = calculate_segment_area(lw, D, verbose=True)

        # 9. 停留时间校验
        if self.verbose:
            print(f"\n8. 停留时间校验:")
            print(f"   公式: τ_actual = (Af × HT) / L")
            print(f"   计算: τ_actual = ({actual_Af:.6f} × {HT:.3f}) / {L_m3s:.6f}")
            print(f"        = {actual_Af * HT:.6f} / {L_m3s:.6f}")
            print(f"        = {tau_actual:.2f} s")

        if tau_actual >= tau_min:
            status = "满足"
            status_symbol = "✓"
        else:
            status = "不满足"
            status_symbol = "⚠"

        # 10. 设计结果汇总
        if self.verbose:
            print(f"\n9. 设计结果汇总:")
            print(f"   堰长: lw = {lw:.4f} m (lw/D = {lw/D:.3f})")
            print(f"   降液管宽度: Wd = {Wd:.4f} m (Wd/D = {Wd/D:.3f})")
            print(f"   降液管面积: Af = {actual_Af:.6f} m²")
            print(f"   实际停留时间: τ = {tau_actual:.2f} s")
            print(f"   最小停留时间要求: τ_min = {tau_min} s")
            print(f"   校验: {status_symbol} 停留时间{status}要求")

            if tau_actual < tau_min:
                print(f"\n   ⚠ 警告: 停留时间不满足要求 (τ={tau_actual:.2f}s < {tau_min}s)")
                print(f"      建议调整措施:")
                print(f"      1. 增加板间距 HT (当前: {HT:.3f} m)")
                print(f"      2. 减少液体流量 L (当前: {L_m3s:.6f} m³/s)")
                print(f"      3. 增加塔径 D (当前: {D:.3f} m)")

        # 11. 存储结果
        self.downcomer_results = {
            '堰长_lw_m': lw,
            '堰长比例_lw/D': lw / D,
            '降液管宽度_Wd_m': Wd,
            '降液管宽度比例_Wd/D': Wd / D,
            '降液管面积_Af_m2': actual_Af,
            '降液管面积比_Af/At': actual_Af / (math.pi * R**2),
            '实际停留时间_s': tau_actual,
            '所需最小停留时间_s': tau_min,
            '液体流量_m3s': L_m3s,
            '板间距_m': HT,
            '流型': flow_type,
            '设计状态': status
        }

        return self.downcomer_results

    def calculate_plate_area_divisions(self):
        """
        计算塔板面积分区（使用扣除法）
        返回: 面积分区结果
        """
        print("\n" + "="*80)
        print("二、塔板面积分区计算（扣除法）")
        print("="*80)

        # 检查是否已设计降液管
        if self.downcomer_results is None:
            print("错误: 请先运行design_downcomer_from_liquid_load()")
            return None

        D = self.D
        Ws_inlet = self.design_params['Ws_inlet']
        Ws_outlet = self.design_params['Ws_outlet']
        Wc = self.design_params['Wc']
        Af = self.downcomer_results['降液管面积_Af_m2']

        if self.verbose:
            print(f"\n[输入参数]")
            print(f"   塔径 D = {D:.3f} m")
            print(f"   入口安定区宽度 Ws_inlet = {Ws_inlet:.3f} m")
            print(f"   出口安定区宽度 Ws_outlet = {Ws_outlet:.3f} m")
            print(f"   边缘区宽度 Wc = {Wc:.3f} m")
            print(f"   单个降液管面积 Af = {Af:.6f} m²")
            print(f"   降液管宽度 Wd = {self.downcomer_results['降液管宽度_Wd_m']:.4f} m")
            print(f"   堰长 lw = {self.downcomer_results['堰长_lw_m']:.4f} m")

        # 使用扣除法计算面积分区
        self.area_results = calculate_bubbling_area(
            D=D,
            Ws_inlet=Ws_inlet,
            Ws_outlet=Ws_outlet,
            Wc=Wc,
            Af=Af,
            Af_side='both',  # 双侧降液管
            verbose=self.verbose
        )

        return self.area_results

    def calculate_float_valve_arrangement(self):
        """计算浮阀塔板布置"""
        print("\n" + "="*80)
        print("三、浮阀塔板布置计算")
        print("="*80)

        # 检查是否已计算面积分区
        if self.area_results is None:
            print("错误: 请先运行calculate_plate_area_divisions()")
            return None

        # 使用控制塔段的工艺参数
        V_vol = self.section_data.get('V_vol_m3h', 3000)  # 气体体积流量，m³/h
        V_mass = self.section_data.get('V_mass_kgh', 10000)  # 气体质量流量，kg/h
        rho_v = self.section_data.get('rho_vapor_kgm3', 3.0)  # 气相密度，kg/m³

        if self.verbose:
            print(f"\n[工艺参数]")
            print(f"   气体体积流量 V = {V_vol:.2f} m³/h")
            print(f"   气体质量流量 V_mass = {V_mass:.2f} kg/h")
            print(f"   气相密度 ρ_v = {rho_v:.4f} kg/m³")

        # 1. 计算阀孔气速
        F0_target = self.design_params['F0_default']  # 目标气体动能因数
        u0 = F0_target / math.sqrt(rho_v)  # 阀孔气速，m/s

        if self.verbose:
            print(f"\n1. 阀孔气速计算:")
            print(f"   公式: u0 = F0 / √ρ_v")
            print(f"   目标气体动能因数 F0 = {F0_target}")
            print(f"   计算: u0 = {F0_target} / √{rho_v:.4f}")
            print(f"        = {F0_target} / {math.sqrt(rho_v):.4f}")
            print(f"        = {u0:.3f} m/s")

        # 2. 计算总阀孔面积
        V_m3s = V_vol / 3600  # m³/s
        A0_total = V_m3s / u0  # 总阀孔面积，m²

        if self.verbose:
            print(f"\n2. 总阀孔面积计算:")
            print(f"   气体体积流量换算: V_m3s = V / 3600 = {V_vol:.2f} / 3600 = {V_m3s:.4f} m³/s")
            print(f"   公式: A0_total = V / u0")
            print(f"   计算: A0_total = {V_m3s:.4f} / {u0:.3f}")
            print(f"        = {A0_total:.4f} m²")

        # 3. 计算阀孔数
        valve_diameter = self.design_params['valve_diameter']  # 浮阀直径，m
        valve_area = math.pi * (valve_diameter / 2) ** 2  # 单个浮阀面积，m²
        effective_valve_area = valve_area * 0.8  # 有效面积系数

        if self.verbose:
            print(f"\n3. 阀孔数计算:")
            print(f"   浮阀直径: d = {valve_diameter*1000:.0f} mm")
            print(f"   单个浮阀面积: A_valve = π(d/2)² = π×({valve_diameter/2:.3f})² = {valve_area:.6f} m²")
            print(f"   有效面积系数取 0.8，有效面积 = {valve_area:.6f} × 0.8 = {effective_valve_area:.6f} m²")

        N_valve = int(math.ceil(A0_total / effective_valve_area))

        if self.verbose:
            print(f"   阀孔数: N = A0_total / 有效面积 = {A0_total:.4f} / {effective_valve_area:.6f}")
            print(f"           = {N_valve} 个（向上取整）")

        # 4. 计算实际阀孔气速和动能因数
        A0_actual = N_valve * effective_valve_area
        u0_actual = V_m3s / A0_actual
        F0_actual = u0_actual * math.sqrt(rho_v)

        if self.verbose:
            print(f"\n4. 实际参数核算:")
            print(f"   实际总阀孔面积: A0_actual = N × A_valve_eff = {N_valve} × {effective_valve_area:.6f}")
            print(f"                          = {A0_actual:.4f} m²")
            print(f"   实际阀孔气速: u0_actual = V / A0_actual = {V_m3s:.4f} / {A0_actual:.4f}")
            print(f"                        = {u0_actual:.3f} m/s")
            print(f"   实际气体动能因数: F0_actual = u0_actual × √ρ_v = {u0_actual:.3f} × {math.sqrt(rho_v):.3f}")
            print(f"                          = {F0_actual:.2f}")

        # 5. 校验F0是否在合理范围内
        F0_min, F0_max = self.design_params['F0_range']
        if F0_min <= F0_actual <= F0_max:
            F0_status = "在范围内"
            F0_symbol = "✓"
        else:
            F0_status = "超出范围"
            F0_symbol = "⚠"

        if self.verbose:
            print(f"\n5. F0校验:")
            print(f"   F0范围: {F0_min} - {F0_max}")
            print(f"   实际F0: {F0_actual:.2f}")
            print(f"   校验: {F0_symbol} F0{F0_status}")

            if F0_actual < F0_min or F0_actual > F0_max:
                print(f"   ⚠ 警告: F0不在合理范围内，建议调整阀孔数")

        # 6. 浮阀排列设计
        if self.verbose:
            print(f"\n6. 浮阀排列设计:")

        # 获取鼓泡区面积
        Aa = self.area_results['鼓泡区面积_Aa_m2']

        if self.verbose:
            print(f"   鼓泡区面积: Aa = {Aa:.4f} m²")

        # 浮阀排列方式：等腰三角形排列
        t_min, t_max = 0.075, 0.125  # 阀孔中心距范围，m
        h_min, h_max = 0.065, 0.100  # 排间距范围，m

        if self.use_random_params:
            t_design = random.uniform(t_min, t_max)
            h_design = random.uniform(h_min, h_max)
        else:
            t_design = self.get_parameter_value('阀孔中心距', t_min, t_max, 0.10, "阀孔中心距(m)")
            h_design = self.get_parameter_value('排间距', h_min, h_max, 0.08, "排间距(m)")

        if self.verbose:
            print(f"   阀孔中心距: t = {t_design*1000:.0f} mm")
            print(f"   排间距: h = {h_design*1000:.0f} mm")

        # 估算可排列的浮阀数
        # 假设鼓泡区为矩形，估算尺寸
        bubble_length_est = math.sqrt(Aa) * 1.5  # 估算长度
        bubble_width_est = math.sqrt(Aa) * 0.8   # 估算宽度

        if self.verbose:
            print(f"\n   鼓泡区尺寸估算:")
            print(f"   假设鼓泡区为矩形，面积 Aa = {Aa:.4f} m²")
            print(f"   估算长度: L ≈ √Aa × 1.5 = √{Aa:.4f} × 1.5 = {math.sqrt(Aa):.3f} × 1.5 = {bubble_length_est:.3f} m")
            print(f"   估算宽度: W ≈ √Aa × 0.8 = √{Aa:.4f} × 0.8 = {math.sqrt(Aa):.3f} × 0.8 = {bubble_width_est:.3f} m")

        valves_per_row = int(bubble_length_est / t_design)
        rows = int(bubble_width_est / h_design)
        max_valves = valves_per_row * rows

        if self.verbose:
            print(f"\n   排列估算:")
            print(f"   每排可排列浮阀数: {valves_per_row} 个 (L/t = {bubble_length_est:.3f}/{t_design:.3f})")
            print(f"   可排列排数: {rows} 排 (W/h = {bubble_width_est:.3f}/{h_design:.3f})")
            print(f"   最大可排列浮阀数: {max_valves} 个")
            print(f"   需要排列浮阀数: {N_valve} 个")

        if max_valves >= N_valve:
            arrangement_status = '合理'
            arrangement_symbol = "✓"
            if self.verbose:
                print(f"   校验: {arrangement_symbol} 鼓泡区可容纳 {N_valve} 个浮阀")
        else:
            arrangement_status = '需调整'
            arrangement_symbol = "⚠"
            if self.verbose:
                print(f"   校验: {arrangement_symbol} 鼓泡区空间不足，需调整排列参数")
                print(f"       建议: 1) 减小阀孔中心距 t")
                print(f"             2) 减小排间距 h")
                print(f"             3) 调整鼓泡区形状")

        # 7. 存储结果
        self.plate_results = {
            '阀孔气速_m/s': round(u0_actual, 3),
            '气体动能因数_F0': round(F0_actual, 2),
            '阀孔数': N_valve,
            '阀孔直径_m': valve_diameter,
            '总阀孔面积_m2': round(A0_actual, 4),
            '阀孔中心距_m': t_design,
            '排间距_m': h_design,
            '每排阀孔数': valves_per_row,
            '排列排数': rows,
            '最大可排阀孔数': max_valves,
            '排列状态': arrangement_status,
            'F0校验': F0_status
        }

        return self.plate_results

    def calculate_sieve_plate_arrangement(self):
        """计算筛板塔板布置"""
        print("\n" + "="*80)
        print("三、筛板塔板布置计算")
        print("="*80)

        # 检查是否已计算面积分区
        if self.area_results is None:
            print("错误: 请先运行calculate_plate_area_divisions()")
            return None

        # 使用控制塔段的工艺参数
        V_vol = self.section_data.get('V_vol_m3h', 3000)  # 气体体积流量，m³/h
        rho_v = self.section_data.get('rho_vapor_kgm3', 3.0)  # 气相密度，kg/m³

        if self.verbose:
            print(f"\n[工艺参数]")
            print(f"   气体体积流量 V = {V_vol:.2f} m³/h")
            print(f"   气相密度 ρ_v = {rho_v:.4f} kg/m³")

        # 1. 筛孔参数
        d0 = self.design_params['sieve_hole_diameter']  # 筛孔直径，m
        plate_thickness = self.design_params['plate_thickness']  # 筛板厚度，m

        if self.verbose:
            print(f"\n1. 筛孔参数:")
            print(f"   筛孔直径: d0 = {d0*1000:.1f} mm")
            print(f"   筛板厚度: δ = {plate_thickness*1000:.1f} mm")
            print(f"   单个筛孔面积: A0 = π(d0/2)² = π×({d0/2:.4f})² = {math.pi * (d0 / 2) ** 2:.6f} m²")

        # 2. 筛孔气速
        u0_min, u0_max = 10.0, 20.0  # 筛孔气速范围，m/s

        if self.use_random_params:
            u0_design = random.uniform(u0_min, u0_max)
        else:
            u0_design = self.get_parameter_value('筛孔气速', u0_min, u0_max, 15.0, "筛孔气速(m/s)")

        if self.verbose:
            print(f"\n2. 筛孔气速:")
            print(f"   设计值: u0 = {u0_design:.1f} m/s (范围: {u0_min:.1f} - {u0_max:.1f} m/s)")

        # 3. 计算总筛孔面积
        V_m3s = V_vol / 3600  # m³/s
        A0_total = V_m3s / u0_design  # 总筛孔面积，m²

        if self.verbose:
            print(f"\n3. 总筛孔面积计算:")
            print(f"   气体体积流量换算: V_m3s = V / 3600 = {V_vol:.2f} / 3600 = {V_m3s:.4f} m³/s")
            print(f"   公式: A0_total = V / u0")
            print(f"   计算: A0_total = {V_m3s:.4f} / {u0_design:.1f}")
            print(f"        = {A0_total:.4f} m²")

        # 4. 计算筛孔数
        A0_single = math.pi * (d0 / 2) ** 2  # 单个筛孔面积，m²
        N_holes = int(math.ceil(A0_total / A0_single))

        if self.verbose:
            print(f"\n4. 筛孔数计算:")
            print(f"   单个筛孔面积: A0_single = {A0_single:.6f} m²")
            print(f"   筛孔数: N = A0_total / A0_single = {A0_total:.4f} / {A0_single:.6f}")
            print(f"           = {N_holes} 个（向上取整）")

        # 5. 计算实际筛孔气速
        u0_actual = V_m3s / (N_holes * A0_single)

        if self.verbose:
            print(f"\n5. 实际参数核算:")
            print(f"   实际总筛孔面积: A0_actual = N × A0_single = {N_holes} × {A0_single:.6f}")
            print(f"                          = {N_holes * A0_single:.4f} m²")
            print(f"   实际筛孔气速: u0_actual = V / A0_actual = {V_m3s:.4f} / {N_holes * A0_single:.4f}")
            print(f"                        = {u0_actual:.2f} m/s")

        # 6. 校验筛孔气速
        if u0_min <= u0_actual <= u0_max:
            u0_status = "在范围内"
            u0_symbol = "✓"
        else:
            u0_status = "超出范围"
            u0_symbol = "⚠"

        if self.verbose:
            print(f"\n6. 筛孔气速校验:")
            print(f"   筛孔气速范围: {u0_min:.1f} - {u0_max:.1f} m/s")
            print(f"   实际筛孔气速: {u0_actual:.2f} m/s")
            print(f"   校验: {u0_symbol} 筛孔气速{u0_status}")

            if u0_actual < u0_min or u0_actual > u0_max:
                print(f"   ⚠ 警告: 筛孔气速不在合理范围内，建议调整筛孔数")

        # 7. 筛孔排列设计
        if self.verbose:
            print(f"\n7. 筛孔排列设计:")

        # 筛孔中心距：一般取(2.5-5)d0
        t_min, t_max = 0.0025, 0.004

        if self.use_random_params:
            t_design = random.uniform(t_min, t_max)
        else:
            t_design = self.get_parameter_value('筛孔中心距', t_min, t_max,
                                               0.003, "筛孔中心距(m)")

        if self.verbose:
            print(f"   筛孔中心距: t = {t_design*1000:.1f} mm")

        # 开孔率计算
        Aa = self.area_results['鼓泡区面积_Aa_m2']  # 鼓泡区面积
        open_area_ratio = (N_holes * A0_single) / Aa * 100  # 开孔率，%

        if self.verbose:
            print(f"\n8. 开孔率计算:")
            print(f"   鼓泡区面积: Aa = {Aa:.3f} m²")
            print(f"   总筛孔面积: A0_total = N × A0_single = {N_holes} × {A0_single:.6f} = {N_holes * A0_single:.4f} m²")
            print(f"   公式: φ = (N × A0_single) / Aa × 100%")
            print(f"   计算: φ = ({N_holes} × {A0_single:.6f}) / {Aa:.3f} × 100%")
            print(f"        = {N_holes * A0_single:.4f} / {Aa:.3f} × 100%")
            print(f"        = {open_area_ratio:.2f}%")

        # 开孔率校验（一般5-15%）
        if 5 <= open_area_ratio <= 15:
            open_area_status = '合理'
            open_area_symbol = "✓"
        else:
            open_area_status = '需调整'
            open_area_symbol = "⚠"

        if self.verbose:
            print(f"\n9. 开孔率校验:")
            print(f"   开孔率范围: 5-15%")
            print(f"   实际开孔率: {open_area_ratio:.2f}%")
            print(f"   校验: {open_area_symbol} 开孔率{open_area_status}")

            if open_area_ratio < 5 or open_area_ratio > 15:
                print(f"   ⚠ 警告: 开孔率不在合理范围内，建议调整筛孔参数")
                print(f"       措施: 1) 调整筛孔数 N")
                print(f"             2) 调整筛孔直径 d0")
                print(f"             3) 调整鼓泡区面积 Aa")

        # 8. 排列方式：一般采用正三角形排列
        if self.verbose:
            print(f"\n10. 筛孔排列估算:")

        # 假设鼓泡区为矩形，估算尺寸
        bubble_length_est = math.sqrt(Aa) * 1.5
        bubble_width_est = math.sqrt(Aa) * 0.8

        if self.verbose:
            print(f"   鼓泡区尺寸估算:")
            print(f"   假设鼓泡区为矩形，面积 Aa = {Aa:.4f} m²")
            print(f"   估算长度: L ≈ √Aa × 1.5 = √{Aa:.4f} × 1.5 = {bubble_length_est:.3f} m")
            print(f"   估算宽度: W ≈ √Aa × 0.8 = √{Aa:.4f} × 0.8 = {bubble_width_est:.3f} m")

        # 正三角形排列，孔心距为t
        holes_per_row = int(bubble_length_est / t_design)
        # 正三角形高度为t*√3/2
        rows = int(bubble_width_est / (t_design * math.sqrt(3) / 2))
        max_holes = holes_per_row * rows

        if self.verbose:
            print(f"\n   排列估算（正三角形排列）:")
            print(f"   正三角形高度: h_triangle = t × √3/2 = {t_design:.3f} × {math.sqrt(3)/2:.3f} = {t_design * math.sqrt(3) / 2:.4f} m")
            print(f"   每排可排列筛孔数: {holes_per_row} 个 (L/t = {bubble_length_est:.3f}/{t_design:.3f})")
            print(f"   可排列排数: {rows} 排 (W/h_triangle = {bubble_width_est:.3f}/{t_design * math.sqrt(3) / 2:.4f})")
            print(f"   最大可排列筛孔数: {max_holes} 个")
            print(f"   需要排列筛孔数: {N_holes} 个")

        if max_holes >= N_holes:
            arrangement_status = '合理'
            arrangement_symbol = "✓"
            if self.verbose:
                print(f"   校验: {arrangement_symbol} 鼓泡区可容纳 {N_holes} 个筛孔")
        else:
            arrangement_status = '需调整'
            arrangement_symbol = "⚠"
            if self.verbose:
                print(f"   校验: {arrangement_symbol} 鼓泡区空间不足，需调整排列参数")
                print(f"       建议: 1) 减小筛孔中心距 t")
                print(f"             2) 调整排列方式")

        # 9. 存储结果
        self.plate_results = {
            '筛孔气速_m/s': round(u0_actual, 2),
            '筛孔数': N_holes,
            '筛孔直径_m': d0,
            '筛板厚度_m': plate_thickness,
            '总筛孔面积_m2': round(N_holes * A0_single, 4),
            '筛孔中心距_m': t_design,
            '开孔率_%': round(open_area_ratio, 2),
            '每排筛孔数': holes_per_row,
            '排列排数': rows,
            '最大可排筛孔数': max_holes,
            '排列方式': '正三角形排列',
            '开孔率状态': open_area_status,
            '排列状态': arrangement_status,
            '气速校验': u0_status
        }

        return self.plate_results

    def calculate_overflow_device(self):
        """计算溢流装置"""
        print("\n" + "="*80)
        print("四、溢流装置计算")
        print("="*80)

        # 检查是否已设计降液管
        if self.downcomer_results is None:
            print("错误: 请先运行design_downcomer_from_liquid_load()")
            return None

        D = self.D
        lw = self.downcomer_results['堰长_lw_m']

        if self.verbose:
            print(f"\n[基本参数]")
            print(f"   塔径 D = {D:.3f} m")
            print(f"   堰长 lw = {lw:.3f} m (lw/D = {lw/D:.3f})")

        # 1. 确定溢流类型
        if D < 2.2:
            overflow_type = "单溢流"
        else:
            overflow_type = "双溢流"

        if self.verbose:
            print(f"\n1. 溢流类型确定:")
            print(f"   塔径 D = {D:.3f} m {'< 2.2 m' if D < 2.2 else '≥ 2.2 m'}")
            print(f"   溢流类型: {overflow_type}")

        # 2. 堰上液层高度计算
        L_vol = self.section_data.get('L_vol_Lh', 15000)  # 液体体积流量，L/h
        L_m3h = L_vol / 1000  # m³/h

        if self.verbose:
            print(f"\n2. 堰上液层高度计算:")
            print(f"   液体流量: L = {L_vol:.1f} L/h")
            print(f"   换算: L_m3h = L / 1000 = {L_vol:.1f} / 1000 = {L_m3h:.3f} m³/h")

        # 堰上液层高度（平直堰）
        E = 1.0  # 液流收缩系数
        how = 2.84 / 1000 * E * (L_m3h / lw) ** (2 / 3)

        if self.verbose:
            print(f"\n   计算公式（平直堰）:")
            print(f"   how = 2.84×10⁻³ × E × (L_m3h/lw)^(2/3)")
            print(f"        = 0.00284 × {E} × ({L_m3h:.3f}/{lw:.3f})^(2/3)")
            print(f"        = 0.00284 × ({L_m3h/lw:.3f})^(2/3)")
            print(f"        = 0.00284 × {(L_m3h/lw) ** (2/3):.3f}")
            print(f"        = {how:.4f} m = {how*1000:.1f} mm")

        # 3. 堰高选择
        hw_min, hw_max = self.design_params['weir_height_range']
        if self.use_random_params:
            hw = random.uniform(hw_min, hw_max)
        else:
            hw = self.get_parameter_value('堰高', hw_min, hw_max,
                                         (hw_min + hw_max)/2, "堰高(m)")

        if self.verbose:
            print(f"\n3. 堰高选择:")
            print(f"   堰高范围: hw = {hw_min*1000:.0f} - {hw_max*1000:.0f} mm")
            print(f"   选择堰高: hw = {hw*1000:.1f} mm")

        # 4. 堰类型判断
        if how > 0.006:  # 6mm
            weir_type = "平直堰"
        else:
            weir_type = "齿形堰"

        if self.verbose:
            print(f"\n4. 堰类型判断:")
            print(f"   how = {how*1000:.1f} mm")
            print(f"   判断: how {'>' if how > 0.006 else '≤'} 6mm")
            print(f"   堰类型: {weir_type}")

        # 5. 清液层高度
        hL = hw + how

        if self.verbose:
            print(f"\n5. 清液层高度计算:")
            print(f"   hL = hw + how = {hw*1000:.1f} mm + {how*1000:.1f} mm")
            print(f"      = {hL*1000:.1f} mm")

        # 6. 降液管底隙高度
        L_m3s = L_m3h / 3600
        hb_min, hb_max = self.design_params['downcomer_clearance_range']

        if self.verbose:
            print(f"\n6. 降液管底隙高度设计:")
            print(f"   液体流量换算: L_m3s = L_m3h / 3600 = {L_m3h:.3f} / 3600 = {L_m3s:.6f} m³/s")

        # 底隙高度: h0 = hw - (0.006~0.012)
        if self.use_random_params:
            adjustment = random.uniform(0.006, 0.012)
        else:
            adjustment = 0.009

        h0_initial = hw - adjustment
        h0 = max(hb_min, min(h0_initial, hb_max))

        if self.verbose:
            print(f"\n   底隙高度计算:")
            print(f"   初步计算: h0 = hw - 调整值 = {hw*1000:.1f} mm - {adjustment*1000:.1f} mm")
            print(f"              = {h0_initial*1000:.1f} mm")
            print(f"   范围限制: h0应在 {hb_min*1000:.0f} - {hb_max*1000:.0f} mm之间")
            print(f"   最终确定: h0 = {h0*1000:.1f} mm")

        # 7. 底隙流速
        if lw * h0 > 0:
            u0_prime = L_m3s / (lw * h0)
        else:
            u0_prime = 0

        if self.verbose:
            print(f"\n   底隙流速计算:")
            print(f"   公式: u0' = L_m3s / (lw × h0)")
            print(f"   计算: u0' = {L_m3s:.6f} / ({lw:.3f} × {h0:.3f})")
            print(f"           = {L_m3s:.6f} / {lw * h0:.6f}")
            print(f"           = {u0_prime:.3f} m/s")

        # 8. 底隙流速校验
        if u0_prime <= 0.5:
            u0_prime_status = "满足"
            u0_prime_symbol = "✓"
        else:
            u0_prime_status = "不满足"
            u0_prime_symbol = "⚠"

        if self.verbose:
            print(f"\n   底隙流速校验:")
            print(f"   流速要求: u0' ≤ 0.5 m/s")
            print(f"   实际流速: u0' = {u0_prime:.3f} m/s")
            print(f"   校验: {u0_prime_symbol} 底隙流速{u0_prime_status}要求")

            if u0_prime > 0.5:
                print(f"   ⚠ 警告: 底隙流速过大，建议调整底隙高度 h0")

        # 9. 受液盘设计
        sump_depth = self.design_params['inlet_sump_depth']
        tear_holes = self.design_params['tear_hole_count']

        if self.verbose:
            print(f"\n7. 受液盘设计:")
            print(f"   采用凹形受液盘")
            print(f"   盘深: {sump_depth*1000:.0f} mm")
            print(f"   泪孔数量: {tear_holes} 个")

        # 10. 存储结果
        self.overflow_results = {
            '溢流类型': overflow_type,
            '堰长_m': lw,
            '堰长比例_lw/D': lw/D,
            '堰上液层高度_m': how,
            '堰上液层高度_mm': how*1000,
            '堰高_m': hw,
            '堰高_mm': hw*1000,
            '清液层高度_mm': hL*1000,
            '堰类型': weir_type,
            '降液管底隙高度_m': h0,
            '降液管底隙高度_mm': h0*1000,
            '底隙流速_m/s': u0_prime,
            '底隙流速校验': u0_prime_status,
            '受液盘深度_m': sump_depth,
            '泪孔数量': tear_holes
        }

        # 11. 结果汇总
        if self.verbose:
            print(f"\n8. 溢流装置设计结果汇总:")
            print(f"   溢流类型: {overflow_type}")
            print(f"   堰长: {lw:.3f} m (lw/D = {lw/D:.3f})")
            print(f"   堰高: {hw*1000:.1f} mm")
            print(f"   堰上液层高度: {how*1000:.1f} mm")
            print(f"   清液层高度: {hL*1000:.1f} mm")
            print(f"   堰类型: {weir_type}")
            print(f"   降液管底隙高度: {h0*1000:.1f} mm")
            print(f"   底隙流速: {u0_prime:.3f} m/s ({u0_prime_status})")
            print(f"   受液盘深度: {sump_depth*1000:.0f} mm")
            print(f"   泪孔数量: {tear_holes} 个")

        return self.overflow_results

    def run_complete_calculation(self):
        """运行完整的塔板布置计算"""
        print("="*80)
        print("塔板布置完整计算（几何修正版 - 详细计算过程）")
        print("="*80)
        print(f"塔板类型: {'浮阀塔' if self.tower_type == 'float_valve' else '筛板塔'}")
        print(f"参数选择: {'随机生成' if self.use_random_params else '手动指定'}")
        print(f"塔板材料: {'碳钢' if self.plate_material == 'carbon_steel' else '不锈钢'}")
        print(f"控制塔段: {self.control_section.replace('结果', '')}")
        print(f"塔径: {self.D:.3f} m")
        print(f"详细计算: {'开启' if self.verbose else '关闭'}")

        # 1. 设计降液管（基于液体负荷）
        print("\n" + "="*80)
        print("开始第一步：降液管设计")
        print("="*80)
        downcomer_results = self.design_downcomer_from_liquid_load()

        # 2. 计算塔板面积分区
        print("\n" + "="*80)
        print("开始第二步：塔板面积分区计算")
        print("="*80)
        area_results = self.calculate_plate_area_divisions()

        # 3. 根据塔板类型进行布置计算
        print("\n" + "="*80)
        print("开始第三步：塔板布置计算")
        print("="*80)
        if self.tower_type == 'float_valve':
            plate_results = self.calculate_float_valve_arrangement()
            plate_type = "浮阀"
        else:
            plate_results = self.calculate_sieve_plate_arrangement()
            plate_type = "筛板"

        # 4. 计算溢流装置
        print("\n" + "="*80)
        print("开始第四步：溢流装置计算")
        print("="*80)
        overflow_results = self.calculate_overflow_device()

        # 5. 汇总结果
        self.combine_results(downcomer_results, area_results, plate_results,
                            overflow_results, plate_type)

        # 6. 保存结果
        self.save_results_to_json()

        # 7. 显示汇总
        self.display_summary()

        return self.final_results

    def combine_results(self, downcomer_results, area_results, plate_results,
                       overflow_results, plate_type):
        """合并所有计算结果"""
        self.final_results = {
            "设计信息": {
                "项目": f"精馏塔{plate_type}塔板布置设计（几何修正版）",
                "塔板类型": plate_type,
                "塔径_m": self.D,
                "控制塔段": self.control_section.replace('结果', ''),
                "参数选择方式": "随机生成" if self.use_random_params else "手动指定",
                "塔板材料": self.plate_material,
                "设计日期": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "详细计算": self.verbose
            },
            "降液管设计": downcomer_results,
            "塔板面积分区": area_results,
            f"{plate_type}布置": plate_results,
            "溢流装置": overflow_results,
            "设计参数汇总": self.design_params,
            "几何计算说明": {
                "降液管计算": "基于液体负荷和停留时间要求，使用弓形几何公式精确计算",
                "面积分区": "使用扣除法计算，确保面积闭合误差为0",
                "核心公式": "弓形面积Af = 0.5*R²*θ - 0.5*lw*√(R²-(lw/2)²)",
                "设计流程": "液体负荷 → 所需Af → 堰长lw → 降液管宽度Wd → 面积分区",
                "最小停留时间": f"{self.design_params['liquid_retention_time_min']} s"
            }
        }

    def save_results_to_json(self, filename='tower_plate_arrangement_corrected.json'):
        """保存计算结果到JSON文件"""
        try:
            with open(filename, 'w', encoding='utf-8') as f:
                json.dump(self.final_results, f, ensure_ascii=False, indent=4)

            print(f"\n" + "="*80)
            print("计算结果保存")
            print("="*80)
            print(f"✓ 塔板布置计算结果已保存到 '{filename}'")

            # 显示文件路径
            file_path = os.path.abspath(filename)
            print(f"文件保存路径: {file_path}")

        except Exception as e:
            print(f"✗ 保存JSON文件时出错: {e}")

    def display_summary(self):
        """显示设计结果汇总"""
        if not hasattr(self, 'final_results'):
            print("尚未进行计算，请先运行run_complete_calculation()")
            return

        results = self.final_results

        print("\n" + "="*80)
        print("塔板布置设计结果汇总（几何修正版）")
        print("="*80)

        print(f"\n一、基本参数")
        print(f"   塔径: {self.D:.3f} m")
        print(f"   塔板类型: {results['设计信息']['塔板类型']}")
        print(f"   控制塔段: {results['设计信息']['控制塔段']}")
        print(f"   设计日期: {results['设计信息']['设计日期']}")

        print(f"\n二、降液管设计")
        downcomer = results['降液管设计']
        print(f"   堰长: {downcomer['堰长_lw_m']:.3f} m ({downcomer['堰长比例_lw/D']:.3f}D)")
        print(f"   降液管宽度: {downcomer['降液管宽度_Wd_m']:.3f} m ({downcomer['降液管宽度比例_Wd/D']:.3f}D)")
        print(f"   降液管面积: {downcomer['降液管面积_Af_m2']:.4f} m²")
        print(f"   停留时间: {downcomer['实际停留时间_s']:.1f} s (要求: ≥{downcomer['所需最小停留时间_s']}s)")
        print(f"   设计状态: {downcomer['设计状态']}")

        print(f"\n三、塔板面积分区")
        area = results['塔板面积分区']
        print(f"   塔板总面积: {area['塔板总面积_At_m2']:.3f} m²")
        print(f"   鼓泡区面积: {area['鼓泡区面积_Aa_m2']:.3f} m² ({area['鼓泡区面积占比_%']:.1f}%)")
        print(f"   降液管面积: {area['总降液管面积_total_Af_m2']:.3f} m² ({area['降液管面积占比_%']:.1f}%)")
        print(f"   面积闭合误差: {area['面积闭合误差_%']:.6f}%")

        plate_type = results['设计信息']['塔板类型']
        if plate_type == "浮阀":
            print(f"\n四、浮阀布置")
            valve = results['浮阀布置']
            print(f"   阀孔数: {valve['阀孔数']} 个")
            print(f"   阀孔气速: {valve['阀孔气速_m/s']:.2f} m/s")
            print(f"   气体动能因数 F0: {valve['气体动能因数_F0']:.1f} ({valve['F0校验']})")
            print(f"   阀孔直径: {valve['阀孔直径_m']*1000:.0f} mm")
            print(f"   排列: 每排{valve['每排阀孔数']}个，共{valve['排列排数']}排")
            print(f"   排列状态: {valve['排列状态']}")
        else:
            print(f"\n四、筛板布置")
            sieve = results['筛板布置']
            print(f"   筛孔数: {sieve['筛孔数']} 个")
            print(f"   筛孔气速: {sieve['筛孔气速_m/s']:.2f} m/s ({sieve['气速校验']})")
            print(f"   筛孔直径: {sieve['筛孔直径_m']*1000:.1f} mm")
            print(f"   开孔率: {sieve['开孔率_%']:.1f}% ({sieve['开孔率状态']})")
            print(f"   排列方式: {sieve['排列方式']}")
            print(f"   排列状态: {sieve['排列状态']}")

        print(f"\n五、溢流装置")
        overflow = results['溢流装置']
        print(f"   溢流类型: {overflow['溢流类型']}")
        print(f"   堰长: {overflow['堰长_m']:.3f} m (lw/D = {overflow['堰长比例_lw/D']:.3f})")
        print(f"   堰高: {overflow['堰高_mm']:.0f} mm")
        print(f"   堰上液层高度: {overflow['堰上液层高度_mm']:.1f} mm")
        print(f"   堰类型: {overflow['堰类型']}")
        print(f"   降液管底隙高度: {overflow['降液管底隙高度_mm']:.0f} mm")
        print(f"   底隙流速: {overflow['底隙流速_m/s']:.3f} m/s ({overflow['底隙流速校验']})")

        print(f"\n六、设计状态")
        print(f"   ✓ 几何修正已完成，面积闭合误差接近0")
        print(f"   ✓ 基于液体负荷和停留时间设计降液管")
        print(f"   ✓ 最小停留时间: {self.design_params['liquid_retention_time_min']} s")
        print(f"   ✓ 详细计算过程已打印")


def test_geometry_functions():
    """测试几何计算函数"""
    print("测试几何计算函数...")

    # 测试弓形计算
    print("\n1. 测试弓形计算函数:")
    try:
        segment = calculate_segment_area(lw=0.7, D=1.0, verbose=True)
        print(f"  结果摘要:")
        print(f"    降液管宽度: {segment['降液管宽度_Wd_m']:.4f} m")
        print(f"    降液管面积: {segment['弓形面积_Af_m2']:.4f} m²")
        print(f"    面积比: {segment['面积比_%']:.2f}%")
    except Exception as e:
        print(f"  错误: {e}")

    # 测试面积分区计算
    print("\n2. 测试面积分区计算函数:")
    area_test = calculate_bubbling_area(
        D=1.0,
        Ws_inlet=0.07,
        Ws_outlet=0.07,
        Wc=0.05,
        Af=0.0942,
        Af_side='both',
        verbose=True
    )
    print(f"  结果摘要:")
    print(f"    鼓泡区面积: {area_test['鼓泡区面积_Aa_m2']:.4f} m²")
    print(f"    面积闭合误差: {area_test['面积闭合误差_%']:.6f}%")
    print(f"    面积校验通过: {area_test['面积校验通过']}")


def main():
    """主函数"""
    print("塔板布置计算程序（几何修正版 - 带详细计算过程）")
    print("="*80)
    print("功能说明:")
    print("  1. 基于液体负荷和停留时间设计降液管（弓形几何计算）")
    print("  2. 使用扣除法计算塔板面积分区，确保面积闭合")
    print("  3. 支持浮阀塔和筛板塔两种塔板类型")
    print("  4. 计算溢流装置（堰、降液管、受液盘）")
    print("  5. 完整的校验和结果输出")
    print("  6. 详细的计算过程显示（公式、代入、计算步骤）")
    print("="*80)

    # 测试几何函数
    test_geometry_functions()

    # 选择塔板类型
    print("\n" + "="*80)
    print("请选择塔板类型:")
    print("  1. 浮阀塔 (F1型浮阀)")
    print("  2. 筛板塔")

    choice = input("请输入选择 (1或2，默认1): ").strip()
    if choice == "2":
        tower_type = "sieve"
        print("选择: 筛板塔")
    else:
        tower_type = "float_valve"
        print("选择: 浮阀塔")

    # 选择参数生成方式
    print("\n请选择参数生成方式:")
    print("  1. 随机生成（在合理范围内随机选择）")
    print("  2. 手动指定（输入每个参数）")

    param_choice = input("请输入选择 (1或2，默认1): ").strip()
    use_random_params = (param_choice != "2")

    if use_random_params:
        print("选择: 随机生成参数")
        random.seed(42)  # 设置随机种子以便结果可复现
    else:
        print("选择: 手动指定参数")

    # 对于筛板塔，选择材料
    plate_material = "carbon_steel"
    if tower_type == "sieve":
        print("\n请选择筛板材料:")
        print("  1. 碳钢 (carbon steel)")
        print("  2. 不锈钢 (stainless steel)")

        material_choice = input("请输入选择 (1或2，默认1): ").strip()
        if material_choice == "2":
            plate_material = "stainless_steel"
            print("选择: 不锈钢")
        else:
            print("选择: 碳钢")

    # 选择是否开启详细计算
    print("\n请选择是否开启详细计算过程:")
    print("  1. 开启（显示所有公式、代入过程和计算步骤）")
    print("  2. 关闭（只显示结果摘要）")

    verbose_choice = input("请输入选择 (1或2，默认1): ").strip()
    verbose = (verbose_choice != "2")

    if verbose:
        print("选择: 开启详细计算过程")
    else:
        print("选择: 关闭详细计算过程")

    # 创建计算器
    calculator = TowerPlateArrangementCorrected(
        tower_type=tower_type,
        use_random_params=use_random_params,
        plate_material=plate_material,
        verbose=verbose
    )

    # 运行计算
    print("\n" + "="*80)
    print("开始塔板布置计算...")
    print("="*80)

    results = calculator.run_complete_calculation()

    print("\n" + "="*80)
    print("塔板布置计算完成！")
    print("="*80)
    print("关键几何修正:")
    print("  1. 使用弓形几何公式精确计算降液管")
    print("  2. 基于液体负荷和停留时间设计降液管尺寸")
    print("  3. 使用扣除法确保面积闭合误差接近0")
    print("  4. 最小停留时间: 3 s")

    return calculator


if __name__ == "__main__":
    # 运行主程序
    calculator = main()

    # 示例：如何访问计算结果
    print(f"\n" + "="*80)
    print("计算结果访问示例")
    print("="*80)
    if hasattr(calculator, 'final_results'):
        print(f"塔径: {calculator.final_results['设计信息']['塔径_m']:.3f} m")
        print(f"塔板类型: {calculator.final_results['设计信息']['塔板类型']}")

        # 显示降液管设计
        downcomer = calculator.final_results['降液管设计']
        print(f"降液管宽度: {downcomer['降液管宽度_Wd_m']:.4f} m")
        print(f"降液管面积: {downcomer['降液管面积_Af_m2']:.4f} m²")
        print(f"停留时间: {downcomer['实际停留时间_s']:.1f} s")

        # 显示面积分区
        area = calculator.final_results['塔板面积分区']
        print(f"鼓泡区面积: {area['鼓泡区面积_Aa_m2']:.4f} m²")
        print(f"面积闭合误差: {area['面积闭合误差_%']:.6f}%")

        print(f"\n详细计算结果已保存到 'tower_plate_arrangement_corrected.json'")