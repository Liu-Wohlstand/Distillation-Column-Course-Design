"""
浮阀塔流体力学验算程序（分别计算精馏段和提馏段）
功能：
1. 塔板压降计算 (干板压降 + 液层压降)
2. 液沫夹带计算
3. 漏液点气速计算
4. 降液管液泛验算
5. 精馏段和提馏段负荷性能图绘制

浮阀塔特点：
1. 无筛孔，改为阀孔和浮阀
2. 干板压降计算采用浮阀塔专用公式
3. 漏液点计算采用浮阀塔专用公式
4. 操作弹性计算考虑浮阀塔特性

优化内容：
1. 增加详细计算过程输出，包括公式、带数过程和单位
2. 优化负荷性能线方程推导过程输出
3. 统一格式，使输出更加清晰易读
4. 解决重复计算和输出问题
"""

import json
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, brentq
import warnings
import matplotlib
from matplotlib.font_manager import FontProperties

times_new_roman = FontProperties(family='Times New Roman')
matplotlib.rcParams['font.family'] = 'SimSun'
matplotlib.rcParams['axes.unicode_minus'] = False

warnings.filterwarnings('ignore')


class FloatValveTowerHydrodynamicCheck:
    """浮阀塔流体力学验算器（分别计算精馏段和提馏段）"""

    def __init__(self, plate_arrangement_file='tower_plate_arrangement_optimized.json',
                 diameter_design_file='tower_diameter_approximate_optimized.json',
                 properties_file='tower_properties_results.json'):
        """
        初始化浮阀塔流体力学验算器

        参数：
        plate_arrangement_file: 塔板布置设计文件
        diameter_design_file: 塔径设计文件
        properties_file: 塔工艺特性文件
        """
        self.load_data(plate_arrangement_file, diameter_design_file, properties_file)

        # 物理常数
        self.g = 9.81  # 重力加速度，m/s²

        # 缓存已计算的值，避免重复计算
        self._cached_values = {
            'C0': None,  # 干板阻力系数
            'how': {},  # 堰上液层高度缓存
            'u_ow': {}  # 漏液点气速缓存
        }

        # 验算结果存储 - 使用统一的中文键名
        self.check_results = {'精馏段': {}, '提馏段': {}}
        self.load_performance_data = {'精馏段': {}, '提馏段': {}}

        # 详细输出开关
        self.verbose = True

    def load_data(self, plate_file, diameter_file, properties_file):
        """加载所有设计数据"""
        try:
            # 加载塔板布置设计结果
            with open(plate_file, 'r', encoding='utf-8') as f:
                self.plate_data = json.load(f)

            # 加载塔径设计结果
            with open(diameter_file, 'r', encoding='utf-8') as f:
                self.diameter_data = json.load(f)

            # 加载工艺特性数据
            with open(properties_file, 'r', encoding='utf-8') as f:
                self.properties_data = json.load(f)

            print("✓ 成功加载所有设计数据")

            # 提取浮阀塔关键参数
            self.extract_float_valve_plate_parameters()

            # 分别提取精馏段和提馏段工艺数据
            self.extract_section_parameters()

        except FileNotFoundError as e:
            print(f"✗ 错误: 未找到数据文件 - {e}")
            print("请确保以下文件存在:")
            print(f"  1. {plate_file}")
            print(f"  2. {diameter_file}")
            print(f"  3. {properties_file}")
            raise

    def extract_float_valve_plate_parameters(self):
        """提取浮阀塔关键设计参数（塔板结构参数）"""
        # 塔板基本信息
        self.D = self.plate_data['设计信息']['塔径_m']  # 塔径，m
        self.control_section = self.plate_data['设计信息']['控制塔段']  # 控制塔段

        # 塔板面积分区参数
        area_data = self.plate_data['塔板面积分区']
        self.A_t = area_data['塔板总面积_m2']  # 塔板总面积，m²
        self.A_a = area_data['鼓泡区面积_m2']  # 鼓泡区面积，m²

        # 注意：塔板面积分区中的降液管面积是投影面积，不用于流体计算
        self.A_f_projected = area_data['降液管面积_m2']  # 降液管投影面积，m²（仅用于参考）

        # 浮阀布置参数（假设设计文件中有这些参数）
        valve_data = self.plate_data.get('浮阀布置', {})

        # 阀孔参数
        self.N_valves = valve_data.get('阀孔数', 100)  # 阀孔数
        self.dv = valve_data.get('阀孔直径_m', 0.039)  # 阀孔直径，m
        self.δ = valve_data.get('塔板厚度_m', 0.003)  # 塔板厚度，m
        self.φ = valve_data.get('开孔率_%', 12) / 100  # 开孔率（小数）

        # 单个阀孔面积
        self.Av_single = math.pi * (self.dv / 2) ** 2  # m²
        self.Av_total = self.N_valves * self.Av_single  # 总阀孔面积，m²

        # 浮阀参数
        self.valve_type = valve_data.get('浮阀类型', 'F1')  # 浮阀类型
        self.valve_weight = valve_data.get('阀重_kg', 0.034)  # 浮阀重量，kg

        # 溢流装置参数
        overflow_data = self.plate_data['溢流装置']
        self.lw = overflow_data['堰长_m']  # 堰长，m
        self.hw = overflow_data['堰高_m']  # 堰高，m
        self.how = overflow_data['堰上液层高度_m']  # 堰上液层高度，m
        self.h0 = overflow_data['降液管底隙高度_m']  # 降液管底隙高度，m

        # 修改：使用溢流装置中的降液管面积（液体流通面积）
        self.A_f = overflow_data['降液管面积_m2']  # 降液管流通面积，m²

        # 板间距（使用最终调整后的板间距）
        if '最终塔径下详细参数' in self.diameter_data:
            # 使用经过调整后的最终板间距
            self.HT_rectifying = self.diameter_data['最终塔径下详细参数']['精馏段']['HT']
            self.HT_stripping = self.diameter_data['最终塔径下详细参数']['提馏段']['HT']
        else:
            # 如果不存在最终参数，则使用设计参数
            self.HT_rectifying = self.diameter_data['设计参数']['板间距']['精馏段_m']
            self.HT_stripping = self.diameter_data['设计参数']['板间距']['提馏段_m']

        # 计算清液层高度（使用结构参数计算）
        self.hL = self.hw + self.how  # 清液层高度，m

        print(f"\n浮阀塔设计参数:")
        print(f"塔径: {self.D:.3f} m")
        print(f"阀孔数: {self.N_valves}")
        print(f"阀孔直径: {self.dv * 1000:.1f} mm")
        print(f"开孔率: {self.φ * 100:.1f}%")
        print(f"浮阀类型: {self.valve_type}")
        print(f"阀重: {self.valve_weight:.3f} kg")
        print(f"堰长: {self.lw:.3f} m, 堰高: {self.hw * 1000:.1f} mm")
        print(f"清液层高度: {self.hL * 1000:.1f} mm")
        print(f"精馏段板间距: {self.HT_rectifying:.3f} m")
        print(f"提馏段板间距: {self.HT_stripping:.3f} m")
        print(f"降液管面积（流通）: {self.A_f:.4f} m²")
        print(f"降液管面积（投影）: {self.A_f_projected:.4f} m²")

    def extract_section_parameters(self):
        """分别提取精馏段和提馏段的工艺参数"""
        # 精馏段数据
        rect_data = self.properties_data['精馏段结果']
        self.rectifying = {
            'V_vol_m3h': rect_data['V_vol_m3h'],
            'V_vol_m3s': rect_data['V_vol_m3h'] / 3600,
            'ρ_v': rect_data['rho_vapor_kgm3'],
            'ρ_L': rect_data['rho_liquid_kgm3'],
            'σ': rect_data['sigma_mNm'] / 1000,
            'L_vol_Lh': rect_data.get('L_vol_Lh', 15000),
            'HT': self.HT_rectifying
        }
        self.rectifying['L_vol_m3h'] = self.rectifying['L_vol_Lh'] / 1000
        self.rectifying['L_vol_m3s'] = self.rectifying['L_vol_m3h'] / 3600

        # 提馏段数据
        strip_data = self.properties_data['提馏段结果']
        self.stripping = {
            'V_vol_m3h': strip_data['V_vol_m3h'],
            'V_vol_m3s': strip_data['V_vol_m3h'] / 3600,
            'ρ_v': strip_data['rho_vapor_kgm3'],
            'ρ_L': strip_data['rho_liquid_kgm3'],
            'σ': strip_data['sigma_mNm'] / 1000,
            'L_vol_Lh': strip_data.get('L_vol_Lh', 20000),  # 提馏段液相流量通常更大
            'HT': self.HT_stripping
        }
        self.stripping['L_vol_m3h'] = self.stripping['L_vol_Lh'] / 1000
        self.stripping['L_vol_m3s'] = self.stripping['L_vol_m3h'] / 3600

        print(f"\n精馏段工艺参数:")
        print(f"  气相流量: {self.rectifying['V_vol_m3s']:.6f} m³/s")
        print(f"  液相流量: {self.rectifying['L_vol_m3s']:.6f} m³/s")
        print(f"  气相密度: {self.rectifying['ρ_v']:.4f} kg/m³")
        print(f"  液相密度: {self.rectifying['ρ_L']:.1f} kg/m³")
        print(f"  表面张力: {self.rectifying['σ'] * 1000:.2f} mN/m")

        print(f"\n提馏段工艺参数:")
        print(f"  气相流量: {self.stripping['V_vol_m3s']:.6f} m³/s")
        print(f"  液相流量: {self.stripping['L_vol_m3s']:.6f} m³/s")
        print(f"  气相密度: {self.stripping['ρ_v']:.4f} kg/m³")
        print(f"  液相密度: {self.stripping['ρ_L']:.1f} kg/m³")
        print(f"  表面张力: {self.stripping['σ'] * 1000:.2f} mN/m")

    def print_formula(self, formula_name, formula):
        """打印公式"""
        if self.verbose:
            print(f"  └ 公式: {formula_name}: {formula}")

    def print_calculation(self, step_name, calculation, result, unit=""):
        """打印计算步骤"""
        if self.verbose:
            print(f"    ├ {step_name}: {calculation} = {result:.4f} {unit}")

    def calculate_dry_plate_resistance_coefficient(self, skip_print=False):
        """
        计算干板阻力系数C0（浮阀塔专用）
        公式: C0 = 5.34 + 2.88/Re^0.5  当阀全开时，C0 ≈ 0.7-0.8
        简化取C0 = 0.75（常见值）
        """
        # 检查缓存
        if self._cached_values['C0'] is not None:
            if not skip_print and self.verbose:
                print(f"\n  使用缓存的干板阻力系数C0: {self._cached_values['C0']:.4f}")
            return self._cached_values['C0']

        if self.verbose and not skip_print:
            print(f"\n  计算干板阻力系数C0（浮阀塔）:")

        # 浮阀塔干板阻力系数通常取0.7-0.8，这里取0.75
        C0 = 0.75

        if self.verbose and not skip_print:
            self.print_formula("干板阻力系数C0", "C0 = 0.75 (浮阀塔经验值)")
            self.print_calculation("C0", "0.75", C0, "")

        # 缓存结果
        self._cached_values['C0'] = C0
        return C0

    def calculate_plate_pressure_drop(self, V_m3s, ρ_v, ρ_L, σ, HT, section_name="精馏段"):
        """
        计算塔板压降 hp = hc + hl （浮阀塔无表面张力压降）

        参数:
        V_m3s: 气相体积流量，m³/s
        ρ_v: 气相密度，kg/m³
        ρ_L: 液相密度，kg/m³
        σ: 表面张力，N/m
        HT: 板间距，m
        section_name: 塔段名称

        返回:
        hp: 总压降，m液柱
        hc: 干板压降，m液柱
        hl: 液层压降，m液柱
        hσ: 表面张力压降，m液柱（浮阀塔为0）
        """
        if self.verbose:
            print(f"\n  计算{section_name}塔板压降（浮阀塔）:")

        # 1. 计算阀孔气速
        u0 = V_m3s / self.Av_total  # 阀孔气速，m/s
        self.print_calculation("阀孔气速u0", f"V/Av_total = {V_m3s:.6f}/{self.Av_total:.6f}", u0, "m/s")

        # 2. 计算干板阻力系数C0
        C0 = self.calculate_dry_plate_resistance_coefficient()

        # 3. 计算干板压降 hc（浮阀塔公式）
        self.print_formula("干板压降hc", "hc = ξ × (ρ_v/ρ_L) × (u0²/(2g))")

        # 浮阀塔阻力系数ξ = 5.34（阀全开）
        ξ = 5.34

        hc = ξ * (ρ_v / ρ_L) * (u0 ** 2 / (2 * self.g))
        self.print_calculation("干板压降hc", f"5.34 × ({ρ_v:.4f}/{ρ_L:.1f}) × ({u0:.4f}²/(2×9.81))", hc, "m液柱")

        # 4. 计算液层压降 hl
        # 首先计算Fa（气体动能因数）
        u_a = V_m3s / self.A_a  # 鼓泡区气速，m/s
        F_a = u_a * math.sqrt(ρ_v)  # 气体动能因数

        self.print_calculation("鼓泡区气速u_a", f"V/A_a = {V_m3s:.6f}/{self.A_a:.4f}", u_a, "m/s")
        self.print_calculation("气体动能因数F_a", f"u_a × √ρ_v = {u_a:.4f} × √{ρ_v:.4f}", F_a, "m/s·√(kg/m³)")

        # 计算充气系数 ε0
        self.print_formula("充气系数ε0", "ε0 = 0.971 - 0.355*F_a + 0.0757*F_a²")

        ε0 = 0.971 - 0.355 * F_a + 0.0757 * F_a ** 2
        if self.verbose:
            print(f"      ε0 = 0.971 - 0.355×{F_a:.4f} + 0.0757×{F_a:.4f}²")
            print(f"         = 0.971 - {0.355 * F_a:.4f} + {0.0757 * F_a ** 2:.4f}")
            print(f"         = {ε0:.4f}")

        # 限制ε0在合理范围内
        ε0 = max(0.4, min(ε0, 0.7))
        if self.verbose:
            print(f"      限制在0.4~0.7之间: ε0 = {ε0:.4f}")

        # 液层压降 hl = ε0 * hL
        self.print_formula("液层压降hl", "hl = ε0 × hL")
        hl = ε0 * self.hL
        self.print_calculation("液层压降hl", f"{ε0:.4f} × {self.hL:.4f}", hl, "m液柱")

        # 5. 表面张力压降 hσ（浮阀塔通常不考虑）
        hσ = 0
        if self.verbose:
            print(f"    表面张力压降hσ: 浮阀塔通常不考虑，hσ = 0 m液柱")

        # 6. 总压降 hp = hc + hl + hσ
        self.print_formula("总压降hp", "hp = hc + hl + hσ")
        hp = hc + hl + hσ
        self.print_calculation("总压降hp", f"{hc:.4f} + {hl:.4f} + {hσ:.4f}", hp, "m液柱")

        return hp, hc, hl, hσ

    def calculate_entrainment(self, V_m3s, L_m3s, ρ_v, ρ_L, σ, HT, section_name="精馏段"):
        """
        计算液沫夹带量（使用亨特经验公式） - 浮阀塔与筛板塔相同

        参数:
        V_m3s: 气相体积流量，m³/s
        L_m3s: 液相体积流量，m³/s
        ρ_v: 气相密度，kg/m³
        ρ_L: 液相密度，kg/m³
        σ: 表面张力，N/m
        HT: 板间距，m
        section_name: 塔段名称

        返回:
        ev: 液沫夹带量，kg液体/kg气体
        """
        if self.verbose:
            print(f"\n  计算{section_name}液沫夹带量:")

        # 计算实际空塔气速（使用投影面积计算塔板有效面积）
        self.print_formula("空塔气速u", "u = V / (A_t - A_f_projected)")
        u = V_m3s / (self.A_t - self.A_f_projected)  # 有效塔截面积，m²
        self.print_calculation("空塔气速u",
                               f"{V_m3s:.6f} / ({self.A_t:.4f} - {self.A_f_projected:.4f})",
                               u, "m/s")

        # 计算液气比参数FLV
        self.print_formula("液气比参数FLV", "FLV = (L/V) × √(ρ_L/ρ_v)")
        FLV = (L_m3s / V_m3s) * math.sqrt(ρ_L / ρ_v) if V_m3s > 0 else 0
        self.print_calculation("FLV",
                               f"({L_m3s:.6f}/{V_m3s:.6f}) × √({ρ_L:.1f}/{ρ_v:.4f})",
                               FLV, "")

        # 计算板间距减去清液层高度
        HT_minus_hL = HT - self.hL
        self.print_calculation("HT - hL", f"{HT:.4f} - {self.hL:.4f}", HT_minus_hL, "m")

        # 计算C20经验系数
        self.print_formula("C20经验系数", "C20 = 0.0162 - 0.0648*X + 0.181*Y + 0.0162*X² - 0.139*X*Y + 0.185*Y²")
        X = FLV
        Y = HT_minus_hL
        C20 = (0.0162 - 0.0648 * X + 0.181 * Y +
               0.0162 * X ** 2 - 0.139 * X * Y + 0.185 * Y ** 2)

        if self.verbose:
            print(f"    代入计算:")
            print(f"      X = FLV = {FLV:.4f}")
            print(f"      Y = HT - hL = {HT_minus_hL:.4f}")
            print(
                f"      C20 = 0.0162 - 0.0648×{FLV:.4f} + 0.181×{HT_minus_hL:.4f} + 0.0162×{FLV:.4f}² - 0.139×{FLV:.4f}×{HT_minus_hL:.4f} + 0.185×{HT_minus_hL:.4f}²")
            print(
                f"          = {0.0162:.4f} - {0.0648 * FLV:.4f} + {0.181 * HT_minus_hL:.4f} + {0.0162 * FLV ** 2:.4f} - {0.139 * FLV * HT_minus_hL:.4f} + {0.185 * HT_minus_hL ** 2:.4f}")
            print(f"          = {C20:.4f}")

        # 表面张力校正
        self.print_formula("表面张力校正C", "C = C20 × (σ/0.02)^0.2")
        C = C20 * (σ * 1000 / 20) ** 0.2  # σ单位转换为mN/m
        self.print_calculation("校正系数C",
                               f"{C20:.4f} × ({σ * 1000:.2f}/20)^0.2 = {C20:.4f} × {((σ * 1000) / 20) ** 0.2:.4f}",
                               C, "")

        # 计算液泛气速
        self.print_formula("液泛气速u_f", "u_f = C × √[(ρ_L - ρ_v)/ρ_v]")
        u_f = C * math.sqrt((ρ_L - ρ_v) / ρ_v)
        self.print_calculation("液泛气速u_f",
                               f"{C:.4f} × √[({ρ_L:.1f} - {ρ_v:.4f})/{ρ_v:.4f}] = {C:.4f} × √{((ρ_L - ρ_v) / ρ_v):.4f}",
                               u_f, "m/s")

        # 亨特经验公式计算液沫夹带量
        σ_mN_m = σ * 1000  # 转换为mN/m

        if HT_minus_hL <= 0:
            ev = 0
            if self.verbose:
                print(f"    HT - hL ≤ 0, 液沫夹带量 ev = 0")
        else:
            self.print_formula("亨特经验公式", "ev = (5.7×10⁻⁶/σ) × (u/(HT-hL))^3.2")
            ev = (5.7e-6 / σ_mN_m) * (u / (HT - self.hL)) ** 3.2
            self.print_calculation("液沫夹带量ev",
                                   f"(5.7×10⁻⁶/{σ_mN_m:.2f}) × ({u:.4f}/{HT_minus_hL:.4f})^3.2",
                                   ev, "kg液体/kg气体")

        # 限制液沫夹带量在合理范围内
        ev = max(0, ev)

        # 验算标准：ev < 0.1 kg液体/kg气体
        ev_check = ev < 0.1

        if self.verbose:
            print(f"    限值: ev < 0.1 kg液体/kg气体")
            print(f"    结果: ev = {ev:.6f} {'<' if ev_check else '≥'} 0.1")
            print(f"    {'✓ 合格' if ev_check else '✗ 不合格'}")

        return ev, ev_check

    def calculate_how(self, L_m3s, skip_print=False):
        """
        计算堰上液层高度（带缓存）
        """
        # 生成缓存键
        cache_key = f"{L_m3s:.8f}"

        if cache_key in self._cached_values['how']:
            if not skip_print and self.verbose:
                print(f"    使用缓存的堰上液层高度: {self._cached_values['how'][cache_key]:.6f} m")
            return self._cached_values['how'][cache_key]

        # 计算堰上液层高度（随液相流量变化）
        if not skip_print and self.verbose:
            self.print_formula("堰上液层高度how", "how = 2.84×10⁻³ × E × (L_h/l_w)^(2/3)")

        E = 1.0  # 液流收缩系数
        L_h = L_m3s * 3600  # 转换为m³/h
        how = 2.84 / 1000 * E * (L_h / self.lw) ** (2 / 3)  # m

        if not skip_print and self.verbose:
            self.print_calculation("堰上液层高度how",
                                   f"2.84×10⁻³ × 1.0 × ({L_h:.1f}/{self.lw:.3f})^(2/3)",
                                   how, "m")

        # 缓存结果
        self._cached_values['how'][cache_key] = how
        return how

    def calculate_leakage_point_velocity(self, L_m3s, ρ_v, ρ_L, σ, section_name="精馏段", skip_print=False):
        """
        计算漏液点气速（浮阀塔专用公式）

        浮阀塔漏液点气速公式：
        u_ow = (4.4/C0) × √[(0.0056 + 0.13hL - hσ) × (ρ_L/ρ_v)]
        或简化为：u_ow = 5.0 × √(ρ_L/ρ_v)

        参数:
        L_m3s: 液相体积流量，m³/s
        ρ_v: 气相密度，kg/m³
        ρ_L: 液相密度，kg/m³
        σ: 表面张力，N/m
        section_name: 塔段名称
        skip_print: 是否跳过打印（用于避免重复输出）

        返回:
        u_ow: 漏液点气速，m/s
        """
        # 生成缓存键
        cache_key = f"{L_m3s:.8f}_{ρ_v:.8f}_{ρ_L:.8f}_{σ:.8f}"

        if cache_key in self._cached_values['u_ow']:
            if not skip_print and self.verbose:
                print(f"    使用缓存的漏液点气速: {self._cached_values['u_ow'][cache_key]:.6f} m/s")
            return self._cached_values['u_ow'][cache_key]

        if not skip_print and self.verbose:
            print(f"\n  计算{section_name}漏液点气速（浮阀塔）:")

        # 计算堰上液层高度
        how = self.calculate_how(L_m3s, skip_print=skip_print)

        # 清液层高度
        if not skip_print and self.verbose:
            self.print_formula("清液层高度hL", "hL = hw + how")
        hL = self.hw + how
        if not skip_print and self.verbose:
            self.print_calculation("清液层高度hL", f"{self.hw:.4f} + {how:.4f}", hL, "m")

        # 浮阀塔漏液点气速简化公式
        if not skip_print and self.verbose:
            self.print_formula("漏液点气速u_ow", "u_ow = 5.0 × √(ρ_L/ρ_v)  (浮阀塔经验公式)")

        u_ow = 5.0 * math.sqrt(ρ_L / ρ_v)

        if not skip_print and self.verbose:
            self.print_calculation("漏液点气速u_ow",
                                   f"5.0 × √({ρ_L:.1f}/{ρ_v:.4f}) = 5.0 × √{ρ_L / ρ_v:.4f}",
                                   u_ow, "m/s")

        # 缓存结果
        self._cached_values['u_ow'][cache_key] = u_ow
        return u_ow

    def calculate_downcomer_flooding(self, V_m3s, L_m3s, ρ_v, ρ_L, σ, HT, section_name="精馏段"):
        """
        计算降液管液泛（浮阀塔与筛板塔相同）

        参数:
        V_m3s: 气相体积流量，m³/s
        L_m3s: 液相体积流量，m³/s
        ρ_v: 气相密度，kg/m³
        ρ_L: 液相密度，kg/m³
        σ: 表面张力，N/m
        HT: 板间距，m
        section_name: 塔段名称

        返回:
        Hd: 降液管内清液层高度，m
        Hd_max: 允许最大清液层高度，m
        flooding_check: 是否发生液泛
        """
        if self.verbose:
            print(f"\n  计算{section_name}降液管液泛:")

        # 1. 计算塔板压降 hp
        hp, hc, hl, hσ = self.calculate_plate_pressure_drop(V_m3s, ρ_v, ρ_L, σ, HT, section_name)

        # 2. 计算清液层高度 hL
        how = self.calculate_how(L_m3s)
        hL = self.hw + how

        if self.verbose:
            self.print_formula("清液层高度hL", "hL = hw + how")
            self.print_calculation("清液层高度hL", f"{self.hw:.4f} + {how:.4f}", hL, "m")

        # 3. 计算液体通过降液管的压头损失 hd
        self.print_formula("降液管压头损失hd", "hd = 0.153 × (L_s/(l_w × h_0))²")
        if L_m3s > 0 and self.lw * self.h0 > 0:
            hd = 0.153 * (L_m3s / (self.lw * self.h0)) ** 2
            self.print_calculation("降液管压头损失hd",
                                   f"0.153 × ({L_m3s:.6f}/({self.lw:.3f}×{self.h0:.3f}))²",
                                   hd, "m")
        else:
            hd = 0
            print(f"    hd = 0 (L_s ≤ 0 或 l_w×h_0 ≤ 0)")

        # 4. 降液管内清液层高度 Hd = hp + hL + hd
        self.print_formula("降液管液层高度Hd", "Hd = hp + hL + hd")
        Hd = hp + hL + hd
        self.print_calculation("降液管液层高度Hd", f"{hp:.4f} + {hL:.4f} + {hd:.4f}", Hd, "m")

        # 5. 计算允许最大清液层高度 Hd_max = Φ * (HT + hw)
        self.print_formula("允许最大液层高度Hd_max", "Hd_max = Φ × (HT + hw)")
        Φ = 0.5  # 一般物系
        Hd_max = Φ * (HT + self.hw)
        self.print_calculation("允许最大液层高度Hd_max",
                               f"{Φ:.1f} × ({HT:.4f} + {self.hw:.4f})",
                               Hd_max, "m")

        # 6. 液泛验算
        flooding_check = Hd <= Hd_max

        if self.verbose:
            print(f"    液泛条件: Hd ≤ Hd_max")
            print(f"    验算: {Hd:.4f} {'≤' if flooding_check else '>'} {Hd_max:.4f}")
            print(f"    {'✓ 合格' if flooding_check else '✗ 不合格'}")

        return Hd, Hd_max, flooding_check

    def calculate_entrainment_limit_line(self, ρ_v, ρ_L, σ, HT, section_name="精馏段"):
        """
        计算过量液沫夹带线（ev = 0.1 kg液体/kg气体）
        返回: 液相流量Ls与气相流量Vs的关系
        """
        print(f"\n  ┌─────────────────────────────────────")
        print(f"  │ 计算{section_name}过量液沫夹带线 (ev = 0.1)")
        print(f"  └─────────────────────────────────────")

        if self.verbose:
            print(f"    方程推导:")
            print(f"      1. 液沫夹带量公式: ev = (5.7×10⁻⁶/σ) × (u/(HT-hL))^3.2")
            print(f"      2. 空塔气速: u = Vs / (A_t - A_f_projected)")
            print(f"      3. 清液层高度: hL = hw + how")
            print(f"      4. 堰上液层高度: how = 2.84×10⁻³ × E × (L_h/l_w)^(2/3)")
            print(f"      5. 令 ev = 0.1，求解 Vs 与 Ls 的关系")

        # 创建液相流量范围 (m³/s)
        if section_name == "精馏段":
            Ls_design = self.rectifying['L_vol_m3s']
        else:
            Ls_design = self.stripping['L_vol_m3s']

        Ls_min = Ls_design * 0.1
        Ls_max = min(Ls_design * 3, 0.02)  # 限制最大流量
        Ls_values = np.linspace(Ls_min, Ls_max, 50)
        Vs_values = []

        if self.verbose:
            print(f"\n    计算过程 (示例点 Ls = {Ls_values[0]:.6f} m³/s):")

        for i, Ls in enumerate(Ls_values):
            # 对于每个Ls，求解满足ev=0.1的Vs
            def ev_equation(Vs):
                # 计算空塔气速（使用投影面积计算塔板有效面积）
                u = Vs / (self.A_t - self.A_f_projected)

                # 亨特经验公式计算ev
                if HT - self.hL <= 0:
                    ev = 0
                else:
                    # 修正的亨特公式
                    σ_mN_m = σ * 1000
                    ev = (5.7e-6 / σ_mN_m) * (u / (HT - self.hL)) ** 3.2

                return ev - 0.1  # 目标: ev = 0.1

            # 使用数值方法求解
            try:
                # 设置Vs搜索范围
                if section_name == "精馏段":
                    Vs_design = self.rectifying['V_vol_m3s']
                else:
                    Vs_design = self.stripping['V_vol_m3s']

                Vs_min = 0.01 * Vs_design
                Vs_max = 3 * Vs_design

                # 检查函数值符号
                f_min = ev_equation(Vs_min)
                f_max = ev_equation(Vs_max)

                if f_min * f_max < 0:
                    # 有根，使用brentq方法
                    Vs_root = brentq(ev_equation, Vs_min, Vs_max)
                    Vs_values.append(Vs_root)

                    if self.verbose and i == 0:
                        print(f"      Ls = {Ls:.6f} m³/s:")
                        print(f"        方程: ev(Vs) - 0.1 = 0")
                        print(f"        搜索区间: [{Vs_min:.6f}, {Vs_max:.6f}] m³/s")
                        print(f"        求解得: Vs = {Vs_root:.6f} m³/s")
                        print(f"        验证: ev = {ev_equation(Vs_root) + 0.1:.6f}")
                else:
                    # 无根，使用近似值
                    Vs_values.append(Vs_max if f_max < 0 else Vs_min)
                    if self.verbose and i == 0:
                        print(f"      Ls = {Ls:.6f} m³/s:")
                        print(f"        无解，使用边界值: Vs = {Vs_values[-1]:.6f} m³/s")
            except:
                # 计算失败，使用设计值
                if section_name == "精馏段":
                    Vs_values.append(self.rectifying['V_vol_m3s'])
                else:
                    Vs_values.append(self.stripping['V_vol_m3s'])

        return Ls_values, np.array(Vs_values)

    def calculate_flooding_limit_line(self, ρ_v, ρ_L, σ, HT, section_name="精馏段"):
        """
        计算溢流液泛线（Hd = Hd_max）
        返回: 气相流量Vs与液相流量Ls的关系
        """
        print(f"\n  ┌─────────────────────────────────────")
        print(f"  │ 计算{section_name}溢流液泛线 (Hd = Hd_max)")
        print(f"  └─────────────────────────────────────")

        if self.verbose:
            print(f"    方程推导:")
            print(f"      1. 降液管液层高度: Hd = hp + hL + hd")
            print(f"      2. 塔板压降: hp = hc + hl + hσ")
            print(f"        其中: hc = ξ × (ρ_v/ρ_L) × (u0²/(2g))  (浮阀塔公式)")
            print(f"              hl = ε0 × hL")
            print(f"              ε0 = 0.971 - 0.355×F_a + 0.0757×F_a²")
            print(f"              hσ = 0 (浮阀塔不考虑)")
            print(f"      3. 降液管压头损失: hd = 0.153 × (L_s/(l_w × h_0))²")
            print(f"      4. 允许最大高度: Hd_max = Φ × (HT + hw)")
            print(f"      5. 令 Hd = Hd_max，求解 Vs 与 Ls 的关系")

        # 允许最大清液层高度
        Φ = 0.5
        Hd_max = Φ * (HT + self.hw)

        if self.verbose:
            print(f"      Hd_max = {Φ:.1f} × ({HT:.3f} + {self.hw:.3f}) = {Hd_max:.4f} m")

        # 创建液相流量范围 (m³/s)
        if section_name == "精馏段":
            Ls_design = self.rectifying['L_vol_m3s']
        else:
            Ls_design = self.stripping['L_vol_m3s']

        Ls_min = max(1e-6, Ls_design * 0.1)
        Ls_max = min(Ls_design * 3, 0.02)
        Ls_values = np.linspace(Ls_min, Ls_max, 50)
        Vs_values = []

        if self.verbose:
            print(f"\n    计算过程 (示例点 Ls = {Ls_values[0]:.6f} m³/s):")

        for i, Ls in enumerate(Ls_values):
            # 对于每个Ls，求解满足Hd=Hd_max的Vs
            def flooding_equation(Vs):
                # 计算塔板压降（浮阀塔公式）
                u0 = Vs / self.Av_total

                # 干板压降（浮阀塔）
                ξ = 5.34
                hc = ξ * (ρ_v / ρ_L) * (u0 ** 2 / (2 * self.g))

                # 鼓泡区气速和Fa
                u_a = Vs / self.A_a
                F_a = u_a * math.sqrt(ρ_v)

                # 充气系数
                ε0 = 0.971 - 0.355 * F_a + 0.0757 * F_a ** 2
                ε0 = max(0.4, min(ε0, 0.7))

                # 堰上液层高度（使用缓存函数）
                how = self.calculate_how(Ls, skip_print=True)
                hL = self.hw + how

                # 液层压降
                hl = ε0 * hL

                # 表面张力压降（浮阀塔不考虑）
                hσ = 0

                # 总压降
                hp = hc + hl + hσ

                # 降液管压头损失
                if Ls > 0 and self.lw * self.h0 > 0:
                    hd = 0.153 * (Ls / (self.lw * self.h0)) ** 2
                else:
                    hd = 0

                # 降液管内清液层高度
                Hd = hp + hL + hd

                return Hd - Hd_max  # 目标: Hd = Hd_max

            # 使用数值方法求解
            try:
                # 设置Vs搜索范围
                if section_name == "精馏段":
                    Vs_design = self.rectifying['V_vol_m3s']
                else:
                    Vs_design = self.stripping['V_vol_m3s']

                Vs_min = 0.01 * Vs_design
                Vs_max = 3 * Vs_design

                # 检查函数值符号
                f_min = flooding_equation(Vs_min)
                f_max = flooding_equation(Vs_max)

                if f_min * f_max < 0:
                    # 有根，使用brentq方法
                    Vs_root = brentq(flooding_equation, Vs_min, Vs_max)
                    Vs_values.append(Vs_root)

                    if self.verbose and i == 0:
                        print(f"      Ls = {Ls:.6f} m³/s:")
                        print(f"        方程: Hd(Vs) - Hd_max = 0")
                        print(f"        搜索区间: [{Vs_min:.6f}, {Vs_max:.6f}] m³/s")
                        print(f"        求解得: Vs = {Vs_root:.6f} m³/s")
                else:
                    # 无根，使用近似值
                    Vs_values.append(Vs_max if f_max < 0 else Vs_min)
                    if self.verbose and i == 0:
                        print(f"      Ls = {Ls:.6f} m³/s:")
                        print(f"        无解，使用边界值: Vs = {Vs_values[-1]:.6f} m³/s")
            except Exception as e:
                # 计算失败，使用近似值
                if section_name == "精馏段":
                    Vs_values.append(self.rectifying['V_vol_m3s'])
                else:
                    Vs_values.append(self.stripping['V_vol_m3s'])

        return Ls_values, np.array(Vs_values)

    def calculate_liquid_upper_limit_line(self):
        """
        计算液相负荷上限线（降液管内停留时间τ = 3~5s）
        返回: 液相流量上限值

        修正：使用流通面积self.A_f而不是投影面积
        """
        print(f"\n  ┌─────────────────────────────────────")
        print(f"  │ 计算液相负荷上限线 (停留时间τ = 5s)")
        print(f"  └─────────────────────────────────────")

        if self.verbose:
            print(f"    方程推导:")
            print(f"      1. 停留时间公式: τ = A_f × HT / L_s")
            print(f"      2. 取最小停留时间 τ_min = 5 s (最严格)")
            print(f"      3. 最大液相流量: L_s,max = A_f × HT / τ_min")
            print(f"      4. 使用较小的板间距以确保安全")

        # 取停留时间τ = 5s（最严格）
        τ_min = 5  # s

        # 计算最大液相流量 Ls,max = Af * HT / τ
        # 注意：这里使用较小的板间距以确保安全（取精馏段和提馏段的最小值）
        HT_min = min(self.HT_rectifying, self.HT_stripping)

        # 修正：使用流通面积self.A_f
        self.print_formula("最大液相流量Ls,max", "Ls,max = A_f × HT_min / τ_min")
        Ls_max = self.A_f * HT_min / τ_min

        self.print_calculation("Ls,max",
                               f"{self.A_f:.4f} × {HT_min:.3f} / {τ_min}",
                               Ls_max, "m³/s")

        # 这是一条垂直线，所有Vs对应同一个Ls_max
        Ls_values = np.array([Ls_max, Ls_max])
        Vs_min = 0
        Vs_max = max(self.rectifying['V_vol_m3s'], self.stripping['V_vol_m3s']) * 3
        Vs_values = np.array([Vs_min, Vs_max])

        return Ls_values, Vs_values

    def calculate_leakage_limit_line(self, ρ_v, ρ_L, σ, section_name="精馏段"):
        """
        计算漏液线（u0 = u_ow）
        返回: 气相流量Vs与液相流量Ls的关系
        """
        print(f"\n  ┌─────────────────────────────────────")
        print(f"  │ 计算{section_name}漏液线 (u0 = u_ow)")
        print(f"  └─────────────────────────────────────")

        if self.verbose:
            print(f"    方程推导:")
            print(f"      1. 阀孔气速: u0 = Vs / Av_total")
            print(f"      2. 漏液点气速: u_ow = 5.0 × √(ρ_L/ρ_v) (浮阀塔经验公式)")
            print(f"      3. 令 u0 = u_ow，则 Vs = u_ow × Av_total")

        # 创建液相流量范围 (m³/s)
        if section_name == "精馏段":
            Ls_design = self.rectifying['L_vol_m3s']
        else:
            Ls_design = self.stripping['L_vol_m3s']

        Ls_min = max(1e-6, Ls_design * 0.1)
        Ls_max = min(Ls_design * 3, 0.02)
        Ls_values = np.linspace(Ls_min, Ls_max, 50)
        Vs_values = []

        if self.verbose:
            print(f"\n    计算过程 (示例点 Ls = {Ls_values[0]:.6f} m³/s):")

        for i, Ls in enumerate(Ls_values):
            # 计算漏液点气速（使用缓存版本）
            u_ow = self.calculate_leakage_point_velocity(Ls, ρ_v, ρ_L, σ, section_name, skip_print=(i > 0))

            # 对应的气相流量 Vs = u_ow * Av_total
            self.print_formula("气相流量Vs", "Vs = u_ow × Av_total")
            Vs = u_ow * self.Av_total
            self.print_calculation("Vs", f"{u_ow:.4f} × {self.Av_total:.6f}", Vs, "m³/s")

            Vs_values.append(Vs)

            if self.verbose and i == 0:
                print(f"      Ls = {Ls:.6f} m³/s:")
                print(f"        漏液点气速 u_ow = {u_ow:.4f} m/s")
                print(f"        气相流量 Vs = {u_ow:.4f} × {self.Av_total:.6f} = {Vs:.6f} m³/s")

        return Ls_values, np.array(Vs_values)

    def calculate_liquid_lower_limit_line(self):
        """
        计算液相负荷下限线（how = 6mm）
        返回: 液相流量下限值
        """
        print(f"\n  ┌─────────────────────────────────────")
        print(f"  │ 计算液相负荷下限线 (how = 6mm)")
        print(f"  └─────────────────────────────────────")

        if self.verbose:
            print(f"    方程推导:")
            print(f"      1. 堰上液层高度公式: how = 2.84×10⁻³ × E × (L_h/l_w)^(2/3)")
            print(f"      2. 堰上液层高度下限 how_min = 0.006 m")
            print(f"      3. 反推最小液相流量: L_h,min = l_w × (how_min×1000/(2.84×E))^(3/2)")

        # 堰上液层高度下限 how_min = 0.006 m
        how_min = 0.006  # m

        # 计算最小液相流量
        self.print_formula("最小液相流量L_h,min", "L_h,min = l_w × [how_min×1000/(2.84×E)]^(3/2)")
        E = 1.0
        Lh_min = self.lw * (how_min * 1000 / (2.84 * E)) ** 1.5  # m³/h

        self.print_calculation("L_h,min",
                               f"{self.lw:.3f} × [{how_min * 1000:.1f}/(2.84×1.0)]^1.5",
                               Lh_min, "m³/h")

        Ls_min = Lh_min / 3600  # m³/s
        self.print_calculation("L_s,min", f"{Lh_min:.2f} / 3600", Ls_min, "m³/s")

        # 这是一条垂直线，所有Vs对应同一个Ls_min
        Ls_values = np.array([Ls_min, Ls_min])
        Vs_min = 0
        Vs_max = max(self.rectifying['V_vol_m3s'], self.stripping['V_vol_m3s']) * 3
        Vs_values = np.array([Vs_min, Vs_max])

        return Ls_values, Vs_values

    def generate_load_performance_diagram_for_section(self, section_name="精馏段"):
        """为指定塔段生成负荷性能图"""
        print(f"\n{'=' * 80}")
        print(f"{section_name}负荷性能图计算")
        print(f"{'=' * 80}")

        # 获取塔段参数
        if section_name == "精馏段":
            ρ_v = self.rectifying['ρ_v']
            ρ_L = self.rectifying['ρ_L']
            σ = self.rectifying['σ']
            HT = self.rectifying['HT']
            Ls_design = self.rectifying['L_vol_m3s']
            Vs_design = self.rectifying['V_vol_m3s']
        else:
            ρ_v = self.stripping['ρ_v']
            ρ_L = self.stripping['ρ_L']
            σ = self.stripping['σ']
            HT = self.stripping['HT']
            Ls_design = self.stripping['L_vol_m3s']
            Vs_design = self.stripping['V_vol_m3s']

        # 计算各条限制线
        print(f"\n计算{section_name}各条限制线...")

        # 1. 过量液沫夹带线
        print(f"\n1. 过量液沫夹带线 (ev = 0.1 kg液体/kg气体):")
        Ls_ent, Vs_ent = self.calculate_entrainment_limit_line(ρ_v, ρ_L, σ, HT, section_name)

        # 2. 溢流液泛线
        print(f"\n2. 溢流液泛线 (Hd = Hd_max):")
        Ls_flood, Vs_flood = self.calculate_flooding_limit_line(ρ_v, ρ_L, σ, HT, section_name)

        # 3. 液相负荷上限线（两条塔段共用）
        print(f"\n3. 液相负荷上限线 (停留时间τ = 5s):")
        Ls_upper, Vs_upper = self.calculate_liquid_upper_limit_line()

        # 4. 漏液线
        print(f"\n4. 漏液线 (u0 = u_ow):")
        Ls_leak, Vs_leak = self.calculate_leakage_limit_line(ρ_v, ρ_L, σ, section_name)

        # 5. 液相负荷下限线（两条塔段共用）
        print(f"\n5. 液相负荷下限线 (how = 6mm):")
        Ls_lower, Vs_lower = self.calculate_liquid_lower_limit_line()

        # 存储负荷性能数据
        self.load_performance_data[section_name] = {
            'entrainment_line': {'Ls': Ls_ent.tolist(), 'Vs': Vs_ent.tolist()},
            'flooding_line': {'Ls': Ls_flood.tolist(), 'Vs': Vs_flood.tolist()},
            'liquid_upper_line': {'Ls': Ls_upper.tolist(), 'Vs': Vs_upper.tolist()},
            'leakage_line': {'Ls': Ls_leak.tolist(), 'Vs': Vs_leak.tolist()},
            'liquid_lower_line': {'Ls': Ls_lower.tolist(), 'Vs': Vs_lower.tolist()},
            'operating_point': {'Ls': float(Ls_design), 'Vs': float(Vs_design)}
        }

        # 绘制负荷性能图
        self.plot_load_performance_diagram_for_section(
            Ls_ent, Vs_ent, Ls_flood, Vs_flood,
            Ls_upper, Vs_upper, Ls_leak, Vs_leak,
            Ls_lower, Vs_lower, Ls_design, Vs_design,
            section_name
        )

        # 计算操作弹性
        operational_elasticity = self.calculate_operational_elasticity_for_section(
            Ls_ent, Vs_ent, Ls_flood, Vs_flood,
            Ls_upper, Vs_upper, Ls_leak, Vs_leak,
            Ls_lower, Vs_lower, Ls_design, Vs_design,
            section_name
        )

        return operational_elasticity

    def plot_load_performance_diagram_for_section(self, Ls_ent, Vs_ent, Ls_flood, Vs_flood,
                                                  Ls_upper, Vs_upper, Ls_leak, Vs_leak,
                                                  Ls_lower, Vs_lower, Ls_design, Vs_design,
                                                  section_name="精馏段"):
        """绘制指定塔段的负荷性能图"""
        plt.figure(figsize=(12, 8))

        # 绘制各条限制线
        plt.plot(Ls_ent, Vs_ent, 'b-', linewidth=2, label='过量液沫夹带线 (ev=0.1)')
        plt.plot(Ls_flood, Vs_flood, 'r-', linewidth=2, label='溢流液泛线')
        plt.plot(Ls_upper, Vs_upper, 'g--', linewidth=2, label='液相负荷上限线 (τ=5s)')
        plt.plot(Ls_leak, Vs_leak, 'm-', linewidth=2, label='漏液线')
        plt.plot(Ls_lower, Vs_lower, 'c--', linewidth=2, label='液相负荷下限线 (how=6mm)')

        # 绘制操作点
        plt.plot(Ls_design, Vs_design, 'ko', markersize=10,
                 markerfacecolor='yellow', label=f'{section_name}操作点', markeredgewidth=2)

        # 标记操作点
        plt.annotate(f'{section_name}操作点\nLs={Ls_design:.6f} m³/s\nVs={Vs_design:.4f} m³/s',
                     xy=(Ls_design, Vs_design),
                     xytext=(Ls_design * 1.1, Vs_design * 1.1),
                     arrowprops=dict(arrowstyle='->', color='black'))

        # 填充稳定操作区（简化的近似方法，避免复杂的插值）
        try:
            # 使用简单的多边形填充方法
            # 收集所有边界点
            all_points = []

            # 液相负荷下限线
            all_points.append((Ls_lower[0], 0))
            all_points.append((Ls_lower[0], Vs_upper[0] * 0.1))  # 接近0的小值

            # 液相负荷上限线
            all_points.append((Ls_upper[0], 0))
            all_points.append((Ls_upper[0], Vs_upper[0] * 0.1))

            # 收集有效的限制线点
            if len(Ls_ent) > 1:
                for i in range(len(Ls_ent)):
                    if Ls_lower[0] <= Ls_ent[i] <= Ls_upper[0]:
                        all_points.append((Ls_ent[i], Vs_ent[i]))

            if len(Ls_flood) > 1:
                for i in range(len(Ls_flood)):
                    if Ls_lower[0] <= Ls_flood[i] <= Ls_upper[0]:
                        all_points.append((Ls_flood[i], Vs_flood[i]))

            if len(Ls_leak) > 1:
                for i in range(len(Ls_leak)):
                    if Ls_lower[0] <= Ls_leak[i] <= Ls_upper[0]:
                        all_points.append((Ls_leak[i], Vs_leak[i]))

            # 按角度排序点
            if len(all_points) >= 3:
                # 计算质心
                cx = sum(p[0] for p in all_points) / len(all_points)
                cy = sum(p[1] for p in all_points) / len(all_points)

                # 按极角排序
                def angle(p):
                    return math.atan2(p[1] - cy, p[0] - cx)

                all_points.sort(key=angle)

                # 提取坐标
                poly_Ls = [p[0] for p in all_points]
                poly_Vs = [p[1] for p in all_points]

                # 填充多边形
                plt.fill(poly_Ls, poly_Vs, alpha=0.2, color='gray', label='稳定操作区')
        except Exception as e:
            print(f"  警告: 无法计算稳定操作区: {e}")
            # 使用简化的矩形区域作为近似
            Ls_fill_min = Ls_lower[0]
            Ls_fill_max = Ls_upper[0]
            if Ls_fill_max > Ls_fill_min:
                # 计算近似的高度
                if len(Vs_ent) > 0 and len(Vs_leak) > 0:
                    Vs_fill_min = min(Vs_leak) * 0.8
                    Vs_fill_max = min(min(Vs_ent), min(Vs_flood)) * 0.8 if len(Vs_ent) > 0 and len(
                        Vs_flood) > 0 else 0.1
                    if Vs_fill_max > Vs_fill_min:
                        plt.fill_between([Ls_fill_min, Ls_fill_max],
                                         [Vs_fill_min, Vs_fill_min],
                                         [Vs_fill_max, Vs_fill_max],
                                         alpha=0.2, color='gray', label='稳定操作区(近似)')

        # 设置图形属性
        plt.xlabel('液相流量 Ls (m³/s)', fontsize=12)
        plt.ylabel('气相流量 Vs (m³/s)', fontsize=12)
        plt.title(f'{section_name}负荷性能图（浮阀塔）', fontsize=14, fontweight='bold')
        plt.grid(True, alpha=0.3)
        plt.legend(loc='best', fontsize=10)
        plt.tight_layout()

        # 设置坐标轴范围
        plt.xlim(0, max(Ls_upper) * 1.2)
        plt.ylim(0, max(Vs_ent) * 1.2 if len(Vs_ent) > 0 else Vs_design * 2)

        # 保存图形
        filename = f'float_valve_tower_load_performance_{section_name}.png'
        plt.savefig(filename, dpi=300)
        print(f"  ✓ {section_name}负荷性能图已保存为 '{filename}'")

        plt.show()

    def calculate_operational_elasticity_for_section(self, Ls_ent, Vs_ent, Ls_flood, Vs_flood,
                                                     Ls_upper, Vs_upper, Ls_leak, Vs_leak,
                                                     Ls_lower, Vs_lower, Ls_design, Vs_design,
                                                     section_name="精馏段"):
        """计算指定塔段的操作弹性"""
        print(f"\n{section_name}操作弹性计算")

        try:
            # 检查是否有足够的数据点
            if len(Ls_ent) < 2 or len(Ls_flood) < 2 or len(Ls_leak) < 2:
                print(f"  警告: {section_name}限制线数据点不足，使用简化计算方法")

                # 使用简化的方法：直接在设计点附近估算
                if section_name == "精馏段":
                    Vs_design_val = self.rectifying['V_vol_m3s']
                else:
                    Vs_design_val = self.stripping['V_vol_m3s']

                # 简化的估算：假设最大气速为设计值的1.5倍，最小气速为设计值的0.5倍
                Vs_max = Vs_design_val * 1.5
                Vs_min = Vs_design_val * 0.5
                operational_elasticity = Vs_max / Vs_min if Vs_min > 0 else 0

                print(f"  使用简化计算方法:")
                print(f"  最大操作气速 Vs,max = {Vs_max:.6f} m³/s (估算值)")
                print(f"  最小操作气速 Vs,min = {Vs_min:.6f} m³/s (估算值)")
                print(f"  操作弹性 = {operational_elasticity:.2f}")

                # 确保check_results字典存在
                if section_name not in self.check_results:
                    self.check_results[section_name] = {}

                # 存储结果
                self.check_results[section_name]['operational_elasticity'] = {
                    'Vs_max_m3s': float(Vs_max),
                    'Vs_min_m3s': float(Vs_min),
                    'operational_elasticity': float(operational_elasticity),
                    'elasticity_rating': "估算值"
                }

                return operational_elasticity

            # 在操作点液相流量处，计算各限制线对应的气相流量
            # 使用线性插值，但先检查数据是否有效
            def interpolate_at_Ls(Ls_values, Vs_values, target_Ls):
                """在target_Ls处线性插值Vs值"""
                if len(Ls_values) < 2:
                    return None

                # 找到最接近的两个点
                idx = np.searchsorted(Ls_values, target_Ls)
                if idx == 0:
                    return Vs_values[0]
                elif idx == len(Ls_values):
                    return Vs_values[-1]
                else:
                    # 线性插值
                    x0, x1 = Ls_values[idx - 1], Ls_values[idx]
                    y0, y1 = Vs_values[idx - 1], Vs_values[idx]
                    return y0 + (y1 - y0) * (target_Ls - x0) / (x1 - x0)

            # 液相负荷上下限
            Ls_upper_limit = Ls_upper[0]
            Ls_lower_limit = Ls_lower[0]

            print(f"\n   在{section_name}操作点液相流量 Ls = {Ls_design:.6f} m³/s 处:")

            # 确保Ls_design在限制线范围内
            if Ls_design < min(Ls_ent) or Ls_design > max(Ls_ent):
                print(f"   警告: 操作点液相流量不在液沫夹带线范围内")
                Vs_ent_at_op = None
            else:
                Vs_ent_at_op = interpolate_at_Ls(Ls_ent, Vs_ent, Ls_design)
                if Vs_ent_at_op is not None:
                    print(f"   过量液沫夹带限制: Vs ≤ {Vs_ent_at_op:.6f} m³/s")

            if Ls_design < min(Ls_flood) or Ls_design > max(Ls_flood):
                print(f"   警告: 操作点液相流量不在液泛线范围内")
                Vs_flood_at_op = None
            else:
                Vs_flood_at_op = interpolate_at_Ls(Ls_flood, Vs_flood, Ls_design)
                if Vs_flood_at_op is not None:
                    print(f"   溢流液泛限制: Vs ≤ {Vs_flood_at_op:.6f} m³/s")

            if Ls_design < min(Ls_leak) or Ls_design > max(Ls_leak):
                print(f"   警告: 操作点液相流量不在漏液线范围内")
                Vs_leak_at_op = None
            else:
                Vs_leak_at_op = interpolate_at_Ls(Ls_leak, Vs_leak, Ls_design)
                if Vs_leak_at_op is not None:
                    print(f"   漏液限制: Vs ≥ {Vs_leak_at_op:.6f} m³/s")

            print(f"   液相负荷上限: Ls ≤ {Ls_upper_limit:.6f} m³/s")
            print(f"   液相负荷下限: Ls ≥ {Ls_lower_limit:.6f} m³/s")

            # 检查操作点是否在液相负荷范围内
            if Ls_design < Ls_lower_limit or Ls_design > Ls_upper_limit:
                print(f"   ⚠️ 警告: {section_name}操作点不在液相负荷范围内!")

            # 确定最大和最小操作气速
            Vs_max_values = []
            if Vs_ent_at_op is not None:
                Vs_max_values.append(Vs_ent_at_op)
            if Vs_flood_at_op is not None:
                Vs_max_values.append(Vs_flood_at_op)

            Vs_min_values = []
            if Vs_leak_at_op is not None:
                Vs_min_values.append(Vs_leak_at_op)

            if not Vs_max_values:
                print(f"   警告: 无法确定最大操作气速，使用设计值的1.5倍作为估算")
                Vs_max = Vs_design * 1.5
            else:
                Vs_max = min(Vs_max_values)

            if not Vs_min_values:
                print(f"   警告: 无法确定最小操作气速，使用设计值的0.5倍作为估算")
                Vs_min = Vs_design * 0.5
            else:
                Vs_min = max(Vs_min_values)

            # 操作弹性 = Vs,max / Vs,min
            if Vs_min > 0:
                operational_elasticity = Vs_max / Vs_min
                print(f"\n   {section_name}最大操作气速 Vs,max = {Vs_max:.6f} m³/s")
                print(f"   {section_name}最小操作气速 Vs,min = {Vs_min:.6f} m³/s")
                print(
                    f"   {section_name}操作弹性 = Vs,max / Vs,min = {Vs_max:.6f} / {Vs_min:.6f} = {operational_elasticity:.2f}")

                # 评价操作弹性
                if operational_elasticity > 3:
                    elasticity_rating = "优秀"
                elif operational_elasticity > 2:
                    elasticity_rating = "良好"
                elif operational_elasticity > 1.5:
                    elasticity_rating = "合格"
                else:
                    elasticity_rating = "不足"

                print(f"   {section_name}操作弹性评价: {elasticity_rating}")

                # 确保check_results字典存在
                if section_name not in self.check_results:
                    self.check_results[section_name] = {}

                # 存储操作弹性结果
                self.check_results[section_name]['operational_elasticity'] = {
                    'Vs_max_m3s': float(Vs_max),
                    'Vs_min_m3s': float(Vs_min),
                    'operational_elasticity': float(operational_elasticity),
                    'elasticity_rating': elasticity_rating
                }

                return operational_elasticity
            else:
                print("   ✗ 无法计算操作弹性（Vs_min ≤ 0）")
                return 0

        except Exception as e:
            print(f"   ✗ 计算操作弹性时出错: {str(e)}")
            print(f"  使用简化计算方法")

            # 使用简化的方法
            if section_name == "精馏段":
                Vs_design_val = self.rectifying['V_vol_m3s']
            else:
                Vs_design_val = self.stripping['V_vol_m3s']

            # 简化的估算
            Vs_max = Vs_design_val * 1.5
            Vs_min = Vs_design_val * 0.5
            operational_elasticity = Vs_max / Vs_min if Vs_min > 0 else 0

            # 确保check_results字典存在
            if section_name not in self.check_results:
                self.check_results[section_name] = {}

            self.check_results[section_name]['operational_elasticity'] = {
                'Vs_max_m3s': float(Vs_max),
                'Vs_min_m3s': float(Vs_min),
                'operational_elasticity': float(operational_elasticity),
                'elasticity_rating': "估算值(出错时备用)"
            }

            return operational_elasticity

    def run_section_check(self, section_name="精馏段"):
        """运行指定塔段的流体力学验算"""
        print(f"\n{'=' * 80}")
        print(f"{section_name}流体力学验算（浮阀塔）")
        print(f"{'=' * 80}")

        # 获取塔段参数
        if section_name == "精馏段":
            ρ_v = self.rectifying['ρ_v']
            ρ_L = self.rectifying['ρ_L']
            σ = self.rectifying['σ']
            HT = self.rectifying['HT']
            V_m3s = self.rectifying['V_vol_m3s']
            L_m3s = self.rectifying['L_vol_m3s']
        else:
            ρ_v = self.stripping['ρ_v']
            ρ_L = self.stripping['ρ_L']
            σ = self.stripping['σ']
            HT = self.stripping['HT']
            V_m3s = self.stripping['V_vol_m3s']
            L_m3s = self.stripping['L_vol_m3s']

        print(f"\n{section_name}工艺参数:")
        print(f"  气相流量: {V_m3s:.6f} m³/s")
        print(f"  液相流量: {L_m3s:.6f} m³/s")
        print(f"  气相密度: {ρ_v:.4f} kg/m³")
        print(f"  液相密度: {ρ_L:.1f} kg/m³")
        print(f"  表面张力: {σ * 1000:.2f} mN/m")
        print(f"  板间距: {HT:.3f} m")

        # 1. 塔板压降验算
        print(f"\n一、{section_name}塔板压降验算")
        hp, hc, hl, hσ = self.calculate_plate_pressure_drop(V_m3s, ρ_v, ρ_L, σ, HT, section_name)
        hp_pa = hp * ρ_L * self.g  # 转换为Pa

        print(f"\n   {section_name}塔板压降:")
        print(f"   干板压降 hc: {hc * 1000:.1f} mm液柱")
        print(f"   液层压降 hl: {hl * 1000:.1f} mm液柱")
        print(f"   表面张力压降 hσ: {hσ * 1000:.2f} mm液柱")
        print(f"   总压降 hp: {hp * 1000:.1f} mm液柱 = {hp_pa:.0f} Pa")

        # 压降评价
        if hp_pa < 400:
            hp_rating = "偏低"
        elif hp_pa <= 800:
            hp_rating = "合理"
        else:
            hp_rating = "偏高"
        print(f"   评价: {hp_rating}")

        # 2. 液沫夹带验算
        print(f"\n二、{section_name}液沫夹带验算")
        ev, ev_check = self.calculate_entrainment(V_m3s, L_m3s, ρ_v, ρ_L, σ, HT, section_name)
        print(f"   液沫夹带量: {ev:.6f} kg液体/kg气体")
        print(f"   限值: < 0.1")
        ev_status = "✓ 合格" if ev_check else "✗ 不合格"
        print(f"   结果: {ev_status}")

        # 3. 漏液和稳定性验算
        print(f"\n三、{section_name}漏液和稳定性验算")
        # 计算阀孔气速
        u0 = V_m3s / self.Av_total
        # 计算漏液点气速（使用缓存版本）
        u_ow = self.calculate_leakage_point_velocity(L_m3s, ρ_v, ρ_L, σ, section_name)
        # 稳定性系数
        K = u0 / u_ow if u_ow > 0 else float('inf')
        print(f"   设计阀孔气速 u0: {u0:.4f} m/s")
        print(f"   漏液点气速 u_ow: {u_ow:.4f} m/s")
        print(f"   稳定性系数 K: {K:.2f}")
        print(f"   要求: > 1.5")
        K_check = K > 1.5
        K_status = "✓ 合格" if K_check else "✗ 不合格"
        print(f"   结果: {K_status}")

        # 4. 降液管液泛验算
        print(f"\n四、{section_name}降液管液泛验算")
        Hd, Hd_max, flooding_check = self.calculate_downcomer_flooding(V_m3s, L_m3s, ρ_v, ρ_L, σ, HT, section_name)
        print(f"   降液管液层高度 Hd: {Hd:.4f} m")
        print(f"   允许最大高度 Hd_max: {Hd_max:.3f} m")
        flooding_status = "✓ 合格" if flooding_check else "✗ 不合格"
        print(f"   结果: {flooding_status}")

        # 5. 负荷性能图
        print(f"\n五、{section_name}负荷性能图")
        try:
            operational_elasticity = self.generate_load_performance_diagram_for_section(section_name)
        except Exception as e:
            print(f"   警告: 生成负荷性能图时出错: {str(e)}")
            operational_elasticity = 0

        # 确保check_results字典存在
        if section_name not in self.check_results:
            self.check_results[section_name] = {}

        # 存储验算结果
        self.check_results[section_name]['summary'] = {
            'plate_pressure_drop': {
                'hp_m': float(hp),
                'hp_mm': float(hp * 1000),
                'hp_pa': float(hp_pa),
                'rating': hp_rating
            },
            'entrainment': {
                'ev_kgkg': float(ev),
                'check': bool(ev_check),
                'limit': 0.1
            },
            'stability': {
                'u0_ms': float(u0),
                'u_ow_ms': float(u_ow),
                'K': float(K),
                'check': bool(K_check),
                'limit': 1.5
            },
            'flooding': {
                'Hd_m': float(Hd),
                'Hd_max_m': float(Hd_max),
                'check': bool(flooding_check)
            },
            'overall_check': bool(ev_check and K_check and flooding_check)
        }

        return self.check_results[section_name]

    def run_complete_check(self):
        """运行完整的浮阀塔流体力学验算（包括精馏段和提馏段）"""
        print("=" * 80)
        print("浮阀塔流体力学验算程序")
        print("=" * 80)

        try:
            # 分别验算精馏段和提馏段
            print("\n一、精馏段验算")
            rect_results = self.run_section_check("精馏段")

            print("\n\n二、提馏段验算")
            strip_results = self.run_section_check("提馏段")

            # 比较两个塔段的操作弹性
            print("\n" + "=" * 80)
            print("精馏段与提馏段对比分析")
            print("=" * 80)

            rect_elasticity = self.check_results['精馏段'].get('operational_elasticity', {}).get(
                'operational_elasticity', 0)
            strip_elasticity = self.check_results['提馏段'].get('operational_elasticity', {}).get(
                'operational_elasticity', 0)

            print(f"\n操作弹性对比:")
            print(f"  精馏段操作弹性: {rect_elasticity:.2f}")
            print(f"  提馏段操作弹性: {strip_elasticity:.2f}")

            # 确定控制塔段（操作弹性较小的塔段）
            if rect_elasticity > 0 and strip_elasticity > 0:
                control_section = "精馏段" if rect_elasticity < strip_elasticity else "提馏段"
                print(f"\n控制塔段: {control_section} (操作弹性较小)")
                print(f"设计依据: 以{control_section}的操作弹性为设计基准")
            else:
                print(f"\n无法确定控制塔段，请检查计算结果")

            # 总体评价
            print(f"\n总体评价:")
            rect_ok = self.check_results['精馏段']['summary']['overall_check']
            strip_ok = self.check_results['提馏段']['summary']['overall_check']

            if rect_ok and strip_ok:
                print("  ✅ 精馏段和提馏段所有验算项目均合格，设计合理")
            else:
                print("  ⚠️ 部分验算项目不合格，需要调整设计")
                if not rect_ok:
                    print("  精馏段存在问题:")
                    summary = self.check_results['精馏段']['summary']
                    if not summary['entrainment']['check']:
                        print("    - 液沫夹带超标")
                    if not summary['stability']['check']:
                        print("    - 稳定性不足")
                    if not summary['flooding']['check']:
                        print("    - 降液管液泛")
                if not strip_ok:
                    print("  提馏段存在问题:")
                    summary = self.check_results['提馏段']['summary']
                    if not summary['entrainment']['check']:
                        print("    - 液沫夹带超标")
                    if not summary['stability']['check']:
                        print("    - 稳定性不足")
                    if not summary['flooding']['check']:
                        print("    - 降液管液泛")

            # 保存验算结果
            self.save_check_results()

            return self.check_results

        except Exception as e:
            print(f"✗ 程序运行出错: {str(e)}")
            import traceback
            traceback.print_exc()
            return None

    def save_check_results(self, filename='float_valve_tower_hydrodynamic_check_complete.json'):
        """保存完整的验算结果到JSON文件"""
        results = {
            "设计信息": {
                "项目": "浮阀塔流体力学验算（精馏段+提馏段）",
                "塔径_m": self.D,
                "控制塔段": self.control_section,
                "验算日期": "2025-12-28",
                "设计依据文件": [
                    "tower_plate_arrangement_optimized.json",
                    "tower_diameter_approximate_optimized.json",
                    "tower_properties_results.json"
                ]
            },
            "塔板结构参数": {
                "阀孔数": self.N_valves,
                "阀孔直径_m": self.dv,
                "阀孔直径_mm": self.dv * 1000,
                "开孔率_%": self.φ * 100,
                "塔板厚度_m": self.δ,
                "浮阀类型": self.valve_type,
                "阀重_kg": self.valve_weight,
                "总阀孔面积_m2": self.Av_total,
                "鼓泡区面积_m2": self.A_a,
                "堰长_m": self.lw,
                "堰高_m": self.hw,
                "降液管底隙高度_m": self.h0,
                "降液管面积_流通_m2": self.A_f,
                "降液管面积_投影_m2": self.A_f_projected,
                "精馏段板间距_m": self.HT_rectifying,
                "提馏段板间距_m": self.HT_stripping
            },
            "工艺参数": {
                "精馏段": {
                    "气相流量_m3s": self.rectifying['V_vol_m3s'],
                    "液相流量_m3s": self.rectifying['L_vol_m3s'],
                    "气相密度_kgm3": self.rectifying['ρ_v'],
                    "液相密度_kgm3": self.rectifying['ρ_L'],
                    "表面张力_Nm": self.rectifying['σ']
                },
                "提馏段": {
                    "气相流量_m3s": self.stripping['V_vol_m3s'],
                    "液相流量_m3s": self.stripping['L_vol_m3s'],
                    "气相密度_kgm3": self.stripping['ρ_v'],
                    "液相密度_kgm3": self.stripping['ρ_L'],
                    "表面张力_Nm": self.stripping['σ']
                }
            },
            "验算结果": self.check_results,
            "负荷性能数据": self.load_performance_data
        }

        try:
            with open(filename, 'w', encoding='utf-8') as f:
                json.dump(results, f, ensure_ascii=False, indent=4)

            print(f"\n✓ 完整验算结果已保存到 '{filename}'")
            print(f"✓ 精馏段负荷性能图已保存为 'float_valve_tower_load_performance_精馏段.png'")
            print(f"✓ 提馏段负荷性能图已保存为 'float_valve_tower_load_performance_提馏段.png'")

        except Exception as e:
            print(f"✗ 保存验算结果时出错: {e}")


def main():
    """主函数"""
    print("浮阀塔流体力学验算程序（精馏段+提馏段）")
    print("=" * 80)
    print("功能说明:")
    print("  1. 分别计算精馏段和提馏段的塔板压降")
    print("  2. 分别计算精馏段和提馏段的液沫夹带")
    print("  3. 分别计算精馏段和提馏段的漏液点和稳定性")
    print("  4. 分别计算精馏段和提馏段的降液管液泛")
    print("  5. 分别绘制精馏段和提馏段的负荷性能图")
    print("  6. 比较两个塔段的操作弹性，确定控制塔段")
    print("=" * 80)

    try:
        # 创建验算器
        checker = FloatValveTowerHydrodynamicCheck(
            plate_arrangement_file='tower_plate_arrangement_optimized.json',
            diameter_design_file='tower_diameter_approximate_optimized.json',
            properties_file='tower_properties_results.json'
        )

        # 运行完整验算
        results = checker.run_complete_check()

        print("\n" + "=" * 80)
        print("浮阀塔流体力学验算完成！")
        print("=" * 80)

        return checker

    except Exception as e:
        print(f"✗ 程序运行出错: {str(e)}")
        import traceback
        traceback.print_exc()
        return None


if __name__ == "__main__":
    checker = main()