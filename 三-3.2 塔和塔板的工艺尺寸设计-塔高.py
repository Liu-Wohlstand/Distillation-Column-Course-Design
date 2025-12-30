"""
塔高设计程序
根据塔板数、板间距、人孔布置等计算塔体总有效高度
"""

import math
import json


class TowerHeightCalculator:
    """塔高设计计算器"""

    def __init__(self, tower_data_file='tower_properties_results.json',
                 calc_data_file='benzene_toluene_results.json', interactive_mode=False):
        """
        初始化塔高计算器

        参数:
        tower_data_file: 塔工艺特性结果文件
        calc_data_file: 计算过程结果文件
        interactive_mode: 是否使用交互模式输入参数
        """
        # 加载塔工艺条件数据
        try:
            with open(tower_data_file, 'r', encoding='utf-8') as f:
                self.tower_props = json.load(f)
        except FileNotFoundError:
            print(f"警告: 未找到塔工艺特性文件 '{tower_data_file}'，将使用默认值。")
            self.tower_props = {}

        # 加载计算过程数据
        try:
            with open(calc_data_file, 'r', encoding='utf-8') as f:
                self.calc_data = json.load(f)
        except FileNotFoundError:
            print(f"警告: 未找到计算过程文件 '{calc_data_file}'，将使用默认值。")
            self.calc_data = {}

        # 设计参数配置
        self.design_config = {}

        if interactive_mode:
            self.get_design_parameters_interactive()
        else:
            self.set_default_design_parameters()

        # 从计算结果获取实际塔板数
        self.N1_actual = math.ceil(self.calc_data.get('N1_actual', 10))  # 默认10块
        self.N2_actual = math.ceil(self.calc_data.get('N2_actual', 8))  # 默认8块
        self.N_total = self.N1_actual + self.N2_actual  # 总实际塔板数

        # 进料板位置（假设在精馏段最后一块板）
        self.feed_stage = self.N1_actual

        print("=" * 70)
        print("塔高设计计算程序")
        print("=" * 70)
        print(f"总实际塔板数: {self.N_total} 块")
        print(f"  精馏段: {self.N1_actual} 块")
        print(f"  提馏段: {self.N2_actual} 块")
        print(f"进料板位置: 第 {self.feed_stage} 块塔板")
        print()

    def get_design_parameters_interactive(self):
        """通过交互方式获取设计参数"""
        print("=" * 70)
        print("塔高设计参数输入（按Enter使用默认值）")
        print("=" * 70)

        # 塔顶空间高度
        print("\n1. 塔顶空间高度（HD）")
        print("-" * 40)
        print("设计要点:")
        print("  • 保证塔顶气液分离，减少雾沫夹带")
        print("  • 通常为1.0-2.0m，塔径大时取大值")
        print("  • 常压塔: 1.0-1.5m，减压塔: 1.5-2.0m")
        hd_min = float(1.0)
        hd_max = float(1.5)
        hd_design = float(input(f"设计塔顶空间高度(m) [建议: 1.2-1.8, 默认{1.2}]: ") or 1.2)

        # 板间距
        print("\n2. 板间距（HT）")
        print("-" * 40)
        print("设计要点:")
        print("  • 影响雾沫夹带、操作弹性、塔高")
        print("  • 常压塔: 0.3-0.6m，减压塔: 0.45-0.8m")
        print("  • 大塔径、高液气比时取大值")
        ht_rect = float(input(f"精馏段板间距(m) [参考: 0.3-0.6, 默认{0.45}]: ") or 0.45)
        ht_strip = float(input(f"提馏段板间距(m) [参考: 0.3-0.6, 默认{0.45}]: ") or 0.45)

        # 塔釜设计
        print("\n3. 塔釜设计参数")
        print("-" * 40)
        print("设计要点:")
        print("  • h1: 保证底层塔板操作稳定，防止液泛")
        print("  • 装填系数: 防止液泛和波动，一般0.6-0.7")
        print("  • 塔釜高度: 保证足够停留时间（通常5-15分钟）")
        h1_min = float( 0.5)
        h1_max = float(0.7)
        h1_design = float(input(f"设计釜液面距底层塔板距离(m) [建议: 0.6-0.8, 默认{0.6}]: ") or 0.6)

        fill_factor = float(input(f"装填系数(0.6-0.7) [参考: 真空塔0.5-0.6, 常压塔0.6-0.7, 默认{0.65}]: ") or 0.65)

        print("\n塔釜高度要求:")
        print("  • 小塔（D<1m）: 1.0-1.5m")
        print("  • 中塔（D=1-2m）: 1.5-2.0m")
        print("  • 大塔（D>2m）: 2.0-2.5m")
        hb_min = float( 1.5)
        hb_max = float( 2.0)

        # 人孔设计
        print("\n4. 人孔设计参数")
        print("-" * 40)
        print("设计要点:")
        print("  • 人孔间距: 一般8-12块板设置一个")
        print("  • 人孔直径: D<1m时用400mm, D=1-2m时用450mm, D>2m时用500mm")
        print("  • 人孔处板间距: 一般比正常板间距大150-200mm")

        manhole_int_min = int( 8)
        manhole_int_max = int(10)
        manhole_interval = int(input(f"设计人孔间距(板数) [建议: 8-10, 默认{9}]: ") or 9)

        print("\n人孔直径选择:")
        print("  • D<1.0m: 0.4m (Φ400mm)")
        print("  • D=1.0-1.8m: 0.45m (Φ450mm)")
        print("  • D>1.8m: 0.5m (Φ500mm)")
        manhole_dia = float(input(f"人孔直径(m) [参考: 0.4-0.5, 默认{0.5}]: ") or 0.5)

        manhole_ht_rect = float(input(f"精馏段人孔处板间距(m) [比正常板间距大0.15-0.2m, 默认{0.6}]: ") or 0.6)
        manhole_ht_strip = float(input(f"提馏段人孔处板间距(m) [比正常板间距大0.15-0.2m, 默认{0.6}]: ") or 0.6)

        # 进料板处空间高度
        print("\n5. 进料板处空间高度")
        print("-" * 40)
        print("设计要点:")
        print("  • 保证进料分布和两相分离")
        print("  • 一般比正常板间距大20%-50%")
        hf_factor = float(input(f"进料板处空间高度系数 [参考: 1.2-1.5, 默认{1.2}]: ") or 1.2)

        # 安全系数
        print("\n6. 安全系数")
        print("-" * 40)
        print("设计要点:")
        print("  • 考虑安装误差、热膨胀等")
        print("  • 一般取5%-15%")
        safety_margin = float(input(f"安全余量(如0.1表示10%) [参考: 0.05-0.15, 默认{0.1}]: ") or 0.1)

        # 塔板数（如果数据文件未找到）
        print("\n7. 塔板数（如果数据文件未找到）")
        print("-" * 40)
        print("设计要点:")
        print("  • 精馏段: 根据分离要求计算")
        print("  • 提馏段: 根据分离要求计算")
        if 'N1_actual' not in self.calc_data:
            n1 = int(input(f"精馏段塔板数 [参考: 5-30块, 默认{10}]: ") or 10)
            n2 = int(input(f"提馏段塔板数 [参考: 5-30块, 默认{8}]: ") or 8)
            self.calc_data['N1_actual'] = n1
            self.calc_data['N2_actual'] = n2

        # 验证输入参数
        self.validate_parameters(hd_design, hd_min, hd_max, "塔顶空间高度")
        self.validate_parameters(h1_design, h1_min, h1_max, "釜液面距底层塔板距离")
        self.validate_parameters(fill_factor, 0.5, 0.7, "装填系数")
        self.validate_parameters(hf_factor, 1.1, 1.6, "进料板处空间系数")
        self.validate_parameters(safety_margin, 0.05, 0.15, "安全余量")

        # 保存到设计配置
        self.design_config = {
            'HD_min': hd_min,
            'HD_max': hd_max,
            'HD_design': hd_design,

            'HT_rectifying': ht_rect,
            'HT_stripping': ht_strip,

            'h1_min': h1_min,
            'h1_max': h1_max,
            'h1_design': h1_design,
            'fill_factor': fill_factor,
            'HB_min': hb_min,
            'HB_max': hb_max,

            'manhole_interval_min': manhole_int_min,
            'manhole_interval_max': manhole_int_max,
            'manhole_interval': manhole_interval,
            'manhole_diameter': manhole_dia,
            'manhole_HT_rectifying': manhole_ht_rect,
            'manhole_HT_stripping': manhole_ht_strip,

            'HF_factor': hf_factor,

            'safety_margin': safety_margin
        }

        print("\n" + "=" * 70)
        print("设计参数输入完成!")
        print("=" * 70)

    def validate_parameters(self, value, min_val, max_val, param_name):
        """验证参数是否在合理范围内"""
        if value < min_val or value > max_val:
            print(f"⚠  警告: {param_name}={value} 超出建议范围({min_val}-{max_val})")
            confirm = input(f"  是否继续使用此值? (y/n) [默认y]: ").lower() or 'y'
            if confirm != 'y':
                raise ValueError(f"{param_name} 超出合理范围")

    def set_default_design_parameters(self):
        """设置默认设计参数"""
        self.design_config = {
            # 塔顶空间高度 (m)
            'HD_min': 1.0,  # 最小，参考：常压塔1.0-1.5m
            'HD_max': 1.5,  # 最大，参考：常压塔1.5-2.0m
            'HD_design': 1.2,  # 设计值，参考：常压塔推荐1.2-1.5m

            # 板间距 (m)
            'HT_rectifying': 0.45,  # 精馏段板间距，参考：常压塔0.3-0.6m
            'HT_stripping': 0.45,  # 提馏段板间距，参考：常压塔0.3-0.6m

            # 塔釜设计
            'h1_min': 0.5,  # 最小釜液面距底层塔板距离，参考：0.5-1.0m
            'h1_max': 0.7,  # 最大，参考：0.7-1.2m
            'h1_design': 0.6,  # 设计值，参考：推荐0.6-0.8m

            'fill_factor': 0.65,  # 装填系数，一般0.6-0.7，真空塔0.5-0.6
            'HB_min': 1.5,  # 最小塔釜高度，参考：中塔1.5-2.0m
            'HB_max': 2.0,  # 最大塔釜高度，参考：中塔1.5-2.5m

            # 人孔设计
            'manhole_interval_min': 8,  # 最小人孔间距 (板数)，参考：6-10
            'manhole_interval_max': 10,  # 最大人孔间距，参考：8-12
            'manhole_interval': 9,  # 设计人孔间距，参考：8-10

            'manhole_diameter': 0.5,  # 人孔直径 (m)，参考：Φ400-500mm

            'manhole_HT_rectifying': 0.6,  # 精馏段人孔处板间距，参考：比正常大0.15-0.2m
            'manhole_HT_stripping': 0.6,  # 提馏段人孔处板间距，参考：比正常大0.15-0.2m

            # 进料板处空间高度
            'HF_factor': 1.2,  # 进料板处空间高度系数（相对于正常板间距），参考：1.2-1.5

            # 安全系数
            'safety_margin': 0.1,  # 10%的安全余量，参考：5%-15%
        }

    def calculate_tower_height(self):
        """计算塔体总有效高度"""
        print("一、塔体总有效高度计算")
        print("=" * 70)

        # 1. 塔顶空间高度
        HD = self.design_config['HD_design']
        print(f"\n1. 塔顶空间高度 HD")
        print("-" * 40)
        print(f"取值范围: {self.design_config['HD_min']:.1f} ~ {self.design_config['HD_max']:.1f} m")
        print(f"设计取值: HD = {HD:.2f} m")
        print(f"设计参考:")
        print(f"  • 常压塔: 1.0-1.5m")
        print(f"  • 减压塔: 1.5-2.0m")
        print(f"  • 作用: 保证气液分离，减少雾沫夹带")

        # 2. 塔釜高度计算
        print(f"\n2. 塔釜高度 HB")
        print("-" * 40)

        # 釜液面距底层塔板距离
        h1 = self.design_config['h1_design']
        print(f"(1) 釜液面距底层塔板距离 h1")
        print(f"    取值范围: {self.design_config['h1_min']:.1f} ~ {self.design_config['h1_max']:.1f} m")
        print(f"    设计取值: h1 = {h1:.2f} m")
        print(f"    设计参考:")
        print(f"      • 范围: 0.5-1.2m")
        print(f"      • 推荐: 0.6-0.8m")
        print(f"      • 作用: 保证底层塔板稳定操作，防止液泛")

        # 计算釜内液层高度 h2
        φ = self.design_config['fill_factor']
        # φ = h2/(h1+h2) => h2 = φ*h1/(1-φ)
        h2 = φ * h1 / (1 - φ)
        print(f"\n(2) 釜内液层高度 h2")
        print(f"    装填系数 φ = {φ:.2f}")
        print(f"    设计参考:")
        print(f"      • 常压塔: 0.6-0.7")
        print(f"      • 减压塔: 0.5-0.6")
        print(f"      • 作用: 防止液泛和液面波动")
        print(f"    计算公式: h2 = φ × h1 / (1 - φ)")
        print(f"             = {φ:.2f} × {h1:.2f} / (1 - {φ:.2f})")
        print(f"             = {h2:.2f} m")

        # 塔釜总高度
        HB = h1 + h2
        print(f"\n(3) 塔釜总高度 HB")
        print(f"    HB = h1 + h2 = {h1:.2f} + {h2:.2f} = {HB:.2f} m")
        print(f"    要求: HB > {self.design_config['HB_min']:.1f} ~ {self.design_config['HB_max']:.1f} m")
        print(f"    设计结果: {'满足' if HB >= self.design_config['HB_min'] else '不满足'}")
        print(f"    设计参考:")
        print(f"      • 停留时间: 通常5-15分钟")
        print(f"      • 小塔(D<1m): 1.0-1.5m")
        print(f"      • 中塔(D=1-2m): 1.5-2.0m")
        print(f"      • 大塔(D>2m): 2.0-2.5m")

        # 3. 人孔设计与布置
        print(f"\n3. 人孔设计与布置")
        print("-" * 40)

        manhole_interval = self.design_config['manhole_interval']
        print(f"人孔设置原则: 每隔 {manhole_interval} 块塔板设一个人孔")
        print(
            f"取值范围: 每隔 {self.design_config['manhole_interval_min']} ~ {self.design_config['manhole_interval_max']} 块塔板设一个人孔")
        print(f"设计参考:")
        print(f"  • 一般间距: 8-12块板")
        print(f"  • 作用: 安装、检修、清洗")

        # 计算总人孔数
        S_total = self.N_total // manhole_interval
        print(f"\n总人孔数 S = N_total // {manhole_interval}")
        print(f"          = {self.N_total} // {manhole_interval}")
        print(f"          = {S_total} 个")

        # 确定人孔具体位置
        manhole_positions = []
        for i in range(1, S_total + 1):
            pos = i * manhole_interval
            if pos <= self.N_total:
                manhole_positions.append(pos)

        print(f"人孔位置（从上到下）: {manhole_positions}")

        # 统计精馏段和提馏段的人孔数
        S_rect = sum(1 for pos in manhole_positions if pos <= self.N1_actual)
        S_strip = S_total - S_rect

        print(f"精馏段人孔数: {S_rect} 个")
        print(f"提馏段人孔数: {S_strip} 个")

        # 人孔处板间距
        HT_manhole_rect = self.design_config['manhole_HT_rectifying']
        HT_manhole_strip = self.design_config['manhole_HT_stripping']
        print(f"精馏段人孔处板间距: {HT_manhole_rect:.2f} m")
        print(f"提馏段人孔处板间距: {HT_manhole_strip:.2f} m")
        print(f"人孔直径: Φ{self.design_config['manhole_diameter'] * 1000:.0f} mm")
        print(f"设计参考:")
        print(f"  • 人孔直径选择:")
        print(f"      D<1.0m: Φ400mm")
        print(f"      D=1.0-1.8m: Φ450mm")
        print(f"      D>1.8m: Φ500mm")
        print(f"  • 板间距增加: 比正常板间距大150-200mm")

        # 4. 进料板处空间高度
        print(f"\n4. 进料板处空间高度")
        print("-" * 40)

        # 进料板处空间高度通常比正常板间距大
        HF_factor = self.design_config['HF_factor']
        HF = HF_factor * self.design_config['HT_rectifying']

        print(f"进料板处空间高度系数: {HF_factor}")
        print(f"正常板间距: {self.design_config['HT_rectifying']:.2f} m")
        print(f"进料板处空间高度 HF = {HF_factor} × {self.design_config['HT_rectifying']:.2f}")
        print(f"                     = {HF:.2f} m")
        print(f"设计参考:")
        print(f"  • 系数范围: 1.2-1.5")
        print(f"  • 作用: 保证进料分布和两相分离")

        # 5. 板间距设置
        print(f"\n5. 板间距设置")
        print("-" * 40)

        HT_rect = self.design_config['HT_rectifying']
        HT_strip = self.design_config['HT_stripping']

        print(f"精馏段正常板间距: HT_rect = {HT_rect:.2f} m")
        print(f"提馏段正常板间距: HT_strip = {HT_strip:.2f} m")
        print(f"人孔处板间距增加量: ΔHT = {HT_manhole_rect - HT_rect:.2f} m")
        print(f"设计参考:")
        print(f"  • 常压塔: 0.3-0.6m")
        print(f"  • 减压塔: 0.45-0.8m")
        print(f"  • 影响因素: 雾沫夹带、操作弹性、塔径")

        print(f"\n6. 塔体总有效高度计算")
        print("-" * 40)

        # 使用您提供的符号公式
        print("塔体总有效高度计算公式（根据分段板间距修正）:")
        print("H = H_D + H_B + n_p × H_p + (n_R - 1 - n_p,R) × H_T,R + (n_S - 1 - n_p,S) × H_T,S + H_F")
        print()

        # 解释符号含义
        print("符号说明:")
        print("  H   : 塔体总有效高度")
        print("  H_D : 塔顶空间高度 = {:.2f} m".format(HD))
        print("  H_B : 塔底空间高度 = {:.2f} m".format(HB))
        print("  n_p : 全塔人孔总数 = {} 个".format(S_total))
        print("  n_p,R: 精馏段人孔数 = {} 个".format(S_rect))
        print("  n_p,S: 提馏段人孔数 = {} 个".format(S_strip))
        print("  H_p : 人孔处塔板间距（精馏段={:.2f}m, 提馏段={:.2f}m）".format(HT_manhole_rect, HT_manhole_strip))
        print("  H_T,R: 精馏段普通板间距 = {:.2f} m".format(HT_rect))
        print("  H_T,S: 提馏段普通板间距 = {:.2f} m".format(HT_strip))
        print("  H_F : 进料板处空间高度 = {:.2f} m".format(HF))
        print("  n_R : 精馏段塔板数 = {} 块".format(self.N1_actual))
        print("  n_S : 提馏段塔板数 = {} 块".format(self.N2_actual))
        print()

        # 详细的分步计算
        print("详细计算步骤:")
        print("-" * 40)

        # 1. 计算精馏段非人孔板总高
        H_rect_normal = (self.N1_actual - 1 - S_rect) * HT_rect
        print("(1) 精馏段非人孔板总高:")
        print(f"    H_rect_normal = (n_R - 1 - n_p,R) × H_T,R")
        print(f"                  = ({self.N1_actual} - 1 - {S_rect}) × {HT_rect:.2f}")
        print(f"                  = {self.N1_actual - 1 - S_rect} × {HT_rect:.2f}")
        print(f"                  = {H_rect_normal:.2f} m")

        # 2. 计算提馏段非人孔板总高
        H_strip_normal = (self.N2_actual - 1 - S_strip) * HT_strip
        print(f"\n(2) 提馏段非人孔板总高:")
        print(f"    H_strip_normal = (n_S - 1 - n_p,S) × H_T,S")
        print(f"                   = ({self.N2_actual} - 1 - {S_strip}) × {HT_strip:.2f}")
        print(f"                   = {self.N2_actual - 1 - S_strip} × {HT_strip:.2f}")
        print(f"                   = {H_strip_normal:.2f} m")

        # 3. 计算人孔总高（注意：精馏段和提馏段人孔高度可能不同）
        H_manhole_total = S_rect * HT_manhole_rect + S_strip * HT_manhole_strip
        print(f"\n(3) 全塔人孔板总高:")
        print(f"    H_manhole_total = n_p,R × H_p,R + n_p,S × H_p,S")
        print(f"                    = {S_rect} × {HT_manhole_rect:.2f} + {S_strip} × {HT_manhole_strip:.2f}")
        print(f"                    = {S_rect * HT_manhole_rect:.2f} + {S_strip * HT_manhole_strip:.2f}")
        print(f"                    = {H_manhole_total:.2f} m")

        # 4. 汇总计算（不含安全余量）
        print(f"\n(4) 塔体总有效高度计算:")
        print(f"    H = H_D + H_B + n_p × H_p + H_rect_normal + H_strip_normal + H_F")
        print(
            f"      = {HD:.2f} + {HB:.2f} + {H_manhole_total:.2f} + {H_rect_normal:.2f} + {H_strip_normal:.2f} + {HF:.2f}")
        print(f"      = {HD:.2f} + {HB:.2f} + {H_manhole_total:.2f} + {H_rect_normal + H_strip_normal:.2f} + {HF:.2f}")
        print(f"      = {HD + HB + H_manhole_total + H_rect_normal + H_strip_normal + HF:.2f} m")

        H_total = HD + HB + H_manhole_total + H_rect_normal + H_strip_normal + HF

        # 5. 考虑安全余量
        safety_margin = self.design_config['safety_margin']
        H_total_with_margin = H_total * (1 + safety_margin)
        print(f"\n(5) 考虑安全余量:")
        print(f"    安全余量: {safety_margin * 100:.0f}%")
        print(f"    H_design = H × (1 + 安全余量)")
        print(f"             = {H_total:.2f} × (1 + {safety_margin:.2f})")
        print(f"             = {H_total_with_margin:.2f} m")
        print(f"    设计参考: 5%-15%，考虑安装误差、热膨胀等")

        # 7. 存储计算结果（更新部分）
        self.height_results = {
            'HD': round(HD, 3),
            'h1': round(h1, 3),
            'h2': round(h2, 3),
            'HB': round(HB, 3),
            'S_total': S_total,
            'S_rectifying': S_rect,
            'S_stripping': S_strip,
            'manhole_positions': manhole_positions,
            'manhole_diameter': self.design_config['manhole_diameter'],
            'HT_rectifying': HT_rect,
            'HT_stripping': HT_strip,
            'HT_manhole_rectifying': HT_manhole_rect,
            'HT_manhole_stripping': HT_manhole_strip,
            'HF': round(HF, 3),
            #'H_rectifying': round(H_rect, 3),
            #'H_stripping': round(H_strip, 3),
            'H_rect_normal': round(H_rect_normal, 3),  # 新增
            'H_strip_normal': round(H_strip_normal, 3),  # 新增
            'H_manhole_total': round(H_manhole_total, 3),  # 新增
            'H_total': round(H_total, 3),
            'H_design': round(H_total_with_margin, 3),
            'safety_margin': safety_margin
        }

        return self.height_results

    def display_design_summary(self):
        """显示设计结果汇总"""
        results = self.height_results

        print("\n" + "=" * 70)
        print("塔高设计结果汇总")
        print("=" * 70)

        print(f"\n一、塔板布置")
        print("-" * 40)
        print(f"总实际塔板数: {self.N_total} 块")
        print(f"  精馏段: {self.N1_actual} 块")
        print(f"  提馏段: {self.N2_actual} 块")
        print(f"进料板位置: 第 {self.feed_stage} 块塔板")

        print(f"\n二、塔体各部分高度")
        print("-" * 40)
        print(f"1. 塔顶空间:")
        print(f"   高度 HD = {results['HD']:.2f} m")
        print(f"   范围: {self.design_config['HD_min']:.1f} ~ {self.design_config['HD_max']:.1f} m")
        print(f"   参考范围: 常压塔1.0-1.5m, 减压塔1.5-2.0m")

        print(f"\n2. 塔釜空间:")
        print(f"   釜液面距底层塔板 h1 = {results['h1']:.2f} m")
        print(f"     范围: {self.design_config['h1_min']:.1f} ~ {self.design_config['h1_max']:.1f} m")
        print(f"     参考范围: 0.5-1.2m, 推荐0.6-0.8m")
        print(f"   釜内液层高度 h2 = {results['h2']:.2f} m")
        print(f"     装填系数 φ = {self.design_config['fill_factor']:.2f}")
        print(f"     参考范围: 常压塔0.6-0.7, 减压塔0.5-0.6")
        print(f"   塔釜总高度 HB = h1 + h2 = {results['HB']:.2f} m")
        print(f"     要求: > {self.design_config['HB_min']:.1f} m")
        print(f"     参考范围: 小塔1.0-1.5m, 中塔1.5-2.0m, 大塔2.0-2.5m")

        print(f"\n3. 人孔布置:")
        print(f"   总人孔数: {results['S_total']} 个")
        print(f"     精馏段: {results['S_rectifying']} 个")
        print(f"     提馏段: {results['S_stripping']} 个")
        print(f"   人孔位置: {results['manhole_positions']}")
        print(f"   人孔直径: Φ{results['manhole_diameter'] * 1000:.0f} mm")
        print(f"   参考选择: D<1m:Φ400, D=1-1.8m:Φ450, D>1.8m:Φ500")
        print(f"   人孔处板间距:")
        print(f"     精馏段: {results['HT_manhole_rectifying']:.2f} m")
        print(f"     提馏段: {results['HT_manhole_stripping']:.2f} m")
        print(f"     参考: 比正常板间距大150-200mm")

        print(f"\n4. 板间距:")
        print(f"   精馏段正常板间距: {results['HT_rectifying']:.2f} m")
        print(f"   提馏段正常板间距: {results['HT_stripping']:.2f} m")
        print(f"   参考范围: 常压塔0.3-0.6m, 减压塔0.45-0.8m")

        print(f"\n5. 进料板处空间:")
        print(f"   高度 HF = {results['HF']:.2f} m")
        print(f"   系数: {self.design_config['HF_factor']} × 正常板间距")
        print(f"   参考范围: 系数1.2-1.5")

        print(f"\n三、塔体总高度计算")
        print("-" * 40)
        #print(f"精馏段高度: {results['H_rectifying']:.2f} m")
        #print(f"提馏段高度: {results['H_stripping']:.2f} m")
        print(f"塔体有效高度: {results['H_total']:.2f} m")
        print(f"考虑安全余量 ({results['safety_margin'] * 100:.0f}%):")
        print(f"设计高度 H_design = {results['H_design']:.2f} m")
        print(f"安全余量参考范围: 5%-15%")

        print(f"\n四、设计参数")
        print("-" * 40)
        print(f"安全余量: {results['safety_margin'] * 100:.0f}%")
        print(f"人孔间距: 每隔 {self.design_config['manhole_interval']} 块塔板设一个人孔")
        print(f"装填系数: {self.design_config['fill_factor']:.2f}")

        return results

    def save_results_to_json(self, filename='tower_height_design.json'):
        """保存设计结果到JSON文件"""
        # 创建完整的设计结果
        design_results = {
            "设计参数参考范围": {
                "塔顶空间高度_m": {
                    "取值": self.height_results['HD'],
                    "范围": f"{self.design_config['HD_min']:.1f}-{self.design_config['HD_max']:.1f}",
                    "参考": "常压塔1.0-1.5m, 减压塔1.5-2.0m"
                },
                "塔釜液面距底层塔板_m": {
                    "取值": self.height_results['h1'],
                    "范围": f"{self.design_config['h1_min']:.1f}-{self.design_config['h1_max']:.1f}",
                    "参考": "0.5-1.2m, 推荐0.6-0.8m"
                },
                "装填系数": {
                    "取值": self.design_config['fill_factor'],
                    "范围": "0.5-0.7",
                    "参考": "常压塔0.6-0.7, 减压塔0.5-0.6"
                },
                "板间距_m": {
                    "取值": f"{self.height_results['HT_rectifying']:.3f}/{self.height_results['HT_stripping']:.3f}",
                    "范围": "0.3-0.8",
                    "参考": "常压塔0.3-0.6m, 减压塔0.45-0.8m"
                },
                "安全余量": {
                    "取值": self.height_results['safety_margin'],
                    "范围": "0.05-0.15",
                    "参考": "5%-15%"
                }
            },
            "塔板布置": {
                "总塔板数": self.N_total,
                "精馏段塔板数": self.N1_actual,
                "提馏段塔板数": self.N2_actual,
                "进料板位置": self.feed_stage
            },
            "塔体各部分高度": {
                "塔顶空间高度_m": self.height_results['HD'],
                "塔釜液面距底层塔板_m": self.height_results['h1'],
                "釜内液层高度_m": self.height_results['h2'],
                "塔釜总高度_m": self.height_results['HB'],
                "装填系数": self.design_config['fill_factor']
            },
            "人孔设计": {
                "总人孔数": self.height_results['S_total'],
                "精馏段人孔数": self.height_results['S_rectifying'],
                "提馏段人孔数": self.height_results['S_stripping'],
                "人孔位置": self.height_results['manhole_positions'],
                "人孔直径_m": self.height_results['manhole_diameter'],
                "人孔处板间距_精馏段_m": self.height_results['HT_manhole_rectifying'],
                "人孔处板间距_提馏段_m": self.height_results['HT_manhole_stripping']
            },
            "板间距设计": {
                "正常板间距_精馏段_m": self.height_results['HT_rectifying'],
                "正常板间距_提馏段_m": self.height_results['HT_stripping']
            },
            "进料板处空间": {
                "高度_m": self.height_results['HF'],
                "系数": self.design_config['HF_factor']
            },
            "高度计算结果": {
                #"精馏段高度_m": self.height_results['H_rectifying'],
                #"提馏段高度_m": self.height_results['H_stripping'],
                "塔体有效高度_m": self.height_results['H_total'],
                "设计高度_m": self.height_results['H_design'],
                "安全余量": self.height_results['safety_margin']
            },
            "设计参数": {
                "安全余量": self.design_config['safety_margin'],
                "人孔间距": self.design_config['manhole_interval'],
                "装填系数": self.design_config['fill_factor']
            }
        }

        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(design_results, f, indent=4, ensure_ascii=False)

        print(f"\n设计结果已保存到: {filename}")

        return filename

    def export_design_formula(self):
        """导出设计计算公式"""
        formula_text = f"""
塔体总有效高度计算公式:

H = HD + HB + S_rect×HT'_rect + S_strip×HT'_strip + 
    (N1 - 1 - S_rect)×HT_rect + (N2 - 1 - S_strip)×HT_strip + HF

其中:
  HD  = 塔顶空间高度 = {self.height_results['HD']:.2f} m
      (参考范围: 常压塔1.0-1.5m, 减压塔1.5-2.0m)

  HB  = 塔釜总高度 = h1 + h2 = {self.height_results['h1']:.2f} + {self.height_results['h2']:.2f} = {self.height_results['HB']:.2f} m
      (h1参考: 0.5-1.2m, 推荐0.6-0.8m)

  S_rect = 精馏段人孔数 = {self.height_results['S_rectifying']} 个
      (参考间距: 8-12块板设一个)

  S_strip = 提馏段人孔数 = {self.height_results['S_stripping']} 个

  HT'_rect = 精馏段人孔处板间距 = {self.height_results['HT_manhole_rectifying']:.2f} m
      (参考: 比正常板间距大150-200mm)

  HT'_strip = 提馏段人孔处板间距 = {self.height_results['HT_manhole_stripping']:.2f} m

  N1 = 精馏段塔板数 = {self.N1_actual} 块

  N2 = 提馏段塔板数 = {self.N2_actual} 块

  HT_rect = 精馏段正常板间距 = {self.height_results['HT_rectifying']:.2f} m
      (参考范围: 常压塔0.3-0.6m, 减压塔0.45-0.8m)

  HT_strip = 提馏段正常板间距 = {self.height_results['HT_stripping']:.2f} m

  HF = 进料板处空间高度 = {self.height_results['HF']:.2f} m
      (参考系数: 1.2-1.5)

代入计算:
H = {self.height_results['HD']:.2f} + {self.height_results['HB']:.2f} + 
    {self.height_results['S_rectifying']}×{self.height_results['HT_manhole_rectifying']:.2f} + 
    {self.height_results['S_stripping']}×{self.height_results['HT_manhole_stripping']:.2f} + 
    ({self.N1_actual} - 1 - {self.height_results['S_rectifying']})×{self.height_results['HT_rectifying']:.2f} + 
    ({self.N2_actual} - 1 - {self.height_results['S_stripping']})×{self.height_results['HT_stripping']:.2f} + 
    {self.height_results['HF']:.2f}

H = {self.height_results['H_total']:.2f} m

考虑安全余量 ({self.height_results['safety_margin'] * 100:.0f}%):
H_design = {self.height_results['H_total']:.2f} × (1 + {self.height_results['safety_margin']:.2f})
         = {self.height_results['H_design']:.2f} m
         (安全余量参考范围: 5%-15%)
"""
        # 保存到文件
        with open('tower_height_formula.txt', 'w', encoding='utf-8') as f:
            f.write(formula_text)

        print(f"设计公式已保存到: tower_height_formula.txt")

        return formula_text


def main():
    """主函数"""
    print("塔高设计计算程序")
    print("=" * 70)
    print("设计参考标准:")
    print("  • 《化工容器设计》")
    print("  • 《塔器设计手册》")
    print("  • 《石油化工设计手册》")
    print("  • HG/T 20592-2009 钢制管法兰、垫片、紧固件")
    print("=" * 70)

    # 询问是否使用交互模式
    use_interactive = input("是否使用交互模式输入设计参数? (y/n) [默认y]: ").lower() or 'y'
    interactive_mode = use_interactive == 'y'

    try:
        # 创建设计计算器
        calculator = TowerHeightCalculator(interactive_mode=interactive_mode)

        # 计算塔体高度
        height_results = calculator.calculate_tower_height()

        # 显示设计结果汇总
        calculator.display_design_summary()

        # 保存设计结果
        calculator.save_results_to_json()

        # 导出设计公式
        calculator.export_design_formula()

        print("\n" + "=" * 70)
        print("塔高设计计算完成!")
        print("=" * 70)

        return calculator

    except FileNotFoundError as e:
        print(f"错误: 未找到数据文件 - {e}")
        print("请确保已运行之前的计算程序并生成所需的数据文件。")
        return None
    except Exception as e:
        print(f"设计过程中发生错误: {e}")
        import traceback
        traceback.print_exc()
        return None


# 如果直接运行此文件
if __name__ == "__main__":
    calculator = main()

    if calculator:
        # 打印关键设计结果
        print("\n关键设计结果:")
        print("-" * 40)
        print(f"塔体总有效高度: {calculator.height_results['H_total']:.2f} m")
        print(f"设计高度(含安全余量): {calculator.height_results['H_design']:.2f} m")
        print(f"塔板数: {calculator.N_total} 块")
        print(f"人孔数: {calculator.height_results['S_total']} 个")
        print(f"塔釜高度: {calculator.height_results['HB']:.2f} m")
        print(f"塔顶空间: {calculator.height_results['HD']:.2f} m")