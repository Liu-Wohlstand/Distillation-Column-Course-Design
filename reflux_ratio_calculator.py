"""
增强版回流比计算模块
考虑过冷回流的影响
"""

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei']
matplotlib.rcParams['axes.unicode_minus'] = False


class EnhancedRefluxCalculator:
    """考虑过冷回流的增强版回流比计算器"""

    def __init__(self, equilibrium_data, Cp_data=None, r_data=None, M_components=None):
        """
        初始化计算器

        参数:
        equilibrium_data: 相平衡数据，字典格式 {x: [y, temperature]}
        Cp_data: 比热容数据，字典格式 {component: {T: Cp}}
        r_data: 汽化潜热数据，字典格式 {component: r}
        M_components: 组分摩尔质量，字典格式 {component: M}
        """
        # 相平衡数据
        self.equilibrium_data = equilibrium_data
        self.x_data = list(equilibrium_data.keys())
        self.y_data = [equilibrium_data[x][0] for x in self.x_data]
        self.t_data = [equilibrium_data[x][1] for x in self.x_data]

        # 排序
        sorted_idx = np.argsort(self.x_data)
        self.x_data = np.array([self.x_data[i] for i in sorted_idx])
        self.y_data = np.array([self.y_data[i] for i in sorted_idx])
        self.t_data = np.array([self.t_data[i] for i in sorted_idx])

        # 创建插值函数
        self.y_from_x = interp1d(self.x_data, self.y_data, kind='cubic', fill_value="extrapolate")
        self.x_from_y = interp1d(self.y_data, self.x_data, kind='cubic', fill_value="extrapolate")
        self.t_from_x = interp1d(self.x_data, self.t_data, kind='linear', fill_value="extrapolate")

        # 物性数据
        self.Cp_data = Cp_data or {}
        self.r_data = r_data or {}
        self.M_components = M_components or {}

        # 设计参数（由主程序设置）
        self.xD = None
        self.xF = None
        self.xW = None
        self.q = None
        self.R_min = None
        self.R_multiplier = None
        self.R_actual = None
        self.R_internal = None
        self.reflux_temp = None
        self.feed_stage = None
        self.total_stages = None
        self.ET = None  # 塔板效率
        self.stage_positions = []  # 阶梯图结果

    def set_design_parameters(self, xD, xF, xW, q, R_min, R_multiplier, 
                             R_actual, R_internal, reflux_temp, 
                             feed_stage, total_stages, ET=None, stage_positions=None):
        """设置主程序的设计参数（完整版）"""
        self.xD = xD
        self.xF = xF
        self.xW = xW
        self.q = q
        self.R_min = R_min
        self.R_multiplier = R_multiplier
        self.R_actual = R_actual
        self.R_internal = R_internal
        self.reflux_temp = reflux_temp
        self.feed_stage = feed_stage
        self.total_stages = total_stages
        self.ET = ET or 0.5  # 默认塔板效率
        self.stage_positions = stage_positions or []

    def calculate_R_min(self, xD, xF, q=1.0, method='graphical'):
        """
        计算最小回流比（理论值，不考虑过冷）
        参数同上一个版本
        """
        # 计算q线与平衡线的交点
        x_eq, y_eq = self._find_q_equilibrium_intersection(xF, q)

        if method == 'analytical':
            if abs(q - 1) < 1e-6:  # q=1，泡点进料
                y_eq = float(self.y_from_x(xF))
                R_min = (xD - y_eq) / (y_eq - xF)
            else:
                R_min = (xD - y_eq) / (y_eq - x_eq)

        elif method == 'graphical':
            if x_eq == xD:
                R_min = float('inf')
            else:
                R_min = (xD - y_eq) / (y_eq - x_eq)

        else:
            raise ValueError("method参数必须为 'graphical' 或 'analytical'")

        return R_min, (x_eq, y_eq)

    def calculate_actual_reflux_ratio(self, R_min, R_multiplier=1.2,
                                      reflux_temp=None, top_temp=None,
                                      x_top=None):
        """
        计算考虑过冷回流的实际回流比

        参数:
        R_min: 最小回流比（理论值）
        R_multiplier: 回流比倍数
        reflux_temp: 回流液温度 (°C)
        top_temp: 塔顶温度 (°C)，如果为None则根据x_top计算
        x_top: 塔顶组成，用于计算塔顶温度

        返回:
        R_actual: 实际外回流比
        R_internal: 实际内回流比
        correction_factor: 修正系数
        """
        # 计算理论操作回流比
        R_theoretical = R_min * R_multiplier

        # 如果没有指定过冷回流温度，直接返回理论值
        if reflux_temp is None or top_temp is None:
            if top_temp is None and x_top is not None:
                top_temp = float(self.t_from_x(x_top))
            else:
                return R_theoretical, R_theoretical, 1.0

        # 计算过冷回流修正
        correction_factor = self._calculate_subcooling_correction(
            reflux_temp, top_temp, x_top
        )

        # 计算内回流比
        R_internal = R_theoretical * correction_factor

        # 实际外回流比（保持不变）
        R_actual = R_theoretical

        return R_actual, R_internal, correction_factor

    def _calculate_subcooling_correction(self, T_reflux, T_top, x_top):
        """
        计算过冷回流修正系数

        修正公式: R' = R * [1 + Cp*(T_top - T_reflux)/ΔH_v]
        其中:
        R' : 内回流比
        R  : 外回流比
        """
        if not self.Cp_data or not self.r_data or not self.M_components:
            print("警告：缺少物性数据，无法计算过冷回流修正")
            return 1.0

        # 计算塔顶组成的物性
        # 假设组分A为轻组分（苯），B为重组分（甲苯）
        component_A = 'benzene'
        component_B = 'toluene'

        # 计算平均温度下的比热容
        T_avg = (T_reflux + T_top) / 2

        # 插值获取比热容
        if component_A in self.Cp_data and component_B in self.Cp_data:
            # 获取比热容 (kJ/(kg·K))
            Cp_A = self._interpolate_property(T_avg, self.Cp_data[component_A])
            Cp_B = self._interpolate_property(T_avg, self.Cp_data[component_B])

            # 转换为摩尔比热容 (kJ/(kmol·K))
            M_A = self.M_components.get(component_A, 78.11)
            M_B = self.M_components.get(component_B, 92.14)
            Cp_A_mol = Cp_A * M_A
            Cp_B_mol = Cp_B * M_B

            # 混合物摩尔比热容
            Cp_mix = x_top * Cp_A_mol + (1 - x_top) * Cp_B_mol
        else:
            # 如果没有数据，使用估算值
            Cp_mix = 150  # kJ/(kmol·K) 的估算值
            print(f"使用估算的比热容: {Cp_mix} kJ/(kmol·K)")

        # 获取汽化潜热
        if component_A in self.r_data and component_B in self.r_data:
            r_A = self.r_data[component_A] * self.M_components.get(component_A, 78.11)
            r_B = self.r_data[component_B] * self.M_components.get(component_B, 92.14)
            r_mix = x_top * r_A + (1 - x_top) * r_B
        else:
            # 如果没有数据，使用估算值
            r_mix = 30000  # kJ/kmol 的估算值
            print(f"使用估算的汽化潜热: {r_mix} kJ/kmol")

        # 计算修正系数
        delta_T = T_top - T_reflux  # K
        correction = 1 + (Cp_mix * delta_T) / r_mix

        print(f"\n过冷回流修正计算:")
        print(f"  回流液温度: {T_reflux} °C")
        print(f"  塔顶温度: {T_top} °C")
        print(f"  温差 ΔT: {delta_T:.1f} K")
        print(f"  混合物比热容 Cp: {Cp_mix:.1f} kJ/(kmol·K)")
        print(f"  混合物汽化潜热 ΔHv: {r_mix:.0f} kJ/kmol")
        print(f"  修正系数: {correction:.3f}")

        return correction

    def _interpolate_property(self, T, property_dict):
        """插值获取物性数据"""
        temperatures = list(property_dict.keys())
        values = list(property_dict.values())

        if T <= min(temperatures):
            return values[0]
        elif T >= max(temperatures):
            return values[-1]
        else:
            # 线性插值
            return np.interp(T, temperatures, values)

    def _find_q_equilibrium_intersection(self, xF, q):
        """寻找q线与平衡线的交点（同上个版本）"""
        if abs(q - 1) < 1e-6:
            y_eq = float(self.y_from_x(xF))
            return xF, y_eq

        if abs(q) < 1e-6:
            try:
                x_eq = float(self.x_from_y(xF))
                return x_eq, xF
            except:
                return self._numerical_intersection(xF, q)

        return self._numerical_intersection(xF, q)

    def _numerical_intersection(self, xF, q):
        """数值求解交点"""

        def q_line_y(x):
            return q / (q - 1) * x - xF / (q - 1)

        x_range = np.linspace(0, 1, 1000)
        min_diff = float('inf')
        best_x, best_y = 0, 0

        for x in x_range:
            if x < 0 or x > 1:
                continue

            try:
                y_eq = float(self.y_from_x(x))
            except:
                continue

            y_q = q_line_y(x)
            diff = abs(y_eq - y_q)

            if diff < min_diff:
                min_diff = diff
                best_x = x
                best_y = (y_eq + y_q) / 2

            if diff < 1e-6:
                break

        return best_x, best_y

    def perform_mccabe_thiele_analysis(self, xD, xF, xW, q, R_internal, 
                                       custom_stage_positions=None):
        """
        执行完整的McCabe-Thiele阶梯图分析
        支持使用自定义的阶梯图结果
        """
        
        # 如果提供了自定义阶梯图结果，直接使用
        if custom_stage_positions and len(custom_stage_positions) > 0:
            stages = []
            for i, stage_info in enumerate(custom_stage_positions):
                if len(stage_info) == 3:
                    section, x_break, y_zhezhong = stage_info
                    stages.append({
                        'stage': i + 1,
                        'x_start': x_break,
                        'y_op': y_zhezhong,
                        'x_eq': self._find_x_from_y(y_zhezhong),
                        'is_feed': (section == '提馏段' and i == 0) or 
                                   (section == '精馏段' and i > 0 and stages[-1].get('is_feed', False))
                    })
                elif len(stage_info) == 4:
                    section, x_break, y_zheqi, y_zhezhong = stage_info
                    stages.append({
                        'stage': i + 1,
                        'x_start': x_break,
                        'y_op': y_zheqi,
                        'x_eq': self._find_x_from_y(y_zheqi),
                        'is_feed': (section == '提馏段' and i == 0) or 
                                   (section == '精馏段' and i > 0 and stages[-1].get('is_feed', False))
                    })
            
            # 计算总理论板数
            total_stages = len(stages)
            
            # 查找进料板位置
            feed_stage = None
            for i, stage in enumerate(stages):
                if stage.get('is_feed', False):
                    feed_stage = i + 1
                    break
            
            # 计算操作线交点
            if abs(q - 1) < 1e-6:  # q=1
                x_intersection = xF
                y_intersection = R_internal/(R_internal+1) * xF + xD/(R_internal+1)
            else:
                # q线与精馏段操作线交点
                A1 = R_internal/(R_internal+1)
                B1 = xD/(R_internal+1)
                A2 = q/(q-1)
                B2 = -xF/(q-1)
                x_intersection = (B2 - B1)/(A1 - A2)
                y_intersection = A1 * x_intersection + B1
            
            return {
                'stages': stages,
                'total_stages': total_stages,
                'feed_stage': feed_stage or total_stages,
                'intersection': (x_intersection, y_intersection)
            }
        
        # 否则使用原有的计算方法
        return self._calculate_mccabe_thiele_analysis(xD, xF, xW, q, R_internal)

    def _find_x_from_y(self, y_val):
        """从y值找到对应的x值（平衡线上）"""
        try:
            return float(self.x_from_y(y_val))
        except:
            # 如果超出插值范围，使用线性插值
            if y_val >= max(self.y_data):
                return 1.0
            elif y_val <= min(self.y_data):
                return 0.0
            else:
                # 线性插值
                idx = np.searchsorted(self.y_data, y_val)
                x1, x2 = self.x_data[idx-1], self.x_data[idx]
                y1, y2 = self.y_data[idx-1], self.y_data[idx]
                return x1 + (x2 - x1) * (y_val - y1) / (y2 - y1)

    def _calculate_mccabe_thiele_analysis(self, xD, xF, xW, q, R_internal):
        """原始的计算方法"""
        # 1. 计算操作线交点
        if abs(q - 1) < 1e-6:  # q=1
            x_intersection = xF
            y_intersection = R_internal/(R_internal+1) * xF + xD/(R_internal+1)
        else:
            # q线与精馏段操作线交点
            A1 = R_internal/(R_internal+1)
            B1 = xD/(R_internal+1)
            A2 = q/(q-1)
            B2 = -xF/(q-1)
            x_intersection = (B2 - B1)/(A1 - A2)
            y_intersection = A1 * x_intersection + B1

        # 2. 阶梯图解法
        stages = []
        x_current = xD
        stage_count = 0

        while x_current > xW + 1e-4 and stage_count < 50:
            stage_count += 1

            # 从操作线到平衡线（垂直）
            if stage_count == 1:
                # 第一块板从塔顶开始
                y_current = xD
            elif x_current >= x_intersection:
                # 在精馏段
                y_current = R_internal/(R_internal+1) * x_current + xD/(R_internal+1)
            else:
                # 在提馏段
                slope = (y_intersection - xW) / (x_intersection - xW)
                y_current = slope * (x_current - xW) + xW

            # 从平衡线到操作线（水平）
            x_next = self._find_x_from_y(y_current)

            stages.append({
                'stage': stage_count,
                'x_start': x_current,
                'y_op': y_current,
                'x_eq': x_next,
                'is_feed': False
            })

            # 检查是否到达进料板
            if x_current > x_intersection and x_next <= x_intersection:
                stages[-1]['is_feed'] = True

            x_current = x_next

        return {
            'stages': stages,
            'total_stages': len(stages),
            'feed_stage': next((i+1 for i, s in enumerate(stages) if s['is_feed']), len(stages)),
            'intersection': (x_intersection, y_intersection)
        }

    def plot_enhanced_diagram(self, xD=None, xF=None, xW=None, q=1.0, R_multiplier=1.2,
                              reflux_temp=None, R_internal=None,
                              custom_stage_positions=None,
                              figsize=(16, 12), save_path=None):
        """
        绘制考虑过冷回流的完整分析图，包括阶梯图
        优先使用已设置的参数，如果没有则使用传入的参数
        """
        # 优先使用已设置的参数
        xD = xD if xD is not None else self.xD
        xF = xF if xF is not None else self.xF
        xW = xW if xW is not None else self.xW
        q = q if q is not None else self.q
        reflux_temp = reflux_temp if reflux_temp is not None else self.reflux_temp
        R_multiplier = R_multiplier if R_multiplier is not None else self.R_multiplier
        
        if xD is None or xF is None or xW is None:
            raise ValueError("必须提供xD, xF, xW参数或先调用set_design_parameters")

        # 计算理论最小回流比
        R_min, (x_eq, y_eq) = self.calculate_R_min(xD, xF, q)

        # 计算塔顶温度
        top_temp = float(self.t_from_x(xD))

        # 计算实际回流比（考虑过冷）
        if R_internal is None:
            if self.R_actual is not None and self.R_internal is not None:
                # 使用已设置的回流比
                R_actual = self.R_actual
                R_internal = self.R_internal
                correction = R_internal / R_actual if R_actual > 0 else 1.0
            else:
                # 计算新的回流比
                R_actual, R_internal, correction = self.calculate_actual_reflux_ratio(
                    R_min, R_multiplier, reflux_temp, top_temp, xD
                )
        else:
            R_actual = R_min * R_multiplier
            correction = R_internal / R_actual if R_actual > 0 else 1.0

        # 执行McCabe-Thiele分析
        mccabe_results = self.perform_mccabe_thiele_analysis(
            xD, xF, xW, q, R_internal, custom_stage_positions
        )

        # 更新设计参数
        self.set_design_parameters(
            xD, xF, xW, q, R_min, R_multiplier, R_actual, R_internal,
            reflux_temp, mccabe_results['feed_stage'], 
            mccabe_results['total_stages'], 
            self.ET, custom_stage_positions
        )

        # 创建图形
        fig, axes = plt.subplots(2, 3, figsize=figsize)
        axes = axes.flatten()

        # 1. 平衡线与操作线图
        ax1 = axes[0]
        self._plot_equilibrium_operating_lines(ax1, xD, xF, xW, q, R_min, R_actual, R_internal,
                                              mccabe_results['intersection'][0],
                                              mccabe_results['intersection'][1])

        # 2. 温度分布图
        ax2 = axes[1]
        self._plot_temperature_distribution(ax2, xD, xF, xW, reflux_temp, top_temp)

        # 3. 回流比影响图
        ax3 = axes[2]
        self._plot_reflux_ratio_effects(ax3, R_min, R_actual, R_internal, correction)

        # 4. 能量分析图
        ax4 = axes[3]
        self._plot_energy_analysis(ax4, R_actual, R_internal, xD, reflux_temp, top_temp)

        # 5. McCabe-Thiele阶梯图
        ax5 = axes[4]
        self._plot_mccabe_thiele_diagram(ax5, xD, xF, xW, q, R_internal, mccabe_results)

        # 6. 组成分布图
        ax6 = axes[5]
        self._plot_composition_profile(ax6, mccabe_results['stages'], xF)

        # 设置总标题
        if reflux_temp:
            title = (f'苯-甲苯精馏塔设计分析\n'
                     f'xD={xD:.3f}, xF={xF:.3f}, xW={xW:.3f}, q={q:.2f}, '
                     f'回流温度={reflux_temp}°C\n'
                     f'R_min={R_min:.3f}, R_实际={R_actual:.3f}, R_内={R_internal:.3f}, '
                     f'理论板数={mccabe_results["total_stages"]}')
        else:
            title = (f'苯-甲苯精馏塔设计分析\n'
                     f'xD={xD:.3f}, xF={xF:.3f}, xW={xW:.3f}, q={q:.2f}, '
                     f'理论板数={mccabe_results["total_stages"]}')

        plt.suptitle(title, fontsize=14, fontweight='bold')
        plt.tight_layout(rect=[0, 0, 1, 0.96])  # 为总标题留出空间

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"图表已保存至: {save_path}")

        return fig, axes, R_min, R_actual, R_internal, mccabe_results

    def _plot_equilibrium_operating_lines(self, ax, xD, xF, xW, q, R_min, R_actual, R_internal, x_intersection, y_intersection):
        """绘制平衡线与操作线"""
        # 平衡线
        x_smooth = np.linspace(0, 1, 200)
        y_smooth = [float(self.y_from_x(x)) for x in x_smooth]

        ax.plot(x_smooth, y_smooth, 'b-', linewidth=2, label='平衡线')
        ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, label='对角线')

        # q线
        if abs(q - 1) < 1e-6:
            ax.axvline(x=xF, color='purple', linestyle='--', linewidth=1.5, label=f'q线 (x={xF:.3f})')
        else:
            x_q = np.linspace(0, 1, 100)
            y_q = q / (q - 1) * x_q - xF / (q - 1)
            valid_idx = (y_q >= 0) & (y_q <= 1)
            ax.plot(x_q[valid_idx], y_q[valid_idx], 'purple', linestyle='--',
                    linewidth=1.5, label=f'q线 (q={q:.2f})')

        # 标记交点
        ax.scatter([x_intersection], [y_intersection], color='red', s=80, zorder=5,
                   label=f'交点 (x={x_intersection:.3f}, y={y_intersection:.3f})')

        # 最小回流比操作线
        if R_min != float('inf'):
            ax.plot([x_intersection, xD], [y_intersection, xD], 'r-', linewidth=1.5,
                    label=f'最小回流比线 (R_min={R_min:.3f})')

        # 实际外回流比操作线（精馏段）
        slope_actual = R_actual / (R_actual + 1)
        intercept_actual = xD / (R_actual + 1)
        x_op_actual = np.linspace(0, xD, 100)
        y_op_actual = slope_actual * x_op_actual + intercept_actual
        ax.plot(x_op_actual, y_op_actual, 'g-', linewidth=2,
                label=f'外回流比线 (R={R_actual:.3f})')

        # 内回流比操作线（精馏段）
        slope_internal = R_internal / (R_internal + 1)
        intercept_internal = xD / (R_internal + 1)
        x_op_internal = np.linspace(0, xD, 100)
        y_op_internal = slope_internal * x_op_internal + intercept_internal
        ax.plot(x_op_internal, y_op_internal, 'orange', linestyle='--', linewidth=2,
                label=f'内回流比线 (R\'={R_internal:.3f})')

        # 提馏段操作线
        slope_stripping = (y_intersection - xW) / (x_intersection - xW)
        x_op_stripping = np.linspace(xW, x_intersection, 100)
        y_op_stripping = slope_stripping * (x_op_stripping - xW) + xW
        ax.plot(x_op_stripping, y_op_stripping, 'brown', linewidth=2,
                label='提馏段操作线')

        # 标记关键点
        ax.scatter([xD, xF, xW], [xD, xF, xW], color='blue', s=60, zorder=5)
        ax.text(xD, xD, f' 塔顶', verticalalignment='bottom', fontsize=9)
        ax.text(xF, xF, f' 进料', verticalalignment='top', fontsize=9)
        ax.text(xW, xW, f' 塔釜', verticalalignment='top', fontsize=9)

        ax.set_xlabel('液相组成 x', fontsize=11)
        ax.set_ylabel('气相组成 y', fontsize=11)
        ax.set_title('相平衡与操作线', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper left', fontsize=8)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

    def _plot_temperature_distribution(self, ax, xD, xF, xW, reflux_temp, top_temp):
        """绘制温度分布图"""
        # 计算塔釜温度
        t_bottom = float(self.t_from_x(xW))

        # 塔内温度分布（简化）
        stages = np.arange(0, 21)  # 假设20块理论板
        # 非线性温度分布（更符合实际情况）
        x_positions = np.linspace(0, 1, len(stages))
        temperatures = top_temp + (t_bottom - top_temp) * (1 - np.exp(-3*x_positions))/(1 - np.exp(-3))

        ax.plot(stages, temperatures, 'b-o', linewidth=2, markersize=4)

        # 标记关键温度点
        ax.axhline(y=top_temp, color='r', linestyle='--', alpha=0.7,
                   label=f'塔顶温度: {top_temp:.1f}°C')

        if reflux_temp:
            ax.axhline(y=reflux_temp, color='g', linestyle='--', alpha=0.7,
                       label=f'回流温度: {reflux_temp}°C')
            # 填充过冷区域
            ax.fill_between(stages, reflux_temp, top_temp,
                            alpha=0.2, color='green', label='过冷区')

        ax.axhline(y=t_bottom, color='orange', linestyle='--', alpha=0.7,
                   label=f'塔釜温度: {t_bottom:.1f}°C')

        ax.set_xlabel('理论板位置 (从上到下)', fontsize=11)
        ax.set_ylabel('温度 (°C)', fontsize=11)
        ax.set_title('塔内温度分布', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='best', fontsize=8)
        ax.set_xlim(0, 20)

    def _plot_reflux_ratio_effects(self, ax, R_min, R_actual, R_internal, correction):
        """绘制回流比影响图"""
        # 回流比倍数范围
        multipliers = np.linspace(1.1, 2.5, 20)
        R_theoretical = R_min * multipliers

        # 计算对应的内回流比（考虑过冷修正）
        R_internal_vals = R_theoretical * correction

        ax.plot(multipliers, R_theoretical, 'b-', linewidth=2, label='外回流比')
        ax.plot(multipliers, R_internal_vals, 'g--', linewidth=2, label='内回流比')

        # 标记当前设计点
        current_multiplier = R_actual / R_min if R_min > 0 else 1.0
        ax.scatter([current_multiplier], [R_actual], color='red', s=100, zorder=5, label='当前设计点')
        ax.scatter([current_multiplier], [R_internal], color='red', s=100, zorder=5)

        ax.axvline(x=current_multiplier, color='r', linestyle=':', alpha=0.7)

        ax.set_xlabel('回流比倍数 (R/R_min)', fontsize=11)
        ax.set_ylabel('回流比', fontsize=11)
        ax.set_title('回流比对设计的影响', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper left', fontsize=8)

    def _plot_energy_analysis(self, ax, R_actual, R_internal, xD, reflux_temp, top_temp):
        """绘制能量分析图"""
        # 假设的汽化潜热（使用苯-甲苯混合物的平均值）
        delta_Hv = 30000  # kJ/kmol

        # 计算热负荷
        D = 100  # 假设馏出液流量为100 kmol/h
        Q_condenser_theoretical = R_actual * D * delta_Hv
        Q_reboiler_theoretical = (R_actual + 1) * D * delta_Hv

        Q_condenser_actual = R_internal * D * delta_Hv
        Q_reboiler_actual = (R_internal + 1) * D * delta_Hv

        # 计算过冷回流额外热负荷
        if reflux_temp and top_temp:
            # 假设比热容
            Cp = 150  # kJ/(kmol·K)
            delta_T = top_temp - reflux_temp
            Q_subcooling = R_actual * D * Cp * delta_T
        else:
            Q_subcooling = 0

        # 创建柱状图
        categories = ['冷凝器', '再沸器']
        theoretical = [Q_condenser_theoretical / 1000, Q_reboiler_theoretical / 1000]  # 转换为MJ/h
        actual = [Q_condenser_actual / 1000, Q_reboiler_actual / 1000]

        x = np.arange(len(categories))
        width = 0.35

        ax.bar(x - width / 2, theoretical, width, label='理论值', color='blue', alpha=0.7)
        ax.bar(x + width / 2, actual, width, label='实际值(过冷)', color='green', alpha=0.7)

        ax.set_xlabel('设备', fontsize=11)
        ax.set_ylabel('热负荷 (MJ/h)', fontsize=11)
        ax.set_title('热负荷分析', fontsize=12)
        ax.set_xticks(x)
        ax.set_xticklabels(categories)
        ax.grid(True, alpha=0.3, axis='y')
        ax.legend(fontsize=8)

        # 添加数值标签
        for i, (t, a) in enumerate(zip(theoretical, actual)):
            ax.text(i - width / 2, t + max(theoretical) * 0.02, f'{t:.0f}',
                    ha='center', va='bottom', fontsize=8)
            ax.text(i + width / 2, a + max(actual) * 0.02, f'{a:.0f}',
                    ha='center', va='bottom', fontsize=8)

        if Q_subcooling > 0:
            ax.text(0.5, 0.95, f'过冷回流额外热负荷: {Q_subcooling / 1000:.0f} MJ/h',
                    transform=ax.transAxes, ha='center', fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.5))

    def _plot_mccabe_thiele_diagram(self, ax, xD, xF, xW, q, R_internal, mccabe_results):
        """绘制McCabe-Thiele阶梯图"""
        # 平衡线
        x_smooth = np.linspace(0, 1, 200)
        y_smooth = [float(self.y_from_x(x)) for x in x_smooth]

        ax.plot(x_smooth, y_smooth, 'b-', linewidth=2, label='平衡线')
        ax.plot([0, 1], [0, 1], 'k--', alpha=0.3)

        # q线
        if abs(q - 1) < 1e-6:
            ax.axvline(x=xF, color='purple', linestyle='--', linewidth=1.5)
        else:
            x_q = np.linspace(0, 1, 100)
            y_q = q / (q - 1) * x_q - xF / (q - 1)
            valid_idx = (y_q >= 0) & (y_q <= 1)
            ax.plot(x_q[valid_idx], y_q[valid_idx], 'purple', linestyle='--', linewidth=1.5)

        # 精馏段操作线
        slope_rect = R_internal / (R_internal + 1)
        intercept_rect = xD / (R_internal + 1)
        x_rect = np.linspace(0, xD, 100)
        y_rect = slope_rect * x_rect + intercept_rect
        ax.plot(x_rect, y_rect, 'g-', linewidth=2)

        # 提馏段操作线
        x_intersection, y_intersection = mccabe_results['intersection']
        slope_strip = (y_intersection - xW) / (x_intersection - xW)
        x_strip = np.linspace(xW, x_intersection, 100)
        y_strip = slope_strip * (x_strip - xW) + xW
        ax.plot(x_strip, y_strip, 'brown', linewidth=2)

        # 绘制阶梯
        stages = mccabe_results['stages']
        for i, stage in enumerate(stages):
            # 垂直线（从操作线到平衡线）
            ax.plot([stage['x_start'], stage['x_start']],
                   [stage['y_op'], stage['y_op']], 'r-', linewidth=1, alpha=0.7)

            # 水平线（从平衡线到操作线）
            ax.plot([stage['x_start'], stage['x_eq']],
                   [stage['y_op'], stage['y_op']], 'r-', linewidth=1, alpha=0.7)

            # 标记进料板
            if stage['is_feed']:
                ax.plot(stage['x_start'], stage['y_op'], 's', color='orange',
                       markersize=8, label='进料板' if i==0 else "")

        # 标记关键点
        ax.scatter([xD, xF, xW], [xD, xF, xW], color='blue', s=40, zorder=5)
        ax.scatter([x_intersection], [y_intersection], color='red', s=40, zorder=5)

        ax.set_xlabel('液相组成 x', fontsize=11)
        ax.set_ylabel('气相组成 y', fontsize=11)
        ax.set_title(f'McCabe-Thiele阶梯图 (理论板数: {len(stages)})', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper left', fontsize=8)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

    def _plot_composition_profile(self, ax, stages, xF):
        """绘制组成分布图"""
        if not stages:
            return

        stage_numbers = [s['stage'] for s in stages]
        x_compositions = [s['x_eq'] for s in stages]
        y_compositions = [s['y_op'] for s in stages]

        ax.plot(stage_numbers, x_compositions, 'b-o', linewidth=2, markersize=4, label='液相组成 x')
        ax.plot(stage_numbers, y_compositions, 'r-s', linewidth=2, markersize=4, label='气相组成 y')

        # 标记进料板位置
        feed_stage_idx = next((i for i, s in enumerate(stages) if s.get('is_feed', False)), -1)
        if feed_stage_idx >= 0:
            feed_stage = stage_numbers[feed_stage_idx]
            ax.axvline(x=feed_stage, color='g', linestyle='--', alpha=0.7,
                      label=f'进料板 (第{feed_stage}块)')
            ax.axhline(y=xF, color='orange', linestyle=':', alpha=0.7, label=f'进料组成 xF={xF:.3f}')

        ax.set_xlabel('理论板序号', fontsize=11)
        ax.set_ylabel('组成', fontsize=11)
        ax.set_title('塔内组成分布', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right', fontsize=8)
        ax.set_xlim(1, len(stages))

    def get_design_summary(self):
        """获取设计汇总信息"""
        if not all([self.xD, self.xF, self.xW, self.q, self.R_actual, self.R_internal]):
            return "设计参数未完全设置"

        summary = f"""
        精馏塔设计参数汇总:
        ===================================
        1. 组成参数:
          塔顶组成 (xD): {self.xD:.4f}
          进料组成 (xF): {self.xF:.4f}
          塔釜组成 (xW): {self.xW:.4f}
        
        2. 操作参数:
          进料热状况 (q): {self.q:.3f}
          最小回流比 (R_min): {self.R_min:.3f}
          回流比倍数: {self.R_multiplier:.2f}
          外回流比 (R): {self.R_actual:.3f}
          内回流比 (R'): {self.R_internal:.3f}
          回流液温度: {self.reflux_temp if self.reflux_temp else '未指定'} °C
        
        3. 塔板参数:
          理论板总数: {self.total_stages if self.total_stages else '未计算'}
          进料板位置: 第{self.feed_stage if self.feed_stage else '未计算'}块
          塔板效率: {self.ET*100:.1f}%
        
        4. 温度参数:
          塔顶温度: {float(self.t_from_x(self.xD)):.1f} °C
          塔釜温度: {float(self.t_from_x(self.xW)):.1f} °C
        ===================================
        """
        return summary


# 示例用法
if __name__ == "__main__":
    # 苯-甲苯相平衡数据
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
               0.6: [0.791, 88.70],
               0.65: [0.825, 87.63],
               0.7: [0.857, 86.52],
               0.75: [0.885, 85.44],
               0.8: [0.912, 84.40],
               0.85: [0.936, 83.33],
               0.9: [0.959, 82.25],
               0.95: [0.98, 81.11],
               0.97: [0.988, 80.66],
               0.99: [0.9961, 80.21],
               1.0: [1.0, 80.01]}

    # 物性数据（苯-甲苯）
    Cp_benzene = {0: 1.570, 20: 1.716, 40: 1.767, 60: 1.828, 80: 1.881, 100: 1.953}
    Cp_toluene = {0: 1.630, 20: 1.681, 40: 1.757, 60: 1.834, 80: 1.902, 100: 1.970, 120: 2.073}

    Cp_data = {
        'benzene': Cp_benzene,
        'toluene': Cp_toluene
    }

    r_data = {
        'benzene': 393.9,  # kJ/kg
        'toluene': 363.0  # kJ/kg
    }

    M_components = {
        'benzene': 78.11,  # g/mol
        'toluene': 92.14  # g/mol
    }

    # 创建增强版计算器
    calculator = EnhancedRefluxCalculator(
        dic_txy,
        Cp_data=Cp_data,
        r_data=r_data,
        M_components=M_components
    )

    # 设计参数（应与主程序一致）
    xD = 0.95  # 塔顶组成
    xF = 0.25  # 进料组成
    xW = 0.05  # 塔釜组成
    q = 1.266  # 进料热状况参数（冷液进料）
    reflux_temp = 60  # 回流液温度 (°C)
    R_multiplier = 1.5  # 回流比倍数

    print("=" * 70)
    print("考虑过冷回流的精馏塔设计分析")
    print("=" * 70)

    # 计算最小回流比
    R_min, intersection = calculator.calculate_R_min(xD, xF, q)
    print(f"最小回流比 R_min = {R_min:.4f}")
    print(f"q线与平衡线交点: x={intersection[0]:.4f}, y={intersection[1]:.4f}")

    # 计算塔顶温度
    top_temp = float(calculator.t_from_x(xD))
    print(f"塔顶温度: {top_temp:.1f} °C")

    # 计算考虑过冷回流的实际回流比
    R_actual, R_internal, correction = calculator.calculate_actual_reflux_ratio(
        R_min, R_multiplier, reflux_temp, top_temp, xD
    )

    print(f"\n回流比计算结果:")
    print(f"理论外回流比 (R_min×{R_multiplier}): {R_actual:.4f}")
    print(f"实际内回流比 (考虑过冷): {R_internal:.4f}")
    print(f"过冷回流修正系数: {correction:.4f}")
    print(f"过冷度: {top_temp - reflux_temp:.1f} °C")

    # 执行McCabe-Thiele分析
    mccabe_results = calculator.perform_mccabe_thiele_analysis(xD, xF, xW, q, R_internal)

    print(f"\nMcCabe-Thiele分析结果:")
    print(f"理论板总数: {mccabe_results['total_stages']}")
    print(f"进料板位置: 第{mccabe_results['feed_stage']}块理论板")
    print(f"操作线交点: x={mccabe_results['intersection'][0]:.4f}, y={mccabe_results['intersection'][1]:.4f}")

    # 设置设计参数
    calculator.set_design_parameters(
        xD=xD, xF=xF, xW=xW, q=q, R_min=R_min, R_multiplier=R_multiplier,
        R_actual=R_actual, R_internal=R_internal,
        reflux_temp=reflux_temp, feed_stage=mccabe_results['feed_stage'],
        total_stages=mccabe_results['total_stages']
    )

    # 输出设计汇总
    print(calculator.get_design_summary())

    # 绘制分析图
    fig, axes, R_min_calc, R_actual_calc, R_internal_calc, mccabe_results = calculator.plot_enhanced_diagram(
        xD, xF, xW, q, R_multiplier, reflux_temp, R_internal, save_path="enhanced_reflux_analysis.png"
    )

    print("\n分析图表已生成")
    plt.show()
