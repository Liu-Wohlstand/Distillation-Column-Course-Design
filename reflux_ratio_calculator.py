"""
增强版回流比计算模块（无绘图，支持 q=1）
"""

import numpy as np
from scipy.interpolate import interp1d


class EnhancedRefluxCalculator:
    """考虑过冷回流的增强版回流比计算器（无绘图）"""

    def __init__(self, equilibrium_data, Cp_data=None, r_data=None, M_components=None):
        # 相平衡数据
        self.equilibrium_data = equilibrium_data
        self.x_data = list(equilibrium_data.keys())
        self.y_data = [equilibrium_data[x][0] for x in self.x_data]
        self.t_data = [equilibrium_data[x][1] for x in self.x_data]

        sorted_idx = np.argsort(self.x_data)
        self.x_data = np.array([self.x_data[i] for i in sorted_idx])
        self.y_data = np.array([self.y_data[i] for i in sorted_idx])
        self.t_data = np.array([self.t_data[i] for i in sorted_idx])

        # 插值函数（限制外推边界，避免负值）
        self.y_from_x = interp1d(self.x_data, self.y_data, kind='cubic',
                                 bounds_error=False, fill_value=(self.y_data[0], self.y_data[-1]))
        self.x_from_y = interp1d(self.y_data, self.x_data, kind='cubic',
                                 bounds_error=False, fill_value=(self.x_data[0], self.x_data[-1]))
        self.t_from_x = interp1d(self.x_data, self.t_data, kind='linear',
                                 bounds_error=False, fill_value=(self.t_data[0], self.t_data[-1]))

        self.Cp_data = Cp_data or {}
        self.r_data = r_data or {}
        self.M_components = M_components or {}

        # 设计参数
        self.xD = self.xF = self.xW = self.q = self.R_min = None
        self.R_multiplier = self.R_actual = self.R_internal = self.reflux_temp = None
        self.feed_stage = self.total_stages = self.ET = None
        self.stage_positions = []

    def set_design_parameters(self, xD, xF, xW, q, R_min, R_multiplier,
                             R_actual, R_internal, reflux_temp,
                             feed_stage, total_stages, ET=None, stage_positions=None):
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
        self.ET = ET or 0.5
        self.stage_positions = stage_positions or []

    def calculate_R_min(self, xD, xF, q=1.0, method='graphical'):
        """计算最小回流比（正确处理 q=1）"""
        # 强制将非常接近 1 的 q 设为精确 1
        if abs(q - 1) < 1e-6:
            q = 1.0

        # 获取 q 线与平衡线的交点
        x_eq, y_eq = self._get_intersection(xF, q)

        if abs(q - 1) < 1e-6:
            # 泡点进料：直接使用解析公式
            y_eq = float(self.y_from_x(xF))
            R_min = (xD - y_eq) / (y_eq - xF)
        else:
            if x_eq == xD:
                R_min = float('inf')
            else:
                R_min = (xD - y_eq) / (y_eq - x_eq)

        # 最终检查
        if R_min < 0:
            print(f"警告：最小回流比 R_min = {R_min:.4f} 仍为负，请检查进料条件或相平衡数据。")
        return R_min, (x_eq, y_eq)

    def _get_intersection(self, xF, q):
        """获取 q 线与平衡线的交点（q=1 时直接返回）"""
        if abs(q - 1) < 1e-6:
            # q=1 泡点进料：垂直线 x = xF
            y_eq = float(self.y_from_x(xF))
            return xF, y_eq

        if abs(q) < 1e-6:
            # q=0 饱和蒸汽进料：水平线 y = xF
            try:
                x_eq = float(self.x_from_y(xF))
                return x_eq, xF
            except:
                return self._numerical_intersection(xF, q)

        return self._numerical_intersection(xF, q)

    def _numerical_intersection(self, xF, q):
        """数值求解交点（用于 q≠1 且 q≠0）"""
        def q_line_y(x):
            return q / (q - 1) * x - xF / (q - 1)

        x_range = np.linspace(0, 1, 1000)
        best_x, best_y = 0, 0
        min_diff = float('inf')
        for x in x_range:
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

    def calculate_actual_reflux_ratio(self, R_min, R_multiplier=1.2,
                                      reflux_temp=None, top_temp=None,
                                      x_top=None):
        """计算考虑过冷回流的实际回流比"""
        R_theoretical = R_min * R_multiplier
        if reflux_temp is None or top_temp is None:
            if top_temp is None and x_top is not None:
                top_temp = float(self.t_from_x(x_top))
            else:
                return R_theoretical, R_theoretical, 1.0

        correction = self._subcooling_correction(reflux_temp, top_temp, x_top)
        R_internal = R_theoretical * correction
        return R_theoretical, R_internal, correction

    def _subcooling_correction(self, T_reflux, T_top, x_top):
        """过冷回流修正系数"""
        if not self.Cp_data or not self.r_data or not self.M_components:
            print("警告：缺少物性数据，修正系数=1.0")
            return 1.0

        T_avg = (T_reflux + T_top) / 2
        # 苯
        Cp_benzene = self._interpolate_property(T_avg, self.Cp_data.get('benzene', {}))
        M_benzene = self.M_components.get('benzene', 78.11)
        Cp_benz_mol = Cp_benzene * M_benzene
        # 甲苯
        Cp_toluene = self._interpolate_property(T_avg, self.Cp_data.get('toluene', {}))
        M_toluene = self.M_components.get('toluene', 92.14)
        Cp_tolu_mol = Cp_toluene * M_toluene

        Cp_mix = x_top * Cp_benz_mol + (1 - x_top) * Cp_tolu_mol

        r_benz = self.r_data.get('benzene', 393.9) * M_benzene
        r_tolu = self.r_data.get('toluene', 363.0) * M_toluene
        r_mix = x_top * r_benz + (1 - x_top) * r_tolu

        delta_T = T_top - T_reflux
        corr = 1 + (Cp_mix * delta_T) / r_mix
        print(f"\n过冷回流修正: ΔT={delta_T:.1f}K, Cp={Cp_mix:.1f}, ΔHv={r_mix:.0f}, 修正={corr:.3f}")
        return corr

    def _interpolate_property(self, T, prop_dict):
        if not prop_dict:
            return 1.8  # 默认值
        temps = sorted(prop_dict.keys())
        vals = [prop_dict[t] for t in temps]
        if T <= temps[0]:
            return vals[0]
        if T >= temps[-1]:
            return vals[-1]
        return np.interp(T, temps, vals)

    def get_design_summary(self):
        if None in [self.xD, self.xF, self.xW, self.q, self.R_actual]:
            return "设计参数未完全设置"
        return f"""
        精馏塔设计参数汇总:
        ===================================
        塔顶 xD={self.xD:.4f}  进料 xF={self.xF:.4f}  塔釜 xW={self.xW:.4f}
        q={self.q:.3f}  R_min={self.R_min:.3f}  R外={self.R_actual:.3f}  R内={self.R_internal:.3f}
        理论板数={self.total_stages}  进料板={self.feed_stage}
        塔板效率={self.ET*100:.1f}%
        塔顶温度={self.t_from_x(self.xD):.1f}°C  塔釜温度={self.t_from_x(self.xW):.1f}°C
        ===================================
        """
