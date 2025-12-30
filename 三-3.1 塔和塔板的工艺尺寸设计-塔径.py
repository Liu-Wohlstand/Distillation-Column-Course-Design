import json
import math
import os
import random
from datetime import datetime


class TowerDiameterApproximateOptimizedCalculator:
    """塔径近似计算器（使用史密斯关联图经验公式）- 带自动优化功能"""

    STANDARD_DIAMETERS = [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                          1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]

    # 空塔气速的经验范围 (m/s)
    U_ALLOWED_RANGE = (0.3, 1.5)
    # 泛点率的合理范围 (%)
    FLOODING_ALLOWED_RANGE = (60, 80)

    # 可调整参数的范围
    HT_RANGE = (0.3, 0.6)  # 板间距范围，m
    HL_RANGE = (0.05, 0.08)  # 清液层高度范围，m
    FLOODING_FRACTION_RANGE = (0.6, 0.8)  # 设计泛点率范围

    def __init__(self, json_file_path='tower_properties_results.json'):
        self.json_file_path = json_file_path
        self.tower_data = None
        self.load_tower_data()

        # 设计参数（初始值）
        self.HT_rectifying = 0.45  # 精馏段板间距，m
        self.HT_stripping = 0.45  # 提馏段板间距，m
        self.hL_rectifying = 0.06  # 精馏段清液层高度，m
        self.hL_stripping = 0.06  # 提馏段清液层高度，m
        self.design_flooding_fraction = 0.7  # 设计泛点率

        # 存储最终计算结果
        self.final_results = None

    def load_tower_data(self):
        """从JSON文件加载塔工艺数据"""
        try:
            with open(self.json_file_path, 'r', encoding='utf-8') as f:
                self.tower_data = json.load(f)
            print(f"✓ 成功从 '{self.json_file_path}' 加载塔工艺数据\n")
        except FileNotFoundError:
            print(f"✗ 错误: 文件 '{self.json_file_path}' 未找到")
            self.tower_data = None

    def calculate_F_LV(self, section_name):
        """计算流动参数F_LV"""
        if not self.tower_data or section_name not in self.tower_data:
            return None

        data = self.tower_data[section_name]

        # 使用质量流量计算
        L_mass = data['L_mass_kgh']  # kg/h
        V_mass = data['V_mass_kgh']  # kg/h
        ρ_L = data['rho_liquid_kgm3']  # kg/m³
        ρ_V = data['rho_vapor_kgm3']  # kg/m³

        mass_ratio = L_mass / V_mass
        density_ratio_sqrt = math.sqrt(ρ_V / ρ_L)
        F_LV_mass = mass_ratio * density_ratio_sqrt

        return F_LV_mass

    def approximate_C20_from_smith_chart(self, F_LV, HT_minus_hL):
        """
        使用史密斯关联图的经验公式计算C20
        公式: C20 = 0.0162 - 0.0648*X + 0.181*Y + 0.0162*X^2 - 0.139*X*Y + 0.185*Y^2
        其中: X = F_LV, Y = HT_minus_hL
        """
        X = F_LV  # 流动参数
        Y = HT_minus_hL  # HT-hL值

        # 使用提供的公式计算C20
        C20 = (0.0162 - 0.0648 * X + 0.181 * Y +
               0.0162 * X ** 2 - 0.139 * X * Y + 0.185 * Y ** 2)

        # 确保在合理范围内（筛板塔C20通常在0.03-0.12之间）
        if C20 < 0.03:
            C20 = 0.03
        elif C20 > 0.15:
            C20 = 0.15

        return C20

    def check_parameters(self, u_actual, flooding_percent):
        """检查空塔气速和泛点率是否在允许范围内"""
        velocity_ok = self.U_ALLOWED_RANGE[0] <= u_actual <= self.U_ALLOWED_RANGE[1]
        flooding_ok = self.FLOODING_ALLOWED_RANGE[0] <= flooding_percent <= self.FLOODING_ALLOWED_RANGE[1]

        return velocity_ok and flooding_ok, velocity_ok, flooding_ok

    def calculate_section_diameter_with_params(self, section_name, HT, hL, design_fraction, verbose=True):
        """使用指定参数计算塔段塔径"""
        if not self.tower_data or section_name not in self.tower_data:
            return None

        data = self.tower_data[section_name]

        # 1. 计算F_LV
        F_LV = self.calculate_F_LV(section_name)

        # 2. 计算HT-hL
        HT_minus_hL = HT - hL

        # 3. 使用经验公式计算C20
        C20 = self.approximate_C20_from_smith_chart(F_LV, HT_minus_hL)

        # 4. 表面张力校正
        σ = data['sigma_mNm']
        C = C20 * (σ / 20) ** 0.2

        # 5. 计算液泛气速
        ρ_L = data['rho_liquid_kgm3']
        ρ_V = data['rho_vapor_kgm3']
        u_f = C * math.sqrt((ρ_L - ρ_V) / ρ_V)

        # 6. 设计空塔气速
        u_design = design_fraction * u_f

        # 7. 计算塔截面积
        V_m3h = data['V_vol_m3h']
        V_m3s = V_m3h / 3600
        A_t = V_m3s / u_design

        # 8. 计算理论塔径
        D_calc = math.sqrt(4 * A_t / math.pi)

        # 9. 塔径圆整
        D_rounded = min(self.STANDARD_DIAMETERS, key=lambda x: abs(x - D_calc))

        # 10. 计算实际参数
        A_actual = math.pi * D_rounded ** 2 / 4
        u_actual = V_m3s / A_actual
        flooding_percent = (u_actual / u_f) * 100

        # 11. 检查参数
        all_ok, velocity_ok, flooding_ok = self.check_parameters(u_actual, flooding_percent)

        results = {
            'section': section_name.replace('结果', ''),
            'F_LV': round(F_LV, 4),
            'HT': round(HT, 2),
            'hL': round(hL, 3),
            'HT_minus_hL': round(HT_minus_hL, 3),
            'C20_approximate': round(C20, 4),
            'C_corrected': round(C, 4),
            'u_f': round(u_f, 4),
            '设计泛点率': design_fraction,
            'u_design': round(u_design, 4),
            'u_actual': round(u_actual, 4),
            'flooding_percent': round(flooding_percent, 1),
            'D_calc': round(D_calc, 4),
            'D_rounded': round(D_rounded, 3),
            'A_actual': round(A_actual, 4),
            'V_flow_m3s': round(V_m3s, 6),
            'parameters_ok': all_ok,
            'velocity_ok': velocity_ok,
            'flooding_ok': flooding_ok
        }

        if verbose:
            self.print_section_results_detailed(results)

        return results

    def print_section_results_detailed(self, results):
        """详细打印塔段计算结果"""
        print(f"\n{'=' * 80}")
        print(f"            {results['section']}塔径详细计算过程")
        print(f"{'=' * 80}")

        print(f"1. 流动参数 F_LV = {results['F_LV']:.4f}")
        print(f"2. 板间距 HT = {results['HT']:.2f} m, 清液层高度 hL = {results['hL']:.3f} m")
        print(f"   HT - hL = {results['HT_minus_hL']:.3f} m")
        print(f"3. 经验公式计算C20 = {results['C20_approximate']:.4f}")
        print(f"4. 表面张力校正后 C = {results['C_corrected']:.4f}")
        print(f"5. 液泛气速 u_f = {results['u_f']:.4f} m/s")
        print(f"6. 设计泛点率 = {results['设计泛点率']:.2f}")
        print(f"   设计空塔气速 u_design = {results['u_design']:.4f} m/s")
        print(f"7. 气体体积流量 V = {results['V_flow_m3s']:.6f} m³/s")
        print(f"8. 理论塔截面积 A_t = {results['A_actual']:.4f} m²")
        print(f"9. 理论塔径 D_calc = {results['D_calc']:.4f} m")
        print(f"10. 圆整塔径 D_rounded = {results['D_rounded']:.3f} m")
        print(f"11. 实际塔截面积 A_actual = {results['A_actual']:.4f} m²")
        print(f"12. 实际空塔气速 u_actual = {results['u_actual']:.4f} m/s")
        print(f"13. 实际泛点率 = {results['flooding_percent']:.1f}%")

        print(f"\n{'=' * 60}")
        print("14. 参数校验:")
        print(f"{'=' * 60}")

        # 空塔气速校验
        if results['velocity_ok']:
            velocity_status = "✅ 合理"
            velocity_suggestion = f"气速在合理范围内"
        else:
            if results['u_actual'] < self.U_ALLOWED_RANGE[0]:
                velocity_status = "❌ 过低"
                velocity_suggestion = f"建议: 考虑减小塔径以提高气速"
            else:
                velocity_status = "❌ 过高"
                velocity_suggestion = f"建议: 考虑增大塔径以降低气速"

        print(f"   空塔气速校验:")
        print(f"     实际值: {results['u_actual']:.4f} m/s")
        print(f"     允许范围: {self.U_ALLOWED_RANGE[0]:.1f} - {self.U_ALLOWED_RANGE[1]:.1f} m/s")
        print(f"     结果: {velocity_status}")
        print(f"     {velocity_suggestion}")

        # 泛点率校验
        if results['flooding_ok']:
            flooding_status = "✅ 合理"
            flooding_suggestion = f"泛点率在合理操作范围内"
        else:
            if results['flooding_percent'] < self.FLOODING_ALLOWED_RANGE[0]:
                flooding_status = "⚠️ 偏低"
                flooding_suggestion = f"建议: 可以考虑适当减小塔径以提高设备利用率"
            else:
                flooding_status = "⚠️ 偏高"
                flooding_suggestion = f"建议: 考虑增大塔径以降低泛点率"

        print(f"\n   泛点率校验:")
        print(f"     实际值: {results['flooding_percent']:.1f}%")
        print(f"     合理范围: {self.FLOODING_ALLOWED_RANGE[0]} - {self.FLOODING_ALLOWED_RANGE[1]}%")
        print(f"     结果: {flooding_status}")
        print(f"     {flooding_suggestion}")

        print(f"\n   总体校验结果: {'✅ 所有参数合理' if results['parameters_ok'] else '❌ 参数需要调整'}")

    def optimize_section_parameters(self, section_name, max_attempts=50):
        """自动优化单个塔段的参数"""
        print(f"\n{'=' * 80}")
        print(f"【自动优化】{section_name.replace('结果', '')} - 尝试寻找满足条件的参数组合")
        print(f"{'=' * 80}")

        best_result = None
        best_score = -float('inf')

        for attempt in range(1, max_attempts + 1):
            # 随机生成参数组合
            HT = random.uniform(self.HT_RANGE[0], self.HT_RANGE[1])
            hL = random.uniform(self.HL_RANGE[0], self.HL_RANGE[1])
            design_fraction = random.uniform(self.FLOODING_FRACTION_RANGE[0], self.FLOODING_FRACTION_RANGE[1])

            # 计算塔径
            result = self.calculate_section_diameter_with_params(section_name, HT, hL, design_fraction, verbose=False)

            if result['parameters_ok']:
                # 计算得分（越接近目标越好）
                velocity_mid = (self.U_ALLOWED_RANGE[0] + self.U_ALLOWED_RANGE[1]) / 2
                flooding_mid = (self.FLOODING_ALLOWED_RANGE[0] + self.FLOODING_ALLOWED_RANGE[1]) / 2

                velocity_score = 1 / (abs(result['u_actual'] - velocity_mid) + 0.1)
                flooding_score = 1 / (abs(result['flooding_percent'] - flooding_mid) + 0.1)
                total_score = velocity_score + flooding_score

                result['score'] = total_score
                result['attempt'] = attempt

                if total_score > best_score:
                    best_score = total_score
                    best_result = result

                print(f"  尝试 {attempt:2d}: HT={HT:.2f}m, hL={hL:.3f}m, 泛点率={design_fraction:.2f}, "
                      f"D={result['D_rounded']:.3f}m, u={result['u_actual']:.4f}m/s, "
                      f"泛点={result['flooding_percent']:.1f}% ✓")
            else:
                print(f"  尝试 {attempt:2d}: 参数不满足要求 ✗")

        if best_result:
            print(f"\n✓ 找到最优参数组合（第{best_result['attempt']}次尝试）:")
            print(f"  板间距 HT = {best_result['HT']} m")
            print(f"  清液层高度 hL = {best_result['hL']} m")
            print(f"  设计泛点率 = {best_result['设计泛点率']:.2f}")
            print(f"  塔径 D = {best_result['D_rounded']} m")
            print(f"  空塔气速 u = {best_result['u_actual']:.4f} m/s")
            print(f"  实际泛点率 = {best_result['flooding_percent']}%")

            # 更新参数
            if "精馏段" in section_name:
                self.HT_rectifying = best_result['HT']
                self.hL_rectifying = best_result['hL']
            else:
                self.HT_stripping = best_result['HT']
                self.hL_stripping = best_result['hL']

            # 详细打印结果
            self.print_section_results_detailed(best_result)

            return best_result
        else:
            print(f"\n✗ 经过 {max_attempts} 次尝试仍未找到满足条件的参数组合")
            return None

    def calculate_final_diameter_parameters(self, final_D, results_rect, results_strip):
        """计算最终塔径下的详细参数（类似精馏段和提馏段的详细输出）"""
        print(f"\n{'=' * 80}")
        print("最终塔径下详细参数核算")
        print(f"{'=' * 80}")

        A_final = math.pi * final_D ** 2 / 4
        print(f"最终塔径: {final_D:.3f} m")
        print(f"最终塔截面积: A = π × {final_D:.3f}² / 4 = {A_final:.4f} m²")

        final_params = {}

        for section_name, original_results, section_label in [
            ('精馏段结果', results_rect, '精馏段'),
            ('提馏段结果', results_strip, '提馏段')
        ]:
            print(f"\n{'=' * 60}")
            print(f"{section_label}在最终塔径 {final_D:.3f} m 下的详细核算:")
            print(f"{'=' * 60}")

            data = self.tower_data[section_name]
            V_m3h = data['V_vol_m3h']
            V_m3s = V_m3h / 3600

            print(f"1. 气体体积流量:")
            print(f"   V = {V_m3h:.2f} m³/h")
            print(f"     = {V_m3h:.2f} / 3600")
            print(f"     = {V_m3s:.6f} m³/s")

            u_actual_final = V_m3s / A_final
            print(f"2. 空塔气速:")
            print(f"   u = V / A = {V_m3s:.6f} / {A_final:.4f}")
            print(f"     = {u_actual_final:.4f} m/s")

            # 重新计算HT-hL和C20（使用优化后的参数）
            if section_label == '精馏段':
                HT = self.HT_rectifying
                hL = self.hL_rectifying
            else:
                HT = self.HT_stripping
                hL = self.hL_stripping

            HT_minus_hL = HT - hL
            F_LV = self.calculate_F_LV(section_name)
            C20 = self.approximate_C20_from_smith_chart(F_LV, HT_minus_hL)
            C = C20 * (data['sigma_mNm'] / 20) ** 0.2

            print(f"3. 重新计算HT-hL和C值:")
            print(f"   HT = {HT:.2f} m, hL = {hL:.3f} m")
            print(f"   HT - hL = {HT_minus_hL:.3f} m")
            print(f"   F_LV = {F_LV:.4f}")
            print(f"   经验公式计算C20 = {C20:.4f}")
            print(f"   表面张力校正后 C = {C:.4f}")

            ρ_L = data['rho_liquid_kgm3']
            ρ_V = data['rho_vapor_kgm3']
            u_f = C * math.sqrt((ρ_L - ρ_V) / ρ_V)
            print(f"4. 液泛气速:")
            print(f"   u_f = {C:.4f} × √[({ρ_L:.1f} - {ρ_V:.4f})/{ρ_V:.4f}]")
            print(f"       = {u_f:.4f} m/s")

            flooding_percent_final = (u_actual_final / u_f) * 100
            print(f"5. 实际泛点率:")
            print(f"   泛点率 = ({u_actual_final:.4f} / {u_f:.4f}) × 100%")
            print(f"          = {flooding_percent_final:.2f}%")

            # 参数校验
            print(f"\n{'=' * 60}")
            print(f"参数校验:")
            print(f"{'=' * 60}")

            all_ok, velocity_ok, flooding_ok = self.check_parameters(u_actual_final, flooding_percent_final)

            # 空塔气速校验
            if velocity_ok:
                velocity_status = "✅ 合理"
                velocity_suggestion = f"气速在合理范围内"
            else:
                if u_actual_final < self.U_ALLOWED_RANGE[0]:
                    velocity_status = "❌ 过低"
                    velocity_suggestion = f"建议: 考虑减小塔径以提高气速"
                else:
                    velocity_status = "❌ 过高"
                    velocity_suggestion = f"建议: 考虑增大塔径以降低气速"

            print(f"   空塔气速校验:")
            print(f"     实际值: {u_actual_final:.4f} m/s")
            print(f"     允许范围: {self.U_ALLOWED_RANGE[0]:.1f} - {self.U_ALLOWED_RANGE[1]:.1f} m/s")
            print(f"     结果: {velocity_status}")
            print(f"     {velocity_suggestion}")

            # 泛点率校验
            if flooding_ok:
                flooding_status = "✅ 合理"
                flooding_suggestion = f"泛点率在合理操作范围内"
            else:
                if flooding_percent_final < self.FLOODING_ALLOWED_RANGE[0]:
                    flooding_status = "⚠️ 偏低"
                    flooding_suggestion = f"建议: 可以考虑适当减小塔径以提高设备利用率"
                else:
                    flooding_status = "⚠️ 偏高"
                    flooding_suggestion = f"建议: 考虑增大塔径以降低泛点率"

            print(f"\n   泛点率校验:")
            print(f"     实际值: {flooding_percent_final:.2f}%")
            print(f"     合理范围: {self.FLOODING_ALLOWED_RANGE[0]} - {self.FLOODING_ALLOWED_RANGE[1]}%")
            print(f"     结果: {flooding_status}")
            print(f"     {flooding_suggestion}")

            print(f"\n   总体校验结果: {'✅ 所有参数合理' if all_ok else '❌ 参数需要调整'}")

            # 存储参数
            final_params[section_label] = {
                'u_actual_final': u_actual_final,
                'flooding_percent_final': flooding_percent_final,
                'velocity_ok': velocity_ok,
                'flooding_ok': flooding_ok,
                'all_ok': all_ok,
                'HT': HT,
                'hL': hL,
                'HT_minus_hL': HT_minus_hL,
                'C20': C20,
                'C': C,
                'u_f': u_f
            }

        return final_params

    def optimize_final_diameter_parameters(self, results_rect, results_strip, max_attempts=100):
        """针对最终塔径优化参数，确保所有参数都在合理范围内"""
        print(f"\n{'=' * 80}")
        print("【最终塔径参数优化】启动自动参数调整...")
        print(f"{'=' * 80}")

        final_D = max(results_rect['D_rounded'], results_strip['D_rounded'])
        controlling_section = "精馏段" if final_D == results_rect['D_rounded'] else "提馏段"

        best_solution = None
        best_score = -float('inf')

        for attempt in range(1, max_attempts + 1):
            # 随机生成参数组合
            HT_rect = random.uniform(self.HT_RANGE[0], self.HT_RANGE[1])
            hL_rect = random.uniform(self.HL_RANGE[0], self.HL_RANGE[1])
            HT_strip = random.uniform(self.HT_RANGE[0], self.HT_RANGE[1])
            hL_strip = random.uniform(self.HL_RANGE[0], self.HL_RANGE[1])

            # 分别计算两段在最终塔径下的参数
            all_ok = True
            score = 0

            # 计算精馏段参数
            data_rect = self.tower_data['精馏段结果']
            HT_minus_hL_rect = HT_rect - hL_rect
            F_LV_rect = self.calculate_F_LV('精馏段结果')
            C20_rect = self.approximate_C20_from_smith_chart(F_LV_rect, HT_minus_hL_rect)
            C_rect = C20_rect * (data_rect['sigma_mNm'] / 20) ** 0.2
            ρ_L_rect = data_rect['rho_liquid_kgm3']
            ρ_V_rect = data_rect['rho_vapor_kgm3']
            u_f_rect = C_rect * math.sqrt((ρ_L_rect - ρ_V_rect) / ρ_V_rect)

            # 最终塔径下的实际参数
            A_final = math.pi * final_D ** 2 / 4
            V_m3h_rect = data_rect['V_vol_m3h']
            V_m3s_rect = V_m3h_rect / 3600
            u_actual_rect_final = V_m3s_rect / A_final
            flooding_percent_rect_final = (u_actual_rect_final / u_f_rect) * 100

            # 计算提馏段参数
            data_strip = self.tower_data['提馏段结果']
            HT_minus_hL_strip = HT_strip - hL_strip
            F_LV_strip = self.calculate_F_LV('提馏段结果')
            C20_strip = self.approximate_C20_from_smith_chart(F_LV_strip, HT_minus_hL_strip)
            C_strip = C20_strip * (data_strip['sigma_mNm'] / 20) ** 0.2
            ρ_L_strip = data_strip['rho_liquid_kgm3']
            ρ_V_strip = data_strip['rho_vapor_kgm3']
            u_f_strip = C_strip * math.sqrt((ρ_L_strip - ρ_V_strip) / ρ_V_strip)

            V_m3h_strip = data_strip['V_vol_m3h']
            V_m3s_strip = V_m3h_strip / 3600
            u_actual_strip_final = V_m3s_strip / A_final
            flooding_percent_strip_final = (u_actual_strip_final / u_f_strip) * 100

            # 检查参数是否合理
            rect_ok, rect_v_ok, rect_f_ok = self.check_parameters(u_actual_rect_final, flooding_percent_rect_final)
            strip_ok, strip_v_ok, strip_f_ok = self.check_parameters(u_actual_strip_final, flooding_percent_strip_final)

            all_ok = rect_ok and strip_ok

            if all_ok:
                # 计算得分（参数越接近范围中间值得分越高）
                velocity_mid = (self.U_ALLOWED_RANGE[0] + self.U_ALLOWED_RANGE[1]) / 2
                flooding_mid = (self.FLOODING_ALLOWED_RANGE[0] + self.FLOODING_ALLOWED_RANGE[1]) / 2

                velocity_score_rect = 1 / (abs(u_actual_rect_final - velocity_mid) + 0.1)
                flooding_score_rect = 1 / (abs(flooding_percent_rect_final - flooding_mid) + 0.1)
                velocity_score_strip = 1 / (abs(u_actual_strip_final - velocity_mid) + 0.1)
                flooding_score_strip = 1 / (abs(flooding_percent_strip_final - flooding_mid) + 0.1)

                total_score = velocity_score_rect + flooding_score_rect + velocity_score_strip + flooding_score_strip

                solution = {
                    'HT_rect': HT_rect,
                    'hL_rect': hL_rect,
                    'HT_strip': HT_strip,
                    'hL_strip': hL_strip,
                    'u_actual_rect': u_actual_rect_final,
                    'flooding_rect': flooding_percent_rect_final,
                    'u_actual_strip': u_actual_strip_final,
                    'flooding_strip': flooding_percent_strip_final,
                    'score': total_score,
                    'attempt': attempt
                }

                if total_score > best_score:
                    best_score = total_score
                    best_solution = solution

                '''print(f"  尝试 {attempt:3d}: HT_rect={HT_rect:.2f}m, hL_rect={hL_rect:.3f}m, "
                      f"u_rect={u_actual_rect_final:.4f}m/s, f_rect={flooding_percent_rect_final:.1f}%, "
                      f"HT_strip={HT_strip:.2f}m, hL_strip={hL_strip:.3f}m, "
                      f"u_strip={u_actual_strip_final:.4f}m/s, f_strip={flooding_percent_strip_final:.1f}% ✓")'''
            else:
                None
                # print(f"  尝试 {attempt:3d}: 参数不满足要求 ✗")

        if best_solution:
            print(f"\n✓ 找到最优参数组合（第{best_solution['attempt']}次尝试）:")
            print(f"  精馏段: HT={best_solution['HT_rect']:.2f}m, hL={best_solution['hL_rect']:.3f}m")
            print(
                f"          空塔气速={best_solution['u_actual_rect']:.4f}m/s, 泛点率={best_solution['flooding_rect']:.1f}%")
            print(f"  提馏段: HT={best_solution['HT_strip']:.2f}m, hL={best_solution['hL_strip']:.3f}m")
            print(
                f"          空塔气速={best_solution['u_actual_strip']:.4f}m/s, 泛点率={best_solution['flooding_strip']:.1f}%")
            print(f"  最终塔径: {final_D:.3f} m")

            # 更新参数
            self.HT_rectifying = best_solution['HT_rect']
            self.hL_rectifying = best_solution['hL_rect']
            self.HT_stripping = best_solution['HT_strip']
            self.hL_stripping = best_solution['hL_strip']

            # 存储优化后的结果
            self.final_results = {
                'final_D': final_D,
                'controlling_section': controlling_section,
                'HT_rectifying': best_solution['HT_rect'],
                'hL_rectifying': best_solution['hL_rect'],
                'HT_stripping': best_solution['HT_strip'],
                'hL_stripping': best_solution['hL_strip'],
                'u_actual_rect_final': best_solution['u_actual_rect'],
                'flooding_percent_rect_final': best_solution['flooding_rect'],
                'u_actual_strip_final': best_solution['u_actual_strip'],
                'flooding_percent_strip_final': best_solution['flooding_strip'],
                'optimized': True
            }

            return best_solution
        else:
            print(f"\n✗ 经过 {max_attempts} 次尝试仍未找到满足条件的参数组合")
            return None

    def print_final_check_detailed_calculation(self, final_D, final_params):
        """
        打印最终校核的详细计算过程
        包含所有公式和代入过程，适合抄写到课程设计报告中
        """
        print(f"\n{'=' * 80}")
        print("第六步：最终校核结果详细计算过程（适用于课程设计报告）")
        print(f"{'=' * 80}")

        print("一、设计基本参数")
        print("-" * 40)
        print(f"最终塔径：D = {final_D:.3f} m")
        print(f"精馏段板间距：HT = {self.HT_rectifying:.2f} m")
        print(f"精馏段清液层高度：hL = {self.hL_rectifying:.3f} m")
        print(f"提馏段板间距：HT = {self.HT_stripping:.2f} m")
        print(f"提馏段清液层高度：hL = {self.hL_stripping:.3f} m")
        print(f"设计泛点率：f = {self.design_flooding_fraction:.2f}")

        # 详细计算过程（分精馏段和提馏段）
        for section_name, data_key, section_label in [
            ('精馏段', '精馏段结果', '精馏段'),
            ('提馏段', '提馏段结果', '提馏段')
        ]:
            data = self.tower_data[data_key]
            params = final_params[section_label]

            print(f"\n二、{section_name}详细校核计算过程")
            print("-" * 50)

            # 步骤1：计算实际塔截面积
            print("1. 实际塔截面积计算：")
            print("   公式：A = π × D² / 4")
            print(f"   代入：A = 3.1416 × ({final_D:.3f})² / 4")
            A_actual = math.pi * final_D ** 2 / 4
            print(f"   计算得：A = {A_actual:.4f} m²")

            # 步骤2：计算实际空塔气速
            print("\n2. 实际空塔气速计算：")
            print("   公式：u_实际 = V_s / A")
            V_m3h = data['V_vol_m3h']
            V_m3s = V_m3h / 3600
            print(f"   气体体积流量：V = {V_m3h:.2f} m³/h = {V_m3h:.2f} / 3600 = {V_m3s:.6f} m³/s")
            u_actual = V_m3s / A_actual
            print(f"   代入：u_实际 = {V_m3s:.6f} / {A_actual:.4f}")
            print(f"   计算得：u_实际 = {u_actual:.4f} m/s")

            # 步骤3：计算流动参数F_LV
            print("\n3. 流动参数F_LV计算：")
            print("   公式：F_LV = (L_mass / V_mass) × √(ρ_V / ρ_L)")
            L_mass = data['L_mass_kgh']
            V_mass = data['V_mass_kgh']
            ρ_L = data['rho_liquid_kgm3']
            ρ_V = data['rho_vapor_kgm3']
            print(f"   液相质量流量：L = {L_mass:.1f} kg/h")
            print(f"   气相质量流量：V = {V_mass:.1f} kg/h")
            print(f"   液相密度：ρ_L = {ρ_L:.1f} kg/m³")
            print(f"   气相密度：ρ_V = {ρ_V:.4f} kg/m³")
            F_LV = (L_mass / V_mass) * math.sqrt(ρ_V / ρ_L)
            print(f"   代入：F_LV = ({L_mass:.1f}/{V_mass:.1f}) × √({ρ_V:.4f}/{ρ_L:.1f})")
            print(f"   计算得：F_LV = {F_LV:.4f}")

            # 步骤4：计算HT-hL
            print("\n4. HT - hL计算：")
            HT = self.HT_rectifying if section_name == '精馏段' else self.HT_stripping
            hL = self.hL_rectifying if section_name == '精馏段' else self.hL_stripping
            HT_minus_hL = HT - hL
            print(f"   板间距：HT = {HT:.2f} m")
            print(f"   清液层高度：hL = {hL:.3f} m")
            print(f"   HT - hL = {HT:.2f} - {hL:.3f} = {HT_minus_hL:.3f} m")

            # 步骤5：计算C20（史密斯经验公式）
            print("\n5. C20计算（史密斯经验公式）：")
            print("   公式：C20 = 0.0162 - 0.0648×X + 0.181×Y + 0.0162×X² - 0.139×X×Y + 0.185×Y²")
            print("   其中：X = F_LV, Y = HT - hL")
            print(f"   代入：X = {F_LV:.4f}, Y = {HT_minus_hL:.3f}")

            X = F_LV
            Y = HT_minus_hL
            C20 = (0.0162 - 0.0648 * X + 0.181 * Y +
                   0.0162 * X ** 2 - 0.139 * X * Y + 0.185 * Y ** 2)

            # 显示每一步计算
            print("   分步计算：")
            print(f"   0.0162 = 0.0162")
            print(f"   -0.0648×{X:.4f} = {(-0.0648 * X):.6f}")
            print(f"   0.181×{Y:.3f} = {(0.181 * Y):.6f}")
            print(f"   0.0162×({X:.4f})² = {(0.0162 * X ** 2):.6f}")
            print(f"   -0.139×{X:.4f}×{Y:.3f} = {(-0.139 * X * Y):.6f}")
            print(f"   0.185×({Y:.3f})² = {(0.185 * Y ** 2):.6f}")

            sum_parts = 0.0162 + (-0.0648 * X) + (0.181 * Y) + (0.0162 * X ** 2) + (-0.139 * X * Y) + (0.185 * Y ** 2)
            print(f"   总和 = {sum_parts:.6f}")

            # 限制在合理范围
            if C20 < 0.03:
                C20 = 0.03
                print(f"   计算值小于0.03，取C20 = 0.0300（最小值限制）")
            elif C20 > 0.15:
                C20 = 0.15
                print(f"   计算值大于0.15，取C20 = 0.1500（最大值限制）")
            else:
                print(f"   计算得：C20 = {C20:.4f}")

            # 步骤6：表面张力校正
            print("\n6. 表面张力校正计算：")
            print("   公式：C = C20 × (σ/20)^0.2")
            σ = data['sigma_mNm']
            print(f"   表面张力：σ = {σ:.2f} mN/m")
            C = C20 * (σ / 20) ** 0.2
            print(f"   代入：C = {C20:.4f} × ({σ:.2f}/20)^0.2")
            print(f"   计算得：C = {C:.4f}")

            # 步骤7：计算液泛气速
            print("\n7. 液泛气速计算：")
            print("   公式：u_f = C × √[(ρ_L - ρ_V) / ρ_V]")
            u_f = C * math.sqrt((ρ_L - ρ_V) / ρ_V)
            print(f"   代入：u_f = {C:.4f} × √[({ρ_L:.1f} - {ρ_V:.4f}) / {ρ_V:.4f}]")
            print(f"   计算得：u_f = {u_f:.4f} m/s")

            # 步骤8：计算实际泛点率
            print("\n8. 实际泛点率计算：")
            print("   公式：实际泛点率 = (u_实际 / u_f) × 100%")
            flooding_percent = (u_actual / u_f) * 100
            print(f"   代入：实际泛点率 = ({u_actual:.4f} / {u_f:.4f}) × 100%")
            print(f"   计算得：实际泛点率 = {flooding_percent:.2f}%")

            # 步骤9：参数校核
            print("\n9. 参数校核：")
            print("   (1) 空塔气速校核：")
            print(f"       实际值：{u_actual:.4f} m/s")
            print(f"       允许范围：{self.U_ALLOWED_RANGE[0]:.1f} ~ {self.U_ALLOWED_RANGE[1]:.1f} m/s")

            if self.U_ALLOWED_RANGE[0] <= u_actual <= self.U_ALLOWED_RANGE[1]:
                print("       √ 在合理范围内")
            else:
                if u_actual < self.U_ALLOWED_RANGE[0]:
                    print(f"       × 偏低（低于最小值{self.U_ALLOWED_RANGE[0]:.1f} m/s）")
                else:
                    print(f"       × 偏高（高于最大值{self.U_ALLOWED_RANGE[1]:.1f} m/s）")

            print("\n   (2) 泛点率校核：")
            print(f"       实际值：{flooding_percent:.2f}%")
            print(f"       合理范围：{self.FLOODING_ALLOWED_RANGE[0]} ~ {self.FLOODING_ALLOWED_RANGE[1]}%")

            if self.FLOODING_ALLOWED_RANGE[0] <= flooding_percent <= self.FLOODING_ALLOWED_RANGE[1]:
                print("       √ 在合理范围内")
            else:
                if flooding_percent < self.FLOODING_ALLOWED_RANGE[0]:
                    print(f"       × 偏低（低于最小值{self.FLOODING_ALLOWED_RANGE[0]}%）")
                else:
                    print(f"       × 偏高（高于最大值{self.FLOODING_ALLOWED_RANGE[1]}%）")

            print("\n" + "-" * 50)

        # 总体评价
        print("\n三、总体评价与结论")
        print("-" * 40)

        all_ok = True
        for section_label, params in final_params.items():
            if not params['all_ok']:
                all_ok = False
                break

        if all_ok:
            print("✅ 所有参数均在合理范围内，设计合格。")
            print("\n设计结论：")
            print(f"1. 最终塔径确定为 D = {final_D:.3f} m")
            print(f"2. 精馏段板间距 HT = {self.HT_rectifying:.2f} m")
            print(f"3. 精馏段清液层高度 hL = {self.hL_rectifying:.3f} m")
            print(f"4. 提馏段板间距 HT = {self.HT_stripping:.2f} m")
            print(f"5. 提馏段清液层高度 hL = {self.hL_stripping:.3f} m")
            print("6. 各段空塔气速和泛点率均在允许范围内")
        else:
            print("⚠️ 部分参数超出允许范围，需要调整。")
            print("\n存在问题及建议：")

            for section_label, params in final_params.items():
                if not params['velocity_ok']:
                    if params['u_actual_final'] < self.U_ALLOWED_RANGE[0]:
                        print(f"   - {section_label}空塔气速({params['u_actual_final']:.4f} m/s)偏低")
                        print(f"     建议：适当减小塔径或提高设计泛点率")
                    else:
                        print(f"   - {section_label}空塔气速({params['u_actual_final']:.4f} m/s)偏高")
                        print(f"     建议：适当增大塔径或降低设计泛点率")

                if not params['flooding_ok']:
                    if params['flooding_percent_final'] < self.FLOODING_ALLOWED_RANGE[0]:
                        print(f"   - {section_label}泛点率({params['flooding_percent_final']:.2f}%)偏低")
                        print(f"     建议：适当减小塔径以提高设备利用率")
                    else:
                        print(f"   - {section_label}泛点率({params['flooding_percent_final']:.2f}%)偏高")
                        print(f"     建议：适当增大塔径以降低泛点率")

        print(f"\n{'=' * 80}")
        print("最终校核完成！")
        print(f"{'=' * 80}")

    # ==================== 新增功能：基于优化参数从头详细计算 ====================
    def print_optimized_parameters_detailed_calculation(self):
        """
        使用优化后的参数自动从头计算塔径（详细步骤）
        在自动优化后自动调用，展示完整的计算过程
        """
        print(f"\n{'=' * 80}")
        print("【优化后参数详细计算】基于优化参数从头计算塔径（完整步骤）")
        print(f"{'=' * 80}")

        print("\n说明：")
        print("1. 使用优化后的板间距和清液层高度")
        print("2. 按照优化得到的设计泛点率从头计算塔径")
        print("3. 详细展示每一步的计算公式、代入数值的过程和中间结果")
        print("4. 最后进行校核验证")
        print("5. 输出格式适合抄写到课程设计报告中")

        # 使用优化后的参数
        HT_rect = self.HT_rectifying
        hL_rect = self.hL_rectifying
        HT_strip = self.HT_stripping
        hL_strip = self.hL_stripping
        flooding_fraction = self.design_flooding_fraction

        print(f"\n使用的优化参数:")
        print(f"  精馏段: HT={HT_rect:.2f}m, hL={hL_rect:.3f}m")
        print(f"  提馏段: HT={HT_strip:.2f}m, hL={hL_strip:.3f}m")
        print(f"  设计泛点率: f={flooding_fraction:.2f}")

        print(f"\n{'=' * 80}")
        print("第一部分：精馏段塔径详细计算")
        print(f"{'=' * 80}")

        # 精馏段详细计算
        rect_results = self._calculate_section_full_detail('精馏段结果', HT_rect, hL_rect, flooding_fraction)

        print(f"\n{'=' * 80}")
        print("第二部分：提馏段塔径详细计算")
        print(f"{'=' * 80}")

        # 提馏段详细计算
        strip_results = self._calculate_section_full_detail('提馏段结果', HT_strip, hL_strip, flooding_fraction)

        print(f"\n{'=' * 80}")
        print("第三部分：确定最终塔径")
        print(f"{'=' * 80}")

        D_rect = rect_results['D_rounded']
        D_strip = strip_results['D_rounded']

        print(f"精馏段计算塔径: D_rect = {D_rect:.3f} m")
        print(f"提馏段计算塔径: D_strip = {D_strip:.3f} m")

        final_D = max(D_rect, D_strip)

        if final_D == D_rect:
            print(f"最终塔径取较大值: D = {final_D:.3f} m (由精馏段控制)")
            controlling_section = "精馏段"
        else:
            print(f"最终塔径取较大值: D = {final_D:.3f} m (由提馏段控制)")
            controlling_section = "提馏段"

        print(f"\n{'=' * 80}")
        print("第四部分：最终塔径下详细校核")
        print(f"{'=' * 80}")

        # 在最终塔径下重新校核
        final_check = self._verify_final_diameter_detail(final_D, HT_rect, hL_rect, HT_strip, hL_strip)

        print(f"\n{'=' * 80}")
        print("第五部分：设计总结")
        print(f"{'=' * 80}")

        all_ok = True
        for section in ['精馏段', '提馏段']:
            if section in final_check and not final_check[section]['all_ok']:
                all_ok = False
                break

        if all_ok:
            print("✅ 设计合格！所有参数均在合理范围内")
        else:
            print("⚠️ 设计存在问题，请检查参数")

        print(f"\n最终设计结果:")
        print(f"  最终塔径: D = {final_D:.3f} m")
        print(f"  控制塔段: {controlling_section}")
        print(f"  精馏段板间距: HT = {HT_rect:.2f} m")
        print(f"  精馏段清液层高度: hL = {hL_rect:.3f} m")
        print(f"  提馏段板间距: HT = {HT_strip:.2f} m")
        print(f"  提馏段清液层高度: hL = {hL_strip:.3f} m")
        print(f"  设计泛点率: f = {flooding_fraction:.2f}")

        # 保存详细报告
        self._save_detailed_report(
            rect_results, strip_results, final_D, controlling_section,
            HT_rect, hL_rect, HT_strip, hL_strip, flooding_fraction, final_check
        )

        return {
            '精馏段结果': rect_results,
            '提馏段结果': strip_results,
            '最终塔径': final_D,
            '控制塔段': controlling_section,
            '最终校核': final_check
        }

    def _calculate_section_full_detail(self, section_name, HT, hL, design_fraction):
        """
        详细计算单个塔段塔径（完整步骤，包含所有公式）
        """
        if not self.tower_data or section_name not in self.tower_data:
            return None

        data = self.tower_data[section_name]
        section_label = section_name.replace('结果', '')

        print(f"\n【{section_label}】详细计算:")
        print("-" * 60)

        # 1. 计算流动参数F_LV
        print("\n1. 计算流动参数F_LV:")
        print("   公式: F_LV = (L_mass / V_mass) × √(ρ_V / ρ_L)")

        L_mass = data['L_mass_kgh']
        V_mass = data['V_mass_kgh']
        ρ_L = data['rho_liquid_kgm3']
        ρ_V = data['rho_vapor_kgm3']

        print(f"   已知:")
        print(f"     L = {L_mass:.1f} kg/h")
        print(f"     V = {V_mass:.1f} kg/h")
        print(f"     ρ_L = {ρ_L:.1f} kg/m³")
        print(f"     ρ_V = {ρ_V:.4f} kg/m³")

        mass_ratio = L_mass / V_mass
        density_ratio_sqrt = math.sqrt(ρ_V / ρ_L)
        F_LV = mass_ratio * density_ratio_sqrt

        print(f"   计算:")
        print(f"     L/V = {L_mass:.1f}/{V_mass:.1f} = {mass_ratio:.4f}")
        print(f"     √(ρ_V/ρ_L) = √({ρ_V:.4f}/{ρ_L:.1f}) = {density_ratio_sqrt:.4f}")
        print(f"     F_LV = {mass_ratio:.4f} × {density_ratio_sqrt:.4f} = {F_LV:.4f}")

        # 2. 计算HT-hL
        print(f"\n2. 计算HT-hL:")
        HT_minus_hL = HT - hL
        print(f"   已知: HT = {HT:.2f} m, hL = {hL:.3f} m")
        print(f"   计算: HT - hL = {HT:.2f} - {hL:.3f} = {HT_minus_hL:.3f} m")

        # 3. 计算C20（史密斯经验公式）
        print(f"\n3. 使用史密斯经验公式计算C20:")
        print("   公式: C20 = 0.0162 - 0.0648×X + 0.181×Y + 0.0162×X² - 0.139×X×Y + 0.185×Y²")
        print(f"   其中: X = F_LV = {F_LV:.4f}, Y = HT - hL = {HT_minus_hL:.3f}")

        X = F_LV
        Y = HT_minus_hL
        C20 = (0.0162 - 0.0648 * X + 0.181 * Y +
               0.0162 * X ** 2 - 0.139 * X * Y + 0.185 * Y ** 2)

        print(f"   逐步计算:")
        term1 = 0.0162
        term2 = -0.0648 * X
        term3 = 0.181 * Y
        term4 = 0.0162 * X ** 2
        term5 = -0.139 * X * Y
        term6 = 0.185 * Y ** 2

        print(f"     0.0162 = {term1:.6f}")
        print(f"     -0.0648×{X:.4f} = {term2:.6f}")
        print(f"     0.181×{Y:.3f} = {term3:.6f}")
        print(f"     0.0162×({X:.4f})² = {term4:.6f}")
        print(f"     -0.139×{X:.4f}×{Y:.3f} = {term5:.6f}")
        print(f"     0.185×({Y:.3f})² = {term6:.6f}")

        sum_parts = term1 + term2 + term3 + term4 + term5 + term6
        print(f"     总和 = {sum_parts:.6f}")

        # 限制在合理范围
        if C20 < 0.03:
            C20 = 0.03
            print(f"   应用最小值限制: C20 = 0.0300")
        elif C20 > 0.15:
            C20 = 0.15
            print(f"   应用最大值限制: C20 = 0.1500")
        else:
            print(f"   计算结果: C20 = {C20:.4f}")

        # 4. 表面张力校正
        print(f"\n4. 表面张力校正:")
        print("   公式: C = C20 × (σ/20)^0.2")
        σ = data['sigma_mNm']
        print(f"   已知表面张力: σ = {σ:.2f} mN/m")
        C = C20 * (σ / 20) ** 0.2
        print(f"   计算: C = {C20:.4f} × ({σ:.2f}/20)^0.2 = {C:.4f}")

        # 5. 计算液泛气速
        print(f"\n5. 计算液泛气速u_f:")
        print("   公式: u_f = C × √[(ρ_L - ρ_V) / ρ_V]")
        u_f = C * math.sqrt((ρ_L - ρ_V) / ρ_V)
        print(f"   计算:")
        density_diff = ρ_L - ρ_V
        density_ratio = density_diff / ρ_V
        print(f"     ρ_L - ρ_V = {ρ_L:.1f} - {ρ_V:.4f} = {density_diff:.4f}")
        print(f"     (ρ_L - ρ_V)/ρ_V = {density_diff:.4f}/{ρ_V:.4f} = {density_ratio:.4f}")
        sqrt_density_ratio = math.sqrt(density_ratio)
        print(f"     √[(ρ_L - ρ_V)/ρ_V] = √{density_ratio:.4f} = {sqrt_density_ratio:.4f}")
        print(f"     u_f = {C:.4f} × {sqrt_density_ratio:.4f} = {u_f:.4f} m/s")

        # 6. 计算设计空塔气速
        print(f"\n6. 计算设计空塔气速u_design:")
        print("   公式: u_design = f × u_f")
        print(f"   设计泛点率: f = {design_fraction:.2f}")
        u_design = design_fraction * u_f
        print(f"   计算: u_design = {design_fraction:.2f} × {u_f:.4f} = {u_design:.4f} m/s")

        # 7. 计算气体体积流量
        print(f"\n7. 计算气体体积流量V:")
        V_m3h = data['V_vol_m3h']
        V_m3s = V_m3h / 3600
        print(f"   已知: V = {V_m3h:.2f} m³/h")
        print(f"   单位换算: V = {V_m3h:.2f} ÷ 3600 = {V_m3s:.6f} m³/s")

        # 8. 计算理论塔截面积
        print(f"\n8. 计算理论塔截面积A_t:")
        print("   公式: A_t = V / u_design")
        A_t = V_m3s / u_design
        print(f"   计算: A_t = {V_m3s:.6f} ÷ {u_design:.4f} = {A_t:.4f} m²")

        # 9. 计算理论塔径
        print(f"\n9. 计算理论塔径D_calc:")
        print("   公式: D_calc = √(4 × A_t / π)")
        D_calc = math.sqrt(4 * A_t / math.pi)
        print(f"   计算:")
        print(f"     4 × A_t / π = 4 × {A_t:.4f} ÷ 3.1416 = {4 * A_t / math.pi:.4f}")
        print(f"     D_calc = √{4 * A_t / math.pi:.4f} = {D_calc:.4f} m")

        # 10. 塔径圆整
        print(f"\n10. 塔径圆整:")
        D_rounded = min(self.STANDARD_DIAMETERS, key=lambda x: abs(x - D_calc))
        print(f"   计算塔径: {D_calc:.4f} m")
        print(f"   圆整为标准塔径: {D_rounded:.3f} m")

        # 11. 计算实际塔截面积
        print(f"\n11. 计算实际塔截面积A_actual:")
        print("   公式: A_actual = π × D_rounded² / 4")
        A_actual = math.pi * D_rounded ** 2 / 4
        print(f"   计算: A_actual = 3.1416 × ({D_rounded:.3f})² ÷ 4 = {A_actual:.4f} m²")

        # 12. 计算实际空塔气速
        print(f"\n12. 计算实际空塔气速u_actual:")
        print("   公式: u_actual = V / A_actual")
        u_actual = V_m3s / A_actual
        print(f"   计算: u_actual = {V_m3s:.6f} ÷ {A_actual:.4f} = {u_actual:.4f} m/s")

        # 13. 计算实际泛点率
        print(f"\n13. 计算实际泛点率:")
        print("   公式: 实际泛点率 = (u_actual / u_f) × 100%")
        flooding_percent = (u_actual / u_f) * 100
        print(f"   计算: 实际泛点率 = ({u_actual:.4f} ÷ {u_f:.4f}) × 100% = {flooding_percent:.2f}%")

        # 14. 参数校验
        print(f"\n14. 参数校验:")
        all_ok, velocity_ok, flooding_ok = self.check_parameters(u_actual, flooding_percent)

        print(f"   (1) 空塔气速校验:")
        print(f"       实际值: {u_actual:.4f} m/s")
        print(f"       允许范围: {self.U_ALLOWED_RANGE[0]:.1f} ~ {self.U_ALLOWED_RANGE[1]:.1f} m/s")
        if velocity_ok:
            print(f"       ✅ 在合理范围内")
        else:
            if u_actual < self.U_ALLOWED_RANGE[0]:
                print(f"       ❌ 偏低（低于最小值{self.U_ALLOWED_RANGE[0]:.1f} m/s）")
            else:
                print(f"       ❌ 偏高（高于最大值{self.U_ALLOWED_RANGE[1]:.1f} m/s）")

        print(f"\n   (2) 泛点率校验:")
        print(f"       实际值: {flooding_percent:.2f}%")
        print(f"       合理范围: {self.FLOODING_ALLOWED_RANGE[0]} ~ {self.FLOODING_ALLOWED_RANGE[1]}%")
        if flooding_ok:
            print(f"       ✅ 在合理范围内")
        else:
            if flooding_percent < self.FLOODING_ALLOWED_RANGE[0]:
                print(f"       ❌ 偏低（低于最小值{self.FLOODING_ALLOWED_RANGE[0]}%）")
            else:
                print(f"       ❌ 偏高（高于最大值{self.FLOODING_ALLOWED_RANGE[1]}%）")

        print(f"\n   总体校验结果: {'✅ 通过' if all_ok else '❌ 未通过'}")

        return {
            'section': section_label,
            'F_LV': round(F_LV, 4),
            'HT': round(HT, 2),
            'hL': round(hL, 3),
            'HT_minus_hL': round(HT_minus_hL, 3),
            'C20': round(C20, 4),
            'C': round(C, 4),
            'u_f': round(u_f, 4),
            'u_design': round(u_design, 4),
            'u_actual': round(u_actual, 4),
            'flooding_percent': round(flooding_percent, 2),
            'D_calc': round(D_calc, 4),
            'D_rounded': round(D_rounded, 3),
            'A_actual': round(A_actual, 4),
            'V_flow_m3s': round(V_m3s, 6),
            'parameters_ok': all_ok,
            'velocity_ok': velocity_ok,
            'flooding_ok': flooding_ok
        }

    def _verify_final_diameter_detail(self, final_D, HT_rect, hL_rect, HT_strip, hL_strip):
        """在给定塔径下进行最终详细校核"""
        print(f"最终塔径: D = {final_D:.3f} m")

        A_final = math.pi * final_D ** 2 / 4
        print(f"塔截面积: A = π × D² / 4 = 3.1416 × ({final_D:.3f})² ÷ 4 = {A_final:.4f} m²")

        final_results = {}

        for section_name, HT, hL, section_label in [
            ('精馏段结果', HT_rect, hL_rect, '精馏段'),
            ('提馏段结果', HT_strip, hL_strip, '提馏段')
        ]:
            print(f"\n【{section_label}】最终校核:")
            print("-" * 40)

            data = self.tower_data[section_name]

            # 气体体积流量
            V_m3h = data['V_vol_m3h']
            V_m3s = V_m3h / 3600
            print(f"1. 气体体积流量:")
            print(f"   V = {V_m3h:.2f} m³/h = {V_m3h:.2f} ÷ 3600 = {V_m3s:.6f} m³/s")

            # 实际空塔气速
            u_actual = V_m3s / A_final
            print(f"2. 实际空塔气速:")
            print(f"   u_实际 = V / A = {V_m3s:.6f} ÷ {A_final:.4f} = {u_actual:.4f} m/s")

            # 重新计算F_LV和C值
            F_LV = self.calculate_F_LV(section_name)
            HT_minus_hL = HT - hL
            C20 = self.approximate_C20_from_smith_chart(F_LV, HT_minus_hL)
            σ = data['sigma_mNm']
            C = C20 * (σ / 20) ** 0.2

            print(f"3. 重新计算流动参数和C值:")
            print(f"   F_LV = {F_LV:.4f}")
            print(f"   HT = {HT:.2f} m, hL = {hL:.3f} m")
            print(f"   HT - hL = {HT_minus_hL:.3f} m")
            print(f"   C20 = {C20:.4f}")
            print(f"   C = {C20:.4f} × ({σ:.2f}/20)^0.2 = {C:.4f}")

            # 液泛气速
            ρ_L = data['rho_liquid_kgm3']
            ρ_V = data['rho_vapor_kgm3']
            u_f = C * math.sqrt((ρ_L - ρ_V) / ρ_V)
            print(f"4. 液泛气速:")
            print(f"   u_f = {C:.4f} × √[({ρ_L:.1f} - {ρ_V:.4f})/{ρ_V:.4f}] = {u_f:.4f} m/s")

            # 实际泛点率
            flooding_percent = (u_actual / u_f) * 100
            print(f"5. 实际泛点率:")
            print(f"   实际泛点率 = ({u_actual:.4f} ÷ {u_f:.4f}) × 100% = {flooding_percent:.2f}%")

            # 参数校验
            all_ok, velocity_ok, flooding_ok = self.check_parameters(u_actual, flooding_percent)

            print(f"\n6. 校核结果:")
            if all_ok:
                print(f"   ✅ {section_label}校核通过")
            else:
                print(f"   ❌ {section_label}校核未通过")

            # 存储结果
            final_results[section_label] = {
                'u_actual': round(u_actual, 4),
                'flooding_percent': round(flooding_percent, 2),
                'velocity_ok': velocity_ok,
                'flooding_ok': flooding_ok,
                'all_ok': all_ok,
                'HT': HT,
                'hL': hL
            }

        return final_results

    def _save_detailed_report(self, rect_results, strip_results, final_D, controlling_section,
                              HT_rect, hL_rect, HT_strip, hL_strip, flooding_fraction, final_check):
        """保存详细计算报告"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"优化后参数详细计算报告_{timestamp}.txt"

        try:
            with open(filename, 'w', encoding='utf-8') as f:
                f.write("=" * 80 + "\n")
                f.write("苯-甲苯精馏塔塔径详细计算报告（优化后参数）\n")
                f.write("=" * 80 + "\n\n")

                f.write("一、优化后参数\n")
                f.write("-" * 40 + "\n")
                f.write(f"精馏段: HT = {HT_rect:.2f}m, hL = {hL_rect:.3f}m\n")
                f.write(f"提馏段: HT = {HT_strip:.2f}m, hL = {hL_strip:.3f}m\n")
                f.write(f"设计泛点率: f = {flooding_fraction:.2f}\n\n")

                f.write("二、精馏段计算结果\n")
                f.write("-" * 40 + "\n")
                self._write_section_summary(f, rect_results, '精馏段')

                f.write("\n三、提馏段计算结果\n")
                f.write("-" * 40 + "\n")
                self._write_section_summary(f, strip_results, '提馏段')

                f.write("\n四、最终塔径\n")
                f.write("-" * 40 + "\n")
                f.write(f"精馏段塔径: {rect_results['D_rounded']:.3f} m\n")
                f.write(f"提馏段塔径: {strip_results['D_rounded']:.3f} m\n")
                f.write(f"最终塔径: {final_D:.3f} m\n")
                f.write(f"控制塔段: {controlling_section}\n\n")

                f.write("五、最终校核结果\n")
                f.write("-" * 40 + "\n")

                all_ok = True
                for section in ['精馏段', '提馏段']:
                    if section in final_check:
                        result = final_check[section]
                        f.write(f"{section}:\n")
                        f.write(f"  空塔气速: {result['u_actual']:.4f} m/s")
                        if result['velocity_ok']:
                            f.write(" (合理)\n")
                        else:
                            f.write(" (不合理)\n")

                        f.write(f"  泛点率: {result['flooding_percent']:.2f}%")
                        if result['flooding_ok']:
                            f.write(" (合理)\n")
                        else:
                            f.write(" (不合理)\n")

                        if not result['all_ok']:
                            all_ok = False

                f.write("\n六、设计结论\n")
                f.write("-" * 40 + "\n")
                if all_ok:
                    f.write("✅ 所有参数均在合理范围内，设计合格。\n")
                else:
                    f.write("❌ 部分参数超出允许范围，需要调整。\n")

                f.write(f"\n报告生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

            print(f"\n✓ 详细计算报告已保存到: '{filename}'")

        except Exception as e:
            print(f"保存报告时出错: {e}")

    def _write_section_summary(self, file, results, section_name):
        """将塔段计算结果摘要写入文件"""
        file.write(f"{section_name}:\n")
        file.write(f"  流动参数 F_LV = {results['F_LV']:.4f}\n")
        file.write(f"  板间距 HT = {results['HT']:.2f} m\n")
        file.write(f"  清液层高度 hL = {results['hL']:.3f} m\n")
        file.write(f"  C20 = {results['C20']:.4f}\n")
        file.write(f"  校正后 C = {results['C']:.4f}\n")
        file.write(f"  液泛气速 u_f = {results['u_f']:.4f} m/s\n")
        file.write(f"  实际空塔气速 u_实际 = {results['u_actual']:.4f} m/s\n")
        file.write(f"  理论塔径 D_calc = {results['D_calc']:.4f} m\n")
        file.write(f"  圆整塔径 D_rounded = {results['D_rounded']:.3f} m\n")
        file.write(f"  实际泛点率 = {results['flooding_percent']:.2f}%\n")
        file.write(f"  校核结果: {'通过' if results['parameters_ok'] else '未通过'}\n")

    def run_complete_calculation(self):
        """运行完整的塔径计算"""
        if not self.tower_data:
            print("✗ 错误: 未加载塔工艺数据，无法进行计算")
            return None

        print("塔径计算程序启动（史密斯经验公式法）")
        print("=" * 80)
        print("说明：使用史密斯关联图经验公式计算塔径")
        print("     优点：无需手动查图，自动计算C20")
        print("=" * 80)

        # 显示校验参数范围
        print(f"\n校验参数范围:")
        print(f"  空塔气速允许范围: {self.U_ALLOWED_RANGE[0]:.1f} - {self.U_ALLOWED_RANGE[1]:.1f} m/s")
        print(f"  泛点率合理范围: {self.FLOODING_ALLOWED_RANGE[0]} - {self.FLOODING_ALLOWED_RANGE[1]}%")
        print(f"  板间距调整范围: {self.HT_RANGE[0]} - {self.HT_RANGE[1]} m")
        print(f"  清液层高度调整范围: {self.HL_RANGE[0]} - {self.HL_RANGE[1]} m")

        # 询问是否修改校验范围
        change_range = 'n'  # input("是否修改校验参数范围？(y/n): ").lower()
        if change_range == 'y':
            try:
                u_min = float(input(f"请输入最小空塔气速 (当前 {self.U_ALLOWED_RANGE[0]:.1f} m/s): "))
                u_max = float(input(f"请输入最大空塔气速 (当前 {self.U_ALLOWED_RANGE[1]:.1f} m/s): "))
                self.U_ALLOWED_RANGE = (u_min, u_max)

                flood_min = float(input(f"请输入最小泛点率 (当前 {self.FLOODING_ALLOWED_RANGE[0]}%): "))
                flood_max = float(input(f"请输入最大泛点率 (当前 {self.FLOODING_ALLOWED_RANGE[1]}%): "))
                self.FLOODING_ALLOWED_RANGE = (flood_min, flood_max)

                ht_min = float(input(f"请输入最小板间距 (当前 {self.HT_RANGE[0]} m): "))
                ht_max = float(input(f"请输入最大板间距 (当前 {self.HT_RANGE[1]} m): "))
                self.HT_RANGE = (ht_min, ht_max)

                hl_min = float(input(f"请输入最小清液层高度 (当前 {self.HL_RANGE[0]} m): "))
                hl_max = float(input(f"请输入最大清液层高度 (当前 {self.HL_RANGE[1]} m): "))
                self.HL_RANGE = (hl_min, hl_max)

                print(f"已更新校验范围:")
                print(f"  空塔气速 {u_min:.1f}-{u_max:.1f} m/s, 泛点率 {flood_min}-{flood_max}%")
                print(f"  板间距 {ht_min:.2f}-{ht_max:.2f} m, 清液层高度 {hl_min:.3f}-{hl_max:.3f} m")
            except ValueError:
                print("输入无效，使用默认范围")

        # 设置设计泛点率
        print(f"\n当前设计泛点率: {self.design_flooding_fraction}")
        change_fraction = input("是否修改设计泛点率？(y/n): ").lower()
        if change_fraction == 'y':
            try:
                fraction = float(
                    input(f"请输入设计泛点率 ({self.FLOODING_FRACTION_RANGE[0]}-{self.FLOODING_FRACTION_RANGE[1]}): "))
                if self.FLOODING_FRACTION_RANGE[0] <= fraction <= self.FLOODING_FRACTION_RANGE[1]:
                    self.design_flooding_fraction = fraction
                else:
                    print("超出合理范围，使用默认值")
            except ValueError:
                print("输入无效，使用默认值")

        # 设置清液层高度
        print(f"\n当前清液层高度设置:")
        print(f"  精馏段 hL = {self.hL_rectifying:.3f} m")
        print(f"  提馏段 hL = {self.hL_stripping:.3f} m")

        change_hL = input("是否修改清液层高度？(y/n): ").lower()
        if change_hL == 'y':
            try:
                hL_rect = float(input(f"请输入精馏段清液层高度hL (m): "))
                hL_strip = float(input(f"请输入提馏段清液层高度hL (m): "))
                self.hL_rectifying = hL_rect
                self.hL_stripping = hL_strip
            except ValueError:
                print("输入无效，使用默认值")

        # 精馏段计算
        print("\n" + "=" * 80)
        print("第一步：精馏段计算")
        print("=" * 80)

        # 使用初始参数计算
        results_rect = self.calculate_section_diameter_with_params(
            '精馏段结果',
            self.HT_rectifying,
            self.hL_rectifying,
            self.design_flooding_fraction
        )

        # 如果参数不合理，启动自动优化
        if not results_rect['parameters_ok']:
            print(f"\n⚠️ 参数校验未通过，启动自动参数调整...")
            optimized_rect = self.optimize_section_parameters('精馏段结果')
            if optimized_rect:
                results_rect = optimized_rect

        # 提馏段计算
        print("\n\n" + "=" * 80)
        print("第二步：提馏段计算")
        print("=" * 80)

        # 使用初始参数计算
        results_strip = self.calculate_section_diameter_with_params(
            '提馏段结果',
            self.HT_stripping,
            self.hL_stripping,
            self.design_flooding_fraction
        )

        # 如果参数不合理，启动自动优化
        if not results_strip['parameters_ok']:
            print(f"\n⚠️ 参数校验未通过，启动自动参数调整...")
            optimized_strip = self.optimize_section_parameters('提馏段结果')
            if optimized_strip:
                results_strip = optimized_strip

        return results_rect, results_strip

    def determine_final_tower_diameter(self, results_rect, results_strip, auto_optimize_final=True):
        """确定最终塔径，并详细核算参数"""
        print("\n" + "=" * 80)
        print("第三步：确定最终塔径")
        print("=" * 80)

        D_rect = results_rect['D_rounded']
        D_strip = results_strip['D_rounded']

        final_D = max(D_rect, D_strip)

        print(f"精馏段塔径: {D_rect:.3f} m")
        print(f"提馏段塔径: {D_strip:.3f} m")
        print(f"最终塔径: {final_D:.3f} m (取较大值)")

        # 确定由哪个塔段控制
        if final_D == D_rect:
            controlling_section = "精馏段"
        else:
            controlling_section = "提馏段"

        print(f"控制塔段: {controlling_section}")

        # 用最终塔径详细核算参数
        print(f"\n" + "=" * 80)
        print("第四步：最终塔径下详细参数核算")
        print("=" * 80)

        # 计算最终塔径下的详细参数
        final_params = self.calculate_final_diameter_parameters(final_D, results_rect, results_strip)

        # 检查是否所有参数都满足要求
        all_ok = True
        warnings = []

        for section_label, params in final_params.items():
            if not params['all_ok']:
                all_ok = False
                warnings.append(f"{section_label}参数不满足要求")
                if not params['velocity_ok']:
                    warnings.append(
                        f"  {section_label}空塔气速({params['u_actual_final']:.4f} m/s)超出范围 {self.U_ALLOWED_RANGE[0]:.1f}-{self.U_ALLOWED_RANGE[1]:.1f} m/s")
                if not params['flooding_ok']:
                    warnings.append(
                        f"  {section_label}泛点率({params['flooding_percent_final']:.2f}%)超出范围 {self.FLOODING_ALLOWED_RANGE[0]}-{self.FLOODING_ALLOWED_RANGE[1]}%")

        # 如果参数不满足且启用了自动优化
        if not all_ok and auto_optimize_final:
            print(f"\n⚠️ 最终塔径下参数不满足要求，启动自动优化...")

            # 询问用户是否要进行自动优化
            optimize_choice = 'y'  # input("是否进行自动参数优化以找到满足条件的参数组合？(y/n): ").lower()
            if optimize_choice == 'y':
                optimized_result = self.optimize_final_diameter_parameters(results_rect, results_strip)

                if optimized_result:
                    print(f"\n✓ 最终塔径参数优化成功！")

                    # 使用优化后的参数重新计算
                    final_params = self.calculate_final_diameter_parameters(final_D, results_rect, results_strip)

                    # 重新检查参数
                    all_ok = True
                    for params in final_params.values():
                        if not params['all_ok']:
                            all_ok = False
                            break

                    warnings = []
                else:
                    print(f"\n⚠️ 自动优化失败，请手动调整参数")
            else:
                print(f"\n⚠️ 用户选择不进行自动优化")

        # 生成设计总结
        print(f"\n" + "=" * 80)
        print("第五步：塔径设计总结")
        print("=" * 80)
        print(f"最终塔径: {final_D:.3f} m")
        print(f"设计方法: 史密斯关联图经验公式法")
        print(f"板间距: 精馏段 {self.HT_rectifying:.2f}m, 提馏段 {self.HT_stripping:.2f}m")
        print(f"清液层高度: 精馏段 {self.hL_rectifying:.3f}m, 提馏段 {self.hL_stripping:.3f}m")

        print(f"\n最终校验结果:")
        if all_ok:
            print("✅ 所有参数均在合理范围内，设计合格")
        else:
            print("⚠️ 注意: 部分参数超出允许范围")
            for warning in warnings:
                print(f"  - {warning}")

            print(f"\n设计建议:")
            if final_D == D_rect and final_D > D_strip:
                print("  - 当前塔径由精馏段控制，可考虑使用异径塔设计")
            elif final_D == D_strip and final_D > D_rect:
                print("  - 当前塔径由提馏段控制，可考虑使用异径塔设计")

            # 针对泛点率的调整建议
            for section_label, params in final_params.items():
                if params['flooding_percent_final'] < self.FLOODING_ALLOWED_RANGE[0]:
                    print(f"  - {section_label}泛点率({params['flooding_percent_final']:.2f}%)偏低，建议增大设计泛点率")
                elif params['flooding_percent_final'] > self.FLOODING_ALLOWED_RANGE[1]:
                    print(f"  - {section_label}泛点率({params['flooding_percent_final']:.2f}%)偏高，建议减小设计泛点率")

            print("  - 如需调整，可重新运行程序修改设计参数或调整校验范围")

        # ============ 新增：自动进行详细计算 ============
        print(f"\n{'=' * 80}")
        print("第六步：基于优化参数从头详细计算（适合课程设计报告）")
        print("=" * 80)

        print("说明：")
        print("1. 使用优化后的板间距和清液层高度")
        print("2. 从头计算塔径，展示完整公式和计算步骤")
        print("3. 适合抄写到课程设计报告中")

        # 自动开始详细计算
        self.print_optimized_parameters_detailed_calculation()
        # ===========================================

        return final_D, controlling_section, final_params

    def save_results_to_json(self, results_rect, results_strip, final_D, controlling_section, final_params,
                             filename='tower_diameter_approximate_optimized.json'):
        """将计算结果保存为JSON文件"""

        results_dict = {
            "设计信息": {
                "项目": "苯-甲苯精馏塔塔径设计（经验公式优化版）",
                "设计方法": "史密斯关联图经验公式法",
                "设计日期": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "特点": "使用经验公式估算C20，带自动优化功能"
            },
            "设计参数": {
                "板间距": {
                    "精馏段_m": round(self.HT_rectifying, 2),
                    "提馏段_m": round(self.HT_stripping, 2)
                },
                "清液层高度": {
                    "精馏段_m": round(self.hL_rectifying, 3),
                    "提馏段_m": round(self.hL_stripping, 3)
                },
                "设计泛点率": self.design_flooding_fraction
            },
            "精馏段计算结果": results_rect,
            "提馏段计算结果": results_strip,
            "最终设计结果": {
                "最终塔径_m": round(final_D, 3),
                "控制塔段": controlling_section,
                "精馏段塔径_m": results_rect['D_rounded'],
                "提馏段塔径_m": results_strip['D_rounded'],
                "塔径确定原则": "取精馏段和提馏段中较大的塔径作为全塔直径"
            },
            "最终塔径下详细参数": final_params,
            "校验参数范围": {
                "空塔气速允许范围_m/s": self.U_ALLOWED_RANGE,
                "泛点率合理范围_%": self.FLOODING_ALLOWED_RANGE,
                "板间距调整范围_m": self.HT_RANGE,
                "清液层高度调整范围_m": self.HL_RANGE,
                "设计泛点率调整范围": self.FLOODING_FRACTION_RANGE
            }
        }

        # 如果有优化后的结果，也保存
        if self.final_results and self.final_results.get('optimized'):
            results_dict["优化后参数"] = {
                "优化状态": "已优化",
                "优化后板间距": {
                    "精馏段_m": round(self.final_results.get('HT_rectifying', self.HT_rectifying), 2),
                    "提馏段_m": round(self.final_results.get('HT_stripping', self.HT_stripping), 2)
                },
                "优化后清液层高度": {
                    "精馏段_m": round(self.final_results.get('hL_rectifying', self.hL_rectifying), 3),
                    "提馏段_m": round(self.final_results.get('hL_stripping', self.hL_stripping), 3)
                }
            }

        try:
            with open(filename, 'w', encoding='utf-8') as f:
                json.dump(results_dict, f, ensure_ascii=False, indent=4)
            print(f"\n✓ 计算结果已保存到 '{filename}'")

            # 显示文件大小
            file_size = os.path.getsize(filename)
            print(f"   文件大小: {file_size} 字节")

        except Exception as e:
            print(f"✗ 保存JSON文件时出错: {e}")


def main():
    """主函数"""
    print("苯-甲苯精馏塔塔径计算程序（经验公式优化版）")
    print("=" * 80)
    print("功能说明:")
    print("  1. 使用史密斯关联图经验公式自动计算C20值")
    print("  2. 当计算结果不满足校验条件时，程序会自动尝试调整参数")
    print("  3. 调整的参数包括：板间距、清液层高度、设计泛点率")
    print("  4. 最终塔径确定后会详细核算参数，与精馏段和提馏段计算一样详细")
    print("  5. 自动输出基于优化参数的详细计算过程（适合课程设计报告）")
    print("=" * 80)

    calculator = TowerDiameterApproximateOptimizedCalculator('tower_properties_results.json')

    results = calculator.run_complete_calculation()

    if results:
        results_rect, results_strip = results
        final_D, controlling_section, final_params = calculator.determine_final_tower_diameter(
            results_rect, results_strip
        )

        # 保存结果为JSON文件
        calculator.save_results_to_json(
            results_rect, results_strip, final_D, controlling_section, final_params,
            'tower_diameter_approximate_optimized.json'
        )

        print("\n" + "=" * 80)
        print("塔径计算完成！")
        print("=" * 80)
        print(f"最终塔径: {final_D:.3f} m")
        print(f"控制塔段: {controlling_section}")
        print(f"设计方法: 史密斯关联图经验公式法")
        print(f"计算结果已保存到 'tower_diameter_approximate_optimized.json'")

        # 返回结果供可能的外部调用
        return {
            'final_diameter': final_D,
            'controlling_section': controlling_section,
            'rectifying_results': results_rect,
            'stripping_results': results_strip,
            'final_params': final_params
        }
    else:
        print("✗ 计算失败，请检查数据和输入")
        return None


if __name__ == "__main__":
    # 运行主程序
    design_results = main()