import json
import numpy as np
import os
import math

def linear_interpolate(value, data_dict):
    """
    通用线性插值函数
    参数:
    value: 要插值的自变量值
    data_dict: 字典，键为自变量，值为因变量
    返回:
    插值结果
    """
    if value in data_dict:
        return data_dict[value]
    keys = sorted(data_dict.keys())
    if value <= keys[0]:
        return data_dict[keys[0]]
    if value >= keys[-1]:
        return data_dict[keys[-1]]
    for i in range(len(keys) - 1):
        if keys[i] <= value <= keys[i + 1]:
            x1, x2 = keys[i], keys[i + 1]
            y1, y2 = data_dict[x1], data_dict[x2]
            result = y1 + (y2 - y1) * (value - x1) / (x2 - x1)
            return result
    return data_dict[keys[-1]]


def get_interpolation_formula(value, data_dict, property_name, component_name):
    """
    获取插值计算过程和公式
    """
    if value in data_dict:
        return f"插值公式: y = y₁ (当x=x₁)\n计算过程: 温度{value:.2f}℃在数据点上，直接取{property_name}={data_dict[value]:.2f}"

    keys = sorted(data_dict.keys())

    if value <= keys[0]:
        x1 = keys[0]
        return f"插值公式: y = y₁ (外推，x≤x₁)\n计算过程: 温度{value:.2f}℃≤{x1}℃，取{property_name}={data_dict[x1]:.2f}"

    if value >= keys[-1]:
        x1 = keys[-1]
        return f"插值公式: y = y₁ (外推，x≥x₁)\n计算过程: 温度{value:.2f}℃≥{x1}℃，取{property_name}={data_dict[x1]:.2f}"

    for i in range(len(keys) - 1):
        if keys[i] <= value <= keys[i + 1]:
            x1, x2 = keys[i], keys[i + 1]
            y1, y2 = data_dict[x1], data_dict[x2]

            formula = f"插值公式: y = y₁ + (y₂ - y₁) × (x - x₁) / (x₂ - x₁)"
            process = f"计算过程: {component_name}{property_name}在{value:.2f}℃\n"
            process += f"         x₁={x1}℃, y₁={y1:.2f}; x₂={x2}℃, y₂={y2:.2f}\n"
            process += f"         = {y1:.2f} + ({y2:.2f} - {y1:.2f}) × ({value:.2f} - {x1}) / ({x2} - {x1})\n"
            process += f"         = {y1:.2f} + {y2 - y1:.2f} × {value - x1:.1f} / {x2 - x1:.1f}\n"
            process += f"         = {y1:.2f} + {(y2 - y1) * (value - x1) / (x2 - x1):.2f}\n"
            process += f"         = {y1 + (y2 - y1) * (value - x1) / (x2 - x1):.2f}"

            return f"{formula}\n{process}"


class TowerPropertiesCalculator:
    """塔工艺条件和物性数据计算器"""

    def __init__(self, data_file='benzene_toluene_results.json'):
        """初始化，加载数据"""
        with open(data_file, 'r', encoding='utf-8') as f:
            self.data = json.load(f)

        # 苯-甲苯物性常数
        self.M_benzene = 78.11  # g/mol
        self.M_toluene = 92.14  # g/mol

        # 苯的密度数据 (温度℃: 密度kg/m³)
        self.rho_benzene_dict = {
            20: 877.4, 40: 857.3, 60: 836.6, 80: 815.0,
            90: 803.9, 100: 792.5, 110: 780.3, 120: 768.9
        }

        # 甲苯的密度数据 (温度℃: 密度kg/m³)
        self.rho_toluene_dict = {
            0: 885.6, 20: 867.0, 40: 848.2, 60: 829.3,
            80: 810.0, 90: 800.2, 100: 790.3, 120: 770.0
        }

        # 苯的表面张力数据 (温度℃: 表面张力mN/m)
        self.gamma_benzene_dict = {
            20: 20.28, 40: 26.25, 60: 23.74, 80: 21.27,
            90: 20.06, 100: 18.85, 110: 17.66, 120: 16.49
        }

        # 甲苯的表面张力数据 (温度℃: 表面张力mN/m)
        self.gamma_toluene_dict = {
            20: 28.54, 40: 26.22, 60: 23.94, 80: 21.69,
            90: 20.59, 100: 19.94, 110: 18.41, 120: 17.31
        }

        # 苯的粘度数据 (温度℃: 粘度mPa·s)
        self.benzene_viscosity_dict = {
            0: 0.742, 20: 0.638, 40: 0.485, 60: 0.381,
            80: 0.308, 100: 0.255, 120: 0.215
        }

        # 甲苯的粘度数据 (温度℃: 粘度mPa·s)
        self.toluene_viscosity_dict = {
            0: 0.758, 20: 0.580, 40: 0.459, 60: 0.373,
            80: 0.311, 100: 0.264, 120: 0.228
        }

        print("=" * 60)
        print("塔工艺条件和物性数据计算")
        print("=" * 60)
        print(f"已从 {data_file} 加载精馏塔计算结果\n")

    def calculate_operating_pressure(self, top_pressure=101.325, section='rectifying'):
        """计算操作压强"""
        pressure_drop_per_tray = 0.7  # kPa/tray

        print("1. 操作压强计算")
        print("-" * 40)

        if section == 'rectifying':
            N_trays =  math.ceil(self.data['N1_actual'])
            avg_pressure = top_pressure + (N_trays / 2) * pressure_drop_per_tray
            print(f"公式: P_avg = P_top + (N1_actual / 2) × ΔP_per_tray")
            print(f"计算: P_avg = {top_pressure} + ({N_trays} / 2) × {pressure_drop_per_tray}")
            print(f"     = {top_pressure} + {N_trays / 2:.1f} × {pressure_drop_per_tray}")
        else:
            N_trays = math.ceil(self.data['N2_actual'])
            feed_stage = self.data['feed_stage_position']
            avg_pressure = top_pressure + (feed_stage + N_trays / 2) * pressure_drop_per_tray
            print(f"公式: P_avg = P_top + (feed_stage + N2_actual/2) × ΔP_per_tray")
            print(f"计算: P_avg = {top_pressure} + ({feed_stage} + {N_trays}/2) × {pressure_drop_per_tray}")
            print(f"     = {top_pressure} + ({feed_stage} + {N_trays / 2:.1f}) × {pressure_drop_per_tray}")

        print(f"结果: P_avg = {avg_pressure:.2f} kPa\n")
        return avg_pressure

    def calculate_section_temperatures(self, section='rectifying'):
        """计算精馏段或提馏段平均温度"""
        print("2. 平均温度计算")
        print("-" * 40)

        if section == 'rectifying':
            t1 = self.data['t_top']
            t2 = self.data['t_feed_bubble']
            print(f"公式: t_avg = (t_top + t_feed_bubble) / 2")
            print(f"计算: t_avg = ({t1:.1f} + {t2:.1f}) / 2")
        else:
            t1 = self.data['t_feed_bubble']
            t2 = self.data['t_bottom']
            print(f"公式: t_avg = (t_feed_bubble + t_bottom) / 2")
            print(f"计算: t_avg = ({t1:.1f} + {t2:.1f}) / 2")

        t_avg = (t1 + t2) / 2
        print(f"     = {t_avg:.1f} °C")
        print(f"温度范围: {t1:.1f} - {t2:.1f} °C\n")
        return t_avg, t1, t2

    def calculate_average_molecular_weight(self, section='rectifying'):
        """计算平均相对分子质量"""
        print("3. 平均分子量计算")
        print("-" * 40)

        if section == 'rectifying':
            x_D = self.data['x_D']
            x_F = self.data['x_F']
            x_avg = (x_D + x_F) / 2
            print(f"公式: x_avg = (x_D + x_F) / 2")
            print(f"计算: x_avg = ({x_D:.4f} + {x_F:.4f}) / 2")
        else:
            x_F = self.data['x_F']
            x_W = self.data['x_W']
            x_avg = (x_F + x_W) / 2
            print(f"公式: x_avg = (x_F + x_W) / 2")
            print(f"计算: x_avg = ({x_F:.4f} + {x_W:.4f}) / 2")

        print(f"     = {x_avg:.4f}")
        print(f"\n公式: M_avg = x_avg × M_苯 + (1 - x_avg) × M_甲苯")
        print(f"计算: M_avg = {x_avg:.4f} × {self.M_benzene} + {1 - x_avg:.4f} × {self.M_toluene}")

        M_avg = x_avg * self.M_benzene + (1 - x_avg) * self.M_toluene
        print(f"结果: M_avg = {M_avg:.2f} g/mol\n")
        return M_avg, x_avg

    def calculate_liquid_density(self, temperature, composition):
        """计算液体密度 (kg/m³)"""
        print("4. 液体密度计算")
        print("-" * 40)

        # 插值获取苯和甲苯的密度
        print("(1) 苯的密度插值:")
        print(get_interpolation_formula(temperature, self.rho_benzene_dict, "密度", "苯"))
        rho_benzene = linear_interpolate(temperature, self.rho_benzene_dict)
        print(f"结果: ρ_苯 = {rho_benzene:.1f} kg/m³\n")

        print("(2) 甲苯的密度插值:")
        print(get_interpolation_formula(temperature, self.rho_toluene_dict, "密度", "甲苯"))
        rho_toluene = linear_interpolate(temperature, self.rho_toluene_dict)
        print(f"结果: ρ_甲苯 = {rho_toluene:.1f} kg/m³\n")

        # 使用摩尔分数和摩尔体积计算混合物密度
        print("(3) 混合物密度计算:")
        print("公式: Vm_i = M_i / ρ_i / 1000")
        print("公式: Vm_mix = x × Vm_苯 + (1-x) × Vm_甲苯")
        print("公式: ρ_mix = M_mix / Vm_mix / 1000")

        # 计算苯和甲苯的摩尔体积 (m³/kmol)
        Vm_benzene = self.M_benzene / rho_benzene / 1000
        Vm_toluene = self.M_toluene / rho_toluene / 1000

        print(f"\n计算苯摩尔体积:")
        print(f"Vm_苯 = {self.M_benzene} / {rho_benzene:.1f} / 1000 = {Vm_benzene:.6f} m³/kmol")

        print(f"\n计算甲苯摩尔体积:")
        print(f"Vm_甲苯 = {self.M_toluene} / {rho_toluene:.1f} / 1000 = {Vm_toluene:.6f} m³/kmol")

        # 混合物的平均摩尔体积
        Vm_mix = composition * Vm_benzene + (1 - composition) * Vm_toluene
        print(f"\n计算混合物摩尔体积:")
        print(f"Vm_mix = {composition:.4f} × {Vm_benzene:.6f} + {1 - composition:.4f} × {Vm_toluene:.6f}")
        print(f"       = {Vm_mix:.6f} m³/kmol")

        # 混合物的平均分子量
        M_mix = composition * self.M_benzene + (1 - composition) * self.M_toluene
        # 混合物密度
        rho_mix = M_mix / Vm_mix / 1000

        print(f"\n计算混合物密度:")
        print(f"ρ_mix = {M_mix:.2f} / {Vm_mix:.6f} / 1000 = {rho_mix:.1f} kg/m³\n")

        return rho_mix, rho_benzene, rho_toluene

    def calculate_vapor_density(self, temperature, pressure, molecular_weight):
        """计算蒸汽密度 (kg/m³)"""
        print("5. 蒸汽密度计算")
        print("-" * 40)

        R = 8.314  # J/(mol·K)
        T_K = temperature + 273.15
        P_Pa = pressure * 1000

        print("公式: ρ_vapor = (P × M) / (R × T)")
        print(f"其中: P = {pressure:.2f} kPa = {P_Pa:.0f} Pa")
        print(f"      M = {molecular_weight:.2f} g/mol = {molecular_weight / 1000:.4f} kg/mol")
        print(f"      R = 8.314 J/(mol·K)")
        print(f"      T = {temperature:.1f} °C = {T_K:.2f} K")

        # 理想气体密度
        rho_vapor = (P_Pa * molecular_weight / 1000) / (R * T_K)

        print(f"\n计算: ρ_vapor = ({P_Pa:.0f} × {molecular_weight / 1000:.4f}) / (8.314 × {T_K:.2f})")
        print(f"            = {rho_vapor:.3f} kg/m³\n")

        return rho_vapor

    def calculate_surface_tension(self, temperature, composition, location=""):
        """计算液体表面张力 (mN/m)"""
        print(f"计算{location}表面张力:")
        print("-" * 30)

        # 苯的表面张力
        print("(1) 苯的表面张力插值:")
        print(get_interpolation_formula(temperature, self.gamma_benzene_dict, "表面张力", "苯"))
        sigma_benzene = linear_interpolate(temperature, self.gamma_benzene_dict)
        print(f"结果: σ_苯 = {sigma_benzene:.2f} mN/m\n")

        # 甲苯的表面张力
        print("(2) 甲苯的表面张力插值:")
        print(get_interpolation_formula(temperature, self.gamma_toluene_dict, "表面张力", "甲苯"))
        sigma_toluene = linear_interpolate(temperature, self.gamma_toluene_dict)
        print(f"结果: σ_甲苯 = {sigma_toluene:.2f} mN/m\n")

        # 使用理想混合规则
        print("(3) 混合物表面张力:")
        print("公式: σ_mix = x × σ_苯 + (1 - x) × σ_甲苯")
        print(f"计算: σ_mix = {composition:.4f} × {sigma_benzene:.2f} + {1 - composition:.4f} × {sigma_toluene:.2f}")
        print(f"          = {composition * sigma_benzene:.2f} + {(1 - composition) * sigma_toluene:.2f}")

        sigma_mix = composition * sigma_benzene + (1 - composition) * sigma_toluene
        print(f"结果: σ_{location} = {sigma_mix:.2f} mN/m\n")

        return sigma_mix, sigma_benzene, sigma_toluene

    def calculate_viscosity(self, temperature, composition, location=""):
        """计算液体粘度 (mPa·s)"""
        print(f"计算{location}粘度:")
        print("-" * 30)

        # 苯的粘度
        print("(1) 苯的粘度插值:")
        print(get_interpolation_formula(temperature, self.benzene_viscosity_dict, "粘度", "苯"))
        mu_benzene = linear_interpolate(temperature, self.benzene_viscosity_dict)
        print(f"结果: μ_苯 = {mu_benzene:.3f} mPa·s\n")

        # 甲苯的粘度
        print("(2) 甲苯的粘度插值:")
        print(get_interpolation_formula(temperature, self.toluene_viscosity_dict, "粘度", "甲苯"))
        mu_toluene = linear_interpolate(temperature, self.toluene_viscosity_dict)
        print(f"结果: μ_甲苯 = {mu_toluene:.3f} mPa·s\n")

        # 混合物粘度（对数混合规则）
        print("(3) 混合物粘度:")
        print("公式: ln(μ_mix) = x × ln(μ_苯) + (1-x) × ln(μ_甲苯)")
        print(
            f"计算: ln(μ_mix) = {composition:.4f} × ln({mu_benzene:.3f}) + {1 - composition:.4f} × ln({mu_toluene:.3f})")
        print(
            f"              = {composition:.4f} × {np.log(mu_benzene):.3f} + {1 - composition:.4f} × {np.log(mu_toluene):.3f}")
        print(f"              = {composition * np.log(mu_benzene):.3f} + {(1 - composition) * np.log(mu_toluene):.3f}")

        ln_mu_mix = composition * np.log(mu_benzene) + (1 - composition) * np.log(mu_toluene)
        mu_mix = np.exp(ln_mu_mix)

        print(f"      ln(μ_mix) = {ln_mu_mix:.3f}")
        print(f"      μ_{location} = exp({ln_mu_mix:.3f}) = {mu_mix:.3f} mPa·s\n")

        return mu_mix, mu_benzene, mu_toluene

    def calculate_average_surface_tension(self, t1, t2, x1, x2, location1="端点1", location2="端点2"):
        """计算两端点间的平均表面张力"""
        print(f"\n计算{location1}和{location2}的平均表面张力:")
        print("-" * 40)

        # 计算端点1的表面张力
        print(f"1) {location1}表面张力计算:")
        sigma1, sigma1_benzene, sigma1_toluene = self.calculate_surface_tension(t1, x1, location1)

        # 计算端点2的表面张力
        print(f"2) {location2}表面张力计算:")
        sigma2, sigma2_benzene, sigma2_toluene = self.calculate_surface_tension(t2, x2, location2)

        # 表面张力使用算术平均
        sigma_avg = (sigma1 + sigma2) / 2

        print(f"3) 平均表面张力计算:")
        print(f"公式: σ_avg = (σ_{location1} + σ_{location2}) / 2")
        print(f"计算: σ_avg = ({sigma1:.2f} + {sigma2:.2f}) / 2")
        print(f"          = {sigma1 + sigma2:.2f} / 2")
        print(f"结果: σ_avg = {sigma_avg:.2f} mN/m\n")

        return sigma_avg, sigma1, sigma2

    def calculate_average_viscosity(self, t1, t2, x1, x2, location1="端点1", location2="端点2"):
        """计算两端点间的平均黏度"""
        print(f"\n计算{location1}和{location2}的平均粘度:")
        print("-" * 40)

        # 计算端点1的黏度
        print(f"1) {location1}粘度计算:")
        mu1, mu1_benzene, mu1_toluene = self.calculate_viscosity(t1, x1, location1)

        # 计算端点2的黏度
        print(f"2) {location2}粘度计算:")
        mu2, mu2_benzene, mu2_toluene = self.calculate_viscosity(t2, x2, location2)

        # 黏度使用对数平均
        if mu1 == mu2:
            mu_avg = mu1
        elif mu1 > 0 and mu2 > 0:
            mu_avg = (mu1 - mu2) / np.log(mu1 / mu2)
        else:
            mu_avg = (mu1 + mu2) / 2

        print(f"3) 平均粘度计算:")
        if mu1 != mu2 and mu1 > 0 and mu2 > 0:
            print(f"公式: μ_avg = (μ_{location1} - μ_{location2}) / ln(μ_{location1}/μ_{location2})")
            print(f"计算: μ_avg = ({mu1:.3f} - {mu2:.3f}) / ln({mu1:.3f}/{mu2:.3f})")
            print(f"          = {mu1 - mu2:.3f} / ln({mu1 / mu2:.3f})")
            print(f"          = {mu1 - mu2:.3f} / {np.log(mu1 / mu2):.3f}")
        else:
            print(f"公式: μ_avg = (μ_{location1} + μ_{location2}) / 2")
            print(f"计算: μ_avg = ({mu1:.3f} + {mu2:.3f}) / 2")
            print(f"          = {mu1 + mu2:.3f} / 2")

        print(f"结果: μ_avg = {mu_avg:.3f} mPa·s\n")

        return mu_avg, mu1, mu2

    def calculate_flow_rates(self, section='rectifying'):
        """计算气液相流量"""
        print("8. 气液相流量计算")
        print("-" * 40)

        D = self.data['D_kmol']
        R = self.data['R']
        F = self.data['F_kmol_per_h']
        q = self.data['q_value']

        if section == 'rectifying':
            # 精馏段
            print(f"公式: L = R × D")
            print(f"计算: L = {R:.4f} × {D:.4f}")
            L = R * D
            print(f"     = {L:.2f} kmol/h\n")

            print(f"公式: V = (R + 1) × D")
            print(f"计算: V = ({R:.4f} + 1) × {D:.4f}")
            V = (R + 1) * D
            print(f"     = {V:.2f} kmol/h")
        else:
            # 提馏段
            # 先计算精馏段流量
            L_rect = R * D
            V_rect = (R + 1) * D

            print(f"首先计算精馏段流量作为基础:")
            print(f"精馏段: L = R × D = {R:.4f} × {D:.4f} = {L_rect:.2f} kmol/h")
            print(f"精馏段: V = (R+1) × D = {R + 1:.4f} × {D:.4f} = {V_rect:.2f} kmol/h\n")

            print(f"提馏段计算公式:")
            print(f"公式: L' = L + q × F")
            print(f"计算: L' = {L_rect:.4f} + {q:.4f} × {F:.4f}")
            L = L_rect + q * F
            print(f"       = {L:.2f} kmol/h\n")

            print(f"公式: V' = V - (1 - q) × F")
            print(f"计算: V' = {V_rect:.4f} - (1 - {q:.4f}) × {F:.4f}")
            V = V_rect - (1 - q) * F
            print(f"       = {V:.2f} kmol/h")

        return L, V

    def run_complete_calculation(self, section='rectifying'):
        """运行完整的工艺条件计算"""
        section_name = '精馏段' if section == 'rectifying' else '提馏段'

        print(f"\n{'=' * 60}")
        print(f"四、工艺条件计算 - {section_name}")
        print(f"{'=' * 60}")

        # 1. 操作压强
        P_avg = self.calculate_operating_pressure(section=section)

        # 2. 温度
        t_avg, t1, t2 = self.calculate_section_temperatures(section=section)

        # 3. 平均分子量
        M_avg, x_avg = self.calculate_average_molecular_weight(section=section)

        # 4. 液体密度
        rho_liquid, rho_b, rho_t = self.calculate_liquid_density(t_avg, x_avg)

        # 5. 蒸汽密度
        rho_vapor = self.calculate_vapor_density(t_avg, P_avg, M_avg)

        # 6. 液体表面张力
        if section == 'rectifying':
            print(f"\n6. 精馏段表面张力计算（塔顶到进料板）")
            sigma_avg, sigma1, sigma2 = self.calculate_average_surface_tension(
                self.data['t_top'], self.data['t_feed_bubble'],
                self.data['x_D'], self.data['x_F'],
                "塔顶", "进料板"
            )
        else:
            print(f"\n6. 提馏段表面张力计算（进料板到塔釜）")
            sigma_avg, sigma1, sigma2 = self.calculate_average_surface_tension(
                self.data['t_feed_bubble'], self.data['t_bottom'],
                self.data['x_F'], self.data['x_W'],
                "进料板", "塔釜"
            )

        # 7. 液体粘度
        if section == 'rectifying':
            print(f"\n7. 精馏段粘度计算（塔顶到进料板）")
            mu_avg, mu1, mu2 = self.calculate_average_viscosity(
                self.data['t_top'], self.data['t_feed_bubble'],
                self.data['x_D'], self.data['x_F'],
                "塔顶", "进料板"
            )
        else:
            print(f"\n7. 提馏段粘度计算（进料板到塔釜）")
            mu_avg, mu1, mu2 = self.calculate_average_viscosity(
                self.data['t_feed_bubble'], self.data['t_bottom'],
                self.data['x_F'], self.data['x_W'],
                "进料板", "塔釜"
            )

        # 8. 流量计算
        print(f"\n{'=' * 60}")
        print(f"五、气液负荷计算 - {section_name}")
        print(f"{'=' * 60}")
        L, V = self.calculate_flow_rates(section=section)

        # 9. 质量流量和体积流量


        print("\n9. 质量流量计算")
        print("-" * 40)
        print(f"公式: L_mass = L × M_avg")
        print(f"计算: L_mass = {L:.4f} × {M_avg:.2f}")
        L_mass = L * M_avg
        print(f"     = {L_mass:.2f} kg/h\n")

        print(f"公式: V_mass = V × M_avg")
        print(f"计算: V_mass = {V:.4f} × {M_avg:.2f}")
        V_mass = V * M_avg
        print(f"     = {V_mass:.2f} kg/h\n")

        print("\n10. 体积流量计算")
        print("-" * 40)
        print(f"公式: L_vol = L_mass / ρ_liquid × 1000")
        print(f"计算: L_vol = {L_mass:.2f} / {rho_liquid:.1f} × 1000")
        L_vol = L_mass / rho_liquid * 1000
        print(f"     = {L_vol:.2f} L/h\n")

        print(f"公式: V_vol = V_mass / ρ_vapor")
        print(f"计算: V_vol = {V_mass:.2f} / {rho_vapor:.4f}")
        V_vol = V_mass / rho_vapor
        print(f"     = {V_vol:.2f} m³/h\n")

        # 汇总结果
        results = {
            'section': section_name,
            'P_avg_kPa': round(P_avg, 2),
            't_avg_C': round(t_avg, 1),
            'M_avg_gmol': round(M_avg, 2),
            'x_avg': round(x_avg, 4),
            'rho_liquid_kgm3': round(rho_liquid, 1),
            'rho_vapor_kgm3': round(rho_vapor, 4),
            'sigma_mNm': round(sigma_avg, 2),
            'mu_mPas': round(mu_avg, 3),
            'sigma_endpoint1_mNm': round(sigma1, 2),
            'sigma_endpoint2_mNm': round(sigma2, 2),
            'mu_endpoint1_mPas': round(mu1, 3),
            'mu_endpoint2_mPas': round(mu2, 3),
            'L_kmolh': round(L, 4),
            'V_kmolh': round(V, 4),
            'L_mass_kgh': round(L_mass, 2),
            'V_mass_kgh': round(V_mass, 2),
            'L_vol_Lh': round(L_vol, 1),
            'V_vol_m3h': round(V_vol, 2)
        }

        return results

    def save_results_to_json(self, rectifying_results, stripping_results):
        """保存计算结果到JSON文件"""
        output_data = {
            '精馏段结果': rectifying_results,
            '提馏段结果': stripping_results
        }

        output_file = 'tower_properties_results.json'
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(output_data, f, indent=2, ensure_ascii=False)

        print(f"\n{'=' * 60}")
        print(f"计算结果已保存到: {output_file}")
        print(f"{'=' * 60}")

        file_path = os.path.abspath(output_file)
        print(f"文件保存路径: {file_path}")

        return output_file


if __name__ == "__main__":
    calculator = TowerPropertiesCalculator('benzene_toluene_results.json')

    print("\n开始计算精馏段...")
    results_rectifying = calculator.run_complete_calculation(section='rectifying')

    print("\n" + "=" * 60)
    print("开始计算提馏段...")
    print("=" * 60)

    results_stripping = calculator.run_complete_calculation(section='stripping')

    output_file = calculator.save_results_to_json(results_rectifying, results_stripping)

    print("\n" + "=" * 60)
    print("计算完成！结果汇总：")
    print("=" * 60)

    for result in [results_rectifying, results_stripping]:
        section_name = result['section']
        print(f"\n{section_name}主要结果:")
        print(f"  平均操作压强: {result['P_avg_kPa']} kPa")
        print(f"  平均温度: {result['t_avg_C']} °C")
        print(f"  平均分子量: {result['M_avg_gmol']} g/mol")
        print(f"  液体密度: {result['rho_liquid_kgm3']} kg/m³")
        print(f"  蒸汽密度: {result['rho_vapor_kgm3']} kg/m³")
        print(f"  液相流量: {result['L_kmolh']:.2f} kmol/h ({result['L_vol_Lh']:.2f} L/h)")
        print(f"  气相流量: {result['V_kmolh']:.2f} kmol/h ({result['V_vol_m3h']:.2f} m³/h)")