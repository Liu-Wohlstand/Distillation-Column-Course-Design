import json
import math
from typing import Dict, Any, Tuple, List
from datetime import datetime

# ======== 以上请勿更改 ==========

# ======== 用户可修改参数范围 ==========
'''
当然这些参数范围和默认值改不改都无所谓，只是为了显示美观和实现傻瓜操作
最终的结果肯定是取决于交互输入的参数值，具体的范围参考课程要求即可
'''
# 1. 冷却水参数范围
COOLING_WATER_IN_TEMP_RANGE = (15.0, 35.0)  # 进口温度范围 (°C)
COOLING_WATER_OUT_TEMP_RANGE = (35.0, 45.0)  # 出口温度范围 (°C)
COOLING_WATER_DELTA_T_DEFAULT = 15.0  # 默认温升 (°C)

# 2. 加热蒸汽参数范围
STEAM_PRESSURE_RANGE = (0.1, 1.0)  # 蒸汽压力范围 (MPa，表压)
STEAM_MIN_TEMP_ABOVE_BOTTOM = 10.0  # 蒸汽温度最低高于塔釜温度 (°C)

# 3. 传热系数范围 (W/(m²·K))
CONDENSER_HEAT_TRANSFER_COEFFICIENT_RANGE = (400, 800)  # 冷凝器传热系数范围
REBOILER_HEAT_TRANSFER_COEFFICIENT_RANGE = (700, 1200)  # 再沸器传热系数范围
CONDENSER_K_DEFAULT = 600  # 冷凝器默认传热系数
REBOILER_K_DEFAULT = 900  # 再沸器默认传热系数

# 4. 管道设计流速范围 (m/s)
PIPE_VELOCITY_RANGES = {
    "进料管(液体)": (0.5, 1.5),
    "回流管(液体)": (0.5, 1.5),
    "塔顶出料管(液体)": (0.5, 1.5),
    "塔釜出料管(液体)": (0.5, 1.5),
    "塔顶蒸汽管(气体)": (15.0, 25.0),
    "加热蒸汽管(气体)": (20.0, 40.0),
    "冷却水管(液体)": (1.0, 2.5)
}

# 5. 安全系数
SAFETY_FACTOR_AREA = 1.2  # 传热面积安全系数
SAFETY_FACTOR_PIPE = 1.1  # 管道尺寸安全系数


# ======== 以下请勿更改 ==========


class TowerAncillaryEquipmentDesign:
    """板式塔附属设备设计计算类"""

    def __init__(self, input_files: Dict[str, str]):
        """
        初始化类，加载输入文件数据

        Args:
            input_files: 字典，包含各个输入文件的数据
        """
        self.data = {}
        for file_name, content in input_files.items():
            if isinstance(content, str):
                self.data[file_name] = json.loads(content)
            else:
                self.data[file_name] = content

        # 提取关键参数
        self.extract_key_parameters()

        # 设置计算精度
        self.precision = 4

        # 初始化加热蒸汽和冷却水流量（将在后续设计中计算）
        self.steam_flow_rate = None
        self.cooling_water_flow_rate = None
        self.steam_pressure = None
        self.steam_temp = None

    def extract_key_parameters(self):
        """从加载的数据中提取关键参数"""
        print("\n" + "=" * 80)
        print("关键参数提取")
        print("=" * 80)

        try:
            # 从苯-甲苯结果文件
            benzene_data = self.data.get('benzene_toluene_results.json', {})
            self.F_mass = benzene_data.get('F_mass', 7500.0)  # kg/h
            self.D_mass = benzene_data.get('D_mass', 1906.78)  # kg/h
            self.W_mass = benzene_data.get('W_mass', 5593.22)  # kg/h
            self.t_top = benzene_data.get('t_top', 81.73)  # °C
            self.t_bottom = benzene_data.get('t_bottom', 108.83)  # °C
            self.x_D = benzene_data.get('x_D', 0.923)  # 苯的摩尔分数
            self.x_W = benzene_data.get('x_W', 0.029)  # 苯的摩尔分数
            self.R = benzene_data.get('R', 2.8525)  # 回流比

            print(f"1. 物料衡算:")
            print(f"   进料流量 F = {self.F_mass:.2f} kg/h")
            print(f"   馏出液流量 D = {self.D_mass:.2f} kg/h")
            print(f"   塔釜液流量 W = {self.W_mass:.2f} kg/h")
            print(f"   塔顶温度 t_top = {self.t_top:.2f} °C")
            print(f"   塔釜温度 t_bottom = {self.t_bottom:.2f} °C")
            print(f"   回流比 R = {self.R:.4f}")

            # 从塔特性结果文件
            tower_props = self.data.get('tower_properties_results.json', {})
            rectifying = tower_props.get('精馏段结果', {})
            stripping = tower_props.get('提馏段结果', {})

            self.V_mass_rect = rectifying.get('V_mass_kgh', 7762.66)  # kg/h 精馏段气相质量流量
            self.L_mass_rect = rectifying.get('L_mass_kgh', 5747.71)  # kg/h 精馏段液相质量流量
            self.rho_v_rect = rectifying.get('rho_vapor_kgm3', 2.9271)  # kg/m³ 精馏段气相密度
            self.rho_l_rect = rectifying.get('rho_liquid_kgm3', 802.6)  # kg/m³ 精馏段液相密度
            self.t_avg_rect = rectifying.get('t_avg_C', 89.7)  # °C 精馏段平均温度

            self.V_mass_strip = stripping.get('V_mass_kgh', 10378.12)  # kg/h 提馏段气相质量流量
            self.L_mass_strip = stripping.get('L_mass_kgh', 15863.17)  # kg/h 提馏段液相质量流量
            self.rho_v_strip = stripping.get('rho_vapor_kgm3', 3.1337)  # kg/m³ 提馏段气相密度
            self.rho_l_strip = stripping.get('rho_liquid_kgm3', 787.2)  # kg/m³ 提馏段液相密度
            self.t_avg_strip = stripping.get('t_avg_C', 103.3)  # °C 提馏段平均温度

            print(f"\n2. 精馏段参数:")
            print(f"   气相流量 V_rect = {self.V_mass_rect:.2f} kg/h")
            print(f"   液相密度 ρ_l_rect = {self.rho_l_rect:.2f} kg/m³")
            print(f"   气相密度 ρ_v_rect = {self.rho_v_rect:.4f} kg/m³")

            print(f"\n3. 提馏段参数:")
            print(f"   气相流量 V_strip = {self.V_mass_strip:.2f} kg/h")
            print(f"   液相密度 ρ_l_strip = {self.rho_l_strip:.2f} kg/m³")
            print(f"   气相密度 ρ_v_strip = {self.rho_v_strip:.4f} kg/m³")

            # 从塔径设计文件
            tower_diameter = self.data.get('tower_diameter_approximate_optimized.json', {})
            final_design = tower_diameter.get('最终设计结果', {})
            self.D_tower = final_design.get('最终塔径_m', 1.0)  # m 塔径

            print(f"\n4. 塔体参数:")
            print(f"   塔径 D = {self.D_tower:.3f} m")

            print("\n关键参数提取完成！")

        except Exception as e:
            print(f"参数提取错误: {e}")
            raise

    def design_condenser(self, water_in_temp: float, water_out_temp: float) -> Dict[str, Any]:
        """
        设计塔顶冷凝器

        Args:
            water_in_temp: 冷却水进口温度 (°C)
            water_out_temp: 冷却水出口温度 (°C)

        Returns:
            冷凝器设计结果字典
        """
        print("\n" + "=" * 80)
        print("塔顶冷凝器设计 - 详细计算过程")
        print("=" * 80)

        # 检查冷却水参数是否在范围内
        if not (COOLING_WATER_IN_TEMP_RANGE[0] <= water_in_temp <= COOLING_WATER_IN_TEMP_RANGE[1]):
            print(f"警告: 冷却水进口温度 {water_in_temp}°C 超出推荐范围 {COOLING_WATER_IN_TEMP_RANGE}")
            water_in_temp = max(COOLING_WATER_IN_TEMP_RANGE[0],
                                min(water_in_temp, COOLING_WATER_IN_TEMP_RANGE[1]))
            print(f"调整为: {water_in_temp}°C")

        if not (COOLING_WATER_OUT_TEMP_RANGE[0] <= water_out_temp <= COOLING_WATER_OUT_TEMP_RANGE[1]):
            print(f"警告: 冷却水出口温度 {water_out_temp}°C 超出推荐范围 {COOLING_WATER_OUT_TEMP_RANGE}")
            water_out_temp = max(COOLING_WATER_OUT_TEMP_RANGE[0],
                                 min(water_out_temp, COOLING_WATER_OUT_TEMP_RANGE[1]))
            print(f"调整为: {water_out_temp}°C")

        # 1. 热负荷计算
        print("\n1. 热负荷计算:")
        print("   公式: V_mass = (R + 1) × D_mass")
        print(f"         = ({self.R:.4f} + 1) × {self.D_mass:.2f} kg/h")
        V_mass = (self.R + 1) * self.D_mass  # kg/h
        print(f"         = {V_mass:.2f} kg/h")

        print("\n   公式: Q = V_mass × r_avg")
        print("   其中 r_avg = x_D × r_benzene + (1 - x_D) × r_toluene")
        print(f"        = {self.x_D:.4f} × 393.9 + (1 - {self.x_D:.4f}) × 363.0")

        # 苯和甲苯的汽化潜热 (kJ/kg)
        latent_heat_benzene = 393.9  # kJ/kg @ 80°C
        latent_heat_toluene = 363.0  # kJ/kg @ 110°C

        # 计算平均汽化潜热 (根据摩尔分数加权)
        avg_latent_heat = (self.x_D * latent_heat_benzene +
                           (1 - self.x_D) * latent_heat_toluene)  # kJ/kg
        print(f"        = {self.x_D * latent_heat_benzene:.2f} + {(1 - self.x_D) * latent_heat_toluene:.2f}")
        print(f"        = {avg_latent_heat:.2f} kJ/kg")

        # 热负荷 Q = V * r
        Q = V_mass * avg_latent_heat  # kJ/h
        print(f"\n   Q = {V_mass:.2f} kg/h × {avg_latent_heat:.2f} kJ/kg")
        print(f"     = {Q:.2f} kJ/h")

        Q_kW = Q / 3600  # kW
        print(f"     = {Q:.2f} kJ/h ÷ 3600 s/h")
        print(f"     = {Q_kW:.2f} kW")

        # 2. 冷却水用量计算
        print("\n2. 冷却水用量计算:")
        print("   公式: W_water = Q ÷ [c_p_water × (t_out - t_in)]")
        cp_water = 4.187  # kJ/(kg·°C)
        delta_T_water = water_out_temp - water_in_temp

        if delta_T_water <= 0:
            print(f"   警告: 冷却水出口温度({water_out_temp}°C)必须高于进口温度({water_in_temp}°C)!")
            delta_T_water = COOLING_WATER_DELTA_T_DEFAULT  # 使用配置的默认值
            water_out_temp = water_in_temp + delta_T_water
            print(f"   自动调整为: 出口温度 = {water_out_temp}°C，温差 = {delta_T_water}°C")

        print(f"   c_p_water = {cp_water} kJ/(kg·°C)")
        print(f"   ΔT_water = {water_out_temp} - {water_in_temp} = {delta_T_water} °C")

        water_flow_rate = Q / (cp_water * delta_T_water)  # kg/h
        print(f"\n   W_water = {Q:.2f} kJ/h ÷ [{cp_water} kJ/(kg·°C) × {delta_T_water} °C]")
        print(f"           = {Q:.2f} ÷ {cp_water * delta_T_water:.2f}")
        print(f"           = {water_flow_rate:.2f} kg/h")

        water_flow_rate_m3h = water_flow_rate / 1000  # m³/h (水的密度为1000 kg/m³)
        print(f"           = {water_flow_rate:.2f} kg/h ÷ 1000 kg/m³")
        print(f"           = {water_flow_rate_m3h:.2f} m³/h")

        # 保存冷却水流量供管道计算使用
        self.cooling_water_flow_rate = water_flow_rate

        # 3. 平均温差计算
        print("\n3. 对数平均温差计算:")
        print("   公式: ΔT_lm = (ΔT1 - ΔT2) ÷ ln(ΔT1/ΔT2)")
        print("   其中: ΔT1 = T_hot_in - T_cold_in")
        print("          ΔT2 = T_hot_out - T_cold_out")

        T_hot = self.t_top  # °C 塔顶蒸汽温度
        T_cold_in = water_in_temp
        T_cold_out = water_out_temp

        delta_T1 = T_hot - T_cold_in
        delta_T2 = T_hot - T_cold_out

        print(f"\n   ΔT1 = {T_hot} - {T_cold_in} = {delta_T1:.2f} °C")
        print(f"   ΔT2 = {T_hot} - {T_cold_out} = {delta_T2:.2f} °C")

        if abs(delta_T1 - delta_T2) < 1e-6:
            delta_T_lm = delta_T1
            print(f"\n   ΔT1 ≈ ΔT2，取 ΔT_lm = ΔT1 = {delta_T_lm:.2f} °C")
        else:
            delta_T_lm = (delta_T1 - delta_T2) / math.log(delta_T1 / delta_T2)
            print(f"\n   ΔT_lm = ({delta_T1:.2f} - {delta_T2:.2f}) ÷ ln({delta_T1:.2f}/{delta_T2:.2f})")
            print(f"         = {delta_T1 - delta_T2:.2f} ÷ {math.log(delta_T1 / delta_T2):.4f}")
            print(f"         = {delta_T_lm:.2f} °C")

        # 4. 传热面积计算
        print("\n4. 传热面积计算:")
        print("   公式: A = Q ÷ (K × ΔT_lm)")
        print(f"   其中 K 为传热系数，推荐范围 {CONDENSER_HEAT_TRANSFER_COEFFICIENT_RANGE} W/(m²·K)")

        K = CONDENSER_K_DEFAULT  # W/(m²·K)
        print(f"   采用默认值 K = {K} W/(m²·K) = {K / 1000:.3f} kW/(m²·K)")

        # 考虑安全系数
        A_calc = Q_W / (K * delta_T_lm)  # m²
        A = A_calc * SAFETY_FACTOR_AREA  # 考虑安全系数
        print(f"\n   Q = {Q_kW:.2f} kW = {Q_W:.0f} W")
        print(f"   ΔT_lm = {delta_T_lm:.2f} K")
        print(f"\n   计算面积: A_calc = {Q_W:.0f} W ÷ [{K} W/(m²·K) × {delta_T_lm:.2f} K]")
        print(f"                = {Q_W:.0f} ÷ {K * delta_T_lm:.2f}")
        print(f"                = {A_calc:.2f} m²")
        print(f"   考虑安全系数 {SAFETY_FACTOR_AREA}: A = {A_calc:.2f} × {SAFETY_FACTOR_AREA} = {A:.2f} m²")

        # 5. 冷凝器选型建议
        print("\n5. 选型建议:")
        print(f"   根据计算传热面积 A = {A:.2f} m²，推荐选用管壳式冷凝器")

        # 返回结果
        result = {
            "设计参数": {
                "热负荷_kW": round(Q_kW, 2),
                "冷却水进口温度_C": water_in_temp,
                "冷却水出口温度_C": water_out_temp,
                "冷却水用量_kg_h": round(water_flow_rate, 2),
                "对数平均温差_C": round(delta_T_lm, 2),
                "计算传热面积_m2": round(A_calc, 2),
                "设计传热面积_m2": round(A, 2),
                "设计传热系数_W_m2K": K,
                "安全系数": SAFETY_FACTOR_AREA
            },
            "计算过程": {
                "塔顶蒸汽流量_kg_h": round(V_mass, 2),
                "平均汽化潜热_kJ_kg": round(avg_latent_heat, 2),
                "热负荷_kJ_h": round(Q, 2),
                "冷却水温差_C": round(delta_T_water, 2),
                "ΔT1_C": round(delta_T1, 2),
                "ΔT2_C": round(delta_T2, 2)
            },
            "选型建议": {
                "推荐类型": "管壳式冷凝器",
                "推荐型号": f"BEM型，换热面积≥{math.ceil(A)}m²",
                "制造商": "兰州兰石换热设备有限责任公司"
            }
        }

        return result

    def design_reboiler(self, steam_pressure: float = 0.3) -> Dict[str, Any]:
        """
        设计塔釜再沸器

        Args:
            steam_pressure: 加热蒸汽压力 (MPa，表压)

        Returns:
            再沸器设计结果字典
        """
        print("\n" + "=" * 80)
        print("塔釜再沸器设计 - 详细计算过程")
        print("=" * 80)

        # 检查蒸汽压力是否在范围内
        if not (STEAM_PRESSURE_RANGE[0] <= steam_pressure <= STEAM_PRESSURE_RANGE[1]):
            print(f"警告: 蒸汽压力 {steam_pressure} MPa 超出推荐范围 {STEAM_PRESSURE_RANGE}")
            steam_pressure = max(STEAM_PRESSURE_RANGE[0],
                                 min(steam_pressure, STEAM_PRESSURE_RANGE[1]))
            print(f"调整为: {steam_pressure} MPa")

        # 保存蒸汽参数供管道计算使用
        self.steam_pressure = steam_pressure

        # 根据蒸汽压力确定蒸汽温度
        # 饱和蒸汽温度与压力的关系（近似公式）
        # 对于水蒸气：t_sat ≈ 100 + 6.14 × (P - 0.1013) / 0.1
        P_abs = steam_pressure + 0.1013  # MPa，绝对压力
        steam_temp = 100 + 6.14 * (P_abs - 0.1013) / 0.1  # °C
        self.steam_temp = steam_temp

        print(f"\n1. 加热蒸汽参数确定:")
        print(f"   输入蒸汽压力: {steam_pressure} MPa (表压)")
        print(f"   绝对压力: P_abs = {steam_pressure} + 0.1013 = {P_abs:.4f} MPa")
        print(f"   饱和蒸汽温度: t_sat = 100 + 6.14 × ({P_abs:.4f} - 0.1013) ÷ 0.1")
        print(f"                  = {steam_temp:.2f} °C")

        if steam_temp <= self.t_bottom + STEAM_MIN_TEMP_ABOVE_BOTTOM:
            print(
                f"\n   警告: 蒸汽饱和温度({steam_temp:.2f}°C)低于要求温度({self.t_bottom + STEAM_MIN_TEMP_ABOVE_BOTTOM:.2f}°C)!")
            print(f"   要求蒸汽温度至少比塔釜温度高 {STEAM_MIN_TEMP_ABOVE_BOTTOM}°C")

            # 计算所需最低蒸汽压力
            required_temp = self.t_bottom + STEAM_MIN_TEMP_ABOVE_BOTTOM
            required_P_abs = 0.1013 + (required_temp - 100) * 0.1 / 6.14
            required_P_gauge = required_P_abs - 0.1013

            print(f"   建议蒸汽压力至少为: {required_P_gauge:.3f} MPa (表压)")
            steam_pressure = max(steam_pressure, required_P_gauge)
            P_abs = steam_pressure + 0.1013
            steam_temp = 100 + 6.14 * (P_abs - 0.1013) / 0.1
            print(f"   采用蒸汽压力: {steam_pressure:.3f} MPa，饱和温度: {steam_temp:.2f} °C")

            # 更新保存的参数
            self.steam_pressure = steam_pressure
            self.steam_temp = steam_temp

        # 1. 热负荷计算
        print("\n2. 热负荷计算:")
        print("   公式: Q = V_mass_strip × r_avg")
        print("   其中 r_avg = x_W × r_benzene + (1 - x_W) × r_toluene")
        print(f"        = {self.x_W:.4f} × 360.0 + (1 - {self.x_W:.4f}) × 335.0")

        # 苯和甲苯的汽化潜热 (kJ/kg) @ 塔釜温度
        latent_heat_benzene = 360.0  # kJ/kg @ 110°C
        latent_heat_toluene = 335.0  # kJ/kg @ 110°C

        avg_latent_heat = (self.x_W * latent_heat_benzene +
                           (1 - self.x_W) * latent_heat_toluene)  # kJ/kg
        print(f"        = {self.x_W * latent_heat_benzene:.2f} + {(1 - self.x_W) * latent_heat_toluene:.2f}")
        print(f"        = {avg_latent_heat:.2f} kJ/kg")

        Q = self.V_mass_strip * avg_latent_heat  # kJ/h
        print(f"\n   Q = {self.V_mass_strip:.2f} kg/h × {avg_latent_heat:.2f} kJ/kg")
        print(f"     = {Q:.2f} kJ/h")

        Q_kW = Q / 3600  # kW
        print(f"     = {Q:.2f} kJ/h ÷ 3600 s/h")
        print(f"     = {Q_kW:.2f} kW")

        # 2. 加热蒸汽用量计算
        print("\n3. 加热蒸汽用量计算:")
        print("   公式: W_steam = Q ÷ r_steam")
        print("   其中 r_steam 为蒸汽的汽化潜热")

        # 蒸汽的汽化潜热随温度变化，近似公式
        r_steam = 2257 - 2.5 * (steam_temp - 100)  # kJ/kg
        print(f"   r_steam = 2257 - 2.5 × ({steam_temp:.1f} - 100)")
        print(f"           = {r_steam:.2f} kJ/kg")

        steam_flow_rate = Q / r_steam  # kg/h
        print(f"\n   W_steam = {Q:.2f} kJ/h ÷ {r_steam:.2f} kJ/kg")
        print(f"           = {steam_flow_rate:.2f} kg/h")

        # 保存蒸汽流量供管道计算使用
        self.steam_flow_rate = steam_flow_rate

        # 3. 平均温差计算
        print("\n4. 对数平均温差计算:")
        print("   公式: ΔT_lm = ΔT (对于蒸汽温度恒定的情况)")
        print("   因为蒸汽在饱和温度下冷凝，温度不变")

        delta_T = steam_temp - self.t_bottom
        print(f"\n   ΔT = 蒸汽饱和温度 - 塔釜温度")
        print(f"      = {steam_temp:.2f} °C - {self.t_bottom:.2f} °C")
        print(f"      = {delta_T:.2f} °C")

        delta_T_lm = delta_T  # 对于恒温蒸汽

        # 4. 传热面积计算
        print("\n5. 传热面积计算:")
        print("   公式: A = Q ÷ (K × ΔT_lm)")
        print(f"   其中 K 为传热系数，推荐范围 {REBOILER_HEAT_TRANSFER_COEFFICIENT_RANGE} W/(m²·K)")

        K = REBOILER_K_DEFAULT  # W/(m²·K)
        print(f"   采用默认值 K = {K} W/(m²·K) = {K / 1000:.3f} kW/(m²·K)")

        # 单位转换
        Q_W = Q_kW * 1000  # W
        print(f"\n   Q = {Q_kW:.2f} kW = {Q_W:.0f} W")
        print(f"   ΔT_lm = {delta_T_lm:.2f} K")

        # 考虑安全系数
        A_calc = Q_W / (K * delta_T_lm)  # m²
        A = A_calc * SAFETY_FACTOR_AREA  # 考虑安全系数
        print(f"\n   计算面积: A_calc = {Q_W:.0f} W ÷ [{K} W/(m²·K) × {delta_T_lm:.2f} K]")
        print(f"                = {Q_W:.0f} ÷ {K * delta_T_lm:.2f}")
        print(f"                = {A_calc:.2f} m²")
        print(f"   考虑安全系数 {SAFETY_FACTOR_AREA}: A = {A_calc:.2f} × {SAFETY_FACTOR_AREA} = {A:.2f} m²")

        # 5. 再沸器选型建议
        print("\n6. 选型建议:")
        print(f"   根据计算传热面积 A = {A:.2f} m²，推荐选用釜式再沸器")

        # 返回结果
        result = {
            "设计参数": {
                "热负荷_kW": round(Q_kW, 2),
                "加热蒸汽压力_MPa": round(steam_pressure, 3),
                "加热蒸汽温度_C": round(steam_temp, 2),
                "加热蒸汽用量_kg_h": round(steam_flow_rate, 2),
                "对数平均温差_C": round(delta_T_lm, 2),
                "计算传热面积_m2": round(A_calc, 2),
                "设计传热面积_m2": round(A, 2),
                "设计传热系数_W_m2K": K,
                "安全系数": SAFETY_FACTOR_AREA
            },
            "计算过程": {
                "蒸汽汽化潜热_kJ_kg": round(r_steam, 2),
                "平均汽化潜热_kJ_kg": round(avg_latent_heat, 2),
                "热负荷_kJ_h": round(Q, 2),
                "温差_C": round(delta_T, 2)
            },
            "选型建议": {
                "推荐类型": "釜式再沸器",
                "推荐型号": f"KET型，换热面积≥{math.ceil(A)}m²",
                "制造商": "上海锅炉厂有限公司",
                "操作压力": f"{steam_pressure:.3f} MPa (表压)"
            }
        }

        return result

    def calculate_pipe_sizes(self, fluid_type: str = "直管") -> Dict[str, Any]:
        """
        计算各接管尺寸

        Args:
            fluid_type: 流体类型，可选"直管"或其他类型

        Returns:
            各接管尺寸计算结果字典
        """
        print("\n" + "=" * 80)
        print(f"接管尺寸计算 - 详细计算过程 ({fluid_type})")
        print("=" * 80)

        # 使用全局定义的流速范围
        velocity_ranges = PIPE_VELOCITY_RANGES

        # 取中间值作为设计流速
        design_velocities = {
            key: (v_range[0] + v_range[1]) / 2 for key, v_range in velocity_ranges.items()
        }

        results = {}

        print(f"\n设计流速范围:")
        for pipe_type, (v_min, v_max) in velocity_ranges.items():
            v_design = design_velocities[pipe_type]
            print(f"  {pipe_type}: {v_min:.1f}-{v_max:.1f} m/s，取 {v_design:.2f} m/s")

        # 1. 进料管
        print("\n1. 进料管计算:")
        print("   公式: d = √[4 × V_vol_s ÷ (π × v)]")

        # 假设进料为泡点进料，密度取平均值
        rho_feed = (800 + 850) / 2  # kg/m³，苯-甲苯混合液的平均密度
        F_vol = self.F_mass / rho_feed  # m³/h
        F_vol_s = F_vol / 3600  # m³/s
        v = design_velocities["进料管(液体)"]

        print(f"\n   进料流量: F = {self.F_mass:.2f} kg/h")
        print(f"   进料密度: ρ = {rho_feed:.0f} kg/m³")
        print(f"   体积流量: V_vol = F ÷ ρ = {self.F_mass:.2f} ÷ {rho_feed:.0f} = {F_vol:.4f} m³/h")
        print(f"              V_vol_s = {F_vol:.4f} ÷ 3600 = {F_vol_s:.6f} m³/s")
        print(f"   设计流速: v = {v:.2f} m/s")

        d = math.sqrt(4 * F_vol_s / (math.pi * v))
        print(f"\n   计算内径: d = √[4 × {F_vol_s:.6f} ÷ (π × {v:.2f})]")
        print(f"             = √[{4 * F_vol_s:.6f} ÷ {math.pi * v:.4f}]")
        print(f"             = √[{4 * F_vol_s / (math.pi * v):.6f}]")
        print(f"             = {d:.4f} m = {d * 1000:.1f} mm")

        dn = self.get_standard_dn(d)
        print(f"   公称直径: DN{dn}")

        results["进料管"] = {
            "质量流量_kg_h": round(self.F_mass, 2),
            "体积流量_m3_s": round(F_vol_s, 6),
            "设计流速_m_s": v,
            "计算内径_m": round(d, 4),
            "计算内径_mm": round(d * 1000, 1),
            "公称直径_DN": dn,
            "推荐规格": f"DN{dn} 不锈钢管"
        }

        # 2. 回流管
        print("\n2. 回流管计算:")
        print("   公式: d = √[4 × V_vol_s ÷ (π × v)]")

        # 回流液量 L = R * D
        L_mass = self.R * self.D_mass  # kg/h
        L_vol = L_mass / self.rho_l_rect  # m³/h
        L_vol_s = L_vol / 3600  # m³/s
        v = design_velocities["回流管(液体)"]

        print(f"\n   回流液量: L = R × D = {self.R:.4f} × {self.D_mass:.2f} = {L_mass:.2f} kg/h")
        print(f"   液相密度: ρ = {self.rho_l_rect:.1f} kg/m³")
        print(f"   体积流量: V_vol = L ÷ ρ = {L_mass:.2f} ÷ {self.rho_l_rect:.1f} = {L_vol:.4f} m³/h")
        print(f"              V_vol_s = {L_vol:.4f} ÷ 3600 = {L_vol_s:.6f} m³/s")
        print(f"   设计流速: v = {v:.2f} m/s")

        d = math.sqrt(4 * L_vol_s / (math.pi * v))
        print(f"\n   计算内径: d = √[4 × {L_vol_s:.6f} ÷ (π × {v:.2f})]")
        print(f"             = {d:.4f} m = {d * 1000:.1f} mm")

        dn = self.get_standard_dn(d)
        print(f"   公称直径: DN{dn}")

        results["回流管"] = {
            "质量流量_kg_h": round(L_mass, 2),
            "体积流量_m3_s": round(L_vol_s, 6),
            "设计流速_m_s": v,
            "计算内径_m": round(d, 4),
            "计算内径_mm": round(d * 1000, 1),
            "公称直径_DN": dn,
            "推荐规格": f"DN{dn} 不锈钢管"
        }

        # 3. 塔顶蒸汽管
        print("\n3. 塔顶蒸汽管计算:")
        print("   公式: d = √[4 × V_vol_s ÷ (π × v)]")

        # 塔顶蒸汽量 V = (R+1)*D
        V_mass = (self.R + 1) * self.D_mass  # kg/h
        V_vol = V_mass / self.rho_v_rect  # m³/h
        V_vol_s = V_vol / 3600  # m³/s
        v = design_velocities["塔顶蒸汽管(气体)"]

        print(f"\n   塔顶蒸汽量: V = (R+1) × D = ({self.R:.4f}+1) × {self.D_mass:.2f} = {V_mass:.2f} kg/h")
        print(f"   气相密度: ρ = {self.rho_v_rect:.4f} kg/m³")
        print(f"   体积流量: V_vol = V ÷ ρ = {V_mass:.2f} ÷ {self.rho_v_rect:.4f} = {V_vol:.4f} m³/h")
        print(f"              V_vol_s = {V_vol:.4f} ÷ 3600 = {V_vol_s:.6f} m³/s")
        print(f"   设计流速: v = {v:.2f} m/s")

        d = math.sqrt(4 * V_vol_s / (math.pi * v))
        print(f"\n   计算内径: d = √[4 × {V_vol_s:.6f} ÷ (π × {v:.2f})]")
        print(f"             = {d:.4f} m = {d * 1000:.1f} mm")

        dn = self.get_standard_dn(d)
        print(f"   公称直径: DN{dn}")

        results["塔顶蒸汽管"] = {
            "质量流量_kg_h": round(V_mass, 2),
            "体积流量_m3_s": round(V_vol_s, 6),
            "设计流速_m_s": v,
            "计算内径_m": round(d, 4),
            "计算内径_mm": round(d * 1000, 1),
            "公称直径_DN": dn,
            "推荐规格": f"DN{dn} 不锈钢管，保温层厚度50mm"
        }

        # 4. 塔顶出料管
        print("\n4. 塔顶出料管计算:")
        print("   公式: d = √[4 × V_vol_s ÷ (π × v)]")

        D_vol = self.D_mass / self.rho_l_rect  # m³/h
        D_vol_s = D_vol / 3600  # m³/s
        v = design_velocities["塔顶出料管(液体)"]

        print(f"\n   塔顶出料量: D = {self.D_mass:.2f} kg/h")
        print(f"   液相密度: ρ = {self.rho_l_rect:.1f} kg/m³")
        print(f"   体积流量: V_vol = D ÷ ρ = {self.D_mass:.2f} ÷ {self.rho_l_rect:.1f} = {D_vol:.4f} m³/h")
        print(f"              V_vol_s = {D_vol:.4f} ÷ 3600 = {D_vol_s:.6f} m³/s")
        print(f"   设计流速: v = {v:.2f} m/s")

        d = math.sqrt(4 * D_vol_s / (math.pi * v))
        print(f"\n   计算内径: d = √[4 × {D_vol_s:.6f} ÷ (π × {v:.2f})]")
        print(f"             = {d:.4f} m = {d * 1000:.1f} mm")

        dn = self.get_standard_dn(d)
        print(f"   公称直径: DN{dn}")

        results["塔顶出料管"] = {
            "质量流量_kg_h": round(self.D_mass, 2),
            "体积流量_m3_s": round(D_vol_s, 6),
            "设计流速_m_s": v,
            "计算内径_m": round(d, 4),
            "计算内径_mm": round(d * 1000, 1),
            "公称直径_DN": dn,
            "推荐规格": f"DN{dn} 不锈钢管"
        }

        # 5. 塔釜出料管
        print("\n5. 塔釜出料管计算:")
        print("   公式: d = √[4 × V_vol_s ÷ (π × v)]")

        W_vol = self.W_mass / self.rho_l_strip  # m³/h
        W_vol_s = W_vol / 3600  # m³/s
        v = design_velocities["塔釜出料管(液体)"]

        print(f"\n   塔釜出料量: W = {self.W_mass:.2f} kg/h")
        print(f"   液相密度: ρ = {self.rho_l_strip:.1f} kg/m³")
        print(f"   体积流量: V_vol = W ÷ ρ = {self.W_mass:.2f} ÷ {self.rho_l_strip:.1f} = {W_vol:.4f} m³/h")
        print(f"              V_vol_s = {W_vol:.4f} ÷ 3600 = {W_vol_s:.6f} m³/s")
        print(f"   设计流速: v = {v:.2f} m/s")

        d = math.sqrt(4 * W_vol_s / (math.pi * v))
        print(f"\n   计算内径: d = √[4 × {W_vol_s:.6f} ÷ (π × {v:.2f})]")
        print(f"             = {d:.4f} m = {d * 1000:.1f} mm")

        dn = self.get_standard_dn(d)
        print(f"   公称直径: DN{dn}")

        results["塔釜出料管"] = {
            "质量流量_kg_h": round(self.W_mass, 2),
            "体积流量_m3_s": round(W_vol_s, 6),
            "设计流速_m_s": v,
            "计算内径_m": round(d, 4),
            "计算内径_mm": round(d * 1000, 1),
            "公称直径_DN": dn,
            "推荐规格": f"DN{dn} 不锈钢管"
        }

        # 6. 加热蒸汽管（新增）
        print("\n6. 加热蒸汽管计算:")
        print("   公式: d = √[4 × V_vol_s ÷ (π × v)]")

        if self.steam_flow_rate is not None:
            # 计算蒸汽密度（根据压力和温度）
            # 对于饱和蒸汽，密度可近似计算：ρ = P × M / (R × T)
            # 其中：P为绝对压力(Pa)，M为水的摩尔质量(0.018 kg/mol)，R为气体常数(8.314 J/(mol·K))，T为绝对温度(K)

            P_abs = (self.steam_pressure + 0.1013) * 1e6  # 转换为Pa
            T_K = self.steam_temp + 273.15  # 转换为K
            M = 0.018  # kg/mol
            R = 8.314  # J/(mol·K)

            rho_steam = (P_abs * M) / (R * T_K)  # kg/m³

            steam_vol = self.steam_flow_rate / rho_steam  # m³/h
            steam_vol_s = steam_vol / 3600  # m³/s
            v = design_velocities["加热蒸汽管(气体)"]

            print(f"\n   加热蒸汽流量: W_steam = {self.steam_flow_rate:.2f} kg/h")
            print(
                f"   蒸汽绝对压力: P_abs = {self.steam_pressure + 0.1013:.4f} MPa = {(self.steam_pressure + 0.1013) * 1e6:.0f} Pa")
            print(f"   蒸汽温度: T = {self.steam_temp:.2f} °C = {T_K:.2f} K")
            print(f"   蒸汽密度: ρ = P × M ÷ (R × T)")
            print(f"              = {P_abs:.0f} × {M} ÷ ({R} × {T_K:.2f})")
            print(f"              = {rho_steam:.4f} kg/m³")
            print(
                f"   体积流量: V_vol = W_steam ÷ ρ = {self.steam_flow_rate:.2f} ÷ {rho_steam:.4f} = {steam_vol:.4f} m³/h")
            print(f"              V_vol_s = {steam_vol:.4f} ÷ 3600 = {steam_vol_s:.6f} m³/s")
            print(f"   设计流速: v = {v:.2f} m/s")

            d = math.sqrt(4 * steam_vol_s / (math.pi * v))
            print(f"\n   计算内径: d = √[4 × {steam_vol_s:.6f} ÷ (π × {v:.2f})]")
            print(f"             = {d:.4f} m = {d * 1000:.1f} mm")

            dn = self.get_standard_dn(d)
            print(f"   公称直径: DN{dn}")

            results["加热蒸汽管"] = {
                "质量流量_kg_h": round(self.steam_flow_rate, 2),
                "蒸汽压力_MPa": round(self.steam_pressure, 3),
                "蒸汽温度_C": round(self.steam_temp, 2),
                "蒸汽密度_kg_m3": round(rho_steam, 4),
                "体积流量_m3_s": round(steam_vol_s, 6),
                "设计流速_m_s": v,
                "计算内径_m": round(d, 4),
                "计算内径_mm": round(d * 1000, 1),
                "公称直径_DN": dn,
                "推荐规格": f"DN{dn} 不锈钢管，保温层厚度100mm，耐压{self.steam_pressure + 0.2:.2f}MPa"
            }
        else:
            print("   警告: 蒸汽流量未计算，请先设计再沸器")
            results["加热蒸汽管"] = {
                "状态": "未计算",
                "原因": "请先设计再沸器以获取蒸汽流量参数"
            }

        # 7. 冷却水管（新增）
        print("\n7. 冷却水管计算:")
        print("   公式: d = √[4 × V_vol_s ÷ (π × v)]")

        if self.cooling_water_flow_rate is not None:
            # 冷却水密度取1000 kg/m³
            rho_water = 1000  # kg/m³
            water_vol = self.cooling_water_flow_rate / rho_water  # m³/h
            water_vol_s = water_vol / 3600  # m³/s
            v = design_velocities["冷却水管(液体)"]

            print(f"\n   冷却水流量: W_water = {self.cooling_water_flow_rate:.2f} kg/h")
            print(f"   水的密度: ρ = {rho_water} kg/m³")
            print(
                f"   体积流量: V_vol = W_water ÷ ρ = {self.cooling_water_flow_rate:.2f} ÷ {rho_water} = {water_vol:.4f} m³/h")
            print(f"              V_vol_s = {water_vol:.4f} ÷ 3600 = {water_vol_s:.6f} m³/s")
            print(f"   设计流速: v = {v:.2f} m/s")

            d = math.sqrt(4 * water_vol_s / (math.pi * v))
            print(f"\n   计算内径: d = √[4 × {water_vol_s:.6f} ÷ (π × {v:.2f})]")
            print(f"             = {d:.4f} m = {d * 1000:.1f} mm")

            dn = self.get_standard_dn(d)
            print(f"   公称直径: DN{dn}")

            results["冷却水管"] = {
                "质量流量_kg_h": round(self.cooling_water_flow_rate, 2),
                "体积流量_m3_h": round(water_vol, 2),
                "体积流量_m3_s": round(water_vol_s, 6),
                "设计流速_m_s": v,
                "计算内径_m": round(d, 4),
                "计算内径_mm": round(d * 1000, 1),
                "公称直径_DN": dn,
                "推荐规格": f"DN{dn} 碳钢管，普通自来水管道"
            }
        else:
            print("   警告: 冷却水流量未计算，请先设计冷凝器")
            results["冷却水管"] = {
                "状态": "未计算",
                "原因": "请先设计冷凝器以获取冷却水流量参数"
            }

        return results

    def get_standard_dn(self, inner_diameter_m: float) -> int:
        """
        根据计算内径获取标准公称直径DN

        Args:
            inner_diameter_m: 计算内径 (m)

        Returns:
            标准公称直径DN
        """
        inner_diameter_mm = inner_diameter_m * 1000

        # 标准DN系列 (mm)
        dn_series = [15, 20, 25, 32, 40, 50, 65, 80, 100, 125, 150,
                     200, 250, 300, 350, 400, 450, 500]

        # 找到最接近且不小于计算值的标准DN（考虑安全系数）
        closest_dn = None
        for dn in sorted(dn_series):
            if dn >= inner_diameter_mm * SAFETY_FACTOR_PIPE:  # 考虑安全系数
                closest_dn = dn
                break

        # 如果所有标准值都小于计算值，返回最大值
        if closest_dn is None:
            closest_dn = max(dn_series)

        return closest_dn

    def design_all_equipment(self,
                             water_in_temp: float = 25.0,
                             water_out_temp: float = 40.0,
                             steam_pressure: float = 0.3,
                             pipe_type: str = "直管") -> Dict[str, Any]:
        """
        设计所有附属设备

        Args:
            water_in_temp: 冷却水进口温度 (°C)
            water_out_temp: 冷却水出口温度 (°C)
            steam_pressure: 加热蒸汽压力 (MPa，表压)
            pipe_type: 管道类型

        Returns:
            所有设计结果字典
        """
        print("=" * 80)
        print("板式塔附属设备设计 - 详细计算报告")
        print("=" * 80)
        print(f"设计日期: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"项目: 苯-甲苯精馏塔附属设备设计")
        print(f"\n使用参数范围:")
        print(f"  冷却水进口温度: {COOLING_WATER_IN_TEMP_RANGE[0]}~{COOLING_WATER_IN_TEMP_RANGE[1]}°C")
        print(f"  冷却水出口温度: {COOLING_WATER_OUT_TEMP_RANGE[0]}~{COOLING_WATER_OUT_TEMP_RANGE[1]}°C")
        print(f"  蒸汽压力范围: {STEAM_PRESSURE_RANGE[0]}~{STEAM_PRESSURE_RANGE[1]} MPa")
        print(f"  安全系数 - 面积: {SAFETY_FACTOR_AREA}, 管道: {SAFETY_FACTOR_PIPE}")

        # 设计冷凝器
        print("\n" + "=" * 80)
        print("第一部分：塔顶冷凝器设计")
        print("=" * 80)
        condenser_result = self.design_condenser(water_in_temp, water_out_temp)

        # 设计再沸器
        print("\n" + "=" * 80)
        print("第二部分：塔釜再沸器设计")
        print("=" * 80)
        reboiler_result = self.design_reboiler(steam_pressure)

        # 计算接管尺寸
        print("\n" + "=" * 80)
        print("第三部分：接管尺寸计算")
        print("=" * 80)
        pipe_result = self.calculate_pipe_sizes(pipe_type)

        # 汇总结果
        all_results = {
            "项目信息": {
                "项目名称": "苯-甲苯精馏塔附属设备设计",
                "设计日期": datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                "设计依据": [
                    "benzene_toluene_results.json",
                    "tower_properties_results.json",
                    "tower_diameter_approximate_optimized.json"
                ],
                "塔体直径_m": self.D_tower,
                "塔顶温度_C": self.t_top,
                "塔釜温度_C": self.t_bottom
            },
            "设计参数范围": {
                "冷却水进口温度范围_C": COOLING_WATER_IN_TEMP_RANGE,
                "冷却水出口温度范围_C": COOLING_WATER_OUT_TEMP_RANGE,
                "冷却水温差默认值_C": COOLING_WATER_DELTA_T_DEFAULT,
                "蒸汽压力范围_MPa": STEAM_PRESSURE_RANGE,
                "蒸汽温度最小高于塔釜_C": STEAM_MIN_TEMP_ABOVE_BOTTOM,
                "冷凝器传热系数范围_W_m2K": CONDENSER_HEAT_TRANSFER_COEFFICIENT_RANGE,
                "再沸器传热系数范围_W_m2K": REBOILER_HEAT_TRANSFER_COEFFICIENT_RANGE,
                "管道流速范围_m_s": PIPE_VELOCITY_RANGES,
                "安全系数": {
                    "传热面积": SAFETY_FACTOR_AREA,
                    "管道尺寸": SAFETY_FACTOR_PIPE
                }
            },
            "塔顶冷凝器": condenser_result,
            "塔釜再沸器": reboiler_result,
            "接管尺寸": pipe_result,
            "设计说明": {
                "冷凝器设计": "采用冷却水作为冷却介质，设计泛点率控制在合理范围内",
                "再沸器设计": "采用饱和蒸汽作为加热介质，保证足够的传热温差",
                "管道设计": f"采用{pipe_type}设计，流速控制在推荐范围内",
                "材料选择": "接触物料的管道和设备部件选用不锈钢材质，其他采用碳钢",
                "安全要求": "蒸汽管道需设置安全阀和疏水器，保温层厚度按标准设计"
            }
        }

        # 保存结果到文件
        self.save_results(all_results)

        # 打印设计摘要
        self.print_summary(all_results)

        return all_results

    def save_results(self, results: Dict[str, Any], filename: str = "tower_ancillary_equipment_results.json"):
        """保存设计结果到JSON文件"""
        try:
            with open(filename, 'w', encoding='utf-8') as f:
                json.dump(results, f, ensure_ascii=False, indent=2)
            print(f"\n✓ 详细设计结果已保存到文件: {filename}")
        except Exception as e:
            print(f"✗ 保存结果时出错: {e}")

    def print_summary(self, results: Dict[str, Any]):
        """打印设计摘要"""
        print("\n" + "=" * 80)
        print("设计摘要")
        print("=" * 80)

        condenser = results.get("塔顶冷凝器", {})
        reboiler = results.get("塔釜再沸器", {})
        pipes = results.get("接管尺寸", {})

        print("\n一、塔顶冷凝器:")
        if "设计参数" in condenser:
            params = condenser["设计参数"]
            print(f"   1. 热负荷: {params.get('热负荷_kW', 0):.2f} kW")
            print(
                f"   2. 传热面积: {params.get('设计传热面积_m2', 0):.2f} m² (计算值: {params.get('计算传热面积_m2', 0):.2f} m²)")
            print(
                f"   3. 冷却水用量: {params.get('冷却水用量_kg_h', 0):.0f} kg/h ({params.get('冷却水用量_kg_h', 0) / 1000:.1f} m³/h)")
            print(f"   4. 冷却水温差: {params.get('冷却水出口温度_C', 0) - params.get('冷却水进口温度_C', 0):.1f} °C")

        print("\n二、塔釜再沸器:")
        if "设计参数" in reboiler:
            params = reboiler["设计参数"]
            print(f"   1. 热负荷: {params.get('热负荷_kW', 0):.2f} kW")
            print(
                f"   2. 传热面积: {params.get('设计传热面积_m2', 0):.2f} m² (计算值: {params.get('计算传热面积_m2', 0):.2f} m²)")
            print(f"   3. 加热蒸汽用量: {params.get('加热蒸汽用量_kg_h', 0):.0f} kg/h")
            print(f"   4. 蒸汽压力: {params.get('加热蒸汽压力_MPa', 0):.3f} MPa")
            print(f"   5. 蒸汽温度: {params.get('加热蒸汽温度_C', 0):.1f} °C")

        print("\n三、主要接管尺寸:")
        important_pipes = ["进料管", "回流管", "塔顶蒸汽管", "塔顶出料管", "塔釜出料管", "加热蒸汽管", "冷却水管"]
        for pipe_name in important_pipes:
            if pipe_name in pipes:
                pipe_info = pipes[pipe_name]
                if "公称直径_DN" in pipe_info:
                    print(
                        f"   {pipe_name}: DN{pipe_info.get('公称直径_DN', 0)} ({pipe_info.get('计算内径_mm', 0):.1f} mm)")
                elif "状态" in pipe_info:
                    print(f"   {pipe_name}: {pipe_info.get('状态', '未计算')}")

        print(f"\n注: 设计中使用安全系数 - 面积: {SAFETY_FACTOR_AREA}, 管道: {SAFETY_FACTOR_PIPE}")
        print("\n" + "=" * 80)
        print("设计完成！")
        print("=" * 80)


def interactive_design():
    """交互式设计函数"""
    print("板式塔附属设备设计程序")
    print("=" * 80)
    print("本程序将计算塔顶冷凝器、塔釜再沸器和各接管尺寸")
    print("=" * 80)

    # 显示当前参数范围
    print("\n当前参数范围设置:")
    print(f"  冷却水进口温度: {COOLING_WATER_IN_TEMP_RANGE[0]}~{COOLING_WATER_IN_TEMP_RANGE[1]}°C")
    print(f"  冷却水出口温度: {COOLING_WATER_OUT_TEMP_RANGE[0]}~{COOLING_WATER_OUT_TEMP_RANGE[1]}°C")
    print(f"  蒸汽压力范围: {STEAM_PRESSURE_RANGE[0]}~{STEAM_PRESSURE_RANGE[1]} MPa")
    print(f"  安全系数 - 面积: {SAFETY_FACTOR_AREA}, 管道: {SAFETY_FACTOR_PIPE}")
    print("\n如需修改参数范围，请在代码开头的'用户可修改参数范围'部分进行修改")

    # 加载文件数据
    input_files = {}

    # 尝试加载所有文件
    file_names = [
        'benzene_toluene_results.json',
        'sieve_tower_hydrodynamic_check_complete.json',
        'tower_diameter_approximate_optimized.json',
        'tower_height_design.json',
        'tower_plate_arrangement_corrected.json',
        'tower_properties_results.json'
    ]

    for file_name in file_names:
        try:
            with open(file_name, 'r', encoding='utf-8') as f:
                input_files[file_name] = json.load(f)
            print(f"✓ 已加载文件: {file_name}")
        except FileNotFoundError:
            print(f"⚠ 文件未找到: {file_name}，将使用默认值")

    # 创建设计实例
    try:
        designer = TowerAncillaryEquipmentDesign(input_files)
    except Exception as e:
        print(f"✗ 初始化失败: {e}")
        return

    # 用户输入参数
    print("\n" + "=" * 80)
    print("请输入设计参数")
    print("=" * 80)

    try:
        print("\n一、冷却水参数:")
        print(f"   推荐范围: 进口温度 {COOLING_WATER_IN_TEMP_RANGE[0]}~{COOLING_WATER_IN_TEMP_RANGE[1]}°C")
        water_in_temp = float(input(f"   冷却水进口温度 (°C) [默认25.0]: ") or "25.0")

        print(f"   推荐范围: 出口温度 {COOLING_WATER_OUT_TEMP_RANGE[0]}~{COOLING_WATER_OUT_TEMP_RANGE[1]}°C")
        water_out_temp = float(input(f"   冷却水出口温度 (°C) [默认40.0]: ") or "40.0")

        print("\n二、加热蒸汽参数:")
        print(f"   注意: 蒸汽温度必须高于塔釜温度({designer.t_bottom:.2f}°C)至少{STEAM_MIN_TEMP_ABOVE_BOTTOM}°C")
        print(f"   推荐压力范围: {STEAM_PRESSURE_RANGE[0]}~{STEAM_PRESSURE_RANGE[1]} MPa")
        steam_pressure = float(input("   加热蒸汽压力 (MPa，表压) [默认0.3]: ") or "0.3")

        print("\n三、管道类型:")
        print("   1. 直管 (默认)")
        print("   2. 弯管")
        print("   3. 变径管")
        pipe_choice = input("   选择管道类型 [默认1]: ") or "1"

        pipe_types = {"1": "直管", "2": "弯管", "3": "变径管"}
        pipe_type = pipe_types.get(pipe_choice, "直管")

        # 执行设计
        print("\n" + "=" * 80)
        print("开始设计计算...")
        print("=" * 80)

        results = designer.design_all_equipment(
            water_in_temp=water_in_temp,
            water_out_temp=water_out_temp,
            steam_pressure=steam_pressure,
            pipe_type=pipe_type
        )

        # 询问是否查看详细结果
        view_details = input("\n是否查看详细JSON格式结果？(y/n) [默认n]: ") or "n"
        if view_details.lower() == 'y':
            print("\n详细设计结果(JSON格式):")
            print(json.dumps(results, ensure_ascii=False, indent=2))

        print("\n" + "=" * 80)
        print("设计程序结束")
        print("=" * 80)

    except ValueError as e:
        print(f"✗ 输入错误: {e}")
    except Exception as e:
        print(f"✗ 设计过程中出错: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    # 运行交互式设计
    interactive_design()