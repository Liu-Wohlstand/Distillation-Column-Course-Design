import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import matplotlib
from matplotlib.font_manager import FontProperties

# ========== 字体设置 ==========
times_new_roman = FontProperties(family='Times New Roman')
matplotlib.rcParams['font.family'] = 'SimSun'
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.use('TkAgg')

# ========== 读取设计数据 ==========
with open('tower_plate_arrangement_corrected.json', 'r', encoding='utf-8') as f:
    plate = json.load(f)
with open('tower_properties_results.json', 'r', encoding='utf-8') as f:
    props = json.load(f)
with open('tower_diameter_approximate_optimized.json', 'r', encoding='utf-8') as f:
    diam = json.load(f)

# ========== 几何参数 ==========
D = diam['最终设计结果']['最终塔径_m']           
At = np.pi * (D/2)**2                            # 塔截面积
Af_single = plate['降液管设计']['降液管面积_Af_m2']   
Af_proj = 2 * Af_single                          # 双侧降液管投影面积
Aa = plate['塔板面积分区']['鼓泡区面积_Aa_m2']      
lw = plate['降液管设计']['堰长_lw_m']            
hw = plate['溢流装置']['堰高_m']                 
h0 = plate['溢流装置']['降液管底隙高度_m']       
A0_total = plate['筛板布置']['总筛孔面积_m2']
d0 = plate['筛板布置']['筛孔直径_m']              
C0 = 0.8          # 孔流系数
phi = 0.5         # 液泛安全系数
g = 9.81

# ========== 物性数据 ==========
rect = props['精馏段结果']
strip = props['提馏段结果']
sections = {
    '精馏段': {
        'rho_L': rect['rho_liquid_kgm3'],
        'rho_V': rect['rho_vapor_kgm3'],
        'sigma': rect['sigma_mNm'],
        'V_vol': rect['V_vol_m3h'] / 3600,
        'L_vol': rect['L_vol_Lh'] / 1000 / 3600,
    },
    '提馏段': {
        'rho_L': strip['rho_liquid_kgm3'],
        'rho_V': strip['rho_vapor_kgm3'],
        'sigma': strip['sigma_mNm'],
        'V_vol': strip['V_vol_m3h'] / 3600,
        'L_vol': strip['L_vol_Lh'] / 1000 / 3600,
    }
}

# ========== 辅助函数 ==========
def calc_h_sigma(sigma, rhoL):
    return 4 * sigma / 1000 / (rhoL * g * d0)

def calc_how(Ls, lw):
    Lh = Ls * 3600
    return 2.84 / 1000 * (Lh / lw) ** (2/3)

def calc_hL(hw, how):
    return hw + how

def calc_hc(Vs, rhoV, rhoL):
    u0 = Vs / A0_total
    term = (u0 / C0)**2 * (rhoV / rhoL) * (1 - (A0_total / Aa)**2)
    return 0.051 * term

def calc_hl(Vs, rhoV, hL):
    ua = Vs / Aa
    Fa = ua * np.sqrt(rhoV)
    eps = 0.971 - 0.355*Fa + 0.0757*Fa**2
    eps = max(0.4, min(0.7, eps))
    return eps * hL

def calc_hd(Ls):
    return 0.153 * (Ls / (lw * h0)) ** 2

def calc_hp(Vs, Ls, rhoV, rhoL, sigma):
    how = calc_how(Ls, lw)
    hL = calc_hL(hw, how)
    hc = calc_hc(Vs, rhoV, rhoL)
    hl = calc_hl(Vs, rhoV, hL)
    hs = calc_h_sigma(sigma, rhoL)
    return hc + hl + hs, hL, how

# ========== 负荷线详细推导输出 ==========
def print_load_line_details(sec_name, data):
    print(f"\n{'='*80}")
    print(f"【{sec_name}】负荷性能图五条线详细推导")
    print('='*80)
    rhoL = data['rho_L']
    rhoV = data['rho_V']
    sigma = data['sigma']
    HT = 0.45

    # 1. 液相下限线 (how = 6 mm)
    print("\n1. 液相下限线 (how = 6 mm)")
    print("   公式：how = 2.84×10⁻³ × (Lh/lw)^{2/3} = 0.006")
    print(f"   已知 lw = {lw:.3f} m")
    Lh_low = ((0.006 / 0.00284) ** (3/2)) * lw
    Ls_low = Lh_low / 3600
    print(f"   解得 Lh = [{0.006/0.00284:.3f}]^{1.5} × lw = {Lh_low:.3f} m³/h")
    print(f"   ∴ Ls = Lh/3600 = {Ls_low:.6f} m³/s")
    print("   该线为垂直线，与纵坐标无关。")

    # 2. 液相上限线 (停留时间 τ = 3 s)
    print("\n2. 液相上限线 (降液管停留时间 τ = 3 s)")
    print("   公式：τ = Af × HT / Ls = 3")
    print(f"   Af = {Af_single:.4f} m², HT = {HT} m")
    Ls_up = Af_single * HT / 3
    print(f"   Ls = Af×HT/3 = {Af_single:.4f}×{HT}/{3} = {Ls_up:.5f} m³/s")
    print("   该线为垂直线。")

    # 3. 漏液线
    print("\n3. 漏液线")
    print("   推导：由漏液点气速公式 u_ow = 4.4 C0 √[ (0.0056+0.13hL - hσ) ρL/ρV ]")
    print("   其中 u_ow = Vs / A0_total，hL = hw + how，how = 2.84×10⁻³ (Lh/lw)^{2/3}，hσ = 4σ/(ρL g d0)")
    print("   整理得 Vs = A0_total × 4.4 C0 √[ (0.0056 + 0.13(hw+how) - hσ) ρL/ρV ]")
    print("   代入已知常数：")
    print(f"   A0_total = {A0_total:.4f} m², C0 = {C0}")
    print(f"   hw = {hw:.4f} m, hσ = {calc_h_sigma(sigma, rhoL):.4f} m")
    print(f"   ρL/ρV = {rhoL/rhoV:.2f}")
    print("   函数关系：Vs = f(Ls) 通过 how(Ls) 表达，绘图用数值点生成。")

    # 4. 液沫夹带线 (ev = 0.1)
    print("\n4. 液沫夹带线 (ev = 0.1 kg液/kg气)")
    print("   亨特公式：ev = (5.7×10⁻⁶ / σ) × [u / (HT - hL)]^{3.2}")
    print("   其中 u = Vs / (At - Af_proj)，hL = hw + how")
    print("   令 ev = 0.1 得：")
    print(f"   0.1 = (5.7e-6 / {sigma/1000:.5f}) × [ Vs / ((At-Af_proj)(HT - hL)) ]^{3.2}")
    print("   即 Vs = (0.1 × σ/5.7e-6)^{1/3.2} × (At-Af_proj) × (HT - hL)")
    print(f"   代入 At = {At:.4f} m², Af_proj = {Af_proj:.4f} m²")
    coeff = (0.1 * (sigma/1000) / 5.7e-6) ** (1/3.2) * (At - Af_proj)
    print(f"   所以 Vs = {coeff:.4f} × (HT - hL)   (无量纲系数)")
    print(f"   其中 hL = hw + how(Ls)，因此 Vs 随 Ls 增大而减小。")

    # 5. 液泛线
    print("\n5. 液泛线 (降液管液泛)")
    print("   液泛条件：Hd = φ (HT + hw)")
    print("   其中 Hd = hp + hL + hd")
    print("   hp = hc + hl + hσ, hc ∝ Vs², hl ∝ ε0 hL, hd ∝ Ls²")
    print("   最终得到关于 Vs 和 Ls 的隐式方程：")
    print("   hc(Vs) + hl(Vs, hL) + hσ + hL + hd(Ls) = φ (HT + hw)")
    print("   本程序采用数值求解，绘图时通过 fsolve 获得 Vs(Ls)。")

# ========== 绘制负荷性能图 ==========
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, (sec_name, data) in zip(axes, sections.items()):
    rhoL = data['rho_L']
    rhoV = data['rho_V']
    sigma = data['sigma']
    HT = 0.45
    Ls_range = np.linspace(0.002, 0.035, 100)

    # 漏液线
    def leak_Vs(Ls):
        how = calc_how(Ls, lw)
        hL = calc_hL(hw, how)
        hs = calc_h_sigma(sigma, rhoL)
        term = (0.0056 + 0.13*hL - hs) * rhoL / rhoV
        if term < 0:
            term = 1e-6
        return A0_total * 4.4 * C0 * np.sqrt(term)
    Vs_leak = [leak_Vs(Ls) for Ls in Ls_range]

    # 液沫夹带线 (ev=0.1)
    def entrain_Vs(Ls):
        how = calc_how(Ls, lw)
        hL = calc_hL(hw, how)
        def residual(Vs):
            u = Vs / (At - Af_proj)
            ev = (5.7e-6 / (sigma/1000)) * (u / (HT - hL)) ** 3.2
            return ev - 0.1
        try:
            Vs_sol = fsolve(residual, 2.0)[0]
            return max(Vs_sol, 0.0)
        except:
            return 0.0
    Vs_entrain = [entrain_Vs(Ls) for Ls in Ls_range]

    # 液泛线
    def flooding_Vs(Ls):
        def func(Vs):
            hp, hL, _ = calc_hp(Vs, Ls, rhoV, rhoL, sigma)
            hd = calc_hd(Ls)
            Hd = hp + hL + hd
            return Hd - 0.5*(HT + hw)
        try:
            Vs_sol = fsolve(func, 2.0)[0]
            return max(Vs_sol, 0.0)
        except:
            return 0.0
    Vs_flood = [flooding_Vs(Ls) for Ls in Ls_range]

    # 液相下限线 (修正单位)
    Lh_low = ((0.006 / 0.00284) ** (3/2)) * lw
    Ls_low = Lh_low / 3600
    # 液相上限线
    Ls_up = Af_single * HT / 3

    # 绘图
    ax.plot(Ls_range, Vs_leak, 'b-', label='漏液线')
    ax.plot(Ls_range, Vs_entrain, 'g-', label='液沫夹带线 (ev=0.1)')
    ax.plot(Ls_range, Vs_flood, 'r-', label='液泛线')
    ax.axvline(Ls_low, color='orange', linestyle='--', label='液相下限 (how=6mm)')
    ax.axvline(Ls_up, color='purple', linestyle='--', label='液相上限 (τ=3s)')

    Ls_op = data['L_vol']
    Vs_op = data['V_vol']
    ax.plot(Ls_op, Vs_op, 'ko', markersize=8, label='操作点')

    ax.set_xlabel('液相流量 $L_s$ (m³/s)')
    ax.set_ylabel('气相流量 $V_s$ (m³/s)')
    ax.set_title(f'{sec_name} 负荷性能图')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right')
    ax.set_xlim(0, 0.035)
    ax.set_ylim(0, 6)

    # 输出该段的负荷线推导
    print_load_line_details(sec_name, data)

plt.tight_layout()
plt.savefig('负荷性能图.png', dpi=150)
plt.show()

# ========== 流体力学验算详细步骤 ==========
print("\n" + "="*80)
print("筛板塔流体力学详细验算步骤（课程设计用）")
print("="*80)

for sec_name, data in sections.items():
    print(f"\n\n【{sec_name}】")
    Ls = data['L_vol']
    Vs = data['V_vol']
    rhoL = data['rho_L']
    rhoV = data['rho_V']
    sigma = data['sigma']
    HT = 0.45

    # 1. 堰上液层高度
    Lh = Ls * 3600
    print("\n1. 堰上液层高度 how")
    print("   公式: how = 2.84×10⁻³ × (Lh/lw)^{2/3}")
    print(f"   已知: Lh = {Lh:.2f} m³/h, lw = {lw:.3f} m")
    how = calc_how(Ls, lw)
    print(f"   计算: how = {how*1000:.1f} mm = {how:.4f} m")

    # 2. 清液层高度
    print("\n2. 清液层高度 hL")
    print("   公式: hL = hw + how")
    print(f"   已知: hw = {hw*1000:.1f} mm, how = {how*1000:.1f} mm")
    hL = calc_hL(hw, how)
    print(f"   计算: hL = {hL*1000:.1f} mm = {hL:.4f} m")

    # 3. 干板压降
    u0 = Vs / A0_total
    print("\n3. 干板压降 hc")
    print("   公式: hc = 0.051 × (u0/C0)² × (ρv/ρL) × [1 - (A0/Aa)²]")
    print(f"   u0 = Vs/A0 = {Vs:.4f}/{A0_total:.4f} = {u0:.2f} m/s")
    print(f"   (u0/C0)² = ({u0:.2f}/{C0})² = { (u0/C0)**2:.2f}")
    print(f"   ρv/ρL = {rhoV:.4f}/{rhoL:.1f} = {rhoV/rhoL:.6f}")
    print(f"   A0/Aa = {A0_total:.4f}/{Aa:.4f} = {A0_total/Aa:.4f}")
    print(f"   1 - (A0/Aa)² = 1 - { (A0_total/Aa)**2:.4f} = {1 - (A0_total/Aa)**2:.4f}")
    hc = calc_hc(Vs, rhoV, rhoL)
    print(f"   计算: hc = {hc*1000:.1f} mm = {hc:.4f} m")

    # 4. 液层压降
    ua = Vs / Aa
    Fa = ua * np.sqrt(rhoV)
    print("\n4. 液层压降 hl")
    print("   公式: hl = ε0 × hL, 其中 ε0 = 0.971 - 0.355×Fa + 0.0757×Fa², Fa = ua×√ρv")
    print(f"   ua = Vs/Aa = {Vs:.4f}/{Aa:.4f} = {ua:.3f} m/s")
    print(f"   Fa = {ua:.3f} × √{rhoV:.4f} = {Fa:.3f}")
    eps = 0.971 - 0.355*Fa + 0.0757*Fa**2
    eps = max(0.4, min(0.7, eps))
    print(f"   ε0 = 0.971 - 0.355×{Fa:.3f} + 0.0757×{Fa:.3f}² = {eps:.4f} (限制在0.4~0.7)")
    hl = calc_hl(Vs, rhoV, hL)
    print(f"   hl = {eps:.4f} × {hL*1000:.1f} mm = {hl*1000:.1f} mm")

    # 5. 表面张力压降
    print("\n5. 表面张力压降 hσ")
    print("   公式: hσ = 4σ / (ρL × g × d0)")
    print(f"   σ = {sigma:.2f} mN/m = {sigma/1000:.4f} N/m")
    hs = calc_h_sigma(sigma, rhoL)
    print(f"   代入: hσ = 4×{sigma/1000:.4f} / ({rhoL:.1f}×9.81×{d0:.4f}) = {hs*1000:.2f} mm")

    # 6. 总压降
    hp = hc + hl + hs
    print("\n6. 塔板总压降 hp")
    print(f"   hp = {hc*1000:.1f} + {hl*1000:.1f} + {hs*1000:.2f} = {hp*1000:.1f} mm液柱")
    print(f"   折合压强 ΔP = hp × ρL × g = {hp:.4f} × {rhoL:.1f} × 9.81 = {hp*rhoL*g:.0f} Pa")

    # 7. 降液管液层高度
    hd = calc_hd(Ls)
    print("\n7. 降液管液层高度 Hd")
    print("   降液管底隙压头损失: hd = 0.153 × (Ls/(lw×h0))²")
    print(f"   Ls/(lw×h0) = {Ls:.4f}/({lw:.3f}×{h0:.3f}) = {Ls/(lw*h0):.3f}")
    print(f"   hd = 0.153 × { (Ls/(lw*h0))**2:.4f} = {hd*1000:.1f} mm")
    Hd = hp + hL + hd
    Hd_max = phi * (HT + hw)
    print(f"   Hd = {hp*1000:.1f} + {hL*1000:.1f} + {hd*1000:.1f} = {Hd*1000:.1f} mm")
    print(f"   允许最大值 Hd,max = Φ×(HT+hw) = {phi}×({HT*1000:.0f}+{hw*1000:.0f}) = {Hd_max*1000:.0f} mm")
    print(f"   判定: {'合格' if Hd<=Hd_max else '不合格'}")

    # 8. 液沫夹带
    u = Vs / (At - Af_proj)
    ev = (5.7e-6 / (sigma/1000)) * (u / (HT - hL)) ** 3.2
    print("\n8. 液沫夹带 ev")
    print("   公式: ev = (5.7×10⁻⁶ / σ) × [u / (HT - hL)]^{3.2}")
    print(f"   u = Vs/(At - Af_proj) = {Vs:.4f} / ({At:.3f} - {Af_proj:.3f}) = {u:.3f} m/s")
    print(f"   HT - hL = {HT:.2f} - {hL:.4f} = {HT - hL:.3f} m")
    print(f"   代入: ev = {5.7e-6:.1e} / {sigma/1000:.4f} × ({u:.3f}/{HT-hL:.3f})^{3.2} = {ev:.6f} kg液/kg气")
    print(f"   限值 0.1 → {'合格' if ev<=0.1 else '不合格'}")

    # 9. 漏液验算
    term = (0.0056 + 0.13*hL - hs) * rhoL / rhoV
    u_ow = 4.4 * C0 * np.sqrt(term)
    K = u0 / u_ow
    print("\n9. 漏液验算")
    print("   漏液点气速公式: u_ow = 4.4×C0 × √[ (0.0056 + 0.13hL - hσ) × (ρL/ρv) ]")
    print(f"   括号内: 0.0056 + 0.13×{hL:.4f} - {hs:.4f} = {0.0056+0.13*hL-hs:.4f}")
    print(f"   ρL/ρv = {rhoL:.1f}/{rhoV:.4f} = {rhoL/rhoV:.2f}")
    print(f"   乘积分 = {term:.4f}")
    print(f"   u_ow = 4.4×{C0} × √{term:.4f} = {u_ow:.2f} m/s")
    print(f"   实际筛孔气速 u0 = {u0:.2f} m/s")
    print(f"   稳定系数 K = u0 / u_ow = {u0:.2f}/{u_ow:.2f} = {K:.2f} (要求 ≥1.5) → {'合格' if K>=1.5 else '不合格'}")

    # 10. 降液管停留时间
    tau = Af_single * HT / Ls
    print("\n10. 降液管停留时间 τ")
    print(f"   τ = Af×HT / Ls = {Af_single:.4f}×{HT:.2f} / {Ls:.4f} = {tau:.2f} s (要求 ≥3s) → {'合格' if tau>=3 else '不合格'}")

print("\n" + "="*80)
print("验算结束")
print("="*80)

# 设计建议
print("\n【设计调整建议】")
print("提馏段降液管液泛不合格 (Hd=303.5mm > 247mm)，建议采取以下措施之一：")
print("  1. 增大降液管面积（增加堰长或塔径）")
print("  2. 增大板间距 HT")
print("  3. 降低堰高 hw（但需保证 how>6mm）")
print("  4. 增大降液管底隙 h0（降低 hd）")
