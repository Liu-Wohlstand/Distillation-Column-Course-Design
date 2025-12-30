import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


def calculate_average_tower_temperature(txy_dict, xD, xW, pressure=None,
                                        method='arithmetic', plot_diagram=True,
                                        diagram_title="t-x-y Diagram"):
    """
    根据t-x-y数据、塔顶和塔釜组成计算塔平均温度

    参数:
    txy_dict: 字典，格式为 {x: [y, temperature]}
    xD: 塔顶液相组成 (馏出液组成)
    xW: 塔釜液相组成 (釜液组成)
    pressure: 操作压力 (仅用于图表标注)
    method: 平均温度计算方法，可选 'arithmetic' 或 'weighted'
    plot_diagram: 是否绘制t-x-y图
    diagram_title: 图表标题

    返回:
    tuple: (t_top, t_bottom, t_avg, fig, ax) 或 (t_top, t_bottom, t_avg)
    """

    # 1. 数据预处理
    x_list = list(txy_dict.keys())
    t_list = [txy_dict[x][1] for x in x_list]
    y_list = [txy_dict[x][0] for x in x_list]

    # 确保数据有序
    sorted_idx = np.argsort(x_list)
    x_sorted = np.array([x_list[i] for i in sorted_idx])
    t_sorted = np.array([t_list[i] for i in sorted_idx])
    y_sorted = np.array([y_list[i] for i in sorted_idx])

    # 2. 插值函数（用于查找任意组成对应的温度）
    def interpolate_temperature(x_val, x_data, t_data):
        """线性插值获取给定组成对应的温度"""
        if x_val <= min(x_data):
            return t_data[0]
        elif x_val >= max(x_data):
            return t_data[-1]
        else:
            # 找到相邻的数据点
            idx = np.searchsorted(x_data, x_val)
            x1, x2 = x_data[idx - 1], x_data[idx]
            t1, t2 = t_data[idx - 1], t_data[idx]
            # 线性插值
            return t1 + (t2 - t1) * (x_val - x1) / (x2 - x1)

    # 3. 计算塔顶和塔釜温度
    try:
        t_top = interpolate_temperature(xD, x_sorted, t_sorted)
        t_bottom = interpolate_temperature(xW, x_sorted, t_sorted)
    except Exception as e:
        raise ValueError(f"插值计算失败: {e}. 请确保xD({xD})和xW({xW})在数据范围内[0, 1]")

    # 4. 计算平均温度
    if method == 'arithmetic':
        # 算术平均
        t_avg = (t_top + t_bottom) / 2
    elif method == 'weighted':
        # 基于组成差的加权平均
        weight_top = 1 - xW
        weight_bottom = xD
        total_weight = weight_top + weight_bottom
        t_avg = (t_top * weight_top + t_bottom * weight_bottom) / total_weight
    else:
        raise ValueError("method参数必须为 'arithmetic' 或 'weighted'")

    # 5. 绘制t-x-y图（可选）
    fig, ax = None, None
    if plot_diagram:
        fig, ax = plt.subplots(figsize=(10, 8))

        # 绘制泡点线和露点线
        ax.plot(x_sorted, t_sorted, 'b-', linewidth=2, label='Bubble point line', marker='o', markersize=5)
        ax.plot(y_sorted, t_sorted, 'r--', linewidth=2, label='Dew point line', marker='s', markersize=5)

        # 标记塔顶和塔釜组成点
        ax.scatter([xD, xW], [t_top, t_bottom], color='green', s=100, zorder=5,
                   label=f'塔操作点\n塔顶(x={xD:.3f}, T={t_top:.1f}°C)\n塔釜(x={xW:.3f}, T={t_bottom:.1f}°C)')

        # 添加水平线表示温度范围
        ax.axhline(y=t_top, color='g', linestyle=':', alpha=0.5, linewidth=1)
        ax.axhline(y=t_bottom, color='g', linestyle=':', alpha=0.5, linewidth=1)
        ax.axhline(y=t_avg, color='orange', linestyle='-', alpha=0.8, linewidth=2,
                   label=f'平均温度 T_avg={t_avg:.1f}°C')

        # 标记平均温度
        ax.text(0.02, t_avg, f' T_avg = {t_avg:.1f}°C',
                verticalalignment='bottom', color='orange', fontweight='bold')

        # 设置图表属性
        ax.set_xlabel('液相摩尔分数 (x) / 气相摩尔分数 (y)', fontsize=12)
        ax.set_ylabel('温度 (°C)', fontsize=12)

        # 标题中添加压力信息（如果提供）
        if pressure:
            title = f'{diagram_title}\n操作压力: {pressure} kPa'
        else:
            title = diagram_title
        ax.set_title(title, fontsize=14, fontweight='bold')

        ax.grid(True, which='both', linestyle='--', alpha=0.6)
        ax.legend(loc='upper right', fontsize=10)

        # 设置坐标轴范围
        ax.set_xlim(-0.05, 1.05)
        min_temp = min(t_sorted) * 0.98
        max_temp = max(t_sorted) * 1.02
        ax.set_ylim(min_temp, max_temp)

        # 添加更精细的网格
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.xaxis.set_minor_locator(MultipleLocator(0.02))
        ax.yaxis.set_major_locator(MultipleLocator(5))
        ax.yaxis.set_minor_locator(MultipleLocator(1))

        plt.tight_layout()

    # 6. 打印计算结果
    print("=" * 50)
    print("精馏塔温度计算结果")
    print("=" * 50)
    print(f"塔顶组成 xD = {xD:.4f}")
    print(f"塔釜组成 xW = {xW:.4f}")
    print(f"塔顶温度 T_top = {t_top:.2f} °C")
    print(f"塔釜温度 T_bottom = {t_bottom:.2f} °C")
    print(f"塔平均温度 T_avg = {t_avg:.2f} °C (计算方法: {method})")
    print(f"温度范围 ΔT = {abs(t_top - t_bottom):.2f} °C")
    print("=" * 50)

    # 返回结果
    if plot_diagram:
        return t_top, t_bottom, t_avg, fig, ax
    else:
        return t_top, t_bottom, t_avg


# 专门用于估算物性查询温度的函数（简洁版）
def get_property_temperature(txy_dict, xD, xW):
    """
    快速获取用于物性查询的代表性温度

    参数:
    txy_dict: t-x-y数据字典
    xD: 塔顶组成
    xW: 塔釜组成

    返回:
    float: 用于物性查询的代表性温度
    """
    # 直接使用算术平均方法
    t_top, t_bottom, t_avg = calculate_average_tower_temperature(
        txy_dict, xD, xW, plot_diagram=False, method='arithmetic'
    )
    return t_avg


# 示例用法
if __name__ == "__main__":
    # 您的t-x-y数据字典
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
               0.6: [0.791, 88.70],  # 注意：这个y值可能有误，原字典中是80.80，但前面是90.11，后面是87.63
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

    # 示例：计算塔平均温度（课程设计常用场景）
    print("示例1: 完整计算并绘图")
    xD = 0.95  # 塔顶高纯度
    xW = 0.05  # 塔釜低纯度
    t_top, t_bottom, t_avg, fig, ax = calculate_average_tower_temperature(
        dic_txy, xD, xW, pressure=101.3, method='arithmetic',
        diagram_title="苯-甲苯体系 t-x-y 图 (101.3 kPa)"
    )

    # 保存图表
    fig.savefig('txy_diagram_with_tower_points.png', dpi=300, bbox_inches='tight')
    plt.show()

    print("\n" + "=" * 50)
    print("示例2: 快速获取物性查询温度（集成到其他计算中）")
    # 这是您可以在实际塔板数计算程序中调用的方式
    property_temp = get_property_temperature(dic_txy, xD, xW)
    print(f"用于物性查询的代表温度: {property_temp:.2f} °C")