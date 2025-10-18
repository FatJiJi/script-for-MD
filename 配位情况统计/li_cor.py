#统计某一帧所有Li原子配位情况的脚本
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.distances import distance_array

def analyze_li_coordination(u, li_resname, frame_number, threshold=2.5):
    # 定位所有Li原子（通过残基名）
    li_atoms = u.select_atoms(f'resname {li_resname}')
    if len(li_atoms) == 0:
        print(f"警告：未找到名为 '{li_resname}' 的残基，请确认名称是否正确。")
        return

    print(f"轨迹中Li原子数量: {len(li_atoms)}")
    # 跳转到指定帧
    u.trajectory[frame_number]

    # 取当前帧所有非Li原子（用于统计配位分子）
    non_li_atoms = u.select_atoms(f'not resname {li_resname}')
    non_li_resids = np.unique(non_li_atoms.resids)

    # 建立分子残基字典，key=残基号，value=AtomGroup
    res_groups = {}
    for resid in non_li_resids:
        res_groups[resid] = non_li_atoms.select_atoms(f'resid {resid}')

    threshold = float(threshold)

    print(f"\n第 {frame_number} 帧中各 Li 原子的配位情况（距离阈值：{threshold} Å）：\n")

    # 记录所有Li配位模式
    coordination_patterns = {}

    for li_atom in li_atoms:
        li_pos = li_atom.position.reshape(1, 3)
        # 统计该Li周围分子类别和数量
        mol_count = {}
        for resid, group in res_groups.items():
            dists = distance_array(li_pos, group.positions, box=u.dimensions)
            if np.any(dists <= threshold):
                resname = group.resnames[0]
                mol_count[resname] = mol_count.get(resname, 0) + 1

        if not mol_count:
            print(f"Li原子 {li_atom.index} 周围 {threshold} Å 内无配位分子")
        else:
            total_mol = sum(mol_count.values())
            print(f"Li原子 {li_atom.index} 周围 {threshold} Å 内有 {total_mol} 个分子：")
            for k, v in mol_count.items():
                print(f"    {k}: {v} 个")

            # 生成配位模式元组 (按字母序排列)
            pattern = tuple(sorted(mol_count.items()))
            coordination_patterns[pattern] = coordination_patterns.get(pattern, 0) + 1

    print("\n最终统计：")
    for pattern, count in sorted(coordination_patterns.items(), key=lambda x: x[1], reverse=True):
        pattern_str = ', '.join([f"{name}×{num}" for name, num in pattern])
        print(f"    配位模式 ({pattern_str}) ： {count} 个Li")

if __name__ == "__main__":
    print("=== 所有 Li 的配位分析 ===")
    li_resname = input("请输入Li的残基名（例如 Li）：").strip()
    frame_number = int(input("请输入要分析的帧编号: ").strip())
    threshold = input("请输入距离阈值（单位 Å，默认 2.5）: ").strip()
    if threshold == '':
        threshold = 2.5
    else:
        threshold = float(threshold)

    # 你的文件名
    topology_file = 'eq2.gro'
    trajectory_file = 'eq2_fixed_10-20ns.xtc'

    u = mda.Universe(topology_file, trajectory_file)
    analyze_li_coordination(u, li_resname, frame_number, threshold)
