import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.nsgrid import FastNS

def count_all_residues_near_li(li_index, start_frame, end_frame, threshold=2.5, output_file='dmc_count_all_residues.txt'):
    topology_file = 'eq2.gro'
    trajectory_file = 'eq2_fixed_10-20ns.xtc'

    u = mda.Universe(topology_file, trajectory_file)
    li_atom = u.atoms[li_index]

    all_resnames = set(u.atoms.residues.resnames)
    if 'Li' in all_resnames:
        all_resnames.remove('Li')
    sorted_resnames = sorted(all_resnames)

    with open(output_file, 'w') as f:
        header = "Frame\t" + "\t".join(sorted_resnames)
        print(header)          # 控制台打印表头
        f.write(header + "\n") # 文件写入表头

        for ts in u.trajectory[start_frame:end_frame + 1]:
            li_pos = li_atom.position.reshape(1, 3)
            counts = {}
            for resname in sorted_resnames:
                target_atoms = u.select_atoms(f'resname {resname}')
                if len(target_atoms) == 0:
                    counts[resname] = 0
                    continue

                ns = FastNS(threshold, target_atoms.positions, box=u.dimensions)
                results = ns.search(li_pos)
                pairs = results.get_pairs()
                if pairs.size == 0:
                    neighbor_indices = np.array([], dtype=int)
                else:
                    neighbor_indices = pairs[:, 1]
                nearby_resids = set(target_atoms[neighbor_indices].resids)
                counts[resname] = len(nearby_resids)

            count_str = "\t".join(str(counts[res]) for res in sorted_resnames)
            print(f"{ts.frame}\t{count_str}")
            f.write(f"{ts.frame}\t{count_str}\n")

if __name__ == "__main__":
    print("=== 分子计数程序（自动所有残基版） ===")
    li_index = int(input("请输入 Li 原子索引（0开始）: ").strip())
    start_frame = int(input("请输入起始帧编号: ").strip())
    end_frame = int(input("请输入终止帧编号: ").strip())
    threshold_input = input("请输入距离阈值（单位 Å，默认2.5）: ").strip()
    threshold = float(threshold_input) if threshold_input else 2.5

    count_all_residues_near_li(
        li_index=li_index,
        start_frame=start_frame,
        end_frame=end_frame,
        threshold=threshold
    )
