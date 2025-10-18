#计算某原子位移的平均值
import MDAnalysis as mda
import numpy as np

def calculate_avg_displacement(gro_file, xtc_file, atom_name, frame1, frame2):
    # 载入系统，包含拓扑gro和轨迹xtc
    u = mda.Universe(gro_file, xtc_file)
    
    # 按原子名选择原子（例如 Li）
    atoms = u.select_atoms(f"name {atom_name}*")#注意名字有的名字比如是Li0128081，所以要加*
    if len(atoms) == 0:
        raise ValueError(f"No atoms with name '{atom_name}' found.")
    
    # 定位到第frame1帧，读取该帧中选中原子的坐标
    u.trajectory[frame1]
    pos1 = atoms.positions.copy()
    
    # 定位到第frame2帧，读取该帧中选中原子的坐标
    u.trajectory[frame2]
    pos2 = atoms.positions.copy()
    
    # 计算两帧间的位移向量，单位由Å转为nm（1Å = 0.1 nm）
    displacements = (pos2 - pos1) * 0.1
    
    # 计算每个原子位移向量的模长（即每个原子移动的距离，单位nm）
    displacement_magnitudes = np.linalg.norm(displacements, axis=1)
    
    # 计算所有选中原子位移距离的平均值（标量平均）
    avg_displacement_magnitude = np.mean(displacement_magnitudes)
    
    # 打印每个原子的移动距离
    for i, disp in enumerate(displacement_magnitudes):
        print(f"Atom {i+1}: displacement = {disp:.5f} nm")
    
    # 打印平均位移大小
    print(f"\nAverage displacement magnitude (nm): {avg_displacement_magnitude:.5f}")

if __name__ == "__main__":
    # 交互式输入，用户输入拓扑和轨迹文件路径
    gro_file = input("请输入结构文件路径（如 topology.gro）：").strip()
    xtc_file = input("请输入轨迹文件路径（如 trajectory.xtc）：").strip()
    
    # 交互式输入，用户输入需要计算的原子名和两帧索引
    atom_name = input("请输入想统计位移的原子名（如 Li）：").strip()
    frame1 = int(input("请输入第一个帧的索引（0开始）："))
    frame2 = int(input("请输入第二个帧的索引（0开始）："))
    
    # 调用函数执行计算和打印
    calculate_avg_displacement(gro_file, xtc_file, atom_name, frame1, frame2)
