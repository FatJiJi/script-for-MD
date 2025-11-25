#计算平均位移大小和差分法估算的速度（有限差分近似速度）,在去除pbc确保分子完整后使用
import MDAnalysis as mda
import numpy as np

def calculate_avg_displacement(gro_file, xtc_file, atom_name, frame1, frame2):
    # 载入系统，包含结构文件（gro）和轨迹文件（xtc）
    u = mda.Universe(gro_file, xtc_file)
    
    # 按原子名选择原子，使用通配符*匹配可能的后缀
    atoms = u.select_atoms(f"name {atom_name}*")
    if len(atoms) == 0:
        raise ValueError(f"No atoms with name starting with '{atom_name}' found.")
    
    # 定位到第frame1帧，读取该帧中选中原子的坐标
    u.trajectory[frame1]
    pos1 = atoms.positions.copy()   # 获取该帧的原子坐标（单位：Å）
    time1 = u.trajectory.time       # 获取该帧对应的模拟时间（单位：ps）
    
    # 定位到第frame2帧，读取该帧中选中原子的坐标
    u.trajectory[frame2]
    pos2 = atoms.positions.copy()
    time2 = u.trajectory.time
    
    # 计算两帧之间的时间间隔（单位：ps）
    delta_t = time2 - time1
    if delta_t <= 0:
        raise ValueError(f"Second frame time ({time2}) must be greater than first frame time ({time1}).")
    
    # 打印时间信息，方便确认
    print(f"Frame {frame1} time: {time1:.3f} ps")
    print(f"Frame {frame2} time: {time2:.3f} ps")
    print(f"Time interval: {delta_t:.3f} ps\n")
    
    # 计算两帧间的位移向量，单位转换：Å → nm （1 Å = 0.1 nm）
    displacements = (pos2 - pos1) * 0.1
    
    # 计算每个原子位移的大小（距离），单位：nm
    displacement_magnitudes = np.linalg.norm(displacements, axis=1)
    
    # 计算每个原子平均速度 = 位移 / 时间间隔，单位：nm/ps
    velocities = displacement_magnitudes / delta_t
    
    # 计算所有选中原子的平均位移和平均速度
    avg_displacement = np.mean(displacement_magnitudes)
    avg_speed = np.mean(velocities)
    
    # 打印每个原子的位移和对应速度
    for i, (disp, vel) in enumerate(zip(displacement_magnitudes, velocities)):
        print(f"Atom {i+1}: displacement = {disp:.5f} nm, speed = {vel:.5f} nm/ps")
    
    # 打印所有原子的平均位移和平均速度
    print(f"\nAverage displacement magnitude: {avg_displacement:.5f} nm")
    print(f"Average speed: {avg_speed:.5f} nm/ps")

if __name__ == "__main__":
    # 交互式输入，用户输入结构文件和轨迹文件路径
    gro_file = input("请输入结构文件路径（如 topology.gro）：").strip()
    xtc_file = input("请输入轨迹文件路径（如 trajectory.xtc）：").strip()
    
    # 输入想统计的原子名，及计算的两帧索引
    atom_name = input("请输入想统计位移的原子名（如 Li）：").strip()
    frame1 = int(input("请输入第一个帧的索引（0开始）："))
    frame2 = int(input("请输入第二个帧的索引（0开始）："))
    
    # 调用函数执行计算和打印
    calculate_avg_displacement(gro_file, xtc_file, atom_name, frame1, frame2)
