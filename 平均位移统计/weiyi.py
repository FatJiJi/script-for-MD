from MDAnalysis.coordinates.PDB import PDBReader
import numpy as np

def read_positions(pdb_file):
    with PDBReader(pdb_file) as reader:
        coords = reader.trajectory[0].positions
    return coords  # 单位为 Å（angstrom）

# 修改对应文件路径 
file1 = "5000_Li.pdb"  # t1 帧
file2 = "5100_Li.pdb"  # t2 帧
delta_t = 100  # ps，两帧之间的时间间隔

# 读取两帧的坐标
pos1 = read_positions(file1)
pos2 = read_positions(file2)

# 计算每个原子的位移向量，单位转换 Å -> nm
displacements = (pos2 - pos1) * 0.1  # 位移向量，单位 nm

# 计算每个原子的位移大小（标量）
displacement_magnitudes = np.linalg.norm(displacements, axis=1)  # 位移大小，单位 nm

# 计算所有原子的平均位移向量（矢量平均）
avg_displacement = np.mean(displacements, axis=0)  # 平均位移向量，单位 nm

# 计算平均位移大小（标量）
avg_displacement_magnitude = np.linalg.norm(avg_displacement)  # 平均位移大小，单位 nm

# 打印每个原子的位移大小
for i, disp in enumerate(displacement_magnitudes):
    print(f"Atom {i+1}: displacement = {disp:.5f} nm")

print(f"\nAverage displacement vector (nm): {avg_displacement}")
print(f"Average displacement magnitude (nm): {avg_displacement_magnitude:.5f}")
