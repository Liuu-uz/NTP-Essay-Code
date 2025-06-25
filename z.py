import argparse
import mdtraj as md
import os
from multiprocessing import Pool


def parse_filename(filename):
    """
    解析文件名获取蛋白质ID和链ID
    例如: 1d4x_A_af.cif -> protein_id='1d4x', chain_id='A'
    """
    basename = os.path.splitext(filename)[0]  # 去掉.cif
    parts = basename.split('_')
    
    if len(parts) >= 2:
        protein_id = parts[0]
        chain_id = parts[1]
        return protein_id, chain_id
    else:
        return basename, None


def convert_cif_to_pdb(args_tuple):
    """
    将单个CIF文件转换为PDB文件
    """
    cif_file_path, output_folder = args_tuple
    
    try:
        # 获取文件路径信息
        filename = os.path.basename(cif_file_path)
        protein_id, chain_id = parse_filename(filename)
        
        # 加载CIF文件
        structure = md.load(cif_file_path)
        
        # 确保输出文件夹存在
        os.makedirs(output_folder, exist_ok=True)
        
        # 构造输出文件名
        if chain_id:
            output_filename = f"{protein_id}_{chain_id}.pdb"
        else:
            output_filename = f"{protein_id}.pdb"
        
        output_path = os.path.join(output_folder, output_filename)
        
        # 保存为PDB格式
        structure.save_pdb(output_path)
        
        print(f"转换完成: {filename} -> {output_filename}")
        return output_path
        
    except Exception as e:
        print(f"转换文件 {cif_file_path} 时出错: {str(e)}")
        return None


def get_cif_files(folder):
    """
    获取文件夹中所有CIF文件的路径列表
    """
    if not os.path.exists(folder):
        raise ValueError(f"文件夹不存在: {folder}")
    
    all_files = sorted(os.listdir(folder), key=lambda x: x.lower())
    cif_files = [
        os.path.join(folder, f) 
        for f in all_files 
        if f.lower().endswith('.cif')
    ]
    return cif_files


def batch_convert_cif_to_pdb(folder, output_folder="pdb_alpha", processes=4):
    """
    批量转换文件夹中的所有CIF文件为PDB格式
    """
    try:
        cif_files = get_cif_files(folder)
        print(f"找到 {len(cif_files)} 个CIF文件")
        
        if not cif_files:
            print("警告: 在指定文件夹中未找到CIF文件")
            return []
            
    except ValueError as e:
        print(f"错误: {e}")
        return []

    # 创建输出文件夹的完整路径
    if not os.path.isabs(output_folder):
        # 如果是相对路径，基于输入文件夹创建
        full_output_folder = os.path.join(folder, output_folder)
    else:
        full_output_folder = output_folder
    
    print(f"PDB文件将保存到: {full_output_folder}")
    
    # 确保输出文件夹存在
    os.makedirs(full_output_folder, exist_ok=True)

    # 准备参数元组
    args_list = [(cif_file, full_output_folder) for cif_file in cif_files]

    # 并行转换
    print(f"开始转换，使用 {processes} 个进程...")
    
    converted_files = []
    
    with Pool(processes=processes) as pool:
        results = pool.imap(convert_cif_to_pdb, args_list)
        
        for i, result in enumerate(results, 1):
            print(f"完成 {i}/{len(cif_files)}")
            if result:
                converted_files.append(result)
    
    print(f"转换完成！成功转换 {len(converted_files)} 个文件")
    print(f"所有PDB文件已保存到: {full_output_folder}")
    return converted_files


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='将CIF文件批量转换为PDB格式')
    parser.add_argument('--folder', required=True, help='包含CIF文件的目录')
    parser.add_argument('--output', default='pdb_alpha', help='输出文件夹名称 (默认: pdb_alpha)')
    parser.add_argument('--processes', type=int, default=4, help='并行进程数 (默认: 4)')
    parser.add_argument('--single', help='转换单个文件（提供完整路径）')
    
    args = parser.parse_args()

    if args.single:
        # 转换单个文件
        input_dir = os.path.dirname(args.single)
        output_folder = os.path.join(input_dir, args.output)
        os.makedirs(output_folder, exist_ok=True)
        
        result = convert_cif_to_pdb((args.single, output_folder))
        if result:
            print(f"文件已转换为: {result}")
    else:
        # 批量转换
        converted_files = batch_convert_cif_to_pdb(args.folder, args.output, args.processes)
        
        if converted_files:
            print("\n转换后的PDB文件:")
            for pdb_file in converted_files:
                print(f"  {os.path.basename(pdb_file)}")