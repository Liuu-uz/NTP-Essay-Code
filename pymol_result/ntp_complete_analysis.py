#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
基于UniProt ID的NTP结合位点分析脚本
从Excel读取UniProt ID，查找对应FASTA序列和PDB结构
"""

import os
import urllib.request
import time
import json
import requests
from pathlib import Path

def read_excel_data(excel_file='NTP_Analysis_Report.xlsx'):
    """
    读取Excel文件，提取NTP processing SF found的数据
    返回第1列（文件夹名）、第3列（UniProt ID）、第7列的信息
    """
    print("📊 读取Excel数据...")
    
    try:
        import pandas as pd
        
        if not os.path.exists(excel_file):
            print(f"❌ Excel文件不存在: {excel_file}")
            return None
            
        # 读取Excel文件
        df = pd.read_excel(excel_file)
        print(f"✅ 读取到 {len(df)} 行数据")
        print(f"📋 列名: {list(df.columns)}")
        
        # 筛选最后一列为 "NTP processing SF found" 的行
        last_col = df.columns[-1]  # 最后一列
        ntp_rows = df[df[last_col] == 'NTP processing SF found']
        print(f"✅ 筛选到 {len(ntp_rows)} 行 NTP processing SF found 数据")
        
        if len(ntp_rows) == 0:
            print("❌ 没有找到 NTP processing SF found 的数据")
            return None
        
        # 提取第1、3、7列数据
        col_1 = df.columns[0]  # 第1列 (文件夹名)
        col_3 = df.columns[2]  # 第3列 (UniProt ID)
        col_7 = df.columns[6]  # 第7列
        
        print(f"📋 提取列:")
        print(f"  第1列='{col_1}' (文件夹名)")
        print(f"  第3列='{col_3}' (UniProt ID)")
        print(f"  第7列='{col_7}'")
        
        uniprot_data = []
        for idx, row in ntp_rows.iterrows():
            folder_name = str(row[col_1]).strip()
            uniprot_id = str(row[col_3]).strip()
            col7_value = str(row[col_7]).strip()
            
            # 跳过空值
            if any(val in ['nan', 'NaN', '', 'N/A'] for val in [folder_name, uniprot_id]):
                continue
                
            uniprot_data.append({
                'folder_name': folder_name,
                'uniprot_id': uniprot_id,
                'col7_info': col7_value,
                'row_index': idx
            })
        
        print(f"✅ 有效数据: {len(uniprot_data)} 条")
        for item in uniprot_data[:5]:  # 显示前5条
            print(f"  • {item['folder_name']}/{item['uniprot_id']} -> {item['col7_info']}")
        
        return uniprot_data
        
    except Exception as e:
        print(f"❌ Excel读取失败: {e}")
        return None

def find_fasta_sequences(uniprot_data, fasta_base_dir='/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/fasta_files'):
    """
    根据文件夹名和UniProt ID查找对应的FASTA序列
    路径结构: fasta_base_dir/第一列文件夹名/第三列UniProt_ID.fasta
    """
    print(f"\n🔍 在 {fasta_base_dir} 中查找FASTA序列...")
    
    if not os.path.exists(fasta_base_dir):
        print(f"❌ FASTA基础目录不存在: {fasta_base_dir}")
        return []
    
    # 先列出基础目录下的所有文件夹
    try:
        available_folders = [f for f in os.listdir(fasta_base_dir) 
                           if os.path.isdir(os.path.join(fasta_base_dir, f))]
        print(f"📁 可用文件夹: {len(available_folders)} 个")
        for folder in available_folders[:5]:  # 显示前5个
            print(f"  • {folder}")
        if len(available_folders) > 5:
            print(f"  ... 还有 {len(available_folders)-5} 个")
    except Exception as e:
        print(f"❌ 无法读取目录: {e}")
        return []
    
    found_sequences = []
    
    for item in uniprot_data:
        folder_name = item['folder_name']
        uniprot_id = item['uniprot_id']
        
        print(f"\n🔍 查找 [{folder_name}] 中的 [{uniprot_id}]...")
        
        # 目标文件夹路径
        target_folder = os.path.join(fasta_base_dir, folder_name)
        
        if not os.path.exists(target_folder):
            print(f"  ❌ 文件夹不存在: {folder_name}")
            # 尝试查找相似的文件夹名
            similar_folders = [f for f in available_folders if folder_name.lower() in f.lower() or f.lower() in folder_name.lower()]
            if similar_folders:
                print(f"  💡 相似文件夹: {', '.join(similar_folders[:3])}")
            continue
        
        # 列出文件夹中的文件
        try:
            files_in_folder = os.listdir(target_folder)
            print(f"  📁 文件夹 [{folder_name}] 包含 {len(files_in_folder)} 个文件")
        except Exception as e:
            print(f"  ❌ 无法读取文件夹: {e}")
            continue
        
        # 构建可能的FASTA文件路径
        possible_files = [
            f"{uniprot_id}.fasta",
            f"{uniprot_id}.fa", 
            f"{uniprot_id}.seq",
            f"{uniprot_id}.txt"
        ]
        
        # 查找精确匹配的文件
        fasta_file = None
        matched_filename = None
        
        for filename in possible_files:
            file_path = os.path.join(target_folder, filename)
            if os.path.exists(file_path):
                fasta_file = file_path
                matched_filename = filename
                break
        
        # 如果精确匹配失败，尝试模糊匹配
        if not fasta_file:
            print(f"  🔍 精确匹配失败，尝试模糊匹配...")
            for file in files_in_folder:
                # 检查文件名是否包含UniProt ID
                if (uniprot_id.lower() in file.lower() and 
                    file.lower().endswith(('.fasta', '.fa', '.seq', '.txt'))):
                    fasta_file = os.path.join(target_folder, file)
                    matched_filename = file
                    print(f"  💡 模糊匹配到: {file}")
                    break
        
        if fasta_file:
            try:
                # 读取FASTA序列
                with open(fasta_file, 'r') as f:
                    content = f.read().strip()
                
                if not content:
                    print(f"  ❌ 文件为空: {matched_filename}")
                    continue
                
                # 解析FASTA格式
                lines = content.split('\n')
                
                # 处理多种FASTA格式
                if lines[0].startswith('>'):
                    # 标准FASTA格式
                    header = lines[0][1:]  # 去掉>
                    sequence = ''.join(line.strip() for line in lines[1:] if not line.startswith('>'))
                else:
                    # 纯序列格式（无header）
                    header = f"{uniprot_id} from {matched_filename}"
                    sequence = ''.join(line.strip() for line in lines if line.strip())
                
                # 验证序列（应该主要包含氨基酸字母）
                if sequence and len(sequence) > 20:  # 至少20个氨基酸
                    # 移除非氨基酸字符
                    clean_sequence = ''.join(c for c in sequence.upper() if c in 'ACDEFGHIKLMNPQRSTVWY')
                    
                    if len(clean_sequence) >= len(sequence) * 0.8:  # 至少80%是有效氨基酸
                        item['fasta_file'] = fasta_file
                        item['fasta_header'] = header
                        item['sequence'] = clean_sequence
                        item['sequence_length'] = len(clean_sequence)
                        item['matched_filename'] = matched_filename
                        
                        found_sequences.append(item)
                        print(f"  ✅ 成功: {matched_filename} ({len(clean_sequence)} aa)")
                    else:
                        print(f"  ❌ 序列质量差: {matched_filename} (有效氨基酸比例: {len(clean_sequence)/len(sequence)*100:.1f}%)")
                else:
                    print(f"  ❌ 序列太短: {matched_filename} ({len(sequence)} 字符)")
                    
            except Exception as e:
                print(f"  ❌ 读取失败: {matched_filename} - {e}")
        else:
            print(f"  ❌ 未找到FASTA文件")
            # 显示文件夹中的前几个文件供参考
            fasta_files = [f for f in files_in_folder if f.lower().endswith(('.fasta', '.fa', '.seq', '.txt'))]
            if fasta_files:
                print(f"  💡 文件夹中的FASTA文件: {', '.join(fasta_files[:3])}")
    
    print(f"\n✅ 成功找到 {len(found_sequences)}/{len(uniprot_data)} 个FASTA序列")
    return found_sequences

def query_uniprot_for_pdb(uniprot_id):
    """
    通过UniProt API查询对应的PDB结构
    """
    try:
        # UniProt API查询PDB信息
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        response = requests.get(url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            
            # 提取PDB信息
            pdb_entries = []
            if 'uniProtKBCrossReferences' in data:
                for ref in data['uniProtKBCrossReferences']:
                    if ref.get('database') == 'PDB':
                        pdb_id = ref.get('id', '').lower()
                        if pdb_id:
                            pdb_entries.append(pdb_id)
            
            return pdb_entries
        else:
            return []
            
    except Exception as e:
        print(f"    ⚠️  UniProt查询失败: {e}")
        return []

def find_pdb_structures_for_sequences(sequence_data):
    """
    为每个UniProt ID查找对应的PDB结构
    """
    print(f"\n🔍 查找PDB结构...")
    
    sequences_with_pdb = []
    
    for item in sequence_data:
        uniprot_id = item['uniprot_id']
        print(f"🔍 查询 {uniprot_id} 的PDB结构...")
        
        # 查询UniProt数据库
        pdb_entries = query_uniprot_for_pdb(uniprot_id)
        
        if pdb_entries:
            item['pdb_entries'] = pdb_entries
            sequences_with_pdb.append(item)
            print(f"  ✅ 找到PDB: {', '.join(pdb_entries)}")
        else:
            print(f"  ❌ 未找到PDB结构")
            item['pdb_entries'] = []
        
        time.sleep(0.5)  # 避免请求过快
    
    print(f"✅ {len(sequences_with_pdb)}/{len(sequence_data)} 个序列有PDB结构")
    return sequences_with_pdb

def download_pdb_files(sequence_data, output_dir='ntp_analysis_output'):
    """
    下载找到的PDB文件
    """
    print(f"\n📥 下载PDB文件...")
    
    pdb_dir = os.path.join(output_dir, 'pdb_files')
    os.makedirs(pdb_dir, exist_ok=True)
    
    successful_downloads = []
    
    for item in sequence_data:
        if not item.get('pdb_entries'):
            continue
            
        uniprot_id = item['uniprot_id']
        pdb_entries = item['pdb_entries']
        
        print(f"📥 下载 {uniprot_id} 的PDB文件...")
        
        downloaded_pdbs = []
        
        for pdb_id in pdb_entries[:3]:  # 限制最多3个PDB文件
            pdb_file = os.path.join(pdb_dir, f"{pdb_id}.pdb")
            
            # 检查文件是否已存在
            if os.path.exists(pdb_file) and os.path.getsize(pdb_file) > 1000:
                print(f"  ✅ {pdb_id} 已存在")
                downloaded_pdbs.append(pdb_id)
                continue
            
            # 下载PDB文件
            try:
                url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
                urllib.request.urlretrieve(url, pdb_file)
                time.sleep(1)
                
                if os.path.exists(pdb_file) and os.path.getsize(pdb_file) > 1000:
                    print(f"  ✅ {pdb_id} 下载成功")
                    downloaded_pdbs.append(pdb_id)
                else:
                    print(f"  ❌ {pdb_id} 下载失败")
                    
            except Exception as e:
                print(f"  ❌ {pdb_id} 下载错误: {e}")
        
        if downloaded_pdbs:
            item['downloaded_pdbs'] = downloaded_pdbs
            successful_downloads.append(item)
    
    print(f"✅ 成功下载: {len(successful_downloads)} 个UniProt条目的PDB文件")
    return successful_downloads

def analyze_pdb_for_atp(pdb_file):
    """
    分析PDB文件中的ATP分子
    """
    atp_types = ['ATP', 'ADP', 'AMP', 'ANP', 'ACP', 'AGS', 'GTP', 'GDP', 'GMP']
    found_ligands = set()
    
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('HETATM'):
                    res_name = line[17:20].strip()
                    if res_name in atp_types:
                        found_ligands.add(res_name)
        return list(found_ligands)
    except:
        return []

def generate_pymol_script(sequence_data, output_dir='ntp_analysis_output'):
    """
    生成PyMOL脚本文件
    """
    print(f"\n📝 生成PyMOL脚本...")
    
    script_file = os.path.join(output_dir, 'pymol_uniprot_analysis.py')
    
    with open(script_file, 'w') as f:
        f.write("#!/usr/bin/env python\n")
        f.write("# -*- coding: utf-8 -*-\n\n")
        f.write("# PyMOL脚本 - 基于UniProt ID的NTP结合位点分析\n")
        f.write("# 在PyMOL中运行: run pymol_uniprot_analysis.py\n\n")
        
        f.write("import pymol\nfrom pymol import cmd\nimport os\n\n")
        f.write("# 启动PyMOL\n")
        f.write("try:\n")
        f.write("    pymol.finish_launching(['pymol', '-qc'])\n")
        f.write("except:\n")
        f.write("    pymol.finish_launching()\n\n")
        
        f.write("# 清理环境\n")
        f.write("cmd.delete('all')\n")
        f.write("cmd.reinitialize()\n\n")
        
        # 定义颜色
        colors = ['blue', 'red', 'green', 'orange', 'yellow', 'purple', 'cyan', 'magenta']
        
        f.write("# 加载结构\n")
        f.write("loaded_structures = []\n")
        f.write("pdb_dir = 'pdb_files'\n\n")
        
        struct_count = 0
        for item in sequence_data:
            if not item.get('downloaded_pdbs'):
                continue
                
            uniprot_id = item['uniprot_id']
            folder_name = item['folder_name']
            
            for pdb_id in item['downloaded_pdbs']:
                color = colors[struct_count % len(colors)]
                struct_name = f"protein_{uniprot_id}_{pdb_id}"
                
                f.write(f"\n# 加载 {uniprot_id} -> {pdb_id}\n")
                f.write(f"pdb_file = os.path.join(pdb_dir, '{pdb_id}.pdb')\n")
                f.write(f"if os.path.exists(pdb_file):\n")
                f.write(f"    cmd.load(pdb_file, '{pdb_id}')\n")
                f.write(f"    cmd.set_name('{pdb_id}', '{struct_name}')\n")
                f.write(f"    cmd.color('{color}', '{struct_name}')\n")
                f.write(f"    loaded_structures.append(('{uniprot_id}', '{pdb_id}', '{struct_name}'))\n")
                f.write(f"    print(f'✅ 加载 {uniprot_id} -> {pdb_id}')\n")
                f.write(f"    \n")
                
                # 分析ATP
                pdb_file_path = os.path.join(output_dir, 'pdb_files', f'{pdb_id}.pdb')
                if os.path.exists(pdb_file_path):
                    found_ligands = analyze_pdb_for_atp(pdb_file_path)
                    if found_ligands:
                        ligand_str = "' or resn '".join(found_ligands)
                        f.write(f"    # ATP分析 for {pdb_id}\n")
                        f.write(f"    cmd.select('atp_{struct_name}', '{struct_name} and (resn {ligand_str})')\n")
                        f.write(f"    if cmd.count_atoms('atp_{struct_name}') > 0:\n")
                        f.write(f"        cmd.select('active_site_{struct_name}', '{struct_name} and (atp_{struct_name} around 5)')\n")
                        f.write(f"        print(f'  ✅ 发现ATP: {'+'.join(found_ligands)}')\n")
                        f.write(f"    \n")
                
                struct_count += 1
                
                if struct_count >= 8:  # 限制结构数量
                    break
            
            if struct_count >= 8:
                break
        
        f.write("print(f'总共加载: {len(loaded_structures)} 个结构')\n\n")
        
        f.write("# 结构对齐\n")
        f.write("if len(loaded_structures) > 1:\n")
        f.write("    reference = loaded_structures[0][2]  # 第一个结构作为参考\n")
        f.write("    print(f'参考结构: {loaded_structures[0][0]} -> {loaded_structures[0][1]}')\n")
        f.write("    \n")
        f.write("    for uniprot_id, pdb_id, struct_name in loaded_structures[1:]:\n")
        f.write("        try:\n")
        f.write("            result = cmd.align(struct_name, reference)\n")
        f.write("            rmsd = result[0]\n")
        f.write("            print(f'{uniprot_id} -> {pdb_id} 对齐 RMSD: {rmsd:.3f} Å')\n")
        f.write("        except Exception as e:\n")
        f.write("            print(f'{uniprot_id} -> {pdb_id} 对齐失败: {e}')\n\n")
        
        f.write("# 设置显示样式\n")
        f.write("for uniprot_id, pdb_id, struct_name in loaded_structures:\n")
        f.write("    cmd.hide('everything', struct_name)\n")
        f.write("    cmd.show('cartoon', struct_name)\n")
        f.write("    cmd.set('cartoon_transparency', 0.3, struct_name)\n")
        f.write("    \n")
        f.write("    # ATP显示\n")
        f.write("    if cmd.count_atoms(f'atp_{struct_name}') > 0:\n")
        f.write("        cmd.show('sticks', f'atp_{struct_name}')\n")
        f.write("        cmd.color('forest', f'atp_{struct_name} and elem C')\n")
        f.write("        cmd.color('red', f'atp_{struct_name} and elem O')\n")
        f.write("        cmd.color('blue', f'atp_{struct_name} and elem N')\n")
        f.write("        cmd.color('orange', f'atp_{struct_name} and elem P')\n")
        f.write("        \n")
        f.write("        cmd.show('sticks', f'active_site_{struct_name}')\n")
        f.write("        cmd.color('yellow', f'active_site_{struct_name} and elem C')\n")
        f.write("        \n")
        f.write("        cmd.distance(f'hbonds_{struct_name}', f'atp_{struct_name}', f'active_site_{struct_name}', cutoff=3.5, mode=2)\n")
        f.write("        cmd.hide('labels', f'hbonds_{struct_name}')\n")
        f.write("        cmd.color('cyan', f'hbonds_{struct_name}')\n\n")
        
        f.write("# 设置视角和渲染\n")
        f.write("cmd.set('cartoon_fancy_helices', 1)\n")
        f.write("cmd.set('stick_radius', 0.2)\n")
        f.write("cmd.set('ray_shadows', 1)\n")
        f.write("cmd.center('all')\n")
        f.write("cmd.zoom('all', buffer=10)\n\n")
        
        f.write("# 保存图像\n")
        f.write("cmd.png('images/UniProt_NTP_structures_overview.png', width=1800, height=1400, dpi=300, ray=1)\n")
        f.write("print('✅ 保存全景图: UniProt_NTP_structures_overview.png')\n\n")
        
        f.write("# 保存会话\n")
        f.write("cmd.save('uniprot_ntp_analysis.pse')\n")
        f.write("print('💾 保存会话: uniprot_ntp_analysis.pse')\n")
        f.write("print('🎉 PyMOL分析完成!')\n")
    
    return script_file

def main():
    """主函数"""
    print("🧬 基于UniProt ID的NTP结合位点分析")
    print("=" * 60)
    
    # 1. 读取Excel数据
    excel_file = 'NTP_Analysis_Report.xlsx'
    uniprot_data = read_excel_data(excel_file)
    
    if not uniprot_data:
        print("❌ 没有有效的UniProt数据，退出分析")
        return
    
    # 2. 查找FASTA序列
    sequence_data = find_fasta_sequences(uniprot_data)
    
    if not sequence_data:
        print("❌ 没有找到FASTA序列，退出分析")
        return
    
    # 3. 查找PDB结构
    sequences_with_pdb = find_pdb_structures_for_sequences(sequence_data)
    
    if not sequences_with_pdb:
        print("❌ 没有找到对应的PDB结构，退出分析")
        return
    
    # 4. 下载PDB文件
    output_dir = 'ntp_analysis_output'
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, 'images'), exist_ok=True)
    
    successful_downloads = download_pdb_files(sequences_with_pdb, output_dir)
    
    if not successful_downloads:
        print("❌ 没有成功下载PDB文件，退出分析")
        return
    
    # 5. 生成PyMOL脚本
    script_file = generate_pymol_script(successful_downloads, output_dir)
    
    # 6. 生成详细报告
    report_file = os.path.join(output_dir, 'uniprot_analysis_report.txt')
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("基于UniProt ID的NTP结合位点分析报告\n")
        f.write("=" * 60 + "\n\n")
        
        f.write(f"Excel数据: {len(uniprot_data)} 条UniProt记录\n")
        f.write(f"找到FASTA: {len(sequence_data)} 个序列\n")
        f.write(f"有PDB结构: {len(sequences_with_pdb)} 个\n")
        f.write(f"成功下载: {len(successful_downloads)} 个\n\n")
        
        f.write("详细信息:\n")
        f.write("-" * 40 + "\n")
        
        for item in successful_downloads:
            f.write(f"UniProt ID: {item['uniprot_id']}\n")
            f.write(f"  文件夹: {item['folder_name']}\n")
            f.write(f"  序列长度: {item['sequence_length']} aa\n")
            f.write(f"  PDB结构: {', '.join(item['downloaded_pdbs'])}\n")
            
            # 分析每个PDB的ATP含量
            for pdb_id in item['downloaded_pdbs']:
                pdb_file = os.path.join(output_dir, 'pdb_files', f'{pdb_id}.pdb')
                if os.path.exists(pdb_file):
                    atp_ligands = analyze_pdb_for_atp(pdb_file)
                    if atp_ligands:
                        f.write(f"    {pdb_id}: ATP配体 = {', '.join(atp_ligands)}\n")
                    else:
                        f.write(f"    {pdb_id}: 无ATP配体\n")
            f.write("\n")
    
    print("\n" + "=" * 60)
    print("🎉 分析完成!")
    print("=" * 60)
    
    print(f"✅ 处理了 {len(uniprot_data)} 个UniProt条目")
    print(f"🔍 找到 {len(sequence_data)} 个FASTA序列")
    print(f"🧪 获得 {len(successful_downloads)} 个有PDB结构的条目")
    
    total_pdbs = sum(len(item['downloaded_pdbs']) for item in successful_downloads)
    print(f"📥 下载了 {total_pdbs} 个PDB文件")
    
    print(f"\n📁 结果保存在: {os.path.abspath(output_dir)}")
    print("📊 主要输出:")
    print(f"  📝 {os.path.basename(script_file)}")
    print("  📄 uniprot_analysis_report.txt")
    print("  📁 pdb_files/ (PDB文件)")
    
    print("\n🚀 下一步:")
    print("1. 安装PyMOL: conda install -c conda-forge pymol-open-source")
    print("2. 启动PyMOL")
    print(f"3. 在PyMOL中运行: run {os.path.basename(script_file)}")

if __name__ == "__main__":
    main()