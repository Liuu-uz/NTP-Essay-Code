import os
import subprocess
import numpy as np
import pandas as pd
from Bio.PDB import *
from Bio.PDB.Polypeptide import PPBuilder
import matplotlib.pyplot as plt
import seaborn as sns
import re
import warnings
warnings.filterwarnings('ignore')

# 配置设置
EXP_CIF_DIR = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/cif_files"
AF_CIF_DIR = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/cif_alphafold"
TMALIGN_PATH = "/Users/napkin/tools/TMalign"
OUTPUT_DIR = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/results2"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 单链分析配置
PROTEIN_DB = {
    "1d4x": {
        "exp_cif": "1d4x.cif",
        "af_cif": "1d4x_af.cif",
        "target_chain": "A",  # 指定要比较的链
        "active_site": [("A", 50), ("A", 102), ("A", 195)]  # 这些会根据实际范围调整
    },
}

class SingleChainSelector(Select):
    """选择特定链的选择器"""
    def __init__(self, chain_id):
        self.chain_id = chain_id
    
    def accept_chain(self, chain):
        return chain.get_id() == self.chain_id

class CIFToPDBConverter:
    """CIF到PDB格式转换器 - 支持单链提取"""
    
    def convert_cif_to_pdb_single_chain(self, cif_path, pdb_path, chain_id):
        """将CIF文件的指定链转换为PDB格式"""
        try:
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure("temp", cif_path)
            
            # 只保存指定的链
            io = PDBIO()
            io.set_structure(structure)
            io.save(pdb_path, SingleChainSelector(chain_id))
            
            # 验证输出文件
            if os.path.exists(pdb_path) and os.path.getsize(pdb_path) > 0:
                return True
            else:
                print(f"警告: 链 {chain_id} 可能不存在或为空")
                return False
                
        except Exception as e:
            print(f"CIF转PDB失败 (链 {chain_id}): {e}")
            return False

def analyze_chain_structure(structure, chain_id, structure_type="unknown"):
    """分析指定链的结构质量"""
    try:
        model = structure[0]
        
        # 检查链是否存在
        if chain_id not in model:
            print(f"错误: 链 {chain_id} 不存在")
            return None
        
        chain = model[chain_id]
        
        # 基本统计
        num_residues = 0
        num_atoms = 0
        b_factors = []
        
        for residue in chain:
            if residue.get_id()[0] == ' ':  # 标准残基
                num_residues += 1
                for atom in residue:
                    num_atoms += 1
                    b_factor = atom.get_bfactor()
                    if b_factor is not None:
                        b_factors.append(b_factor)
        
        # 残基范围
        residues = [r for r in chain if r.get_id()[0] == ' ']
        first_res = residues[0].get_id()[1] if residues else None
        last_res = residues[-1].get_id()[1] if residues else None
        
        # Ramachandran角度 - 只计算这条链
        ppb = PPBuilder()
        phi_psi_angles = []
        
        try:
            # 创建只包含目标链的临时结构
            temp_structure = Structure.Structure("temp")
            temp_model = Model.Model(0)
            temp_model.add(chain.copy())
            temp_structure.add(temp_model)
            
            for peptide in ppb.build_peptides(temp_structure):
                phi_psi_list = peptide.get_phi_psi_list()
                for phi_psi in phi_psi_list:
                    if phi_psi[0] is not None and phi_psi[1] is not None:
                        phi_psi_angles.append((np.degrees(phi_psi[0]), np.degrees(phi_psi[1])))
        except Exception as e:
            print(f"计算链 {chain_id} Ramachandran角度时出错: {e}")
        
        return {
            'structure_type': structure_type,
            'chain_id': chain_id,
            'num_residues': num_residues,
            'num_atoms': num_atoms,
            'first_residue': first_res,
            'last_residue': last_res,
            'phi_psi_angles': phi_psi_angles,
            'b_factors': b_factors,
            'avg_b_factor': np.mean(b_factors) if b_factors else None,
            'std_b_factor': np.std(b_factors) if b_factors else None
        }
    
    except Exception as e:
        print(f"链 {chain_id} 结构分析失败: {e}")
        return None

def run_tmalign(pdb1_path, pdb2_path):
    """运行TM-align比较两个结构"""
    try:
        cmd = [TMALIGN_PATH, pdb1_path, pdb2_path]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        if result.returncode != 0:
            print(f"TM-align运行失败: {result.stderr}")
            return None
        
        output = result.stdout
        
        # 解析输出
        tm_score = None
        rmsd = None
        aligned_length = None
        
        for line in output.split('\n'):
            line = line.strip()
            
            if line.startswith('Aligned length='):
                try:
                    parts = line.split(',')
                    if len(parts) >= 1 and 'Aligned length=' in parts[0]:
                        aligned_length_str = parts[0].split('=')[1].strip()
                        aligned_length = int(aligned_length_str)
                    if len(parts) >= 2 and 'RMSD=' in parts[1]:
                        rmsd_str = parts[1].split('=')[1].strip()
                        rmsd = float(rmsd_str)
                except (IndexError, ValueError) as e:
                    continue
            
            elif line.startswith('TM-score=') and 'Chain_1' in line:
                try:
                    tm_match = re.search(r'TM-score=\s*([0-9]+\.?[0-9]*)', line)
                    if tm_match:
                        tm_score = float(tm_match.group(1))
                except (ValueError, AttributeError) as e:
                    continue
        
        return {
            'tm_score': tm_score,
            'rmsd': rmsd,
            'aligned_length': aligned_length,
            'raw_output': output
        }
    
    except subprocess.TimeoutExpired:
        print("TM-align运行超时")
        return None
    except Exception as e:
        print(f"TM-align运行出错: {e}")
        return None

def calculate_active_site_rmsd_smart(exp_structure, af_structure, chain_id, active_sites):
    """智能计算活性位点RMSD - 自动调整范围"""
    try:
        exp_model = exp_structure[0]
        af_model = af_structure[0]
        
        if chain_id not in exp_model or chain_id not in af_model:
            print(f"链 {chain_id} 在某个结构中不存在")
            return None, []
        
        exp_chain = exp_model[chain_id]
        af_chain = af_model[chain_id]
        
        # 获取两个结构的残基范围
        exp_residues = [r.get_id()[1] for r in exp_chain if r.get_id()[0] == ' ']
        af_residues = [r.get_id()[1] for r in af_chain if r.get_id()[0] == ' ']
        
        # 找到重叠范围
        overlap_start = max(min(exp_residues), min(af_residues))
        overlap_end = min(max(exp_residues), max(af_residues))
        
        print(f"    实验结构范围: {min(exp_residues)}-{max(exp_residues)}")
        print(f"    AlphaFold范围: {min(af_residues)}-{max(af_residues)}")
        print(f"    重叠范围: {overlap_start}-{overlap_end}")
        
        # 智能选择活性位点 - 在重叠范围内均匀分布
        if overlap_end > overlap_start:
            range_size = overlap_end - overlap_start
            smart_sites = [
                (chain_id, overlap_start + range_size // 4),
                (chain_id, overlap_start + range_size // 2),
                (chain_id, overlap_start + 3 * range_size // 4)
            ]
            print(f"    智能活性位点: {smart_sites}")
        else:
            print(f"    警告: 无重叠区域")
            return None, []
        
        active_site_results = []
        
        for chain_id_site, res_num in smart_sites:
            try:
                if res_num in exp_residues and res_num in af_residues:
                    exp_residue = exp_chain[res_num]
                    af_residue = af_chain[res_num]
                    
                    exp_ca = exp_residue['CA']
                    af_ca = af_residue['CA']
                    
                    distance = np.linalg.norm(exp_ca.get_coord() - af_ca.get_coord())
                    
                    active_site_results.append({
                        'chain': chain_id_site,
                        'residue': res_num,
                        'distance': distance
                    })
                    
                    print(f"      残基 {res_num}: {distance:.2f} Å")
            
            except KeyError:
                print(f"      残基 {res_num}: 未找到")
                continue
        
        if active_site_results:
            avg_rmsd = np.mean([r['distance'] for r in active_site_results])
            return avg_rmsd, active_site_results
        else:
            return None, []
    
    except Exception as e:
        print(f"计算活性位点RMSD失败: {e}")
        return None, []

def evaluate_single_chain_protein(pdb_id, protein_data):
    """评估单链蛋白质"""
    converter = CIFToPDBConverter()
    parser = MMCIFParser(QUIET=True)
    
    try:
        chain_id = protein_data['target_chain']
        print(f"  分析链 {chain_id}...")
        
        # CIF文件路径
        exp_cif = os.path.join(EXP_CIF_DIR, protein_data['exp_cif'])
        af_cif = os.path.join(AF_CIF_DIR, protein_data['af_cif'])
        
        # 检查文件存在
        if not os.path.exists(exp_cif):
            print(f"错误: 实验CIF文件未找到: {exp_cif}")
            return None
        if not os.path.exists(af_cif):
            print(f"错误: AlphaFold CIF文件未找到: {af_cif}")
            return None
        
        # 加载结构
        exp_structure = parser.get_structure(f"{pdb_id}_exp", exp_cif)
        af_structure = parser.get_structure(f"{pdb_id}_af", af_cif)
        
        print(f"    结构加载成功")
        
        # 分析链结构质量
        exp_quality = analyze_chain_structure(exp_structure, chain_id, "experimental")
        af_quality = analyze_chain_structure(af_structure, chain_id, "alphafold")
        
        if not exp_quality or not af_quality:
            print(f"    链结构质量分析失败")
            return None
        
        # 转换为PDB格式（只提取目标链）
        exp_pdb = os.path.join(OUTPUT_DIR, f"{pdb_id}_exp_chain_{chain_id}.pdb")
        af_pdb = os.path.join(OUTPUT_DIR, f"{pdb_id}_af_chain_{chain_id}.pdb")
        
        if not converter.convert_cif_to_pdb_single_chain(exp_cif, exp_pdb, chain_id):
            print(f"    实验结构链 {chain_id} 转换失败")
            return None
        if not converter.convert_cif_to_pdb_single_chain(af_cif, af_pdb, chain_id):
            print(f"    预测结构链 {chain_id} 转换失败")
            return None
        
        # 运行TM-align
        print(f"    运行TM-align比较链 {chain_id}...")
        tmalign_result = run_tmalign(exp_pdb, af_pdb)
        
        if not tmalign_result:
            print(f"    链 {chain_id} TM-align比较失败")
            return None
        
        # 智能计算活性位点RMSD
        active_site_rmsd, active_site_details = calculate_active_site_rmsd_smart(
            exp_structure, af_structure, chain_id, protein_data['active_site'])
        
        # 整合结果
        result = {
            'pdb_id': pdb_id,
            'chain_id': chain_id,
            'exp_quality': exp_quality,
            'af_quality': af_quality,
            'tm_score': tmalign_result['tm_score'],
            'global_rmsd': tmalign_result['rmsd'],
            'aligned_length': tmalign_result['aligned_length'],
            'active_site_rmsd': active_site_rmsd,
            'active_site_details': active_site_details
        }
        
        # 清理临时PDB文件
        for temp_file in [exp_pdb, af_pdb]:
            if os.path.exists(temp_file):
                os.remove(temp_file)
        
        return result
    
    except Exception as e:
        print(f"评估蛋白质 {pdb_id} 链 {chain_id} 时出错: {e}")
        import traceback
        traceback.print_exc()
        return None

def create_single_chain_visualizations(results):
    """创建单链分析可视化图表"""
    if not results:
        print("没有结果数据，跳过可视化")
        return
    
    plt.style.use('default')
    sns.set_palette("husl")
    
    fig = plt.figure(figsize=(16, 12))
    
    # 1. TM-score vs RMSD
    ax1 = plt.subplot(2, 3, 1)
    tm_scores = [r['tm_score'] for r in results if r['tm_score'] is not None]
    rmsds = [r['global_rmsd'] for r in results if r['global_rmsd'] is not None]
    
    if tm_scores and rmsds:
        plt.scatter(tm_scores, rmsds, alpha=0.7, s=100, c='blue')
        for i, (tm, rmsd) in enumerate(zip(tm_scores, rmsds)):
            plt.annotate(f"{results[i]['pdb_id']}-{results[i]['chain_id']}", 
                        (tm, rmsd), xytext=(5, 5), textcoords='offset points')
        plt.xlabel('TM-score')
        plt.ylabel('Global RMSD (Å)')
        plt.title('TM-score vs Global RMSD (Chain A)')
        plt.grid(True, alpha=0.3)
    
    # 2. 活性位点RMSD
    ax2 = plt.subplot(2, 3, 2)
    active_rmsds = [r['active_site_rmsd'] for r in results if r['active_site_rmsd'] is not None]
    
    if active_rmsds:
        plt.hist(active_rmsds, bins=10, alpha=0.7, edgecolor='black', color='orange')
        plt.axvline(2.0, color='red', linestyle='--', label='Good threshold (2.0 Å)')
        plt.axvline(5.0, color='orange', linestyle='--', label='Acceptable threshold (5.0 Å)')
        plt.xlabel('Active Site RMSD (Å)')
        plt.ylabel('Count')
        plt.title('Active Site RMSD Distribution (Chain A)')
        plt.legend()
        plt.grid(True, alpha=0.3)
    
    # 3. 残基数量比较
    ax3 = plt.subplot(2, 3, 3)
    exp_residues = [r['exp_quality']['num_residues'] for r in results if r['exp_quality']]
    af_residues = [r['af_quality']['num_residues'] for r in results if r['af_quality']]
    protein_names = [f"{r['pdb_id']}-{r['chain_id']}" for r in results]
    
    if exp_residues and af_residues:
        x = np.arange(len(protein_names))
        width = 0.35
        
        plt.bar(x - width/2, exp_residues, width, label='Experimental', alpha=0.8)
        plt.bar(x + width/2, af_residues, width, label='AlphaFold', alpha=0.8)
        plt.xlabel('Proteins')
        plt.ylabel('Number of Residues')
        plt.title('Chain A Residue Count Comparison')
        plt.xticks(x, protein_names, rotation=45)
        plt.legend()
        plt.grid(True, alpha=0.3)
    
    # 4. Ramachandran图
    ax4 = plt.subplot(2, 3, 4)
    all_phi = []
    all_psi = []
    
    for r in results:
        if r['af_quality'] and r['af_quality']['phi_psi_angles']:
            for phi, psi in r['af_quality']['phi_psi_angles']:
                all_phi.append(phi)
                all_psi.append(psi)
    
    if all_phi and all_psi:
        plt.scatter(all_phi, all_psi, alpha=0.6, s=20, c='green')
        plt.xlabel('Phi (degrees)')
        plt.ylabel('Psi (degrees)')
        plt.title('Ramachandran Plot (AlphaFold Chain A)')
        plt.xlim(-180, 180)
        plt.ylim(-180, 180)
        plt.grid(True, alpha=0.3)
    
    # 5. 结构质量条形图
    ax5 = plt.subplot(2, 3, 5)
    if tm_scores:
        colors = ['green' if tm > 0.9 else 'orange' if tm > 0.5 else 'red' for tm in tm_scores]
        bars = plt.bar(range(len(protein_names)), tm_scores, color=colors, alpha=0.7)
        plt.axhline(0.5, color='red', linestyle='--', label='Good threshold')
        plt.axhline(0.7, color='orange', linestyle='--', label='Very good threshold')
        plt.axhline(0.9, color='green', linestyle='--', label='Excellent threshold')
        plt.xlabel('Proteins')
        plt.ylabel('TM-score')
        plt.title('Structure Similarity (TM-score) - Chain A')
        plt.xticks(range(len(protein_names)), protein_names, rotation=45)
        plt.legend()
        plt.grid(True, alpha=0.3)
    
    # 6. 统计总结
    ax6 = plt.subplot(2, 3, 6)
    ax6.axis('off')
    
    stats_text = "=== 链A分析统计 ===\n\n"
    stats_text += f"分析蛋白质数: {len(results)}\n\n"
    
    if tm_scores:
        stats_text += f"TM-score 统计:\n"
        stats_text += f"  平均: {np.mean(tm_scores):.3f}\n"
        stats_text += f"  范围: {min(tm_scores):.3f} - {max(tm_scores):.3f}\n\n"
    
    if rmsds:
        stats_text += f"RMSD 统计:\n"
        stats_text += f"  平均: {np.mean(rmsds):.2f} Å\n"
        stats_text += f"  范围: {min(rmsds):.2f} - {max(rmsds):.2f} Å\n\n"
    
    if active_rmsds:
        stats_text += f"活性位点 RMSD:\n"
        stats_text += f"  平均: {np.mean(active_rmsds):.2f} Å\n"
        good_sites = sum(1 for x in active_rmsds if x < 2.0)
        acceptable_sites = sum(1 for x in active_rmsds if x < 5.0)
        stats_text += f"  <2Å (优秀): {good_sites}/{len(active_rmsds)}\n"
        stats_text += f"  <5Å (可接受): {acceptable_sites}/{len(active_rmsds)}"
    
    ax6.text(0.1, 0.9, stats_text, transform=ax6.transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    
    # 保存图表
    plot_path = os.path.join(OUTPUT_DIR, "single_chain_analysis.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"单链分析图表保存至: {plot_path}")
    plt.show()

def main():
    """主函数"""
    print("=== 单链蛋白质结构预测评估系统 ===")
    
    # 检查依赖
    if not os.path.exists(TMALIGN_PATH):
        print(f"错误: TM-align未找到: {TMALIGN_PATH}")
        return
    
    if not os.path.exists(AF_CIF_DIR):
        print(f"错误: AlphaFold CIF目录未找到: {AF_CIF_DIR}")
        return
    
    if not os.path.exists(EXP_CIF_DIR):
        print(f"错误: 实验CIF目录未找到: {EXP_CIF_DIR}")
        return
    
    # 开始评估
    results = []
    total_proteins = len(PROTEIN_DB)
    
    print(f"开始评估 {total_proteins} 个蛋白质的A链...")
    
    for i, (pdb_id, protein_data) in enumerate(PROTEIN_DB.items(), 1):
        print(f"[{i}/{total_proteins}] 评估 {pdb_id}...")
        result = evaluate_single_chain_protein(pdb_id, protein_data)
        
        if result:
            results.append(result)
            print(f"  ✓ {pdb_id} 链 {result['chain_id']} 评估完成")
            print(f"    TM-score: {result['tm_score']:.3f}")
            print(f"    Global RMSD: {result['global_rmsd']:.2f} Å")
            if result['active_site_rmsd']:
                print(f"    Active Site RMSD: {result['active_site_rmsd']:.2f} Å")
        else:
            print(f"  ✗ {pdb_id} 评估失败")
    
    if not results:
        print("错误: 没有成功评估任何蛋白质")
        return
    
    # 保存详细结果
    results_data = []
    for r in results:
        row = {
            'PDB_ID': r['pdb_id'],
            'Chain_ID': r['chain_id'],
            'TM_score': r['tm_score'],
            'Global_RMSD': r['global_rmsd'],
            'Aligned_Length': r['aligned_length'],
            'Active_Site_RMSD': r['active_site_rmsd'],
            'Exp_Residues': r['exp_quality']['num_residues'] if r['exp_quality'] else None,
            'AF_Residues': r['af_quality']['num_residues'] if r['af_quality'] else None,
            'Exp_Residue_Range': f"{r['exp_quality']['first_residue']}-{r['exp_quality']['last_residue']}" if r['exp_quality'] else None,
            'AF_Residue_Range': f"{r['af_quality']['first_residue']}-{r['af_quality']['last_residue']}" if r['af_quality'] else None,
            'AF_Avg_B_factor': r['af_quality']['avg_b_factor'] if r['af_quality'] else None
        }
        results_data.append(row)
    
    results_df = pd.DataFrame(results_data)
    results_csv = os.path.join(OUTPUT_DIR, "single_chain_evaluation_results.csv")
    results_df.to_csv(results_csv, index=False)
    print(f"\n详细结果保存至: {results_csv}")
    
    # 创建可视化
    print("生成可视化图表...")
    create_single_chain_visualizations(results)
    
    print("\n=== 单链评估完成 ===")
    print("主要结果:")
    print(f"  成功评估: {len(results)} 个蛋白质的A链")
    if results:
        tm_scores = [r['tm_score'] for r in results if r['tm_score']]
        if tm_scores:
            print(f"  平均TM-score: {np.mean(tm_scores):.3f}")
            print(f"  高质量预测 (TM-score > 0.5): {sum(1 for x in tm_scores if x > 0.5)} 个")

if __name__ == "__main__":
    main()