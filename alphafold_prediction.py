#!/usr/bin/env python3
"""
AlphaFold 批处理脚本 - 序列长度过滤: 200-1000 AA
"""

import os
import json
import subprocess
import time
from pathlib import Path
import argparse

def check_system_environment():
    """检查系统环境"""
    print("\n💻 系统环境检查:")
    
    # 检查Python版本
    import sys
    print(f"   Python版本: {sys.version.split()[0]}")
    
    # 检查CPU信息
    try:
        import multiprocessing
        cpu_count = multiprocessing.cpu_count()
        print(f"   CPU核心数: {cpu_count}")
    except:
        print("   CPU信息: 无法获取")
    
    # 检查内存信息
    try:
        with open('/proc/meminfo', 'r') as f:
            meminfo = f.read()
            for line in meminfo.split('\n'):
                if 'MemTotal' in line:
                    mem_total = int(line.split()[1]) // 1024  # 转换为MB
                    print(f"   系统内存: {mem_total//1024:.1f} GB")
                    break
                elif 'MemAvailable' in line:
                    mem_available = int(line.split()[1]) // 1024  # 转换为MB
                    print(f"   可用内存: {mem_available//1024:.1f} GB")
                    break
    except:
        print("   内存信息: 无法获取")

def convert_cif_to_pdb(output_dir, sequence_name, global_pdb_dir="./all_pdb_files"):
    """自动将CIF文件转换为PDB文件"""
    output_path = Path(output_dir)
    
    # 查找CIF文件
    sample_0_cif = None
    all_cif_files = list(output_path.glob("**/*.cif"))
    
    for cif_file in all_cif_files:
        if "sample_0" in cif_file.name or len(all_cif_files) == 1:
            sample_0_cif = cif_file
            break
    
    if not sample_0_cif:
        print(f"⚠️  没有找到CIF文件")
        return 0
    
    print(f"🔄 发现CIF文件，开始转换: {sample_0_cif.name}")
    
    # 找到转换脚本
    script_dir = Path(__file__).parent
    converter_script = script_dir / "convert_cif_to_pdb.py"
    
    if not converter_script.exists():
        print(f"⚠️  转换脚本不存在: {converter_script}")
        return 0
    
    # 创建全局PDB目录
    global_pdb_path = Path(global_pdb_dir)
    global_pdb_path.mkdir(exist_ok=True)
    
    try:
        # 生成PDB文件名
        pdb_filename = f"{sequence_name}.pdb"
        pdb_file = global_pdb_path / pdb_filename
        
        # 调用转换脚本
        cmd = ['python3', str(converter_script), str(sample_0_cif), str(pdb_file)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        if result.returncode == 0:
            print(f"   ✅ {sample_0_cif.name} → all_pdb_files/{pdb_filename}")
            
            # 验证PDB文件内容
            if pdb_file.exists() and pdb_file.stat().st_size > 0:
                with open(pdb_file, 'r') as f:
                    atom_lines = [line for line in f if line.startswith('ATOM')]
                    print(f"   📊 生成PDB包含 {len(atom_lines)} 个原子记录")
                return 1
            else:
                print(f"   ❌ PDB文件为空")
                return 0
        else:
            print(f"   ❌ 转换失败: {result.stderr.strip()}")
            return 0
            
    except Exception as e:
        print(f"   💥 转换异常: {e}")
        return 0

def cleanup_intermediate_files(output_dir):
    """清理中间文件"""
    output_path = Path(output_dir)
    cleaned_files = 0
    saved_space = 0
    
    keep_extensions = {'.pdb', '.json', '.cif'}
    keep_patterns = ['confidence', 'summary', 'result']
    protected_dirs = {'pdb_files'}
    
    try:
        for file_path in output_path.rglob('*'):
            if file_path.is_file():
                is_protected = any(protected_dir in file_path.parts for protected_dir in protected_dirs)
                
                if is_protected:
                    continue
                
                should_keep = (
                    file_path.suffix.lower() in keep_extensions or
                    any(pattern in file_path.name.lower() for pattern in keep_patterns)
                )
                
                if not should_keep:
                    try:
                        file_size = file_path.stat().st_size
                        file_path.unlink()
                        cleaned_files += 1
                        saved_space += file_size
                    except Exception as e:
                        print(f"⚠️  无法删除 {file_path}: {e}")
        
        if cleaned_files > 0:
            print(f"🧹 清理完成: 删除 {cleaned_files} 个中间文件, 节省空间 {saved_space/1024/1024:.1f} MB")
    
    except Exception as e:
        print(f"⚠️  清理过程出错: {e}")

def get_sequence_length(json_file):
    """获取序列长度"""
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        if isinstance(data, list) and len(data) > 0:
            sequence = data[0]['sequences'][0]['proteinChain']['sequence']
        else:
            sequence = data['sequences'][0]['proteinChain']['sequence']
        
        return len(sequence)
    except:
        return 0



def run_alphafold_prediction(json_file, output_base_dir="./predictions", cleanup_intermediate=True, global_pdb_dir="./all_pdb_files"):
    """运行单个AlphaFold预测"""
    sequence_name = json_file.stem
    output_dir = Path(output_base_dir) / f"{sequence_name}_output"
    
    # 获取序列长度
    seq_length = get_sequence_length(json_file)
    
    print(f"\n🔬 预测: {json_file.name}")
    print(f"📏 序列长度: {seq_length} AA")
    
    # 构建命令 - 使用MSA服务器
    cmd = [
        'protenix', 'predict',
        '--input', str(json_file),
        '--out_dir', str(output_dir),
        '--seeds', '101',
        '--use_msa_server'
    ]
    
    print(f"🚀 命令: {' '.join(cmd)}")
    
    start_time = time.time()
    
    try:
        # 运行预测
        #result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)  # 2小时超时
        result = subprocess.run(cmd, text=True)
        
        duration = time.time() - start_time
        
        if result.returncode == 0:
            # 检查输出文件
            pdb_files = list(output_dir.glob("**/*.pdb"))
            cif_files = list(output_dir.glob("**/*.cif"))
            structure_files = list(output_dir.glob("**/*structure*"))
            model_files = list(output_dir.glob("**/*model*"))
            all_files = list(output_dir.glob("**/*"))
            
            print(f"📁 输出文件检查:")
            print(f"   PDB文件: {len(pdb_files)}个")
            print(f"   CIF文件: {len(cif_files)}个") 
            print(f"   结构文件: {len(structure_files)}个")
            print(f"   模型文件: {len(model_files)}个")
            print(f"   总文件数: {len([f for f in all_files if f.is_file()])}个")
            
            # 转换CIF文件
            converted_count = 0
            if cif_files:
                converted_count = convert_cif_to_pdb(output_dir, sequence_name, global_pdb_dir)
            else:
                print("   没有发现CIF文件，跳过转换")
            
            total_structure_files = len(pdb_files) + len(cif_files) + len(structure_files) + len(model_files)
            
            # 清理中间文件
            if cleanup_intermediate and total_structure_files > 0:
                print(f"🧹 开始清理中间文件...")
                cleanup_intermediate_files(output_dir)
            
            print(f"✅ 成功! 耗时: {duration/60:.1f}分钟, 结构文件: {total_structure_files}个")
            
            return {
                'status': 'success',
                'duration': duration,
                'pdb_count': len(pdb_files),
                'cif_count': len(cif_files),
                'structure_count': len(structure_files),
                'model_count': len(model_files),
                'converted_count': converted_count,
                'total_files': len([f for f in all_files if f.is_file()]),
                'sequence_length': seq_length,
                'cleaned_up': cleanup_intermediate and total_structure_files > 0
            }
        else:
            print(f"❌ 失败! 耗时: {duration/60:.1f}分钟")
            
            # 显示错误信息
            stderr_output = result.stderr.strip()
            stdout_output = result.stdout.strip()
            
            print(f"📤 返回码: {result.returncode}")
            
            if stdout_output:
                print(f"📋 标准输出:")
                print(stdout_output[:500] + "..." if len(stdout_output) > 500 else stdout_output)
            
            if stderr_output:
                print(f"⚠️  标准错误:")
                print(stderr_output[:500] + "..." if len(stderr_output) > 500 else stderr_output)
            
            return {
                'status': 'failed',
                'duration': duration,
                'error': stderr_output,
                'stdout': stdout_output,
                'sequence_length': seq_length
            }
    
    except subprocess.TimeoutExpired:
        print(f"⏰ 超时! (2小时)")
        return {
            'status': 'timeout',
            'duration': 7200,
            'sequence_length': seq_length
        }
    
    except Exception as e:
        print(f"💥 异常: {e}")
        return {
            'status': 'error',
            'error': str(e),
            'sequence_length': seq_length
        }

def process_batch(input_dir, start_from=0, limit=None, output_dir="./predictions", cleanup_intermediate=True, global_pdb_dir="./all_pdb_files"):
    """批处理所有JSON文件 - 序列长度200-1000 AA"""
    input_path = Path(input_dir)
    json_files = sorted(list(input_path.glob("*.json")))
    
    if not json_files:
        print(f"❌ 在 {input_dir} 中没有找到JSON文件")
        return
    
    print(f"📁 找到 {len(json_files)} 个JSON文件")
    
    # 按序列长度过滤和排序 (200-1000 AA)
    print("📏 按序列长度排序并过滤（只处理200-1000 AA）...")
    files_with_length = []
    too_short_count = 0
    too_long_count = 0
    
    for json_file in json_files:
        seq_length = get_sequence_length(json_file)
        if seq_length > 0:
            if 258 <= seq_length <= 1000:
                files_with_length.append((json_file, seq_length))
            elif seq_length < 258:
                too_short_count += 1
                if too_short_count <= 3:
                    print(f"⏩ 跳过过短序列: {json_file.name} ({seq_length} AA)")
            else:  # seq_length > 1000
                too_long_count += 1
                if too_long_count <= 3:
                    print(f"⏩ 跳过超长序列: {json_file.name} ({seq_length} AA)")
        else:
            print(f"⚠️  无法获取序列长度: {json_file.name}")
    
    if too_short_count > 3:
        print(f"⏩ ... 还有 {too_short_count - 3} 个过短序列被跳过")
    if too_long_count > 3:
        print(f"⏩ ... 还有 {too_long_count - 3} 个超长序列被跳过")
    
    print(f"🔍 过滤结果:")
    print(f"  - 跳过过短序列: {too_short_count} 个 (<200 AA)")
    print(f"  - 跳过超长序列: {too_long_count} 个 (>1000 AA)")
    print(f"  - 处理队列: {len(files_with_length)} 个序列 (200-1000 AA)")
    
    if not files_with_length:
        print("❌ 没有找到符合条件的序列")
        return
    
    # 按长度排序（从短到长）
    files_with_length.sort(key=lambda x: x[1])
    sorted_files = [item[0] for item in files_with_length]
    
    # 应用起始位置和限制
    if start_from > 0:
        sorted_files = sorted_files[start_from:]
        print(f"⏩ 从第 {start_from + 1} 个文件开始")
    
    if limit:
        sorted_files = sorted_files[:limit]
        print(f"🔢 限制处理 {limit} 个文件")
    
    total_files = len(sorted_files)
    if total_files == 0:
        print("❌ 没有文件需要处理")
        return
    
    # 创建输出目录
    Path(output_dir).mkdir(exist_ok=True)
    
    print(f"\n🚀 开始批处理 {total_files} 个序列...")
    print("=" * 80)
    
    results = []
    success_count = 0
    total_duration = 0
    start_time = time.time()
    
    for i, json_file in enumerate(sorted_files, 1):
        current_length = get_sequence_length(json_file)
        
        print(f"\n📊 进度: [{i}/{total_files}] - {json_file.name} ({current_length} AA)")
        
        result = run_alphafold_prediction(json_file, output_dir, cleanup_intermediate, global_pdb_dir)
        result['file'] = json_file.name
        results.append(result)
        
        if result['status'] == 'success':
            success_count += 1
            total_duration += result['duration']
            avg_duration = total_duration / success_count
            
            remaining_files = total_files - i
            estimated_remaining_time = (remaining_files * avg_duration) / 3600
            
            print(f"📈 统计: 成功 {success_count}/{i}, 平均 {avg_duration/60:.1f}分钟/个")
            print(f"⏱️  预计剩余: {estimated_remaining_time:.1f}小时")
        
        print("-" * 50)
    
    # 最终统计
    batch_duration = time.time() - start_time
    failed_count = total_files - success_count
    
    print(f"\n🏁 批处理完成!")
    print("=" * 80)
    print(f"📊 总结:")
    print(f"  - 总文件数: {total_files}")
    print(f"  - 成功: {success_count} ({success_count/total_files*100:.1f}%)")
    print(f"  - 失败: {failed_count} ({failed_count/total_files*100:.1f}%)")
    print(f"  - 总耗时: {batch_duration/3600:.1f} 小时")
    if success_count > 0:
        print(f"  - 平均预测时间: {total_duration/success_count/60:.1f} 分钟")
    
    # 保存结果
    save_results(results, output_dir)

def save_results(results, output_dir):
    """保存结果报告"""
    timestamp = time.strftime('%Y%m%d_%H%M%S')
    report_file = Path(output_dir) / f"batch_results_{timestamp}.json"
    
    success_results = [r for r in results if r['status'] == 'success']
    
    report = {
        'metadata': {
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'total_files': len(results),
            'success_count': len(success_results),
            'failed_count': len(results) - len(success_results),
            'length_filter': '200-1000 AA only'
        },
        'statistics': {
            'avg_duration': sum(r.get('duration', 0) for r in success_results) / len(success_results) if success_results else 0,
            'total_prediction_time': sum(r.get('duration', 0) for r in success_results),
            'sequence_length_stats': {
                'min': min(r.get('sequence_length', 0) for r in results if r.get('sequence_length', 0) > 0) if results else 0,
                'max': max(r.get('sequence_length', 0) for r in results),
                'avg': sum(r.get('sequence_length', 0) for r in results) / len(results) if results else 0
            }
        },
        'results': results
    }
    
    with open(report_file, 'w', encoding='utf-8') as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    
    print(f"📄 详细报告已保存: {report_file}")

def main():
    parser = argparse.ArgumentParser(description='AlphaFold批处理脚本 (200-1000 AA)')
    parser.add_argument('--input_dir', '-i', default='/home/webserver/student/students_webserver/zhijing/input_jsons', help='JSON文件输入目录')
    parser.add_argument('--output_dir', '-o', default='./alphafold_predictions', help='输出目录')
    parser.add_argument('--start_from', '-f', type=int, default=0, help='从第几个文件开始')
    parser.add_argument('--limit', '-l', type=int, help='处理文件数量限制')
    parser.add_argument('--no_cleanup', action='store_true', help='不清理中间文件')
    
    args = parser.parse_args()
    
    print("💻 AlphaFold批处理器 (200-1000 AA)")
    print("📏 序列长度过滤: 200-1000 氨基酸")
    print("=" * 60)
    
    # 检查系统环境
    check_system_environment()
    
    print("🚀 准备开始批量处理...")
    
    if not os.path.exists(args.input_dir):
        print(f"❌ 输入目录不存在: {args.input_dir}")
        return
    
    print(f"\n📁 输入目录: {args.input_dir}")
    print(f"📤 输出目录: {args.output_dir}")
    print(f"🔢 处理范围: 从第{args.start_from + 1}个开始" + (f", 限制{args.limit}个" if args.limit else ", 无限制"))
    print(f"🧹 自动清理: {'禁用' if args.no_cleanup else '启用'}")
    print("=" * 60)
    
    # 开始批处理
    process_batch(
        input_dir=args.input_dir,
        start_from=args.start_from,
        limit=args.limit,
        output_dir=args.output_dir,
        cleanup_intermediate=not args.no_cleanup,
        global_pdb_dir="./all_pdb_files"
    )

if __name__ == "__main__":
    main()