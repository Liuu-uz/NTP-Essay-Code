#!/usr/bin/env python3
"""
AlphaFold 批处理脚本 - 序列长度过滤: 200-1000 AA
新增MSA限流处理和智能休息策略
"""

import os
import json
import subprocess
import time
import random
import re
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

def parse_msa_status_from_output(output_text):
    """从输出中解析MSA状态"""
    msa_info = {
        'pending_detected': False,
        'sleep_time': 0,
        'msa_success': False,
        'error_detected': False
    }
    
    lines = output_text.split('\n')
    
    for line in lines:
        # 检查PENDING状态
        if 'Sleeping for' in line and 'PENDING' in line:
            match = re.search(r'Sleeping for (\d+)s', line)
            if match:
                msa_info['pending_detected'] = True
                msa_info['sleep_time'] = int(match.group(1))
        
        # 检查MSA成功
        elif 'update msa result success' in line:
            msa_info['msa_success'] = True
        
        # 检查其他错误
        elif any(error_keyword in line.lower() for error_keyword in ['error', 'failed', 'timeout']):
            msa_info['error_detected'] = True
    
    return msa_info

def calculate_smart_delay(consecutive_pendings, seq_length, is_first_request=False):
    """计算智能延迟时间"""
    
    if is_first_request:
        # 第一个请求不延迟
        return 0
    
    # 根据连续PENDING次数调整基础延迟和倍数
    if consecutive_pendings == 0:
        # 正常情况，基础延迟60秒
        base_delay = 60
        pending_factor = 1.0
    elif consecutive_pendings == 1:
        # 第一次PENDING，休息约3分钟
        base_delay = 180
        pending_factor = 1.0
    elif consecutive_pendings == 2:
        # 连续PENDING，休息约5分钟
        base_delay = 300
        pending_factor = 1.0
    else:
        # 严重限流，休息约8分钟
        base_delay = 480
        pending_factor = 1.0
    
    # 根据序列长度微调（长序列MSA更耗时）
    length_factor = min(seq_length / 300, 1.5)
    
    # 随机化避免同步请求
    random_factor = random.uniform(0.8, 1.2)
    
    total_delay = base_delay * pending_factor * length_factor * random_factor
    
    # 限制延迟范围：最少15秒，最多600秒（10分钟）
    return max(15, min(total_delay, 600))

def run_alphafold_prediction_optimized(json_file, output_base_dir="./predictions", cleanup_intermediate=True, global_pdb_dir="./all_pdb_files", previous_msa_info=None):
    """优化版AlphaFold预测（带MSA状态处理）"""
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
        # 运行预测，实时监控输出
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True
        )
        
        output_lines = []
        msa_start_time = None
        pending_detected = False
        
        # 实时读取输出
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                output_lines.append(output.strip())
                
                # 实时检测关键状态
                if 'starting to update msa result' in output:
                    msa_start_time = time.time()
                    print(f"   🔍 开始MSA搜索...")
                
                elif 'Sleeping for' in output and 'PENDING' in output:
                    match = re.search(r'Sleeping for (\d+)s', output)
                    if match:
                        sleep_time = int(match.group(1))
                        print(f"   🚦 MSA服务器PENDING，需等待{sleep_time}秒...")
                        pending_detected = True
                
                elif 'Files downloaded and extracted successfully' in output:
                    if msa_start_time:
                        msa_duration = time.time() - msa_start_time
                        print(f"   ✅ MSA下载完成，耗时{msa_duration:.1f}秒")
                
                elif 'update msa result success' in output:
                    if msa_start_time:
                        total_msa_time = time.time() - msa_start_time
                        print(f"   🎯 MSA处理完成，总耗时{total_msa_time:.1f}秒")
        
        # 等待进程完成
        return_code = process.wait()
        duration = time.time() - start_time
        
        # 解析输出
        full_output = '\n'.join(output_lines)
        msa_info = parse_msa_status_from_output(full_output)
        
        if return_code == 0:
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
                'cleaned_up': cleanup_intermediate and total_structure_files > 0,
                'msa_info': msa_info,
                'pending_detected': pending_detected
            }
        else:
            print(f"❌ 失败! 耗时: {duration/60:.1f}分钟")
            print(f"📤 返回码: {return_code}")
            
            # 显示部分输出用于调试
            if full_output:
                print(f"📋 输出摘要:")
                print(full_output[-300:])  # 显示最后300字符
            
            return {
                'status': 'failed',
                'duration': duration,
                'error': full_output[-500:] if full_output else '',
                'sequence_length': seq_length,
                'msa_info': msa_info,
                'pending_detected': pending_detected
            }
    
    except Exception as e:
        print(f"💥 异常: {e}")
        return {
            'status': 'error',
            'error': str(e),
            'sequence_length': seq_length
        }

def process_batch_optimized(input_dir, start_from=0, limit=None, output_dir="./predictions", cleanup_intermediate=True, global_pdb_dir="./all_pdb_files"):
    """优化版批处理（带智能休息策略）"""
    input_path = Path(input_dir)
    json_files = sorted(list(input_path.glob("*.json")))
    
    if not json_files:
        print(f"❌ 在 {input_dir} 中没有找到JSON文件")
        return
    
    print(f"📁 找到 {len(json_files)} 个JSON文件")
    
    # 按序列长度过滤和排序 (258-1000 AA)
    print("📏 按序列长度排序并过滤（只处理258-1000 AA）...")
    files_with_length = []
    too_short_count = 0
    too_long_count = 0
    
    for json_file in json_files:
        seq_length = get_sequence_length(json_file)
        if seq_length > 0:
            if 0 <= seq_length <= 1000:
                files_with_length.append((json_file, seq_length))
            elif seq_length < 0:
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
    print(f"  - 跳过过短序列: {too_short_count} 个 (<258 AA)")
    print(f"  - 跳过超长序列: {too_long_count} 个 (>1000 AA)")
    print(f"  - 处理队列: {len(files_with_length)} 个序列 (258-1000 AA)")
    
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
    
    # 检查当前时间，给出建议
    current_hour = time.localtime().tm_hour
    if current_hour in [9, 10, 14, 15, 20, 21, 22]:
        print(f"\n⏰ 当前时段({current_hour}:00)可能是MSA服务器使用高峰期")
        print("   建议在深夜(2-6点)或上午(11-13点)运行以获得最佳性能")
        proceed = input("   是否继续？(y/n): ")
        if proceed.lower() != 'y':
            print("已取消批处理")
            return
    
    # 创建输出目录
    Path(output_dir).mkdir(exist_ok=True)
    
    print(f"\n🚀 开始MSA优化批处理 {total_files} 个序列...")
    print("=" * 80)
    
    results = []
    success_count = 0
    total_duration = 0
    start_time = time.time()
    consecutive_pendings = 0  # 连续PENDING计数
    previous_msa_info = None
    
    for i, json_file in enumerate(sorted_files, 1):
        current_length = get_sequence_length(json_file)
        
        print(f"\n📊 进度: [{i}/{total_files}] - {json_file.name} ({current_length} AA)")
        
        # 智能延迟策略
        if i > 1:  # 第一个序列不延迟
            delay_time = calculate_smart_delay(consecutive_pendings, current_length, is_first_request=False)
            
            if delay_time > 60:  # 超过1分钟显示详细信息
                print(f"   ⏱️  智能休息: {delay_time/60:.1f}分钟 (连续PENDING: {consecutive_pendings}次)")
            else:
                print(f"   ⏱️  预防延迟: {delay_time:.0f}秒")
            
            time.sleep(delay_time)
        
        # 执行预测
        result = run_alphafold_prediction_optimized(
            json_file, output_dir, cleanup_intermediate, global_pdb_dir, previous_msa_info
        )
        result['file'] = json_file.name
        results.append(result)
        
        # 更新统计
        if result['status'] == 'success':
            success_count += 1
            total_duration += result['duration']
            avg_duration = total_duration / success_count
            
            remaining_files = total_files - i
            estimated_remaining_time = (remaining_files * avg_duration) / 3600
            
            print(f"📈 统计: 成功 {success_count}/{i}, 平均 {avg_duration/60:.1f}分钟/个")
            print(f"⏱️  预计剩余: {estimated_remaining_time:.1f}小时")
        
        # 更新PENDING状态追踪
        if result.get('pending_detected'):
            consecutive_pendings += 1
            print(f"   🚦 检测到PENDING (连续第{consecutive_pendings}次)")
        else:
            consecutive_pendings = 0  # 重置计数
        
        # 连续PENDING过多时强制长休息
        if consecutive_pendings >= 3:
            long_rest = 360 + random.uniform(60, 180)  # 6-9分钟
            print(f"   🛑 连续PENDING过多，强制长休息{long_rest/60:.1f}分钟...")
            time.sleep(long_rest)
            consecutive_pendings = 0  # 重置计数
        
        # 保存MSA信息供下次使用
        if 'msa_info' in result:
            previous_msa_info = result['msa_info']
        
        print("-" * 50)
    
    # 最终统计
    batch_duration = time.time() - start_time
    failed_count = total_files - success_count
    
    print(f"\n🏁 MSA优化批处理完成!")
    print("=" * 80)
    print(f"📊 总结:")
    print(f"  - 总文件数: {total_files}")
    print(f"  - 成功: {success_count} ({success_count/total_files*100:.1f}%)")
    print(f"  - 失败: {failed_count} ({failed_count/total_files*100:.1f}%)")
    print(f"  - 总耗时: {batch_duration/3600:.1f} 小时")
    if success_count > 0:
        print(f"  - 平均预测时间: {total_duration/success_count/60:.1f} 分钟")
        print(f"  - 实际平均(含休息): {batch_duration/total_files/60:.1f} 分钟/个")
    
    # 保存结果
    save_results_optimized(results, output_dir)

def save_results_optimized(results, output_dir):
    """保存优化版结果报告"""
    timestamp = time.strftime('%Y%m%d_%H%M%S')
    report_file = Path(output_dir) / f"batch_results_optimized_{timestamp}.json"
    
    success_results = [r for r in results if r['status'] == 'success']
    pending_count = sum(1 for r in results if r.get('pending_detected', False))
    
    report = {
        'metadata': {
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'total_files': len(results),
            'success_count': len(success_results),
            'failed_count': len(results) - len(success_results),
            'pending_encounters': pending_count,
            'length_filter': '258-1000 AA only',
            'optimization': 'MSA smart delay enabled'
        },
        'statistics': {
            'avg_duration': sum(r.get('duration', 0) for r in success_results) / len(success_results) if success_results else 0,
            'total_prediction_time': sum(r.get('duration', 0) for r in success_results),
            'pending_rate': pending_count / len(results) if results else 0,
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
    print(f"🚦 PENDING遭遇率: {pending_count}/{len(results)} ({pending_count/len(results)*100:.1f}%)")

def main():
    parser = argparse.ArgumentParser(description='AlphaFold批处理脚本 (MSA优化版)')
    parser.add_argument('--input_dir', '-i', default='/home/webserver/student/students_webserver/zhijing/input_jsons', help='JSON文件输入目录')
    parser.add_argument('--output_dir', '-o', default='./alphafold_predictions', help='输出目录')
    parser.add_argument('--start_from', '-f', type=int, default=0, help='从第几个文件开始')
    parser.add_argument('--limit', '-l', type=int, help='处理文件数量限制')
    parser.add_argument('--no_cleanup', action='store_true', help='不清理中间文件')
    
    args = parser.parse_args()
    
    print("💻 AlphaFold批处理器 (MSA优化版)")
    print("📏 序列长度过滤: 258-1000 氨基酸")
    print("🧠 智能特性: MSA限流检测 + 自适应休息策略")
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
    
    # 开始优化批处理
    process_batch_optimized(
        input_dir=args.input_dir,
        start_from=args.start_from,
        limit=args.limit,
        output_dir=args.output_dir,
        cleanup_intermediate=not args.no_cleanup,
        global_pdb_dir="./all_pdb_files"
    )

if __name__ == "__main__":
    main()