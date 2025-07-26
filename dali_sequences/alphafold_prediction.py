#!/usr/bin/env python3
"""
RTX A4000 简化批处理脚本 - 专为单GPU优化
"""

import os
import json
import subprocess
import time
from pathlib import Path
import argparse

def setup_gpu_environment():
    """设置RTX A4000优化环境"""
    os.environ['CUDA_VISIBLE_DEVICES'] = '0'
    os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'max_split_size_mb:512'
    os.environ['CUDA_LAUNCH_BLOCKING'] = '0'
    print("🎮 GPU环境已配置: RTX A4000")

def cleanup_intermediate_files(output_dir):
    """清理中间文件，只保留PDB文件和重要结果"""
    output_path = Path(output_dir)
    cleaned_files = 0
    saved_space = 0
    
    # 要保留的文件类型
    keep_extensions = {'.pdb', '.json'}  # PDB文件和JSON结果文件
    keep_patterns = ['confidence', 'summary', 'result']  # 包含这些关键词的文件
    
    try:
        for file_path in output_path.rglob('*'):
            if file_path.is_file():
                # 检查是否应该保留
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

def get_optimal_params(sequence_length):
    """根据RTX A4000的16GB显存优化参数"""
    if sequence_length < 100:
        # 短序列 - 可以用更高精度
        return {
            'cycle': 4,
            'step': 150,
            'sample': 1,
            'seeds': '101,102',
            'mode': '平衡模式'
        }
    elif sequence_length < 300:
        # 中等序列 - 快速模式
        return {
            'cycle': 3,
            'step': 100,
            'sample': 1,
            'seeds': '101,102',
            'mode': '快速模式'
        }
    elif sequence_length < 500:
        # 长序列 - 超快速模式
        return {
            'cycle': 2,
            'step': 75,
            'sample': 1,
            'seeds': '101',
            'mode': '超快速模式'
        }
    else:
        # 超长序列 - 极速模式
        return {
            'cycle': 1,
            'step': 50,
            'sample': 1,
            'seeds': '101',
            'mode': '极速模式'
        }

def run_alphafold_prediction(json_file, output_base_dir="./predictions", cleanup_intermediate=True):
    """运行单个AlphaFold预测"""
    sequence_name = json_file.stem
    output_dir = Path(output_base_dir) / f"{sequence_name}_output"
    
    # 获取序列长度和优化参数
    seq_length = get_sequence_length(json_file)
    params = get_optimal_params(seq_length)
    
    print(f"\n🔬 预测: {json_file.name}")
    print(f"📏 序列长度: {seq_length} AA")
    print(f"⚡ 模式: {params['mode']}")
    
    # 构建命令
    cmd = [
        'protenix', 'predict',
        '--input', str(json_file),
        '--out_dir', str(output_dir),
        '--seeds', params['seeds'],
        '--cycle', str(params['cycle']),
        '--step', str(params['step']),
        '--sample', str(params['sample']),
        '--use_msa_server'
    ]
    
    print(f"🚀 命令: {' '.join(cmd)}")
    
    start_time = time.time()
    
    try:
        # 运行预测
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)  # 30分钟超时
        
        duration = time.time() - start_time
        
        if result.returncode == 0:
            # 检查输出文件
            pdb_files = list(output_dir.glob("**/*.pdb"))
            
            # 清理中间文件
            if cleanup_intermediate and pdb_files:
                cleanup_intermediate_files(output_dir)
            
            print(f"✅ 成功! 耗时: {duration:.1f}秒, PDB文件: {len(pdb_files)}个")
            
            return {
                'status': 'success',
                'duration': duration,
                'pdb_count': len(pdb_files),
                'mode': params['mode'],
                'sequence_length': seq_length,
                'cleaned_up': cleanup_intermediate
            }
        else:
            print(f"❌ 失败! 耗时: {duration:.1f}秒")
            print(f"错误: {result.stderr.strip()[:200]}...")
            
            return {
                'status': 'failed',
                'duration': duration,
                'error': result.stderr.strip(),
                'mode': params['mode'],
                'sequence_length': seq_length
            }
    
    except subprocess.TimeoutExpired:
        print(f"⏰ 超时! (30分钟)")
        return {
            'status': 'timeout',
            'duration': 1800,
            'mode': params['mode'],
            'sequence_length': seq_length
        }
    
    except Exception as e:
        print(f"💥 异常: {e}")
        return {
            'status': 'error',
            'error': str(e),
            'mode': params['mode'],
            'sequence_length': seq_length
        }

def process_batch(input_dir, start_from=0, limit=None, output_dir="./predictions", cleanup_intermediate=True):
    """批处理所有JSON文件"""
    input_path = Path(input_dir)
    json_files = sorted(list(input_path.glob("*.json")))
    
    if not json_files:
        print(f"❌ 在 {input_dir} 中没有找到JSON文件")
        return
    
    print(f"📁 找到 {len(json_files)} 个JSON文件")
    
    # 应用起始位置和限制
    if start_from > 0:
        json_files = json_files[start_from:]
        print(f"⏩ 从第 {start_from + 1} 个文件开始")
    
    if limit:
        json_files = json_files[:limit]
        print(f"🔢 限制处理 {limit} 个文件")
    
    total_files = len(json_files)
    if total_files == 0:
        print("❌ 没有文件需要处理")
        return
    
    # 创建输出目录
    Path(output_dir).mkdir(exist_ok=True)
    
    print(f"\n🚀 开始批处理 {total_files} 个序列...")
    print("=" * 60)
    
    results = []
    success_count = 0
    total_duration = 0
    start_time = time.time()
    
    for i, json_file in enumerate(json_files, 1):
        print(f"\n📊 进度: [{i}/{total_files}] - {json_file.name}")
        
        result = run_alphafold_prediction(json_file, output_dir, cleanup_intermediate)
        result['file'] = json_file.name
        results.append(result)
        
        if result['status'] == 'success':
            success_count += 1
            total_duration += result['duration']
            avg_duration = total_duration / success_count
            
            # 估算剩余时间
            remaining = total_files - i
            eta_minutes = (remaining * avg_duration) / 60
            
            print(f"📈 统计: 成功 {success_count}/{i}, 平均 {avg_duration:.1f}s/个, 预计剩余 {eta_minutes:.1f}分钟")
        
        print("-" * 40)
    
    # 最终统计
    batch_duration = time.time() - start_time
    failed_count = total_files - success_count
    
    print(f"\n🏁 批处理完成!")
    print("=" * 60)
    print(f"📊 总结:")
    print(f"  - 总文件数: {total_files}")
    print(f"  - 成功: {success_count} ({success_count/total_files*100:.1f}%)")
    print(f"  - 失败: {failed_count} ({failed_count/total_files*100:.1f}%)")
    print(f"  - 总耗时: {batch_duration/60:.1f} 分钟")
    if success_count > 0:
        print(f"  - 平均预测时间: {total_duration/success_count:.1f} 秒")
        print(f"  - 预测总时长: {total_duration/60:.1f} 分钟")
        print(f"  - RTX A4000 利用率: {(total_duration/batch_duration)*100:.1f}%")
    
    # 保存结果
    save_results(results, output_dir)
    
    # 显示失败详情
    failed_results = [r for r in results if r['status'] != 'success']
    if failed_results:
        print(f"\n❌ 失败文件详情:")
        for result in failed_results[:5]:
            print(f"  - {result['file']}: {result.get('error', 'Unknown error')[:100]}...")
        if len(failed_results) > 5:
            print(f"  ... 还有 {len(failed_results) - 5} 个失败")

def save_results(results, output_dir):
    """保存结果报告"""
    timestamp = time.strftime('%Y%m%d_%H%M%S')
    report_file = Path(output_dir) / f"rtx4000_batch_results_{timestamp}.json"
    
    # 统计信息
    success_results = [r for r in results if r['status'] == 'success']
    
    report = {
        'metadata': {
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'gpu_model': 'RTX A4000',
            'total_files': len(results),
            'success_count': len(success_results),
            'failed_count': len(results) - len(success_results)
        },
        'statistics': {
            'avg_duration': sum(r.get('duration', 0) for r in success_results) / len(success_results) if success_results else 0,
            'total_prediction_time': sum(r.get('duration', 0) for r in success_results),
            'mode_distribution': {},
            'sequence_length_stats': {
                'min': min(r.get('sequence_length', 0) for r in results if r.get('sequence_length', 0) > 0) if results else 0,
                'max': max(r.get('sequence_length', 0) for r in results),
                'avg': sum(r.get('sequence_length', 0) for r in results) / len(results) if results else 0
            }
        },
        'results': results
    }
    
    # 统计模式分布
    for result in results:
        mode = result.get('mode', 'unknown')
        report['statistics']['mode_distribution'][mode] = report['statistics']['mode_distribution'].get(mode, 0) + 1
    
    with open(report_file, 'w', encoding='utf-8') as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    
    print(f"📄 详细报告已保存: {report_file}")

def main():
    parser = argparse.ArgumentParser(description='RTX A4000 AlphaFold批处理脚本')
    parser.add_argument('--input_dir', '-i', default='/home/webserver/student/students_webserver/zhijing/input_jsons', help='JSON文件输入目录')
    parser.add_argument('--output_dir', '-o', default='./alphafold_predictions', help='输出目录')
    parser.add_argument('--start_from', '-f', type=int, default=0, help='从第几个文件开始')
    parser.add_argument('--limit', '-l', type=int, help='处理文件数量限制')
    parser.add_argument('--no_cleanup', action='store_true', help='不清理中间文件')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input_dir):
        print(f"❌ 输入目录不存在: {args.input_dir}")
        return
    
    print("🎮 RTX A4000 AlphaFold批处理器")
    print("=" * 50)
    print(f"📁 输入目录: {args.input_dir}")
    print(f"📤 输出目录: {args.output_dir}")
    print(f"🔢 处理范围: 从第{args.start_from + 1}个开始" + (f", 限制{args.limit}个" if args.limit else ", 无限制"))
    print(f"🧹 自动清理: {'禁用' if args.no_cleanup else '启用'}")
    print("=" * 50)
    
    # 设置GPU环境
    setup_gpu_environment()
    
    # 开始批处理
    process_batch(
        input_dir=args.input_dir,
        start_from=args.start_from,
        limit=args.limit,
        output_dir=args.output_dir,
        cleanup_intermediate=not args.no_cleanup
    )

if __name__ == "__main__":
    main()