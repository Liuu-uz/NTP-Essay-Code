#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FASTA文件生成器
功能：将蛋白质序列转换为标准FASTA格式文件
"""

import os
import re
from datetime import datetime

def create_fasta_file(sequence, sequence_id=None, description="", output_dir="fasta_files"):
    """
    创建FASTA格式文件
    
    参数:
    - sequence: 蛋白质序列
    - sequence_id: 序列ID（可选）
    - description: 序列描述（可选）
    - output_dir: 输出目录
    
    返回:
    - fasta_path: 生成的FASTA文件路径
    - clean_seq_id: 清理后的序列ID
    """
    
    print("📄 开始生成FASTA文件...")
    
    # 1. 清理序列（移除空格、换行符、非氨基酸字符）
    print("🧹 清理序列...")
    clean_sequence = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', sequence.upper())
    
    if len(clean_sequence) == 0:
        raise ValueError("❌ 序列为空或不包含有效的氨基酸")
    
    print(f"   原始序列长度: {len(sequence)} 字符")
    print(f"   清理后序列长度: {len(clean_sequence)} 氨基酸")
    
    # 2. 生成序列ID
    if sequence_id is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        sequence_id = f"seq_{timestamp}"
        print(f"🏷️  自动生成序列ID: {sequence_id}")
    else:
        print(f"🏷️  使用提供的序列ID: {sequence_id}")
    
    # 3. 清理序列ID（只保留字母数字和下划线，符合superfamily.pl要求）
    clean_seq_id = re.sub(r'[^a-zA-Z0-9_]', '_', sequence_id)
    if clean_seq_id != sequence_id:
        print(f"🔧 序列ID已清理: {sequence_id} → {clean_seq_id}")
    
    # 4. 创建输出目录
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"📁 创建输出目录: {output_dir}")
    
    # 5. 创建FASTA内容
    fasta_content = f">{clean_seq_id}"
    if description:
        fasta_content += f" {description}"
    fasta_content += f"\n{clean_sequence}\n"
    
    # 6. 保存到文件（使用.fa扩展名，符合superfamily.pl要求）
    fasta_filename = f"{clean_seq_id}.fa"
    fasta_path = os.path.join(output_dir, fasta_filename)
    
    with open(fasta_path, 'w') as f:
        f.write(fasta_content)
    
    print(f"✅ FASTA文件已创建: {fasta_path}")
    print(f"   序列ID: {clean_seq_id}")
    print(f"   序列长度: {len(clean_sequence)} 氨基酸")
    print(f"   文件大小: {os.path.getsize(fasta_path)} 字节")
    
    return fasta_path, clean_seq_id

def batch_create_fasta(sequences_list, output_dir="fasta_files"):
    """
    批量创建FASTA文件
    
    参数:
    - sequences_list: 序列信息列表，每个元素包含 sequence, sequence_id, description
    - output_dir: 输出目录
    
    返回:
    - created_files: 创建的文件列表
    """
    
    print(f"🚀 开始批量生成 {len(sequences_list)} 个FASTA文件...")
    
    created_files = []
    
    for i, seq_info in enumerate(sequences_list, 1):
        print(f"\n📄 [{i}/{len(sequences_list)}] 处理序列: {seq_info.get('sequence_id', 'unknown')}")
        
        try:
            fasta_path, clean_seq_id = create_fasta_file(
                sequence=seq_info['sequence'],
                sequence_id=seq_info.get('sequence_id'),
                description=seq_info.get('description', ''),
                output_dir=output_dir
            )
            
            created_files.append({
                'sequence_id': clean_seq_id,
                'fasta_path': fasta_path,
                'original_info': seq_info
            })
            
        except Exception as e:
            print(f"❌ 处理序列失败: {e}")
            continue
    
    print(f"\n✅ 批量生成完成! 成功创建 {len(created_files)} 个文件")
    return created_files

def validate_fasta_file(fasta_path):
    """
    验证FASTA文件格式
    """
    print(f"🔍 验证FASTA文件: {fasta_path}")
    
    if not os.path.exists(fasta_path):
        print("❌ 文件不存在")
        return False
    
    try:
        with open(fasta_path, 'r') as f:
            content = f.read()
        
        lines = content.strip().split('\n')
        
        # 检查第一行是否以>开头
        if not lines[0].startswith('>'):
            print("❌ 第一行不是有效的FASTA标题行")
            return False
        
        # 检查序列行
        sequence_lines = lines[1:]
        sequence = ''.join(sequence_lines)
        
        # 检查是否包含有效氨基酸
        valid_aa = re.match(r'^[ACDEFGHIKLMNPQRSTVWY]+$', sequence)
        if not valid_aa:
            print("❌ 序列包含无效字符")
            return False
        
        print(f"✅ FASTA文件验证通过")
        print(f"   标题: {lines[0]}")
        print(f"   序列长度: {len(sequence)} 氨基酸")
        return True
        
    except Exception as e:
        print(f"❌ 验证失败: {e}")
        return False

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FASTA文件生成器 - 自动运行版本
在下面的SEQUENCES_TO_PROCESS列表中添加你的序列，直接运行即可
"""

import os
import re
from datetime import datetime

# ==================== 序列配置区域 ====================
# 在这里添加你要生成FASTA的序列
SEQUENCES_TO_PROCESS = [
    {
        "sequence": """
        MMKKIDVKILDPRVGKEFPLPTYATSGSAGLDLRACLNDAVELAPGDTTLVPTGLAIHIA
        DPSLAAMMLPRSGLGHKHGIVLGNLVGLIDSDYQGQLMISVWNRGQDSFTIQPGERIAQM
        IFVPVVQAEFNLVEDFDFRQ
        """,
        "sequence_id": "test",
        "description": "test"
    },
    
    # 在这里添加你的序列 - 复制上面的格式
    {
        "sequence": """
        把你的蛋白质序列粘贴到这里
        可以包含换行符和空格，程序会自动清理
        """,
        "sequence_id": "my_protein_001",
        "description": "我的测试蛋白质"
    }
    
    # 继续添加更多序列...
]

# 输出目录设置
OUTPUT_DIR = "fasta_files"

# ==================== 处理函数 ====================

def create_fasta_file(sequence, sequence_id=None, description="", output_dir=OUTPUT_DIR):
    """
    创建FASTA格式文件
    """
    
    # 1. 清理序列（移除空格、换行符、非氨基酸字符）
    clean_sequence = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', sequence.upper())
    
    if len(clean_sequence) == 0:
        raise ValueError("❌ 序列为空或不包含有效的氨基酸")
    
    # 2. 生成序列ID
    if sequence_id is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        sequence_id = f"seq_{timestamp}"
    
    # 3. 清理序列ID（只保留字母数字和下划线，符合superfamily.pl要求）
    clean_seq_id = re.sub(r'[^a-zA-Z0-9_]', '_', sequence_id)
    
    # 4. 创建输出目录
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 5. 创建FASTA内容
    fasta_content = f">{clean_seq_id}"
    if description:
        fasta_content += f" {description}"
    fasta_content += f"\n{clean_sequence}\n"
    
    # 6. 保存到文件（使用.fa扩展名，符合superfamily.pl要求）
    fasta_filename = f"{clean_seq_id}.fa"
    fasta_path = os.path.join(output_dir, fasta_filename)
    
    with open(fasta_path, 'w') as f:
        f.write(fasta_content)
    
    return fasta_path, clean_seq_id, len(clean_sequence)

def auto_generate_fasta_files():
    """
    自动生成所有配置的FASTA文件
    """
    print("🧬 FASTA文件自动生成器")
    print("=" * 50)
    
    # 过滤有效序列
    valid_sequences = []
    for seq_info in SEQUENCES_TO_PROCESS:
        if seq_info['sequence'].strip() and "把你的" not in seq_info['sequence']:
            valid_sequences.append(seq_info)
        else:
            print(f"⚠️  跳过空序列或模板: {seq_info.get('sequence_id', 'unknown')}")
    
    if not valid_sequences:
        print("❌ 没有找到有效的序列")
        print("💡 请在脚本顶部的SEQUENCES_TO_PROCESS中添加你的序列")
        return []
    
    print(f"📄 开始生成 {len(valid_sequences)} 个FASTA文件...")
    
    created_files = []
    
    for i, seq_info in enumerate(valid_sequences, 1):
        print(f"\n[{i}/{len(valid_sequences)}] 处理序列: {seq_info['sequence_id']}")
        
        try:
            fasta_path, clean_seq_id, seq_length = create_fasta_file(
                sequence=seq_info['sequence'],
                sequence_id=seq_info['sequence_id'],
                description=seq_info.get('description', '')
            )
            
            print(f"✅ 生成成功: {fasta_path}")
            print(f"   序列ID: {clean_seq_id}")
            print(f"   序列长度: {seq_length} 氨基酸")
            
            created_files.append({
                'sequence_id': clean_seq_id,
                'fasta_path': fasta_path,
                'length': seq_length
            })
            
        except Exception as e:
            print(f"❌ 生成失败: {e}")
            continue
    
    # 总结报告
    print(f"\n{'='*50}")
    print("📊 生成完成统计")
    print(f"{'='*50}")
    print(f"配置序列数: {len(SEQUENCES_TO_PROCESS)}")
    print(f"有效序列数: {len(valid_sequences)}")
    print(f"成功生成数: {len(created_files)}")
    print(f"输出目录: {OUTPUT_DIR}")
    
    if created_files:
        print(f"\n📋 生成的文件列表:")
        for file_info in created_files:
            print(f"  - {file_info['fasta_path']} ({file_info['length']} aa)")
        
        print(f"\n🚀 下一步:")
        print(f"   使用 remote_superfamily_runner.py 上传并分析这些文件")
    
    return created_files

def validate_all_files():
    """
    验证生成的所有FASTA文件
    """
    print(f"\n🔍 验证 {OUTPUT_DIR} 目录中的FASTA文件...")
    
    if not os.path.exists(OUTPUT_DIR):
        print(f"❌ 目录不存在: {OUTPUT_DIR}")
        return
    
    fasta_files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith('.fa')]
    
    if not fasta_files:
        print(f"❌ 目录中没有找到.fa文件")
        return
    
    print(f"找到 {len(fasta_files)} 个FASTA文件:")
    
    for fasta_file in fasta_files:
        fasta_path = os.path.join(OUTPUT_DIR, fasta_file)
        
        try:
            with open(fasta_path, 'r') as f:
                content = f.read().strip()
            
            lines = content.split('\n')
            header = lines[0] if lines else ""
            sequence = ''.join(lines[1:]) if len(lines) > 1 else ""
            
            if header.startswith('>') and sequence:
                valid_aa = re.match(r'^[ACDEFGHIKLMNPQRSTVWY]+$', sequence)
                if valid_aa:
                    print(f"  ✅ {fasta_file} - {len(sequence)} aa")
                else:
                    print(f"  ⚠️  {fasta_file} - 包含无效字符")
            else:
                print(f"  ❌ {fasta_file} - 格式错误")
                
        except Exception as e:
            print(f"  ❌ {fasta_file} - 读取失败: {e}")

def list_generated_files():
    """
    列出已生成的文件，供其他脚本调用
    """
    if not os.path.exists(OUTPUT_DIR):
        return []
    
    fasta_files = []
    for filename in os.listdir(OUTPUT_DIR):
        if filename.endswith('.fa'):
            fasta_path = os.path.join(OUTPUT_DIR, filename)
            sequence_id = filename.replace('.fa', '')
            fasta_files.append({
                'sequence_id': sequence_id,
                'fasta_path': fasta_path,
                'filename': filename
            })
    
    return fasta_files

def main():
    """
    主函数 - 自动运行
    """
    # 自动生成FASTA文件
    created_files = auto_generate_fasta_files()
    
    # 验证生成的文件
    if created_files:
        validate_all_files()
    
    return created_files

if __name__ == "__main__":
    # 直接运行，不需要用户交互
    main()