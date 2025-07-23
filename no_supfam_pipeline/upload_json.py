#!/usr/bin/env python3
"""
批量上传脚本 - 上传DALI生成的所有JSON文件到远程服务器
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path
import glob

class BatchFileUploader:
    def __init__(self):
        # 固定路径 - DALI生成的JSON文件夹
        self.json_folder = "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/dali_sequences"
        
        # 远程服务器配置
        self.remote_server = "webserver@coulomb.phys.ucl.ac.uk"
        self.remote_input_dir = "~/student/students_webserver/zhijing/input_jsons"
    
    def get_json_files(self):
        """获取所有JSON文件"""
        if not os.path.exists(self.json_folder):
            print(f"❌ JSON文件夹不存在: {self.json_folder}")
            return []
        
        json_files = glob.glob(os.path.join(self.json_folder, "*.json"))
        
        if not json_files:
            print(f"❌ 在文件夹中没有找到JSON文件: {self.json_folder}")
            return []
        
        print(f"📁 找到 {len(json_files)} 个JSON文件")
        return json_files
    
    def verify_json_file(self, json_file):
        """验证单个JSON文件"""
        if not os.path.exists(json_file):
            print(f"❌ 文件不存在: {json_file}")
            return False
        
        if os.path.getsize(json_file) == 0:
            print(f"❌ 文件为空: {json_file}")
            return False
        
        return True
    
    def upload_single_file(self, json_file):
        """上传单个JSON文件"""
        filename = os.path.basename(json_file)
        
        print(f"📤 上传: {filename}")
        
        cmd = f"scp '{json_file}' {self.remote_server}:{self.remote_input_dir}/"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"❌ 上传失败 {filename}: {result.stderr}")
            return False
        
        print(f"✅ 上传成功: {filename}")
        return True
    
    def upload_all_files(self):
        """批量上传所有JSON文件"""
        print("🚀 开始批量上传DALI JSON文件...")
        print(f"📂 源文件夹: {self.json_folder}")
        print(f"🖥️  目标服务器: {self.remote_server}")
        print(f"📁 目标路径: {self.remote_input_dir}")
        print("="*60)
        
        # 获取所有JSON文件
        json_files = self.get_json_files()
        if not json_files:
            return False
        
        # 统计变量
        total_files = len(json_files)
        successful_uploads = 0
        failed_uploads = 0
        
        # 逐个上传文件
        for i, json_file in enumerate(json_files, 1):
            filename = os.path.basename(json_file)
            print(f"[{i}/{total_files}] 处理: {filename}")
            
            # 验证文件
            if not self.verify_json_file(json_file):
                failed_uploads += 1
                continue
            
            # 上传文件
            if self.upload_single_file(json_file):
                successful_uploads += 1
            else:
                failed_uploads += 1
            
            print()  # 空行分隔
        
        # 输出统计结果
        print("="*60)
        print("📊 上传统计:")
        print(f"   总文件数: {total_files}")
        print(f"   成功上传: {successful_uploads}")
        print(f"   失败数量: {failed_uploads}")
        print("="*60)
        
        if successful_uploads > 0:
            print("🎉 批量上传完成!")
            print(f"📝 后续步骤:")
            print(f"   1. SSH登录远程服务器: ssh {self.remote_server}")
            print(f"   2. 切换到工作目录: cd ~/student/students_webserver/zhijing")
            print(f"   3. 批量运行远程脚本 (示例):")
            
            # 显示几个示例命令
            example_files = [os.path.basename(f) for f in json_files[:3]]
            for filename in example_files:
                print(f"      python remote_pipeline.py {filename}")
            
            if len(json_files) > 3:
                print(f"      ... (还有 {len(json_files)-3} 个文件)")
        
        return successful_uploads > 0

class SingleFileUploader:
    def __init__(self, json_file_path):
        self.json_file_path = json_file_path
        self.json_filename = os.path.basename(json_file_path)
        self.output_name = os.path.splitext(self.json_filename)[0]
        self.remote_server = "webserver@coulomb.phys.ucl.ac.uk"
        self.remote_input_dir = "~/student/students_webserver/zhijing/input_jsons"
    
    def verify_json_file(self):
        """验证JSON文件是否存在"""
        if not os.path.exists(self.json_file_path):
            print(f"❌ JSON文件不存在: {self.json_file_path}")
            return False
        
        print(f"✅ JSON文件验证成功: {self.json_file_path}")
        return True
    
    def upload_json_to_remote(self):
        """上传JSON文件到远程服务器"""
        print(f"📤 上传JSON文件到远程服务器...")
        
        cmd = f"scp '{self.json_file_path}' {self.remote_server}:{self.remote_input_dir}/"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"❌ 上传失败: {result.stderr}")
            return False
        
        print(f"✅ 文件上传成功: {self.json_filename}")
        return True
    
    def upload_file(self):
        """上传单个文件的主函数"""
        print(f"🚀 开始文件上传...")
        print(f"📋 JSON文件: {self.json_file_path}")
        print(f"📋 文件名: {self.json_filename}")
        print(f"📋 输出名: {self.output_name}")
        print("="*50)
        
        # 验证文件
        if not self.verify_json_file():
            return False
        
        # 上传文件
        if not self.upload_json_to_remote():
            return False
        
        print("="*50)
        print("🎉 文件上传完成!")
        print(f"📝 后续步骤:")
        print(f"   1. SSH登录远程服务器: ssh {self.remote_server}")
        print(f"   2. 切换到工作目录: cd ~/student/students_webserver/zhijing")
        print(f"   3. 运行远程脚本: python remote_pipeline.py {self.json_filename}")
        print("="*50)
        
        return True

def main():
    parser = argparse.ArgumentParser(
        description='上传JSON文件到远程服务器',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  python upload_json.py --batch                    # 批量上传所有DALI JSON文件
  python upload_json.py single_file.json          # 上传单个文件
  python upload_json.py --file single_file.json   # 上传单个文件(参数形式)
        """
    )
    
    parser.add_argument('json_file', nargs='?', help='单个JSON文件路径')
    parser.add_argument('--file', help='单个JSON文件路径(参数形式)')
    parser.add_argument('--batch', action='store_true', 
                       help='批量上传所有DALI生成的JSON文件')
    
    args = parser.parse_args()
    
    if args.batch:
        # 批量上传模式
        uploader = BatchFileUploader()
        success = uploader.upload_all_files()
    else:
        # 单文件上传模式
        json_file = args.json_file or args.file
        if not json_file:
            print("❌ 错误: 请指定JSON文件路径或使用 --batch 进行批量上传")
            print("使用方法: python upload_json.py --batch")
            print("或者:     python upload_json.py your_file.json")
            sys.exit(1)
        
        uploader = SingleFileUploader(json_file)
        success = uploader.upload_file()
    
    if success:
        sys.exit(0)
    else:
        print("💥 上传失败!")
        sys.exit(1)

if __name__ == "__main__":
    main()