#!/usr/bin/env python3
"""
本地上传脚本 - 只负责上传JSON文件到远程服务器
"""

import os
import sys
import subprocess
import argparse

class FileUploader:
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
        
        cmd = f"scp {self.json_file_path} {self.remote_server}:{self.remote_input_dir}/"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"❌ 上传失败: {result.stderr}")
            return False
        
        print(f"✅ 文件上传成功: {self.json_filename}")
        return True
    
    def upload_file(self):
        """上传文件的主函数"""
        print(f"🚀 开始上传文件...")
        print(f"📋 JSON文件: {self.json_file_path}")
        print(f"📋 文件名: {self.json_filename}")
        print(f"📋 输出名称: {self.output_name}")
        print("="*50)
        
        # 验证文件
        if not self.verify_json_file():
            return False
        
        # 上传文件
        if not self.upload_json_to_remote():
            return False
        
        print("="*50)
        print("🎉 文件上传完成!")
        print(f"📝 接下来请:")
        print(f"   1. SSH登录到远程服务器: ssh {self.remote_server}")
        print(f"   2. 切换到工作目录: cd ~/student/students_webserver/zhijing")
        print(f"   3. 运行远程脚本: python remote_pipeline.py {self.json_filename}")
        print("="*50)
        
        return True

def main():
    parser = argparse.ArgumentParser(description='上传JSON文件到远程服务器')
    parser.add_argument('json_file', help='JSON文件路径')
    
    args = parser.parse_args()
    
    uploader = FileUploader(args.json_file)
    success = uploader.upload_file()
    
    if success:
        sys.exit(0)
    else:
        print("💥 上传失败!")
        sys.exit(1)

if __name__ == "__main__":
    main()