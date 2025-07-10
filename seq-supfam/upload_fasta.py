#!/usr/bin/env python3
"""
交互式FASTA上传器
上传文件后，打开交互式SSH会话让你手动运行分析
"""

import os
import sys
import subprocess

class InteractiveUploader:
    def __init__(self, fasta_file_path):
        self.fasta_file_path = fasta_file_path
        self.fasta_filename = os.path.basename(fasta_file_path)
        self.sequence_id = os.path.splitext(self.fasta_filename)[0]
        
        # 固定配置
        self.remote_server = "webserver@coulomb.phys.ucl.ac.uk"
        self.upload_dir = "~/student/students_webserver/zhijing/fasta_sequence"
        self.superfamily_dir = "/mnt/data2/supfam/supfam"
        
    def verify_and_upload(self):
        """验证并上传文件"""
        if not os.path.exists(self.fasta_file_path):
            print(f"❌ FASTA文件不存在: {self.fasta_file_path}")
            return False
        
        if not self.fasta_filename.endswith('.fa'):
            print(f"❌ 文件必须以.fa结尾: {self.fasta_filename}")
            return False
        
        print(f"✅ FASTA文件验证成功: {self.fasta_file_path}")
        
        # 上传文件
        print(f"📤 上传FASTA文件到远程服务器...")
        cmd = f"scp {self.fasta_file_path} {self.remote_server}:{self.upload_dir}/"
        result = subprocess.run(cmd, shell=True)
        
        if result.returncode != 0:
            print(f"❌ 上传失败")
            return False
        
        print(f"✅ 文件上传成功: {self.fasta_filename}")
        return True
    
    def copy_file_remote(self):
        """在远程服务器上复制文件到SUPERFAMILY目录"""
        print(f"📋 复制文件到SUPERFAMILY目录...")
        
        copy_cmd = f"cp {self.upload_dir}/{self.fasta_filename} {self.superfamily_dir}/"
        ssh_cmd = f"ssh {self.remote_server} '{copy_cmd}'"
        
        result = subprocess.run(ssh_cmd, shell=True)
        
        if result.returncode != 0:
            print(f"❌ 复制文件失败")
            return False
        
        print(f"✅ 文件已复制到SUPERFAMILY目录")
        return True
    
    def start_interactive_session(self):
        """启动交互式SSH会话"""
        print(f"\n🚀 启动交互式SSH会话...")
        print("="*60)
        print(f"📋 文件已准备好: {self.superfamily_dir}/{self.fasta_filename}")
        print(f"📋 序列ID: {self.sequence_id}")
        print("="*60)
        print("🔄 现在将打开SSH会话，请在远程服务器上运行:")
        print(f"   cd {self.superfamily_dir}")
        print(f"   ./superfamily.pl {self.fasta_filename}")
        print()
        print("⚠️  分析完成后，输入 'exit' 退出SSH会话")
        print("   脚本将自动下载结果文件")
        print("="*60)
        
        input("按回车键开始SSH会话...")
        
        # 启动交互式SSH，直接进入SUPERFAMILY目录
        ssh_cmd = f"ssh -t {self.remote_server} 'cd {self.superfamily_dir}; bash'"
        
        print(f"\n🔗 正在连接到远程服务器...")
        subprocess.run(ssh_cmd, shell=True)
        
        print(f"\n✅ SSH会话已结束")
        return True
    
    def download_results(self):
        """下载HTML结果文件"""
        print(f"\n📥 开始下载HTML结果文件...")
        
        # 创建HTML结果目录
        html_dir = "html_results"
        if not os.path.exists(html_dir):
            os.makedirs(html_dir)
            print(f"📁 创建HTML结果目录: {html_dir}")
        
        # 只下载HTML文件
        html_files_to_try = [
            f"{self.sequence_id}.html",
            f"{self.sequence_id}_html"
        ]
        
        downloaded_files = []
        
        for file_name in html_files_to_try:
            remote_path = f"{self.remote_server}:{self.superfamily_dir}/{file_name}"
            local_path = os.path.join(html_dir, f"{self.sequence_id}.html")  # 统一命名为 .html
            
            print(f"   尝试下载: {file_name}")
            cmd = f"scp {remote_path} {local_path}"
            result = subprocess.run(cmd, shell=True, capture_output=True)
            
            if result.returncode == 0:
                downloaded_files.append(local_path)
                print(f"   ✅ HTML文件下载成功")
                print(f"      保存位置: {local_path}")
                
                # 检查文件大小
                if os.path.exists(local_path):
                    file_size = os.path.getsize(local_path)
                    print(f"      文件大小: {file_size} 字节")
                    
                    if file_size > 500:
                        print(f"      ✅ HTML文件看起来完整")
                    else:
                        print(f"      ⚠️  HTML文件可能不完整")
                
                break  # 成功下载一个就停止尝试
            else:
                print(f"   ❌ {file_name} 下载失败 (可能不存在)")
        
        if not downloaded_files:
            print(f"   ⚠️  没有找到HTML结果文件")
        
        return downloaded_files
    
    def process_file(self):
        """处理文件的主函数"""
        print(f"🧬 交互式FASTA处理器")
        print(f"📋 文件: {self.fasta_filename}")
        print(f"📋 序列ID: {self.sequence_id}")
        print("="*50)
        
        # 1. 验证并上传文件
        if not self.verify_and_upload():
            return False
        
        # 2. 复制文件到SUPERFAMILY目录
        if not self.copy_file_remote():
            return False
        
        # 3. 启动交互式会话
        if not self.start_interactive_session():
            return False
        
        # 4. 下载结果
        downloaded_files = self.download_results()
        
        # 5. 总结
        print(f"\n{'='*60}")
        print("🎉 处理完成!")
        
        if downloaded_files:
            print(f"📋 HTML结果文件:")
            for file_path in downloaded_files:
                if os.path.exists(file_path):
                    file_size = os.path.getsize(file_path)
                    print(f"   📄 {file_path} ({file_size} 字节)")
                    print(f"   🌐 可以在浏览器中打开查看结果")
        else:
            print("⚠️  没有下载到HTML结果文件")
            print("💡 请检查分析是否成功完成")
        
        print("="*60)
        return True

def process_all_files():
    """批量处理所有文件"""
    fasta_dir = "fasta_files"
    
    if not os.path.exists(fasta_dir):
        print(f"❌ 目录不存在: {fasta_dir}")
        return False
    
    # 查找所有.fa文件
    fasta_files = []
    for filename in os.listdir(fasta_dir):
        if filename.endswith('.fa'):
            fasta_path = os.path.join(fasta_dir, filename)
            fasta_files.append(fasta_path)
    
    if not fasta_files:
        print(f"❌ 在 {fasta_dir} 目录中没有找到.fa文件")
        return False
    
    print(f"🧬 批量交互式处理器")
    print("="*50)
    print(f"📁 找到 {len(fasta_files)} 个FASTA文件:")
    for fasta_path in fasta_files:
        filename = os.path.basename(fasta_path)
        print(f"   - {filename}")
    
    # 确认处理
    confirm = input(f"\n逐个处理这些文件? (y/N): ").strip().lower()
    if confirm not in ['y', 'yes']:
        print("👋 操作已取消")
        return False
    
    # 逐个处理
    for i, fasta_path in enumerate(fasta_files, 1):
        print(f"\n{'='*60}")
        print(f"处理第 {i}/{len(fasta_files)} 个文件")
        print(f"{'='*60}")
        
        uploader = InteractiveUploader(fasta_path)
        uploader.process_file()
        
        if i < len(fasta_files):
            input(f"\n按回车键继续处理下一个文件...")
    
    print(f"\n🎉 所有文件处理完成!")
    return True

def main():
    if len(sys.argv) == 1:
        # 批量模式
        process_all_files()
    elif len(sys.argv) == 2:
        # 单文件模式
        fasta_file = sys.argv[1]
        uploader = InteractiveUploader(fasta_file)
        uploader.process_file()
    else:
        print("使用方法:")
        print("  python interactive_uploader.py                # 批量处理所有文件")
        print("  python interactive_uploader.py <file.fa>      # 处理单个文件")
        print()
        print("示例:")
        print("  python interactive_uploader.py fasta_files/test.fa")
        print("  python interactive_uploader.py")

if __name__ == "__main__":
    main()