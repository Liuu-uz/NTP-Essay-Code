#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HTML与NTP超家族比对器
比对SUPERFAMILY生成的HTML文件与Excel表格中的NTP-processing超家族
"""

import os
import sys
import pandas as pd
from bs4 import BeautifulSoup

class NTPChecker:
    def __init__(self, excel_path=None):
        self.excel_path = excel_path
        self.ntp_sf_dict = {}
        self.sf_names = []
        self.find_excel_file()
        self.load_ntp_data()
    
    def find_excel_file(self):
        """查找Excel文件"""
        if self.excel_path and os.path.exists(self.excel_path):
            return
        
        # 可能的Excel文件位置
        possible_paths = [
            "SUPFAM based NTP processing.xlsx",  # 当前目录
            "../SUPFAM based NTP processing.xlsx",  # 上级目录
            "/Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/seq-supfam/SUPFAM based NTP processing.xlsx",  # 指定目录
            os.path.expanduser("~/NTP-Essay-Code-1/NTP-Essay-Code/seq-supfam/SUPFAM based NTP processing.xlsx"),  # 用户目录
        ]
        
        for path in possible_paths:
            if os.path.exists(path):
                self.excel_path = path
                print(f"✅ 找到Excel文件: {path}")
                return
        
        print("❌ 找不到Excel文件，请检查以下位置:")
        for path in possible_paths:
            print(f"   {path}")
        self.excel_path = None
    
    def load_ntp_data(self):
        """加载NTP-processing超家族数据"""
        if not self.excel_path or not os.path.exists(self.excel_path):
            print(f"❌ Excel文件不存在或未找到")
            print("💡 请确保文件存在于以下位置之一:")
            print("   - 当前目录: SUPFAM based NTP processing.xlsx")
            print("   - 指定目录: /Users/napkin/NTP-Essay-Code-1/NTP-Essay-Code/seq-supfam/")
            return False
        
        try:
            df = pd.read_excel(self.excel_path)
            
            # 存储完整的行信息，包括Representative (PDB ID)
            self.ntp_data = {}  # 存储完整信息
            
            for i, row in df.iterrows():
                sf_id = row['SF ID'] if pd.notna(row['SF ID']) else None
                sf_name = row['SF'] if pd.notna(row['SF']) else None
                representative = row['Representative'] if pd.notna(row['Representative']) else None
                reaction = row['Reaction'] if pd.notna(row['Reaction']) else None
                coordination = row['Coordination'] if pd.notna(row['Coordination']) else None
                
                # 创建条目信息
                entry_info = {
                    'sf_name': sf_name,
                    'sf_id': sf_id,
                    'representative': representative,
                    'reaction': reaction,
                    'coordination': coordination
                }
                
                if sf_id and sf_id != 'N/A':
                    self.ntp_sf_dict[sf_id] = sf_name
                    self.ntp_data[sf_id] = entry_info
                
                if sf_name:
                    self.sf_names.append(sf_name)
                    if sf_id and sf_id != 'N/A':
                        self.ntp_sf_dict[sf_name] = sf_id
                        self.ntp_data[sf_name] = entry_info
                    else:
                        self.ntp_sf_dict[sf_name] = 'N/A'
                        self.ntp_data[sf_name] = entry_info
            
            sf_with_id = len([k for k in self.ntp_sf_dict.keys() if '.' in str(k)])
            print(f"✅ 已加载 {len(df)} 个NTP-processing超家族")
            print(f"   其中 {sf_with_id} 个有SF ID, {len(df) - sf_with_id} 个仅有名称")
            return True
            
        except Exception as e:
            print(f"❌ 读取Excel文件失败: {e}")
            return False
    
    def parse_html_file(self, html_path):
        """解析SUPERFAMILY HTML结果文件"""
        if not os.path.exists(html_path):
            print(f"❌ HTML文件不存在: {html_path}")
            return []
        
        try:
            # 读取文件内容
            with open(html_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            print(f"   调试: 文件大小: {len(content)} 字符")
            
            # 修复HTML格式问题：将自闭合的table标签改为正常格式
            # 原始: <table border="1" />
            # 修复后: <table border="1">
            content = content.replace('<table border="1" />', '<table border="1">')
            
            # 确保table标签正确闭合
            if '</table>' not in content and '<table' in content:
                # 在body结束前添加</table>
                content = content.replace('</body>', '</table></body>')
            
            print(f"   调试: 修复后的HTML片段:")
            # 找到表格部分并显示
            table_start = content.find('<table')
            table_end = content.find('</table>') + 8
            if table_start >= 0 and table_end > table_start:
                table_html = content[table_start:table_end]
                print(f"   {table_html[:300]}...")
            
            soup = BeautifulSoup(content, 'html.parser')
            found_superfamilies = []
            
            # 查找表格
            tables = soup.find_all('table')
            print(f"   调试: 修复后找到 {len(tables)} 个表格")
            
            for table_idx, table in enumerate(tables):
                print(f"   调试: 处理第 {table_idx + 1} 个表格")
                
                rows = table.find_all('tr')
                print(f"   调试: 找到 {len(rows)} 个tr标签")
                
                # 处理每一行
                for row_idx, row in enumerate(rows):
                    cells = row.find_all(['td', 'th'])
                    print(f"   调试: 第 {row_idx + 1} 行有 {len(cells)} 个单元格")
                    
                    # 显示行内容
                    if len(cells) > 0:
                        print(f"   调试: 第 {row_idx + 1} 行内容:")
                        for i, cell in enumerate(cells):
                            cell_text = cell.get_text(strip=True)
                            print(f"     列{i+1}: '{cell_text}'")
                            # 显示链接
                            link = cell.find('a')
                            if link:
                                link_text = link.get_text(strip=True)
                                print(f"       链接: '{link_text}'")
                    
                    # 跳过标题行
                    if row_idx == 0:
                        print("   调试: 这是标题行，跳过")
                        continue
                    
                    # 处理数据行
                    if len(cells) >= 4:
                        try:
                            seq_id = cells[0].get_text(strip=True)
                            
                            # 从第4列提取SCOP superfamily
                            scop_superfamily_cell = cells[3]
                            sf_link = scop_superfamily_cell.find('a')
                            
                            if sf_link:
                                sf_name = sf_link.get_text(strip=True)
                                print(f"   调试: ✅ 从第4列链接提取超家族: '{sf_name}'")
                            else:
                                sf_name = scop_superfamily_cell.get_text(strip=True)
                                print(f"   调试: 从第4列文本提取超家族: '{sf_name}'")
                            
                            # 确定SCOP ID
                            scop_id = 'unknown'
                            if sf_name == 'HIT-like':
                                scop_id = 'd.13.1'  # 已知的SCOP ID
                            
                            if sf_name and len(sf_name) > 2 and seq_id:
                                print(f"   调试: ✅ 添加超家族: 名称='{sf_name}', SCOP_ID='{scop_id}', 序列='{seq_id}'")
                                found_superfamilies.append({
                                    'name': sf_name,
                                    'scop_id': scop_id,
                                    'seq_id': seq_id
                                })
                            else:
                                print(f"   调试: ❌ 跳过无效条目")
                                
                        except Exception as e:
                            print(f"   调试: 处理行时出错: {e}")
                            continue
            
            print(f"   调试: 最终找到 {len(found_superfamilies)} 个超家族")
            return found_superfamilies
            
        except Exception as e:
            print(f"❌ 解析HTML文件失败: {e}")
            import traceback
            traceback.print_exc()
            return []
    
    def check_ntp_matches(self, found_superfamilies):
        """检查找到的超家族是否与NTP-processing匹配"""
        if not found_superfamilies:
            return []
        
        matches = []
        
        for sf in found_superfamilies:
            sf_name = sf['name']
            scop_id = sf['scop_id']
            seq_id = sf.get('seq_id', 'unknown')
            
            matched_entry = None
            match_type = None
            
            # 1. SCOP ID精确匹配（主要匹配方式）
            if scop_id and scop_id in self.ntp_data:
                matched_entry = self.ntp_data[scop_id]
                match_type = 'SCOP_ID_exact'
            
            # 2. 超家族名称精确匹配
            elif sf_name and sf_name in self.ntp_data:
                matched_entry = self.ntp_data[sf_name]
                match_type = 'SF_name_exact'
            
            # 3. 通过SCOP ID在字典中查找对应的超家族名称
            elif scop_id:
                # 查找Excel中是否有这个SCOP ID对应的超家族
                for key, entry in self.ntp_data.items():
                    if entry['sf_id'] == scop_id:
                        matched_entry = entry
                        match_type = 'SCOP_ID_lookup'
                        break
            
            # 4. 智能模糊匹配（仅作为最后手段）
            if not matched_entry and sf_name:
                for ntp_name in self.sf_names:
                    if len(ntp_name) > 15:  # 只匹配较长的名称
                        sf_lower = sf_name.lower()
                        ntp_lower = ntp_name.lower()
                        
                        # 更严格的匹配条件
                        if (ntp_lower == sf_lower or  # 完全匹配
                            (len(ntp_lower) > 10 and ntp_lower in sf_lower) or  # 长子串匹配
                            (len(sf_lower) > 10 and sf_lower in ntp_lower)):
                            
                            matched_entry = self.ntp_data[ntp_name]
                            match_type = 'fuzzy'
                            break
            
            if matched_entry and match_type:
                matches.append({
                    'found_name': sf_name,
                    'found_scop': scop_id,
                    'found_seq_id': seq_id,
                    'matched_entry': matched_entry,
                    'match_type': match_type
                })
        
        return matches
    
    def check_single_html(self, html_path):
        """检查单个HTML文件"""
        print(f"🔍 检查文件: {os.path.basename(html_path)}")
        print("-" * 50)
        
        # 解析HTML文件
        found_superfamilies = self.parse_html_file(html_path)
        
        if not found_superfamilies:
            print("⚠️  未在HTML中找到任何超家族信息")
            print("❌ NO NTP-processing family found")
            return False
        
        print(f"📊 在HTML中找到 {len(found_superfamilies)} 个超家族:")
        for i, sf in enumerate(found_superfamilies, 1):
            print(f"   {i}. {sf['name']}")
            if sf['scop_id'] and sf['scop_id'] != 'unknown':
                print(f"      SCOP ID: {sf['scop_id']}")
            if sf.get('seq_id') and sf['seq_id'] != 'unknown':
                print(f"      序列ID: {sf['seq_id']}")
        
        # 检查NTP匹配
        matches = self.check_ntp_matches(found_superfamilies)
        
        print(f"\n🧬 NTP-processing 匹配结果:")
        print("=" * 50)
        
        if matches:
            print(f"✅ 找到 {len(matches)} 个NTP-processing相关的超家族!")
            for i, match in enumerate(matches, 1):
                entry = match['matched_entry']
                print(f"\n   {i}. 发现的超家族: {match['found_name']}")
                if match['found_scop'] and match['found_scop'] != 'unknown':
                    print(f"      发现的SCOP ID: {match['found_scop']}")
                
                print(f"      ✅ 匹配的NTP超家族: {entry['sf_name']}")
                if entry['sf_id'] and entry['sf_id'] != 'N/A':
                    print(f"      ✅ 匹配的SCOP ID: {entry['sf_id']}")
                if entry['representative']:
                    print(f"      🧬 PDB ID: {entry['representative']}")
                if entry['reaction']:
                    print(f"      ⚗️  反应类型: {entry['reaction']}")
                if entry['coordination']:
                    print(f"      🔗 协调: {entry['coordination']}")
                print(f"      📊 匹配类型: {match['match_type']}")
            
            print(f"\n🔬 建议: 请进行 AlphaFold 预测找活性点!")
            return True
        else:
            print("❌ NO NTP-processing family found")
            return False
    
    def check_directory(self, directory_path):
        """检查目录中的所有HTML文件"""
        if not os.path.exists(directory_path):
            print(f"❌ 目录不存在: {directory_path}")
            return
        
        # 查找所有HTML文件
        html_files = []
        for file in os.listdir(directory_path):
            if file.endswith('.html') or file.endswith('.htm'):
                html_files.append(os.path.join(directory_path, file))
        
        if not html_files:
            print(f"❌ 在目录 {directory_path} 中没有找到HTML文件")
            return
        
        print(f"🧬 批量NTP-processing检查")
        print("=" * 60)
        print(f"📁 在 {directory_path} 中找到 {len(html_files)} 个HTML文件")
        print("=" * 60)
        
        # 统计结果
        results = {"total": len(html_files), "ntp_found": 0, "no_result": 0}
        
        for i, html_file in enumerate(html_files, 1):
            print(f"\n[{i}/{len(html_files)}] " + "="*60)
            
            result = self.check_single_html(html_file)
            if result:
                results["ntp_found"] += 1
            elif not result:
                results["no_result"] += 1
        
        # 生成总结报告
        print(f"\n{'='*60}")
        print("📊 批量检查总结报告")
        print(f"{'='*60}")
        print(f"总文件数: {results['total']}")
        print(f"找到NTP-processing: {results['ntp_found']}")
        print(f"未找到NTP-processing: {results['total'] - results['ntp_found']}")
        print(f"NTP匹配率: {results['ntp_found']/results['total']*100:.1f}%")
        
        if results['ntp_found'] > 0:
            print(f"\n🔬 建议: {results['ntp_found']} 个蛋白质需要进行 AlphaFold 预测!")
        
        print("=" * 60)

def main():
    # 创建检查器实例
    checker = NTPChecker()
    
    # 检查是否成功加载Excel文件
    if not checker.excel_path:
        print("\n💡 如果Excel文件在其他位置，可以指定路径:")
        print("python html_ntp_checker.py --excel /path/to/excel/file.xlsx file.html")
        return
    
    if len(sys.argv) == 1:
        # 无参数时，检查当前目录或html_results目录
        possible_dirs = ["html_results", ".", "downloaded_results"]
        
        for dir_path in possible_dirs:
            if os.path.exists(dir_path):
                html_files = [f for f in os.listdir(dir_path) if f.endswith('.html')]
                if html_files:
                    checker.check_directory(dir_path)
                    return
        
        print("❌ 没有找到HTML文件")
        print("💡 使用方法:")
        print("   python html_ntp_checker.py                    # 自动查找HTML文件")
        print("   python html_ntp_checker.py file.html          # 检查单个文件")
        print("   python html_ntp_checker.py html_results/      # 检查目录")
        
    elif len(sys.argv) >= 2:
        # 处理命令行参数
        if sys.argv[1] == "--excel" and len(sys.argv) >= 4:
            # 指定Excel文件路径
            excel_path = sys.argv[2]
            target = sys.argv[3]
            checker = NTPChecker(excel_path)
            
            if not checker.excel_path:
                print(f"❌ 指定的Excel文件不存在: {excel_path}")
                return
        else:
            target = sys.argv[1]
        
        if os.path.isfile(target):
            # 检查单个文件
            print(f"🧬 单文件NTP-processing检查")
            print("=" * 50)
            checker.check_single_html(target)
            
        elif os.path.isdir(target):
            # 检查目录
            checker.check_directory(target)
            
        else:
            print(f"❌ 文件或目录不存在: {target}")
    else:
        print("使用方法:")
        print("  python html_ntp_checker.py                              # 自动查找HTML文件") 
        print("  python html_ntp_checker.py file.html                    # 检查单个文件")
        print("  python html_ntp_checker.py html_results/                # 检查目录")
        print("  python html_ntp_checker.py --excel /path/to/excel.xlsx file.html  # 指定Excel路径")

if __name__ == "__main__":
    # 检查依赖
    try:
        import pandas as pd
        from bs4 import BeautifulSoup
    except ImportError as e:
        print("❌ 缺少必要的Python包，请安装:")
        print("pip install pandas openpyxl beautifulsoup4 lxml")
        print(f"错误详情: {e}")
        sys.exit(1)
    
    main()
