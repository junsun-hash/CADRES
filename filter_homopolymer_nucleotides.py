#!/usr/bin/env python
import argparse
import sys
import os
import subprocess
import tempfile

def check_homopolymer(sequence, edit_nuc):
    """
    检查变异是否位于同聚物区域
    序列长度为9bp，变异位点在索引4的位置
    原Perl脚本检查5个窗口，每个窗口4个碱基
    """
    s = sequence.upper()
    edit_nuc = edit_nuc.upper()
    
    # 窗口1: -4, -3, -2, -1 (索引0,1,2,3)
    if all(base == edit_nuc for base in (s[0], s[1], s[2], s[3])):
        return True
    # 窗口2: -3, -2, -1, +1 (索引1,2,3,5)
    if all(base == edit_nuc for base in (s[1], s[2], s[3], s[5])):
        return True
    # 窗口3: -2, -1, +1, +2 (索引2,3,5,6)
    if all(base == edit_nuc for base in (s[2], s[3], s[5], s[6])):
        return True
    # 窗口4: -1, +1, +2, +3 (索引3,5,6,7)
    if all(base == edit_nuc for base in (s[3], s[5], s[6], s[7])):
        return True
    # 窗口5: +1, +2, +3, +4 (索引5,6,7,8)
    if all(base == edit_nuc for base in (s[5], s[6], s[7], s[8])):
        return True
        
    return False


def filter_homopolymers(infile, outfile, refgenome_path):
    """
    过滤位于同聚物区域的变异位点
    Python替换filter_homopolymer_nucleotides.pl，使用fastaFromBed工具
    """
    left_buffer = 4
    right_buffer = 4
    
    # 创建临时文件
    tmp_bed = infile + '.nuctmp'
    tmp_fa = infile + '.tmpnuctmp'
    
    try:
        with open(infile, 'r') as sites_file, \
             open(outfile, 'w') as output_file, \
             open(f"{outfile}_failed", 'w') as failed_file:
            
            for line in sites_file:
                line = line.strip()
                if not line:
                    continue
                
                fields = line.split()
                chrom = fields[0]
                pos = int(fields[1])
                edit_nuc = fields[4]
                
                # 计算BED格式的起始和结束位置
                # BED格式是0-based，半开区间 [start, end)
                # Perl脚本: startpos = position - leftbuffer, endpos = position + rightbuffer + 1
                startpos = pos - left_buffer
                endpos = pos + right_buffer + 1
                
                # 写入临时BED文件
                with open(tmp_bed, 'w') as bed_file:
                    bed_file.write(f"{chrom}\t{startpos}\t{endpos}\n")
                
                # 使用fastaFromBed提取序列
                cmd = [
                    "fastaFromBed",
                    "-fi", refgenome_path,
                    "-bed", tmp_bed,
                    "-fo", tmp_fa
                ]
                
                try:
                    subprocess.run(cmd, check=True, capture_output=True, text=True)
                except subprocess.CalledProcessError as e:
                    sys.stderr.write(f"Error running fastaFromBed: {e.stderr}\n")
                    continue
                except FileNotFoundError:
                    sys.stderr.write("Error: fastaFromBed not found. Please ensure BEDTools is installed.\n")
                    sys.exit(1)
                
                # 读取提取的序列
                seq = ''
                with open(tmp_fa, 'r') as fa_file:
                    for fa_line in fa_file:
                        fa_line = fa_line.strip()
                        if not fa_line.startswith('>'):
                            seq += fa_line
                
                # 检查同聚物
                is_homopolymer = check_homopolymer(seq, edit_nuc)
                
                if not is_homopolymer:
                    output_file.write(line + '\n')
                else:
                    failed_file.write(line + '\n')
                
                # 删除临时文件
                try:
                    os.remove(tmp_fa)
                    os.remove(tmp_bed)
                except OSError:
                    pass
    
    except IOError as e:
        sys.stderr.write(f"Error opening files: {e}\n")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="同聚物过滤器。Python替换filter_homopolymer_nucleotides.pl。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-i", "--infile", required=True, help="包含待过滤变异位点的文件")
    parser.add_argument("-o", "--outfile", required=True, help="过滤后变异位点的输出文件")
    parser.add_argument("-r", "--refgenome", required=True, help="参考基因组FASTA格式文件")
    
    args = parser.parse_args()
    filter_homopolymers(args.infile, args.outfile, args.refgenome)
    print("同聚物过滤完成。")

if __name__ == '__main__':
    main()