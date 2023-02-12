#!/usr/bin/env python3.6
# -*- coding: utf-8 -*-

#导入模块，初始传递命令、变量等
import argparse
import re

parser = argparse.ArgumentParser(description = '该脚本用于统计 fasta 文件中所含序列的基本信息', add_help = False, usage = '\npython3.6 fasta_stat.py -i [input.fasta] -o [output.file]\npython3.6 fasta_stat.py --input [input.fasta] --output [output.file]')
required = parser.add_argument_group()
optional = parser.add_argument_group()
required.add_argument('-i', '--input', metavar = '[input.fasta]', help = '输入 fasta 文件', required = True)
required.add_argument('-o', '--output', metavar = '[output.file]', help = '输出统计结果', required = True)
optional.add_argument('-h', '--help', action = 'help', help = '帮助信息')
args = parser.parse_args()

dict = {}; dict_refer = {}
seqs_sum = 0
base_sum = 0
GC_sum = 0
N_sum = 0
Ln = 0; base_sum_n = 0
N50 = 0; N75 = 0; N90 = 0
n1000000 = 0; n1000000_len = 0
n100000 = 0; n100000_len = 0
n10000 = 0; n10000_len = 0
n1000 = 0; n1000_len = 0

#读入 fasta 文件，统计 scaffolds（每条 + 所有）总数、长度、GC 含量等基本信息
with open(args.input, 'r') as read_fas:
	for line in read_fas:
		if line[0] == '>':
			key = line.strip('[ >\n]')
			dict[key] = 0
			seqs_sum += 1
		else:
			value = line.strip()
			seqs_len = len(value)
			base_sum += seqs_len
			dict[key] += seqs_len
			GC_sum += len(re.findall('[gcGC]', value))
			N_sum += len(re.findall('[nN]', value))

read_fas.close()

#统计 N50、N75、N90 以及各长度区间 scaffolds 数量及长度
dict_sort = sorted(dict.items(), key = lambda x:x[1], reverse = True)
base_sum_n50 = 5 * base_sum/10
base_sum_n75 = 7.5 * base_sum/10
base_sum_n90 = 9 * base_sum/10

for value in dict_sort:
	base_sum_n += value[1]
	Ln += 1
	if N50:
		if N75:
			if base_sum_n >= base_sum_n90:
				N90 = value[1]
				L90 = Ln
				break
		else:
			if base_sum_n >= base_sum_n75:
				N75 = value[1]
				L75 = Ln
	else:
		if base_sum_n >= base_sum_n50:
			N50 = value[1]
			L50 = Ln

for value in dict_sort:
	if value[1] >= 1000000:
		n1000000 += 1
		n1000000_len += value[1]
	if value[1] >= 100000:
		n100000 += 1
		n100000_len += value[1]
	if value[1] >= 10000:
		n10000 += 1
		n10000_len += value[1]
	if value[1] >= 1000:
		n1000 += 1
		n1000_len += value[1]
	if value[1] < 1000:
		break

#统计总览，并输出结果
basic_stat = open(args.output, 'w')
print(f'Total scaffolds\t{seqs_sum}', file = basic_stat)
print(f'Total base (bp)\t{base_sum}', file = basic_stat)
print(f'Total N (bp)\t{N_sum}', file = basic_stat)
print(f'Average length (bp)\t{int(base_sum/seqs_sum)}', file = basic_stat)
print(f'Longest scaffold (bp)\t{dict_sort[0][1]}', file = basic_stat)
print(f'Shortest scaffold (bp)\t{dict_sort[-1][1]}', file = basic_stat)
print(f'L50\t{L50}', file = basic_stat)
print(f'N50\t{N50}', file = basic_stat)
print(f'L75\t{L75}', file = basic_stat)
print(f'N75\t{N75}', file = basic_stat)
print(f'L90\t{L90}', file = basic_stat)
print(f'N90\t{N90}', file = basic_stat)
print(f'GC (%)\t{round(100 * GC_sum / base_sum, 2)}', file = basic_stat)
print(f'Total scaffolds (>= 1000000 bp)\t{n1000000}', file = basic_stat)
print(f'Total scaffolds (>= 100000 bp)\t{n100000}', file = basic_stat)
print(f'Total scaffolds (>= 10000 bp)\t{n10000}', file = basic_stat)
print(f'Total scaffolds (>= 1000 bp)\t{n1000}', file = basic_stat)
print(f'Total length (>= 1000000 bp)\t{n1000000_len}', file = basic_stat)
print(f'Total length (>= 100000 bp)\t{n100000_len}', file = basic_stat)
print(f'Total length (>= 10000 bp)\t{n10000_len}', file = basic_stat)
print(f'Total length (>= 1000 bp)\t{n1000_len}', file = basic_stat)
basic_stat.close()

