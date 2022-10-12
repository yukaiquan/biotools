#!/use/bin/env python
# -*- coding: utf-8 -*-
'''
Copyright [yukaiquan 1962568272@qq.com]
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
    http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.


This script is used to convert kegg_k_map.json to kegg_k_map.tsv
usage: python kegg_k_map.py ko00001.json kegg_k_map.tsv
output_file = `
09100 Metabolism	09101 Carbohydrate metabolism	ko00010	Glycolysis / Gluconeogenesis	K00844	HK; hexokinase	[EC:2.7.1.1]
09100 Metabolism	09101 Carbohydrate metabolism	ko00010	Glycolysis / Gluconeogenesis	K12407	GCK; glucokinase	[EC:2.7.1.2]
09100 Metabolism	09101 Carbohydrate metabolism	ko00010	Glycolysis / Gluconeogenesis	K00845	glk; glucokinase	[EC:2.7.1.2]
09100 Metabolism	09101 Carbohydrate metabolism	ko00010	Glycolysis / Gluconeogenesis	K25026	glk; glucokinase	[EC:2.7.1.2]
`
'''


import json
import re
import sys


def json_to_list(ko_map_data: dict) -> list:
    """Convert json data to list"""
    ko_map_list: list = []
    for data in ko_map_data['children']:
        leve_1 = data['name']
        for data2 in data['children']:
            leve_2 = data2['name']
            for data3 in data2['children']:
                leve_3 = data3['name']
                # {'name': '01000 Enzymes [BR:ko01000]'} has no 'children' so we need to use 'if'
                if 'children' in data3:
                    for data4 in data3['children']:
                        K_data = data4['name']
                        ko_map_list.append([leve_1, leve_2, leve_3, K_data])
                else:
                    ko_map_list.append([leve_1, leve_2, leve_3, ''])

    return ko_map_list


def list_to_tsv(ko_map_list: list, output_file: str) -> bool:
    """
    :param ko_map_list: list
    :param output_file: str
    :return: str
    """
    output_tsv = open(output_file, 'w')
    for l in ko_map_list:
        level1 = l[0]
        level2 = l[1]
        level3 = l[2]
        # test_re = '00010 Glycolysis / Gluconeogenesis [PATH:ko00010]'
        # 01000 Enzymes [BR:ko01000] has no K number
        if '01000 Enzymes [BR:ko01000]' == level3:
            continue
        try:
            path = 'ko' + \
                re.search(r'(\d+) ([\w\W ]+) (\[\w+:\w+\d+\])', l[2]).group(1)
        except AttributeError:
            if '99980' in level3:
                path = 'k099980'
            elif '09113' in level3:
                path = 'ko09113'
        try:
            path_ann = re.search(
                r'(\d+) ([\w\W ]+) (\[\w+:\w+\d+\])', l[2]).group(2)
        except AttributeError:
            if '99980' in level3:
                path_ann = 'Enzymes with EC numbers'
            elif '09113' in level3:
                path_ann = 'Global maps only'
        if '' == l[3]:
            k_id = ''
            k_ann = ''
            ec_id = ''
        else:
            if 'K86952' in l[3]:
                k_id = 'K86952'
                k_ann = ''
                ec_id = ''
                continue
            else:
                try:
                    k_id = re.search(
                        r'(K\d{5})  (.*) (\[EC:.*\])', l[3]).group(1)
                    k_ann = re.search(
                        r'(K\d{5})  (.*) (\[EC:.*\])', l[3]).group(2)
                    ec_id = re.search(
                        r'(K\d{5})  (.*) (\[EC:.*\])', l[3]).group(3)
                except AttributeError:
                    k_id = re.search(r'(K\d{5})  (.*)', l[3]).group(1)
                    k_ann = re.search(r'(K\d{5})  (.*)', l[3]).group(2)
                    ec_id = ''
        output_tsv.write(
            f'{level1}\t{level2}\t{path}\t{path_ann}\t{k_id}\t{k_ann}\t{ec_id}\n')
    output_tsv.close()
    return True


# test_re = 'K01810  GPI, pgi; glucose-6-phosphate isomerase [EC:5.3.1.9]'
# test_re = 'K00246  frdC; succinate dehydrogenase subunit C'

if __name__ == '__main__':
    try:
        inpuy_file = sys.argv[1]
        output_file = sys.argv[2]
        ko_map_data: dict = json.load(open(inpuy_file, 'r'))
        ko_map_list: list = json_to_list(ko_map_data)
        list_to_tsv(ko_map_list, output_file)
    except KeyboardInterrupt:
        print('KeyboardInterrupt')
        sys.exit(0)
