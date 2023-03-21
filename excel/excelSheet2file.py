import xlrd
import pandas as pd
import sys
import os

input_excel = sys.argv[1]


def excel2csv(excel_file):
    # 打开excel文件
    workbook = xlrd.open_workbook(excel_file)
    # 获取所有sheet名字
    sheet_names = workbook.sheet_names()
    for worksheet_name in sheet_names:
        # 遍历每个sheet并用Pandas读取
        data_xls = pd.read_excel(excel_file, worksheet_name, index_col=None)
        # 获取excel当前目录
        dir_path = os.path.abspath(os.path.dirname(excel_file))
        # 转换成csv并保存到excel所在目录下的csv文件夹中
        csv_path = dir_path+'\\csv\\'
        if not os.path.exists(csv_path):
            os.mkdir(csv_path)
        data_xls.to_csv(csv_path+worksheet_name+'.csv',
                        index=None, encoding='utf-8')


excel2csv(input_excel)
