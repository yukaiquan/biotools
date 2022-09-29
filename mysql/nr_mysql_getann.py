#!/use/bin/env python3
# -*- coding: utf-8 -*-
# Author: 2022-2050, yukaiquan
# Created Time: 2022-09-29 15:00:00
# Description: 从mysqld的nr中获取注释信息


from pymysql_comm import UsingMysql


def get_count(cursor):
    cursor.execute(
        "select nr_id,annotation from nr_annotation where nr_id='A0A059TC02.1' ")

    # 使用 fetchone() 方法获取单条数据.
    data = cursor.fetchone()

    print("-- 当前数量: %s " % data['annotation'])


with UsingMysql(log_time=True) as um:
    get_count(um.cursor)
