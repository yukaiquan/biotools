#!/use/bin/env python3
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
input_file = `
Penicillin binding protein 3
Histidine decarboxylase
Nitric oxide synthase, inducible
`
output_file = `
Penicillin binding protein 3 'P33317', 'A0A2T2NDB6', 'A0A2N0V3U4', 'A0A8I2JQW5'
Histidine decarboxylase ['P0A9K8', 'A0A2T2NDB6', 'A0A2N0V3U4', 'A0A8I2JQW5',...]
Nitric oxide synthase, inducible ['P0A9K8', 'A0A2T2NDB6', 'A0A2N0V3U4', 'A0A8I2JQW5',...]
`
'''

import requests
from urllib import parse
import sys
import time


CONFIG = {
    "header_user_agent": "Mozilla/5.0 (iPhone; CPU iPhone OS 16_0 like Mac OS X) AppleWebKit/605.1.15 (KHTML, like Gecko) Mobile/15E148 MicroMessenger/8.0.27(0x18001b36) NetType/WIFI Language/zh_CN",
    "accept": "application/json",
}

# url = r"https://rest.uniprot.org/uniprotkb/search?size=1&query=Vasopressin%20V1a%20receptor"
# https://rest.uniprot.org/uniprotkb/search?size=5&query=Penicillin%20binding%20protein%203


def get_json_uniprot_id(url_list: list, CONFIG: dict) -> dict:
    uniport_dict: dict = {}
    uniport_list: list = []
    for url in url_list:
        # print(url)
        res_htlm = requests.get(url, headers=CONFIG, timeout=10)
        # 逆向解析
        name = parse.unquote(url.split('query=')[1])
        if res_htlm.status_code == 200:
            try:
                json_data = res_htlm.json()

                for uniport_id in json_data['results']:
                    # just one best match to dict
                    if name not in uniport_dict:
                        uniport_dict[name] = uniport_id['primaryAccession']
                    else:
                        break
                    # uniport_list.append(uniport_id['primaryAccession'])
            except:
                print(name, "not found in uniprot 404")

        else:
            print(
                f"Failed to get Uniport ID{name}, status code: {res_htlm.status_code}")
        # uniport_dict[name] = uniport_list
        time.sleep(0.01)
    return uniport_dict


def generator_url(uniportname_list: list) -> list:
    url_list: list = []
    for uniportname in uniportname_list:
        # 正向解析
        name_filds = parse.quote(uniportname)
        url_list.append(
            f"https://rest.uniprot.org/uniprotkb/search?size=2&query={name_filds}")
    return url_list


def main(input_id: str, output_file: str):
    start_time = time.time()
    uniportname_list: list = []
    with open(input_id, 'r') as f:
        for line in f.readlines():
            uniportname_list.append(line.strip())
    url_list = generator_url(uniportname_list)
    uniport_dict = get_json_uniprot_id(url_list, CONFIG)
    with open(output_file, 'w') as f:
        for key, value in uniport_dict.items():
            f.write(f"{key}\t{value}\n")
    stop_time = time.time()
    print(f"Total time: {stop_time-start_time}")


if __name__ == "__main__":
    try:
        input_id = sys.argv[1]
        output_file = sys.argv[2]
        main(input_id, output_file)
    except KeyboardInterrupt:
        print("KeyboardInterrupt")
        sys.exit(0)
