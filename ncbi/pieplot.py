# -*- coding: utf-8 -*-
from palettable.cartocolors.qualitative import Vivid_9
import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')


def draw_pie(labels: list, quants: list, title: str, colors: list = ["coral", "blue", "red", "green",
                                                                     "yellow", "orange"], output: str = "pie.png"):
    # make a square figure
    plt.figure(1, figsize=(15, 10))
    # For China, make the piece explode a bit
    expl = [0.1, 0, 0, 0, 0, 0, 0, 0]  # 第二块即China离开圆心0.1
    # Colors used. Recycle if not enough.
    # Pie Plot
    # autopct: format of "percent" string;百分数格式
    plt.pie(quants, explode=expl, colors=colors, labels=labels,
            autopct='%1.1f%%', pctdistance=0.8, shadow=True)
    plt.title(title, bbox={'facecolor': '0.8', 'pad': 5})
    # plt.show()
    plt.savefig(output)
    plt.close()
    return True


def main():
    input_data = np.loadtxt(sys.argv[1], dtype=str, delimiter="\t")
    title = sys.argv[2]
    output = sys.argv[3]
    input_nd = input_data[:8]
    # print(input_nd)
    labels = input_nd[:, 1]
    quants = input_nd[:, 0].astype(float)
    draw_state = draw_pie(labels=labels, quants=quants,
                          title=title, colors=Vivid_9.hex_colors, output=output)
    if draw_state:
        print("Draw pie plot successfully!")
    else:
        print("Draw pie plot failed!")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("User Interrupted!")
        exit(1)
