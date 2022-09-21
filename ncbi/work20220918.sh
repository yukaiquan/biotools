bash pieplot.sh
cat *.ann_detail.tsv > Avena_oat_all_ann_detail.tsv
python3 group_by_sum.py Avena_oat_all_ann_detail.tsv 1 0 Avena_oat_all_ann_detail_sum.tsv
