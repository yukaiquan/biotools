cp target/x86_64-unknown-linux-musl/ example/
cp target/x86_64-unknown-linux-musl/release/merge_chr_pos example/


merge_chr_pos.exe -i input.txt -o output_txt.txt
merge_chr_pos -i input.txt -o output_txt.txt