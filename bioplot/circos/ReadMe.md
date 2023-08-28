## Circos Plot

### 1. install

#### Try conda

```bash
conda install -c bioconda circos
```

#### Try source code

```bash
# install perl
......
# install circos
......
failed
```

i try to install circos by code, but failed.
i try to install circos by source code, but failed.

#### Try singularity

```bash
singularity pull --arch amd64 library://yukaiquan/bioplot/circos:0.69.9
```

### 2. run

#### 1. mount file path

you must devide the path, for example:

```bash
$ pwd
/mnt/c/Workplace/soft/img
```

```bash
$ ls
circos.sif
example/
circos.png
circos.svg
```

#### 2. run

```bash
singularity exec -B /mnt/c/Workplace/soft/img/:/root/ circos.sif circos -conf example/etc/circos.conf
```

You will get the circos.png and circos.svg in the /mnt/c/Workplace/soft/img/ directory.
