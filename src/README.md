# 说明

## 一、安装

```bash
podman machine init
podman machine start


unzip docker-images-tar.zip
tar -zxvf x86-64-images.tar.gz
podman load < nginx:latest-amd64.tar
podman load < redis:latest-amd64.tar
podman load < mysql:8.0.33-amd64.tar

podman run -itd --name cdesk_nginx -v /Users/chengchenzhao/CDesk/src/nginx.conf:/etc/nginx/nginx.conf -v /Users/chengchenzhao/CDesk/src/html:/usr/share/nginx/html --network host nginx
podman run -itd --name cdesk_redis -v /Users/chengchenzhao/CDesk/src/redis.conf:/etc/redis/redis.conf --network host redis
podman run -itd --name cdesk_mysql -e MYSQL_ROOT_PASSWORD=11223344 --network host mysql:8.0.33

podman exec -it cdesk_mysql /bin/bash

mysql -u root -p
alter user `root`@'%' identified with mysql_native_password by '11223344';
flush privileges;

exit
exit

curl https://raw.githubusercontent.com/mamba-org/micromamba-releases/main/install.sh | bash
# curl micro.mamba.pm/install.sh | bash
/Users/chengchenzhao/.local/bin/micromamba shell init -s bash -p ~/micromamba

micromamba activate

micromamba install -y python=3.9 -c conda-forge

~/micromamba/bin/python3 -m pip install markupsafe==2.0.1 cython fastapi python-multipart pandas pymysql celery redis uvicorn gunicorn paramiko DBUtils pycryptodome matplotlib seaborn pyg2plot scikit-learn umap-learn chardet bs4 requests adjustText imageio twobitreader pyg2plot PyPDF2==2.9 khmer -i https://mirrors.aliyun.com/pypi/simple

# CDesk_cli
sudo ln -s `pwd`/CDesk_cli /usr/local/bin/

# # ubuntu
# ufw allow 8888

export PATH=~/micromamba/bin:$PATH
```

<!-- ``` bash
# Effective Genome Size
unique-kmers.py hg38.fa
```
 -->

## 二、启动

```bash
cd /Users/chengchenzhao/CDesk/src

~/micromamba/bin/python3 start.py
```
