#!/bin/bash
 
PATH=/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin:/usr/local/sbin:~/bin
export PATH
 
#判断要创建的文件是否存在，如果文件名不存在
if [ ! "$1" ]
	then
		echo '请输入要新建的文件名称，例如(sh AutoHead.sh Test.sh)'
		exit 1
fi 
#如果文件已经创建，直接用vim打开
if [ -f "$1" ]
	then
		vi "$1"
		exit 2
fi 
#创建定义的文件
touch "$1" 
#添加注释信息 
echo "#!/bin/bash">>"$1"
echo "">>"$1" 
echo "# Description : ">>"$1"
echo "# Author      : Chengchen Zhao">>"$1"
echo "# Version     : 1.0">>"$1"
echo "# Time        : `date "+%Y-%m-%d %H:%M:%S"`">>"$1" 
# echo "PATH=/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin:/usr/local/sbin:~/bin">>"$1"
# echo "export PATH">>"$1"
echo "">>"$1" 
echo "###################################################################################################">>"$1"
echo "################################     PARAMETER/FUNCTION DEFINE     ################################">>"$1"
echo "###################################################################################################">>"$1"
echo "">>"$1" 
#打开定义的文件
