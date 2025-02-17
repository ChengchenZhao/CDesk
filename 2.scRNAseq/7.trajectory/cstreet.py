import argparse
import warnings
import subprocess

warnings.filterwarnings('ignore')

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run CStreet analysis")
    
    parser.add_argument("-i", type=str, help="Input expression matrixes, separate by ',' ")
    parser.add_argument("-s", type=str, help="Input cell state file, separate by ',' ")
    parser.add_argument("-n", type=str,default='',help="Output directory name")
    parser.add_argument("-o", type=str, help="Output directory")
    return parser.parse_args()

def main():
    # 解析命令行参数
    args = parse_arguments()
    
    expMatrix = args.i.split(',')
    cellState = args.s.split(',')
    
    cmd = ['CStreet','-i',*expMatrix,'-s',*cellState,'-n',args.n,'-o',args.o]
    subprocess.run(cmd, check=True)

if __name__ == "__main__":

    main()
