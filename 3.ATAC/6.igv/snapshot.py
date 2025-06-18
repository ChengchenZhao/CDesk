import argparse
import subprocess
import os
import tempfile
import sys
import json

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run IGV snapshots visualization")
    parser.add_argument("-i",required=True, type=str, help="Specify the input for IGV snapshots, separate by ,")
    parser.add_argument("-o",required=True,type=str, help="Output directory")
    parser.add_argument("--genome", required=True, type=str, help="Select the genome id")
    parser.add_argument("--display_mode",default='collapse',help="Display mode", choices=["collapse","expand","squish"])
    parser.add_argument("--region",required=True,help="The region file")
    parser.add_argument("--type", default='png', type=str,choices=['svg','png'] ,help="svg or png output")
    parser.add_argument("config_path", help="Path to configuration file")
    return parser.parse_args()

def main():
    # 解析命令行参数
    args = parse_arguments()
    os.makedirs(args.o, exist_ok=True)

    igv_script = f"""
    new
    genome {args.genome}
    snapshotDirectory {args.o}
    """

    for input_file in (args.i).split(','):
        input_file = input_file.strip()
        if os.path.exists(input_file):
            if input_file.endswith('.bam'):
                bam_index_file = input_file + ".bai"
                if not os.path.exists(bam_index_file):
                    print('没有找到对应的bai文件，请检查')
                    sys.exit(1)
            igv_script += f"load {input_file}\n    "
        else:
            print(f"{input_file}不存在")
            sys.exit(1)

    with open(args.region, 'r') as file:
        lines = [line.strip() for line in file if line.strip()]
    regions = ','.join(lines)  

    for region in regions.split(','):
        region_save = region.replace(':','_')
        region_save = region_save.replace('-','_')
        igv_script += f"""
    goto {region}
    {args.display_mode}
    sort position
    snapshot {region_save+'.'+args.type}
    """
    igv_script += """
    exit
    """

    with open(args.config_path, "r") as f:
        config = json.load(f)
    igv_path = os.path.join(config['software']['igv'],'igv.sh')
    
    script_file = os.path.join(f'{args.o}/igv_batch_script.txt')
    try:
        with open(script_file, "w") as f:
            f.write(igv_script)
    except Exception as e:
        print(f"Error writing IGV script to file: {e}")
        return

    xvfb_run_path = config['software']['xvfb-run']
    xvfb_cmd = [
        xvfb_run_path,
        "-a",
        igv_path,
        "-b", script_file
    ]

    try:
        subprocess.run(xvfb_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error during IGV execution: {e}")
        return
    except Exception as e:
        print(f"Unexpected error: {e}")
        return

    print("Snapshot generation completed successfully.")

if __name__ == "__main__":
    main()
