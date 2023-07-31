import os

def parse_filelist(filelist):
    try:
        with open(filelist, 'r') as f:
            # Check if the file is empty
            if not f.read().strip():
                raise ValueError("Error: The file is empty.")

            # Move the file pointer back to the beginning of the file
            f.seek(0)

            # Process the file content and create a dictionary
            file_dict = {}
            paired_flag = None
            for line in f.readlines():
                parts = line.strip().split('\t')
                if len(parts) != 2:
                    raise ValueError(f"Error: Incorrect format in line: {line.strip()}")

                # Check if the second column can be split by ";"
                second_col_parts = parts[1].split(';')
                if len(second_col_parts) == 2:
                    paired_flag = True
                elif len(second_col_parts) == 1:
                    paired_flag = False
                else:
                    raise ValueError(f"Error: Incorrect number of fields in the second column: {parts[1]}")

                file_dict[parts[0]] = second_col_parts

            # Check if all lines have consistent splitting results
            if paired_flag is not None:
                for line in f.readlines():
                    parts = line.strip().split('\t')
                    second_col_parts = parts[1].split(';')
                    if (len(second_col_parts) == 2 and not paired_flag) or (len(second_col_parts) == 1 and paired_flag):
                        raise ValueError("Error: Inconsistent splitting results in the second column.")
        
        for key, file_list in file_dict.items():
            for file_path in file_list:
                if not os.path.exists(file_path):
                    raise FileNotFoundError(f"file not found: {file_path}")
                elif os.path.getsize(file_path) == 0:
                    raise ValueError(f"file empty: {file_path}")
                
            return file_dict, paired_flag

    except FileNotFoundError:
        raise FileNotFoundError(f"Error: File '{filelist}' not found.")

def get_chromosomes(xg_file):
    cmd = 'vg paths --list -x '+ xg_file
    chromosomes = [line.strip() for line in os.popen(cmd).readlines()]
    return chromosomes

