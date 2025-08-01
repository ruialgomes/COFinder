import re
import os
import stat
import shutil

# Function to format time and print
def format_time(time_in_sec):
    minutes = int(time_in_sec // 60)
    seconds = int(time_in_sec % 60)
    if minutes == 0:
        return f"{seconds} seconds"
    elif seconds == 0:
        return f"{minutes} minutes"
    else:
        return f"{minutes} minutes and {seconds} seconds"


# write to an output file
def write_output(path, output):
    with open(path, "a+") as outfile:
        outfile.write(output)


# Function to replace forbidden characters in a identifier
def change_chars(input_text):
    searched_chars = r'[\\/:"*?<>|,\n\t^.;|@_# ]'
    return re.sub(searched_chars, '-', input_text)


# Function to check sequence length
def check_len(tested_seq_len, fixed_len, max_prot_len_perc):
    lower_bound = fixed_len * (1 - max_prot_len_perc)
    upper_bound = fixed_len * (1 + max_prot_len_perc)
    return lower_bound <= tested_seq_len <= upper_bound

def ensure_executable(path):
    mode = os.stat(path).st_mode
    os.chmod(path, mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

def create_dirs(output_path,output_orfs_path):
    if os.path.exists(output_orfs_path):
        shutil.rmtree(output_orfs_path)
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.makedirs(output_orfs_path, exist_ok=True)
    os.makedirs(output_path, exist_ok=True)
