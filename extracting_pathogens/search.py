#!/data/hb-austick/work/miniforge3/bin/python3

import os
import argparse

def search_words_in_files(file_list, words_file, output_file):
    # Clear the output file if it exists
    with open(output_file, 'w') as f:
        pass

    # Read the list of files
    with open(file_list, 'r') as f:
        files = [line.strip() for line in f]

    # Read the list of words
    with open(words_file, 'r') as f:
        words = [line.strip() for line in f]

    # Search for each word in each file
    with open(output_file, 'a') as out_f:
        for word in words:
            for file in files:
                if os.path.isfile(file):
                    with open(file, 'r') as in_f:
                        for line_num, line in enumerate(in_f, 1):
                            if word in line:
                                out_f.write(f"{file}:{line_num}:{line.strip()}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Search for words in a list of files.")
    parser.add_argument("file_list", help="File containing the list of .txt files")
    parser.add_argument("words_file", help="File containing the 4000 words")
    parser.add_argument("output_file", help="Output file to write the results")

    args = parser.parse_args()

    search_words_in_files(args.file_list, args.words_file, args.output_file)
