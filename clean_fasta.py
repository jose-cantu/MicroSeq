import argparse 

def clean_fasta(input_file, output_file):
    """
    Remove any extra whitespace in sequence lines, keep headers intact. 
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                outfile.write(line.strip() + "\n")
            else:
                clean_line = "".join(line.split())
                outfile.write(clean_line + "\n")

def main():
    parser = argparse.ArgumentParser(description="Cleaning a FASTA file.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA path")
    parser.add_argument("-i", "--output", required=True, help="Output FASTA path")
    args = parser.parse_args() 

    clean_fasta(args.input, args.output)

if __name__ == "__main__":
    main() 
