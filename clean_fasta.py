import argparse 

def clean_fasta(input_file, output_file):
    """
    Function that reads and write an input file and an output FASTA file(s) and from there removes trailing and leading whitespace in header and removes all whitespace for sequenece lines. 
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile: 
        for line in infile:
            if line.startswith('>'):
                outfile.write(line.strip() + "\n")
            else:
                cleaned_line = "".join(line.split())
                outfile.write(cleaned_line + "\n")

def main():
    parser = argparse.ArgumentParser(description="Clean a FASTA file.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA path")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA path")
    args = parser.parse_args()

    clean_fasta(args.input, args.output)

if __name__ == "__main__":
    main() 

