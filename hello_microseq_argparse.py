import argparse 

def main():
    parser = argparse.ArgumentParser(description="Creating a test program with arguments involved.")
    parser.add_argument("-n", "--name", default="MicroSeq!", help="Name to greet")
    args = parser.parse_args() 
    print(f"Hello, {args.name}") 

if __name__ == "__main__":
    main() 

