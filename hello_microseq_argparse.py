import argparse 

def main():
    parser = argparse.ArgumentParser(description="A test program for test arguments.")
    parser.add_argument("-n", "--name", default="MicroSeq!", help="Name to greet the person")
    args = parser.parse_args()
    print(f"Hello, {args.name}!")

if __name__ == "__main__":
    main() 
