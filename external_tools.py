import subprocess 

def call_echo(message):
    """
    Example function that calls the 'echo' command,
    simulating a real eternal tool call.
    """
    cmd = ['echo', message]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return result.stdout.strip() 

if __name__ == "__main__":
    output = call_echo("Testing echo command")
    print("output from echo:", output) 
