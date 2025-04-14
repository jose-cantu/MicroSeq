import subprocess

def call_echo(message):
    """
    Example function that calls the 'echo' command simulating a real external tool call such as Trimmomatic or blastn or cap3.
    """
    cmd = ['echo', message]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True) 
    # capture_out=True tells Python to capture the commands standard output (stdout) and standard error (stderr) instead of printing them on your own screen. End up as result.stdout and result.stderr. 

    #text=True output treated as a string rather than raw bytes so result.stdout and result.stderr are nromal Python strings. 

    #check=True if echo or any command you run exits with a non-zero status indicating an error, Python raises a CalledProcessError exception. Good for error handling since you dont have to manually check the exit code. 
    return result.stdout.strip() 

if __name__ == "__main__":
    output = call_echo("Testing echo command")
    print("output from echo:", output) 
