import argparse

def print_help():
    """Print help message."""
    print ("\nUsage: python -m Module_3 --config-name=[config]")
    print ("Example: python -m Module_3 --config-name=config1")
    print ("config file should be in config/ directory\n")

if __name__ == "__main__":
    import sys
    from .main_infer import infer
    
    """Main function."""
    ########## preparation ##########
    # Read config file
    if len(sys.argv) != 2:
        print_help()
        raise Exception("Wrong number of arguments")

    infer()