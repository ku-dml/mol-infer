def print_help():
    """Print help message."""
    # print("Usage: python3 infer.py [config.yaml]")
    print ("Usage: python -m Module_3 [config.yaml]")
    print ("Example: python -m Module_3 config/config.yaml")

if __name__ == "__main__":
    import sys
    from .main_infer import infer
    
    """Main function."""
    ########## preparation ##########
    # Read config file
    if len(sys.argv) != 2:
        print_help()
        raise Exception("Wrong number of arguments")
    # Read config file
    config_filename: str = sys.argv[1]
    infer(config_filename)