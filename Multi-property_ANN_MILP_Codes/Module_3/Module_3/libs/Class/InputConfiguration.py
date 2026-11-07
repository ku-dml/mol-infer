import sys
import yaml

class Config:
    # constructor
    def __init__(self, config_path):
        self.config_path = config_path
        # read config file
        try :
            with open(self.config_path, "r", encoding="utf-8") as file:
                self.config_data = yaml.safe_load(file)
        except Exception as ex:
            print(f"Error: {ex}")
            sys.exit(1)
        # get values
        self.__get_value("output_prefix")
        self.__get_value("instance_file")
        self.__get_value("fringe_tree_file")
        self.__get_value("input_data" , verbose=False)
            
    def __get_value(self, key, verbose=True):
        # set self.key = config_data[key] and save the value to self.key
        value = self.config_data.get(key)
        if value is None:
            print(f"Error: Please specify ({key}) in config.yaml")
            sys.exit(1)
        self.__dict__[key] = value
        if verbose:
            print(f"{key}: {value}")
    
    def get_with_index(self, index, key, verbose=True):
        # Get config.input_data[index][key]
        value = self.input_data[index].get(key)
        if value is None:
            print(f"\tError: Please specify ({key}) in config.yaml")
            sys.exit(1)
        if verbose:
            print(f"\t{key}: {value}")
        return value
