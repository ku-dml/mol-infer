import argparse
import json
import sys
from dataclasses import dataclass, asdict
from typing import List, Union, Any, Dict, Type
import yaml
import os

# Define Dataclasses

@dataclass
class CommonConfig:
    output_prefix: str
    instance_file: str
    fringe_tree_file: str

@dataclass
class BasePropertyConfig:
    property: str
    model: str
    desc_filename: str
    fringe_filename: str
    values_filename: str
    target_value_lower_bound: float
    target_value_upper_bound: float
    sdf_filename: str

@dataclass
class ANNConfig(BasePropertyConfig):
    biases_filename: str
    weights_filename: str
    desc_norm_selected_filename: str

@dataclass
class LRConfig(BasePropertyConfig):
    LR_filename: str
    desc_norm_filename: str

@dataclass
class RFConfig(BasePropertyConfig):
    rf_filename: str
    desc_norm_selected_filename: str

@dataclass
class Config:
    common: CommonConfig
    input_data: List[BasePropertyConfig]

# Helper Functions for Dataclass Conversion

CONFIG_MODEL_MAP: Dict[str, Type[BasePropertyConfig]] = {
    'ANN': ANNConfig,
    'LR': LRConfig,
    'RF': RFConfig
}

def dict_to_common_config(data: Dict[str, Any]) -> CommonConfig:
    return CommonConfig(
        output_prefix=data.get('output_prefix'),
        instance_file=data.get('instance_file'),
        fringe_tree_file=data.get('fringe_tree_file')
    )

def dict_to_property_config(data: Dict[str, Any]) -> BasePropertyConfig:
    model = data.get('model')
    
    if not model:
        raise ValueError(f"model is not specified in property configuration: {data}")
    
    config_class = CONFIG_MODEL_MAP.get(model)

    if not config_class:
        raise ValueError(f"Unsupported model type: {model}")
    return config_class(**data)

def config_to_dict(config: Config) -> Dict[str, Any]:
    return {
        "output_prefix": config.common.output_prefix,
        "instance_file": config.common.instance_file,
        "fringe_tree_file": config.common.fringe_tree_file,
        'input_data': [asdict(prop) for prop in config.input_data]
    }

def dict_to_config(data: Dict[str, Any]) -> Config:
    common = dict_to_common_config(data)
    input_data = [dict_to_property_config(prop) for prop in data.get('input_data', [])]
    return Config(common=common, input_data=input_data)

def load_config(input_yaml_file: str) -> Config:
    with open(input_yaml_file, 'r') as f:
        config_dict = yaml.safe_load(f)
    return dict_to_config(config_dict)

def save_config(config: Config, output_file: str):
    with open(output_file, 'w') as f:
        yaml.dump(config_to_dict(config), f, sort_keys=False)
    print(f'Config file saved at {output_file}')

# File Checking Functions

def validate_file_extension(filename: str, ext: str) -> str:
    if not filename.lower().endswith(ext):
        raise argparse.ArgumentTypeError(f'File must have a {ext} extension')
    return filename

def yaml_file(filename: str) -> str:
    return validate_file_extension(filename, '.yaml')

# Configuration Manipulation Functions

def make_sample_config(output_file: str):
    sample_common = CommonConfig(
        output_prefix='Module_3/samples/output/test',
        instance_file='Module_3/samples/input/instance_files/instance_a_2LMM.txt',
        fringe_tree_file='Module_3/samples/input/fringe_files/ins_a_fringe_2LMM.txt'
    )

    sample_properties: List[BasePropertyConfig] = [
        LRConfig(
            property='Mp',
            model='LR',
            LR_filename='Module_3/samples/input/data/Mp_small_LR.txt',
            desc_filename='Module_3/samples/input/data/Mp_small_desc.csv',
            desc_norm_filename='Module_3/samples/input/data/Mp_small_desc_norm.csv',
            fringe_filename='Module_3/samples/input/data/Mp_small_fringe.txt',
            values_filename='Module_3/samples/input/data/Mp_small_values.txt',
            target_value_lower_bound=-2000.0,
            target_value_upper_bound=1000.0,
            sdf_filename='Module_3/samples/input/data/Mp_small.sdf'
        ),
        RFConfig(
            property='Bp',
            model='RF',
            rf_filename='Module_3/samples/input/data/Bp_small_RF.txt',
            desc_filename='Module_3/samples/input/data/Bp_small_desc.csv',
            desc_norm_selected_filename='Module_3/samples/input/data/Bp_small_desc_norm_selected.csv',
            fringe_filename='Module_3/samples/input/data/Bp_small_fringe.txt',
            values_filename='Module_3/samples/input/data/Bp_small_values.txt',
            target_value_lower_bound=0.0,
            target_value_upper_bound=4500.0,
            sdf_filename='Module_3/samples/input/data/Bp_small.sdf'
        ),
        ANNConfig(
            property='Fp',
            model='ANN',
            biases_filename='Module_3/samples/input/data/Fp_small_biases.txt',
            weights_filename='Module_3/samples/input/data/Fp_small_weights.txt',
            desc_filename='Module_3/samples/input/data/Fp_small_desc.csv',
            desc_norm_selected_filename='Module_3/samples/input/data/Fp_small_desc_norm_selected.csv',
            fringe_filename='Module_3/samples/input/data/Fp_small_fringe.txt',
            values_filename='Module_3/samples/input/data/Fp_small_values.txt',
            target_value_lower_bound=-100.0,
            target_value_upper_bound=0.0,
            sdf_filename='Module_3/samples/input/data/Fp_small.sdf'
        )
    ]

    sample_config = Config(
        common=sample_common,
        input_data=sample_properties
    )

    save_config(sample_config, output_file)
    print(f'Sample config file created at {output_file}')

def parse_key_value(arg: str) -> tuple:
    if '=' not in arg:
        raise argparse.ArgumentTypeError("Arguments must be in key=value format")
    key, value = arg.split('=', 1)
    return key, parse_value(value)

def parse_value(value: str) -> Any:
    """
    Converts strings to appropriate types.
    Handles numbers, boolean values, and null.
    """
    value = value.strip()
    if value.lower() == 'null':
        return None
    elif value.lower() == 'true':
        return True
    elif value.lower() == 'false':
        return False
    else:
        try:
            if '.' in value:
                return float(value)
            else:
                return int(value)
        except ValueError:
            return value  # Return as string if conversion fails


def make_common(output_file: str, key_values: List[tuple]):
    data_dict = dict(key_values)
    try:
        common_config = dict_to_common_config(data_dict)
    except TypeError as e:
        print(f'Invalid common parameters: {e}')
        sys.exit(1)

    config = Config(common=common_config, input_data=[])
    save_config(config, output_file)
    print(f'Common parameters created at {output_file}')

def edit_common(input_file: str, key_values: List[tuple]):
    config = load_config(input_file)
    data_dict = dict(key_values)

    for key, value in data_dict.items():
        if hasattr(config.common, key):
            setattr(config.common, key, value)
            print(f'Updated common parameter: {key} = {value}')
        else:
            print(f'Invalid key "{key}" for common configuration.')
            sys.exit(1)

    save_config(config, input_file)
    print(f'Common parameters updated in {input_file}')

def edit_property(input_file: str, property_name: str, key_values: List[tuple]):
    config = load_config(input_file)
    data_dict = dict(key_values)

    for prop in config.input_data:
        if prop.property == property_name:
            for key, value in data_dict.items():
                if hasattr(prop, key):
                    setattr(prop, key, value)
                    print(f'Updated property "{property_name}": {key} = {value}')
                else:
                    print(f'Invalid key "{key}" for property "{property_name}".')
                    sys.exit(1)
            save_config(config, input_file)
            print(f'Property "{property_name}" updated in {input_file}')
            return

    print(f'Property "{property_name}" not found.')
    sys.exit(1)

def add_property(input_file: str, key_values: List[tuple]):
    config = load_config(input_file)
    data_dict = dict(key_values)

    try:
        new_property = dict_to_property_config(data_dict)
    except (ValueError, TypeError) as e:
        print(f'Error adding property: {e}')
        sys.exit(1)

    if any(prop.property == new_property.property for prop in config.input_data):
        print(f'Property "{new_property.property}" already exists.')
        sys.exit(1)

    config.input_data.append(new_property)
    save_config(config, input_file)
    print(f'Property "{new_property.property}" added to {input_file}')

def delete_property(input_file: str, property_name: str):
    config = load_config(input_file)
    original_length = len(config.input_data)
    config.input_data = [prop for prop in config.input_data if prop.property != property_name]

    if len(config.input_data) == original_length:
        print(f'Property "{property_name}" not found.')
        sys.exit(1)

    save_config(config, input_file)
    print(f'Property "{property_name}" deleted from {input_file}')


def setup_args() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='Tool to edit config file by script')
    subparsers = parser.add_subparsers(dest='command', required=True, help='Sub-command help')

    # Subcommand 'make_sample'
    subparsers.add_parser('make_sample', help='Make sample config file').add_argument(
        '-o', '--output_yaml_file', type=yaml_file, required=True, help='Output YAML file'
    )

    # Subcommand 'make_common'
    make_common_parser = subparsers.add_parser('make_common', help='Create config file with common parameters')
    make_common_parser.add_argument('-o', '--output_yaml_file', type=yaml_file, required=True, help='Output YAML file')
    make_common_parser.add_argument('--set', type=parse_key_value, nargs='+', help='Set key=value pairs', required=True)

    # Subcommand 'edit_common'
    edit_common_parser = subparsers.add_parser('edit_common', help='Edit common parameters')
    edit_common_parser.add_argument('-i', '--input_yaml_file', type=yaml_file, required=True, help='Target YAML file')
    edit_common_parser.add_argument('--set', type=parse_key_value, nargs='+', help='Set key=value pairs', required=True)

    # Subcommand 'edit_property'
    edit_property_parser = subparsers.add_parser('edit_property', help='Edit property parameters')
    edit_property_parser.add_argument('-i', '--input_yaml_file', type=yaml_file, required=True, help='Target YAML file')
    edit_property_parser.add_argument('-p', '--property', type=str, required=True, help='Property name to edit')
    edit_property_parser.add_argument('--set', type=parse_key_value, nargs='+', help='Set key=value pairs', required=True)

    # Subcommand 'add_property'
    add_property_parser = subparsers.add_parser('add_property', help='Add property parameters')
    add_property_parser.add_argument('-i', '--input_yaml_file', type=yaml_file, required=True, help='Target YAML file')
    add_property_parser.add_argument('--set', type=parse_key_value, nargs='+', help='Set key=value pairs', required=True)

    # Subcommand 'delete_property'
    delete_property_parser = subparsers.add_parser('delete_property', help='delete property parameters')
    delete_property_parser.add_argument('-i', '--input_yaml_file', type=yaml_file, required=True, help='Target YAML file')
    delete_property_parser.add_argument('-p', '--property', type=str, required=True, help='Property name to delete')

    return parser

def main():
    parser = setup_args()
    args = parser.parse_args()

    command = args.command
    if command == 'make_sample':
        make_sample_config(args.output_yaml_file)
    elif command == 'make_common':
        make_common(args.output_yaml_file, args.set)
    elif command == 'edit_common':
        edit_common(args.input_yaml_file, args.set)
    elif command == 'edit_property':
        edit_property(args.input_yaml_file, args.property, args.set)
    elif command == 'add_property':
        add_property(args.input_yaml_file, args.set)
    elif command == 'delete_property':
        delete_property(args.input_yaml_file, args.property)
    else:
        parser.print_help()

if __name__ == '__main__':
    main()
