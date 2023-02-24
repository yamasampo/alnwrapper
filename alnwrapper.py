import os
import sys

from typing import Dict, List, Callable, Union, Any

ALIGNER = {'clustalo', 'muscle'}
MODE = {'pw', 'pw_batch'}

def main(
        aligner, 
        mode, 
        input1='', 
        input2='', 
        output_folder='', 
        file_exist_ok=False, 
        **aligner_args):
    """Runs an alignment software """
    # Check if the given aligner is supported
    check_aligner(aligner)
    
    # Check if the given mode is supported
    check_mode(mode)

    # # Check output
    # check_output(file_exist_ok, args['output'])

    # Generate shell command from the given parameters
    commands = get_command(
        aligner, mode, input1, input2, output_folder, **aligner_args)

    # Run processes
    run_process(commands)
    
def run_process(
        commands, 
        # thread_num NOTE: here we can set multi-thread option. Currently I 
        # made it disable because I want to make sure that equivalent option 
        # is not specified in aligner_args. 
        ):
    for command in commands:
        os.system(command)

def get_command(aligner, mode, input1='', input2='', output_folder='', **aligner_args):
    if mode == 'pw':
        args_str = get_string_of_aligner_args(aligner_args)
        return [f'{aligner} {args_str}']

    if mode == 'pw_batch':
        if input1 == '':
            raise ValueError('Empty input.')
        if input2 == '':
            raise ValueError('Empty input.')

        input_item_list = make_pair_fasta(input1, input2, output_folder)
        args_str = get_string_of_aligner_args(aligner_args)
        commands = []

        for input_item in input_item_list:
            input_file = input_item + '.fa'
            output_file = input_item + '.aln'

            if aligner == 'clustalo':
                command = f'{aligner} -i {input_file} -o {output_file} {args_str}'
                
            elif aligner == 'muscle':
                command = f'{aligner} -align {input_file} -output {output_file} {args_str}'

            commands.append(command)

    return commands

def get_string_of_aligner_args(aligner_args):
    args_str = ' '.join([
        '{k} {v}'
        for k, v in aligner_args.items()
        if k not in {'input', 'output'}
    ])

    return args_str

def make_pair_fasta(input1, input2, output_folder):
    fasta1 = read_2D_list(input1, '>', join_value_lines=True, skip_empty_lines=True)
    fasta2 = read_2D_list(input2, '>', join_value_lines=True, skip_empty_lines=True)
    fasta_item_list = []

    for key1, value1 in fasta1.items():
        key1_short = key1.split()[0]

        for key2, value2 in fasta2.items():
            key2_short = key2.split()[0]

            out_path = os.path.join(output_folder, f'{key1_short}_{key2_short}.fa')
            if os.path.isfile(out_path):
                raise FileExistsError(out_path)

            out_fasta = f'>{key1}\n{value1}\n>{key2}\n{value2}'
            with open(out_path, 'w') as f:
                print(out_fasta, file=f)

            fasta_item_list.append(out_path[:-3])

    return fasta_item_list

# ============= File Reader [start] ============= #

def read_2D_list(
        file_path       : str, 
        item_divisor    : str = '>',
        comments        : List[str] = ['/*', '#'], 
        value_parser    : Callable[[str], Any] = None, 
        key_parser      : Callable[[str], Any] = None, 
        join_value_lines: bool = False,
        skip_empty_lines: bool = True,
        # cut_inline_comment = False # TODO: Implement in the future
        ) -> Dict[str, Union[str, List[Any]]]:
    """Read a plain-text file containing 2D list data. 

    This function requires a string specifies a division of items (item_divisor). 
    itemnum can be included. 

    Parameters
    ----------
    file_path: str
        Path to input file.
    item_divisor: str, optional (default: '>')
        String specifies a division of items. If a line starts with item_divisor, 
        then the line is considered a start of a new item, therefore the last 
        item is stored as a separate element in an output dictionary. 
    comments: List[str], optional (default: ['/*', '#'])
        Strings that indicate comment lines.
    apply_func: Callable[[str], Any], optional (default: text.do_nothing)
        Function applied to each of value lines.
    join_value_lines: bool, optional (default: False)
        Whether to combine elements of a value into one string. For example, 
        give True when reading a multi-FASTA file. 
    skip_empty_lines: bool, optional (default: True)
        Whether to skip empty lines. 

    Return
    ------
    Dict[str, Union[str, List[Any]]]
        Dictionary of items. A string following to item_divisor is a key and 
        lines until the next item_divisor is stored as value. 
    """
    
    # Initialize a list that will be returned from this function
    items: Dict[str, List[str]] = {}
    key = ''
    value = []
    exp_itemnum = None

    if value_parser == None:
        value_parser = do_nothing
    if key_parser == None:
        key_parser = do_nothing

    # Check if item_divisor is a string and is not empty. 
    check_item_divisor(item_divisor)

    # Open the input file by a read mode
    with open(file_path, 'r') as f:
        # For each line
        for l in f:
            # Remove empty characters (e.g., space, tab or next line) on both 
            # left and right ends
            line = l.strip()

            # If a line is empty
            if line == '':
                # If skip_empty_lines option is True
                if skip_empty_lines:
                    # Go to the next line
                    continue

            # If a line starts with 'itemnum:'
            if line.startswith('itemnum:'):
                # Get the expected number of items
                exp_itemnum = get_itemnum(line)
                # Go to the next line
                continue

            # Check if this a comment line
            full_line_comment = [
                comment # Character indicating that this is a comment line 
                for comment in comments if line.startswith(comment)
            ]
            # If a line starts with one of the comment, 
            if len(full_line_comment) > 0:
                # Go to the next line
                continue

            # If a line starts with the item_divisor, 
            if line.startswith(item_divisor):
                # If this is not empty (meaning that this is not the first item)
                if key != '':
                    # Assign key and value to the items dictionary
                    if join_value_lines:
                        items[key] = ''.join(value)
                    else:
                        items[key] = value

                key = key_parser(line.split(item_divisor)[1])
                assert key != '', 'Empty key for 2D list is not supported.'
                value = []

            else:
                value.append(value_parser(line))
    
    if key != '':
        # Assign key and value to the items dictionary
        if join_value_lines:
            items[key] = ''.join(value)
        else:
            items[key] = value

    if exp_itemnum != None:
        # Raise Assertion error if the observed number of items is different 
        # from the expected. 
        check_item_number(len(items), exp_itemnum)

    return items

def do_nothing(input_str: str) -> str:
    """Returns an input argument as it is. This can be used as a default value 
    of a Callable object. 
    """
    return input_str

def get_itemnum(itemnum_str: str) -> int:
    """Reads itemnum string and returns integer. """
    return int(itemnum_str.split('itemnum:')[1].strip())

# ============= File Reader [end] ============= #

# ============= Parameter Checker [start] ============= #

def check_item_number(obs_itemnum, exp_itemnum):
    msg = f'Wrong number of items found: {obs_itemnum}. {exp_itemnum} was expected.'
    assert obs_itemnum == exp_itemnum, msg
        
def check_item_divisor(item_divisor):
    assert isinstance(item_divisor, str)
    assert item_divisor != ''

def check_output(output, file_exist_ok):
    """Raise FileExistsError if a given file/folder already exists and 
    file_exist_ok=False; otherwise this function will return True.

    """
    if file_exist_ok:
        return True

    if os.path.isfile(output):
        raise FileExistsError(output)

    if os.path.isdir(output):
        raise FileExistsError(output)

    return True

def check_mode(mode):
    if mode not in MODE:
        msg = f'Unknown mode is found: {mode}. Please choose from {MODE}.'
        raise ValueError(msg)

def check_aligner(aligner):
    if aligner not in ALIGNER:
        msg = f'Unknown aligner is found: {aligner}. Please choose from {ALIGNER}.'
        raise ValueError(msg)

# ============= Parameter Checker [end] ============= #

if __name__ == '__main__':
    file_exist_ok = True
    # print(sys.argv)
    main(
        sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], 
        file_exist_ok, 
        # sys.argv[5:]
    )

