def parse_function(file_name):
    with open(file_name, 'r') as file:
        return eval(f'lambda x: {file.readline()}')