
# csv.field_size_limit(1000000000)


def process_file(input_file, output_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    string_N80 = 'N' * 80  + '\n'
    # Remover linhas que começam com '>'
    # lines = [line for line in lines if line.startswith('>') line = string_N80 + string_N80 else line ]

    output_lines = []
    for line in lines:
        line = line.strip()
        if len(output_lines) == 0:            
            output_lines.append(line + '\n')
        elif line.startswith('>'):
            output_lines.append(string_N80)
            output_lines.append(string_N80)
        else:
            line = line.upper()
            if len(line) < 80:
                line = line + 'N' * (80-len(line))
            output_lines.append(line + '\n')

    lines = output_lines
    output_lines = []
    count_N80_line = 0

    line_count = 0
    for line in lines:
        line_count += 1
        if line == string_N80:
            count_N80_line += 1
            if count_N80_line > 2:
                continue
            # else:
            #     print(f"Entered N line but didnt skip: {line_count}")
        else:
            count_N80_line = 0
        output_lines.append(line)

    # Escrever resultado no arquivo de saída
    with open(output_file, 'w') as file:
        file.writelines(output_lines)

# Exemplo de uso
input_file = 'data/human_genome.fna'
output_file = 'data/modified_human_genome.fna'
process_file(input_file, output_file)
