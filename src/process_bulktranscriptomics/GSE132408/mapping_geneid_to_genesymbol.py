with open("path_to_geneid_mapping.tab") as mapping_table:
    mapping_table.readline()
    geneid_genesymbol = {}
    for line in mapping_table:
        line = line.strip().split("\t")
        line[1] = line[1].split(" ")[0]
        if line[2] not in geneid_genesymbol:
            geneid_genesymbol[line[2]] = set()
        geneid_genesymbol[line[2]].add(line[1])

with open("path_to_geneid_list") as geneid_list:
    with open("output_file", "w") as output_file:
        for line in geneid_list:
            line = line.strip().split("\t")
            if line[0] in geneid_genesymbol:
                for gs in geneid_genesymbol[line[0]]:
                    output_file.write(gs + "\t"+ "\t".join(line) + "\n")
