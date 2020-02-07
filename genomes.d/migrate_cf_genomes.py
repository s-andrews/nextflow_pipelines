#!/bi/apps/python/3.7.3/bin/python3
from glob import glob
import re


def process_file (file):
    print(f"Processing {file}")

    with open(file) as fh:
        name = None
        species = None
        attributes = {}
        for line in fh:
            line = line[11:]
            sections = re.split("\s+",line)

            if len(sections) < 3:
                continue

            if name is None:
                print(f"Got name {sections[1]} from {line}")
                name = sections[1]
                species = sections[3]

            attributes[sections[0]] = sections[2]

    outfile = f"{name}.genome"

    with open(outfile,"w") as out:
        out.write(f"name\t{name}\n")
        out.write(f"species\t{species}\n")

        for attribute in attributes.keys():
            out.write(f"{attribute}\t{attributes[attribute]}\n")

def main():
    cf_directory = "/bi/apps/clusterflow/0.6_devel/genomes.d/"

    cf_files = glob(cf_directory+"/*.config")

    for file in cf_files:
        process_file(file)


if __name__ == "__main__":
    main()

