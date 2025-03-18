# Discover_DNMs
Leveraging Larger Family to detect DNMs


python3 Discover_DNMs.py -h
usage: deNOV_M.py [-h] [--depth_min DEPTH_MIN] [--depth_max DEPTH_MAX] [--input_file INPUT_FILE]
                  [--child [CHILD [CHILD ...]]] [--output_file OUTPUT_FILE] [--out_dir OUT_DIR] [--parent1_id PARENT1_ID]
                  [--parent2_id PARENT2_ID]

De Novo Mutations Discovery in the Offspring, Maher ALnajjar, 2025, MATE-GBI-Genetics and Genomics Dep. Dr.Barta's LAB,
Gödöllő, HUNGARY

optional arguments:
  -h, --help            show this help message and exit
  --depth_min DEPTH_MIN, -d_min DEPTH_MIN
                        Define the MINIMUM Filtering Depth, default=15
  --depth_max DEPTH_MAX, -d_max DEPTH_MAX
                        Define the MAXIMUM Filtering Depth, default=45
  --input_file INPUT_FILE, -i INPUT_FILE
                        The file is a standard VCF file coming from GATK for example
  --child [CHILD [CHILD ...]], -c [CHILD [CHILD ...]]
                        list of children sample id(s) separated by spaces, otherwise all samples in the vcf will be
                        considered
  --output_file OUTPUT_FILE, -o OUTPUT_FILE
                        output file
  --out_dir OUT_DIR, -o_dir OUT_DIR
                        indicate an output directory
  --parent1_id PARENT1_ID, -p1 PARENT1_ID
                        One of the parents id
  --parent2_id PARENT2_ID, -p2 PARENT2_ID
                        The other parent id
