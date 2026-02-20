import os, sys
import random

MOTIF_BED = sys.argv[1]
directory = 'BED'
outdirectory='results'
# for each variant calculate the closest

os.system(f'mkdir -p {outdirectory}')

UTR3_BED = os.path.join(directory, '3UTR.bed')
SNP_UTR3_BED = os.path.join(directory, 'SNP_3UTR.bed')

random_number = random.randint(1,10000)

def motif_3utr(UTR3_BED, MOTIF_BED, outfile):
    command = ['bedtools', 'intersect',
               '-a', UTR3_BED,
               '-b', f'{random_number}.tmp.SNP_3UTR.bed',
               '-u',
               '-s',
               '>', f'{random_number}.tmp.3UTR.bed'
              ]

    os.system(' '.join(command))

    command = ['bedtools', 'intersect',
               '-a', MOTIF_BED,
               '-b', f'{random_number}.tmp.3UTR.bed',
               '-s',
               '>', f'{random_number}.tmp.motif.bed'
              ]

    os.system(' '.join(command))

    command = ['bedtools', 'closest',
               '-a', f'{random_number}.tmp.SNP_3UTR.bed',
               '-b', f'{random_number}.tmp.motif.bed',
               '-D', 'ref',
               '-s',
               '>>', outfile
              ]

    os.system(' '.join(command))

    return

def get_motif_distances(SNP_UTR3_BED, UTR3_BED, MOTIF_BED, outfile='test.tsv'):
    n = 0

    # make sure outfile is empty before writing distances
    with open(outfile, 'w') as fout:
        fout.close()

    with open(SNP_UTR3_BED, 'r') as bed_file:
            for line in bed_file.readlines():
                with open(f'{random_number}.tmp.SNP_3UTR.bed', 'w') as tmp_bed:
                    tmp_bed.write(line)
                    tmp_bed.close()

                motif_3utr(UTR3_BED, MOTIF_BED, outfile)
    print(f'{outfile} made')

outfile = os.path.basename(MOTIF_BED).replace('bed', 'results.tsv')
outfile = os.path.join(outdirectory, outfile)
get_motif_distances(SNP_UTR3_BED, UTR3_BED, MOTIF_BED, outfile=outfile)

os.system(f'rm -f {random_number}.tmp.*')
