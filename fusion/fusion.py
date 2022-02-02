import argparse
import glob
import os
import logging
import shutil
import subprocess
from pathlib import Path

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

parser = argparse.ArgumentParser()

parser.add_argument('sample', type=str,
    help='Sample id')

parser.add_argument('fq_1', type=str,
    help='RNA-seq fastq 1')

parser.add_argument('fq_2', type=str,
    help='RNA-seq fastq 2')

parser.add_argument('--cpu', type=str, default=1,
    help='Number of CPUs to run on')

parser.add_argument('--genome-lib-dir', type=str,
    help='GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir directory location')

parser.add_argument('--genome-db', type=str,
    help='ericscript_db_homosapiens_ensembl84 directory location')

parser.add_argument('--bwts', type=str,
    help='bwts location')

parser.add_argument('--integrate-fasta', type=str,
    help='fasta used during integration step')

parser.add_argument('--integrate-annotations', type=str,
    help='location of annot.ensembl.GRCh38.txt used during integration step')

parser.add_argument('--combine-call-script', type=str,
    help='location of combine_call.pl')

parser.add_argument('--filter-script', type=str,
    help='location of filter.pl')

args = parser.parse_args()


FUSION_OUT = 'STAR_FUSION'
ERICSCRIPT_OUT = 'ERICSCRIPT'
ERICSCRIPT = 'ericscript.pl'
STAR_OUT = 'STAR'
INTEGRATE_OUT = 'INTEGRATE'
MERGE_FUSIONS_OUT = 'Merged_Fusions'


def run_star_fusion():
    logging.info('running star fusion')
    cmd = f'STAR-Fusion --left_fq {args.fq_1} --right_fq {args.fq_2} --CPU {args.cpu} --examine_coding_effect -O {FUSION_OUT} --genome_lib_dir {args.genome_lib_dir} --verbose_level 2'
    logging.info(f'executing command: {cmd}')
    output = subprocess.check_output(cmd, shell=True)
    logging.info('tool output:')
    logging.info(output)

    logging.info('cleaning large fusion files')
    if os.path.exists(f'{FUSION_OUT}/Aligned.out.bam'):
        os.remove(f'{FUSION_OUT}/Aligned.out.bam')

    dirs = glob.glob(f'{FUSION_OUT}/_*')
    for d in dirs:
        if os.path.exists(d):
            shutil.rmtree(d)

    if os.path.exists(f'{FUSION_OUT}/star-fusion.preliminary'):
        shutil.rmtree(f'{FUSION_OUT}/star-fusion.preliminary')


def run_ericscript():
    logging.info('running ericscript')
    cmd = f'{ERICSCRIPT} -o {ERICSCRIPT_OUT} --remove -ntrim 0 --refid homo_sapiens -db {args.genome_db} -p {args.cpu} -name {args.sample} {args.fq_1} {args.fq_2}'
    logging.info(f'executing command: {cmd}')
    output = subprocess.check_output(cmd, shell=True)
    logging.info('tool output:')
    logging.info(output)


def run_star():
    logging.info('running star')
    Path(STAR_OUT).mkdir(parents=True, exist_ok=True)
    star_genome_dir = os.path.join(args.genome_lib_dir, 'ref_genome.fa.star.idx')
    cmd = f'STAR --runThreadN 80 --genomeDir {star_genome_dir} --readFilesCommand zcat --readFilesIn {args.fq_1} {args.fq_2} --outSAMtype BAM Unsorted --outFileNamePrefix {STAR_OUT}/ --chimSegmentMin 18 --chimOutType WithinBAM --outSAMunmapped Within'
    logging.info(f'executing command: {cmd}')
    output = subprocess.check_output(cmd, shell=True)
    logging.info('tool output:')
    logging.info(output)

    logging.info('sorting bam')
    cmd = f'samtools sort {STAR_OUT}/Aligned.out.bam {STAR_OUT}/Aligned.out.sorted'
    logging.info(f'executing command: {cmd}')
    output = subprocess.check_output(cmd, shell=True)

    logging.info('indexing bam')
    cmd = f'samtools index {STAR_OUT}/Aligned.out.sorted.bam'
    logging.info(f'executing command: {cmd}')
    output = subprocess.check_output(cmd, shell=True)


def run_integrate():
    logging.info('running integrate')
    Path(INTEGRATE_OUT).mkdir(parents=True, exist_ok=True)
    cmd = f'Integrate fusion -reads {INTEGRATE_OUT}/reads.txt -sum {INTEGRATE_OUT}/summary.tsv -ex {INTEGRATE_OUT}/exons.tsv -bk {INTEGRATE_OUT}/breakpoints.tsv -vcf {INTEGRATE_OUT}/bk_sv.vcf -bedpe {INTEGRATE_OUT}/fusions.bedpe {INTEGRATE_OUT} {args.integrate_annotations} {args.bwts} {STAR_OUT}/Aligned.out.sorted.bam {STAR_OUT}/Aligned.out.sorted.bam'
    logging.info(f'executing command: {cmd}')
    output = subprocess.check_output(cmd, shell=True)
    logging.info('tool output:')
    logging.info(output)

    logging.info('clean STAR outputs')
    for fp in [f'{STAR_OUT}/Aligned.out.bam',
               f'{STAR_OUT}/Aligned.out.sorted.bam',
               f'{STAR_OUT}/Aligned.out.sorted.bam.bai']:
        if os.path.exists(fp):
            os.remove(fp)

    for fp in [f'{STAR_OUT}/star_STARtmp']:
        if os.path.exists(fp):
            shutil.rmtree(fp)


def run_merge_fusion():
    logging.info('running merge fusions')
    Path(MERGE_FUSIONS_OUT).mkdir(parents=True, exist_ok=True)
    cmd = f'perl {args.combine_call_script} {args.sample} {FUSION_OUT}/star-fusion.fusion_predictions.abridged.coding_effect.tsv {ERICSCRIPT_OUT}/{args.sample}.results.total.tsv {INTEGRATE_OUT}/summary.tsv {INTEGRATE_OUT}/breakpoints.tsv {MERGE_FUSIONS_OUT}'
    logging.info(f'executing command: {cmd}')
    output = subprocess.check_output(cmd, shell=True)
    logging.info('tool output:')
    logging.info(output)

    logging.info('running merge fusions filtering')
    cmd = f'perl {args.filter_script} {args.merge_fusions_out}'
    logging.info(f'executing command: {cmd}')
    output = subprocess.check_output(cmd, shell=True)


def main():
    logging.info('starting fusion pipeline')
    run_star_fusion()
    run_ericscript()
    run_star()
    run_integrate()
    run_merge_fusion()
    logging.info('fusion pipeline finished')


if __name__ == '__main__':
    main()
