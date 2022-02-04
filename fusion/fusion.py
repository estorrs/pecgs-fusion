import argparse
import glob
import os
import logging
import shutil
import subprocess
# we are python2 :(
# from pathlib import Path

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

parser.add_argument('--filter-database', type=str,
    help='Filter database location. ....../FilterDatabase')

parser.add_argument('--fusion-annotator-dir', type=str,
    help='fusion annotation location. ....../FusionAnnotater')

parser.add_argument('--integrate-executable', type=str,
    help='Location of Integrate executable')

parser.add_argument('--integrate-fasta', type=str,
    help='fasta used during integration step')

parser.add_argument('--integrate-annotations', type=str,
    help='location of annot.ensembl.GRCh38.txt used during integration step')

parser.add_argument('--combine-call-script', type=str,
    help='location of combine_call.pl')

parser.add_argument('--filter-script', type=str,
    help='location of filter.pl')

print('here')

args = parser.parse_args()

print('here')


FUSION_OUT = 'STAR_FUSION'
ERICSCRIPT_OUT = 'ERICSCRIPT'
ERICSCRIPT = 'ericscript.pl'
STAR_OUT = 'STAR'
INTEGRATE_OUT = 'INTEGRATE'
MERGE_FUSIONS_OUT = 'Merged_Fusions'


def run_star_fusion():
    logging.info('running star fusion')
    # ugh cant use f strings bc of python2
    #cmd = f'STAR-Fusion --left_fq {args.fq_1} --right_fq {args.fq_2} --CPU {args.cpu} --examine_coding_effect -O {FUSION_OUT} --genome_lib_dir {args.genome_lib_dir} --verbose_level 2'
    cmd = 'STAR-Fusion --left_fq {fq_1} --right_fq {fq_2} --CPU {cpu} --examine_coding_effect -O {fo} --genome_lib_dir {genome_lib_dir} --verbose_level 2'.format(
        fq_1=args.fq_1, fq_2=args.fq_2, cpu=args.cpu, fo=FUSION_OUT, genome_lib_dir=args.genome_lib_dir, )
    logging.info('executing command: {cmd}'.format(cmd=cmd))
    output = subprocess.check_output(cmd, shell=True)
    logging.info('tool output:')
    logging.info(output)

    logging.info('cleaning large fusion files')
    if os.path.exists('{fo}/Aligned.out.bam'.format(fo=FUSION_OUT)):
        os.remove('{fo}/Aligned.out.bam'.format(fo=FUSION_OUT))

    dirs = glob.glob('{fo}/_*'.format(fo=FUSION_OUT))
    for d in dirs:
        if os.path.exists(d):
            shutil.rmtree(d)

    if os.path.exists('{fo}/star-fusion.preliminary'.format(fo=FUSION_OUT)):
        shutil.rmtree('{fo}/star-fusion.preliminary'.format(fo=FUSION_OUT))


def run_ericscript():
    logging.info('running ericscript')
    # cmd = f'{ERICSCRIPT} -o {ERICSCRIPT_OUT} --remove -ntrim 0 --refid homo_sapiens -db {args.genome_db} -p {args.cpu} -name {args.sample} {args.fq_1} {args.fq_2}'
    cmd = '{es} -o {eo} --remove -ntrim 0 --refid homo_sapiens -db {genome_db} -p {cpu} -name {sample} {fq_1} {fq_2}'.format(
        es=ERICSCRIPT, eo=ERICSCRIPT_OUT, genome_db=args.genome_db, cpu=args.cpu, sample=args.sample, fq_1=args.fq_1, fq_2=args.fq_2)
    logging.info('executing command: {cmd}'.format(cmd=cmd))
    output = subprocess.check_output(cmd, shell=True)
    logging.info('tool output:')
    logging.info(output)


def run_star():
    logging.info('running star')
    # Path(STAR_OUT).mkdir(parents=True, exist_ok=True)
    if not os.path.exists(STAR_OUT):
        os.mkdir(STAR_OUT)
    star_genome_dir = os.path.join(args.genome_lib_dir, 'ref_genome.fa.star.idx')
    # cmd = f'STAR --runThreadN 80 --genomeDir {star_genome_dir} --readFilesCommand zcat --readFilesIn {args.fq_1} {args.fq_2} --outSAMtype BAM Unsorted --outFileNamePrefix {STAR_OUT}/ --chimSegmentMin 18 --chimOutType WithinBAM --outSAMunmapped Within'
    cmd = 'STAR --runThreadN 80 --genomeDir {star_genome_dir} --readFilesCommand zcat --readFilesIn {fq_1} {fq_2} --outSAMtype BAM Unsorted --outFileNamePrefix {so}/ --chimSegmentMin 18 --chimOutType WithinBAM --outSAMunmapped Within'.format(
        star_genome_dir=star_genome_dir, fq_1=args.fq_1, fq_2=args.fq_2, so=STAR_OUT)
    logging.info('executing command: {cmd}'.format(cmd=cmd))
    output = subprocess.check_output(cmd, shell=True)
    logging.info('tool output:')
    logging.info(output)

    logging.info('sorting bam')
    cmd = 'samtools sort {so}/Aligned.out.bam {so}/Aligned.out.sorted'.format(
        so=STAR_OUT)
    logging.info('executing command: {cmd}'.format(cmd=cmd))
    output = subprocess.check_output(cmd, shell=True)

    logging.info('indexing bam')
    cmd = 'samtools index {so}/Aligned.out.sorted.bam'.format(
        so=STAR_OUT)
    logging.info('executing command: {cmd}'.format(cmd=cmd))
    output = subprocess.check_output(cmd, shell=True)


def run_integrate():
    logging.info('running integrate')
    # Path(INTEGRATE_OUT).mkdir(parents=True, exist_ok=True)
    if not os.path.exists(INTEGRATE_OUT):
        os.mkdir(INTEGRATE_OUT)
    # cmd = f'Integrate fusion -reads {INTEGRATE_OUT}/reads.txt -sum {INTEGRATE_OUT}/summary.tsv -ex {INTEGRATE_OUT}/exons.tsv -bk {INTEGRATE_OUT}/breakpoints.tsv -vcf {INTEGRATE_OUT}/bk_sv.vcf -bedpe {INTEGRATE_OUT}/fusions.bedpe {INTEGRATE_OUT} {args.integrate_annotations} {args.bwts} {STAR_OUT}/Aligned.out.sorted.bam {STAR_OUT}/Aligned.out.sorted.bam'
    cmd = '{ie} fusion -reads {io}/reads.txt -sum {io}/summary.tsv -ex {io}/exons.tsv -bk {io}/breakpoints.tsv -vcf {io}/bk_sv.vcf -bedpe {io}/fusions.bedpe {integrate_fasta} {integrate_annotations} {bwts} {so}/Aligned.out.sorted.bam {so}/Aligned.out.sorted.bam'.format(
        ie=args.integrate_executable, io=INTEGRATE_OUT, integrate_fasta=args.integrate_fasta, integrate_annotations=args.integrate_annotations, bwts=args.bwts, so=STAR_OUT)
    logging.info('executing command: {cmd}'.format(cmd=cmd))
    output = subprocess.check_output(cmd, shell=True)
    logging.info('tool output:')
    logging.info(output)

    logging.info('clean STAR outputs')
    for fp in ['{so}/Aligned.out.bam'.format(so=STAR_OUT),
               '{so}/Aligned.out.sorted.bam'.format(so=STAR_OUT),
               '{so}/Aligned.out.sorted.bam.bai'.format(so=STAR_OUT)]:
        if os.path.exists(fp):
            os.remove(fp)

    for fp in ['{so}/star_STARtmp'.format(so=STAR_OUT)]:
        if os.path.exists(fp):
            shutil.rmtree(fp)


def run_merge_fusion():
    logging.info('running merge fusions')
    # Path(MERGE_FUSIONS_OUT).mkdir(parents=True, exist_ok=True)
    if not os.path.exists(MERGE_FUSIONS_OUT):
        os.mkdir(MERGE_FUSIONS_OUT)
    # cmd = f'perl {args.combine_call_script} {args.sample} {FUSION_OUT}/star-fusion.fusion_predictions.abridged.coding_effect.tsv {ERICSCRIPT_OUT}/{args.sample}.results.total.tsv {INTEGRATE_OUT}/summary.tsv {INTEGRATE_OUT}/breakpoints.tsv {MERGE_FUSIONS_OUT}'
    cmd = 'perl {combine_call_script} {sample} {fo}/star-fusion.fusion_predictions.abridged.coding_effect.tsv {eo}/{sample}.results.total.tsv {io}/summary.tsv {io}/breakpoints.tsv {mf} {fd} {sb} {fa}'.format(
        combine_call_script=args.combine_call_script, sample=args.sample, fo=FUSION_OUT, eo=ERICSCRIPT_OUT, io=INTEGRATE_OUT, mf=MERGE_FUSIONS_OUT, fd=args.filter_database, sb=args.genome_lib_dir, fa=args.fusion_annotator_dir)
    logging.info('executing command: {cmd}'.format(cmd=cmd))
    output = subprocess.check_output(cmd, shell=True)
    logging.info('tool output:')
    logging.info(output)

    logging.info('running merge fusions filtering')
    cmd = 'perl {filter_script} {mf} {sample} {df}'.format(
        filter_script=args.filter_script, mf=MERGE_FUSIONS_OUT, sample=args.sample, fd=args.filter_database)
    logging.info('executing command: {cmd}'.format(cmd=cmd))
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
