$namespaces:
  sbg: https://www.sevenbridges.com/
arguments:
- position: 0
  prefix: --combine-call-script
  valueFrom: /pecgs-fusion/fusion/combine_call.pl
- position: 0
  prefix: --filter-script
  valueFrom: /pecgs-fusion/fusion/filter.pl
baseCommand:
- /pecgs-fusion/fusion/fusion.py
class: CommandLineTool
cwlVersion: v1.0
id: fusion
inputs:
- id: sample
  inputBinding:
    position: '1'
  type: string
- id: fq_1
  inputBinding:
    position: '2'
  type: File
- id: fq2
  inputBinding:
    position: '3'
  type: File
- id: genome_lib_dir
  inputBinding:
    position: '0'
    prefix: --genome-lib-dir
  type: Directory
- id: genome_db
  inputBinding:
    position: '0'
    prefix: --genome-db
  type: Directory
- id: bwts
  inputBinding:
    position: '0'
    prefix: --bwts
  type: Directory
- id: integrate_fasta
  inputBinding:
    position: '0'
    prefix: --integrate-fasta
  type: File
- id: integrate_annotations
  inputBinding:
    position: '0'
    prefix: --integrate-annotations
  type: File
- id: cpu
  inputBinding:
    position: '0'
    prefix: -cpu
  type: string?
label: fusion
outputs:
- id: output_summary
  outputBinding:
    glob: output
  type: File
- id: output_dis
  outputBinding:
    glob: output_dis
  type: File
- id: output_germline
  outputBinding:
    glob: output_germline
  type: File
- id: output_somatic
  outputBinding:
    glob: output_somatic
  type: File
requirements:
- class: DockerRequirement
  dockerPull: estorrs/pecgs_fusion:0.0.1
- class: ResourceRequirement
  ramMin: 50000
