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
- python
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
- id: fq_2
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
- id: integrate_executable
  inputBinding:
    position: '0'
    prefix: --integrate-executable
  type: File
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
- id: filter_database
  inputBinding:
    position: '0'
    prefix: --filter-database
  type: Directory
- id: fusion_annotator_dir
  inputBinding:
    position: '0'
    prefix: --fusion-annotator-dir
  type: Directory
- default: 1
  id: cpu
  inputBinding:
    position: '0'
    prefix: --cpu
  type: int?
- default: /miniconda/envs/fusion/bin:$PATH
  id: environ_PATH
  type: string?
- default: /miniconda/envs/fusion/lib/:$LD_LIBRARY_PATH
  id: environ_LD_LIBRARY_PATH
  type: string?
label: fusion
outputs:
- id: filtered_fusions
  outputBinding:
    glob: Merged_Fusions/Filtered*.tsv
  type: File
- id: total_fusions
  outputBinding:
    glob: Merged_Fusions/Total*.tsv
  type: File
requirements:
- class: DockerRequirement
  dockerPull: estorrs/pecgs_fusion:0.0.2
- class: ResourceRequirement
  coresMin: $(inputs.cpu)
  ramMin: 50000
- class: EnvVarRequirement
  envDef:
    LD_LIBRARY_PATH: $(inputs.environ_LD_LIBRARY_PATH)
    PATH: $(inputs.environ_PATH)
