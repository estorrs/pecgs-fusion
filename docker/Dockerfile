FROM ubuntu:18.04

RUN apt-get update && apt-get install -y vim wget curl git

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.3-Linux-x86_64.sh -O ~/miniconda.sh
RUN bash ~/miniconda.sh -b -p /miniconda
ENV PATH="/miniconda/bin:$PATH"

# get env file seperately so it doesn't reinstall every time
COPY ./env.yaml /pecgs-fusion/env.yaml
RUN conda env create --file /pecgs-fusion/env.yaml 

COPY . /pecgs-fusion

ENV PATH="/miniconda/envs/fusion/bin:$PATH"
ENV LD_LIBRARY_PATH="/miniconda/envs/fusion/lib/:$LD_LIBRARY_PATH"

# Needs this to work on compute1 LSF bc persmissions are weird with how docker is launched
RUN chmod -R +777 /miniconda/envs/fusion/share/ericscript-0.5.5-5/
# need to make this executable for the merge script
RUN chmod +x /miniconda/envs/fusion/lib/STAR-Fusion/FusionAnnotator/util/fusion_to_coding_region_effect.pl

CMD /bin/bash
