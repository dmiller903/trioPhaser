FROM lifebitai/shapeit4
  
LABEL \
  author="Dustin Miller" \
  description="Tools needed to phase a trio with SHAPEIT4 and Mendelian inheritance logic" \
  maintainer="dustin_miller@byu.edu"

RUN \
# Install dependencies and other useful programs
apt-get update \
  && apt-get install -y --no-install-recommends wget \
  vim=2:8.0.1453-1ubuntu1.4 \
  unzip=6.0-21ubuntu1.1 \
  zip=3.0-11build1 \
  tar=1.29b-2ubuntu0.2 \
  git=1:2.17.1-1ubuntu0.8 \
  # Install miniconda2 and config channels
  && wget --no-check-certificate https://repo.anaconda.com/miniconda/Miniconda3-py38_4.9.2-Linux-x86_64.sh \
  && bash Miniconda3-py38_4.9.2-Linux-x86_64.sh -b -p \
  && export PATH="/root/miniconda3/bin:$PATH" \
  && /root/miniconda3/bin/conda init \
  && conda config --append channels bioconda \
  # Install GATK4
  && conda install -c bioconda gatk4=4.0.5.1 \
  # Install bcftools
  && conda install -c bioconda bcftools=1.9 \
  # Install samtools
  && conda install -c bioconda samtools=1.9 \
  # Download reference files
  && wget --no-check-certificate https://files.osf.io/v1/resources/3znuj/providers/osfstorage/5d9f54d2a7bc73000ee99fd6/?zip= -O /fasta_references.zip \
  # Clone trioPhaser from github and place python scripts at root
  && git config --global http.sslverify "false" \
  && git clone https://github.com/dmiller903/trioPhaser.git \
  && cp trioPhaser/trio_phaser.py

ENV \
  PATH=${PATH}:/root/miniconda3/bin/