### Quick Install Guide Meant for First Time Usage Only 

```bash 
# 1 Clone the repo:
git clone https://github.com/your-org/microseq_tests.git
cd microseq_tests 

# 2 Create the Conda env located in env yaml 
conda env create -f config/environment.yml 
conda activate MicroSeq 

# 3 Install & download databases 
pip install -e .  # This sets the current version of MicroSeq 
microseq-setup    # downloads 16S Databases GG2, SILVA, NCBI16S, updates ~/.bashrc 
source ~/.bashrc  # reloads BLASTDB exports (or you can open a new shell next time?) 

# Sanity check you can run to make sure databses work with blasting  

printf ">demo\nACGT\n" > tiny.fasta
for d in gg2 silva ncbi16s; do
  microseq blast -i tiny.fasta -d $d -o /tmp/$d.tsv \
                 --identity 50 --qcov 10 \
      && echo "[✓] $d OK"
done
rm tiny.fasta /tmp/*.tsv
``` 












MicroSeq's GUI uses PySide6 / Qt 6, licensed under the LGPL v3. See LICENSES/Qt_LGPLv3.txt for details. 
