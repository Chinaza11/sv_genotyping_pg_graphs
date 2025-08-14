#! /bin/bash

#SBATCH --partition=Orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --time=10-00:00:00
#SBATCH --mem=200GB
#SBATCH --mail-user=cnnamdi@charlotte.edu
#SBATCH --mail-type=ALL

echo "======================================================"
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID : $SLURM_JOB_NAME"
echo "Node List   : $SLURM_JOB_NODELIST"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo ""

module load parallel/20240422
pjl="parallel_jobs_log"

#=================================================================================
#
#=================================================================================

map_reads2() {
        echo ""
        echo "--------------------------------------"
        read_name=$(basename $1 | awk -F'_' '{print $1}')

        if [[ "$1" == *_1.fastq || "$1" == *_1.fastq.gz ]]; then
                read_name=$(basename $1 | awk -F'_' '{print $1}')
                read_2=$(echo $1 | sed 's/_1/_2/')
        elif [[ "$1" == *_R1.fastq || "$1" == *_R1.fastq.gz ]]; then
                read_name=$(basename $1 | awk -F'.' '{print $1}' | sed 's/_R1//')
                read_2=$(echo "$1" | sed 's/_R1/_R2/')
        else
                echo "Unrecognized read 1 format: $1" >&2
                return 1
        fi

        echo $read_name
        echo $1
        echo $read_2
	echo $2
        echo "--------------------------------------"
}


#reads_path="/projects/cooper_research1/Wild_Sorghum_WGS"
#map_reads2 $reads_path/Grif16309_1.fastq "mappings"

#reads_path="/projects/cooper_research1/TERRA_RAW"
#map_reads2 $reads_path/PI_673843_TCCGCGAA_R1.fastq.gz

#=================================================================================
#
#=================================================================================

vg_call() {
        echo ""
        echo "--------------------------------------"
        sample=$(basename $1 | awk -F'.' '{print $1}')

        echo $sample
        echo $1
	echo $2
	echo $2/$sample.vcf
        echo "--------------------------------------"

}

export -f vg_call

#find mappings/*gam | parallel -j 14 --joblog $pjl/sandbox_vg_call.txt vg_call {} "mappings"

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

for file in vg_call_wild_sorghum_wgs/*.filtered.vcf;
do
        echo ""
        echo "--------------------------------------"
        sample=$(basename $file | awk -F'.' '{print $1}')

        echo $sample
	grep -E "0/1|1/0" $file | wc -l
done

COMMENT

# slurm-4761593.out

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

for file in vg_call_wild_sorghum_wgs/*.filtered.vcf;
do
        echo ""
        echo "--------------------------------------"
        sample=$(basename $file | awk -F'.' '{print $1}')

        echo $sample
        grep -E "1/2" $file | wc -l
done

COMMENT

# slurm-4761613.out

#=================================================================================
#
#=================================================================================

#files=$(ls vg_call_wild_sorghum_wgs/*.filtered.vcf | paste -sd ' ' -)
#echo $files

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

awk '!/^#/ 
	{ for (i=10; i<=NF; i++) 
		{ split($i, a, ":"); if (a[1] != "./." && a[1] != ".") print a[1] 
		} 
	}' vg_call_wild_sorghum_wgs/merged.vcf | sort -u

COMMENT

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

grep -v '^#' vg_call_terra_raw/*.filtered.vcf | \
cut -f10- | \
tr '\t' '\n' | \
cut -d ':' -f1 | \
sort | \
uniq | \
paste -sd ',' -

COMMENT

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

awk 'BEGIN {i=0}
        /^#/ {print; next} {
        i++; $2 = $2 "_" i; print
        }' OFS='\t' vg_call_terra_raw/merged.filtered.vcf > vg_call_terra_raw/merged.filtered2.vcf

COMMENT

#=================================================================================
# TERRA_RAW: deduplicating
#=================================================================================

: << 'COMMENT'

awk '
  BEGIN { OFS = "\t" }
  /^#/ { print; next }  # Print header lines unchanged
  {
    key = $1 ":" $2      # Combine CHROM and POS
    if (!(key in seen)) {
      print
      seen[key] = 1
    }
  }
' vg_call_terra_raw/merged.filtered.vcf > vg_call_terra_raw/merged.filtered3.vcf

COMMENT

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

module load plink/2.00

plink2 \
        --vcf vg_call_terra_raw/merged.filtered.vcf \
        --chr-set 20 no-xy no-mt \
        --threads 30 \
	--pca \
	--geno 0.5 \
        --out sandbox/plink

COMMENT

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env2

read_1="/projects/cooper_research1/TERRA_RAW/PI145626_TCTCGCGC-TATAGCCT_R1.fastq.gz"
read_2="/projects/cooper_research1/TERRA_RAW/PI145626_TCTCGCGC-TATAGCCT_R2.fastq.gz"
sample_name="PI145626_TCTCGCGC-TATAGCCT"

vg giraffe -Z sorghum-pg/sorghum-pg.d2.gbz \
	-m sorghum-pg/sorghum-pg.d2.shortread.withzip.min \
	-z sorghum-pg/sorghum-pg.d2.shortread.zipcodes \
	-d sorghum-pg/sorghum-pg.d2.dist \
	--progress \
	-t 30 \
	--log-reads \
	-f $read_1 -f $read_2 >sandbox/$sample_name.mapped.gam

COMMENT

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

cd  ~/packages/varsim/tests/quickstart_test/ 
./quickstart.sh

# slurm-5130429.out

COMMENT

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

cd /projects/cooper_research2/chinaza/packages/varsim/tests/quickstart_test/
./quickstart.sh

# slurm-5135830.out

COMMENT

#=================================================================================
#
#=================================================================================

#: << 'COMMENT'

cd /projects/cooper_research2/chinaza/packages/varsim/tests/quickstart_test_2/
./quickstart.sh

# slurm-___.out

#COMMENT

#=================================================================================
#
#=================================================================================

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
