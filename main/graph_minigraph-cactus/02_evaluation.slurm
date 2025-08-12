#! /bin/bash

#SBATCH --partition=Orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=10-00:00:00
#SBATCH --mem=100GB
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
# simulating structural variation
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate survivor

cd evaluation/sandbox/

# SURVIVOR simSV parameter_file

SURVIVOR simSV riouncc.fa parameter_file 0 0 simulated

# SURVIVOR eval simulated.vcf simulated.bed 10 eval_res

COMMENT

# slurm-4875594.out

#=================================================================================
# function for simulating illumina short paired-end reads from a genome
#=================================================================================

simulate_short_reads() {
	echo ""
        echo "--------------------------------------"
        sample=$(basename $1 | awk -F'.' '{print $1}')

        echo $1
        echo $sample
        echo "--------------------------------------"

	art_illumina -ss HS25 -p -sam -i $1 \
		-l 150 -f 45 -m 500 -s 50 \
		-o $2/simulated_reads/$sample
	}

export -f simulate_short_reads

#=================================================================================
# simulate illumina short paired-end reads
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate art

find evaluation/sandbox/ -type f \( -name "riouncc.fa" -o -name "simulated.fasta" \) \
	| parallel -j 2 --joblog $pjl/sandbox_simulate_short_reads.txt simulate_short_reads {} "evaluation/sandbox"

COMMENT

# slurm-4876026.out

#=================================================================================
# pangenome graph construction using minigraph-cactus
#=================================================================================

: << 'COMMENT'

cactus_dir="/projects/cooper_research2/chinaza/packages/cactus/cactus-bin-v2.9.3/"
source "${cactus_dir}venv-cactus-v2.9.3/bin/activate"
#cactus --version

#toil clean jobstore

cactus-pangenome jobstore \
	evaluation/sandbox/pangenome/seq_file.txt \
        --outDir evaluation/sandbox/pangenome \
        --outName sorghum-pg \
        --reference riouncc \
        --vcf --giraffe --gfa --gbz
# slurm-4875865.out

halStats evaluation/sandbox/pangenome/sorghum-pg.full.hal
# slurm-4998990.out

COMMENT

#=================================================================================
# activate conda environment
# activate ./01_graph_minigraph-cactus.slurm script
# create base directory
#=================================================================================

module load anaconda3/2023.09
conda activate vg_env2

source ./01_graph_minigraph-cactus.slurm

base="evaluation/sandbox"
export base

#=================================================================================
# recreate minimizer index (that includes zipcodes) and a new zipcodes file
# map some reads
#=================================================================================

: << 'COMMENT'

map_reads "$base/simulated_reads/riouncc1.fq" "$base/pangenome" "$base/mappings"

COMMENT

# slurm-4999007.out

#=================================================================================
# map other reads
#=================================================================================

: << 'COMMENT'

map_reads2 "$base/simulated_reads/simulated1.fq" "$base/pangenome" "$base/mappings"

COMMENT

# slurm-4999866.out

#=================================================================================
# stats of alignments (gam file)
#=================================================================================

: << 'COMMENT'

find $base/mappings/*gam \
	| parallel -j 2 --joblog $pjl/vg_stats_eval_sandbox.txt vg_stats {}

COMMENT

# slurm-5012567.out

#=================================================================================
# SV calling and genotyping
#=================================================================================

# compute the snarls
#vg snarls $base/pangenome/sorghum-pg.d2.gbz > $base/pangenome/sorghum-pg.d2.snarls
# slurm-5012576.out


: << 'COMMENT'

find $base/mappings/*gam \
	| parallel -j 2 --joblog $pjl/vg_call_eval_sandbox.txt \
	vg_call {} "$base/vg_call" "$base/pangenome"

COMMENT

#slurm-5012581.out

#=================================================================================
# filter VCF file, exclude rows with lowdepth and lowad and rows with SVs less than 50bp
#=================================================================================

: << 'COMMENT'

#rm $base/vg_call/*filtered.vcf

find $base/vg_call/*vcf \
       | parallel -j 2 --joblog $pjl/filter_vcf_eval_sandbox.txt filter_vcf {} "$base/vg_call/"

#slurm-5026737.out

COMMENT

#=================================================================================
# map simulated reads to original pangenome that was constructed without the simulated genome
# genotype mapped reads, filter it
#=================================================================================

#map_reads2 "$base/simulated_reads/simulated1.fq" "sorghum-pg" "$base/mappings/mapped_to_original_pangenome"
# ==> slurm-5026730.out

#vg_call "$base/mappings/mapped_to_original_pangenome/simulated.mapped.gam" \
#	"$base/vg_call/mapped_to_original_pangenome" "sorghum-pg"
# ==> slurm-5027338.out

#rm $base/vg_call/mapped_to_original_pangenome/*.filtered.vcf
#filter_vcf $base/vg_call/mapped_to_original_pangenome/simulated.vcf "$base/vg_call/mapped_to_original_pangenome"

#=================================================================================
#
#=================================================================================

: << 'COMMENT'



COMMENT

#=================================================================================
#
#=================================================================================

extract_fasta_region() {
  local fasta_file=$1
  local chromosome=$2
  local start=$3
  local end=$4

  awk -v target="$chromosome" -v start="$start" -v end="$end" '
  BEGIN {
      pos = 0
      capture = 0
  }
  /^>/ {
      curr_chr = substr($0, 2)
      capture = (curr_chr == target)
      pos = 0
      next
  }
  capture {
      line = $0
      len = length(line)
      if (pos + len >= start) {
          for (i = 1; i <= len; i++) {
              pos++
              if (pos >= start && pos <= end) {
                  seq = seq substr(line, i, 1)
              }
              if (pos > end) {
                  print seq
                  exit
              }
          }
      } else {
          pos += len
      }
  }
  ' "$fasta_file"
}

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

cd $base
{
	awk '$5 == "DEL"' simulated.bed | sort -k1,1 -k2,2n | \
	while read chr1 start chr2 end type; do
		echo ">$chr1:$start-$end"
		extract_fasta_region riouncc.fa "$chr1" "$start" "$end"
	done
} > simulated.deletions.fa

#slurm-5199075.out

COMMENT

#=================================================================================
#
#=================================================================================

process_vcf_variants() {
  local vcf_file="$1"
  local outdir="$2"

  mkdir -p "$outdir"

  if [[ ! -f "$vcf_file" ]]; then
    echo "VCF file not found: $vcf_file"
    return 1
  fi

  awk -v outdir="$outdir" '
       BEGIN { FS="\t" }
       !/^#/ {
           ref_len = length($4)
           alt_len = length($5)

           # Skip multi-allelic ALT or equal-length variants
           if ($5 ~ /,/) next
           if (alt_len == ref_len) next

           if (alt_len > ref_len) {
               type = "INSERTION"
               seq = $5
           } else {
               type = "DELETION"
               seq = $4
           }

           chr = $1
           pos = $2
           out_file = outdir "/" chr "_" pos "_" type ".txt"

           print ">" chr ":" pos ":" type > out_file
           print seq >> out_file
       }' "$vcf_file"
}

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

cd $base
process_vcf_variants fn.vcf "fn"

COMMENT

#=================================================================================
#
#=================================================================================

split_fasta_to_files() {
	local fasta_file="$1"
	local output_dir="$2"

	mkdir -p "$output_dir"

	awk -v outdir="$output_dir" '
		/^>/ {
			if (seq) {
				print seq > filename
				seq = ""
			}
			# Clean the header to create a valid filename
			header = substr($0, 2)
			gsub(/[ \t:]/, "_", header)
			filename = outdir "/" header ".fa"
			print $0 > filename
			next
		}
		{
			seq = seq $0
		}
		END {
			if (seq) {
				print seq > filename
			}
		}
	' "$fasta_file"
}

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

cd $base
split_fasta_to_files simulated.insertions.fa simulated.insertions
split_fasta_to_files simulated.deletions.fa simulated.deletions

COMMENT

#=================================================================================
#
#=================================================================================

variant_alignment() {
	local filepath="$1"
	local filename=$(basename "$filepath")
	local name="${filename%.txt}"

	# Split name by underscore
	IFS='_' read -r chr pos type <<< "$name"

	# Set simulated directory based on type
	if [[ "$type" == "DELETION" ]]; then
		sim_dir="simulated.deletions"
	elif [[ "$type" == "INSERTION" ]]; then
		sim_dir="simulated.insertions"
	else
		echo "Unknown variant type: $type"
		return 1
	fi

	# Create output directory using only the filename (no path)
	mkdir -p "$base/$2/$name"

	# Loop through matching simulated variant files
	for sim_file in "$base/$sim_dir/$chr"_*; do
		if [[ -f "$sim_file" ]]; then
			sim_name=$(basename "$sim_file" .fa)
			output_file="$base/$2/$name/$sim_name.txt"

			# Run Smith-Waterman alignment
			water -asequence "$filepath" \
			      -bsequence "$sim_file" \
			      -gapopen 10 \
			      -gapextend 0.5 \
			      -outfile "$output_file"
		fi
	done

	# Parse alignment output and collect stats
	stats_output="$base/$2/$name/alignment_summary.tsv"
	echo -e "Asequence\tAsequence_len\tBsequence\tBsequence_len\tIdentity_raw\tIdentity_pct\tSimilarity_raw\tSimilarity_pct\tGaps_pct" > "$stats_output"

	for result_file in "$base/$2/$name/"*; do
		if [[ -f "$result_file" ]]; then
			aseq=$(grep 'asequence' "$result_file" | awk -F'/' '{print $NF}')
                        bseq=$(grep 'bsequence' "$result_file" | awk -F'/' '{print $NF}')

			aseq_len=$(grep 'asequence' $result_file | awk -F'asequence ' '{print $2}' | xargs cat | grep -v '^>' | tr -d '\n' | wc -c)
			bseq_len=$(grep 'bsequence' $result_file | awk -F'bsequence ' '{print $2}' | xargs cat | grep -v '^>' | tr -d '\n' | wc -c)

			identity=$(grep 'Identity:' "$result_file" | awk '{print $3, $4}' | sed 's/[()]//g')
			similarity=$(grep 'Similarity:' "$result_file" | awk '{print $3, $4}' | sed 's/[()]//g')
			gaps=$(grep 'Gaps:' "$result_file" | awk '{print $3, $4}')

			# Extract raw and percent components
			identity_raw=$(echo "$identity" | awk '{print $1}')
			identity_pct=$(echo "$identity" | awk '{print $2}')
			sim_raw=$(echo "$similarity" | awk '{print $1}')
			sim_pct=$(echo "$similarity" | awk '{print $2}')

			echo -e "${aseq}\t${aseq_len}\t${bseq}\t${bseq_len}\t${identity_raw}\t${identity_pct}\t${sim_raw}\t${sim_pct}\t${gaps}" >> "$stats_output"
		fi
	done

	awk 'NR==1{print; next} {print | "sort -k6,6nr -k8,8nr"}' "$stats_output" > "${stats_output}.tmp" && mv "${stats_output}.tmp" "$stats_output"
}

export -f variant_alignment

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

module load emboss/6.6.0

#variant_alignment $base/fn/chr1_36797141_DELETION.txt "fn"

find $base/fn/*.txt | parallel -j 60 --joblog $pjl/variant_alignment.txt variant_alignment {} "fn"

COMMENT

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

for file in $base/fn/*txt; do
	folder_name=$(basename $file .txt)
	echo ""
	echo $folder_name
	head -n 6 $base/fn/$folder_name/alignment_summary.tsv
done

COMMENT

# slurm-5208619.out

#=================================================================================
#
#=================================================================================

: << 'COMMENT'



COMMENT

#=================================================================================
#
#=================================================================================

: << 'COMMENT'



COMMENT

#=================================================================================
#
#=================================================================================

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
