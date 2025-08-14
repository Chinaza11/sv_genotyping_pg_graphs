#! /bin/bash

#SBATCH --partition=Orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
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
# Plan for evaluating the SVs genotyped by the pan-genome graph: Plan1
#=================================================================================

: << 'COMMENT'

- Pick a genotype: riouncc.
- Simulate SVs using SURVIVOR. This outputs a new altered genome of the riouncc.
- Rebuild the pan-genome graph and include this newly altered genome.
- Simulate short reads from the altered genome and the original genome using ART.
- Map simulated short reads (from the altered and the original genome) to the new pan-genome graph and genotype variants.
- Compare SVs called from both simulated short reads data (from the altered and the original genome). Tools: Truvari. Ideal outcome: 
	- All SVs from the simulated short read data from the original genome is present in the SVs from the altered genome.
	- The difference between the SVs in the original genome and the altered genome is the artificial SVs created.
	- Extraneous SVs from the altered genome simulated short reads are minimal or non-existent.

- At the end of my analysis I realized that Plan1 will not allow me achieve my goal 
  because simulated short reads data from original genome and altered genome might 
  cover different regions of the pan-genome.

COMMENT

#=================================================================================
# simulating structural variants
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate survivor

cd evaluation/sandbox/

# SURVIVOR simSV parameter_file
#parameter_file was edited so SURVIVOR simulates 2000 indels of length 50-2000

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
# Plan for evaluating the SVs genotyped by the pan-genome graph: Plan2
#=================================================================================

: << 'COMMENT'

- Pick a genotype: riouncc.
- Simulate SVs using SURVIVOR. This outputs a new altered genome of the riouncc.
- Rebuild the pan-genome graph and include this newly altered genome.
- Simulate short reads from the altered genome using ART.
- Map simulated short reads to the original pan-genome graph (that does not have the altered genome) and genotype variants.
- Map simulated short reads to the new pan-genome graph (that has the altered genome) and genotype variants.
- Compare SVs genotyped when the original pan-genome graph is used vs when the new pan-genome graph is used. Tools: Truvari. Ideal outcome:
	- All SVs genotyped when the original pan-genome graph is used is present in the SVs genotyped when the new pan-genome graph is used.
	- The difference between the two groups of SVs is the artificial SVs created.

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
# compare SVs genotyped when the original pan-genome graph is used vs when the new pan-genome graph is used
#=================================================================================

: << 'COMMENT'

# --------------------------------------------------------------------------------
# copy SVs to be compared to truvari folder
cp $base/vg_call/simulated.filtered.vcf $base/truvari/sv_from_mapping_to_new_pangenome.vcf
cp $base/vg_call/mapped_to_original_pangenome/simulated.filtered.vcf $base/truvari/sv_from_mapping_to_original_pangenome.vcf

# --------------------------------------------------------------------------------

cd $base/truvari

# --------------------------------------------------------------------------------
# removing riouncc#0# from chromosome name so it doesn't throw an error
sed -i 's/riouncc#0#//g' sv_from_mapping_to_new_pangenome.vcf
sed -i 's/riouncc#0#//g' sv_from_mapping_to_original_pangenome.vcf

# --------------------------------------------------------------------------------
# remove multiallelic loci from vcf files so truvari doesn't throw an error

remove_multiallelic_loci() {

        local sample=$(basename $file | awk -F'.' '{print $1}')
        echo "sample:" $sample

        echo ""
        local unique_genotypes=$(grep -v '^#' $1 | cut -f10- | \
                                        tr '\t' '\n' | cut -d ':' -f1 | \
                                        sort | uniq | paste -sd '|' -)

        echo "unique_genotypes:" $unique_genotypes

        grep -Ev '1/2' $1 > $sample.filtered2.vcf

        echo ""
        local unique_genotypes=$(grep -v '^#' $sample.filtered2.vcf | cut -f10- | \
                                        tr '\t' '\n' | cut -d ':' -f1 | \
                                        sort | uniq | paste -sd '|' -)

        echo "unique_genotypes:" $unique_genotypes

}

for file in *.vcf;
do
        echo "--------------------------------------------------------------"
        echo $file
        remove_multiallelic_loci $file
done


# --------------------------------------------------------------------------------

# gzip vcf files, index them and run truvari

module load samtools/1.19

bgzip *filtered2.vcf
for file in *.gz; do echo $file; tabix -p vcf $file; done

truvari bench -b sv_from_mapping_to_new_pangenome.filtered2.vcf.gz \
	-c sv_from_mapping_to_original_pangenome.filtered2.vcf.gz -o output_dir

COMMENT

#slurm-5717474.out

#=================================================================================
# comments from truvari result
#=================================================================================

: << 'COMMENT'

- Based on the ideal outcome in Plan2, 128 SVs (FN) should be the newly simulated SVs.
	- Below code will investigate this.
- Expected less non-matching calls from sv_from_mapping_to_original_pangenome.filtered2.vcf,
  110 SVs (FP) is high.
	- Need to look into Truvari parameters some more.

COMMENT

#=================================================================================
# function for extracting the simulated deletions
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
# extracting the simulated deletions from riouncc that gave rise to the altered genome
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
# function for extracting the deletions and insertions from a VCF file into individual files
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
# extracting the deletions and insertions in fn.vcf into individual files
#=================================================================================

: << 'COMMENT'

cd $base
mkdir -p fn
cp truvari/output_dir/fn.vcf.gz fn

cd fn
gunzip fn.vcf.gz

process_vcf_variants fn.vcf "individual_fn"
 
# slurm-5731284.out

COMMENT

#=================================================================================
# function for extracting sequences in a fasta file into individual files
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
# extracting the sequences in the simulated insertions and deletions fasta files
# into individual files
#=================================================================================

: << 'COMMENT'

cd $base
split_fasta_to_files simulated.insertions.fa simulated.insertions
split_fasta_to_files simulated.deletions.fa simulated.deletions

COMMENT

#=================================================================================
# function for comparing the fn variants with the simulated variants using the 
# smith-waterman algorithm
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
# compare the fn variants with the simulated variants using the smith-waterman algorithm
#=================================================================================

: << 'COMMENT'

module load emboss/6.6.0

#variant_alignment $base/fn/chr1_36797141_DELETION.txt "fn"

find $base/fn/individual_fn/*.txt | parallel -j 60 --joblog $pjl/variant_alignment.txt variant_alignment {} "fn/individual_fn"

# slurm-5731421.out

COMMENT

#=================================================================================
# get the top 5 result of the alignment of each fn variants to the simulated variants 
#=================================================================================

: << 'COMMENT'

{
	for file in $base/fn/individual_fn/*txt; do
		folder_name=$(basename $file .txt)
		echo ""
		echo $folder_name
		head -n 6 $base/fn/individual_fn/$folder_name/alignment_summary.tsv
	done
} > $base/fn/top_5_result.txt

COMMENT

#=================================================================================
# examining the exact matches/close enough matches between the fn variants and the
# simulated variants
#=================================================================================

: << 'COMMENT'

for file in $base/fn/individual_fn/*txt; do
    folder_name=$(basename $file .txt)

    result=$(head -n 11 "$base/fn/individual_fn/$folder_name/alignment_summary.tsv" \
        | tail -n +2 \
        | awk '{
            a_len = $2
            b_len = $4

            split($5, id_raw, "/")
            id_pct = gensub(/%/, "", "g", $6)

            smallest = (a_len < b_len) ? a_len : b_len
            largest = (a_len > b_len) ? a_len : b_len

            if (id_raw[2] >= 0.75 * largest &&
                ((a_len - b_len) < 0 ? (b_len - a_len) : (a_len - b_len)) <= 0.25 * smallest &&
                id_pct > 80) {
                print
            }
        }')

    if [[ -n $result ]]; then
        echo ""
        echo "$folder_name"
        echo "$result"
    fi
done

COMMENT

# slurm-5777606.out: no exact matches or close enough matches

#=================================================================================
# examining simulated variants localized in a subsection of the fn variants and vice versa
#=================================================================================

: << 'COMMENT'

count=0

for file in $base/fn/individual_fn/*txt; do
    folder_name=$(basename $file .txt)

    result=$(head -n 11 "$base/fn/individual_fn/$folder_name/alignment_summary.tsv" \
        | tail -n +2 \
        | awk '{
            a_len = $2
            b_len = $4

            split($5, id_raw, "/")
            id_pct = gensub(/%/, "", "g", $6)

            smallest = (a_len < b_len) ? a_len : b_len
            largest = (a_len > b_len) ? a_len : b_len

            if ((id_raw[2] >= 0.75 * smallest) && id_pct > 90) {
                print
            }
        }')

    if [[ -n $result ]]; then
        ((count++))
	echo ""
        echo "$folder_name"
        echo "$result"
    fi
done

echo ""
echo "Total matches: $count"

COMMENT

# slurm-5783563.out: 18 variants (out of the 128). Would be less than 18 after removal
# 		     of the TATA regions.

#=================================================================================
#
#=================================================================================

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
