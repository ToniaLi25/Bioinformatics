#!/bin/bash
# Run only Step 6: MACS3 and SEACR peak calling for the CUT&Tag pipeline
# Usage:
#   bash run_step6.sh           # run for all samples defined below
#   bash run_step6.sh SAMPLE    # run only for SAMPLE (e.g. WT-K27-2)

set -euo pipefail
IFS=$'\n\t'

# ---------- Config (edit if needed) ----------
project_dir="/mnt/nas/litongtong/Mat2a/MAT2A_H3K4_K9_K27me3_Cuttag"
SEACR_script="/media/liuhb/cell/ltt/software/SEACR-master/SEACR_1.3.sh"
SEACR_MODE="stringent"   # "stringent" or "non" (or "relaxed" depending on SEACR)
macs_out_dir="${project_dir}/05.callpeaks/MACS3"
seacr_out_dir="${project_dir}/05.callpeaks/SEACR"
default_frag_length=150

mkdir -p "${macs_out_dir}" "${seacr_out_dir}"

# ---------- Sample pairs (treatment -> control) ----------
declare -A sample_pairs=(
  ["WT-K27-1"]="WT-IgG-3"
  ["WT-K27-2"]="WT-IgG-6"
  ["WT-K4-1"]="WT-IgG-1"
  ["WT-K4-2"]="WT-IgG-4"
  ["WT-K9-1"]="WT-IgG-2"
  ["WT-K9-2"]="WT-IgG-5"
  ["KO-K27-1"]="KO-IgG-3"
  ["KO-K27-2"]="KO-IgG-6"
  ["KO-K4-1"]="KO-IgG-1"
  ["KO-K4-2"]="KO-IgG-4"
  ["KO-K9-1"]="KO-IgG-2"
  ["KO-K9-2"]="KO-IgG-5"
)

# Optional single-sample argument
SAMPLE_FILTER=""
if [ "$#" -ge 1 ]; then
  SAMPLE_FILTER="$1"
  echo "Running only for sample: ${SAMPLE_FILTER}"
fi

# ---------- Functions ----------
function run_macs3() {
  local treatment="$1"
  local control="$2"
  local predict_log="${macs_out_dir}/${treatment}_predictd.log"
  local out_prefix="${treatment}_vs_${control}"

  echo "[MACS3] Starting ${out_prefix} ..."

  # Try to extract frag_length from existing predictd log; if missing, warn and use default
  if [ -s "${predict_log}" ]; then
    frag_length=$(grep -m1 'predicted fragment length is' "${predict_log}" | awk '{print $(NF-1)}' || true)
    if ! [[ "${frag_length}" =~ ^[0-9]+$ ]]; then
      echo "[MACS3] Warning: invalid frag_length ('${frag_length}') in ${predict_log}; using default ${default_frag_length}"
      frag_length=${default_frag_length}
    fi
  else
    echo "[MACS3] predictd log not found for ${treatment} at ${predict_log}; using default frag_length=${default_frag_length}"
    frag_length=${default_frag_length}
  fi

  # Run MACS3 callpeak (paired treatment/control)
  macs3 callpeak \
    -t "${project_dir}/04.bam/${treatment}_bowtie2.sorted.bam" \
    -c "${project_dir}/04.bam/${control}_bowtie2.sorted.bam" \
    -q 0.01 -g mm --nomodel --extsize "${frag_length}" \
    -n "${out_prefix}" --outdir "${macs_out_dir}" &> "${macs_out_dir}/${treatment}_macs3.log" || echo "[MACS3] callpeak error for ${treatment}; see ${macs_out_dir}/${treatment}_macs3.log" >&2

  # Check outputs
  if compgen -G "${macs_out_dir}/${out_prefix}"* > /dev/null; then
    echo "[MACS3] Completed ${out_prefix}"
  else
    echo "[MACS3] No outputs found for ${out_prefix}; check ${macs_out_dir}/${treatment}_macs3.log" >&2
  fi
}

function run_seacr() {
  local treatment="$1"
  local control="$2"
  local t_bed="${project_dir}/04.bedgraph/${treatment}_bowtie2.normalized.bedgraph"
  local c_bed="${project_dir}/04.bedgraph/${control}_bowtie2.normalized.bedgraph"
  local out_prefix="${seacr_out_dir}/${treatment}_vs_${control}_seacr"
  local logf="${seacr_out_dir}/${treatment}_seacr.log"

  echo "[SEACR] Starting ${treatment} vs ${control} ..."

  if [ ! -s "${t_bed}" ]; then
    echo "[SEACR] Skipping: treatment bedgraph missing or empty: ${t_bed}" >&2
    return
  fi

  if [ -s "${c_bed}" ]; then
    echo "[SEACR] Using control bedgraph: ${c_bed}"
    bash "${SEACR_script}" "${t_bed}" "${c_bed}" "${SEACR_MODE}" "${out_prefix}" &> "${logf}" || echo "[SEACR] failed for ${treatment}; see ${logf}" >&2
  else
    echo "[SEACR] Control bedgraph not found; running SEACR single-sample mode (using 'non' as control token)"
    bash "${SEACR_script}" "${t_bed}" "non" "${SEACR_MODE}" "${out_prefix}" &> "${logf}" || echo "[SEACR] failed for ${treatment}; see ${logf}" >&2
  fi

  # Check outputs
  if compgen -G "${out_prefix}"* > /dev/null; then
    echo "[SEACR] Completed ${treatment} (output prefix ${out_prefix})"
  else
    echo "[SEACR] No outputs for ${treatment}; check ${logf}" >&2
  fi
}

# ---------- Main loop ----------
for treatment in "${!sample_pairs[@]}"; do
  if [ -n "${SAMPLE_FILTER}" ] && [ "${treatment}" != "${SAMPLE_FILTER}" ]; then
    continue
  fi
  control=${sample_pairs[$treatment]}

  # Run MACS3 and SEACR for this pair
  run_macs3 "${treatment}" "${control}"
  run_seacr "${treatment}" "${control}"
done

echo "Step 6 finished. Check ${macs_out_dir} and ${seacr_out_dir} for results and logs."
