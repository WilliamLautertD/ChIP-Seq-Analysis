from __future__ import annotations

from pathlib import Path
import subprocess
import textwrap

import pandas as pd
import streamlit as st

st.set_page_config(page_title="ChIP-seq Pipeline Dashboard", layout="wide")

SAMPLE_COLUMNS = [
    "sample",
    "fastq_r1",
    "fastq_r2",
    "condition",
    "replicate",
    "read_group_library",
    "read_group_unit",
    "control_sample",
]

DEFAULT_ROWS = [
    {
        "sample": "Sample_Chip_01",
        "fastq_r1": "/path/to/Sample_Chip_01_R1.fastq.gz",
        "fastq_r2": "/path/to/Sample_Chip_01_R2.fastq.gz",
        "condition": "H3K27ac",
        "replicate": "1",
        "read_group_library": "lib001",
        "read_group_unit": "unitA",
        "control_sample": "Sample_Input_01",
    },
    {
        "sample": "Sample_Input_01",
        "fastq_r1": "/path/to/Sample_Input_01_R1.fastq.gz",
        "fastq_r2": "/path/to/Sample_Input_01_R2.fastq.gz",
        "condition": "Input",
        "replicate": "1",
        "read_group_library": "lib001",
        "read_group_unit": "unitA",
        "control_sample": "",
    },
]


def load_samples(path: str) -> pd.DataFrame:
    sample_path = Path(path)
    if sample_path.exists():
        try:
            df = pd.read_csv(sample_path, sep="\t")
            for col in SAMPLE_COLUMNS:
                if col not in df.columns:
                    df[col] = ""
            return df[SAMPLE_COLUMNS].fillna("")
        except Exception as exc:  # noqa: BLE001
            st.warning(f"Could not load existing TSV from {path}: {exc}")
    return pd.DataFrame(DEFAULT_ROWS, columns=SAMPLE_COLUMNS)


def validate_samples(df: pd.DataFrame) -> list[str]:
    errors: list[str] = []
    blank_samples = df[df["sample"].astype(str).str.strip() == ""]
    if not blank_samples.empty:
        errors.append("Every row must have a non-empty sample name.")

    missing_fastq = df[
        (df["fastq_r1"].astype(str).str.strip() == "")
        | (df["fastq_r2"].astype(str).str.strip() == "")
    ]
    if not missing_fastq.empty:
        errors.append("Every row must include fastq_r1 and fastq_r2 paths.")

    input_mask = df["condition"].astype(str).str.lower() == "input"
    missing_control = df[~input_mask & (df["control_sample"].astype(str).str.strip() == "")]
    if not missing_control.empty:
        errors.append(
            "Non-Input samples must have a control_sample that references an Input sample ID."
        )

    sample_ids = set(df["sample"].astype(str).str.strip())
    bad_control_refs = [
        ref
        for ref in df["control_sample"].astype(str).str.strip()
        if ref and ref not in sample_ids
    ]
    if bad_control_refs:
        errors.append(
            "Some control_sample values are not present in the sample column: "
            + ", ".join(sorted(set(bad_control_refs)))
        )

    return errors


def results_snapshot(results_root: Path) -> pd.DataFrame:
    rows = []
    folders = {
        "Raw FASTQC HTML": results_root / "qc" / "raw",
        "Trimmed FASTQC HTML": results_root / "qc" / "trimmed",
        "Fastp HTML": results_root / "fastq_trimmed",
        "Filtered BAM": results_root / "mapping" / "bam",
        "BigWig": results_root / "coverage" / "bigwig",
        "bamCompare BigWig": results_root / "coverage" / "bamcompare",
    }

    for label, folder in folders.items():
        if folder.exists():
            if "FASTQC" in label:
                count = len(list(folder.rglob("*_fastqc.html")))
            elif label == "Fastp HTML":
                count = len(list(folder.glob("*.fastp.html")))
            elif label == "Filtered BAM":
                count = len(list(folder.glob("*.filtered.bam")))
            elif label == "BigWig":
                count = len(list(folder.glob("*.bw")))
            else:
                count = len(list(folder.glob("*.ratio.bw")))
            status = "Available" if count > 0 else "Folder exists (no files found yet)"
        else:
            count = 0
            status = "Not found"

        rows.append({"Artifact": label, "Count": count, "Status": status, "Path": str(folder)})

    return pd.DataFrame(rows)


st.title("ChIP-seq Nextflow Dashboard")
st.caption(
    "Build sample sheets + config safely, run the pipeline, and monitor QC/processed outputs "
    "without editing files by hand in the terminal."
)

run_tab, qc_tab = st.tabs(["Run Setup", "QC & Processed Files"])

with run_tab:
    st.subheader("1) HPC login and working location")
    hpc_col1, hpc_col2 = st.columns(2)
    with hpc_col1:
        hpc_username = st.text_input("HPC username", value="")
        hpc_host = st.text_input("HPC host", value="hpc.example.edu")
    with hpc_col2:
        hpc_workdir = st.text_input(
            "HPC project/work directory",
            value="/path/on/hpc/ChIP-Seq-Analysis",
            help="Directory on HPC where this repo and Nextflow pipeline are available.",
        )
        use_ssh_run = st.checkbox(
            "Run on HPC via SSH command",
            value=True,
            help="Uses local ssh command with your configured key-based login.",
        )

    if hpc_username.strip() and hpc_host.strip():
        st.caption(f"HPC session target: `{hpc_username.strip()}@{hpc_host.strip()}`")
    else:
        st.info("Enter HPC username + host so users can run the pipeline on HPC without manual command editing.")

    st.subheader("2) Pipeline locations")
    col_a, col_b = st.columns(2)
    with col_a:
        pipeline_file = st.text_input("Nextflow pipeline file", value="chipseq_main.nf")
        config_file = st.text_input("Nextflow config file", value="chipseq_nextflow.config")
        samples_file = st.text_input("Samples TSV output path", value="config/chipseq_samples.tsv")
    with col_b:
        outdir = st.text_input("Results output directory", value="results/chipseq_pipeline")
        profile = st.text_input("Execution profile", value="slurm")

    st.subheader("3) Reference and runtime settings")
    ref_col1, ref_col2 = st.columns(2)
    with ref_col1:
        reference_fasta = st.text_input("Reference FASTA path", value="/path/to/reference/genome.fa")
        bwa_index_prefix = st.text_input(
            "BWA index prefix (basename)",
            value="/path/to/reference/index_basename",
            help="Equivalent to params.bwa_index_prefix in chipseq_nextflow.config",
        )
    with ref_col2:
        effective_genome_size = st.number_input(
            "Effective genome size", min_value=1, value=2913022398, step=1
        )
        min_mapq = st.number_input("Minimum MAPQ", min_value=0, value=30, step=1)
        exclude_flags = st.number_input("Exclude flags", min_value=0, value=3332, step=1)

    make_bigwig = st.checkbox("Generate per-sample bigWig", value=True)
    make_bamcompare = st.checkbox("Generate ChIP/Input bamCompare bigWig", value=True)

    st.subheader("4) Sample sheet editor")
    current_df = load_samples(samples_file)
    edited_df = st.data_editor(
        current_df,
        num_rows="dynamic",
        use_container_width=True,
        column_config={
            "sample": st.column_config.TextColumn(required=True),
            "fastq_r1": st.column_config.TextColumn(required=True),
            "fastq_r2": st.column_config.TextColumn(required=True),
            "control_sample": st.column_config.TextColumn(
                help="For non-Input rows, set to sample ID of Input control"
            ),
        },
        key="sample_editor",
    )

    errors = validate_samples(edited_df.fillna(""))
    if errors:
        for err in errors:
            st.error(err)
    else:
        st.success("Sample sheet validation passed.")

    write_col1, write_col2 = st.columns(2)
    with write_col1:
        if st.button("Save samples TSV", type="primary", use_container_width=True):
            Path(samples_file).parent.mkdir(parents=True, exist_ok=True)
            edited_df.fillna("").to_csv(samples_file, sep="\t", index=False)
            st.success(f"Saved: {samples_file}")

    generated_config = textwrap.dedent(
        f"""
        params {{
            samples = '{samples_file}'
            outdir = '{outdir}'
            reference_fasta = '{reference_fasta}'
            bwa_index_prefix = '{bwa_index_prefix}'
            min_mapq = {int(min_mapq)}
            exclude_flags = {int(exclude_flags)}
            make_bigwig = {'true' if make_bigwig else 'false'}
            make_bamcompare = {'true' if make_bamcompare else 'false'}
            effective_genome_size = {int(effective_genome_size)}
        }}
        """
    ).strip()

    with write_col2:
        if st.button("Write dashboard config override", use_container_width=True):
            override_path = Path("config/chipseq_dashboard_override.config")
            override_path.write_text(generated_config + "\n", encoding="utf-8")
            st.success(f"Saved: {override_path}")

    st.subheader("5) Run command")
    local_run_cmd = (
        f"nextflow run {pipeline_file} -c {config_file} -c config/chipseq_dashboard_override.config "
        f"-profile {profile} -resume"
    )
    remote_run_cmd = (
        f"ssh {hpc_username.strip()}@{hpc_host.strip()} "
        f"\"cd {hpc_workdir} && {local_run_cmd}\""
        if hpc_username.strip() and hpc_host.strip()
        else ""
    )

    st.markdown("**Local command (run from current machine):**")
    st.code(local_run_cmd, language="bash")
    st.markdown("**HPC command (recommended):**")
    st.code(
        remote_run_cmd if remote_run_cmd else "Provide HPC username and host to generate SSH command.",
        language="bash",
    )

    if st.button("Run pipeline from dashboard"):
        cmd_to_run = remote_run_cmd if use_ssh_run and remote_run_cmd else local_run_cmd
        with st.spinner("Running Nextflow. This can take a long time..."):
            proc = subprocess.run(
                cmd_to_run,
                shell=True,
                check=False,
                capture_output=True,
                text=True,
            )
        st.write("Exit code:", proc.returncode)
        if proc.stdout:
            st.text_area("Stdout", proc.stdout[-15000:], height=240)
        if proc.stderr:
            st.text_area("Stderr", proc.stderr[-15000:], height=240)
        if proc.returncode == 0:
            st.success("Pipeline command finished successfully.")
        else:
            st.error("Pipeline command failed. Check stderr above.")

with qc_tab:
    st.subheader("QC and processed files overview")
    results_root = Path(outdir)
    snapshot_df = results_snapshot(results_root)
    st.dataframe(snapshot_df, use_container_width=True)

    multiqc_report = results_root / "qc" / "multiqc_report.html"
    if multiqc_report.exists():
        st.success(f"MultiQC report found: {multiqc_report}")
    else:
        st.info(f"MultiQC report not found yet at: {multiqc_report}")

    st.markdown(
        """
        **Expected final data flow**

        - QC and process outputs are reviewed here.
        - bigWig (`.bw`) files remain on HPC storage.
        - Users can download `.bw` later via their usual HPC connection tools (SCP/SFTP/Globus).
        """
    )

    if hpc_username.strip() and hpc_host.strip():
        st.markdown("**Optional `.bw` download command template (outside dashboard):**")
        st.code(
            (
                f"scp {hpc_username.strip()}@{hpc_host.strip()}:{outdir}/coverage/bigwig/*.bw ./"
            ),
            language="bash",
        )
