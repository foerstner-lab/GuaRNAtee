# GuaRNAtee

## What is GuaRNAtee?
GuaRNAtee is a genome-wide computational tool for the suggestive identification of full-length small RNA transcripts in prokaryotes. It integrates differential RNA-Seq (dRNA-Seq) and Term-Seq coverage signals to detect coherent 5′ and 3′ transcript extremities for each condition. It relies on a cutoff-free, outlier-based approach applied in sliding windows over transcriptome coverage to extract extremities based solely on the local signal behavior. These paired coverage anomalies are nominated as potential sRNA loci, enabling the robust detection of transcript boundaries even under low signal-to-noise conditions that are typical for bacterial sRNAs. GuaRNAtee therefore prioritizes candidates for downstream validation and genome annotation without overclaiming final certainty.

In benchmarking against existing annotations, GuaRNAtee was able to recall almost all previously annotated small RNAs, while additionally suggesting numerous novel candidates with rankings. By directly leveraging sequencing-derived terminus evidence, GuaRNAtee supports hypothesis generation in small RNA biology and facilitates the discovery of previously unannotated RNA elements, including those embedded within coding sequences or UTRs.

---

## Installation

### Prerequisites
- Python == 3.12
- Dependencies listed in [requirements.txt](requirements.txt)

### Install Dependencies
```bash
pip install -r requirements.txt
```

---

## Usage

### Quick Start

```bash
python guarnatee.py \
  --gffs genome_annotation.gff known_sRNAs.gff \
  --fastas genome_sequence.fasta \
  --wigs data.wig:condition:replicate:strand:library_mode \
         data.wig:condition:replicate:strand:library_mode \
  --out_dir output/
```

### Command Line Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--gffs` | Yes | GFF3 annotation files (space separated) |
| `--fastas` | Yes | FASTA genome files (space separated) |
| `--wigs` | Yes | Wiggle files with annotations (see format below) |
| `--out_dir` | Yes | Output directory for results |
| `--config_file` | No | Configuration file (default: config.cfg) |
| `--verbose` | No | Enable verbose logging for debugging |
| `--version` | No | Show version and exit |

### Wiggle File Format

The `--wigs` argument takes this format:
```
<file_path>:<condition>:<replicate>:<strand>:<lib_mode>
```

**Fields:**
- `file_path` - Path to wiggle file
- `condition` - Experimental condition (e.g., OD_0.2, LB_0.4)
- `replicate` - Replicate number (e.g., 1, 2, 3)
- `strand` - Strand: `f` (forward) or `r` (reverse)
- `lib_mode` - Library mode (see below)

**Library Modes:**
| Mode | Description | Use Case |
|------|-------------|----------|
| `FL` | Full-length | Both 5' and 3' ends in same file |
| `P1` | Paired 5' end | Paired-end sequencing (5' end) |
| `P2` | Paired 3' end | Paired-end sequencing (3' end) |
| `5E` | 5' enriched | Separate 5' end enrichment |
| `3E` | 3' enriched | Separate 3' end enrichment |
| `d5E` | Differential 5' | Differential 5' end enrichment |

### Example Usage

#### Example 1: Full-length library
```bash
python guarnatee.py \
  --gffs genome.gff known_sRNAs.gff \
  --fastas genome.fasta \
  --wigs data_forward.wig:treatment:1:f:FL \
         data_reverse.wig:treatment:1:r:FL \
  --out_dir results/
```

#### Example 2: Paired-end library
```bash
python guarnatee.py \
  --gffs genome.gff known_sRNAs.gff \
  --fastas genome.fasta \
  --wigs p1_fwd.wig:OD_0.2:1:f:P1 \
  --wigs p2_fwd.wig:OD_0.2:1:f:P2 \
  --wigs p1_rev.wig:OD_0.2:1:r:P1 \
         p2_rev.wig:OD_0.2:1:r:P2 \
  --out_dir output/
```

#### Example 3: Multiple conditions with custom config
```bash
python guarnatee.py \
  --gffs annotations.gff known_sRNAs.gff \
  --fastas reference.fasta \
  --wigs cond1_rep1.wig:condition1:1:f:FL \
         cond1_rep2.wig:condition1:2:f:FL \
         cond2_rep1.wig:condition2:1:f:FL \
         cond2_rep2.wig:condition2:2:f:FL \
  --config_file custom_config.cfg \
  --threshold 2 \
  --out_dir results/ \
  --verbose
```

---

## Output Files

GuaRNAtee generates the following output files:

### GFF Files
```
<organism>_candidates.gff
```
RNA candidates in GFF3 format, one file per organism.

### Excel Files
```
<organism>_candidates.xlsx
```
Organized workbook with sheets:
- **intergenic** - Intergenic and ncRNA candidates
- **ORF_int** - Candidates internal to ORFs
- **others** - Antisense and cross-feature candidates

### Statistics File
```
stats.tsv
```
Processing statistics including peak counts and connections.

---

## Configuration

Default configuration is in `config.cfg`. Key parameters:

```ini
[ALL_CONFIGURATIONS]
min_height=10          # Minimum peak height
min_len=40            # Minimum candidate length
max_len=400           # Maximum candidate length
min_step_factor=1.2   # Minimum step factor for peaks
min_distance=10       # Minimum distance between peaks
```

See [config.cfg](config.cfg) for all options.

---

## License
This software is licensed under [ISC License](LICENSE)

---

## Citation
If you use this software, please cite:

```bibtex
@software{guarnatee,
  title = {GuaRNAtee: Genome-wide RNA annotation tool},
  author = {Elhossary, Muhammad and Förstner, Konrad U.},
  year = {2021},
  url = {https://github.com/foerstner-lab/GuaRNAtee}
}
```

See [CITATION.cff](CITATION.cff) for details.

---

## Contributing

We welcome contributions! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

---

## Support

- **Issues:** [GitHub Issues](https://github.com/foerstner-lab/GuaRNAtee/issues)
- **Email:** elhossary@zbmed.de
- **Documentation:** See docs listed above

---

## Authors

- **Muhammad Elhossary** - *Lead Developer* - elhossary@zbmed.de
- **Konrad U. Förstner** - *Principal Investigator* - konrad@foerstner.org
