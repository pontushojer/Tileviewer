# Tileviewer

Tileviewer is an interactive viewer for spatial mapping data from an Illumina
TruPath experiment. In TruPath, long DNA molecules are used to generate
clusters directly on a flowcell, so clusters from the same original molecule
appear spatially close to each other. Tileviewer plots each cluster at its
flowcell tile position and colors it by where its read pair maps in the
genome, making it easy to see whether neighboring clusters map to nearby
genomic positions.

Tileviewer needs two inputs:

- A sparse matrix (`.npz`) of cluster positions on one tile.
- The FASTA index (`.fai`) of the reference used for mapping.

[Section 1](#1-install-and-run-tileviewer) covers running the viewer from existing inputs. [Section 2](#2-create-the-input-files-from-fastq) at the end,
shows how to create them from FASTQ files.

For example input files, see the [example](example) directory.

## 1. Install and run Tileviewer

### Install

Install the locked Pixi environment from the project directory (containing
`pixi.toml` and `pixi.lock`):

```bash
pixi install
```

### Launch

```bash
pixi run tileviewer
```

Open the URL printed by Streamlit, then upload the `.npz` matrix and the
matching `.fai` in the sidebar. The `.fai` must come from the same reference
used for mapping.

### Input format

- The `.npz` is a SciPy CSR matrix with integer values. Row index is X and
  column index is Y.
- A mapped cluster stores its 1-based absolute genome coordinate, where
  chromosomes are concatenated in `.fai` order.
- A cluster that is unmapped, filtered, or has an incomplete pair is stored as
  an explicit `0`. 

### Using the viewer

- **Crop:** use the sidebar sliders, or "Enter exact bounds" for exact
  coordinates. Bounds are zero-based and the upper bound is exclusive.
- **Dot view:** crops containing at most 100,000 clusters show each cluster as
  an individually sized dot (set the diameter in the sidebar). Larger crops
  show a raster overview; crop further for dots and full hover coverage.
- **Colors:** mapped clusters are colored by absolute genome position. Unmapped
  clusters are shown as magenta crosses (color is adjustable).
- **Chromosome + interval:** in dot view, choose a chromosome and interval to
  spread the colormap across that region and grey out clusters outside it. This
  makes it easier to spot spatially adjacent clusters mapping to nearby
  positions.
- **Hover:** shows grid coordinates, absolute position, and
  chromosome:position.
- **Summary:** the panel reports grid cells, mapped and unmapped clusters,
  empty grid cells, mapped density, and a chromosome distribution for the crop.

## 2. Create the input files from FASTQ

These standalone command-line examples extract one tile's paired reads, map
them to a reference, and convert the alignment to the `.npz` matrix using the
included `sam_to_npz.py`. Snakemake is not required.

### Requirements

- `seqkit` for FASTQ extraction, `minibwa` for indexing and mapping, and
  `samtools` for creating the reference index. Install or load these separately
  if they are not already on your `PATH`.
- The paired FASTQ files for the run and the reference FASTA used for mapping.

Create the reference index, which is also the `.fai` uploaded to the viewer:

```bash
samtools faidx genome.fa
```

### 2.1 Extract the reads for one tile

Set the paths and identifiers used in the read names. This example assumes
names of the form `instrument:run:flowcell:lane:tile:x:y`; adjust the pattern
to match your data.

```bash
FASTQ_DIR=/path/to/fastqs
FLOWCELL=44KDCAGV1
LANE=1
TILE=1101

R1_IN="$FASTQ_DIR/sample_R1.fastq.gz"
R2_IN="$FASTQ_DIR/sample_R2.fastq.gz"

seqkit grep -r -n \
  -p "LH00207:420:${FLOWCELL}:${LANE}:${TILE}:.*" \
  -o "${TILE}_R1.fastq.gz" "$R1_IN"
seqkit grep -r -n \
  -p "LH00207:420:${FLOWCELL}:${LANE}:${TILE}:.*" \
  -o "${TILE}_R2.fastq.gz" "$R2_IN"
```

The same pattern is applied to both mates so the pairs stay consistent. Confirm
both outputs contain reads before mapping.

### 2.2 Index the reference and map the tile

Build the `minibwa` index once per reference, then map the pair of FASTQs:

```bash
THREADS=20
minibwa index -t"$THREADS" genome.fa
minibwa map \
  -o "${TILE}.sam" \
  -t"$THREADS" \
  genome.fa \
  "${TILE}_R1.fastq.gz" \
  "${TILE}_R2.fastq.gz"
```

The SAM `QNAME` must keep the tile X and Y coordinates as its final two
colon-separated fields; the converter uses them to place each cluster.

### 2.3 Convert the SAM to NPZ

```bash
mkdir -p matrix
pixi run python sam_to_npz.py \
  --sam "${TILE}.sam" \
  --fai genome.fa.fai \
  --out "${TILE}.npz"
```

Defaults require MAPQ >= 40 and a maximum difference of 2,000 between the two
reads' SAM `POS` values. Change them with `--min-mapq` and `--max-insert-size`:

```bash
pixi run python sam_to_npz.py \
  --sam "${TILE}.sam" \
  --fai genome.fa.fai \
  --out "${TILE}.npz" \
  --min-mapq 30 \
  --max-insert-size 3000
```

Notes:

- The insert-size filter uses the absolute difference between the two `POS`
  values, not the CIGAR-based fragment span or template length. The pair
  midpoint is likewise computed from `POS` values.
- The helper prints mapped and unmapped counts, incomplete pairs, references
  absent from the FAI, and duplicate X/Y coordinates. For a duplicate
  coordinate it keeps one value, preferring a mapped one; review that count to
  catch coordinate collisions or duplicated records.
- The helper reads SAM text and holds pair data in memory before building the
  matrix, so process one tile at a time. BAM input is not supported directly.

Then follow section 1 to view `${TILE}.npz` with `genome.fa.fai`.
