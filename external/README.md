# External repositories adapted for HiRepEat framework

We build functionality on top of the following repositories from Gene Myers and Richard Durbin for LTR-RT detection framework, as described below. 

---
## alntools (https://github.com/richarddurbin/alntools) 


### `taco.c` — generalization

Generalized `TACO` idea for compressing tandem repeats, and supporting coordinate transformations as described below (SALSA):

```
taco compress    [-o <prefix>] [-f] <input.1aln> <seqFile>
taco collapse    [-o <prefix>] [-f] <spec.tsv> <input.fa>
taco reconstruct [-o <output.fa>] <input.1taco> <input.1seq>
taco chain       [-o <prefix>] <a.1taco> <a.1seq> <b.1taco> <b.1seq> [...]
taco merge       [-o <prefix>] <a.1taco> <a.1seq> <b.1taco> <b.1seq>
taco lift        [-r] <annotations.bed> <map.1taco> [map2.1taco ...]
taco info        <input.1taco>
taco windows     [-f <flank>] <input.1taco>
```

`compress` outputs **2** files — `.1taco` (coordinate map)
plus `.1seq` (the compressed sequence)

Our clustering and iterative collapse approach for nested LTR-RTs depends on `compress`, `collapse`, `reconstruct`, `lift`, `info`
and `windows` — see [engine/pipeline.py](../engine/pipeline.py) and
[seq/info.py](../seq/info.py). `windows` is used for iterative collapse at specific regions, as described below.

### `salsa.{c,h}` the secret sauce for `taco`
**S**orted **A**rray for **L**iftover of **S**equence **A**ddresses: the coordinate map between original and taco-compressed sequences.  One `salsa` per
sequence, holding a sorted array of per-compression-event records plus a
cumulative-shift array for O(log n) coordinate transforms.  We define the `.1taco`
OneCode file type (`c` sequence / `I` id / `L` lift event / `S` removed DNA).
This schema is also in `alntools.h` so `taco` can read.

### `-r <file.bed>` — addition to FastLTR, FasTAN

For the iterative collapse approach
`-r` takes a BED file specifiying
`name beg end` (interval for restriction).


---
## FASTAN (https://github.com/thegenemyers/FASTAN) 

- **Sub-8bp blocks** (`FasTAN.c`, `FastLTR.c`): `spectrum_block` now returns −1
  immediately for `len < 8`, which holds no 8-mer, handling an edge case in our iterative LTR-RT detection.
- **`GDB.c`**: `scount`/`Icount`/`ncount` are zeroed before the `oneStats` calls, otherwise being read uninitialized 
- **`FastLTR.c`**: `DIAG_MIN` 3000 → 2000, to detect smaller but still full-length LTR-RTs that are otherwise missed.
- **`FastLTR.c, FASTAN.c`**: support for the `-r` option to restrict to specific regions, used in iterative collapse.

## ONEcode (https://github.com/thegenemyers/ONEcode)

Mostly identical to original, with some patches around writing `TMPDIR` temp files; and making `int`→`I64` in `vcEncode`/`vcDecode`/`vcAddToTable`/
`Compress_DNA`/`Uncompress_DNA`) to handle >4Gb chromosome lengths.

