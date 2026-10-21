# Two-axis material source laws

Manifest v9 and material source v3 can supply one scalar property as a function
of two independent physical coordinates. The ordinary `matdb-pack` compiler
normalizes the coordinates and values into a `PropertyValue::Surface2D` in
FSMATPK v8. Material-card and ordered-interface packs reuse the same payload.
Older manifest and source versions retain their existing interpretation.

## Source grammar

The new `surface2d` record has twelve tab-separated fields:

| Position | Field |
| --- | --- |
| 1 | Literal `surface2d` |
| 2–4 | Claim id, observation references, property name |
| 5–6 | First source axis name and unit |
| 7–8 | Second source axis name and unit |
| 9 | Property-value unit |
| 10 | Comma-separated `first-coordinate:second-coordinate:value` samples |
| 11 | Comma-separated `first-index:second-index:third-index` triangles, or `-` |
| 12 | `linear` or `tabulated` |

Indices are canonical unsigned decimal integers, starting at zero, and refer to
the declared sample order. The compiler preserves that order and connectivity.
Both axes require explicit manifest `axis` rows; `dimensional` declares that
only dimensions are asserted. Existing property mappings, validity,
uncertainty, frame and provenance records apply unchanged. Both source validity
intervals must contain every normalized sample.

This synthetic example declares a temperature/pressure response. It demonstrates
the format and is not a measured material dataset.

`manifest.tsv`:

```text
frankensim.matdb-manifest.v9
pack_id	synthetic-two-axis
redistribution	permitted	Synthetic example, CC0-1.0
citation	Authored format example; no empirical material claim
license	CC0-1.0
source	primary	source.tsv	material-tsv-v1
axis	temp	T	absolute-temperature
axis	pressure	p	dimensional
```

`source.tsv`:

```text
frankensim.matdb-source.v3
observation	test	synthetic specimen	authored two-axis law	no empirical claim
surface2d	response	test	response	temp	degC	pressure	kPa	GPa	0:100:10,100:100:20,0:200:40	0:1:2	linear
uncertainty	response	unstated	-	-	-	-
validity	response	temp	0	100	degC
validity	response	pressure	100	200	kPa
```

Compile with the existing command:

```bash
cargo run -p xtask -- matdb-pack --manifest manifest.tsv --out example.fsmatpk
```

## Support and units

`linear` permits barycentric interpolation only inside declared triangles.
Disconnected patches and holes remain disconnected or unsupported. The
rectangular validity bounds do not fill the missing upper-right corner of the
example's triangle. `tabulated` admits only exact stored coordinate pairs;
`-` supplies no connecting triangles. Isolated exact samples can coexist with
triangle patches. No triangulation, grid completion or extrapolation is inferred.

The shared runtime rejects degenerate/overlapping triangles, invalid indices,
duplicate coordinates and unsupported source geometry. Runtime sample and
triangle budgets also bound compilation. Coordinate normalization uses the
same declared quantity kinds and unit transforms as validity bounds, including
absolute-temperature offsets and frequency conventions. Each of the three
sample components retains its own source-unit normalization receipt.

Joint-statistics references are `claim:surface-x:sample`,
`claim:surface-y:sample`, and `claim:surface-value:sample`. Here `x` and `y`
mean the first and second declared coordinates, not spatial directions.
Covariance uses the existing product of component scales; affine coordinate
offsets never enter covariance. Interpolation transports the declared nominal
law; it does not supply measured accuracy or infer uncertainty correlation.
