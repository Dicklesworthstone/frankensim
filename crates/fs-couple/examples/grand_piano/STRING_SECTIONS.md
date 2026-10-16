# Material-driven bonded speaking strings

The `grand_piano` example now accepts **`--string-sections FILE`**. This is an
end-to-end path from source-declared density and Young's modulus through
`fs-matdb` material claims, `fs-material::section::ConcentricElasticSection`,
the physical speaking-string mass and bending coefficients, and the existing
contact / string / bridge / soundboard / acoustic output solver. The physics
kernel does not know the material's name.

Input is UTF-8 with exactly these two header records:

```text
frankensim-string-sections-v1
layer,key,outer_radius_m,density_kg_m3,young_pa,source,license
layer,69,0.0005,7800,200000000000,illustrative synthetic steel-like data,CC0-1.0
```

The above file is valid for a **one-key** scale containing MIDI key 69.
For a full piano, declare at least one ordered `layer` for **every** key in
the supplied or preset scale, including unplayed keys. A separate row with
the same key and a strictly larger radius adds another perfectly bonded,
concentric annular layer. The first layer always starts at the centre
(radius zero); hollow strings, omissions, nonpositive/NaN parameters,
reordered radii, duplicate columns and unadmitted keys refuse atomically.
Up to 128 layers per key are admitted. Source and license fields must be
explicit and nonempty; importing them does not verify their claims.

For a single key, first export an existing one-key scale CSV with the
ordinary `--scale` interface, or provide a one-key scale using the documented
SI schema. Then:

```bash
cargo run --release -p fs-couple --example grand_piano -- \
  --scale one-key-scale.csv --string-sections a4-bonded.fsss \
  --render output.wav --note 69 --duration 3
```

Substitute a second file with different density/modulus and the **same**
string geometry and applied tension. The new density changes modal inertia
and frequency; the new elastic modulus changes flexural rigidity and
inharmonicity. The string's strike projection, bridge forcing, acoustic
radiation, and pressure WAV are all downstream of the same `Course`.
Tension remains physical prestress rather than a secret compensating
detune; explicit `--concert-pitch` or `--stretch-tuning` retunes *after*
the material substitution, if requested. `--dump-scale` exports the
effective coefficients actually fed to the solver.

**Physical boundary.** This is the *perfectly bonded, linear-elastic,
concentric rod* reduction `mu=sum(rho A)`, `EI=sum(E I)`. It is appropriate
to homogeneous or genuinely bonded circular composites, but **not**
a mechanically loose helically wound piano bass string. Wound strings
require winding pitch, wrapping inertia, and an independently justified
effective bending/slip model; this interface explicitly makes no claim
that a copper or silver winding is a bonded rod. It does not infer
temperature-dependent material states, axial EA for nonlinear stretching,
viscoelastic bending loss, fatigue, or contact properties from two scalars.
The `--string-stretching` and source-specific `--rt0425-string-damping`
options therefore cannot be combined with this override. This is
parameter-driven physical propagation, not experimental validation of
gold/silver alloy strings.

The source model deliberately keeps material claims and scientific
evidence distinct: property declarations without specimen observations
are unvalidated assertions, not certified measurements. The existing
energy, mesh and audio-band limitations still apply.
