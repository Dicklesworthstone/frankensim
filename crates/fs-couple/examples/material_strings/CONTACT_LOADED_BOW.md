# Bow traction from solved normal contact

`--stop-profile` can pair an existing bow tangent with its corresponding normal
contact. Add `bow-contact BOW_INDEX STOP_ID` after that stop is declared and
before the optional terminal motion block. Bow indices follow canonical
schedule-track order. Both sides must occupy the same string and station;
one normal may supply only one scalar tangential law.

Example, using the existing synthetic material packs from README.md:

```bash
cargo run --release -p fs-couple --example material_strings -- \
  /path/to/contact-bow.material-strings /path/to/new-contact-bow.wav \
  --stop-profile crates/fs-couple/examples/material_strings/contact-bow.stops
```

Copy `contact-bow.material-strings` beside its compiled `wire.fsmatpk` and
`body.fsmatpk`. The output path must not exist. The example contact law is an
authored synthetic fixture, NOT measured hair, rosin or finger tissue.

```text
frankensim-string-stops-v1
stop bow_surface 0 0.2 normal 0 -1 fixed 0 5000 1 0.2 10 0.001
bow-contact 0 bow_surface
motion 0.001 0.1 3
knot 480 -0.000003
knot 960 -0.000003
knot 1440 0.00003
```

For a contact-loaded bow, every scheduled normal-force value MUST be zero:
there is no additional prescribed load. Zero does not disable its friction.
Its normal reaction is solved from the moving support, string and receiver
motion; that SAME held reaction sets the friction capacity in the joint step.
When the reaction vanishes, so does traction. Tangential speed remains authored.
Unmapped bows retain the original prescribed-load/release semantics. Conflicting
initial or future normal commands refuse before material loading.

Two transverse planes are required. The explicit normal must be perpendicular
to the bow tangent. With torsion selected, `normal x tangent = +1` must match the
existing positive-radius torque convention; the opposite surface orientation
refuses rather than silently changing the contact map. Thus normal `(0,-1)`
works with tangent `(1,0)`. The material adapter rechecks the actual geometry.
A selected `bridge:STRING` receiver still has to measure the declared normal.
It does not invent a tangential degree of freedom for that body.

Material density, elasticity, bending/torsion losses and gas continue through
the existing construction and full joint response. Normal actuation work stays
separate from tangential work and dissipation. The original interface card,
when supplied through the Rust API, checks the actual resolved pressure and
all three slip rates. File friction coefficients still have authored authority;
bulk solid properties do not by themselves determine an interface friction law.

This is one tangent per normal, fixed contact stations and prescribed surface
motion. It does not solve a shared bow-hair structure, evolving contact area,
rosin temperature/wear, or friction in two tangential directions. The original
bounded coordinate solve can refuse near strongly coupled/nonmonotone regimes.
Held-reaction, contact, modal and time-step refinement remain necessary; this
is not measured violin/cello/bass fidelity. No pitch or gain is assigned.

Focused checks (added, native execution unavailable in the authoring session):

```bash
cargo test -p fs-couple --lib multiple::normal_load -- --test-threads=1
cargo test -p fs-couple --test material_bow_interface contact_load:: -- --test-threads=1
cargo test -p fs-couple --example material_strings contact_loaded_bow:: -- --test-threads=1
```
