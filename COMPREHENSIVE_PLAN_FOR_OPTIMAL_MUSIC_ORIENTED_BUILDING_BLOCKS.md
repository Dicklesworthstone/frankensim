# Comprehensive Plan for Optimal Music-Oriented Building Blocks

> Nameless physics. Emergent sound. **Menus, not winners.**
> Each filling (woodwind, brass, jet-pipe, string, piano, voice,
> circuit, shell) carries several *orthogonal images* of the
> same cards. They do not share state layouts. They couple
> only at power ports. A bake-off on named QoIs decides
> what may claim what — several may remain legal at once.
>
> Companion to `COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md` and
> new-domains §5.11. Does not amend the Decalogue.
>
> Status: **v8**. v1–v4 each crowned one discretization.
> v5 gave each filling one stack. v6 made the **menu** the
> design. v7 spelled out the ingest law (mesh + cards →
> charts → images). **v8 binds the plan to the tree**:
> every image names its owning crate or extension seam
> (D23), the claims matrix becomes a machine-checked
> registry (§8), real-time becomes a measured budget law
> (§9), the product output surface gets a track (T-Out),
> known duplications get owners (§10), and the plan is
> FROZEN — further revisions require executed beads (§17).
>
> Tags: `[S]` `[F]` `[M]`. Date: 2026-08-14.

---

## 0. How to read this

§1.4 is the rule. §3 is exactness classes. §5 is the menus —
now with an **Owner** column (D23). §6 is the selector.
**§7 is the ingest law**: a licensed mesh plus `fs-matdb` /
`GasState` / `Uniaxial` / visco cards become charts, then
authority images, then performance images. §7.2 is the
trumpet. §7.5 is the piano (felt on the *existing*
hysteresis / visco stack). §7.8 is live vs not — corrected
against the actual tree. **§8 is the claims registry — the
product.** §9 is the budget law and the output surface.
§10 is the reuse ledger (duplications and their owners).
§11 is governance. §12 is parallel tracks with bead IDs.

---

## 1. Intent

### 1.1 Hear physics, play physics

Geometry + materials + gesture + air. No samples. No
instrument crates. Real-time on a named core *when the
selected image's measured budget row says so* (§9) — not
when a mythical universal stepper does.

### 1.2 Success is a matrix, not a scalar

A filling **F** with image **I** succeeds on claim **C**
when the registry gate for `(F, I, C)` is green (§8).
Example:

- Clarinet × TMM × “Ernoult peak cents” can pass while
  clarinet × 6-cell pHS × “brightness to 6 kHz” is
  **refused** (Nyquist). That is not a conflict.
- Piano × modal ZOH × “partials” can pass while
  piano × characteristic air × anything is **not on
  the menu**.

We do not require one image to win every claim.

### 1.3 Why one stack per family was still dumb

A sax **phrase** wants spatial characteristics + `Y(σ)`.
A sax **held note** wants a 12-state VFIT of TMM `R(ω)`
(cheaper, often more accurate in-band than a long split
if you only care about the loop). A sax **design
question** (“where is the 3rd peak?”) wants TMM, not
time stepping. Those three are **orthogonal**. Shipping
only the phrase image makes sustain expensive. Shipping
only VFIT makes slurs a sample bank. Shipping only TMM
makes attacks impossible.

A guitar **pluck** wants modal ZOH. A guitar **wolf /
body filter** can be a VFIT of the plate 1-port *or* a
modal bake — bake them off. A guitar **electric** path
should not touch either air model.

v5 listed one “RT” column. v6+ lists a **menu** and a
**selector**, not a winner.

### 1.4 Design rule

```text
cards  →  several images {I_k}, each with:
            operator limit, exactness class,
            QoI set, budget, failure mode, OWNER
       →  selector(gesture, claim, CPU)
       →  one or more kernels this block
```

Images of the same filling **may run in the same
process** (HB offline “will it speak?” + TD attack) or
**never share a binary path** (piano never links the
tube stepper). Non-overlapping is a feature.

---

## 2. Shared doctrine (still thin)

D1–D13 as before (nameless, no samples, lanes/rates,
reduce don’t invent constitutives, correctness before
*that kernel’s* fusion, Franken-only, certificates,
named budgets, psycho judges, flags, passive LTV ports,
no synth crate, calibrate cards).

**D14.** Newton only nonlinear islands.

**D15.** Coarsen inside a class, toward *that image’s*
QoIs.

**D16.** Multi-rate by island stiffness.

**D17.** Lift when *that* image’s coordinates change.

**D18.** Electrical-length split is **one optional image
of 1-D air**, not a law for plates or circuits.

**D19.** Claims are `(filling, image, QoI)` — and they
live in the registry (§8), not in prose.

**D20.** Interiors couple only at ports. New ports are
**`fs-couple` `PortSchema` v2 descriptors**; the raw
legacy `Port` container is refused by `interconnect`
and must never be re-minted.

**D21. Menus, not winners.** `[S]`
A filling’s CONTRACT lists every legal image. Deleting
an image because another “won” a bake-off is forbidden
unless its failure mode is empty *and* its budget is
strictly worse on every claim. Usually both stay. In
the registry an image row may move to `refused`; it
may never silently vanish.

**D22. Complementary duals.** `[S]`
Frequency-domain oracles (TMM, HB, eigensolve) do not
have to time-step. Time-domain images do not have to
replace them. Ask lock questions in `ω`; ask attacks
in `t`.

**D23. Bind before you build.** `[S]`
Every image on every menu names its **owner**: either
an existing `crate::API` or “NEW: extends `<crate>` at
`<seam>`”. A new crate needs a displacement argument —
which existing crate was tried, and why extension fails.
This plan is the binding record; §7.8 and §10 are its
inventory. The failure mode this rule exists for is
real: the tree already carries Hertz three times and
Hunt–Crossley twice (§10).

**D24. Consolidate before new physics.** `[S]`
When a menu image touches a surface that exists more
than once (contact laws, ZOH runtimes, WAV encoders,
plate assemblers, port containers, steppers), the image
binds to the owner named in §10. Adding an Nth copy is
a refusal, not a style choice.

**D25. The registry is the claim.** `[S]`
D19/D21 are enforced by `instrument-claims.json` plus an
xtask gate (§8), mirroring `capability-maturity.json`
practice. A ◆ live-default without a green gate *and* a
measured budget row is a policy violation the gate
refuses.

---

## 3. Exactness classes

Borrowed from the mega-fused *tier* idea, applied to
physics — not to bit-identity of audio.

| Class | Meaning | Example |
| --- | --- | --- |
| **X-Exact** | Discrete solution of the named operator, up to rounding | Modal ZOH of a linear oscillator; shift of lossless `w^±` |
| **X-Consist** | Converges under documented refinement | Fine LC / fine multimodal TMM |
| **X-Struct** | Passive/structure-preserving approximation with residual | Foster ZK; VFIT `R(s)`; IRKA 1-port |
| **X-Est** | Estimate-only (`fs-io` / adapter / uncertified lab) | Raw mesh skeleton; unscored LBM jet |

A performance kernel may be X-Exact for its *limit*
(lossless delay) and X-Struct for the rest (ZK Foster).
The gate states both. We never call a 6-cell ladder
X-Exact for a distributed bore.

---

## 4. Shared kit — with owners

Ports, cards, extractors, TMM (air oracle), eigensolvers
(vibration oracle), rational fitting, psycho judges,
fingerprint ledger, gesture signals. Kernels are
**separate**. No `k_everything`. Per D23, the kernels
now have names *in the tree*:

| Kernel | Owner today |
| --- | --- |
| `k_tmm` | `fs-duct` (ZK/Bessel losses, cones at local radius, tone-hole T-junctions, Levine–Schwinger / flanged mouths) — live |
| `k_char` | `fs_vfit::DelayedFilter` (exact-FIR scattering port, passivity-enforced) ⊕ `fs_couple::driving_point::characteristic_line` (FIR *is* the TMM `R(ω)`) — live |
| `k_lc` | `fs_phs` ladders / `acoustic_chain*` + `fs_phs::step` — live |
| `k_modal` | `fs-modal` (inertia-certified `slice_window`) + `fs_couple::modal_acoustic_time` (exact-ZOH audio runtime) — live |
| `k_gonzalez_red` | `fs-nlmodal` (SOS quartic storage on `fs_phs::step`) — live |
| `k_vfit` | `fs-vfit` (relaxed VF ⊕ Loewner cross-check, Hamiltonian-exact passivity certificate, convex residue repair) — live |
| `k_valve` | `fs_phs::{bernoulli_volume_flow, quasistatic_aperture_opening}` + 1/2-DOF `mass_spring_damper` — live |
| `k_felt` | `fs-material` (`Uniaxial` + `FractionalZener`→Prony) ⊕ `fs-contact` `FiniteGapPoint` geometry ⊕ island on `fs_phs::step` — trait/visco/gap live; wool-felt implementor NEW (bead 87zbd) |
| `k_fric` | `fs-tribo` (Stribeck, partial slip) + `fs_couple::stribeck_friction` — live |
| `k_jet` | `fs-aeroac` jet-labium lab over `fs-lbm` — live in 2-D (tonal only, typed broadband refusal); 3-D operator staged, no sweep |
| `k_ckt` | `fs_phs::DescriptorPortHamiltonian` + `step_descriptor` — live core. Constrained Kirchhoff DAE is fs-phs’s *recorded deferral with a named trigger* (“first consumer needing constraints”); T-EM **is** that consumer. Device cards NEW |
| `k_hb` | NEW. Home is the periodic-orbit facility already reserved by `fs-vmanifest::i09` (shooting / collocation / harmonic balance) — it also serves flutter and thermoacoustics, so it is **not** a music-local kernel |

Also owned, and previously uncredited:

- **Measurement loop**: `fs-modalid` (FRF → `(f, η, φ)`,
  stabilization, MAC) is the simulate-vs-measure half of
  every realism gate (§8).
- **Psycho judges**: `fs-psycho` (ISO 532-1 loudness,
  DIN sharpness, Daniel–Weber roughness, ECMA tonality,
  fluctuation strength) under its `LISTENING_LAW`.
- **Gesture**: typed schedules extend `fs-scenario`
  (its `acoustic.rs` already carries the aperture
  scenario). Gesture is geometry in time, not MIDI.
- **Convention seam** (pinned, load-bearing): acoustics
  crates (`fs-duct`, `fs-bem`) use `e^{−iωt}`; `fs-vfit`
  fits `e^{+iωt}`. Conjugate before fitting — the
  clarinet casebook pins the failure (no stable rational
  model fits unconjugated data).

---

## 5. Menus per filling

Each filling: **cards**, **images** (orthogonal),
**selector hints**, **failure modes**, **owner** (D23).
Images marked ◆ are the usual live default, not the only
legal one. “live” = implemented and tested in the tree
today; “NEW” = greenfield with the named extension seam.

### 5.1 Single-reed / double-reed winds

**Cards:** `A(x)` or cylinder list, hole table, cane/lay,
gas, mouth flange.

| Image | Class | Good for | Failure mode | Owner |
| --- | --- | --- | --- | --- |
| **TMM** | X-Consist / X-Struct (ZK) | Peaks, design, Ernoult, R(ω) oracle | No attack, no live `σ(t)` | `fs-duct` — live |
| **Char + lumped holes + ZK Foster** ◆ | X-Exact delay + X-Struct loss/holes | Played phrases, slurs | Bad if junctions don’t match TMM R | `fs_vfit::DelayedFilter` ⊕ `fs_couple::{driving_point, reed_bore}`; holes via `fs_phs` side-hole lengths — live |
| **Fine LC/pHS + `step`** | X-Consist | Short chambers; TD oracle vs Gonzalez | HF dispersion if `Δx` large (disclose `f_Ny`) | `fs_phs::{lc_ladder, acoustic_chain*}` + `fs_phs::step` — live |
| **VFIT/IRKA 1-port of TMM R** | X-Struct | Held fingering, cheap sustain | Live `σ`; out-of-band | `fs-vfit` — live (clarinet casebook) |
| **HB + TMM Z(nω)** | X-Struct | “Will this reed speak?” | Transients, multiphonics | NEW: i09 periodic-orbit facility |

**Selector:** `dσ/dt > ε` → char; settled → may hop VFIT
(lift). Design CLI → TMM. Studio bounce may use fine
pHS if `f_Ny` covers the claim.

**Do not:** FIR keymap; Newton the bore; treat 6-cell
`f_Ny` as 8 kHz truth.

### 5.2 Brass

**Cards:** flare `A(x)`, mouthpiece, lip pair, `Z_L` table.

| Image | Class | Good for | Failure | Owner |
| --- | --- | --- | --- | --- |
| **Multimodal TMM** | X-Consist | Impedance, cutoff, design | No lips in time | NEW: `fs-duct` MM expansion — its own recorded deferral; bead zolja |
| **HB / describing function** ◆ lock | X-Struct | Slot, pitch lock | Attack, lipping | NEW: i09 periodic-orbit facility |
| **MM characteristic lines + lip island** ◆ play | X-Exact per mode + island | Attack, slurs, slide | If plane-wave only, dull / won’t lock | lines NEW (zolja); island = `fs_phs::bernoulli_volume_flow` + msd — live |
| **BEM-baked `Z_L(ω)`** | X-Struct / X-Est | Unflanged mouth | In-loop BEM | solver live (`fs_bem::helmholtz` Burton–Miller, batch solves, convention-matched); NEW: sweep driver + `fs-duct Termination::Tabulated` (zolja) |
| **Plane-wave char (5.1)** | — | Debug / mute | **Not** a trumpet claim | live (5.1 stack) |

Plane-wave and multimodal **both stay on the menu**.
The trumpet *claim* requires MM. A brass mute study
might use plane-wave + table. Honest status: MM does
not exist yet, so **no trumpet claim is currently
possible** — that is exactly what bead zolja unblocks.

### 5.3 Jet-pipes (flute, recorder, flue, organ rank)

**Cards:** pipe `A(x)`, chimney, jet-card from lab.

| Image | Class | Good for | Failure | Owner |
| --- | --- | --- | --- | --- |
| **LBM/aeroac lab** | X-Est→X-Struct when residual issued | Mint jet card | Not RT | `fs-lbm` + `fs-aeroac::jetlab` — live 2-D (tonal; broadband is a typed refusal); D3Q19 operator staged, no 3-D sweep |
| **Jet island + char pipe** ◆ | X-Struct + X-Exact delay | Flute / recorder phrase | No card → refuse | pipe = live char stack; jet-card format NEW |
| **TMM pipe only** | X-Consist | Cut/tone-hole design | No source | `fs-duct` — live |
| **N independent 5.3 voices** | — | Organ rank | One fat 3-D flue | composition of live parts |

Organ and flute share kernels, **not** a shared “wind
stepper” with reeds. Reed images are off this menu.

### 5.4 Linear plucked / struck strings + body

**Cards:** μ, T, EI, plate DKT/orthotropy, cavity.

| Image | Class | Good for | Failure | Owner |
| --- | --- | --- | --- | --- |
| **Modal ZOH string** ◆ | X-Exact (truncated basis) | Pluck, harp, ordinary guitar | Missing high modes / tension change | `fs_couple::modal_acoustic_time` — live |
| **Eigensolve plate + modal ZOH** ◆ | X-Exact truncated | Body radiation | Attack thump if N too small | `fs-plate` (DKT, prestress, stiffeners) → `fs-modal::slice_window` — live |
| **VFIT of plate driving-point** | X-Struct | Cheap body as 1-port on the bridge | Nonlinear top | `fs-vfit` — live |
| **Fine pHS string (many modes)** | X-Consist | Authority vs ZOH truncation | RT polyphony | `fs_phs::modal_bank` — live |
| **Dispersive waveguide** `[F]` | X-Struct | Alternative string (stiffness on the delay) | Bake-off vs modal; don’t delete modal | NEW: extends `fs_vfit::DelayedFilter` (dispersive allpass) |

Modal and (optional) dispersive waveguide are
**orthogonal string images**. Bake-off on inharmonicity
and a pluck spectrogram; **keep both** if both pass and
budgets differ (polyphony vs one hero string).

**Off menu:** air-column characteristics, 48 kHz DKT.

### 5.5 Nonlinear string / plate

| Image | Class | Good for | Failure | Owner |
| --- | --- | --- | --- | --- |
| **KC / von Karman Gonzalez on N modes** ◆ | X-Consist on the reduced pHS | Loud guitar, gong-ish | N too small; mesh-scale folds | `fs-nlmodal` — live (FE-sampled modes injected via `fs_couple::thin_plate`) |
| **Linear modal (5.4)** | X-Exact truncated | Soft dynamics | Large-amplitude pitch glide | live (5.4) |
| **Full-mesh NL** `[M]` | X-Consist | Offline crash | RT | NEW offline lab |

Selector: `|Δℓ/ℓ|` or energy in cubic terms vs a
threshold. Soft notes may stay on 5.4 in the same
voice (two images, one filling).

### 5.6 Bowed

| Image | Class | Good for | Failure | Owner |
| --- | --- | --- | --- | --- |
| **Modal + Stribeck island** ◆ | X-Exact modes + X-Struct friction | Violin-class mechanism | Not a measured rosin unless card says so | `fs-tribo` + `fs_couple::stribeck_friction` — live |
| **KC + friction** | P-NLH + island | Loud / high bow | Cost | compose `fs-nlmodal` × `fs-tribo` (parts live) |
| **Forced Helmholtz oscillator** | — | **Forbidden** (not physics) | Always | — |

Determinism note: `fs-tribo` declares **no cross-ISA
bit-stability** (platform `hypot`). Bowed-image registry
rows must not claim cross-ISA replay until that is
routed or re-declared (§8).

### 5.7 Piano-class

**Cards:** string `(μ, T, EI)` from steel + tension as
state; soundboard / rim from spruce `OrthotropicElastic`
(already in `fs-material`); hammer **felt** from the
existing hysteresis stack, not a new constitutive
family — see §7.5.

| Image | Class | Good for | Failure | Owner |
| --- | --- | --- | --- | --- |
| **Modal strings + Uniaxial/visco felt island + modal board** ◆ | X-Exact + X-Struct felt | One note / polyphony | Felt card wrong or out of visco band | trait + visco + `FiniteGapPoint` live; NEW: wool-felt `Uniaxial` implementor + felt pack (bead 87zbd) |
| **3-string unison coupling** | X-Consist | Chorus, aftersound | Extra states | compose (parts live) |
| **Linear hammer spring** | X-Exact | Debug | No hammer spectral tilt | trivial — live |
| **Hunt–Crossley hammer** | X-Struct | Debug / bake-off vs hysteresis | Residual-strain / tilt of real felt | `fs_contact::NormalPatchLaw::HuntCrossley*` — live. **Reuse this one; a third HC implementation is refused (D24)** |
| **Duplex as extra modal segments** | X-Exact | Aftersound | — | live modal segments |
| **JCA/Biot felt-as-absorber** | X-Struct | Case linings, under-string cloth | **Not** the hammer contact | NEW — one taxonomy row exists (`docs/MATERIAL_PROPERTY_TAXONOMY.md`), **zero code, zero packs** |
| **Air-column / char tube** | — | **Off menu** | Always for this filling | — |

The live default felt image **consumes physics already
in the tree**:

1. `fs-material::Uniaxial` — the hysteretic fiber
   contract (committed reversal state, consistent
   tangent, cycle dissipation > 0, gated by the
   existing mt-001 / mt-004 conformance tests). Piano
   felt is a **new implementor** of that trait
   (wool-felt envelope + residual crush), not Mander
   concrete and not Menegotto–Pinto steel wearing a
   felt hat.
2. `fs-material::visco` — `FractionalZener` is the
   fitting-side law the crate already names for
   **wood, polymers, and felt** (near-constant loss
   factor). `lower_to_prony` emits a
   `GeneralizedMaxwell` with a certified band; the
   island steps that Prony series (the step is already
   `expm1`-guarded for audio-rate `dt` against long
   relaxation times). Out of band refuses
   (`FS-MAT-VISCO-OUT-OF-BAND`).
3. `fs-contact` — the hammer–string *geometry* of the
   contact: `FiniteGapPoint` sampled from the felt
   thickness field via the existing
   `sample_finite_gap_from_chart`. Hunt–Crossley is
   the **debug / second image**, not the felt law —
   and it is fs-contact’s existing implementation.
4. Taxonomy row (JCA / Biot: flow resistivity,
   porosity, tortuosity, characteristic lengths) —
   felt as a **porous absorber** (linings). A
   different claim from hammer compression, and
   currently pure greenfield.

Hunt–Crossley stays on the menu for bake-off against
a hammer-velocity coupon. It does **not** replace
the hysteresis stack. A linear spring is debug only.

### 5.8 Bars / plates / bells

| Image | Class | Good for | Failure | Owner |
| --- | --- | --- | --- | --- |
| **1-D beam modal** ◆ | X-Exact truncated | Xylophone, bar | 3-D bell | compose: analytic ω (`fs_nlmodal::prestressed_beam_omega`) + ZOH runtime (parts live) |
| **DKT plate modal** ◆ | X-Exact truncated | Gongs, soundboards | Deep shell | `fs-plate` — live (Leissa-gated) |
| **Axisymmetric shell** `[F]` | X-Consist | Bells | Cost, mesh | NEW: the flat-facet-shell follow-up `fs-plate`’s own CONTRACT records — **not** a new crate, and **not** `fs-solid::shell` (that is an estimate-only Euler-disc surrogate, §10) |
| **Strike island** | X-Struct | Attack | — | `fs-contact` / `fs-dcontact` — live |
| **Char air** | — | Off menu | — | — |

Beam vs plate vs shell are **different charts**, not
refinements of one stepper.

### 5.9 Voice (ooh / aah, then more)

| Image | Class | Good for | Failure | Owner |
| --- | --- | --- | --- | --- |
| **TMM of `A(x)`** | X-Consist | Formant oracle | No glottis in time | `fs-duct` — live |
| **KL / char tract + glottal island** ◆ | X-Exact delay + island | Sung vowels, articulation | 3-D mouth, fricatives | substrate live: char stack ⊕ `bernoulli_volume_flow` ⊕ `fs-dcontact` fold collision; articulation foothold `DelayedFilter::modulate_delay`. Voice *product* is NEW — zero glottis/formant code exists |
| **VFIT 1-port of tract R** | X-Struct | Held vowel cheap | Moving tongue | `fs-vfit` — live |
| **Two-mass glottis** vs **1-DOF** | two P-Valve images | Bake-off on spectra | — | 1-DOF parts live (`fs_phs` msd + aperture); two-mass NEW composition |
| **WallPin tissue** | X-Struct | Bandwidth / extra formants | — | `fs_phs::WallPin` — live |
| **Jet card in tract** | 5.3 add-on | Fricatives later | Don’t block /u//a/ | NEW (5.3 card) |
| **Formant filter / vocoder / WAV** | — | **Forbidden** | Always | — |

1-DOF and two-mass **both stay**. Articulated KL and
held VFIT **both stay**.

### 5.10 Electric (pickup → circuit → speaker)

| Image | Class | Good for | Failure | Owner |
| --- | --- | --- | --- | --- |
| **Lumped Faraday** ◆ | X-Struct | Volt from `v_string` | Magnet geometry | NEW: EM port module in `fs-couple` sitting on a 5.4/5.5 port — zero pickup code exists today |
| **Magnetostatic bake** `[F]` | X-Struct | Position-dependent B | Cost | NEW offline lab |
| **Split circuit (linear factor + device)** ◆ | X-Consist linear / X-Struct device | Distortion as constitutive | Stiff device without island | core live: `fs_phs::DescriptorPortHamiltonian` + `step_descriptor`. Kirchhoff-law DAE interconnection is fs-phs’s **named deferral trigger**, and this track is the first consumer. Device (triode/diode) cards NEW — zero device code exists |
| **Thiele–Small / plate cabinet** | X-Struct | Speaker + box | Room | compose `fs_phs` msd / transformer / radiation pieces + `fs-plate` (parts live) |
| **Cabinet IR pack** | — | **Forbidden** | Always | — |

String images are **whatever 5.4/5.5 selected**. The
electric path does not re-discretize the string.

### 5.11 Cross-filling graphs

String×plate×cavity: pick 5.4/5.5 string image × 5.4/5.8
body image. Reed×plate: 5.1 char/VFIT × 5.4 modal.
Pickup sits on a 5.4/5.5 port. **No fused interior.**
All joins go through `fs-couple` (`join_port`,
`PortSchema` v2, energy/window audits) — the crate that
already owns the composed instrument loop.

---

## 6. Selector and bake-off

### 6.1 Selector (runtime)

```text
select(F, gesture, claim, budget) → {I_k}
```

The selector v1 is a **module in `fs-couple`** (the only
crate that already composes images into audible loops),
consuming the claims registry (§8) and the measured
budget rows (§9). It is not a new crate and not physics.

Typical policies (not laws):

- Phrase / `dθ/dt` large → spatial or island images.
- Held / sustain → VFIT or modal 1-port if gated.
- “Where are the peaks?” → TMM / eigensolve, no audio.
- “Will it speak?” → HB/describing function.
- CPU panic → drop to a *pre-gated* cheaper image of
  the **same filling**, never to a forbidden one.

### 6.2 Bake-off (offline, required to add an image)

Same cards, same fixture, same QoI. Report residual,
budget, failure mode. Outcomes: **keep both**, keep
one for a subset of claims, or refuse the newcomer.
“Keep both” is the default when they are orthogonal
(TD vs FD, phrase vs hold, 1-DOF vs two-mass). Every
bake-off outcome is a **receipt referenced by the
registry row** (§8) — not a paragraph in a doc.

### 6.3 Hopping mid-note

Allowed if D17 lift exists and a G7-like fixture
shows no click. Example: char slur → VFIT sustain.
Forbidden: hop to an image that failed that filling’s
gate.

---

## 7. From 3-D geometry and materials to sound

The 3-D trumpet is not what we time-step. It is the
**source of every number** the physics is allowed to
play.

We do **not** drop a photogrammetry mesh into a 48 kHz
Navier–Stokes solver, paint it “brass,” and hope a Bb
comes out. That is computationally impossible on an
audio thread and the wrong reduction. A real trumpet
is already a reduced machine: a long, slowly flaring
air column, a couple of electrically short cavities,
three switchable extra tubes, a radiating bell, and a
pair of soft orifices that lock to that column. The
job of the mesh and the material cards is to **mint
those numbers honestly**, then let the named operators
in §5 produce the sound.

The same sentence, with different charts, is how we
get every other instrument.

### 7.1 The ingest law

```text
licensed mesh(es) + region labels
        + matdb / GasState cards
        + gesture map
              │  fs-io quarantine (Estimate-only)
              ▼
     extractors  →  one or more *charts*
              │  certificates (volume, A(x)>0, …)
              ▼
     authority images (TMM, eigensolve, BEM table, lip FEM reduce)
              │  bake-off (§6)
              ▼
     performance images (char, modal ZOH, VFIT, island)
              │  selector(gesture, claim, CPU)
              ▼
           pascals
```

Three things never enter the audio callback:

1. The triangle mesh.
2. A 3-D lip / reed / felt FEM.
3. A “trumpet brightness” EQ, or any other invented
   timbre.

What *does* enter is a **chart**: a certified,
low-dimensional description of the same object.
`A(x)`, valve lengths, plate mid-surface, lip reduced
masses, a `Z_L(ω)` table, a felt `Uniaxial` + Prony
card. If the extractor cannot certify the chart, that
path is X-Est or refused. We do not invent numbers to
keep the note playing.

**Extractor owner and authority colors** (corrected in
v8). The mesh → `A(x)` extractor is **NEW code in
`fs-query`** — nothing closer exists. Its footholds are
already there: `medial_poles` (interior pole cloud with
radii — the medial-axis approximation), `thickness_at`,
and the rigorous `geometric_moments` enclosures. Honest
labels per fs-query’s own contracts: the centerline /
inscribed-radius / thickness machinery is
**Estimate-class by construction** (the implicit marcher
has no no-tunneling theorem); only the moment/volume
enclosures are certified. So an `A(x)` receipt is
**Estimate, cross-checked by pole radii and thickness,
closed by a certified volume enclosure** — never sold as
a certified chart.

Materials are **cards on labeled regions**, not a
vibe. They come from `fs-matdb` receipts into
`fs-material` laws that already exist. Credit where
due (v8 correction): the seed corpus already holds
**≈99 provenance-complete packs**, including yellow
brass with measured damping, four spruce packs (two
with measured loss factors), maple / rosewood /
mahogany tonewoods, phosphor bronze, **music wire**,
and an instrument grain-axis convention. What is
actually missing is narrow: **felt, cane, tissue,
gut**. Absent data is a named refusal, and that
refusal is the population signal.

| Region | Card (examples) | What it is *not* |
| --- | --- | --- |
| Air inside / outside | `GasState` from `GasSpec` + T, p (`fs-material::gas`, already live). **Humidity is NOT an input**: `c = sqrt(γRT)`; RH enters only ISO 9613 absorption. A moist-air `GasSpec` mixture rule (variable `M`, `γ`) is NEW work | “trumpet air” |
| Yellow-brass body | isotropic `E, ν, ρ` via `resolve_isotropic_elastic_state_point` (packs exist: C26000 / yellow-brass damping) | a brightness knob |
| Lip / fold tissue | tissue `ρ, E` (or Lamé), loss factor — **pack NEW** | cane reed, or a mass typed by ear |
| Cane reed | orthotropic cane — **pack NEW** | a generic “reed” |
| Spruce top / soundboard | `OrthotropicElastic` via the existing `resolve_orthotropic_elastic_state_point` (spruce packs exist, two with measured loss) | a body IR |
| Piano-hammer felt | **new `Uniaxial` implementor** + `FractionalZener` → `GeneralizedMaxwell` (the visco crate already names felt); **coupon pack NEW** | Hunt–Crossley as the felt; Mander concrete in a felt hat |
| Felt lining / cloth | JCA / Biot porous-absorber row — **one taxonomy row, zero code today** | the hammer contact |
| Steel string | `ρ, E`, and **tension as state** (music-wire packs exist) | a sampled string |

Gesture is **geometry in time**, not a MIDI note
number: valve down → extra tube length; key down →
hole `σ`; lip pressure → orifice and pre-stress; bow
`v` and `F_n`; hammer velocity into the felt card.
Typed gesture schedules extend `fs-scenario`.

### 7.2 Worked filling: a brass trumpet

There is **no `Trumpet` type**. The instrument is a
composition of nameless charts.

**What we take in**

1. **A real interior.** A cited, licensed watertight
   mesh of the *air*: mouthpiece cup → throat →
   leadpipe → valve cluster → yards → flare → bell.
   Separately, the **exterior** of the bell (that is
   what radiates). Valve pistons in both states, or at
   least the **crook centerlines** from the CAD, so we
   know how much extra tube each piston inserts.

   The valve *casing*, the engraving, the lacquer, the
   serial number, and the photoreal shading do **not**
   enter the wave field unless we make an explicit
   leak / rattle claim. Acoustically we need the lumen
   and the bell exterior. Looks-like-a-trumpet for
   LUMEN is a different pipeline.

2. **A brass card.** Density, Young’s modulus, Poisson
   ratio of the actual alloy (the packs exist). Used
   for:
   - optional wall yield (`WallPin` in `fs-phs`:
     resistance, mass, stiffness per area, then
     `Y' = 2πa · slant / Z'` in the duct),
   - optional bell-shell modes `[F]`,
   - **not** the speed of sound in the air. That is
     `GasState`. Brass does not make the air
     brass-colored.

3. **Lips that are shaped like lips and made of lip.**
   Two inputs, both serious:
   - a **lip-shaped mesh** (rest gap, thickness,
     width, the actual orifice),
   - a **tissue card** (density, elasticity, loss —
     lip, not cane, not rubber; pack NEW).

   Those two go into an **offline authority lab**:
   `fs-solid::linear3` (3-D small-strain FEM whose
   output is the ordinary `(K, M)` pencil) →
   `fs-modal` → a NEW reduction recipe. That lab
   *mints a valve card*: reduced masses, stiffnesses,
   damping, and an orifice law whose rest shape came
   from the mesh. Same moral status as the LBM jet:
   it is a 3-D experiment that produces a card. The
   triangle mesh does **not** run at 48 kHz.

4. **Air.** `GasState` from T and p. Hot stage vs cold
   rehearsal is `c(T)`, which warps every delay and
   every cutoff. (Humid breath needs the NEW moist-air
   mixture rule — do not pretend the current card
   models it.) Lungs are an upstream compliance plus a
   blowing-pressure gesture, not a “note-on.”

5. **Player gesture.** Three bits (1st / 2nd / 3rd
   valve) plus continuous lip pre-stress / aperture
   plus blow pressure. Optional tuning slide as a
   continuous length.

**What the extractor does (not a `Trumpet` type)**

| 3-D thing | Chart it becomes |
| --- | --- |
| Interior lumen | Centerline + inscribed radius → `A(x)` samples (fs-query extension; Estimate + certified volume closure) |
| Mouthpiece cup | Electrically short → **lumped** sections / Helmholtz-ish cavity (D18) |
| Cylindrical yards | **Characteristic** delays of length `ℓ/c(T)` |
| Flare | `A(x)` + **local cutoff** of higher modes `f_c(s) ∼ c j_{mn} / (2π r(s))` → P-MMWave lines that exist only where they propagate |
| Valve *up* | Bypass length `L_i` |
| Valve *down* | Bypass + **crook** length `L_i + ΔL_i` taken from the CAD of that slide |
| Bell exterior | Offline `fs-bem::helmholtz` sweep → `Z_L(ω)` table → `fs-duct Termination::Tabulated` (bead zolja) |
| Brass wall | Optional `WallPin` from brass `σ, K, r` and local `a` |
| Lip mesh + tissue | Offline `fs-solid`→`fs-modal` reduce → `(m1,m2,k,c, orifice)` valve card |

Certificates we actually check, with their honest
colors: reconstructed air volume vs the mesh
(**certified enclosure** via `geometric_moments`);
`A(x) > 0` everywhere and centerline sanity
(**Estimate**, cross-checked by medial poles /
thickness); each valve `ΔL` matches the CAD within a
band; `Z_L` is passive (**fs-vfit certificate class
recorded** — Hamiltonian-exact or grid-only, named).
Fail a certificate and that image is Estimate-only or
refused.

**The playing loop (this is where the sound comes from)**

1. Lungs hold a blowing pressure behind the lips.
2. The lips (tissue card, lip-shaped rest orifice)
   form a pressure-controlled valve. Bernoulli through
   the gap. The gap opens and closes because the
   tissue has mass and stiffness.
3. Downstream, the air column presents a reflection
   `R(ω)` at the mouthpiece. That `R` is a function of
   `A(x)`, the valves, the flare, the bell `Z_L`, the
   wall card, and `GasState`.
4. The lip oscillator **locks** to a harmonic of that
   column — the one whose impedance peak sits where
   the lips can still close. That is why a trumpet
   plays a harmonic series, why you can lip up/down
   inside a slot, and why a valve that lengthens the
   column drops the whole series.
5. What you hear is the **radiated** field at the
   bell, from the same `Z_L` that terminated the
   column. Not a sampled trumpet.

HB / describing function (menu 5.2) can predict the
lock offline. The time-domain image produces the
attack. They are orthogonal, not rivals.

**How depressing a “key” makes a different note**

A trumpet “key” is a **piston that switches which
brass tube the air sees**. It is not a Boolean synth
parameter and it is not a MIDI note number.

When the first valve goes down, the extractor (or the
already-baked crook table) **inserts ~16 cm of extra
tube** into `A(x)`. The authority TMM rebuilds the
input impedance the lips see. The performance image
updates three **control-rate** lengths (or three
bypass junctions) and lifts state (D17) so there is
no click. Pitch is not assigned. Pitch **emerges**
because `R(ω)` moved and the lips locked to a new
harmonic of that column.

**How it sounds like a trumpet rather than a cone**

- **Cup + leadpipe geometry** → the high-frequency
  `R(ω)` the lips see. That is attack and slot, and
  it is why a Bach 1-1/2C and a Schilke 14A4a are
  different instruments on the same body.
- **Flare + higher modes**, cutoffs from `r(x)`, not
  a brightness knob. Higher-mode lines exist only
  where they propagate.
- **Brass** is in the wall / optional bell shell. It
  is a small, honest yield and a possible shell
  color. It is not “add even harmonics.”
- **Lips** are tissue-card dynamics plus Bernoulli
  through an orifice whose rest shape came from a
  lip mesh.
- **Air** is `GasState`. Change temperature, every
  delay and every cutoff moves.
- **Valves** change length, hence peaks, hence lock.

A plane-wave column plus lips can still produce *a*
note. We **do not claim** that note is a trumpet.
The trumpet claim requires the flare’s higher modes
and a mouthpiece cup that is actually that cup.
Plane-wave stays on the menu as debug / mute
(menu 5.2). It is not the instrument.

**What we refuse (honesty)**

- 48 kHz 3-D CFD of the mouth and bell.
- A lip triangle mesh in the audio thread.
- “Brass = add even harmonics.”
- MIDI note number as the physical pitch.
- Inventing a 3-D jet or a mean-flow `k/(1−M²)` hack
  to make it “more 3-D.”

### 7.3 The same law for every filling

Always: **mesh + region labels + material cards →
charts → menu images → selector.** The *chart* changes.

| Filling | 3-D / geometry in | Material cards | Gesture → geometry | What emerges |
| --- | --- | --- | --- | --- |
| **Trumpet / horn** | Interior lumen, crooks, exterior bell, lip mesh (offline) | Brass, tissue, air | Valves = extra `ΔL`; lip stress; blow `p` | Lock pitch, flare brightness |
| **Trombone** | Same + slide yards | Brass, tissue, air | Slide = continuous `L(t)` + lift | Gliss, lock |
| **Clarinet / sax** | Bore + mouthpiece + tone holes + pads; reed/lay mesh offline | Cane or synthetic, plated body, air | Key = hole `σ` / chimney; reed blow `p` | Fingering peaks, reed lock |
| **Flute / flue** | Headjoint + chimney + body holes; labium edge | Metal/wood, air | Jet speed / angle; keys = `σ` | Edge-tone lock (needs jet card) |
| **Voice /u/ /a/** | MRI or published `A(x)`; optional face mesh for radiation bake | Tissue `WallPin`, air | Tongue/lips = `A(x,t)` | Formants; sung interaction via `R` at glottis |
| **Guitar (acoustic)** | Top/back/sides mesh; strings as 1-D; f-holes; cavity volume | Spruce ortho, steel/nylon, air in box | Fret = shorten `L` + obstacle height; pluck | Partials, body mobility |
| **Guitar (electric)** | String 1-D + pickup pose + **schematic** (electrical geometry) + cab mesh | Steel, magnet bake, device cards, cone | Fret, pluck; knobs = circuit params | Volt then cone `p` |
| **Violin** | Corpus mesh → plate/cavity charts; string 1-D; bridge as transformer | Maple/spruce, gut/steel, rosin constitutive | Bow `v`, `F_n`, station; finger = `L` + obstacle | Stick-slip + body |
| **Piano** | Soundboard / rim mesh; hammer felt thickness field; string lengths | `Uniaxial`+Prony felt, steel, spruce `OrthotropicElastic` | Key `v` into felt island; pedal = coupling | Spectral tilt, aftersound |
| **Bell / bar** | Bell scan or bar CAD | Bronze / rosewood | Strike position + `v` | Partial ratios from eigensolve |
| **Organ rank** | Many 5.3 pipes from a chest drawing | Wood/metal, air | Stop = which pipes get the jet card | Chorus of cheap pipes |

One physical object can emit **several charts** (guitar:
strings *and* top *and* cavity). Each chart feeds only
the images on **its** menu.

### 7.4 Lips, reeds, felt, fingers — “shaped like X, made of X”

We take the 3-D shape and the tissue/cane/felt card
**seriously in the authority lane**, then **reduce**:

| Player bit | Authority (3-D + material) | Card that plays |
| --- | --- | --- |
| Lips | `fs-solid::linear3` FEM / beam of a lip-shaped mesh, tissue `E,ρ,ζ` → `fs-modal` reduce; orifice from rest gap | Two-mass / beam + Bernoulli |
| Cane reed | Shell/plate of a reed blank, cane ortho; lay as fs-contact patch law from lay mesh | 1-DOF or measured MSD + lay |
| Vocal folds | Same idea as lips, glottal gap from mesh | 1-DOF and/or two-mass (both on menu) |
| Piano felt | Thickness field + wool-felt coupon on the **existing** `Uniaxial` + `FractionalZener`/`GeneralizedMaxwell` stack | Substepped hysteretic island (Hunt–Crossley is debug only) |
| Fingertip / fret | Digitizer height field + soft contact card | Obstacle + optional Coulomb |
| Bow hair | Optional 1-D hair `[F]`; else station + Stribeck card | Friction island |

If we skip the 3-D reduce and type masses by hand, that
is **Estimate-only** until a coupon or a mesh-reduce
receipt exists. We still do not put the lip FEM — or
the felt continuum — in the callback.

### 7.5 Worked filling: a piano, felt on physics we already have

Same ingest law as the trumpet. Different charts.
Still no `Piano` type.

**What we take in**

1. **Soundboard + rim mesh**, labeled spruce / maple /
   rib regions. Card: `OrthotropicElastic` already in
   `fs-material` (nine engineering constants resolved
   atomically by the existing
   `resolve_orthotropic_elastic_state_point`; grain
   axes owned by the geometry consumer — the matdb
   axis-convention file exists). Authority: `fs-plate`
   eigensolve → modal chart. Performance: modal ZOH of
   those modes (menu 5.7 / 5.4).
2. **Strings as 1-D.** Lengths, speaking points,
   duplex segments from the plate/bridge geometry.
   Card: steel `ρ, E` (music-wire packs exist), and
   **tension as state** (not a frozen MIDI frequency).
   Image: modal ZOH (menu 5.4); 3-string unison as a
   second image.
3. **Hammer felt that is felt.** A thickness field
   (possibly from a scan or a manufacturer section)
   plus a **wool-felt coupon**. This is where we
   refuse to invent a new constitutive family.

   FrankenSim already captured the physics felt
   needs:

   | Already in the tree | Role for the hammer |
   | --- | --- |
   | `fs-material::Uniaxial` (hysteretic fiber contract: reversal state, consistent tangent, cycle dissipation, FD-gated by mt-001 / mt-004) | Felt compression is uniaxial with residual crush. New implementor: wool-felt envelope, *not* `ManderConcrete` or `MenegottoPintoSteel`. Same trait, same thermodynamic gates. |
   | `fs-material::visco::FractionalZener` | Fitting-side law the crate **already names for felt** (near-constant loss factor over decades). |
   | `lower_to_prony` → `GeneralizedMaxwell` | Runtime island: exact-exponential Prony step (already `expm1`-guarded for audio-rate `dt`), dissipation ledger, out-of-band refusal. |
   | `fs-contact` `FiniteGapPoint` + `sample_finite_gap_from_chart` | Geometry of the hammer–string interface from the felt thickness field — the sampler exists. |
   | `fs_contact::NormalPatchLaw::HuntCrossley{Sphere, EllipticParaboloid}` | **Debug / bake-off image**, not the felt — reuse this implementation (D24). |
   | Taxonomy JCA / Biot row (flow resistivity, porosity, tortuosity, Λ, Λ′) | Felt as **porous absorber** — case linings, under-string cloth. A different claim, and zero code today. |

   The 3-D felt continuum (if we ever run one) is an
   **authority lab** that mints the `Uniaxial` + Prony
   card, same moral status as the lip FEM. It does
   not run at 48 kHz.

4. **Air in the box** is a cavity chart coupled at
   ports to the board, not a characteristic tube
   (off menu 5.7).
5. **Player gesture.** Key velocity → hammer
   incoming `v` into the felt island. Sustain pedal
   couples / uncouples string terminations. Soft
   pedal is a hammer-shift / sample-count change of
   the *same* cards, not a different instrument.

**The playing loop**

1. Key velocity sets hammer kinetic energy.
2. The felt island (`Uniaxial` stress + Prony
   internal variables) meets the string modal
   coordinates through the `fs-contact` gap. The
   force is hysteretic: loading and unloading are
   different paths; residual crush stays in state;
   the dissipation ledger must stay non-negative.
3. That contact force is the only thing that
   *tilts* the string’s partials. A linear spring
   cannot claim a piano attack.
4. The strings drive the board at the bridge port.
   The board radiates. Pedal state changes which
   string terminations are live.

Pitch is the string’s tension and speaking length,
not a MIDI `ω`. Timbre of the attack is the felt
card. Aftersound is unison coupling + duplex +
board. None of those are samples.

**What we refuse**

- Inventing a second felt constitutive next to
  `Uniaxial` / `FractionalZener` because “piano is
  special.”
- Using `ManderConcrete` or `MenegottoPintoSteel`
  as the hammer (wrong envelope, wrong units of
  meaning).
- Treating Hunt–Crossley as the piano-felt law —
  or implementing it a third time.
- JCA/Biot as a substitute for hammer compression
  (wrong physics: that row is absorption).
- MIDI note number as `ω`.
- A 48 kHz 3-D felt FEM in the callback.

### 7.6 Keys, valves, holes, slides, pedals

These are **time-varying charts**, not note-on:

| Player action | Chart change |
| --- | --- |
| Trumpet piston | Insert/remove crook length `ΔL_i` from CAD |
| Trombone slide | Continuous centerline length |
| Woodwind key | Hole `σ ∈ [0,1]` + chimney + pad `C` from the hole’s 3-D |
| Clarinet register | A small vent at a station (another hole) |
| Guitar fret | Speaking length `L` + obstacle at fret height |
| Piano key | Hammer velocity; not a MIDI note as `ω` |
| Sustain pedal | Couple / uncouple string terminations |
| Mute / hand in bell | Change `Z_L` bake (another exterior chart) |

The performance image updates these at **control rate**
and lifts state (D17). The sound changes because `R(ω)`
or the modal lengths changed — not because a sample
switched.

A deliberate boundary: gestures stay **chart
schedules**, not mechanism simulations. `fs-mbd`
explicitly disclaims joints, constraints, and impacts,
and the `fs-kinematics` lane it defers to does not
exist. Nobody “fixes” a valve gesture by simulating
the piston mechanism.

### 7.7 What “looks like a real trumpet” does *not* require

Photoreal shading is LUMEN, not this plan. Acoustically
we need the **interior lumen and the bell exterior**,
not the valve *casing* as a solid in the wave field
(unless we claim leak / rattle — then it is another
contact/cavity chart). Engraving, lacquer, and serial
numbers do not enter `A(x)` unless a measurement says
the wall card changed.

A piano that “looks like a piano” needs the
soundboard mid-surface, the speaking lengths, and
the felt thickness field. Ivory, case veneer, and
the music desk do not enter the hammer island.

### 7.8 What is live today vs what this still needs

Live in the tree, verified against the code (v8):

- `fs-duct` — 1-D viscothermal TMM: ZK / all-regime /
  Bessel losses, cones cascaded at local radius,
  tone-hole T-junctions with Dalmont matching and
  chimney wall law, Levine–Schwinger (0.6133 a) and
  flanged (0.8216 a) mouths, `ka ≤ 1` named refusal
- `fs-phs` — Dirac algebra + Gonzalez discrete-gradient
  stepping with the exact energy ledger; the acoustic
  element zoo (cylinders, spherical cones, chains,
  side-hole lengths, `WallPin` / `ViscothermalPin`,
  Foster fitting, modal banks, Bernoulli aperture,
  `DescriptorPortHamiltonian`)
- `fs-vfit` — VF + Loewner, two-tier passivity
  certificate + convex repair, biquads / state space /
  `DelayedFilter` scattering runtime
- `fs-modal` / `fs-plate` / `fs-nlmodal` — certified
  eigenproblem; orthotropic DKT plates with prestress
  and stiffeners (Leissa-gated); von Kármán / KC on
  FE-sampled modes
- `fs-couple` — the composed loops: `reed_bore` (locks
  near the quarter-wave), string–plate–duct joins,
  thin-plate wiring, exact-ZOH modal audio runtime,
  Aitken added-mass iteration, energy/window audits,
  `PortSchema` v2, physically scaled WAV bytes
  (test-only)
- `fs-bem::helmholtz` — Burton–Miller exterior
  radiation, radiation-impedance matrices, batch
  solves, SH directivity, deterministic casebook
- `fs-contact` — Hertz / Hunt–Crossley / FiniteGapPoint
  patch laws with receipts; `fs-dcontact` — distributed
  unilateral lays/rattles as pHS storage; `fs-tribo` —
  Stribeck / partial slip
- `fs-material` — `GasState` (T, p), `Uniaxial`
  hysteresis with mt-gates, `FractionalZener` →
  `lower_to_prony` (names felt), orthotropic resolvers
- `fs-matdb` — ≈99 seed packs incl. brass damping,
  spruce loss factors, music wire, axis convention
- `fs-aeroac` over `fs-lbm` — the 2-D jet-labium lab
  with hysteresis-ramp protocol and the typed 2-D
  broadband refusal
- `fs-psycho` — ISO/DIN/ECMA metrics under the
  listening law; `fs-modalid` — FRF identification
- `Region` / `Chart` in `fs-geom`; `medial_poles`,
  thickness, and certified moments in `fs-query`

Not yet, with owners, and we will not fake them:

- Mesh → `A(x)` extractor — NEW in `fs-query`
  (centerline chaining + plane-section area; Estimate
  + certified volume closure). No centerline or
  cross-section-slicing code exists anywhere today
- Valve-crook `ΔL` from CAD — NEW (extractor sibling)
- Offline lip / reed reduce — NEW recipe over
  `fs-solid::linear3` → `fs-modal`
- Wool-felt `Uniaxial` implementor + cited felt coupon
  pack — NEW (bead 87zbd); felt / cane / tissue / gut
  packs are the missing matdb families
- Multimodal horn — NEW in `fs-duct` (its recorded
  deferral; the trumpet-claim gate; bead zolja)
- `Z_L(ω)` bake plumbing — NEW: sweep driver over
  `solve_radiation_batch` + `Termination::Tabulated`
  (bead zolja). The BEM physics is done
- HB / describing function — NEW (i09 periodic-orbit
  facility; zero HB code exists)
- Voice product — NEW (zero glottis/formant code;
  substrate ports live)
- Electric path — NEW (zero Faraday / Thiele–Small /
  device code; `DescriptorPortHamiltonian` core live)
- Moist-air `GasSpec` mixture rule — NEW
- JCA / Biot poroacoustics — NEW (one doc row only)
- 3-D jet card — blocked on a 3-D sweep (2-D is
  pinned tonal; the D3Q19 operator is staged)
- Honest mean flow — NEW (not a `k/(1−M²)` hack)
- Product output surface — NEW (bead ib15w, §9)

---

## 8. The claims registry — the product

D19/D21/D25 as code, mirroring how this repository
already enforces claims (`capability-maturity.json`,
`schema-policy.json`, vv-scorecard). Bead **mc31g**.

**Artifact.** Tracked `instrument-claims.json`. One row
per `(filling, image, QoI)`:

```text
{ filling, image, qoi,
  owner_crates[], exactness_class,
  gate: ungated | green | refused,
  evidence,            # test path / package / receipt
  determinism_class,   # per-row, honest
  budget_row,          # §9 measured samples/sec ref
  corpus_refs[] }      # §8 validation data
```

**Gate.** `xtask check-instrument-claims` joins
check-all: a ◆ live-default without a green gate *and*
a measured budget row refuses; an image row may move to
`refused` but may never vanish (D21). Bake-off outcomes
(§6.2) are receipts referenced by rows.

**Validation data.** Reuse `fs-vvreg`’s Level-A/B/C
corpus discipline — do not invent a parallel music
registry. Corpora to register (licensing resolved
*before* ingestion — the same law as “Licensed A(x)”):

- Ernoult-class woodwind resonance / peak-cents data
- Published input-impedance datasets (clarinet,
  trumpet) for TMM and char gates
- Leissa plate tables (already fs-plate’s gate) and,
  when licensable, guitar-top holography (fs-plate’s
  CONTRACT records its absence honestly)
- Felt force-compression coupons (hammer-velocity
  spectral-tilt gate)
- Published glottal-flow waveforms for the voice gate

**Measured side.** `fs-modalid` is the loop: measured
FRFs → `(f, η, φ)` + MAC against simulated modes. Its
recorded no-claim — no license-compatible published FRF
*trace* found yet — becomes a corpus-hunt registry row,
not a fabricated benchmark.

**Listening law.** Carried verbatim from `fs-psycho`:
the metrics are never a substitute for human listening.
A registry row may cite psycho metrics as evidence; the
“sounds right” adjudication is a **recorded human
listening receipt** (who, when, what fixture).

**Determinism truth (2026-08-14).** The entire music
stack is one-host bit-deterministic; **zero cross-ISA
goldens exist**; `fs-tribo` explicitly declares no
cross-ISA claim (platform `hypot`) and sits inside the
bowed path. Rows must say so. G0 passivity remains
mandatory for every admitted image.

---

## 9. Budget law and the product surface

### 9.1 No ◆ without a measured budget

Real-time was a founding ask; today the plan carries
zero measured numbers, and the implementations are
correctness-first (`fs_phs::step` is Newton with an FD
Jacobian per step at 1e-13 — the right first move, and
nowhere near a callback budget). Law: a ◆ live-default
requires a **measured samples/sec row** (state count,
machine fingerprint, headroom at 48 kHz) produced by a
roofline-style lane (bead ib15w), recorded in the
registry. “Real-time” is measured headroom, never
prose — the repo-wide performance doctrine applies.

### 9.2 Fusion comes after the gate, kernel by kernel

Per D5 (correctness before *that kernel’s* fusion) and
the mega-fused-kernel doctrine: fusion candidates, in
likely order of payoff — `DelayedFilter` push/incoming
(FIR), the modal-bank ZOH update, the Prony island
step, the aperture solve (bisection → island Newton
with an analytic Jacobian), `fs-dcontact` probes.
Each fusion must reproduce its gated image bit-for-bit
in strict mode or ship as a declared fast mode.

### 9.3 One product surface (T-Out, bead ib15w)

Today pascals→PCM exists **twice and disjointly**:
`fs_couple::pcm_wav` (mono, physically scaled, never
peak-normalizes — and test-only: no binary writes it)
and the fs-euler-disc-e2e cinematic stack (stereo,
48 kHz-pinned, receipt-hashed, disk-written). Law:

- **One owner** for pascals→PCM. Converge at the
  existing h7xu5.7.8 bridge seam; keep the
  never-peak-normalize rule; no third writer.
- `fs-couple` exposes a **block render API**:
  audio-rate blocks, control-rate chart/gesture
  updates between blocks (gesture as fs-scenario
  schedules), no allocation inside a block. The
  real-time callback is a consumer of that API.
- A CLI renders a gated assembly to WAV with declared
  full-scale pascals.
- The euler-disc exact-ZOH runtime stays cinematic;
  `fs_couple::modal_acoustic_time` is the music owner
  (§10).

---

## 10. Reuse ledger — duplications and their owners

The tree already contains the duplications D23/D24
exist to stop. Recorded here so no image multiplies
them:

| Surface | Situation today | Owner ruling for this plan |
| --- | --- | --- |
| Hertz algebra | Exists 3×: `fs-tribo`, `fs-contact` patch laws, `fs-euler-disc-e2e` | Music images bind to `fs-contact` patch laws; no fourth |
| Hunt–Crossley | 2 shaped laws: `fs-contact` (rate factor on true Hertz) and `fs-dcontact` (χ on the penalty force) — algebraically different, both legitimate | Point-patch debug/bake-off = `fs-contact`; distributed lay/rattle loss = `fs-dcontact` χ. Both stay, named. No third |
| Contact `(K, α)` data | `fs-dcontact` provenance field is a recorded placeholder | Packs go to `fs-matdb` (the wiring point its CONTRACT names) |
| Flat-plate assembly | `fs-plate` (validated DKT) vs `fs-solid::shell` (estimate-only Euler-disc surrogate, self-disclaimed) | Music bodies = `fs-plate`; shells = fs-plate’s recorded flat-facet follow-up. `fs-solid::shell` is never a music owner |
| Exact-ZOH modal runtime | 2×: `fs_couple::modal_acoustic_time` and the euler-disc cinematic runtime | Music owner = fs-couple; convergence only through h7xu5.7.8 |
| WAV encoders | 2×: `fs_couple::pcm_wav`, euler-disc `audio_artifact` | §9.3 — one owner, no third |
| Port containers | `fs-couple` raw `Port` is legacy; `interconnect` refuses it | `PortSchema` v2 only (D20) |
| Steppers | Exactly three, by design: `fs_phs::step` (Gonzalez, L2), fs-couple exact ZOH (L3 audio), `GeneralizedMaxwell::step` (Prony internal). `fs-time` (gen-α, IMEX, exp) is deliberately unused here — L2 `fs-phs` cannot depend on L3 `fs-time`, and the audio steppers’ exactness classes are the point | A fourth stepper needs a displacement argument in this plan |
| `TravelingWaveLine` | Dead: zero in-tree consumers; fs-couple’s CONTRACT negates it twice | No music image may cite it; the char image is the exact-FIR port |
| `fs-acoustics` | The new-domains §5.11 crate was never built; the lane landed as fs-duct / fs-phs / fs-bem / fs-couple / fs-vfit / fs-psycho instead | Name retired. Bead 0ja4 carries the repoint comment (§11); residual §5.11 rungs extend existing crates |

---

## 11. Governance and hygiene

- **New-domains reconciliation.** §5.11’s `fs-acoustics`
  is retired as a crate name (§10). A status note now
  sits under that section; bead
  `frankensim-ext-acoustics-core-0ja4` carries the
  repoint comment listing which rungs are satisfied by
  which existing crates and which remain (interior
  volume acoustics, both-domain power closure, convected
  Helmholtz/LEE, NAFEMS R0083 / NASA CAA decks).
- **Constitution.** The flagship plan gains one
  companion-pointer paragraph; sound is otherwise
  invisible in it, which mislabels a 12-crate lane as
  unchartered work.
- **Beads.** The musical-acoustics epic (tx6k8) closed
  2026-08-13 while physics kept landing as comments on
  the closed epic — new work under this plan gets new
  beads. Filed 2026-08-14 under labels `music-plan`,
  `musical-acoustics-gap`:
  - `frankensim-music-t-brass-zl-mm-zolja`
  - `frankensim-music-t-piano-felt-87zbd`
  - `frankensim-music-t-out-render-ib15w`
  - `frankensim-music-claims-registry-mc31g`
- **Hygiene debts (recorded, not expanded here).**
  `fs-duct` and `fs-plate` have no `tests/` directories
  (inline-only batteries); the acoustics crates are not
  enrolled in the strict-Clippy rollup; cross-ISA
  determinism goldens are absent across the stack (a
  G5 audit away, and unclaimed until run).

---

## 12. Parallel tracks

Tracks own their menus. **No shared M1.**

| Track | First slice (does not block others) | Bead |
| --- | --- | --- |
| T-Brass | `Z_L(ω)` bake driver → `Termination::Tabulated` → MM lines (the trumpet-claim unblock); MM TMM + HB follow | **zolja** |
| T-Piano | Wool-felt `Uniaxial` + felt pack + `FiniteGapPoint` island + 1 string + 1 board; Hunt–Crossley (fs-contact’s) is the debug image | **87zbd** |
| T-Out | Block render API + one pascals→PCM owner + the budget lane | **ib15w** |
| T-Claims | `instrument-claims.json` + xtask gate + corpus rows | **mc31g** |
| T-Wind | Char+holes+valve **and** TMM R gate; VFIT hold as second image (mostly live — needs registry rows) | next |
| T-String | Modal ZOH gate (exists); optional dispersive-WG bake-off | next |
| T-NL | KC on N modes (exists — needs registry rows) | next |
| T-Bow | Friction island (exists — determinism row honest re fs-tribo) | next |
| T-Voice | Two `A(x)` + 1-DOF **and** TMM formants; two-mass as second valve (greenfield product on live ports) | later |
| T-EM | `DescriptorPortHamiltonian` Kirchhoff DAE (the named fs-phs trigger) + Faraday port + device cards (greenfield) | later |
| T-Jet | 3-D sweep → lab card + one pipe (2-D refuses broadband by type) | later |
| T-Shell | Bar modal; fs-plate flat-facet shell later | later |
| T-Sched | Selector module in fs-couple; calls kernels; no physics | later |

Fuse **after that track’s bake-off**, per kernel (§9.2).

---

## 13. Never

Samples; IR cabinets as truth; instrument crates; one
stepper for all fillings; **one image per filling as
policy**; Newton on linear pHS; coarse ladder as HF air
truth; FIR keymaps; formant EQ; forced Helmholtz saw;
in-loop BEM/LBM; deleting a passing orthogonal image
because another is newer; Hunt–Crossley as the
piano-felt law; `ManderConcrete` / `MenegottoPintoSteel`
as a hammer; a parallel “music felt” constitutive next
to `Uniaxial` / `FractionalZener`.

New in v8: resurrecting `TravelingWaveLine`; a second
`Port` container next to `PortSchema` v2; a fourth
stepper without a displacement argument; a parallel
`fs-acoustics` crate; a fourth Hertz or third
Hunt–Crossley implementation; a ◆ live-default without
a green gate and a measured budget row; deleting a
registry image row instead of marking it `refused`;
claiming humidity effects from a card that does not
model them; simulating key/valve mechanisms instead of
chart schedules.

---

## 14. Open questions (local to a menu)

Thiran order (T-Wind char). When to hop char→VFIT.
Wool-felt `Uniaxial` envelope (which coupon); Prony
band vs hammer spectrum. Hunt–Crossley remains
**debug** until a bake-off against the hysteresis
island says otherwise. JCA/Biot linings are a
different claim, not an OQ on the hammer.
Lip pair vs beam (both on brass menu). Licensed `A(x)`
**and licensed impedance/FRF corpora** (fs-modalid found
no license-compatible published FRF trace — the hunt is
a registry row). Jet-card residual. Whether
dispersive-WG string earns a keep-both vs modal.
Moist-air `GasSpec` mixture rule (variable `M`, `γ`).
Scattering-form `|S| ≤ 1` fitting for reflectance ports
(fs-vfit’s recorded follow-up; today’s law is
impedance-form positive realness). MM line-count policy
per flare. Where the felt coupon comes from.

**No track blocked by another track’s OQ.**

---

## 15. Advice condensed

1. **Menus.** Several orthogonal images per filling;
   bake-off; keep both when they serve different claims.
2. **FD for design/lock, TD for attack/play.** Don’t
   make TMM sing and don’t make a delay line design a
   bore.
3. **Long air = delay+junctions *or* VFIT hold *or*
   TMM.** Not a coarse ladder claiming 8 kHz.
4. **Vibration = modal ZOH *or* reduced Gonzalez.**
   Not characteristics.
5. **Contact / valve / device = small islands**,
   possibly several constitutives on the menu.
   **Reuse the constitutive stack we already have**
   (`Uniaxial`, `FractionalZener`/`GeneralizedMaxwell`,
   `fs-contact` patches). Do not mint a parallel
   “music felt” law.
6. **Brass ≠ clarinet with more states.** Separate
   menu (MM + HB + lips).
7. **Voice ≠ woodwind with flesh Q.** Separate menu
   (`A(x)` + glottis); may *reuse* `k_char`, not
   `k_valve` defaults.
8. **Electric ≠ acoustic with a magnet.** Circuit
   split + Faraday — on the descriptor-pHS core that
   already exists.
9. **Selector is the product.** Gesture × claim ×
   budget → image set — reading the registry, not
   prose.
10. **Fuse each kernel after its own gate.** Parallel
    calendars.
11. **Bind before you build** (D23). Every image names
    its owner; a new crate needs a displacement
    argument.
12. **A real-time claim is a measured budget row**
    (§9). Never prose.

---

## 16. Diff from v6 / v7

| v6 | v7 | v8 |
| --- | --- | --- |
| Menu of orthogonal images | Same menus; **ingest law** spelled out | Menus bound to **owners** (D23); live/NEW status corrected against the tree |
| Bake-off; keep-both default | Unchanged | Bake-offs emit receipts referenced by the **claims registry** (§8) |
| D18 is *one* air image | Unchanged | Unchanged |
| One-page geometry stub | Full mesh → cards → playing loop (trumpet §7.2, piano §7.5) | Extractor owner (`fs-query`) + honest authority colors; humidity / JCA / matdb facts corrected |
| Felt: two P-Contact images, equal | Felt **is** `Uniaxial` + `FractionalZener`/Prony; Hunt–Crossley is debug | Unchanged, plus: the HC debug image is `fs-contact`’s existing law (no third implementation) |
| — | — | NEW: budget law + product surface (§9), reuse ledger (§10), governance (§11), beads filed (§12), **frozen** (§17) |

v4/v5 physics (d’Alembert, don’t Newton linear pHS,
Schur island, R(ω) for reeds) remain **entries on
menus**, not the constitution of every filling.

---

## 17. Planning-workflow next — and the freeze

The v7 request (“beadify menus”) is **done**: zolja
(T-Brass), 87zbd (T-Piano), ib15w (T-Out), mc31g
(T-Claims), labels `music-plan` /
`musical-acoustics-gap`. Do not bead “the” instrument
stepper.

**Freeze rule.** This plan went v3→v8 in one day. The
next unit of value is a bead executed against it, not
v9. A future revision must cite landed commits (which
bead, which gate went green, which registry rows
changed) — doc-only churn is refused.

---

*End of v8. Constitution:
`COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md`. Shared ports.
Orthogonal interiors. **Menus, not winners.**
3-D + cards mint charts; charts play — through the
crates we already have.*
