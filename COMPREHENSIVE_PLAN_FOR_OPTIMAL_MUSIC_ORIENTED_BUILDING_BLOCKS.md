# Comprehensive Plan for Optimal Music-Oriented Building Blocks

> Nameless physics. Emergent sound. The *same PDE*, split the
> way it wants to be split: **characteristics for the wave,
> Cauer/Foster for electrically short bits and relaxation,
> a tiny implicit island for the valve/felt/bow.**
>
> TMM is the frequency oracle of that PDE, not a source of
> fitted delay lines. A 6-cell LC ladder is not “the”
> physics of a 34 cm bore — it is already a reduction, and
> a dispersive one.
>
> Companion to `COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md` and
> the new-domains plan §5.11. Does not amend the Decalogue.
>
> Status: **v4**. v1 FIR bake. v2 playable scatter + loop.
> v3 don’t Newton a linear pHS. v4 don’t replace d’Alembert
> with a coarse ladder and call that accuracy.
>
> Tags: `[S]` `[F]` `[M]`. Date: 2026-08-14.

---

## 0. How to read this

The argument is §1.4–§1.5 and §6. Doctrine §2. Order of
work §12. Family names are fillings, never crates.

---

## 1. Intent

### 1.1 Hear physics, play physics

Materials + geometry + gesture + air. Reed, lips, jet,
string, felt, shell, vocal `A(x)`, pickup/circuit/cone.
No samples. No instrument types.

### 1.2 Success

Honesty, **same PDE**, playability (time-varying geometry
without rebake-per-key), perceptual adequacy, named 48 kHz
budget on a *phrase*. Fail the budget → renderer. Fail the
PDE → synth. Fail playability → sampler of our own notes.

### 1.3 Workflows

Author → authority (TMM + well-resolved TD) → **split
realization** (not a second theory) → play → ear-weighted
gate on **R(ω) at the valve** and radiated `p` → escalate
→ optional adjoint calibration of *cards*.

### 1.4 Three category errors (v1–v3)

**E1. Frozen `Z_in` FIR as the instrument.** (v1) A
fingering phrase is not `2^N` impulse responses.

**E2. Newton on a linear bore.** (v3) Quadratic `H`,
state-independent `(J,R)` ⇒ linear midpoint or ZOH, one
factor per control tick.

**E3. “The LC graph we already wrote is the physics.”**
(v3’s remaining miss.) `acoustic_chain` with ~6 cells on
a speaking length `L` has `Δx ~ L/6` and a spatial
Nyquist `c/(2Δx) ~ 3c/L` — about 3 kHz on a 34 cm air
column. Brightness, reed locking on upper harmonics, and
brass flare live *above* that. TMM does not have that
dispersion. A bilinear WDF of those 6 cells **inherits
it**. Making Newton cheap on a wrong spatial operator
does not make a better clarinet.

The one-pole `TravelingWaveLine` was wrong about *loss*.
It was not wrong that **1-D lossless acoustics is a
delay**. That is d’Alembert, not a hack.

### 1.5 The statement that is both more accurate and cheaper

The linear duct is a **singularly perturbed hyperbolic
system**: a principal wave part plus a stiff, spatially
local viscothermal relaxation (ZK / Foster).

```text
∂t u + A(x) ∂x u  =  B(x) u  +  (memory / Foster)
     \_________/     \___________________________/
      exact shift         local ODE / Cauer
      (characteristics)   (electrically short)
```

**Optimal discretization** (operator split, `[S]`):

| Piece | Electrical length | Realization | Cost |
| --- | --- | --- | --- |
| Lossless principal part | `kℓ ≳ 0.2` | Characteristics: delay of `ψ` or `(p±Zc U)` | O(1) per section, **no HF ladder dispersion** |
| Area jump, sphere factor `1/x` | interface | Scattering / transformer (already `x` on ψ) | O(1) |
| Hole, pad, mouth, wall, ZK | `kℓ ≪ 1` or local in `ω` | Cauer/Foster / T-junction, bilinear or ZOH | O(n_br) |
| Short chamber, mouthpiece, compact neck | `kℓ ≪ 1` | Lumped LC (already the chimney rule) | O(1) |
| Valve, felt, bow, jet, KC `H` | — | Small implicit island | O(1)–O(10) |

TMM is the **symbol** of the left-hand side in `ω`. A
split stepper is correct when its driving-point
**R(ω)** (what the reed sees) and radiated **p/U** match
TMM inside an ear-weighted band.

This is Smith/Bilbao *as a splitting of our PDE*, not as
a pin-frequency fit, and not as “throw away the delay
and keep 6 inductors.”

We already half-do this: compact chimneys lump; long
necks become a line; cones use `ψ = xp` (the spherical
characteristic). v4 **makes that the global rule.**

---

## 2. Doctrine

D1–D13 as in v3 (nameless, no samples, two lanes / three
rates, reduce don’t invent, correctness first, Franken-
only, certificates, named RT budget, psycho judges, flags,
passive time-varying ports, no synth crate, calibrate
cards).

**D14.** Do not Newton a linear pHS.

**D15.** Child charts are coarsenings with a QoI receipt —
but **do not coarsen a long hyperbolic section into few
LC cells** to “save states.” That spends states on
numerical dispersion. Spend them on Foster/junctions.

**D16.** Multi-rate by stiffness (felt µs, audio, control).

**D17.** Lift state when `σ`, delay length, or reduced
basis changes.

**D18. Split by electrical length.** `[S]`
`k_max ℓ < ε` (ε ~ 0.2, same spirit as `tap_is_line`) ⇒
lump. Else ⇒ characteristics (+ local relaxation).
Cone long-wave ⇒ characteristics on `ψ`, not a stack of
cylinders unless multimodal.

**D19. The valve-loop QoI is R(ω), not only Z_in.** `[S]`
Playability is the reflection the island sees. Radiated
timbre is a second QoI (`p_rad / U_mouth` or baffled `p`).
Gates and DWR mark cells/branches against **those**, not
L2(`p`) in the bore.

**D20. Strang (or WDF) pieces must each be passive; the
composition needs a receipt.** `[S]`
Midpoint of a passive linear pHS is passive. WDF of a
passive adaptor is passive. Naive split of shift + loss
can drift. G10 below.

---

## 3. Relation to existing plans

Constitution wins. `fs-acoustics` / `fs-bem` supply
exterior tables and shells. `fs-phs` Galerkin stays the
tool for **modal** islands (strings, plates, frozen
1-ports), not the default for a long played bore.
Mega-fused work waits until the split 1-port×island
exists.

---

## 4. Current state, reread again

TMM + ODE + reed + string + plate are real. Ernoult
frozen peaks are real. Helmholtz `k` is real.

What we do not have: a **characteristic image** of the
duct that TMM would recognize above the ladder Nyquist;
a **linear** stepper for the LTI piece; a slur that is
not Newton-on-80-states; R(ω)-weighted gates; lip/jet/
felt/pickup physics.

`TravelingWaveLine` stays dead (one-pole loss). Its delay
is rehabilitated as d’Alembert of the principal part.

---

## 5. Atlas

Media, meshes, valves, strings, plates, felt, pickup,
circuit, speaker, stochastic sources: as before.

### 5.1 Linear acoustics — three equivalent writings

1. **TMM** — exact (to the section model) in `ω`. Oracle
   for R, Z, peaks. Bessel ZK, cones, holes, walls.
2. **Well-resolved pHS** — Cauer of the same telegraph
   system with `Δx` set by `k_max`, *or* characteristics
   + Foster. Time-domain authority for transients and
   nonlinear coupling.
3. **Split performance image** — D18. Same operators,
   cheap.

6-cell `acoustic_chain` is a **debug / low-band** child,
not the gold standard. If we keep it, say so: “valid
below `f_Ny`.” Do not gate brightness against it.

### 5.2 Characteristic sections `[S]`

Wave variables (cylinder): `w^± = (p ± Z0 U)/2` with
`Z0 = ρc/S` of the **inviscid** principal part. Advance:
shift by `Δt` corresponding to `ℓ/c`. Fractional delay:
Thiran/allpass, error in the ear band gated.

Sphere: same on `ψ = x p`, `U` reconstructed as we
already do (near-field shunt is a **lump** on `ψ`, D18).

Viscothermal: do **not** put `c_eff(ω)` into the delay
every sample. Keep delay at `c`. Put ZK in a **local
passive operator** on `w^±` (Foster of the boundary
layer / `F(r_v)` remainder). TMM checks the composite.

### 5.3 Junctions `[S]`

Area jump: lossless scatter from `S_L, S_R` (exact).
Hole: existing T-junction + `Y(σ,ω)` → low-order Foster
in `σ` (evaluate at control rate). Mouth: mass + Foster
`R(ω)` or table `Z_L`. Wall: existing `WallPin` as shunt
Foster. Mutual hole series: extra `t_s` as now.

Time-varying `σ`: interpolate **Foster parameters** under
D11, lift junction state (D17).

### 5.4 When modal / VFIT 1-port wins

Frozen (or slowly varying) linear map as seen by one
port: **IRKA / `fs-vfit` / Galerkin onto peaks** of TMM
`Z_in` or `R`. Fewest states per ear residual. Attack
or slur uses the spatial split; after `σ` settles, **drop
to the 1-port** with a lift (two-phase note). Sustain of
a clarinet can be ~10–20 passive states + delay memory,
not 80 LC cells.

### 5.5 Nonlinear island + Schur 1-port

```text
(w^-, x_F)  --linear split-->  p_face
u_face = F(p_up − p_face, y)     # Bernoulli / jet
y' = f(y, p_face)                # MSD / two-mass / felt
```

Implicit in `F,y` only. Linear piece: shift + ZOH Foster
+ scatter. Gate vs Gonzalez on a **short** fixture and
vs TMM `R` on the long fixture.

### 5.6 Goal-oriented coarsening

Adjoint of **R(ω)** at the reed and of **p_rad**. Delete
Foster branches / extra slices that do not move those
QoIs by more than the ear budget. Never delete the
optical length (the delay).

### 5.7 Describing function / HB `[F]`

For lock: valve describing function + `Z(nω)` from TMM
(Fletcher). Brass/voice “will this speak?” without a
transient. Attacks still TD.

### 5.8 Multi-rate

Felt substep. Board radiation may run at block rate if
G9 says so. Control: `σ`, `A(x)`, `T`, slide.

---

## 6. Architecture

```
 geometry + gas + gesture
           │
           ▼
    TMM oracle  (R, Z, peaks, p_rad/U)
           │
    split realization (D18)
      delays ‖ lumps ‖ Foster
           │
    control tick: update σ, A, c(T); lift; D11
           │
    audio: shift + scatter + ZOH Foster
         + 1–3 dim island
           │
    optional: settle → VFIT 1-port sustain
           │
        observer (piston / table / ISO)
```

### 6.1 Why this is more accurate *and* faster

| Approach | HF accuracy | Playable `σ(t)` | Cost / sample |
| --- | --- | --- | --- |
| 6-cell Gonzalez | Poor (dispersion) | Yes, expensive | Newton ~ n²–n³ |
| 6-cell linear midpoint (v3) | Still poor | Yes | Factor at control, O(n) audio |
| Fine LC (80 cells) midpoint | Good | Yes | Heavy memory + factor |
| Pin-fit scatter (v2) | Accidental | Yes | Cheap, wrong Zc(ω) |
| **Split D18** | Good if TMM-gated | Yes | **O(sections + n_br)** , mostly shift |
| Frozen VFIT 1-port | Good in-band | No | O(n_red) |

### 6.2 Control tick

Rebuild junction Foster from `σ`. Warp delays by `c(T)`.
Lift `w^±` and Foster states. D11 ramp test in CI, not
every note.

### 6.3 Audio tick

The mega-fused target is **shift + 2×2 scatter + island
solve**. Not a sparse 80×80.

### 6.4 Quality knob

`{TMM-offline, split-fine, split-coarse, VFIT-sustain}`.
Same PDE, different Δx / n_br / 1-port switch.

### 6.5 Inverse calibration

Fit `A(x)`, wall `r`, lip rest, felt — so **TMM + island**
match a cited fixture. Do not ship the fixture WAV.

---

## 7. Geometry

Cited mesh/table → `fs-io` → centerline/`A(x)`/shell →
**section list tagged lump vs characteristic** (D18) →
TMM + split image. MRI vowels are already `A(x)`.

---

## 8. Fillings

- **Reed winds:** island + characteristic bore + lumped
  holes. MVP: Ernoult *slur* with R(ω) vs TMM, not a
  brighter 6-cell Newton.
- **Brass:** lip island + multimodal characteristics
  (higher modes are extra delay lines above cutoff) +
  `Z_L` table + HB/describing-function lock.
- **Jet:** lab card or refuse.
- **Strings:** modal/KC as now; linear strings are
  already characteristics in mode space (ZOH).
- **Piano:** D16 felt + modal strings/board.
- **Voice:** glottal island + `A(x)` characteristics
  (classic Kelly–Lochbaum **as D18 of the tract PDE**).
  Ooh/aah = two cards. Interaction = R(ω) of the tract.
- **Electric:** Faraday + circuit island (split RLC vs
  device) + speaker 1-port.

---

## 9. Validation

| Gate | Question |
| --- | --- |
| G0 | pHS / passivity / D11 |
| G2 | Ernoult, published Z, Fant `A(x)` |
| G6 | Ear-weighted **R** and **p_rad** vs TMM |
| G7 | Slur / vibrato: no click, supply-rate |
| G8 | Island+split vs Gonzalez on a *short* nonlinear fixture |
| G9 | Child vs fine split (n_br, Δx) on R and p_rad |
| **G10** | Split vs TMM: max ear-weighted \|R−R_TMM\|; composition passivity |
| **G11** | Spatial Nyquist disclosed; no brightness claim above it |

Reference audio: test-only.

---

## 10. Never

Samples, IR cabinets as truth, instrument crates, in-loop
CFD, fake mean flow, one-pole *loss*, Newton on linear
pHS, **claiming a coarse ladder matches TMM in the
partials**, `2^N` FIR keymaps, pin-frequency delay fits
that ignore `Zc(ω)`.

---

## 11. Open questions

| ID | Question | Blocks |
| --- | --- | --- |
| OQ-1 | ZK on characteristics: per-junction Foster vs distributed splitting along the delay? Start junction-local (matches TMM stations). | M1 |
| OQ-2 | Fractional delay: Thiran order vs ear-band group-delay error. | M1 |
| OQ-3 | When to hop spatial split → VFIT sustain (dwell time, residual). | M1 |
| OQ-4 | Licensed `A(x)`. | M2 |
| OQ-5 | Lip pair vs beam. | M5 |
| OQ-6 | Felt paper. | M4 |
| OQ-7 | How many multimodal lines for a flare before the ear saturates. | M5 |

---

## 12. Phases

**M0** — D14 `step_linear` for true LTI pHS (plates
lumps, short necks). D18 tagging on sections. G11
Nyquist printed on every `acoustic_chain` used as a
claim. TMM remains the linear oracle.

**M1 — Load-bearing** `[S]`
Characteristic image of a **cylinder + holes + reed**:

1. Delays + lumped T-junctions + mouth Foster (physics
   we already have, new assembly).
2. ZK as junction Foster, G10 vs TMM `R`.
3. Valve island, G8 vs Gonzalez on a 1-section fixture.
4. `σ(t)` + lift + G7 slur on Ernoult.
5. Optional VFIT hop on a held fingering (OQ-3).

**Do not** implement WDF-of-6-cells as M1. That repeats
E3.

**M2** — `A(x)` characteristics (vowels, cones as `ψ`
delays).

**M3** — Pickup + circuit split + speaker.

**M4** — Felt D16 + board modal.

**M5** — Multimodal delays + lips + `Z_L` + describing
function / HB.

**M6** — Jet card or refuse.

**M7** — Shell modal.

**M8** — Fuse shift×scatter×island; SIMD on VFIT/modal;
scheduler. After G7+G10.

**M9** — Gesture skin.

```
D18 cylinder image → G10 vs TMM → island → G7 slur = M1
        ↘ A(x)/ψ → M2
        ↘ circuit → M3
        ↘ felt → M4
multimodal + HB → M5
jet → M6
shell → M7
G7∧G10 → M8
```

---

## 13. M1 work packages

- Section tagger: lump vs characteristic (shared ε with
  `tap_is_line`).
- Characteristic stepper: `w^±` shift, fractional delay,
  area-jump scatter.
- Hole/mouth/ZK as existing Foster/T-junction **on the
  wave variables**.
- G10: R(ω) vs `input_impedance` / `input_impedance_wall`.
- Valve island + G8.
- Lift + D11 + G7 on Ernoult `xxxx→xxxo`.
- Ear-weighted R residual helper.
- Disclose `f_Ny` on ladder children (G11).

---

## 14. Advice condensed

1. **A long bore is a delay plus local operators**, not
   a handful of inductors. That is more accurate at 3–8 kHz
   and mostly a memmove.
2. **TMM judges R(ω) and p_rad.** The time stepper is
   wrong until those match, not until energy on a 6-cell
   toy matches Gonzalez.
3. **Newton only the island.** Linear midpoint for lumps
   that really are lumps.
4. **Playable holes are local Y(σ), not a new FIR.**
5. **Sustain may collapse to a VFIT 1-port** after the
   slur; lift into it.
6. **Coarsen Foster, never optical length.**
7. **Voice is Kelly–Lochbaum as D18 of A(x)**, not a
   formant filter, not 6 tract cells.
8. **Brass is extra modal delays + lips + HB**, not a
   brighter cone ladder.
9. **Fuse the shift×island product last.**
10. **Calibrate cards; never ship the take.**

---

## 15. Diff from v3

| v3 | v4 |
| --- | --- |
| Same LC graph, cheaper integrator | Same *PDE*; LC only if electrically short |
| WDF = bilinear of 6 cells | That WDF still lies above ~c/(2Δx) |
| Scatter identified or inherited | Scatter = characteristics + TMM junctions |
| QoI ≈ ear-weighted Z_in | **R at the valve** + radiated p |
| M1 = step_linear on acoustic_chain | M1 = characteristic cylinder + TMM G10 + slur |
| Delay lines as a possible image | Delay is the *principal-part solution* |

v3’s D14–D17 stay. They apply to lumps, plates, Foster,
and islands. They do not justify a coarse ladder as a
trumpet.

---

## 16. Planning-workflow next

Review should try to break D18 (when does a ladder beat
characteristics? very short, very nonlinear, 3-D).
Beadify **only** M0–M1: tagger, `w^±` stepper, G10, island,
Ernoult slur. If G10 cannot meet TMM R in the ear band,
the split is wrong — fix the ZK junction, do not add
cells.

---

*End of v4. Constitution:
`COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md`. The cheap, correct
1-D wave is d’Alembert; Cauer is for what is not a wave.*
