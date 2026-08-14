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
> Status: **v7**. v1–v4 each crowned one discretization.
> v5 gave each filling one stack. That was still one-size
> *inside* a family. v6 makes the **menu** the design.
> v7 is how a real 3-D mesh and named material cards
> become those menus (trumpet playing loop, every filling,
> piano felt on the *existing* hysteresis / visco stack).
>
> Tags: `[S]` `[F]` `[M]`. Date: 2026-08-14.

---

## 0. How to read this

§1.4 is the rule. §4 is exactness classes. §5 is the menus.
§6 is the selector. **§7 is the ingest law**: a licensed
mesh plus `fs-matdb` / `GasState` / `Uniaxial` / visco
cards become charts, then authority images, then
performance images. §7.2 is the trumpet (lips, brass,
valves, emergent lock). §7.5 is the piano (felt on the
hysteresis stack we already have). §7.8 is live vs not.
§8 is parallel tracks.

---

## 1. Intent

### 1.1 Hear physics, play physics

Geometry + materials + gesture + air. No samples. No
instrument crates. Real-time on a named core *when the
selected image’s budget says so* — not when a mythical
universal stepper does.

### 1.2 Success is a matrix, not a scalar

A filling **F** with image **I** succeeds on claim **C**
when the gate for `(F, I, C)` is green. Example:

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

v5 listed one “RT” column. v6 lists a **menu** and a
**selector**, not a winner.

### 1.4 Design rule

```text
cards  →  several images {I_k}, each with:
            operator limit, exactness class,
            QoI set, budget, failure mode
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

**D19.** Claims are `(filling, image, QoI)`.

**D20.** Interiors couple only at ports.

**D21. Menus, not winners.** `[S]`
A filling’s CONTRACT lists every legal image. Deleting
an image because another “won” a bake-off is forbidden
unless its failure mode is empty *and* its budget is
strictly worse on every claim. Usually both stay.

**D22. Complementary duals.** `[S]`
Frequency-domain oracles (TMM, HB, eigensolve) do not
have to time-step. Time-domain images do not have to
replace them. Ask lock questions in `ω`; ask attacks
in `t`.

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

## 4. Shared kit (still not a mega-stepper)

Ports, cards, extractors, TMM (air oracle), eigensolvers
(vibration oracle), `fs-vfit`/IRKA, `fs-psycho`,
fingerprint ledger, gesture signals.

Kernels are **separate**: `k_tmm`, `k_char`, `k_lc`,
`k_modal`, `k_gonzalez_red`, `k_vfit`, `k_valve`,
`k_felt`, `k_fric`, `k_jet`, `k_ckt`, `k_hb`. No
`k_everything`.

---

## 5. Menus per filling

Each filling: **cards**, **images** (orthogonal),
**selector hints**, **failure modes**. Images marked
◆ are the usual live default, not the only legal one.

### 5.1 Single-reed / double-reed winds

**Cards:** `A(x)` or cylinder list, hole table, cane/lay,
gas, mouth flange.

| Image | Class | Good for | Failure mode |
| --- | --- | --- | --- |
| **TMM** | X-Consist / X-Struct (ZK) | Peaks, design, Ernoult, R(ω) oracle | No attack, no live `σ(t)` |
| **Char + lumped holes + ZK Foster** ◆ | X-Exact delay + X-Struct loss/holes | Played phrases, slurs | Bad if junctions don’t match TMM R |
| **Fine LC/pHS + `step_linear`** | X-Consist | Short chambers; TD oracle vs Gonzalez | HF dispersion if `Δx` large (disclose `f_Ny`) |
| **VFIT/IRKA 1-port of TMM R** | X-Struct | Held fingering, cheap sustain | Live `σ`; out-of-band |
| **HB + TMM Z(nω)** | X-Struct | “Will this reed speak?” | Transients, multiphonics |

**Selector:** `dσ/dt > ε` → char; settled → may hop VFIT
(lift). Design CLI → TMM. Studio bounce may use fine
pHS if `f_Ny` covers the claim.

**Do not:** FIR keymap; Newton the bore; treat 6-cell
`f_Ny` as 8 kHz truth.

### 5.2 Brass

**Cards:** flare `A(x)`, mouthpiece, lip pair, `Z_L` table.

| Image | Class | Good for | Failure |
| --- | --- | --- | --- |
| **Multimodal TMM** | X-Consist | Impedance, cutoff, design | No lips in time |
| **HB / describing function** ◆ lock | X-Struct | Slot, pitch lock | Attack, lipping |
| **MM characteristic lines + lip island** ◆ play | X-Exact per mode + island | Attack, slurs, slide | If plane-wave only, dull / won’t lock |
| **BEM-baked `Z_L(ω)`** | X-Struct / X-Est | Unflanged mouth | In-loop BEM |
| **Plane-wave char (5.1)** | — | Debug / mute | **Not** a trumpet claim |

Plane-wave and multimodal **both stay on the menu**.
The trumpet *claim* requires MM. A brass mute study
might use plane-wave + table.

### 5.3 Jet-pipes (flute, recorder, flue, organ rank)

**Cards:** pipe `A(x)`, chimney, jet-card from lab.

| Image | Class | Good for | Failure |
| --- | --- | --- | --- |
| **LBM/aeroac lab** | X-Est→X-Struct when residual issued | Mint jet card | Not RT |
| **Jet island + char pipe** ◆ | X-Struct + X-Exact delay | Flute / recorder phrase | No card → refuse |
| **TMM pipe only** | X-Consist | Cut/tone-hole design | No source |
| **N independent 5.3 voices** | — | Organ rank | One fat 3-D flue |

Organ and flute share kernels, **not** a shared “wind
stepper” with reeds. Reed images are off this menu.

### 5.4 Linear plucked / struck strings + body

**Cards:** μ, T, EI, plate DKT/orthotropy, cavity.

| Image | Class | Good for | Failure |
| --- | --- | --- | --- |
| **Modal ZOH string** ◆ | X-Exact (truncated basis) | Pluck, harp, ordinary guitar | Missing high modes / tension change |
| **Eigensolve plate + modal ZOH** ◆ | X-Exact truncated | Body radiation | Attack thump if N too small |
| **VFIT of plate driving-point** | X-Struct | Cheap body as 1-port on the bridge | Nonlinear top |
| **Fine pHS string (many modes)** | X-Consist | Authority vs ZOH truncation | RT polyphony |
| **Dispersive waveguide** `[F]` | X-Struct | Alternative string (stiffness on the delay) | Bake-off vs modal; don’t delete modal |

Modal and (optional) dispersive waveguide are
**orthogonal string images**. Bake-off on inharmonicity
and a pluck spectrogram; **keep both** if both pass and
budgets differ (polyphony vs one hero string).

**Off menu:** air-column characteristics, 48 kHz DKT.

### 5.5 Nonlinear string / plate

| Image | Class | Good for | Failure |
| --- | --- | --- | --- |
| **KC / von Karman Gonzalez on N modes** ◆ | X-Consist on the reduced pHS | Loud guitar, gong-ish | N too small; mesh-scale folds |
| **Linear modal (5.4)** | X-Exact truncated | Soft dynamics | Large-amplitude pitch glide |
| **Full-mesh NL** `[M]` | X-Consist | Offline crash | RT |

Selector: `|Δℓ/ℓ|` or energy in cubic terms vs a
threshold. Soft notes may stay on 5.4 in the same
voice (two images, one filling).

### 5.6 Bowed

| Image | Class | Good for | Failure |
| --- | --- | --- | --- |
| **Modal + Stribeck island** ◆ | X-Exact modes + X-Struct friction | Violin-class mechanism | Not a measured rosin unless card says so |
| **KC + friction** | P-NLH + island | Loud / high bow | Cost |
| **Forced Helmholtz oscillator** | — | **Forbidden** (not physics) | Always |

### 5.7 Piano-class

**Cards:** string `(μ, T, EI)` from steel + tension as
state; soundboard / rim from spruce `OrthotropicElastic`
(already in `fs-material`); hammer **felt** from the
existing hysteresis stack, not a new constitutive
family — see §7.5.

| Image | Class | Good for | Failure |
| --- | --- | --- | --- |
| **Modal strings + Uniaxial/visco felt island + modal board** ◆ | X-Exact + X-Struct felt | One note / polyphony | Felt card wrong or out of visco band |
| **3-string unison coupling** | X-Consist | Chorus, aftersound | Extra states |
| **Linear hammer spring** | X-Exact | Debug | No hammer spectral tilt |
| **Hunt–Crossley hammer** | X-Struct | Debug / bake-off vs hysteresis | Residual-strain / tilt of real felt |
| **Duplex as extra modal segments** | X-Exact | Aftersound | — |
| **JCA/Biot felt-as-absorber** | X-Struct | Case linings, under-string cloth | **Not** the hammer contact |
| **Air-column / char tube** | — | **Off menu** | Always for this filling |

The live default felt image **consumes physics already
in the tree**:

1. `fs-material::Uniaxial` — the hysteretic fiber
   contract (committed reversal state, consistent
   tangent, cycle dissipation > 0). Piano felt is a
   **new implementor** of that trait (wool-felt
   envelope + residual crush), not Mander concrete
   and not Menegotto–Pinto steel wearing a felt hat.
2. `fs-material::visco` — `FractionalZener` is the
   fitting-side law the crate already names for
   **wood, polymers, and felt** (near-constant loss
   factor). `lower_to_prony` emits a
   `GeneralizedMaxwell` with a certified band; the
   island steps that Prony series. Out of band
   refuses (`FS-MAT-VISCO-OUT-OF-BAND`).
3. `fs-contact` — the hammer–string *geometry* of the
   contact (FiniteGapPoint from the felt thickness
   field, or a Hertz patch). Hunt–Crossley is the
   **debug / second image**, not the felt law.
4. Taxonomy row (JCA / Biot: flow resistivity,
   porosity, tortuosity, characteristic lengths) —
   felt as a **porous absorber** (linings). A
   different claim from hammer compression.

Hunt–Crossley stays on the menu for bake-off against
a hammer-velocity coupon. It does **not** replace
the hysteresis stack. A linear spring is debug only.

### 5.8 Bars / plates / bells

| Image | Class | Good for | Failure |
| --- | --- | --- | --- |
| **1-D beam modal** ◆ | X-Exact truncated | Xylophone, bar | 3-D bell |
| **DKT plate modal** ◆ | X-Exact truncated | Gongs, soundboards | Deep shell |
| **Axisymmetric shell** `[F]` | X-Consist | Bells | Cost, mesh |
| **Strike island** | X-Struct | Attack | — |
| **Char air** | — | Off menu | — |

Beam vs plate vs shell are **different charts**, not
refinements of one stepper.

### 5.9 Voice (ooh / aah, then more)

| Image | Class | Good for | Failure |
| --- | --- | --- | --- |
| **TMM of `A(x)`** | X-Consist | Formant oracle | No glottis in time |
| **KL / char tract + glottal island** ◆ | X-Exact delay + island | Sung vowels, articulation | 3-D mouth, fricatives |
| **VFIT 1-port of tract R** | X-Struct | Held vowel cheap | Moving tongue |
| **Two-mass glottis** vs **1-DOF** | two P-Valve images | Bake-off on spectra | — |
| **WallPin tissue** | X-Struct | Bandwidth / extra formants | — |
| **Jet card in tract** | 5.3 add-on | Fricatives later | Don’t block /u//a/ |
| **Formant filter / vocoder / WAV** | — | **Forbidden** | Always |

1-DOF and two-mass **both stay**. Articulated KL and
held VFIT **both stay**.

### 5.10 Electric (pickup → circuit → speaker)

| Image | Class | Good for | Failure |
| --- | --- | --- | --- |
| **Lumped Faraday** ◆ | X-Struct | Volt from `v_string` | Magnet geometry |
| **Magnetostatic bake** `[F]` | X-Struct | Position-dependent B | Cost |
| **Split circuit (linear factor + device)** ◆ | X-Consist linear / X-Struct device | Distortion as constitutive | Stiff device without island |
| **Thiele–Small / plate cabinet** | X-Struct | Speaker + box | Room |
| **Cabinet IR pack** | — | **Forbidden** | Always |

String images are **whatever 5.4/5.5 selected**. The
electric path does not re-discretize the string.

### 5.11 Cross-filling graphs

String×plate×cavity: pick 5.4/5.5 string image × 5.4/5.8
body image. Reed×plate: 5.1 char/VFIT × 5.4 modal.
Pickup sits on a 5.4/5.5 port. **No fused interior.**

---

## 6. Selector and bake-off

### 6.1 Selector (runtime)

```text
select(F, gesture, claim, budget) → {I_k}
```

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
(TD vs FD, phrase vs hold, 1-DOF vs two-mass).

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

Materials are **cards on labeled regions**, not a
vibe. They come from `fs-matdb` receipts into
`fs-material` laws that already exist:

| Region | Card (examples) | What it is *not* |
| --- | --- | --- |
| Air inside / outside | `GasState` from `GasSpec` + T, p, humidity (`fs-material::gas`, already live) | “trumpet air” |
| Yellow-brass body | isotropic `E, ν, ρ` via `resolve_isotropic_elastic_state_point` (70/30 CuZn from matdb / a cited coupon) | a brightness knob |
| Lip / fold tissue | tissue `ρ, E` (or Lamé), loss factor | cane reed, or a mass typed by ear |
| Cane reed | orthotropic cane | a generic “reed” |
| Spruce top / soundboard | `OrthotropicElastic` (already in `fs-material`) | a body IR |
| Piano-hammer felt | **new `Uniaxial` implementor** + `FractionalZener` → `GeneralizedMaxwell` (the visco crate already names felt) | Hunt–Crossley as the felt; Mander concrete in a felt hat |
| Felt lining / cloth | JCA / Biot porous-absorber row (taxonomy §6) | the hammer contact |
| Steel string | `ρ, E`, and **tension as state** | a sampled string |

Gesture is **geometry in time**, not a MIDI note
number: valve down → extra tube length; key down →
hole `σ`; lip pressure → orifice and pre-stress; bow
`v` and `F_n`; hammer velocity into the felt card.

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
   ratio of the actual alloy. Used for:
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
     lip, not cane, not rubber).

   Those two go into an **offline authority lab**
   (FEM / beam of the mesh + tissue law). That lab
   *mints a valve card*: reduced masses, stiffnesses,
   damping, and an orifice law whose rest shape came
   from the mesh. Same moral status as the LBM jet:
   it is a 3-D experiment that produces a card. The
   triangle mesh does **not** run at 48 kHz.

4. **Air.** `GasState`. Hot stage vs cold rehearsal is
   `c(T)`, which warps every delay and every cutoff.
   Lungs are an upstream compliance plus a
   blowing-pressure gesture, not a “note-on.”

5. **Player gesture.** Three bits (1st / 2nd / 3rd
   valve) plus continuous lip pre-stress / aperture
   plus blow pressure. Optional tuning slide as a
   continuous length.

**What the extractor does (not a `Trumpet` type)**

| 3-D thing | Chart it becomes |
| --- | --- |
| Interior lumen | Centerline + inscribed radius → `A(x)` samples |
| Mouthpiece cup | Electrically short → **lumped** sections / Helmholtz-ish cavity (D18) |
| Cylindrical yards | **Characteristic** delays of length `ℓ/c(T)` |
| Flare | `A(x)` + **local cutoff** of higher modes `f_c(s) ∼ c j_{mn} / (2π r(s))` → P-MMWave lines that exist only where they propagate |
| Valve *up* | Bypass length `L_i` |
| Valve *down* | Bypass + **crook** length `L_i + ΔL_i` taken from the CAD of that slide |
| Bell exterior | Offline BEM or piston → `Z_L(ω)` table (P-Rad bake) |
| Brass wall | Optional `WallPin` from brass `σ, K, r` and local `a` |
| Lip mesh + tissue | Offline reduce → `(m1,m2,k,c, orifice)` valve card |

Certificates we actually check: reconstructed air
volume vs the mesh; `A(x) > 0` everywhere; each valve
`ΔL` matches the CAD within a band; `Z_L` is passive.
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
| Lips | FEM/beam of a lip-shaped mesh, tissue `E,ρ,ζ`; orifice from rest gap | Two-mass / beam + Bernoulli |
| Cane reed | Shell/plate of a reed blank, cane ortho; lay as Hunt–Crossley from lay mesh | 1-DOF or measured MSD + lay |
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
   `fs-material` (nine engineering constants, grain
   axes owned by the geometry consumer). Authority:
   plate eigensolve → modal chart. Performance:
   modal ZOH of those modes (menu 5.7 / 5.4).
2. **Strings as 1-D.** Lengths, speaking points,
   duplex segments from the plate/bridge geometry.
   Card: steel `ρ, E`, and **tension as state** (not
   a frozen MIDI frequency). Image: modal ZOH
   (menu 5.4); 3-string unison as a second image.
3. **Hammer felt that is felt.** A thickness field
   (possibly from a scan or a manufacturer section)
   plus a **wool-felt coupon**. This is where we
   refuse to invent a new constitutive family.

   FrankenSim already captured the physics felt
   needs:

   | Already in the tree | Role for the hammer |
   | --- | --- |
   | `fs-material::Uniaxial` (hysteretic fiber contract: reversal state, consistent tangent, cycle dissipation, FD-gated) | Felt compression is uniaxial with residual crush. New implementor: wool-felt envelope, *not* `ManderConcrete` or `MenegottoPintoSteel`. Same trait, same thermodynamic gates (mt-001 / mt-004). |
   | `fs-material::visco::FractionalZener` | Fitting-side law the crate **already names for felt** (near-constant loss factor over decades). |
   | `lower_to_prony` → `GeneralizedMaxwell` | Runtime island: exact-exponential Prony step, dissipation ledger, out-of-band refusal. |
   | `fs-contact` FiniteGapPoint / Hertz patch | Geometry of the hammer–string interface from the felt thickness field. |
   | `HuntCrossleySphere` / `HuntCrossleyEllipticParaboloid` | **Debug / bake-off image**, not the felt. |
   | Taxonomy JCA / Biot (flow resistivity, porosity, tortuosity, Λ, Λ′) | Felt as **porous absorber** — case linings, under-string cloth. A different claim. |

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
- Treating Hunt–Crossley as the piano-felt law.
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

Already in the tree, and already the right kinds of
objects:

- `GasState` / `GasSpec` in `fs-material`
- `IsotropicElastic` / `OrthotropicElastic` +
  `resolve_*_state_point` from `fs-matdb` receipts
- `Uniaxial` hysteretic contract (`MenegottoPintoSteel`,
  `ManderConcrete`) with FD-gated tangents and
  cycle-dissipation fixtures
- `FractionalZener`, `GeneralizedMaxwell`,
  `lower_to_prony` — the crate text already names
  **felt** as a target
- `WallPin` in `fs-phs`; slant / Bessel / chimney
  in `fs-duct`
- `Region` / `Chart` in `fs-geom` (the right
  abstraction; not yet a lumen extractor)
- `fs-contact` Hertz / Hunt–Crossley / FiniteGapPoint
- `fs-vfit`, plates, couple ports, BEM crate, LBM crate
- `fs-matdb` cards with provenance (not yet a
  brass–tissue–felt music pack)

Not yet, and we will not fake them:

- A certified mesh → `A(x)` extractor (centerline +
  inscribed radius + volume certificate)
- Valve-crook `ΔL` from CAD
- Offline lip / reed reduce that mints a valve card
  from a mesh + tissue card
- Wool-felt `Uniaxial` implementor + a cited felt
  coupon lowered through `FractionalZener`
- Multimodal horn as a default (plane-wave is not
  a trumpet claim)
- Unflanged 3-D exterior BEM in the live loop
  (bake a `Z_L` table offline)
- 3-D jet (unminted; flute/organ wait on a lab card)
- Honest mean flow (not a `k/(1−M²)` hack)

---

## 8. Parallel tracks

Tracks own their menus. **No shared M1.**

| Track | First slice (does not block others) |
| --- | --- |
| T-Wind | Char+holes+valve **and** TMM R gate; VFIT hold as second image |
| T-String | Modal ZOH gate (mostly exists); optional dispersive-WG bake-off |
| T-NL | KC on N modes |
| T-Bow | Friction island |
| T-Piano | `Uniaxial`+Prony felt island + 1 string + 1 board; Hunt–Crossley is the debug image |
| T-Voice | Two `A(x)` + 1-DOF **and** TMM formants; two-mass as second valve |
| T-EM | Faraday + split circuit |
| T-Brass | MM TMM + HB; MM char+lips as second image |
| T-Jet | Lab card + one pipe |
| T-Shell | Bar modal; shell later |
| T-Sched | Calls kernels; no physics |

Fuse **after that track’s bake-off**, per kernel.

---

## 9. Validation matrix (sketch)

Rows: fillings. Columns: images. Cells: QoIs or “off
menu.” Empty legal cell = not yet gated (cannot claim).
Refused cell = tried, failed, stays off the claim.

G0 passivity still applies to every *admitted* image.

---

## 10. Never

Samples; IR cabinets as truth; instrument crates; one
stepper for all fillings; **one image per filling as
policy**; Newton on linear pHS; coarse ladder as HF air
truth; FIR keymaps; formant EQ; forced Helmholtz saw;
in-loop BEM/LBM; deleting a passing orthogonal image
because another is newer; Hunt–Crossley as the
piano-felt law; `ManderConcrete` / `MenegottoPintoSteel`
as a hammer; a parallel “music felt” constitutive next
to `Uniaxial` / `FractionalZener`.

---

## 11. Open questions (local to a menu)

Thiran order (T-Wind char). When to hop char→VFIT.
Wool-felt `Uniaxial` envelope (which coupon); Prony
band vs hammer spectrum. Hunt–Crossley remains
**debug** until a bake-off against the hysteresis
island says otherwise. JCA/Biot linings are a
different claim, not an OQ on the hammer.
Lip pair vs beam (both on brass menu). Licensed `A(x)`.
Jet-card residual. Whether dispersive-WG string earns
a keep-both vs modal.

**No track blocked by another track’s OQ.**

---

## 12. Advice condensed

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
   split + Faraday.
9. **Selector is the product.** Gesture × claim ×
   budget → image set.
10. **Fuse each kernel after its own gate.** Parallel
    calendars.

---

## 13. Diff from v5 / v6

| v5 | v6 | v7 |
| --- | --- | --- |
| One RT stack per filling | Menu of orthogonal images | Same menus; **ingest law** spelled out |
| Implied winner | Bake-off; keep-both default | Unchanged |
| D18 as the wind law | D18 is *one* air image | Unchanged |
| “Extractors emit charts” | One-page geometry stub | Full mesh → cards → playing loop (trumpet §7.2, piano §7.5) |
| Felt = Hunt–Crossley vs “paper” | Two P-Contact images, equal | Felt **is** `Uniaxial` + `FractionalZener`/`Prony`; Hunt–Crossley is debug |
| T-Wind M1 = the char path | T-Wind ships TMM **and** char | Unchanged |

v4/v5 physics (d’Alembert, don’t Newton linear pHS,
Schur island, R(ω) for reeds) remain **entries on
menus**, not the constitution of every filling.

---

## 14. Planning-workflow next

Review should add a filling we still over-merged
(organ vs flute is already split by cardinality;
electric bass vs guitar is cards not tracks). Beadify
**menus**: T-Wind (TMM + char + R gate + slur),
T-String (modal gate + optional WG bake-off), and
T-Piano (`Uniaxial`+Prony felt island, not a new
constitutive crate) in parallel. Do not bead “the”
instrument stepper.

---

*End of v7. Constitution:
`COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md`. Shared ports.
Orthogonal interiors. **Menus, not winners.**
3-D + cards mint charts; charts play.*
