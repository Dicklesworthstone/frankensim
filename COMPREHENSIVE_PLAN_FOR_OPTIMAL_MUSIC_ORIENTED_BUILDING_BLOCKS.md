# Comprehensive Plan for Optimal Music-Oriented Building Blocks

> Nameless physics, emergent sound, two lanes, one set of
> equations. Real-time by *splitting integrators and
> structure-preserving reduction* — not by fitting a second
> invented waveguide, and not by running Gonzalez Newton on a
> linear bore.
>
> Companion to `COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md` and
> `COMPREHENSIVE_PLAN_TO_EXTEND_FRANKENSIM_TO_NEW_DOMAINS.md`
> §5.11. Does **not** amend the Decalogue.
>
> Status: **v3**. v1 = doctrine + naïve FIR bake. v2 = playable
> scatter + valve loop + ear-weighted gate. v3 = the same
> physics as the authority, cheaper *because* we respect
> linearity, pHS structure, and multiple time scales.
>
> Tags: `[S]` solid, `[F]` frontier, `[M]` moonshot.
> Date: 2026-08-14.

---

## 0. How to read this

Load-bearing: §1.4 (the category error), §2, §6 (split
stepper + 1-port Schur + parameterized reduction), §12
(M1 is now “linear duct cheap, valve small”).

If a section names “trumpet” or “ooh,” it is a *filling*,
never a crate.

---

## 1. Intent

### 1.1 What we want to hear

Notes whose pitch, timbre, attack, and *playing response*
come from materials + geometry + gesture + air: reed and
holes, lips and flare, jet and pipe, string and plate,
felt and soundboard, shell strike, vocal `A(x)`, pickup
and circuit and cone. No samples, no instrument types.

### 1.2 Success (all five)

1. **Honesty** — named primitive, material, geometry, tier
   or no-claim.
2. **Same equations** — the audible stepper is a
   discretization or Galerkin image of the authority
   network, not a separately invented delay line.
3. **Playability** — holes, `A(x)`, slide, bow force change
   without a rebake-per-key and without energy clicks.
4. **Perceptual adequacy** — ear-weighted residual +
   `fs-psycho` + listening. Metrics do not replace ears.
5. **Budget** — 48 kHz on one named Mac core for a *played
   phrase*, zero heap in the audio loop. Authority may be
   slow.

Fail (5) → scientific renderer. Fail (1)–(2) → synth.
Fail (3) → sampler of our own notes. Fail (4) → paper.

### 1.3 Workflows

W1 author → W2 authority solve → W3 bake (parameterized
*realization of the same network*, plus optional modal
image) → W4 play (split stepper) → W5 gate → W6 escalate
→ W7 `[F]` calibrate cards with adjoints (never ship the
recording).

### 1.4 The category error (why v2 was still not optimal)

We already write the bore as a **port-Hamiltonian LC
ladder** with linear Foster wall laws and a **tiny
nonlinear valve**. Then we time-step the *whole* thing
with Gonzalez approximate-Newton, as if the duct were a
von Karman plate.

For quadratic `H` and linear `(J,R,G)`, Gonzalez **is**
implicit midpoint: **one linear system per step**, or
**one factorization per control-rate geometry**. There is
no Jacobian to difference. Treating a linear bore as a
generic nonlinear DAE is why couple tests take minutes
and why a synth seems impossible.

v2’s “identify a scatter-chain from TMM at a pin
frequency” is the next-order version of
`TravelingWaveLine`: a *fit* that forgets
frequency-dependent `Zc` and then patches Foster on the
side. The honest scatter/WDF object is the **bilinear
transform of the same LC graph we already have**, not a
curve-fit.

**Optimal and more correct is the same statement:**

> Linear (or linear-at-control-rate) pHS → structure-
> preserving *linear* stepper or Galerkin. Nonlinear
> islands (valve, felt, bow, jet, KC, von Karman) stay
> small and implicit. Spatial/model order is
> *goal-oriented* toward the listener, not uniform.

That is physically the telegraph + ZK + Bernoulli system.
It is also what real-time circuit simulators do.

---

## 2. Doctrine

D1 nameless primitives. D2 materials+geometry, recordings
test-only. D3 two lanes, three rates (gesture / control /
audio). D4 reduce, do not invent a second physics; WDF
and modal banks are *images* of our network. D5
correctness before fusion. D6 Franken-only + `fs-io`
quarantine. D7 certificates. D8 named RT budget. D9
`fs-psycho` judges, does not sing. D10 moonshots flagged.
D11 time-varying ports stay passive. D12 no `fs-bake` /
`fs-synth` crate. D13 calibrate cards, do not clone WAVs.

**D14. Do not Newton a linear pHS.** `[S]`
If `H` is quadratic and `(J,R)` do not depend on `x`,
the stepper is a linear solve (or a closed ZOH in a
modal basis). Finite-difference Newton on that system is
a bug.

**D15. Same graph, coarsened — not a new graph.** `[S]`
Performance spatial meshes are child charts of the
authority mesh (fewer cells, fewer Foster branches),
with a goal-oriented residual (first impedance peaks,
ERB-weighted `Z`, radiated `p`). A Kelly–Lochbaum chain
is legal only as the bilinear/WDF of that child LC
graph.

**D16. Multi-rate by stiffness, not by folklore.** `[S]`
Felt impact (~µs), string/bore (48 kHz), board radiation
(can be coarser), gesture (100 Hz). Substep only the
stiff island. One rate for everything is how pianos miss
RT *and* how impacts go unstable.

**D17. State must lift when the image changes.** `[S]`
Changing `σ` or the reduced basis is a change of
coordinates. Define `x_r⁺ = Π(σ⁺) x` (or `V(σ⁺)⁺ x`)
so energy and the audible wave do not jump. Coefficient
interpolation without a lift is how slurs click.

---

## 3. Relation to existing plans

Constitution wins. `fs-acoustics` in the new-domains plan
is the 3-D/NVH stack we *consume* (exterior tables,
shells). Live `fs-phs` already has Gonzalez, Galerkin
reduction, Bernoulli, LC/ψ ducts — v3 **uses that
algebra**, it does not replace it with a synth engine.
Mega-fused doctrine applies after the split stepper
exists, on the remaining hot inner products.

---

## 4. Current state, reread

We have the *right objects* (pHS duct, TMM, reed, string,
plate) and the *wrong default integrator* for the linear
majority of the state. Bessel `k` is now Helmholtz. Ernoult
frozen peaks work. We cannot slur a hole cheaply. FIR `R(ω)`
is a frozen-note renderer. `reduce_galerkin` exists and is
not the audio path. `fs-dwr` exists and is not used to
coarsen bores. `fs-vfit` already passivity-repairs radiation
filters.

Missing physics (still): multimodal tube, lip pair, jet
card, felt, pickup+circuit+speaker, `A(x)` chart, thin
shell, unflanged `Z_L` table, mean flow.

---

## 5. Atlas

Media, geometry charts, valves, strings, plates, contact,
pickup/circuit/speaker, radiation, stochastic physics:
as in v2. Changes below are the ones that matter.

### 5.1 The linear acoustic network (authority)

The 1-D duct **is** the telegraph system

```text
∂p/∂x = −Z'(ω) U ,   ∂U/∂x = −Y'(ω) p
```

with ZK `Z',Y'` (Bessel in frequency, Foster/Cauer in
time), plus T-junctions, pad `C`, mouth `Z_L`, walls.
ODE `acoustic_chain` is a spatial finite-volume / Cauer
of that system. TMM is the same system in `ω`.

**These two must stay certified against each other**
(ear-weighted `Z_in`, peak cents). Today they can drift
(convention bugs already taught us that). One composition
picks **one** time-domain authority (pHS) and uses TMM
as the frequency oracle, or the reverse for frozen linear
maps.

### 5.2 Child discretizations (performance, same graph)

| Image | What it is | When |
| --- | --- | --- |
| **Fine LC + linear midpoint** | Same pHS, D14 stepper | Authority-quality TD, or RT if `n` small |
| **Coarse LC / WDF** | Bilinear of the child graph | Time-varying holes, vowels, slides |
| **Modal / Krylov / IRKA** | Galerkin or pH-IRKA of the *linear* operator | Frozen or slowly varying geometry; strings; plates; shells |
| **Parameterized MOR `[F]`** | `V(σ)`, interpolatory pMOR | Fast `σ` *and* few global modes |
| **Driving-point 1-port** | Schur complement of the linear network at the valve | Always, as the loop interface |
| **FIR** | IFFT of a *frozen* 1-port | Held note, radiation, pickup |

Scatter-chain in v2 **is** the WDF/bilinear row, once it
is derived from the LC graph rather than identified at
one `ω_pin`.

### 5.3 Nonlinear islands

Bernoulli, MSD, two-mass, Hunt–Crossley, hysteretic felt,
Stribeck, `JetOscillator`, KC `H`, von Karman `H`.
These keep Gonzalez or a small Newton / K-method.
They never swallow the linear duct’s Jacobian.

### 5.4 `ImplicitPortLoop` = Schur + small island

At each audio sample (geometry frozen on the control
interval):

```text
p = C x + D u          # linear 1-port of the duct
ẋ = A x + B u          # or scatter step; A,B,C,D from factorization
u = F(p_up − p, y)     # Bernoulli / jet
ẏ = f_island(y, p)     # reed/lip/fold ODE
```

`A,B,C,D` (or the WDF adaptor list) are **rebuilt at
control rate**, not identified from a pin frequency.
The nonlinear solve is 1–3 scalars, not 80.

If `F` is monotone dissipative, a scalar Newton or a
WDF adaptor is enough. Gate the discrete loop against
full Gonzalez on a tiny fixture (one cell + valve)
*and* against the fine network on Ernoult.

### 5.5 Goal-oriented coarsening `[S/F]`

Use the listener (or the first `N` impedance peaks) as
the QoI. `fs-dwr` / adjoint of `Z_in` or radiated `p`
marks cells and Foster branches that do not pay. Coarsen
there. This is how a 40-slice cone becomes 8 *for the
ear*, with a receipt, instead of a vibe.

### 5.6 Multi-rate islands

| Island | Rate |
| --- | --- |
| Felt / impact | sub-audio (fixed substeps or event) |
| Valve loop + bore / string | 48 kHz |
| Plate/shell radiation, ISO path | block or 2× downsample if gated |
| `σ`, `A(x)`, gas, slide | control (0.5–2 kHz) |
| Score / MIDI | gesture |

### 5.7 Harmonic-balance authority `[F]`

For *periodic* brass/reed lock, frequency-domain HB of
the valve + multimodal `Z(nω)` is the right *existence*
oracle (does this lip + this flare lock at this pitch?).
Cheaper than a long Gonzalez transient. Attacks still
need TD. Two authorities, one physics.

---

## 6. Architecture

```
        geometry + materials + gesture
                     │
                     ▼
           pHS network + TMM oracle
           (fine LC, Foster, valve)
                     │
         child mesh / Galerkin / pMOR
                     │
         ┌───────────┴───────────┐
         │  control tick         │
         │  factor A(σ), lift x  │
         └───────────┬───────────┘
                     │
         audio: 1-port linear step
              + 1–3 dim island
                     │
                  pascals
```

### 6.1 Control tick (the real “bake update”)

On `σ`, `A(x)`, `T`, `μ` change:

1. Rebuild child network or `V(σ)`.
2. **Lift state** (D17).
3. Factor the linear midpoint / scatter map, or refresh
   modal `e^{A Δt}`.
4. Check D11 on the discrete energy.

This is O(n) to O(n³) at 1 kHz for n ~ 20–80 — cheap
compared to Newton-every-sample on the same `n`.

### 6.2 Audio tick

Linear 1-port step (matrix–vector or WDF sweep) +
scalar/small island. Fuse those two (mega-fused target).
No heap. SoA. SIMD on modal/FIR inner products only
after this is green.

### 6.3 Which image when

| Geometry in time | Prefer |
| --- | --- |
| Holes, vowels, slide | Child LC / WDF + lift |
| Held note, string, plate, bell | Modal / IRKA 1-port |
| Fast `σ` *and* few global peaks | pMOR `[F]` |
| Studio | Fine LC + linear midpoint, Gonzalez on islands |

### 6.4 Quality knob

Same gesture. `{fine, child, modal_lo}`. Each is a
gated image, not a different constitutive law.

### 6.5 Polyphony scheduler

Steal decaying modal tails; drop multimodal lines;
never drop the foreground valve loop; substep felt only
on attacks.

### 6.6 Inverse calibration `[F]`

Adjoint/DFO on cards (`A(x)`, wall `r`, felt, lip rest
length) so *authority* matches a cited Z(ω) or held
note. Runtime still plays the network.

### 6.7 Hardware truth (audio, not GEMM)

Delay lines and modal states are **bandwidth**. Keep
one voice’s delays contiguous. Do not parallelize inside
a 3 ms line. Unified memory: avoid extra transposes.
Fingerprint cycles/sample; ledger or shut up.

---

## 7. Geometry

Cited meshes/tables → `fs-io` → extractors → charts →
**fine network** → **child** under DWR/ear gate.
`A(x)` papers skip meshing. No scraped CAD.

---

## 8. Family fillings

Unchanged morally from v2; the *stepper* underneath
changes.

- **Reed winds:** valve island + child LC/WDF of the
  bore + `σ(t)`. Mouthpiece = sections. MVP = slur on
  Ernoult geometry with D14/D17, not a new FIR.
- **Brass:** lip island + multimodal child + HB lock
  oracle + `Z_L` table. Slide = delay length at control
  rate + lift.
- **Jet:** lab-minted island or refuse.
- **Strings:** already modal/KC; KC keeps Gonzalez (true
  nonlinear `H`); linear strings are ZOH (D14).
- **Bow:** Stribeck island + string.
- **Piano:** felt island *substepped* (D16) + modal
  strings + modal board.
- **Bells:** modal shell + strike island.
- **Voice:** glottal island + `A(x)` child network;
  source–filter interaction is the 1-port. Tissue =
  `WallPin`. Ooh/aah = two `A(x)` cards.
- **Electric:** Faraday + circuit island + speaker
  1-port. Circuit may be stiff → same split (linear RLC
  factored, nonlinear device as island).

---

## 9. Validation

G0–G6 as v2 (passivity, Ernoult, ear-weighted, psycho).

**G7 played phrase** (slur / vibrato / crescendo): no
click, no extra ring, supply-rate holds.

**G8 split-vs-monolithic:** on a fixture, linear-midpoint
+ island versus full Gonzalez: energy and audible peaks
inside a declared band. This *is* the proof that D14 did
not change the physics.

**G9 child-vs-fine:** ear-weighted `Z` and G7 on the
coarsened graph.

Reference WAVs remain test-only.

---

## 10. Never

Samples, cabinet IRs as truth, instrument crates, in-loop
CFD/BEM, brightness EQ, one-pole bores, fake mean flow,
`2^N` FIR fingering banks, Newton on linear pHS, a second
waveguide theory fitted at one frequency.

---

## 11. Open questions

| ID | Question | Blocks |
| --- | --- | --- |
| OQ-1 | WDF sweep vs factored implicit midpoint vs modal 1-port as the *first* M1 linear image? Implement midpoint first (it is already Gonzalez for quadratic H). Compare WDF as a child. | M1 |
| OQ-2 | Control rate vs factor cost on the Mac fingerprint. | M1 |
| OQ-3 | Lift operator: energy projection vs interpolation of distributed `p,U`. | M1 |
| OQ-4 | Licensed `A(x)` sets. | M2 |
| OQ-5 | Two-mass vs beam lips. | M5 |
| OQ-6 | Felt constitutive paper. | M4 |
| OQ-7 | HB brass oracle: how many harmonics. | M5 |
| OQ-8 | Listening budgets, then tighten. | M1 |

---

## 12. Phases

**M0** — Document D14–D17. Add a linear-pHS fast path
*detection* (quadratic H, state-independent J,R) in
`fs-phs::step` or a sibling `step_linear`. This is
product, not ceremony: the same API, a cheaper correct
integrator.

**M1 — Load-bearing** `[S]`
Reed + cylinder + live vents:

1. Linear midpoint (or factored implicit) on the existing
   LC+Foster graph (D14).
2. Valve as the only Newton/K-method island.
3. Control-rate `σ` + state lift (D17) + D11.
4. G7 slur + G8 vs Gonzalez + Ernoult frozen G2.

Optional child: bilinear WDF of the *same* cells, G9.

If M1 is late, stop adding Foster branches for “RT.”

**M2** — `A(x)` + vowels + cones on the same stepper.

**M3** — Pickup + circuit island + speaker (circuit uses
the same split).

**M4** — Felt island + D16 substep + board modal.

**M5** — Multimodal child + lips + `Z_L` table + HB
oracle.

**M6** — Jet island from the lab, or refuse.

**M7** — Shell modal.

**M8** — Fuse the 1-port×island inner loop; SIMD modal/
FIR; scheduler; fingerprint ledger. **After G7/G8.**

**M9** — Gesture skin (`σ(t)`, pressure). Still no
instrument crate.

```
D14 step_linear → island valve → lift σ → G7/G8 = M1
        ↘ A(x) → M2
        ↘ circuit split → M3
        ↘ felt substep → M4
multimodal + HB → M5
jet card → M6
shell → M7
G7 green → M8 fusion
```

---

## 13. Work packages (M1 first)

- `step_linear` / detect quadratic+LTI pHS; test = bit or
  tight residual vs Gonzalez on a lossless cylinder.
- Schur 1-port of `acoustic_chain` at the reed face.
- Valve island (K-method or 1–2 Newton) + G8.
- Control-rate `σ` + lift + D11 ramp + G7 slur.
- Ear-weighted `Z` residual helper (ERB + peak align).
- Child WDF only after midpoint path is green (G9).
- DWR coarsen cells against first peaks `[F]` after M1.

Do not bead M5–M7 until M1 has slurred.

---

## 14. Advice condensed

1. **The bore is linear. Stop Newtoning it.** That is the
   highest-leverage correctness *and* speed move in the
   repo today.
2. **WDF/scatter is the bilinear of our LC graph**, not
   a pin-frequency fit. Fitting is how we re-invent
   `TravelingWaveLine`.
3. **Galerkin / IRKA for frozen geometry; spatial child
   for moving holes.** pMOR later if both are true.
4. **Condense the duct to a 1-port; implicit only the
   valve/felt/bow.** Schur + island.
5. **Lift state when `σ` or `V` changes.** Passive
   coefficients are not enough.
6. **Coarsen with an adjoint/ear QoI**, not a slice
   count guessed in a comment.
7. **Substep stiff contact.** One rate is neither
   accurate nor fast.
8. **HB for brass lock, TD for attacks.**
9. **Meshes → extractors → fine network → child.**
   Not 48 kHz FEM, not handmade formants.
10. **Fuse last**, on the 1-port×island product, after
    G7/G8.

---

## 15. Diff from v2 (reviewers)

| v2 | v3 |
| --- | --- |
| Identify scatter at `ω_pin` | Bilinear/WDF of the *same* LC; or don’t scatter yet |
| Gonzalez vs “cheap image” | Gonzalez for nonlinear `H`; linear midpoint for the bore |
| Bake as a new object family | Child chart + factored linear map + small island |
| Valve loop as a glue primitive | Schur 1-port + island (the loop is the interface) |
| Passive coefficient interp | That **plus** state lift (D17) |
| Uniform then ERB residual | Goal-oriented spatial coarsening too |
| One audio rate | Multi-rate by stiffness |
| Brass = multimodal scatter | Plus HB lock authority |
| M1 = scatter + slur | M1 = D14 linear duct + island + lift + slur |

v2 was right that frozen FIR cannot play, and that the
valve–bore loop is the instrument. It still proposed a
*second* waveguide. v3 says: **play the network we
already wrote, with the integrator it always deserved.**

---

## 16. Planning-workflow next

v3 is the first plan that is simultaneously more
accurate and more likely to be real-time. Next review
should attack OQ-1 (midpoint vs WDF first) and the
lift operator (OQ-3). Beadify **only** M0–M1
(`step_linear`, 1-port, island, lift, G7/G8).

---

*End of v3. Constitution:
`COMPREHENSIVE_PLAN_FOR_FRANKENSIM.md`. The cheap stepper
is a theorem about linear pHS, not a new instrument
engine.*
