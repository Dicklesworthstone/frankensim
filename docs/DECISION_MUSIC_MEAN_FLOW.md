# Decision Record: Mean Flow in the Music Program — Escalation Ladder

| | |
|---|---|
| **Record** | `DECISION_MUSIC_MEAN_FLOW` |
| **Status** | Ratified (bead `frankensim-music-v8-root-3ez8g.10.5`, owner-delegated execution) |
| **Date** | 2026-08-23 |
| **Scope** | `fs-duct`, `fs-phs` acoustic surfaces (music program); T-Jet child `.10.5` |
| **Supersedes** | nothing |
| **Feeds** | future convected-model card lane (rung B); fs-vvreg amplitude fixtures |

## Context

Real wind instruments carry mean flow: the player blows through them, and at
forte the quiescent-medium assumption behind every current `fs-duct`
transmission-line/modal claim and the `fs-phs` acoustic pHS components
(Helmholtz resonator, `acoustic_chain` LC ladders) strains. The plan's NEVER
list records the standing temptation: the `k/(1-M^2)` convected-wavenumber
substitution — cheap, plausible-looking, and wrong outside its narrow
derivation class. Deadline pressure is exactly when such a hack gets committed,
so the escalation ladder is written down BEFORE anyone needs it.

This record quantifies where mean flow actually matters, fixes the ladder,
names the falsifier that forces escalation on evidence, and patches the two
crates' no-claims sections. It is a DECISION RECORD: no solver work lands here.

## 1. Where mean flow matters, and how much (all Estimate)

Three distinct regimes, quantified by scaling laws (source classes named;
licensing-first: no transcribed tables, graphs, or numeric datasets —
magnitudes are order estimates derived from the cited scalings):

**(a) Through-flow (steady blowing component in the bore).**
Woodwind/brass through-flow Mach numbers at forte are typically
`M̄ ≈ 0.01–0.08` (blowing-pressure classes in the wind-literature range
translate to these in-bore values; the jet itself lives outside the bore
proper in flutes and is the card lane's business). Convected duct acoustics
shifts resonance frequencies by first order in `M̄` only when the termination
is asymmetric (radiating end): fractional shift `Δf/f = O(M̄)` per half-wave
accumulated phase; symmetric closed configurations see `O(M̄²)`. Concretely:
at `M̄ = 0.05` a third-mode bore resonance moves single-digit Hz (a few to a
few tens of cents depending on Q), and input-impedance MAGNITUDE minima
(loss-dominated anti-resonances) shift by percent-to-tens-of-percent, while
maxima barely move. Below `M̄ ≈ 0.01` all effects sit under typical
measurement scatter. [Source class: standard uniform-flow duct acoustics —
convected wavenumber pairs `k± = k/(1∓M̄)`, Morse & Ingard-style treatments;
Eversman's duct-acoustics review literature.]

**(b) Amplitude-dependent (acoustic Mach) effects.**
Peak acoustic Mach in bores at forte reaches `M_a ≈ 0.02–0.15` locally near
fipples, register/valve galleries, and tone holes. Two mechanisms dominate
long before waveform steepening matters (`L_shock = c/(β ω M_a)` stays
multi-meter for these values): (i) vortex-shedding-dominated end corrections
at tone holes and mouth windows, which go nonlinear once the acoustic
displacement approaches the hole radius — i.e. `ξ = M_a c/ω ~ mm` around
`M_a ≈ 0.02–0.06` in the 300–1000 Hz band; (ii) jet-drive interaction at the
embouchure/labium, ALREADY owned by the card lane. Effect size: impedance
MINIMA deepen/shift by percent-to-tens-of-percent; resonance peaks move
little. [Source class: wind-instrument aeroacoustics reviews
(Hirschberg-school), Fletcher & Rossing-scale textbook treatments.]

**(c) Vena contracta at tone holes at high amplitude.** Same physical family
as (b)(i); magnitude as above; modeling home is the tone-hole card lane, not
the transmission line.

Everything else in the program's current claims (small-amplitude impedance
surveys, pianissimo–mezzoforte synthesis, transfer-function archiving) sits
one to two orders of magnitude below regime onset.

## 2. The escalation ladder

**Rung A — quiescent medium with AMPLITUDE NO-CLAIMS (current state, default).**
Every `fs-duct` and `fs-phs` acoustic claim states: quiescent medium
(`M̄ = 0`), valid for peak local acoustic Mach `< 0.01` unless the owning
surface states a different bound. This is now WRITTEN into both crates'
no-claims sections (patched by this record), not merely implied.

**Rung B — convected propagation as a NAMED model (admitted, not built here).**
Uniform, subsonic, constant `M̄`; rigid smooth walls; plane or `m = 0` modes;
passivity-preserving formulation with upstream/downstream wavenumber pair
`k± = k/(1∓M̄)` AND its matching impedance scaling, entering ONLY as an
explicit, named model choice (`ConvectedUniformV1`-class card) with its
derivation class printed in the card. It must never appear as a silent
multiplier inside existing quiescent functions. Losses couple through the
convection factorization chosen by that card; the card owns its own
validity envelope and refusals.

**Rung C — beyond uniform flow ([M], feature-gated, out of near-term scope).**
Shear flow (LEE/APE-class) or direct CAA. Requires its own frozen model card,
registered corpus, and ratification before any implementation. Nothing in
this program schedules Rung C work.

Escalation moves ONE rung at a time and only on the falsifier below.

## 3. The falsifier that forces escalation

A named measured-deviation class, checked against the registered corpus
machinery (`fs-vvreg` acoustic rows), NOT against enthusiasm:

> Fixture: a registered bore/input-impedance corpus row measured at TWO
> authored drive levels spanning at least a factor of five in level.
> Trigger: the quiescent model's deviation from measurement exceeds its
> authored amplitude envelope at the HIGH level while remaining inside it
> at the LOW level, reproduced on an independent rerun of the same fixture.
> That monotone-in-level signature distinguishes genuine weakly-nonlinear /
> convected physics from ordinary model error (which does not track level).

On trigger: escalate the affected surface one rung (A→B) with a model card;
the affected claims' envelopes shrink to the passing region meanwhile. No
trigger, no machinery — even if profiling makes someone eager.

## 4. CONTRACT language patched by this record

- `crates/fs-duct/CONTRACT.md` — the straight-smooth-walls and modal-module
  bullets now state the quiescent/amplitude boundary explicitly and point
  here instead of leaving "no mean flow" as bare absence.
- `crates/fs-phs/CONTRACT.md` — new no-claims bullet for the acoustic pHS
  components (lumped Helmholtz resonator, `AcousticTap`s, `acoustic_chain`
  inviscid ladders): same quiescent + `M_a < 0.01` boundary, same pointer.

## NEVER list entry

**FORBIDDEN**: substituting `k → k/(1-M²)` (or any equivalent single-factor
rescale of `k` or `Z`) inside quiescent models as a "mean-flow fix".
Reasons, each independently disqualifying:

1. It is the symmetric-closure special case posing as the general result:
   real instruments have an asymmetric radiating end, so the leading shift is
   `O(M̄)` with endpoint-dependent sign — the substitution gets both wrong.
2. It double-counts: wall-loss and end-correction terms were linearized about
   quiescence; rescaling `k` underneath them corrupts their first-order
   structure instead of convection-correcting it.
3. It breaks the passivity/energy bookkeeping the `fs-phs` layer guarantees
   by construction; there is no port-Hamiltonian realization of the hack.
4. Deployed without a class statement it has no refusal surface — it cannot
   be falsified, which violates the program's certificates-over-vibes law.

Any PR introducing it is refused in review regardless of numerical plausibility.

## No-claim boundaries of THIS record

Literature magnitudes are order-of-effect estimates tied to named scaling
laws, not transcribed measurements (licensing-first). No convected model is
implemented or claimed here. Rung B admission changes nothing until a card
lands and a corpus falsifier fires. This record certifies DECISIONS and
boundaries, not numbers on instruments.
