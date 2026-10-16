# Measured Cycle-Time Baseline Dossier — Incumbent Electronics-Cooling Workflow

<!-- marker: fs-wedge-cycle-time-baseline-dossier-v1 -->

Assembled: 2026-07-22 · Bead: frankensim-extreal-program-f85xj.1.3
Supersedes: the `baseline_days = 5.0` scoping placeholder in `CHT_BASELINE`.

This dossier is the provenance record behind the `MeasuredCycleTimeBaseline`
constant in `fs-wedge`. Every load-bearing number below was traced to a
reachable primary or near-primary source during assembly, with the sentence
containing the number quoted. Numbers that could not be verified are listed in
the Exclusions section and are NOT load-bearing. The baseline is an ESTIMATED
quantity assembled from published sources — not a measurement we executed —
and the Rust record labels it accordingly.

---

## 1. Protocol (what makes a figure admissible)

1. **Task frame.** The representative incumbent task is: thermal analysis of an
   electronics assembly (PCB + components + heat sink, forced convection /
   conjugate heat transfer) taken from dirty CAD to an engineering report in a
   commercial tool (Simcenter Flotherm/FLOEFD, Ansys Icepak/Fluent, Cadence
   Celsius, 6SigmaET class), one design iteration.
2. **Step decomposition.** Six steps: CAD preparation (import/cleanup/
   defeaturing), meshing, solver setup (materials/BCs/power maps), solve,
   post-processing, report assembly.
3. **Admissible sources**, in descending authority: (a) an incumbent-workflow
   run executed and timed by us (none yet — see §6); (b) government-lab or
   peer-reviewed time studies; (c) practitioner surveys (self-report; vendor
   sponsorship flagged); (d) vendor case studies (admissible only with the
   marketing-bias caveat attached); (e) vendor marketing claims (context only,
   never load-bearing).
4. **Every admitted figure carries**: the figure and units, what exactly was
   measured (engineer-time vs wall-clock, tool, problem size), full citation
   with URL, the verbatim sentence containing the number, and a bias class.
5. **Envelope derivation.** Per-step low/high bounds are set from the
   admissible figures for a task of the representative frame's scale;
   defense-lab-scale figures (Sandia DART) bound context but do not set the
   representative high unless no closer figure exists. The total envelope is
   the sum of step bounds; the low total is a simultaneous-best-case and is
   deliberately conservative as a kill-criterion denominator.
6. **No false precision.** The record stores a range, never a point value.
   The retired 5.0-day placeholder happens to lie inside the range; that is a
   consistency observation, not a validation of the placeholder.

## 2. Per-step evidence

### 2.1 CAD preparation (import, cleanup, defeaturing) — low 1 h, high 16 h

| Figure | Source | Measured | Bias class |
|---|---|---|---|
| 39% of thermal engineers spend >1 h importing CAD; 24% <10 min | 6SigmaET (Future Facilities) "State of Thermal" survey, n>170, via engineering.com, 2018-01-25, <https://www.engineering.com/70-of-engineers-left-out-in-the-cold-by-thermal-simulations/> — "39 percent spend over an hour importing CAD data, while 24 percent claimed they can run the same task in under 10 minutes." | Self-reported engineer-time per model, tool-agnostic | Vendor-run survey (6SigmaET differentiates on this pain point) |
| 30% spend >1 day building models; 37% <1 h | Same source — "only 37 percent of engineers spend less than an hour building their models, while 30 percent spend more than a day." | Self-reported model build incl. geometry idealization | Vendor-run survey |
| Preprocessing = 38% of total simulation time | Tech-Clarity, "Addressing the Bottlenecks of FEA Simulation," 2016-07-19, <https://tech-clarity.com/simulation_bottlenecks/5467> — "Preprocessing is the most time consuming part of the simulation process, taking up 38% of total simulation time." n>160 manufacturers. | Self-reported share, all-FEA (not thermal-specific) | Analyst survey, typically vendor-sponsored |
| Defeaturing 48.8 h + decomposition 36.8 h initial engineer-hours (heat-transfer problems) | Sandia SAND2005-4647 "DART System Analysis," 2005, Table 5, <https://www.osti.gov/servlets/purl/876325> | Structured time study, 23 analysts / 34 datasets; heat-transfer n=5 | Gov-lab study; defense-scale problems — context/upper bound, not representative-task scale |

Bound rationale: 1 h floor for clean-CAD import plus minimal cleanup
(6SigmaET fast quartile); 16 h (2 working days) ceiling for a dirty assembly,
consistent with the 30%->1-day model-build row while staying below
defense-lab scale.

### 2.2 Meshing — low 0.5 h, high 24 h

| Figure | Source | Measured | Bias class |
|---|---|---|---|
| 41% spend >1 h gridding; 10% >1 day | 6SigmaET 2018 (as above) — "41 percent also said that they typically spend over an hour gridding their designs, while 10 percent spend more than a day." | Self-reported per design | Vendor-run survey |
| 1.9M-cell mesh in <3 h wall-clock (FpBGA 208 package + test board, Flotherm XT, 12-core Xeon) | Thales Corporate Engineering case via ETS Solution Asia (Siemens reseller), <https://www.etssolution-asia.com/blog/simcenter-flotherm-accurate-simulations-for-electronics-thermal-design> — "1.9-million-cell mesh in less than three hours of computation on a 12-core Intel Xeon processor." | Wall-clock mesh compute, one detailed package | Vendor-channel case study (favorable framing) |
| Meshing 19.6 h + mesh manipulation 11.3 h initial engineer-hours (heat transfer) | Sandia SAND2005-4647, Table 5 | Time study (defense scale) | Gov-lab; context/upper |
| Surface meshing "minutes to more than two weeks"; volume "minutes to days"; no participant meshed in one attempt | GMGW-1 summary, AIAA 2018-0128, <https://ntrs.nasa.gov/api/citations/20180006182/downloads/20180006182.pdf> | Participant questionnaires, n=51, clean aerospace research geometry | Peer-reviewed workshop; aerospace, context only |

Bound rationale: 0.5 h floor for automated hex-dominant gridding in a
CAD-embedded tool; 24 h ceiling covers the 10%->1-day survey tail for
electronics-scale models.

### 2.3 Solver setup (materials, BCs, power maps) — low 0.5 h, high 16 h

| Figure | Source | Measured | Bias class |
|---|---|---|---|
| 39% spend >1 h defining material properties | 6SigmaET 2018 — "39 percent spend more than an hour defining properties" | Self-reported | Vendor-run survey |
| 23% spend >1 h on boundary conditions | 6SigmaET 2018 — "23 percent must spend over an hour setting boundary conditions." | Self-reported | Vendor-run survey |
| Parameter assignment 19.4 h; simulation-model assembly 55.4 h (heat transfer, initial) | Sandia SAND2005-4647, Table 5 | Time study (defense scale) | Gov-lab; context/upper |
| Package thermal model: "days of effort" in CAD vs minutes with a generator | Siemens Simcenter blog, 2023-04-11, <https://blogs.sw.siemens.com/simcenter/creating-a-cpu-from-scratch/> | Anecdote, one component | Vendor blog; asserted baseline |

Bound rationale: 0.5 h floor when material/component libraries exist; 16 h
ceiling for a board-level power-map + interface-card setup assembled by hand.

### 2.4 Solve — low 0.5 h, high 48 h (wall-clock)

| Figure | Source | Measured | Bias class |
|---|---|---|---|
| 66% spend up to a day or more per solve; 14% <30 min | 6SigmaET 2018 — "66 percent of thermal engineers spend up to a day or more solving their simulations—and yet 14 percent said that they can perform this task in under 30 minutes." | Self-reported wall-clock | Vendor-run survey |
| <4.5 h to full convergence, 1.9M cells, 12-core Xeon, 10 GB RAM | Thales / Flotherm XT case (as §2.2) — "The powerful solver needed less than 4.5 hours to reach full convergence on the same processor using 10 Gb of memory." | Wall-clock, one detailed package + board | Vendor-channel case |
| Icepak transient runs 730–4800 s (~12–80 min), 0.3–3M cells | Padmanabhan (ZF), SAE 2025-01-5073, author MS <https://arxiv.org/pdf/2606.11226>, Table 3 | Wall-clock, module scale, laptop-class host | Peer-reviewed |
| Data-hall CHT steady solve 0.44–5.95 h (32→1 cores), 10M cells, 6SigmaDCX | Wang et al., "Kalibre," ACM BuildSys '20, <https://arxiv.org/pdf/2001.10681>, Table 2 | Wall-clock per steady solve | Peer-reviewed; data-center scale |
| RNG k-ε chassis CHT: "two days of continuous runs" (900k cells, 2008 hardware) | Öztürk & Tari, IEEE TCAPT 31(3), 2008, <https://users.metu.edu.tr/itari/OzturkTariIEEETCAPT2008.pdf> | Wall-clock | Peer-reviewed; dated hardware |

Bound rationale: 0.5 h floor (module-scale steady case on current hardware);
48 h ceiling (turbulence-model-heavy chassis CHT; counted 1:1 into cycle time
because a running solve blocks the iteration regardless of engineer
attention).

### 2.5 Post-processing — low 1 h, high 16 h

| Figure | Source | Measured | Bias class |
|---|---|---|---|
| Post-processing 123.4 h initial / 245.3 h expected = 37% of total engineer-time — the LARGEST single step for heat-transfer problems | Sandia SAND2005-4647, Tables 5/7/8 | Time study (defense scale) | Gov-lab; strongest available signal that post-processing, not meshing, dominates thermal analysis |
| Thermal-to-structural results handoff: half a day of manual scripting per transfer | Siemens Mynaric case study, <https://resources.sw.siemens.com/en-US/case-study-mynaric/> | Incumbent manual baseline in a vendor case | Vendor case study |

Bound rationale: 1 h floor for standard temperature-map review; 16 h ceiling
scaled down from the DART fraction to representative-task size.

### 2.6 Report assembly — low 1 h, high 8 h

Least-documented step in public sources; no independent electronics-specific
figure was found (explicit evidence gap). DART folds report generation into
its post-processing step I. The bounds are set from the DART post-processing
share net of visualization review, and are flagged as the weakest row in the
envelope. Improving this row is exactly the kind of update the protocol's
review cadence exists for.

## 3. Derived envelope

- Per-step sums: **low 4.5 engineer-hours**, **high 128 hours** (solve counted
  as wall-clock).
- At the record's explicit 8 h/working-day conversion: **low ≈ 0.6 working
  days, high = 16 working days** for one first-pass iteration.
- Re-iterations (design change → re-solve) are cheaper than first passes
  because mesh/setup are partially reusable; the envelope covers the first
  full pass. Sources on re-iteration cost are vendor-cases only (Thales,
  GOLDTek) and were not used to narrow the envelope.
- The retired placeholder (5.0 days) lies inside [0.6, 16]; the range, not the
  point, is the citable object.
- Kill-criterion consequence (target ≥3× reduction): the conservative verdict
  uses the LOW bound as denominator — FrankenSim's per-iteration cycle time
  must be ≤ 0.2 working days (~1.5 h) to claim "met" against even the fastest
  incumbents; ≤ 5.3 days only clears the HIGH bound and yields
  "indeterminate," not "met."

## 4. Genealogy note: the "meshing is 80%" folklore

The widely-repeated "geometry-to-mesh is 70–80% of analysis time" claim traces
to exactly one measurement: Sandia's DART study (SAND2005-4647, 2005), whose
community model gives geometry 57% / meshing 20% / everything else 23%
(steps A–E = 77%). Its most-cited appearance — Hughes, Cottrell & Bazilevs,
CMAME 194 (2005) p. 4136, "It is estimated that about 80% of overall analysis
time is devoted to mesh generation…" — carries **no citation** in the primary.
The Sandia pie entered the isogeometric-analysis literature in 2010 as a
courtesy figure credited to Blacker/Hardwick/Clay, i.e., the DART data
transmitted as personal communication. NASA CFD Vision 2030 (NASA/CR-2014-
218178) states meshing is a "principal bottleneck" and "the dominant cost in
terms of human intervention" but **never states a percentage**. Crucially for
this vertical, DART itself says heat transfer is the exception: "Thermal
analysts therefore spend less time on geometry and meshing" (geo+mesh 41%,
post-processing 37%). Any FrankenSim marketing that repeats "80% meshing" for
the cooling vertical would be citing folklore against our own evidence.

## 5. Source-quality summary

| Source | Year | Type | Verified | Load-bearing |
|---|---|---|---|---|
| Sandia SAND2005-4647 (DART) | 2005 | Gov-lab time study | Full PDF read | Yes (context/upper bounds; fractions) |
| 6SigmaET "State of Thermal" via engineering.com | 2018 | Vendor-run survey, n>170 | Article fetched, quotes extracted | Yes (only electronics-thermal per-step survey) |
| Tech-Clarity FEA bottlenecks | 2016 | Analyst survey, n>160 | Page fetched | Yes (preprocessing share) |
| Thales / Flotherm XT via ETS Solution Asia | ~2018 | Vendor-channel case | Page fetched | Yes (concrete mesh+solve wall-clock) |
| GMGW-1, AIAA 2018-0128 | 2018 | Peer-reviewed workshop | NTRS PDF read | Context only (aerospace) |
| Wang et al. Kalibre, BuildSys '20 | 2020 | Peer-reviewed | arXiv PDF read | Yes (solve scaling) |
| Öztürk & Tari, IEEE TCAPT | 2008 | Peer-reviewed | Author PDF read | Yes (solve upper bound) |
| Padmanabhan, SAE 2025-01-5073 | 2025 | Peer-reviewed | Author MS read | Yes (solve lower range) |
| Siemens Mynaric / GOLDTek / Flotrend cases | 2020s | Vendor cases | Pages fetched | Caveated context |
| Hughes et al. CMAME 194; Bazilevs et al. CMAME 199; Hughes & Evans ICES 10-18 | 2005–2010 | Peer-reviewed | PDFs read | Genealogy only |
| NASA CFD Vision 2030 | 2014 | NASA CR | NTRS PDF read | Qualitative context |

## 6. Exclusions (not load-bearing)

- Pointwise/Cadence "75% of CFD analysis in meshing" — originals unreachable
  (301/403); UNVERIFIED wording; folklore layer.
- Ansys Icepak "3–5× quicker setup, 80% fewer clicks" — search-snippet only;
  ansys.com fetch timed out; vendor self-comparison.
- Siemens FLOEFD "65–75% simulation-time reduction" / Cadence Celsius "10×
  faster than legacy" / Flotherm BCI-ROM "40,000× faster" — vendor marketing,
  no methodology; context only.
- Aberdeen-via-Mentor "thermal verification −33%; re-spins −500%" — secondhand,
  original unreachable, arithmetically garbled as transcribed.
- NAFEMS "How Analysis Engineers Spend their Time" (Benchmark, Jan 2020) —
  exists, paywalled; percentages unverified. Candidate for acquisition at the
  next review.
- NASA Glenn "70% of total analysis time" (DRAGON-grid lineage) — primary not
  located.
- 6SigmaET "Thermal Focus: IT Report" whitepaper — URL now redirects to a 403
  page; unrecoverable.

## 7. Review cadence and falsifier hooks

- **Review cadence:** quarterly, alongside the ratification record's review
  date (bead .1.4). A review either re-affirms the envelope, narrows it with
  new admissible sources (NAFEMS 2020 acquisition; an executed incumbent run),
  or widens it with contrary evidence.
- **Upgrade path:** an actually executed and timed incumbent-workflow run on
  the reference cooling task supersedes every survey row at the moment it is
  recorded (protocol source class (a)); until then the record stays
  published-source-derived and ESTIMATED.
- **Falsifier hook for the kill criterion:** if a future measured FrankenSim
  iteration time cannot clear 3× against the LOW bound, the claim reverts to
  "indeterminate" and marketing may not state the reduction factor as met.
