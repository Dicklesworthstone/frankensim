//! IRP M5.5-W2 (pexfp.16): the plan §1.4 simplification-equivalence
//! comparison on the M5.5-W1 seam (`fs_ladder::axis`). The comparator is
//! qualified on closed-form cylinder pairs, then run on a tone hole with and
//! without an undercut through fs-duct's 1-D profiled chimney. Gauntlet G2
//! (closed form) and G3 (metamorphic reduction). Each comparison prints one
//! JSON line.
//!
//! Synthetic declared fixtures (§1.3), no specimen claimed: dry air at
//! 293.15 K; a clarinet-like 0.6 m bore of radius 7.5 mm, 0.42 m above and
//! 0.18 m below one side hole of radius 3.5 mm, flanged open end. The
//! equivalence tolerances, 0.1 dB on raw complex outputs and 1 cent on
//! resonance frequencies, are declared here before any evaluation.
//!
//! The band runs from 50 Hz to the bore's first non-plane mode,
//! `1.8412 c / (2 pi a)` (about 13.4 kHz), the validity limit of the
//! plane-wave chain; a single hole has no lattice cutoff, and the hole's own
//! `k b = 1` point (about 15.6 kHz at b = 3.5 mm) lies above it, so the bore
//! limit binds. Its 2177 frequencies
//! are uniform in `sqrt(f)`, so the spacing follows the resonance width
//! (which grows like `sqrt(f)` under boundary-layer loss): about 1.2 Hz at
//! the first resonance (143 Hz, half-power width about 5 Hz) and 12 Hz at the
//! top (width about 50 Hz).

use std::sync::Arc;
use std::time::Instant;

use fs_blake3::{ContentHash, hash_domain};
use fs_duct::{
    Duct, HoleState, LossModel, PROFILED_CHIMNEY_LIMITATION, RadiusProfile, RiccatiSteps, Segment,
    Termination, input_impedance, input_impedance_load, profiled_tone_hole_junction,
    profiled_tone_hole_shunt, segment_wave, tone_hole_shunt,
};
use fs_ladder::axis::{
    AnalysisPolicy, AxisVerdict, BoundSource, BoundedSolution, Coordinate, Coverage, Dimension,
    EnsembleDecl, EquivalenceDomain, EquivalenceReport, EquivalenceVerdict, ErrorModel,
    LadderReport, Metric, Observable, RefinableAxis, RungRecord, SensitivityFixture,
    compare_equivalence, evaluate_axis, qoi_bound,
};
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_math::det;

/// Declared before evaluation: 0.1 dB on |Z_in| and on |p_out/U_in|.
const TOLERANCE_DB: f64 = 0.1;
/// Declared before evaluation: 1 cent on every resonance frequency.
const TOLERANCE_CENTS: f64 = 1.0;
const Z_UNIT: &str = "Pa s/m^3";
const P_UNIT: &str = "Pa/(m^3/s)";
const F_UNIT: &str = "Hz";
const STEPS_AXIS: &str = "riccati_steps";
/// Riccati step rungs of every chimney profile (fourth-order Magnus). The
/// closed hole's low-frequency outputs reach roundoff by 32 steps, so the
/// ladder stops at 16.
const STEP_RUNGS: [usize; 4] = [2, 4, 8, 16];
const BORE: f64 = 0.0075;
const UPSTREAM: f64 = 0.42;
const DOWNSTREAM: f64 = 0.18;
const HOLE: f64 = 0.0035;
/// First zero of `J1'`: the first non-plane mode's `k a`.
const J1_PRIME_ZERO: f64 = 1.841_183_781_340_659;
const BAND_LOW_HZ: f64 = 50.0;
/// Frequency grid: 2177 points uniform in `sqrt(f)`, sampled on nested
/// rungs of 273, 545, 1089 and 2177 points.
const GRID_POINTS: usize = 2177;
const GRID_STRIDES: [usize; 4] = [8, 4, 2, 1];
const STATES: [HoleState; 2] = [HoleState::Open, HoleState::Closed];

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
}
fn tau(f: f64) -> f64 {
    core::f64::consts::TAU * f
}
/// The bore's plane-wave cutoff, the top of the band [Hz].
fn band_high_hz() -> f64 {
    J1_PRIME_ZERO * air().sound_speed / (core::f64::consts::TAU * BORE)
}
fn freqs() -> Vec<f64> {
    let (u0, u1) = (det::sqrt(BAND_LOW_HZ), det::sqrt(band_high_hz()));
    (0..GRID_POINTS)
        .map(|j| {
            let u = u0 + (u1 - u0) * j as f64 / (GRID_POINTS - 1) as f64;
            u * u
        })
        .collect()
}

/// Location of the maximum of a unimodal `eval` on `[a, b]` by golden-section
/// search, to the floating-point resolution of the bracket.
fn golden_max(eval: &dyn Fn(f64) -> f64, mut a: f64, mut b: f64) -> f64 {
    let g = 0.5 * (det::sqrt(5.0) - 1.0);
    let (mut c, mut d) = (b - g * (b - a), a + g * (b - a));
    let (mut fc, mut fd) = (eval(c), eval(d));
    for _ in 0..80 {
        if fc > fd {
            (b, d, fd) = (d, c, fc);
            c = b - g * (b - a);
            fc = eval(c);
        } else {
            (a, c, fc) = (c, d, fd);
            d = a + g * (b - a);
            fd = eval(d);
        }
    }
    0.5 * (a + b)
}

/// A frozen chimney description: an optional undercut cone (bore-side
/// radius, depth) under a cylinder of the hole radius, total height `height`.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Chimney {
    undercut: Option<(f64, f64)>,
    height: f64,
}

impl Chimney {
    fn identity(&self) -> ContentHash {
        let text =
            format!("bore {BORE} up {UPSTREAM} down {DOWNSTREAM} hole {HOLE} flanged; {self:?}");
        hash_domain("org.frankensim.test.tone-hole-description", text.as_bytes())
    }
    /// Sections from the bore outward. Every profile section carries the
    /// Riccati step axis; a plain chimney is a constant-radius profile so
    /// both sides of a comparison are resolved on the same axis.
    fn sections(&self, steps: usize) -> Vec<Segment> {
        let profile = |p: RadiusProfile| Segment::Profile {
            profile: Arc::new(p),
            steps: RiccatiSteps(steps),
        };
        match self.undercut {
            None => vec![profile(
                RadiusProfile::power(HOLE, 1.0, 0.0, self.height).unwrap(),
            )],
            Some((radius, depth)) => {
                let slope = (HOLE - radius) / depth;
                let cone =
                    RadiusProfile::spline(vec![(0.0, radius), (depth, HOLE)], Some([slope, slope]))
                        .unwrap();
                vec![
                    profile(cone),
                    Segment::Cylinder {
                        radius: HOLE,
                        length: self.height - depth,
                    },
                ]
            }
        }
    }
}

/// Raw outputs at one frequency and hole state: the input impedance and the
/// pressure at the open end's plane per unit inlet volume flow. The hole's
/// 2-port joins the two bore cylinders: `Z_a = (J0 Z_b + J1)/(J2 Z_b + J3)`,
/// `U_b/U_a = 1/(J2 Z_b + J3)`.
fn outputs(chimney: &Chimney, steps: usize, f: f64, hole: HoleState) -> (C64, C64) {
    let (gas, omega, loss) = (air(), tau(f), LossModel::Bessel);
    let cylinder = |length| Duct {
        segments: vec![Segment::Cylinder {
            radius: BORE,
            length,
        }],
    };
    let below = input_impedance(
        &cylinder(DOWNSTREAM),
        &gas,
        omega,
        loss,
        Termination::FlangedOpen,
    )
    .unwrap();
    let j = profiled_tone_hole_junction(
        &gas,
        &chimney.sections(steps),
        hole,
        omega,
        loss,
        BORE,
        None,
    )
    .unwrap();
    let denominator = j[2] * below.impedance + j[3];
    let z_hole_plane = (j[0] * below.impedance + j[1]) * denominator.recip();
    let above =
        input_impedance_load(&cylinder(UPSTREAM), &gas, omega, loss, z_hole_plane, None).unwrap();
    let p_out = below.p_mouth_over_u_in * denominator.recip() * above.u_mouth_over_u_in;
    (above.impedance, p_out)
}

fn policy() -> AnalysisPolicy {
    AnalysisPolicy {
        maximum_order: 16.0,
        order_tolerance: 0.2,
        relative_fit_tolerance: 0.01,
        asymptotic_tolerance: 0.1,
        minimum_refinement_ratio: 1.5,
    }
}

fn observable(name: String, value: f64, unit: &str) -> Observable {
    Observable {
        name,
        value,
        unit: unit.into(),
    }
}

/// Axis liveness on a separate fixture: the 1.2 b undercut chimney's open
/// shunt at 2 kHz moves between 1 and 64 Riccati steps.
fn sensitivity() -> SensitivityFixture {
    let chimney = Chimney {
        undercut: Some((1.2 * HOLE, 0.0015)),
        height: 0.004,
    };
    let shunt = |steps: usize| {
        profiled_tone_hole_shunt(
            &air(),
            &chimney.sections(steps),
            HoleState::Open,
            tau(2000.0),
            LossModel::Bessel,
            BORE,
            None,
        )
        .unwrap()
        .re
    };
    let (coarse, fine) = (shunt(1), shunt(64));
    SensitivityFixture {
        fixture: "1.2 b undercut chimney, open shunt at 2 kHz, Riccati steps 1 vs 64".into(),
        realized: [1.0, 64.0],
        qoi: [
            observable("shunt.re".into(), coarse, Z_UNIT),
            observable("shunt.re".into(), fine, Z_UNIT),
        ],
        minimum_change: 1e-9 * fine.abs(),
    }
}

/// The M5.5-W1 report of one description: Re/Im of both outputs at every
/// (state, frequency) point, on every step rung.
fn ladder(chimney: &Chimney) -> LadderReport {
    let f = freqs();
    let rungs: Vec<RungRecord> = STEP_RUNGS
        .iter()
        .map(|&steps| {
            let start = Instant::now();
            let mut qoi = Vec::with_capacity(4 * STATES.len() * f.len());
            for (s, &hole) in STATES.iter().enumerate() {
                for (j, &hz) in f.iter().enumerate() {
                    let i = s * f.len() + j;
                    let (z, p) = outputs(chimney, steps, hz, hole);
                    qoi.push(observable(format!("zin.re[{i}]"), z.re, Z_UNIT));
                    qoi.push(observable(format!("zin.im[{i}]"), z.im, Z_UNIT));
                    qoi.push(observable(format!("pout.re[{i}]"), p.re, P_UNIT));
                    qoi.push(observable(format!("pout.im[{i}]"), p.im, P_UNIT));
                }
            }
            RungRecord {
                requested: steps as f64,
                dimensions: vec![Dimension {
                    name: STEPS_AXIS.into(),
                    value: steps as f64,
                    unit: "1".into(),
                }],
                qoi,
                bounds: vec![],
                regime: "plane-wave chain, hole open and closed".into(),
                solver_residual: 0.0,
                wall_s: start.elapsed().as_secs_f64(),
            }
        })
        .collect();
    finish("chimney-ladder", chimney, rungs)
}

/// Every |Z_in| maximum of each hole state: bracketed by the finest grid
/// rung of `grid` (a sample above its left and not below its right
/// neighbour), then refined by golden-section search on the continuous
/// frequency at every Riccati-step rung. Returns the M5.5-W1 report and the
/// observable names, `fres[state.k]` in order.
fn resonance_ladder(chimney: &Chimney, grid: &LadderReport) -> (LadderReport, Vec<String>) {
    let f = freqs();
    let n = f.len();
    let finest = grid.joint.rungs.last().expect("grid ladder rungs");
    let mut brackets = Vec::new();
    for s in 0..STATES.len() {
        let magnitude = |j: usize| {
            let i = 4 * (s * n + j);
            C64::new(finest.qoi[i].value, finest.qoi[i + 1].value).abs()
        };
        let mut k = 0;
        for j in 1..n - 1 {
            if magnitude(j) > magnitude(j - 1) && magnitude(j) >= magnitude(j + 1) {
                brackets.push((format!("fres[{s}.{k}]"), s, f[j - 1], f[j + 1]));
                k += 1;
            }
        }
    }
    let rungs: Vec<RungRecord> = STEP_RUNGS
        .iter()
        .map(|&steps| {
            let start = Instant::now();
            let qoi = brackets
                .iter()
                .map(|(name, s, lo, hi)| {
                    let hole = STATES[*s];
                    let peak =
                        golden_max(&|hz| outputs(chimney, steps, hz, hole).0.abs(), *lo, *hi);
                    observable(name.clone(), peak, F_UNIT)
                })
                .collect();
            RungRecord {
                requested: steps as f64,
                dimensions: vec![Dimension {
                    name: STEPS_AXIS.into(),
                    value: steps as f64,
                    unit: "1".into(),
                }],
                qoi,
                bounds: vec![],
                regime: "plane-wave chain, |Z_in| maxima, hole open and closed".into(),
                solver_residual: 0.0,
                wall_s: start.elapsed().as_secs_f64(),
            }
        })
        .collect();
    let names = brackets.into_iter().map(|b| b.0).collect();
    (finish("resonance-ladder", chimney, rungs), names)
}

/// The Riccati-step M5.5-W1 report of one set of rungs, with a receipt line
/// counting the per-observable verdicts.
fn finish(case: &str, chimney: &Chimney, rungs: Vec<RungRecord>) -> LadderReport {
    let axis = RefinableAxis {
        name: STEPS_AXIS.into(),
        unit: "1".into(),
        model: ErrorModel::Algebraic { formal_order: 4.0 },
        group: None,
        coordinate: Coordinate::Reciprocal,
        bound_targets: vec![],
    };
    let report = evaluate_axis(&axis, &rungs, &sensitivity(), &policy());
    let mut counts = std::collections::BTreeMap::new();
    for q in &report.qois {
        let kind = match &q.verdict {
            AxisVerdict::Asymptotic { gci: Some(_), .. } => "asymptotic",
            AxisVerdict::Asymptotic { gci: None, .. } => "asymptotic-no-gci",
            AxisVerdict::PreAsymptotic { .. } => "pre-asymptotic",
            AxisVerdict::Scattered => "scattered",
            AxisVerdict::Oscillatory => "oscillatory",
            AxisVerdict::Bounded { .. } => "bounded",
            AxisVerdict::Unresolved => "unresolved",
            AxisVerdict::ObservedOnly { .. } => "observed-only",
            AxisVerdict::InsufficientRungs => "insufficient",
            AxisVerdict::Refused { .. } => "refused",
        };
        *counts.entry(kind).or_insert(0_usize) += 1;
    }
    let unbounded: Vec<String> = report
        .qois
        .iter()
        .enumerate()
        .filter(|(_, q)| {
            !matches!(
                q.verdict,
                AxisVerdict::Asymptotic { gci: Some(_), .. } | AxisVerdict::Unresolved
            )
        })
        .take(4)
        .map(|(j, q)| {
            let values: Vec<f64> = rungs.iter().map(|r| r.qoi[j].value).collect();
            format!("{} {:?} {values:?}", q.name, q.verdict)
        })
        .collect();
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"{case}\",\"inputs\":\"{chimney:?}\",\
         \"axis\":\"{STEPS_AXIS}\",\"rungs\":{STEP_RUNGS:?},\"qois\":{},\"verdicts\":{counts:?},\
         \"unbounded_examples\":{unbounded:?},\"wall_s\":{:?}}}",
        report.qois.len(),
        rungs.iter().map(|r| r.wall_s).collect::<Vec<_>>()
    );
    LadderReport::new(
        policy(),
        vec![report],
        vec![],
        rungs,
        EnsembleDecl::Deterministic {
            reason: "passive frequency-domain chain".into(),
        },
    )
    .expect("complete M5.5-W1 ladder")
}

fn domain(limitations: Vec<String>) -> EquivalenceDomain {
    let n = GRID_POINTS;
    EquivalenceDomain {
        description: format!(
            "{BAND_LOW_HZ}-{:.1} Hz (bore plane-wave cutoff), {n} points uniform in sqrt(f); \
             holes open and closed; receivers: inlet impedance and open-end plane pressure per \
             inlet volume flow",
            band_high_hz()
        ),
        coverage: Coverage::Sampled {
            rungs: GRID_STRIDES
                .iter()
                .map(|&stride| {
                    (0..STATES.len())
                        .flat_map(|s| (0..n).step_by(stride).map(move |j| s * n + j))
                        .collect()
                })
                .collect(),
            lines: (0..STATES.len())
                .map(|s| (s * n..(s + 1) * n).collect())
                .collect(),
            coordinate: STATES.iter().flat_map(|_| freqs()).collect(),
        },
        ensemble: EnsembleDecl::Deterministic {
            reason: "passive frequency-domain chain".into(),
        },
        limitations,
    }
}

/// A resonance set's values and bounds, from its own report.
fn resonance_side(chimney: &Chimney, report: &LadderReport, names: &[String]) -> BoundedSolution {
    let finest = report.joint.rungs.last().expect("resonance rungs");
    let mut side = BoundedSolution {
        description: chimney.identity(),
        ladder: report.hash,
        unit: F_UNIT.into(),
        values: vec![],
        bounds: vec![],
        sources: vec![],
    };
    for name in names {
        let q = finest
            .qoi
            .iter()
            .find(|q| &q.name == name)
            .expect("resonance");
        let (bound, source) =
            qoi_bound(report, name).unwrap_or_else(|| panic!("{name}: no ladder bound"));
        side.values.push((q.value, 0.0));
        side.bounds.push(bound);
        side.sources.push(source);
    }
    side
}

fn receipt(case: &str, kind: &str, report: &EquivalenceReport, start: Instant) {
    let (max_lower, max_upper) = (
        report.lower.iter().copied().fold(0.0, f64::max),
        report.sampled_maxima.last().unwrap() + report.sampling_term,
    );
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"{case}\",\"receiver\":\"{kind}\",\
         \"detailed_hash\":\"{}\",\"simplified_hash\":\"{}\",\"domain_points\":{},\
         \"metric\":\"{:?}\",\"tolerance\":{},\"max_upper\":{max_upper:.6e},\
         \"max_lower\":{max_lower:.6e},\"sampled_maxima\":{:?},\"sampling_terms\":{:?},\
         \"sampling_term\":{:.3e},\"spread_bounded_points\":{:?},\
         \"verdict\":\"{}\",\"limitations\":{:?},\"ladder_report_hashes\":[\"{}\",\"{}\"],\
         \"elapsed_s\":{:.3}}}",
        report.descriptions[0].to_hex(),
        report.descriptions[1].to_hex(),
        report.upper.len(),
        report.metric,
        report.tolerance,
        report.sampled_maxima,
        report.sampling_terms,
        report.sampling_term,
        report.spread_bounded,
        verdict_name(&report.verdict),
        report.limitations,
        report.ladders[0].to_hex(),
        report.ladders[1].to_hex(),
        start.elapsed().as_secs_f64()
    );
}

/// A detailed vs simplified pair, each resolved on its own ladders, with the
/// 1-D chimney limitation carried: input impedance and receiver transfer in
/// dB over the sampled band, then every resonance frequency in cents.
fn compare(detailed: &Chimney, simplified: &Chimney, case: &str) -> [EquivalenceReport; 3] {
    let start = Instant::now();
    let reports = [ladder(detailed), ladder(simplified)];
    let points = |kind: &str| -> Vec<(String, String)> {
        (0..STATES.len() * GRID_POINTS)
            .map(|i| (format!("{kind}.re[{i}]"), format!("{kind}.im[{i}]")))
            .collect()
    };
    let limitation = vec![PROFILED_CHIMNEY_LIMITATION.to_string()];
    let [zin, pout] = [("zin", Z_UNIT), ("pout", P_UNIT)].map(|(kind, unit)| {
        let side = |c: &Chimney, r: &LadderReport| {
            BoundedSolution::from_ladder(c.identity(), r, unit, &points(kind))
                .unwrap_or_else(|e| panic!("{case} {kind}: {e}"))
        };
        let report = compare_equivalence(
            &side(detailed, &reports[0]),
            &side(simplified, &reports[1]),
            &domain(limitation.clone()),
            Metric::LevelDb,
            TOLERANCE_DB,
        )
        .unwrap();
        receipt(case, kind, &report, start);
        report
    });
    let (detailed_peaks, names) = resonance_ladder(detailed, &reports[0]);
    let (simplified_peaks, simplified_names) = resonance_ladder(simplified, &reports[1]);
    assert_eq!(
        names, simplified_names,
        "{case}: the two descriptions have different resonance sets in the band"
    );
    let resonances = compare_equivalence(
        &resonance_side(detailed, &detailed_peaks, &names),
        &resonance_side(simplified, &simplified_peaks, &names),
        &EquivalenceDomain {
            description: format!(
                "every |Z_in| maximum in {BAND_LOW_HZ}-{:.1} Hz, holes open and closed, found on \
                 the {GRID_POINTS}-point grid and refined: {} resonances",
                band_high_hz(),
                names.len()
            ),
            coverage: Coverage::Enumerated,
            ensemble: EnsembleDecl::Deterministic {
                reason: "passive frequency-domain chain".into(),
            },
            limitations: limitation,
        },
        Metric::Cents,
        TOLERANCE_CENTS,
    )
    .unwrap();
    receipt(case, "resonance-frequency", &resonances, start);
    [zin, pout, resonances]
}

fn verdict_name(v: &EquivalenceVerdict) -> &'static str {
    match v {
        EquivalenceVerdict::Equivalent { .. } => "Equivalent",
        EquivalenceVerdict::NotEquivalent { .. } => "NotEquivalent",
        EquivalenceVerdict::Inconclusive { .. } => "Inconclusive",
    }
}

const PLAIN: Chimney = Chimney {
    undercut: None,
    height: 0.004,
};
const UNDERCUT: Chimney = Chimney {
    undercut: Some((1.2 * HOLE, 0.0015)),
    height: 0.004,
};

#[test]
fn profiled_chimney_reduces_to_cylindrical_tone_hole() {
    // G3: one cylindrical section is the cylindrical tone hole, bit for bit
    // in the shunt; a constant-radius profile matches it to roundoff; and
    // the hole's 2-port between two cylinders matches fs-duct's own chain
    // with a `Segment::ToneHole`.
    let gas = air();
    let (mut shunt_err, mut profile_err, mut chain_err) = (0.0_f64, 0.0_f64, 0.0_f64);
    let cylinder = [Segment::Cylinder {
        radius: HOLE,
        length: 0.004,
    }];
    for hz in [50.0, 700.0, 2500.0, 4000.0] {
        for hole in [HoleState::Open, HoleState::Closed, HoleState::Vent(0.4)] {
            for loss in [LossModel::Lossless, LossModel::WideTube, LossModel::Bessel] {
                let reference =
                    tone_hole_shunt(&gas, HOLE, 0.004, hole, tau(hz), loss, BORE).unwrap();
                let profiled =
                    profiled_tone_hole_shunt(&gas, &cylinder, hole, tau(hz), loss, BORE, None)
                        .unwrap();
                shunt_err = shunt_err.max((profiled - reference).abs() / reference.abs());
                let smooth = profiled_tone_hole_shunt(
                    &gas,
                    &PLAIN.sections(64),
                    hole,
                    tau(hz),
                    loss,
                    BORE,
                    None,
                )
                .unwrap();
                profile_err = profile_err.max((smooth - reference).abs() / reference.abs());
            }
            let chain = input_impedance(
                &Duct {
                    segments: vec![
                        Segment::Cylinder {
                            radius: BORE,
                            length: UPSTREAM,
                        },
                        Segment::ToneHole {
                            hole_radius: HOLE,
                            chimney_height: 0.004,
                            bore_radius: BORE,
                            state: hole,
                        },
                        Segment::Cylinder {
                            radius: BORE,
                            length: DOWNSTREAM,
                        },
                    ],
                },
                &gas,
                tau(hz),
                LossModel::Bessel,
                Termination::FlangedOpen,
            )
            .unwrap();
            // `outputs` uses the constant-radius profile; compare the chain
            // against the exact cylindrical section through the same 2-port.
            let j = profiled_tone_hole_junction(
                &gas,
                &cylinder,
                hole,
                tau(hz),
                LossModel::Bessel,
                BORE,
                None,
            )
            .unwrap();
            let below = input_impedance(
                &Duct {
                    segments: vec![Segment::Cylinder {
                        radius: BORE,
                        length: DOWNSTREAM,
                    }],
                },
                &gas,
                tau(hz),
                LossModel::Bessel,
                Termination::FlangedOpen,
            )
            .unwrap();
            let z_a = (j[0] * below.impedance + j[1]) * (j[2] * below.impedance + j[3]).recip();
            let above = input_impedance_load(
                &Duct {
                    segments: vec![Segment::Cylinder {
                        radius: BORE,
                        length: UPSTREAM,
                    }],
                },
                &gas,
                tau(hz),
                LossModel::Bessel,
                z_a,
                None,
            )
            .unwrap();
            chain_err = chain_err
                .max((above.impedance - chain.impedance).abs() / chain.impedance.abs())
                .max(
                    (outputs(&PLAIN, 64, hz, hole).0 - chain.impedance).abs()
                        / chain.impedance.abs(),
                );
        }
    }
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"profiled-chimney-reduction\",\"observed\":{{\
         \"shunt_rel\":{shunt_err:e},\"constant_profile_rel\":{profile_err:e},\
         \"chain_rel\":{chain_err:e}}},\"tol\":1e-12}}"
    );
    assert!(
        shunt_err == 0.0,
        "one cylindrical section vs tone_hole_shunt {shunt_err:e}"
    );
    assert!(
        profile_err < 1e-12,
        "constant-radius profile {profile_err:e}"
    );
    assert!(
        chain_err < 1e-12,
        "2-port vs Segment::ToneHole chain {chain_err:e}"
    );
    // Refusals: no sections, a tone-hole section, a bore-side radius at the bore.
    let call = |sections: &[Segment]| {
        profiled_tone_hole_shunt(
            &gas,
            sections,
            HoleState::Open,
            tau(500.0),
            LossModel::Bessel,
            BORE,
            None,
        )
    };
    assert!(call(&[]).is_err());
    assert!(
        call(&[Segment::ToneHole {
            hole_radius: HOLE,
            chimney_height: 0.004,
            bore_radius: BORE,
            state: HoleState::Open,
        }])
        .is_err()
    );
    assert!(
        call(&[Segment::Cylinder {
            radius: BORE,
            length: 0.004,
        }])
        .is_err()
    );
    assert!(
        call(&[Segment::Cylinder {
            radius: HOLE,
            length: -0.004,
        }])
        .is_err()
    );
}

/// Closed-form input impedance of a closed cylinder, `i Zc cot(kL)` under
/// `e^{-i omega t}`, with the exact Bessel (or lossless) `k` and `Zc`.
fn closed_cylinder(length: f64, f: f64, loss: LossModel) -> C64 {
    let wave = segment_wave(&air(), BORE, tau(f), loss).unwrap();
    let kl = wave.wavenumber.scale(length);
    let exp = |z: C64| {
        let (s, c) = det::sin_cos(z.im);
        C64::new(c, s).scale(det::exp(z.re))
    };
    let (e_pos, e_neg) = (exp(C64::new(-kl.im, kl.re)), exp(C64::new(kl.im, -kl.re)));
    let i = C64::new(0.0, 1.0);
    i * wave.characteristic_impedance * (i * (e_pos + e_neg) * (e_pos - e_neg).recip())
}

/// Two closed cylinders differing in length by 1 mm, solved by fs-duct's
/// exact transfer matrix (no discretization axis): the declared bound is a
/// roundoff allowance of 1e-12 |Z|. `n` frequencies over `band`.
fn length_pair(
    band: (f64, f64),
    n: usize,
    loss: LossModel,
) -> (BoundedSolution, BoundedSolution, Vec<f64>) {
    let f: Vec<f64> = (0..n)
        .map(|j| band.0 + (band.1 - band.0) * j as f64 / (n - 1) as f64)
        .collect();
    let solve = |length: f64, tag: &str| {
        let values: Vec<(f64, f64)> = f
            .iter()
            .map(|&hz| {
                let z = input_impedance(
                    &Duct {
                        segments: vec![Segment::Cylinder {
                            radius: BORE,
                            length,
                        }],
                    },
                    &air(),
                    tau(hz),
                    loss,
                    Termination::Closed,
                )
                .unwrap()
                .impedance;
                (z.re, z.im)
            })
            .collect();
        BoundedSolution {
            description: hash_domain("org.frankensim.test.cylinder", tag.as_bytes()),
            ladder: hash_domain(
                "org.frankensim.test.cylinder",
                b"closed-form transfer matrix: no discretization axis; 1e-12 |Z| roundoff allowance",
            ),
            unit: Z_UNIT.into(),
            bounds: values
                .iter()
                .map(|&(re, im)| 1e-12 * det::sqrt(re * re + im * im))
                .collect(),
            sources: vec![BoundSource::Declared; values.len()],
            values,
        }
    };
    (
        solve(0.300, "closed 0.300 m"),
        solve(0.301, "closed 0.301 m"),
        f,
    )
}

/// The last `count` nested rungs over `n = 2^m + 1` points (strides
/// `2^(m-k)`), on one line with the point index as its coordinate.
fn plain_domain(n: usize, count: usize) -> EquivalenceDomain {
    let mut strides = vec![];
    let mut stride = n - 1;
    while stride >= 1 {
        strides.push(stride);
        stride /= 2;
    }
    let kept = &strides[strides.len() - count..];
    EquivalenceDomain {
        description: format!("closed cylinders, {n} frequencies"),
        coverage: Coverage::Sampled {
            rungs: kept.iter().map(|&s| (0..n).step_by(s).collect()).collect(),
            lines: vec![(0..n).collect()],
            coordinate: (0..n).map(|j| j as f64).collect(),
        },
        ensemble: EnsembleDecl::Deterministic {
            reason: "passive frequency-domain chain".into(),
        },
        limitations: vec![],
    }
}

#[test]
fn analytic_length_pair_bound_contains_truth() {
    // G2: the exact |dZ| of a 1 mm length change lies inside every point's
    // [lower, upper], lossless (below the first cot pole) and Bessel-lossy
    // (across the first resonance); and the reported band maximum
    // U_K + T is not below the true supremum (65 x denser exact grid).
    for (loss, band) in [
        (LossModel::Lossless, (50.0, 500.0)),
        (LossModel::Bessel, (400.0, 750.0)),
    ] {
        let (long, short, f) = length_pair(band, 257, loss);
        let report =
            compare_equivalence(&short, &long, &plain_domain(257, 4), Metric::Absolute, 1e9)
                .unwrap();
        let truth =
            |hz: f64| (closed_cylinder(0.301, hz, loss) - closed_cylinder(0.300, hz, loss)).abs();
        let mut inside = 0;
        for (i, &hz) in f.iter().enumerate() {
            let t = truth(hz);
            assert!(
                report.lower[i] <= t && t <= report.upper[i],
                "{loss:?} {hz} Hz: |dZ| {t:e} outside [{:e}, {:e}]",
                report.lower[i],
                report.upper[i]
            );
            inside += 1;
        }
        let dense = 16 * 256 + 1;
        let sup = (0..dense)
            .map(|j| truth(band.0 + (band.1 - band.0) * j as f64 / (dense - 1) as f64))
            .fold(0.0, f64::max);
        let reported = report.sampled_maxima.last().unwrap() + report.sampling_term;
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"analytic-length-pair\",\"inputs\":{{\"loss\":\"{loss:?}\",\
             \"band_hz\":[{},{}]}},\"observed\":{{\"points_inside\":{inside},\"sampled_maxima\":{:?},\
             \"sampling_term\":{:e},\"reported_sup\":{reported:e},\"dense_sup\":{sup:e}}}}}",
            band.0, band.1, report.sampled_maxima, report.sampling_term
        );
        assert!(
            report.sampling_term.is_finite(),
            "{loss:?}: the sampled maximum settles"
        );
        assert!(
            reported >= sup,
            "{loss:?}: reported {reported:e} < true sup {sup:e}"
        );
    }
}

#[test]
fn sampling_term_converges() {
    // G1: the lossy pair across its resonance (half-power width about 10 Hz)
    // on 257 frequencies, sampled on five nested rungs of 17 to 257 points.
    // The per-rung curvature terms are printed with the gap to the dense
    // supremum; on the finest rungs the term shrinks, stays finite and
    // covers the gap.
    let band = (400.0, 750.0);
    let (long, short, _) = length_pair(band, 257, LossModel::Bessel);
    let report =
        compare_equivalence(&short, &long, &plain_domain(257, 5), Metric::Absolute, 1e9).unwrap();
    let truth = |hz: f64| {
        (closed_cylinder(0.301, hz, LossModel::Bessel)
            - closed_cylinder(0.300, hz, LossModel::Bessel))
        .abs()
    };
    let dense = 16 * 256 + 1;
    let sup = (0..dense)
        .map(|j| truth(band.0 + (band.1 - band.0) * j as f64 / (dense - 1) as f64))
        .fold(0.0, f64::max);
    for (k, (u, t)) in report
        .sampled_maxima
        .iter()
        .zip(&report.sampling_terms)
        .enumerate()
    {
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"sampling-term\",\"axis\":\"grid\",\"rung\":{},\
             \"observed\":{{\"sampled_max\":{u:e},\"sampling_term\":{t:e},\"gap_to_sup\":{:e}}}}}",
            [17, 33, 65, 129, 257][k],
            sup - u
        );
    }
    let terms = &report.sampling_terms;
    let last = terms.len() - 1;
    assert!(report.sampling_term.is_finite(), "{terms:?}");
    assert!(terms[last] < terms[last - 1], "the term shrinks: {terms:?}");
    assert!(
        terms[last] > 0.0,
        "liveness: no grid samples the peak exactly"
    );
    for k in last - 1..=last {
        assert!(
            report.sampled_maxima[k] + terms[k] >= sup,
            "rung {k}: term {:e} misses gap {:e}",
            terms[k],
            sup - report.sampled_maxima[k]
        );
    }
}

#[test]
fn straddling_bound_is_inconclusive() {
    // An injected numerical uncertainty of 20 % of the tolerance on both
    // sides, with the tolerance set at the sampled maximum, straddles it.
    let (long, short, _) = length_pair((50.0, 500.0), 65, LossModel::Lossless);
    let exact =
        compare_equivalence(&short, &long, &plain_domain(65, 3), Metric::Absolute, 1e9).unwrap();
    let tolerance = *exact.sampled_maxima.last().unwrap();
    let inject = |s: &BoundedSolution| BoundedSolution {
        bounds: s.bounds.iter().map(|b| b + 0.2 * tolerance).collect(),
        ..s.clone()
    };
    let report = compare_equivalence(
        &inject(&short),
        &inject(&long),
        &plain_domain(65, 3),
        Metric::Absolute,
        tolerance,
    )
    .unwrap();
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"straddling\",\"tolerance\":{tolerance:e},\"verdict\":\"{}\"}}",
        verdict_name(&report.verdict)
    );
    assert!(
        matches!(report.verdict, EquivalenceVerdict::Inconclusive { ref straddling } if !straddling.is_empty()),
        "{:?}",
        report.verdict
    );
    assert!(!report.authorizes_omission());
}

#[test]
fn undercut_pair_reports_bound_with_1d_limitation() {
    // The 1.2 b undercut (1.5 mm deep) vs its omission, both resolved on
    // their own Riccati ladders. Whatever the verdict, the report carries
    // the 1-D limitation and never authorizes omitting a real undercut.
    for report in compare(&UNDERCUT, &PLAIN, "undercut-1.2b-vs-plain") {
        assert_eq!(
            report.limitations,
            vec![PROFILED_CHIMNEY_LIMITATION.to_string()]
        );
        assert!(!report.authorizes_omission());
        if !matches!(report.metric, Metric::Cents) {
            assert_eq!(report.upper.len(), STATES.len() * GRID_POINTS);
        }
        assert!(report.sampled_maxima.windows(2).all(|w| w[1] >= w[0]));
        assert!(
            report
                .lower
                .iter()
                .zip(&report.upper)
                .all(|(lo, hi)| 0.0 <= *lo && lo <= hi)
        );
    }
}

#[test]
fn large_simplification_not_equivalent() {
    // Omitting a 1.5 b undercut must be caught.
    let large = Chimney {
        undercut: Some((1.5 * HOLE, 0.0015)),
        height: 0.004,
    };
    for report in compare(&large, &PLAIN, "undercut-1.5b-vs-plain") {
        assert!(
            matches!(report.verdict, EquivalenceVerdict::NotEquivalent { .. }),
            "{:?}",
            report.verdict
        );
    }
}

#[test]
fn combined_simplifications_tested_jointly() {
    // Undercut omission and chimney-height rounding (4.4 mm -> 4 mm), each
    // alone and both together; the joint pair is its own comparison, never
    // inferred from the single ones.
    let detailed = Chimney {
        undercut: Some((1.2 * HOLE, 0.0015)),
        height: 0.0044,
    };
    let no_undercut = Chimney {
        undercut: None,
        height: 0.0044,
    };
    let rounded = Chimney {
        undercut: Some((1.2 * HOLE, 0.0015)),
        height: 0.004,
    };
    let both = Chimney {
        undercut: None,
        height: 0.004,
    };
    let singles = [
        compare(&detailed, &no_undercut, "combined:undercut-omitted"),
        compare(&detailed, &rounded, "combined:height-rounded"),
    ];
    let joint = compare(&detailed, &both, "combined:both");
    for (k, report) in joint.iter().enumerate() {
        assert_eq!(report.descriptions, [detailed.identity(), both.identity()]);
        assert!(singles.iter().all(|s| s[k].hash != report.hash));
        // The joint upper bound is not the sum of the singles; report it.
        let single_sum = singles[0][k]
            .upper
            .iter()
            .zip(&singles[1][k].upper)
            .map(|(a, b)| a + b);
        let interaction = report
            .upper
            .iter()
            .zip(single_sum)
            .map(|(j, s)| (j - s).abs())
            .fold(0.0, f64::max);
        println!(
            "{{\"suite\":\"fs-couple\",\"case\":\"combined-interaction\",\"receiver\":{k},\
             \"max_abs_joint_minus_single_sum_db\":{interaction:e}}}"
        );
        assert!(interaction.is_finite());
    }
}
