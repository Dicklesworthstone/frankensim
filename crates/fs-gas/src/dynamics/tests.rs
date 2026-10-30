use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Length, MolarMass, Pressure, Temperature, Time};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K, required_nasa9_parameter_dims,
};

fn source(
    cp_linear: f64,
    formation: f64,
    minimum: f64,
    maximum: f64,
) -> FrozenIdealGasMixtureModelV1 {
    // Manufactured source: R_specific=1 J/(kg K), Cp/R=3.5+cp_linear*T.
    // This gives a dimensional SI fixture, not real-gas material data.
    let mut a = [0.0; 9];
    a[2] = 3.5;
    a[3] = cp_linear;
    a[7] = formation;
    let mut parameters = std::collections::BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(a) {
        parameters.insert(name.to_owned(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 1.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId("nasa9-standard-state".into()),
        law_version: 1,
        parameters,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", minimum, maximum),
        sources: Vec::new(),
        provenance: Provenance {
            source: "Manufactured Euler caloric test".into(),
            license: "test fixture".into(),
            artifact: None,
        },
    };
    let convention = StandardStateConventionV1::new(
        StandardStatePhaseV1::Gas,
        ReferenceEquationOfStateV1::IdealGas,
        Pressure::new(1.0),
        ElementalReferenceIdV1::new("euler-fixture").unwrap(),
    )
    .unwrap();
    let species = Nasa9StandardStateModelV1::new(
        SpeciesId::new("test-gas").unwrap(),
        MolarMass::new(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K),
        convention,
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap();
    FrozenIdealGasMixtureModelV1::new(
        vec![species],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
    )
    .unwrap()
}
pub(super) fn gas(cp_linear: f64, formation: f64) -> FrozenGas {
    FrozenGas::new(
        source(cp_linear, formation, 0.001, 50.0),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap()
}
fn field(gas: FrozenGas, n: usize, boundaries: Boundaries, f: impl Fn(f64) -> [f64; 3]) -> Euler1d {
    let q = (0..n)
        .map(|i| {
            let p = f((i as f64 + 0.5) / n as f64);
            gas.conserved(p[0], p[1], p[2]).unwrap()
        })
        .collect();
    Euler1d::new(gas, Length::new(1.0 / n as f64), 0.7, q, boundaries, || {
        false
    })
    .unwrap()
}
fn run(flow: &mut Euler1d, end: f64) -> usize {
    let mut steps = 0;
    while flow.time().value() < end {
        let report = flow
            .step(
                Time::new(end - flow.time().value()),
                StepControl::default(),
                || false,
            )
            .unwrap();
        assert!(report.maximum_cfl <= 0.40000000000001);
        assert!(report.conservation_residual.mass.abs() < 3e-13);
        assert!(report.conservation_residual.momentum.abs() < 3e-12);
        assert!(report.conservation_residual.energy.abs() < 3e-11);
        steps += 1;
        assert!(steps < 10_000);
    }
    steps
}
const OPEN: Boundaries = Boundaries::Ends {
    left: Boundary::Transmissive,
    right: Boundary::Transmissive,
};

#[path = "geometry_tests.rs"]
mod geometry_tests;

#[path = "normal_flux_tests.rs"]
mod normal_flux_tests;

#[path = "spatial_flux_tests.rs"]
mod spatial_flux_tests;

#[test]
fn g0_source_energy_gauge_variable_cv_and_supersonic_round_trip() {
    let ordinary = gas(0.4, 0.0);
    for reference in [-1e18, 0.0, 1e18] {
        let shifted = gas(0.4, reference);
        for t in [0.1_f64, 1.0, 3.0, 20.0] {
            for velocity in [-30.0, 0.0, 45.0] {
                let q = shifted.conserved(0.8, velocity, 0.8 * t).unwrap();
                assert_eq!(q, ordinary.conserved(0.8, velocity, 0.8 * t).unwrap());
                let p = shifted.primitive(q, || false).unwrap();
                assert!((p.temperature - t).abs() < 2e-9);
                let gamma = (3.5 + 0.4 * t) / (2.5 + 0.4 * t);
                assert!((p.sound_speed - (gamma * t).sqrt()).abs() < 1e-8);
                let expected = 0.8
                    * (2.5 * (t - 0.001)
                        + 0.2 * (t * t - 0.001 * 0.001)
                        + 0.5 * velocity * velocity);
                assert!((q.energy - expected).abs() < 1e-10);
                let thermo = shifted.thermodynamics(q, || false).unwrap();
                assert_eq!(
                    thermo.receipt().component_receipts()[0].card_identity(),
                    shifted.source().components()[0].regions()[0].card_identity()
                );
            }
        }
    }
}

#[test]
fn g0_hllc_contact_consistency_and_rarefaction_fallback() {
    let g = gas(0.0, 1e14);
    let l = g.conserved(1.0, 0.0, 1.0).unwrap();
    let r = g.conserved(4.0, 0.0, 1.0).unwrap();
    let f = g.hllc_flux(l, r, || false).unwrap();
    assert_eq!(f.mass, 0.0);
    assert_eq!(f.energy, 0.0);
    assert!((f.momentum - 1.0).abs() < 1e-14);
    assert!(!f.used_hll_fallback);
    for u in [-8.0, 0.0, 10.0] {
        let q = g.conserved(1.2, u, 1.7).unwrap();
        let f = g.hllc_flux(q, q, || false).unwrap();
        assert!((f.mass - q.momentum).abs() < 1e-13);
        assert!((f.momentum - (q.momentum * u + 1.7)).abs() < 1e-12);
        assert!((f.energy - (q.energy + 1.7) * u).abs() < 1e-10);
    }
    let f = g
        .hllc_flux(
            g.conserved(1.0, -2.0, 0.4).unwrap(),
            g.conserved(1.0, 2.0, 0.4).unwrap(),
            || false,
        )
        .unwrap();
    assert!(f.used_hll_fallback);
    assert!(f.signal_speed > 2.0);
}

#[test]
fn g0_periodic_uniform_and_closed_wall_conservation() {
    let mut uniform = field(gas(0.4, -1e16), 32, Boundaries::Periodic, |_| {
        [1.1, 2.0, 1.7]
    });
    let before = uniform.states().to_vec();
    run(&mut uniform, 0.1);
    assert_eq!(uniform.states(), before);
    let walls = Boundaries::Ends {
        left: Boundary::Reflecting,
        right: Boundary::Reflecting,
    };
    let mut closed = field(gas(0.4, -4.0), 64, walls, |_| [1.0, 0.6, 1.0]);
    let initial = closed.totals();
    for _ in 0..30 {
        let report = closed
            .step(Time::new(0.002), StepControl::default(), || false)
            .unwrap();
        assert_eq!(report.boundary_change.mass, 0.0);
        assert_eq!(report.boundary_change.energy, 0.0);
    }
    assert!((closed.totals().mass - initial.mass).abs() < 2e-14);
    assert!((closed.totals().energy - initial.energy).abs() < 1e-13);
}

#[test]
fn g1_smooth_contact_muscl_refines_and_conserves() {
    let mut errors = Vec::new();
    for n in [40, 80, 160] {
        let mut flow = field(gas(0.0, 0.0), n, Boundaries::Periodic, |x| {
            [1.0 + 0.2 * (2.0 * std::f64::consts::PI * x).sin(), 1.0, 1.0]
        });
        let initial = flow.totals();
        run(&mut flow, 0.2);
        let error = flow
            .states()
            .iter()
            .enumerate()
            .map(|(i, q)| {
                let x = (i as f64 + 0.5) / n as f64;
                (q.density - (1.0 + 0.2 * (2.0 * std::f64::consts::PI * (x - 0.2)).sin())).abs()
            })
            .sum::<f64>()
            / n as f64;
        errors.push(error);
        assert!((flow.totals().mass - initial.mass).abs() < 5e-14);
        assert!((flow.totals().momentum - initial.momentum).abs() < 5e-14);
        assert!((flow.totals().energy - initial.energy).abs() < 5e-13);
    }
    eprintln!("smooth contact L1 errors: {errors:?}");
    assert!(errors[0] / errors[1] > 2.7 && errors[1] / errors[2] > 2.7);
}

#[test]
fn g2_sod_and_lax_shock_tubes_match_independent_exact_riemann_solution() {
    for (left, right, time, tolerance) in [
        ([1.0, 0.0, 1.0], [0.125, 0.0, 0.1], 0.2, 0.012),
        ([0.445, 0.698, 3.528], [0.5, 0.0, 0.571], 0.12, 0.03),
    ] {
        let mut errors = Vec::new();
        for n in [120, 240] {
            let mut flow = field(
                gas(0.0, 0.0),
                n,
                OPEN,
                |x| if x < 0.5 { left } else { right },
            );
            run(&mut flow, time);
            let exact = ExactRiemann::new(left, right);
            let error = flow
                .states()
                .iter()
                .enumerate()
                .map(|(i, q)| {
                    let xi = ((i as f64 + 0.5) / n as f64 - 0.5) / time;
                    (q.density - exact.density(xi)).abs()
                })
                .sum::<f64>()
                / n as f64;
            errors.push(error);
        }
        eprintln!("shock left={left:?} density L1 errors: {errors:?}");
        assert!(errors[1] < tolerance, "{errors:?}");
        assert!(errors[0] / errors[1] > 1.3, "{errors:?}");
    }
}

#[test]
fn g1_variable_cp_evolution_carries_the_original_acoustic_speed() {
    let cp_linear = 0.4_f64;
    let a = ((3.5 + cp_linear) / (2.5 + cp_linear)).sqrt();
    let epsilon = 1e-5;
    let mut flow = field(gas(cp_linear, -1e15), 128, Boundaries::Periodic, |x| {
        let wave = epsilon * (2.0 * std::f64::consts::PI * x).sin();
        [1.0 + wave, a * wave, 1.0 + a * a * wave]
    });
    run(&mut flow, 0.2);
    let (mut sine, mut cosine) = (0.0, 0.0);
    for (i, q) in flow.states().iter().enumerate() {
        let phase = 2.0 * std::f64::consts::PI * (i as f64 + 0.5) / 128.0;
        sine += (q.density - 1.0) * phase.sin();
        cosine += (q.density - 1.0) * phase.cos();
    }
    let measured = (-cosine).atan2(sine) / (2.0 * std::f64::consts::PI * 0.2);
    eprintln!("variable-Cp acoustic speed: measured={measured}, source={a}");
    assert!((measured / a - 1.0).abs() < 0.002);
    assert!((measured / (1.4_f64.sqrt()) - 1.0).abs() > 0.01);
    // A nonlinear variable-Cp shock evolution also closes every conserved flux.
    let mut shock = field(gas(cp_linear, 0.0), 96, OPEN, |x| {
        if x < 0.5 {
            [1.0, 0.0, 4.0]
        } else {
            [0.2, 0.0, 0.2]
        }
    });
    run(&mut shock, 0.08);
    assert!(shock.states().iter().any(|q| q.momentum > 0.1));
    for &q in shock.states() {
        let p = shock.gas().primitive(q, || false).unwrap();
        assert!(p.temperature > 0.0 && p.pressure > 0.0);
    }
}

#[test]
fn g0_cancellation_invalid_controls_and_shock_source_refusal_are_atomic() {
    let make = || {
        field(gas(0.0, 0.0), 32, OPEN, |x| {
            if x < 0.5 {
                [1.0, 0.0, 1.0]
            } else {
                [0.125, 0.0, 0.1]
            }
        })
    };
    let mut measured = make();
    let mut calls = 0;
    measured
        .step(Time::new(0.001), StepControl::default(), || {
            calls += 1;
            false
        })
        .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut flow = make();
        let before = flow.states().to_vec();
        let mut seen = 0;
        assert_eq!(
            flow.step(Time::new(0.001), StepControl::default(), || {
                seen += 1;
                seen >= stop
            }),
            Err(Error::Cancelled)
        );
        assert_eq!(flow.states(), before);
        assert_eq!(flow.time().value(), 0.0);
    }
    let mut flow = make();
    let before = flow.states().to_vec();
    assert!(
        flow.step(
            Time::new(0.1),
            StepControl {
                cfl: 0.6,
                ..StepControl::default()
            },
            || false
        )
        .is_err()
    );
    assert_eq!(flow.states(), before);
    let narrow = FrozenGas::new(
        source(0.0, 0.0, 0.001, 0.6),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(0.6)),
        || false,
    )
    .unwrap();
    let mut collision = field(narrow, 32, OPEN, |x| {
        [1.0, if x < 0.5 { 5.0 } else { -5.0 }, 0.6 - 1e-12]
    });
    let before = collision.states().to_vec();
    let result = collision.step(
        Time::new(0.1),
        StepControl {
            maximum_halvings: 3,
            ..StepControl::default()
        },
        || false,
    );
    assert!(
        matches!(result, Err(Error::StepRetries { .. })),
        "{result:?}"
    );
    assert_eq!(collision.states(), before);
    assert_eq!(collision.time().value(), 0.0);
}

#[test]
fn g0_rarefaction_remains_positive_and_bad_caloric_path_refuses() {
    let mut flow = field(gas(0.0, 0.0), 128, OPEN, |x| {
        [1.0, if x < 0.5 { -1.5 } else { 1.5 }, 0.4]
    });
    run(&mut flow, 0.06);
    for &q in flow.states() {
        let p = flow.gas().primitive(q, || false).unwrap();
        assert!(p.density > 0.0 && p.pressure > 0.0 && p.temperature >= 0.001);
    }
    assert!(
        FrozenGas::new(
            source(-1.0, 0.0, 0.001, 50.0),
            CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
            || false
        )
        .is_err()
    );
    assert_eq!(
        FrozenGas::new(
            source(0.0, 0.0, 0.001, 50.0),
            CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
            || true
        ),
        Err(Error::Cancelled)
    );
}

// Independent exact gamma-law Riemann oracle from shock Hugoniot and rarefaction
// invariants; no HLLC flux, source inversion or finite-volume code is reused.
struct ExactRiemann {
    left: [f64; 3],
    right: [f64; 3],
    pressure: f64,
    velocity: f64,
}
impl ExactRiemann {
    fn new(left: [f64; 3], right: [f64; 3]) -> Self {
        let f = |p: f64, q: [f64; 3]| {
            if p > q[2] {
                (p - q[2]) * (2.0 / (2.4 * q[0]) / (p + q[2] / 6.0)).sqrt()
            } else {
                2.0 * (1.4 * q[2] / q[0]).sqrt() / 0.4 * ((p / q[2]).powf(1.0 / 7.0) - 1.0)
            }
        };
        let residual = |p| f(p, left) + f(p, right) + right[1] - left[1];
        let (mut lo, mut hi) = (1e-14, left[2].max(right[2]));
        for _ in 0..30 {
            if residual(hi) >= 0.0 {
                break;
            }
            hi *= 2.0;
        }
        assert!(residual(lo) < 0.0 && residual(hi) >= 0.0);
        for _ in 0..100 {
            let p = 0.5 * (lo + hi);
            if residual(p) < 0.0 {
                lo = p;
            } else {
                hi = p;
            }
        }
        let pressure = 0.5 * (lo + hi);
        let velocity = 0.5 * (left[1] + right[1] + f(pressure, right) - f(pressure, left));
        if left == [1.0, 0.0, 1.0] && right == [0.125, 0.0, 0.1] {
            assert!((pressure - 0.30313017805064685).abs() < 1e-13);
            assert!((velocity - 0.9274526200489499).abs() < 1e-13);
        }
        Self {
            left,
            right,
            pressure,
            velocity,
        }
    }
    fn density(&self, xi: f64) -> f64 {
        let is_left = xi <= self.velocity;
        let q = if is_left { self.left } else { self.right };
        let direction = if is_left { -1.0 } else { 1.0 };
        let a = (1.4 * q[2] / q[0]).sqrt();
        let ratio = self.pressure / q[2];
        if ratio > 1.0 {
            let s = q[1] + direction * a * ((6.0 / 7.0) * ratio + 1.0 / 7.0).sqrt();
            if direction * (xi - s) > 0.0 {
                q[0]
            } else {
                q[0] * (ratio + 1.0 / 6.0) / (ratio / 6.0 + 1.0)
            }
        } else {
            let astar = a * ratio.powf(1.0 / 7.0);
            let head = q[1] + direction * a;
            let tail = self.velocity + direction * astar;
            if direction * (xi - head) >= 0.0 {
                q[0]
            } else if direction * (xi - tail) <= 0.0 {
                q[0] * ratio.powf(1.0 / 1.4)
            } else {
                let local_a = (2.0 / 2.4) * (a + direction * 0.2 * (xi - q[1]));
                q[0] * (local_a / a).powf(5.0)
            }
        }
    }
}
