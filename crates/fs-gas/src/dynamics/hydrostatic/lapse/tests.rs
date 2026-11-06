//! G1/G3/G4: independent hydrostatics, true cell storage and live gravity.
use super::*;
use crate::dynamics::{Boundary, Reconstruction, StepControl};
use fs_qty::Time;

const CLOSED: Boundaries = Boundaries::Ends {
    left: Boundary::Reflecting,
    right: Boundary::Reflecting,
};

fn reference(g: f64, lapse: f64) -> LapseRateAtmosphere {
    LapseRateAtmosphere {
        acceleration: Acceleration::new(g),
        reference_temperature: Temperature::new(2.0),
        reference_pressure: Pressure::new(3.0),
        reference_position: Length::new(0.0),
        temperature_gradient_k_m: lapse,
    }
}
fn column(g: f64, lapse: f64, formation: f64, n: usize) -> Euler1d {
    Euler1d::new_lapse_rate_atmosphere(
        crate::dynamics::tests::gas(0.4, formation),
        Mesh1d::uniform(n, Length::new(1.0 / n as f64), 0.7).unwrap(),
        reference(g, lapse),
        AtmosphereQuadrature::default(),
        CLOSED,
        || false,
    )
    .unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} versus {b:e}");
}

#[test]
fn g1_lapse_integrals_match_independent_variable_cp_antiderivatives() {
    for g in [-3.0, 3.0] {
        let lapse = -0.4;
        let mesh = Mesh1d::new(vec![0.0, 0.25, 1.0], vec![0.7; 3], vec![0.175, 0.525]).unwrap();
        let flow = Euler1d::new_lapse_rate_atmosphere(
            crate::dynamics::tests::gas(0.4, 1e18),
            mesh,
            reference(g, lapse),
            AtmosphereQuadrature::default(),
            CLOSED,
            || false,
        )
        .unwrap();
        let exponent = g / lapse;
        let power = |x: f64, shift: f64| (1.0 + lapse * x / 2.0).powf(exponent + shift);
        for i in 0..flow.mesh().cell_count() {
            let left = flow.mesh().positions()[i];
            let right = flow.mesh().positions()[i + 1];
            // Integrals of rho, rho*T and rho*T². The source has R=1,
            // Cv(T)=2.5+0.4*T and reference sensible energy e(0.001)=0.
            let mass = 3.0 * (power(right, 0.0) - power(left, 0.0)) / g;
            let first = 6.0 * (power(right, 1.0) - power(left, 1.0)) / (lapse * (exponent + 1.0));
            let second = 12.0 * (power(right, 2.0) - power(left, 2.0)) / (lapse * (exponent + 2.0));
            let energy = 2.5 * first + 0.2 * second - (0.0025 + 0.0000002) * mass;
            close(flow.states()[i].density * (right - left), mass, 3e-12);
            close(flow.states()[i].energy * (right - left), energy, 2e-11);
            let (_, face) = flow.hydrostatic_reference_face(i + 1).unwrap();
            close(face.temperature, 2.0 + lapse * right, 3e-14);
            close(face.pressure, 3.0 * power(right, 0.0), 2e-12);
        }
        assert_eq!(
            flow.hydrostatic_reference(),
            Some(HydrostaticReference::LapseRate(reference(g, lapse)))
        );
        assert_eq!(flow.atmosphere(), None);
    }
}

#[test]
fn g1_zero_gravity_isobaric_and_zero_lapse_isothermal_limits() {
    let flow = column(0.0, -0.4, 0.0, 8);
    // With zero gravity the density still varies inversely with temperature.
    close(
        flow.totals().mass,
        0.7 * (3.0 / -0.4) * (1.6_f64 / 2.0).ln(),
        3e-12,
    );
    for i in 0..=8 {
        close(
            flow.hydrostatic_reference_face(i).unwrap().1.pressure,
            3.0,
            3e-14,
        );
    }
    for g in [-3.0, 0.0, 3.0] {
        let lapse = column(g, 0.0, 0.0, 8);
        let isothermal = Euler1d::new_isothermal_atmosphere(
            crate::dynamics::tests::gas(0.4, 0.0),
            Mesh1d::uniform(8, Length::new(0.125), 0.7).unwrap(),
            IsothermalAtmosphere {
                acceleration: Acceleration::new(g),
                temperature: Temperature::new(2.0),
                reference_pressure: Pressure::new(3.0),
                reference_position: Length::new(0.0),
            },
            CLOSED,
            || false,
        )
        .unwrap();
        assert_eq!(lapse.states(), isothermal.states());
        for tiny in [-1e-14, 1e-14] {
            let nearly = column(g, tiny, 0.0, 8);
            for i in 0..=8 {
                close(
                    nearly.hydrostatic_reference_face(i).unwrap().1.pressure,
                    isothermal.hydrostatic_reference_face(i).unwrap().1.pressure,
                    3e-12,
                );
            }
        }
    }
}

#[test]
fn g1_stable_and_unstable_lapses_keep_source_cp_and_exact_mechanical_rest() {
    let gas = crate::dynamics::tests::gas(0.4, 0.0);
    for lapse in [-0.4, -1.2] {
        let r = reference(-3.0, lapse);
        let x = 0.3;
        let t = 2.0 + lapse * x;
        let expected = (-3.0 / t) * (-3.0 / (3.5 + 0.4 * t) - lapse);
        let n2 = r
            .buoyancy_frequency_squared(&gas, Length::new(x), || false)
            .unwrap();
        close(n2, expected, 2e-14);
        assert_eq!(n2 > 0.0, lapse == -0.4);
        let frozen_cp = (-3.0 / t) * (-3.0 / 4.3 - lapse);
        assert!((n2 - frozen_cp).abs() > 1e-3);
        for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
            let mut flow = column(-3.0, lapse, 0.0, 16);
            let initial = flow.states().to_vec();
            let mass = flow.totals().mass;
            for _ in 0..24 {
                let report = flow
                    .step(
                        Time::new(0.02),
                        StepControl {
                            reconstruction,
                            ..StepControl::default()
                        },
                        || false,
                    )
                    .unwrap();
                assert_eq!(report.gravity_work, 0.0);
                close(
                    report.gravity_momentum_impulse,
                    -3.0 * mass * report.elapsed.value(),
                    3e-11,
                );
                assert_eq!(report.halvings, 0);
            }
            assert_eq!(flow.states(), initial);
        }
    }
}

#[test]
fn g3_perturbations_evolve_with_original_energy_gauge_and_gravity_work() {
    let mut first = column(-3.0, -0.4, 0.0, 24);
    let mut shifted = column(-3.0, -0.4, 1e18, 24);
    assert_eq!(first.states(), shifted.states());
    for flow in [&mut first, &mut shifted] {
        let gas = flow.gas().clone();
        flow.apply_cell_sources(
            |i, q, p| {
                gas.conserved(
                    q.density,
                    if (9..15).contains(&i) { 0.02 } else { 0.0 },
                    p.pressure,
                )
            },
            || false,
        )
        .unwrap();
        let initial = flow.states().to_vec();
        let mass = flow.totals().mass;
        let energy = flow.totals().energy + flow.gravitational_potential_energy();
        for _ in 0..24 {
            let report = flow
                .step(Time::new(0.005), StepControl::default(), || false)
                .unwrap();
            assert!(report.conservation_residual.energy.abs() < 3e-13);
        }
        assert_ne!(flow.states(), initial);
        close(flow.totals().mass, mass, 3e-13);
        close(
            flow.totals().energy + flow.gravitational_potential_energy(),
            energy,
            3e-12,
        );
    }
    assert_eq!(first.states(), shifted.states());
}

#[test]
fn g0_source_domain_endpoint_quadrature_and_boundary_refusals() {
    let make = |r, q, b| {
        Euler1d::new_lapse_rate_atmosphere(
            crate::dynamics::tests::gas(0.4, 0.0),
            Mesh1d::uniform(2, Length::new(0.5), 0.7).unwrap(),
            r,
            q,
            b,
            || false,
        )
    };
    for r in [
        reference(-3.0, -2.0),
        reference(-3.0, 100.0),
        reference(f64::INFINITY, -0.4),
        reference(-3.0, f64::NAN),
        LapseRateAtmosphere {
            reference_pressure: Pressure::new(0.0),
            ..reference(-3.0, -0.4)
        },
        LapseRateAtmosphere {
            reference_position: Length::new(f64::NAN),
            ..reference(-3.0, -0.4)
        },
    ] {
        assert!(make(r, AtmosphereQuadrature::default(), CLOSED).is_err());
    }
    assert!(
        make(
            reference(-3.0, -0.4),
            AtmosphereQuadrature::default(),
            Boundaries::Periodic
        )
        .is_err()
    );
    for q in [
        AtmosphereQuadrature {
            maximum_panels: 1,
            ..AtmosphereQuadrature::default()
        },
        AtmosphereQuadrature {
            maximum_panels: 129,
            ..AtmosphereQuadrature::default()
        },
        AtmosphereQuadrature {
            relative_tolerance: f64::NAN,
            ..AtmosphereQuadrature::default()
        },
    ] {
        assert!(make(reference(-3.0, -0.4), q, CLOSED).is_err());
    }
    assert!(
        make(
            reference(-3.0, -1.99),
            AtmosphereQuadrature {
                maximum_panels: 2,
                relative_tolerance: 1e-16
            },
            CLOSED
        )
        .is_err()
    );
}

#[test]
fn g4_preparation_and_evolution_cancellation_preserve_accepted_state() {
    let mut polls = 0;
    let failed = Euler1d::new_lapse_rate_atmosphere(
        crate::dynamics::tests::gas(0.4, 0.0),
        Mesh1d::uniform(8, Length::new(0.125), 0.7).unwrap(),
        reference(-3.0, -0.4),
        AtmosphereQuadrature::default(),
        CLOSED,
        || {
            polls += 1;
            polls == 30
        },
    );
    assert!(matches!(failed, Err(Error::Cancelled)));
    let mut interrupted = column(-3.0, -0.4, 0.0, 8);
    let mut clean = column(-3.0, -0.4, 0.0, 8);
    let before = interrupted.states().to_vec();
    let time = interrupted.time();
    let mut polls = 0;
    assert_eq!(
        interrupted.step(Time::new(0.01), StepControl::default(), || {
            polls += 1;
            polls == 30
        }),
        Err(Error::Cancelled)
    );
    assert_eq!(interrupted.states(), before);
    assert_eq!(interrupted.time(), time);
    let a = interrupted
        .step(Time::new(0.01), StepControl::default(), || false)
        .unwrap();
    let b = clean
        .step(Time::new(0.01), StepControl::default(), || false)
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(interrupted.states(), clean.states());
}
