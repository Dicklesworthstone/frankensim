//! G1/G3/G4: source isentropes, independent volume integrals and evolved balance.

use super::*;
use crate::dynamics::{Boundary, Reconstruction, StepControl};
use fs_qty::Time;

const CLOSED: Boundaries = Boundaries::Ends {
    left: Boundary::Reflecting,
    right: Boundary::Reflecting,
};

fn reference(g: f64) -> IsentropicAtmosphere {
    IsentropicAtmosphere {
        acceleration: Acceleration::new(g),
        reference_temperature: Temperature::new(2.0),
        reference_pressure: Pressure::new(3.0),
        reference_position: Length::new(0.0),
    }
}

fn column(cp_linear: f64, g: f64, formation: f64, n: usize) -> Euler1d {
    Euler1d::new_isentropic_atmosphere(
        crate::dynamics::tests::gas(cp_linear, formation),
        Mesh1d::uniform(n, Length::new(1.0 / n as f64), 0.7).unwrap(),
        reference(g),
        AtmosphereQuadrature::default(),
        CLOSED,
        || false,
    )
    .unwrap()
}

#[test]
fn g1_constant_cp_column_matches_analytic_polytrope_and_energy_integrals() {
    for g in [-3.0, 3.0] {
        let flow = column(0.0, g, 0.0, 8);
        let balance = flow.hydrostatic.as_ref().unwrap();
        for (i, face) in balance.faces.iter().enumerate() {
            let x = flow.mesh().positions()[i];
            let t = 2.0 + g * x / 3.5;
            let p = 3.0 * (t / 2.0).powf(3.5);
            assert!((face.1.temperature - t).abs() < 3e-13);
            assert!((face.1.pressure - p).abs() < 3e-12);
        }
        for (i, q) in flow.states().iter().enumerate() {
            let left = flow.mesh().positions()[i];
            let right = flow.mesh().positions()[i + 1];
            let tl = 2.0 + g * left / 3.5;
            let tr = 2.0 + g * right / 3.5;
            let mass_per_area = 3.0 * ((tr / 2.0).powf(3.5) - (tl / 2.0).powf(3.5)) / g;
            let pressure_integral =
                3.0 * 2.0 / ((g / 3.5) * 4.5) * ((tr / 2.0).powf(4.5) - (tl / 2.0).powf(4.5));
            let energy_per_area = 2.5 * (pressure_integral - 0.001 * mass_per_area);
            assert!((q.density * (right - left) - mass_per_area).abs() < 2e-12);
            assert!((q.energy * (right - left) - energy_per_area).abs() < 1e-11);
        }
    }
}

#[test]
fn g1_variable_cp_column_matches_independent_high_precision_cell_integrals() {
    // Independent oracle: H(T)=3.5T+.2T²=7.8-3x,
    // p=3(T/2)^3.5 exp[.4(T-2)], rho=p/T. Integrate rho and
    // rho[2.5(T-.001)+.2(T²-.001²)] after dx=(3.5+.4T)dT/(-3).
    // A 70-digit Decimal power-series integral in T and independent SciPy
    // quadrature in x agreed to 4.45e-16. Neither uses the production solver.
    let mesh = Mesh1d::new(vec![0.0, 0.25, 1.0], vec![0.7; 3], vec![0.175, 0.525]).unwrap();
    let flow = Euler1d::new_isentropic_atmosphere(
        crate::dynamics::tests::gas(0.4, 1e18),
        mesh,
        reference(-3.0),
        AtmosphereQuadrature::default(),
        CLOSED,
        || false,
    )
    .unwrap();
    let expected = [
        (1.298_458_824_526_370_7, 7.173_365_842_652_091),
        (0.692_120_249_840_788_5, 3.128_895_909_685_369_6),
    ];
    for (q, (density, energy)) in flow.states().iter().zip(expected) {
        assert!((q.density - density).abs() < 2e-12, "{q:?}");
        assert!((q.energy - energy).abs() < 1e-11, "{q:?}");
        assert_eq!(q.momentum, 0.0);
    }
    assert!((flow.totals().mass - 0.590_593_425_458_528_9).abs() < 2e-12);
    assert!((flow.totals().energy - 2.898_009_375_048_935).abs() < 1e-11);
    let balance = flow.hydrostatic.as_ref().unwrap();
    for (face, pressure) in
        balance
            .faces
            .iter()
            .zip([3.0, 2.026_155_881_605_222, 0.468_885_319_463_447_8])
    {
        assert!((face.1.pressure - pressure).abs() < 3e-12);
        let entropy = flow
            .gas()
            .source()
            .specific_entropy_change(
                Temperature::new(2.0),
                Pressure::new(3.0),
                Temperature::new(face.1.temperature),
                Pressure::new(face.1.pressure),
            )
            .unwrap()
            .value();
        assert!(entropy.abs() < 2e-12);
    }
    // Cp(T), not Cp(T_ref), controls the lapse and final pressure.
    let frozen_cp_pressure = 3.0_f64 * (1.0_f64 - 3.0 / (4.3 * 2.0)).powf(4.3);
    assert!((balance.faces[2].1.pressure - frozen_cp_pressure).abs() > 0.001);
    assert_eq!(
        flow.hydrostatic_reference(),
        Some(HydrostaticReference::Isentropic(reference(-3.0)))
    );
    assert_eq!(flow.atmosphere(), None);
}

#[test]
fn g1_isentropic_rest_is_preserved_and_reports_physical_gravitational_impulse() {
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        for g in [-3.0, 0.0, 3.0] {
            let mut flow = column(0.4, g, 0.0, 16);
            let initial = flow.states().to_vec();
            let mass = flow.totals().mass;
            for _ in 0..30 {
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
                assert_eq!(report.halvings, 0);
                let expected = g * mass * report.elapsed.value();
                assert!((report.gravity_momentum_impulse - expected).abs() < 3e-12);
            }
            assert_eq!(flow.states(), initial);
        }
    }
}

#[test]
fn g3_isentropic_disturbances_evolve_with_conserved_total_energy_and_source_gauge() {
    let mut a = column(0.4, -3.0, 0.0, 24);
    let mut b = column(0.4, -3.0, 1e18, 24);
    assert_eq!(a.states(), b.states());
    for flow in [&mut a, &mut b] {
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
        for _ in 0..30 {
            let report = flow
                .step(Time::new(0.005), StepControl::default(), || false)
                .unwrap();
            assert!(report.conservation_residual.energy.abs() < 3e-13);
        }
        assert_ne!(flow.states(), initial);
        assert!((flow.totals().mass - mass).abs() < 3e-13);
        assert!(
            (flow.totals().energy + flow.gravitational_potential_energy() - energy).abs() < 3e-12
        );
    }
    assert_eq!(a.states(), b.states());
}

#[test]
fn g0_isentropic_domain_quadrature_and_unbound_network_refusals() {
    let mesh = || Mesh1d::uniform(2, Length::new(0.5), 0.7).unwrap();
    let make = |reference, quadrature, boundaries| {
        Euler1d::new_isentropic_atmosphere(
            crate::dynamics::tests::gas(0.4, 0.0),
            mesh(),
            reference,
            quadrature,
            boundaries,
            || false,
        )
    };
    assert!(matches!(
        make(reference(-20.0), AtmosphereQuadrature::default(), CLOSED),
        Err(Error::CaloricPath(_))
    ));
    assert!(
        make(
            reference(-3.0),
            AtmosphereQuadrature::default(),
            Boundaries::Periodic
        )
        .is_err()
    );
    for quadrature in [
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
        AtmosphereQuadrature {
            relative_tolerance: 1e-16,
            maximum_panels: 2,
        },
    ] {
        assert!(make(reference(-3.0), quadrature, CLOSED).is_err());
    }
    let flow = column(0.4, -3.0, 0.0, 4);
    assert!(
        crate::dynamics::network::GasNetwork1d::new(vec![flow], vec![], vec![], || false).is_err()
    );
}

#[test]
fn g4_isentropic_preparation_and_time_step_are_cancellable_and_restartable() {
    let mut polls = 0;
    let result = Euler1d::new_isentropic_atmosphere(
        crate::dynamics::tests::gas(0.4, 0.0),
        Mesh1d::uniform(8, Length::new(0.125), 0.7).unwrap(),
        reference(-3.0),
        AtmosphereQuadrature::default(),
        CLOSED,
        || {
            polls += 1;
            polls == 50
        },
    );
    assert!(matches!(result, Err(Error::Cancelled)));
    let mut interrupted = column(0.4, -3.0, 0.0, 8);
    let mut uninterrupted = column(0.4, -3.0, 0.0, 8);
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
    let b = uninterrupted
        .step(Time::new(0.01), StepControl::default(), || false)
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(interrupted.states(), uninterrupted.states());
}
