//! G1/G3/G4: analytic atmospheres, real perturbation evolution and conservation.

use super::*;
use crate::dynamics::{Boundary, Reconstruction, StepControl};
use fs_qty::Time;

const CLOSED: Boundaries = Boundaries::Ends {
    left: Boundary::Reflecting,
    right: Boundary::Reflecting,
};

fn reference(g: f64) -> IsothermalAtmosphere {
    IsothermalAtmosphere {
        acceleration: Acceleration::new(g),
        temperature: Temperature::new(2.0),
        reference_pressure: Pressure::new(3.0),
        reference_position: Length::new(0.0),
    }
}

fn stretched_mesh(n: usize, shift: f64) -> Mesh1d {
    let positions: Vec<_> = (0..=n)
        .map(|i| shift + (i as f64 / n as f64).powi(2))
        .collect();
    let volumes = positions.windows(2).map(|p| 0.7 * (p[1] - p[0])).collect();
    Mesh1d::new(positions, vec![0.7; n + 1], volumes).unwrap()
}

fn atmosphere(n: usize, g: f64, formation: f64) -> Euler1d {
    Euler1d::new_isothermal_atmosphere(
        crate::dynamics::tests::gas(0.4, formation),
        stretched_mesh(n, 0.0),
        reference(g),
        CLOSED,
        || false,
    )
    .unwrap()
}

#[test]
fn g0_identical_and_stationary_contact_fluxes_are_exact() {
    let gas = crate::dynamics::tests::gas(0.4, 0.0);
    for velocity in [0.0, -0.17, 0.31] {
        let q = gas.conserved(1.7, velocity, 2.9).unwrap();
        let p = gas.primitive(q, || false).unwrap();
        let f = gas.hllc_flux(q, q, || false).unwrap();
        assert_eq!(f.mass, q.momentum);
        assert_eq!(f.momentum, q.momentum * p.velocity + p.pressure);
        assert_eq!(f.energy, (q.energy + p.pressure) * p.velocity);
        assert!(!f.used_hll_fallback);
    }
    let f = gas
        .hllc_normal_flux_from_primitive(
            crate::dynamics::NormalPrimitive {
                density: 1.25,
                normal_velocity: 0.0,
                transverse_velocity: 0.4,
                pressure: 2.5,
            },
            crate::dynamics::NormalPrimitive {
                density: 2.5,
                normal_velocity: 0.0,
                transverse_velocity: -0.9,
                pressure: 2.5,
            },
            || false,
        )
        .unwrap();
    assert_eq!(f.mass, 0.0);
    assert_eq!(f.normal_momentum, 2.5);
    assert_eq!(f.transverse_momentum, 0.0);
    assert_eq!(f.energy, 0.0);
}

#[test]
fn g1_exact_cell_averages_and_supported_thermodynamics() {
    for g in [-4.0, 4.0, 0.0, 1e-9] {
        let flow = atmosphere(19, g, 1e18);
        let b = g / 2.0;
        let expected_mass = 0.7 * 1.5 * if b == 0.0 { 1.0 } else { b.exp_m1() / b };
        assert!((flow.totals().mass - expected_mass).abs() < 3e-14 * expected_mass);
        for (i, q) in flow.states().iter().enumerate() {
            let l = flow.mesh().positions()[i];
            let r = flow.mesh().positions()[i + 1];
            let average = if b == 0.0 {
                1.5
            } else {
                1.5 * (b * l).exp() * (b * (r - l)).exp_m1() / (b * (r - l))
            };
            assert!((q.density - average).abs() < 3e-14 * average);
            let p = flow.gas().primitive(*q, || false).unwrap();
            assert!((p.temperature - 2.0).abs() < 1e-11);
        }
    }
}

#[test]
fn g1_resting_stratified_atmosphere_survives_both_reconstructions() {
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        for g in [-4.0, 4.0] {
            let mut flow = atmosphere(32, g, 0.0);
            let before = flow.states().to_vec();
            for _ in 0..40 {
                let r = flow
                    .step(
                        Time::new(0.1),
                        StepControl {
                            reconstruction,
                            ..StepControl::default()
                        },
                        || false,
                    )
                    .unwrap();
                assert_eq!(r.halvings, 0);
                assert_eq!(r.gravity_work, 0.0);
                assert!(r.conservation_residual.momentum.abs() < 2e-17);
            }
            assert_eq!(flow.states(), before);
        }
    }
}

#[test]
fn g1_disturbed_atmosphere_evolves_and_conserves_gas_plus_potential_energy() {
    let mut flow = atmosphere(40, -3.0, 0.0);
    let gas = flow.gas().clone();
    flow.apply_cell_sources(
        |i, q, p| {
            let bump = if (15..25).contains(&i) { 1.02 } else { 1.0 };
            gas.conserved(q.density, 0.0, bump * p.pressure)
        },
        || false,
    )
    .unwrap();
    let initial = flow.states().to_vec();
    let mass = flow.totals().mass;
    let energy = flow.totals().energy + flow.gravitational_potential_energy();
    let mut work = 0.0;
    let gas_energy = flow.totals().energy;
    for _ in 0..80 {
        let r = flow
            .step(Time::new(0.02), StepControl::default(), || false)
            .unwrap();
        work += r.gravity_work;
        assert!(r.conservation_residual.mass.abs() < 1e-14);
        assert!(r.conservation_residual.energy.abs() < 3e-13);
        assert!(r.conservation_residual.momentum.abs() < 1e-14);
    }
    assert_ne!(flow.states(), initial);
    assert!(flow.states().iter().any(|q| q.momentum.abs() > 1e-5));
    assert!((flow.totals().mass - mass).abs() < 2e-14);
    assert!((flow.totals().energy + flow.gravitational_potential_energy() - energy).abs() < 3e-12);
    assert!((flow.totals().energy - gas_energy - work).abs() < 3e-12);
}

#[test]
fn g3_formation_reference_and_coordinate_translation_leave_motion_unchanged() {
    let mut a = atmosphere(16, -2.0, 0.0);
    let mut b = Euler1d::new_isothermal_atmosphere(
        crate::dynamics::tests::gas(0.4, 1e18),
        stretched_mesh(16, 8.0),
        IsothermalAtmosphere {
            reference_position: Length::new(8.0),
            ..reference(-2.0)
        },
        CLOSED,
        || false,
    )
    .unwrap();
    for field in [&mut a, &mut b] {
        let gas = field.gas().clone();
        field
            .apply_cell_sources(
                |i, q, p| gas.conserved(q.density, if i == 7 { 0.03 } else { 0.0 }, p.pressure),
                || false,
            )
            .unwrap();
        for _ in 0..12 {
            field
                .step(Time::new(0.001), StepControl::default(), || false)
                .unwrap();
        }
    }
    assert_eq!(a.states(), b.states());
    assert_eq!(
        a.gravitational_potential_energy(),
        b.gravitational_potential_energy()
    );
}

#[test]
fn g4_cancellation_keeps_state_and_clock_and_restart_is_deterministic() {
    let mut flow = atmosphere(12, -2.0, 0.0);
    let before = flow.states().to_vec();
    let mut polls = 0;
    assert_eq!(
        flow.step(Time::new(0.01), StepControl::default(), || {
            polls += 1;
            polls == 70
        }),
        Err(Error::Cancelled)
    );
    assert_eq!(flow.states(), before);
    assert_eq!(flow.time().value(), 0.0);
    let mut twin = atmosphere(12, -2.0, 0.0);
    let a = flow
        .step(Time::new(0.01), StepControl::default(), || false)
        .unwrap();
    let b = twin
        .step(Time::new(0.01), StepControl::default(), || false)
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(flow.states(), twin.states());
}

#[test]
fn g1_nonequilibrium_uniform_gas_converges_to_analytic_free_fall() {
    let mut errors = Vec::new();
    for n in [32, 64, 128] {
        let gas = crate::dynamics::tests::gas(0.0, 0.0);
        let mut flow = Euler1d::new_isothermal_atmosphere(
            gas.clone(),
            Mesh1d::uniform(n, Length::new(1.0 / n as f64), 0.7).unwrap(),
            reference(-2.0),
            Boundaries::Ends {
                left: Boundary::Transmissive,
                right: Boundary::Transmissive,
            },
            || false,
        )
        .unwrap();
        // This constant-pressure state is deliberately NOT hydrostatic. In
        // the interior, uniform gravity must accelerate it with u=g*t.
        let uniform = gas.conserved(1.5, 0.0, 3.0).unwrap();
        flow.apply_cell_sources(|_, _, _| Ok(uniform), || false)
            .unwrap();
        let end = 0.005;
        while flow.time().value() < end {
            flow.step(
                Time::new(end - flow.time().value()),
                StepControl::default(),
                || false,
            )
            .unwrap();
        }
        let mut error = 0.0;
        for q in &flow.states()[n / 4..3 * n / 4] {
            let p = gas.primitive(*q, || false).unwrap();
            error += (p.velocity - (-2.0 * end)).abs();
            assert!((p.temperature - 2.0).abs() < 2e-5);
        }
        errors.push(error / (n / 2) as f64);
    }
    assert!(errors[2] < 2e-6, "free-fall errors: {errors:?}");
    assert!(
        errors[0] > 3.0 * errors[1] && errors[1] > 3.0 * errors[2],
        "free-fall errors: {errors:?}"
    );
}

#[test]
fn g0_incompatible_geometry_sources_and_unbound_ports_refuse() {
    let gas = crate::dynamics::tests::gas(0.0, 0.0);
    for invalid in [
        IsothermalAtmosphere {
            acceleration: Acceleration::new(f64::NAN),
            ..reference(-2.0)
        },
        IsothermalAtmosphere {
            temperature: Temperature::new(100.0),
            ..reference(-2.0)
        },
        IsothermalAtmosphere {
            reference_pressure: Pressure::new(-1.0),
            ..reference(-2.0)
        },
        reference(-1e6),
    ] {
        assert!(
            Euler1d::new_isothermal_atmosphere(
                gas.clone(),
                stretched_mesh(4, 0.0),
                invalid,
                CLOSED,
                || false
            )
            .is_err()
        );
    }
    let bad_volume = Mesh1d::new(vec![0.0, 0.5, 1.0], vec![1.0; 3], vec![0.4; 2]).unwrap();
    assert!(
        Euler1d::new_isothermal_atmosphere(
            gas.clone(),
            bad_volume,
            reference(-2.0),
            CLOSED,
            || false
        )
        .is_err()
    );
    assert!(
        Euler1d::new_isothermal_atmosphere(
            gas,
            stretched_mesh(4, 0.0),
            reference(-2.0),
            Boundaries::Periodic,
            || false
        )
        .is_err()
    );
    let mut flow = atmosphere(4, -2.0, 0.0);
    assert!(flow.set_boundaries(Boundaries::Periodic, || false).is_err());
    assert_eq!(flow.boundaries(), CLOSED);
    assert!(
        flow.step_with_ports(
            Time::new(0.01),
            StepControl::default(),
            |_| Ok([None; 2]),
            || false
        )
        .is_err()
    );
    assert!(
        super::super::network::GasNetwork1d::new(vec![flow], vec![], vec![], || false).is_err()
    );
}
