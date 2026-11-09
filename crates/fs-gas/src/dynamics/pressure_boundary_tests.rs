//! G0/G1/G3/G4 pressure-boundary integration through the original RK owner.

use super::*;
use crate::dynamics::tests::gas;

fn ends(left: f64, right: f64) -> Boundaries {
    Boundaries::Ends {
        left: Boundary::PressureOutlet { pressure_pa: left },
        right: Boundary::PressureOutlet { pressure_pa: right },
    }
}

fn uniform(boundaries: Boundaries) -> Euler1d {
    let source = gas(0.2, 0.0);
    let state = source.conserved(1.0, 0.0, 1.0).unwrap();
    Euler1d::new(
        source,
        Length::new(0.05),
        0.7,
        vec![state; 20],
        boundaries,
        || false,
    )
    .unwrap()
}

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-11 * (1.0 + a.abs().max(b.abs())),
        "{a} != {b}"
    );
}

#[test]
fn equal_pressure_outlets_preserve_rest_and_their_opposed_support_impulses() {
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        let mut flow = uniform(ends(1.0, 1.0));
        let original = flow.states().to_vec();
        for _ in 0..8 {
            let report = flow
                .step(
                    Time::new(0.001),
                    StepControl {
                        reconstruction,
                        ..StepControl::default()
                    },
                    || false,
                )
                .unwrap();
            for (state, reference) in flow.states().iter().zip(&original) {
                close(state.density, reference.density);
                close(state.momentum, reference.momentum);
                close(state.energy, reference.energy);
            }
            close(report.boundary_change.mass, 0.0);
            close(report.boundary_change.momentum, 0.0);
            close(report.boundary_change.energy, 0.0);
            close(
                report.left_boundary_change.momentum,
                0.7 * report.elapsed.value(),
            );
            close(
                report.right_boundary_change.momentum,
                -0.7 * report.elapsed.value(),
            );
            close(report.conservation_residual.mass, 0.0);
            close(report.conservation_residual.momentum, 0.0);
            close(report.conservation_residual.energy, 0.0);
        }
    }
}

#[test]
fn declared_outlets_match_public_characteristic_ports_at_both_rk_stages() {
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        let source = gas(0.2, 0.0);
        let states: Vec<_> = (0..20)
            .map(|i| {
                let x = (i as f64 + 0.5) / 20.0;
                source
                    .conserved(1.0 + 0.1 * x, 0.4 * (x - 0.5), 1.0 + 0.1 * x)
                    .unwrap()
            })
            .collect();
        let mut declared = Euler1d::new(
            source.clone(),
            Length::new(0.05),
            0.7,
            states.clone(),
            ends(0.8, 0.7),
            || false,
        )
        .unwrap();
        let mut explicit = Euler1d::new(
            source.clone(),
            Length::new(0.05),
            0.7,
            states,
            Boundaries::Ends {
                left: Boundary::Transmissive,
                right: Boundary::Transmissive,
            },
            || false,
        )
        .unwrap();
        let controls = StepControl {
            reconstruction,
            ..StepControl::default()
        };
        for _ in 0..12 {
            let a = declared.step(Time::new(0.001), controls, || false).unwrap();
            let mut stages = 0;
            let b = explicit
                .step_with_ports(
                    Time::new(0.001),
                    controls,
                    |stage| {
                        stages += 1;
                        Ok([
                            Some(source.pressure_outlet_flux(stage.left.0, -1.0, 0.8, || false)?),
                            Some(source.pressure_outlet_flux(stage.right.0, 1.0, 0.7, || false)?),
                        ])
                    },
                    || false,
                )
                .unwrap();
            assert_eq!(stages, 2);
            assert_eq!(a.elapsed, b.elapsed);
            assert!(a.boundary_change.mass < 0.0);
            for (x, y) in declared.states().iter().zip(explicit.states()) {
                close(x.density, y.density);
                close(x.momentum, y.momentum);
                close(x.energy, y.energy);
            }
            close(a.boundary_change.mass, b.boundary_change.mass);
            close(a.boundary_change.momentum, b.boundary_change.momentum);
            close(a.boundary_change.energy, b.boundary_change.energy);
            close(a.conservation_residual.mass, 0.0);
            close(a.conservation_residual.momentum, 0.0);
            close(a.conservation_residual.energy, 0.0);
        }
    }
}

#[test]
fn vacuum_ends_match_choked_low_pressure_discharge_with_conservative_rk_steps() {
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        let mut vacuum = uniform(ends(0.0, 0.0));
        let mut low_pressure = uniform(ends(1e-12, 1e-12));
        let controls = StepControl { reconstruction, ..StepControl::default() };
        for _ in 0..4 {
            let report = vacuum.step(Time::new(0.001), controls, || false).unwrap();
            let low = low_pressure.step(Time::new(0.001), controls, || false).unwrap();
            assert_eq!(vacuum.states(), low_pressure.states());
            assert_eq!(report.boundary_change, low.boundary_change);
            assert!(report.boundary_change.mass < 0.0);
            assert!(report.boundary_change.energy < 0.0);
            close(report.conservation_residual.mass, 0.0);
            close(report.conservation_residual.momentum, 0.0);
            close(report.conservation_residual.energy, 0.0);
        }
    }
}

#[test]
fn invalid_pressure_backflow_and_cancelled_changes_publish_nothing() {
    let mut flow = uniform(ends(1.0, 1.0));
    let states = flow.states().to_vec();
    let boundaries = flow.boundaries();
    for pressure in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(flow.set_boundaries(ends(pressure, 1.0), || false).is_err());
        assert!(flow.set_boundaries(ends(1.0, pressure), || false).is_err());
        assert_eq!(flow.boundaries(), boundaries);
    }
    assert_eq!(
        flow.set_boundaries(ends(0.8, 0.7), || true),
        Err(Error::Cancelled)
    );
    assert_eq!(flow.boundaries(), boundaries);
    flow.set_boundaries(ends(2.0, 2.0), || false).unwrap();
    let controls = StepControl {
        maximum_halvings: 0,
        ..StepControl::default()
    };
    assert!(flow.step(Time::new(0.001), controls, || false).is_err());
    assert_eq!(flow.states(), states);
    assert_eq!(flow.time(), Time::new(0.0));
    flow.set_boundaries(ends(0.8, 0.7), || false).unwrap();
    assert_eq!(
        flow.step(Time::new(0.001), controls, || true),
        Err(Error::Cancelled)
    );
    assert_eq!(flow.states(), states);
    assert_eq!(flow.time(), Time::new(0.0));
    assert!(
        flow.step(Time::new(0.001), controls, || false)
            .unwrap()
            .boundary_change
            .mass
            < 0.0
    );
}
