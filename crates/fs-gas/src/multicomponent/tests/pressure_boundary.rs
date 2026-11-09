//! G0/G3/G4: pressure discharge carries the accepted species and energy state.

use super::*;

fn ends(pressure: f64) -> Boundaries {
    Boundaries::Ends {
        left: Boundary::PressureOutlet {
            pressure_pa: pressure,
        },
        right: Boundary::PressureOutlet {
            pressure_pa: pressure,
        },
    }
}

fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 2e-11 * (1.0 + a.abs().max(b.abs())),
        "{a} != {b}"
    );
}

#[test]
fn variable_caloric_pressure_outlets_share_reservoir_outflow_and_species_balances() {
    for reconstruction in [Reconstruction::FirstOrder, Reconstruction::Muscl] {
        let source = gas(vec![
            source("A", R, 3.5, 0.1, 1e18),
            source("B", 2.0 * R, 4.2, 0.2, -1e18),
        ]);
        let mut declared = flow(source.clone(), 20, ends(0.7), |x| {
            ([0.4, 0.6], 0.4 * (x - 0.5), 1.0)
        });
        let mut reservoir = flow(
            source.clone(),
            20,
            Boundaries::Ends {
                left: Boundary::Transmissive,
                right: Boundary::Transmissive,
            },
            |x| ([0.4, 0.6], 0.4 * (x - 0.5), 1.0),
        );
        // Outgoing composition comes from the field. The reservoir deliberately
        // has different species fractions: its state cannot contaminate outflow.
        let exterior = source.conserved(&[0.7, 0.0], 0.0, 0.7).unwrap();
        let controls = StepControl {
            reconstruction,
            ..StepControl::default()
        };
        for _ in 0..10 {
            let a = declared.step(Time::new(0.001), controls, || false).unwrap();
            let b = reservoir
                .step_with_reservoirs(
                    Time::new(0.001),
                    controls,
                    [Some(exterior.as_cell()); 2],
                    || false,
                )
                .unwrap();
            assert_eq!(a.elapsed, b.elapsed);
            assert!(a.boundary_change.mass() < 0.0);
            for i in 0..declared.cell_count() {
                let x = declared.cell(i).unwrap();
                let y = reservoir.cell(i).unwrap();
                for (&dx, &dy) in x.partial_densities.iter().zip(y.partial_densities) {
                    near(dx, dy);
                }
                near(x.momentum, y.momentum);
                near(x.energy, y.energy);
                near(x.partial_densities[0] / x.density(), 0.4);
            }
            for i in 0..source.species_count() {
                near(
                    a.boundary_change.species_masses[i],
                    b.boundary_change.species_masses[i],
                );
                near(a.conservation_residual.species_masses[i], 0.0);
            }
            near(a.boundary_change.energy, b.boundary_change.energy);
            near(a.conservation_residual.momentum, 0.0);
            near(a.conservation_residual.energy, 0.0);
        }
    }
}

#[test]
fn vacuum_discharge_preserves_each_species_and_matches_low_backpressure() {
    let source = standard(0.0);
    let mut vacuum = flow(source.clone(), 12, ends(0.0), |_| ([0.4, 0.6], 0.0, 1.0));
    let mut low_pressure = flow(source, 12, ends(1e-12), |_| ([0.4, 0.6], 0.0, 1.0));
    for _ in 0..4 {
        let report = vacuum
            .step(Time::new(0.001), StepControl::default(), || false)
            .unwrap();
        low_pressure
            .step(Time::new(0.001), StepControl::default(), || false)
            .unwrap();
        assert_eq!(vacuum.states, low_pressure.states);
        let loss = report.boundary_change.mass();
        assert!(loss < 0.0);
        for (i, fraction) in [0.4, 0.6].into_iter().enumerate() {
            near(report.boundary_change.species_masses[i], fraction * loss);
            near(report.conservation_residual.species_masses[i], 0.0);
        }
        near(report.conservation_residual.energy, 0.0);
    }
}

#[test]
fn mixture_pressure_history_and_refusals_preserve_the_accepted_field() {
    let mut gas_flow = flow(standard(0.0), 12, ends(1.0), |_| ([0.4, 0.6], 0.0, 1.0));
    let original = gas_flow.states.clone();
    let controls = StepControl {
        maximum_halvings: 0,
        ..StepControl::default()
    };
    let rest = gas_flow.step(Time::new(0.001), controls, || false).unwrap();
    for (&state, &reference) in gas_flow.states.iter().zip(&original) {
        near(state, reference);
    }
    near(rest.boundary_change.mass(), 0.0);
    near(rest.boundary_change.energy, 0.0);
    let original = gas_flow.states.clone();
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(gas_flow.set_boundaries(ends(invalid), || false).is_err());
        assert_eq!(gas_flow.boundaries(), ends(1.0));
    }
    assert_eq!(
        gas_flow.set_boundaries(ends(0.7), || true),
        Err(Error::Cancelled)
    );
    assert_eq!(gas_flow.boundaries(), ends(1.0));
    gas_flow.set_boundaries(ends(2.0), || false).unwrap();
    let time = gas_flow.time();
    assert!(gas_flow.step(Time::new(0.001), controls, || false).is_err());
    assert_eq!(gas_flow.states, original);
    assert_eq!(gas_flow.time(), time);
    gas_flow.set_boundaries(ends(0.7), || false).unwrap();
    assert!(
        gas_flow
            .step(Time::new(0.001), controls, || false)
            .unwrap()
            .boundary_change
            .mass()
            < 0.0
    );
}
