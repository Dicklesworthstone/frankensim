//! Independent Cartesian diffusion symbols, reservoir efforts and publication.
use super::*;
use crate::compressible::mixture::molecular::StabilizedDissipationControl;

fn controls(step: f64) -> StabilizedDissipationControl {
    StabilizedDissipationControl {
        transport: DissipationControl {
            maximum_step_s: step,
            ..DissipationControl::default()
        },
        ..StabilizedDissipationControl::default()
    }
}

fn prepared(
    flow: &MixtureEuler2d,
    law: ConstantTransport,
) -> MolecularTransport2d<ConstantTransport> {
    let source = PlanarTransport::new(flow.gas().clone(), law, || false).unwrap();
    let budget = MolecularTransport2d::required_stabilized_workspace_bytes(flow, &source).unwrap();
    MolecularTransport2d::new(flow, source, budget, || false).unwrap()
}

// Evaluate 1-b_s+b_s P_s(1+w_1*z) by the ordinary Legendre polynomial.
// No original stage coefficients or production recurrence enter the oracle.
fn amplification(stages: usize, z: f64) -> f64 {
    let s = stages as f64;
    let b = (s * s + s - 2.0) / (2.0 * s * (s + 1.0));
    let x = 1.0 + 4.0 * z / (s * s + s - 2.0);
    let (mut older, mut previous) = (1.0, x);
    for j in 2..=stages {
        let jf = j as f64;
        let next = ((2.0 * jf - 1.0) * x * previous - (jf - 1.0) * older) / jf;
        older = previous;
        previous = next;
    }
    1.0 - b + b * previous
}

#[test]
fn g1_stabilized_oblique_species_mode_matches_two_axis_legendre_symbol() {
    let grid = CartesianGrid::new(16, 12, 0.05, 0.08, 0.3).unwrap();
    let diffusion = 0.02;
    let duration = 0.5;
    let wave = |x: [f64; 2]| (2.0 * std::f64::consts::PI * (x[0] / 0.8 + 2.0 * x[1] / 0.96)).sin();
    let make = || {
        field(equal_mass_gas(), grid, PERIODIC, |x| {
            let a = 0.5 + 0.1 * wave(x);
            ([a, 1.0 - a, 0.0], [0.0; 2], 1.0)
        })
    };
    let mut flow = make();
    let mut op = prepared(&flow, law(flow.gas(), diffusion, 0.0, 0.0, 0.0));
    let report = op
        .advance_stabilized(&mut flow, Time::new(duration), controls(duration), || false)
        .unwrap();
    assert_eq!(report.transport.accepted_steps, 1);
    assert_eq!(report.transport.rejected_steps, 0);
    assert!(report.maximum_stages_used > 8);
    assert_eq!(report.operator_evaluations, report.maximum_stages_used + 1);
    assert!(report.maximum_scaled_parabolic_number <= 0.40000000000001);
    assert!(report.transport.maximum_parabolic_number > 10.0);
    let pi = std::f64::consts::PI;
    let eigenvalue = 4.0
        * diffusion
        * ((pi / 16.0).sin().powi(2) / 0.05_f64.powi(2)
            + (2.0 * pi / 12.0).sin().powi(2) / 0.08_f64.powi(2));
    let decay = amplification(report.maximum_stages_used, -eigenvalue * duration);
    for index in 0..flow.cell_count() {
        let q = flow.cell(index).unwrap();
        let expected = 0.5 + 0.1 * decay * wave(grid.cell_center(index).unwrap());
        close(q.partial_densities[0], expected, 3e-12);
        close(q.partial_densities[1], 1.0 - expected, 3e-12);
        assert_eq!(q.partial_densities[2], 0.0);
        let p = q.primitive(flow.gas(), || false).unwrap();
        close(p.density, 1.0, 3e-12);
        // Unequal Cp and original formation references still retain the exact
        // isothermal mode only when the original partial enthalpies are carried.
        close(p.temperature, 1.0, 3e-10);
    }
    for residual in report.transport.conservation_residual.species_masses {
        close(residual, 0.0, 3e-12);
    }
    close(report.transport.absolute_energy_residual_j, 0.0, 3e-10);
    assert_eq!(flow.time().value(), 0.0);
    assert_eq!(op.last_step(), Some(&report.transport));
    let mut explicit = make();
    let mut explicit_op = prepare(&explicit, law(explicit.gas(), diffusion, 0.0, 0.0, 0.0));
    let baseline = explicit_op
        .advance(
            &mut explicit,
            Time::new(duration),
            controls(duration).transport,
            || false,
        )
        .unwrap();
    // Two RHS calls per accepted SSP step in this constant linear fixture.
    // Compare stability work, not equal temporal error at this large macrostep.
    assert!(4 * report.operator_evaluations < 2 * baseline.accepted_steps);
}

fn thermal_error(step: f64) -> f64 {
    let grid = CartesianGrid::new(12, 8, 0.1, 0.125, 0.3).unwrap();
    let wave = |x: [f64; 2]| (2.0 * std::f64::consts::PI * (x[0] / 1.2 + x[1])).sin();
    let mut flow = field(equal_mass_gas(), grid, PERIODIC, |x| {
        ([1.0, 0.0, 0.0], [0.0; 2], 1.0 + 0.05 * wave(x))
    });
    let conductivity = 0.05;
    let duration = 0.2;
    let mut op = prepared(&flow, law(flow.gas(), 1e-6, 0.0, 0.0, conductivity));
    let report = op
        .advance_stabilized(&mut flow, Time::new(duration), controls(step), || false)
        .unwrap();
    let pi = std::f64::consts::PI;
    let rate = 4.0
        * (conductivity / 2.5)
        * ((pi / 12.0).sin().powi(2) / 0.1_f64.powi(2)
            + (pi / 8.0).sin().powi(2) / 0.125_f64.powi(2));
    let mut error = 0.0_f64;
    for index in 0..flow.cell_count() {
        let q = flow.cell(index).unwrap();
        assert_eq!(q.partial_densities, [1.0, 0.0, 0.0]);
        let expected =
            1.0 + 0.05 * (-rate * duration).exp() * wave(grid.cell_center(index).unwrap());
        let p = q.primitive(flow.gas(), || false).unwrap();
        error = error.max((p.temperature - expected).abs());
    }
    close(report.transport.conservation_residual.energy, 0.0, 3e-12);
    error
}

#[test]
fn g1_stabilized_fourier_mode_refines_to_independent_semidiscrete_exponential() {
    let coarse = thermal_error(0.02);
    let fine = thermal_error(0.01);
    assert!(fine < 1e-6, "fine={fine:e}");
    assert!(fine < 0.3 * coarse, "coarse={coarse:e}, fine={fine:e}");
}

#[test]
fn g1_stabilized_baths_integrate_signed_original_species_and_heat_efforts() {
    let topology = Boundaries2d {
        x: OPEN,
        y: AxisBoundaries::Periodic,
    };
    let grid = CartesianGrid::new(1, 1, 0.1, 0.2, 0.3).unwrap();
    let inputs = [
        ReservoirBoundary {
            side: BoundarySide::XMin,
            partial_densities: &[0.8, 0.2, 0.0],
            velocity: [0.0; 2],
            pressure_pa: 1.0,
        },
        ReservoirBoundary {
            side: BoundarySide::XMax,
            partial_densities: &[0.6, 0.4, 0.0],
            velocity: [0.0; 2],
            pressure_pa: 1.0,
        },
    ];
    let mut flow = super::super::reservoir::bind(
        field(equal_mass_gas(), grid, topology, |_| {
            ([0.2, 0.8, 0.0], [0.0; 2], 1.0)
        }),
        &inputs,
    );
    let duration = 0.4;
    let diffusion = 0.03;
    let mut op = prepared(&flow, law(flow.gas(), diffusion, 0.0, 0.0, 0.0));
    let report = op
        .advance_stabilized(&mut flow, Time::new(duration), controls(duration), || false)
        .unwrap();
    assert_eq!(report.transport.accepted_steps, 1);
    assert!(report.maximum_stages_used > 2);
    let rate = 4.0 * diffusion / grid.spacing[0].powi(2);
    let expected = 0.7 + (0.2 - 0.7) * amplification(report.maximum_stages_used, -rate * duration);
    let q = flow.cell(0).unwrap();
    close(q.partial_densities[0], expected, 3e-12);
    close(
        q.primitive(flow.gas(), || false).unwrap().temperature,
        1.0,
        3e-10,
    );
    let change = grid.cell_volume() * (expected - 0.2);
    let conductance = 2.0 * grid.spacing[1] * grid.depth * diffusion / grid.spacing[0];
    let difference = duration * conductance * (0.8 - 0.6);
    for (side, inward) in [
        (BoundarySide::XMin, 0.5 * (change + difference)),
        (BoundarySide::XMax, 0.5 * (change - difference)),
    ] {
        let exchange = op.boundary_impulses(side).next().unwrap();
        close(exchange.species_masses[0], -inward, 3e-13);
        close(exchange.species_masses[1], inward, 3e-13);
        assert_eq!(exchange.species_masses[2], 0.0);
    }
    close(
        report.transport.boundary_change.species_masses[0],
        change,
        3e-13,
    );
    close(report.transport.boundary_change.mass(), 0.0, 3e-13);
    close(report.transport.absolute_energy_residual_j, 0.0, 3e-10);

    // Pure A with hotter baths isolates Fourier storage; species pressure
    // driving vanishes exactly. The same matrix eigenvalue applies to T.
    let hot = inputs.map(|input| ReservoirBoundary {
        partial_densities: &[1.0, 0.0, 0.0],
        pressure_pa: 3.0,
        ..input
    });
    let mut thermal = super::super::reservoir::bind(
        field(equal_mass_gas(), grid, topology, |_| {
            ([1.0, 0.0, 0.0], [0.0; 2], 2.0)
        }),
        &hot,
    );
    let mut heat = prepared(
        &thermal,
        law(thermal.gas(), 1e-6, 0.0, 0.0, 2.5 * diffusion),
    );
    let report = heat
        .advance_stabilized(
            &mut thermal,
            Time::new(duration),
            controls(duration),
            || false,
        )
        .unwrap();
    let temperature = 3.0 - amplification(report.maximum_stages_used, -rate * duration);
    close(
        thermal
            .cell(0)
            .unwrap()
            .primitive(thermal.gas(), || false)
            .unwrap()
            .temperature,
        temperature,
        3e-11,
    );
    let energy = 2.5 * grid.cell_volume() * (temperature - 2.0);
    close(report.transport.boundary_change.energy, energy, 3e-12);
    for side in [BoundarySide::XMin, BoundarySide::XMax] {
        let exchange = heat.boundary_impulses(side).next().unwrap();
        close(exchange.energy, -0.5 * energy, 3e-12);
        assert_eq!(exchange.species_masses, [0.0; 3]);
    }
}

#[test]
fn g4_stabilized_workspace_work_and_late_cancellation_preserve_accepted_fields() {
    let create = || {
        let mut flow = wall_field();
        // Excess caller storage remains represented after every O(1) swap.
        let mut states = Vec::with_capacity(flow.states.len() + 37);
        states.extend_from_slice(&flow.states);
        flow.workspace_bytes += (states.capacity() - flow.states.capacity()) * size_of::<f64>();
        flow.states = states;
        let mut op = prepared(&flow, law(flow.gas(), 0.02, 0.03, 0.01, 0.02));
        op.advance_stabilized(&mut flow, Time::new(0.001), controls(0.001), || false)
            .unwrap();
        flow.step(Time::new(1e-4), StepControl::default(), || false)
            .unwrap();
        (flow, op)
    };
    let (mut reference, mut reference_op) = create();
    let mut polls = 0;
    let expected = reference_op
        .advance_stabilized(&mut reference, Time::new(0.006), controls(0.002), || {
            polls += 1;
            false
        })
        .unwrap();
    for work_refusal in [false, true] {
        let (mut flow, mut op) = create();
        let before = flow.states().to_vec();
        let time = flow.time();
        let wave_report = flow.last_step().copied();
        let report = op.last_step().copied();
        let loads = impulses(&op);
        let bytes = flow.workspace_bytes() + op.workspace_bytes();
        let mut control = controls(0.002);
        if work_refusal {
            control.maximum_operator_evaluations = 4;
        }
        let mut visits = 0;
        let result = op.advance_stabilized(&mut flow, Time::new(0.006), control, || {
            visits += 1;
            !work_refusal && visits >= polls
        });
        if work_refusal {
            assert!(matches!(
                result,
                Err(DissipationError::OperatorWorkLimit {
                    operator_evaluations: 4,
                    ..
                })
            ));
        } else {
            assert_eq!(result, Err(DissipationError::Cancelled));
        }
        assert_eq!(flow.states(), before);
        assert_eq!(flow.time(), time);
        assert_eq!(flow.last_step().copied(), wave_report);
        assert_eq!(op.last_step().copied(), report);
        assert_eq!(impulses(&op), loads);
        assert_eq!(flow.workspace_bytes() + op.workspace_bytes(), bytes);
        let replay = op
            .advance_stabilized(&mut flow, Time::new(0.006), controls(0.002), || false)
            .unwrap();
        assert_eq!(replay.operator_evaluations, expected.operator_evaluations);
        for (actual, expected) in flow.states().iter().zip(reference.states()) {
            close(*actual, *expected, 3e-11);
        }
        assert_eq!(flow.workspace_bytes() + op.workspace_bytes(), bytes);
    }

    let mut flow = wall_field();
    let source = PlanarTransport::new(
        flow.gas().clone(),
        law(flow.gas(), 0.02, 0.03, 0.01, 0.02),
        || false,
    )
    .unwrap();
    let required =
        MolecularTransport2d::required_stabilized_workspace_bytes(&flow, &source).unwrap();
    let mut op = MolecularTransport2d::new(&flow, source, required - 1, || false).unwrap();
    let before = flow.states().to_vec();
    let bytes = flow.workspace_bytes() + op.workspace_bytes();
    assert!(matches!(
        op.advance_stabilized(&mut flow, Time::new(0.001), controls(0.001), || false),
        Err(DissipationError::Invalid(_))
    ));
    assert_eq!(flow.states(), before);
    assert_eq!(
        flow.workspace_bytes() + op.workspace_bytes(),
        bytes,
        "insufficient budget refused before lazy reservation"
    );
    assert!(op.last_step().is_none());
}
