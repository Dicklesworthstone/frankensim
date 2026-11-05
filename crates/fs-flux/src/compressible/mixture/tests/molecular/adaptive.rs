//! Local accuracy against independent modes, physical face efforts and rollback.
use super::*;
use crate::compressible::mixture::molecular::{
    AdaptiveDissipationControl, AdaptiveMolecularReport,
};

fn controls(step: f64, tolerance: f64) -> AdaptiveDissipationControl {
    let mut controls = AdaptiveDissipationControl::default();
    controls.stabilized.transport.maximum_step_s = step;
    controls.relative_tolerance = 0.0;
    controls.absolute_partial_density_tolerance_kg_m3 = tolerance;
    controls.absolute_temperature_tolerance_k = tolerance;
    controls.absolute_velocity_tolerance_m_s = tolerance;
    controls
}

fn prepared(
    flow: &MixtureEuler2d,
    law: ConstantTransport,
) -> MolecularTransport2d<ConstantTransport> {
    let source = PlanarTransport::new(flow.gas().clone(), law, || false).unwrap();
    let bytes = MolecularTransport2d::required_adaptive_workspace_bytes(flow, &source).unwrap();
    MolecularTransport2d::new(flow, source, bytes, || false).unwrap()
}

// Independent Legendre-polynomial evaluation, not the production recurrence.
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
fn g1_adaptive_planar_halves_share_the_full_step_scheme_and_native_rhs() {
    let grid = CartesianGrid::new(16, 12, 0.05, 0.08, 0.3).unwrap();
    let diffusion = 0.02;
    let duration = 0.1;
    let wave = |x: [f64; 2]| (std::f64::consts::TAU * (x[0] / 0.8 + 2.0 * x[1] / 0.96)).sin();
    let mut flow = field(equal_mass_gas(), grid, PERIODIC, |x| {
        let a = 0.5 + 0.1 * wave(x);
        ([a, 1.0 - a, 0.0], [0.0; 2], 1.0)
    });
    let source = PlanarTransport::new(
        flow.gas().clone(),
        law(flow.gas(), diffusion, 0.0, 0.0, 0.0),
        || false,
    )
    .unwrap();
    let bytes = MolecularTransport2d::required_adaptive_workspace_bytes(&flow, &source).unwrap();
    let mut op = MolecularTransport2d::new(&flow, source, bytes, || false).unwrap();
    let report = op
        .advance_adaptive(
            &mut flow,
            Time::new(duration),
            controls(duration, 1.0),
            || false,
        )
        .unwrap();
    assert_eq!(report.stabilized.transport.accepted_steps, 1);
    assert_eq!(report.stabilized.transport.rejected_steps, 0);
    assert_eq!(report.temporal_rejections, 0);
    let stages = report.stabilized.maximum_stages_used;
    assert!(stages > 2);
    assert_eq!(report.stabilized.operator_evaluations, 3 * stages + 1);
    assert!(report.maximum_estimated_local_error <= 1.0);
    assert!(report.stabilized.maximum_scaled_parabolic_number <= 0.40000000000001);
    assert_eq!(flow.workspace_bytes() + op.workspace_bytes(), bytes);
    let pi = std::f64::consts::PI;
    let rate = 4.0
        * diffusion
        * ((pi / 16.0).sin().powi(2) / 0.05_f64.powi(2)
            + (2.0 * pi / 12.0).sin().powi(2) / 0.08_f64.powi(2));
    let decay = amplification(stages, -0.5 * rate * duration).powi(2);
    for index in 0..flow.cell_count() {
        let q = flow.cell(index).unwrap();
        let a = 0.5 + 0.1 * decay * wave(grid.cell_center(index).unwrap());
        close(q.partial_densities[0], a, 3e-12);
        close(q.partial_densities[1], 1.0 - a, 3e-12);
        assert_eq!(q.partial_densities[2], 0.0);
        close(
            q.primitive(flow.gas(), || false).unwrap().temperature,
            1.0,
            3e-10,
        );
    }
    close(
        report.stabilized.transport.absolute_energy_residual_j,
        0.0,
        3e-10,
    );
    assert_eq!(flow.time().value(), 0.0);
    assert_eq!(op.last_step(), Some(&report.stabilized.transport));
}

fn relaxation(tolerance: f64, thermal: bool) -> (f64, AdaptiveMolecularReport) {
    let grid = CartesianGrid::new(12, 8, 0.1, 0.125, 0.3).unwrap();
    let wave = |x: [f64; 2]| (std::f64::consts::TAU * (x[0] / 1.2 + x[1])).sin();
    let diffusivity = 0.04;
    let duration = 0.4;
    let mut flow = field(equal_mass_gas(), grid, PERIODIC, |x| {
        if thermal {
            ([1.0, 0.0, 0.0], [0.0; 2], 1.0 + 0.1 * wave(x))
        } else {
            let a = 0.5 + 0.1 * wave(x);
            ([a, 1.0 - a, 0.0], [0.0; 2], 1.0)
        }
    });
    let mut op = prepared(
        &flow,
        law(
            flow.gas(),
            diffusivity,
            0.0,
            0.0,
            if thermal { 2.5 * diffusivity } else { 0.0 },
        ),
    );
    let report = op
        .advance_adaptive(
            &mut flow,
            Time::new(duration),
            controls(duration, tolerance),
            || false,
        )
        .unwrap();
    let pi = std::f64::consts::PI;
    let rate = 4.0
        * diffusivity
        * ((pi / 12.0).sin().powi(2) / 0.1_f64.powi(2)
            + (pi / 8.0).sin().powi(2) / 0.125_f64.powi(2));
    let decay = (-rate * duration).exp();
    let mut error = 0.0_f64;
    for index in 0..flow.cell_count() {
        let q = flow.cell(index).unwrap();
        assert_eq!(q.partial_densities[2], 0.0);
        let exact =
            if thermal { 1.0 } else { 0.5 } + 0.1 * decay * wave(grid.cell_center(index).unwrap());
        let value = if thermal {
            q.primitive(flow.gas(), || false).unwrap().temperature
        } else {
            q.partial_densities[0]
        };
        error = error.max((value - exact).abs());
    }
    assert!(report.maximum_estimated_local_error <= 1.0);
    close(
        report.stabilized.transport.absolute_energy_residual_j,
        0.0,
        3e-9,
    );
    (error, report)
}

#[test]
fn g1_physical_tolerances_refine_species_and_fourier_relaxation() {
    for thermal in [false, true] {
        let (coarse, loose) = relaxation(1e-4, thermal);
        let (fine, tight) = relaxation(1e-7, thermal);
        assert!(fine < 2e-5, "thermal={thermal}, fine={fine:e}");
        assert!(
            fine < 0.25 * coarse,
            "thermal={thermal}, coarse={coarse:e}, fine={fine:e}"
        );
        assert!(tight.temporal_rejections > 0);
        assert!(tight.stabilized.operator_evaluations > loose.stabilized.operator_evaluations);
    }
}

#[test]
fn g1_g3_adaptive_bath_quadrature_retains_opposed_efforts_and_formation_gauge() {
    let grid = CartesianGrid::new(1, 1, 0.1, 0.2, 0.3).unwrap();
    let diffusion = 0.03;
    let duration = 1.0;
    let mut reference = None;
    for formation in [0.0, 1e8] {
        let gas = gas(vec![
            source("A", R, 3.5, 0.0, 40.0 + formation),
            source("B", R, 4.5, 0.0, -50.0 - 2.0 * formation),
            source("absent-C", R, 3.5, 0.0, 2.0 + 3.0 * formation),
        ]);
        let mut flow = super::super::reservoir::bind(
            field(
                gas,
                grid,
                Boundaries2d {
                    x: OPEN,
                    y: AxisBoundaries::Periodic,
                },
                |_| ([0.2, 0.8, 0.0], [0.0; 2], 1.0),
            ),
            &[
                ReservoirBoundary {
                    side: BoundarySide::XMin,
                    partial_densities: &[0.9, 0.1, 0.0],
                    velocity: [0.0; 2],
                    pressure_pa: 1.0,
                },
                ReservoirBoundary {
                    side: BoundarySide::XMax,
                    partial_densities: &[0.7, 0.3, 0.0],
                    velocity: [0.0; 2],
                    pressure_pa: 1.0,
                },
            ],
        );
        let mut op = prepared(&flow, law(flow.gas(), diffusion, 0.0, 0.0, 0.0));
        let report = op
            .advance_adaptive(
                &mut flow,
                Time::new(duration),
                controls(duration, 1e-8),
                || false,
            )
            .unwrap();
        let rate = 4.0 * diffusion / grid.spacing()[0].powi(2);
        let expected = 0.8 + (0.2 - 0.8) * (-rate * duration).exp();
        close(flow.cell(0).unwrap().partial_densities[0], expected, 2e-6);
        let total = grid.cell_volume() * (expected - 0.2);
        let conductance = 2.0 * grid.face_area(BoundarySide::XMin) * diffusion / grid.spacing()[0];
        let difference = duration * conductance * (0.9 - 0.7);
        let left = op.boundary_impulses(BoundarySide::XMin).next().unwrap();
        let right = op.boundary_impulses(BoundarySide::XMax).next().unwrap();
        assert!(left.species_masses[0] < 0.0 && right.species_masses[0] > 0.0);
        for (face, inward) in [
            (left, 0.5 * (total + difference)),
            (right, 0.5 * (total - difference)),
        ] {
            close(face.species_masses[0], -inward, 1e-8);
            close(face.species_masses[1], inward, 1e-8);
            assert_eq!(face.species_masses[2], 0.0);
        }
        assert!(report.maximum_estimated_local_error <= 1.0);
        let signature = (
            flow.states().to_vec(),
            impulses(&op),
            report.stabilized.operator_evaluations,
            report.temporal_rejections,
        );
        if let Some(expected) = &reference {
            assert_eq!(&signature, expected);
        } else {
            reference = Some(signature);
        }
    }
}

#[test]
fn g4_adaptive_complete_work_storage_and_late_cancellation_preserve_publication() {
    let create = || {
        let mut flow = wall_field();
        let mut states = Vec::with_capacity(flow.states.len() + 37);
        states.extend_from_slice(&flow.states);
        flow.workspace_bytes += (states.capacity() - flow.states.capacity()) * size_of::<f64>();
        flow.states = states;
        let mut op = prepared(&flow, law(flow.gas(), 0.02, 0.03, 0.01, 0.02));
        op.advance_adaptive(&mut flow, Time::new(0.001), controls(0.001, 1e-8), || false)
            .unwrap();
        flow.step(Time::new(1e-4), StepControl::default(), || false)
            .unwrap();
        (flow, op)
    };
    let (mut reference, mut reference_op) = create();
    let mut polls = 0;
    let expected = reference_op
        .advance_adaptive(
            &mut reference,
            Time::new(0.006),
            controls(0.002, 1e-8),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for work_refusal in [false, true] {
        let (mut flow, mut op) = create();
        let before = flow.states().to_vec();
        let time = flow.time();
        let wave_report = flow.last_step().copied();
        let report = op.last_step().copied();
        let loads = impulses(&op);
        let bytes = flow.workspace_bytes() + op.workspace_bytes();
        let mut control = controls(0.002, 1e-8);
        if work_refusal {
            control.stabilized.maximum_operator_evaluations =
                expected.stabilized.operator_evaluations - 1;
        }
        let mut visits = 0;
        let result = op.advance_adaptive(&mut flow, Time::new(0.006), control, || {
            visits += 1;
            !work_refusal && visits >= polls
        });
        if work_refusal {
            assert!(
                matches!(result, Err(DissipationError::OperatorWorkLimit { operator_evaluations, .. })
                if operator_evaluations == control.stabilized.maximum_operator_evaluations)
            );
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
            .advance_adaptive(&mut flow, Time::new(0.006), controls(0.002, 1e-8), || false)
            .unwrap();
        assert_eq!(replay, expected);
        assert_eq!(flow.states(), reference.states());
        assert_eq!(impulses(&op), impulses(&reference_op));
    }
    let mut flow = wall_field();
    let source = PlanarTransport::new(
        flow.gas().clone(),
        law(flow.gas(), 0.02, 0.03, 0.01, 0.02),
        || false,
    )
    .unwrap();
    let required = MolecularTransport2d::required_adaptive_workspace_bytes(&flow, &source).unwrap();
    let mut op = MolecularTransport2d::new(&flow, source, required - 1, || false).unwrap();
    let before = flow.states().to_vec();
    let bytes = flow.workspace_bytes() + op.workspace_bytes();
    assert!(matches!(
        op.advance_adaptive(&mut flow, Time::new(0.001), controls(0.001, 1e-8), || false),
        Err(DissipationError::Invalid(_))
    ));
    assert_eq!(flow.states(), before);
    assert_eq!(
        flow.workspace_bytes() + op.workspace_bytes(),
        bytes,
        "all pending groups admitted before lazy allocation"
    );
    assert!(op.last_step().is_none());
}
