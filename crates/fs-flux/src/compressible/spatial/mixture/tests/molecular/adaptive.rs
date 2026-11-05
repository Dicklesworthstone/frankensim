//! G1/G3/G4: physical temporal accuracy and transactional volumetric diffusion.
use super::*;
use crate::compressible::spatial::mixture::ReservoirFaceState;
use crate::compressible::spatial::mixture::molecular::{
    AdaptiveDissipationControl, AdaptiveMolecularReport, InertialWall3d, WallInertiaControl3d,
};
use fs_qty::Mass;

fn accurate(
    flow: &MixtureEuler3d,
    law: ConstantTransport,
    extra: usize,
) -> (MolecularTransport3d<ConstantTransport>, usize) {
    let source = SpatialTransport::new(flow.gas().clone(), law, || false).unwrap();
    let bytes =
        MolecularTransport3d::required_adaptive_workspace_bytes(flow, &source).unwrap() + extra;
    (
        MolecularTransport3d::new(flow, source, bytes, || false).unwrap(),
        bytes,
    )
}

fn accuracy(tolerance: f64) -> AdaptiveDissipationControl {
    AdaptiveDissipationControl {
        relative_tolerance: 0.0,
        absolute_partial_density_tolerance_kg_m3: tolerance,
        absolute_temperature_tolerance_k: tolerance,
        absolute_velocity_tolerance_m_s: tolerance,
        ..AdaptiveDissipationControl::default()
    }
}

fn balanced(report: AdaptiveMolecularReport) {
    assert!(report.maximum_estimated_local_error <= 1.0);
    assert!(report.stabilized.maximum_scaled_parabolic_number <= 0.40000000000001);
    for value in report
        .stabilized
        .transport
        .conservation_residual
        .species_masses
    {
        close(value, 0.0, 4e-11);
    }
    for value in report.stabilized.transport.conservation_residual.momentum {
        close(value, 0.0, 4e-11);
    }
    close(
        report.stabilized.transport.conservation_residual.energy,
        0.0,
        4e-11,
    );
    close(
        report.stabilized.transport.absolute_energy_residual_j,
        0.0,
        4e-9,
    );
    assert!(report.stabilized.transport.wall_inertia.is_none());
}

#[test]
fn g1_adaptive_three_axis_species_and_heat_refine_against_independent_decay() {
    let n = 8;
    let h = 1.0 / n as f64;
    let diffusivity = 0.03;
    let duration = 0.25;
    let grid = CartesianGrid::new(n, n, n, h, h, h).unwrap();
    let mode = |p: [f64; 3]| (2.0 * std::f64::consts::PI * (p[0] + p[1] + p[2])).sin();
    // Independent finite-volume Fourier symbol; comparison removes spatial
    // error and isolates time evolution of the original molecular operators.
    let decay = (-12.0 * diffusivity * (std::f64::consts::PI * h).sin().powi(2) / h.powi(2)
        * duration)
        .exp();
    for species_mode in [true, false] {
        let mut errors = Vec::new();
        let mut work = Vec::new();
        for tolerance in [1e-3, 1e-6] {
            let mut flow = field(equal_mass_gas(), grid, PERIODIC, |p| {
                let wave = 0.15 * mode(p);
                if species_mode {
                    ([0.5 + wave, 0.5 - wave, 0.0], [0.0; 3], 1.0)
                } else {
                    ([1.0, 0.0, 0.0], [0.0; 3], 1.0 + wave)
                }
            });
            let native_law = law(
                flow.gas(),
                diffusivity,
                0.0,
                0.0,
                if species_mode { 0.0 } else { 2.5 * diffusivity },
            );
            let (mut operator, bytes) = accurate(&flow, native_law, 0);
            let report = operator
                .advance_adaptive(&mut flow, Time::new(duration), accuracy(tolerance), || {
                    false
                })
                .unwrap();
            balanced(report);
            assert_eq!(operator.workspace_bytes() + flow.workspace_bytes(), bytes);
            assert_eq!(flow.time().value(), 0.0);
            assert!(flow.last_step().is_none());
            if tolerance == 1e-3 {
                assert!(report.stabilized.maximum_stages_used > 2);
            } else {
                assert!(report.temporal_rejections > 0);
            }
            let mut error = 0.0_f64;
            for i in 0..flow.cell_count() {
                let q = flow.cell(i).unwrap();
                let p = q.primitive(flow.gas(), || false).unwrap();
                let wave = 0.15 * decay * mode(grid.cell_center(i).unwrap());
                error = error.max(if species_mode {
                    close(p.temperature, 1.0, 4e-10);
                    close(p.pressure, 1.0, 4e-10);
                    (q.partial_densities[0] - (0.5 + wave)).abs()
                } else {
                    assert_eq!(q.partial_densities, [1.0, 0.0, 0.0]);
                    (p.temperature - (1.0 + wave)).abs()
                });
                assert_eq!(q.partial_densities[2], 0.0);
                close(p.density, 1.0, 4e-12);
                for value in p.velocity {
                    close(value, 0.0, 4e-12);
                }
            }
            errors.push(error);
            work.push(report.stabilized.operator_evaluations);
        }
        assert!(errors[1] < 0.2 * errors[0], "temporal errors {errors:?}");
        assert!(errors[1] < 5e-5, "temporal errors {errors:?}");
        assert!(work[1] > work[0]);
    }
}

#[test]
fn g1_adaptive_six_profile_faces_retain_counterdiffusion_and_sensible_enthalpy() {
    let grid = CartesianGrid::new(4, 3, 2, 0.2, 0.3, 0.4).unwrap();
    let slope = [0.1, 0.2, 0.15];
    let fraction = |p: [f64; 3]| 0.2 + (0..3).map(|d| slope[d] * p[d]).sum::<f64>();
    let diffusion = 0.02;
    let duration = 0.1;
    let mut flow = field(
        equal_mass_gas(),
        grid,
        Boundaries3d {
            x: OPEN,
            y: OPEN,
            z: OPEN,
        },
        |p| {
            let a = fraction(p);
            ([a, 1.0 - a, 0.0], [0.0; 3], 1.0)
        },
    );
    let bytes = flow
        .required_reservoir_profile_workspace_bytes(&BoundarySide::ALL)
        .unwrap();
    flow.prepare_reservoir_profiles(
        &BoundarySide::ALL,
        bytes,
        |_, _, point, densities| {
            let a = fraction(point);
            densities.copy_from_slice(&[a, 1.0 - a, 0.0]);
            Ok(ReservoirFaceState {
                velocity: [0.0; 3],
                pressure_pa: 1.0,
            })
        },
        || false,
    )
    .unwrap();
    let original = flow.states().to_vec();
    let (mut operator, bytes) = accurate(&flow, law(flow.gas(), diffusion, 0.0, 0.0, 0.0), 0);
    let report = operator
        .advance_adaptive(&mut flow, Time::new(duration), accuracy(1e-7), || false)
        .unwrap();
    balanced(report);
    assert_eq!(report.stabilized.transport.accepted_steps, 1);
    assert_eq!(report.stabilized.transport.rejected_steps, 0);
    // The full/first-half anchor and qualified midpoint RHS are reused.
    assert_eq!(
        report.stabilized.operator_evaluations,
        3 * report.stabilized.maximum_stages_used + 1
    );
    assert_eq!(operator.workspace_bytes() + flow.workspace_bytes(), bytes);
    for (actual, expected) in flow.states().iter().zip(original) {
        close(*actual, expected, 3e-11);
    }
    for side in BoundarySide::ALL {
        let flux = -diffusion * slope[side.axis()];
        let scale = side.sign() * grid.face_area(side) * duration;
        for transfer in operator.boundary_impulses(side) {
            close(transfer.species_masses[0] / scale, flux, 2e-10);
            close(transfer.species_masses[1] / scale, -flux, 2e-10);
            assert_eq!(transfer.species_masses[2], 0.0);
            // Original Cp_A-Cp_B=-1 and the original Tmin=0.001 chart.
            close(transfer.energy / scale, -0.999 * flux, 3e-10);
        }
    }
}

fn channel() -> (MixtureEuler3d, MolecularTransport3d<ConstantTransport>) {
    let flow = field(
        equal_mass_gas(),
        CartesianGrid::new(3, 4, 2, 1.0 / 3.0, 0.25, 0.5).unwrap(),
        Boundaries3d {
            y: WALLS,
            ..PERIODIC
        },
        |_| ([1.0, 0.0, 0.0], [0.05, 0.0, -0.03], 1.0),
    );
    let (operator, _) = accurate(&flow, law(flow.gas(), 0.002, 0.03, 0.01, 0.04), 0);
    let operator = operator
        .with_wall(
            BoundarySide::YMax,
            MolecularWall3d::NoSlip {
                velocity: [0.2, 0.0, 0.1],
                thermal: WallThermalCondition::Isothermal(Temperature::new(2.0)),
            },
        )
        .unwrap();
    (flow, operator)
}

#[test]
fn g1_g4_adaptive_heated_wall_refinement_and_late_refusals_preserve_both_owners() {
    let controls = accuracy(1e-7);
    let duration = 0.2;
    let (mut reference, mut reference_operator) = channel();
    // Retain a real earlier publication before each failed complete interval.
    reference_operator
        .advance_adaptive(&mut reference, Time::new(0.02), controls, || false)
        .unwrap();
    let mut polls = 0usize;
    let expected = reference_operator
        .advance_adaptive(&mut reference, Time::new(duration), controls, || {
            polls += 1;
            false
        })
        .unwrap();
    balanced(expected);
    assert!(expected.stabilized.transport.boundary_change.energy > 0.0);
    assert!(expected.stabilized.transport.boundary_change.momentum[0] > 0.0);
    assert!(expected.stabilized.transport.boundary_change.momentum[2] > 0.0);
    assert!(expected.stabilized.operator_evaluations > 7 && polls > 10);

    let (mut small, mut small_operator) = channel();
    small_operator
        .advance_adaptive(&mut small, Time::new(0.02), controls, || false)
        .unwrap();
    let explicit = small_operator
        .advance(
            &mut small,
            Time::new(duration),
            DissipationControl {
                maximum_step_s: 0.0005,
                ..DissipationControl::default()
            },
            || false,
        )
        .unwrap();
    for i in 0..small.cell_count() {
        let a = small
            .cell(i)
            .unwrap()
            .primitive(small.gas(), || false)
            .unwrap();
        let b = reference
            .cell(i)
            .unwrap()
            .primitive(reference.gas(), || false)
            .unwrap();
        close(a.temperature, b.temperature, 2e-5);
        for axis in 0..3 {
            close(a.velocity[axis], b.velocity[axis], 2e-5);
        }
    }
    close(
        explicit.boundary_change.energy,
        expected.stabilized.transport.boundary_change.energy,
        2e-5,
    );

    for cancel in [false, true] {
        let (mut flow, mut operator) = channel();
        operator
            .advance_adaptive(&mut flow, Time::new(0.02), controls, || false)
            .unwrap();
        let state = flow.states().to_vec();
        let report = operator.last_step().copied();
        let faces = impulses(&operator);
        let walls = BoundarySide::ALL.map(|side| operator.wall(side));
        let bytes = operator.workspace_bytes();
        let mut attempt = controls;
        if !cancel {
            attempt.stabilized.maximum_operator_evaluations =
                expected.stabilized.operator_evaluations - 1;
        }
        let mut observed = 0usize;
        let refusal = operator
            .advance_adaptive(&mut flow, Time::new(duration), attempt, || {
                observed += 1;
                cancel && observed == polls
            })
            .unwrap_err();
        if cancel {
            assert!(matches!(refusal, DissipationError::Cancelled));
        } else {
            assert!(
                matches!(refusal, DissipationError::OperatorWorkLimit { operator_evaluations, .. }
                if operator_evaluations == expected.stabilized.operator_evaluations - 1)
            );
        }
        assert_eq!(flow.states(), state);
        assert_eq!(flow.time().value(), 0.0);
        assert!(flow.last_step().is_none());
        assert_eq!(operator.last_step().copied(), report);
        assert_eq!(impulses(&operator), faces);
        assert_eq!(BoundarySide::ALL.map(|side| operator.wall(side)), walls);
        assert_eq!(operator.workspace_bytes(), bytes);
        let replay = operator
            .advance_adaptive(&mut flow, Time::new(duration), controls, || false)
            .unwrap();
        assert_eq!(replay, expected);
        assert_eq!(flow.states(), reference.states());
        assert_eq!(impulses(&operator), impulses(&reference_operator));
    }
}

#[test]
fn g0_adaptive_storage_and_finite_inertia_refuse_before_publication() {
    let (mut flow, _) = channel();
    let make_source = || {
        SpatialTransport::new(
            flow.gas().clone(),
            law(flow.gas(), 0.002, 0.03, 0.01, 0.04),
            || false,
        )
        .unwrap()
    };
    let source = make_source();
    let bytes = MolecularTransport3d::required_adaptive_workspace_bytes(&flow, &source).unwrap();
    let mut limited = MolecularTransport3d::new(&flow, source, bytes - 1, || false).unwrap();
    let source = make_source();
    let initial = flow.states().to_vec();
    let before_bytes = limited.workspace_bytes();
    assert!(matches!(
        limited.advance_adaptive(&mut flow, Time::new(0.1), accuracy(1e-7), || false),
        Err(DissipationError::Invalid(_))
    ));
    assert_eq!(limited.workspace_bytes(), before_bytes);
    assert!(limited.last_step().is_none());
    assert_eq!(flow.states(), initial);

    let mut finite = MolecularTransport3d::new(
        &flow,
        source,
        bytes + WallInertiaControl3d::solver_workspace_bytes(),
        || false,
    )
    .unwrap()
    .with_wall(
        BoundarySide::YMax,
        MolecularWall3d::NoSlip {
            velocity: [0.2, 0.0, 0.1],
            thermal: WallThermalCondition::Adiabatic,
        },
    )
    .unwrap()
    .with_wall_inertia(
        &[InertialWall3d {
            side: BoundarySide::YMax,
            mass: Mass::new(1.0),
        }],
        WallInertiaControl3d::default(),
    )
    .unwrap();
    let wall = finite.wall(BoundarySide::YMax);
    assert!(matches!(
        finite.advance_adaptive(&mut flow, Time::new(0.1), accuracy(1e-7), || false),
        Err(DissipationError::Invalid(
            "stabilized molecular transport requires prescribed walls without finite inertia"
        ))
    ));
    assert_eq!(flow.states(), initial);
    assert_eq!(finite.wall(BoundarySide::YMax), wall);
    assert!(finite.last_step().is_none());
    finite
        .advance(
            &mut flow,
            Time::new(0.001),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    assert_ne!(finite.wall(BoundarySide::YMax), wall);
}
