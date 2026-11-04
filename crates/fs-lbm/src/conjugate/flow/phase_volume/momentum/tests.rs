//! Independent phase-volume, pressure-impulse and material-momentum balances.

use super::*;
use crate::conjugate::flow::{
    CoexistenceVolumeControls, DensityForceControls, DensityMomentumReconstructionControls,
    DensityProjectionControls, MaterialAdvectionControls, ThermalVolumeControls,
};
use fs_blake3::ContentHash;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumPhaseState};

fn chart() -> EquilibriumEnthalpyPhaseCurve {
    // Original latent interval: v(h)=1/2+h/4 for 0<=h<=2.
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([181; 32]),
        [
            (-1.0, 300.0, 0.0, 2.0),
            (0.0, 400.0, 0.0, 2.0),
            (2.0, 400.0, 1.0, 1.0),
            (3.0, 500.0, 1.0, 1.0),
        ]
        .into_iter()
        .map(|(h, t, liquid, rho)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: liquid,
            bulk_density_kg_m3: rho,
        })
        .collect(),
    )
    .unwrap()
}

pub(super) fn controls() -> ThermalVolumeMomentumControls {
    let projection = DensityProjectionControls {
        duration_s: 0.1,
        relative_tolerance: 1e-12,
        maximum_volume_defect_m3_s: 1e-10,
        maximum_energy_defect_j: 1e-8,
    };
    ThermalVolumeMomentumControls {
        volume: PiecewiseThermalVolumeControls {
            thermal: ThermalVolumeControls {
                volume: CoexistenceVolumeControls {
                    transport: MaterialAdvectionControls {
                        projection,
                        max_cells: 100,
                        max_connections: 1000,
                        max_iterations: 1000,
                        relative_tolerance: 1e-12,
                        maximum_mass_defect_kg: 1e-10,
                        maximum_enthalpy_defect_j: 1e-9,
                    },
                    maximum_cell_volume_defect_m3: 1e-10,
                },
                max_volume_iterations: 20,
                source_relaxation: 1.0,
            },
            max_mixing_intervals: 1000,
        },
        mechanics: MaterialDensityMomentumControls {
            reconstruction: DensityMomentumReconstructionControls {
                max_cells: 100,
                maximum_mass_defect_kg: 1e-10,
                maximum_momentum_defect_kg_m_s: 1e-9,
                maximum_energy_defect_j: 1e-8,
            },
            forcing: DensityForceControls {
                max_cells: 100,
                maximum_impulse_defect_kg_m_s: 1e-9,
                maximum_energy_defect_j: 1e-8,
            },
            projection,
            maximum_feedback_momentum_defect_kg_m_s: 1e-9,
            maximum_feedback_angular_impulse_defect_kg_m2_s: 1e-9,
            maximum_feedback_energy_defect_j: 1e-8,
        },
        momentum: MaterialMomentumControls {
            max_iterations: 1000,
            relative_tolerance: 1e-12,
            maximum_momentum_defect_kg_m_s: 1e-9,
            maximum_energy_defect_j: 1e-8,
        },
    }
}

fn input<'a>(
    mass: &'a [f64],
    states: &'a [Option<EquilibriumPhaseState>],
    increments: &'a [f64],
    momentum: &'a [[f64; 3]],
    force: &'a [[f64; 3]],
) -> ThermalVolumeMomentumInput<'a> {
    ThermalVolumeMomentumInput {
        volume: CoexistenceVolumeInput {
            mass_kg: mass,
            phase_states: states,
            net_enthalpy_increment_j: increments,
            boundary_inflow: [None; 6],
        },
        momentum_kg_m_s: momentum,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: force,
        boundary_velocity_m_s: [None; 6],
    }
}

pub(super) fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 4e-10 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn latent_expansion_and_contraction_drive_one_conservative_momentum_step_on_every_axis() {
    let curve = chart();
    let mass = [2.0, 1.0];
    let states = [
        Some(curve.state_at_specific_enthalpy(0.0).unwrap()),
        Some(curve.state_at_specific_enthalpy(2.0).unwrap()),
    ];
    let increments = [0.2, -0.2];
    let zero = [[0.0; 3]; 2];
    for axis in 0..3 {
        let mut shape = [1; 3];
        shape[axis] = 2;
        let domain = VoxelDomain::new(shape[0], shape[1], shape[2], 1.0).unwrap();
        let field = FlowField::quiescent(&domain);
        let result = field
            .advance_piecewise_thermal_volume_momentum(
                &domain,
                [FlowFace::Wall; 6],
                &curve,
                input(&mass, &states, &increments, &zero, &zero),
                controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        let volume = &result.thermal.volume;
        let transport = &volume.transport;
        // delta_V=(1/4)*delta_H=+/-1/20, so Q01=1/2. Backward
        // donor balance gives M0'=2/(1+1/20), M1'=3-M0'.
        close(volume.volume_source_m3_s[0], 0.5);
        close(volume.volume_source_m3_s[1], -0.5);
        close(
            transport.projection.field.outward(&domain, 0, 2 * axis + 1),
            0.5,
        );
        close(transport.mass_kg[0], 40.0 / 21.0);
        close(transport.mass_kg[1], 23.0 / 21.0);
        close(transport.specific_enthalpy_j_kg[0], 0.1);
        close(transport.specific_enthalpy_j_kg[1], 38.0 / 23.0);
        close(transport.mass_kg.iter().sum(), 3.0);
        close(
            transport
                .mass_kg
                .iter()
                .zip(&transport.specific_enthalpy_j_kg)
                .map(|(m, h)| m * h)
                .sum(),
            2.0,
        );
        close(transport.report.incoming_mass_kg, 0.0);
        close(transport.report.outgoing_mass_kg, 0.0);
        // Actual dual mass is (2+1)/2. J=Mdual*Q=3/4, hence dp=7.5.
        // Its cell shares produce common velocity 1/4 before advection.
        close(
            transport.projection.pressure_correction_pa[0]
                - transport.projection.pressure_correction_pa[1],
            7.5,
        );
        for component in 0..3 {
            let expected = if component == axis { 0.75 } else { 0.0 };
            close(
                result.force_projection.pressure_impulse_kg_m_s[component],
                expected,
            );
            close(
                result.feedback.cell_momentum_change_kg_m_s[component],
                expected,
            );
            close(
                result.momentum.report.momentum_after_kg_m_s[component],
                expected,
            );
        }
        close(result.momentum.momentum_kg_m_s[0][axis], 10.0 / 21.0);
        close(result.momentum.momentum_kg_m_s[1][axis], 23.0 / 84.0);
        for (cell, phase) in volume.phase_states.iter().enumerate() {
            let phase = phase.unwrap();
            close(transport.mass_kg[cell] / phase.bulk_density_kg_m3(), 1.0);
            close(result.momentum.velocity_m_s[cell][axis], 0.25);
            close(phase.temperature_k(), 400.0);
            assert_eq!(phase.phase_curve_identity(), curve.identity());
        }
        // Source work=3/8, split into physical returned K and separately
        // reported projection/scatter losses. No part is silently added to h.
        close(volume.report.required_source_pressure_work_j, 3.0 / 8.0);
        close(result.momentum.report.kinetic_energy_after_j, 3.0 / 32.0);
        close(result.feedback.numerical_scatter_loss_j, 3.0 / 32.0);
        close(transport.projection.report.projection_loss_j, 3.0 / 16.0);
        close(result.momentum.report.numerical_transport_loss_j, 0.0);
        close(volume.report.material_enthalpy_change_j, 0.0);
        assert_eq!(result.momentum.report.represented_impulse_kg_m_s, [0.0; 3]);
        assert_eq!(field, FlowField::quiescent(&domain));
    }
}

#[test]
fn continuation_reconstructs_returned_mass_and_momentum_not_old_transport_velocity() {
    let curve = chart();
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let faces = [FlowFace::Wall; 6];
    let field = FlowField::quiescent(&domain);
    let mass = [2.0, 1.0];
    let states = [
        Some(curve.state_at_specific_enthalpy(0.0).unwrap()),
        Some(curve.state_at_specific_enthalpy(2.0).unwrap()),
    ];
    let zero = [[0.0; 3]; 2];
    let first = field
        .advance_piecewise_thermal_volume_momentum(
            &domain,
            faces,
            &curve,
            input(&mass, &states, &[0.2, -0.2], &zero, &zero),
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let volume = &first.thermal.volume;
    let next = input(
        &volume.transport.mass_kg,
        &volume.phase_states,
        &[0.0; 2],
        &first.momentum.momentum_kg_m_s,
        &zero,
    );
    let continued = volume
        .transport
        .projection
        .field
        .advance_piecewise_thermal_volume_momentum(
            &domain,
            faces,
            &curve,
            next,
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    // Returned cell speed is 1/4, while the old interior transport speed is
    // 1/2. With zero thermal source in a sealed line, the new pressure and
    // held-wall reactions stop exactly the actual supplied cell momentum.
    close(
        continued
            .thermal
            .volume
            .transport
            .projection
            .field
            .flux_x(1, 0, 0),
        0.0,
    );
    close(continued.feedback.cell_momentum_change_kg_m_s[0], -0.75);
    for cell in 0..2 {
        close(
            continued.thermal.volume.transport.mass_kg[cell],
            volume.transport.mass_kg[cell],
        );
        close(
            continued.thermal.volume.transport.specific_enthalpy_j_kg[cell],
            volume.transport.specific_enthalpy_j_kg[cell],
        );
        for component in 0..3 {
            close(continued.momentum.momentum_kg_m_s[cell][component], 0.0);
        }
    }
    let mut unrelated = FlowField::quiescent(&domain);
    *unrelated.face_mut(0, 1) = 17.0;
    let same = unrelated
        .advance_piecewise_thermal_volume_momentum(
            &domain,
            faces,
            &curve,
            next,
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(same, continued);
}

#[test]
fn contraction_pressure_gauge_and_late_incoming_momentum_refusal_preserve_retry() {
    let curve = chart();
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let original = field.clone();
    let mass = [4.0 / 3.0];
    let states = [Some(curve.state_at_specific_enthalpy(1.0).unwrap())];
    let increments = [-0.2];
    let zero = [[0.0; 3]];
    let mut data = input(&mass, &states, &increments, &zero, &zero);
    data.volume.boundary_inflow[1] = states[0];
    data.boundary_velocity_m_s[1] = Some([-0.5, 0.0, 0.0]);
    let mut faces = [FlowFace::Wall; 6];
    faces[1] = FlowFace::Free;
    for reference in [0.0, 1e12] {
        let x = [None, Some(reference)];
        let transverse = [None; 2];
        let profile = StaggeredBoundaryPressure {
            x_pa: &x,
            y_pa: &transverse,
            z_pa: &transverse,
        };
        let run = |data, control, gate: &CancelGate| {
            field.advance_piecewise_thermal_volume_momentum_with_pressure_profile(
                &domain, faces, profile, &curve, data, control, gate,
            )
        };
        let expected = run(data, controls(), &CancelGate::new_clock_free()).unwrap();
        let volume = &expected.thermal.volume;
        // delta_V=-1/20; inflow mass=(1/20)*(4/3)=1/15. Independent
        // M/Mh balances give M'=7/5 and h'=6/7, again filling one m3.
        close(volume.volume_source_m3_s[0], -0.5);
        close(volume.transport.mass_kg[0], 7.0 / 5.0);
        close(volume.transport.specific_enthalpy_j_kg[0], 6.0 / 7.0);
        close(volume.transport.report.incoming_mass_kg, 1.0 / 15.0);
        close(
            volume.transport.mass_kg[0] / volume.phase_states[0].unwrap().bulk_density_kg_m3(),
            1.0,
        );
        // One open half-cell contributes J=(M/2)*(-1/2)=-1/3.
        // Incoming momentum is (1/15)*(-1/2); no second pressure kick.
        close(expected.feedback.cell_momentum_change_kg_m_s[0], -1.0 / 3.0);
        close(expected.momentum.momentum_kg_m_s[0][0], -11.0 / 30.0);
        close(expected.momentum.velocity_m_s[0][0], -11.0 / 42.0);
        close(volume.report.required_source_pressure_work_j, 1.0 / 6.0);
        close(expected.feedback.composed_energy_defect_j, 0.0);

        let mut missing = data;
        missing.boundary_velocity_m_s = [None; 6];
        // This refusal occurs after thermal-volume, pressure and cell-feedback
        // acceptance, when the original accepted mass route needs momentum.
        assert!(matches!(
            run(missing, controls(), &CancelGate::new_clock_free()),
            Err(ChtError::InvalidInput { field: "flow.material_advection", reason })
                if reason.contains("incoming mass route lacks an authored velocity")
        ));
        let mut wrong_clock = controls();
        wrong_clock.mechanics.projection.duration_s *= 2.0;
        assert!(run(data, wrong_clock, &CancelGate::new_clock_free()).is_err());
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(matches!(
            run(data, controls(), &cancelled),
            Err(ChtError::Cancelled)
        ));
        assert_eq!(field, original);
        assert_eq!(mass, [4.0 / 3.0]);
        assert_eq!(
            states,
            [Some(curve.state_at_specific_enthalpy(1.0).unwrap())]
        );
        assert_eq!(increments, [-0.2]);
        assert_eq!(zero, [[0.0; 3]]);
        assert_eq!(
            run(data, controls(), &CancelGate::new_clock_free()).unwrap(),
            expected
        );
    }
}
