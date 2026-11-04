//! Molecular heat must produce the original chart's occupied-volume response.

use super::super::tests::{close, controls as momentum_controls};
use super::*;
use crate::conjugate::FvBoundary;
use crate::conjugate::flow::{
    CoexistenceVolumeInput, MaterialAdvectionInput, MaterialMomentumAdvanceControls,
    MaterialMomentumAdvanceInput, MaterialVolumeSource,
};
use fs_blake3::ContentHash;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};

fn controls() -> ThermalVolumeMomentumControls {
    let mut value = momentum_controls();
    value.mechanics.projection.duration_s = 0.5;
    value.volume.thermal.volume.transport.projection = value.mechanics.projection;
    value
}

fn viscosity_controls() -> MaterialViscosityControls {
    MaterialViscosityControls {
        max_cells: 100,
        max_sweeps: 100,
        relative_tolerance: 1e-12,
        maximum_impulse_defect_kg_m_s: 1e-9,
        maximum_angular_impulse_defect_kg_m2_s: 1e-9,
        maximum_energy_defect_j: 1e-8,
    }
}

fn liquid_curve(changing_density: bool) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([197; 32]),
        SolidLiquidPhase::Liquid,
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 300.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 2.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 2.0,
                temperature_k: 302.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: if changing_density { 1.0 } else { 2.0 },
            },
        ],
    )
    .unwrap()
}

fn mouths(axis: usize) -> ([FlowFace; 6], [FvBoundary; 6]) {
    let mut faces = [FlowFace::Wall; 6];
    let mut boundaries = [FvBoundary::Symmetry; 6];
    for side in 0..2 {
        faces[2 * axis + side] = FlowFace::Free;
        boundaries[2 * axis + side] = FvBoundary::Outlet;
        boundaries[2 * ((axis + 1) % 3) + side] = FvBoundary::Wall { velocity: [0.0; 3] };
    }
    (faces, boundaries)
}

#[test]
fn molecular_heat_drives_exact_density_volume_and_momentum_on_every_axis() {
    let curve = liquid_curve(true);
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let state = curve.state_at_specific_enthalpy(0.0).unwrap();
    let mass = [2.0];
    let phases = [Some(state)];
    let source = FluidTransportState::try_new(2.0, 1.0, 1.0, 1.0).unwrap();
    let sources: [Option<&dyn FluidTransport>; 1] = [Some(&source)];
    for axis in 0..3 {
        let (faces, boundaries) = mouths(axis);
        let mut momentum = [[0.0; 3]];
        momentum[0][axis] = 2.0;
        let mut incoming_velocity = [0.0; 3];
        incoming_velocity[axis] = 0.5;
        let mut input = ThermalVolumeMomentumInput {
            volume: CoexistenceVolumeInput {
                mass_kg: &mass,
                phase_states: &phases,
                net_enthalpy_increment_j: &[0.0],
                boundary_inflow: [None; 6],
            },
            momentum_kg_m_s: &momentum,
            acceleration_m_s2: [0.0; 3],
            force_density_n_m3: &[[0.0; 3]],
            boundary_velocity_m_s: [None; 6],
        };
        input.volume.boundary_inflow[2 * axis] = Some(state);
        input.boundary_velocity_m_s[2 * axis] = Some(incoming_velocity);
        for reference in [0.0, 1e12] {
            let mut planes = [[None; 2]; 3];
            planes[axis] = [Some(reference); 2];
            let result = field
                .advance_piecewise_thermal_volume_momentum_with_viscosity_and_pressure_profile(
                    &domain,
                    faces,
                    StaggeredBoundaryPressure {
                        x_pa: &planes[0],
                        y_pa: &planes[1],
                        z_pa: &planes[2],
                    },
                    &curve,
                    input,
                    MaterialMomentumViscosityInput {
                        sources: &sources,
                        boundaries,
                    },
                    controls(),
                    viscosity_controls(),
                    &CancelGate::new_clock_free(),
                )
                .unwrap();
            // M=2, u=1, and the two stationary transverse walls give
            // 2*(u*-1)/dt=-4*mu*u*. With dt=1/2 and mu=1, u*=1/2.
            // Physical heat=dt*4*mu*u*^2=1/2; BE loss=M*(du)^2/2=1/4.
            close(
                result.viscosity_reconstruction.cell_momentum_kg_m_s[axis],
                2.0,
            );
            close(result.viscosity.cell_momentum_change_kg_m_s[axis], -1.0);
            close(result.viscosity.cell_kinetic_energy_change_j, -0.75);
            close(result.viscosity.molecular_heat_j, 0.5);
            close(result.viscosity.cell_molecular_heat_j[0], 0.5);
            close(result.viscosity.backward_euler_loss_j, 0.25);
            close(result.viscosity.numerical_scatter_loss_j, 0.0);
            close(result.viscosity.wall_work_j, 0.0);
            close(result.viscosity.energy_defect_j, 0.0);
            let evolution = &result.evolution;
            close(evolution.reconstruction.cell_momentum_kg_m_s[axis], 1.0);
            let volume = &evolution.thermal.volume;
            let transport = &volume.transport;
            // v(h)=1/2+h/4. Q_mu=1/2 changes the initial occupied volume
            // by 1/8, hence S=1/4. Equal reservoir heads preserve mean u*=1/2,
            // giving face fluxes [3/8,5/8] and zero net pressure impulse.
            close(volume.report.requested_enthalpy_increment_j, 0.5);
            close(volume.volume_source_m3_s[0], 0.25);
            close(
                transport.projection.field.outward(&domain, 0, 2 * axis),
                -3.0 / 8.0,
            );
            close(
                transport.projection.field.outward(&domain, 0, 2 * axis + 1),
                5.0 / 8.0,
            );
            close(
                evolution.force_projection.pressure_impulse_kg_m_s[axis],
                0.0,
            );
            close(evolution.feedback.cell_momentum_change_kg_m_s[axis], 0.0);
            close(volume.report.required_source_pressure_work_j, 1.0 / 32.0);
            close(transport.projection.report.projection_loss_j, 1.0 / 64.0);
            close(evolution.feedback.numerical_scatter_loss_j, 1.0 / 64.0);
            // Inlet mass=3/8 and implicit outlet mass=(5/16)*M_new.
            // M_new=(2+3/8)/(1+5/16)=38/21; H_new=(1/2)/(1+5/16).
            close(transport.mass_kg[0], 38.0 / 21.0);
            close(transport.specific_enthalpy_j_kg[0], 4.0 / 19.0);
            close(transport.report.incoming_mass_kg, 3.0 / 8.0);
            close(transport.report.outgoing_mass_kg, 95.0 / 168.0);
            close(
                transport.mass_kg[0] / volume.phase_states[0].unwrap().bulk_density_kg_m3(),
                1.0,
            );
            close(volume.report.material_enthalpy_change_j, 8.0 / 21.0);
            close(evolution.momentum.momentum_kg_m_s[0][axis], 19.0 / 21.0);
            close(evolution.momentum.velocity_m_s[0][axis], 0.5);
            close(evolution.momentum.report.numerical_transport_loss_j, 0.0);
            assert_eq!(
                evolution.momentum.report.represented_impulse_kg_m_s,
                [0.0; 3]
            );
            close(result.absolute_heat_assembly_defect_j, 0.0);
            close(result.absolute_enthalpy_defect_j, 0.0);
        }
    }
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(mass, [2.0]);
    assert_eq!(phases, [Some(state)]);
}

#[test]
fn source_density_and_late_inflow_refusals_preserve_the_original_state_and_retry() {
    let curve = liquid_curve(true);
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let state = curve.state_at_specific_enthalpy(0.0).unwrap();
    let mass = [2.0];
    let phases = [Some(state)];
    let momentum = [[2.0, 0.0, 0.0]];
    let (faces, boundaries) = mouths(0);
    let source = FluidTransportState::try_new(2.0, 1.0, 1.0, 1.0).unwrap();
    let sources: [Option<&dyn FluidTransport>; 1] = [Some(&source)];
    let viscosity = MaterialMomentumViscosityInput {
        sources: &sources,
        boundaries,
    };
    let mut input = ThermalVolumeMomentumInput {
        volume: CoexistenceVolumeInput {
            mass_kg: &mass,
            phase_states: &phases,
            net_enthalpy_increment_j: &[0.0],
            boundary_inflow: [None; 6],
        },
        momentum_kg_m_s: &momentum,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: &[[0.0; 3]],
        boundary_velocity_m_s: [None; 6],
    };
    input.volume.boundary_inflow[0] = Some(state);
    input.boundary_velocity_m_s[0] = Some([0.5, 0.0, 0.0]);
    let x = [Some(1e12); 2];
    let other = [None; 2];
    let run = |input, viscosity, gate: &CancelGate| {
        field.advance_piecewise_thermal_volume_momentum_with_viscosity_and_pressure_profile(
            &domain,
            faces,
            StaggeredBoundaryPressure {
                x_pa: &x,
                y_pa: &other,
                z_pa: &other,
            },
            &curve,
            input,
            viscosity,
            controls(),
            viscosity_controls(),
            gate,
        )
    };
    let expected = run(input, viscosity, &CancelGate::new_clock_free()).unwrap();
    // This API owns a chart-changing-density evolution: even one ULP of
    // source/chart rho disagreement refuses before molecular heat is used.
    let wrong =
        FluidTransportState::try_new(f64::from_bits(2.0_f64.to_bits() + 1), 1.0, 1.0, 1.0).unwrap();
    let wrong_sources: [Option<&dyn FluidTransport>; 1] = [Some(&wrong)];
    assert!(matches!(
        run(input, MaterialMomentumViscosityInput { sources: &wrong_sources, ..viscosity }, &CancelGate::new_clock_free()),
        Err(ChtError::InvalidInput { reason, .. }) if reason.contains("density")
    ));
    let mut missing_velocity = input;
    missing_velocity.boundary_velocity_m_s = [None; 6];
    assert!(matches!(
        run(missing_velocity, viscosity, &CancelGate::new_clock_free()),
        Err(ChtError::InvalidInput { reason, .. }) if reason.contains("incoming mass route lacks an authored velocity")
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        run(input, viscosity, &cancelled),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(mass, [2.0]);
    assert_eq!(phases, [Some(state)]);
    assert_eq!(momentum, [[2.0, 0.0, 0.0]]);
    assert_eq!(
        run(input, viscosity, &CancelGate::new_clock_free()).unwrap(),
        expected
    );
}

#[test]
fn constant_density_limit_retains_the_original_closed_viscous_material_evolution() {
    let curve = liquid_curve(false);
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let phases = [Some(curve.state_at_specific_enthalpy(0.0).unwrap()); 2];
    let mass = [2.0; 2];
    let momentum = [[2.0, 0.0, 0.0], [0.0; 3]];
    let source = FluidTransportState::try_new(2.0, 1.0, 1.0, 1.0).unwrap();
    let sources: [Option<&dyn FluidTransport>; 2] = [Some(&source); 2];
    let viscosity = MaterialMomentumViscosityInput {
        sources: &sources,
        boundaries: [FvBoundary::Wall { velocity: [0.0; 3] }; 6],
    };
    let input = ThermalVolumeMomentumInput {
        volume: CoexistenceVolumeInput {
            mass_kg: &mass,
            phase_states: &phases,
            net_enthalpy_increment_j: &[0.0; 2],
            boundary_inflow: [None; 6],
        },
        momentum_kg_m_s: &momentum,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: &[[0.0; 3]; 2],
        boundary_velocity_m_s: [None; 6],
    };
    let controls = controls();
    let gate = CancelGate::new_clock_free();
    let coupled = field
        .advance_piecewise_thermal_volume_momentum_with_viscosity(
            &domain,
            &curve,
            input,
            viscosity,
            controls,
            viscosity_controls(),
            &gate,
        )
        .unwrap();
    let original = field
        .advance_material_momentum_with_viscosity(
            &domain,
            MaterialMomentumAdvanceInput {
                material: MaterialAdvectionInput {
                    mass_kg: &mass,
                    specific_enthalpy_j_kg: &[0.0; 2],
                    volume_sources: &[MaterialVolumeSource::None; 2],
                    boundary_inflow: [None; 6],
                },
                momentum_kg_m_s: &momentum,
                acceleration_m_s2: [0.0; 3],
                force_density_n_m3: &[[0.0; 3]; 2],
                source_velocity_m_s: &[None; 2],
                boundary_velocity_m_s: [None; 6],
            },
            viscosity,
            MaterialMomentumAdvanceControls {
                mechanics: controls.mechanics,
                material: controls.volume.thermal.volume.transport,
                momentum: controls.momentum,
            },
            viscosity_controls(),
            &gate,
        )
        .unwrap();
    assert!(coupled.viscosity.molecular_heat_j > 0.0);
    assert_eq!(coupled.viscosity, original.viscosity);
    assert_eq!(
        coupled.viscosity_reconstruction,
        original.viscosity_reconstruction
    );
    assert_eq!(coupled.evolution.momentum, original.advance.momentum);
    assert_eq!(
        coupled.evolution.thermal.volume.transport,
        original.advance.transport
    );
    assert_eq!(
        coupled.evolution.thermal.volume.volume_source_m3_s,
        [0.0; 2]
    );
    close(
        coupled
            .evolution
            .thermal
            .volume
            .report
            .material_enthalpy_change_j,
        coupled.viscosity.molecular_heat_j,
    );
    assert_eq!(field, FlowField::quiescent(&domain));
}
