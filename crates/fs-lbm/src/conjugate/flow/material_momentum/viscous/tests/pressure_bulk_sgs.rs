//! Independent molecular volume response with reservoir traction and SGS.

use super::*;
use crate::sgs3::{SgsModel, SgsUse};

fn model(model: SgsModel) -> FvSgsConfig {
    FvSgsConfig {
        model,
        filter_width_factor: 1.0,
        maximum_eddy_to_molecular: 10.0,
        usage: SgsUse::Exploration,
    }
}

fn longitudinal_profile() -> StaggeredBoundaryPressure<'static> {
    StaggeredBoundaryPressure {
        x_pa: &[Some(0.0), None, Some(0.0)],
        y_pa: &[None; 4],
        z_pa: &[None; 4],
    }
}

fn reservoirs<'a>(
    mass: &'a [f64; 2],
    h: &'a [f64; 2],
    p: &'a [[f64; 3]; 2],
    incoming_h: f64,
    incoming_x: [f64; 2],
) -> MaterialMomentumAdvanceInput<'a> {
    let mut data = input(mass, h, p);
    for (side, speed) in incoming_x.into_iter().enumerate() {
        data.material.boundary_inflow[side] = Some(AdvectedMaterial {
            density_kg_m3: 2.0,
            specific_enthalpy_j_kg: incoming_h,
        });
        data.boundary_velocity_m_s[side] = Some([speed, 0.0, 0.0]);
    }
    data
}

fn balance(report: &MaterialViscosityReport) {
    close(
        report.native_kinetic_energy_change_j
            + report.molecular_heat_j
            + report.sgs_transfer_j
            + report.backward_euler_loss_j
            - report.wall_work_j,
        0.0,
    );
    close(
        report.cell_kinetic_energy_change_j
            + report.molecular_heat_j
            + report.sgs_transfer_j
            + report.backward_euler_loss_j
            + report.numerical_scatter_loss_j
            - report.wall_work_j
            - report.representation_work_j,
        report.energy_defect_j,
    );
    close(report.energy_defect_j, 0.0);
}

#[test]
fn reservoir_normal_mode_preserves_authored_bulk_heat_and_actual_density_sgs() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [2.0; 2];
    let h = [10.0; 2];
    let p = [[2.0, 0.0, 0.0], [-2.0, 0.0, 0.0]];
    let cs: f64 = 0.2;
    let low_source = FluidTransportState::try_new(0.001, 2.0, 1.0, 1.0).unwrap();
    let high_source = FluidTransportState::try_new(50_000.0, 2.0, 1.0, 1.0).unwrap();
    for bulk in [0.0_f64, 0.7, 3.0] {
        // Native velocities [a,0,-a], with physical mouth masses 1 kg:
        // 2(a-1) + (4mu/3+zeta)*a + (4rho/3)*Cs^2*sqrt(2)*a^2 = 0.
        // Molecular bulk enters only the linear stiffness and physical heat.
        let longitudinal = 4.0 / 3.0 + bulk;
        let modeled = (8.0 / 3.0) * cs.powi(2) * 2.0_f64.sqrt();
        let linear = 2.0 + longitudinal;
        let speed = 4.0 / (linear + (linear.powi(2) + 8.0 * modeled).sqrt());
        let heat = longitudinal * speed.powi(2);
        let transfer = modeled * speed.powi(3);
        let heated_h = 10.0 + heat / 4.0;
        let data = reservoirs(&mass, &h, &p, heated_h, [0.0; 2]);
        let run = |source: &dyn FluidTransport| {
            field
                .advance_material_momentum_with_bulk_viscosity_sgs_and_pressure_profile(
                    &domain,
                    reservoir_faces(),
                    longitudinal_profile(),
                    data,
                    MaterialMomentumViscosityInput {
                        sources: &[Some(source); 2],
                        boundaries: reservoir_boundaries(),
                    },
                    &[bulk; 2],
                    model(SgsModel::Smagorinsky { cs }),
                    controls(),
                    viscous_controls(),
                    &CancelGate::new_clock_free(),
                )
                .unwrap()
        };
        let next = run(&low_source);
        let report = &next.viscosity;
        close(report.molecular_heat_j, heat);
        close(report.sgs_transfer_j, transfer);
        close(report.native_kinetic_energy_change_j, speed.powi(2) - 1.0);
        close(report.backward_euler_loss_j, (speed - 1.0).powi(2));
        close(report.numerical_scatter_loss_j, 0.5 * (speed - 1.0).powi(2));
        close(
            report.cell_kinetic_energy_change_j,
            0.5 * (1.0 + speed).powi(2) - 2.0,
        );
        close(report.representation_work_j, 0.0);
        close(report.wall_work_j, 0.0);
        close(next.represented_molecular_heat_j, heat);
        balance(report);
        for cell in 0..2 {
            close(report.cell_molecular_heat_j[cell], heat / 2.0);
            close(report.cell_sgs_transfer_j[cell], transfer / 2.0);
            close(next.advance.transport.mass_kg[cell], 2.0);
            close(
                next.advance.transport.specific_enthalpy_j_kg[cell],
                heated_h,
            );
            // The conservative feedback preserves original interior half-cell
            // counterflow even when its two half momenta cancel at the face.
            let sign = if cell == 0 { 1.0 } else { -1.0 };
            close(
                next.advance.momentum.momentum_kg_m_s[cell][0],
                sign * (1.0 + speed) / 2.0,
            );
        }
        let projection = &next.advance.transport.projection;
        for face in 0..3 {
            close(projection.field.flux_x(face, 0, 0), 0.0);
        }
        close(projection.report.open_boundary_pressure_work_j, 0.0);
        assert_eq!(run(&high_source), next);
    }
}

#[test]
fn independent_bulk_does_not_change_solenoidal_mouth_shear_or_heat_transport() {
    let domain = VoxelDomain::new(1, 2, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [2.0; 2];
    let h = [10.0; 2];
    let p = [[2.0, 0.0, 0.0], [-2.0, 0.0, 0.0]];
    let source = liquid(1.0);
    let sources: [Option<&dyn FluidTransport>; 2] = [Some(&source); 2];
    let cs: f64 = 0.2;
    let speed = 4.0 / (3.0 + (9.0 + 32.0 * cs.powi(2)).sqrt());
    let heat = 2.0 * speed.powi(2);
    let transfer = 8.0 * cs.powi(2) * speed.powi(3);
    let heated_h = 10.0 + heat / 4.0;
    let data = reservoirs(&mass, &h, &p, heated_h, [speed, -speed]);
    for bulk in [0.0, 0.7, 3.0] {
        let next = field
            .advance_material_momentum_with_bulk_viscosity_sgs_and_pressure_profile(
                &domain,
                reservoir_faces(),
                StaggeredBoundaryPressure {
                    x_pa: &[Some(0.0); 4],
                    y_pa: &[None; 3],
                    z_pa: &[None; 4],
                },
                data,
                MaterialMomentumViscosityInput {
                    sources: &sources,
                    boundaries: reservoir_boundaries(),
                },
                &[bulk; 2],
                model(SgsModel::Smagorinsky { cs }),
                controls(),
                viscous_controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        close(next.viscosity.molecular_heat_j, heat);
        close(next.viscosity.sgs_transfer_j, transfer);
        close(
            next.viscosity.native_kinetic_energy_change_j,
            2.0 * (speed.powi(2) - 1.0),
        );
        close(
            next.viscosity.backward_euler_loss_j,
            2.0 * (speed - 1.0).powi(2),
        );
        close(next.viscosity.numerical_scatter_loss_j, 0.0);
        balance(&next.viscosity);
        for cell in 0..2 {
            let sign = if cell == 0 { 1.0 } else { -1.0 };
            close(next.advance.transport.mass_kg[cell], 2.0);
            close(
                next.advance.transport.specific_enthalpy_j_kg[cell],
                heated_h,
            );
            close(
                next.advance.momentum.momentum_kg_m_s[cell][0],
                sign * 2.0 * speed,
            );
            for side in 0..2 {
                close(
                    next.advance
                        .transport
                        .projection
                        .field
                        .flux_x(side, cell, 0),
                    sign * speed,
                );
            }
        }
    }
}

#[test]
fn zero_closures_retain_bulk_pressure_path_and_bulk_cannot_hide_the_shear_cap() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [2.0; 2];
    let h = [10.0; 2];
    let p = [[2.0, 0.0, 0.0], [-2.0, 0.0, 0.0]];
    let source = liquid(1.0);
    let sources: [Option<&dyn FluidTransport>; 2] = [Some(&source); 2];
    let viscosity = MaterialMomentumViscosityInput {
        sources: &sources,
        boundaries: reservoir_boundaries(),
    };
    let data = reservoirs(&mass, &h, &p, 10.0, [0.0; 2]);
    let gate = CancelGate::new_clock_free();
    let bulk = [0.7; 2];
    let original = field
        .advance_material_momentum_with_bulk_viscosity_and_pressure_profile(
            &domain,
            reservoir_faces(),
            longitudinal_profile(),
            data,
            viscosity,
            &bulk,
            controls(),
            viscous_controls(),
            &gate,
        )
        .unwrap();
    let run = |bulk: &[f64], sgs| {
        field.advance_material_momentum_with_bulk_viscosity_sgs_and_pressure_profile(
            &domain,
            reservoir_faces(),
            longitudinal_profile(),
            data,
            viscosity,
            bulk,
            sgs,
            controls(),
            viscous_controls(),
            &gate,
        )
    };
    for selected in [
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
    ] {
        let next = run(
            &bulk,
            FvSgsConfig {
                maximum_eddy_to_molecular: 0.0,
                ..model(selected)
            },
        )
        .unwrap();
        assert_eq!(next, original);
        assert!(next.viscosity.cell_sgs_transfer_j.is_empty());
    }
    let active = model(SgsModel::Smagorinsky { cs: 0.2 });
    let accepted = run(&bulk, active).unwrap();
    for invalid in [&[0.7][..], &[-1.0, 0.7], &[0.7, f64::NAN]] {
        assert!(run(invalid, active).is_err());
        assert!(run(invalid, model(SgsModel::Smagorinsky { cs: 0.0 })).is_err());
    }
    assert!(matches!(
        run(
            &[1e6; 2],
            FvSgsConfig {
                maximum_eddy_to_molecular: 0.01,
                ..active
            }
        ),
        Err(ChtError::InvalidInput {
            field: "simple.sgs",
            ..
        })
    ));
    assert_eq!(run(&bulk, active).unwrap(), accepted);
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(mass, [2.0; 2]);
    assert_eq!(h, [10.0; 2]);
    assert_eq!(p, [[2.0, 0.0, 0.0], [-2.0, 0.0, 0.0]]);
}
