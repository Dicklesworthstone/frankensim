//! Reservoir half-cell inertia, boundary SGS shear and original pressure work.

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

fn pressure(x_pa: &[Option<f64>; 4]) -> StaggeredBoundaryPressure<'_> {
    StaggeredBoundaryPressure {
        x_pa,
        y_pa: &[None; 3],
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

/// Two transverse cells, each one cell deep between the x reservoirs, carry
/// initial u_x=[+1,-1]. At either physical mouth half-edge gamma=-2a,
/// Delta=1 and nu_t=2 Cs^2 a. Each normal face has physical mass 1 kg.
/// Its independent BE equation is 2(a-1)+(1+4 Cs^2 a)*a=0.
fn mouth_speed(cs: f64) -> f64 {
    4.0 / (3.0 + (9.0 + 32.0 * cs.powi(2)).sqrt())
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
fn opposed_mouth_shear_uses_physical_half_inertia_and_separate_sgs_transfer() {
    let domain = VoxelDomain::new(1, 2, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [2.0; 2];
    let h = [10.0; 2];
    let p = [[2.0, 0.0, 0.0], [-2.0, 0.0, 0.0]];
    let source = FluidTransportState::try_new(999.0, 2.0, 1.0, 1.0).unwrap();
    let sources: [Option<&dyn FluidTransport>; 2] = [Some(&source); 2];
    let cs: f64 = 0.2;
    let speed = mouth_speed(cs);
    let heat = 2.0 * speed.powi(2);
    let transfer = 8.0 * cs.powi(2) * speed.powi(3);
    let heated_h = 10.0 + heat / 4.0;
    let data = reservoirs(&mass, &h, &p, heated_h, [speed, -speed]);
    let run = |sources| {
        field
            .advance_material_momentum_with_viscosity_sgs_and_pressure_profile(
                &domain,
                reservoir_faces(),
                pressure(&[Some(0.0); 4]),
                data,
                MaterialMomentumViscosityInput {
                    sources,
                    boundaries: reservoir_boundaries(),
                },
                model(SgsModel::Smagorinsky { cs }),
                controls(),
                viscous_controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap()
    };
    let next = run(&sources);
    let report = &next.viscosity;
    // There are no interior four-cell shear edges. This loss comes entirely
    // from the two reservoir half-edges, whose native gradient is twice the
    // cell-centred gradient. Averaging cell nu_t cannot satisfy this root.
    close(
        report.native_kinetic_energy_change_j,
        2.0 * (speed.powi(2) - 1.0),
    );
    close(report.backward_euler_loss_j, 2.0 * (speed - 1.0).powi(2));
    close(report.molecular_heat_j, heat);
    close(report.sgs_transfer_j, transfer);
    close(report.numerical_scatter_loss_j, 0.0);
    close(report.representation_work_j, 0.0);
    close(report.wall_work_j, 0.0);
    close(report.face_angular_impulse_kg_m2_s[2], 2.0 * (speed - 1.0));
    close(next.represented_molecular_heat_j, heat);
    balance(report);
    for cell in 0..2 {
        let sign = if cell == 0 { 1.0 } else { -1.0 };
        close(report.cell_molecular_heat_j[cell], heat / 2.0);
        close(report.cell_sgs_transfer_j[cell], transfer / 2.0);
        close(next.advance.transport.mass_kg[cell], mass[cell]);
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
    close(next.advance.transport.projection.field.flux_y(0, 1, 0), 0.0);
    close(next.advance.transport.report.incoming_mass_kg, 2.0 * speed);
    close(next.advance.transport.report.outgoing_mass_kg, 2.0 * speed);
    close(
        next.advance
            .transport
            .projection
            .report
            .open_boundary_pressure_work_j,
        0.0,
    );
    close(
        next.advance.transport.projection.report.projection_loss_j,
        0.0,
    );
    // This is a bidirectional reservoir flow; the same accepted material
    // routes preserve the authored incoming specific h and normal velocity.
    // Source reference rho changes neither mouth inertia nor dynamic SGS.
    let other = FluidTransportState::try_new(0.001, 2.0, 1.0, 1.0).unwrap();
    let other_sources: [Option<&dyn FluidTransport>; 2] = [Some(&other); 2];
    assert_eq!(run(&other_sources), next);
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(mass, [2.0; 2]);
    assert_eq!(h, [10.0; 2]);
    assert_eq!(p, [[2.0, 0.0, 0.0], [-2.0, 0.0, 0.0]]);
}

#[test]
fn spatial_reservoir_pressure_is_applied_once_after_the_sgs_predictor() {
    let domain = VoxelDomain::new(1, 2, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [2.0; 2];
    let h = [10.0; 2];
    let p = [[2.0, 0.0, 0.0], [-2.0, 0.0, 0.0]];
    let source = liquid(1.0);
    let sources: [Option<&dyn FluidTransport>; 2] = [Some(&source); 2];
    let speed = mouth_speed(0.2);
    let heated_h = 10.0 + 0.5 * speed.powi(2);
    let data = reservoirs(&mass, &h, &p, heated_h, [0.0; 2]);
    let run = |profile: StaggeredBoundaryPressure<'_>| {
        field
            .advance_material_momentum_with_viscosity_sgs_and_pressure_profile(
                &domain,
                reservoir_faces(),
                profile,
                data,
                MaterialMomentumViscosityInput {
                    sources: &sources,
                    boundaries: reservoir_boundaries(),
                },
                model(SgsModel::Smagorinsky { cs: 0.2 }),
                controls(),
                viscous_controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap()
    };
    let unforced = run(pressure(&[Some(0.0); 4]));
    for reference in [0.0, 1e12] {
        let x = [
            Some(reference + 2.0),
            Some(reference),
            Some(reference + 6.0),
            Some(reference + 2.0),
        ];
        let next = run(pressure(&x));
        assert_eq!(next.viscosity, unforced.viscosity);
        close(next.represented_molecular_heat_j, 2.0 * speed.powi(2));
        // The physical density projection has two pressure rows:
        // 5p0-p1=4 and -p0+5p1=16, hence p0=1.5, p1=3.5 Pa
        // relative to the datum. Normal mouth masses are 1, interior y mass
        // is 2. This independently fixes all five nonzero face velocities.
        let projection = &next.advance.transport.projection;
        for (cell, expected) in [[speed + 0.25, speed + 0.75], [-speed + 1.25, -speed + 0.75]]
            .into_iter()
            .enumerate()
        {
            for (side, value) in expected.into_iter().enumerate() {
                close(projection.field.flux_x(side, cell, 0), value);
            }
        }
        close(projection.field.flux_y(0, 1, 0), -0.5);
        close(
            next.advance.force_projection.pressure_impulse_kg_m_s[0],
            3.0,
        );
        close(
            projection.report.open_boundary_pressure_work_j,
            speed - 3.25,
        );
        close(projection.report.kinetic_energy_change_j, 1.625 - speed);
        close(projection.report.projection_loss_j, 1.625);
        close(projection.report.energy_defect_j, 0.0);
        for cell in 0..2 {
            close(next.advance.transport.mass_kg[cell], 2.0);
            close(
                next.advance.transport.specific_enthalpy_j_kg[cell],
                heated_h,
            );
        }
    }
}

#[test]
fn every_zero_model_retains_the_exact_molecular_reservoir_advance() {
    let domain = VoxelDomain::new(1, 2, 1, 1.0).unwrap();
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
    let x = [Some(2.0), Some(0.0), Some(6.0), Some(2.0)];
    let gate = CancelGate::new_clock_free();
    let original = field
        .advance_material_momentum_with_viscosity_and_pressure_profile(
            &domain,
            reservoir_faces(),
            pressure(&x),
            data,
            viscosity,
            controls(),
            viscous_controls(),
            &gate,
        )
        .unwrap();
    for selected in [
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
    ] {
        let next = field
            .advance_material_momentum_with_viscosity_sgs_and_pressure_profile(
                &domain,
                reservoir_faces(),
                pressure(&x),
                data,
                viscosity,
                FvSgsConfig {
                    maximum_eddy_to_molecular: 0.0,
                    ..model(selected)
                },
                controls(),
                viscous_controls(),
                &gate,
            )
            .unwrap();
        assert_eq!(next, original);
        assert!(next.viscosity.cell_sgs_transfer_j.is_empty());
        assert_eq!(next.viscosity.sgs_transfer_j, 0.0);
    }
}

#[test]
fn mouth_edge_cap_and_late_refusals_leave_original_material_replayable() {
    let domain = VoxelDomain::new(1, 2, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [2.0; 2];
    let h = [10.0; 2];
    let p = [[2.0, 0.0, 0.0], [-2.0, 0.0, 0.0]];
    let source = FluidTransportState::try_new(0.001, 2.0, 1.0, 1.0).unwrap();
    let sources: [Option<&dyn FluidTransport>; 2] = [Some(&source); 2];
    let data = reservoirs(&mass, &h, &p, 10.0, [0.0; 2]);
    let selected = model(SgsModel::Smagorinsky { cs: 0.2 });
    let x = [Some(0.0); 4];
    let run = |data, boundaries, profile, selected, gate: &CancelGate| {
        field.advance_material_momentum_with_viscosity_sgs_and_pressure_profile(
            &domain,
            reservoir_faces(),
            profile,
            data,
            MaterialMomentumViscosityInput {
                sources: &sources,
                boundaries,
            },
            selected,
            controls(),
            viscous_controls(),
            gate,
        )
    };
    let gate = CancelGate::new_clock_free();
    let accepted = run(data, reservoir_boundaries(), pressure(&x), selected, &gate).unwrap();
    // Initial cell nu_t=.04 passes a .1 molecular ratio cap at actual rho2.
    // Mouth nu_t=.08 gives ratio .16 and must refuse before any solve; tiny
    // source reference rho must not hide this native boundary-edge violation.
    assert!(
        run(
            data,
            reservoir_boundaries(),
            pressure(&x),
            FvSgsConfig {
                maximum_eddy_to_molecular: 0.1,
                ..selected
            },
            &gate
        )
        .is_err()
    );
    let mut nonzero_trace = data;
    nonzero_trace.boundary_velocity_m_s[0] = Some([0.0, 1.0, 0.0]);
    assert!(
        run(
            nonzero_trace,
            reservoir_boundaries(),
            pressure(&x),
            selected,
            &gate
        )
        .is_err()
    );
    let mut conflicting = reservoir_boundaries();
    conflicting[0] = FvBoundary::Pressure { pressure_pa: 0.0 };
    assert!(run(data, conflicting, pressure(&x), selected, &gate).is_err());
    assert!(
        run(
            data,
            reservoir_boundaries(),
            pressure(&x),
            FvSgsConfig {
                usage: SgsUse::AcceptanceRender,
                ..selected
            },
            &gate
        )
        .is_err()
    );
    // These fail after preparing the viscous candidate: an invalid original
    // pressure profile and unrepresentable physical molecular heating.
    let short_x = [Some(0.0); 3];
    let bad_profile = StaggeredBoundaryPressure {
        x_pa: &short_x,
        ..pressure(&x)
    };
    assert!(run(data, reservoir_boundaries(), bad_profile, selected, &gate).is_err());
    let huge_h = [1e30; 2];
    let caloric = reservoirs(&mass, &huge_h, &p, 10.0, [0.0; 2]);
    assert!(
        run(
            caloric,
            reservoir_boundaries(),
            pressure(&x),
            selected,
            &gate
        )
        .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        run(
            data,
            reservoir_boundaries(),
            pressure(&x),
            selected,
            &cancelled
        ),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(
        run(data, reservoir_boundaries(), pressure(&x), selected, &gate).unwrap(),
        accepted
    );
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(mass, [2.0; 2]);
    assert_eq!(h, [10.0; 2]);
    assert_eq!(p, [[2.0, 0.0, 0.0], [-2.0, 0.0, 0.0]]);
}
