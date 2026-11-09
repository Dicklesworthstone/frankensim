//! Independent molecular bulk response composed with deviatoric material SGS.

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
fn graded_longitudinal_root_retains_authored_bulk_heat_and_separate_sgs_transfer() {
    let domain = VoxelDomain::graded([vec![1.0, 3.0], vec![2.0], vec![1.0]]).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [4.0, 36.0];
    let h = [10.0, 30.0];
    let p = [[4.0, 0.0, 0.0], [36.0, 0.0, 0.0]];
    let slow = liquid(1.0);
    let fast = liquid(3.0);
    let sources: [Option<&dyn FluidTransport>; 2] = [Some(&slow), Some(&fast)];
    let volume: [f64; 2] = [2.0, 6.0];
    let width: [f64; 2] = [1.0, 3.0];
    let rho: [f64; 2] = [2.0, 6.0];
    let mu: [f64; 2] = [1.0, 3.0];
    let cs: f64 = 0.2;
    let (dt, area, face_mass) = (0.5, 2.0, 20.0);
    for bulk in [[0.0, 0.0], [0.2, 4.0], [4.0, 20.0]] {
        let next = field
            .advance_material_momentum_with_bulk_viscosity_and_sgs(
                &domain,
                input(&mass, &h, &p),
                MaterialMomentumViscosityInput {
                    sources: &sources,
                    boundaries: [FvBoundary::Symmetry; 6],
                },
                &bulk,
                model(SgsModel::Smagorinsky { cs }),
                controls(),
                viscous_controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        // Independent one-DOF equation:
        // M(u-1)/dt + K_mu,zeta*u + C_sgs*u^2 = 0.
        // Authored zeta changes only K; there is no modeled bulk stress.
        let longitudinal: [f64; 2] = core::array::from_fn(|c| (4.0 / 3.0) * mu[c] + bulk[c]);
        let molecular = (0..2)
            .map(|c| area * longitudinal[c] / width[c])
            .sum::<f64>();
        let modeled = (0..2)
            .map(|c| {
                (4.0 / 3.0) * area * rho[c] * cs.powi(2) * volume[c].cbrt().powi(2) * 2.0_f64.sqrt()
                    / width[c].powi(2)
            })
            .sum::<f64>();
        let linear = face_mass + dt * molecular;
        let speed =
            2.0 * face_mass / (linear + (linear.powi(2) + 4.0 * dt * modeled * face_mass).sqrt());
        close(
            next.viscosity.face_impulse_kg_m_s[0],
            face_mass * (speed - 1.0),
        );
        close(
            next.viscosity.native_kinetic_energy_change_j,
            0.5 * face_mass * (speed.powi(2) - 1.0),
        );
        let mut heat_sum = 0.0;
        let mut transfer_sum = 0.0;
        for c in 0..2 {
            let rate = speed / width[c];
            let nu_t = cs.powi(2) * volume[c].cbrt().powi(2) * 2.0_f64.sqrt() * rate;
            let heat = dt * longitudinal[c] * volume[c] * rate.powi(2);
            let transfer = dt * (4.0 / 3.0) * rho[c] * nu_t * volume[c] * rate.powi(2);
            heat_sum += heat;
            transfer_sum += transfer;
            close(next.viscosity.cell_molecular_heat_j[c], heat);
            close(next.viscosity.cell_sgs_transfer_j[c], transfer);
            close(next.advance.transport.mass_kg[c], mass[c]);
            close(
                next.advance.transport.specific_enthalpy_j_kg[c],
                h[c] + heat / mass[c],
            );
        }
        close(next.viscosity.molecular_heat_j, heat_sum);
        close(next.viscosity.sgs_transfer_j, transfer_sum);
        close(next.represented_molecular_heat_j, heat_sum);
        assert!(transfer_sum > 0.0);
        balance(&next.viscosity);
    }
}

#[test]
fn solenoidal_native_mode_and_material_transport_are_independent_of_bulk_viscosity() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let (rho, mu, amplitude, dt, cs): (f64, f64, f64, f64, f64) = (1.7, 0.8, 0.2, 0.2, 0.3);
    let mass = [rho; 4];
    let h = [10.0; 4];
    // Each cell has du_x/dx = -du_y/dy. The paired tangential
    // derivatives also cancel, so div(u)=0 and SGS sees |S|=2a.
    let p: [[f64; 3]; 4] = core::array::from_fn(|cell| {
        let at = domain.coords(cell);
        [
            rho * if at[1] == 0 { amplitude } else { -amplitude },
            rho * if at[0] == 0 { -amplitude } else { amplitude },
            0.0,
        ]
    });
    let source = FluidTransportState::try_new(900.0, 2.0, 1.0, mu).unwrap();
    let sources: [Option<&dyn FluidTransport>; 4] = [Some(&source); 4];
    let data = MaterialMomentumAdvanceInput {
        material: MaterialAdvectionInput {
            mass_kg: &mass,
            specific_enthalpy_j_kg: &h,
            volume_sources: &[MaterialVolumeSource::None; 4],
            boundary_inflow: [None; 6],
        },
        momentum_kg_m_s: &p,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: &[[0.0; 3]; 4],
        source_velocity_m_s: &[None; 4],
        boundary_velocity_m_s: [None; 6],
    };
    let mut advance_controls = controls();
    advance_controls.mechanics.projection.duration_s = dt;
    advance_controls.material.projection = advance_controls.mechanics.projection;
    // Four equivalent moving faces: rho(a-a0)/dt + 4mu*a + 8rho*Cs^2*a^2=0.
    let linear = rho + 4.0 * dt * mu;
    let quadratic = 8.0 * dt * rho * cs.powi(2);
    let accepted = 2.0 * rho * amplitude
        / (linear + (linear.powi(2) + 4.0 * quadratic * rho * amplitude).sqrt());
    let heat = 16.0 * dt * mu * accepted.powi(2);
    let transfer = 32.0 * dt * rho * cs.powi(2) * accepted.powi(3);
    let mut reference: Option<ViscousMaterialMomentumAdvance> = None;
    for bulk in [0.0, 0.7, 3.0] {
        let next = field
            .advance_material_momentum_with_bulk_viscosity_and_sgs(
                &domain,
                data,
                MaterialMomentumViscosityInput {
                    sources: &sources,
                    boundaries: [FvBoundary::Symmetry; 6],
                },
                &[bulk; 4],
                model(SgsModel::Smagorinsky { cs }),
                advance_controls,
                viscous_controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        close(next.viscosity.molecular_heat_j, heat);
        close(next.viscosity.sgs_transfer_j, transfer);
        close(next.represented_molecular_heat_j, heat);
        close(
            next.viscosity.native_kinetic_energy_change_j,
            2.0 * rho * (accepted.powi(2) - amplitude.powi(2)),
        );
        close(
            next.viscosity.backward_euler_loss_j,
            2.0 * rho * (accepted - amplitude).powi(2),
        );
        for c in 0..4 {
            close(next.viscosity.cell_molecular_heat_j[c], heat / 4.0);
            close(next.viscosity.cell_sgs_transfer_j[c], transfer / 4.0);
            close(next.advance.transport.mass_kg[c], rho);
            close(
                next.advance.transport.specific_enthalpy_j_kg[c],
                h[c] + heat / (4.0 * rho),
            );
        }
        balance(&next.viscosity);
        if let Some(first) = &reference {
            for c in 0..4 {
                for a in 0..3 {
                    close(
                        next.advance.momentum.momentum_kg_m_s[c][a],
                        first.advance.momentum.momentum_kg_m_s[c][a],
                    );
                }
            }
            for (actual, expected) in next
                .advance
                .transport
                .projection
                .field
                .face_arrays()
                .into_iter()
                .zip(first.advance.transport.projection.field.face_arrays())
            {
                for (&a, &e) in actual.iter().zip(expected) {
                    close(a, e);
                }
            }
        } else {
            reference = Some(next);
        }
    }
}

#[test]
fn zero_sgs_preserves_original_bulk_and_invalid_bulk_cannot_relax_the_shear_cap() {
    let domain = VoxelDomain::graded([vec![1.0, 3.0], vec![2.0], vec![1.0]]).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [4.0, 36.0];
    let h = [10.0, 30.0];
    let p = [[4.0, 0.0, 0.0], [36.0, 0.0, 0.0]];
    let source = liquid(1.0);
    let sources: [Option<&dyn FluidTransport>; 2] = [Some(&source); 2];
    let viscosity = MaterialMomentumViscosityInput {
        sources: &sources,
        boundaries: [FvBoundary::Symmetry; 6],
    };
    let bulk = [0.2, 4.0];
    let gate = CancelGate::new_clock_free();
    let run = |bulk: &[f64], sgs| {
        field.advance_material_momentum_with_bulk_viscosity_and_sgs(
            &domain,
            input(&mass, &h, &p),
            viscosity,
            bulk,
            sgs,
            controls(),
            viscous_controls(),
            &gate,
        )
    };
    let original = field
        .advance_material_momentum_with_bulk_viscosity(
            &domain,
            input(&mass, &h, &p),
            viscosity,
            &bulk,
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
        let next = run(&bulk, model(selected)).unwrap();
        assert_eq!(next, original);
        assert_eq!(next.viscosity.sgs_transfer_j, 0.0);
        assert!(next.viscosity.cell_sgs_transfer_j.is_empty());
    }
    let selected = model(SgsModel::Smagorinsky { cs: 0.2 });
    let accepted = run(&bulk, selected).unwrap();
    for invalid in [
        &[0.0][..],
        &[-1.0, 0.0],
        &[0.0, f64::NAN],
        &[f64::INFINITY, 0.0],
    ] {
        assert!(run(invalid, selected).is_err());
        assert!(run(invalid, model(SgsModel::Smagorinsky { cs: 0.0 })).is_err());
    }
    // The cap compares rho*nu_t to the original SHEAR mu, regardless of how
    // large an independently authored bulk damping coefficient becomes.
    assert!(matches!(
        run(
            &[1e6; 2],
            FvSgsConfig {
                maximum_eddy_to_molecular: 0.01,
                ..selected
            }
        ),
        Err(ChtError::InvalidInput {
            field: "simple.sgs",
            ..
        }),
    ));
    assert_eq!(run(&bulk, selected).unwrap(), accepted);
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(mass, [4.0, 36.0]);
    assert_eq!(h, [10.0, 30.0]);
    assert_eq!(p, [[4.0, 0.0, 0.0], [36.0, 0.0, 0.0]]);
    assert_eq!(bulk, [0.2, 4.0]);
}
