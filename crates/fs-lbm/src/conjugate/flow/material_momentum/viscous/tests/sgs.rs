//! Actual-density SGS composition on the original molecular/pressure/transport owner.

use super::*;
use crate::sgs3::{SgsModel, SgsUse};
use fs_material::MaterialError;
use std::cell::Cell;

fn model(model: SgsModel) -> FvSgsConfig {
    FvSgsConfig {
        model,
        filter_width_factor: 1.0,
        maximum_eddy_to_molecular: 10.0,
        usage: SgsUse::Exploration,
    }
}

/// Independent positive solution of M(u-u0)/dt + K_mu*u + C_sgs*u^2=0.
/// Only the interior x face moves. The SGS normal stress is 4/3 mu_t du/dx,
/// while the original molecular normal stress remains 2 mu du/dx.
fn graded_speed(mu: [f64; 2], cs: f64) -> f64 {
    let volume: [f64; 2] = [2.0, 6.0];
    let width: [f64; 2] = [1.0, 3.0];
    let rho: [f64; 2] = [2.0, 6.0];
    let (face_mass, area, dt) = (20.0, 2.0, 0.5);
    let molecular = (0..2).map(|c| 2.0 * area * mu[c] / width[c]).sum::<f64>();
    let modeled = (0..2)
        .map(|c| {
            (4.0 / 3.0) * area * rho[c] * cs.powi(2) * volume[c].cbrt().powi(2) * 2.0_f64.sqrt()
                / width[c].powi(2)
        })
        .sum::<f64>();
    let linear = face_mass + dt * molecular;
    2.0 * face_mass / (linear + (linear.powi(2) + 4.0 * dt * modeled * face_mass).sqrt())
}

fn assert_energy_balance(report: &MaterialViscosityReport) {
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
fn graded_actual_density_uses_deviatoric_sgs_and_heats_only_molecular_loss() {
    let domain = VoxelDomain::graded([vec![1.0, 3.0], vec![2.0], vec![1.0]]).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [4.0, 36.0];
    let h = [10.0, 30.0];
    let p = [[4.0, 0.0, 0.0], [36.0, 0.0, 0.0]];
    let slow = liquid(1.0);
    let fast = liquid(3.0);
    let sources: [Option<&dyn FluidTransport>; 2] = [Some(&slow), Some(&fast)];
    let cs = 0.2;
    let run = |sources| {
        field
            .advance_material_momentum_with_viscosity_and_sgs(
                &domain,
                input(&mass, &h, &p),
                MaterialMomentumViscosityInput {
                    sources,
                    boundaries: [FvBoundary::Symmetry; 6],
                },
                model(SgsModel::Smagorinsky { cs }),
                controls(),
                viscous_controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap()
    };
    let next = run(&sources);
    let speed = graded_speed([1.0, 3.0], cs);
    close(next.viscosity.face_impulse_kg_m_s[0], 20.0 * (speed - 1.0));
    close(
        next.viscosity.native_kinetic_energy_change_j,
        10.0 * (speed.powi(2) - 1.0),
    );
    close(
        next.viscosity.backward_euler_loss_j,
        10.0 * (speed - 1.0).powi(2),
    );
    assert_eq!(next.viscosity.cell_sgs_transfer_j.len(), 2);
    let mut molecular_heat = 0.0;
    let mut modeled_transfer = 0.0;
    for c in 0..2 {
        let width: f64 = if c == 0 { 1.0 } else { 3.0 };
        let volume = domain.volume(c);
        let rho = mass[c] / volume;
        let mu = if c == 0 { 1.0 } else { 3.0 };
        let rate = speed / width;
        let nu_t = cs.powi(2) * volume.cbrt().powi(2) * 2.0_f64.sqrt() * rate;
        let heat = 0.5 * 2.0 * mu * volume * rate.powi(2);
        let transfer = 0.5 * (4.0 / 3.0) * rho * nu_t * volume * rate.powi(2);
        molecular_heat += heat;
        modeled_transfer += transfer;
        close(next.viscosity.cell_molecular_heat_j[c], heat);
        close(next.viscosity.cell_sgs_transfer_j[c], transfer);
        close(next.advance.transport.mass_kg[c], mass[c]);
        close(
            next.advance.transport.specific_enthalpy_j_kg[c],
            h[c] + heat / mass[c],
        );
        close(next.advance.momentum.momentum_kg_m_s[c][0], 0.0);
    }
    close(next.viscosity.molecular_heat_j, molecular_heat);
    close(next.viscosity.sgs_transfer_j, modeled_transfer);
    close(next.represented_molecular_heat_j, molecular_heat);
    assert!(modeled_transfer > 0.1);
    assert_energy_balance(&next.viscosity);
    // The original closed projection removes the remaining throughflow once,
    // and its zero accepted mass routes leave both inventories in place.
    close(next.advance.transport.projection.field.flux_x(1, 0, 0), 0.0);
    assert_eq!(
        next.advance.momentum.report.represented_impulse_kg_m_s,
        [0.0; 3]
    );

    // Original source rho may differ arbitrarily from actual M/V. It supplies
    // neither face inertia, dynamic SGS conversion nor viscosity-cap density.
    let low_reference = FluidTransportState::try_new(0.001, 2.0, 1.0, 1.0).unwrap();
    let high_reference = FluidTransportState::try_new(20_000.0, 2.0, 1.0, 3.0).unwrap();
    let other_sources: [Option<&dyn FluidTransport>; 2] =
        [Some(&low_reference), Some(&high_reference)];
    assert_eq!(run(&other_sources), next);
}

#[test]
fn every_zero_model_is_the_exact_original_material_advance() {
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
    let gate = CancelGate::new_clock_free();
    let original = field
        .advance_material_momentum_with_viscosity(
            &domain,
            input(&mass, &h, &p),
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
            .advance_material_momentum_with_viscosity_and_sgs(
                &domain,
                input(&mass, &h, &p),
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

struct Source {
    state: Cell<FluidTransportState>,
    calls: Cell<usize>,
    refuse: Cell<bool>,
}

impl FluidTransport for Source {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.calls.set(self.calls.get() + 1);
        if self.refuse.get() {
            Err(MaterialError::State {
                what: "source outside its supported state".into(),
            })
        } else {
            Ok(self.state.get())
        }
    }
}

#[test]
fn original_sources_refresh_once_per_step_without_reusing_sgs_or_reference_density() {
    let domain = VoxelDomain::graded([vec![1.0, 3.0], vec![2.0], vec![1.0]]).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [4.0, 36.0];
    let h = [10.0, 30.0];
    let p = [[4.0, 0.0, 0.0], [36.0, 0.0, 0.0]];
    let source = |mu| Source {
        state: Cell::new(FluidTransportState::try_new(1000.0, 2.0, 1.0, mu).unwrap()),
        calls: Cell::new(0),
        refuse: Cell::new(false),
    };
    let low = source(1.0);
    let high = source(3.0);
    let sources: [Option<&dyn FluidTransport>; 2] = [Some(&low), Some(&high)];
    let run = || {
        field.advance_material_momentum_with_viscosity_and_sgs(
            &domain,
            input(&mass, &h, &p),
            MaterialMomentumViscosityInput {
                sources: &sources,
                boundaries: [FvBoundary::Symmetry; 6],
            },
            model(SgsModel::Smagorinsky { cs: 0.2 }),
            controls(),
            viscous_controls(),
            &CancelGate::new_clock_free(),
        )
    };
    let first = run().unwrap();
    assert_eq!((low.calls.get(), high.calls.get()), (1, 1));
    low.state
        .set(FluidTransportState::try_new(0.001, 2.0, 1.0, 4.0).unwrap());
    let second = run().unwrap();
    assert_eq!((low.calls.get(), high.calls.get()), (2, 2));
    let speed = graded_speed([4.0, 3.0], 0.2);
    close(
        second.viscosity.face_impulse_kg_m_s[0],
        20.0 * (speed - 1.0),
    );
    assert!(second.viscosity.face_impulse_kg_m_s[0] < first.viscosity.face_impulse_kg_m_s[0]);
    assert_energy_balance(&second.viscosity);
    high.refuse.set(true);
    assert!(matches!(
        run(),
        Err(ChtError::InvalidInput {
            field: "flow.viscosity_source",
            ..
        })
    ));
    assert_eq!((low.calls.get(), high.calls.get()), (3, 3));
    high.refuse.set(false);
    assert_eq!(run().unwrap(), second);
    assert_eq!((low.calls.get(), high.calls.get()), (4, 4));
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(mass, [4.0, 36.0]);
    assert_eq!(h, [10.0, 30.0]);
    assert_eq!(p, [[4.0, 0.0, 0.0], [36.0, 0.0, 0.0]]);
}

#[test]
fn sgs_cap_source_boundary_work_and_cancellation_refuse_without_publication() {
    let domain = VoxelDomain::graded([vec![1.0, 3.0], vec![2.0], vec![1.0]]).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [4.0, 36.0];
    let h = [10.0, 30.0];
    let p = [[4.0, 0.0, 0.0], [36.0, 0.0, 0.0]];
    let source = FluidTransportState::try_new(0.001, 2.0, 1.0, 1.0).unwrap();
    let sources: [Option<&dyn FluidTransport>; 2] = [Some(&source); 2];
    let viscosity = MaterialMomentumViscosityInput {
        sources: &sources,
        boundaries: [FvBoundary::Symmetry; 6],
    };
    let selected = model(SgsModel::Smagorinsky { cs: 0.2 });
    let run = |data, viscosity, selected, work, gate: &CancelGate| {
        field.advance_material_momentum_with_viscosity_and_sgs(
            &domain,
            data,
            viscosity,
            selected,
            controls(),
            work,
            gate,
        )
    };
    let gate = CancelGate::new_clock_free();
    let accepted = run(
        input(&mass, &h, &p),
        viscosity,
        selected,
        viscous_controls(),
        &gate,
    )
    .unwrap();
    for invalid in [
        FvSgsConfig {
            maximum_eddy_to_molecular: 0.01,
            ..selected
        },
        FvSgsConfig {
            filter_width_factor: 0.0,
            ..selected
        },
        FvSgsConfig {
            model: SgsModel::Smagorinsky { cs: f64::NAN },
            ..selected
        },
        FvSgsConfig {
            usage: SgsUse::AcceptanceRender,
            ..selected
        },
    ] {
        assert!(
            run(
                input(&mass, &h, &p),
                viscosity,
                invalid,
                viscous_controls(),
                &gate
            )
            .is_err()
        );
    }
    for work in [
        MaterialViscosityControls {
            max_cells: 1,
            ..viscous_controls()
        },
        MaterialViscosityControls {
            max_sweeps: 0,
            ..viscous_controls()
        },
        MaterialViscosityControls {
            max_sweeps: 1,
            ..viscous_controls()
        },
    ] {
        assert!(run(input(&mass, &h, &p), viscosity, selected, work, &gate).is_err());
    }
    let missing: [Option<&dyn FluidTransport>; 2] = [Some(&source), None];
    assert!(
        run(
            input(&mass, &h, &p),
            MaterialMomentumViscosityInput {
                sources: &missing,
                ..viscosity
            },
            selected,
            viscous_controls(),
            &gate,
        )
        .is_err()
    );
    let mut open = viscosity;
    open.boundaries[0] = FvBoundary::Outlet;
    assert!(
        run(
            input(&mass, &h, &p),
            open,
            selected,
            viscous_controls(),
            &gate
        )
        .is_err()
    );
    // This refusal occurs after the accepted native molecular/SGS endpoint:
    // unrepresentable physical heat cannot escape as momentum-only success.
    assert!(
        run(
            input(&mass, &[1e30; 2], &p),
            viscosity,
            selected,
            viscous_controls(),
            &gate
        )
        .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        run(
            input(&mass, &h, &p),
            viscosity,
            selected,
            viscous_controls(),
            &cancelled
        ),
        Err(ChtError::Cancelled),
    ));
    struct Cancelling<'a>(&'a CancelGate, FluidTransportState);
    impl FluidTransport for Cancelling<'_> {
        fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
            self.0.request();
            Ok(self.1)
        }
    }
    let during_source = CancelGate::new_clock_free();
    let cancels = Cancelling(&during_source, source);
    let cancel_sources: [Option<&dyn FluidTransport>; 2] = [Some(&cancels), Some(&source)];
    assert!(matches!(
        run(
            input(&mass, &h, &p),
            MaterialMomentumViscosityInput {
                sources: &cancel_sources,
                ..viscosity
            },
            selected,
            viscous_controls(),
            &during_source,
        ),
        Err(ChtError::Cancelled),
    ));
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(mass, [4.0, 36.0]);
    assert_eq!(h, [10.0, 30.0]);
    assert_eq!(p, [[4.0, 0.0, 0.0], [36.0, 0.0, 0.0]]);
    assert_eq!(
        run(
            input(&mass, &h, &p),
            viscosity,
            selected,
            viscous_controls(),
            &gate
        )
        .unwrap(),
        accepted
    );
}
