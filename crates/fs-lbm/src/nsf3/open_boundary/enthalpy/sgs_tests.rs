//! Coupled SGS heat contacts use one original face and conservative history.
use super::tests::{grid, thermal};
use super::*;
use crate::d3q27::D3q27;
use crate::nsf3::{NsfBoundaryDrive, NsfCharacteristicPatch, NsfSgsConfig};
use crate::sgs3::{SgsConfig, SgsError, SgsModel, SgsUse};

fn config(model: SgsModel) -> NsfSgsConfig {
    NsfSgsConfig {
        closure: SgsConfig {
            model,
            filter_width_lu: 1.0,
            maximum_eddy_to_molecular: 100.0,
            usage: SgsUse::Exploration,
        },
        turbulent_prandtl: 0.9,
    }
}

fn contact_grid(
    axis: usize,
    high: bool,
    workers: usize,
    model: SgsModel,
    sgs_first: bool,
    walls: bool,
) -> NsfGrid3<D3q27> {
    let template = grid(axis, high);
    let mut spec = template.spec.clone();
    spec.workers = workers;
    if walls {
        spec.periodic[(axis + 1) % 3] = false;
    }
    let dims = spec.dims;
    let unit = template.law.dx_m / template.law.dt_s;
    let mut result = NsfGrid3::new(
        spec,
        |_| false,
        |x| {
            let phase = std::array::from_fn::<_, 3, _>(|a| {
                core::f64::consts::TAU * (x[a] as f64 + 0.5) / dims[a] as f64
            });
            let a = axis;
            let b = (axis + 1) % 3;
            let temperature = 300.4 + 0.05 * phase[a].cos();
            let mut velocity = [0.0; 3];
            velocity[a] =
                (if high { -2.0 } else { 2.0 }) / unit + 0.0004 * phase[a].sin() * phase[b].cos();
            velocity[b] = 0.25 / unit - 0.0004 * phase[a].cos() * phase[b].sin();
            (300.0 / temperature, velocity, temperature)
        },
    )
    .unwrap()
    .with_characteristic_axis(*template.characteristic_axis().unwrap())
    .unwrap();
    if sgs_first {
        result = result.with_sgs(config(model)).unwrap();
    }
    result = result
        .with_characteristic_thermal_policy([NsfCharacteristicThermalPolicy::EnthalpyFlux; 2])
        .unwrap();
    let mut drives = [thermal(-3.0, 308.0); 2];
    drives[usize::from(high)] = thermal(8.0, 301.0);
    result.set_boundary_drives(drives).unwrap();
    if !sgs_first {
        result = result.with_sgs(config(model)).unwrap();
    }
    result
}

fn check_flux(grid: &NsfGrid3<D3q27>, axis: usize, high: bool, upstream: f64) {
    let (rho, _) = grid.macro_fields(&grid.f);
    // Independent stress values verify that mechanical traction work stays
    // separate from prescribed thermal enthalpy; it is not added twice.
    let sig = vec![[1e-5, 2e-5, -0.5e-5, 0.7e-5, -0.2e-5, 0.3e-5]; grid.n];
    let mut active = false;
    for side in [false, true] {
        let direction = if side { 1 } else { -1 };
        for trace in grid.open_boundary_traces(side).unwrap() {
            let i = grid.index(trace.cell);
            let face = grid
                .characteristic_face(grid.open_faces.as_ref(), i, axis, direction)
                .unwrap();
            let normal = direction as f64 * face.u[axis];
            let (total, _, _) = grid
                .open_energy_flux(
                    face,
                    i,
                    axis,
                    direction,
                    &rho,
                    &grid.theta,
                    &sig,
                    &grid.sgs_transport,
                )
                .unwrap();
            let upstream = grid.characteristic_patch[usize::from(side)]
                .filter(|patch| patch.shape.contains(trace.cell, axis))
                .and_then(|patch| match patch.drive {
                    NsfBoundaryDrive::PressureTemperature {
                        inflow_temperature_k,
                        ..
                    } => Some(inflow_temperature_k),
                    _ => None,
                })
                .unwrap_or(upstream);
            let enthalpy = if side == high {
                assert!(normal < 0.0);
                let pi = grid.law.pressure_pa(rho[i], grid.theta[i]);
                let molecular = grid
                    .enthalpy_flux_theta(
                        grid.theta[i],
                        pi,
                        trace.pressure_pa,
                        normal,
                        upstream,
                        0.0,
                    )
                    .unwrap();
                if grid.sgs_transport[i].turbulent_conductivity_lu > 1e-9 {
                    active = true;
                    assert!((face.theta - grid.theta[i]).abs() < (molecular - grid.theta[i]).abs());
                }
                grid.contact_enthalpy(grid.law.theta_of_kelvin(upstream), trace.pressure_pa)
                    .unwrap()
                    .1
            } else {
                assert!(normal > 0.0);
                assert_eq!(face.theta.to_bits(), grid.theta[i].to_bits());
                grid.contact_enthalpy(grid.theta[i], trace.pressure_pa)
                    .unwrap()
                    .1
            };
            let row = match axis {
                0 => [sig[i][0], sig[i][3], sig[i][4]],
                1 => [sig[i][3], sig[i][1], sig[i][5]],
                _ => [sig[i][4], sig[i][5], sig[i][2]],
            };
            let work: f64 = row.iter().zip(face.u).map(|(s, u)| s * u).sum();
            let kinetic = 0.5 * face.u.iter().map(|v| v * v).sum::<f64>();
            let expected = face.rho * normal * (enthalpy + kinetic);
            assert!(
                (direction as f64 * (total + work) - expected).abs()
                    < 3e-13 * expected.abs().max(1e-10)
            );
        }
    }
    assert!(
        active,
        "the incoming thermal contact must exercise modeled heat"
    );
}

#[test]
fn g1_coupled_enthalpy_uses_total_heat_on_all_axes_sides_models_and_attachment_orders() {
    for axis in 0..3 {
        for high in [false, true] {
            for model in [
                SgsModel::Smagorinsky { cs: 0.17 },
                SgsModel::Wale { cw: 0.5 },
                SgsModel::Vreman { cv: 0.07 },
            ] {
                let mut first = contact_grid(axis, high, 1, model, true, false);
                let second = contact_grid(axis, high, 1, model, false, false);
                assert_eq!(first.total_energy_j(), second.total_energy_j());
                for side in [false, true] {
                    for (a, b) in first
                        .open_boundary_traces(side)
                        .unwrap()
                        .iter()
                        .zip(second.open_boundary_traces(side).unwrap())
                    {
                        assert!((a.temperature_k - b.temperature_k).abs() < 2e-10);
                        assert!(
                            (a.outward_normal_velocity_m_s - b.outward_normal_velocity_m_s).abs()
                                < 2e-10
                        );
                    }
                }
                check_flux(&first, axis, high, 301.0);
                first = first
                    .with_characteristic_patch(
                        high,
                        NsfCharacteristicPatch {
                            start: [1, 1],
                            end: [2, 2],
                        },
                        thermal(8.0, 300.8),
                    )
                    .unwrap();
                check_flux(&first, axis, high, 301.0);
            }
        }
    }
}

#[test]
fn g1_force_predictor_closes_velocity_and_matches_the_face_sgs_coefficient() {
    let grid = contact_grid(0, false, 1, SgsModel::Wale { cw: 0.5 }, true, false);
    let (rho, momentum) = grid.macro_fields(&grid.f);
    let physical = grid.velocity_field(&rho, &momentum, &grid.force);
    let mut guess = physical.clone();
    for velocity in &mut guess {
        velocity[0] += 0.0001;
    }
    let result = grid
        .sgs_contact_predictor(&rho, &momentum, &grid.theta, &grid.material, &guess)
        .unwrap();
    assert_eq!(
        result.velocity,
        grid.velocity_field(&rho, &momentum, &result.force)
    );
    let from_physical = grid
        .sgs_field(
            &rho,
            &result.velocity,
            &grid.theta,
            &grid.material,
            result.faces.as_ref(),
        )
        .unwrap();
    for (a, b) in result.transport.iter().zip(from_physical) {
        assert!((a.turbulent_conductivity_lu - b.turbulent_conductivity_lu).abs() < 2e-14);
    }
    for (a, b) in physical.iter().zip(result.velocity) {
        for (a, b) in a.iter().zip(b) {
            assert!((a - b).abs() < 3e-14);
        }
    }
    for trace in grid.open_boundary_traces(false).unwrap() {
        let i = grid.index(trace.cell);
        let face = grid
            .characteristic_face(result.faces.as_ref(), i, 0, -1)
            .unwrap();
        let (total, _, _) = grid
            .open_energy_flux(
                face,
                i,
                0,
                -1,
                &rho,
                &grid.theta,
                &vec![[0.0; 6]; grid.n],
                &result.transport,
            )
            .unwrap();
        let enthalpy = grid
            .contact_enthalpy(grid.law.theta_of_kelvin(301.0), trace.pressure_pa)
            .unwrap()
            .1;
        let expected =
            face.rho * face.u[0] * (enthalpy + 0.5 * face.u.iter().map(|v| v * v).sum::<f64>());
        assert!((total - expected).abs() < 3e-15);
    }
}

#[test]
fn g0_g4_g5_coupled_enthalpy_trials_conserve_replay_and_drop_every_modeled_history() {
    for model in [
        SgsModel::Smagorinsky { cs: 0.17 },
        SgsModel::Wale { cw: 0.5 },
        SgsModel::Vreman { cv: 0.07 },
    ] {
        let mut serial = contact_grid(0, false, 1, model, true, false);
        let mut parallel = contact_grid(0, false, 3, model, true, false);
        let before = (
            serial.state_hash(),
            serial.force.clone(),
            serial.div.clone(),
            serial.sgs_transport.clone(),
            serial.open_boundary_traces(false),
            serial.cumulative_ledger(),
            serial.boundary_drives(),
        );
        let drives = [thermal(8.0, 300.9), thermal(-3.0, 308.0)];
        let work = 2 * serial.cell_count();
        let trial = serial
            .preview_boundary_interval(drives, 2, work, || false)
            .unwrap();
        let mass: f64 = trial
            .boundary_exchange()
            .iter()
            .map(|x| x.mass_out_kg)
            .sum();
        let energy: f64 = trial
            .boundary_exchange()
            .iter()
            .map(|x| x.total_energy_out_j)
            .sum();
        assert!(
            (trial.grid().total_mass_kg() - trial.previous_mass_kg() + mass).abs()
                < 5e-12 * trial.previous_mass_kg()
        );
        assert!(
            (trial.grid().total_energy_j() - trial.previous_energy_j() + energy).abs()
                < 5e-12 * trial.previous_energy_j()
        );
        let expected = (
            trial.grid().state_hash(),
            *trial.energy_ledger(),
            *trial.boundary_exchange(),
        );
        drop(trial);
        assert_eq!(
            (
                serial.state_hash(),
                serial.force.clone(),
                serial.div.clone(),
                serial.sgs_transport.clone(),
                serial.open_boundary_traces(false),
                serial.cumulative_ledger(),
                serial.boundary_drives()
            ),
            before
        );
        let calls = std::cell::Cell::new(0);
        assert!(
            serial
                .preview_boundary_interval(drives, 2, work, || {
                    calls.set(calls.get() + 1);
                    calls.get() == 3
                })
                .is_err()
        );
        assert_eq!(serial.state_hash(), before.0);
        let retry = serial
            .preview_boundary_interval(drives, 2, work, || false)
            .unwrap();
        assert_eq!(
            (
                retry.grid().state_hash(),
                *retry.energy_ledger(),
                *retry.boundary_exchange()
            ),
            expected
        );
        retry.accept();
        parallel
            .preview_boundary_interval(drives, 2, work, || false)
            .unwrap()
            .accept();
        assert_eq!(serial.state_hash(), parallel.state_hash());
        assert_eq!(serial.sgs_transport, parallel.sgs_transport);
        assert_eq!(serial.cumulative_ledger(), parallel.cumulative_ledger());
    }
}

#[test]
fn g0_rejected_enthalpy_sgs_rebind_retains_original_force_flux_and_controls() {
    let template = grid(0, false);
    let mut selected =
        NsfGrid3::<D3q27>::new(template.spec.clone(), |_| false, |_| (1.0, [0.0; 3], 300.0))
            .unwrap()
            .with_characteristic_axis(*template.characteristic_axis().unwrap())
            .unwrap()
            .with_characteristic_thermal_policy([NsfCharacteristicThermalPolicy::EnthalpyFlux; 2])
            .unwrap();
    let mut settings = config(SgsModel::Smagorinsky { cs: 0.17 });
    settings.closure.maximum_eddy_to_molecular = 1e-5;
    selected = selected.with_sgs(settings).unwrap();
    let before = (
        selected.state_hash(),
        selected.force.clone(),
        selected.div.clone(),
        selected.sgs_transport.clone(),
        selected.open_boundary_traces(false),
        selected.boundary_drives(),
    );
    let result =
        selected.set_boundary_drives([thermal(10.0, 301.0), NsfBoundaryDrive::PressureOffset(0.0)]);
    assert!(matches!(
        result,
        Err(NsfError::Sgs(SgsError::ViscosityRatioExceeded { .. }))
    ));
    assert_eq!(
        (
            selected.state_hash(),
            selected.force.clone(),
            selected.div.clone(),
            selected.sgs_transport.clone(),
            selected.open_boundary_traces(false),
            selected.boundary_drives()
        ),
        before
    );
}

#[test]
fn g1_transverse_wall_rebind_updates_the_coupled_contact_without_changing_physical_energy() {
    let mut selected = contact_grid(0, false, 1, SgsModel::Smagorinsky { cs: 0.17 }, true, true);
    let physical: Vec<_> = (0..selected.n)
        .map(|i| selected.cell(selected.at(i)))
        .collect();
    let energy = selected.total_energy_j();
    let transport = selected.sgs_transport.clone();
    selected
        .set_edge_velocity(1, false, [0.001, 0.0, 0.0])
        .unwrap();
    assert_eq!(selected.total_energy_j(), energy);
    assert_ne!(selected.sgs_transport, transport);
    for (i, old) in physical.iter().enumerate() {
        let new = selected.cell(selected.at(i));
        assert!((new.0 - old.0).abs() < 2e-14);
        assert_eq!(new.2, old.2);
        for (a, b) in new.1.into_iter().zip(old.1) {
            assert!((a - b).abs() < 2e-14);
        }
    }
    check_flux(&selected, 0, false, 301.0);
}

#[test]
fn g0_zero_sgs_enthalpy_contacts_preserve_molecular_trajectories() {
    for model in [
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
    ] {
        let make = || {
            grid(0, false)
                .with_characteristic_thermal_policy(
                    [NsfCharacteristicThermalPolicy::EnthalpyFlux; 2],
                )
                .unwrap()
        };
        let mut molecular = make();
        let mut zero = make().with_sgs(config(model)).unwrap();
        for solver in [&mut molecular, &mut zero] {
            solver
                .set_boundary_drives([thermal(8.0, 301.0), thermal(-3.0, 308.0)])
                .unwrap();
        }
        for _ in 0..3 {
            assert_eq!(molecular.state_hash(), zero.state_hash());
            assert_eq!(molecular.step().unwrap(), zero.step().unwrap());
        }
    }
}
