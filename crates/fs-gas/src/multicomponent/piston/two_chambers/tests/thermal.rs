//! Original three-body calorics, moving pressure work and atomic publication.

use super::*;
use crate::multicomponent::wall_exchange::{
    MaterialWallCaloric, SharedWallControl, SharedWallStorage, WallContact, WallControl, WallDrag,
    WallThermal,
};
use fs_thermal::enthalpy::{
    Boundary as ThermalBoundary, CellGeometry, CellId, DepositedCell, DepositionError, MaterialId,
    ThermalLaw, ThermalStorageBasis,
};

mod suspension;

fn heat_controls() -> WallControl {
    WallControl {
        absolute_energy_tolerance_j: 1e-13,
        relative_energy_tolerance: 1e-13,
        ..WallControl::default()
    }
}

fn fields(count: [usize; 2], datum: f64) -> (MixtureEuler1d, MixtureEuler1d) {
    (
        chamber(left_gas(datum), 0.0, 1.0, count[0], AREA, 1.0),
        chamber(right_gas(datum), 1.0, 1.0, count[1], AREA, 1.0),
    )
}

fn finite(
    left: &MixtureEuler1d,
    right: &MixtureEuler1d,
    temperature: f64,
    conductance: [f64; 2],
) -> SharedPistonThermal<'static> {
    SharedPistonThermal::new(
        left,
        right,
        SharedWallStorage {
            thermal: WallThermal::Lumped {
                temperature_k: temperature,
                heat_capacity_j_per_k: 0.6,
                conductance_w_per_k: 0.0,
            },
            heat_storage: None,
        },
        conductance,
        SharedWallControl::default(),
        || false,
    )
    .unwrap()
}

fn wall_temperature(wall: &SharedPistonThermal<'_>) -> f64 {
    match wall.storage().thermal {
        WallThermal::Lumped { temperature_k, .. } => temperature_k,
        thermal => wall
            .storage()
            .heat_storage
            .as_ref()
            .unwrap()
            .temperature_k(thermal)
            .unwrap(),
    }
}

fn advance(
    left: &mut MixtureEuler1d,
    right: &mut MixtureEuler1d,
    thermal: &mut SharedPistonThermal<'_>,
    piston: &mut Piston,
    duration: f64,
    cancelled: impl FnMut() -> bool,
) -> Result<PistonStep<SharedPistonThermalStepReport>, SharedPistonThermalError> {
    left.step_with_shared_piston_and_heat(
        right,
        thermal,
        piston,
        Time::new(duration),
        StepControl::default(),
        heat_controls(),
        0.0,
        controls(),
        cancelled,
    )
}

// Independent constant-capacity three-body backward Euler solution. Each
// b_i accounts for its ORIGINAL gas amount; there is only one Cw and Tw.
fn linear_heat(t: [f64; 2], capacity: [f64; 2], tw: f64, g: [f64; 2], dt: f64) -> (f64, [f64; 2]) {
    let b = [
        capacity[0] * (g[0] * dt) / (capacity[0] + g[0] * dt),
        capacity[1] * (g[1] * dt) / (capacity[1] + g[1] * dt),
    ];
    let endpoint = (0.6 * tw + b[0] * t[0] + b[1] * t[1]) / (0.6 + b[0] + b[1]);
    (
        endpoint,
        [b[0] * (endpoint - t[0]), b[1] * (endpoint - t[1])],
    )
}

fn end_temperature(flow: &MixtureEuler1d, index: usize, remove_heat: f64) -> f64 {
    let cell = flow.cell(index).unwrap();
    flow.gas()
        .primitive(
            crate::multicomponent::Cell {
                energy: cell.energy - remove_heat / flow.mesh().volumes()[index],
                ..cell
            },
            || false,
        )
        .unwrap()
        .temperature
}

#[test]
fn g1_one_shared_heat_inventory_drives_both_pressures_with_old_and_moved_volume_balances() {
    let g = [1.7, 0.4];
    for initial_wall in [0.4, 3.0] {
        let (mut left, mut right) = fields([2, 2], 1e18);
        let mut wall = finite(&left, &right, initial_wall, g);
        let mut piston = Piston::new(0.4, 0.0).unwrap();
        let original = [left.totals(), right.totals()];
        let mut old_t = [1.0, 2.0];
        let mut old_wall = initial_wall;
        for _ in 0..6 {
            let capacity = [
                2.5 * left.cell(1).unwrap().density() * left.mesh().volumes()[1],
                1.75 * right.cell(0).unwrap().density() * right.mesh().volumes()[0],
            ];
            let report = advance(&mut left, &mut right, &mut wall, &mut piston, 0.01, || {
                false
            })
            .unwrap();
            assert!(report.admitted(controls()));
            let gas = report.gas;
            let hdt = 0.5 * gas.transport.left.elapsed.value();
            let (middle_wall, q1) = linear_heat(old_t, capacity, old_wall, g, hdt);
            let actual_q2 = [
                gas.after.left.energy - gas.transport.left.after.energy,
                gas.after.right.energy - gas.transport.right.after.energy,
            ];
            let middle_t = [
                end_temperature(&left, 1, actual_q2[0]),
                end_temperature(&right, 0, actual_q2[1]),
            ];
            let moved_capacity = [
                2.5 * left.cell(1).unwrap().density() * left.mesh().volumes()[1],
                1.75 * right.cell(0).unwrap().density() * right.mesh().volumes()[0],
            ];
            let (end_wall, q2) = linear_heat(middle_t, moved_capacity, middle_wall, g, hdt);
            for (side, (transport, before, after)) in [
                (gas.transport.left, gas.before.left, gas.after.left),
                (gas.transport.right, gas.before.right, gas.after.right),
            ]
            .into_iter()
            .enumerate()
            {
                close(transport.before.energy - before.energy, q1[side], 3e-12);
                close(after.energy - transport.after.energy, q2[side], 3e-12);
                close(gas.heat.heat_to_gas_j[side], q1[side] + q2[side], 5e-12);
                for k in 0..before.species_count {
                    close(before.species_masses[k], after.species_masses[k], 3e-14);
                }
            }
            close(wall_temperature(&wall), end_wall, 1e-11);
            close(
                gas.heat.wall_energy_change_j,
                0.6 * (end_wall - old_wall),
                1e-11,
            );
            close(gas.heat.energy_residual_j, 0.0, 1e-12);
            assert_eq!(wall.last_report(), Some(gas.heat));
            for half in gas.heat.stages {
                assert_eq!(half.elapsed.value(), hdt);
                assert_eq!(half.wall_impulse_ns, 0.0);
                assert_eq!(half.mechanical_work_to_gas_j, 0.0);
                assert_eq!(half.drag_heating_j, 0.0);
                assert_eq!(half.external_heat_to_gas_j, 0.0);
            }
            close(
                gas.transport.left.right_boundary_change.energy
                    + gas.transport.right.left_boundary_change.energy,
                -report.midpoint_velocity_m_s * report.gas_impulse_n_s,
                2e-14,
            );
            assert_eq!(left.cell(0).unwrap().partial_densities[1], 0.0);
            assert_eq!(left.time(), right.time());
            assert_eq!(left.mesh().positions()[2], right.mesh().positions()[0]);
            old_t = [
                end_temperature(&left, 1, 0.0),
                end_temperature(&right, 0, 0.0),
            ];
            old_wall = end_wall;
        }
        let sign = (initial_wall - 1.0).signum();
        assert!(
            sign * piston.velocity_m_s() > 0.0,
            "heat must drive the real pressure root"
        );
        assert!(sign * (pressure(&left, 1) - pressure(&right, 0)) > 0.0);
        close(
            sum([
                left.totals().energy - original[0].energy,
                right.totals().energy - original[1].energy,
                0.6 * (wall_temperature(&wall) - initial_wall),
                0.5 * piston.mass_kg() * piston.velocity_m_s().powi(2),
            ]),
            0.0,
            1.3e-10,
        );
    }
}

#[test]
fn g0_zero_conductance_retains_native_pressure_work_and_stiff_shared_heat_is_implicit() {
    let (mut left, mut right) = fields([2, 2], 0.0);
    let (mut native_left, mut native_right) = fields([2, 2], 0.0);
    let mut wall = finite(&left, &right, 3.0, [0.0; 2]);
    let mut piston = Piston::new(0.4, 0.2).unwrap();
    let mut native_piston = piston;
    let native = native_left
        .step_with_shared_piston(
            &mut native_right,
            &mut native_piston,
            Time::new(0.01),
            StepControl::default(),
            0.0,
            controls(),
            || false,
        )
        .unwrap();
    let report = advance(&mut left, &mut right, &mut wall, &mut piston, 0.01, || {
        false
    })
    .unwrap();
    assert_eq!(report.gas.transport, native.gas);
    assert_eq!(left.states, native_left.states);
    assert_eq!(right.states, native_right.states);
    assert_eq!(piston, native_piston);
    assert_eq!(wall_temperature(&wall), 3.0);
    assert_eq!(report.gas.heat.heat_to_gas_j, [0.0; 2]);

    let (mut left, mut right) = fields([2, 2], 0.0);
    let conductance = [1e5, 2e5];
    let mut wall = finite(&left, &right, 3.0, conductance);
    let mut piston = Piston::new(0.4, 0.0).unwrap();
    let mut mechanical = controls();
    mechanical.energy_tolerance_j = 1e-8;
    let report = left
        .step_with_shared_piston_and_heat(
            &mut right,
            &mut wall,
            &mut piston,
            Time::new(0.01),
            StepControl::default(),
            WallControl::default(),
            0.0,
            mechanical,
            || false,
        )
        .unwrap();
    assert_eq!(report.gas.restarts, 0);
    assert_eq!(report.gas.transport.left.elapsed.value(), 0.01);
    let (middle_wall, heat) = linear_heat(
        [1.0, 2.0],
        [2.5 * AREA / 2.0, 1.75 * AREA / 2.0],
        3.0,
        conductance,
        0.005,
    );
    close(
        report.gas.transport.left.before.energy - report.gas.before.left.energy,
        heat[0],
        2e-9,
    );
    close(
        report.gas.transport.right.before.energy - report.gas.before.right.energy,
        heat[1],
        2e-9,
    );
    close(
        report.gas.heat.stages[0].wall_energy_change_j,
        0.6 * (middle_wall - 3.0),
        2e-9,
    );
    assert!((1.0..3.0).contains(&wall_temperature(&wall)));
    assert!(report.admitted(mechanical));
}

// A declared original constitutive test source, not measured material data.
// u(T)=offset+a*T²/2, so its capacity changes with actual wall temperature.
struct QuadraticStore {
    offset: f64,
    basis: ThermalStorageBasis,
}
impl ThermalLaw for QuadraticStore {
    fn storage_basis(&self) -> ThermalStorageBasis {
        self.basis
    }
    fn temperature_k(&self, u: f64) -> Result<f64, DepositionError> {
        if u <= self.offset || !u.is_finite() {
            return Err(DepositionError::MaterialDomain);
        }
        Ok((5.0 * (u - self.offset)).sqrt())
    }
    fn maximum_temperature_slope(&self) -> f64 {
        100.0
    }
    fn validate_cell(&self, cell: &DepositedCell) -> Result<(), DepositionError> {
        if cell.mass_kg != 0.4 || cell.geometry.volume_m3 != 0.1 {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(())
    }
}

fn source_storage(offset: f64, basis: ThermalStorageBasis) -> SharedWallStorage<'static> {
    let source = MaterialWallCaloric::from_thermal_law(
        QuadraticStore { offset, basis },
        DepositedCell {
            id: CellId(17),
            segment_id: 3,
            body_id: 8,
            material_id: MaterialId(29),
            mass_kg: 0.4,
            specific_h_j_kg: offset + 1.8,
            geometry: CellGeometry {
                start_m: [0.0; 3],
                end_m: [1.0, 0.0, 0.0],
                volume_m3: 0.1,
                surface_area_m2: 1.0,
            },
            boundary: ThermalBoundary {
                ambient_k: 1.0,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w: 0.0,
            },
        },
        [offset + 0.01, offset + 20.0],
    )
    .unwrap();
    SharedWallStorage {
        thermal: source.initial_contact(0.0, WallDrag::None).unwrap().thermal,
        heat_storage: Some(source),
    }
}

#[test]
fn g1_original_nonlinear_wall_storage_closes_without_a_tangent_capacity_or_hidden_enthalpy_work() {
    let mut results = Vec::new();
    for datum in [0.0, 1024.0] {
        let (mut left, mut right) = fields([2, 2], datum * 1e12);
        let g = [1.7, 0.4];
        let mut wall = SharedPistonThermal::new(
            &left,
            &right,
            source_storage(datum, ThermalStorageBasis::InternalEnergy),
            g,
            SharedWallControl::default(),
            || false,
        )
        .unwrap();
        let mut piston = Piston::new(0.4, 0.0).unwrap();
        let initial = wall.storage().thermal;
        let report = advance(&mut left, &mut right, &mut wall, &mut piston, 0.01, || {
            false
        })
        .unwrap();
        let hdt = 0.5 * report.gas.transport.left.elapsed.value();
        let capacity = [2.5 * AREA / 2.0, 1.75 * AREA / 2.0];
        let b = [
            capacity[0] * g[0] * hdt / (capacity[0] + g[0] * hdt),
            capacity[1] * g[1] * hdt / (capacity[1] + g[1] * hdt),
        ];
        let a = 0.08; // mass*a_source/2, counted once for the shared solid.
        let c = -a * 9.0 - b[0] - 2.0 * b[1];
        let middle_wall = -2.0 * c / (b[0] + b[1] + ((b[0] + b[1]).powi(2) - 4.0 * a * c).sqrt());
        close(
            report.gas.transport.left.before.energy - report.gas.before.left.energy,
            b[0] * (middle_wall - 1.0),
            8e-12,
        );
        close(
            report.gas.transport.right.before.energy - report.gas.before.right.energy,
            b[1] * (middle_wall - 2.0),
            8e-12,
        );
        close(
            report.gas.heat.stages[0].wall_energy_change_j,
            a * (middle_wall.powi(2) - 9.0),
            8e-12,
        );
        let source = wall.storage().heat_storage.as_ref().unwrap();
        assert_eq!(source.initial_cell().unwrap().material_id, MaterialId(29));
        close(
            report.gas.heat.wall_energy_change_j,
            source
                .energy_change_j(
                    WallContact {
                        drag: WallDrag::None,
                        thermal: initial,
                    },
                    WallContact {
                        drag: WallDrag::None,
                        thermal: wall.storage().thermal,
                    },
                )
                .unwrap(),
            0.0,
        );
        assert!(report.admitted(controls()));
        results.push((
            left.states,
            right.states,
            piston.velocity_m_s(),
            wall_temperature(&wall),
        ));
    }
    for (a, b) in results[0]
        .0
        .iter()
        .chain(&results[0].1)
        .zip(results[1].0.iter().chain(&results[1].1))
    {
        close(*a, *b, 1e-11);
    }
    close(results[0].2, results[1].2, 1e-11);
    close(results[0].3, results[1].3, 1e-11);
    let (left, right) = fields([2, 2], 0.0);
    assert!(
        SharedPistonThermal::new(
            &left,
            &right,
            source_storage(0.0, ThermalStorageBasis::Enthalpy),
            [1.0; 2],
            SharedWallControl::default(),
            || false
        )
        .is_err()
    );
}

#[test]
fn g4_heat_restarts_on_acoustic_shortening_and_late_cancellation_preserves_all_accepted_state() {
    let (mut left, mut right) = fields([12, 20], 0.0);
    let mut wall = finite(&left, &right, 10.0, [100.0; 2]);
    let mut piston = Piston::new(0.4, 0.0).unwrap();
    let mut polls = 0;
    let report = advance(&mut left, &mut right, &mut wall, &mut piston, 1.0, || {
        polls += 1;
        false
    })
    .unwrap();
    assert!(report.gas.restarts > 0);
    let (mut replay_left, mut replay_right) = fields([12, 20], 0.0);
    let mut replay_wall = finite(&replay_left, &replay_right, 10.0, [100.0; 2]);
    let mut replay_piston = Piston::new(0.4, 0.0).unwrap();
    let replay = advance(
        &mut replay_left,
        &mut replay_right,
        &mut replay_wall,
        &mut replay_piston,
        report.gas.transport.left.elapsed.value(),
        || false,
    )
    .unwrap();
    assert_eq!(left.states, replay_left.states);
    assert_eq!(right.states, replay_right.states);
    assert_eq!(piston, replay_piston);
    assert_eq!(wall.storage(), replay_wall.storage());
    assert_eq!(report.gas.heat, replay.gas.heat);

    let (mut failed_left, mut failed_right) = fields([12, 20], 0.0);
    let mut failed_wall = finite(&failed_left, &failed_right, 10.0, [100.0; 2]);
    let mut failed_piston = Piston::new(0.4, 0.0).unwrap();
    let initial_piston = failed_piston;
    let states = [failed_left.states.clone(), failed_right.states.clone()];
    let meshes = [failed_left.mesh.clone(), failed_right.mesh.clone()];
    let storage = failed_wall.storage().clone();
    let mut calls = 0;
    let error = advance(
        &mut failed_left,
        &mut failed_right,
        &mut failed_wall,
        &mut failed_piston,
        1.0,
        || {
            calls += 1;
            calls == polls
        },
    )
    .unwrap_err();
    assert_eq!(error, SharedPistonThermalError::Cancelled);
    assert_eq!(failed_left.states, states[0]);
    assert_eq!(failed_right.states, states[1]);
    assert_eq!(failed_left.mesh, meshes[0]);
    assert_eq!(failed_right.mesh, meshes[1]);
    assert_eq!(failed_left.time().value(), 0.0);
    assert_eq!(failed_right.time().value(), 0.0);
    assert_eq!(failed_wall.storage(), &storage);
    assert_eq!(failed_wall.last_report(), None);
    assert_eq!(failed_piston, initial_piston);
    let retried = advance(
        &mut failed_left,
        &mut failed_right,
        &mut failed_wall,
        &mut failed_piston,
        1.0,
        || false,
    )
    .unwrap();
    assert_eq!(failed_left.states, left.states);
    assert_eq!(failed_right.states, right.states);
    assert_eq!(failed_wall.storage(), wall.storage());
    assert_eq!(retried, report);

    let (mut left, mut right) = fields([2, 2], 0.0);
    let mut wall = SharedPistonThermal::new(
        &left,
        &right,
        storage,
        [1.0; 2],
        SharedWallControl {
            maximum_source_evaluations: 1,
        },
        || false,
    )
    .unwrap();
    let original = [left.states.clone(), right.states.clone()];
    let mut piston = initial_piston;
    assert!(
        left.step_with_shared_piston_and_heat(
            &mut right,
            &mut wall,
            &mut piston,
            Time::new(0.01),
            StepControl {
                maximum_halvings: 0,
                ..StepControl::default()
            },
            heat_controls(),
            0.0,
            controls(),
            || false
        )
        .is_err()
    );
    assert_eq!(left.states, original[0]);
    assert_eq!(right.states, original[1]);
    assert_eq!(piston, initial_piston);
    assert_eq!(wall.last_report(), None);
}

fn contacts() -> (MixtureEuler1d, MixtureEuler1d) {
    let fields: Vec<_> = [64, 80]
        .into_iter()
        .enumerate()
        .map(|(side, count)| {
            let g = gas(vec![
                source("A", R, 3.5, 0.4, 1e18),
                source("B", 2.0 * R, 4.5, 0.2, -2e18),
            ]);
            let start = side as f64;
            let positions: Vec<_> = (0..=count)
                .map(|i| start + i as f64 / count as f64)
                .collect();
            let volumes = positions.windows(2).map(|x| AREA * (x[1] - x[0])).collect();
            let mesh = Mesh1d::new(positions, vec![AREA; count + 1], volumes).unwrap();
            let states = (0..count)
                .map(|i| {
                    let d = if i < count / 2 {
                        [1.0, 0.0]
                    } else {
                        [0.0, 1.0]
                    };
                    g.conserved(&d, 0.0, [1.2, 0.8][side]).unwrap()
                })
                .collect();
            MixtureEuler1d::new_on_mesh(g, mesh, states, CLOSED, || false).unwrap()
        })
        .collect();
    let mut it = fields.into_iter();
    (it.next().unwrap(), it.next().unwrap())
}

#[test]
fn g0_heat_changes_wall_pressure_while_material_contacts_keep_their_independent_pressure() {
    let (mut left, mut right) = contacts();
    let (mut native_left, mut native_right) = contacts();
    let mut wall = finite(&left, &right, 3.0, [1.7, 0.4]);
    let mut native_wall = finite(&native_left, &native_right, 3.0, [1.7, 0.4]);
    let mut piston = Piston::new(0.4, 0.25).unwrap();
    let mut native_piston = piston;
    let mut workspaces = [64, 80].map(|n| {
        ContactWorkspace::new(
            n,
            ContactWorkspace::required_workspace_bytes(n).unwrap(),
            || false,
        )
        .unwrap()
    });
    let equilibrium = [ContactEquilibriumControl {
        relative_tolerance: 1e-10,
        maximum_energy_defect_j: 1.0,
    }; 2];
    let mut defect = 0.0;
    for _ in 0..8 {
        advance(
            &mut native_left,
            &mut native_right,
            &mut native_wall,
            &mut native_piston,
            0.001,
            || false,
        )
        .unwrap();
        let [a, b] = &mut workspaces;
        let report = left
            .step_with_shared_piston_heat_and_contact_equilibrium(
                &mut right,
                &mut wall,
                &mut piston,
                Time::new(0.001),
                StepControl::default(),
                heat_controls(),
                0.0,
                controls(),
                equilibrium,
                [a, b],
                || false,
            )
            .unwrap();
        let correction = sum([
            report.gas.transport.left.correction.signed_energy_j,
            report.gas.transport.right.correction.signed_energy_j,
        ]);
        defect += correction.abs();
        close(
            sum([
                report.gas.after.left.energy - report.gas.before.left.energy,
                report.gas.after.right.energy - report.gas.before.right.energy,
                report.gas.heat.wall_energy_change_j,
                report.kinetic_energy_change_j,
                -correction,
            ]),
            report.total_energy_residual_j,
            3e-14,
        );
        assert!(report.admitted(controls()));
        assert!(report.gas.heat.heat_to_gas_j[0] > 0.0);
    }
    assert!(defect > 1e-9);
    for (i, flow, native) in [(0, &left, &native_left), (1, &right, &native_right)] {
        let mut native_error = 0.0_f64;
        for cell in flow.cell_count() / 2 - 3..flow.cell_count() / 2 + 3 {
            close(pressure(flow, cell), [1.2, 0.8][i], 3e-10);
            native_error = native_error.max((pressure(native, cell) - [1.2, 0.8][i]).abs());
        }
        assert!(native_error > 1e-5);
    }
    assert!(
        (pressure(&left, left.cell_count() - 1) - 1.2).abs() > 1e-4,
        "physical heating/compression at the piston must not be projected away"
    );
}

#[test]
fn g2_heat_driven_pressure_and_piston_velocity_converge_under_time_refinement() {
    let mut endpoints = Vec::new();
    for steps in [2, 4, 8, 32] {
        let (mut left, mut right) = fields([2, 2], 0.0);
        let mut wall = finite(&left, &right, 3.0, [30.0, 7.0]);
        let mut piston = Piston::new(0.4, 0.0).unwrap();
        let dt = 0.04 / steps as f64;
        for _ in 0..steps {
            let report =
                advance(&mut left, &mut right, &mut wall, &mut piston, dt, || false).unwrap();
            assert_eq!(report.gas.transport.left.elapsed.value(), dt);
            assert_eq!(report.gas.restarts, 0);
        }
        close(left.time().value(), 0.04, 2e-16);
        endpoints.push([
            piston.velocity_m_s(),
            pressure(&left, 1),
            pressure(&right, 0),
        ]);
    }
    // This compares the composed time discretization at fixed spatial
    // resolution; it neither claims second order nor validates a material.
    for component in 0..3 {
        let error = |i: usize| (endpoints[i][component] - endpoints[3][component]).abs();
        assert!(
            error(1) < 0.8 * error(0),
            "component {component}: {endpoints:?}"
        );
        assert!(
            error(2) < 0.8 * error(1),
            "component {component}: {endpoints:?}"
        );
    }
}
