//! G1 original shared-wall calorimetry and G4 coupled network publication.
use super::*;
use crate::multicomponent::network::dissipating::{DissipatingNetwork1d, DissipatingStepControl};
use crate::multicomponent::network::dissipation::tests::law;
use crate::multicomponent::network::wall_exchange::{
    SharedWall, SharedWallControl, SharedWallGas, SharedWallLink, SharedWallPort, SharedWallStorage,
};
use crate::multicomponent::wall_exchange::{
    MaterialWallCaloric, SharedRecoilControl, SlidingSupport,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionError, MaterialId, ThermalLaw,
    ThermalStorageBasis,
};

const WALL_MASS: f64 = 0.6;

struct Storage {
    basis: ThermalStorageBasis,
    plateau: bool,
}
impl ThermalLaw for Storage {
    fn storage_basis(&self) -> ThermalStorageBasis {
        self.basis
    }
    fn temperature_k(&self, u: f64) -> Result<f64, DepositionError> {
        if !u.is_finite() || u <= 0.0 {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if !self.plateau {
            // Original manufactured u(T)=0.2*T², hence Cv(T)=0.4*T.
            (5.0 * u).sqrt()
        } else if u < 3.0 {
            2.0 + u
        } else if u <= 8.0 {
            5.0
        } else {
            5.0 + 0.5 * (u - 8.0)
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        100.0
    }
}

fn source(
    basis: ThermalStorageBasis,
    plateau: bool,
    range: [f64; 2],
    id: u64,
) -> MaterialWallCaloric<'static> {
    source_amount(basis, plateau, range, id, WALL_MASS)
}
fn source_amount(
    basis: ThermalStorageBasis,
    plateau: bool,
    range: [f64; 2],
    id: u64,
    mass: f64,
) -> MaterialWallCaloric<'static> {
    MaterialWallCaloric::from_thermal_law(
        Storage { basis, plateau },
        DepositedCell {
            id: CellId(id),
            segment_id: id,
            body_id: id,
            material_id: MaterialId(29),
            // This mass belongs to the common support inertia, once.
            mass_kg: mass,
            specific_h_j_kg: if plateau { 6.0 } else { 5.0 },
            geometry: CellGeometry {
                start_m: [0.0; 3],
                end_m: [1.0, 0.0, 0.0],
                volume_m3: 0.06,
                surface_area_m2: 1.0,
            },
            boundary: Boundary {
                ambient_k: 2.0,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w: 0.0,
            },
        },
        range,
    )
    .unwrap()
}

fn support() -> SlidingSupport {
    SlidingSupport {
        mass_kg: WALL_MASS,
        velocity_m_s: 0.5,
    }
}
fn base(
    n: &MixtureNetwork1d,
    supports: Vec<Option<SlidingSupport>>,
) -> NetworkWallOperator<'static> {
    uniform_operator(
        n,
        WallContact {
            drag: WallDrag::Linear { rate_per_s: 3.0 },
            thermal: WallThermal::Adiabatic,
        },
    )
    .with_shared_branch_supports(supports, SharedRecoilControl::default(), 1e-11, || false)
    .unwrap()
}
fn group(
    source: MaterialWallCaloric<'static>,
    branches: &[usize],
    films: [f64; 2],
) -> SharedWall<'static> {
    SharedWall {
        storage: SharedWallStorage {
            thermal: source.initial_contact(0.0, WallDrag::None).unwrap().thermal,
            heat_storage: Some(source),
        },
        ports: branches
            .iter()
            .flat_map(|&branch| {
                (0..2).map(move |cell| SharedWallPort {
                    gas: SharedWallGas::Branch { branch, cell },
                    conductance_w_per_k: films[cell],
                })
            })
            .collect(),
    }
}
fn bind(
    n: &MixtureNetwork1d,
    source: MaterialWallCaloric<'static>,
    films: [f64; 2],
    fractions: Option<[f64; 2]>,
    limit: usize,
) -> Result<NetworkWallOperator<'static>, WallExchangeError> {
    let wall = base(n, vec![Some(support())]).with_shared_walls(
        n,
        vec![group(source, &[0], films)],
        SharedWallControl {
            maximum_source_evaluations: limit,
        },
        || false,
    )?;
    match fractions {
        Some(f) => wall.with_shared_branch_heat_partition(vec![f.to_vec()], || false),
        None => Ok(wall),
    }
}
fn stored(wall: &NetworkWallOperator<'_>, index: usize) -> f64 {
    let WallThermal::Caloric {
        mass_kg,
        specific_energy_j_kg,
        ..
    } = wall.shared_wall(index).unwrap().storage.thermal
    else {
        panic!("original material energy coordinate");
    };
    mass_kg * specific_energy_j_kg
}
fn temperature(wall: &NetworkWallOperator<'_>, index: usize) -> f64 {
    let state = &wall.shared_wall(index).unwrap().storage;
    state
        .heat_storage
        .as_ref()
        .unwrap()
        .temperature_k(state.thermal)
        .unwrap()
}

// Equal initial velocities reduce the gas to m=0.4 kg against M=0.6 kg.
// Its slip obeys s'= -3*(1+m/M)*s; total thermalized loss is mu*(s0²-s²)/2.
fn slip_loss(dt: f64) -> f64 {
    0.5 * (0.4 * WALL_MASS / (0.4 + WALL_MASS)) * 2.25 * -(-10.0 * dt).exp_m1()
}
fn quadratic(a: f64, b: f64, c: f64) -> f64 {
    -2.0 * c / (b + (b * b - 4.0 * a * c).sqrt())
}
fn thermal_oracle(fractions: [f64; 2], films: [f64; 2], dt: f64) -> ([f64; 2], f64) {
    let loss = slip_loss(dt);
    let wall_heat = loss * (1.0 - 0.5 * (fractions[0] + fractions[1]));
    // The original NASA fixture has Cv=2.5+0.1*T, with 0.2 kg per cell.
    let gases = |tw| {
        core::array::from_fn(|i| {
            quadratic(
                0.01,
                0.5 + dt * films[i],
                -1.04 - 0.5 * fractions[i] * loss - dt * films[i] * tw,
            )
        })
    };
    let hot_wall = (25.0 + wall_heat / (0.2 * WALL_MASS)).sqrt();
    let mut low = 2.0_f64;
    let mut high = hot_wall.max(gases(hot_wall)[0]).max(gases(hot_wall)[1]);
    // Independent analytic storage equations; no production gas or thermal
    // query participates in this monotone one-wall reference calculation.
    for _ in 0..80 {
        let mid = 0.5 * (low + high);
        let gas = gases(mid);
        let residual = 0.2 * WALL_MASS * (mid * mid - 25.0) - wall_heat
            + dt * (films[0] * (mid - gas[0]) + films[1] * (mid - gas[1]));
        if residual < 0.0 {
            low = mid;
        } else {
            high = mid;
        }
    }
    let wall = 0.5 * (low + high);
    (gases(wall), wall)
}

#[test]
fn g1_shared_partition_matches_original_variable_capacity_and_one_wall_inventory() {
    let g = gas(1, true, 1e6);
    let dt = 0.2_f64;
    let slip = 1.5 * (-5.0 * dt).exp();
    let gas_velocity = 1.1 + 0.6 * slip;
    let wall_velocity = 1.1 - 0.4 * slip;
    let kinetic = 0.5 * WALL_MASS * (wall_velocity - 0.5) * (wall_velocity + 0.5);
    for fractions in [[0.0; 2], [0.35; 2], [1.0; 2], [0.2, 0.8]] {
        for film in [0.0, 0.3] {
            let films = [film, 2.0 * film];
            let mut n = network(&g, &[2.0], false);
            let before = fields(&n);
            let original = source(ThermalStorageBasis::InternalEnergy, false, [0.01, 50.0], 17);
            let mut wall = bind(&n, original.clone(), films, Some(fractions), 65_536).unwrap();
            let report = wall
                .advance_for(&mut n, Time::new(dt), controls(), || false)
                .unwrap_or_else(|error| {
                    panic!("fractions={fractions:?}, films={films:?}: {error:?}")
                });
            let (tg, tw) = thermal_oracle(fractions, films, dt);
            let loss = slip_loss(dt);
            let wall_heat = loss * (1.0 - 0.5 * (fractions[0] + fractions[1]));
            let mut heat = 0.0;
            let mut gas_change = 0.0;
            for i in 0..2 {
                let cell = n.branches()[0].cell(i).unwrap();
                let actual = g.primitive(cell, || false).unwrap();
                close(actual.velocity, gas_velocity, 3e-12);
                close(actual.temperature, tg[i], 3e-10);
                assert_eq!(cell.partial_densities, before.0[0][i].partial_densities);
                let q = dt * films[i] * (tw - tg[i]);
                let sensible = 0.2 * (tg[i] - 2.0) * (2.5 + 0.05 * (tg[i] + 2.0));
                close(sensible, 0.5 * fractions[i] * loss + q, 3e-12);
                let change = 0.2 * (cell.energy - before.0[0][i].energy);
                close(
                    change,
                    0.1 * (gas_velocity - 2.0) * (gas_velocity + 2.0) + sensible,
                    3e-11,
                );
                gas_change += change;
                heat += q;
                assert_eq!(
                    wall.branch_contacts(0).unwrap()[i].thermal,
                    WallThermal::Adiabatic
                );
            }
            let receipt = wall.accepted_branch_shared_support_exchange(0).unwrap();
            close(
                wall.branch_shared_support(0).unwrap().velocity_m_s,
                wall_velocity,
                3e-12,
            );
            close(receipt.drag_heating_j, loss, 3e-12);
            close(receipt.wall_slip_heating_j, wall_heat, 3e-12);
            close(receipt.gas_energy_change_j, -kinetic - wall_heat, 3e-12);
            close(receipt.kinetic_energy_change_j, kinetic, 3e-12);
            close(receipt.energy_residual_j, 0.0, 3e-12);
            close(temperature(&wall, 0), tw, 3e-10);
            close(stored(&wall, 0), 0.2 * WALL_MASS * tw * tw, 3e-11);
            close(
                gas_change + kinetic + stored(&wall, 0) - WALL_MASS * 5.0,
                0.0,
                4e-11,
            );
            close(report.wall_energy_change_j, wall_heat - heat, 4e-11);
            close(report.heat_to_gas_j, heat - wall_heat, 4e-11);
            close(report.mechanical_work_to_gas_j, -kinetic, 3e-12);
            close(report.energy_residual_j, 0.0, 4e-11);
            assert_eq!(report.external_heat_to_gas_j, 0.0);
            assert!(wall.shared_wall(1).is_none());
            assert_eq!(n.time(), Time::new(0.0));
            if fractions == [1.0; 2] {
                let mut reference = network(&g, &[2.0], false);
                let mut default = bind(&reference, original, films, None, 65_536).unwrap();
                assert_eq!(
                    default
                        .advance_for(&mut reference, Time::new(dt), controls(), || false)
                        .unwrap(),
                    report
                );
                assert_eq!(fields(&reference), fields(&n));
                assert_eq!(
                    default.shared_wall(0).unwrap().storage.thermal,
                    wall.shared_wall(0).unwrap().storage.thermal
                );
                assert_eq!(
                    default.accepted_branch_shared_support_exchange(0),
                    wall.accepted_branch_shared_support_exchange(0)
                );
                assert_eq!(
                    default.last_shared_wall_report(),
                    wall.last_shared_wall_report()
                );
            }
        }
    }
}

#[test]
fn g1_latent_inventory_and_shared_inertia_continue_through_both_network_source_halves() {
    let g = gas(1, true, 0.0);
    for film in [0.0, 0.3] {
        let n = network(&g, &[2.0], false);
        let original = source(ThermalStorageBasis::InternalEnergy, true, [0.01, 50.0], 17);
        let mut wall = bind(&n, original, [film; 2], Some([0.35; 2]), 65_536).unwrap();
        let mut flow =
            DissipatingNetwork1d::new(n, law(&g, 0.004, 0.0, 0.02, 0.0), vec![], || false).unwrap();
        let initial_gas = flow.totals().energy;
        let initial_wall = stored(&wall, 0);
        let mut assigned_heat = 0.0;
        for _ in 0..5 {
            let report = flow
                .step_with_wall(
                    &mut wall,
                    Time::new(0.01),
                    DissipatingStepControl::default(),
                    controls(),
                    || false,
                )
                .unwrap();
            let halves = report.wall.unwrap();
            assert_eq!(halves.first.elapsed.value(), 0.5 * report.elapsed.value());
            assert_eq!(halves.second.elapsed.value(), 0.5 * report.elapsed.value());
            let velocity = wall.branch_shared_support(0).unwrap().velocity_m_s;
            let kinetic = 0.5 * WALL_MASS * (velocity - 0.5) * (velocity + 0.5);
            close(
                flow.totals().energy - initial_gas + stored(&wall, 0) - initial_wall + kinetic,
                0.0,
                7e-11,
            );
            assert_eq!(temperature(&wall, 0), 5.0);
            close(report.absolute_energy_residual_j, 0.0, 5e-11);
            let receipt = wall.accepted_branch_shared_support_exchange(0).unwrap();
            assigned_heat += receipt.wall_slip_heating_j;
            close(receipt.energy_residual_j, 0.0, 3e-12);
            assert_eq!(wall.last_report().unwrap().external_heat_to_gas_j, 0.0);
        }
        close(assigned_heat, 0.65 * slip_loss(flow.time().value()), 3e-11);
        assert!((stored(&wall, 0) - initial_wall).abs() > 1e-4);
        if film == 0.0 {
            close(stored(&wall, 0) - initial_wall, assigned_heat, 3e-11);
        }
    }
}

#[test]
fn g4_shared_heat_refusal_preserves_accepted_state_and_exact_retry() {
    let g = gas(1, true, 0.0);
    let snapshot = |wall: &NetworkWallOperator<'_>| {
        (
            wall.shared_wall(0).unwrap().storage.thermal,
            wall.branch_shared_support(0),
            wall.accepted_branch_shared_support_exchange(0),
            wall.last_report(),
            wall.last_shared_wall_report(),
            wall.branch_contacts(0).unwrap().to_vec(),
        )
    };
    for mode in 0..3 {
        let first_dt = if mode == 2 { 0.0 } else { 1e-8 };
        let make = || {
            let mut n = network(&g, &[2.0], false);
            let range = if mode == 1 {
                [4.99999, 5.00001]
            } else {
                [0.01, 50.0]
            };
            let original = source(ThermalStorageBasis::InternalEnergy, false, range, 17);
            let mut wall = bind(
                &n,
                original,
                [0.3; 2],
                Some([0.35; 2]),
                if mode == 2 { 64 } else { 65_536 },
            )
            .unwrap();
            wall.advance_for(&mut n, Time::new(first_dt), controls(), || false)
                .unwrap();
            (n, wall)
        };
        let (mut n, mut wall) = make();
        let (mut expected, mut reference) = make();
        let retry_dt = if mode == 0 { 0.2 } else { first_dt };
        let mut polls = 0;
        let expected_report = reference
            .advance_for(&mut expected, Time::new(retry_dt), controls(), || {
                polls += 1;
                false
            })
            .unwrap();
        let before = fields(&n);
        let accepted = snapshot(&wall);
        let mut visits = 0;
        let refused = wall.advance_for(&mut n, Time::new(0.2), controls(), || {
            visits += 1;
            mode == 0 && visits == polls
        });
        match mode {
            0 => assert_eq!(refused, Err(WallExchangeError::Cancelled)),
            1 => assert!(refused.is_err()),
            _ => assert_eq!(
                refused,
                Err(WallExchangeError::Invalid(
                    "shared-wall nested source work budget"
                ))
            ),
        }
        assert_eq!(fields(&n), before);
        assert_eq!(snapshot(&wall), accepted);
        assert_eq!(
            wall.advance_for(&mut n, Time::new(retry_dt), controls(), || false)
                .unwrap(),
            expected_report
        );
        assert_eq!(fields(&n), fields(&expected));
        assert_eq!(snapshot(&wall), snapshot(&reference));
    }
}

#[test]
fn g0_shared_partition_admits_one_original_thermal_owner_for_each_heated_support() {
    let g = gas(1, true, 0.0);
    let n = network(&g, &[2.0], false);
    let unpartitioned = |basis| {
        bind(
            &n,
            source(basis, false, [0.01, 50.0], 17),
            [0.3; 2],
            None,
            65_536,
        )
        .unwrap()
    };
    for fractions in [
        vec![],
        vec![vec![0.35]],
        vec![vec![-0.1; 2]],
        vec![vec![1.1; 2]],
        vec![vec![f64::NAN; 2]],
        vec![vec![f64::INFINITY; 2]],
    ] {
        assert!(
            unpartitioned(ThermalStorageBasis::InternalEnergy)
                .with_shared_branch_heat_partition(fractions, || false)
                .is_err()
        );
    }
    for basis in [
        ThermalStorageBasis::Enthalpy,
        ThermalStorageBasis::DeclaredSensibleStorage,
    ] {
        assert!(
            unpartitioned(basis)
                .with_shared_branch_heat_partition(vec![vec![0.35; 2]], || false)
                .is_err()
        );
    }
    assert!(
        base(&n, vec![Some(support())])
            .with_shared_branch_heat_partition(vec![vec![0.35; 2]], || false)
            .is_err()
    );
    assert!(
        base(&n, vec![Some(support())])
            .with_shared_branch_heat_partition(vec![vec![1.0; 2]], || false)
            .is_ok()
    );
    let unsupported = base(&n, vec![None])
        .with_shared_walls(
            &n,
            vec![group(
                source(ThermalStorageBasis::InternalEnergy, false, [0.01, 50.0], 17),
                &[0],
                [0.3; 2],
            )],
            SharedWallControl::default(),
            || false,
        )
        .unwrap();
    assert!(
        unsupported
            .with_shared_branch_heat_partition(vec![vec![0.35; 2]], || false)
            .is_err()
    );
    assert!(
        unpartitioned(ThermalStorageBasis::InternalEnergy)
            .with_shared_branch_heat_partition(vec![vec![0.35; 2]], || false)
            .unwrap()
            .with_shared_branch_heat_partition(vec![vec![0.5; 2]], || false)
            .is_err()
    );

    let mut two = network(&g, &[2.0, -0.4], false);
    let shared = base(&two, vec![Some(support()); 2])
        .with_shared_walls(
            &two,
            vec![group(
                source(ThermalStorageBasis::InternalEnergy, false, [0.01, 50.0], 17),
                &[0, 1],
                [0.3; 2],
            )],
            SharedWallControl::default(),
            || false,
        )
        .unwrap();
    assert!(
        shared
            .with_shared_branch_heat_partition(vec![vec![0.35; 2]; 2], || false)
            .is_err()
    );
    let mut separate = base(&two, vec![Some(support()); 2])
        .with_shared_walls(
            &two,
            (0..2)
                .map(|branch| {
                    group(
                        source(
                            ThermalStorageBasis::InternalEnergy,
                            false,
                            [0.01, 50.0],
                            17 + branch as u64,
                        ),
                        &[branch],
                        [0.3; 2],
                    )
                })
                .collect(),
            SharedWallControl::default(),
            || false,
        )
        .unwrap()
        .with_shared_branch_heat_partition(vec![vec![0.35; 2]; 2], || false)
        .unwrap();
    separate
        .advance_for(&mut two, Time::new(0.1), controls(), || false)
        .unwrap();
    assert!(
        separate
            .accepted_branch_shared_support_exchange(0)
            .unwrap()
            .wall_slip_heating_j
            > 0.0
    );
    assert!(
        separate
            .accepted_branch_shared_support_exchange(1)
            .unwrap()
            .wall_slip_heating_j
            > 0.0
    );

    // One mechanical body may own distinct thermal regions. Their masses
    // sum to the one support inertia, and an existing link transfers heat
    // between the actual regions AFTER unequal frictional deposition.
    let mut regional = network(&g, &[2.0], false);
    let before_gas = regional.totals().energy;
    let regions = (0..2)
        .map(|cell| {
            let mut region = group(
                source_amount(
                    ThermalStorageBasis::InternalEnergy,
                    false,
                    [0.01, 50.0],
                    31 + cell as u64,
                    0.5 * WALL_MASS,
                ),
                &[0],
                [0.0; 2],
            );
            region.ports = vec![region.ports[cell]];
            region
        })
        .collect();
    let mut linked = base(&regional, vec![Some(support())])
        .with_shared_walls(&regional, regions, SharedWallControl::default(), || false)
        .unwrap()
        .with_shared_branch_heat_partition(vec![vec![0.0, 1.0]], || false)
        .unwrap()
        .with_shared_wall_links(
            vec![SharedWallLink {
                walls: [0, 1],
                conductance_w_per_k: 1.0,
            }],
            || false,
        )
        .unwrap();
    let report = linked
        .advance_for(&mut regional, Time::new(0.2), controls(), || false)
        .unwrap();
    let transferred = linked.shared_wall_link_reports().unwrap()[0].heat_a_to_b_j;
    let assigned = 0.5 * slip_loss(0.2);
    assert!(transferred > 0.0);
    close(
        transferred,
        0.2 * (temperature(&linked, 0) - temperature(&linked, 1)),
        3e-11,
    );
    close(
        stored(&linked, 0),
        0.5 * WALL_MASS * 5.0 + assigned - transferred,
        3e-11,
    );
    close(
        stored(&linked, 1),
        0.5 * WALL_MASS * 5.0 + transferred,
        3e-11,
    );
    close(report.wall_energy_change_j, assigned, 4e-11);
    close(report.heat_to_gas_j, -assigned, 4e-11);
    let receipt = linked.accepted_branch_shared_support_exchange(0).unwrap();
    close(
        regional.totals().energy - before_gas + stored(&linked, 0) + stored(&linked, 1)
            - WALL_MASS * 5.0
            + receipt.kinetic_energy_change_j,
        0.0,
        5e-11,
    );
}
