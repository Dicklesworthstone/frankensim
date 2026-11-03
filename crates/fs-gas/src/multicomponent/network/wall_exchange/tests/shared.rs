//! G1 source calorimetry and G3/G4 continuing composition with one wall inventory.
use super::*;
use crate::chemistry::IsochoricKineticsSolve;
use crate::multicomponent::network::dissipation::tests::owned;
use crate::multicomponent::network::reactive::{ReactingNetwork1d, ReactiveNetworkError};
use crate::multicomponent::network::wall_exchange::{
    SharedWall, SharedWallControl, SharedWallGas, SharedWallPort, SharedWallStorage,
};
use crate::multicomponent::reactive::ReactiveStepControl;
use crate::multicomponent::tests::{bank, source};
use crate::multicomponent::wall_exchange::MaterialWallCaloric;
use crate::multicomponent::{CaloricDomain, Conserved};
use fs_qty::Temperature;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionError, MaterialId, ThermalLaw,
    ThermalStorageBasis,
};
use fs_thermochem::kinetics::{
    ArrheniusRateV1, ElementaryReactionV1, IdealGasKineticsV1, ReactionProgressRateQuantityV1,
};
use fs_thermochem::{
    MolarEnergyQuantityV1, ReactionId, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};

struct Plateau;
impl ThermalLaw for Plateau {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(0.0..=20.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if h < 3.0 {
            2.0 + h
        } else if h <= 8.0 {
            5.0
        } else {
            5.0 + (h - 8.0) / 2.0
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}
struct UnitInternalStorage;
impl ThermalLaw for UnitInternalStorage {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::InternalEnergy
    }
    fn temperature_k(&self, u: f64) -> Result<f64, DepositionError> {
        if u.is_finite() && u > 0.0 {
            Ok(u)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}
fn specimen(e: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(71),
        segment_id: 2,
        body_id: 9,
        material_id: MaterialId(38),
        mass_kg: 1.0,
        specific_h_j_kg: e,
        geometry: CellGeometry {
            start_m: [0.0; 3],
            end_m: [1.0, 0.0, 0.0],
            volume_m3: 0.1,
            surface_area_m2: 1.0,
        },
        boundary: Boundary {
            ambient_k: 2.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    }
}
fn shared_storage(
    source: MaterialWallCaloric<'static>,
    ports: Vec<SharedWallPort>,
) -> SharedWall<'static> {
    SharedWall {
        storage: SharedWallStorage {
            thermal: source.initial_contact(0.0, WallDrag::None).unwrap().thermal,
            heat_storage: Some(source),
        },
        ports,
    }
}
fn wall_energy(wall: &NetworkWallOperator) -> f64 {
    let WallThermal::Caloric {
        mass_kg,
        specific_energy_j_kg,
        ..
    } = wall.shared_wall(0).unwrap().storage.thermal
    else {
        panic!("original energy coordinate");
    };
    mass_kg * specific_energy_j_kg
}
fn adiabatic(n: &MixtureNetwork1d) -> NetworkWallOperator<'static> {
    uniform_operator(
        n,
        WallContact {
            drag: WallDrag::None,
            thermal: WallThermal::Adiabatic,
        },
    )
}

#[test]
fn unequal_current_mixtures_share_one_original_latent_inventory_and_energy_gauges_cancel() {
    let densities = [[0.6, 0.4], [0.2, 0.8], [1.2, 0.3]];
    let temperatures = [2.0, 8.0, 3.0];
    let conductances = [0.3, 0.7, 0.4];
    let volumes = [0.2, 0.2, 0.1];
    let dt = 0.2;
    let mut gauge_results = Vec::new();
    for formation in [0.0, 1e6] {
        let original = gas(2, true, formation);
        // Fund caloric inversion below the tighter coupled heat allowance,
        // including equal-temperature endpoints of the passive root bracket.
        let mut domain = original.domain();
        domain.absolute_energy_tolerance = 1e-13;
        domain.relative_energy_tolerance = 1e-14;
        let g = MixtureGas::new(original.source().clone(), domain, || false).unwrap();
        let states: Vec<_> = densities
            .iter()
            .zip(temperatures)
            .map(|(d, t)| g.conserved(d, 0.0, (d[0] + 0.5 * d[1]) * t).unwrap())
            .collect();
        let branch = MixtureEuler1d::new(
            g.clone(),
            Length::new(0.5),
            0.4,
            states[..2].to_vec(),
            Boundaries::Periodic,
            || false,
        )
        .unwrap();
        let mut n = MixtureNetwork1d::new(
            vec![branch],
            vec![Plenum::new(0.1, states[2].clone()).unwrap()],
            vec![],
            || false,
        )
        .unwrap();
        let ports: Vec<_> = [
            SharedWallGas::Branch { branch: 0, cell: 0 },
            SharedWallGas::Branch { branch: 0, cell: 1 },
            SharedWallGas::Plenum(0),
        ]
        .into_iter()
        .zip(conductances)
        .map(|(gas, conductance_w_per_k)| SharedWallPort {
            gas,
            conductance_w_per_k,
        })
        .collect();
        let original = shared_storage(
            MaterialWallCaloric::from_thermal_law(Plateau, specimen(6.0), [0.0, 20.0]).unwrap(),
            ports,
        );
        let mut duplicate = original.clone();
        duplicate.ports.push(duplicate.ports[0]);
        assert!(matches!(
            adiabatic(&n).with_shared_walls(
                &n,
                vec![duplicate],
                SharedWallControl::default(),
                || false
            ),
            Err(WallExchangeError::Invalid(_))
        ));
        let mut wall = adiabatic(&n)
            .with_shared_walls(&n, vec![original], SharedWallControl::default(), || false)
            .unwrap();
        let before = fields(&n);
        let report = wall
            .advance_for(&mut n, Time::new(dt), controls(), || false)
            .unwrap_or_else(|error| panic!("formation datum {formation}: {error:?}"));
        let after = fields(&n);
        let mut heat = 0.0;
        for i in 0..3 {
            // Independent integration of the two original species polynomials:
            // rho*Cv(T) = 2.5*rho_A + 1.65*rho_B + (0.1*rho_A+0.025*rho_B)*T.
            let [a, b] = densities[i];
            let c0 = 2.5 * a + 1.65 * b;
            let c1 = 0.1 * a + 0.025 * b;
            let aa = 0.5 * volumes[i] * c1;
            let bb = volumes[i] * c0 + dt * conductances[i];
            let cc = -volumes[i] * (c0 * temperatures[i] + 0.5 * c1 * temperatures[i].powi(2))
                - dt * conductances[i] * 5.0;
            let expected_t = -2.0 * cc / (bb + (bb * bb - 4.0 * aa * cc).sqrt());
            let expected_q = dt * conductances[i] * (5.0 - expected_t);
            let (old, new) = if i < 2 {
                (&before.0[0][i], &after.0[0][i])
            } else {
                (&before.1[0], &after.1[0])
            };
            assert_eq!(new.partial_densities, old.partial_densities);
            assert_eq!(new.momentum, old.momentum);
            close(
                g.primitive(new.as_cell(), || false).unwrap().temperature,
                expected_t,
                2e-10,
            );
            close(volumes[i] * (new.energy - old.energy), expected_q, 2e-11);
            heat += expected_q;
        }
        let accepted = &wall.shared_wall(0).unwrap().storage;
        let source = accepted.heat_storage.as_ref().unwrap();
        assert_eq!(source.storage_basis(), ThermalStorageBasis::Enthalpy);
        assert_eq!(source.initial_cell().unwrap().material_id, MaterialId(38));
        assert_eq!(source.temperature_k(accepted.thermal).unwrap(), 5.0);
        close(wall_energy(&wall), 6.0 - heat, 2e-11);
        close(report.heat_to_gas_j, heat, 2e-11);
        close(report.wall_energy_change_j, -heat, 2e-11);
        close(report.energy_residual_j, 0.0, 2e-11);
        assert_eq!(report.external_heat_to_gas_j, 0.0);
        assert_eq!(n.time(), Time::new(0.0));
        gauge_results.push((after, wall_energy(&wall)));
    }
    assert_eq!(gauge_results[0], gauge_results[1]);
}

fn reacting(limit: usize, range: [f64; 2]) -> (ReactingNetwork1d, NetworkWallOperator<'static>) {
    let thermodynamics = bank(vec![
        source("A", R, 3.5, 0.0, 0.0),
        source("B", R, 5.0, 0.0, -10.0),
    ]);
    let reaction = ElementaryReactionV1::new(
        ReactionId::new("shared-wall-isomers").unwrap(),
        vec![(SpeciesId::new("A").unwrap(), 1)],
        vec![(SpeciesId::new("B").unwrap(), 1)],
        ArrheniusRateV1::new(
            Temperature::new(1.0),
            ReactionProgressRateQuantityV1::new(100.0),
            0.0,
            MolarEnergyQuantityV1::new(0.0),
        )
        .unwrap(),
    )
    .unwrap();
    let mechanism = IdealGasKineticsV1::new(
        thermodynamics.clone(),
        "shared-wall-isomers",
        1,
        vec![reaction],
        1e-12,
    )
    .unwrap();
    let g = MixtureGas::new(
        thermodynamics,
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    let branch = MixtureEuler1d::new(
        g.clone(),
        Length::new(0.5),
        0.4,
        vec![g.conserved(&[1.0, 0.0], 0.0, 2.0).unwrap(); 2],
        Boundaries::Periodic,
        || false,
    )
    .unwrap();
    let n = MixtureNetwork1d::new(vec![branch], vec![], vec![], || false).unwrap();
    let storage =
        MaterialWallCaloric::from_thermal_law(UnitInternalStorage, specimen(5.0), range).unwrap();
    let shared = shared_storage(
        storage,
        (0..2)
            .map(|cell| SharedWallPort {
                gas: SharedWallGas::Branch { branch: 0, cell },
                conductance_w_per_k: 100.0,
            })
            .collect(),
    );
    let wall = adiabatic(&n)
        .with_shared_walls(
            &n,
            vec![shared],
            SharedWallControl {
                maximum_source_evaluations: limit,
            },
            || false,
        )
        .unwrap();
    let chemistry = IsochoricKineticsSolve {
        minimum_temperature_k: 0.001,
        maximum_temperature_k: 50.0,
        energy_tolerance_j: 1e-10,
        initial_step_s: 1e-5,
        ..IsochoricKineticsSolve::default()
    };
    (
        ReactingNetwork1d::new(n, mechanism, chemistry, || false).unwrap(),
        wall,
    )
}
fn reactive_controls() -> ReactiveStepControl {
    ReactiveStepControl {
        maximum_coupling_retries: 0,
        ..ReactiveStepControl::default()
    }
}
fn reactive_states(flow: &ReactingNetwork1d) -> Vec<Conserved> {
    (0..2)
        .map(|i| owned(flow.branches()[0].cell(i).unwrap()))
        .collect()
}

#[test]
fn reaction_between_shared_wall_halves_refreshes_original_calorics_and_one_candidate_inventory() {
    let (mut flow, mut wall) = reacting(65_536, [0.01, 50.0]);
    let initial_gas = flow.totals().absolute_energy(flow.gas()).unwrap();
    let report = flow
        .step_with_wall(
            &mut wall,
            Time::new(0.001),
            reactive_controls(),
            controls(),
            || false,
        )
        .unwrap();
    let half = report.elapsed.value() * 0.5;
    let k = 100.0 * half;
    // Two identical gas volumes contact ONE Cwall=1 inventory. The first
    // source half's exact three-body equation gives each individual Q1.
    let q1 = k * (5.0 - 2.0) / (1.0 + k / 0.5 + 2.0 * k);
    let tw1 = 5.0 - 2.0 * q1;
    let state = flow.branches()[0].cell(0).unwrap();
    let [a, b] = [state.partial_densities[0], state.partial_densities[1]];
    assert!(
        b > 1e-3,
        "reaction must create a calorically distinct product"
    );
    close(a + b, 1.0, 1e-11);
    let cv = 2.5 * a + 4.0 * b;
    // Formation internal energy is -10*rho_B. Both chemical halves conserve
    // original absolute energy, and this uniform periodic branch has no flux.
    let before_second = (5.0 + q1 / 0.2 + 10.0 * b) / cv;
    let q2 = k * (tw1 - before_second) / (1.0 + k / (0.2 * cv) + 2.0 * k);
    let halves = report.wall.unwrap();
    close(halves.first.heat_to_gas_j, 2.0 * q1, 2e-10);
    close(halves.second.heat_to_gas_j, 2.0 * q2, 2e-9);
    close(wall_energy(&wall), tw1 - 2.0 * q2, 2e-9);
    close(
        flow.gas().primitive(state, || false).unwrap().temperature,
        before_second + q2 / (0.2 * cv),
        2e-9,
    );
    close(
        flow.totals().absolute_energy(flow.gas()).unwrap() - initial_gas + wall_energy(&wall) - 5.0,
        0.0,
        3e-9,
    );
    close(report.absolute_energy_residual_j, 0.0, 2e-10);
    assert_eq!(wall.last_report().unwrap().external_heat_to_gas_j, 0.0);
    assert_eq!(
        wall.last_shared_wall_report().unwrap().elapsed,
        report.elapsed
    );
    assert_eq!(
        wall.shared_wall(0)
            .unwrap()
            .storage
            .heat_storage
            .as_ref()
            .unwrap()
            .storage_basis(),
        ThermalStorageBasis::InternalEnergy
    );
}

#[test]
fn shared_wall_work_limits_late_source_refusal_and_cancellation_preserve_the_complete_trial() {
    let (mut reference, mut reference_wall) = reacting(65_536, [0.01, 50.0]);
    let mut polls = 0;
    let expected = reference
        .step_with_wall(
            &mut reference_wall,
            Time::new(0.001),
            reactive_controls(),
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let work = reference_wall
        .last_shared_wall_report()
        .unwrap()
        .source_evaluations;
    assert!(work > 1 && polls > 1);
    // The budget spans BOTH halves, and the exact observed work is sufficient.
    let (mut exact, mut exact_wall) = reacting(work, [0.01, 50.0]);
    assert_eq!(
        exact
            .step_with_wall(
                &mut exact_wall,
                Time::new(0.001),
                reactive_controls(),
                controls(),
                || false
            )
            .unwrap(),
        expected
    );
    assert_eq!(reactive_states(&exact), reactive_states(&reference));
    for failure in 0..3 {
        let range = if failure == 1 {
            [4.74, 6.0]
        } else {
            [0.01, 50.0]
        };
        let (mut flow, mut wall) = reacting(if failure == 0 { work - 1 } else { 65_536 }, range);
        let before = reactive_states(&flow);
        let inventory = wall.shared_wall(0).unwrap().clone();
        let mut seen = 0;
        let error = flow
            .step_with_wall(
                &mut wall,
                Time::new(0.001),
                reactive_controls(),
                controls(),
                || {
                    seen += 1;
                    failure == 2 && seen >= polls
                },
            )
            .unwrap_err();
        fn wall_refusal(error: &ReactiveNetworkError) -> bool {
            match error {
                ReactiveNetworkError::Wall(_) => true,
                ReactiveNetworkError::RetryLimit { last, .. } => wall_refusal(last),
                _ => false,
            }
        }
        if failure == 2 {
            assert_eq!(error, ReactiveNetworkError::Cancelled);
        } else {
            assert!(wall_refusal(&error), "{error:?}");
        }
        assert_eq!(reactive_states(&flow), before);
        assert_eq!(flow.time(), Time::new(0.0));
        assert!(flow.last_report().is_none());
        assert_eq!(wall.shared_wall(0).unwrap(), &inventory);
        assert!(wall.last_report().is_none());
        assert!(wall.last_shared_wall_report().is_none());
        if failure == 2 {
            assert_eq!(
                flow.step_with_wall(
                    &mut wall,
                    Time::new(0.001),
                    reactive_controls(),
                    controls(),
                    || false
                )
                .unwrap(),
                expected
            );
            assert_eq!(reactive_states(&flow), reactive_states(&reference));
            // Independent constructors retain different immutable source Arcs;
            // replay requires identical physical state and original ports.
            let actual = wall.shared_wall(0).unwrap();
            let expected_wall = reference_wall.shared_wall(0).unwrap();
            assert_eq!(actual.storage.thermal, expected_wall.storage.thermal);
            assert_eq!(actual.ports, expected_wall.ports);
            assert_eq!(
                wall.last_shared_wall_report(),
                reference_wall.last_shared_wall_report()
            );
        }
    }
}
