//! Original Fourier/film equations and one nonlinear inventory across gas ports.

use super::conduction::source;
use super::*;
use crate::multicomponent::network::wall_exchange::{
    MaterialWallConduction, SharedWall, SharedWallControl, SharedWallGas, SharedWallLink,
    SharedWallPort, SharedWallStorage,
};
use fs_material::thermal::conduction::ReservoirLayerGeometry;
use fs_thermal::enthalpy::DepositionError;

const DT: f64 = 0.2;
const WALL: f64 = 5.0;
const GAS: [f64; 3] = [3.0, 7.0, 4.0];
const FACE: [f64; 3] = [4.0, 6.0, 4.5];
const DENSITY: [f64; 3] = [1.0, 2.0, 1.5];
const VOLUME: [f64; 3] = [0.2, 0.2, 0.1];

fn energy(temperature: f64) -> f64 {
    // The original one-kilogram source has Cwall(T)=1+0.2T.
    temperature + 0.1 * temperature * temperature
}

fn manufactured(
    gas: &MixtureGas,
    geometry: ReservoirLayerGeometry,
    multiplier: f64,
    span: [f64; 2],
    maximum_source_evaluations: usize,
) -> (
    MixtureNetwork1d,
    NetworkWallOperator<'static>,
    [f64; 3],
    f64,
) {
    // Manufacture the ENDPOINT equations, not the numerical root iteration.
    // Fourier: P_i=C integral(Ts_i..Tw,T dT). Film: P_i=G_i(Ts_i-Tg_i).
    let power = FACE.map(|face| 0.5 * multiplier * (WALL * WALL - face * face));
    let heat = power.map(|power| DT * power);
    let film = std::array::from_fn::<_, 3, _>(|i| power[i] / (FACE[i] - GAS[i]));
    let initial = std::array::from_fn::<_, 3, _>(|i| {
        let capacity = 2.5 * DENSITY[i] * VOLUME[i];
        let temperature = GAS[i] - heat[i] / capacity;
        // The original gas has R_specific=1; conserved() takes pressure.
        gas.conserved(&[DENSITY[i]], 0.0, DENSITY[i] * temperature)
            .unwrap()
    });
    let branch = MixtureEuler1d::new(
        gas.clone(),
        Length::new(0.5),
        0.4,
        initial[..2].to_vec(),
        Boundaries::Periodic,
        || false,
    )
    .unwrap();
    let plenum = Plenum::new(VOLUME[2], initial[2].clone()).unwrap();
    let network = MixtureNetwork1d::new(vec![branch], vec![plenum], vec![], || false).unwrap();
    let initial_energy = energy(WALL) + heat.iter().sum::<f64>();
    let initial_wall = 2.0 * initial_energy / (1.0 + (1.0 + 0.4 * initial_energy).sqrt());
    let (transport, caloric) = source([0.0, 1.0], [1.0, 0.2], span);
    let layer = MaterialWallConduction::from_geometry(transport, geometry)
        .unwrap()
        .with_integrated_profile()
        .unwrap();
    let ports = [
        SharedWallGas::Branch { branch: 0, cell: 0 },
        SharedWallGas::Branch { branch: 0, cell: 1 },
        SharedWallGas::Plenum(0),
    ]
    .into_iter()
    .zip(film)
    .map(|(gas, conductance_w_per_k)| SharedWallPort {
        gas,
        conductance_w_per_k,
    })
    .collect();
    let storage = SharedWall {
        storage: SharedWallStorage {
            thermal: caloric
                .contact(initial_wall, 0.0, WallDrag::None)
                .unwrap()
                .thermal,
            heat_storage: Some(caloric),
        },
        ports,
    };
    let wall = uniform_operator(
        &network,
        WallContact {
            drag: WallDrag::None,
            thermal: WallThermal::Adiabatic,
        },
    )
    .with_shared_walls(
        &network,
        vec![storage],
        SharedWallControl {
            maximum_source_evaluations,
        },
        || false,
    )
    .unwrap()
    .with_shared_wall_conduction(vec![vec![Some(layer); 3]], || false)
    .unwrap();
    (network, wall, heat, initial_wall)
}

#[test]
fn shared_integrated_profiles_close_opposing_films_and_one_original_wall_inventory() {
    let g = gas(1, false, 0.0);
    // Independent resistance geometry: A/L and 2*pi*length/ln(r_outer/r_inner).
    // The cylinder exposes its OUTER face, whose area is not the multiplier.
    let geometries = [
        (ReservoirLayerGeometry::planar(0.1, 1.0).unwrap(), 0.1),
        (
            ReservoirLayerGeometry::cylindrical_face(0.2, 0.1, 0.01).unwrap(),
            std::f64::consts::TAU * 0.01 / 2.0_f64.ln(),
        ),
    ];
    for (geometry, multiplier, span) in geometries.into_iter().flat_map(|(geometry, multiplier)| {
        [[1.0, 10.0], [3.5, 10.0], [3.91, 10.0]].map(|span| (geometry, multiplier, span))
    }) {
        let (mut n, mut wall, heat, initial_wall) =
            manufactured(&g, geometry, multiplier, span, 65_536);
        let before = fields(&n);
        let initial_gas_energy = n.totals().absolute_energy(&g).unwrap();
        let original = wall.shared_wall(0).unwrap().clone();
        let layers = wall.shared_wall_conduction(0).unwrap().to_vec();
        let report = wall
            .advance_for(&mut n, Time::new(DT), controls(), || false)
            .unwrap();
        let accepted = wall.shared_wall(0).unwrap();
        let caloric = accepted.storage.heat_storage.as_ref().unwrap();
        let tw = caloric.temperature_k(accepted.storage.thermal).unwrap();
        close(tw, WALL, 3e-10);
        assert!(wall.shared_wall(1).is_none(), "exactly one physical wall");
        assert_eq!(accepted.ports, original.ports, "films must not compound");
        assert_eq!(accepted.storage.heat_storage, original.storage.heat_storage);
        assert_eq!(wall.shared_wall_conduction(0).unwrap(), layers);
        if span[0] > GAS[0] {
            // Only the actual solid profile must be supported. In particular,
            // unused outer root brackets must not query k at the colder gas.
            let transport = layers[0].as_ref().unwrap().source().unwrap();
            assert!(transport.at_temperature(GAS[0]).is_err());
            assert!(transport.at_temperature(FACE[0]).is_ok());
        }
        let after = fields(&n);
        let mut total_heat = 0.0;
        for i in 0..3 {
            let (old, new) = if i < 2 {
                (&before.0[0][i], &after.0[0][i])
            } else {
                (&before.1[0], &after.1[0])
            };
            let tg = g.primitive(new.as_cell(), || false).unwrap().temperature;
            let q = VOLUME[i] * (new.energy - old.energy);
            let film = accepted.ports[i].conductance_w_per_k;
            let face = tg + q / (DT * film);
            close(tg, GAS[i], 3e-10);
            close(q, heat[i], 3e-10);
            close(face, FACE[i], 3e-10);
            close(q, DT * 0.5 * multiplier * (tw * tw - face * face), 3e-10);
            close(q, DT * film * (FACE[i] - tg), 3e-10);
            assert_eq!(new.partial_densities, old.partial_densities);
            assert_eq!(new.momentum, old.momentum);
            // A wall-endpoint conductivity sample would conserve energy while
            // satisfying a DIFFERENT constitutive law; these checks reject it.
            let sampled = 1.0 / (1.0 / film + 1.0 / (multiplier * WALL));
            assert!((DT * sampled * (WALL - GAS[i]) - heat[i]).abs() > 1e-4);
            total_heat += q;
        }
        assert!(heat[0] > 0.0 && heat[1] < 0.0 && heat[2] > 0.0);
        let wall_change = energy(tw) - energy(initial_wall);
        close(wall_change, -heat.iter().sum::<f64>(), 3e-10);
        close(total_heat + wall_change, 0.0, 3e-10);
        close(report.heat_to_gas_j, total_heat, 3e-10);
        close(report.wall_energy_change_j, wall_change, 3e-10);
        close(
            n.totals().absolute_energy(&g).unwrap() - initial_gas_energy,
            total_heat,
            3e-10,
        );
        close(report.energy_residual_j, 0.0, 3e-10);
        close(report.maximum_thermal_residual_j, 0.0, 3e-10);
        assert_eq!(report.external_heat_to_gas_j, 0.0);
        let shared = wall.last_shared_wall_report().unwrap();
        close(shared.energy_residual_j, 0.0, 3e-10);
        close(shared.maximum_contact_residual_j, 0.0, 3e-10);
        assert_eq!(n.time(), Time::new(0.0));
    }
}

#[test]
fn shared_integrated_profile_refuses_an_unsupported_face_without_publishing() {
    let g = gas(1, false, 0.0);
    let (mut n, mut wall, _, initial_wall) = manufactured(
        &g,
        ReservoirLayerGeometry::planar(0.1, 1.0).unwrap(),
        0.1,
        [4.1, 10.0],
        65_536,
    );
    // Both wall states are supported; the cold port's manufactured Ts=4 is not.
    // Sampling only k(Twall) would silently accept this unsupported profile.
    assert!(initial_wall > 4.1 && WALL > 4.1 && FACE[0] < 4.1);
    let before = fields(&n);
    let original = wall.shared_wall(0).unwrap().clone();
    let layers = wall.shared_wall_conduction(0).unwrap().to_vec();
    assert!(matches!(
        wall.advance_for(&mut n, Time::new(DT), controls(), || false),
        Err(WallExchangeError::Thermal(DepositionError::MaterialDomain))
    ));
    assert_eq!(fields(&n), before);
    assert_eq!(wall.shared_wall(0), Some(&original));
    assert_eq!(wall.shared_wall_conduction(0).unwrap(), layers);
    assert!(wall.last_report().is_none());
    assert!(wall.last_shared_wall_report().is_none());
    assert_eq!(n.time(), Time::new(0.0));
}

#[test]
fn shared_integrated_profiles_meter_original_work_and_replay_after_late_cancellation() {
    let g = gas(1, false, 0.0);
    let build = |limit| {
        manufactured(
            &g,
            ReservoirLayerGeometry::planar(0.1, 1.0).unwrap(),
            0.1,
            [1.0, 10.0],
            limit,
        )
    };
    let (mut n, mut wall, _, _) = build(65_536);
    let first = wall
        .advance_for(&mut n, Time::new(DT), controls(), || false)
        .unwrap();
    let work = wall.last_shared_wall_report().unwrap().source_evaluations;
    assert!(
        work > 3,
        "original gas, storage and profile queries are metered"
    );
    let (mut exact, mut exact_wall, _, _) = build(work);
    assert_eq!(
        exact_wall
            .advance_for(&mut exact, Time::new(DT), controls(), || false)
            .unwrap(),
        first
    );
    assert_eq!(fields(&exact), fields(&n));
    assert_eq!(exact_wall.shared_wall(0), wall.shared_wall(0));
    let (mut denied, mut denied_wall, _, _) = build(work - 1);
    let before = fields(&denied);
    let original = denied_wall.shared_wall(0).unwrap().clone();
    assert_eq!(
        denied_wall
            .advance_for(&mut denied, Time::new(DT), controls(), || false)
            .unwrap_err(),
        WallExchangeError::Invalid("shared-wall nested source work budget")
    );
    assert_eq!(fields(&denied), before);
    assert_eq!(denied_wall.shared_wall(0), Some(&original));
    assert!(denied_wall.last_report().is_none());
    assert!(denied_wall.last_shared_wall_report().is_none());

    let (mut retry, mut retry_wall, _, _) = build(65_536);
    assert_eq!(
        retry_wall
            .advance_for(&mut retry, Time::new(DT), controls(), || false)
            .unwrap(),
        first
    );
    let accepted = fields(&retry);
    let storage = retry_wall.shared_wall(0).unwrap().clone();
    let prior_shared_report = retry_wall.last_shared_wall_report();
    let mut polls = 0;
    let next = wall
        .advance_for(&mut n, Time::new(DT), controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut seen = 0;
    assert_eq!(
        retry_wall
            .advance_for(&mut retry, Time::new(DT), controls(), || {
                seen += 1;
                seen >= polls
            })
            .unwrap_err(),
        WallExchangeError::Cancelled
    );
    assert_eq!(fields(&retry), accepted);
    assert_eq!(retry_wall.shared_wall(0), Some(&storage));
    assert_eq!(retry_wall.last_report(), Some(first));
    assert_eq!(retry_wall.last_shared_wall_report(), prior_shared_report);
    assert_eq!(retry.time(), Time::new(0.0));
    assert_eq!(
        retry_wall
            .advance_for(&mut retry, Time::new(DT), controls(), || false)
            .unwrap(),
        next
    );
    assert_eq!(fields(&retry), fields(&n));
    assert_eq!(retry_wall.shared_wall(0), wall.shared_wall(0));
    assert_eq!(
        retry_wall.last_shared_wall_report(),
        wall.last_shared_wall_report()
    );
}

#[test]
fn a_later_wall_link_must_admit_the_original_material_face_before_publication() {
    let g = gas(1, false, 0.0);
    let lower = 4.61;
    let (transport, hot) = source([0.0, 1.0], [1.0, 0.0], [lower, 10.0]);
    let (_, cold) = source([0.0, 1.0], [1.0, 0.0], [1.0, 10.0]);
    let layer = MaterialWallConduction::from_geometry(
        transport,
        ReservoirLayerGeometry::planar(0.1, 1.0).unwrap(),
    )
    .unwrap()
    .with_integrated_profile()
    .unwrap();
    let walls: Vec<_> = [(hot, 5.5, 0.01), (cold, 3.0, 0.0)]
        .into_iter()
        .enumerate()
        .map(|(cell, (caloric, temperature, film))| SharedWall {
            storage: SharedWallStorage {
                thermal: caloric
                    .contact(temperature, 0.0, WallDrag::None)
                    .unwrap()
                    .thermal,
                heat_storage: Some(caloric),
            },
            ports: vec![SharedWallPort {
                gas: SharedWallGas::Branch { branch: 0, cell },
                conductance_w_per_k: film,
            }],
        })
        .collect();
    let bind = |n: &MixtureNetwork1d| {
        uniform_operator(
            n,
            WallContact {
                drag: WallDrag::None,
                thermal: WallThermal::Adiabatic,
            },
        )
        .with_shared_walls(n, walls.clone(), SharedWallControl::default(), || false)
        .unwrap()
        .with_shared_wall_conduction(vec![vec![Some(layer.clone())], vec![None]], || false)
        .unwrap()
    };
    // The existing owner solves gas contacts, then the full-duration wall link.
    // CpA=CpB=1, Cgas=0.5, dt=1, film=.01, C=.1 give
    // q=.05*((5.5-q)^2-(2+102q)^2), independent of either scalar solver.
    let q = 2.0 * 1.3125 / (21.95 + (21.95_f64.powi(2) + 4.0 * 520.15 * 1.3125).sqrt());
    let gas_endpoint = 2.0 + 2.0 * q;
    let wall_after_gas = 5.5 - q;
    let face_after_gas = 2.0 + 102.0 * q;
    close(
        q,
        0.05 * (wall_after_gas.powi(2) - face_after_gas.powi(2)),
        3e-13,
    );
    // A G=1, dt=1 link between two unit capacities transfers (Ta-Tb)/3.
    let final_wall = wall_after_gas - (wall_after_gas - 3.0) / 3.0;
    let final_face = (0.01 + final_wall * final_wall + 0.2 * gas_endpoint).sqrt() - 0.1;
    assert!(face_after_gas > lower);
    assert!(final_wall > lower && final_face < lower);

    let mut gas_only = network(&g, &[0.0], false);
    let mut unlinked = bind(&gas_only);
    unlinked
        .advance_for(&mut gas_only, Time::new(1.0), controls(), || false)
        .unwrap();
    let accepted = &unlinked.shared_wall(0).unwrap().storage;
    close(
        accepted
            .heat_storage
            .as_ref()
            .unwrap()
            .temperature_k(accepted.thermal)
            .unwrap(),
        wall_after_gas,
        3e-10,
    );
    close(
        g.primitive(gas_only.branches()[0].cell(0).unwrap(), || false)
            .unwrap()
            .temperature,
        gas_endpoint,
        3e-10,
    );

    let mut n = network(&g, &[0.0], false);
    let mut linked = bind(&n)
        .with_shared_wall_links(
            vec![SharedWallLink {
                walls: [0, 1],
                conductance_w_per_k: 1.0,
            }],
            || false,
        )
        .unwrap();
    let before = fields(&n);
    let accepted: Vec<_> = (0..2)
        .map(|i| linked.shared_wall(i).unwrap().clone())
        .collect();
    assert!(matches!(
        linked.advance_for(&mut n, Time::new(1.0), controls(), || false),
        Err(WallExchangeError::Thermal(DepositionError::MaterialDomain))
    ));
    assert_eq!(fields(&n), before);
    for (i, original) in accepted.iter().enumerate() {
        assert_eq!(linked.shared_wall(i), Some(original));
    }
    assert!(linked.last_report().is_none());
    assert!(linked.last_shared_wall_report().is_none());
    assert!(linked.shared_wall_link_reports().is_none());
    assert_eq!(n.time(), Time::new(0.0));
}
