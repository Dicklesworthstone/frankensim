//! Original two-sided conduction, nonlinear latent transfer and stage ownership.

use super::super::conduction::source as temperature_source;
use super::super::phase_conduction::source as phase_source;
use super::material_conduction::{
    frozen_at, frozen_walls, mixture_at, phase_layer, phase_wall, port,
};
use super::*;
use crate::multicomponent::network::wall_exchange::{
    MaterialWallConduction, SharedWallLink, WallConductionEvaluation,
};
use fs_qty::Area;

fn link(a: usize, b: usize, interface: f64) -> SharedWallLink {
    SharedWallLink {
        walls: [a, b],
        conductance_w_per_k: interface,
    }
}
fn material_wall(
    t: f64,
    capacity: f64,
    k: [f64; 2],
    geometry: [f64; 2],
    cell: usize,
) -> (SharedWall<'static>, MaterialWallConduction) {
    let (transport, caloric) = temperature_source(k, [capacity, 0.0], [1.0, 10.0]);
    let wall = SharedWall {
        storage: SharedWallStorage {
            thermal: caloric.contact(t, 0.0, WallDrag::None).unwrap().thermal,
            heat_storage: Some(caloric),
        },
        ports: vec![port(cell, 0.0)],
    };
    let layer =
        MaterialWallConduction::new(transport, Length::new(geometry[0]), Area::new(geometry[1]))
            .unwrap();
    (wall, layer)
}
fn temperature(wall: &SharedWall<'_>) -> f64 {
    let storage = &wall.storage;
    storage
        .heat_storage
        .as_ref()
        .unwrap()
        .temperature_k(storage.thermal)
        .unwrap()
}
fn enthalpy(wall: &SharedWall<'_>) -> f64 {
    let WallThermal::Caloric {
        mass_kg,
        specific_energy_j_kg,
        ..
    } = wall.storage.thermal
    else {
        panic!("original phase coordinate")
    };
    assert_eq!(mass_kg, 1.0);
    specific_energy_j_kg
}
fn bound(
    n: &MixtureNetwork1d,
    walls: Vec<SharedWall<'static>>,
    links: Vec<SharedWallLink>,
    layers: Vec<[Option<MaterialWallConduction>; 2]>,
    limit: usize,
) -> NetworkWallOperator<'static> {
    adiabatic(n)
        .with_shared_walls(
            n,
            walls,
            SharedWallControl {
                maximum_source_evaluations: limit,
            },
            || false,
        )
        .unwrap()
        .with_shared_wall_links(links, || false)
        .unwrap()
        .with_shared_wall_link_conduction(layers, || false)
        .unwrap()
}

#[test]
fn two_original_material_states_and_geometry_match_oriented_calorimetry_in_both_facades() {
    let g = gas(1, false, 0.0);
    for scale in [1.0, 3.0] {
        for reverse in [false, true] {
            let (a, ka) = material_wall(8.0, 0.7, [0.2, 0.1], [0.2, 0.3], 0);
            let (b, kb) = material_wall(2.0, 1.1, [scale, 0.2 * scale], [0.1, 0.4], 1);
            let orientation = if reverse { [1, 0] } else { [0, 1] };
            let layers = if reverse {
                [Some(kb), Some(ka)]
            } else {
                [Some(ka), Some(kb)]
            };
            let port = SharedWallLink {
                walls: orientation,
                conductance_w_per_k: 2.5,
            };
            let mut mixture = mixture_at(&g, 5.0);
            let mut frozen = frozen_at(&g, 5.0);
            let mut mw = bound(
                &mixture,
                vec![a.clone(), b.clone()],
                vec![port],
                vec![layers.clone()],
                65_536,
            );
            let mut fw = frozen_walls(&frozen, vec![a, b])
                .with_shared_wall_links(vec![port], || false)
                .unwrap()
                .with_shared_wall_link_conduction(vec![layers.clone()], || false)
                .unwrap();
            let original_gas = fields(&mixture);
            let (mut ta, mut tb) = (8.0_f64, 2.0_f64);
            for _ in 0..3 {
                let dt = 0.01;
                // Each material is queried at its OWN entering temperature.
                let effective = 1.0
                    / (1.0 / 2.5
                        + 0.2 / (0.3 * (0.2 + 0.1 * ta))
                        + 0.1 / (0.4 * scale * (1.0 + 0.2 * tb)));
                let q =
                    dt * effective * (ta - tb) / (1.0 + dt * effective * (1.0 / 0.7 + 1.0 / 1.1));
                let mr = mw
                    .advance_for(&mut mixture, Time::new(dt), controls(), || false)
                    .unwrap();
                let fr = frozen
                    .step_with_walls(
                        Time::new(dt),
                        StepControl::default(),
                        &mut fw,
                        controls(),
                        || false,
                    )
                    .unwrap();
                assert_eq!(fr.elapsed, Time::new(dt));
                ta -= q / 0.7;
                tb += q / 1.1;
                for (i, expected) in [ta, tb].into_iter().enumerate() {
                    close(temperature(mw.shared_wall(i).unwrap()), expected, 2e-11);
                    close(temperature(fw.shared_wall(i).unwrap()), expected, 2e-11);
                }
                let signed = if reverse { -q } else { q };
                for receipt in [
                    mw.shared_wall_link_reports().unwrap()[0],
                    fw.shared_wall_link_reports().unwrap()[0],
                ] {
                    close(receipt.heat_a_to_b_j, signed, 2e-11);
                    close(receipt.wall_energy_change_j[0], -signed, 2e-11);
                    close(receipt.wall_energy_change_j[1], signed, 2e-11);
                    close(receipt.flux_residual_j, 0.0, 2e-11);
                    close(receipt.energy_residual_j, 0.0, 2e-11);
                }
                close(0.7 * ta + 1.1 * tb, 0.7 * 8.0 + 1.1 * 2.0, 2e-11);
                assert_eq!(fields(&mixture), original_gas);
                assert_eq!(mr.heat_to_gas_j, 0.0);
                assert_eq!(mr.external_heat_to_gas_j, 0.0);
                close(
                    fr.thermal.unwrap().gas_wall_conservation_residual_j,
                    0.0,
                    2e-11,
                );
                assert_eq!(mw.shared_wall_links().unwrap(), [port]);
                assert_eq!(mw.shared_wall_link_conduction(0), Some(&layers));
                assert_eq!(fw.shared_wall_link_conduction(0), Some(&layers));
            }
        }
    }
}

#[test]
fn growing_phase_conductivity_transfers_more_than_the_initial_flux_bound() {
    let g = gas(1, false, 0.0);
    let source = phase_source(1.0, 1.0, 4.0);
    for reverse in [false, true] {
        let (hot, hot_layer) = material_wall(8.0, 1.0, [0.2, 0.1], [0.1, 0.25], 0);
        let cold = phase_wall(&source, 4.0, vec![port(1, 0.0)]);
        let cold_layer = phase_layer(&source, 0.1, WallConductionEvaluation::ImplicitEndpoint);
        let (connection, layers) = if reverse {
            (link(1, 0, 0.7), [Some(cold_layer), Some(hot_layer)])
        } else {
            (link(0, 1, 0.7), [Some(hot_layer), Some(cold_layer)])
        };
        let mut n = mixture_at(&g, 5.0);
        let mut wall = bound(&n, vec![hot, cold], vec![connection], vec![layers], 65_536);
        let (mut ta, mut h) = (8.0_f64, 4.0_f64);
        for interval in 0..2 {
            let dt = 1.0;
            // The hot side's original k(T) is lagged; cold-side phase k grows on
            // its T=5 plateau. Gsolid_cold=a+bQ and Thot_new=Thot_old-Q/C, C=1.
            // Eliminating the series resistances gives A Q²+B Q-C=0.
            let interface_and_a = 1.0 / (1.0 / 0.7 + 0.1 / (0.25 * (0.2 + 0.1 * ta)));
            let a = 2.5 * (0.1 + 0.19 * (h - 4.0));
            let b = 2.5 * 0.19;
            let difference = ta - 5.0;
            let aa = b * (1.0 + dt * interface_and_a);
            let bb = interface_and_a + a - dt * interface_and_a * (b * difference - a);
            let cc = dt * interface_and_a * a * difference;
            let q = 2.0 * cc / (bb + (bb * bb + 4.0 * aa * cc).sqrt());
            let old_flux_bound = dt * difference / (1.0 / interface_and_a + 1.0 / a);
            if interval == 0 {
                assert!(q > 1.2 * old_flux_bound);
            }
            let report = wall
                .advance_for(&mut n, Time::new(dt), controls(), || false)
                .unwrap();
            ta -= q;
            h += q;
            close(temperature(wall.shared_wall(0).unwrap()), ta, 2e-11);
            assert_eq!(temperature(wall.shared_wall(1).unwrap()), 5.0);
            close(enthalpy(wall.shared_wall(1).unwrap()), h, 2e-11);
            close(ta + h, 12.0, 2e-11);
            let receipt = wall.shared_wall_link_reports().unwrap()[0];
            close(receipt.heat_a_to_b_j, if reverse { -q } else { q }, 2e-11);
            close(receipt.energy_residual_j, 0.0, 2e-11);
            close(receipt.flux_residual_j, 0.0, 2e-11);
            assert_eq!(report.heat_to_gas_j, 0.0);
            assert_eq!(report.external_heat_to_gas_j, 0.0);
            assert!(std::sync::Arc::ptr_eq(
                wall.shared_wall(1)
                    .unwrap()
                    .storage
                    .heat_storage
                    .as_ref()
                    .unwrap()
                    .phase_curve()
                    .unwrap(),
                source.curve()
            ));
        }
    }
    // A source from a different original phase chart cannot bind this wall.
    let other = phase_source(2.0, 1.0, 4.0);
    let cold = phase_wall(&source, 4.0, vec![port(1, 0.0)]);
    let (hot, _) = material_wall(8.0, 1.0, [1.0, 0.0], [0.1, 0.25], 0);
    let n = mixture_at(&g, 5.0);
    assert!(
        adiabatic(&n)
            .with_shared_walls(&n, vec![hot, cold], SharedWallControl::default(), || false)
            .unwrap()
            .with_shared_wall_links(vec![link(0, 1, 0.7)], || false)
            .unwrap()
            .with_shared_wall_link_conduction(
                vec![[
                    None,
                    Some(phase_layer(
                        &other,
                        0.1,
                        WallConductionEvaluation::ImplicitEndpoint
                    ))
                ]],
                || false
            )
            .is_err()
    );
}

#[test]
fn lagged_links_refresh_after_preceding_links_and_share_work_and_rollback() {
    let g = gas(1, false, 0.0);
    let (a, _) = material_wall(8.0, 1.0, [1.0, 0.0], [0.1, 0.25], 0);
    let (b, layer) = material_wall(2.0, 1.0, [0.2, 0.1], [0.1, 0.25], 1);
    let build = |limit, material, coefficient| {
        let n = mixture_at(&g, 5.0);
        let wall = adiabatic(&n)
            .with_shared_walls(
                &n,
                vec![a.clone(), b.clone()],
                SharedWallControl {
                    maximum_source_evaluations: limit,
                },
                || false,
            )
            .unwrap()
            .with_shared_wall_links(vec![link(0, 1, 0.7), link(1, 0, coefficient)], || false)
            .unwrap();
        let wall = if material {
            wall.with_shared_wall_link_conduction(
                vec![[None, None], [Some(layer.clone()), None]],
                || false,
            )
            .unwrap()
        } else {
            wall
        };
        (n, wall)
    };
    let (mut n, mut wall) = build(65_536, true, 0.4);
    let first = wall
        .advance_for(&mut n, Time::new(0.2), controls(), || false)
        .unwrap();
    let q1 = 0.2 * 0.7 * 6.0 / (1.0 + 2.0 * 0.2 * 0.7);
    let middle = [8.0 - q1, 2.0 + q1];
    let effective = 1.0 / (1.0 / 0.4 + 0.1 / (0.25 * (0.2 + 0.1 * middle[1])));
    let q2 = 0.2 * effective * (middle[1] - middle[0]) / (1.0 + 2.0 * 0.2 * effective);
    let receipts = wall.shared_wall_link_reports().unwrap();
    assert_eq!(
        receipts[1].initial_state,
        [receipts[0].endpoint_state[1], receipts[0].endpoint_state[0]]
    );
    close(receipts[0].heat_a_to_b_j, q1, 2e-11);
    close(receipts[1].heat_a_to_b_j, q2, 2e-11);
    close(
        temperature(wall.shared_wall(0).unwrap()),
        middle[0] + q2,
        2e-11,
    );
    close(
        temperature(wall.shared_wall(1).unwrap()),
        middle[1] - q2,
        2e-11,
    );
    // Compare identical physical coefficients to isolate counted k queries.
    let solid = layer
        .source()
        .unwrap()
        .at_temperature(receipts[1].initial_temperature_k[0])
        .unwrap()
        .conductivity_w_m_k()
        / (layer.length().value() / layer.area().value());
    let same_coefficient = 0.4_f64.min(solid) / (1.0 + 0.4_f64.min(solid) / 0.4_f64.max(solid));
    let work = wall.last_shared_wall_report().unwrap().source_evaluations;
    let (mut constant, mut cw) = build(65_536, false, same_coefficient);
    cw.advance_for(&mut constant, Time::new(0.2), controls(), || false)
        .unwrap();
    assert_eq!(
        work,
        cw.last_shared_wall_report().unwrap().source_evaluations + 3,
        "initial k, this stage's endpoint, and the complete sequence endpoint"
    );
    assert_eq!(fields(&n), fields(&constant));
    assert_eq!(wall.shared_wall(0), cw.shared_wall(0));
    assert_eq!(wall.shared_wall(1), cw.shared_wall(1));
    let (mut exact, mut ew) = build(work, true, 0.4);
    assert_eq!(
        ew.advance_for(&mut exact, Time::new(0.2), controls(), || false)
            .unwrap(),
        first
    );
    let (mut denied, mut dw) = build(work - 1, true, 0.4);
    let before = fields(&denied);
    assert!(
        dw.advance_for(&mut denied, Time::new(0.2), controls(), || false)
            .is_err()
    );
    assert_eq!(fields(&denied), before);
    assert_eq!(dw.shared_wall(0), Some(&a));
    assert_eq!(dw.shared_wall(1), Some(&b));
    assert!(dw.shared_wall_link_reports().is_none());
    let (mut retry, mut rw) = build(65_536, true, 0.4);
    rw.advance_for(&mut retry, Time::new(0.2), controls(), || false)
        .unwrap();
    let saved_fields = fields(&retry);
    let saved_walls = [
        rw.shared_wall(0).unwrap().clone(),
        rw.shared_wall(1).unwrap().clone(),
    ];
    let saved_receipts = rw.shared_wall_link_reports().unwrap().to_vec();
    let saved_report = rw.last_shared_wall_report();
    let mut polls = 0;
    let next = wall
        .advance_for(&mut n, Time::new(0.2), controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut seen = 0;
    assert_eq!(
        rw.advance_for(&mut retry, Time::new(0.2), controls(), || {
            seen += 1;
            seen >= polls
        })
        .unwrap_err(),
        WallExchangeError::Cancelled
    );
    assert_eq!(fields(&retry), saved_fields);
    for (i, old) in saved_walls.iter().enumerate() {
        assert_eq!(rw.shared_wall(i), Some(old));
    }
    assert_eq!(rw.shared_wall_link_reports().unwrap(), saved_receipts);
    assert_eq!(rw.last_shared_wall_report(), saved_report);
    assert_eq!(rw.last_report(), Some(first));
    assert_eq!(retry.time(), Time::new(0.0));
    assert_eq!(
        rw.advance_for(&mut retry, Time::new(0.2), controls(), || false)
            .unwrap(),
        next
    );
    assert_eq!(
        rw.shared_wall_link_reports(),
        wall.shared_wall_link_reports()
    );
}

#[test]
fn later_links_cannot_leave_an_earlier_pinned_material_outside_its_phase_domain() {
    let g = gas(1, false, 0.0);
    for evaluation in [
        WallConductionEvaluation::InitialState,
        WallConductionEvaluation::ImplicitEndpoint,
    ] {
        let source = phase_source(1.0, 0.001, 4.0);
        let cold = phase_wall(&source, 4.0, vec![port(0, 0.0)]);
        let (rest, _) = material_wall(5.0, 1.0, [1.0, 0.0], [0.1, 0.25], 1);
        let (hot, _) = material_wall(8.0, 1.0, [1.0, 0.0], [0.1, 0.25], 2);
        let storage = vec![cold, rest, hot];
        let layer = phase_layer(&source, 0.1, evaluation);
        let n3 = || {
            let branch = MixtureEuler1d::new(
                g.clone(),
                Length::new(0.5),
                0.4,
                vec![g.conserved(&[1.0], 0.0, 5.0).unwrap(); 3],
                Boundaries::Periodic,
                || false,
            )
            .unwrap();
            MixtureNetwork1d::new(vec![branch], vec![], vec![], || false).unwrap()
        };
        let mut supported = n3();
        let mut first = bound(
            &supported,
            storage.clone(),
            vec![link(0, 1, 0.5)],
            vec![[Some(layer.clone()), None]],
            65_536,
        );
        first
            .advance_for(&mut supported, Time::new(0.01), controls(), || false)
            .unwrap();
        assert_eq!(enthalpy(first.shared_wall(0).unwrap()), 4.0);
        let mut n = n3();
        let before = fields(&n);
        let links = vec![link(0, 1, 0.5), link(2, 0, 1.0)];
        let layers = vec![[Some(layer), None], [None, None]];
        let mut walls = bound(&n, storage.clone(), links.clone(), layers.clone(), 65_536);
        // Link0 is at rest and fully supported. Link1 then melts beyond the
        // narrow claim; its broad alternative must not replace the pinned k.
        assert!(matches!(
            walls.advance_for(&mut n, Time::new(0.01), controls(), || false),
            Err(WallExchangeError::PhaseConduction(_))
        ));
        assert_eq!(fields(&n), before);
        assert_eq!(n.time(), Time::new(0.0));
        for (i, old) in storage.iter().enumerate() {
            assert_eq!(walls.shared_wall(i), Some(old));
        }
        assert!(walls.last_report().is_none());
        assert!(walls.last_shared_wall_report().is_none());
        assert!(walls.shared_wall_link_reports().is_none());
        let fg = frozen_at(&g, 5.0).gas().clone();
        let state = fg.conserved(1.0, 0.0, 5.0).unwrap();
        let branch = crate::dynamics::Euler1d::new(
            fg,
            Length::new(0.5),
            0.4,
            vec![state; 3],
            Boundaries::Periodic,
            || false,
        )
        .unwrap();
        let mut frozen =
            crate::dynamics::network::GasNetwork1d::new(vec![branch], vec![], vec![], || false)
                .unwrap();
        let before = frozen.branches()[0].states().to_vec();
        let mut fw = frozen_walls(&frozen, storage.clone())
            .with_shared_wall_links(links, || false)
            .unwrap()
            .with_shared_wall_link_conduction(layers, || false)
            .unwrap();
        let error = frozen
            .step_with_walls(
                Time::new(0.01),
                StepControl {
                    maximum_halvings: 0,
                    ..StepControl::default()
                },
                &mut fw,
                controls(),
                || false,
            )
            .unwrap_err();
        let crate::dynamics::Error::StepRetries { last, .. } = error else {
            panic!("{error:?}")
        };
        assert!(matches!(*last, crate::dynamics::Error::Thermal(ref e)
            if matches!(**e, WallExchangeError::PhaseConduction(_))));
        assert_eq!(frozen.branches()[0].states(), before);
        assert_eq!(frozen.time(), Time::new(0.0));
        assert!(frozen.last_report().is_none());
        for (i, old) in storage.iter().enumerate() {
            assert_eq!(fw.shared_wall(i), Some(old));
        }
        assert!(fw.last_report().is_none());
        assert!(fw.last_shared_wall_report().is_none());
        assert!(fw.shared_wall_link_reports().is_none());
    }
}
