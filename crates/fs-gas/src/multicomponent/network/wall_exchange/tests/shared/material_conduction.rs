//! Original k(T,phi), one finite inventory, and atomic material-domain refusal.

use super::super::conduction::source as temperature_source;
use super::super::phase_conduction::{cell as phase_cell, source as phase_source};
use super::*;
use crate::dynamics::network::{
    BranchWallContact, GasNetwork1d, NetworkWallOperator as FrozenWalls,
};
use crate::dynamics::{Euler1d, FrozenGas};
use crate::multicomponent::network::wall_exchange::{
    MaterialWallConduction, PhaseWallConductionSource, SharedWallLink, WallConductionEvaluation,
};
use fs_qty::Area;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{Composition, CompositionBasis};
use std::sync::Arc;

fn port(cell: usize, film: f64) -> SharedWallPort {
    SharedWallPort {
        gas: SharedWallGas::Branch { branch: 0, cell },
        conductance_w_per_k: film,
    }
}
fn mixture_at(g: &MixtureGas, temperature: f64) -> MixtureNetwork1d {
    let branch = MixtureEuler1d::new(
        g.clone(),
        Length::new(0.5),
        0.4,
        vec![g.conserved(&[1.0], 0.0, temperature).unwrap(); 2],
        Boundaries::Periodic,
        || false,
    )
    .unwrap();
    MixtureNetwork1d::new(vec![branch], vec![], vec![], || false).unwrap()
}
fn frozen_at(g: &MixtureGas, temperature: f64) -> GasNetwork1d {
    let source = FrozenIdealGasMixtureModelV1::new(
        g.source().models().to_vec(),
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
    )
    .unwrap();
    let gas = FrozenGas::new(source, g.domain(), || false).unwrap();
    let state = gas.conserved(1.0, 0.0, temperature).unwrap();
    let branch = Euler1d::new(
        gas,
        Length::new(0.5),
        0.4,
        vec![state; 2],
        Boundaries::Periodic,
        || false,
    )
    .unwrap();
    GasNetwork1d::new(vec![branch], vec![], vec![], || false).unwrap()
}
fn frozen_walls(n: &GasNetwork1d, walls: Vec<SharedWall<'static>>) -> FrozenWalls<'static> {
    FrozenWalls::new(n, vec![], || false)
        .unwrap()
        .with_branch_contacts(
            n,
            vec![Some(vec![
                BranchWallContact {
                    contact: WallContact {
                        drag: WallDrag::None,
                        thermal: WallThermal::Adiabatic
                    },
                    heat_storage: None,
                };
                2
            ])],
            || false,
        )
        .unwrap()
        .with_shared_walls(n, walls, SharedWallControl::default(), || false)
        .unwrap()
}
fn phase_wall(
    source: &Arc<PhaseWallConductionSource>,
    h: f64,
    ports: Vec<SharedWallPort>,
) -> SharedWall<'static> {
    shared_storage(
        MaterialWallCaloric::from_phase_curve(source.curve().clone(), phase_cell(h), [0.0, 19.0])
            .unwrap(),
        ports,
    )
}
fn phase_layer(
    source: &Arc<PhaseWallConductionSource>,
    length: f64,
    evaluation: WallConductionEvaluation,
) -> MaterialWallConduction {
    MaterialWallConduction::from_phase_source(source.clone(), Length::new(length), Area::new(0.25))
        .unwrap()
        .with_evaluation(evaluation)
}

#[test]
fn implicit_shared_phase_conductivity_matches_independent_roots_and_source_geometry_swaps() {
    let g = gas(1, false, 0.0);
    for (scale, length) in [(1.0, 0.1), (3.0, 0.1), (1.0, 0.2)] {
        for (h0, tg0) in [(10.0, 2.0), (4.0, 8.0)] {
            let source = phase_source(scale, 1.0, h0);
            let layer = phase_layer(&source, length, WallConductionEvaluation::ImplicitEndpoint);
            let mut n = mixture_at(&g, tg0);
            let storage = phase_wall(&source, h0, vec![port(0, 0.3), port(1, 0.3)]);
            let mut wall = adiabatic(&n)
                .with_shared_walls(&n, vec![storage], SharedWallControl::default(), || false)
                .unwrap()
                .with_shared_wall_conduction(vec![vec![Some(layer.clone()); 2]], || false)
                .unwrap();
            let initial_energy = n.totals().absolute_energy(&g).unwrap() + h0;
            let (mut h, mut tg) = (h0, tg0);
            for _ in 0..2 {
                // Tw=5, two identical Cgas=1/2 volumes, one kilogram of wall.
                // h_new=h_old-2Q, so Gsolid=a-bQ. Eliminating both series
                // resistances gives -A Q²+B Q-C=0, independent of the solver.
                let dt = 0.2;
                let a = 0.25 / length * scale * (0.1 + 0.19 * (h - 4.0));
                let b = 0.25 / length * scale * 0.19 * 2.0;
                let difference = 5.0 - tg;
                let qa = b * (1.0 + dt * 0.3 / 0.5);
                let qb = 0.3 + a + dt * 0.3 * (a / 0.5 + b * difference);
                let qc = dt * 0.3 * a * difference;
                let q = 2.0 * qc / (qb + (qb * qb - 4.0 * qa * qc).sqrt());
                let before = fields(&n);
                let report = wall
                    .advance_for(&mut n, Time::new(dt), controls(), || false)
                    .unwrap();
                h -= 2.0 * q;
                tg += q / 0.5;
                for (i, old) in before.0[0].iter().enumerate() {
                    let new = n.branches()[0].cell(i).unwrap();
                    assert_eq!(new.partial_densities, old.partial_densities);
                    assert_eq!(new.momentum, old.momentum);
                    close(g.primitive(new, || false).unwrap().temperature, tg, 2e-10);
                    close(0.2 * (new.energy - old.energy), q, 2e-11);
                }
                close(wall_energy(&wall), h, 2e-11);
                close(report.heat_to_gas_j, 2.0 * q, 2e-11);
                close(report.wall_energy_change_j, -2.0 * q, 2e-11);
                close(report.energy_residual_j, 0.0, 2e-11);
                assert_eq!(report.external_heat_to_gas_j, 0.0);
                let accepted = wall.shared_wall(0).unwrap();
                let caloric = accepted.storage.heat_storage.as_ref().unwrap();
                assert!(Arc::ptr_eq(caloric.phase_curve().unwrap(), source.curve()));
                assert_eq!(
                    caloric.temperature_k(accepted.storage.thermal).unwrap(),
                    5.0
                );
                assert_eq!(accepted.ports, vec![port(0, 0.3), port(1, 0.3)]);
                assert_eq!(
                    wall.shared_wall_conduction(0).unwrap()[0],
                    Some(layer.clone())
                );
                close(
                    n.totals().absolute_energy(&g).unwrap() + wall_energy(&wall),
                    initial_energy,
                    3e-11,
                );
            }
            // A conductivity source from a different original chart cannot be
            // substituted for the caloric owner, even with the same geometry.
            let alien = phase_source(scale + 1.0, 1.0, h0);
            let n = mixture_at(&g, tg0);
            assert!(
                adiabatic(&n)
                    .with_shared_walls(
                        &n,
                        vec![phase_wall(&source, h0, vec![port(0, 0.3)])],
                        SharedWallControl::default(),
                        || false,
                    )
                    .unwrap()
                    .with_shared_wall_conduction(
                        vec![vec![Some(phase_layer(
                            &alien,
                            length,
                            WallConductionEvaluation::ImplicitEndpoint,
                        ))]],
                        || false
                    )
                    .is_err()
            );
        }
    }
}

#[test]
fn lagged_original_conductivity_and_capacity_refresh_with_frozen_mixture_parity() {
    let g = gas(1, false, 0.0);
    let mut mixture = mixture_at(&g, 2.0);
    let mut frozen = frozen_at(&g, 2.0);
    let (transport, caloric) = temperature_source([0.2, 0.1], [1.0, 0.2], [1.0, 10.0]);
    let layer = MaterialWallConduction::new(transport, Length::new(0.1), Area::new(0.25)).unwrap();
    let storage = SharedWall {
        storage: SharedWallStorage {
            thermal: caloric.contact(5.0, 0.0, WallDrag::None).unwrap().thermal,
            heat_storage: Some(caloric),
        },
        ports: vec![port(0, 0.3), port(1, 0.3)],
    };
    let bindings = vec![vec![Some(layer.clone()); 2]];
    let mut fw = frozen_walls(&frozen, vec![storage.clone()])
        .with_shared_wall_conduction(bindings.clone(), || false)
        .unwrap();
    let mut mw = adiabatic(&mixture)
        .with_shared_walls(
            &mixture,
            vec![storage],
            SharedWallControl::default(),
            || false,
        )
        .unwrap()
        .with_shared_wall_conduction(bindings, || false)
        .unwrap();
    let (mut tw, mut tg) = (5.0_f64, 2.0_f64);
    for _ in 0..3 {
        let dt = 0.01;
        let k = 0.2 + 0.1 * tw;
        let effective = 1.0 / (1.0 / 0.3 + 0.1 / (0.25 * k));
        let response = dt * effective / (1.0 + dt * effective / 0.5);
        // Original integrated Cv(T)=1+0.2T, one inventory and two gas ports.
        let linear = 1.0 + 2.0 * response;
        let rhs = tw + 0.1 * tw * tw + 2.0 * response * tg;
        let next = 2.0 * rhs / (linear + (linear * linear + 0.4 * rhs).sqrt());
        let q = response * (next - tg);
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
        tw = next;
        tg += q / 0.5;
        let ms = &mw.shared_wall(0).unwrap().storage;
        let fs = &fw.shared_wall(0).unwrap().storage;
        close(
            ms.heat_storage
                .as_ref()
                .unwrap()
                .temperature_k(ms.thermal)
                .unwrap(),
            tw,
            2e-10,
        );
        close(
            fs.heat_storage
                .as_ref()
                .unwrap()
                .temperature_k(fs.thermal)
                .unwrap(),
            tw,
            2e-10,
        );
        for i in 0..2 {
            let m = mixture.branches()[0].cell(i).unwrap();
            let f = frozen.branches()[0].states()[i];
            close(g.primitive(m, || false).unwrap().temperature, tg, 2e-10);
            close(
                frozen.gas().primitive(f, || false).unwrap().temperature,
                tg,
                2e-10,
            );
            close(m.energy, f.energy, 3e-11);
        }
        let thermal = fr.thermal.unwrap();
        close(mr.heat_to_gas_j, 2.0 * q, 2e-11);
        close(thermal.source.heat_to_gas_j, 2.0 * q, 2e-11);
        close(mr.wall_energy_change_j, -2.0 * q, 2e-11);
        close(thermal.gas_wall_conservation_residual_j, 0.0, 3e-11);
        assert_eq!(
            mw.shared_wall(0).unwrap().ports,
            vec![port(0, 0.3), port(1, 0.3)]
        );
        assert_eq!(
            fw.shared_wall_conduction(0).unwrap()[0],
            Some(layer.clone())
        );
    }
}

#[test]
fn conductivity_queries_share_the_original_budget_and_late_cancellation_rolls_back() {
    let g = gas(1, false, 0.0);
    let source = phase_source(1.0, 1.0, 4.0);
    let layer = phase_layer(&source, 0.1, WallConductionEvaluation::InitialState);
    let build = |limit, material| {
        let n = mixture_at(&g, 8.0);
        let film = if material {
            0.3
        } else {
            0.25 / (1.0 + 0.25 / 0.3)
        };
        let wall = adiabatic(&n)
            .with_shared_walls(
                &n,
                vec![phase_wall(&source, 4.0, vec![port(0, film), port(1, film)])],
                SharedWallControl {
                    maximum_source_evaluations: limit,
                },
                || false,
            )
            .unwrap();
        let wall = if material {
            wall.with_shared_wall_conduction(vec![vec![Some(layer.clone()); 2]], || false)
                .unwrap()
        } else {
            wall
        };
        (n, wall)
    };
    let (mut n, mut wall) = build(65_536, true);
    let report = wall
        .advance_for(&mut n, Time::new(0.2), controls(), || false)
        .unwrap();
    let work = wall.last_shared_wall_report().unwrap().source_evaluations;
    let (mut constant, mut cw) = build(65_536, false);
    cw.advance_for(&mut constant, Time::new(0.2), controls(), || false)
        .unwrap();
    assert_eq!(
        work,
        cw.last_shared_wall_report().unwrap().source_evaluations + 4,
        "two original initial k queries plus two endpoint domain admissions"
    );
    assert_eq!(fields(&n), fields(&constant));
    assert_eq!(wall_energy(&wall), wall_energy(&cw));
    let (mut exact, mut ew) = build(work, true);
    assert_eq!(
        ew.advance_for(&mut exact, Time::new(0.2), controls(), || false)
            .unwrap(),
        report
    );
    let (mut denied, mut dw) = build(work - 1, true);
    let before = fields(&denied);
    let original = dw.shared_wall(0).unwrap().clone();
    assert!(
        dw.advance_for(&mut denied, Time::new(0.2), controls(), || false)
            .is_err()
    );
    assert_eq!(fields(&denied), before);
    assert_eq!(dw.shared_wall(0).unwrap(), &original);
    assert!(dw.last_report().is_none());
    assert!(dw.last_shared_wall_report().is_none());
    // A later accepted interval may have a different inversion cost. Cancel
    // at the reference interval's final publication poll with ample work.
    let (mut retry, mut rw) = build(65_536, true);
    rw.advance_for(&mut retry, Time::new(0.2), controls(), || false)
        .unwrap();
    let prior = fields(&retry);
    let accepted = rw.shared_wall(0).unwrap().clone();
    let prior_report = rw.last_shared_wall_report();
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
    assert_eq!(fields(&retry), prior);
    assert_eq!(rw.shared_wall(0).unwrap(), &accepted);
    assert_eq!(rw.last_shared_wall_report(), prior_report);
    assert_eq!(rw.last_report(), Some(report));
    assert_eq!(retry.time(), Time::new(0.0));
    assert_eq!(
        rw.advance_for(&mut retry, Time::new(0.2), controls(), || false)
            .unwrap(),
        next
    );
    assert_eq!(fields(&retry), fields(&n));
    assert_eq!(rw.shared_wall(0), wall.shared_wall(0));
}

#[test]
fn post_link_phase_domain_refusal_preserves_both_networks_without_reselecting_material() {
    let g = gas(1, false, 0.0);
    for evaluation in [
        WallConductionEvaluation::InitialState,
        WallConductionEvaluation::ImplicitEndpoint,
    ] {
        let source = phase_source(1.0, 0.001, 4.0);
        let layer = phase_layer(&source, 0.1, evaluation);
        let storage = vec![
            phase_wall(&source, 4.0, vec![port(0, 0.3)]),
            SharedWall {
                storage: SharedWallStorage {
                    thermal: WallThermal::Lumped {
                        temperature_k: 8.0,
                        heat_capacity_j_per_k: 1.0,
                        conductance_w_per_k: 0.0,
                    },
                    heat_storage: None,
                },
                ports: vec![port(1, 0.0)],
            },
        ];
        let layers = vec![vec![Some(layer)], vec![None]];
        let link = SharedWallLink {
            walls: [1, 0],
            conductance_w_per_k: 1.0,
        };
        let mut n = mixture_at(&g, 5.0);
        let bound = |n: &MixtureNetwork1d| {
            adiabatic(n)
                .with_shared_walls(n, storage.clone(), SharedWallControl::default(), || false)
                .unwrap()
                .with_shared_wall_conduction(layers.clone(), || false)
                .unwrap()
        };
        // The gas stage alone remains at h=4 and is supported. Only the
        // subsequent wall link advances phi beyond the pinned narrow claim.
        let mut unlinked = bound(&n);
        unlinked
            .advance_for(&mut n, Time::new(0.01), controls(), || false)
            .unwrap();
        assert_eq!(wall_energy(&unlinked), 4.0);
        let mut n = mixture_at(&g, 5.0);
        let before = fields(&n);
        let mut linked = bound(&n)
            .with_shared_wall_links(vec![link], || false)
            .unwrap();
        assert!(matches!(
            linked.advance_for(&mut n, Time::new(0.01), controls(), || false),
            Err(WallExchangeError::PhaseConduction(_))
        ));
        assert_eq!(fields(&n), before);
        assert_eq!(n.time(), Time::new(0.0));
        for (i, old) in storage.iter().enumerate() {
            assert_eq!(linked.shared_wall(i), Some(old));
        }
        assert!(linked.last_report().is_none());
        assert!(linked.last_shared_wall_report().is_none());
        assert!(linked.shared_wall_link_reports().is_none());

        let mut frozen = frozen_at(&g, 5.0);
        let before = frozen.branches()[0].states().to_vec();
        let mut walls = frozen_walls(&frozen, storage.clone())
            .with_shared_wall_conduction(layers, || false)
            .unwrap()
            .with_shared_wall_links(vec![link], || false)
            .unwrap();
        // Fix this attempted physical interval; automatic step shortening
        // would eventually restore domain support and legitimately accept.
        let error = frozen
            .step_with_walls(
                Time::new(0.01),
                StepControl {
                    maximum_halvings: 0,
                    ..StepControl::default()
                },
                &mut walls,
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
            assert_eq!(walls.shared_wall(i), Some(old));
        }
        assert!(walls.last_report().is_none());
        assert!(walls.last_shared_wall_report().is_none());
        assert!(walls.shared_wall_link_reports().is_none());
    }
}
