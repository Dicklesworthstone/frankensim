//! Constant support rest geometry on the same native contact/preload owners.
use super::*;
use fs_couple::render::schedule::force::coupled::contact::ModalContact;
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::MultiFrictionModalSystem;
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::surfaces::geometry::FrictionGeometry;

fn offset(card: &fs_couple::bowed_string::BowedStringCard, station: f64, rest: f64) -> f64 {
    let k = stiffness(card);
    let g = [
        0.6 * omitted(card, station, 0.8),
        omitted(card, station, 0.67),
    ];
    (k[0] * g[0] + k[1] * g[1]) * rest + (k[1] * g[0] + k[2] * g[1]) * (-0.7 * rest)
}
fn gap(law: &Obstacle, clearance: f64) -> Obstacle {
    Obstacle::new(
        law.collocation().to_vec(),
        1,
        1,
        vec![clearance],
        law.weights().to_vec(),
        law.stiffness(),
        law.alpha(),
        law.provenance().to_owned(),
    )
    .unwrap()
    .with_internal_loss(law.internal_loss())
    .unwrap()
}
// These fixtures use fixed opposite surfaces. The existing complete columns
// therefore expose their support rows directly, without rewriting their maps.
fn supports(system: &MultiFrictionModalSystem, column: &[f64]) -> Vec<ModalAttachment> {
    let mut offset = system.components()[0].modes().len();
    system
        .components()
        .iter()
        .enumerate()
        .skip(1)
        .map(|(component, body)| {
            let end = offset + body.modes().len();
            let row = ModalAttachment {
                component,
                shapes: column[offset..end].to_vec(),
            };
            offset = end;
            row
        })
        .collect()
}

#[test]
fn g1_affine_stop_clearance_retains_initial_energy_and_source_law_then_preloads_and_releases() {
    with_cx(|cx, gate| {
        for rest in [-2e-6, 2e-6] {
            let base = coupled_bank(cx, gate, true, false, rest, 8).unwrap();
            let bridge_energy = base.system().total_energy_j().unwrap();
            let shift = offset(&base.cards()[0], 0.43, rest);
            assert!(shift != 0.0);
            let mut authored = stop(TransverseDirection::try_new([0.6, 0.8]).unwrap());
            // Deliberately close the original physical contact at q=0. This
            // requires a negative effective gap, which must retain its energy.
            let overlap = 2e-7;
            let source_gap = shift - overlap;
            authored.law = gap(&authored.law, source_gap);
            let automatic = multi(base, gate)
                .with_string_stops(vec![authored.clone()], stop_limits(), gate)
                .unwrap();
            let kept = &automatic.source().stops()[0].law;
            assert_eq!(kept.gaps()[0].to_bits(), source_gap.to_bits());
            assert_eq!(kept.provenance(), authored.law.provenance());
            let lowered = &automatic.system().obstacle(0).unwrap().law;
            near(lowered.gaps()[0], -overlap, 1e-15);
            assert!(lowered.gaps()[0] < 0.0);
            assert_eq!(lowered.internal_loss(), authored.law.internal_loss());
            near(
                automatic.system().total_energy_j().unwrap(),
                bridge_energy + 0.5 * authored.law.stiffness() * overlap * overlap,
                1e-17,
            );

            // Independently calculated finite-tail offset, lowered manually
            // onto the existing normal owner. Reuse already-tested reciprocal
            // maps here to isolate affine geometry, energy and continuation.
            let original = automatic.system().obstacle(0).unwrap();
            let reference_contact = ModalContact {
                left: original.left.clone(),
                right: original.right.clone(),
                law: gap(&authored.law, source_gap - shift),
            };
            let column = automatic
                .system()
                .obstacle_normal_column(0)
                .unwrap()
                .to_vec();
            let corrections = supports(automatic.system(), &column);
            let compliance = automatic
                .system()
                .obstacle_residual_compliance()
                .unwrap()
                .to_vec();
            let (_, reference) =
                multi(coupled_bank(cx, gate, true, false, rest, 8).unwrap(), gate).into_parts();
            let mut reference = reference
                .with_obstacle_residual_ports_and_corrections(
                    vec![(reference_contact, authored.limits)],
                    vec![corrections],
                    compliance,
                    stop_limits(),
                    gate,
                )
                .unwrap();
            let (_, mut automatic) = automatic.into_parts();
            let force: Vec<_> = column.iter().map(|b| 0.01 * b).collect();
            let zero = vec![0.0; force.len()];
            let a = automatic
                .initialize_released_equilibrium(&force, gate)
                .unwrap();
            let b = reference
                .initialize_released_equilibrium(&force, gate)
                .unwrap();
            near(a, b, 1e-16);
            let released = [FrictionDrive {
                speed_m_s: 0.0,
                normal_force_n: 0.0,
            }];
            for step in 0..64 {
                let external = if step < 16 { &force } else { &zero };
                let a = automatic.step(external, &released).unwrap().clone();
                let b = reference.step(external, &released).unwrap().clone();
                near(a.observer_pressure_pa, b.observer_pressure_pa, 1e-10);
                near(a.external_work_j, b.external_work_j, 1e-15);
                near(a.network_energy_j, b.network_energy_j, 1e-15);
                near(
                    a.obstacle_residual_energy_j,
                    b.obstacle_residual_energy_j,
                    1e-15,
                );
                for frame in [&a, &b] {
                    assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
                }
                let a = automatic.last_obstacle_frames().unwrap()[0];
                let b = reference.last_obstacle_frames().unwrap()[0];
                near(a.normal_force_n, b.normal_force_n, 1e-9);
                near(a.stored_energy_j, b.stored_energy_j, 1e-15);
                near(a.dissipation_j, b.dissipation_j, 1e-15);
                assert_eq!(a.actuator_work_j, 0.0);
                assert_eq!(a.support_offset_m, 0.0);
                if step == 0 {
                    assert!(a.normal_force_n > 0.0);
                }
            }
        }
    });
}

#[test]
fn g4_bow_rest_offset_cancels_from_incremental_slip_without_adding_work_or_storage() {
    with_cx(|cx, gate| {
        let rest = 2e-6;
        let automatic = multi(coupled_bank(cx, gate, true, false, rest, 8).unwrap(), gate)
            .with_bow_residual_flexibility(gate)
            .unwrap();
        let original = automatic.system().port(0).unwrap();
        let geometry = FrictionGeometry {
            left: original.left.clone(),
            right: original.right.clone(),
        };
        let corrections = supports(
            automatic.system(),
            automatic.system().friction_port_column(0).unwrap(),
        );
        let compliance = automatic
            .system()
            .friction_residual_compliance()
            .unwrap()
            .to_vec();
        let (_, reference) =
            multi(coupled_bank(cx, gate, true, false, rest, 8).unwrap(), gate).into_parts();
        let mut reference = reference
            .with_friction_residual_ports_and_corrections(
                vec![geometry],
                vec![corrections],
                compliance,
                gate,
            )
            .unwrap();
        let (_, mut automatic) = automatic.into_parts();
        assert_eq!(
            automatic.total_energy_j().unwrap(),
            reference.total_energy_j().unwrap()
        );
        let zero = vec![0.0; automatic.mode_count()];
        for _ in 0..64 {
            let a = automatic.step(&zero, &[drive()]).unwrap().clone();
            let b = reference.step(&zero, &[drive()]).unwrap().clone();
            assert_eq!(a, b);
            assert!(a.energy_residual_j.abs() <= a.energy_tolerance_j);
        }
    });
}
