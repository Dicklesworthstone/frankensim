//! The block heat cache against the plain exact stage (Gauntlet G2, pexfp.270 layer 1). Two identical island networks
//! are stepped side by side: one with the blocks' own `HeatCache`, one with every block's cache disabled, which
//! rebuilds the positive-heat factor on every call exactly as the code did before the cache existed. Every output
//! is compared with `to_bits`, including after a change of a few ulps planted in one block's resistance.
use super::*;
use crate::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_math::c64::C64;

/// The block that receives the planted resistance change. Block 1 holds the heavily damped modes (zeta 0.2): its
/// component loss dominates the component loss summed over blocks, so a few ulps planted there are not rounded away
/// by that sum. Block 0 (zeta 0.01) contributes roughly 1e-4 of it, and a planted change there never reaches the sum.
const PLANTED: usize = 1;

fn bits(a: &[f64]) -> Vec<u64> {
    a.iter().map(|v| v.to_bits()).collect()
}

fn four_models(rate: u32, omega: f64) -> Vec<ModalAcousticTimeModel> {
    let base: Vec<_> = [(omega, 0.01), (1.7 * omega, 0.2)]
        .into_iter()
        .enumerate()
        .map(|(i, (w, zeta))| {
            let mut model = ModalAcousticTimeModel::try_new(
                rate,
                vec![ModalAcousticMode {
                    angular_frequency_rad_s: w,
                    damping_ratio: zeta,
                    pressure_per_modal_velocity: C64::ZERO,
                }],
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap();
            model
                .restore_states(&[ModalAcousticState {
                    displacement_m_sqrt_kg: if i == 0 { 1e-8 } else { -2e-8 },
                    velocity_m_sqrt_kg_per_s: 0.0,
                }])
                .unwrap();
            model
        })
        .collect();
    base.iter().cloned().cycle().take(4).collect()
}

fn assert_same_record(a: &StepRecord, b: &StepRecord, k: usize) {
    assert_eq!(bits(&a.x), bits(&b.x), "x at step {k}");
    assert_eq!(bits(&a.y), bits(&b.y), "y at step {k}");
    assert_eq!(a.delta_h.to_bits(), b.delta_h.to_bits(), "delta_h at step {k}");
    assert_eq!(a.dissipated.to_bits(), b.dissipated.to_bits(), "dissipated at step {k}");
    assert_eq!(a.supplied.to_bits(), b.supplied.to_bits(), "supplied at step {k}");
    assert_eq!(a.newton_iters, b.newton_iters, "newton_iters at step {k}");
    assert_eq!(
        a.solver_residual.to_bits(),
        b.solver_residual.to_bits(),
        "solver_residual at step {k}"
    );
}

fn records_differ(a: &StepRecord, b: &StepRecord) -> bool {
    bits(&a.x) != bits(&b.x)
        || bits(&a.y) != bits(&b.y)
        || a.delta_h.to_bits() != b.delta_h.to_bits()
        || a.dissipated.to_bits() != b.dissipated.to_bits()
        || a.supplied.to_bits() != b.supplied.to_bits()
}

/// The held external force of step `k`: it changes every step for 20 steps (every one a transition miss), is held to
/// step 40 (hits), then the dominant block's coordinate 1 changes by `fulps` ulps and is held to step 50, and is zero
/// after (the rest force only).
fn force(k: usize, fulps: u64) -> [f64; 4] {
    let varying = |k: usize| ((k * 7 % 13) as f64) / 13.0 - 0.5;
    let c = if k < 20 {
        varying(k)
    } else if k < 50 {
        varying(19)
    } else {
        0.0
    };
    let mut f = [1e-4 * c, -2e-4 * c, 3e-4 * c, -4e-4 * c];
    if (40..50).contains(&k) {
        f[1] = f64::from_bits(f[1].to_bits() + fulps);
    }
    f
}

#[test]
fn the_block_heat_cache_replays_the_plain_exact_stage_bit_for_bit() {
    let rate = 48_000;
    let mut models = four_models(rate, 1_000.0);
    let columns = [vec![1.0, 0.0, -0.7, 0.0], vec![0.0, -0.3, 0.0, 1.0]];
    let links = [
        ConnectionLaw {
            stiffness_n_m: 2e5,
            damping_n_s_m: 70.0,
            rest_extension_m: 1e-9,
            constant_loss_reference_omega_rad_s: None,
        },
        ConnectionLaw {
            stiffness_n_m: 3e5,
            damping_n_s_m: 90.0,
            rest_extension_m: -2e-9,
            constant_loss_reference_omega_rad_s: None,
        },
    ];
    let gate = CancelGate::new();
    let mut cached = ExactNetwork::new(&models, &columns, &links, 4, &gate).unwrap();
    let mut plain = ExactNetwork::new(&models, &columns, &links, 4, &gate).unwrap();
    // The control is never planted: it shows the planted change really moves a loss, so the negative can fail.
    let mut control = ExactNetwork::new(&models, &columns, &links, 4, &gate).unwrap();
    for network in [&mut plain, &mut control] {
        for (_, block) in &mut network.blocks {
            block.heat = HeatCache::disabled();
        }
    }
    assert_eq!(cached.blocks.len(), 2, "two noncontiguous islands");
    let dt = 1.0 / f64::from(rate);
    // The square roots that build a factor can absorb a single ulp, so find the smallest planted change (in ulps on
    // every nonzero entry of the planted block's component resistance) that really moves the component loss at the starting
    // state. A negative that cannot fail would show nothing.
    let start_force = [1e-4, -2e-4, 3e-4, -4e-4];
    let (start_record, start_x, start_u) =
        control.advance(&models, &[], &[], &start_force, dt).unwrap();
    let start_loss = control
        .losses(&start_x, &start_u, dt, start_record.dissipated)
        .unwrap()
        .0;
    let ulps = (1..=64_u64)
        .find(|&ulps| {
            let mut probe = ExactNetwork::new(&models, &columns, &links, 4, &gate).unwrap();
            for v in probe.blocks[PLANTED].1.component_r.iter_mut().filter(|v| **v != 0.0) {
                *v = f64::from_bits(v.to_bits() + ulps);
            }
            let loss = probe
                .losses(&start_x, &start_u, dt, start_record.dissipated)
                .unwrap()
                .0;
            loss.to_bits() != start_loss.to_bits()
        })
        .expect("some planted ulp count moves the component loss");
    // The same for the force: the smallest ulp count on the dominant block's coordinate-1 force that moves the
    // record at the starting state, so that serving the previous force's transition would be caught.
    let held = force(25, 0);
    let (held_record, _, _) = control.advance(&models, &[], &[], &held, dt).unwrap();
    let fulps = (1..=64_u64)
        .find(|&fulps| {
            let mut changed = held;
            changed[1] = f64::from_bits(held[1].to_bits() + fulps);
            let (record, _, _) = control.advance(&models, &[], &[], &changed, dt).unwrap();
            records_differ(&held_record, &record)
        })
        .expect("some planted force ulp count moves the record");
    let contact = [1.0, -2.0, 0.3, 0.4];
    let steps = 60;
    let plant_at = 30;
    let mut builds_before_plant = 0;
    let mut moved = 0;
    let mut force_moved = 0;
    for k in 0..steps {
        if k == plant_at {
            builds_before_plant = cached.blocks[PLANTED].1.heat.builds();
            // The planted negative: `ulps` ulps on every nonzero entry of one block's component resistance, on both
            // twins. The cached twin must rebuild that one factor and keep following the plain twin.
            for network in [&mut cached, &mut plain] {
                let resistance = &mut network.blocks[PLANTED].1.component_r;
                assert!(resistance.iter().any(|v| *v != 0.0), "component loss");
                for v in resistance.iter_mut().filter(|v| **v != 0.0) {
                    *v = f64::from_bits(v.to_bits() + ulps);
                }
            }
        }
        // The held force changes every step for 40 steps, then is zero (the rest force only).
        let external = force(k, fulps);
        let (rc, xc, uc) = cached.advance(&models, &[], &[], &external, dt).unwrap();
        let (rp, xp, up) = plain.advance(&models, &[], &[], &external, dt).unwrap();
        assert_same_record(&rc, &rp, k);
        if (40..50).contains(&k) {
            // The stale-key negative: the plain stage fed the unchanged force at this same state is what a cache that
            // ignored the force would have served. It must differ on some step, or the key proves nothing.
            let (stale, _, _) = control.advance(&models, &[], &[], &force(39, fulps), dt).unwrap();
            if records_differ(&rp, &stale) {
                force_moved += 1;
            }
        }
        assert_eq!(bits(&xc), bits(&xp), "start state at step {k}");
        assert_eq!(bits(&uc), bits(&up), "held input at step {k}");
        let lc = cached.losses(&xc, &uc, dt, rc.dissipated).unwrap();
        let lp = plain.losses(&xp, &up, dt, rp.dissipated).unwrap();
        for (a, b) in [(lc.0, lp.0), (lc.1, lp.1), (lc.2, lp.2), (lc.3, lp.3)] {
            assert_eq!(a.to_bits(), b.to_bits(), "losses at step {k}");
        }
        let lk = control.losses(&xp, &up, dt, rp.dissipated).unwrap();
        if k >= plant_at && lp.0.to_bits() != lk.0.to_bits() {
            moved += 1;
        }
        let wc = cached
            .connection_work(&xc, &uc, &rc.x, dt, &columns, &links, lc.1, Some(&gate))
            .unwrap();
        let wp = plain
            .connection_work(&xp, &up, &rp.x, dt, &columns, &links, lp.1, Some(&gate))
            .unwrap();
        assert_eq!(wc.len(), wp.len());
        for (a, b) in wc.iter().zip(&wp) {
            assert_eq!(a.energy_j.to_bits(), b.energy_j.to_bits(), "connection work at step {k}");
        }
        let ec = cached.response(&contact, dt).unwrap();
        let ep = plain.response(&contact, dt).unwrap();
        assert_eq!(bits(&ec), bits(&ep), "contact response at step {k}");
        for (i, model) in models.iter_mut().enumerate() {
            model
                .restore_states(&[ModalAcousticState {
                    displacement_m_sqrt_kg: rc.x[2 * i],
                    velocity_m_sqrt_kg_per_s: rc.x[2 * i + 1],
                }])
                .unwrap();
        }
    }
    // The negatives can fail: the planted ulps moved the component loss against the unplanted control, and the
    // planted force change moved the record against the previous force's.
    assert!(moved > 0, "the planted {ulps}-ulp change never moved the component loss");
    assert!(force_moved > 0, "the planted {fulps}-ulp force change never moved the record");
    // Changing the resistance rebuilt exactly one factor of the planted block (its component channel), once.
    assert_eq!(cached.blocks[PLANTED].1.heat.builds(), builds_before_plant + 1);
    for (_, block) in &cached.blocks {
        // The distinct keys of a block are its system R, the component channel (twice in the planted block), the
        // connection channel and one dashpot row per link: six at most here. Each is built once and then reused on
        // every other call, so builds stay within the cache capacity and reuses dominate them.
        assert!(block.heat.builds() <= 8, "builds {}", block.heat.builds());
        assert!(block.heat.len() <= 8, "entries {}", block.heat.len());
        assert!(
            block.heat.reuses() > 20 * block.heat.builds(),
            "builds {} reuses {}",
            block.heat.builds(),
            block.heat.reuses()
        );
    }
    for (_, block) in &cached.blocks {
        // The transitions: a build for each changed force (and the contact-response column, which competes for the
        // 8 slots while the force changes), then every later call of a held force reuses. Reuses dominate.
        assert!(block.heat.transition_builds() <= 60, "transition builds {}", block.heat.transition_builds());
        assert!(block.heat.transition_len() <= 8, "transition entries {}", block.heat.transition_len());
        assert!(
            block.heat.transition_reuses() > 2 * block.heat.transition_builds(),
            "transition builds {} reuses {}",
            block.heat.transition_builds(),
            block.heat.transition_reuses()
        );
    }
    for (_, block) in &plain.blocks {
        assert_eq!((block.heat.reuses(), block.heat.len()), (0, 0));
        assert_eq!((block.heat.transition_reuses(), block.heat.transition_len()), (0, 0));
        assert!(block.heat.builds() > 100);
        assert!(block.heat.transition_builds() > 100);
    }
}

/// (checks, reuses) of each block's component and connection verdicts.
fn verdict_counts(network: &ExactNetwork) -> Vec<[(u64, u64); 2]> {
    network
        .blocks
        .iter()
        .map(|(_, block)| {
            let v = &block.verdicts;
            [
                (v.component.checks(), v.component.reuses()),
                (v.connection.checks(), v.connection.reuses()),
            ]
        })
        .collect()
}

// frankensim-3sw19: the blocks' PSD verdicts against verdicts disabled, which run the eigendecomposition on every call
// as the code did before the slots existed. Every record and all four losses are compared with to_bits. A one-ulp
// change planted in one block's component resistance re-checks that channel once. An indefinite component resistance
// planted in the warm slot is refused exactly as the disabled twin refuses it, on every attempt.
#[test]
fn the_block_psd_verdicts_replay_the_plain_loss_channels_bit_for_bit_and_never_admit_a_refusal() {
    let rate = 48_000;
    let mut models = four_models(rate, 1_000.0);
    let columns = [vec![1.0, 0.0, -0.7, 0.0], vec![0.0, -0.3, 0.0, 1.0]];
    let links = [
        ConnectionLaw {
            stiffness_n_m: 2e5,
            damping_n_s_m: 70.0,
            rest_extension_m: 1e-9,
            constant_loss_reference_omega_rad_s: None,
        },
        ConnectionLaw {
            stiffness_n_m: 3e5,
            damping_n_s_m: 90.0,
            rest_extension_m: -2e-9,
            constant_loss_reference_omega_rad_s: None,
        },
    ];
    let gate = CancelGate::new();
    let mut slots = ExactNetwork::new(&models, &columns, &links, 4, &gate).unwrap();
    let mut plain = ExactNetwork::new(&models, &columns, &links, 4, &gate).unwrap();
    for (_, block) in &mut plain.blocks {
        block.verdicts = ChannelVerdicts::disabled();
    }
    assert_eq!(slots.blocks.len(), 2, "two noncontiguous islands");
    let dt = 1.0 / f64::from(rate);
    let (steps, plant_at) = (60_u64, 30_u64);
    let start = verdict_counts(&slots);
    let restore = |models: &mut [ModalAcousticTimeModel], x: &[f64]| {
        for (i, model) in models.iter_mut().enumerate() {
            model
                .restore_states(&[ModalAcousticState {
                    displacement_m_sqrt_kg: x[2 * i],
                    velocity_m_sqrt_kg_per_s: x[2 * i + 1],
                }])
                .unwrap();
        }
    };
    for k in 0..steps {
        if k == plant_at {
            for network in [&mut slots, &mut plain] {
                let resistance = &mut network.blocks[PLANTED].1.component_r;
                assert!(resistance.iter().any(|v| *v != 0.0), "component loss");
                for v in resistance.iter_mut().filter(|v| **v != 0.0) {
                    *v = f64::from_bits(v.to_bits() + 1);
                }
            }
        }
        let external = force(k as usize, 1);
        let (rs, xs, us) = slots.advance(&models, &[], &[], &external, dt).unwrap();
        let (rp, xp, up) = plain.advance(&models, &[], &[], &external, dt).unwrap();
        assert_same_record(&rs, &rp, k as usize);
        let ls = slots.losses(&xs, &us, dt, rs.dissipated).unwrap();
        let lp = plain.losses(&xp, &up, dt, rp.dissipated).unwrap();
        for (a, b) in [(ls.0, lp.0), (ls.1, lp.1), (ls.2, lp.2), (ls.3, lp.3)] {
            assert_eq!(a.to_bits(), b.to_bits(), "losses at step {k}");
        }
        restore(&mut models, &rs.x);
    }
    // One loss evaluation per block per step: the planted block's component channel checked once more after the plant,
    // every other call admitted from its slot. The disabled twin checked on every call.
    for (i, (now, then)) in verdict_counts(&slots).iter().zip(&start).enumerate() {
        let component_checks = if i == PLANTED { 2 } else { 1 };
        let component = (now[0].0 - then[0].0, now[0].1 - then[0].1);
        let connection = (now[1].0 - then[1].0, now[1].1 - then[1].1);
        assert_eq!(
            component,
            (component_checks, steps - component_checks),
            "block {i} component"
        );
        assert_eq!(connection, (1, steps - 1), "block {i} connection");
    }
    for (_, block) in &plain.blocks {
        assert_eq!(
            block.verdicts.component.reuses() + block.verdicts.connection.reuses(),
            0
        );
        assert!(
            block.verdicts.component.checks() >= steps
                && block.verdicts.connection.checks() >= steps
        );
    }
    // Planted negative: the warm slot's component resistance made indefinite (its largest diagonal entry negated).
    let external = force(steps as usize, 1);
    let (rs, xs, us) = slots.advance(&models, &[], &[], &external, dt).unwrap();
    let (rp, xp, up) = plain.advance(&models, &[], &[], &external, dt).unwrap();
    for network in [&mut slots, &mut plain] {
        let resistance = &mut network.blocks[PLANTED].1.component_r;
        let n = (1..=resistance.len())
            .find(|n| n * n == resistance.len())
            .unwrap();
        let i = (0..n)
            .max_by(|&a, &b| resistance[a * n + a].total_cmp(&resistance[b * n + b]))
            .unwrap();
        assert!(
            resistance[i * n + i] > 0.0,
            "a positive diagonal entry to negate"
        );
        resistance[i * n + i] = -resistance[i * n + i];
    }
    let held = verdict_counts(&slots)[PLANTED][0];
    for attempt in 0..2 {
        let refused = slots.losses(&xs, &us, dt, rs.dissipated).unwrap_err();
        let expected = plain.losses(&xp, &up, dt, rp.dissipated).unwrap_err();
        assert!(
            matches!(
                refused,
                ModalCouplingError::Propagation(fs_phs::PhsError::NotPsd {
                    what: "loss channel R"
                })
            ),
            "attempt {attempt}: {refused:?}"
        );
        assert_eq!(
            format!("{refused:?}"),
            format!("{expected:?}"),
            "attempt {attempt}"
        );
    }
    assert_eq!(
        verdict_counts(&slots)[PLANTED][0],
        (held.0 + 2, held.1),
        "each refusal ran the eigendecomposition; none was admitted from the slot"
    );
}
