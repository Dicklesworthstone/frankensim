//! G1/G3/G4: exact momentum jumps and continuing coupled motion, not a held pulse.
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::plate::impact::radiation::{Model, Pole};
use fs_couple::render::schedule::ScheduledRenderer;
use fs_couple::render::schedule::force::coupled::render::CoupledModalVoice;
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalColumnConnection, ModalConnection,
    ModalCouplingConfig, ModalCouplingError, PreparedModalBudget,
};
use fs_couple::render::schedule::force::impulse::ports::{PhysicalImpulseBudget, PortImpulse};
use fs_couple::render::schedule::force::impulse::{ImpulseBudget, ImpulseRenderer, ModalImpulse};
use fs_couple::render::schedule::force::{
    ForceInitialization, ForceRenderConfig, ModalForceEvent, ModalForceVoice,
};
use fs_couple::render::{
    GatedRenderOutcome, ModalStringVoice, RenderContext, RenderError, RenderVoice,
};
use fs_exec::CancelGate;
use fs_math::c64::C64;

fn config() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 8,
        max_connections: 8,
        max_setup_terms: 100_000,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 100.0,
        maximum_abs_pressure_pa: 1e5,
        maximum_abs_connection_force_n: 1e5,
        solve_relative_tolerance: 1e-11,
        energy_absolute_tolerance_j: 1e-12,
        energy_relative_tolerance: 1e-10,
    }
}

fn prepared() -> PreparedModalBudget {
    PreparedModalBudget {
        max_modes: 8,
        max_dense_bytes: 1 << 24,
        max_multiply_terms: 1_000_000_000,
    }
}

fn model(
    rate: u32,
    omega: f64,
    q: f64,
    v: f64,
    budget: ModalAcousticTimeBudget,
) -> ModalAcousticTimeModel {
    let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
        rate,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: omega,
            damping_ratio: 0.0,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        budget,
    )
    .unwrap();
    model
        .restore_states(&[ModalAcousticState {
            displacement_m_sqrt_kg: q,
            velocity_m_sqrt_kg_per_s: v,
        }])
        .unwrap();
    model
}

fn states(system: &CoupledModalSystem) -> Vec<ModalAcousticState> {
    system
        .components()
        .iter()
        .flat_map(|m| m.states().iter().copied())
        .collect()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

fn counterwork_model(large_velocity: f64) -> ModalAcousticTimeModel {
    let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
        2000,
        vec![
            ModalAcousticMode {
                angular_frequency_rad_s: 0.0,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::ZERO,
            };
            3
        ],
        ModalAcousticTimeBudget {
            maximum_abs_displacement_m_sqrt_kg: 1e20,
            maximum_abs_velocity_m_sqrt_kg_per_s: 1e20,
            maximum_total_energy_j: 1e40,
            ..ModalAcousticTimeBudget::audible_reference()
        },
    )
    .unwrap();
    model
        .restore_states(
            &[large_velocity, 0.0, large_velocity].map(|velocity_m_sqrt_kg_per_s| {
                ModalAcousticState {
                    displacement_m_sqrt_kg: 0.0,
                    velocity_m_sqrt_kg_per_s,
                }
            }),
        )
        .unwrap();
    model
}

#[test]
fn coupled_work_preserves_quadratic_energy_before_midpoint_and_product_rounding() {
    let gate = CancelGate::new_clock_free();
    // First case: the first endpoint mean falls halfway between f64 values.
    // Second case: both large means are exact, but multiplying either by six
    // rounds away two joules. Compensating only the final sum fixes neither.
    for (large, kick, expected_work) in [
        ((1u64 << 55) as f64, 8.0, 64.5),
        ((1u64 << 52) as f64, 6.0, 36.5),
    ] {
        for use_prepared in [false, true] {
            let mut cfg = config();
            cfg.maximum_total_energy_j = 1e40;
            let models = vec![counterwork_model(large)];
            let mut system = if use_prepared {
                CoupledModalSystem::new_prepared(models, vec![], cfg, prepared(), &gate)
            } else {
                CoupledModalSystem::new(models, vec![], cfg, &gate)
            }
            .unwrap();
            let before = states(&system);
            let impulse = [kick, 1.0, -kick];
            let receipt = system.apply_modal_impulse(&impulse, &gate).unwrap();
            // Unit modal masses: v dot J cancels, leaving |J|²/2 exactly.
            assert_eq!(receipt.kinetic_energy_change_j, expected_work);
            assert_eq!(receipt.external_work_j, expected_work);
            assert_eq!(receipt.work_residual_j, 0.0);
            assert_eq!(receipt.maximum_relative_momentum_residual, 0.0);
            for ((after, old), input) in states(&system).iter().zip(&before).zip(impulse) {
                assert_eq!(after.displacement_m_sqrt_kg, old.displacement_m_sqrt_kg);
                assert_eq!(
                    after.velocity_m_sqrt_kg_per_s - old.velocity_m_sqrt_kg_per_s,
                    input
                );
            }
            // Reversing the actual transfer extracts the same energy rather
            // than manufacturing another positive quadratic-energy deposit.
            let reverse = system
                .apply_modal_impulse(&impulse.map(|j| -j), &gate)
                .unwrap();
            assert_eq!(reverse.kinetic_energy_change_j, -expected_work);
            assert_eq!(reverse.external_work_j, -expected_work);
            assert_eq!(reverse.work_residual_j, 0.0);
            assert_eq!(states(&system), before);
            assert_eq!(system.samples_rendered(), 0);
            assert!(system.last_frame().is_none());
        }
    }
}

#[test]
fn independent_scheduled_work_preserves_midpoint_and_product_tails() {
    for (large, kick, expected_work) in [
        ((1u64 << 55) as f64, 8.0, 64.5),
        ((1u64 << 52) as f64, 6.0, 36.5),
    ] {
        let voice = RenderVoice::ModalString(
            ModalStringVoice::new(counterwork_model(large), vec![0.0; 3]).unwrap(),
        );
        let events = [1.0, -1.0]
            .into_iter()
            .enumerate()
            .flat_map(|(sample, sign)| {
                [kick, 1.0, -kick]
                    .into_iter()
                    .enumerate()
                    .map(move |(mode, input)| ModalImpulse {
                        sample: sample as u64,
                        voice: 0,
                        mode,
                        impulse_n_s_per_sqrt_kg: sign * input,
                    })
            })
            .collect();
        let mut render = ScheduledRenderer::new(RenderContext::new(vec![voice], 2), vec![], 0)
            .unwrap()
            .with_modal_impulses(
                events,
                ImpulseBudget {
                    max_events: 6,
                    max_staged_modes: 6,
                },
            )
            .unwrap();
        render.block(&mut [0.0; 2]).unwrap();
        assert_eq!(
            render
                .receipts()
                .iter()
                .map(|receipt| receipt.kinetic_energy_change_j)
                .collect::<Vec<_>>(),
            [expected_work, -expected_work]
        );
        assert_eq!(
            render
                .modal_states(0)
                .unwrap()
                .iter()
                .map(|state| state.velocity_m_sqrt_kg_per_s)
                .collect::<Vec<_>>(),
            [large, 0.0, large]
        );
    }
}

#[test]
fn reciprocal_physical_impulse_matches_two_mass_solution_at_two_sample_rates() {
    let gate = CancelGate::new_clock_free();
    for rate in [2000, 4000] {
        for use_prepared in [false, true] {
            let models = vec![4.0, 1.0]
                .into_iter()
                .map(|mass| {
                    ModalAcousticTimeModel::try_free_mass(
                        rate,
                        mass,
                        0.0,
                        0.5,
                        ModalAcousticTimeBudget::audible_reference(),
                    )
                    .unwrap()
                })
                .collect();
            let links = vec![ModalConnection {
                left: ModalAttachment {
                    component: 0,
                    shapes: vec![0.5],
                },
                right: ModalAttachment {
                    component: 1,
                    shapes: vec![1.0],
                },
                stiffness_n_m: 400.0,
                damping_n_s_m: 0.0,
                rest_extension_m: 0.0,
            }];
            let mut system = if use_prepared {
                CoupledModalSystem::new_prepared(models, links, config(), prepared(), &gate)
            } else {
                CoupledModalSystem::new(models, links, config(), &gate)
            }
            .unwrap();
            // +0.4 N s on the 4 kg body, -0.4 N s on the 1 kg body.
            let receipt = system.apply_modal_impulse(&[0.2, -0.4], &gate).unwrap();
            close(receipt.external_work_j, 0.1);
            close(receipt.kinetic_energy_change_j, 0.1);
            assert_eq!(receipt.sample, 0);
            assert_eq!(system.samples_rendered(), 0);
            assert!(
                states(&system)
                    .iter()
                    .all(|s| s.displacement_m_sqrt_kg == 0.0)
            );
            assert!(system.last_frame().is_none());
            let omega = 500.0_f64.sqrt();
            for tick in 1..=rate / 50 {
                system.step(&[0.0, 0.0]).unwrap();
                let t = f64::from(tick) / f64::from(rate);
                let relative = 0.5 * (omega * t).sin() / omega;
                let relative_rate = 0.5 * (omega * t).cos();
                let s = states(&system);
                close(s[0].displacement_m_sqrt_kg / 2.0, 0.5 * t + 0.2 * relative);
                close(s[1].displacement_m_sqrt_kg, 0.5 * t - 0.8 * relative);
                close(
                    s[0].velocity_m_sqrt_kg_per_s / 2.0,
                    0.5 + 0.2 * relative_rate,
                );
                close(s[1].velocity_m_sqrt_kg_per_s, 0.5 - 0.8 * relative_rate);
                close(
                    2.0 * s[0].velocity_m_sqrt_kg_per_s + s[1].velocity_m_sqrt_kg_per_s,
                    2.5,
                );
                close(system.total_energy_j().unwrap(), 0.725);
            }
        }
    }
}

#[test]
fn continuing_impulse_preserves_material_acoustic_and_interval_histories() {
    let gate = CancelGate::new_clock_free();
    for use_prepared in [false, true] {
        let models = vec![model(
            2000,
            40.0,
            0.01,
            0.1,
            ModalAcousticTimeBudget::audible_reference(),
        )];
        let system = if use_prepared {
            CoupledModalSystem::new_prepared(models, vec![], config(), prepared(), &gate)
        } else {
            CoupledModalSystem::new(models, vec![], config(), &gate)
        }
        .unwrap();
        let mut system = system
            .with_material_relaxation(
                vec![fs_phs::RelaxationBranch {
                    projection: vec![1.0, 0.0],
                    stiffness: 200.0,
                    relaxation_time_s: 0.05,
                }],
                Some(vec![0.02]),
                &gate,
            )
            .unwrap()
            .with_radiation_load(
                &Model {
                    ports: 1,
                    poles: vec![Pole {
                        omega: 60.0,
                        zeta: 0.1,
                        coupling: vec![10.0],
                    }],
                },
                &[0],
                1,
                &gate,
            )
            .unwrap();
        for _ in 0..3 {
            system.step(&[0.0]).unwrap();
        }
        let before = states(&system);
        let memory = system.material_memory().to_vec();
        let radiation = system.radiation_state().to_vec();
        let frame = system.last_frame().unwrap().clone();
        let energy = system.total_energy_j().unwrap();
        assert!(radiation.iter().any(|x| *x != 0.0));
        let receipt = system.apply_modal_impulse(&[0.125], &gate).unwrap();
        assert_eq!(receipt.sample, 3);
        assert_eq!(system.samples_rendered(), 3);
        assert_eq!(
            states(&system)[0].displacement_m_sqrt_kg,
            before[0].displacement_m_sqrt_kg
        );
        assert_eq!(
            states(&system)[0].velocity_m_sqrt_kg_per_s,
            before[0].velocity_m_sqrt_kg_per_s + 0.125
        );
        assert_eq!(system.material_memory(), memory);
        assert_eq!(system.radiation_state(), radiation);
        assert_eq!(system.last_frame(), Some(&frame));
        close(
            receipt.external_work_j,
            0.125 * before[0].velocity_m_sqrt_kg_per_s + 0.5 * 0.125_f64.powi(2),
        );
        close(
            system.total_energy_j().unwrap() - energy,
            receipt.external_work_j,
        );
        let next = system.step(&[0.0]).unwrap();
        assert_eq!(next.sample, 4);
        close(
            next.total_energy_j() + next.total_dissipation_j() - energy,
            receipt.external_work_j,
        );
    }
}

#[test]
fn component_refusal_and_cancellation_publish_no_partial_impulse() {
    let gate = CancelGate::new_clock_free();
    let create = || {
        let normal = ModalAcousticTimeBudget::audible_reference();
        let bounded = ModalAcousticTimeBudget {
            maximum_abs_velocity_m_sqrt_kg_per_s: 1.0,
            ..normal
        };
        CoupledModalSystem::new(
            vec![
                model(2000, 0.0, 0.0, 0.0, normal),
                model(2000, 0.0, 0.0, 0.0, bounded),
            ],
            vec![],
            config(),
            &gate,
        )
        .unwrap()
    };
    let mut actual = create();
    let before = states(&actual);
    assert!(matches!(
        actual.apply_modal_impulse(&[0.5, 2.0], &gate),
        Err(ModalCouplingError::Component { component: 1, .. })
    ));
    assert_eq!(states(&actual), before);
    for input in [&[0.25][..], &[0.25, f64::NAN][..]] {
        assert!(actual.apply_modal_impulse(input, &gate).is_err());
        assert_eq!(states(&actual), before);
    }
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        actual.apply_modal_impulse(&[0.25, 0.5], &cancelled),
        Err(ModalCouplingError::Cancelled)
    ));
    assert_eq!(states(&actual), before);
    let mut reference = create();
    assert_eq!(
        actual.apply_modal_impulse(&[0.25, 0.5], &gate).unwrap(),
        reference.apply_modal_impulse(&[0.25, 0.5], &gate).unwrap()
    );
    assert_eq!(
        actual.step(&[0.0, 0.0]).unwrap(),
        reference.step(&[0.0, 0.0]).unwrap()
    );
}

#[test]
fn total_energy_pressure_and_connection_force_limits_apply_at_the_jump() {
    let gate = CancelGate::new_clock_free();
    for limit_kind in 0..3 {
        let mut mode = model(
            2000,
            0.0,
            0.0,
            0.0,
            ModalAcousticTimeBudget::audible_reference(),
        );
        if limit_kind == 1 {
            mode = ModalAcousticTimeModel::try_new_with_free_coordinates(
                2000,
                vec![ModalAcousticMode {
                    angular_frequency_rad_s: 40.0,
                    damping_ratio: 0.0,
                    pressure_per_modal_velocity: C64::new(1.0, 0.0),
                }],
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap();
        }
        let mut cfg = config();
        match limit_kind {
            0 => cfg.maximum_total_energy_j = 0.01,
            1 => cfg.maximum_abs_pressure_pa = 0.1,
            _ => cfg.maximum_abs_connection_force_n = 0.1,
        }
        let mut system = CoupledModalSystem::new_with_columns(
            vec![mode],
            vec![ModalColumnConnection {
                column: vec![1.0],
                stiffness_n_m: 0.0,
                damping_n_s_m: 1.0,
                rest_extension_m: 0.0,
            }],
            cfg,
            &gate,
        )
        .unwrap();
        let before = states(&system);
        assert!(matches!(
            system.apply_modal_impulse(&[0.25], &gate),
            Err(ModalCouplingError::Budget { .. })
        ));
        assert_eq!(states(&system), before);
        assert_eq!(system.samples_rendered(), 0);
    }
}

#[test]
fn momentum_loss_cannot_hide_behind_large_energy_or_cancel_across_modes() {
    let gate = CancelGate::new_clock_free();
    let wide = ModalAcousticTimeBudget {
        maximum_abs_velocity_m_sqrt_kg_per_s: 1e20,
        maximum_total_energy_j: 1e40,
        ..ModalAcousticTimeBudget::audible_reference()
    };
    let mut cfg = config();
    cfg.maximum_total_energy_j = 1e40;
    let mut system = CoupledModalSystem::new(
        vec![
            model(2000, 0.0, 0.0, 1e16, wide),
            model(2000, 0.0, 0.0, -1e16, wide),
            model(2000, 0.0, 0.0, 0.0, wide),
        ],
        vec![],
        cfg,
        &gate,
    )
    .unwrap();
    let before = states(&system);
    assert!(matches!(
        system.apply_modal_impulse(&[1.0, -1.0, 0.0], &gate),
        Err(ModalCouplingError::ImpulseResidual { mode: 0, .. })
    ));
    assert_eq!(states(&system), before);
    let receipt = system.apply_modal_impulse(&[0.0, 0.0, 1.0], &gate).unwrap();
    assert_eq!(receipt.kinetic_energy_change_j, 0.5);
    assert_eq!(receipt.external_work_j, 0.5);
    assert_eq!(receipt.work_residual_j, 0.0);
}

fn render_models() -> Vec<ModalAcousticTimeModel> {
    [40.0, 70.0]
        .into_iter()
        .map(|omega| {
            ModalAcousticTimeModel::try_new(
                2000,
                vec![ModalAcousticMode {
                    angular_frequency_rad_s: omega,
                    damping_ratio: 0.01,
                    pressure_per_modal_velocity: C64::new(1.0, 0.1),
                }],
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap()
        })
        .collect()
}

fn render_links() -> Vec<ModalConnection> {
    vec![ModalConnection {
        left: ModalAttachment {
            component: 0,
            shapes: vec![0.5],
        },
        right: ModalAttachment {
            component: 1,
            shapes: vec![1.0],
        },
        stiffness_n_m: 500.0,
        damping_n_s_m: 0.1,
        rest_extension_m: 0.0,
    }]
}

fn render_inputs() -> Vec<PortImpulse> {
    [
        (37, 0, -0.05),
        (0, 0, 0.25),
        (37, 1, 0.025),
        (37, 0, 0.01),
        (127, 1, -0.03),
    ]
    .into_iter()
    .map(|(sample, voice, impulse_n_s)| PortImpulse {
        sample,
        voice,
        port: 0,
        impulse_n_s,
    })
    .collect()
}

fn performance(max_staged_modes: usize) -> Result<ImpulseRenderer, RenderError> {
    let voices = render_models()
        .into_iter()
        .enumerate()
        .map(|(component, model)| {
            ModalForceVoice::new(
                model,
                vec![vec![if component == 0 { 0.5 } else { 1.0 }]],
                vec![if component == 0 { 0.02 } else { -0.01 }],
                ForceInitialization::RetainState,
            )
            .unwrap()
        })
        .collect();
    ScheduledRenderer::from_coupled_modal_forces_and_impulses(
        voices,
        vec![
            ModalForceEvent {
                sample: 37,
                voice: 0,
                port: 0,
                force_n: 0.0,
            },
            ModalForceEvent {
                sample: 37,
                voice: 1,
                port: 0,
                force_n: 0.1,
            },
        ],
        render_inputs(),
        ForceRenderConfig {
            sample_rate_hz: 2000,
            max_block: 128,
            max_events: 2,
            max_controls: 2,
            max_projection_terms: 4,
        },
        PhysicalImpulseBudget {
            max_events: 5,
            max_projection_terms: 4,
            max_staged_modes,
        },
        render_links(),
        config(),
        &CancelGate::new_clock_free(),
    )
}

#[test]
fn physical_coupled_impulses_retain_held_forces_and_callback_boundary_ownership() {
    let gate = CancelGate::new_clock_free();
    let mut direct =
        CoupledModalSystem::new(render_models(), render_links(), config(), &gate).unwrap();
    let mut expected = Vec::new();
    let mut work = Vec::new();
    for sample in 0..128 {
        let mut physical = [0.0; 2];
        for input in render_inputs()
            .into_iter()
            .filter(|event| event.sample == sample)
        {
            physical[input.voice] += input.impulse_n_s;
        }
        if physical != [0.0; 2] {
            work.push(
                direct
                    .apply_modal_impulse(&[0.5 * physical[0], physical[1]], &gate)
                    .unwrap()
                    .kinetic_energy_change_j,
            );
        }
        let force = if sample < 37 {
            [0.01, -0.01]
        } else {
            [0.0, 0.1]
        };
        expected.push(direct.step(&force).unwrap().observer_pressure_pa);
    }
    assert!(expected.iter().any(|x| x.abs() > 0.01));
    for block in [1, 37, 128] {
        let mut render = performance(6).unwrap();
        let mut actual = [0.0; 128];
        for chunk in actual.chunks_mut(block) {
            render.block(chunk).unwrap();
        }
        assert_eq!(actual.as_slice(), expected, "callback length {block}");
        assert_eq!(
            states(render.coupled_modal_system(0).unwrap()),
            states(&direct)
        );
        assert_eq!(
            render
                .receipts()
                .iter()
                .map(|r| r.kinetic_energy_change_j)
                .collect::<Vec<_>>(),
            work
        );
        assert_eq!(
            render
                .receipts()
                .iter()
                .map(|r| r.sample)
                .collect::<Vec<_>>(),
            [0, 37, 127]
        );
    }
    // Four projected scalar entries fit, but three samples require SIX modes
    // of complete native-network staging. No component is silently skipped.
    assert!(matches!(performance(5), Err(RenderError::Sizing { .. })));
    let mut resumed = performance(6).unwrap();
    let mut actual = [0.0; 128];
    resumed.block(&mut actual[..37]).unwrap();
    assert_eq!(resumed.receipts().len(), 1);
    assert_eq!(resumed.renderer().pending_controls().len(), 2);
    let boundary = states(resumed.coupled_modal_system(0).unwrap());
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    let mut untouched = [123.0; 91];
    assert_eq!(
        resumed
            .render_under_gate(&cancelled, &mut untouched, 91, 1)
            .unwrap(),
        GatedRenderOutcome::Cancelled { blocks: 0 }
    );
    assert_eq!(untouched, [123.0; 91]);
    assert_eq!(states(resumed.coupled_modal_system(0).unwrap()), boundary);
    assert_eq!(resumed.receipts().len(), 1);
    resumed.block(&mut actual[37..]).unwrap();
    assert_eq!(actual.as_slice(), expected);
}

#[test]
fn same_time_failure_is_atomic_across_coupled_and_independent_voices() {
    let gate = CancelGate::new_clock_free();
    let system = CoupledModalSystem::new(render_models(), render_links(), config(), &gate).unwrap();
    let before = states(&system);
    let bounded = model(
        2000,
        0.0,
        0.0,
        0.0,
        ModalAcousticTimeBudget {
            maximum_abs_velocity_m_sqrt_kg_per_s: 0.1,
            ..ModalAcousticTimeBudget::audible_reference()
        },
    );
    let voices = vec![
        RenderVoice::CoupledModal(Box::new(
            CoupledModalVoice::new(system, vec![0.0; 2]).unwrap(),
        )),
        RenderVoice::ModalString(ModalStringVoice::new(bounded, vec![0.0]).unwrap()),
    ];
    let mut render = ScheduledRenderer::new(RenderContext::new(voices, 32), vec![], 0)
        .unwrap()
        .with_modal_impulses(
            vec![
                ModalImpulse {
                    sample: 0,
                    voice: 0,
                    mode: 0,
                    impulse_n_s_per_sqrt_kg: 0.25,
                },
                ModalImpulse {
                    sample: 0,
                    voice: 1,
                    mode: 0,
                    impulse_n_s_per_sqrt_kg: 0.25,
                },
            ],
            ImpulseBudget {
                max_events: 2,
                max_staged_modes: 3,
            },
        )
        .unwrap();
    let mut output = [123.0];
    assert!(matches!(
        render.block(&mut output),
        Err(RenderError::Modal(_))
    ));
    assert_eq!(output, [123.0]);
    assert_eq!(render.samples_rendered(), 0);
    assert_eq!(states(render.coupled_modal_system(0).unwrap()), before);
    assert_eq!(
        render.modal_states(1).unwrap()[0].velocity_m_sqrt_kg_per_s,
        0.0
    );
    assert!(render.receipts().is_empty());
    assert!(matches!(
        render.block(&mut output),
        Err(RenderError::Poisoned)
    ));
}

#[test]
fn independent_signed_work_retains_small_energy_between_opposite_modal_works() {
    let large = 2.0_f64.powi(28);
    for (small_before, small_impulse, expected_work) in [(0.0, 1.0, 0.5), (1.0, -1.0, -0.5)] {
        let mut source = ModalAcousticTimeModel::try_new_with_free_coordinates(
            2000,
            vec![
                ModalAcousticMode {
                    angular_frequency_rad_s: 0.0,
                    damping_ratio: 0.0,
                    pressure_per_modal_velocity: C64::ZERO,
                };
                3
            ],
            ModalAcousticTimeBudget {
                maximum_abs_displacement_m_sqrt_kg: 1e12,
                maximum_abs_velocity_m_sqrt_kg_per_s: 1e20,
                maximum_total_energy_j: 1e40,
                ..ModalAcousticTimeBudget::audible_reference()
            },
        )
        .unwrap();
        source
            .restore_states(&[0.0, small_before, large].map(|velocity_m_sqrt_kg_per_s| {
                ModalAcousticState {
                    displacement_m_sqrt_kg: 0.0,
                    velocity_m_sqrt_kg_per_s,
                }
            }))
            .unwrap();
        let voice = RenderVoice::ModalString(ModalStringVoice::new(source, vec![0.0; 3]).unwrap());
        let mut render = ScheduledRenderer::new(RenderContext::new(vec![voice], 1), vec![], 0)
            .unwrap()
            .with_modal_impulses(
                [large, small_impulse, -large]
                    .into_iter()
                    .enumerate()
                    .map(|(mode, impulse_n_s_per_sqrt_kg)| ModalImpulse {
                        sample: 0,
                        voice: 0,
                        mode,
                        impulse_n_s_per_sqrt_kg,
                    })
                    .collect(),
                ImpulseBudget {
                    max_events: 3,
                    max_staged_modes: 3,
                },
            )
            .unwrap();
        render.block(&mut [0.0]).unwrap();
        // Per-mode works in basis order are +2^55, +/-0.5, -2^55 J.
        // Subtracting endpoint total energies or a plain sum would return zero.
        assert_eq!(render.receipts()[0].kinetic_energy_change_j, expected_work);
        assert_eq!(
            render
                .modal_states(0)
                .unwrap()
                .iter()
                .map(|state| state.velocity_m_sqrt_kg_per_s)
                .collect::<Vec<_>>(),
            [large, small_before + small_impulse, 0.0]
        );
    }
}

#[test]
fn unrepresented_independent_impulse_refuses_before_any_same_time_publication() {
    let gate = CancelGate::new_clock_free();
    for include_coupled in [false, true] {
        let mut voices = Vec::new();
        let mut coupled_before = None;
        let mut inputs = Vec::new();
        if include_coupled {
            let system =
                CoupledModalSystem::new(render_models(), render_links(), config(), &gate).unwrap();
            coupled_before = Some(states(&system));
            voices.push(RenderVoice::CoupledModal(Box::new(
                CoupledModalVoice::new(system, vec![0.0; 2]).unwrap(),
            )));
            inputs.push(ModalImpulse {
                sample: 0,
                voice: 0,
                mode: 0,
                impulse_n_s_per_sqrt_kg: 0.25,
            });
        }
        let independent_slot = voices.len();
        let source = model(
            2000,
            0.0,
            0.0,
            1e16,
            ModalAcousticTimeBudget {
                maximum_abs_velocity_m_sqrt_kg_per_s: 1e20,
                maximum_total_energy_j: 1e40,
                ..ModalAcousticTimeBudget::audible_reference()
            },
        );
        let independent_before = source.states().to_vec();
        voices.push(RenderVoice::ModalString(
            ModalStringVoice::new(source, vec![0.0]).unwrap(),
        ));
        inputs.push(ModalImpulse {
            sample: 0,
            voice: independent_slot,
            mode: 0,
            impulse_n_s_per_sqrt_kg: 1.0,
        });
        let mut render = ScheduledRenderer::new(RenderContext::new(voices, 1), vec![], 0)
            .unwrap()
            .with_modal_impulses(
                inputs,
                ImpulseBudget {
                    max_events: 2,
                    max_staged_modes: 3,
                },
            )
            .unwrap();
        let mut output = [123.0];
        assert!(matches!(
            render.block(&mut output),
            Err(RenderError::Control { .. })
        ));
        assert_eq!(output, [123.0]);
        assert_eq!(render.samples_rendered(), 0);
        assert_eq!(
            render.modal_states(independent_slot).unwrap(),
            independent_before
        );
        if let Some(before) = coupled_before {
            assert_eq!(states(render.coupled_modal_system(0).unwrap()), before);
        }
        assert!(render.receipts().is_empty());
        assert!(matches!(
            render.block(&mut output),
            Err(RenderError::Poisoned)
        ));
    }
}
