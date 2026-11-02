use super::*;
use fs_couple::render::plate::impact::material::heating::SqueezeFilmLayerHeatRoute;
use fs_couple::render::plate::impact::material::heating::layered_film::{
    GapViscosityProfile, ViscosityLayer,
};
use fs_couple::render::plate::impact::material::heating::liquid_film::LiquidFilmLayerBinding;

fn profile(interface: f64, low: f64, high: f64) -> GapViscosityProfile {
    GapViscosityProfile::new(vec![
        ViscosityLayer {
            upper_fraction: interface,
            dynamic_viscosity_pa_s: low,
        },
        ViscosityLayer {
            upper_fraction: 1.0,
            dynamic_viscosity_pa_s: high,
        },
    ])
    .unwrap()
}

fn layer_bindings(card: &MaterialCard) -> Vec<LiquidFilmLayerBinding<'_>> {
    // Authoring order must not change physical layer identity.
    [(1, "sink"), (0, "fluid")]
        .into_iter()
        .map(|(layer, reservoir)| LiquidFilmLayerBinding {
            film: 0,
            channel: 0,
            layer,
            reservoir: reservoir.into(),
            law: law(card),
            maximum_relative_viscosity_change: 0.1,
        })
        .collect()
}

fn layered_inner<'m>(
    card: &MaterialCard,
    c: &'m EquilibriumEnthalpyPhaseCurve,
    interface: f64,
) -> HeatedImpactSystem<'m> {
    HeatedImpactSystem::new_with_layer_losses(
        mechanics(law(card).reference_state(), 0.0)
            .prepare_analytic()
            .unwrap(),
        stores(c),
        vec![],
        vec![],
        vec![],
        vec![SqueezeFilmLayerHeatRoute {
            film: 0,
            channel: 0,
            profile: profile(interface, 0.01, 0.01),
            layer_destinations: vec!["fluid".into(), "sink".into()],
        }],
        vec![],
        1e-13,
    )
    .unwrap()
}

fn layered<'m>(
    card: &'m MaterialCard,
    c: &'m EquilibriumEnthalpyPhaseCurve,
    interface: f64,
) -> LiquidFilmImpactSystem<'m> {
    LiquidFilmImpactSystem::new_with_layers(
        layered_inner(card, c, interface),
        vec![],
        vec![],
        layer_bindings(card),
    )
    .unwrap()
}

/// Independent exact polynomial integrals, well conditioned for these fixtures.
fn moments(interface: f64, low: f64, high: f64) -> (f64, [f64; 2]) {
    let a = interface / low + (1.0 - interface) / high;
    let b = 0.5 * (interface.powi(2) / low + (1.0 - interface.powi(2)) / high);
    let z = b / a;
    let d0 = (interface.powi(3) / 3.0 - z * interface.powi(2) + z * z * interface) / low;
    let d1 = ((1.0 - interface.powi(3)) / 3.0 - z * (1.0 - interface.powi(2))
        + z * z * (1.0 - interface))
        / high;
    (1.0 / (12.0 * (d0 + d1)), [d0 / (d0 + d1), d1 / (d0 + d1)])
}

#[test]
fn source_layers_change_actual_motion_and_partition_only_the_completed_profiles_work() {
    // G1/G3: independent Poiseuille flux and through-gap work integrals at each
    // OLD accepted profile, followed by actual native material/thermal evolution.
    let card = card(0.2, 0.001);
    let c = curve(card.content_hash(), 0.2, 1.0);
    let law = law(&card);
    for analytic in [false, true] {
        let mut live = layered(&card, &c, 0.25);
        let mut frozen = layered_inner(&card, &c, 0.25);
        live.set_analytic_newton(analytic);
        frozen.set_analytic_newton(analytic);
        let initial = live.system().mechanics().stored_energy_j();
        assert!(live.sources().is_empty() && live.channel_sources().is_empty());
        assert_eq!(
            live.layer_bindings()
                .iter()
                .map(|b| b.layer)
                .collect::<Vec<_>>(),
            [0, 1]
        );
        for _ in 0..64 {
            let old = live.system().mechanics().state().to_vec();
            let old_mu = live
                .layer_sources()
                .iter()
                .map(ResolvedLiquidState::dynamic_viscosity_pa_s)
                .collect::<Vec<_>>();
            let (equivalent, fractions) = moments(0.25, old_mu[0], old_mu[1]);
            let received = live
                .system()
                .reservoirs()
                .iter()
                .map(CaloricReservoir::received_heat_j)
                .collect::<Vec<_>>();
            let frame = live.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
            frozen
                .step(&[0.0; 2], &CancelGate::new_clock_free())
                .unwrap();
            let next = live.system().mechanics().state();
            let closure =
                0.5 * ((old[0] + next[0]) / 0.2_f64.sqrt() - (old[2] + next[2]) / 0.4_f64.sqrt());
            let speed =
                0.5 * ((old[1] + next[1]) / 0.2_f64.sqrt() - (old[3] + next[3]) / 0.4_f64.sqrt());
            let g = WIDTH * (GAP - closure).powi(3) / (12.0 * equivalent * LENGTH);
            let heat = DT * (AREA * speed).powi(2) / g;
            near(frame.thermal.impact.deposited_heat_j, heat, 2e-20);
            let profile = live.system().layer_heat_profile(0, 0).unwrap();
            assert_eq!(profile.layers()[0].upper_fraction, 0.25);
            assert_eq!(profile.layers()[1].upper_fraction, 1.0);
            for layer in 0..2 {
                let store = &live.system().reservoirs()[layer];
                near(
                    store.received_heat_j() - received[layer],
                    heat * fractions[layer],
                    2e-20,
                );
                let source = law
                    .at_temperature(Temperature::new(store.state().temperature_k()))
                    .unwrap();
                assert_eq!(&live.layer_sources()[layer], &source);
                assert_eq!(
                    profile.layers()[layer].dynamic_viscosity_pa_s,
                    source.dynamic_viscosity_pa_s()
                );
            }
            let (next_mu, next_fraction) = moments(
                0.25,
                live.layer_sources()[0].dynamic_viscosity_pa_s(),
                live.layer_sources()[1].dynamic_viscosity_pa_s(),
            );
            near(
                profile.hydraulic_equivalent_viscosity_pa_s(),
                next_mu,
                2e-16,
            );
            for (actual, expected) in profile
                .layer_pressure_work_fractions()
                .iter()
                .zip(next_fraction)
            {
                near(*actual, expected, 2e-14);
            }
            assert_eq!(
                live.system()
                    .mechanics()
                    .squeeze_film(0)
                    .unwrap()
                    .channel_viscosities_pa_s(),
                &[profile.hydraulic_equivalent_viscosity_pa_s()]
            );
        }
        assert_ne!(
            live.layer_sources()[0].dynamic_viscosity_pa_s(),
            live.layer_sources()[1].dynamic_viscosity_pa_s()
        );
        assert_ne!(
            live.system()
                .layer_heat_profile(0, 0)
                .unwrap()
                .layer_pressure_work_fractions(),
            frozen
                .layer_heat_profile(0, 0)
                .unwrap()
                .layer_pressure_work_fractions()
        );
        assert_ne!(
            live.system().mechanics().state(),
            frozen.mechanics().state()
        );
        near(
            live.system().mechanics().stored_energy_j()
                + live
                    .system()
                    .reservoirs()
                    .iter()
                    .map(CaloricReservoir::storage_change_j)
                    .sum::<f64>(),
            initial,
            1e-13,
        );
    }
}

#[test]
fn one_layer_recovers_uniform_feedback_and_rebinding_preserves_the_accepted_clock() {
    let card = card(0.2, 0.001);
    let c = curve(card.content_hash(), 0.2, 1.0);
    for analytic in [false, true] {
        let inner = HeatedImpactSystem::new_with_layer_losses(
            mechanics(law(&card).reference_state(), 0.0)
                .prepare_analytic()
                .unwrap(),
            stores(&c),
            vec![],
            vec![],
            vec![],
            vec![SqueezeFilmLayerHeatRoute {
                film: 0,
                channel: 0,
                profile: GapViscosityProfile::new(vec![ViscosityLayer {
                    upper_fraction: 1.0,
                    dynamic_viscosity_pa_s: 0.01,
                }])
                .unwrap(),
                layer_destinations: vec!["fluid".into()],
            }],
            vec![],
            1e-13,
        )
        .unwrap();
        let mut live = LiquidFilmImpactSystem::new_with_layers(
            inner,
            vec![],
            vec![],
            vec![LiquidFilmLayerBinding {
                film: 0,
                channel: 0,
                layer: 0,
                reservoir: "fluid".into(),
                law: law(&card),
                maximum_relative_viscosity_change: 0.1,
            }],
        )
        .unwrap();
        let mut uniform = coupled(&card, &c);
        live.set_analytic_newton(analytic);
        uniform.set_analytic_newton(analytic);
        for _ in 0..16 {
            live.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
            uniform
                .step(&[0.0; 2], &CancelGate::new_clock_free())
                .unwrap();
            assert_eq!(
                live.system().mechanics().state(),
                uniform.system().mechanics().state()
            );
            assert_eq!(live.system().reservoirs(), uniform.system().reservoirs());
            assert_eq!(live.layer_sources(), uniform.sources());
        }
        let (inner, whole, channels, layers) = live.into_layered_parts();
        let samples = inner.mechanics().samples();
        let mut rebound =
            LiquidFilmImpactSystem::new_with_layers(inner, whole, channels, layers).unwrap();
        assert_eq!(rebound.system().mechanics().samples(), samples);
        rebound
            .step(&[0.0; 2], &CancelGate::new_clock_free())
            .unwrap();
        uniform
            .step(&[0.0; 2], &CancelGate::new_clock_free())
            .unwrap();
        assert_eq!(
            rebound.system().mechanics().state(),
            uniform.system().mechanics().state()
        );
        assert_eq!(rebound.system().reservoirs(), uniform.system().reservoirs());
        assert_eq!(rebound.layer_sources(), uniform.sources());
    }
}

#[test]
fn layer_sources_follow_actual_network_temperatures_and_conserve_total_energy() {
    let card = card(0.2, 0.001);
    let c = curve(card.content_hash(), 0.2, 1.0);
    let law = law(&card);
    let mut live = layered(&card, &c, 0.25);
    let mut thermal = network(&c);
    let initial = live.system().mechanics().stored_energy_j();
    for _ in 0..32 {
        let frame = live
            .step_with_thermal_network(
                &mut thermal,
                &mut (),
                &local_network_bindings(),
                StepConfig {
                    max_dt_s: DT / 4.0,
                    energy_tolerance_j: 1e-13,
                },
                &[0.0; 2],
                &CancelGate::new_clock_free(),
                |_, _| Ok(None),
            )
            .unwrap();
        assert!(frame.network.substeps >= 4);
        near(frame.external_heat_j, 0.0, 1e-19);
        assert_eq!(
            thermal.time_s(),
            live.system().mechanics().samples() as f64 * DT
        );
        for (layer, binding) in local_network_bindings().iter().enumerate() {
            let state = c
                .state_at_specific_enthalpy(thermal.cells()[&binding.cell].specific_h_j_kg)
                .unwrap();
            assert_eq!(live.system().reservoirs()[layer].state(), state);
            assert_eq!(
                live.layer_sources()[layer],
                law.at_temperature(Temperature::new(state.temperature_k()))
                    .unwrap()
            );
        }
    }
    let stores = live.system().reservoirs();
    assert!(stores[0].storage_change_j() > stores[0].received_heat_j());
    assert!(stores[1].storage_change_j() < stores[1].received_heat_j());
    near(
        live.system().mechanics().stored_energy_j()
            + stores
                .iter()
                .map(CaloricReservoir::storage_change_j)
                .sum::<f64>(),
        initial,
        1e-13,
    );
}

#[test]
fn late_layer_refusals_and_cancellation_preserve_profiles_sources_and_network() {
    // G4: the second, very thin layer can change far beyond its admitted
    // temporal allowance while its effect on hydraulic viscosity is tiny.
    let card = card(0.2, 0.001);
    let c = curve(card.content_hash(), 0.2, 1.0);
    let interface = 0.9999;
    let hypothetical = profile(interface, 0.01, 0.0082).hydraulic_equivalent_viscosity_pa_s();
    assert!((hypothetical - 0.01).abs() / 0.01 < 1e-3);
    for (rise_k, cancel) in [(0.002, false), (0.0002, false), (0.0, true)] {
        let mut live = layered(&card, &c, interface);
        let mut thermal = network(&c);
        let gate = CancelGate::new_clock_free();
        let mut observer = LocalHeater {
            power_w: MASS * 0.2 * rise_k / DT,
            steps: 0,
            cancel: cancel.then_some(&gate),
        };
        let motion = live.system().mechanics().state().to_vec();
        let stores = live.system().reservoirs().to_vec();
        let sources = live.layer_sources().to_vec();
        let profile = live.system().layer_heat_profile(0, 0).unwrap().clone();
        let cells = thermal.cells().clone();
        let result = live.step_with_thermal_network(
            &mut thermal,
            &mut observer,
            &local_network_bindings(),
            StepConfig {
                max_dt_s: DT / 4.0,
                energy_tolerance_j: 1e-13,
            },
            &[0.0; 2],
            &gate,
            |_, _| Ok(None),
        );
        use fs_couple::render::plate::impact::ImpactError;
        let error = result.unwrap_err();
        if cancel {
            assert!(matches!(error, ImpactError::Cancelled));
        } else if rise_k > 0.001 {
            assert!(matches!(error, ImpactError::Owner(_)), "{error}");
        } else {
            assert!(
                matches!(error, ImpactError::Invalid(ref reason)
                if reason.contains("viscosity changed too far")),
                "{error}"
            );
        }
        assert_eq!(live.system().mechanics().state(), motion);
        assert_eq!(live.system().mechanics().samples(), 0);
        assert_eq!(live.system().reservoirs(), stores);
        assert_eq!(live.layer_sources(), sources);
        assert_eq!(live.system().layer_heat_profile(0, 0), Some(&profile));
        assert_eq!(
            live.system()
                .mechanics()
                .squeeze_film(0)
                .unwrap()
                .channel_viscosities_pa_s(),
            &[0.01]
        );
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(thermal.cells(), &cells);
        assert_eq!(observer.steps, 0);
        observer.power_w = 0.0;
        observer.cancel = None;
        let mut fresh = layered(&card, &c, interface);
        let mut fresh_thermal = network(&c);
        live.step_with_thermal_network(
            &mut thermal,
            &mut observer,
            &local_network_bindings(),
            StepConfig {
                max_dt_s: DT / 4.0,
                energy_tolerance_j: 1e-13,
            },
            &[0.0; 2],
            &CancelGate::new_clock_free(),
            |_, _| Ok(None),
        )
        .unwrap();
        fresh
            .step_with_thermal_network(
                &mut fresh_thermal,
                &mut (),
                &local_network_bindings(),
                StepConfig {
                    max_dt_s: DT / 4.0,
                    energy_tolerance_j: 1e-13,
                },
                &[0.0; 2],
                &CancelGate::new_clock_free(),
                |_, _| Ok(None),
            )
            .unwrap();
        assert_eq!(
            live.system().mechanics().state(),
            fresh.system().mechanics().state()
        );
        assert_eq!(live.system().reservoirs(), fresh.system().reservoirs());
        assert_eq!(live.layer_sources(), fresh.layer_sources());
        assert_eq!(
            live.system().layer_heat_profile(0, 0),
            fresh.system().layer_heat_profile(0, 0)
        );
        assert_eq!(thermal.cells(), fresh_thermal.cells());
    }
}

#[test]
fn layer_bindings_require_complete_unique_source_consistent_coverage() {
    let card = card(0.2, 0.001);
    let c = curve(card.content_hash(), 0.2, 1.0);
    let original = layer_bindings(&card);
    let mut outside = original.clone();
    outside[0].layer = 2;
    let mut wrong_store = original.clone();
    wrong_store[0].reservoir = "fluid".into();
    let mut unknown_channel = original.clone();
    unknown_channel[0].channel = 1;
    let mut bad_limit = original.clone();
    bad_limit[0].maximum_relative_viscosity_change = f64::NAN;
    for layers in [
        vec![],
        vec![original[0].clone()],
        vec![original[0].clone(), original[0].clone()],
        outside,
        wrong_store,
        unknown_channel,
        bad_limit,
    ] {
        assert!(
            LiquidFilmImpactSystem::new_with_layers(
                layered_inner(&card, &c, 0.25),
                vec![],
                vec![],
                layers,
            )
            .is_err()
        );
    }
    assert!(
        LiquidFilmImpactSystem::new_with_layers(
            layered_inner(&card, &c, 0.25),
            vec![LiquidFilmBinding {
                film: 0,
                reservoir: "fluid".into(),
                law: law(&card),
                maximum_relative_viscosity_change: 0.1,
            }],
            vec![],
            original.clone(),
        )
        .is_err()
    );
    assert!(
        LiquidFilmImpactSystem::new_with_layers(
            heated(law(&card).reference_state(), &c, 0.0),
            vec![],
            vec![],
            original,
        )
        .is_err()
    );
    let other = super::card(0.3, 0.001);
    let mut mismatch = layer_bindings(&card);
    mismatch[0].law = law(&other);
    assert!(
        LiquidFilmImpactSystem::new_with_layers(
            layered_inner(&card, &c, 0.25),
            vec![],
            vec![],
            mismatch,
        )
        .is_err()
    );
}

#[test]
fn mixed_uniform_channel_and_layer_islands_keep_distinct_source_segments() {
    let card = card(0.2, 0.001);
    let c = curve(card.content_hash(), 0.2, 1.0);
    let law = law(&card);
    // Three independent films share the same coordinates; use a smaller
    // initial speed to stay inside the original mechanical/pressure envelope.
    let (a, _) = ImpactBody::free_mass(0.2, 0.0, 0.0001).unwrap();
    let (b, _) = ImpactBody::free_mass(0.4, 0.0, 0.0).unwrap();
    let mut motion = ImpactSystem::new(
        vec![a, b],
        vec![],
        vec![],
        vec![],
        ImpactConfig {
            dt_s: DT,
            max_steps: 512,
            maximum_energy_j: 1.0,
            energy_absolute_tolerance_j: 1e-14,
            energy_relative_tolerance: 1e-7,
            maximum_generalized_force: 1.0,
        },
    )
    .unwrap();
    for _ in 0..3 {
        motion = film(law.reference_state())
            .unwrap()
            .attach(motion)
            .unwrap()
            .0;
    }
    let inner = HeatedImpactSystem::new_with_layer_losses(
        motion.prepare_analytic().unwrap(),
        stores(&c),
        vec![],
        vec![route()],
        vec![SqueezeFilmChannelHeatRoute {
            film: 1,
            channel: 0,
            to_a: "fluid".into(),
            to_b: "fluid".into(),
            fraction_to_a: 1.0,
        }],
        vec![SqueezeFilmLayerHeatRoute {
            film: 2,
            channel: 0,
            profile: profile(0.25, 0.01, 0.01),
            layer_destinations: vec!["sink".into(), "fluid".into()],
        }],
        vec![],
        1e-13,
    )
    .unwrap();
    let mut live = LiquidFilmImpactSystem::new_with_layers(
        inner,
        vec![LiquidFilmBinding {
            film: 0,
            reservoir: "fluid".into(),
            law: law.clone(),
            maximum_relative_viscosity_change: 0.1,
        }],
        vec![LiquidFilmChannelBinding {
            film: 1,
            channel: 0,
            reservoir: "fluid".into(),
            law: law.clone(),
            maximum_relative_viscosity_change: 0.1,
        }],
        [(1, "fluid"), (0, "sink")]
            .into_iter()
            .map(|(layer, reservoir)| LiquidFilmLayerBinding {
                film: 2,
                channel: 0,
                layer,
                reservoir: reservoir.into(),
                law: law.clone(),
                maximum_relative_viscosity_change: 0.1,
            })
            .collect(),
    )
    .unwrap();
    for _ in 0..8 {
        live.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
        assert_eq!(live.sources().len(), 1);
        assert_eq!(live.channel_sources().len(), 1);
        assert_eq!(live.layer_sources().len(), 2);
        assert_eq!(live.sources(), live.channel_sources());
        assert_eq!(&live.layer_sources()[1], &live.sources()[0]);
        assert_ne!(&live.layer_sources()[0], &live.sources()[0]);
        let lower = live.layer_sources()[0].dynamic_viscosity_pa_s();
        let upper = live.layer_sources()[1].dynamic_viscosity_pa_s();
        let expected = moments(0.25, lower, upper).0;
        near(
            live.system()
                .mechanics()
                .squeeze_film(2)
                .unwrap()
                .channel_viscosities_pa_s()[0],
            expected,
            2e-16,
        );
        for film in [0, 1] {
            assert_eq!(
                live.system()
                    .mechanics()
                    .squeeze_film(film)
                    .unwrap()
                    .channel_viscosities_pa_s(),
                &[upper]
            );
        }
    }
}
