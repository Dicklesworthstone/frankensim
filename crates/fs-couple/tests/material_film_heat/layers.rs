use super::*;
use fs_couple::render::plate::impact::material::heating::{
    SqueezeFilmChannelHeatRoute, SqueezeFilmLayerHeatRoute,
};
use fs_couple::render::plate::impact::material::heating::layered_film::{
    GapViscosityProfile, ViscosityLayer,
};

fn profile(low: f64, high: f64) -> GapViscosityProfile {
    GapViscosityProfile::new(vec![
        ViscosityLayer { upper_fraction: 0.5, dynamic_viscosity_pa_s: low },
        ViscosityLayer { upper_fraction: 1.0, dynamic_viscosity_pa_s: high },
    ]).unwrap()
}

fn layer_route(profile: GapViscosityProfile) -> SqueezeFilmLayerHeatRoute {
    SqueezeFilmLayerHeatRoute {
        film: 0, channel: 0, profile,
        layer_destinations: vec!["fluid".into(), "sink".into()],
    }
}

fn layered<'m>(
    source: &ResolvedLiquidState, c: &'m EquilibriumEnthalpyPhaseCurve,
    profile: GapViscosityProfile,
) -> HeatedImpactSystem<'m> {
    HeatedImpactSystem::new_with_layer_losses(
        mechanics(source, 0.0).prepare_analytic().unwrap(), stores(c),
        vec![], vec![], vec![], vec![layer_route(profile)], vec![], 1e-13,
    ).unwrap()
}

#[test]
fn retained_layers_receive_independent_poiseuille_work_without_changing_motion() {
    let profile = profile(0.005, 0.02);
    let source = liquid(profile.hydraulic_equivalent_viscosity_pa_s(), 0.2);
    let c = curve(source.material().card_identity(), 0.2, 1.0);
    // Independently integrate (z - 0.35)^2 / mu over the two half gaps.
    let fractions = [148.0 / 365.0, 217.0 / 365.0];
    for analytic in [false, true] {
        let mut system = layered(&source, &c, profile.clone());
        let mut unsplit = heated(&source, &c, 0.0);
        system.set_analytic_newton(analytic);
        unsplit.set_analytic_newton(analytic);
        assert_eq!(system.layer_heat_profile(0, 0), Some(&profile));
        let initial = system.mechanics().stored_energy_j();
        let mut heat = 0.0;
        for _ in 0..32 {
            let old = system.mechanics().state().to_vec();
            let frame = system.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
            let other = unsplit.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
            let next = system.mechanics().state();
            assert_eq!(next, unsplit.mechanics().state());
            assert_eq!(frame.deposited_heat_j, other.deposited_heat_j);
            let closure = 0.5 * ((old[0] + next[0]) / 0.2_f64.sqrt()
                - (old[2] + next[2]) / 0.4_f64.sqrt());
            let speed = 0.5 * ((old[1] + next[1]) / 0.2_f64.sqrt()
                - (old[3] + next[3]) / 0.4_f64.sqrt());
            let conductance = WIDTH * (GAP - closure).powi(3)
                / (12.0 * source.dynamic_viscosity_pa_s() * LENGTH);
            near(frame.deposited_heat_j, DT * (AREA * speed).powi(2) / conductance, 2e-20);
            near(frame.remaining_mechanical_loss_j, 0.0, 2e-20);
            heat += frame.deposited_heat_j;
        }
        for (store, fraction) in system.reservoirs().iter().zip(fractions) {
            near(store.received_heat_j(), heat * fraction, 2e-20);
            near(store.storage_change_j(), store.received_heat_j(), 2e-20);
        }
        near(system.mechanics().stored_energy_j()
            + system.reservoirs().iter().map(CaloricReservoir::storage_change_j).sum::<f64>(),
            initial, 1e-13);
    }
}

#[test]
fn one_layer_recovers_the_original_heat_route_and_tiny_layers_keep_positive_work() {
    let uniform = GapViscosityProfile::new(vec![ViscosityLayer {
        upper_fraction: 1.0, dynamic_viscosity_pa_s: 0.01,
    }]).unwrap();
    let source = liquid(0.01, 0.2);
    let c = curve(source.material().card_identity(), 0.2, 1.0);
    for analytic in [false, true] {
        let mut system = HeatedImpactSystem::new_with_layer_losses(
            mechanics(&source, 0.0).prepare_analytic().unwrap(), stores(&c),
            vec![], vec![], vec![], vec![SqueezeFilmLayerHeatRoute {
                film: 0, channel: 0, profile: uniform.clone(),
                layer_destinations: vec!["fluid".into()],
            }], vec![], 1e-13,
        ).unwrap();
        let mut plain = heated(&source, &c, 0.0);
        system.set_analytic_newton(analytic);
        plain.set_analytic_newton(analytic);
        for _ in 0..8 {
            let frame = system.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
            let other = plain.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
            assert_eq!(frame.deposited_heat_j, other.deposited_heat_j);
            assert_eq!(system.mechanics().state(), plain.mechanics().state());
            assert_eq!(system.reservoirs(), plain.reservoirs());
        }
    }
    let disparate = profile(0.001, 1e20);
    let source = liquid(disparate.hydraulic_equivalent_viscosity_pa_s(), 0.2);
    let c = curve(source.material().card_identity(), 0.2, 1.0);
    let expected = disparate.layer_pressure_work_fractions()[1];
    assert!(expected > 0.0 && expected < 1e-18);
    let mut system = layered(&source, &c, disparate);
    let frame = system.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
    let tiny = system.reservoirs()[1].received_heat_j();
    assert!(tiny > 0.0);
    near(tiny / frame.deposited_heat_j, expected, expected * 1e-14);
    assert!(system.reservoirs()[1].state().specific_enthalpy_j_kg() > 0.0);
}

#[test]
fn layer_routes_refuse_ambiguous_coverage_and_failed_steps_publish_nothing() {
    let profile = profile(0.005, 0.02);
    let source = liquid(profile.hydraulic_equivalent_viscosity_pa_s(), 0.2);
    let c = curve(source.material().card_identity(), 0.2, 1.0);
    let good = layer_route(profile.clone());
    let construct = |routes, whole, channels| HeatedImpactSystem::new_with_layer_losses(
        mechanics(&source, 0.0).prepare_analytic().unwrap(), stores(&c),
        vec![], whole, channels, routes, vec![], 1e-13,
    );
    let mut short = good.clone(); short.layer_destinations.pop();
    let mut missing = good.clone(); missing.layer_destinations[1] = "absent".into();
    let mut bad_channel = good.clone(); bad_channel.channel = 1;
    let mut wrong = good.clone(); wrong.profile = self::profile(0.01, 0.02);
    for routes in [vec![short], vec![missing], vec![bad_channel], vec![wrong],
        vec![good.clone(), good.clone()]] {
        assert!(construct(routes, vec![], vec![]).is_err());
    }
    assert!(construct(vec![good.clone()], vec![route()], vec![]).is_err());
    assert!(construct(vec![good], vec![], vec![SqueezeFilmChannelHeatRoute {
        film: 0, channel: 0, to_a: "fluid".into(), to_b: "sink".into(), fraction_to_a: 0.5,
    }]).is_err());

    let narrow = curve(source.material().card_identity(), 0.2, 1e-12);
    let mut rejected = layered(&source, &narrow, profile.clone());
    let old_motion = rejected.mechanics().state().to_vec();
    let old_stores = rejected.reservoirs().to_vec();
    assert!(rejected.step(&[0.0; 2], &CancelGate::new_clock_free()).is_err());
    assert_eq!(rejected.mechanics().state(), old_motion);
    assert_eq!(rejected.reservoirs(), old_stores);
    assert_eq!(rejected.mechanics().samples(), 0);

    let mut system = layered(&source, &c, profile.clone());
    let old_stores = system.reservoirs().to_vec();
    let stopped = CancelGate::new_clock_free(); stopped.request();
    assert!(system.step(&[0.0; 2], &stopped).is_err());
    assert_eq!(system.mechanics().samples(), 0);
    assert_eq!(system.reservoirs(), old_stores);
    let mut clean = layered(&source, &c, profile);
    system.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
    clean.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
    assert_eq!(system.mechanics().state(), clean.mechanics().state());
    assert_eq!(system.reservoirs(), clean.reservoirs());
}

#[test]
fn layer_work_enters_the_original_thermal_network_once() {
    let profile = profile(0.005, 0.02);
    let source = liquid(profile.hydraulic_equivalent_viscosity_pa_s(), 0.2);
    let c = curve(source.material().card_identity(), 0.2, 1.0);
    let mut system = layered(&source, &c, profile);
    let mut thermal = network(&c);
    let initial = system.mechanics().stored_energy_j();
    let bindings = [
        ReservoirCell { reservoir: "fluid".into(), cell: CellId(1) },
        ReservoirCell { reservoir: "sink".into(), cell: CellId(2) },
    ];
    let mut heat = 0.0;
    for _ in 0..16 {
        let frame = system.step_with_thermal_network(
            &mut thermal, &mut (), &bindings,
            StepConfig { max_dt_s: DT/4.0, energy_tolerance_j: 1e-13 },
            &[0.0; 2], &CancelGate::new_clock_free(), |_, _| Ok(None),
        ).unwrap();
        heat += frame.impact.thermal.impact.deposited_heat_j;
        near(frame.external_heat_j, 0.0, 1e-19);
        for (store, binding) in system.reservoirs().iter().zip(&bindings) {
            assert_eq!(store.state().specific_enthalpy_j_kg(),
                thermal.cells()[&binding.cell].specific_h_j_kg);
        }
        assert!(thermal.cells().values().all(|cell| cell.boundary.source_w == 0.0));
    }
    near(system.reservoirs().iter().map(CaloricReservoir::received_heat_j).sum(),
        heat, 2e-20);
    near(system.mechanics().stored_energy_j()
        + system.reservoirs().iter().map(CaloricReservoir::storage_change_j).sum::<f64>(),
        initial, 1e-13);
}
