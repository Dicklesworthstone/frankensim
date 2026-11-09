//! Fixed local emission survives scan-to-SDIRK continuation and stage sampling.
use super::*;
use fs_conduction::PointwiseAmbientRadiation;
use fs_conduction::transient::enthalpy::heterogeneous::EnthalpyRadiation;
use fs_conduction::transient::enthalpy::sdirk::{
    EnthalpySdirk2Input, EnthalpySdirk2Stage, EnthalpySdirk2Step, SDIRK2_GAMMA as G,
};

fn scan<'h, 'm, 'c>(cx: &Cx<'_>, history: &'h ThermalBoundaryHistory<'m>,
    material: &'h ConductivityModel, curves: &'c [EquilibriumEnthalpyPhaseCurve; 2],
    card: &'c InterfaceSystemCard, source: MovingHeatSource,
    patches: Vec<RadiationRegionHistory>) -> ScanningEnthalpyRun<'h, 'm, 'c>
{
    let initial = history.at(cx, 0.).unwrap();
    let contacts = TemperatureContacts::new(cx, history.mesh(), &initial, vec![],
        vec![law(history.mesh(), card)], 1).unwrap();
    ScanningEnthalpyRun::new(cx, base(cx, history, material, curves, 3),
        vec![source], 1, None).unwrap()
        .with_contacts_newton(cx, contacts, coupled_newton()).unwrap()
        .with_pointwise_ambient_radiation(cx, patches).unwrap()
}

// Exact simplex moment of the P1 trace: integral(T_h^4) = A/15 *
// sum_{a+b+c=4} T0^a T1^b T2^c. This does not call production quadrature
// or replace the nonuniform trace with its mean temperature.
fn radiation_w(temperature: &[f64], ambient: f64) -> f64 {
    let t = [temperature[0], temperature[1], temperature[3]];
    let mut moment = 0.;
    for a in 0..=4 {
        for b in 0..=4 - a {
            moment += t[0].powi(a) * t[1].powi(b) * t[2].powi(4 - a - b);
        }
    }
    0.2 * STEFAN_BOLTZMANN_W_M2_K4 * 0.5 * (moment / 15. - ambient.powi(4))
}

fn admit_stage(stage: &EnthalpySdirk2Stage, bath: f64, watts: f64,
    card: &InterfaceSystemCard)
{
    let report = stage.pointwise_radiation.as_ref().unwrap();
    assert_eq!(report.patches.len(), 1);
    close(report.patches[0].patch.ambient_temperature_k(), bath, 1e-12);
    assert!(report.patches[0].material_states.is_empty());
    close(report.outward_heat_w, radiation_w(&stage.temperature, bath), 2e-9);
    close(stage.radiation_out_w, report.outward_heat_w, 1e-12);
    close(stage.source_w, watts, 2e-11);
    assert!(stage.pointwise_enclosure_radiation.is_none());
    assert_eq!(stage.physical_contacts.as_ref().unwrap().states()[0]
        .resistance.card_identity(), card.content_hash());
    assert!(stage.stage_residual_norm_j <= coupled_newton().residual_tolerance_j);
    close(stage.stage_energy_residual_j, 0., coupled_newton().step.energy_tolerance_j);
}

fn admit_interval(step: &EnthalpySdirk2Step, old: &[f64], baths: [f64; 2]) {
    let radiation = step.duration_s * ((1. - G)
        * radiation_w(&step.first_stage.temperature, baths[0])
        + G * radiation_w(&step.second_stage.temperature, baths[1]));
    close(step.radiation_out_j, radiation, 3e-9);
    close(step.stored_energy_change_j,
        stored(old, &step.second_stage.specific_enthalpy_j_kg), 2e-9);
    close(step.stored_energy_change_j + radiation - step.source_j, 0., 2e-8);
    close(step.energy_residual_j, 0., coupled_newton().step.energy_tolerance_j);
}

#[test]
fn g1_pointwise_scan_handoff_preserves_latent_state_and_matches_both_native_sdirk_stages() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh(); let history = control_history(cx, &mesh); let card = card();
        let curves = [chart(false), chart(true)];
        let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(280., 1.6), (380., 3.6)]).unwrap());
        let source = heater(TimeSignal::Ramp { t_start: 0., t_end: 0.2,
            from: QtyAny::new(-0.3, LENGTH), to: QtyAny::new(0.3, LENGTH) },
            TimeSignal::Ramp { t_start: 0., t_end: 0.2,
                from: QtyAny::new(3., COMPONENT_POWER_DIMS),
                to: QtyAny::new(5., COMPONENT_POWER_DIMS) });
        let mut patches = radiation().patches;
        patches[0].ambient_temperature = TimeSignal::Ramp { t_start: 0., t_end: 0.2,
            from: QtyAny::new(300., TEMPERATURE_DIMS), to: QtyAny::new(320., TEMPERATURE_DIMS) };
        let mut evolving = scan(cx, &history, &material, &curves, &card, source.clone(), patches.clone());
        let accepted = evolving.advance_to(cx, 0.1, coupled_newton().step).unwrap();
        let old = accepted.solution().specific_enthalpy_j_kg.clone();
        let temperature = accepted.solution().temperature.clone();
        let fractions = accepted.solution().liquid_mass_fraction.clone();
        assert!(fractions[..3].iter().any(|&f| f > 0. && f < 1.));
        let mut continued = evolving.into_sdirk2().unwrap();
        assert_eq!(continued.conduction().specific_enthalpy_j_kg(), old.as_slice());
        assert_eq!(continued.conduction().temperature(), temperature.as_slice());
        assert_eq!(continued.conduction().liquid_mass_fraction(), fractions.as_slice());
        assert_eq!(continued.conduction().time_s(), 0.1);
        assert_eq!(continued.conduction().accepted_steps(), 1);
        assert_eq!(continued.moving_sources(), &[source.clone()]);
        assert_eq!(continued.pointwise_ambient_patches(), Some(patches.as_slice()));
        assert!(continued.material_ambient_patches().is_none());
        assert!(continued.material_enclosure_patches().is_none());

        let end = 0.15;
        let dt = end - continued.conduction().time_s();
        let times = [G.mul_add(dt, 0.1), end];
        let boundaries = times.map(|time| history.at(cx, time).unwrap());
        let sources = times.map(|time| {
            let sample = source.at(cx, &mesh, time).unwrap();
            ScalarField::Nodal((0..8).map(|v| 2. + sample.source.at(v)).collect())
        });
        let sampled = times.map(|time| vec![patches[0].at(time).unwrap()]);
        let first = PointwiseAmbientRadiation::new(cx, &mesh, &boundaries[0], &sampled[0]).unwrap();
        let second = PointwiseAmbientRadiation::new(cx, &mesh, &boundaries[1], &sampled[1]).unwrap();
        let owners = [&first, &second];
        let inputs = [0, 1].map(|i| EnthalpySdirk2Input {
            problem: ConductionProblem { mesh: &mesh, boundary: &boundaries[i], material: &material,
                element_materials: None, source: &sources[i] },
            interfaces: None, radiation: EnthalpyRadiation::PointwiseAmbient(owners[i]),
            contacts: Some((continued.contacts().unwrap(), coupled_newton())),
        });
        let native = continued.conduction().storage().advance_sdirk2(
            cx, inputs, &old, dt, coupled_newton().step).unwrap();
        let step = continued.advance_to(cx, end, coupled_newton().step).unwrap();
        for (actual, reference) in [(&step.first_stage, &native.first_stage),
            (&step.second_stage, &native.second_stage)]
        {
            assert_eq!(actual.specific_enthalpy_j_kg, reference.specific_enthalpy_j_kg);
            assert_eq!(actual.temperature, reference.temperature);
            assert_eq!(actual.liquid_mass_fraction, reference.liquid_mass_fraction);
        }
        let baths = times.map(|time| patches[0].at(time).unwrap().ambient_temperature_k());
        for (stage, time, bath) in [(&step.first_stage, times[0], baths[0]),
            (&step.second_stage, times[1], baths[1])]
        {
            admit_stage(stage, bath, 3. + 10. * time + 2. / 3., &card);
        }
        assert_ne!(baths[0], baths[1]);
        admit_interval(&step, &old, baths);
        close(step.source_j, dt * (3. + 10. * (0.1 + end) / 2. + 2. / 3.), 2e-11);
        assert_eq!(continued.conduction().accepted_steps(), 2);
        assert_eq!(continued.conduction().time_s(), end);
        assert_eq!(continued.conduction().specific_enthalpy_j_kg(),
            step.second_stage.specific_enthalpy_j_kg.as_slice());
    });
}

#[test]
fn g3_g4_pointwise_sdirk_handoff_retains_bath_event_sides_and_atomic_resume() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh(); let history = control_history(cx, &mesh); let card = card();
        let curves = [chart(false), chart(true)];
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let source = heater(constant(0.1, LENGTH), constant(3., COMPONENT_POWER_DIMS));
        let mut patches = radiation().patches;
        patches[0].ambient_temperature = TimeSignal::Table {
            times: vec![0., 0.12, 0.18], values: vec![300., 330., 290.],
            dims: TEMPERATURE_DIMS, interp: Interp::Hold,
        };
        let mut evolving = scan(cx, &history, &material, &curves, &card, source, patches.clone());
        evolving.advance_to(cx, 0.04, coupled_newton().step).unwrap();
        let mut run = evolving.into_sdirk2().unwrap();
        let initial = run.conduction().specific_enthalpy_j_kg().to_vec();
        let temperature = run.conduction().temperature().to_vec();
        let fractions = run.conduction().liquid_mass_fraction().to_vec();
        let gate = CancelGate::new_clock_free(); gate.request();
        context(&gate, |cancelled| {
            assert!(run.advance_event_step(cancelled, 0.24, 1., coupled_newton().step).is_err());
        });
        // The retained Newton policy also remains authoritative after transfer.
        assert!(run.advance_to(cx, 0.08, config()).is_err());
        assert_eq!(run.conduction().specific_enthalpy_j_kg(), initial.as_slice());
        assert_eq!(run.conduction().temperature(), temperature.as_slice());
        assert_eq!(run.conduction().liquid_mass_fraction(), fractions.as_slice());
        assert_eq!(run.conduction().time_s(), 0.04);
        assert_eq!(run.conduction().accepted_steps(), 1);
        let mut old = initial;
        let mut from = 0.04;
        for (to, bath) in [(0.12, 300.), (0.18, 330.)] {
            let interval = run.advance_event_step(cx, 0.24, 1., coupled_newton().step).unwrap();
            assert_eq!(interval.start_time_s, from);
            assert_eq!(interval.end_time_s, to);
            admit_stage(&interval.report.first_stage, bath, 3. + 2. / 3., &card);
            admit_stage(&interval.report.second_stage, bath, 3. + 2. / 3., &card);
            admit_interval(&interval.report, &old, [bath; 2]);
            close(interval.report.source_j, (to - from) * (3. + 2. / 3.), 2e-11);
            old = interval.report.second_stage.specific_enthalpy_j_kg;
            from = to;
        }
        let temperature = run.conduction().temperature().to_vec();
        let fractions = run.conduction().liquid_mass_fraction().to_vec();
        assert!(run.advance_event_step(cx, 0.24, 1., coupled_newton().step).is_err());
        assert_eq!(run.conduction().specific_enthalpy_j_kg(), old.as_slice());
        assert_eq!(run.conduction().temperature(), temperature.as_slice());
        assert_eq!(run.conduction().liquid_mass_fraction(), fractions.as_slice());
        assert_eq!(run.conduction().time_s(), 0.18);
        assert_eq!(run.conduction().accepted_steps(), 3);
        assert_eq!(run.pointwise_ambient_patches(), Some(patches.as_slice()));
        run.extend_step_budget(4).unwrap();
        let resumed = run.advance_event_step(cx, 0.24, 1., coupled_newton().step).unwrap();
        assert_eq!(resumed.end_time_s, 0.24);
        admit_stage(&resumed.report.first_stage, 290., 3. + 2. / 3., &card);
        admit_stage(&resumed.report.second_stage, 290., 3. + 2. / 3., &card);
        admit_interval(&resumed.report, &old, [290.; 2]);
        assert_eq!(run.conduction().accepted_steps(), 4);
    });
}

#[test]
fn g4_pointwise_sdirk_attachment_uses_native_natural_traces_and_exclusive_ownership() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = separated_mesh(); let curves = [chart(false), chart(true)];
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let history = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&mesh)
            .region("exterior", |face| face.centroid[2] < 0. && face.centroid[1] == 0.,
                ThermalBc::neumann(2.).unwrap()).unwrap().adiabatic_remainder(),
            &Environment::earth_lab(), Vec::new(), 0.).unwrap();
        let create = || base(cx, &history, &material, &curves, 3).into_sdirk2();
        let patches = radiation().patches;
        let mut run = create().with_pointwise_ambient(cx, patches.clone()).unwrap();
        let step = run.advance_to(cx, 0.001, config()).unwrap();
        assert!(step.first_stage.pointwise_radiation.is_some());
        assert!(step.second_stage.pointwise_radiation.is_some());
        close(step.neumann_out_j, 0.001, 1e-12); // 2 W/m^2 on one area-1/2 face.
        close(step.stored_energy_change_j + step.neumann_out_j
            + step.radiation_out_j - step.source_j, 0., config().energy_tolerance_j);
        assert!(run.with_pointwise_ambient(cx, patches.clone()).is_err());
        assert!(create().with_pointwise_ambient(cx, patches.clone()).unwrap()
            .with_pointwise_ambient(cx, patches.clone()).is_err());
        assert!(create().with_pointwise_ambient(cx, vec![]).is_err());
        assert!(create().with_pointwise_ambient(cx, vec![patches[0].clone(), patches[0].clone()]).is_err());
        let gate = CancelGate::new_clock_free(); gate.request();
        context(&gate, |cancelled| {
            assert!(create().with_pointwise_ambient(cancelled, patches.clone()).is_err());
        });
        let mut advanced = create();
        advanced.advance_to(cx, 0.001, config()).unwrap();
        assert!(advanced.with_pointwise_ambient(cx, patches).is_err());
    });
}
