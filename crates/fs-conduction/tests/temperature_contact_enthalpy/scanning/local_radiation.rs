//! Local P1 radiation remains physical throughout the original scan trajectory.
use super::*;

// Exact barycentric integral: the average of a linear trace's fourth power is
// 1/15 times the complete homogeneous polynomial of degree four in its nodes.
// This does not call the production quadrature or replace T_h^4 by mean(T)^4.
fn fourth_moment(t: [f64; 3]) -> f64 {
    let mut total = 0.0;
    for i in 0..=4 {
        for j in 0..=4 - i {
            total += t[0].powi(i) * t[1].powi(j) * t[2].powi(4 - i - j);
        }
    }
    total / 15.0
}

#[test]
fn g1_scan_pointwise_radiation_keeps_natural_flux_and_shared_bath_reactions() {
    context(&CancelGate::new_clock_free(), |cx| {
        let card = card();
        let curves = [chart(false), chart(true)];
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let source = heater(TimeSignal::Ramp { t_start: 0., t_end: 0.1,
            from: QtyAny::new(-0.4, LENGTH), to: QtyAny::new(0.6, LENGTH) },
            constant(2., COMPONENT_POWER_DIMS));
        for mode in 0..3 {
            let mesh = if mode == 0 { separated_mesh() } else { mesh() };
            let history = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&mesh)
                .region("bath", |face| face.centroid[2] < 0. && face.centroid[0] == 0.,
                    ThermalBc::dirichlet(305.).unwrap()).unwrap()
                .region("exterior", |face| face.centroid[2] < 0. && face.centroid[1] == 0.,
                    ThermalBc::neumann(4.).unwrap()).unwrap()
                .adiabatic_remainder(), &Environment::earth_lab(), Vec::new(), 0.).unwrap();
            let initial = history.at(cx, 0.).unwrap();
            let scan = ScanningEnthalpyRun::new(cx, base(cx, &history, &material, &curves, 1),
                vec![source.clone()], 1, None).unwrap();
            let mut scan = if mode == 0 {
                scan.with_pointwise_ambient_radiation(cx, radiation().patches).unwrap()
            } else {
                let contacts = TemperatureContacts::new(cx, &mesh, &initial, vec![],
                    vec![law(&mesh, &card)], 1).unwrap();
                if mode == 1 {
                    scan.with_contacts(cx, contacts, coupling(), 1).unwrap()
                        .with_pointwise_ambient_radiation(cx, radiation().patches).unwrap()
                } else {
                    // Either attachment order keeps the phase-aware Newton owner.
                    scan.with_pointwise_ambient_radiation(cx, radiation().patches).unwrap()
                        .with_contacts_newton(cx, contacts, coupled_newton()).unwrap()
                }
            };
            assert_eq!(scan.pointwise_ambient_patches().unwrap()[0].region, "exterior");
            let controls = if mode == 2 { coupled_newton().step } else { config() };
            let step = scan.advance_to(cx, 0.04, controls).unwrap();
            let solution = step.solution();
            let local = solution.pointwise_radiation.as_ref().expect("local radiation retained");
            assert!(solution.pointwise_enclosure_radiation.is_none());
            let trace = [solution.temperature[0], solution.temperature[1], solution.temperature[3]];
            let power = 0.2 * STEFAN_BOLTZMANN_W_M2_K4 * 0.5
                * (fourth_moment(trace) - 300_f64.powi(4));
            let mean = trace.iter().sum::<f64>() / 3.;
            let mean_power = 0.2 * STEFAN_BOLTZMANN_W_M2_K4 * 0.5
                * (mean.powi(4) - 300_f64.powi(4));
            assert!((power - mean_power).abs() > 1e-8,
                "nonuniform pointwise emission collapsed to an area-mean load");
            close(local.outward_heat_w, power, 1e-8);
            close(solution.neumann_out_w, 2., 1e-12); // Area 1/2 times 4 W/m2.
            close(solution.robin_out_w, 0., 1e-12);
            close(solution.source_w, 2. + 2. / 3., 1e-11);
            close(step.deposition[0].center_m[2], 0., 1e-14);
            for vertex in [0, 2, 3] { close(solution.temperature[vertex], 305., 1e-12); }
            let expected = 0.04 * (solution.source_w + solution.dirichlet_in_w
                - solution.neumann_out_w - power);
            close(stored(&OLD, &solution.specific_enthalpy_j_kg), expected, 3e-8);
            close(solution.energy_residual_j, 0., controls.energy_tolerance_j);
            if mode != 0 {
                let contact = coupled(&step);
                assert!(contact.radiation.is_none() && contact.combined_boundary.is_none());
                close(contact.physical.radiation_out_w, power, 1e-8);
                close(stored(&OLD, &solution.specific_enthalpy_j_kg),
                    0.04 * (solution.source_w + contact.physical.dirichlet_in_w
                        - solution.neumann_out_w - power), 3e-8);
                close(contact.physical.energy_residual_j, 0., controls.energy_tolerance_j);
                assert!(contact.physical.residual_norm_j <= coupling().residual_tolerance_j);
            }
        }
    });
}

#[test]
fn g3_g4_pointwise_scan_bath_events_adaptation_and_refusals_preserve_history() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let card = card();
        let curves = [chart(false), chart(true)];
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let history = ThermalBoundaryHistory::new(cx, ThermalBoundaryBuilder::new(&mesh)
            .region("exterior", |face| face.centroid[2] < 0. && face.centroid[1] == 0.,
                ThermalBc::neumann(0.).unwrap()).unwrap()
            .adiabatic_remainder(), &Environment::earth_lab(), Vec::new(), 0.).unwrap();
        let initial = history.at(cx, 0.).unwrap();
        let mut patches = radiation().patches;
        patches[0].ambient_temperature = TimeSignal::Table {
            times: vec![0., 0.03, 0.06, 0.09], values: vec![300., 330., 310., -1.],
            dims: TEMPERATURE_DIMS, interp: Interp::Hold,
        };
        let source = heater(constant(0.4, LENGTH), TimeSignal::Table {
            times: vec![0., 0.03, 0.06], values: vec![1., 3., 0.],
            dims: COMPONENT_POWER_DIMS, interp: Interp::Hold,
        });
        let create = || {
            let contacts = TemperatureContacts::new(cx, &mesh, &initial, vec![],
                vec![law(&mesh, &card)], 1).unwrap();
            ScanningEnthalpyRun::new(cx, base(cx, &history, &material, &curves, 4),
                vec![source.clone()], 1, None).unwrap()
                .with_contacts_newton(cx, contacts, coupled_newton()).unwrap()
                .with_pointwise_ambient_radiation(cx, patches.clone()).unwrap()
        };
        let mut scan = create();
        let mut manual = create();
        for (end, bath, absorbed) in [(0.03, 300., 1.), (0.06, 330., 3.)] {
            let before = scan.conduction().specific_enthalpy_j_kg().to_vec();
            let cancelled = CancelGate::new_clock_free(); cancelled.request();
            context(&cancelled, |cancelled_cx| {
                assert!(scan.advance_event_step(cancelled_cx, 0.06, 1., coupled_newton().step).is_err());
            });
            assert_eq!(scan.conduction().specific_enthalpy_j_kg(), before.as_slice());
            let accepted = scan.advance_event_step(cx, 0.06, 1., coupled_newton().step).unwrap();
            let reference = manual.advance_event_step(cx, 0.06, 1., coupled_newton().step).unwrap();
            close(accepted.end_time_s, end, 0.);
            close(accepted.report.deposition[0].absorbed_power_w, absorbed, 0.);
            let radiation = accepted.report.solution().pointwise_radiation.as_ref().unwrap();
            close(radiation.patches[0].patch.ambient_temperature_k(), bath, 0.);
            let physical = &coupled(&accepted.report).physical;
            close(stored(&before, &accepted.report.solution().specific_enthalpy_j_kg),
                (accepted.end_time_s - accepted.start_time_s)
                    * (absorbed + 2. / 3. - radiation.outward_heat_w), 3e-8);
            close(physical.energy_residual_j, 0., coupled_newton().step.energy_tolerance_j);
            assert_eq!(accepted.report.solution().specific_enthalpy_j_kg,
                reference.report.solution().specific_enthalpy_j_kg);
        }
        let before = scan.conduction().specific_enthalpy_j_kg().to_vec();
        let pair = scan.advance_adaptive_step(cx, 0.08, 0.02, coupled_newton().step,
            AdaptiveHistoryConfig { absolute_temperature_k: 1000., relative_temperature: 0.,
                min_step_s: 1e-8, max_step_s: 0.02, max_attempts: 2 },
            PhaseErrorTolerance { absolute_enthalpy_j_kg: 1000., relative_enthalpy: 0.,
                absolute_liquid_fraction: 1. }).unwrap();
        assert_eq!(pair.rejected_attempts, 0);
        let mut old = before.clone();
        let mut from = pair.start_time_s;
        let mut integrated = 0.;
        for (end, half) in [(pair.midpoint_time_s, &pair.first_half), (pair.end_time_s, &pair.second_half)] {
            let reference = manual.advance_to(cx, end, coupled_newton().step).unwrap();
            assert_eq!(half.solution().specific_enthalpy_j_kg, reference.solution().specific_enthalpy_j_kg);
            let radiation = half.solution().pointwise_radiation.as_ref().unwrap();
            close(radiation.patches[0].patch.ambient_temperature_k(), 310., 0.);
            let input = (end - from) * (2. / 3. - radiation.outward_heat_w);
            close(stored(&old, &half.solution().specific_enthalpy_j_kg), input, 3e-8);
            integrated += input;
            old = half.solution().specific_enthalpy_j_kg.clone();
            from = end;
        }
        close(stored(&before, scan.conduction().specific_enthalpy_j_kg()), integrated, 6e-8);
        assert_eq!(scan.conduction().accepted_steps(), 4);
        scan.extend_step_budget(5).unwrap();
        let before = scan.conduction().specific_enthalpy_j_kg().to_vec();
        let temperature = scan.conduction().temperature().to_vec();
        let phase = scan.conduction().liquid_mass_fraction().to_vec();
        assert!(scan.advance_to(cx, 0.1, coupled_newton().step).is_err());
        assert_eq!(scan.conduction().specific_enthalpy_j_kg(), before.as_slice());
        assert_eq!(scan.conduction().temperature(), temperature.as_slice());
        assert_eq!(scan.conduction().liquid_mass_fraction(), phase.as_slice());
        close(scan.conduction().time_s(), pair.end_time_s, 0.);
        assert_eq!(scan.conduction().accepted_steps(), 4);
        let resumed = scan.advance_event_step(cx, 0.1, 1., coupled_newton().step).unwrap();
        close(resumed.end_time_s, 0.09, 0.);
        close(resumed.report.solution().pointwise_radiation.as_ref().unwrap().patches[0]
            .patch.ambient_temperature_k(), 310., 0.);
        assert_eq!(scan.conduction().accepted_steps(), 5);
    });
}
