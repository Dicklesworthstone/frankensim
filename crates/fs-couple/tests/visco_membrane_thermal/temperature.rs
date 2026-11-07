//! Actual-temperature kinetics through the original native and caloric owners.
use super::*;
use fs_couple::membrane_thermal::visco::temperature::{
    ViscoMembraneKinetics, ViscoMembraneTemperatureError, advance_candidate as advance_temperature,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use fs_solid::membrane::visco::owned::OwnedFiniteViscoMembraneMotion;

fn calibration() -> TemperatureShift {
    // Explicit synthetic kinetic data, independent of the sheet's elastic card.
    TemperatureShift::new(
        TemperatureShiftModel::WilliamsLandelFerry { c1: 2.0, c2_k: 2.0 },
        300.0,
        [299.0, 301.0],
    )
    .unwrap()
}

fn point_step(
    accepted: &ViscoMembraneState,
    law: &ViscoHyperelastic,
    dt_s: f64,
) -> (ViscoMembraneState, f64) {
    let ViscoMembraneState::Incompressible(accepted) = accepted else {
        unreachable!()
    };
    let f = accepted.deformation_gradient;
    let response = law
        .step_incompressible_plane_stress_response(
            accepted,
            &[f[0], f[1], f[3], f[4]],
            dt_s,
            None,
            local(),
        )
        .unwrap();
    let heat = AREA_M2 * THICKNESS_M * response.material.dissipated_pa;
    (
        ViscoMembraneState::Incompressible(response.material.state),
        heat,
    )
}

fn same_point(actual: &ViscoMembraneState, expected: &ViscoMembraneState) {
    let (ViscoMembraneState::Incompressible(actual), ViscoMembraneState::Incompressible(expected)) =
        (actual, expected)
    else {
        unreachable!()
    };
    assert_eq!(
        actual.inverse_viscous_cauchy_green.len(),
        expected.inverse_viscous_cauchy_green.len()
    );
    for (a, b) in actual
        .deformation_gradient
        .iter()
        .chain(actual.inverse_viscous_cauchy_green.iter().flatten())
        .zip(
            expected
                .deformation_gradient
                .iter()
                .chain(expected.inverse_viscous_cauchy_green.iter().flatten()),
        )
    {
        near(*a, *b, 2e-10);
    }
    for (a, b) in [
        actual.work,
        actual.reference_work,
        actual.stored,
        actual.dissipated,
    ]
    .into_iter()
    .zip([
        expected.work,
        expected.reference_work,
        expected.stored,
        expected.dissipated,
    ]) {
        near(a, b, 2e-9);
    }
}

fn law_at(motion: &OwnedFiniteViscoMembraneMotion, element: usize) -> &ViscoHyperelastic {
    match motion.native_laws() {
        ViscoelasticLawField::Uniform(law) => law,
        ViscoelasticLawField::PerElement(laws) => &laws[element],
    }
}

#[test]
fn held_actual_temperature_matches_native_shifted_and_reduced_time_returns() {
    let (patch, laws, caloric) = fixture(false, &[[0.8; 3]], &[0.08], 0);
    let reference = motion(&patch, &laws, true);
    let original = OwnedFiniteViscoMembraneMotion::from_native(&reference);
    let shift = calibration();
    let kinetics = ViscoMembraneKinetics::try_new(&reference, &[shift], || false).unwrap();
    let temperature = 300.5;
    let dt = 0.015625;
    let shifted = shift.shift_hyperelastic(&laws[0], temperature).unwrap();
    let accepted = &reference.membrane().triangle_states()[0];
    let (expected, heat_j) = point_step(accepted, shifted.law(), dt);
    let (reduced_time, reduced_heat) = point_step(
        accepted,
        &laws[0],
        dt / shift.evaluate(temperature).unwrap().factor(),
    );
    same_point(&expected, &reduced_time);
    near(heat_j, reduced_heat, 2e-12);
    assert!(heat_j > 0.0);
    let mass = AREA_M2 * 0.8;
    // Cp=1 here. Preheat without moving either clock, then remove exactly the
    // independently predicted heat through an authored boundary sink.
    let heat = thermal(
        &caloric,
        &[(30, mass)],
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Boundary(CellId(30), boundary(-heat_j / dt)),
        }],
    )
    .apply_heat_j(
        &BTreeMap::from([(CellId(30), mass * 0.5)]),
        1e-12,
        1e-12,
        || Ok(()),
    )
    .unwrap()
    .0;
    let saved = snapshot(&original.as_native());
    let next = advance_temperature(
        &original,
        &heat,
        &(),
        dt,
        &vec![[0.0; 3]; 3],
        &held(&reference),
        &cells(&caloric, &[30]),
        &kinetics,
        mechanics(),
        config(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    same_point(
        &next.motion().as_native().membrane().triangle_states()[0],
        &expected,
    );
    near(
        next.thermal().temperature_k(CellId(30)).unwrap(),
        temperature,
        2e-10,
    );
    near(next.report().material_heat_j, heat_j, 2e-12);
    near(next.report().deposited_heat_j, heat_j, 2e-12);
    near(next.report().mechanical_energy_change_j, -heat_j, 2e-12);
    near(next.report().thermal.stored_energy_change_j, 0.0, 2e-12);
    near(next.report().external_thermal_heat_j, -heat_j, 2e-12);
    near(next.report().combined_residual_j, 0.0, 2e-12);
    let final_shift = shift
        .shift_hyperelastic(&laws[0], next.thermal().temperature_k(CellId(30)).unwrap())
        .unwrap();
    assert_eq!(law_at(next.motion(), 0), final_shift.law());
    assert_eq!(snapshot(&original.as_native()), saved);
    assert_eq!(kinetics.original_laws(), laws);
    assert_eq!(heat.time_s(), 0.0);
}

#[test]
fn self_heating_accelerates_relaxation_and_refines_with_the_physical_timestep() {
    let (patch, laws, caloric) = fixture(true, &[[0.8; 3]], &[0.1], 0);
    let reference = motion(&patch, &laws, true);
    let original = OwnedFiniteViscoMembraneMotion::from_native(&reference);
    let shift = calibration();
    let kinetics = ViscoMembraneKinetics::try_new(&reference, &[shift], || false).unwrap();
    let heat = thermal(&caloric, &[(30, AREA_M2 * 0.8)], vec![]);
    let bindings = cells(&caloric, &[30]);
    let fixed = held(&reference);
    let load = vec![[0.0; 3]; 3];
    let target = 0.0625;
    let timesteps = [0.015625, 0.0078125, 0.00390625];
    let mut totals = Vec::new();
    for dt in timesteps {
        let controls = ViscoMembraneThermalControls {
            maximum_step_s: dt,
            ..mechanics()
        };
        let next = advance_temperature(
            &original,
            &heat,
            &(),
            target,
            &load,
            &fixed,
            &bindings,
            &kinetics,
            controls,
            config(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
        let actual_t = next.thermal().temperature_k(CellId(30)).unwrap();
        assert!(actual_t > 300.0 && actual_t < 301.0);
        let expected_law = shift.shift_hyperelastic(&laws[0], actual_t).unwrap();
        assert_eq!(law_at(next.motion(), 0), expected_law.law());
        assert!(law_at(next.motion(), 0).branches[0].tau_s < laws[0].branches[0].tau_s);
        near(
            next.report().mechanical_energy_change_j,
            -next.report().material_heat_j,
            2e-12,
        );
        near(
            next.report().thermal.stored_energy_change_j,
            next.report().material_heat_j,
            2e-12,
        );
        near(next.report().external_thermal_heat_j, 0.0, 2e-12);
        near(next.report().combined_residual_j, 0.0, 2e-12);
        totals.push(next.report().material_heat_j);
    }
    let coarse_change = (totals[1] - totals[0]).abs();
    let fine_change = (totals[2] - totals[1]).abs();
    assert!(coarse_change > 1e-11);
    assert!(
        fine_change < coarse_change,
        "heat refinement did not contract: {totals:?}"
    );
    let frozen = advance_candidate(
        &reference,
        &heat,
        &(),
        target,
        &load,
        &fixed,
        &bindings,
        ViscoMembraneThermalControls {
            maximum_step_s: timesteps[2],
            ..mechanics()
        },
        config(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert!(totals[2] > frozen.report().material_heat_j + 1e-9);
    assert_eq!(kinetics.original_laws(), laws);
    assert_eq!(law_at(&original, 0), &laws[0]);
}

#[test]
fn regional_actual_temperatures_and_continuation_keep_the_original_spectra() {
    let density = [[0.8; 3], [1.0, 1.6, 2.2]];
    let (patch, laws, caloric) = fixture(false, &density, &[0.12, 0.3], 0);
    let reference = motion(&patch, &laws, true);
    let original = OwnedFiniteViscoMembraneMotion::from_native(&reference);
    let shifts = [
        calibration(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 120_000.0,
            },
            300.0,
            [299.0, 301.0],
        )
        .unwrap(),
    ];
    let kinetics = ViscoMembraneKinetics::try_new(&reference, &shifts, || false).unwrap();
    let ids = [90, 10];
    let masses = density.map(|rho| AREA_M2 * rho.iter().sum::<f64>() / 3.0);
    let dt = 0.015625;
    let rates = [0.001, -0.002];
    let initial_heat = thermal(
        &caloric,
        &[(ids[0], masses[0]), (ids[1], masses[1])],
        ids.into_iter()
            .zip(rates)
            .map(|(id, watts)| TimedEvent {
                time_s: dt,
                event: Event::Boundary(CellId(id), boundary(watts)),
            })
            .collect(),
    )
    .apply_heat_j(
        &BTreeMap::from([
            (CellId(ids[0]), -0.5 * masses[0]),
            (CellId(ids[1]), 0.5 * masses[1]),
        ]),
        1e-12,
        1e-12,
        || Ok(()),
    )
    .unwrap()
    .0;
    near(
        initial_heat.temperature_k(CellId(ids[0])).unwrap(),
        299.5,
        2e-10,
    );
    near(
        initial_heat.temperature_k(CellId(ids[1])).unwrap(),
        300.5,
        2e-10,
    );
    let mut accepted = original.clone();
    let mut heat = initial_heat.clone();
    let mut expected_states = reference.membrane().triangle_states().to_vec();
    let mut total_heat = 0.0;
    let fixed = held(&reference);
    let load = vec![[0.0; 3]; 6];
    let bindings = cells(&caloric, &ids);
    let controls = ViscoMembraneThermalControls {
        maximum_step_s: dt,
        ..mechanics()
    };
    for interval in 0..3 {
        let start = heat.time_s();
        let mut expected_heat = [0.0; 2];
        for e in 0..2 {
            let temperature = heat.temperature_k(CellId(ids[e])).unwrap();
            let shifted = shifts[e].shift_hyperelastic(&laws[e], temperature).unwrap();
            let (state, loss) = point_step(&expected_states[e], shifted.law(), dt);
            expected_states[e] = state;
            expected_heat[e] = loss;
        }
        let next = advance_temperature(
            &accepted,
            &heat,
            &(),
            start + dt,
            &load,
            &fixed,
            &bindings,
            &kinetics,
            controls,
            config(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
        let native = next.motion().as_native();
        for e in 0..2 {
            same_point(&native.membrane().triangle_states()[e], &expected_states[e]);
            near(next.report().element_heat_j[e], expected_heat[e], 2e-12);
            let external = if interval == 0 { 0.0 } else { rates[e] * dt };
            near(
                next.thermal().cells()[&CellId(ids[e])].specific_h_j_kg
                    - heat.cells()[&CellId(ids[e])].specific_h_j_kg,
                (expected_heat[e] + external) / masses[e],
                2e-10,
            );
            let endpoint_t = next.thermal().temperature_k(CellId(ids[e])).unwrap();
            let endpoint_law = shifts[e].shift_hyperelastic(&laws[e], endpoint_t).unwrap();
            // This checks the published endpoint law, not just the law used at
            // the interval start, and would reject repeatedly shifting tau.
            assert_eq!(law_at(next.motion(), e), endpoint_law.law());
        }
        near(
            next.report().external_thermal_heat_j,
            if interval == 0 {
                0.0
            } else {
                rates.iter().sum::<f64>() * dt
            },
            2e-12,
        );
        total_heat += next.report().material_heat_j;
        let (m, t, _, _) = next.into_parts();
        accepted = m;
        heat = t;
    }
    // One call and three continued calls traverse the same actual intervals,
    // including the scheduled boundary changes at the first endpoint.
    let whole = advance_temperature(
        &original,
        &initial_heat,
        &(),
        3.0 * dt,
        &load,
        &fixed,
        &bindings,
        &kinetics,
        controls,
        config(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    assert_eq!(
        snapshot(&whole.motion().as_native()),
        snapshot(&accepted.as_native())
    );
    assert_eq!(whole.thermal().cells(), heat.cells());
    near(whole.report().material_heat_j, total_heat, 2e-12);
    assert_eq!(kinetics.original_laws(), laws);
    assert_eq!(kinetics.shifts(), shifts);
}

#[test]
fn calibration_law_and_late_thermal_refusals_preserve_accepted_state() {
    let (patch, laws, caloric) = fixture(false, &[[0.8; 3]], &[0.05], 0);
    let reference = motion(&patch, &laws, true);
    let original = OwnedFiniteViscoMembraneMotion::from_native(&reference);
    let shift = calibration();
    assert!(matches!(
        ViscoMembraneKinetics::try_new(&reference, &[], || false),
        Err(ViscoMembraneTemperatureError::Invalid(_))
    ));
    let kinetics = ViscoMembraneKinetics::try_new(&reference, &[shift], || false).unwrap();
    let mass = AREA_M2 * 0.8;
    let heat = thermal(&caloric, &[(30, mass)], vec![]);
    let bindings = cells(&caloric, &[30]);
    let fixed = held(&reference);
    let load = vec![[0.0; 3]; 3];
    let dt = 0.015625;
    let controls = ViscoMembraneThermalControls {
        maximum_step_s: dt,
        ..mechanics()
    };
    let gate = CancelGate::new_clock_free();
    let saved = snapshot(&original.as_native());
    let saved_heat = heat.cells().clone();
    let mut foreign_laws = laws.clone();
    foreign_laws[0].branches[0].energy = BranchEnergy::NeoHookean { mu_pa: 24.0 };
    let foreign_native = motion(&patch, &foreign_laws, true);
    let foreign = OwnedFiniteViscoMembraneMotion::from_native(&foreign_native);
    let saved_foreign = snapshot(&foreign.as_native());
    assert!(matches!(
        advance_temperature(
            &foreign,
            &heat,
            &(),
            dt,
            &load,
            &fixed,
            &bindings,
            &kinetics,
            controls,
            config(),
            &gate,
        ),
        Err(MembraneThermalError::Viscoelastic(_))
    ));
    assert_eq!(snapshot(&foreign.as_native()), saved_foreign);
    let failing = Observer {
        fail_after: Some(0.02),
        ..Observer::default()
    };
    let saved_observer = failing.clone();
    assert!(matches!(
        advance_temperature(
            &original,
            &heat,
            &failing,
            2.0 * dt,
            &load,
            &fixed,
            &bindings,
            &kinetics,
            controls,
            config(),
            &gate,
        ),
        Err(MembraneThermalError::Observer(
            DepositionError::InvalidInput("synthetic late observer refusal")
        ))
    ));
    assert_eq!(failing, saved_observer);
    // Pick a kinetic band whose endpoint lies strictly between the first two
    // physical thermal endpoints. The first interval is demonstrably admissible;
    // the second must fail when its newly heated endpoint is synchronized.
    let (first_state, first_heat) =
        point_step(&reference.membrane().triangle_states()[0], &laws[0], dt);
    let first_t = 300.0 + first_heat / mass;
    let second_law = shift.shift_hyperelastic(&laws[0], first_t).unwrap();
    let (_, second_heat) = point_step(&first_state, second_law.law(), dt);
    assert!(first_heat > 0.0 && second_heat > 0.0);
    let limit = first_t + 0.5 * second_heat / mass;
    let narrow = TemperatureShift::new(shift.model(), 300.0, [299.0, limit]).unwrap();
    let limited = ViscoMembraneKinetics::try_new(&reference, &[narrow], || false).unwrap();
    let prefix = advance_temperature(
        &original,
        &heat,
        &(),
        dt,
        &load,
        &fixed,
        &bindings,
        &limited,
        controls,
        config(),
        &gate,
    )
    .unwrap();
    assert!(prefix.thermal().temperature_k(CellId(30)).unwrap() < limit);
    assert!(matches!(
        advance_temperature(
            &original,
            &heat,
            &(),
            2.0 * dt,
            &load,
            &fixed,
            &bindings,
            &limited,
            controls,
            config(),
            &gate,
        ),
        Err(MembraneThermalError::Temperature(
            ViscoMembraneTemperatureError::Shift { .. }
        ))
    ));
    assert_eq!(snapshot(&original.as_native()), saved);
    assert_eq!(heat.cells(), &saved_heat);
    assert_eq!(heat.time_s(), 0.0);
    assert_eq!(law_at(&original, 0), &laws[0]);
    assert_eq!(kinetics.original_laws(), laws);
    let retry = advance_temperature(
        &original,
        &heat,
        &(),
        2.0 * dt,
        &load,
        &fixed,
        &bindings,
        &kinetics,
        controls,
        config(),
        &gate,
    )
    .unwrap();
    assert!(retry.report().deposited_heat_j > first_heat);
}
