use super::*;

fn controls(lo: f64, hi: f64) -> IsochoricWaterControls {
    IsochoricWaterControls {
        minimum_temperature_k: lo,
        maximum_temperature_k: hi,
        minimum_pressure_pa: 100.0,
        maximum_pressure_pa: 100e6,
        volume_tolerance_m3_kg: 1e-14,
        internal_energy_tolerance_j_kg: 1e-6,
        maximum_pressure_iterations: 96,
        maximum_temperature_iterations: 96,
    }
}
fn pressure() -> DensePressureSolve {
    DensePressureSolve {
        pressure_tolerance_pa: 1e-4,
        density_tolerance_kg_m3: 1e-8,
        maximum_iterations: 128,
    }
}
fn pair(t: f64) -> If97FluidSaturationPair {
    If97FluidSaturationPair::at_temperature(Temperature::new(t), pressure(), || false).unwrap()
}
fn mixed_volume(pair: If97FluidSaturationPair, quality: f64) -> f64 {
    let vl = pair.liquid().specific_volume().value();
    vl + quality * (pair.vapor().specific_volume().value() - vl)
}
fn law(v: f64, lo: f64, hi: f64) -> IsochoricWater {
    IsochoricWater::new_regional_coexistence(
        MassSpecificVolume::new(v),
        controls(lo, hi),
        pressure(),
        || false,
    )
    .unwrap()
}
fn at(law: &IsochoricWater, t: f64) -> IsochoricWaterState {
    law.state_at_temperature(Temperature::new(t), || false)
        .unwrap()
}
fn near(actual: f64, expected: f64, absolute: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= absolute,
        "{actual:e} != {expected:e}, tolerance={absolute:e}"
    );
}
fn original_u_at_fixed_volume(t: f64, v: f64) -> f64 {
    // Independent finite-displacement oracle: use only original u and v,
    // never Cp, thermal expansion, the production quality derivative or Cv.
    let pair = pair(t);
    let l = pair.liquid();
    let g = pair.vapor();
    let x = (v - l.specific_volume().value())
        / (g.specific_volume().value() - l.specific_volume().value());
    l.specific_internal_energy().value()
        + x * (g.specific_internal_energy().value() - l.specific_internal_energy().value())
}

#[test]
fn g3_gibbs_coexistence_mode_preserves_original_isochore_values() {
    let cfg = controls(350.0, 400.0);
    let legacy = IsochoricWater::new(MassSpecificVolume::new(0.1), cfg, || false).unwrap();
    let regional = law(0.1, 350.0, 400.0);
    assert_eq!(legacy.regional_coexistence_controls(), None);
    assert_eq!(regional.regional_coexistence_controls(), Some(pressure()));
    for t in [350.0, 375.0, 400.0] {
        let a = at(&legacy, t);
        let b = at(&regional, t);
        assert_eq!(a.saturation_pair(), b.saturation_pair());
        assert_eq!(a.regional_saturation_pair(), b.regional_saturation_pair());
        assert_eq!(a.vapor_quality(), b.vapor_quality());
        for (x, y) in [
            (a.pressure_pa(), b.pressure_pa()),
            (a.specific_volume_m3_kg(), b.specific_volume_m3_kg()),
            (a.bulk_density_kg_m3(), b.bulk_density_kg_m3()),
            (
                a.specific_internal_energy_j_kg(),
                b.specific_internal_energy_j_kg(),
            ),
            (a.specific_enthalpy_j_kg(), b.specific_enthalpy_j_kg()),
            (a.specific_entropy_j_kg_k(), b.specific_entropy_j_kg_k()),
            (
                a.heat_capacity_at_volume_j_kg_k(),
                b.heat_capacity_at_volume_j_kg_k(),
            ),
        ] {
            assert_eq!(x.to_bits(), y.to_bits());
        }
        assert_eq!(b.source_pressure_residuals_pa(), [0.0; 2]);
        assert_eq!(b.pressure_iterations(), 0);
    }
}

#[test]
fn g1_dense_isochore_keeps_original_pair_density_and_caloric_tangent() {
    for t in [630.0, 640.0, 646.0] {
        let original = pair(t);
        let v = mixed_volume(original, 0.2);
        let law = law(v, t - 0.1, t + 0.1);
        let state = at(&law, t);
        assert_eq!(state.regional_saturation_pair(), Some(original));
        assert_eq!(state.saturation_pair(), None);
        assert_eq!(state.single_phase_source(), None);
        assert_eq!(state.pressure_pa(), original.pressure().value());
        assert!(state.pressure_iterations() > 0);
        let x = state.vapor_quality().unwrap();
        near(x, 0.2, 5e-15);
        let mass = 0.37;
        let volumes = [
            mass * (1.0 - x) * original.liquid().specific_volume().value(),
            mass * x * original.vapor().specific_volume().value(),
        ];
        near(volumes[0] + volumes[1], mass * v, 3e-18);
        near(
            state.bulk_density_kg_m3(),
            mass / (volumes[0] + volumes[1]),
            2e-12,
        );
        let residuals = [
            original.liquid().pressure().value() - original.pressure().value(),
            original.vapor().pressure().value() - original.pressure().value(),
        ];
        assert_eq!(state.source_pressure_residuals_pa(), residuals);
        let phase_work_defect =
            (1.0 - x) * residuals[0] * original.liquid().specific_volume().value()
                + x * residuals[1] * original.vapor().specific_volume().value();
        near(
            state.pressure_closure_defect_j_kg(),
            phase_work_defect,
            32.0 * f64::EPSILON * state.specific_enthalpy_j_kg().abs(),
        );
        let difference = |h| {
            (original_u_at_fixed_volume(t + h, v) - original_u_at_fixed_volume(t - h, v))
                / (2.0 * h)
        };
        let richardson = (4.0 * difference(0.01) - difference(0.02)) / 3.0;
        let cv = state.heat_capacity_at_volume_j_kg_k();
        assert!(cv > 0.0);
        near(cv, richardson, 2e-6 * cv);
        let inverse = law
            .state_at_specific_internal_energy(state.specific_internal_energy_j_kg(), || false)
            .unwrap();
        near(inverse.state().temperature_k(), t, 1e-8);
        assert!(inverse.internal_energy_residual_j_kg().abs() <= 1e-6);
        assert!(inverse.pressure_iterations() > 0);
    }
}

#[test]
fn g1_dense_finite_inventory_heats_with_changing_pressure_and_phase_masses() {
    let v = mixed_volume(pair(640.0), 0.2);
    let law = law(v, 638.0, 642.0);
    let initial = at(&law, 639.0);
    let target = at(&law, 641.0);
    let mass = 0.25;
    let heat =
        mass * (target.specific_internal_energy_j_kg() - initial.specific_internal_energy_j_kg());
    let step = law
        .advance_heat(initial, mass, heat, 5e-5, || false)
        .unwrap();
    near(step.state.temperature_k(), 641.0, 1e-8);
    near(step.internal_energy_change_j, heat, 5e-5);
    assert!(step.state.pressure_pa() > initial.pressure_pa());
    assert_ne!(step.initial_phase_masses_kg, step.phase_masses_kg);
    near(step.phase_masses_kg.iter().sum(), mass, 2e-16);
    near(
        step.state.bulk_density_kg_m3(),
        initial.bulk_density_kg_m3(),
        2e-12,
    );
    assert_eq!(step.fixed_volume_m3, mass * v);
    assert_eq!(step.boundary_work_j, 0.0);
    near(
        step.enthalpy_change_j - step.pressure_enthalpy_change_j,
        heat,
        5e-5,
    );
    assert_eq!(
        step.endpoint_pressure_closure_defects_j,
        [
            mass * initial.pressure_closure_defect_j_kg(),
            mass * step.state.pressure_closure_defect_j_kg(),
        ]
    );
    let cooled = law
        .advance_heat(step.state, mass, -heat, 5e-5, || false)
        .unwrap();
    near(cooled.state.temperature_k(), initial.temperature_k(), 2e-8);
    near(
        cooled.state.vapor_quality().unwrap(),
        initial.vapor_quality().unwrap(),
        2e-10,
    );
}

#[test]
fn g3_rebinding_volume_preserves_mode_controls_and_actual_phase_density() {
    let v = mixed_volume(pair(640.0), 0.2);
    let original = law(v, 639.0, 641.0);
    let mut replacement = *original.controls();
    replacement.internal_energy_tolerance_j_kg = 2e-6;
    let changed = original
        .rebind_specific_volume(MassSpecificVolume::new(1.002 * v), replacement, || false)
        .unwrap();
    assert_eq!(changed.regional_coexistence_controls(), Some(pressure()));
    assert_eq!(changed.controls(), &replacement);
    let a = at(&original, 640.0);
    let b = at(&changed, 640.0);
    assert_eq!(a.regional_saturation_pair(), b.regional_saturation_pair());
    assert!(b.vapor_mass_fraction() > a.vapor_mass_fraction());
    near(
        a.bulk_density_kg_m3() / b.bulk_density_kg_m3(),
        1.002,
        2e-15,
    );
    let legacy = IsochoricWater::new(MassSpecificVolume::new(0.1), controls(350.0, 400.0), || {
        false
    })
    .unwrap();
    let rebound = legacy
        .rebind_specific_volume(MassSpecificVolume::new(0.1002), *legacy.controls(), || {
            false
        })
        .unwrap();
    assert_eq!(rebound.regional_coexistence_controls(), None);
}

#[test]
fn g0_exact_endpoints_remain_pairs_and_unsupported_paths_or_cancellation_refuse() {
    let p = pair(640.0);
    for (x, v) in [
        (0.0, p.liquid().specific_volume().value()),
        (1.0, p.vapor().specific_volume().value()),
    ] {
        let law = law(v, 639.0, 640.0);
        let state = at(&law, 640.0);
        assert_eq!(state.vapor_quality(), Some(x));
        assert_eq!(state.regional_saturation_pair(), Some(p));
        assert_eq!(state.single_phase_source(), None);
    }
    for v in [
        0.99 * p.liquid().specific_volume().value(),
        1.01 * p.vapor().specific_volume().value(),
    ] {
        assert!(matches!(
            IsochoricWater::new_regional_coexistence(
                MassSpecificVolume::new(v),
                controls(639.0, 640.0),
                pressure(),
                || false,
            ),
            Err(IsochoricWaterError::VolumeOutsideDomain { .. })
        ));
    }
    let v = mixed_volume(p, 0.2);
    for cfg in [
        controls(646.0, CRITICAL_TEMPERATURE_K),
        controls(623.0, 624.0),
        controls(623.15, 624.0),
    ] {
        assert!(matches!(
            IsochoricWater::new_regional_coexistence(
                MassSpecificVolume::new(v),
                cfg,
                pressure(),
                || false,
            ),
            Err(IsochoricWaterError::Invalid(_))
        ));
    }
    assert!(
        IsochoricWater::new(MassSpecificVolume::new(v), controls(639.0, 641.0), || false).is_err()
    );
    let mut bad_pressure = pressure();
    bad_pressure.maximum_iterations = 0;
    assert!(
        IsochoricWater::new_regional_coexistence(
            MassSpecificVolume::new(0.1),
            controls(350.0, 400.0),
            bad_pressure,
            || false,
        )
        .is_err()
    );
    let mut polls = 0;
    let cancelled = IsochoricWater::new_regional_coexistence(
        MassSpecificVolume::new(v),
        controls(639.0, 641.0),
        pressure(),
        || {
            polls += 1;
            polls >= 8
        },
    );
    assert_eq!(cancelled, Err(IsochoricWaterError::Cancelled));
    assert_eq!(polls, 8);
}
