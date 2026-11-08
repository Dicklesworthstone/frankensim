use super::*;

// Equal baths and identical cells keep the entire slab on the fusion plateau.
// Internal flux vanishes, leaving an independent scalar physical oracle at
// each external face. Series-mixture k/v is exactly 1/D(x), with D affine.
fn symmetric_plateau<'a>(
    caloric: &'a FreezingWaterThermalLaw,
    transport: &'a FreezingWaterConductivity,
    fraction: f64,
) -> FreezingSlab<'a> {
    let temperature = caloric.chart().melting().temperature().value();
    FreezingSlab::new(
        caloric,
        transport,
        FreezingSlabGeometry {
            thickness_m: 0.004,
            area_m2: 0.001,
            cells: 2,
        },
        caloric
            .chart()
            .state_at_liquid_fraction(fraction, || false)
            .unwrap(),
        [Some(FreezingSlabReservoir {
            temperature_k: temperature + 10.0,
            surface_resistance_m2_k_w: 0.0,
        }); 2],
        &[],
        budget(2),
        || Ok(()),
    )
    .unwrap()
}

fn resistance_coefficients(
    caloric: &FreezingWaterThermalLaw,
    transport: &FreezingWaterConductivity,
) -> [f64; 2] {
    let phase = |fraction| {
        let state = caloric
            .chart()
            .state_at_liquid_fraction(fraction, || false)
            .unwrap();
        let conductivity = transport
            .state_at_material_state(state, || false)
            .unwrap()
            .conductivity_w_m_k();
        state.specific_volume_m3_kg() / conductivity
    };
    [phase(0.0), phase(1.0)]
}

fn assert_plateau_balance(
    initial: &FreezingSlab<'_>,
    next: &FreezingSlab<'_>,
    audit: FreezingSlabAudit,
    initial_fraction: f64,
    expected_fraction: f64,
) {
    let latent = initial.caloric.chart().melting().latent_enthalpy_j_kg();
    let heat = initial.mass_kg() * latent * (expected_fraction - initial_fraction);
    for cell in next.profile(|| Ok(())).unwrap() {
        assert!((cell.material.liquid_mass_fraction() - expected_fraction).abs() < 1e-9);
        assert_eq!(
            next.thermal.cells()[&cell.cell_id].geometry,
            initial.thermal.cells()[&cell.cell_id].geometry,
        );
    }
    assert_eq!(next.mass_kg(), initial.mass_kg());
    assert!(
        next.inventory(|| Ok(())).unwrap().free_volume_m3()
            < initial.inventory(|| Ok(())).unwrap().free_volume_m3()
    );
    assert!((audit.thermal.external_heat_j - heat).abs() < 2e-6);
    assert!(audit.pressure_work_j < 0.0);
    assert!(audit.source_enthalpy_residual_j.abs() < 2e-6);
    assert!(audit.internal_energy_residual_j.abs() < 2e-6);
}

#[test]
fn implicit_expanding_layers_match_the_original_phase_volume_quadratic() {
    let (caloric, transport) = laws();
    let fraction = 0.4;
    let fixed = symmetric_plateau(&caloric, &transport, fraction);
    let initial = fixed.with_material_density_layers(|| Ok(())).unwrap();
    let [ice, liquid] = resistance_coefficients(&caloric, &transport);
    let resistance = (1.0 - fraction) * ice + fraction * liquid;
    let slope = liquid - ice;
    assert!(slope > 0.0);
    let dt = 5.0;
    let mass = initial.mass_kg() / 2.0;
    let area = initial.geometry.area_m2;
    let latent = caloric.chart().melting().latent_enthalpy_j_kg();
    let load = 2.0 * dt * area * area * 10.0 / (mass * mass * latent);
    let change = 2.0 * load / (resistance + (resistance * resistance + 4.0 * slope * load).sqrt());
    let expected = fraction + change;
    assert!(expected < 1.0);
    let control = implicit_controls(dt);
    let mut work = ImplicitConductionWorkspace::new(2).unwrap();
    let (next, audit) = initial
        .advance_implicit_to(dt, control, &mut work, || Ok(()))
        .unwrap();
    assert_plateau_balance(&initial, &next, audit.physical, fraction, expected);
    // Frozen reference length has a measurably different physical resistance.
    let (fixed_next, _) = fixed
        .advance_implicit_to(dt, control, &mut work, || Ok(()))
        .unwrap();
    let fixed_fraction = fixed_next.profile(|| Ok(())).unwrap()[0]
        .material
        .liquid_mass_fraction();
    assert!(expected - fixed_fraction > 1e-5);
    // Continuation retains the selected density geometry without rebinding it.
    let resistance = resistance + slope * change;
    let second = 2.0 * load / (resistance + (resistance * resistance + 4.0 * slope * load).sqrt());
    let (continued, second_audit) = next
        .advance_implicit_to(2.0 * dt, control, &mut work, || Ok(()))
        .unwrap();
    assert_plateau_balance(
        &next,
        &continued,
        second_audit.physical,
        expected,
        expected + second,
    );
}

#[test]
fn explicit_density_layers_refresh_each_step_and_cancel_without_changing_the_slab() {
    let (caloric, transport) = laws();
    let fraction = 0.4;
    let fixed = symmetric_plateau(&caloric, &transport, fraction);
    let original = fixed.profile(|| Ok(())).unwrap();
    assert!(matches!(
        fixed.with_material_density_layers(|| Err(DepositionError::Cancelled)),
        Err(FreezingSlabError::Thermal(DepositionError::Cancelled)),
    ));
    assert_eq!(fixed.profile(|| Ok(())).unwrap(), original);
    let initial = fixed.with_material_density_layers(|| Ok(())).unwrap();
    let [ice, liquid] = resistance_coefficients(&caloric, &transport);
    let mass = initial.mass_kg() / 2.0;
    let area = initial.geometry.area_m2;
    let latent = caloric.chart().melting().latent_enthalpy_j_kg();
    let dt = 0.01;
    let load = 2.0 * dt * area * area * 10.0 / (mass * mass * latent);
    let first = fraction + load / ((1.0 - fraction) * ice + fraction * liquid);
    let expected = first + load / ((1.0 - first) * ice + first * liquid);
    let (next, audit) = initial
        .advance_to(2.0 * dt, controls(dt), || Ok(()))
        .unwrap();
    assert_plateau_balance(&initial, &next, audit, fraction, expected);
}
