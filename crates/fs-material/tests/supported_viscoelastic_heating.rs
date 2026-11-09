use fs_blake3::ContentHash;
use fs_material::elastic::acoustic::viscoelastic::{
    ViscoelasticSlabSource,
    layers::{
        SupportedViscoelasticLayers, ViscoelasticBacking, ViscoelasticLayer,
        ViscoelasticLayerBudget, ViscoelasticLayers, ViscoelasticLayersPower,
    },
    thermal::{
        AcousticBackingThermalCell, AcousticHeatingBudget, AcousticThermalBin,
        SupportedViscoelasticHeating, ViscoelasticHeating,
    },
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase,
};
use fs_material::visco::GeneralizedMaxwell;
use fs_math::c64::C64;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation,
    Event, MaterialId, StepConfig, ThermalLaw, TimedEvent,
};

fn caloric(identity: u8, density: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([identity; 32]),
        SolidLiquidPhase::Solid,
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 300.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: density,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 400.0,
                temperature_k: 400.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: density,
            },
        ],
    )
    .unwrap()
}

fn cell(id: u64, material: u64, density: f64, volume: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: 5,
        body_id: 8,
        material_id: MaterialId(material),
        mass_kg: density * volume,
        specific_h_j_kg: 0.0,
        geometry: CellGeometry {
            start_m: [0.0; 3],
            end_m: [volume, 0.0, 0.0],
            volume_m3: volume,
            surface_area_m2: 1.0,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    }
}

fn config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1.0e-11,
    }
}

fn simulation<'a>(
    cells: Vec<DepositedCell>,
    materials: &[(MaterialId, &'a dyn ThermalLaw)],
) -> DepositionSimulation<'a> {
    let events = cells
        .into_iter()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        })
        .collect();
    DepositionSimulation::new(
        events,
        materials,
        DepositionBudget {
            max_cells: 8,
            max_links: 0,
            max_events: 8,
            max_substeps: 64,
        },
    )
    .unwrap()
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0
}

fn layer<'a>(
    bulk: &'a GeneralizedMaxwell,
    shear: &'a GeneralizedMaxwell,
    density: f64,
    identity: u8,
    thickness: f64,
) -> ViscoelasticLayer<'a> {
    ViscoelasticLayer {
        source: ViscoelasticSlabSource {
            bulk,
            shear,
            density_kg_m3: density,
            frequency_band_rad_s: (10.0, 100.0),
            maximum_absolute_strain: 0.1,
            source_identity: ContentHash([identity; 32]),
        },
        thickness_m: thickness,
    }
}

fn layer_budget(count: usize) -> ViscoelasticLayerBudget {
    ViscoelasticLayerBudget {
        maximum_layers: count,
        maximum_total_branches: 2 * count,
    }
}

fn heating_budget() -> AcousticHeatingBudget {
    AcousticHeatingBudget {
        maximum_bins: 8,
        maximum_power_segments: 8,
        maximum_source_evaluations: 64,
        energy_absolute_tolerance_j: 1.0e-11,
    }
}

fn bin(layer: usize, cell: u64, end: f64) -> AcousticThermalBin {
    AcousticThermalBin {
        layer,
        depths_m: [0.0, end],
        cell: CellId(cell),
        temperature_band_k: [300.0, 310.0],
    }
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5.0e-12 * expected.abs().max(1.0e-300),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn g1_stiff_backing_heat_survives_rounded_clamp_reflection_and_total_heat() {
    let bulk = [
        GeneralizedMaxwell::new(60.0, vec![(20.0, 0.05)]).unwrap(),
        GeneralizedMaxwell::new(130.0, vec![(15.0, 0.02)]).unwrap(),
    ];
    let shear = [
        GeneralizedMaxwell::new(30.0, vec![(15.0, 0.02)]).unwrap(),
        GeneralizedMaxwell::new(45.0, vec![(27.0, 0.03)]).unwrap(),
    ];
    let layers = [
        layer(&bulk[0], &shear[0], 4.0, 3, 0.2),
        layer(&bulk[1], &shear[1], 7.0, 8, 0.3),
    ];
    let fluid = FluidAcousticState::try_new(2.0, 10.0).unwrap();
    let stack = SupportedViscoelasticLayers::try_new(
        &layers,
        &fluid,
        ViscoelasticBacking::Impedance(C64::new(13.0e30, 7.0e30)),
        20.0,
        layer_budget(2),
        &mut || false,
    )
    .unwrap();
    let caloric = [caloric(41, 4.0), caloric(42, 7.0), caloric(43, 2.0)];
    let thermal = simulation(
        vec![cell(1, 7, 4.0, 0.2), cell(2, 8, 7.0, 0.3), cell(3, 9, 2.0, 0.5)],
        &[
            (MaterialId(7), &caloric[0]),
            (MaterialId(8), &caloric[1]),
            (MaterialId(9), &caloric[2]),
        ],
    );
    let bins = [bin(0, 1, 0.2), bin(1, 2, 0.3)];
    let base = SupportedViscoelasticHeating::try_new(
        &stack, C64::ONE, 1.0, &bins, &thermal, heating_budget(), || Ok(()),
    )
    .unwrap();
    // Independent physical transfer/positive-work result, retained even though
    // the exterior reflection rounds to the exact clamped answer.
    near(base.peak_power().backing_absorbed_power_w, 2.3433341559909816e-31);
    let receiver = AcousticBackingThermalCell {
        cell: CellId(3),
        temperature_band_k: [300.0, 310.0],
    };
    let heating = base
        .clone()
        .with_backing_receiver(&thermal, receiver, || Ok(()))
        .unwrap();
    assert_eq!(heating.backing_receiver(), Some(receiver));
    assert!(core::ptr::eq(heating.original_stack(), &stack));
    let (heated, accepted, _) = heating
        .advance_to(&thermal, 2.0, config(0.5), || Ok(()))
        .unwrap();
    let expected = 2.0 * base.peak_power().backing_absorbed_power_w;
    near(accepted.energy().backing_absorbed_j, expected);
    near(accepted.energy().backing_thermal_input_j, expected);
    assert_eq!(accepted.energy().exported_backing_j, 0.0);
    let stored = heated.cells()[&CellId(3)].mass_kg
        * heated.cells()[&CellId(3)].specific_h_j_kg;
    assert!(stored > 0.0);
    near(stored, expected);
    assert_eq!(heated.temperature_k(CellId(3)).unwrap(), 300.0);
    assert!(accepted.energy().bulk_absorbed_j > 0.0);
    assert!(accepted.energy().shear_absorbed_j > 0.0);

    let (exported_thermal, exported, _) = base
        .advance_to(&thermal, 2.0, config(0.5), || Ok(()))
        .unwrap();
    near(exported.energy().exported_backing_j, expected);
    assert_eq!(exported.energy().backing_thermal_input_j, 0.0);
    assert_eq!(exported_thermal.cells()[&CellId(3)].specific_h_j_kg, 0.0);
    assert_eq!(base.time_s(), 0.0);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn g1_sub_ulp_intrinsic_absorption_enters_native_enthalpy_without_temperature_fabrication() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![(1.0e-14, 0.05)]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![(1.0e-14, 0.02)]).unwrap();
    let layers = [layer(&bulk, &shear, 4.0, 9, 1.0e-8)];
    let fluid = FluidAcousticState::try_new(2.0, 10.0).unwrap();
    let stack = SupportedViscoelasticLayers::try_new(
        &layers,
        &fluid,
        ViscoelasticBacking::Clamped,
        20.0,
        layer_budget(1),
        &mut || false,
    )
    .unwrap();
    let caloric = caloric(44, 4.0);
    let thermal = simulation(
        vec![cell(1, 7, 4.0, 1.0e-8)],
        &[(MaterialId(7), &caloric)],
    );
    let heating = SupportedViscoelasticHeating::try_new(
        &stack,
        C64::ONE,
        1.0,
        &[bin(0, 1, 1.0e-8)],
        &thermal,
        heating_budget(),
        || Ok(()),
    )
    .unwrap();
    let power = heating.peak_power();
    let watts = power.bulk_absorbed_power_w + power.shear_absorbed_power_w;
    assert!(watts > 0.0 && watts < f64::EPSILON * power.incident_power_w);
    let (heated, accepted, audit) = heating
        .advance_to(&thermal, 1.0, config(1.0), || Ok(()))
        .unwrap();
    let cell = &heated.cells()[&CellId(1)];
    assert!(cell.specific_h_j_kg > 0.0);
    near(cell.mass_kg * cell.specific_h_j_kg, watts);
    near(accepted.energy().thermal_input_j, watts);
    near(audit.external_heat_j, watts);
    near(accepted.energy().bulk_absorbed_j, power.bulk_absorbed_power_w);
    near(accepted.energy().shear_absorbed_j, power.shear_absorbed_power_w);
    assert_eq!(accepted.energy().backing_absorbed_j, 0.0);
    assert_eq!(accepted.energy().backing_thermal_input_j, 0.0);
    assert_eq!(accepted.energy().exported_backing_j, 0.0);
    // The original caloric chart rounds this tiny temperature rise; enthalpy
    // still retains the absorbed energy, with no invented visible increment.
    assert_eq!(heated.temperature_k(CellId(1)).unwrap(), 300.0);
}

#[test]
fn g0_original_two_fluid_api_and_supported_receiver_source_guards_remain_binding() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![(20.0, 0.05)]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![(15.0, 0.02)]).unwrap();
    let layers = [layer(&bulk, &shear, 4.0, 3, 0.2)];
    let fluid = FluidAcousticState::try_new(2.0, 10.0).unwrap();
    let two_fluid = ViscoelasticLayers::try_new(
        &layers, [&fluid; 2], 20.0, layer_budget(1), &mut || false,
    )
    .unwrap();
    let supported = SupportedViscoelasticLayers::try_new(
        &layers,
        &fluid,
        ViscoelasticBacking::Clamped,
        20.0,
        layer_budget(1),
        &mut || false,
    )
    .unwrap();
    let caloric = caloric(45, 4.0);
    let thermal = simulation(vec![cell(1, 7, 4.0, 0.2)], &[(MaterialId(7), &caloric)]);
    let inputs = [C64::ONE, C64::new(0.2, 0.1)];
    let original: ViscoelasticHeating<'_, '_> = ViscoelasticHeating::try_new(
        &two_fluid,
        inputs,
        1.0,
        &[bin(0, 1, 0.2)],
        &thermal,
        heating_budget(),
        || Ok(()),
    )
    .unwrap();
    let copied = original.clone(); // The original physical stack is not Clone.
    assert!(core::ptr::eq(copied.original_stack(), &two_fluid));
    assert_eq!(copied.incident_peak_pressure_pa(), inputs);
    let native_power: ViscoelasticLayersPower = copied.peak_power();
    assert_eq!(
        native_power,
        two_fluid.drive(inputs, 1.0, &mut || false).unwrap().power()
    );
    let (_, accepted, _) = copied
        .advance_to(&thermal, 0.1, config(0.1), || Ok(()))
        .unwrap();
    assert_eq!(accepted.energy().backing_absorbed_j, 0.0);
    assert_eq!(accepted.energy().backing_thermal_input_j, 0.0);
    assert_eq!(accepted.energy().exported_backing_j, 0.0);

    let heating = SupportedViscoelasticHeating::try_new(
        &supported, C64::ONE, 1.0, &[bin(0, 1, 0.2)], &thermal, heating_budget(), || Ok(()),
    )
    .unwrap();
    for cell in [CellId(1), CellId(99)] {
        assert!(
            heating
                .clone()
                .with_backing_receiver(
                    &thermal,
                    AcousticBackingThermalCell {
                        cell,
                        temperature_band_k: [300.0, 310.0],
                    },
                    || Ok(()),
                )
                .is_err()
        );
    }
    let replacement = self::caloric(46, 4.0);
    let changed = simulation(vec![cell(1, 7, 4.0, 0.2)], &[(MaterialId(7), &replacement)]);
    assert!(heating.advance_to(&changed, 0.1, config(0.1), || Ok(())).is_err());
    assert_eq!(heating.time_s(), 0.0);
    assert_eq!(changed.time_s(), 0.0);

    let zero = [ViscoelasticLayer { thickness_m: 0.0, ..layers[0] }];
    let zero_supported = SupportedViscoelasticLayers::try_new(
        &zero,
        &fluid,
        ViscoelasticBacking::Impedance(C64::from_re(13.0)),
        20.0,
        layer_budget(1),
        &mut || false,
    )
    .unwrap();
    let zero_heating = SupportedViscoelasticHeating::try_new(
        &zero_supported, C64::ONE, 1.0, &[], &thermal, heating_budget(), || Ok(()),
    )
    .unwrap();
    assert_eq!(zero_heating.bin_count(), 0);
    assert_eq!(zero_heating.peak_power().bulk_absorbed_power_w, 0.0);
    assert_eq!(zero_heating.peak_power().shear_absorbed_power_w, 0.0);
    let zero_two_fluid = ViscoelasticLayers::try_new(
        &zero, [&fluid; 2], 20.0, layer_budget(1), &mut || false,
    )
    .unwrap();
    assert!(
        ViscoelasticHeating::try_new(
            &zero_two_fluid, inputs, 1.0, &[], &thermal, heating_budget(), || Ok(()),
        )
        .is_err()
    );
}
