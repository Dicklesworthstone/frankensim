//! Original oscillator-backed layer heat on the native accepted thermal clock.

use super::close;
use fs_blake3::ContentHash;
use fs_couple::elastic_reflection::OscillatorBackedViscoelasticLayers;
use fs_material::elastic::acoustic::viscoelastic::{
    ViscoelasticSlabSource,
    layers::{ViscoelasticLayer, ViscoelasticLayerBudget},
    thermal::{
        AcousticBackingThermalCell, AcousticHeatingBudget, AcousticHeatingError,
        AcousticPowerSegment, AcousticThermalBin,
    },
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_material::visco::GeneralizedMaxwell;
use fs_math::c64::C64;
use fs_phs::WallPin;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, TimedEvent,
};

fn pin(resistance: f64) -> WallPin {
    WallPin {
        surface_density: 1.0,
        stiffness_per_area: 260.0,
        resistance,
    }
}

fn spectra() -> ([GeneralizedMaxwell; 2], [GeneralizedMaxwell; 2]) {
    (
        [
            GeneralizedMaxwell::new(60.0, vec![(20.0, 0.05)]).unwrap(),
            GeneralizedMaxwell::new(130.0, vec![(15.0, 0.02)]).unwrap(),
        ],
        [
            GeneralizedMaxwell::new(30.0, vec![(15.0, 0.02)]).unwrap(),
            GeneralizedMaxwell::new(45.0, vec![(27.0, 0.03)]).unwrap(),
        ],
    )
}

fn mechanics<'a>(
    bulk: &'a [GeneralizedMaxwell; 2],
    shear: &'a [GeneralizedMaxwell; 2],
    fluid: &'a FluidAcousticState,
) -> OscillatorBackedViscoelasticLayers<'a> {
    let layers = std::array::from_fn::<_, 2, _>(|i| ViscoelasticLayer {
        source: ViscoelasticSlabSource {
            bulk: &bulk[i],
            shear: &shear[i],
            density_kg_m3: [4.0, 7.0][i],
            frequency_band_rad_s: (10.0, 100.0),
            maximum_absolute_strain: 0.2,
            source_identity: ContentHash([[3; 32], [8; 32]][i]),
        },
        thickness_m: [0.2, 0.3][i],
    });
    OscillatorBackedViscoelasticLayers::prepare(
        &layers,
        fluid,
        pin(13.0),
        20.0,
        ViscoelasticLayerBudget { maximum_layers: 2, maximum_total_branches: 4 },
    )
    .unwrap()
}

// Original native caloric charts, with explicit source identity, density,
// energy datum and all-solid phase; no replacement test ThermalLaw.
fn caloric(density: f64, identity: u8) -> EquilibriumEnthalpyPhaseCurve {
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

fn boundary(source_w: f64) -> Boundary {
    Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0, reservoir: None, source_w }
}

fn cell(id: u64, material: u64, density: f64, area: f64, depths: [f64; 2]) -> DepositedCell {
    let volume = area * (depths[1] - depths[0]);
    DepositedCell {
        id: CellId(id),
        segment_id: material,
        body_id: 8,
        material_id: MaterialId(material),
        mass_kg: density * volume,
        specific_h_j_kg: 0.0,
        geometry: CellGeometry {
            start_m: [0.0, 0.0, depths[0]],
            end_m: [0.0, 0.0, depths[1]],
            volume_m3: volume,
            surface_area_m2: 10.0,
        },
        boundary: boundary(0.0),
    }
}

fn cells() -> Vec<DepositedCell> {
    vec![
        cell(11, 1, 4.0, 1.0, [0.0, 0.1]),
        cell(12, 1, 4.0, 1.0, [0.1, 0.2]),
        cell(21, 2, 7.0, 1.0, [0.2, 0.35]),
        cell(22, 2, 7.0, 1.0, [0.35, 0.5]),
    ]
}

fn bins() -> [AcousticThermalBin; 4] {
    [(0, 11, [0.0, 0.1]), (0, 12, [0.1, 0.2]),
     (1, 21, [0.0, 0.15]), (1, 22, [0.15, 0.3])]
        .map(|(layer, id, depths_m)| AcousticThermalBin {
            layer,
            depths_m,
            cell: CellId(id),
            temperature_band_k: [300.0, 350.0],
        })
}

fn budget() -> AcousticHeatingBudget {
    AcousticHeatingBudget {
        maximum_bins: 16,
        maximum_power_segments: 8,
        maximum_source_evaluations: 4096,
        energy_absolute_tolerance_j: 1e-9,
    }
}

fn config(dt: f64) -> StepConfig {
    StepConfig { max_dt_s: dt, energy_tolerance_j: 1e-9 }
}

fn simulation<'a>(
    cells: Vec<DepositedCell>,
    laws: &'a [EquilibriumEnthalpyPhaseCurve],
    extra: Vec<TimedEvent>,
) -> DepositionSimulation<'a> {
    let mut events: Vec<_> = cells.into_iter().map(|cell| TimedEvent {
        time_s: 0.0,
        event: Event::Birth(cell),
    }).collect();
    events.extend(extra);
    events.sort_by(|a, b| a.time_s.total_cmp(&b.time_s));
    let materials: Vec<_> = laws.iter().enumerate().map(|(i, law)| {
        (MaterialId(i as u64 + 1), law as &dyn ThermalLaw)
    }).collect();
    DepositionSimulation::new(
        events,
        &materials,
        DepositionBudget { max_cells: 16, max_links: 16, max_events: 32, max_substeps: 4096 },
    )
    .unwrap()
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0
}

#[test]
fn g1_original_spatial_layer_heat_reaches_native_cells_while_backing_work_is_exported() {
    let (bulk, shear) = spectra();
    let fluid = FluidAcousticState::try_new(2.0, 10.0).unwrap();
    let original = mechanics(&bulk, &shear, &fluid);
    let laws = [caloric(4.0, 91), caloric(7.0, 92)];
    let initial = simulation(cells(), &laws, vec![]);
    let mut partition = bins();
    partition.reverse();
    let source = original.prepare_heating(
        C64::ONE, 1.0, &partition, &initial, budget(), || Ok(()),
    ).unwrap();
    assert!(core::ptr::eq(source.mechanics(), &original));
    assert!(core::ptr::eq(source.heating().original_stack(), original.layers()));
    assert_eq!(source.heating().area_m2(), 1.0);
    assert_eq!(source.heating().bin_count(), 4);

    // Independent physical p/v transfer and exact exponential strain integrals.
    // Equal-depth cells do not get equal shares of a standing-wave heat field.
    let expected = [
        [0.001558058059547432, 0.0010745227996878841],
        [0.0018173013797158368, 0.0012533112963557494],
        [0.00035255936035097075, 0.00108256462413651],
        [0.00014418430028242965, 0.0004427306161613427],
    ];
    let (heated, accepted, audit) = source.advance_to(
        &initial, 2.0, config(0.125), || Ok(()),
    ).unwrap();
    for (i, expected) in expected.into_iter().enumerate() {
        let native = accepted.heating().bin_power(i).unwrap();
        close(native.bulk_absorbed_power_w, expected[0], 3e-14);
        close(native.shear_absorbed_power_w, expected[1], 3e-14);
        let before = &initial.cells()[&native.binding.cell];
        let after = &heated.cells()[&native.binding.cell];
        let heat = 2.0 * (expected[0] + expected[1]);
        close(after.specific_h_j_kg, heat / before.mass_kg, 3e-13);
        close(heated.temperature_k(after.id).unwrap(), 300.0 + heat / (4.0 * before.mass_kg), 2e-12);
        assert_eq!(before.geometry, after.geometry);
        assert_eq!(before.mass_kg, after.mass_kg);
        assert_eq!(before.material_id, after.material_id);
        assert_eq!(native.source_identity, ContentHash(if native.binding.layer == 0 { [3; 32] } else { [8; 32] }));
    }
    let e = accepted.heating().energy();
    close(e.incident_j[0], 0.05, 3e-14);
    close(e.outgoing_j[0], 2.0 * 0.006691669902970649, 3e-14);
    close(e.backing_absorbed_j, 2.0 * 0.010583097660791197, 3e-14);
    close(e.exported_backing_j, e.backing_absorbed_j, 3e-15);
    assert_eq!(e.backing_thermal_input_j, 0.0);
    close(e.thermal_input_j, 2.0 * 0.007725232436238154, 3e-14);
    close(e.bulk_absorbed_j + e.shear_absorbed_j, e.thermal_input_j, 3e-14);
    close(audit.external_heat_j, e.thermal_input_j, 3e-14);
    close(heated.energy_j(), e.thermal_input_j, 3e-14);
    close(e.balance_residual_j, 0.0, 3e-14);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(source.heating().time_s(), 0.0);
}

#[test]
fn g1_explicit_backing_receiver_heats_its_own_mass_without_thermalizing_pin_storage() {
    let fluid = FluidAcousticState::try_new(2.0, 10.0).unwrap();
    let laws = [caloric(9.0, 93)];
    let input = C64::new(1.0, 0.2);
    let area = 0.7;
    for resistance in [13.0, 0.0] {
        let original = OscillatorBackedViscoelasticLayers::prepare(
            &[], &fluid, pin(resistance), 20.0,
            ViscoelasticLayerBudget { maximum_layers: 0, maximum_total_branches: 0 },
        ).unwrap();
        let initial = simulation(vec![cell(90, 1, 9.0, area, [0.0, 1.0])], &laws, vec![]);
        assert_ne!(initial.cells()[&CellId(90)].mass_kg, original.backing().surface_density * area);
        let source = original.prepare_heating(
            input, area, &[], &initial, budget(), || Ok(()),
        ).unwrap().with_backing_receiver(
            &initial,
            AcousticBackingThermalCell { cell: CellId(90), temperature_band_k: [300.0, 350.0] },
            || Ok(()),
        ).unwrap().with_power_envelope(
            &initial,
            &[
                AcousticPowerSegment { start_s: 0.0, end_s: 1.0, start_fraction: 0.0, end_fraction: 1.0 },
                AcousticPowerSegment { start_s: 1.0, end_s: 1.5, start_fraction: 1.0, end_fraction: 1.0 },
            ],
            || Ok(()),
        ).unwrap();
        let peak = source.peak_response(&mut || false).unwrap();
        let expected_v = input.scale(2.0) / C64::new(20.0 + resistance, 7.0);
        close((peak.layers.support_normal_velocity_m_s() - expected_v).abs(), 0.0, 3e-15);
        let expected_storage = area * 0.25 * expected_v.norm_sq() * (1.0 + 260.0 / 400.0);
        close(peak.mean_support_stored_energy_j, expected_storage, 3e-16);
        assert!(expected_storage > 0.0);
        let expected_heat = 0.5 * area * resistance * expected_v.norm_sq();
        let (heated, accepted, audit) = source.advance_to(
            &initial, 2.0, config(0.37), || Ok(()),
        ).unwrap();
        let e = accepted.heating().energy();
        close(e.exposure_s, 1.0, 3e-15);
        close(e.backing_absorbed_j, expected_heat, 3e-14);
        close(e.backing_thermal_input_j, expected_heat, 3e-14);
        close(e.thermal_input_j, expected_heat, 3e-14);
        assert_eq!(e.exported_backing_j, 0.0);
        assert_eq!(e.bulk_absorbed_j, 0.0);
        assert_eq!(e.shear_absorbed_j, 0.0);
        close(audit.external_heat_j, expected_heat, 3e-14);
        close(heated.energy_j(), expected_heat, 3e-14);
        assert_eq!(heated.cells()[&CellId(90)].mass_kg, initial.cells()[&CellId(90)].mass_kg);
        assert_eq!(accepted.heating().last_power_fraction(), 0.0);
        // Peak diagnostics retain the authored carrier, even after a dark tail.
        close(accepted.peak_response(&mut || false).unwrap().mean_support_stored_energy_j,
              expected_storage, 3e-16);
        if resistance == 0.0 {
            assert_eq!(heated.energy_j(), 0.0);
            assert_eq!(e.backing_absorbed_j, 0.0);
        }
    }
}

#[test]
fn g3_native_event_cuts_and_split_advances_preserve_exposure_and_original_sources() {
    let (bulk, shear) = spectra();
    let fluid = FluidAcousticState::try_new(2.0, 10.0).unwrap();
    let original = mechanics(&bulk, &shear, &fluid);
    let laws = [caloric(4.0, 91), caloric(7.0, 92)];
    let initial = simulation(cells(), &laws, vec![
        TimedEvent { time_s: 0.35, event: Event::Boundary(CellId(11), boundary(0.003)) },
        TimedEvent { time_s: 0.8, event: Event::Boundary(CellId(11), boundary(0.0)) },
    ]);
    let source = original.prepare_heating(
        C64::ONE, 1.0, &bins(), &initial, budget(), || Ok(()),
    ).unwrap().with_power_envelope(
        &initial,
        &[AcousticPowerSegment { start_s: 0.1, end_s: 0.9, start_fraction: 0.2, end_fraction: 0.8 }],
        || Ok(()),
    ).unwrap();
    let (whole, whole_heat, whole_audit) = source.advance_to(
        &initial, 1.2, config(0.27), || Ok(()),
    ).unwrap();
    let (middle, middle_heat, first_audit) = source.advance_to(
        &initial, 0.43, config(0.27), || Ok(()),
    ).unwrap();
    let (split, split_heat, second_audit) = middle_heat.advance_to(
        &middle, 1.2, config(0.27), || Ok(()),
    ).unwrap();
    for id in [11, 12, 21, 22].map(CellId) {
        close(whole.cells()[&id].specific_h_j_kg, split.cells()[&id].specific_h_j_kg, 3e-13);
    }
    close(whole_heat.heating().energy().exposure_s, 0.4, 3e-15);
    close(split_heat.heating().energy().exposure_s, 0.4, 3e-15);
    let intrinsic = 0.007725232436238154 * 0.4;
    let authored = 0.003 * (0.8 - 0.35);
    close(whole_heat.heating().energy().thermal_input_j, intrinsic, 3e-14);
    close(whole_audit.external_heat_j, intrinsic + authored, 3e-14);
    close(first_audit.external_heat_j + second_audit.external_heat_j,
          whole_audit.external_heat_j, 3e-14);
    close(whole.energy_j(), intrinsic + authored, 3e-14);
    close(split.energy_j(), whole.energy_j(), 3e-14);
    close(split_heat.heating().energy().exported_backing_j,
          0.010583097660791197 * 0.4, 3e-14);
    assert_eq!(whole_audit.events, 2);
    assert_eq!(whole.cells()[&CellId(11)].boundary.source_w, 0.0);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(source.heating().time_s(), 0.0);
}

#[test]
fn g0_g4_native_source_receiver_temperature_budget_and_cancellation_refusals_are_atomic() {
    let (bulk, shear) = spectra();
    let fluid = FluidAcousticState::try_new(2.0, 10.0).unwrap();
    let original = mechanics(&bulk, &shear, &fluid);
    let laws = [caloric(4.0, 91), caloric(7.0, 92)];
    let initial = simulation(cells(), &laws, vec![]);
    assert!(original.prepare_heating(
        C64::ONE, 2.0, &bins(), &initial, budget(), || Ok(()),
    ).is_err());
    let mut wrong_mass = cells();
    wrong_mass[0].mass_kg *= 1.1;
    let wrong_mass = simulation(wrong_mass, &laws, vec![]);
    assert!(original.prepare_heating(
        C64::ONE, 1.0, &bins(), &wrong_mass, budget(), || Ok(()),
    ).is_err());
    let source = original.prepare_heating(
        C64::ONE, 1.0, &bins(), &initial, budget(), || Ok(()),
    ).unwrap();
    assert!(original.prepare_heating(
        C64::ONE, 1.0, &bins(), &initial, budget(), || Ok(()),
    ).unwrap().with_backing_receiver(
        &initial,
        AcousticBackingThermalCell { cell: CellId(11), temperature_band_k: [300.0, 350.0] },
        || Ok(()),
    ).is_err());
    let changed_laws = [caloric(4.0, 94), caloric(7.0, 92)];
    let changed_source = simulation(cells(), &changed_laws, vec![]);
    assert!(source.advance_to(&changed_source, 1.0, config(0.1), || Ok(())).is_err());

    let mut limited_bins = bins();
    for bin in &mut limited_bins {
        bin.temperature_band_k = [300.0, 300.0];
    }
    let limited = original.prepare_heating(
        C64::ONE, 1.0, &limited_bins, &initial, budget(), || Ok(()),
    ).unwrap();
    assert!(matches!(
        limited.advance_to(&initial, 1.0, config(0.1), || Ok(())),
        Err(AcousticHeatingError::TemperatureSupport { .. })
    ));
    let mut calls = 0;
    assert!(matches!(
        source.advance_to(&initial, 1.0, config(0.01), || {
            calls += 1;
            if calls > 100 { Err(DepositionError::Cancelled) } else { Ok(()) }
        }),
        Err(AcousticHeatingError::Thermal(DepositionError::Cancelled))
    ));
    assert!(calls > 100);
    let mut small_budget = budget();
    small_budget.maximum_source_evaluations = 1;
    let budgeted = original.prepare_heating(
        C64::ONE, 1.0, &bins(), &initial, small_budget, || Ok(()),
    ).unwrap();
    assert!(matches!(
        budgeted.advance_to(&initial, 1.0, config(0.1), || Ok(())),
        Err(AcousticHeatingError::Budget)
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.energy_j(), 0.0);
    assert_eq!(source.heating().time_s(), 0.0);
    assert_eq!(source.heating().energy().thermal_input_j, 0.0);
    assert_eq!(limited.heating().time_s(), 0.0);
    assert_eq!(budgeted.heating().time_s(), 0.0);
    let (valid, accepted, _) = source.advance_to(
        &initial, 0.1, config(0.1), || Ok(()),
    ).unwrap();
    assert!(valid.energy_j() > 0.0);
    assert_eq!(accepted.heating().time_s(), 0.1);
    assert!(core::ptr::eq(accepted.mechanics(), &original));
}
