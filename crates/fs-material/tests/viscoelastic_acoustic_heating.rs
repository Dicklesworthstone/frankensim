use fs_blake3::ContentHash;
use fs_material::MaterialError;
use fs_material::elastic::acoustic::viscoelastic::{
    ViscoelasticSlabSource,
    layers::{ViscoelasticLayer, ViscoelasticLayerBudget, ViscoelasticLayers},
    thermal::{
        AcousticHeatingBudget, AcousticHeatingError, AcousticPowerSegment, AcousticThermalBin,
        ViscoelasticHeating,
    },
};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::visco::GeneralizedMaxwell;
use fs_math::c64::C64;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionObserver, DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw,
    ThermalLink, TimedEvent,
};
use std::cell::Cell;

struct CountedFluid(Cell<usize>);
impl FluidAcoustics for CountedFluid {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        self.0.set(self.0.get() + 1);
        FluidAcousticState::try_new(2.0, 10.0)
    }
}

struct Sensible;
impl ThermalLaw for Sensible {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([41; 32])
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && (0.0..=4000.0).contains(&h) {
            Ok(300.0 + h / 4.0)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.25
    }
}
static SENSIBLE: Sensible = Sensible;

struct Plateau;
impl ThermalLaw for Plateau {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([42; 32])
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !h.is_finite() || !(0.0..=4000.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if h <= 4.0 {
            300.0 + h / 4.0
        } else if h <= 12.0 {
            301.0
        } else {
            301.0 + (h - 12.0) / 4.0
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.25
    }
}
static PLATEAU: Plateau = Plateau;

fn boundary(source_w: f64) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w,
    }
}
fn cell(id: u64, density: f64, area: f64, depths: [f64; 2], h: f64) -> DepositedCell {
    let volume = area * (depths[1] - depths[0]);
    DepositedCell {
        id: CellId(id),
        segment_id: 5,
        body_id: 8,
        material_id: MaterialId(7),
        mass_kg: density * volume,
        specific_h_j_kg: h,
        geometry: CellGeometry {
            start_m: [depths[0], 0.0, 0.0],
            end_m: [depths[1], 0.0, 0.0],
            volume_m3: volume,
            surface_area_m2: 10.0,
        },
        boundary: boundary(0.0),
    }
}
fn config(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-9,
    }
}
fn simulation(
    cells: Vec<DepositedCell>,
    law: &'static dyn ThermalLaw,
    extra: Vec<TimedEvent>,
) -> DepositionSimulation<'static> {
    let mut events: Vec<_> = cells
        .into_iter()
        .map(|c| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(c),
        })
        .collect();
    events.extend(extra);
    events.sort_by(|a, b| a.time_s.total_cmp(&b.time_s));
    DepositionSimulation::new(
        events,
        &[(MaterialId(7), law)],
        DepositionBudget {
            max_cells: 32,
            max_links: 32,
            max_events: 64,
            max_substeps: 4096,
        },
    )
    .unwrap()
    .advance_to(0.0, config(1.0), || Ok(()))
    .unwrap()
    .0
}
fn source<'a>(
    bulk: &'a GeneralizedMaxwell,
    shear: &'a GeneralizedMaxwell,
    density: f64,
    id: u8,
    thickness: f64,
) -> ViscoelasticLayer<'a> {
    ViscoelasticLayer {
        source: ViscoelasticSlabSource {
            bulk,
            shear,
            density_kg_m3: density,
            frequency_band_rad_s: (10.0, 100.0),
            maximum_absolute_strain: 0.2,
            source_identity: ContentHash([id; 32]),
        },
        thickness_m: thickness,
    }
}
fn stack<'a>(layers: &[ViscoelasticLayer<'a>], fluid: &'a CountedFluid) -> ViscoelasticLayers<'a> {
    ViscoelasticLayers::try_new(
        layers,
        [fluid; 2],
        20.0,
        ViscoelasticLayerBudget {
            maximum_layers: layers.len(),
            maximum_total_branches: 2 * layers.len(),
        },
        &mut || false,
    )
    .unwrap()
}
fn bin(layer: usize, id: u64, a: f64, b: f64) -> AcousticThermalBin {
    AcousticThermalBin {
        layer,
        depths_m: [a, b],
        cell: CellId(id),
        temperature_band_k: [300.0, 400.0],
    }
}
fn budget() -> AcousticHeatingBudget {
    AcousticHeatingBudget {
        maximum_bins: 32,
        maximum_power_segments: 8,
        maximum_source_evaluations: 4096,
        energy_absolute_tolerance_j: 1e-9,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn total(power: &fs_material::elastic::acoustic::viscoelastic::thermal::AcousticBinPower) -> f64 {
    power.bulk_absorbed_power_w + power.shear_absorbed_power_w
}

#[test]
fn g1_spatial_acoustic_work_heats_original_heterogeneous_cells_without_new_interfaces() {
    let bulk = [
        GeneralizedMaxwell::new(60.0, vec![(20.0, 0.05)]).unwrap(),
        GeneralizedMaxwell::new(130.0, vec![(15.0, 0.02)]).unwrap(),
    ];
    let shear = [
        GeneralizedMaxwell::new(30.0, vec![(15.0, 0.02)]).unwrap(),
        GeneralizedMaxwell::new(45.0, vec![(27.0, 0.03)]).unwrap(),
    ];
    let fluid = CountedFluid(Cell::new(0));
    let layers = [
        source(&bulk[0], &shear[0], 4.0, 3, 0.6),
        source(&bulk[1], &shear[1], 7.0, 8, 0.8),
    ];
    let stack = stack(&layers, &fluid);
    let mut bins = vec![
        bin(0, 11, 0.0, 0.1),
        bin(0, 12, 0.1, 0.6),
        bin(1, 21, 0.0, 0.3),
        bin(1, 22, 0.3, 0.8),
    ];
    let cells = bins
        .iter()
        .map(|b| {
            cell(
                b.cell.0,
                layers[b.layer].source.density_kg_m3,
                0.7,
                b.depths_m,
                0.0,
            )
        })
        .collect();
    let thermal = simulation(cells, &SENSIBLE, vec![]);
    bins.reverse(); // Caller order cannot change the original physical partition.
    let inputs = [C64::new(1.0, 0.2), C64::new(-0.1, 0.3)];
    let heating =
        ViscoelasticHeating::try_new(&stack, inputs, 0.7, &bins, &thermal, budget(), || Ok(()))
            .unwrap();
    assert!(core::ptr::eq(heating.original_stack(), &stack));
    assert_eq!(heating.incident_peak_pressure_pa(), inputs);
    assert_eq!(heating.area_m2(), 0.7);
    assert_eq!(heating.bin_count(), 4);
    assert!(heating.bin_power(4).is_none());
    assert_eq!(fluid.0.get(), 1);
    let (heated, accepted, audit) = heating
        .advance_to(&thermal, 2.0, config(0.125), || Ok(()))
        .unwrap();
    let mut accounted = 0.0;
    for i in 0..heating.bin_count() {
        let p = heating.bin_power(i).unwrap();
        let initial = &thermal.cells()[&p.binding.cell];
        let expected_j = total(p) * 2.0;
        near(
            heated.cells()[&p.binding.cell].specific_h_j_kg,
            expected_j / initial.mass_kg,
            2e-13,
        );
        assert_eq!(heated.cells()[&p.binding.cell].geometry, initial.geometry);
        assert_eq!(heated.cells()[&p.binding.cell].mass_kg, initial.mass_kg);
        assert_eq!(
            p.source_identity,
            layers[p.binding.layer].source.source_identity
        );
        accounted += expected_j;
    }
    near(accepted.energy().thermal_input_j, accounted, 2e-13);
    near(audit.external_heat_j, accounted, 2e-13);
    near(heated.energy_j(), accounted, 2e-13);
    near(
        accepted.energy().bulk_absorbed_j + accepted.energy().shear_absorbed_j,
        accounted,
        2e-13,
    );
    assert!(accepted.energy().balance_residual_j.abs() < 1e-11);
    assert!(accepted.energy().credit_residual_j.abs() < 1e-11);
    assert_eq!(accepted.source_evaluations(), audit.substeps);
    assert_eq!(fluid.0.get(), 1);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(heating.energy().thermal_input_j, 0.0);
}

#[test]
fn g1_g3_power_ramps_dark_gaps_native_events_and_restart_preserve_single_heat_credit() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![(20.0, 0.05)]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![(15.0, 0.02)]).unwrap();
    let fluid = CountedFluid(Cell::new(0));
    let stack = stack(&[source(&bulk, &shear, 4.0, 3, 1.0)], &fluid);
    let bins = [bin(0, 11, 0.0, 0.5), bin(0, 22, 0.5, 1.0)];
    let thermal = simulation(
        bins.iter()
            .map(|b| cell(b.cell.0, 4.0, 1.0, b.depths_m, 0.0))
            .collect(),
        &SENSIBLE,
        vec![
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 9,
                    a: CellId(11),
                    b: CellId(22),
                    area_m2: 1.0,
                    conductance_w_k: 0.2,
                }),
            },
            TimedEvent {
                time_s: 1.75,
                event: Event::Boundary(CellId(22), boundary(1.0)),
            },
        ],
    );
    let envelope = [
        AcousticPowerSegment {
            start_s: 0.5,
            end_s: 1.5,
            start_fraction: 0.0,
            end_fraction: 1.0,
        },
        AcousticPowerSegment {
            start_s: 2.0,
            end_s: 3.0,
            start_fraction: 1.0,
            end_fraction: 0.5,
        },
    ];
    let heating = ViscoelasticHeating::try_new(
        &stack,
        [C64::new(2.0, 0.1), C64::new(0.0, 0.2)],
        1.0,
        &bins,
        &thermal,
        budget(),
        || Ok(()),
    )
    .unwrap()
    .with_power_envelope(&thermal, &envelope, || Ok(()))
    .unwrap();
    assert!(heating.requires_step_interval());
    assert_eq!(heating.power_envelope().unwrap(), envelope);
    assert_eq!(heating.last_power_fraction(), 0.0);
    let peak = heating.peak_power();
    let expected_acoustic = 1.25 * (peak.bulk_absorbed_power_w + peak.shear_absorbed_power_w);
    let (whole, full, audit) = heating
        .advance_to(&thermal, 4.0, config(0.25), || Ok(()))
        .unwrap();
    near(full.energy().exposure_s, 1.25, 1e-14);
    near(full.energy().duration_s, 4.0, 1e-14);
    near(full.energy().thermal_input_j, expected_acoustic, 2e-12);
    near(audit.external_heat_j, expected_acoustic + 2.25, 2e-12);
    near(whole.energy_j(), audit.external_heat_j, 2e-12);
    assert_eq!(audit.events, 1);
    assert_eq!(full.last_power_fraction(), 0.0);
    for side in 0..2 {
        near(
            full.energy().incident_j[side],
            1.25 * peak.incident_power_w[side],
            2e-12,
        );
        near(
            full.energy().outgoing_j[side],
            1.25 * peak.outgoing_power_w[side],
            2e-12,
        );
    }
    let (middle, first, _) = heating
        .advance_to(&thermal, 2.0, config(0.25), || Ok(()))
        .unwrap();
    let (split, resumed, _) = first
        .advance_to(&middle, 4.0, config(0.25), || Ok(()))
        .unwrap();
    assert_eq!(split.cells(), whole.cells());
    assert_eq!(resumed.energy(), full.energy());
    assert_eq!(resumed.time_s(), 4.0);
    assert!(resumed.source_evaluations() < full.source_evaluations());
    assert_eq!(fluid.0.get(), 1);
}

#[test]
fn g1_original_latent_storage_absorbs_sound_without_a_fictitious_temperature_increment() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![(20.0, 0.05)]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![(15.0, 0.02)]).unwrap();
    let fluid = CountedFluid(Cell::new(0));
    let stack = stack(&[source(&bulk, &shear, 4.0, 3, 1.0)], &fluid);
    let thermal = simulation(vec![cell(11, 4.0, 1.0, [0.0, 1.0], 4.0)], &PLATEAU, vec![]);
    let heating = ViscoelasticHeating::try_new(
        &stack,
        [C64::new(2.0, 0.0), C64::ZERO],
        1.0,
        &[bin(0, 11, 0.0, 1.0)],
        &thermal,
        budget(),
        || Ok(()),
    )
    .unwrap();
    let q = total(heating.bin_power(0).unwrap());
    let duration = 4.0 * thermal.cells()[&CellId(11)].mass_kg / q;
    let (heated, accepted, audit) = heating
        .advance_to(&thermal, duration, config(duration), || Ok(()))
        .unwrap();
    near(heated.cells()[&CellId(11)].specific_h_j_kg, 8.0, 1e-12);
    assert_eq!(thermal.temperature_k(CellId(11)).unwrap(), 301.0);
    assert_eq!(heated.temperature_k(CellId(11)).unwrap(), 301.0);
    near(accepted.energy().thermal_input_j, 16.0, 1e-11);
    near(
        heated.energy_j() - thermal.energy_j(),
        audit.external_heat_j,
        1e-11,
    );
}

#[test]
fn g1_spatial_acoustic_heating_and_native_conduction_converge_to_independent_two_cell_solution() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![(20.0, 0.05)]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![(15.0, 0.02)]).unwrap();
    let fluid = CountedFluid(Cell::new(0));
    let stack = stack(&[source(&bulk, &shear, 4.0, 3, 1.0)], &fluid);
    let bins = [bin(0, 11, 0.0, 0.5), bin(0, 22, 0.5, 1.0)];
    let g = 4.0;
    let thermal = simulation(
        bins.iter()
            .map(|b| cell(b.cell.0, 4.0, 1.0, b.depths_m, 0.0))
            .collect(),
        &SENSIBLE,
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 9,
                a: CellId(11),
                b: CellId(22),
                area_m2: 1.0,
                conductance_w_k: g,
            }),
        }],
    );
    let heating = ViscoelasticHeating::try_new(
        &stack,
        [C64::new(4.0, 0.0), C64::ZERO],
        1.0,
        &bins,
        &thermal,
        budget(),
        || Ok(()),
    )
    .unwrap();
    let q = [
        total(heating.bin_power(0).unwrap()),
        total(heating.bin_power(1).unwrap()),
    ];
    assert!((q[0] - q[1]).abs() > 1e-5);
    let capacity = 4.0 * thermal.cells()[&CellId(11)].mass_kg;
    let end = 1.0;
    let mean = 300.0 + (q[0] + q[1]) * end / (2.0 * capacity);
    let difference = (q[0] - q[1]) / (2.0 * g) * -(-2.0 * g * end / capacity).exp_m1();
    let exact = [mean + 0.5 * difference, mean - 0.5 * difference];
    let mut errors = Vec::new();
    for dt in [0.1, 0.05, 0.025] {
        let (heated, accepted, audit) = heating
            .advance_to(&thermal, end, config(dt), || Ok(()))
            .unwrap();
        let error = [CellId(11), CellId(22)]
            .into_iter()
            .zip(exact)
            .map(|(id, t)| (heated.temperature_k(id).unwrap() - t).abs())
            .fold(0.0_f64, f64::max);
        errors.push(error);
        near(audit.external_heat_j, (q[0] + q[1]) * end, 1e-11);
        near(
            accepted.energy().thermal_input_j,
            audit.external_heat_j,
            1e-11,
        );
    }
    assert!((1.8..2.2).contains(&(errors[0] / errors[1])), "{errors:?}");
    assert!((1.8..2.2).contains(&(errors[1] / errors[2])), "{errors:?}");
}

#[test]
fn g0_physical_partition_inventory_and_original_caloric_bindings_refuse_mismatches() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![(20.0, 0.05)]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![(15.0, 0.02)]).unwrap();
    let fluid = CountedFluid(Cell::new(0));
    let stack = stack(&[source(&bulk, &shear, 4.0, 3, 1.0)], &fluid);
    let bins = [bin(0, 11, 0.0, 0.5), bin(0, 22, 0.5, 1.0)];
    let cells: Vec<_> = bins
        .iter()
        .map(|b| cell(b.cell.0, 4.0, 1.0, b.depths_m, 0.0))
        .collect();
    let thermal = simulation(cells.clone(), &SENSIBLE, vec![]);
    let inputs = [C64::ONE, C64::ZERO];
    let create = |bins: &[AcousticThermalBin]| {
        ViscoelasticHeating::try_new(&stack, inputs, 1.0, bins, &thermal, budget(), || Ok(()))
    };
    assert!(create(&[]).is_err());
    assert!(create(&bins[..1]).is_err());
    for changed in [
        AcousticThermalBin {
            depths_m: [0.6, 1.0],
            ..bins[1]
        },
        AcousticThermalBin {
            depths_m: [0.4, 1.0],
            ..bins[1]
        },
        AcousticThermalBin {
            depths_m: [0.5, 0.5],
            ..bins[1]
        },
        AcousticThermalBin {
            cell: CellId(11),
            ..bins[1]
        },
        AcousticThermalBin {
            layer: 1,
            ..bins[1]
        },
        AcousticThermalBin {
            temperature_band_k: [301.0, 400.0],
            ..bins[1]
        },
    ] {
        assert!(create(&[bins[0], changed]).is_err());
    }
    let heating = create(&bins).unwrap();
    let other_caloric = simulation(cells.clone(), &PLATEAU, vec![]);
    assert!(
        heating
            .advance_to(&other_caloric, 1.0, config(0.5), || Ok(()))
            .is_err()
    );
    let mut moved_mass = cells.clone();
    moved_mass[0].mass_kg *= 1.01;
    let wrong_inventory = simulation(moved_mass, &SENSIBLE, vec![]);
    assert!(
        ViscoelasticHeating::try_new(
            &stack,
            inputs,
            1.0,
            &bins,
            &wrong_inventory,
            budget(),
            || Ok(()),
        )
        .is_err()
    );
    assert!(
        heating
            .advance_to(&wrong_inventory, 1.0, config(0.5), || Ok(()))
            .is_err()
    );
    assert!(
        ViscoelasticHeating::try_new(
            &stack,
            inputs,
            1.0,
            &bins,
            &thermal,
            AcousticHeatingBudget {
                maximum_bins: 1,
                ..budget()
            },
            || Ok(()),
        )
        .is_err()
    );
    let invalid_envelopes = [
        AcousticPowerSegment {
            start_s: 0.0,
            end_s: 1.0,
            start_fraction: -1.0,
            end_fraction: 1.0,
        },
        AcousticPowerSegment {
            start_s: 0.0,
            end_s: 1.0,
            start_fraction: 1.0,
            end_fraction: 1.01,
        },
        AcousticPowerSegment {
            start_s: 1.0,
            end_s: 1.0,
            start_fraction: 0.0,
            end_fraction: 1.0,
        },
    ];
    for segment in invalid_envelopes {
        assert!(
            heating
                .clone()
                .with_power_envelope(&thermal, &[segment], || Ok(()))
                .is_err()
        );
    }
}

#[test]
fn g4_endpoint_domains_work_limits_and_cancellation_leave_native_history_and_heat_unpublished() {
    let bulk = GeneralizedMaxwell::new(60.0, vec![(20.0, 0.05)]).unwrap();
    let shear = GeneralizedMaxwell::new(30.0, vec![(15.0, 0.02)]).unwrap();
    let fluid = CountedFluid(Cell::new(0));
    let stack = stack(&[source(&bulk, &shear, 4.0, 3, 1.0)], &fluid);
    let thermal = simulation(vec![cell(11, 4.0, 1.0, [0.0, 1.0], 0.0)], &SENSIBLE, vec![]);
    let binding = bin(0, 11, 0.0, 1.0);
    let limited_temperature = AcousticThermalBin {
        temperature_band_k: [300.0, 300.000001],
        ..binding
    };
    let heating = ViscoelasticHeating::try_new(
        &stack,
        [C64::new(4.0, 0.0), C64::ZERO],
        1.0,
        &[limited_temperature],
        &thermal,
        budget(),
        || Ok(()),
    )
    .unwrap();
    assert!(matches!(
        heating.advance_to(&thermal, 1.0, config(0.5), || Ok(())),
        Err(AcousticHeatingError::TemperatureSupport { .. })
    ));
    let bounded = ViscoelasticHeating::try_new(
        &stack,
        [C64::new(4.0, 0.0), C64::ZERO],
        1.0,
        &[binding],
        &thermal,
        AcousticHeatingBudget {
            maximum_source_evaluations: 1,
            ..budget()
        },
        || Ok(()),
    )
    .unwrap();
    // The first candidate interval succeeds before the second exhausts work.
    assert!(matches!(
        bounded.advance_to(&thermal, 1.0, config(0.5), || Ok(())),
        Err(AcousticHeatingError::Budget)
    ));
    let sound = ViscoelasticHeating::try_new(
        &stack,
        [C64::new(4.0, 0.0), C64::ZERO],
        1.0,
        &[binding],
        &thermal,
        budget(),
        || Ok(()),
    )
    .unwrap();
    let mut calls = 0;
    assert!(matches!(
        sound.advance_to(&thermal, 10.0, config(0.1), || {
            calls += 1;
            if calls == 50 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(AcousticHeatingError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(calls, 50);
    let envelope = [AcousticPowerSegment {
        start_s: 0.0,
        end_s: 1.0,
        start_fraction: 0.0,
        end_fraction: 1.0,
    }];
    let mut scheduled = sound
        .clone()
        .with_power_envelope(&thermal, &envelope, || Ok(()))
        .unwrap();
    assert!(
        scheduled
            .step_sources(&thermal, 0.0, 1.5, 1.5, &mut || Ok(()))
            .is_err()
    );
    assert_eq!(scheduled.energy().thermal_input_j, 0.0);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.energy_j(), 0.0);
    assert_eq!(heating.energy().thermal_input_j, 0.0);
    assert_eq!(bounded.energy().thermal_input_j, 0.0);
    assert_eq!(sound.energy().thermal_input_j, 0.0);
    let (next, accepted, _) = sound
        .advance_to(&thermal, 0.1, config(0.1), || Ok(()))
        .unwrap();
    assert!(next.energy_j() > 0.0 && accepted.energy().thermal_input_j > 0.0);
    assert!(matches!(
        sound.advance_to(&next, 0.2, config(0.1), || Ok(())),
        Err(AcousticHeatingError::ClockMismatch { .. })
    ));
    assert!(
        accepted
            .with_power_envelope(&next, &envelope, || Ok(()))
            .is_err()
    );
}
