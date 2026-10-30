//! G0/G3: original named gas laws reach coupled flow, ducts and source controls.
//! These bounded, unsaturated mixtures are numerical fixtures, not breath data.
use super::*;
use fs_material::gas_species::resolve_gas_mixture;
use fs_qty::semantic::CompositionBasis;

const MIXTURE: &str =
    "mixture 300 100000 mole 5 N2 0.74 O2 0.2 CO2 0.04 He 0.015 H2O 0.005";

fn input(declaration: &str) -> String {
    let record = read_gas_record(&mut declaration.split_ascii_whitespace(), None, None).unwrap();
    let gas = *record.state();
    let dx = (3.0 * gas.specific_gas_constant * gas.temperature).sqrt() / f64::from(RATE);
    let radius = (9.0 * dx * dx / core::f64::consts::PI).sqrt();
    fixture(gas, declaration)
        .replace("fluid_limits 108 250 600", "fluid_limits 108 290 320")
        .replace("coupling 64 6912", "coupling 128 13824")
        .replace("initial rest", "reservoir_aperture 1 low reference 1\nreservoir_cell 5 0 1\ninitial rest")
        .replace("observation inlet-pressure", &format!(
            "upstream_dg dg_tube 0.03 {radius:.17e} matched 1 2 32 50 2000 0.01 1e-11 512\nobservation inlet-pressure"))
        .replace("boundary 1\nhold 32 -0.1",
            "boundary reservoir 2\nhold 16 100002 300\nhold 32 100003 300.01")
}

#[test]
fn g0_g3_named_records_share_one_material_across_flow_ducts_and_reservoir_replay() {
    for declaration in [MIXTURE.to_owned(), MIXTURE.replace("mole", "mass")] {
        let original = read_gas_record(&mut declaration.split_ascii_whitespace(), None, None).unwrap();
        let GasRecord::Mixture(_, original_components) = original else {
            panic!("shared reader must retain named composition");
        };
        let text = input(&declaration);
        let file = read(&text, 8);
        let [GasRecord::Mixture(gas, components)] = file.gas_records() else {
            panic!("native playback must retain its named source record");
        };
        assert_eq!(components, &original_components);
        let system = file.renderer().system();
        let zone = system.zone();
        let NsfGasMedium::NamedMixture(source) = zone.medium() else {
            panic!("named source must not become a fitted GasSpec");
        };
        assert_eq!(source.species(), components.iter().map(|x| x.0).collect::<Vec<_>>());
        assert_eq!(source.composition().basis(), CompositionBasis::MoleFraction);
        assert_eq!(source.composition().fractions(), components.iter().map(|x| x.1).collect::<Vec<_>>());
        assert_eq!(source.reference_state(), gas);
        assert_eq!(system.reference_gas(), gas);
        assert_eq!(system.load().duct().inlet_gas(), gas);
        assert_eq!(system.upstream_load().unwrap().duct().inlet_gas(), gas);
        let reservoir = &zone.reservoir_aperture_at(0).unwrap().state;
        assert_eq!(reservoir.density_kg_m3, gas.density);
        assert_eq!(reservoir.temperature_k, gas.temperature);
        let warmed = resolve_gas_mixture(source.species(), source.composition(), 300.01, 100003.0).unwrap();
        let control = zone.reservoir_aperture_state_for(0, 100003.0, 300.01).unwrap();
        assert_eq!(control.density_kg_m3, warmed.density);
        assert_eq!(control.temperature_k, warmed.temperature);
        assert_eq!(file.renderer().drive_history()[0], NsfDgDrive::Reservoir {
            pressure_pa: 100002.0, temperature_k: 300.0,
        });
        assert_eq!(file.renderer().drive_history()[16], NsfDgDrive::Reservoir {
            pressure_pa: 100003.0, temperature_k: 300.01,
        });

        let (whole, owner) = render(file, 8);
        let (parts, replay) = render(read(&text, 3), 3);
        assert_eq!(whole, parts);
        assert_eq!(owner.system().zone().state_hash(), replay.system().zone().state_hash());
        assert_eq!(owner.system().load().state(), replay.system().load().state());
        assert_eq!(owner.system().upstream_load().unwrap().state(), replay.system().upstream_load().unwrap().state());
        assert_eq!(owner.system().zone().steps(), SAMPLES);
        let zone = owner.system().zone();
        let NsfGasMedium::NamedMixture(source) = zone.medium() else { unreachable!() };
        for cell in [[0, 1, 1], [5, 0, 1], [11, 1, 1]] {
            let local = zone.gas_state(cell).unwrap();
            assert_eq!(local, resolve_gas_mixture(source.species(), source.composition(),
                local.temperature, local.pressure).unwrap());
        }
    }
}

#[test]
fn g3_named_composition_and_source_temperature_change_pressure_at_fixed_geometry() {
    let text = input(MIXTURE);
    // Only composition changes: preserve the physical grid, ducts and clock.
    let helium = text.replace("N2 0.74", "N2 0.64").replace("He 0.015", "He 0.115");
    let baseline = read(&text, 8);
    let enriched = read(&helium, 8);
    assert_eq!(baseline.renderer().system().zone().law().dx_m, enriched.renderer().system().zone().law().dx_m);
    assert_ne!(baseline.renderer().system().zone().law().dt_s, enriched.renderer().system().zone().law().dt_s);
    let gas = baseline.renderer().system().reference_gas();
    let changed_gas = enriched.renderer().system().reference_gas();
    assert_ne!(gas.sound_speed, changed_gas.sound_speed);
    assert_ne!(gas.density, changed_gas.density);
    assert_ne!(gas.dynamic_viscosity, changed_gas.dynamic_viscosity);
    assert_ne!(gas.thermal_conductivity, changed_gas.thermal_conductivity);
    let (pressure, owner) = render(baseline, 8);
    let (changed, _) = render(enriched, 8);
    assert_ne!(pressure, changed);
    let unheated = text.replace("hold 32 100003 300.01", "hold 32 100003 300");
    let (cooler, cool_owner) = render(read(&unheated, 8), 8);
    assert_ne!(pressure, cooler);
    assert_ne!(owner.system().zone().state_hash(), cool_owner.system().zone().state_hash());
    assert_eq!(owner.system().reference_gas(), cool_owner.system().reference_gas());
}

#[test]
fn g0_native_mixture_controls_retain_source_domain_and_condensation_refusals() {
    let text = input(MIXTURE);
    for invalid in [
        text.replace("fluid_limits 108 290 320", "fluid_limits 108 250 600"),
        text.replace("mixture 300 100000", "mixture 324 100000"),
        text.replace("mixture 300 100000", "mixture 300 200001"),
        text.replace("hold 32 100003 300.01", "hold 32 9999 300.01"),
        text.replace("hold 32 100003 300.01", "hold 32 200001 300.01"),
    ] {
        assert!(ResolvedDgFilePerformance::from_bytes(invalid.as_bytes(), 8,
            Arc::new(CancelGate::new_clock_free())).is_err());
    }
    // This wetter reference is admitted across its full declared T window.
    // Its later pressure control crosses the original liquid-water limit.
    let wet = text.replace("N2 0.74", "N2 0.72")
        .replace("H2O 0.005", "H2O 0.025")
        .replace("fluid_limits 108 290 320", "fluid_limits 108 295 320");
    let admitted = read(&wet, 8);
    let zone = admitted.renderer().system().zone();
    assert!(zone.reservoir_aperture_state_for(0, 100003.0, 300.01).is_ok());
    assert!(zone.reservoir_aperture_state_for(0, 200000.0, 300.01).is_err());
    let condensing = wet.replace("hold 32 100003 300.01", "hold 32 200000 300.01");
    assert!(ResolvedDgFilePerformance::from_bytes(condensing.as_bytes(), 8,
        Arc::new(CancelGate::new_clock_free())).is_err());
}
