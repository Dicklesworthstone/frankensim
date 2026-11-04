//! Original source laws and independent energy balances; synthetic gas data.
use fs_bem::{
    helmholtz::laplace::complex::ComplexOptions,
    observation_transfers::{
        HarmonicField, HarmonicVolumeQuadrature, MolecularAbsorptionLimits,
        MolecularBemObservation, MolecularResponseLimits, PreparedMolecularResponse,
        ReceiverSpec, SourceId,
    },
    panel3d::SpherePanels,
};
use fs_couple::molecular_heating::{
    MolecularGasExcitation, MolecularGasHeating, MolecularHeatingError,
};
use fs_matdb::QueryPoint;
use fs_material::{
    relaxation::{RelaxationAbsorption, RelaxingMedium, input::read_relaxation},
    thermochemical_gas::{
        energy::{GasEnergySolve, IsochoricGasParcel},
        packs::from_packs,
    },
};
use fs_math::c64::C64;
use fs_qty::{
    QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{Composition, CompositionBasis};

#[path = "../src/render/schedule/reed/field_input/tests/support.rs"]
mod gas_support;

const VOLUME: f64 = 0.002;
const OMEGA: f64 = 1000.0;

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature, ValueForm::Static,
            )),
            temperature,
        ).unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::Pressure, ValueForm::Static,
            )),
            101325.0,
        ).unwrap()
}

fn gas(temperature: f64, volume: f64) -> IsochoricGasParcel {
    // The reused source fixture has cp/R = 3.5 + 0.001*T, not constant Cv.
    let packs = gas_support::packs("A", 0.028, 3.5, 200.0);
    let model = from_packs(
        &[packs.binding()],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        None,
    ).unwrap();
    IsochoricGasParcel::new(model, point(temperature), volume).unwrap()
}

fn law(relaxing: bool, time_scale: f64, upper_temperature: f64) -> RelaxationAbsorption {
    let data = if relaxing {
        format!(
            "declared 1 A 1 vibration 1 100 1 power {} 300 0 200 {}",
            101.325 * time_scale, upper_temperature,
        )
    } else {
        "declared 1 A 0".into()
    };
    read_relaxation(&mut data.split_ascii_whitespace()).unwrap()
}

fn response(gas: &IsochoricGasParcel, law: &RelaxationAbsorption) -> PreparedMolecularResponse {
    let surface = SpherePanels::icosphere(0.1, 0).unwrap();
    let panels = surface.areas().len();
    let source = law.evaluate(&RelaxingMedium::from_thermochemical(gas.state()).unwrap()).unwrap();
    let observation = MolecularBemObservation::new(
        &surface,
        source,
        &[
            SourceId::StructuralMode("first radial source".into()),
            SourceId::StructuralMode("second radial source".into()),
        ],
        &[vec![1.0; panels], vec![1.0; panels]],
        &[
            ReceiverSpec { id: "first gas sample".into(), position_m: [0.0, 0.0, 0.6] },
            ReceiverSpec { id: "second gas sample".into(), position_m: [0.0, 0.0, 1.0] },
        ],
        0.01,
        [100.0, 2000.0],
        ComplexOptions {
            max_panels: panels,
            max_bytes: 32 * panels * panels,
            max_work: 512 * (panels * panels) as u64,
        },
    ).unwrap();
    observation.prepare_response(
        OMEGA,
        MolecularResponseLimits {
            maximum_retained_entries: 4096,
            maximum_panel_products: 10000,
        },
        || false,
    ).unwrap()
}

fn quadrature() -> [HarmonicVolumeQuadrature; 2] {
    [
        HarmonicVolumeQuadrature { receiver: 0, volume_m3: 0.3 * VOLUME },
        HarmonicVolumeQuadrature { receiver: 1, volume_m3: 0.7 * VOLUME },
    ]
}

fn limits() -> MolecularAbsorptionLimits {
    MolecularAbsorptionLimits {
        maximum_quadrature_points: 2,
        maximum_coherence_modes: 2,
        maximum_products: 1000,
        maximum_workspace_bytes: 1 << 20,
    }
}

fn controls() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 450.0,
        energy_tolerance_j_kg: 1e-8,
        maximum_iterations: 64,
    }
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-30),
        "{actual} != {expected}",
    );
}

fn energy_for_increment(gas: &IsochoricGasParcel, increment: f64) -> f64 {
    let state = gas.state().state();
    gas.mass_kg() * state.specific_gas_constant
        * ((2.5 + 0.001 * state.temperature) * increment + 0.0005 * increment * increment)
}

#[test]
fn g1_molecular_heat_uses_original_nonlinear_internal_energy_and_eos() {
    let gas = gas(300.0, VOLUME);
    let original = gas.clone();
    let law = law(true, 1.0, 500.0);
    let response = response(&gas, &law);
    let quadrature = quadrature();
    let rates = [C64::new(0.4, -0.2), C64::new(-0.1, 0.3)];
    let mut fields = [HarmonicField::default(); 2];
    response.harmonic().evaluate_coherent(&rates, &mut fields, || false).unwrap();
    let compliance = response.source().caloric_response(OMEGA).unwrap().bulk_modulus_pa.recip();
    let expected_power: f64 = quadrature.iter().map(|point| {
        0.5 * OMEGA * compliance.im
            * fields[point.receiver].pressure_pa.norm_sq() * point.volume_m3
    }).sum();
    assert!(expected_power > 0.0);
    let expected_heat = energy_for_increment(&gas, 4.0);
    let duration = expected_heat / expected_power;
    let absorption = response.prepare_absorption(&quadrature, limits(), || false).unwrap();
    let mut heating = MolecularGasHeating::new(absorption, &gas, &law, || false).unwrap();
    let load = heating.heat_load(MolecularGasExcitation::Coherent(&rates), duration, || false).unwrap();
    close(load.absorbed_power_w(), expected_power, 2e-12);
    close(load.absorbed_energy_j(), expected_heat, 2e-12);
    assert_eq!(load.response().harmonic().producer_hash(), response.harmonic().producer_hash());
    assert_eq!(load.response().sources(), response.sources());
    assert_eq!(load.quadrature(), &quadrature);
    close(load.volume_m3(), VOLUME, 2e-15);
    let step = load.heat_gas(controls(), 5.0, || false).unwrap();
    let actual = step.gas().state().state();
    assert!((actual.temperature - 304.0).abs() < 1e-8);
    close(actual.pressure, 101325.0 * 304.0 / 300.0, 1e-10);
    close(actual.density, gas.state().state().density, 1e-13);
    assert_eq!(step.gas().mass_kg(), gas.mass_kg());
    assert_eq!(step.gas().volume_m3(), gas.volume_m3());
    assert_eq!(step.gas().model(), gas.model());
    close(step.energy_balance().stored_energy_change_j, expected_heat, 1e-10);
    assert!(
        step.energy_balance().residual_j.abs()
            <= gas.mass_kg() * controls().energy_tolerance_j_kg
    );
    let wrong_cp_increment = expected_heat / (gas.mass_kg() * gas.state().state().specific_heat_cp);
    assert!((wrong_cp_increment - 4.0).abs() > 0.5, "fixture must distinguish Cp from Cv");
    assert_eq!(
        step.next_relaxation(),
        &law.evaluate(&RelaxingMedium::from_thermochemical(step.gas().state()).unwrap()).unwrap(),
    );
    assert_ne!(step.next_relaxation(), response.source());
    assert_ne!(
        step.next_relaxation().processes()[0].relaxation_time_s,
        response.source().processes()[0].relaxation_time_s,
    );
    // A heated gas cannot be rebound to the old harmonic field as if feedback
    // had already been solved. The original accepted gas remains unchanged.
    let stale = response.prepare_absorption(&quadrature, limits(), || false).unwrap();
    assert!(MolecularGasHeating::new(stale, step.gas(), &law, || false).is_err());
    assert_eq!(gas, original);
}

#[test]
fn g1_coherence_and_lossless_limits_reach_native_gas_storage() {
    let gas = gas(300.0, VOLUME);
    let quadrature = quadrature();
    let law = law(true, 1.0, 500.0);
    let response = response(&gas, &law);
    let absorption = response.prepare_absorption(&quadrature, limits(), || false).unwrap();
    let mut heating = MolecularGasHeating::new(absorption, &gas, &law, || false).unwrap();
    let single = heating.heat_load(
        MolecularGasExcitation::Coherent(&[C64::ONE, C64::ZERO]), 2.0, || false,
    ).unwrap().absorbed_energy_j();
    let coherent = heating.heat_load(
        MolecularGasExcitation::Coherent(&[C64::ONE, C64::ONE]), 2.0, || false,
    ).unwrap().absorbed_energy_j();
    let independent = heating.heat_load(
        MolecularGasExcitation::IndependentSources(&[1.0, 1.0]), 2.0, || false,
    ).unwrap().absorbed_energy_j();
    let modes = heating.heat_load(
        MolecularGasExcitation::CoherenceModes {
            factors: &[C64::ONE, C64::ZERO, C64::ZERO, C64::ONE], modes: 2,
        }, 2.0, || false,
    ).unwrap().absorbed_energy_j();
    close(coherent, 4.0 * single, 2e-12);
    close(independent, 2.0 * single, 2e-12);
    close(modes, independent, 2e-12);
    let cancelled_field = heating.heat_load(
        MolecularGasExcitation::Coherent(&[C64::ONE, C64::from_re(-1.0)]), 2.0, || false,
    ).unwrap();
    assert_eq!(cancelled_field.absorbed_energy_j(), 0.0);
    let unchanged = cancelled_field.heat_gas(controls(), 1.0, || false).unwrap();
    assert_eq!(unchanged.gas(), &gas);
    assert_eq!(unchanged.energy_balance().stored_energy_change_j, 0.0);

    let empty_law = self::law(false, 1.0, 500.0);
    let lossless = self::response(&gas, &empty_law);
    let absorption = lossless.prepare_absorption(&quadrature, limits(), || false).unwrap();
    let mut heating = MolecularGasHeating::new(absorption, &gas, &empty_law, || false).unwrap();
    let load = heating.heat_load(
        MolecularGasExcitation::Coherent(&[C64::ONE, C64::ONE]), 1e6, || false,
    ).unwrap();
    assert_eq!(load.absorbed_power_w(), 0.0);
    assert_eq!(load.heat_gas(controls(), 1.0, || false).unwrap().gas(), &gas);
}

#[test]
fn g3_binding_requires_complete_medium_reservoirs_and_physical_volume() {
    let gas = gas(300.0, VOLUME);
    let original = gas.clone();
    let law = law(true, 1.0, 500.0);
    let response = response(&gas, &law);
    let quadrature = quadrature();
    let wrong_volume = self::gas(300.0, 2.0 * VOLUME);
    let wrong_temperature = self::gas(301.0, VOLUME);
    let wrong_relaxation = self::law(true, 2.0, 500.0);
    for (parcel, relaxation) in [
        (&wrong_volume, &law),
        (&wrong_temperature, &law),
        (&gas, &wrong_relaxation),
    ] {
        let absorption = response.prepare_absorption(&quadrature, limits(), || false).unwrap();
        assert!(MolecularGasHeating::new(absorption, parcel, relaxation, || false).is_err());
    }
    let absorption = response.prepare_absorption(&quadrature, limits(), || false).unwrap();
    assert!(matches!(
        MolecularGasHeating::new(absorption, &gas, &law, || true),
        Err(MolecularHeatingError::Cancelled),
    ));
    assert_eq!(gas, original);
}

#[test]
fn g3_temperature_domain_and_late_cancellation_refuse_atomically() {
    let gas = gas(300.0, VOLUME);
    let original = gas.clone();
    let law = law(true, 1.0, 500.0);
    let response = response(&gas, &law);
    let quadrature = quadrature();
    let absorption = response.prepare_absorption(&quadrature, limits(), || false).unwrap();
    let mut heating = MolecularGasHeating::new(absorption, &gas, &law, || false).unwrap();
    let rates = [C64::ONE, C64::ZERO];
    for duration in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(heating.heat_load(
            MolecularGasExcitation::Coherent(&rates), duration, || false,
        ).is_err());
    }
    let watts = heating.heat_load(
        MolecularGasExcitation::Coherent(&rates), 1.0, || false,
    ).unwrap().absorbed_power_w();
    let duration = energy_for_increment(&gas, 4.0) / watts;
    let load = heating.heat_load(
        MolecularGasExcitation::Coherent(&rates), duration, || false,
    ).unwrap();
    assert!(matches!(
        load.heat_gas(controls(), 3.0, || false),
        Err(MolecularHeatingError::TemperatureIncrement { .. }),
    ));
    let mut polls = 0;
    load.heat_gas(controls(), 5.0, || { polls += 1; false }).unwrap();
    assert!(polls > 3);
    for stop in [1, polls] {
        let mut calls = 0;
        assert!(load.heat_gas(controls(), 5.0, || {
            calls += 1;
            calls == stop
        }).is_err());
        assert_eq!(gas, original);
    }
    let mut insufficient = controls();
    insufficient.maximum_temperature_k = 303.0;
    assert!(load.heat_gas(insufficient, 5.0, || false).is_err());
    assert_eq!(gas, original);

    // The gas source can support 304 K while the original collision law cannot.
    // Refusal happens after caloric inversion, without publishing a hot parcel.
    let narrow_law = self::law(true, 1.0, 303.0);
    let narrow_response = self::response(&gas, &narrow_law);
    let absorption = narrow_response.prepare_absorption(&quadrature, limits(), || false).unwrap();
    let mut narrow = MolecularGasHeating::new(absorption, &gas, &narrow_law, || false).unwrap();
    let load = narrow.heat_load(
        MolecularGasExcitation::Coherent(&rates), duration, || false,
    ).unwrap();
    assert!(matches!(
        load.heat_gas(controls(), 5.0, || false),
        Err(MolecularHeatingError::Material(_)),
    ));
    assert_eq!(gas, original);
}
