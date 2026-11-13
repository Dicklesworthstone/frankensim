//! Original source laws and independent energy balances; synthetic gas data.
#[path = "molecular_heating/finite_window.rs"]
mod finite_window;

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
    spectrum::{MAX_MOLECULAR_SPECTRUM_BINS, MolecularGasSpectrumHeating},
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

#[path = "../../src/render/schedule/reed/field_input/tests/support.rs"]
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
    response_at(gas, law, OMEGA)
}

fn response_at(
    gas: &IsochoricGasParcel,
    law: &RelaxationAbsorption,
    omega: f64,
) -> PreparedMolecularResponse {
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
        omega,
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

fn spectrum<'a>(
    responses: &'a [PreparedMolecularResponse],
    quadrature: &'a [HarmonicVolumeQuadrature],
    gas: &'a IsochoricGasParcel,
    law: &'a RelaxationAbsorption,
) -> MolecularGasSpectrumHeating<'a> {
    let absorptions = responses.iter().map(|response| {
        response.prepare_absorption(quadrature, limits(), || false).unwrap()
    }).collect();
    MolecularGasSpectrumHeating::new(absorptions, gas, law, responses.len(), || false).unwrap()
}

fn direct_power(
    response: &PreparedMolecularResponse,
    quadrature: &[HarmonicVolumeQuadrature],
    rates: &[C64],
) -> f64 {
    let mut fields = [HarmonicField::default(); 2];
    response.harmonic().evaluate_coherent(rates, &mut fields, || false).unwrap();
    let omega = response.harmonic().angular_frequency_rad_s();
    let compliance = response.source().caloric_response(omega).unwrap().bulk_modulus_pa.recip();
    quadrature.iter().map(|point| {
        0.5 * omega * compliance.im
            * fields[point.receiver].pressure_pa.norm_sq() * point.volume_m3
    }).sum()
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

#[test]
fn g1_spectrum_sums_original_frequency_powers_before_nonlinear_gas_heating() {
    let gas = gas(300.0, VOLUME);
    let original = gas.clone();
    let law = law(true, 1.0, 500.0);
    // Preserve input order even when it differs from the summation order.
    let responses = [response_at(&gas, &law, 1500.0), response_at(&gas, &law, 500.0)];
    let quadrature = quadrature();
    let rates = [
        [C64::new(0.4, -0.2), C64::new(-0.1, 0.3)],
        [C64::new(-0.6, 0.1), C64::new(0.2, 0.4)],
    ];
    let powers = [
        direct_power(&responses[0], &quadrature, &rates[0]),
        direct_power(&responses[1], &quadrature, &rates[1]),
    ];
    assert!(powers.iter().all(|power| *power > 0.0));
    let expected_heat = energy_for_increment(&gas, 20.0);
    let duration = expected_heat / (powers[0] + powers[1]);
    let mut heating = spectrum(&responses, &quadrature, &gas, &law);
    let load = heating.heat_load(&[
        MolecularGasExcitation::Coherent(&rates[0]),
        MolecularGasExcitation::Coherent(&rates[1]),
    ], duration, || false).unwrap();
    close(load.absorbed_power_w(), powers[0] + powers[1], 2e-12);
    close(load.absorbed_energy_j(), expected_heat, 2e-12);
    assert_eq!(load.volume_m3(), VOLUME);
    assert_eq!(load.duration_s(), duration);
    for (index, bin) in load.bins().iter().enumerate() {
        assert!(core::ptr::eq(bin.response(), &responses[index]));
        assert_eq!(bin.quadrature(), &quadrature);
        close(bin.absorbed_power_w(), powers[index], 2e-12);
    }
    let step = load.heat_gas(controls(), 21.0, || false).unwrap();
    assert!((step.gas().state().state().temperature - 320.0).abs() < 1e-8);
    close(step.gas().state().state().pressure, 101325.0 * 320.0 / 300.0, 1e-10);
    close(step.energy_balance().stored_energy_change_j, expected_heat, 1e-10);
    assert!(step.energy_balance().residual_j.abs()
        <= gas.mass_kg() * controls().energy_tolerance_j_kg);
    assert_eq!(step.gas().mass_kg(), gas.mass_kg());
    assert_eq!(step.gas().volume_m3(), VOLUME);
    assert_eq!(step.gas().model(), gas.model());
    assert_eq!(step.next_relaxation(), &law.evaluate(
        &RelaxingMedium::from_thermochemical(step.gas().state()).unwrap(),
    ).unwrap());
    assert_ne!(step.next_relaxation(), responses[0].source());
    let wrong_linear_cv = expected_heat / (gas.mass_kg()
        * (gas.state().state().specific_heat_cp - gas.state().state().specific_gas_constant));
    assert!((wrong_linear_cv - 20.0).abs() > 0.05);
    assert_eq!(gas, original);
}

#[test]
fn g1_spectrum_preserves_coherence_within_each_frequency_only() {
    let gas = gas(300.0, VOLUME);
    let law = law(true, 1.0, 500.0);
    let responses = [response_at(&gas, &law, 500.0), response_at(&gas, &law, 1500.0)];
    let quadrature = quadrature();
    let mut heating = spectrum(&responses, &quadrature, &gas, &law);
    let plus = [C64::ONE, C64::ZERO];
    let minus = [C64::from_re(-1.0), C64::ZERO];
    let expected = direct_power(&responses[0], &quadrature, &plus)
        + direct_power(&responses[1], &quadrature, &minus);
    let opposite = heating.heat_load(&[
        MolecularGasExcitation::Coherent(&plus),
        MolecularGasExcitation::Coherent(&minus),
    ], 2.0, || false).unwrap().absorbed_power_w();
    assert!(opposite > 0.0, "opposite phases at different frequencies cannot cancel mean heat");
    close(opposite, expected, 2e-12);
    let two = [C64::ONE, C64::ONE];
    let coherent = heating.heat_load(&[MolecularGasExcitation::Coherent(&two); 2],
        2.0, || false).unwrap().absorbed_power_w();
    close(coherent, 4.0 * expected, 2e-12);
    let independent = heating.heat_load(&[
        MolecularGasExcitation::IndependentSources(&[1.0, 1.0]); 2
    ], 2.0, || false).unwrap().absorbed_power_w();
    let factors = [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE];
    let modes = heating.heat_load(&[
        MolecularGasExcitation::CoherenceModes { factors: &factors, modes: 2 }; 2
    ], 2.0, || false).unwrap().absorbed_power_w();
    close(independent, 2.0 * expected, 2e-12);
    close(modes, independent, 2e-12);
    let cancelling = [C64::ONE, C64::from_re(-1.0)];
    let load = heating.heat_load(&[MolecularGasExcitation::Coherent(&cancelling); 2],
        2.0, || false).unwrap();
    assert_eq!(load.absorbed_energy_j(), 0.0);
    assert_eq!(load.heat_gas(controls(), 1.0, || false).unwrap().gas(), &gas);
}

#[test]
fn g3_spectrum_rejects_duplicate_frequencies_and_retries_cancelled_queries() {
    let gas = gas(300.0, VOLUME);
    let original = gas.clone();
    let law = law(true, 1.0, 500.0);
    let responses = [response_at(&gas, &law, 500.0), response_at(&gas, &law, 1500.0)];
    let quadrature = quadrature();
    let prepare = |index: usize| {
        responses[index].prepare_absorption(&quadrature, limits(), || false).unwrap()
    };
    assert!(MolecularGasSpectrumHeating::new(
        vec![prepare(0), prepare(0)], &gas, &law, 2, || false,
    ).is_err());
    for cap in [0, 1, MAX_MOLECULAR_SPECTRUM_BINS + 1] {
        assert!(MolecularGasSpectrumHeating::new(
            vec![prepare(0), prepare(1)], &gas, &law, cap, || false,
        ).is_err());
    }
    assert!(MolecularGasSpectrumHeating::new(Vec::new(), &gas, &law, 2, || false).is_err());
    let wrong_gas = self::gas(301.0, VOLUME);
    let wrong_response = response_at(&wrong_gas, &law, 1500.0);
    assert!(MolecularGasSpectrumHeating::new(vec![
        prepare(0), wrong_response.prepare_absorption(&quadrature, limits(), || false).unwrap(),
    ], &gas, &law, 2, || false).is_err());
    let mut short_volume = quadrature;
    short_volume[1].volume_m3 *= 0.5;
    assert!(MolecularGasSpectrumHeating::new(vec![
        prepare(0),
        responses[1].prepare_absorption(&short_volume, limits(), || false).unwrap(),
    ], &gas, &law, 2, || false).is_err());

    let mut heating = spectrum(&responses, &quadrature, &gas, &law);
    let rates = [C64::ONE, C64::ZERO];
    let excitations = [MolecularGasExcitation::Coherent(&rates); 2];
    assert!(heating.heat_load(&excitations[..1], 1.0, || false).is_err());
    assert!(heating.heat_load(&excitations, 0.0, || false).is_err());
    let power = heating.heat_load(&excitations, 1.0, || false).unwrap().absorbed_power_w();
    let energy = energy_for_increment(&gas, 4.0);
    let duration = energy / power;
    let mut polls = 0;
    heating.heat_load(&excitations, duration, || { polls += 1; false }).unwrap();
    for stop in [1, polls / 2, polls] {
        let mut calls = 0;
        assert!(heating.heat_load(&excitations, duration, || {
            calls += 1;
            calls == stop
        }).is_err());
        assert_eq!(gas, original);
        close(heating.heat_load(&excitations, duration, || false).unwrap().absorbed_energy_j(),
            energy, 2e-12);
    }
    let load = heating.heat_load(&excitations, duration, || false).unwrap();
    assert!(matches!(load.heat_gas(controls(), 3.0, || false),
        Err(MolecularHeatingError::TemperatureIncrement { .. })));
    let mut polls = 0;
    load.heat_gas(controls(), 5.0, || { polls += 1; false }).unwrap();
    let mut calls = 0;
    assert!(load.heat_gas(controls(), 5.0, || {
        calls += 1;
        calls == polls
    }).is_err());
    assert!((load.heat_gas(controls(), 5.0, || false).unwrap()
        .gas().state().state().temperature - 304.0).abs() < 1e-8);
    let narrow_law = self::law(true, 1.0, 303.0);
    let narrow_responses = [
        response_at(&gas, &narrow_law, 500.0), response_at(&gas, &narrow_law, 1500.0),
    ];
    let mut narrow = spectrum(&narrow_responses, &quadrature, &gas, &narrow_law);
    let load = narrow.heat_load(&excitations, duration, || false).unwrap();
    assert!(matches!(load.heat_gas(controls(), 5.0, || false),
        Err(MolecularHeatingError::Material(_))));
    assert_eq!(gas, original);
}

#[test]
fn g1_spectrum_and_finite_body_exchange_credit_absorption_once() {
    use fs_blake3::ContentHash;
    use fs_conduction::lumped::{LumpedThermalEnvironment, LumpedThermalTransport};
    use fs_couple::cavity_phs::piston::thermal::HeatExchangeSolve;
    use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
    use fs_material::phase::{
        EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase, UniformEnthalpyStepInput,
    };

    let gas = gas(300.0, VOLUME);
    let original = gas.clone();
    let law = law(true, 1.0, 500.0);
    let responses = [response_at(&gas, &law, 500.0), response_at(&gas, &law, 1500.0)];
    let quadrature = quadrature();
    let mut heating = spectrum(&responses, &quadrature, &gas, &law);
    let rates = [C64::ONE, C64::ZERO];
    let excitations = [MolecularGasExcitation::Coherent(&rates); 2];
    let power = heating.heat_load(&excitations, 1.0, || false).unwrap().absorbed_power_w();
    let energy = energy_for_increment(&gas, 4.0);
    let duration = energy / power;
    let load = heating.heat_load(&excitations, duration, || false).unwrap();
    let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([69; 32]), SolidLiquidPhase::Solid,
        [280.0, 340.0].map(|t| EnthalpyPhaseKnot {
            temperature_k: t, specific_enthalpy_j_kg: 1000.0 * (t - 300.0),
            liquid_mass_fraction: 0.0, bulk_density_kg_m3: 1000.0,
        }).to_vec(),
    ).unwrap();
    let body = LumpedThermalEnvironment {
        temperature_k: 999.0, // Replaced by the actual candidate gas endpoint.
        radiation_temperature_k: 300.0,
        convection_w_per_m2_k: 0.5 / (0.01 * duration),
        transport: LumpedThermalTransport::try_declared(100.0, 0.0).unwrap(),
        enthalpy_tolerance_j_kg: 1e-8, maximum_thermal_residual_j: 1e-9,
    };
    let input = UniformEnthalpyStepInput {
        curve: &curve, initial: curve.state_at_specific_enthalpy(-2000.0).unwrap(),
        mass_kg: 0.002, volume_m3: 2e-6, surface_area_m2: 0.01,
        internal_heat_j: 0.0, duration_s: duration,
    };
    let environment = load.finite_body_environment(&body, controls(), HeatExchangeSolve {
        heat_bracket_j: [-20.0, 20.0], source_tolerance_j: 1e-9,
        energy_tolerance_j: 1e-8, maximum_iterations: 64,
    }, 20.0).unwrap();
    assert!(core::ptr::eq(environment.gas, &gas));
    close(environment.gas_heat_j, energy, 2e-12);
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena,
            StreamKey { seed: 0, kernel_id: 0x4d4753, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let result = environment.advance(&cx, input).unwrap();
        // Eliminate the 2 J/K body's endpoint from Q = dt*h*A*(Tb - Tg).
        // The remaining quadratic uses the independently integrated gas Cv(T).
        let a = 0.5 / (1.0 + 0.5 / 2.0);
        let mr = gas.mass_kg() * gas.state().state().specific_gas_constant;
        let b = 2.8 * mr + a;
        let c = energy - 2.0 * a;
        let increment = 2.0 * c / (b + (b * b + 4.0 * 0.0005 * mr * c).sqrt());
        assert!((result.report.next_gas.state().state().temperature - 300.0 - increment).abs() < 1e-7);
        let exchange = -a * (2.0 + increment);
        assert!((result.report.solid_to_gas_heat_j - exchange).abs() < 1e-8);
        assert!((result.state.temperature_k() - 298.0 + exchange / 2.0).abs() < 1e-7);
        close(result.report.external_gas_heat_j, energy, 2e-12);
        let gas_storage = result.report.gas.stored_energy_change_j;
        let body_storage = input.mass_kg * (result.state.specific_enthalpy_j_kg()
            - input.initial.specific_enthalpy_j_kg());
        assert!((gas_storage + body_storage - energy).abs() < 1e-8);
        assert!(result.report.combined_energy_residual_j.abs() < 1e-8);
    });
    assert_eq!(gas, original);
}
