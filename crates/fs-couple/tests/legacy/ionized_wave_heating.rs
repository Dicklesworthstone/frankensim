#![cfg(feature = "ionized-wave-heating")]
//! Actual Maxwell absorption, native gas calorics and original-source refresh.
#[path = "../../../fs-material/tests/common/ionized_gas.rs"]
mod fixture;

use fs_couple::ionized_wave_heating::{IonizedWaveHeatingError, MagnetizedIonizedGasHeating};
use fs_em::plane_wave::{DielectricResponse, plasma::magnetized::MagnetizedIonizedGasInterface};
use fs_evidence::ValidityDomain;
use fs_gas::chemistry::energy::ElectronHeavyEnergySolve;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS,
        MAGNETIC_PERMEABILITY_PROPERTY,
    },
    ionized_gas::{IonizedGasConductivity, IonizedGasMagneticLimits},
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_math::c64::C64;
use fs_qty::{
    QuantitySpec, Temperature, Time,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::electron_collision::MaxwellianRateControlsV1;

#[path = "ionized_wave_heating/reacting.rs"]
mod reacting;

#[path = "ionized_wave_heating/layers.rs"]
mod layers;

fn magnetic_limits() -> IonizedGasMagneticLimits {
    IonizedGasMagneticLimits {
        maximum_cyclotron_to_thermal_energy_ratio: 0.01,
    }
}
fn controls() -> ElectronHeavyEnergySolve {
    let mut controls = ElectronHeavyEnergySolve::new(
        Temperature::new(2000.0),
        Temperature::new(20000.0),
        0.2,
        0.1,
        MaxwellianRateControlsV1 {
            relative_tolerance: 1e-10,
            maximum_subdivisions: 1024,
        },
    );
    controls.energy_tolerance_j = 1e-8;
    controls.temperature_tolerance_k = 1e-7;
    controls.relative_tolerance = 1e-6;
    controls.initial_step_s = 1e-9;
    controls.maximum_step_s = 1e-7;
    controls
}
fn background(source: &IonizedGasConductivity) -> (ResolvedMaterialStatePoint, MaterialCard) {
    let t = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let p = QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
    let validity = ValidityDomain::unconstrained()
        .with_quantity("T_h", t, 2000.0, 20000.0)
        .with_quantity("T_e", t, 2000.0, 20000.0)
        .with_quantity("p", p, 1.0, 1e8);
    let query = QueryPoint::new()
        .with_quantity("T_h", t, source.state().heavy_temperature().value())
        .unwrap()
        .with_quantity("T_e", t, source.state().electron_temperature().value())
        .unwrap()
        .with_quantity("p", p, source.state().pressure().value())
        .unwrap();
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in [
        (
            ELECTRIC_PERMITTIVITY_PROPERTY,
            ELECTRIC_PERMITTIVITY_DIMS,
            1e-11,
        ),
        (
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            1e-6,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: validity.clone(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured bound-charge background for plasma heating".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
                .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "original plasma background".into(),
            phase: "gas".into(),
            process: "finite energy test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    (
        resolve_material_state_point(
            &card,
            &query,
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
        card,
    )
}

fn with_heating(
    coulomb: bool,
    th: f64,
    te: f64,
    excited: bool,
    test: impl FnOnce(&MagnetizedIonizedGasHeating<'_, '_>, &MaterialCard, &MaterialCard),
) {
    let state = fixture::frozen_state(0.001, th, te, 1.0);
    let (ion, ion_card) = fixture::collision_source_and_card(
        "H+",
        &[("T_h", th), ("T_e", te)],
        state.pressure().value(),
        1e9,
    );
    let mut source = IonizedGasConductivity::from_state(
        state.clone(),
        vec![ion, fixture::electron_relaxation(&state, 2e-20)],
        fixture::limits(),
        || false,
    )
    .unwrap();
    if coulomb {
        source = source
            .with_coulomb_drag(fixture::coulomb_model(&state), || false)
            .unwrap();
    }
    let harmonic = source.prepare_harmonic(|| false).unwrap();
    let omega = 1e10;
    let transport = harmonic
        .prepare_magnetized(
            omega,
            [0.0, 0.0, 0.4 * omega * fixture::ME / fixture::E],
            magnetic_limits(),
            || false,
        )
        .unwrap();
    let (background, background_card) = background(&source);
    let interface = MagnetizedIonizedGasInterface::try_new(
        DielectricResponse::Lossless(&background),
        &transport,
        &background,
        [0.0, 0.0, 1.0],
        || false,
    )
    .unwrap();
    let unit = interface
        .prepare_wave([C64::ONE, C64::ZERO, C64::ZERO], || false)
        .unwrap();
    let mut amplitude = f64::INFINITY;
    for (carrier, velocity) in source
        .carriers()
        .iter()
        .zip(&unit.surface_response().carrier_velocity_m_s)
    {
        let squared: f64 = velocity.iter().map(|v| v.norm_sq()).sum();
        let bilinear = velocity.iter().fold(C64::ZERO, |sum, &v| sum + v * v);
        let peak = ((squared + bilinear.abs()) * 0.5).sqrt();
        amplitude = amplitude.min(
            0.1 * source.limits().maximum_drift_thermal_ratio
                * (fixture::KB * carrier.temperature_k() / carrier.mass_kg()).sqrt()
                / peak,
        );
    }
    if !excited {
        amplitude = 0.0;
    }
    let wave = interface
        .prepare_wave([C64::from_re(amplitude), C64::ZERO, C64::ZERO], || false)
        .unwrap();
    let alpha = interface
        .circular_modes()
        .iter()
        .map(|m| -m.transmitted_normal_wavenumber_rad_m.im)
        .fold(0.0_f64, f64::max);
    let width = 0.7 / alpha;
    let slab = wave
        .slab_absorption(state.volume().value() / width, 0.0, width)
        .unwrap();
    let heating = MagnetizedIonizedGasHeating::try_new(slab, || false).unwrap();
    test(&heating, &ion_card, &background_card);
}
fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-100),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn wave_work_advances_original_gas_and_next_sources_feed_a_new_maxwell_solve() {
    with_heating(
        false,
        8000.0,
        8000.0,
        true,
        |heating, ion_card, background_card| {
            let duration = Time::new(1e-6);
            let step = heating
                .advance(
                    duration,
                    controls(),
                    200.0,
                    &[ion_card],
                    background_card,
                    || false,
                )
                .unwrap();
            let initial = heating.source().state();
            let next = step.energy().state();
            let balance = step.energy().balance();
            assert_eq!(next.amounts(), initial.amounts());
            assert_eq!(next.volume(), initial.volume());
            assert_eq!(next.thermodynamics(), initial.thermodynamics());
            assert_eq!(step.energy().initial_state(), initial);
            assert_eq!(step.energy().momentum_sources(), heating.momentum_sources());
            assert!(next.electron_temperature().value() > initial.electron_temperature().value());
            assert!(next.heavy_temperature().value() > initial.heavy_temperature().value());
            near(
                balance.electron_deposited_energy_j,
                duration.value() * heating.absorption().electron_power_w(),
                2e-14,
            );
            near(
                balance.heavy_deposited_energy_j,
                duration.value() * heating.absorption().heavy_power_w(),
                2e-14,
            );
            near(
                balance.electron_deposited_energy_j + balance.heavy_deposited_energy_j,
                step.absorbed_energy_j(),
                2e-11,
            );
            // Independent monatomic Cv closure; formation references cancel once.
            let electron_change = initial.electron_heat_capacity_j_k()
                * (next.electron_temperature().value() - initial.electron_temperature().value());
            let heavy_change = initial.heavy_heat_capacity_j_k()
                * (next.heavy_temperature().value() - initial.heavy_temperature().value());
            assert!(
                (electron_change - balance.electron_deposited_energy_j
                    + balance.electron_to_heavy_energy_j)
                    .abs()
                    < 3e-8
            );
            assert!(
                (heavy_change
                    - balance.heavy_deposited_energy_j
                    - balance.electron_to_heavy_energy_j)
                    .abs()
                    < 3e-8
            );
            assert!((electron_change + heavy_change - step.absorbed_energy_j()).abs() < 4e-8);
            assert_eq!(step.refreshed_transport().state(), next);
            assert_ne!(
                step.refreshed_transport().conductivity_s_m(),
                heating.source().conductivity_s_m()
            );
            assert_eq!(
                step.refreshed_background().card_identity(),
                background_card.content_hash()
            );
            assert_eq!(
                step.refreshed_background().state_coordinate("T_e"),
                Some(next.electron_temperature().value())
            );
            assert_eq!(
                step.refreshed_background().state_coordinate("p"),
                Some(next.pressure().value())
            );
            let previous = heating.absorption().wave();
            let harmonic = step
                .refreshed_transport()
                .prepare_harmonic(|| false)
                .unwrap();
            let transport = harmonic
                .prepare_magnetized(
                    previous.interface().transport().omega_rad_s(),
                    previous.interface().transport().magnetic_field_t(),
                    magnetic_limits(),
                    || false,
                )
                .unwrap();
            let interface = MagnetizedIonizedGasInterface::try_new(
                previous.interface().incident_response(),
                &transport,
                step.refreshed_background(),
                previous.interface().outgoing_normal(),
                || false,
            )
            .unwrap();
            let next_wave = interface
                .prepare_wave(previous.incident_electric_v_m(), || false)
                .unwrap();
            assert_ne!(
                next_wave.transmitted_power_w_m2(),
                previous.transmitted_power_w_m2()
            );
            assert!(next_wave.power_balance_relative_residual() < 1e-10);
        },
    );
}

#[test]
fn signed_coulomb_drift_heat_is_supplied_once_and_baseline_exchange_remains_internal() {
    with_heating(
        true,
        6000.0,
        12000.0,
        true,
        |heating, ion_card, background_card| {
            let duration = Time::new(1e-10);
            let step = heating
                .advance(
                    duration,
                    controls(),
                    200.0,
                    &[ion_card],
                    background_card,
                    || false,
                )
                .unwrap();
            let balance = step.energy().balance();
            assert_ne!(
                heating.absorption().exchange_correction_to_electrons_w(),
                0.0
            );
            near(
                balance.electron_deposited_energy_j,
                duration.value()
                    * (heating.absorption().electron_friction_power_w()
                        + heating.absorption().exchange_correction_to_electrons_w()),
                3e-13,
            );
            near(
                balance.heavy_deposited_energy_j,
                duration.value()
                    * (heating.absorption().heavy_friction_power_w()
                        - heating.absorption().exchange_correction_to_electrons_w()),
                3e-13,
            );
            near(
                balance.electron_deposited_energy_j + balance.heavy_deposited_energy_j,
                step.absorbed_energy_j(),
                3e-11,
            );
            assert!(balance.electron_to_heavy_coulomb_energy_j > 0.0);
            assert!(balance.electron_to_heavy_neutral_energy_j > 0.0);
            near(
                balance.electron_to_heavy_energy_j,
                balance.electron_to_heavy_coulomb_energy_j
                    + balance.electron_to_heavy_neutral_energy_j,
                2e-13,
            );
            assert_eq!(
                step.energy().coulomb_model(),
                heating.source().coulomb_model()
            );
            assert_eq!(
                step.refreshed_transport().coulomb_model(),
                heating.source().coulomb_model()
            );
            assert!(step.energy().state().electron_temperature().value() < 12000.0);
            assert!(step.energy().state().heavy_temperature().value() > 6000.0);
        },
    );
}

#[test]
fn zero_wave_energy_allows_original_thermal_exchange_without_spurious_field_heat() {
    with_heating(
        false,
        6000.0,
        12000.0,
        false,
        |heating, ion_card, background_card| {
            let step = heating
                .advance(
                    Time::new(1e-10),
                    controls(),
                    200.0,
                    &[ion_card],
                    background_card,
                    || false,
                )
                .unwrap();
            let initial = heating.source().state();
            let next = step.energy().state();
            assert_eq!(step.absorbed_energy_j(), 0.0);
            assert_eq!(step.energy().balance().electron_deposited_energy_j, 0.0);
            assert_eq!(step.energy().balance().heavy_deposited_energy_j, 0.0);
            assert!(step.energy().exchanged_energy().value() > 0.0);
            let change = initial.electron_heat_capacity_j_k()
                * (next.electron_temperature().value() - initial.electron_temperature().value())
                + initial.heavy_heat_capacity_j_k()
                    * (next.heavy_temperature().value() - initial.heavy_temperature().value());
            assert!(change.abs() < 4e-8);
        },
    );
}

#[test]
fn volume_source_increment_domain_and_cancellation_refusals_preserve_the_old_wave() {
    with_heating(
        false,
        8000.0,
        8000.0,
        true,
        |heating, ion_card, background_card| {
            let saved = heating.source().clone();
            let before = heating.absorption().absorbed_power_w();
            let wave = heating.absorption().wave();
            let wrong_volume = wave
                .slab_absorption(
                    2.0 * heating.absorption().area_m2(),
                    0.0,
                    heating.absorption().end_depth_m(),
                )
                .unwrap();
            assert!(MagnetizedIonizedGasHeating::try_new(wrong_volume, || false).is_err());
            assert!(matches!(
                heating.advance(
                    Time::new(1e-7),
                    controls(),
                    1e-12,
                    &[ion_card],
                    background_card,
                    || false
                ),
                Err(IonizedWaveHeatingError::TemperatureIncrement { .. })
            ));
            assert!(
                heating
                    .advance(
                        Time::new(1e-10),
                        controls(),
                        200.0,
                        &[],
                        background_card,
                        || false
                    )
                    .is_err()
            );
            assert!(
                heating
                    .advance(
                        Time::new(1e-10),
                        controls(),
                        200.0,
                        &[ion_card],
                        ion_card,
                        || false
                    )
                    .is_err()
            );
            let mut wide = controls();
            wide.maximum_target_thermal_speed_ratio = 0.2;
            assert!(
                heating
                    .advance(
                        Time::new(1e-10),
                        wide,
                        200.0,
                        &[ion_card],
                        background_card,
                        || false
                    )
                    .is_err()
            );
            assert!(matches!(
                heating.advance(
                    Time::new(1e-10),
                    controls(),
                    200.0,
                    &[ion_card],
                    background_card,
                    || true
                ),
                Err(IonizedWaveHeatingError::Cancelled)
            ));
            let mut polls = 0;
            assert!(matches!(
                heating.advance(
                    Time::new(1e-10),
                    controls(),
                    200.0,
                    &[ion_card],
                    background_card,
                    || {
                        polls += 1;
                        polls > 8
                    }
                ),
                Err(IonizedWaveHeatingError::Cancelled)
            ));
            assert!(
                heating
                    .advance(
                        Time::new(0.0),
                        controls(),
                        200.0,
                        &[ion_card],
                        background_card,
                        || false
                    )
                    .is_err()
            );
            assert_eq!(heating.source(), &saved);
            assert_eq!(heating.absorption().absorbed_power_w(), before);
            assert!(
                heating
                    .advance(
                        Time::new(1e-10),
                        controls(),
                        200.0,
                        &[ion_card],
                        background_card,
                        || false
                    )
                    .is_ok()
            );
        },
    );
}

#[test]
fn scalar_electron_drag_cannot_invent_neutral_target_thermal_exchange() {
    let state = fixture::frozen_state(0.001, 8000.0, 8000.0, 1.0);
    let source = IonizedGasConductivity::from_state(
        state.clone(),
        fixture::frozen_rates(&state),
        fixture::limits(),
        || false,
    )
    .unwrap();
    let harmonic = source.prepare_harmonic(|| false).unwrap();
    let transport = harmonic
        .prepare_magnetized(1e10, [0.0; 3], magnetic_limits(), || false)
        .unwrap();
    let (background, _) = background(&source);
    let interface = MagnetizedIonizedGasInterface::try_new(
        DielectricResponse::Lossless(&background),
        &transport,
        &background,
        [0.0, 0.0, 1.0],
        || false,
    )
    .unwrap();
    let wave = interface.prepare_wave([C64::ZERO; 3], || false).unwrap();
    let slab = wave.slab_absorption(1.0, 0.0, 1.0).unwrap();
    assert!(matches!(
        MagnetizedIonizedGasHeating::try_new(slab, || false),
        Err(IonizedWaveHeatingError::Invalid(_))
    ));
}
