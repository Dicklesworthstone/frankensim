//! Synthetic liquid laws test source/physics composition, not qualification of
//! any real saturated specimen. Thermal storage is the authored lumped fixture.

use super::*;
use fs_couple::porous_heating::feedback::covariance::{
    PorousCovarianceBin, PorousCrossSpectrumHeating,
};
use fs_couple::porous_heating::feedback::spectrum::{PorousSpectrumHeating, PorousSpectrumOptions};
use fs_duct::porous::PressureCovariance;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluid;
use fs_material::liquid::density::ResolvedLiquidDensityState;
use fs_material::liquid::{
    KINEMATIC_VISCOSITY_DIMS, KINEMATIC_VISCOSITY_PROPERTY, resolve_liquid_state_from_kinematic,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use std::cell::Cell;
use std::rc::Rc;

const RHO_P: f64 = 1.0 / 1_048_576.0;
const PRESSURE: f64 = 131_072.0;

fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", typed(QuantityKind::AbsoluteTemperature), t)
        .unwrap()
        .with_quantity("p", typed(QuantityKind::Pressure), PRESSURE)
        .unwrap()
        .with_quantity("salinity", QuantitySpec::dimensional(Dims([0; 6])), 0.25)
        .unwrap()
}
fn surface(dims: Dims, value: impl Fn(f64, f64) -> f64) -> PropertyValue {
    PropertyValue::Surface2D {
        axes: ["T".into(), "p".into()],
        axis_dims: [Temperature::DIMS, Pressure::DIMS],
        samples: [
            [256.0, 65_536.0],
            [384.0, 65_536.0],
            [384.0, 196_608.0],
            [256.0, 196_608.0],
        ]
        .map(|[t, p]| ([t, p], value(t, p)))
        .to_vec(),
        triangles: vec![[0, 1, 2], [0, 2, 3]],
        dims,
    }
}
fn liquid_card(rho_t: f64) -> MaterialCard {
    let properties = [
        (
            DENSITY_PROPERTY,
            QuantitySpec::dimensional(Density::DIMS),
            surface(Density::DIMS, |t, p| {
                1000.0 + rho_t * (t - 256.0) + RHO_P * (p - 65_536.0)
            }),
        ),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            QuantitySpec::dimensional(SPECIFIC_HEAT_CAPACITY_DIMS),
            surface(SPECIFIC_HEAT_CAPACITY_DIMS, |t, p| {
                4000.0 + 2.0 * (t - 256.0) + (p - 65_536.0) / 4096.0
            }),
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            typed(QuantityKind::ThermalConductivity),
            surface(THERMAL_CONDUCTIVITY_DIMS, |t, p| {
                0.5 + (t - 256.0) / 1024.0 + (p - 65_536.0) / 1_048_576.0
            }),
        ),
        (
            KINEMATIC_VISCOSITY_PROPERTY,
            QuantitySpec::dimensional(KINEMATIC_VISCOSITY_DIMS),
            surface(KINEMATIC_VISCOSITY_DIMS, |t, _| {
                1.0 / 1_048_576.0 - (t - 256.0) / 1_073_741_824.0
            }),
        ),
    ];
    let mut claims = ClaimSet::new();
    for (name, quantity, value) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value,
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", typed(QuantityKind::AbsoluteTemperature), 256.0, 384.0)
                    .with_quantity("p", typed(QuantityKind::Pressure), 65_536.0, 196_608.0)
                    .with_quantity(
                        "salinity",
                        QuantitySpec::dimensional(Dims([0; 6])),
                        0.0,
                        1.0,
                    ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic liquid-saturated absorption fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored liquid density fixture".into(),
            phase: "liquid".into(),
            process: "numerical coupling test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn reference(card: &MaterialCard) -> ResolvedLiquidDensityState {
    let liquid = resolve_liquid_state_from_kinematic(
        card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    ResolvedLiquidDensityState::try_from_liquid(card, liquid, "T", "p").unwrap()
}
fn liquid_drive(reference: &ResolvedLiquidDensityState, cosine: f64) -> (PorousHarmonicDrive, f64) {
    let state = reference.acoustic_state().unwrap();
    let mut result = drive();
    result.reference_impedance_pa_s_m3 =
        state.density_kg_m3() * state.sound_speed_m_s() / (AREA * cosine);
    result.incident_pressure_pa = [C64::from_re(20_000.0), C64::new(5000.0, -3000.0)];
    let tangential =
        result.angular_frequency_rad_s / state.sound_speed_m_s() * (1.0 - cosine * cosine).sqrt();
    (result, tangential)
}
fn liquid_source<'a>(
    card: &'a MaterialCard,
    reference: &'a ResolvedLiquidDensityState,
) -> IsobaricPorousLiquid<'a> {
    IsobaricPorousLiquid::try_new(card, reference, [280.0, 370.0]).unwrap()
}

#[test]
fn liquid_gradient_keeps_original_thermodynamics_and_coherent_oblique_loss() {
    let card = liquid_card(-0.25);
    let reference = reference(&card);
    let source = liquid_source(&card, &reference);
    let materials = materials();
    let initial = thermal(events([0.0, 40.0 * CP], [0.0; 2]), [&SENSIBLE; 2]);
    let (drive, kt) = liquid_drive(&reference, 0.8);
    let heating = TemperaturePorousHeating::try_new_oblique(
        &initial,
        &cells(&materials),
        source,
        drive,
        kt,
        options(),
        || Ok(()),
    )
    .unwrap();
    let original_states: Vec<_> = [300.0, 340.0]
        .map(|t| reference.at_query_point(&card, &point(t)).unwrap())
        .into_iter()
        .collect();
    let original_layers: Vec<_> = (0..2)
        .map(|i| {
            PorousLayer::try_new_oblique(
                &materials[i],
                &original_states[i],
                drive.angular_frequency_rad_s,
                [0.03, 0.05][i],
                AREA,
                drive.reference_impedance_pa_s_m3,
                kt,
            )
            .unwrap()
        })
        .collect();
    let stack = PorousLayerStack::try_new(&original_layers, options(), || false).unwrap();
    let excitation = stack
        .prepare_excitation(drive.incident_pressure_pa, || false)
        .unwrap();
    for (i, expected) in original_states.iter().enumerate() {
        let actual = heating.fluid_states()[i].source();
        assert_eq!(actual.material(), expected.material());
        assert_eq!(actual.response(), expected.response());
        assert_eq!(actual.material().query().axes()["p"], PRESSURE);
        assert_eq!(actual.material().query().axes()["salinity"], 0.25);
        let t = initial.temperature_k(id(i)).unwrap();
        let rho = 1000.0 - 0.25 * (t - 256.0) + RHO_P * (PRESSURE - 65_536.0);
        let cp = 4000.0 + 2.0 * (t - 256.0) + (PRESSURE - 65_536.0) / 4096.0;
        let cv = cp - t * 0.25_f64.powi(2) / (rho * rho * RHO_P);
        let state = actual.thermoviscous_state().unwrap();
        close(state.specific_heat_cv_j_kg_k(), cv, 2e-10);
        close(
            state.acoustic().sound_speed_m_s().powi(2),
            cp / (cv * RHO_P),
            2e-8,
        );
        close(
            state.transport().dynamic_viscosity_pa_s(),
            rho * (1.0 / 1_048_576.0 - (t - 256.0) / 1_073_741_824.0),
            2e-18,
        );
        assert_eq!(
            *heating.acoustic_layers()[i].material_response(),
            *original_layers[i].material_response(),
        );
        assert_eq!(
            heating.cells()[i].dissipation,
            excitation.layer_excitations()[i]
                .dissipation_between(0.0, original_layers[i].thickness_m())
                .unwrap(),
        );
    }
    assert_eq!(heating.dissipation(), excitation.dissipation());
    assert!(heating.dissipation().viscous_power_w > 0.0);
    assert!(heating.dissipation().thermal_power_w > 0.0);
    close(
        heating.power().incident_power_w - heating.power().outgoing_power_w,
        heating.dissipation().absorbed_power_w,
        2e-12,
    );
}

#[test]
fn liquid_accepted_start_heat_refreshes_response_without_replacing_native_mass() {
    let card = liquid_card(-0.25);
    let reference = reference(&card);
    let materials = materials();
    let authored = [0.3, 0.7];
    let initial = thermal(events([0.0, 30.0 * CP], authored), [&SENSIBLE; 2]);
    let (drive, kt) = liquid_drive(&reference, 0.8);
    let heating = TemperaturePorousHeating::try_new_oblique(
        &initial,
        &cells(&materials),
        liquid_source(&card, &reference),
        drive,
        kt,
        options(),
        || Ok(()),
    )
    .unwrap();
    let dt = 0.25;
    let result = heating
        .advance_to(&initial, dt, config(dt), || Ok(()))
        .unwrap();
    assert_eq!(result.work.fluid_evaluations, 2);
    assert_eq!(result.work.gas_evaluations, 0);
    assert_eq!(result.work.stack_solves, 1);
    close(
        result.step.acoustic.deposited_j,
        dt * heating.dissipation().absorbed_power_w,
        2e-14,
    );
    for i in 0..2 {
        let before = &initial.cells()[&id(i)];
        let after = &result.step.thermal.cells()[&id(i)];
        let heat = heating.cells()[i].dissipation.absorbed_power_w;
        close(
            after.specific_h_j_kg,
            before.specific_h_j_kg + dt * (heat + authored[i]) / before.mass_kg,
            2e-11,
        );
        assert_eq!(after.mass_kg, before.mass_kg);
        assert_eq!(after.geometry, before.geometry);
        assert_eq!(after.material_id, before.material_id);
        let temperature = result.step.thermal.temperature_k(id(i)).unwrap();
        assert_eq!(
            result.heating.fluid_states()[i]
                .source()
                .response()
                .temperature_k(),
            temperature
        );
        assert_ne!(
            result.heating.fluid_states()[i]
                .source()
                .liquid()
                .dynamic_viscosity_pa_s(),
            heating.fluid_states()[i]
                .source()
                .liquid()
                .dynamic_viscosity_pa_s(),
        );
    }
    assert_eq!(
        result.heating.drive().reference_impedance_pa_s_m3,
        drive.reference_impedance_pa_s_m3
    );
    assert_eq!(result.heating.fluid().pressure_pa(), PRESSURE);
    close(result.step.acoustic.residual_j, 0.0, 2e-12);
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn zero_expansion_liquid_has_exactly_zero_thermal_loss_and_reuses_plateau_cache() {
    let card = liquid_card(0.0);
    let reference = reference(&card);
    let materials = materials();
    let initial = thermal(events([0.0; 2], [0.0; 2]), [&PLATEAU; 2]);
    let (drive, kt) = liquid_drive(&reference, 0.7);
    let heating = TemperaturePorousHeating::try_new_oblique(
        &initial,
        &cells(&materials),
        liquid_source(&card, &reference),
        drive,
        kt,
        options(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(heating.dissipation().thermal_power_w, 0.0);
    assert!(heating.dissipation().viscous_power_w > 0.0);
    let result = heating
        .advance_to(&initial, 0.5, config(0.125), || Ok(()))
        .unwrap();
    assert_eq!(result.work, PorousFeedbackWork::default());
    assert_eq!(result.step.acoustic.thermal_j, 0.0);
    close(
        result.step.acoustic.deposited_j,
        0.5 * heating.dissipation().viscous_power_w,
        2e-14,
    );
    for i in 0..2 {
        assert!(std::ptr::eq(
            heating.fluid_states()[i].source(),
            result.heating.fluid_states()[i].source()
        ));
    }
}

#[derive(Clone)]
struct CountedLiquid<'a> {
    source: IsobaricPorousLiquid<'a>,
    evaluations: Rc<Cell<usize>>,
}
impl PorousTemperatureFluid for CountedLiquid<'_> {
    type State = PorousLiquidState;
    fn evaluate(&self, t: f64) -> Result<Self::State, PorousHeatingError> {
        self.evaluations.set(self.evaluations.get() + 1);
        self.source.evaluate(t)
    }
    fn temperature_k(&self, state: &Self::State) -> f64 {
        self.source.temperature_k(state)
    }
}

#[test]
fn liquid_spectral_and_covariance_modes_share_material_evaluations_and_native_heat() {
    let card = liquid_card(-0.25);
    let reference = reference(&card);
    let materials = materials();
    let initial = thermal(events([0.0, 30.0 * CP], [0.3, 0.7]), [&SENSIBLE; 2]);
    let evaluations = Rc::new(Cell::new(0));
    let source = CountedLiquid {
        source: liquid_source(&card, &reference),
        evaluations: Rc::clone(&evaluations),
    };
    let (base, base_kt) = liquid_drive(&reference, 0.8);
    let bins = [1800.0, 500.0, 1100.0].map(|hz| PorousHarmonicDrive {
        angular_frequency_rad_s: core::f64::consts::TAU * hz,
        ..base
    });
    let kt = bins.map(|bin| base_kt * bin.angular_frequency_rad_s / base.angular_frequency_rad_s);
    let options = PorousSpectrumOptions {
        maximum_bins: 3,
        maximum_cell_bins: 6,
        composition: options(),
    };
    let spectrum = PorousSpectrumHeating::try_new_oblique(
        &initial,
        &cells(&materials),
        source.clone(),
        &bins,
        &kt,
        options,
        || Ok(()),
    )
    .unwrap();
    assert_eq!(evaluations.get(), 2);
    for bin in spectrum.bins() {
        for i in 0..2 {
            assert!(std::ptr::eq(
                bin.fluid_states()[i].source(),
                spectrum.fluid_states()[i].source()
            ));
        }
    }
    let covariance_bins = bins.map(|bin| PorousCovarianceBin {
        angular_frequency_rad_s: bin.angular_frequency_rad_s,
        area_m2: bin.area_m2,
        reference_impedance_pa_s_m3: bin.reference_impedance_pa_s_m3,
        covariance: PressureCovariance::from_peak_pressure(bin.incident_pressure_pa).unwrap(),
    });
    let covariance = PorousCrossSpectrumHeating::try_new_oblique(
        &initial,
        &cells(&materials),
        source,
        &covariance_bins,
        &kt,
        options,
        || Ok(()),
    )
    .unwrap();
    assert_eq!(evaluations.get(), 4);
    evaluations.set(0);
    let spectral_step = spectrum
        .advance_to(&initial, 0.25, config(0.25), || Ok(()))
        .unwrap();
    assert_eq!(evaluations.get(), 2);
    let covariance_step = covariance
        .advance_to(&initial, 0.25, config(0.25), || Ok(()))
        .unwrap();
    assert_eq!(evaluations.get(), 4);
    for work in [spectral_step.work, covariance_step.work] {
        assert_eq!(work.fluid_evaluations, 2);
        assert_eq!(work.gas_evaluations, 0);
        assert_eq!(work.stack_solves, 3);
    }
    close(
        spectral_step.step.acoustic.deposited_j,
        covariance_step.step.acoustic.deposited_j,
        2e-12,
    );
    for i in 0..2 {
        close(
            spectral_step.step.thermal.cells()[&id(i)].specific_h_j_kg,
            covariance_step.step.thermal.cells()[&id(i)].specific_h_j_kg,
            2e-10,
        );
    }
}

#[test]
fn liquid_source_and_late_cancellation_refuse_without_publishing_either_candidate() {
    let card = liquid_card(-0.25);
    let reference = reference(&card);
    let wrong_card = liquid_card(-0.5);
    assert!(matches!(
        IsobaricPorousLiquid::try_new(&wrong_card, &reference, [280.0, 370.0]),
        Err(PorousHeatingError::LiquidSource { .. })
    ));
    let materials = materials();
    let initial = thermal(events([0.0; 2], [10.0; 2]), [&SENSIBLE; 2]);
    let (drive, _) = liquid_drive(&reference, 1.0);
    let bounded = IsobaricPorousLiquid::try_new(&card, &reference, [299.0, 300.001]).unwrap();
    let heating = TemperaturePorousHeating::try_new(
        &initial,
        &cells(&materials),
        bounded,
        drive,
        options(),
        || Ok(()),
    )
    .unwrap();
    match heating.advance_to(&initial, 0.25, config(0.25), || Ok(())) {
        Err(PorousHeatingError::LiquidSource {
            card_identity,
            reference_query,
            temperature_k,
            ..
        }) => {
            assert_eq!(card_identity, card.content_hash());
            assert_eq!(reference_query, point(300.0));
            assert!(temperature_k > 300.001);
        }
        _ => panic!("unsupported liquid endpoint must refuse with its original source context"),
    }
    let evaluations = Rc::new(Cell::new(0));
    let source = CountedLiquid {
        source: liquid_source(&card, &reference),
        evaluations: Rc::clone(&evaluations),
    };
    let admitted = TemperaturePorousHeating::try_new(
        &initial,
        &cells(&materials),
        source,
        drive,
        options(),
        || Ok(()),
    )
    .unwrap();
    evaluations.set(0);
    assert!(matches!(
        admitted.advance_to(&initial, 0.25, config(0.25), || {
            if evaluations.get() > 0 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(PorousHeatingError::Thermal(DepositionError::Cancelled)),
    ));
    assert_eq!(initial.time_s(), 0.0);
    for i in 0..2 {
        assert_eq!(initial.cells()[&id(i)].specific_h_j_kg, 0.0);
        assert_eq!(
            admitted.fluid_states()[i]
                .source()
                .response()
                .temperature_k(),
            300.0
        );
    }
    let retry = admitted
        .advance_to(&initial, 0.25, config(0.25), || Ok(()))
        .unwrap();
    assert!(retry.step.acoustic.deposited_j > 0.0);
    assert_eq!(retry.work.fluid_evaluations, 2);
}
