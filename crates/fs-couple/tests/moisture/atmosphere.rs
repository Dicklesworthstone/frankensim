//! Independent finite-reservoir limits and native material/plate propagation.
//! All species, sorption and solid coefficients here are synthetic test data.
use super::*;
use fs_couple::moisture::atmosphere::*;
use fs_material::thermochemical_gas::{energy::parcel::ClosedGasParcel, packs::from_packs};
use fs_qty::chemistry::SpeciesId;
use fs_thermochem::{Composition, CompositionBasis, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

#[path = "../../src/render/schedule/reed/field_input/tests/support.rs"]
mod gas_support;

#[path = "atmosphere/shared.rs"]
mod shared;

const TEMPERATURE: f64 = 293.15;
const PRESSURE: f64 = 101_325.0;
const WATER_MOLAR_MASS: f64 = 0.018_015_28;
const CARRIER_MOLAR_MASS: f64 = 0.028;
const SORPTION_PROPERTY: &str = "equilibrium_dry_basis_moisture";

fn saturation_pressure() -> f64 {
    fs_material::gas::saturation_pressure_water_pa(TEMPERATURE).unwrap()
}
pub(super) fn atmosphere(rh: f64, volume: f64) -> FiniteMoistureAtmosphere {
    let carrier = gas_support::packs("carrier", CARRIER_MOLAR_MASS, 3.5, 200.0);
    let water = gas_support::packs("water", WATER_MOLAR_MASS, 4.5, 200.0);
    let water_fraction = rh * saturation_pressure() / PRESSURE;
    let model = from_packs(
        &[carrier.binding(), water.binding()],
        Composition::new(
            CompositionBasis::MoleFraction,
            vec![1.0 - water_fraction, water_fraction],
        )
        .unwrap(),
        Some(&SpeciesId::new("water").unwrap()),
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            qty(QuantityKind::AbsoluteTemperature),
            TEMPERATURE,
        )
        .unwrap()
        .with_quantity("pressure", qty(QuantityKind::Pressure), PRESSURE)
        .unwrap();
    FiniteMoistureAtmosphere::new(ClosedGasParcel::new(model, point, volume).unwrap()).unwrap()
}
fn sorption_point(rh: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            qty(QuantityKind::AbsoluteTemperature),
            TEMPERATURE,
        )
        .unwrap()
        .with_quantity(
            "relative_humidity",
            qty(QuantityKind::Moisture(MoistureBasis::RelativeHumidity)),
            rh,
        )
        .unwrap()
}
fn sorption_card(slope: f64, output: MoistureBasis, upper: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::with_quantity(SORPTION_PROPERTY, qty(QuantityKind::Moisture(output))),
            value: PropertyValue::Curve {
                abscissa: "relative_humidity".into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, 0.0), (1.0, slope)],
                dims: Dims::NONE,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "relative_humidity",
                    qty(QuantityKind::Moisture(MoistureBasis::RelativeHumidity)),
                    0.0,
                    upper,
                )
                .with_quantity(
                    "temperature",
                    qty(QuantityKind::AbsoluteTemperature),
                    TEMPERATURE,
                    TEMPERATURE,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "authored linear dry-basis isotherm M_eq=a*RH".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored sorption boundary".into(),
            phase: "solid".into(),
            process: "verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
pub(super) fn sorption(slope: f64) -> DryBasisSorption {
    DryBasisSorption::new(
        sorption_card(slope, MoistureBasis::DryMass, 1.0),
        sorption_point(0.5),
        SORPTION_PROPERTY,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
pub(super) fn controls() -> AtmosphereMoistureSolve {
    AtmosphereMoistureSolve {
        relative_humidity_bracket: [0.0, 1.0],
        maximum_iterations: 80,
        mass_tolerance_kg: 1e-13,
    }
}
fn carrier_mass(atmosphere: &FiniteMoistureAtmosphere) -> f64 {
    let gas = atmosphere.gas().state().state();
    gas.pressure
        * (1.0 - gas.water_mole_fraction)
        * atmosphere.gas().volume_m3()
        * CARRIER_MOLAR_MASS
        / (R * TEMPERATURE)
}

#[test]
fn one_cell_finite_water_balance_matches_adsorption_desorption_and_continuation() {
    with_cx(|cx| {
        let slope = 0.25;
        let boundaries = [Some(sorption(slope)), Some(sorption(slope))];
        let source_identity = sorption_card(slope, MoistureBasis::DryMass, 1.0).content_hash();
        let diffusion = 0.01;
        let dt = 2.0;
        let alpha = 4.0 * dt * diffusion; // Unit diffusion length, two half-cell faces.
        for initial_rh in [0.8, 0.2] {
            let mut atmosphere = atmosphere(initial_rh, 1.0);
            let mut moisture = slab(card([diffusion; 2], 0.3), 1.0, vec![0.1], cx);
            let dry_mass = moisture.dry_mass_kg();
            let saturation = saturation_pressure() * WATER_MOLAR_MASS / (R * TEMPERATURE);
            let original_carrier = carrier_mass(&atmosphere);
            let packs = atmosphere.gas().model().source_packs().to_vec();
            let mut previous_transport = None;
            for interval in 0..2 {
                let r0 = atmosphere.relative_humidity().unwrap();
                let m0 = moisture.moisture()[0];
                let expected_rh = (saturation * r0 + dry_mass * alpha * m0 / (1.0 + alpha))
                    / (saturation + dry_mass * alpha * slope / (1.0 + alpha));
                let expected_m = (m0 + alpha * slope * expected_rh) / (1.0 + alpha);
                let prior_gas = *atmosphere.gas().state().state();
                let initial_mass = atmosphere.gas().mass_kg() + moisture.wet_mass_kg();
                let candidate = atmosphere
                    .advance(&moisture, dt, &boundaries, config(), controls(), cx)
                    .unwrap();
                let next_gas = candidate.atmosphere().gas().state().state();
                let report = candidate.step();
                near(
                    candidate.atmosphere().relative_humidity().unwrap(),
                    expected_rh,
                    2e-11,
                );
                near(candidate.moisture().moisture()[0], expected_m, 2e-12);
                near(
                    next_gas.pressure,
                    prior_gas.pressure + saturation_pressure() * (expected_rh - r0),
                    1e-7,
                );
                near(
                    report.water_to_solid_kg,
                    dry_mass * (expected_m - m0),
                    2e-13,
                );
                assert_eq!(
                    report.water_to_solid_kg.is_sign_positive(),
                    initial_rh > 0.5
                );
                near(report.water_balance_residual_kg, 0.0, 1e-13);
                // The authored water NASA-9 law has cp/R=4.5+0.001*T and zero
                // reference offset. At fixed T/carrier inventory, only vapor
                // removal changes gas storage; this is not sorption heat.
                let vapor_u =
                    R * (3.5 * TEMPERATURE + 0.0005 * TEMPERATURE * TEMPERATURE) / WATER_MOLAR_MASS;
                near(
                    report.gas_internal_energy_change_j,
                    -report.water_to_solid_kg * vapor_u,
                    5e-7,
                );
                near(
                    carrier_mass(candidate.atmosphere()),
                    original_carrier,
                    2e-14,
                );
                near(
                    candidate.atmosphere().gas().mass_kg() + candidate.moisture().wet_mass_kg(),
                    initial_mass,
                    2e-13,
                );
                near(
                    report.transport.boundary_water_kg.iter().sum(),
                    report.water_to_solid_kg,
                    2e-13,
                );
                assert_eq!(next_gas.temperature, TEMPERATURE);
                assert_ne!(next_gas.sound_speed, prior_gas.sound_speed);
                assert_ne!(next_gas.dynamic_viscosity, prior_gas.dynamic_viscosity);
                assert_ne!(
                    next_gas.thermal_conductivity,
                    prior_gas.thermal_conductivity
                );
                assert_eq!(
                    candidate.atmosphere().gas().state().source_packs(),
                    packs.as_slice()
                );
                for resolved in report.sorption.iter().flatten() {
                    assert_eq!(resolved.card_identity(), source_identity);
                    near(
                        resolved.property(SORPTION_PROPERTY).unwrap().value_si(),
                        slope * expected_rh,
                        1e-11,
                    );
                    assert_eq!(
                        resolved
                            .property(SORPTION_PROPERTY)
                            .unwrap()
                            .requirement()
                            .quantity(),
                        qty(QuantityKind::Moisture(MoistureBasis::DryMass))
                    );
                }
                if let Some(previous) = &previous_transport {
                    assert_eq!(&report.transport.transport_start, previous);
                }
                previous_transport = Some(report.transport.transport_end.clone());
                assert_eq!(candidate.moisture().time_s(), (interval + 1) as f64 * dt);
                let (next_atmosphere, next_moisture, _) = candidate.into_parts();
                atmosphere = next_atmosphere;
                moisture = next_moisture;
            }
        }
    });
}

#[test]
fn finite_atmosphere_and_diffusivity_data_change_native_plate_stiffness_mass_and_modes() {
    with_cx(|cx| {
        let boundaries = [Some(sorption(0.25)), Some(sorption(0.25))];
        let mut wet_masses = Vec::new();
        let mut eigenvalues = Vec::new();
        for (diffusion, volume) in [(0.01, 0.2), (0.01, 2.0), (0.04, 2.0)] {
            let initial = slab(card([diffusion; 2], 0.3), 1.0, vec![0.1; 4], cx);
            let binding = plate(&initial, PlateMaterialModel::Isotropic, cx);
            let dry = binding.resolve(&initial, cx).unwrap();
            let before = modes(dry.model(), (20.0, 90.0), &SliceOptions::default()).unwrap();
            assert_eq!(before.modes.len(), 1);
            let atmosphere = atmosphere(0.8, volume);
            let (candidate, plate) = binding
                .advance_finite_atmosphere(
                    &atmosphere,
                    &initial,
                    2.0,
                    &boundaries,
                    config(),
                    controls(),
                    cx,
                )
                .unwrap();
            let after = modes(plate.model(), (15.0, 90.0), &SliceOptions::default()).unwrap();
            assert_eq!(after.modes.len(), 1);
            assert!(after.modes[0].lambda < before.modes[0].lambda);
            assert_eq!(plate.time_s(), candidate.moisture().time_s());
            near(
                plate.wet_mass_kg(),
                candidate.moisture().wet_mass_kg(),
                1e-13,
            );
            near(
                plate.wet_mass_kg() - dry.wet_mass_kg(),
                candidate.step().water_to_solid_kg,
                1e-13,
            );
            for (source, &moisture) in plate
                .elastic_states()
                .iter()
                .zip(candidate.moisture().moisture())
            {
                near(
                    source.property("young_modulus").unwrap().value_si(),
                    1e9 * (1.0 - moisture),
                    1e-6,
                );
            }
            assert!(
                plate.elastic_states()[0]
                    .property("young_modulus")
                    .unwrap()
                    .value_si()
                    < plate.elastic_states()[1]
                        .property("young_modulus")
                        .unwrap()
                        .value_si()
            );
            wet_masses.push(plate.wet_mass_kg());
            eigenvalues.push(after.modes[0].lambda);
        }
        assert!(wet_masses[1] > wet_masses[0]);
        assert!(wet_masses[2] > wet_masses[1]);
        assert!(eigenvalues[1] < eigenvalues[0]);
        assert!(eigenvalues[2] < eigenvalues[1]);
    });
}

#[test]
fn sealed_diffusion_and_exact_sorption_equilibrium_retain_the_original_gas() {
    with_cx(|cx| {
        let atmosphere = atmosphere(0.6, 1.0);
        let initial = slab(card([0.01; 2], 0.3), 1.0, vec![0.2, 0.1, 0.05, 0.05], cx);
        let sealed = atmosphere
            .advance(&initial, 1.0, &[None, None], config(), controls(), cx)
            .unwrap();
        assert_eq!(sealed.atmosphere(), &atmosphere);
        assert_eq!(
            sealed.step().vapor_mass_kg[0],
            sealed.step().vapor_mass_kg[1]
        );
        assert_eq!(sealed.step().gas_internal_energy_change_j, 0.0);
        assert_eq!(sealed.step().transport.boundary_water_kg, [0.0; 2]);
        near(
            sealed.moisture().water_mass_kg(),
            initial.water_mass_kg(),
            1e-13,
        );
        assert!(sealed.moisture().moisture()[0] < initial.moisture()[0]);
        assert!(sealed.moisture().moisture()[2] > initial.moisture()[2]);

        let slope = 0.25;
        let m = slope * atmosphere.relative_humidity().unwrap();
        let uniform = slab(card([0.01; 2], 0.3), 1.0, vec![m], cx);
        let boundary = [Some(sorption(slope)), Some(sorption(slope))];
        let unchanged = atmosphere
            .advance(&uniform, 1.0, &boundary, config(), controls(), cx)
            .unwrap();
        assert_eq!(unchanged.atmosphere(), &atmosphere);
        assert_eq!(unchanged.step().gas_internal_energy_change_j, 0.0);
        near(unchanged.moisture().moisture()[0], m, 1e-14);
        near(unchanged.step().water_to_solid_kg, 0.0, 1e-14);
        assert_eq!(unchanged.moisture().time_s(), 1.0);
    });
}

#[test]
fn source_semantics_cancellation_and_late_mechanical_refusal_keep_all_accepted_owners() {
    let valid = sorption_card(0.25, MoistureBasis::DryMass, 1.0);
    let alias = sorption_point(0.5)
        .with_quantity(
            "humidity",
            qty(QuantityKind::Moisture(MoistureBasis::RelativeHumidity)),
            0.1,
        )
        .unwrap();
    let duplicate_temperature = sorption_point(0.5).with("T", 290.0).unwrap();
    let wrong_coordinate = sorption_point(0.5)
        .with_quantity(
            "relative_humidity",
            qty(QuantityKind::Moisture(MoistureBasis::DryMass)),
            0.5,
        )
        .unwrap();
    for point in [alias, duplicate_temperature, wrong_coordinate] {
        assert!(
            DryBasisSorption::new(
                valid.clone(),
                point,
                SORPTION_PROPERTY,
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
    }
    assert!(
        DryBasisSorption::new(
            sorption_card(0.25, MoistureBasis::WetMass, 1.0),
            sorption_point(0.5),
            SORPTION_PROPERTY,
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );

    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let atmosphere = atmosphere(0.8, 2.0);
        let saved = atmosphere.clone();
        let initial = slab(card([0.01; 2], 0.105), 1.0, vec![0.1; 4], &cx);
        let binding = plate(&initial, PlateMaterialModel::Isotropic, &cx);
        let boundaries = [Some(sorption(0.25)), Some(sorption(0.25))];
        let transported = atmosphere
            .advance(&initial, 2.0, &boundaries, config(), controls(), &cx)
            .unwrap();
        assert!(transported.moisture().moisture()[0] > 0.105);
        assert!(
            binding
                .advance_finite_atmosphere(
                    &atmosphere,
                    &initial,
                    2.0,
                    &boundaries,
                    config(),
                    controls(),
                    &cx
                )
                .is_err()
        );
        // The original native diffusion accepts this water endpoint; only the
        // actual mechanical source domain rejects the coupled publication.
        let unsupported = DryBasisSorption::new(
            sorption_card(0.25, MoistureBasis::DryMass, 0.6),
            sorption_point(0.5),
            SORPTION_PROPERTY,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert!(
            atmosphere
                .advance(
                    &initial,
                    2.0,
                    &[Some(unsupported), None],
                    config(),
                    controls(),
                    &cx
                )
                .is_err()
        );
        assert_eq!(atmosphere, saved);
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(initial.moisture(), &[0.1; 4]);
        binding.resolve(&initial, &cx).unwrap();
        gate.request();
        assert!(matches!(
            atmosphere.advance(&initial, 2.0, &boundaries, config(), controls(), &cx),
            Err(AtmosphereMoistureError::Moisture(MoistureError::Cancelled))
        ));
        assert!(
            binding
                .advance_finite_atmosphere(
                    &atmosphere,
                    &initial,
                    2.0,
                    &boundaries,
                    config(),
                    controls(),
                    &cx
                )
                .is_err()
        );
        assert_eq!(atmosphere, saved);
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(initial.moisture(), &[0.1; 4]);
    });
}
