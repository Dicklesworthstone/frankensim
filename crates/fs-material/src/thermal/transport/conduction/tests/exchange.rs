//! Exact native energy exchange, using independent affine-capacity oracles.

use super::*;
use crate::thermal::transport::{
    expansion::{
        FreeExpansionSourceQuery, FreeExpansionStrainSources, PreparedFreeExpansionDensity,
    },
    isobaric::expansion::{
        IsobaricFreeExpansion, conduction::IsobaricExpansionConductionLaw,
    },
    storage::{
        HeatStorageAmount, HeatStorageComponentQuery, HeatStorageThermalLaw, PreparedHeatStorage,
        conduction::HeatStorageConductionLaw,
    },
    volumetric::{
        PreparedVolumetricThermalTransport, VOLUMETRIC_HEAT_CAPACITY_DIMS,
        VolumetricThermalPathQuery, timeline::FixedVolumeThermalLaw,
    },
};
use fs_qty::Dims;

// Independent primitive of C(T)/(bath-T) for C(T)=a+b*T.
fn exposure(initial: f64, target: f64, bath: f64, a: f64, b: f64) -> f64 {
    (a + b * bath) * ((bath - initial) / (bath - target)).ln() - b * (target - initial)
}

#[test]
fn g1_original_cp_and_cv_exchange_match_signed_analytic_energy_and_mean() {
    let (cp_source, cv_source) = (Fixture::new(false), Fixture::new(true));
    let (cp, cv) = (cp_source.cp(), cv_source.cv());
    let cp_k = IsobaricConductionLaw::try_new(
        &cp_source.card, &cp, &cp_source.conductivity, SELECTION,
    ).unwrap();
    let cv_k = IsochoricConductionLaw::try_new(
        &cv_source.card, &cv, &cv_source.conductivity, SELECTION,
    ).unwrap();
    let cp_layer = CaloricConductionLayer::from_isobaric(&cp_k).unwrap();
    let cv_layer = CaloricConductionLayer::from_isochoric(&cv_k).unwrap();
    for (layer, a, b, side, basis) in [
        (&cp_layer, 400.0, 1.0, 1, ThermalStorageBasis::Enthalpy),
        (&cv_layer, 200.0, 0.5, 0, ThermalStorageBasis::InternalEnergy),
    ] {
        // Both paths cross original capacity and independent conductivity knots.
        for (initial, target, bath) in [(350.0, 650.0, 900.0), (650.0, 350.0, 100.0)] {
            let x = exposure(initial, target, bath, a, b);
            let step = layer.linear_heat_exchange(energy(side, initial), bath, x).unwrap();
            let q = (a + 0.5 * b * (initial + target)) * (target - initial);
            assert_eq!(step.storage_basis, basis);
            near(step.temperature_k, target, 2e-10);
            near(step.specific_heat_j_kg, q, 2e-7);
            near(step.final_specific_energy_j_kg, energy(side, target), 2e-7);
            near(step.mean_temperature_k, bath - q / x, 2e-10);
            near(step.inlet_energy_gain, (bath - target) / (bath - initial), 2e-12);
            near(step.inlet_energy_gain + step.inlet_energy_loss, 1.0, 2e-14);
            assert_eq!(
                step.heat_residual_j_kg,
                step.specific_heat_j_kg - step.represented_energy_change_j_kg
            );
            assert_eq!(
                layer.caloric_law().temperature_k(step.final_specific_energy_j_kg).unwrap(),
                step.temperature_k
            );
        }
    }
}

#[test]
fn g1_native_energy_keeps_heat_below_one_temperature_ulp_and_reports_storage_rounding() {
    let source = Fixture::new(false);
    let cp = source.cp();
    let conduction = IsobaricConductionLaw::try_new(
        &source.card, &cp, &source.conductivity, SELECTION,
    ).unwrap();
    let layer = CaloricConductionLayer::from_isobaric(&conduction).unwrap();
    let x = 1e-20;
    for initial in [300.0, 500.0] {
        let input = cp.specific_enthalpy_j_kg(initial).unwrap();
        let step = layer.linear_heat_exchange(input, 600.0, x).unwrap();
        let q = x * (600.0 - initial);
        near(step.specific_heat_j_kg, q, 8.0 * f64::EPSILON * q);
        assert_eq!(step.temperature_k, initial);
        assert_eq!(step.mean_temperature_k, initial);
        assert!(step.inlet_energy_loss > 0.0);
        assert!(step.mean_reservoir_derivative > 0.0);
        if input == 0.0 {
            assert!(step.final_specific_energy_j_kg > 0.0);
            assert_eq!(step.final_specific_energy_j_kg, step.specific_heat_j_kg);
            assert_eq!(step.heat_residual_j_kg, 0.0);
        } else {
            assert_eq!(step.final_specific_energy_j_kg, input);
            assert_eq!(step.represented_energy_change_j_kg, 0.0);
            assert_eq!(step.heat_residual_j_kg, step.specific_heat_j_kg);
        }
        // A temperature-only caller still cannot discard this unresolved heat.
        assert!(cp.linear_heat_exchange(initial, 600.0, x).is_err());
        assert!(cp.temperature_after_specific_heat(initial, q).is_err());
    }
}

#[test]
fn g1_exchange_retains_original_datum_and_composes_in_native_energy() {
    let source = Fixture::new(false);
    let first = source.cp();
    let second = IsobaricHeatCapacityLaw::try_new(
        &source.card,
        IsobaricHeatCapacityQuery {
            specific_heat_cp: &source.capacity, temperature_axis: "T", pressure_axis: "p",
            lower: &source.lower, upper: &source.upper, selection: SELECTION,
        },
        475.0,
    ).unwrap();
    let ka = IsobaricConductionLaw::try_new(
        &source.card, &first, &source.conductivity, SELECTION,
    ).unwrap();
    let kb = IsobaricConductionLaw::try_new(
        &source.card, &second, &source.conductivity, SELECTION,
    ).unwrap();
    let a = CaloricConductionLayer::from_isobaric(&ka).unwrap();
    let b = CaloricConductionLayer::from_isobaric(&kb).unwrap();
    let x = exposure(350.0, 650.0, 900.0, 400.0, 1.0);
    let whole = a.linear_heat_exchange(first.specific_enthalpy_j_kg(350.0).unwrap(), 900.0, x).unwrap();
    let shifted = b.linear_heat_exchange(second.specific_enthalpy_j_kg(350.0).unwrap(), 900.0, x).unwrap();
    near(shifted.temperature_k, whole.temperature_k, 2e-10);
    near(shifted.specific_heat_j_kg, whole.specific_heat_j_kg, 2e-7);
    near(
        shifted.final_specific_energy_j_kg,
        second.specific_enthalpy_j_kg(650.0).unwrap(),
        2e-7,
    );
    let part = a.linear_heat_exchange(whole.initial_specific_energy_j_kg, 900.0, x / 3.0).unwrap();
    let rest = a.linear_heat_exchange(part.final_specific_energy_j_kg, 900.0, 2.0 * x / 3.0).unwrap();
    near(rest.final_specific_energy_j_kg, whole.final_specific_energy_j_kg, 2e-7);
    near(part.specific_heat_j_kg + rest.specific_heat_j_kg, whole.specific_heat_j_kg, 2e-7);
}

#[test]
fn g3_exchange_tangents_match_native_energy_and_bath_perturbations() {
    let source = Fixture::new(false);
    let cp = source.cp();
    let conduction = IsobaricConductionLaw::try_new(
        &source.card, &cp, &source.conductivity, SELECTION,
    ).unwrap();
    let layer = CaloricConductionLayer::from_isobaric(&conduction).unwrap();
    let (input, bath, x) = (energy(1, 500.0), 660.0, 250.0);
    let step = layer.linear_heat_exchange(input, bath, x).unwrap();
    let [lo, hi] = [-0.1, 0.1].map(|d| layer.linear_heat_exchange(input + d, bath, x).unwrap());
    near((hi.final_specific_energy_j_kg - lo.final_specific_energy_j_kg) / 0.2, step.inlet_energy_gain, 1e-7);
    let [lo, hi] = [-0.001, 0.001].map(|d| layer.linear_heat_exchange(input, bath + d, x).unwrap());
    near(
        (hi.specific_heat_j_kg - lo.specific_heat_j_kg) / 0.002,
        step.reservoir_energy_derivative_j_kg_k, 2e-5,
    );
    near((hi.mean_temperature_k - lo.mean_temperature_k) / 0.002, step.mean_reservoir_derivative, 2e-7);
    let [lo, hi] = [-0.01, 0.01].map(|d| layer.linear_heat_exchange(input, bath, x + d).unwrap());
    near(x * (hi.mean_temperature_k - lo.mean_temperature_k) / 0.02, step.mean_log_exposure_derivative_k, 2e-6);
    let unchanged = layer.linear_heat_exchange(input, bath, 0.0).unwrap();
    assert_eq!(unchanged.final_specific_energy_j_kg, input);
    assert_eq!(unchanged.specific_heat_j_kg, 0.0);
    assert_eq!(unchanged.inlet_energy_gain, 1.0);
    assert_eq!(unchanged.inlet_energy_loss, 0.0);
    assert_eq!(unchanged.reservoir_energy_derivative_j_kg_k, 0.0);
}

// Extend the same synthetic Cp/k card with independently declared volume,
// density and free-strain sources. None are inferred from its heat capacity.
fn additional_sources() -> (Fixture, [PropertyKey; 3]) {
    let mut source = Fixture::new(false);
    let keys = [
        PropertyKey::new("native_volumetric_capacity", VOLUMETRIC_HEAT_CAPACITY_DIMS),
        PropertyKey::new("reference_density", DENSITY),
        PropertyKey::new("isotropic_free_strain", Dims::NONE),
    ];
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature, ValueForm::Static,
    ));
    let validity = ValidityDomain::unconstrained()
        .with_quantity("T", absolute, 300.0, 700.0)
        .with_quantity("p", QuantitySpec::dimensional(Pressure::DIMS), 100000.0, 100000.0);
    let mut claims = source.card.claims().clone();
    for (key, knots) in [
        (&keys[0], vec![(300.0, 750.0 * 700.0), (475.0, 750.0 * 875.0), (700.0, 750.0 * 1100.0)]),
        (&keys[1], vec![(300.0, 750.0), (700.0, 750.0)]),
        (&keys[2], vec![(300.0, 0.0), (700.0, 0.02)]),
    ] {
        claims.insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots, dims: key.quantity().dims(),
            },
            validity: validity.clone(), interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated, observations: vec![],
            provenance: Provenance {
                source: "synthetic independent caloric and expansion comparison".into(),
                license: "test-only".into(), artifact: None,
            },
        }).unwrap();
    }
    source.card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "equivalent declared storage fixtures".into(),
            phase: "single phase".into(), process: "numerical comparison".into(), revision: 0,
        },
        claims, vec![],
    ).unwrap();
    (source, keys)
}

#[test]
fn g1_volume_assembly_and_free_expansion_keep_original_capacity_and_energy_basis() {
    let (source, keys) = additional_sources();
    let cp = source.cp();
    let cp_k = IsobaricConductionLaw::try_new(
        &source.card, &cp, &source.conductivity, SELECTION,
    ).unwrap();
    let cp_layer = CaloricConductionLayer::from_isobaric(&cp_k).unwrap();
    let volume = PreparedVolumetricThermalTransport::try_new(
        &source.card,
        VolumetricThermalPathQuery {
            volumetric_heat_capacity: &keys[0], conductivity: &source.conductivity,
            temperature_axis: "T", lower: &source.lower, upper: &source.upper,
            selections: [SELECTION; 2],
        },
    ).unwrap();
    let volume_law = FixedVolumeThermalLaw::try_new(&volume, 750.0, 300.0).unwrap();
    let volume_layer = CaloricConductionLayer::from_fixed_volume(&volume_law).unwrap();
    // 1.5 kg of Cp plus 0.002 m3 of Cvol gives C_total=3*Cp on a 3 kg basis.
    let storage = PreparedHeatStorage::try_new(&[
        HeatStorageComponentQuery {
            id: 1, card: &source.card, capacity: &source.capacity,
            amount: HeatStorageAmount::MassKg(1.5), temperature_axis: "T",
            lower: &source.lower, upper: &source.upper, selection: SELECTION,
        },
        HeatStorageComponentQuery {
            id: 2, card: &source.card, capacity: &keys[0],
            amount: HeatStorageAmount::ReferenceVolumeM3(0.002), temperature_axis: "T",
            lower: &source.lower, upper: &source.upper, selection: SELECTION,
        },
    ]).unwrap();
    let assembled = HeatStorageThermalLaw::try_new(&storage, 3.0, 300.0).unwrap();
    let assembled_k = HeatStorageConductionLaw::try_new(
        &source.card, &assembled, &source.conductivity, "T", &source.lower, &source.upper, SELECTION,
    ).unwrap();
    let assembled_layer = CaloricConductionLayer::from_heat_storage(&assembled_k).unwrap();
    let density = PreparedFreeExpansionDensity::try_new(
        &source.card,
        FreeExpansionSourceQuery {
            density: &keys[1], density_selection: SELECTION, reference: &source.lower,
            strains: FreeExpansionStrainSources::Isotropic { strain: &keys[2], selection: SELECTION },
            temperature_axis: "T", lower: &source.lower, upper: &source.upper,
        },
    ).unwrap();
    let expansion = IsobaricFreeExpansion::try_new(&source.card, &cp, &density).unwrap();
    let expanding_k = IsobaricExpansionConductionLaw::try_new(&expansion, &cp_k).unwrap();
    let expanding_layer = CaloricConductionLayer::from_expanding_isobaric(&expanding_k).unwrap();
    for (initial, target, bath) in [(350.0, 650.0, 900.0), (650.0, 350.0, 100.0)] {
        let input = energy(1, initial);
        let x = exposure(initial, target, bath, 400.0, 1.0);
        let expected = cp_layer.linear_heat_exchange(input, bath, x).unwrap();
        for (layer, basis) in [
            (&volume_layer, ThermalStorageBasis::DeclaredSensibleStorage),
            (&assembled_layer, ThermalStorageBasis::DeclaredSensibleStorage),
            (&expanding_layer, ThermalStorageBasis::Enthalpy),
        ] {
            let step = layer.linear_heat_exchange(input, bath, x).unwrap();
            assert_eq!(step.storage_basis, basis);
            near(step.temperature_k, target, 2e-10);
            near(step.specific_heat_j_kg, expected.specific_heat_j_kg, 2e-7);
            near(step.final_specific_energy_j_kg, expected.final_specific_energy_j_kg, 2e-7);
            near(step.mean_temperature_k, expected.mean_temperature_k, 2e-10);
            near(step.reservoir_energy_derivative_j_kg_k, expected.reservoir_energy_derivative_j_kg_k, 2e-8);
            near(step.inlet_energy_gain, expected.inlet_energy_gain, 2e-12);
            // Source power and conductance exposure convert in the SAME original
            // mass/volume basis, including the derivative with respect to Q0.
            let drive = energy(1, target) - input;
            for feedback in [0.0, 100.0] {
                let driven = layer.driven_heat_exchange(input, drive, feedback).unwrap();
                let reference = cp_layer.driven_heat_exchange(input, drive, feedback).unwrap();
                assert_eq!(driven.storage_basis, basis);
                near(driven.temperature_k, reference.temperature_k, 2e-10);
                near(driven.specific_heat_j_kg, reference.specific_heat_j_kg, 2e-7);
                near(driven.mean_temperature_change_k, reference.mean_temperature_change_k, 2e-10);
                near(driven.source_heat_gain, reference.source_heat_gain, 2e-12);
                near(driven.mean_source_derivative_kg_k_j, reference.mean_source_derivative_kg_k_j, 2e-14);
                if feedback == 0.0 {
                    assert_eq!(driven.specific_heat_j_kg, drive);
                    assert_eq!(driven.source_heat_gain, 1.0);
                }
            }
            // A native energy increment at the datum also survives the basis conversion.
            let weak = layer.linear_heat_exchange(0.0, 600.0, 1e-20).unwrap();
            assert_eq!(weak.temperature_k, 300.0);
            assert!(weak.final_specific_energy_j_kg > 0.0);
            assert_eq!(weak.heat_residual_j_kg, 0.0);
        }
    }
    assert!(
        expansion.bulk_density_kg_m3(energy(1, 650.0)).unwrap().unwrap()
            < expansion.bulk_density_kg_m3(energy(1, 350.0)).unwrap().unwrap()
    );
}

#[test]
fn g4_exchange_keeps_canonical_endpoints_and_refuses_unsupported_paths() {
    let source = Fixture::new(false);
    let cp = source.cp();
    let conduction = IsobaricConductionLaw::try_new(
        &source.card, &cp, &source.conductivity, SELECTION,
    ).unwrap();
    let layer = CaloricConductionLayer::from_isobaric(&conduction).unwrap();
    let input = energy(1, 500.0);
    let bounds = layer.energy_bounds_j_kg();
    for (bath, target) in [(300.0, bounds[0]), (700.0, bounds[1])] {
        let step = layer.linear_heat_exchange(input, bath, 1e12).unwrap();
        assert_eq!(step.final_specific_energy_j_kg, target);
        assert_eq!(step.temperature_k, bath);
        assert_eq!(step.heat_residual_j_kg, step.specific_heat_j_kg - (target - input));
    }
    // A bath outside support is allowed only while the complete material path is supported.
    assert!(layer.linear_heat_exchange(input, 900.0, 1.0).is_ok());
    assert!(layer.linear_heat_exchange(input, 900.0, 1e6).is_err());
    for x in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(layer.linear_heat_exchange(input, 600.0, x).is_err());
    }
    for bath in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(layer.linear_heat_exchange(input, bath, 1.0).is_err());
    }
    for e in [bounds[0].next_down(), bounds[1].next_up(), f64::NAN, f64::INFINITY] {
        assert!(layer.linear_heat_exchange(e, 600.0, 0.0).is_err());
    }
}


#[test]
fn g1_native_cp_cv_source_exchange_retains_original_energy_and_time_mean() {
    let (cp_source, cv_source) = (Fixture::new(false), Fixture::new(true));
    let (cp, cv) = (cp_source.cp(), cv_source.cv());
    let cp_k = IsobaricConductionLaw::try_new(
        &cp_source.card, &cp, &cp_source.conductivity, SELECTION,
    ).unwrap();
    let cv_k = IsochoricConductionLaw::try_new(
        &cv_source.card, &cv, &cv_source.conductivity, SELECTION,
    ).unwrap();
    let cp_layer = CaloricConductionLayer::from_isobaric(&cp_k).unwrap();
    let cv_layer = CaloricConductionLayer::from_isochoric(&cv_k).unwrap();
    for (layer, side, basis) in [
        (&cp_layer, 1, ThermalStorageBasis::Enthalpy),
        (&cv_layer, 0, ThermalStorageBasis::InternalEnergy),
    ] {
        for (initial, target) in [(350.0, 650.0), (650.0, 350.0)] {
            let input = energy(side, initial);
            let drive = energy(side, target) - input;
            let response = layer.driven_heat_exchange(input, drive, 0.0).unwrap();
            assert_eq!(response.storage_basis, basis);
            assert_eq!(response.specific_heat_j_kg, drive);
            assert_eq!(response.source_heat_gain, 1.0);
            near(response.temperature_k, target, 2e-10);
            near(response.final_specific_energy_j_kg, energy(side, target), 2e-7);
            near(response.mean_temperature_k, 508.3333333333333, 2e-10);
            near(response.mean_source_derivative_kg_k_j, (target - response.mean_temperature_k) / drive, 2e-14);
            assert_eq!(response.heat_residual_j_kg,
                response.specific_heat_j_kg - response.represented_energy_change_j_kg);
        }
        let weak = layer.driven_heat_exchange(0.0, 1e-18, 5.0).unwrap();
        assert_eq!(weak.temperature_k, 300.0);
        assert_eq!(weak.mean_temperature_k, 300.0);
        assert!(weak.temperature_change_k > 0.0 && weak.mean_temperature_change_k > 0.0);
        assert!(weak.final_specific_energy_j_kg > 0.0);
        assert_eq!(weak.final_specific_energy_j_kg, weak.specific_heat_j_kg);
        assert_eq!(weak.heat_residual_j_kg, 0.0);
    }
}

