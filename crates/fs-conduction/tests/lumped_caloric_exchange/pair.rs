//! Reciprocal contact must retain each original caloric coordinate and amount.
//!
//! References below are independent finite-body solutions. No fixed bath or
//! representative heat capacity replaces either participating cell.

use super::*;
use fs_conduction::interface::{
    ResistanceUncertainty, SeriesThermalResistance, ThermalResistanceTerm,
};
use fs_conduction::lumped::exchange::caloric::pair::{
    LumpedCaloricPairExchange, LumpedCaloricPairExchangeConfig,
    solve_lumped_caloric_pair_exchange as solve_pair,
};
use fs_material::thermal::transport::{
    conduction::pair::FiniteCaloricExchangeControls,
    volumetric::{
        PreparedVolumetricThermalTransport, VOLUMETRIC_HEAT_CAPACITY_DIMS,
        VolumetricThermalPathQuery, timeline::FixedVolumeThermalLaw,
    },
};

fn resistance(value: f64) -> SeriesThermalResistance {
    SeriesThermalResistance::new(vec![
        ThermalResistanceTerm::declared(
            "finite native contact",
            value,
            ResistanceUncertainty::Unstated,
            "synthetic fixed contact resistance",
        )
        .unwrap(),
    ])
    .unwrap()
}

fn pair_config(duration_s: f64) -> LumpedCaloricPairExchangeConfig {
    LumpedCaloricPairExchangeConfig {
        duration_s,
        characteristic_lengths_m: [0.001; 2],
        contact_areas_m2: [0.25; 2],
        energy_tolerance_j: 1e-7,
        integration: FiniteCaloricExchangeControls::default(),
    }
}

fn pair_cell(id: u64, specific_energy: f64, mass: f64, density: f64) -> DepositedCell {
    let mut body = cell(specific_energy, 300.0);
    body.id = CellId(id);
    body.mass_kg = mass;
    body.geometry.volume_m3 = mass / density;
    body.boundary.htc_w_m2_k = 0.0;
    if id == 43 {
        body.geometry.start_m[0] += 0.02;
        body.geometry.end_m[0] += 0.02;
    }
    body
}

fn add_capacity(
    source: &mut Fixture,
    isochoric: bool,
    key: &PropertyKey,
    knots: Vec<(f64, f64)>,
) {
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature, ValueForm::Static,
    ));
    let (axis, quantity, fixed) = if isochoric {
        ("rho", QuantitySpec::dimensional(DENSITY), RHO)
    } else {
        ("p", QuantitySpec::dimensional(Pressure::DIMS), 100000.0)
    };
    let validity = ValidityDomain::unconstrained()
        .with_quantity("T", absolute, 300.0, 700.0)
        .with_quantity(axis, quantity, fixed, fixed);
    let mut claims = source.card.claims().clone();
    claims.insert_claim(PropertyClaim {
        key: key.clone(),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots,
            dims: key.quantity().dims(),
        },
        validity,
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic independent finite-body capacity".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }).unwrap();
    source.card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "finite-body capacity fixture".into(),
            phase: "single phase".into(),
            process: "independent caloric comparison".into(),
            revision: 0,
        },
        claims,
        vec![],
    ).unwrap();
}

fn constant_fixture(isochoric: bool, capacity: f64) -> Fixture {
    let mut source = Fixture::new(isochoric);
    let key = PropertyKey::new("constant_original_capacity", SPECIFIC_HEAT);
    add_capacity(
        &mut source, isochoric, &key,
        vec![(300.0, capacity), (700.0, capacity)],
    );
    source.capacity = key;
    source
}

fn volume_fixture() -> (Fixture, PropertyKey) {
    let mut source = Fixture::new(false);
    let key = PropertyKey::new("original_volume_capacity", VOLUMETRIC_HEAT_CAPACITY_DIMS);
    add_capacity(
        &mut source, false, &key,
        vec![
            (300.0, 750.0 * 700.0),
            (475.0, 750.0 * 875.0),
            (700.0, 750.0 * 1100.0),
        ],
    );
    (source, key)
}

fn assert_accounts(
    result: &LumpedCaloricPairExchange,
    first: &DepositedCell,
    second: &DepositedCell,
    tolerance: f64,
) {
    let q = result.caloric.heat_first_to_second_j;
    for (side, cell) in [first, second].into_iter().enumerate() {
        let heat = if side == 0 { -q } else { q };
        let stored = cell.mass_kg
            * (result.caloric.final_specific_energies_j_kg[side] - cell.specific_h_j_kg);
        near(result.caloric.represented_energy_change_j[side], stored, tolerance);
        near(stored, heat, tolerance);
        near(result.caloric.energy_residual_j[side], stored - heat, tolerance);
        assert!(result.caloric.energy_residual_j[side].abs() <= tolerance);
    }
    near(
        result.caloric.represented_energy_change_j.iter().sum(),
        0.0,
        tolerance,
    );
    assert!(result.total_energy_residual_j.abs() <= tolerance);
}

// Both total capacities are K*(T+400). Writing y=T+400 preserves
// S=y_a^2+y_b^2; with a=sqrt(S/2) and s=y_a+y_b, exact G*dt is
// K/2 * [a*atanh(s/(2a))-s]_initial^final.
fn affine_pair_oracle(total_capacity_slope: f64, conductance: f64) -> (f64, [f64; 2]) {
    let invariant = 1050.0_f64.powi(2) + 750.0_f64.powi(2);
    let target = [540.0, (invariant - 940.0_f64.powi(2)).sqrt() - 400.0];
    let a = (invariant / 2.0).sqrt();
    let initial_sum = 1800.0;
    let final_sum = target[0] + target[1] + 800.0;
    let exposure = 0.5 * total_capacity_slope
        * (a * ((final_sum / (2.0 * a)).atanh() - (initial_sum / (2.0 * a)).atanh())
            - (final_sum - initial_sum));
    (exposure / conductance, target)
}

#[test]
fn g1_constant_unequal_native_capacities_follow_the_reciprocal_exponential() {
    let (cv_source, cp_source) = (constant_fixture(true, 350.0), constant_fixture(false, 700.0));
    let (cv, cp) = (cv_source.cv(), cp_source.cp());
    let cv_k = IsochoricConductionLaw::try_new(
        &cv_source.card, &cv, &cv_source.conductivity, SELECTION,
    ).unwrap();
    let cp_k = IsobaricConductionLaw::try_new(
        &cp_source.card, &cp, &cp_source.conductivity, SELECTION,
    ).unwrap();
    let cv_layer = CaloricConductionLayer::from_isochoric(&cv_k).unwrap();
    let cp_layer = CaloricConductionLayer::from_isobaric(&cp_k).unwrap();
    let first = pair_cell(42, 350.0 * (650.0 - 300.0), 0.4, RHO);
    let second = pair_cell(43, 700.0 * (350.0 - 300.0), 0.6, RHO);
    let originals = [first.clone(), second.clone()];
    let r = resistance(0.5);
    // Ca=140 J/K, Cb=420 J/K, equilibrium425 K and tau=52.5 s.
    with_cx(|cx| {
        for duration in [0.0, 52.5 * 4.0_f64.ln(), 1e6] {
            let cfg = pair_config(duration);
            let result = solve_pair(
                cx, &first, &cv_layer, &second, &cp_layer, &r, GATE, cfg,
            ).unwrap();
            let remaining = (-duration / 52.5).exp();
            let heat = -31_500.0 * (-duration / 52.5).exp_m1();
            near(result.caloric.temperatures_k[0], 425.0 + 225.0 * remaining, 2e-10);
            near(result.caloric.temperatures_k[1], 425.0 - 75.0 * remaining, 2e-10);
            near(result.caloric.heat_first_to_second_j, heat, 3e-8);
            assert_eq!(
                result.caloric.storage_bases,
                [ThermalStorageBasis::InternalEnergy, ThermalStorageBasis::Enthalpy],
            );
            assert_eq!(result.cell_ids, [42, 43]);
            assert_eq!(result.resistance_k_per_w, 0.5);
            assert_eq!(
                result.caloric_identities,
                [
                    cv_layer.caloric_law().caloric_identity().unwrap(),
                    cp_layer.caloric_law().caloric_identity().unwrap(),
                ],
            );
            assert!(result.caloric.temperatures_k[0] >= 425.0 - 1e-10);
            assert!(result.caloric.temperatures_k[1] <= 425.0 + 1e-10);
            assert_accounts(&result, &first, &second, cfg.energy_tolerance_j);
        }
        // Different native h/u values at the same T are already at equilibrium.
        let mut a = first.clone();
        let mut b = second.clone();
        a.specific_h_j_kg = 350.0 * 200.0;
        b.specific_h_j_kg = 700.0 * 200.0;
        let steady = solve_pair(
            cx, &a, &cv_layer, &b, &cp_layer, &r, GATE, pair_config(100.0),
        ).unwrap();
        assert_eq!(steady.caloric.heat_first_to_second_j, 0.0);
        assert_eq!(steady.caloric.final_specific_energies_j_kg, [a.specific_h_j_kg, b.specific_h_j_kg]);
    });
    assert_eq!([first, second], originals);
}

#[test]
fn g1_original_variable_cv_and_cp_match_closed_form_reversal_and_composition() {
    let (cv_source, cp_source) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (cv_source.cv(), cp_source.cp());
    let cv_k = IsochoricConductionLaw::try_new(
        &cv_source.card, &cv, &cv_source.conductivity, SELECTION,
    ).unwrap();
    let cp_k = IsobaricConductionLaw::try_new(
        &cp_source.card, &cp, &cp_source.conductivity, SELECTION,
    ).unwrap();
    let cv_layer = CaloricConductionLayer::from_isochoric(&cv_k).unwrap();
    let cp_layer = CaloricConductionLayer::from_isobaric(&cp_k).unwrap();
    let first = pair_cell(42, energy(true, 650.0), 0.8, RHO);
    let second = pair_cell(43, energy(false, 350.0), 0.4, RHO);
    let (duration, target) = affine_pair_oracle(0.4, 2.0);
    let r = resistance(0.5);
    let cfg = LumpedCaloricPairExchangeConfig {
        characteristic_lengths_m: [0.002, 0.001],
        contact_areas_m2: [0.125, 0.0625],
        ..pair_config(duration)
    };
    with_cx(|cx| {
        let result = solve_pair(
            cx, &first, &cv_layer, &second, &cp_layer, &r, GATE, cfg,
        ).unwrap();
        near(result.caloric.heat_first_to_second_j, 43_780.0, 3e-5);
        for side in 0..2 {
            near(result.caloric.temperatures_k[side], target[side], 2e-7);
            near(
                result.caloric.final_specific_energies_j_kg[side],
                energy(side == 0, target[side]),
                2e-4,
            );
        }
        assert_accounts(&result, &first, &second, cfg.energy_tolerance_j);
        near(result.biot_numbers[0], 2.0 * 0.002 / (0.125 * 12.0), 1e-15);
        near(result.biot_numbers[1], 2.0 * 0.001 / (0.0625 * 6.0), 1e-15);

        // Initial-capacity exponentials conserve a different caloric model.
        let frozen_a = 525.0 + 125.0 * (-duration / 87.5).exp();
        assert!((frozen_a - result.caloric.temperatures_k[0]).abs() > 5.0);

        let reversed_cfg = LumpedCaloricPairExchangeConfig {
            characteristic_lengths_m: [cfg.characteristic_lengths_m[1], cfg.characteristic_lengths_m[0]],
            contact_areas_m2: [cfg.contact_areas_m2[1], cfg.contact_areas_m2[0]],
            ..cfg
        };
        let reversed = solve_pair(
            cx, &second, &cp_layer, &first, &cv_layer, &r, GATE, reversed_cfg,
        ).unwrap();
        near(reversed.caloric.heat_first_to_second_j, -result.caloric.heat_first_to_second_j, 3e-5);
        for side in 0..2 {
            near(reversed.caloric.temperatures_k[side], result.caloric.temperatures_k[1 - side], 2e-7);
            assert_eq!(reversed.biot_numbers[side], result.biot_numbers[1 - side]);
        }
        assert_eq!(
            reversed.caloric.storage_bases,
            [ThermalStorageBasis::Enthalpy, ThermalStorageBasis::InternalEnergy],
        );
        assert_accounts(&reversed, &second, &first, cfg.energy_tolerance_j);

        let first_half = solve_pair(
            cx, &first, &cv_layer, &second, &cp_layer, &r, GATE,
            LumpedCaloricPairExchangeConfig { duration_s: duration * 0.5, ..cfg },
        ).unwrap();
        let mut a = first.clone();
        let mut b = second.clone();
        a.specific_h_j_kg = first_half.caloric.final_specific_energies_j_kg[0];
        b.specific_h_j_kg = first_half.caloric.final_specific_energies_j_kg[1];
        let second_half = solve_pair(
            cx, &a, &cv_layer, &b, &cp_layer, &r, GATE,
            LumpedCaloricPairExchangeConfig { duration_s: duration * 0.5, ..cfg },
        ).unwrap();
        for side in 0..2 {
            near(second_half.caloric.temperatures_k[side], target[side], 3e-7);
        }
        near(
            first_half.caloric.heat_first_to_second_j + second_half.caloric.heat_first_to_second_j,
            result.caloric.heat_first_to_second_j,
            5e-5,
        );
    });
}

#[test]
fn g1_reference_volume_storage_uses_original_density_and_finite_source_support() {
    let (source, volume_key) = volume_fixture();
    let path = PreparedVolumetricThermalTransport::try_new(
        &source.card,
        VolumetricThermalPathQuery {
            volumetric_heat_capacity: &volume_key,
            conductivity: &source.conductivity,
            temperature_axis: "T",
            lower: &source.lower,
            upper: &source.upper,
            selections: [SELECTION; 2],
        },
    ).unwrap();
    let volume_law = FixedVolumeThermalLaw::try_new(&path, 750.0, 300.0).unwrap();
    let volume_layer = CaloricConductionLayer::from_fixed_volume(&volume_law).unwrap();
    let cv_source = Fixture::new(true);
    let cv = cv_source.cv();
    let cv_k = IsochoricConductionLaw::try_new(
        &cv_source.card, &cv, &cv_source.conductivity, SELECTION,
    ).unwrap();
    let cv_layer = CaloricConductionLayer::from_isochoric(&cv_k).unwrap();
    // Cvol/rho_ref=Cp numerically, but its declared storage basis stays distinct.
    let first = pair_cell(42, energy(false, 650.0), 0.4, 750.0);
    let second = pair_cell(43, energy(true, 350.0), 0.8, RHO);
    let (duration, target) = affine_pair_oracle(0.4, 2.0);
    let r = resistance(0.5);
    with_cx(|cx| {
        let cfg = pair_config(duration);
        let result = solve_pair(
            cx, &first, &volume_layer, &second, &cv_layer, &r, GATE, cfg,
        ).unwrap();
        assert_eq!(
            result.caloric.storage_bases,
            [ThermalStorageBasis::DeclaredSensibleStorage, ThermalStorageBasis::InternalEnergy],
        );
        near(result.caloric.heat_first_to_second_j, 43_780.0, 3e-5);
        for side in 0..2 {
            near(result.caloric.temperatures_k[side], target[side], 2e-7);
        }
        assert_accounts(&result, &first, &second, cfg.energy_tolerance_j);
        let mut wrong_volume = first.clone();
        wrong_volume.geometry.volume_m3 *= 1.01;
        assert!(solve_pair(
            cx, &wrong_volume, &volume_layer, &second, &cv_layer, &r, GATE, cfg,
        ).is_err());
    });

    // The material card extends to700 K, but this prepared original path was
    // admitted only through400 K. Finite contact may not extrapolate that path.
    let upper = QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static)),
            400.0,
        ).unwrap()
        .with_quantity("p", QuantitySpec::dimensional(Pressure::DIMS), 100000.0).unwrap();
    let bounded_path = PreparedVolumetricThermalTransport::try_new(
        &source.card,
        VolumetricThermalPathQuery {
            volumetric_heat_capacity: &volume_key,
            conductivity: &source.conductivity,
            temperature_axis: "T",
            lower: &source.lower,
            upper: &upper,
            selections: [SELECTION; 2],
        },
    ).unwrap();
    let bounded_law = FixedVolumeThermalLaw::try_new(&bounded_path, 750.0, 300.0).unwrap();
    let bounded_layer = CaloricConductionLayer::from_fixed_volume(&bounded_law).unwrap();
    let hot = pair_cell(42, energy(true, 650.0), 0.8, RHO);
    let cold = pair_cell(43, energy(false, 350.0), 0.4, 750.0);
    with_cx(|cx| {
        let short = solve_pair(
            cx, &hot, &cv_layer, &cold, &bounded_layer, &r, GATE, pair_config(1.0),
        ).unwrap();
        assert!((350.0..400.0).contains(&short.caloric.temperatures_k[1]));
        assert!(solve_pair(
            cx, &hot, &cv_layer, &cold, &bounded_layer, &r, GATE, pair_config(duration),
        ).is_err());
    });
}

#[test]
fn g3_mixed_mass_volume_inventory_scales_contact_joules_once() {
    let (source, volume_key) = volume_fixture();
    let cv_source = Fixture::new(true);
    let cv = cv_source.cv();
    let cv_k = IsochoricConductionLaw::try_new(
        &cv_source.card, &cv, &cv_source.conductivity, SELECTION,
    ).unwrap();
    let cv_layer = CaloricConductionLayer::from_isochoric(&cv_k).unwrap();
    for scale in [1.0, 4.0] {
        // .25 kg of Cp + .001 m3 of Cvol is exactly one kilogram's Cp.
        // Neither the volume component nor the reference mass is counted twice.
        let storage = PreparedHeatStorage::try_new(&[
            HeatStorageComponentQuery {
                id: 1, card: &source.card, capacity: &source.capacity,
                amount: HeatStorageAmount::MassKg(0.25 * scale), temperature_axis: "T",
                lower: &source.lower, upper: &source.upper, selection: SELECTION,
            },
            HeatStorageComponentQuery {
                id: 2, card: &source.card, capacity: &volume_key,
                amount: HeatStorageAmount::ReferenceVolumeM3(0.001 * scale), temperature_axis: "T",
                lower: &source.lower, upper: &source.upper, selection: SELECTION,
            },
        ]).unwrap();
        let caloric = HeatStorageThermalLaw::try_new(&storage, scale, 300.0).unwrap();
        let conduction = HeatStorageConductionLaw::try_new(
            &source.card, &caloric, &source.conductivity,
            "T", &source.lower, &source.upper, SELECTION,
        ).unwrap();
        let layer = CaloricConductionLayer::from_heat_storage(&conduction).unwrap();
        let mut first = pair_cell(42, energy(false, 650.0), scale, RHO);
        let mut second = pair_cell(43, energy(true, 350.0), 2.0 * scale, RHO);
        first.geometry.surface_area_m2 *= scale;
        second.geometry.surface_area_m2 *= scale;
        let conductance = 5.0 * scale;
        let r = resistance(1.0 / conductance);
        let (duration, target) = affine_pair_oracle(scale, conductance);
        let cfg = LumpedCaloricPairExchangeConfig {
            contact_areas_m2: [0.25 * scale; 2],
            ..pair_config(duration)
        };
        with_cx(|cx| {
            let result = solve_pair(
                cx, &first, &layer, &second, &cv_layer, &r, GATE, cfg,
            ).unwrap();
            assert_eq!(result.caloric.storage_bases[0], ThermalStorageBasis::DeclaredSensibleStorage);
            near(result.caloric.heat_first_to_second_j, 109_450.0 * scale, 6e-5 * scale);
            for side in 0..2 {
                near(result.caloric.temperatures_k[side], target[side], 2e-7);
            }
            assert_accounts(&result, &first, &second, cfg.energy_tolerance_j);
            let mut wrong_amount = first.clone();
            wrong_amount.mass_kg *= 1.01;
            wrong_amount.geometry.volume_m3 *= 1.01;
            assert!(solve_pair(
                cx, &wrong_amount, &layer, &second, &cv_layer, &r, GATE, cfg,
            ).is_err());
        });
    }
}

#[test]
fn g0_contact_keeps_original_inventory_isolation_biot_and_cancellation_admission() {
    let source = Fixture::new(true);
    let cv = source.cv();
    let k = IsochoricConductionLaw::try_new(
        &source.card, &cv, &source.conductivity, SELECTION,
    ).unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&k).unwrap();
    let first = pair_cell(42, energy(true, 650.0), 0.4, RHO);
    let second = pair_cell(43, energy(true, 350.0), 0.4, RHO);
    let originals = [first.clone(), second.clone()];
    let r = resistance(0.5);
    let cfg = pair_config(10.0);
    with_cx(|cx| {
        for case in 0..5 {
            let mut wrong = first.clone();
            match case {
                0 => wrong.mass_kg *= 1.1,
                1 => wrong.boundary.htc_w_m2_k = 1.0,
                2 => wrong.boundary.source_w = 1.0,
                3 => wrong.boundary.reservoir = Some(Reservoir {
                    area_m2: 0.01, conductance_w_k: 0.1, temperature_k: 500.0,
                }),
                _ => wrong.id = second.id,
            }
            assert!(solve_pair(
                cx, &wrong, &layer, &second, &layer, &r, GATE, cfg,
            ).is_err());
        }
        for refused in [
            LumpedCaloricPairExchangeConfig { contact_areas_m2: [0.251, 0.25], ..cfg },
            // Original k_min=12 gives Bi=.1333, although k(650)=26 would pass.
            LumpedCaloricPairExchangeConfig { characteristic_lengths_m: [0.2, 0.001], ..cfg },
        ] {
            assert!(solve_pair(
                cx, &first, &layer, &second, &layer, &r, GATE, refused,
            ).is_err());
        }
        let mut covered = first.clone();
        covered.boundary.htc_w_m2_k = 8.0;
        covered.boundary.reservoir = Some(Reservoir {
            area_m2: covered.geometry.surface_area_m2,
            conductance_w_k: 0.0,
            temperature_k: 900.0,
        });
        // An entirely covered inactive film supplies no external heat.
        let accepted = solve_pair(
            cx, &covered, &layer, &second, &layer, &r, GATE, cfg,
        ).unwrap();
        assert_accounts(&accepted, &covered, &second, cfg.energy_tolerance_j);
    });
    with_cancelled_cx(|cx| {
        assert!(matches!(
            solve_pair(cx, &first, &layer, &second, &layer, &r, GATE, cfg),
            Err(ConductionError::Cancelled { .. }),
        ));
    });
    assert_eq!([first, second], originals);
}
