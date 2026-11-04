//! G1/G3/G4 thermal absorption through original gas and material source owners.
//! Synthetic source curves test equations, not a calibrated gas/solvent pair.
#[path = "../material_exchange/gas.rs"]
mod gas_fixture;

use fs_couple::cavitation::dissolved_gas::{
    IncompressibleSolvent, ThermalDissolutionControls, ThermalDissolutionError,
    advance_thermal_mixture,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::cavitation::dissolved_gas::{
    DIFFUSIVITY_DIMS, DISSOLVED_DIFFUSIVITY_PROPERTY, DissolutionFilm, DissolvedGasCalorics,
    DissolvedGasControls, DissolvedGasError, DissolvedGasInventory, GAS_FILM_DIFFUSIVITY_PROPERTY,
    GAS_FILM_MOLE_FRACTION_AXIS, GasFilmSource, HENRY_SOLUBILITY_DIMS, HENRY_SOLUBILITY_PROPERTY,
    HenryDissolution, HenrySource, MolarPermeance,
};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{
        GasEnergySolve,
        parcel::{ClosedGasParcel, transfer::species::boundary::GasSpeciesBoundarySolve},
    },
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec, Temperature};
use fs_thermochem::{
    Composition, CompositionBasis, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    mixture::FrozenIdealGasMixtureModelV1,
};

const T0: f64 = 300.0;
const P0: f64 = 100_000.0;
const VL: f64 = 0.002;
const VG: f64 = 0.001;
const M: f64 = 0.028;
const H: f64 = 2e-5;
const D: f64 = 2e-9;
const CP_SOLVENT: f64 = 1000.0;
const CP_DISSOLVED: f64 = 500.0;
const SOLUTION_HEAT: f64 = 6000.0;

fn q(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("temperature", q(QuantityKind::AbsoluteTemperature), t)
        .unwrap()
        .with_quantity("pressure", q(QuantityKind::Pressure), P0)
        .unwrap()
        .with_quantity(
            "dissolved_concentration",
            q(QuantityKind::AmountConcentration),
            0.0,
        )
        .unwrap()
}
fn claim(name: &str, dims: Dims, values: [f64; 2], maximum_temperature: f64) -> PropertyClaim {
    let constant = values[0] == values[1];
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: if constant {
            PropertyValue::Scalar {
                value: values[0],
                dims,
            }
        } else {
            PropertyValue::Curve {
                abscissa: "temperature".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(200.0, values[0]), (500.0, values[1])],
                dims,
            }
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                "temperature",
                q(QuantityKind::AbsoluteTemperature),
                200.0,
                maximum_temperature,
            )
            .with_quantity("pressure", q(QuantityKind::Pressure), 1000.0, 1e7)
            .with_quantity(
                "dissolved_concentration",
                q(QuantityKind::AmountConcentration),
                0.0,
                100.0,
            ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: if constant {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        },
        observations: vec![],
        provenance: Provenance {
            source: "synthetic finite thermal absorption".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn card(name: &str, rows: Vec<PropertyClaim>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for row in rows {
        claims.insert_claim(row).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "liquid".into(),
            process: "analytic fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn caloric(cp: f64) -> ResolvedHeatCapacityCurve {
    let dims = Dims([2, 0, -2, -1, 0, 0]);
    let source = card(
        "explicit pure-solvent or partial solute caloric chart",
        vec![
            claim("cp", dims, [cp; 2], 500.0),
            claim("rho", Density::DIMS, [1000.0; 2], 500.0),
        ],
    );
    ResolvedHeatCapacityCurve::try_new(
        &source,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("cp", dims),
            density: &PropertyKey::new("rho", Density::DIMS),
            temperature_axis: "temperature",
            lower: &point(200.0),
            upper: &point(500.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Liquid,
        0.0,
        1e-8,
    )
    .unwrap()
}
fn gas() -> ClosedGasParcel {
    gas_with_fractions([0.8, 0.2])
}
fn gas_with_fractions(fractions: [f64; 2]) -> ClosedGasParcel {
    let (a, at) = gas_fixture::component("A", M, 0.0);
    let (b, bt) = gas_fixture::component("B", M, 0.0);
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![a, b],
            Composition::new(CompositionBasis::MoleFraction, fractions.to_vec()).unwrap(),
        )
        .unwrap(),
        vec![at, bt],
        None,
    )
    .unwrap();
    ClosedGasParcel::new(model, point(T0), VG).unwrap()
}
fn law(
    gas: &ClosedGasParcel,
    name: &str,
    heat: f64,
    varying: bool,
    max_t: f64,
) -> HenryDissolution {
    let partial = caloric(CP_DISSOLVED);
    let gas_h = gas
        .model()
        .species_specific_enthalpy_j_kg(&id(name), T0)
        .unwrap();
    let dissolved_h = partial
        .state_at_temperature(T0)
        .unwrap()
        .specific_enthalpy_j_kg();
    HenryDissolution::new(
        gas.model(),
        id(name),
        HenrySource {
            material: card(
                "specified solute in one common solvent",
                vec![
                    claim(
                        HENRY_SOLUBILITY_PROPERTY,
                        HENRY_SOLUBILITY_DIMS,
                        if varying { [1.4 * H, 0.2 * H] } else { [H; 2] },
                        max_t,
                    ),
                    claim(
                        DISSOLVED_DIFFUSIVITY_PROPERTY,
                        DIFFUSIVITY_DIMS,
                        [D; 2],
                        max_t,
                    ),
                ],
            ),
            point: point(T0),
            selection: MaterialPropertySelection::SingleClaimOnly,
        },
        DissolvedGasCalorics {
            partial_enthalpy: partial,
            enthalpy_offset_j_kg: gas_h - dissolved_h - heat / M,
        },
    )
    .unwrap()
}
fn solvent() -> IncompressibleSolvent {
    IncompressibleSolvent::new(caloric(CP_SOLVENT), VL).unwrap()
}
fn film() -> DissolutionFilm {
    DissolutionFilm {
        area_m2: 0.01,
        thickness_m: 1e-5,
        gas_side_transfer: Some(MolarPermeance::try_new(H * D / 1e-5).unwrap()),
    }
}
fn sourced_gas_film(law: HenryDissolution, maximum_temperature: f64) -> (HenryDissolution, f64) {
    let diffusivity = 2e-10;
    let fraction = q(QuantityKind::Composition(
        fs_qty::semantic::CompositionBasis::MoleFraction,
    ));
    let mut property = claim(
        GAS_FILM_DIFFUSIVITY_PROPERTY,
        DIFFUSIVITY_DIMS,
        [diffusivity * 200.0 / T0, diffusivity * 500.0 / T0],
        maximum_temperature,
    );
    property.validity =
        property
            .validity
            .with_quantity(GAS_FILM_MOLE_FRACTION_AXIS, fraction, 0.0, 0.05);
    let source = GasFilmSource {
        species: id("A"),
        material: card("synthetic dilute A gas-film diffusivity", vec![property]),
        point: point(T0)
            .with_quantity(GAS_FILM_MOLE_FRACTION_AXIS, fraction, 0.01)
            .unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly,
        maximum_mole_fraction: 0.05,
    };
    let thickness = diffusivity * film().thickness_m / (H * R * T0 * D);
    (law.with_gas_film(source, thickness).unwrap(), thickness)
}
fn controls() -> ThermalDissolutionControls {
    ThermalDissolutionControls {
        exchange: DissolvedGasControls {
            gas: GasSpeciesBoundarySolve {
                gas: GasEnergySolve {
                    minimum_temperature_k: 200.0,
                    maximum_temperature_k: 500.0,
                    energy_tolerance_j_kg: 1e-7,
                    maximum_iterations: 96,
                },
                mass_tolerance_kg: 1e-14,
                energy_tolerance_j: 1e-8,
            },
            maximum_step_s: 1000.0,
            maximum_concentration_mol_m3: 20.0,
            maximum_film_relaxation_ratio: 0.01,
        },
        temperature_bracket_k: [250.0, 450.0],
        maximum_temperature_increment_k: 50.0,
        energy_tolerance_j: 1e-7,
        maximum_iterations: 96,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tol={tolerance:e}"
    );
}

#[test]
fn g1_absorption_and_exsolution_close_independent_finite_capacity_energy() {
    for (initial_moles, solution_heat) in [
        (0.0, SOLUTION_HEAT),
        (0.004, SOLUTION_HEAT),
        (0.0, -SOLUTION_HEAT),
    ] {
        let mut gas = gas();
        let carrier = gas.species_mass_kg(&id("B")).unwrap();
        let mut laws = [law(&gas, "A", solution_heat, false, 500.0)];
        let mut liquids = [DissolvedGasInventory::new(VL, initial_moles, 0.0).unwrap()];
        let films = [film()];
        let solvent = solvent();
        let dt = 300.0;
        let ng = P0 * VG / (R * T0);
        let capacitance = 1.0 / VL + H * R * T0 / VG;
        let conductance = films[0].area_m2
            / (films[0].thickness_m / D + H / films[0].gas_side_transfer.unwrap().value_si());
        let dn = (H * 0.8 * P0 - initial_moles / VL) / capacitance
            * (1.0 - (-conductance * capacitance * dt).exp());
        let capacity = solvent.mass_kg() * CP_SOLVENT
            + (ng - dn) * 2.5 * R
            + (initial_moles + dn) * M * CP_DISSOLVED;
        let expected_t = T0 + dn * (solution_heat - R * T0) / capacity;
        let report = advance_thermal_mixture(
            &mut laws,
            &mut gas,
            &mut liquids,
            &films,
            &solvent,
            dt,
            0.0,
            controls(),
            || false,
        )
        .unwrap();
        near(liquids[0].dissolved_moles(), initial_moles + dn, 2e-15);
        near(report.final_temperature_k, expected_t, 2e-10);
        near(
            gas.state().state().pressure,
            (ng - dn) * R * expected_t / VG,
            2e-7,
        );
        near(gas.species_mass_kg(&id("B")).unwrap(), carrier, 2e-17);
        near(
            gas.species_mass_kg(&id("A")).unwrap() / M + liquids[0].dissolved_moles(),
            0.8 * ng + initial_moles,
            2e-14,
        );
        near(report.energy_residual_j, 0.0, controls().energy_tolerance_j);
        assert_eq!(laws[0].temperature_k(), gas.state().state().temperature);
        assert_eq!(liquids[0].time_s(), dt);
        assert_eq!(report.endpoint_sources.len(), 1);
    }
}

#[test]
fn g1_sealed_mass_path_heats_native_inventories_and_refreshes_original_henry_source() {
    let mut gas = gas();
    let initial = gas.clone();
    let mut laws = [law(&gas, "A", SOLUTION_HEAT, true, 500.0)];
    let mut liquids = [DissolvedGasInventory::new(VL, 0.001, 0.0).unwrap()];
    let solvent = solvent();
    let mut sealed = film();
    sealed.area_m2 = 0.0;
    let heat = 1000.0;
    let ng = P0 * VG / (R * T0);
    let expected_t =
        T0 + heat / (solvent.mass_kg() * CP_SOLVENT + ng * 2.5 * R + 0.001 * M * CP_DISSOLVED);
    let report = advance_thermal_mixture(
        &mut laws,
        &mut gas,
        &mut liquids,
        &[sealed],
        &solvent,
        1.0,
        heat,
        controls(),
        || false,
    )
    .unwrap();
    near(report.final_temperature_k, expected_t, 1e-10);
    near(
        laws[0].henry_mol_m3_pa(),
        H * (1.0 - 0.004 * (expected_t - T0)),
        2e-20,
    );
    assert_eq!(liquids[0].dissolved_moles(), 0.001);
    assert_eq!(gas.mass_kg().to_bits(), initial.mass_kg().to_bits());
    assert_eq!(gas.model(), initial.model());
    let (at, du) = initial.at_temperature(350.0).unwrap();
    assert_eq!(at.mass_kg().to_bits(), initial.mass_kg().to_bits());
    assert_eq!(at.volume_m3(), initial.volume_m3());
    near(du, ng * 2.5 * R * 50.0, 1e-10);
    near(at.state().state().pressure, P0 * 350.0 / T0, 1e-8);
}

#[test]
fn g1_thermal_feedback_trajectory_refines_against_independent_energy_constrained_rk4() {
    // Independent single-coordinate oracle. Constant original capacities make
    // T(n,t) an explicit consequence of TOTAL U, while kinetics still depend
    // on current H(T), partial pressure and both film resistances.
    let total_gas = P0 * VG / (R * T0);
    let total_solute = 0.8 * total_gas;
    let end = 1000.0;
    let power = 20.0;
    let solvent = solvent();
    let films = [film()];
    let temperature = |n: f64, t: f64| {
        T0 + (power * t + n * (SOLUTION_HEAT - R * T0))
            / (solvent.mass_kg() * CP_SOLVENT + (total_gas - n) * 2.5 * R + n * M * CP_DISSOLVED)
    };
    let rate = |n: f64, t: f64| {
        let t_k = temperature(n, t);
        let h = H * (1.0 - 0.004 * (t_k - T0));
        let conductance = films[0].area_m2
            / (films[0].thickness_m / D + h / films[0].gas_side_transfer.unwrap().value_si());
        conductance * (h * (total_solute - n) * R * t_k / VG - n / VL)
    };
    let mut oracle = 0.0;
    let oracle_dt = end / 4096.0;
    for i in 0..4096 {
        let t = i as f64 * oracle_dt;
        let k1 = rate(oracle, t);
        let k2 = rate(oracle + 0.5 * oracle_dt * k1, t + 0.5 * oracle_dt);
        let k3 = rate(oracle + 0.5 * oracle_dt * k2, t + 0.5 * oracle_dt);
        let k4 = rate(oracle + oracle_dt * k3, t + oracle_dt);
        oracle += (oracle_dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }
    let mut errors = Vec::new();
    for steps in [16, 32, 64] {
        let mut gas = gas();
        let mut laws = [law(&gas, "A", SOLUTION_HEAT, true, 500.0)];
        let mut liquids = [DissolvedGasInventory::new(VL, 0.0, 0.0).unwrap()];
        let dt = end / steps as f64;
        for i in 1..=steps {
            let report = advance_thermal_mixture(
                &mut laws,
                &mut gas,
                &mut liquids,
                &films,
                &solvent,
                i as f64 * dt,
                power * dt,
                controls(),
                || false,
            )
            .unwrap();
            near(report.energy_residual_j, 0.0, controls().energy_tolerance_j);
        }
        let n = liquids[0].dissolved_moles();
        near(gas.state().state().temperature, temperature(n, end), 2e-8);
        errors.push((n - oracle).abs());
    }
    assert!(errors[0] > errors[1] && errors[1] > errors[2], "{errors:?}");
    for pair in errors.windows(2) {
        assert!((1.7..2.3).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn g3_opposite_species_share_one_solvent_heat_account_and_canonical_result() {
    let original = gas();
    let laws = [
        law(&original, "A", SOLUTION_HEAT, false, 500.0),
        law(&original, "B", 2.0 * SOLUTION_HEAT, false, 500.0),
    ];
    let liquids = [
        DissolvedGasInventory::new(VL, 0.0, 0.0).unwrap(),
        DissolvedGasInventory::new(VL, 0.003, 0.0).unwrap(),
    ];
    let mut first_gas = original.clone();
    let mut first_laws = laws.clone();
    let mut first_liquids = liquids;
    let mut second_gas = original;
    let mut second_laws = [laws[1].clone(), laws[0].clone()];
    let mut second_liquids = [liquids[1], liquids[0]];
    let a = advance_thermal_mixture(
        &mut first_laws,
        &mut first_gas,
        &mut first_liquids,
        &[film(); 2],
        &solvent(),
        300.0,
        0.0,
        controls(),
        || false,
    )
    .unwrap();
    let b = advance_thermal_mixture(
        &mut second_laws,
        &mut second_gas,
        &mut second_liquids,
        &[film(); 2],
        &solvent(),
        300.0,
        0.0,
        controls(),
        || false,
    )
    .unwrap();
    assert!(a.isothermal_exchange.exchanges[0].dissolved_moles_change > 0.0);
    assert!(a.isothermal_exchange.exchanges[1].dissolved_moles_change < 0.0);
    assert_eq!(a, b);
    assert_eq!(first_gas, second_gas);
    assert_eq!(first_liquids, [second_liquids[1], second_liquids[0]]);
}

#[test]
fn g4_late_source_or_cancellation_refusal_keeps_all_temperatures_sources_and_clocks() {
    let original_gas = gas();
    let original_laws = [law(&original_gas, "A", SOLUTION_HEAT, true, 500.0)];
    let original_liquids = [DissolvedGasInventory::new(VL, 0.0, 0.0).unwrap()];
    let mut calls = 0;
    let mut gas = original_gas.clone();
    let mut laws = original_laws.clone();
    let mut liquids = original_liquids;
    let expected = advance_thermal_mixture(
        &mut laws,
        &mut gas,
        &mut liquids,
        &[film()],
        &solvent(),
        300.0,
        10.0,
        controls(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    let mut trial = 0;
    let mut gas = original_gas.clone();
    let mut laws = original_laws.clone();
    let mut liquids = original_liquids;
    assert!(matches!(
        advance_thermal_mixture(
            &mut laws,
            &mut gas,
            &mut liquids,
            &[film()],
            &solvent(),
            300.0,
            10.0,
            controls(),
            || {
                trial += 1;
                trial == calls
            }
        ),
        Err(ThermalDissolutionError::Cancelled)
    ));
    assert_eq!(gas, original_gas);
    assert_eq!(laws, original_laws);
    assert_eq!(liquids, original_liquids);
    assert_eq!(
        advance_thermal_mixture(
            &mut laws,
            &mut gas,
            &mut liquids,
            &[film()],
            &solvent(),
            300.0,
            10.0,
            controls(),
            || false
        )
        .unwrap(),
        expected
    );
    let mut gas = original_gas.clone();
    let mut laws = [law(&gas, "A", SOLUTION_HEAT, false, 310.0)];
    let saved_laws = laws.clone();
    let mut liquids = original_liquids;
    assert!(
        advance_thermal_mixture(
            &mut laws,
            &mut gas,
            &mut liquids,
            &[film()],
            &solvent(),
            300.0,
            50_000.0,
            controls(),
            || false
        )
        .is_err()
    );
    assert_eq!(gas, original_gas);
    assert_eq!(laws, saved_laws);
    assert_eq!(liquids, original_liquids);
}

#[test]
fn g1_sourced_gas_film_continues_its_temperature_and_both_interface_fluxes() {
    let mut gas = gas_with_fractions([0.01, 0.99]);
    let (law, thickness) = sourced_gas_film(law(&gas, "A", SOLUTION_HEAT, true, 500.0), 500.0);
    let mut laws = [law];
    let mut liquids = [DissolvedGasInventory::new(VL, 0.0, 0.0).unwrap()];
    let film = DissolutionFilm {
        gas_side_transfer: None,
        ..film()
    };
    let report = advance_thermal_mixture(
        &mut laws,
        &mut gas,
        &mut liquids,
        &[film],
        &solvent(),
        300.0,
        20_000.0,
        controls(),
        || false,
    )
    .unwrap();
    let t = report.final_temperature_k;
    assert!(t > 309.0);
    near(
        laws[0].gas_film_diffusivity_m2_s().unwrap(),
        2e-10 * t / T0,
        1e-24,
    );
    assert_eq!(
        laws[0].gas_film_source().unwrap().point.axes()["temperature"],
        t
    );
    let sample = &report.endpoint_sources[0];
    near(
        sample.molar_flux_mol_m2_s,
        D / film.thickness_m
            * (sample.interface_concentration_mol_m3 - sample.concentration_mol_m3),
        1e-18,
    );
    near(
        sample.molar_flux_mol_m2_s,
        (2e-10 * t / T0) / (thickness * R * t)
            * (sample.partial_pressure_pa - sample.interface_partial_pressure_pa),
        1e-18,
    );
    for material in sample.gas_film_material.as_ref().unwrap() {
        assert!(
            material
                .query_point()
                .iter()
                .any(|(axis, value)| axis == "temperature" && *value == t)
        );
    }
    near(report.energy_residual_j, 0.0, controls().energy_tolerance_j);
}

#[test]
fn g4_thermal_continuation_keeps_gas_film_source_and_relaxation_refusals() {
    let original = gas_with_fractions([0.01, 0.99]);
    let film = DissolutionFilm {
        gas_side_transfer: None,
        ..film()
    };
    let liquid = DissolvedGasInventory::new(VL, 0.0, 0.0).unwrap();
    let (supported, _) = sourced_gas_film(law(&original, "A", SOLUTION_HEAT, false, 500.0), 500.0);
    let mut strict = controls().exchange;
    strict.maximum_film_relaxation_ratio = 1e-4;
    assert!(matches!(
        supported.at_conditions(&original, liquid, film, strict),
        Err(DissolvedGasError::GasFilmRelaxation { .. })
    ));

    let (bounded, _) = sourced_gas_film(law(&original, "A", SOLUTION_HEAT, false, 500.0), 305.0);
    let mut gas = original.clone();
    let mut laws = [bounded.clone()];
    let mut liquids = [liquid];
    assert!(matches!(
        advance_thermal_mixture(
            &mut laws,
            &mut gas,
            &mut liquids,
            &[film],
            &solvent(),
            300.0,
            20_000.0,
            controls(),
            || false,
        ),
        Err(ThermalDissolutionError::Dissolution(
            DissolvedGasError::Material(_)
        ))
    ));
    assert_eq!(gas, original);
    assert_eq!(laws, [bounded]);
    assert_eq!(liquids, [liquid]);
}
