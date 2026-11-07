//! Original affine Cp/rho sources and immutable prepared liquid endpoints.
use super::*;
use crate::liquid_temperature::{
    LiquidViscosityBasis,
    tests::{card, point},
};
use crate::phase::{HeatCapacityKnot, HeatCapacitySourceQuery, ResolvedHeatCapacityCurve};
use crate::state_point::MaterialPropertySelection;
use fs_matdb::{MaterialCard, PropertyKey, SelectionPolicy};
use fs_qty::{Density, Dims};

fn original_curve(card: &MaterialCard, datum: f64) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new_exact(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", Dims([2, 0, -2, -1, 0, 0])),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(),
            upper: &point().with("T", 400.0).unwrap(),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Liquid,
        datum,
    )
    .unwrap()
}

fn law(card: &MaterialCard, basis: LiquidViscosityBasis) -> LiquidTemperatureLaw<'_> {
    LiquidTemperatureLaw::try_new(
        card,
        &point(),
        "T",
        basis,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn reanchoring_retains_original_density_expansion_path_and_future_source_updates() {
    let card = card(LiquidViscosityBasis::Dynamic, false);
    let original = law(&card, LiquidViscosityBasis::Dynamic)
        .with_density_expansion(Temperature::new(300.0), Temperature::new(400.0))
        .unwrap();
    let warm = original
        .reanchored_at_temperature(Temperature::new(350.0))
        .unwrap();
    assert!(core::ptr::eq(warm.card(), &card));
    assert_eq!(warm.reference_point(), &point().with("T", 350.0).unwrap());
    assert_eq!(warm.reference_state().density_kg_m3(), 900.0);
    near(
        warm.reference_state().volumetric_expansion_per_k().unwrap(),
        2.0 / 900.0,
        1e-17,
    );
    let hot = warm.at_temperature(Temperature::new(375.0)).unwrap();
    near(
        hot.volumetric_expansion_per_k().unwrap(),
        2.0 / 850.0,
        1e-17,
    );
    assert_eq!(
        hot,
        original.at_temperature(Temperature::new(375.0)).unwrap()
    );
    assert_eq!(
        warm.reference_state(),
        &original.at_temperature(Temperature::new(350.0)).unwrap()
    );
    assert_eq!(original.reference_state().density_kg_m3(), 1000.0);
    assert_eq!(original.reference_point(), &point());
}

#[test]
fn original_affine_cp_heat_resolves_one_consistent_owned_density_endpoint() {
    for basis in [
        LiquidViscosityBasis::Dynamic,
        LiquidViscosityBasis::Kinematic,
    ] {
        let card = card(basis, false);
        let curve = original_curve(&card, 1_000_000_000.0);
        let initial = curve
            .curve()
            .state_at_specific_enthalpy(1_000_000_000.0)
            .unwrap();
        // The endpoint borrows the card, not the old law. Cp(T)=2000+20*(T-300)
        // gives dh=125000 at T=350, not the constant-initial-Cp value T=362.5.
        let endpoint = {
            let original = law(&card, basis)
                .with_density_expansion(Temperature::new(300.0), Temperature::new(400.0))
                .unwrap();
            original
                .prepare_caloric_endpoint(curve.curve(), initial, 125_000.0, 1e-10)
                .unwrap()
        };
        assert_eq!(endpoint.initial(), initial);
        assert_eq!(
            endpoint.state().phase_curve_identity(),
            curve.curve().identity()
        );
        assert_eq!(endpoint.state().temperature_k(), 350.0);
        assert_eq!(endpoint.state().specific_enthalpy_j_kg(), 1_000_125_000.0);
        assert_eq!(endpoint.state().bulk_density_kg_m3(), 900.0);
        assert_eq!(endpoint.source().density_kg_m3(), 900.0);
        assert_eq!(endpoint.source().viscosity_basis(), basis);
        assert_eq!(
            endpoint.source().capillarity(),
            LiquidCapillarity::Neglected
        );
        assert_eq!(
            endpoint.source().material(),
            endpoint.law().reference_state().material()
        );
        assert_eq!(endpoint.net_specific_enthalpy_increment_j_kg(), 125_000.0);
        assert_eq!(
            endpoint.represented_specific_enthalpy_increment_j_kg(),
            125_000.0
        );
        assert_eq!(endpoint.maximum_source_density_defect_kg_m3(), 0.0);
        near(
            endpoint
                .law()
                .reference_state()
                .volumetric_expansion_per_k()
                .unwrap(),
            2.0 / 900.0,
            1e-17,
        );
        near(
            endpoint.source().dynamic_viscosity_pa_s(),
            if basis == LiquidViscosityBasis::Dynamic {
                0.0015
            } else {
                900.0 * 1.5e-6
            },
            1e-17,
        );
        for receipt in curve.receipts() {
            assert_eq!(
                endpoint
                    .source()
                    .material()
                    .property(&receipt.property)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
                receipt.selected,
            );
        }
    }
}

#[test]
fn heating_then_cooling_keeps_original_datum_sources_and_zero_step() {
    let card = card(LiquidViscosityBasis::Dynamic, false);
    let curve = original_curve(&card, -50_000.0);
    let original = law(&card, LiquidViscosityBasis::Dynamic);
    let initial = curve.curve().state_at_specific_enthalpy(-50_000.0).unwrap();
    let warm = original
        .prepare_caloric_endpoint(curve.curve(), initial, 125_000.0, 1e-10)
        .unwrap();
    let cool = warm
        .law()
        .prepare_caloric_endpoint(curve.curve(), warm.state(), -104_000.0, 1e-10)
        .unwrap();
    near(cool.state().temperature_k(), 310.0, 1e-12);
    near(cool.state().bulk_density_kg_m3(), 980.0, 1e-10);
    assert_eq!(cool.state().specific_enthalpy_j_kg(), -29_000.0);
    assert_eq!(
        cool.state().phase_curve_identity(),
        initial.phase_curve_identity()
    );
    assert_eq!(
        cool.law().reference_point(),
        &point().with("T", 310.0).unwrap()
    );
    let unchanged = cool
        .law()
        .prepare_caloric_endpoint(curve.curve(), cool.state(), 0.0, 1e-10)
        .unwrap();
    assert_eq!(unchanged.state(), cool.state());
    assert_eq!(unchanged.source(), cool.source());
    assert_eq!(unchanged.net_specific_enthalpy_increment_j_kg(), 0.0);
    assert_eq!(
        unchanged.represented_specific_enthalpy_increment_j_kg(),
        0.0
    );
}

#[test]
fn original_enthalpy_already_accounts_for_declared_ambient_pressure_work() {
    let card = card(LiquidViscosityBasis::Dynamic, false);
    let curve = original_curve(&card, 0.0);
    let original = law(&card, LiquidViscosityBasis::Dynamic);
    let initial = curve.curve().state_at_specific_enthalpy(0.0).unwrap();
    let endpoint = original
        .prepare_caloric_endpoint(curve.curve(), initial, 125_000.0, 1e-10)
        .unwrap();
    let heat = curve
        .curve()
        .advance_isobaric_heat(initial, 2.0, 100_000.0, 250_000.0)
        .unwrap();
    assert_eq!(endpoint.state(), heat.state);
    near(heat.volume_change_m3, 2.0 / 900.0 - 2.0 / 1000.0, 1e-18);
    assert!(heat.boundary_work_j > 0.0);
    near(
        heat.internal_energy_change_j + heat.boundary_work_j,
        250_000.0,
        1e-10,
    );
    assert_eq!(
        2.0 * endpoint.represented_specific_enthalpy_increment_j_kg(),
        250_000.0
    );
    assert_eq!(
        endpoint
            .source()
            .material()
            .query_point()
            .iter()
            .find(|(name, _)| name == "pressure")
            .unwrap()
            .1,
        100_000.0
    );
}

#[test]
fn accepted_initial_temperature_need_not_equal_original_reference_temperature() {
    let card = card(LiquidViscosityBasis::Dynamic, false);
    let curve = original_curve(&card, 0.0);
    let original = law(&card, LiquidViscosityBasis::Dynamic);
    let initial = curve.curve().state_at_specific_enthalpy(56_250.0).unwrap(); // 325 K
    let endpoint = original
        .prepare_caloric_endpoint(curve.curve(), initial, 68_750.0, 1e-10)
        .unwrap();
    assert_eq!(initial.temperature_k(), 325.0);
    assert_eq!(endpoint.state().temperature_k(), 350.0);
    assert_eq!(endpoint.maximum_source_density_defect_kg_m3(), 0.0);
    assert_eq!(original.reference_point(), &point());
}

#[test]
fn explicit_density_accuracy_never_overwrites_the_original_caloric_chart() {
    let card = card(LiquidViscosityBasis::Dynamic, false);
    let original = law(&card, LiquidViscosityBasis::Dynamic);
    // An explicitly authored chart has rho(350)=950; the original transport
    // source has rho(350)=900. Admission must expose this distinction.
    let curve = EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity_exact(
        card.content_hash(),
        SolidLiquidPhase::Liquid,
        0.0,
        &[
            HeatCapacityKnot {
                temperature_k: 300.0,
                specific_heat_capacity_j_kg_k: 2000.0,
                bulk_density_kg_m3: 1000.0,
            },
            HeatCapacityKnot {
                temperature_k: 400.0,
                specific_heat_capacity_j_kg_k: 4000.0,
                bulk_density_kg_m3: 900.0,
            },
        ],
    )
    .unwrap();
    let initial = curve.state_at_specific_enthalpy(0.0).unwrap();
    assert!(matches!(
        original.prepare_caloric_endpoint(&curve, initial, 125_000.0, 49.0),
        Err(LiquidCaloricError::DensityAgreement {
            observed_kg_m3: 50.0,
            maximum_allowed_kg_m3: 49.0
        }),
    ));
    let endpoint = original
        .prepare_caloric_endpoint(&curve, initial, 125_000.0, 50.0)
        .unwrap();
    assert_eq!(endpoint.state().bulk_density_kg_m3(), 950.0);
    assert_eq!(endpoint.source().density_kg_m3(), 900.0);
    assert_eq!(endpoint.maximum_source_density_defect_kg_m3(), 50.0);
    assert_eq!(endpoint.state().phase_curve_identity(), curve.identity());
}

#[test]
fn endpoint_refuses_new_context_wrong_chart_nonliquid_and_changed_source_support() {
    let card = card(LiquidViscosityBasis::Dynamic, false);
    let curve = original_curve(&card, 0.0);
    let original = law(&card, LiquidViscosityBasis::Dynamic);
    let initial = curve.curve().state_at_specific_enthalpy(0.0).unwrap();
    let changed_context = LiquidTemperatureLaw::try_new(
        &card,
        &point().with("composition", 0.2).unwrap(),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(matches!(
        changed_context.prepare_caloric_endpoint(curve.curve(), initial, 125_000.0, 1e-10),
        Err(LiquidCaloricError::InvalidBinding { .. }),
    ));
    let other = original_curve(&card, 1.0);
    assert!(matches!(
        original.prepare_caloric_endpoint(other.curve(), initial, 1.0, 1e-10),
        Err(LiquidCaloricError::Caloric(
            PhaseStateError::StateCurveMismatch
        )),
    ));
    let solid = EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity_exact(
        card.content_hash(),
        SolidLiquidPhase::Solid,
        0.0,
        &[
            HeatCapacityKnot {
                temperature_k: 300.0,
                specific_heat_capacity_j_kg_k: 2000.0,
                bulk_density_kg_m3: 1000.0,
            },
            HeatCapacityKnot {
                temperature_k: 400.0,
                specific_heat_capacity_j_kg_k: 4000.0,
                bulk_density_kg_m3: 800.0,
            },
        ],
    )
    .unwrap();
    assert!(matches!(
        original.prepare_caloric_endpoint(
            &solid,
            solid.state_at_specific_enthalpy(0.0).unwrap(),
            1.0,
            1e-10
        ),
        Err(LiquidCaloricError::InvalidBinding { .. }),
    ));
    for budget in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            original
                .prepare_caloric_endpoint(curve.curve(), initial, 1.0, budget)
                .is_err()
        );
    }
    // The same transport claim remains pinned even though another formulation
    // exists above 350 K. The caloric sources themselves remain continuous.
    let split = crate::liquid_temperature::tests::card(LiquidViscosityBasis::Dynamic, true);
    let split_law = law(&split, LiquidViscosityBasis::Dynamic);
    let split_curve = original_curve(&split, 0.0);
    let start = split_curve.curve().state_at_specific_enthalpy(0.0).unwrap();
    assert!(matches!(
        split_law.prepare_caloric_endpoint(split_curve.curve(), start, 206_250.0, 1e-10),
        Err(LiquidCaloricError::Material(_)),
    ));
    assert_eq!(split_law.reference_point(), &point());
}
