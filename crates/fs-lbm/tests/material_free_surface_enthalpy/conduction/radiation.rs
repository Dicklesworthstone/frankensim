//! G1/G3/G4: nonlinear radiation reaches native liquid enthalpy and momentum.

use super::*;
use fs_lbm::material::free_surface::enthalpy::conduction::{
    FreeSurfaceEmissivity, FreeSurfaceRadiationFace,
};
use fs_material::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_material::thermal::source::{EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY};

pub(super) fn emissive_card(scale: f64, maximum_temperature_k: f64) -> MaterialCard {
    let base = conducting_card(0.0005);
    let mut claims = base.claims().clone();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, 0.2 * scale), (500.0, 0.8 * scale)],
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with("T", 300.0, maximum_temperature_k)
                .with("pressure", 100000.0, 100000.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic original liquid emissivity".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap()
}

pub(super) fn emissive_law(card: &MaterialCard) -> LiquidTemperatureLaw<'_> {
    let base = LiquidTemperatureLaw::try_new(
        card,
        &point(400.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let mut requirements: Vec<_> = base
        .reference_state()
        .material()
        .properties()
        .iter()
        .map(|property| property.requirement().clone())
        .collect();
    requirements.push(
        ScalarPropertyRequirement::try_new(
            SURFACE_EMISSIVITY_PROPERTY,
            EMISSIVITY_DIMS,
            ScalarAdmissibility::NonNegative,
        )
        .unwrap(),
    );
    let material = resolve_material_state_point(
        card,
        &point(400.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    LiquidTemperatureLaw::try_from_material(card, material, "T", LiquidViscosityBasis::Dynamic)
        .unwrap()
}

fn radiating_face(surroundings_k: f64) -> FreeSurfaceRadiationFace {
    FreeSurfaceRadiationFace {
        liquid: EXPOSED,
        gas: [4, 3, 3],
        area_m2: 1e-6,
        emissivity: 0.8.into(),
        surroundings_k,
    }
}

#[test]
fn g1_g3_original_emissivity_data_changes_finite_heat_and_native_motion() {
    let run = |scale| {
        let card = emissive_card(scale, 500.0);
        let law = emissive_law(&card);
        let chart = caloric(&card);
        let original = law
            .reference_state()
            .material()
            .property(SURFACE_EMISSIVITY_PROPERTY)
            .unwrap()
            .answer()
            .receipt
            .selected;
        with_exposed_law(&card, true, &chart, 200000.0, &law, |thermal, gate| {
            let before = thermal.total_enthalpy_j().unwrap();
            let patch = FreeSurfaceRadiationFace {
                emissivity: FreeSurfaceEmissivity::Material,
                ..radiating_face(1500.0)
            };
            let report = thermal
                .step_with_radiative_environment(
                    &[],
                    &[],
                    controls(),
                    &mut [],
                    &[],
                    &[],
                    &[],
                    &[patch],
                    false,
                    gate,
                )
                .unwrap()
                .conduction;
            let exchange = report.radiation.as_ref().unwrap();
            let (mut temperature, mut heat) = (400.0_f64, 0.0);
            for _ in 0..4 {
                let eps = scale * (0.2 + 0.003 * (temperature - 300.0));
                let q = 0.00025
                    * eps
                    * 5.670374419e-8
                    * 1e-6
                    * (1500.0_f64.powi(4) - temperature.powi(4));
                temperature += q / 0.002;
                heat += q;
            }
            assert_eq!(report.thermal.substeps, 4);
            let endpoint = &exchange.thermal_endpoint[0];
            near(endpoint.liquid_temperature_k, temperature, 3e-13);
            near(
                endpoint.emissivity,
                scale * (0.2 + 0.003 * (temperature - 300.0)),
                3e-13,
            );
            assert_eq!(
                endpoint
                    .emissivity_source
                    .as_ref()
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
                original
            );
            near(exchange.heat_to_liquid_j, heat, 3e-12);
            near(report.flow.external_heat_j, heat, 3e-12);
            near(thermal.total_enthalpy_j().unwrap(), before + heat, 3e-12);
            let material = thermal.material_at(EXPOSED).unwrap();
            assert_eq!(material.material().card_identity(), card.content_hash());
            assert_eq!(
                material
                    .material()
                    .property(SURFACE_EMISSIVITY_PROPERTY)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
                original
            );
            assert_eq!(material.material().query().axes()["pressure"], 100000.0);
            (
                thermal.flow().lattice().grid().state_bits(),
                material.dynamic_viscosity_pa_s(),
                heat,
            )
        })
    };
    let reflecting = run(0.0);
    let dim = run(0.5);
    let bright = run(1.0);
    assert_eq!(reflecting.2, 0.0);
    assert!(bright.2 > dim.2 && bright.1 < dim.1);
    assert_ne!(bright.0, dim.0);
}

#[test]
fn g4_missing_original_or_late_unsupported_emissivity_cannot_advance_any_owner() {
    for retain_source in [false, true] {
        let card = emissive_card(1.0, 400.001);
        let chart = caloric(&card);
        let law = if retain_source {
            emissive_law(&card)
        } else {
            LiquidTemperatureLaw::try_new(
                &card,
                &point(400.0),
                "T",
                LiquidViscosityBasis::Dynamic,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap()
        };
        with_exposed_law(&card, true, &chart, 200000.0, &law, |thermal, gate| {
            let before = thermal.total_enthalpy_j().unwrap();
            let bits = thermal.flow().lattice().grid().state_bits();
            let patch = FreeSurfaceRadiationFace {
                emissivity: FreeSurfaceEmissivity::Material,
                ..radiating_face(1500.0)
            };
            assert!(
                thermal
                    .step_with_radiative_environment(
                        &[],
                        &[],
                        controls(),
                        &mut [],
                        &[],
                        &[],
                        &[],
                        &[patch],
                        false,
                        gate
                    )
                    .is_err()
            );
            assert_eq!(thermal.total_enthalpy_j().unwrap(), before);
            assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
        });
    }
}

#[test]
fn g1_radiation_and_gas_heat_share_one_interval_without_double_counting() {
    let liquid = conducting_card(1.0);
    let gas = gas_source::model(0.0, 250000.0, 100.0);
    let context = gas_source::point(450.0, 100000.0);
    with_exposed_liquid(&liquid, true, |thermal, gate| {
        let before = thermal.total_enthalpy_j().unwrap();
        let report = thermal
            .step_with_radiative_environment(
                &[],
                &[],
                controls(),
                &mut [],
                &[],
                &[],
                &[gas_contact(&gas, &context)],
                &[radiating_face(350.0)],
                false,
                gate,
            )
            .unwrap()
            .conduction;
        let gas_report = report.environment.unwrap();
        let radiation = report.radiation.unwrap();
        let kg = 0.026 * 100.0 * (0.5 + 99000.0 / 249000.0);
        let mut t: f64 = 400.0;
        let mut gas_heat = 0.0;
        let mut radiant_heat = 0.0;
        // Independent fourth powers and two-layer Fourier resistance. The
        // production radiation uses a factored law, not this subtraction.
        for _ in 0..4 {
            let q_gas = 0.00025 * 1e-6 * (450.0 - t) / (0.0005 / (t - 200.0) + 0.0005 / kg);
            let q_rad = 0.00025 * 0.8 * 5.670374419e-8 * 1e-6 * (350.0_f64.powi(4) - t.powi(4));
            t += (q_gas + q_rad) / 0.002;
            gas_heat += q_gas;
            radiant_heat += q_rad;
        }
        assert_eq!(report.thermal.substeps, 4);
        near(radiation.thermal_endpoint[0].liquid_temperature_k, t, 3e-13);
        near(radiation.heat_to_liquid_j, radiant_heat, 3e-12);
        near(gas_report.heat_to_liquid_j, gas_heat, 3e-12);
        assert!(radiation.heat_to_liquid_j < 0.0 && gas_report.heat_to_liquid_j > 0.0);
        near(
            report.thermal.external_heat_j,
            gas_heat + radiant_heat,
            3e-12,
        );
        near(report.flow.external_heat_j, gas_heat + radiant_heat, 3e-12);
        near(
            thermal.total_enthalpy_j().unwrap(),
            before + gas_heat + radiant_heat,
            3e-12,
        );
        assert!(report.flow.balance_defect_j.abs() < 1e-9);
    });
}

fn low_capacity_liquid() -> (MaterialCard, EquilibriumEnthalpyPhaseCurve) {
    let base = card(0.02, 1000.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == "specific_heat_capacity" {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, 0.01), (500.0, 0.01)],
                dims: Dims([2, 0, -2, -1, 0, 0]),
            };
        }
        claims.insert_claim(claim).unwrap();
    }
    let liquid = MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap();
    // Synthetic small capacity makes the radiation time scale comparable to
    // one native step, independently of momentum scaling. It is not a real fluid.
    let chart = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        liquid.content_hash(),
        SolidLiquidPhase::Liquid,
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 300.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 3.0,
                temperature_k: 600.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
        ],
    )
    .unwrap();
    (liquid, chart)
}

#[test]
fn g1_stiff_radiation_subcycles_and_converges_to_nonlinear_cooling_solution() {
    let (liquid, chart) = low_capacity_liquid();
    let run = |max_dt| {
        with_exposed_caloric_liquid(&liquid, false, &chart, 1.0, |thermal, gate| {
            let report = thermal
                .step_with_radiative_environment(
                    &[],
                    &[],
                    FreeSurfaceConductionControls {
                        maximum_dt_s: max_dt,
                        ..controls()
                    },
                    &mut [],
                    &[],
                    &[],
                    &[],
                    &[radiating_face(350.0)],
                    false,
                    gate,
                )
                .unwrap()
                .conduction;
            let radiation = report.radiation.unwrap();
            let t = radiation.thermal_endpoint[0].liquid_temperature_k;
            near(radiation.heat_to_liquid_j, 1e-8 * (t - 400.0), 2e-11);
            assert!(
                t > 350.0 && t < 400.0,
                "radiation must not overshoot its bath: {t}"
            );
            (t, report.thermal.substeps)
        })
    };
    assert!(
        run(0.001).1 > 1,
        "radiative stiffness must limit the shared step"
    );
    // Analytic integral of C dT/dt = eps sigma A (Ts^4 - T^4):
    // F(T/Ts)-F(T0/Ts) = eps sigma A Ts^3 t/C,
    // F(x) = log((x+1)/(x-1))/4 + atan(x)/2, for this cooling branch x>1.
    let primitive = |x: f64| 0.25 * ((x + 1.0) / (x - 1.0)).ln() + 0.5 * x.atan();
    let target =
        primitive(400.0 / 350.0) + 0.8 * 5.670374419e-8 * 1e-6 * 350.0_f64.powi(3) * 0.001 / 1e-8;
    let (mut lo, mut hi) = (350.0, 400.0);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if primitive(mid / 350.0) > target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let exact = 0.5 * (lo + hi);
    let errors = [16.0, 32.0, 64.0].map(|n| (run(0.001 / n).0 - exact).abs());
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 1.9 && pair[0] / pair[1] < 2.2,
            "{errors:?}"
        );
    }
}

#[test]
fn g3_emissivity_changes_actual_heat_viscosity_and_native_momentum() {
    let liquid = conducting_card(1.0);
    let run = |emissivity: f64| {
        with_exposed_liquid(&liquid, true, |thermal, gate| {
            let mut patch = radiating_face(1500.0);
            patch.emissivity = emissivity.into();
            let report = thermal
                .step_with_radiative_environment(
                    &[],
                    &[],
                    controls(),
                    &mut [],
                    &[],
                    &[],
                    &[],
                    &[patch],
                    false,
                    gate,
                )
                .unwrap()
                .conduction;
            (
                thermal.flow().lattice().grid().state_bits(),
                thermal
                    .material_at(EXPOSED)
                    .unwrap()
                    .dynamic_viscosity_pa_s(),
                report.radiation.unwrap().heat_to_liquid_j,
            )
        })
    };
    let insulated = run(0.0);
    let radiating = run(0.8);
    assert_eq!(insulated.2, 0.0);
    assert!(radiating.2 > 0.0 && radiating.1 < insulated.1);
    assert_ne!(radiating.0, insulated.0);
}

#[test]
fn g4_radiation_geometry_source_and_budget_refusals_are_atomic() {
    let liquid = conducting_card(1.0);
    with_exposed_liquid(&liquid, true, |thermal, gate| {
        let before = thermal.total_enthalpy_j().unwrap();
        let bits = thermal.flow().lattice().grid().state_bits();
        let patch = radiating_face(350.0);
        let mut stale = patch;
        stale.gas = [2, 3, 3];
        let mut invalid = patch;
        invalid.emissivity = 1.01.into();
        let mut hot = patch;
        hot.surroundings_k = 10000.0;
        for patches in [vec![stale], vec![invalid], vec![patch, patch], vec![hot]] {
            assert!(
                thermal
                    .step_with_radiative_environment(
                        &[],
                        &[],
                        controls(),
                        &mut [],
                        &[],
                        &[],
                        &[],
                        &patches,
                        false,
                        gate,
                    )
                    .is_err()
            );
            assert_eq!(thermal.total_enthalpy_j().unwrap(), before);
            assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
        }
        assert!(
            thermal
                .step_with_radiative_environment(
                    &[],
                    &[],
                    FreeSurfaceConductionControls {
                        maximum_substeps: 1,
                        ..controls()
                    },
                    &mut [],
                    &[],
                    &[],
                    &[],
                    &[patch],
                    false,
                    gate,
                )
                .is_err()
        );
        assert_eq!(thermal.total_enthalpy_j().unwrap(), before);
        assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(
            thermal
                .step_with_radiative_environment(
                    &[],
                    &[],
                    controls(),
                    &mut [],
                    &[],
                    &[],
                    &[],
                    &[patch],
                    false,
                    &cancelled,
                )
                .is_err()
        );
        assert_eq!(thermal.total_enthalpy_j().unwrap(), before);
        assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
    });
}
