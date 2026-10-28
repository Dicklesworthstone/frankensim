//! Synthetic source data, native J2 work, and actual lumped conduction feedback.
#[path = "thermoplastic_feedback/expansion_curve.rs"]
mod expansion_curve;
#[path = "thermoplastic_feedback/modulus_curve.rs"]
mod modulus_curve;
#[path = "thermoplastic_feedback/kinematic.rs"]
mod kinematic;

use fs_conduction::lumped::{LumpedThermalEnvironment, LumpedThermalTransport};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::plastic::source::{HardeningPropertyBasis, resolve_j2_plasticity};
use fs_material::plastic::thermal::{CaloricPoint, FeedbackConfig, FeedbackError, SourceQuery};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY,
    YOUNG_MODULUS_PROPERTY,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Frequency, Pressure, QuantitySpec, Temperature};

const CP: &str = "specific_heat_capacity";
const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);
const BASIS: HardeningPropertyBasis = HardeningPropertyBasis::PlasticStrainModulus;
const VOLUME: f64 = 1e-6;
const SURFACE: f64 = 1e-3;

fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained(),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic thermoplastic specimen".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn curve(property: &mut PropertyClaim, axis: &str, samples: Vec<(f64, f64)>) {
    let dims = property.key.quantity().dims();
    property.validity =
        ValidityDomain::unconstrained().with(axis, samples[0].0, samples[samples.len() - 1].0);
    property.value = PropertyValue::Curve {
        abscissa: axis.into(),
        abscissa_dims: Temperature::DIMS,
        knots: samples,
        dims,
    };
    property.interpolation = InterpolationPolicy::LinearInside;
}
fn properties(cp: f64, variable_cp: bool) -> Vec<PropertyClaim> {
    let mut result = vec![
        scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1.2e9),
        scalar(POISSON_RATIO_PROPERTY, Dims([0; 6]), 0.2),
        scalar(YIELD_STRESS_PROPERTY, Pressure::DIMS, 3e6),
        scalar(BASIS.property_name(), Pressure::DIMS, 1e8),
        scalar(CP, CP_DIMS, cp),
        scalar(DENSITY_PROPERTY, Density::DIMS, 1000.0),
    ];
    curve(
        &mut result[2],
        "T",
        vec![(280.0, 4e6), (300.0, 3e6), (320.0, 1e6)],
    );
    if variable_cp {
        curve(
            &mut result[4],
            "T",
            vec![(280.0, cp), (300.0, cp), (320.0, 2.0 * cp)],
        );
    }
    result
}

struct Fixture {
    card: MaterialCard,
    caloric: ResolvedHeatCapacityCurve,
    context: QueryPoint,
    axis: String,
}
impl Fixture {
    fn from_properties(
        properties: Vec<PropertyClaim>,
        reference_cp: f64,
        context: QueryPoint,
        axis: &str,
        phase: SolidLiquidPhase,
    ) -> Self {
        let mut claims = ClaimSet::new();
        for property in properties {
            claims.insert_claim(property).unwrap();
        }
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic J2 thermal solid".into(),
                phase: "fixture".into(),
                process: "numerical fixture".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        let lower = context.clone().with(axis, 280.0).unwrap();
        let upper = context.clone().with(axis, 320.0).unwrap();
        let cp_key = PropertyKey::new(CP, CP_DIMS);
        let density_key = PropertyKey::new(DENSITY_PROPERTY, Density::DIMS);
        let caloric = ResolvedHeatCapacityCurve::try_new(
            &card,
            HeatCapacitySourceQuery {
                heat_capacity: &cp_key,
                density: &density_key,
                temperature_axis: axis,
                lower: &lower,
                upper: &upper,
                selection: SelectionPolicy::SingleClaimOnly,
            },
            phase,
            -20.0 * reference_cp,
            1e-4,
        )
        .unwrap();
        Self {
            card,
            caloric,
            context,
            axis: axis.into(),
        }
    }
    fn new(cp: f64, variable_cp: bool) -> Self {
        Self::from_properties(
            properties(cp, variable_cp),
            cp,
            QueryPoint::new().with("T", 300.0).unwrap(),
            "T",
            SolidLiquidPhase::Solid,
        )
    }
    fn query(&self) -> SourceQuery<'_> {
        SourceQuery {
            card: &self.card,
            caloric: &self.caloric,
            context: self.context.clone(),
            temperature_axis: self.axis.clone(),
            selection: MaterialPropertySelection::SingleClaimOnly,
            hardening_basis: BASIS,
            thermal_expansion: None,
            strain_limit: 0.04,
        }
    }
    fn point(&self) -> CaloricPoint<'_> {
        CaloricPoint::try_new(self.query(), 0.0, VOLUME, SURFACE).unwrap()
    }
}

fn expansion_fixture(cp: f64, alpha: f64) -> Fixture {
    let mut claims = properties(cp, false);
    claims.push(scalar(
        fs_material::state_point::LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        fs_material::state_point::INVERSE_TEMPERATURE_DIMS,
        alpha,
    ));
    Fixture::from_properties(
        claims,
        cp,
        QueryPoint::new().with("T", 300.0).unwrap(),
        "T",
        SolidLiquidPhase::Solid,
    )
}

#[test]
fn g1_sourced_expansion_gives_adiabatic_elastic_cooling_and_history_consistent_unloading() {
    context(|gate, cx| {
        let alpha = 1e-5;
        let cp = 500.0;
        let fixture = expansion_fixture(cp, alpha);
        let mut query = fixture.query();
        query.thermal_expansion = Some(MaterialPropertySelection::SingleClaimOnly);
        let point = CaloricPoint::try_new(query, 0.0, VOLUME, SURFACE).unwrap();
        let ambient = environment(300.0, 0.0);
        let strain = [0.001, 0.001, 0.001, 0.0, 0.0, 0.0];
        let loaded = point
            .advance(
                strain,
                1.0,
                config(),
                |input| ambient.advance(cx, input),
                || gate.is_requested(),
            )
            .unwrap();
        // tr(sigma)=9K(eps-alpha*theta); rho Cp theta=-T0 alpha tr(sigma).
        let nine_bulk = 6e9;
        let theta = -300.0 * alpha * nine_bulk * strain[0]
            / (1000.0 * cp - 300.0 * alpha * alpha * nine_bulk);
        close(
            loaded.next.caloric_state().temperature_k(),
            300.0 + theta,
            2e-9,
        );
        assert_eq!(*loaded.next.strain(), strain);
        assert_eq!(loaded.physical_heat_j, 0.0);
        assert!(loaded.reversible_heat_j < 0.0);
        close(loaded.reversible_heat_j, 1000.0 * VOLUME * cp * theta, 1e-9);
        close(loaded.combined_residual_j, 0.0, 1e-8);
        assert!(loaded.thermoelastic_storage_change_j > 0.0);
        assert_eq!(loaded.next.resolved().material().properties().len(), 7);

        let old_t = loaded.next.caloric_state().temperature_k();
        let unloading_theta = old_t * alpha * nine_bulk * strain[0]
            / (1000.0 * cp - old_t * alpha * alpha * nine_bulk);
        let unloaded = loaded
            .next
            .advance(
                [0.0; 6],
                1.0,
                config(),
                |input| ambient.advance(cx, input),
                || gate.is_requested(),
            )
            .unwrap();
        close(
            unloaded.next.caloric_state().temperature_k(),
            old_t + unloading_theta,
            2e-9,
        );
        assert!(unloaded.reversible_heat_j > 0.0);
        assert_eq!(unloaded.next.history().alpha, 0.0);
        close(unloaded.combined_residual_j, 0.0, 1e-8);
        assert_eq!(point.time_s(), 0.0);
    });
}

#[test]
fn g4_expansion_requires_a_supported_source_and_positive_isochoric_capacity() {
    for alpha in [1e-3, -1e-3] {
        let fixture = expansion_fixture(500.0, alpha);
        let mut query = fixture.query();
        query.thermal_expansion = Some(MaterialPropertySelection::SingleClaimOnly);
        assert!(CaloricPoint::try_new(query, 0.0, VOLUME, SURFACE).is_err());
    }
    let mut claims = properties(500.0, false);
    let mut expansion = scalar(
        fs_material::state_point::LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        fs_material::state_point::INVERSE_TEMPERATURE_DIMS,
        1e-5,
    );
    curve(&mut expansion, "T", vec![(280.0, 1e-5), (320.0, 2e-5)]);
    claims.push(expansion);
    let fixture = Fixture::from_properties(
        claims,
        500.0,
        QueryPoint::new().with("T", 300.0).unwrap(),
        "T",
        SolidLiquidPhase::Solid,
    );
    let mut query = fixture.query();
    query.thermal_expansion = Some(MaterialPropertySelection::SingleClaimOnly);
    // A variable coefficient is now admitted with its exact primitive and
    // reversible caloric terms; the source/capacity refusals remain.
    assert!(CaloricPoint::try_new(query, 0.0, VOLUME, SURFACE).is_ok());
    let missing = Fixture::new(500.0, false);
    let mut query = missing.query();
    query.thermal_expansion = Some(MaterialPropertySelection::SingleClaimOnly);
    assert!(CaloricPoint::try_new(query, 0.0, VOLUME, SURFACE).is_err());
}

fn config() -> FeedbackConfig {
    FeedbackConfig {
        max_iterations: 100,
        temperature_tolerance_k: 1e-8,
        energy_tolerance_j: 1e-8,
        yield_tolerance_pa: 1e-3,
        relaxation: 1.0,
    }
}
fn environment(temperature_k: f64, htc: f64) -> LumpedThermalEnvironment {
    LumpedThermalEnvironment {
        temperature_k,
        radiation_temperature_k: temperature_k,
        convection_w_per_m2_k: htc,
        // Explicit independent boundary data; no measured transport provenance
        // is invented for this synthetic source-to-caloric fixture.
        transport: LumpedThermalTransport::try_declared(10.0, 0.0).unwrap(),
        enthalpy_tolerance_j_kg: 1e-9,
        maximum_thermal_residual_j: 1e-10,
    }
}
fn context<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0,
                kernel_id: 0x4a325448,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&gate, &cx)
    })
}
fn shear(xy: f64) -> [f64; 6] {
    [0.0, 0.0, 0.0, xy, 0.0, 0.0]
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}
fn insulated_oracle(cp: f64) -> f64 {
    // Pure shear: q_trial=sqrt(3)*10 MPa, delta_alpha=(q_trial-sigma_y)/1.6 GPa.
    // At T>=300 K, sigma_y=3 MPa-(0.1 MPa/K)*theta. Combining this with
    // rho*Cp*theta=sigma_y*delta_alpha gives this independent quadratic root.
    let q_trial = 3.0_f64.sqrt() * 1e7;
    let slope = 1e5;
    let yield_zero = 3e6;
    let b = 1.6e9 * 1000.0 * cp + slope * (q_trial - 2.0 * yield_zero);
    let k = yield_zero * (q_trial - yield_zero);
    300.0 + 2.0 * k / (b + (b * b + 4.0 * slope * slope * k).sqrt())
}

#[test]
fn material_cp_swap_changes_native_yielding_and_matches_independent_thermal_root() {
    context(|gate, cx| {
        let low = Fixture::new(5.0, false);
        let high = Fixture::new(10.0, false);
        let ambient = environment(300.0, 0.0);
        let low_point = low.point();
        let high_point = high.point();
        let a = low_point
            .advance(
                shear(0.01),
                1.0,
                config(),
                |input| ambient.advance(cx, input),
                || gate.is_requested(),
            )
            .unwrap();
        let b = high_point
            .advance(
                shear(0.01),
                1.0,
                config(),
                |input| ambient.advance(cx, input),
                || gate.is_requested(),
            )
            .unwrap();
        close(
            a.next.caloric_state().temperature_k(),
            insulated_oracle(5.0),
            2e-7,
        );
        close(
            b.next.caloric_state().temperature_k(),
            insulated_oracle(10.0),
            2e-7,
        );
        assert!(
            a.next.caloric_state().temperature_k() > b.next.caloric_state().temperature_k() + 1.0
        );
        assert!(a.next.history().alpha > b.next.history().alpha);
        assert!(a.next.stress_pa()[3] < b.next.stress_pa()[3]);
        for step in [&a, &b] {
            let temperature = step.next.caloric_state().temperature_k();
            let yield_pa = 3e6 - 1e5 * (temperature - 300.0);
            let alpha = (3.0_f64.sqrt() * 1e7 - yield_pa) / 1.6e9;
            close(step.next.history().alpha, alpha, 1e-11);
            close(
                step.next.stress_pa()[3],
                (yield_pa + 1e8 * alpha) / 3.0_f64.sqrt(),
                0.1,
            );
            close(step.physical_heat_j, VOLUME * yield_pa * alpha, 1e-10);
            assert!(step.mechanical.energy.final_state.hardening_j_m3 > 0.0);
            assert!(step.mechanical.energy.numerical_dissipation_j_m3 > 0.0);
            assert_eq!(step.transport.external_heat_j, 0.0);
            assert!(step.combined_residual_j.abs() <= config().energy_tolerance_j);
            close(step.next.mass_kg(), 0.001, 1e-15);
            let bundle = step.next.resolved().material();
            assert_eq!(bundle.properties().len(), 6);
            for property in bundle.properties() {
                let t = property
                    .answer()
                    .receipt
                    .query_point
                    .iter()
                    .find(|(name, _)| name == "T")
                    .unwrap()
                    .1;
                assert_eq!(t, temperature);
            }
        }
        assert_eq!(low_point.history().alpha, 0.0);
        assert_eq!(low_point.caloric_state().temperature_k(), 300.0);
        assert_eq!(low_point.time_s(), 0.0);
    });
}

#[test]
fn cp_curve_boundary_cooling_and_fixed_strain_heating_reach_the_same_history() {
    context(|gate, cx| {
        let constant = Fixture::new(5.0, false);
        let varying = Fixture::new(5.0, true);
        let insulated = environment(300.0, 0.0);
        let cooled = environment(300.0, 10.0);
        let original = constant.point();
        let load = original
            .advance(
                shear(0.01),
                1.0,
                config(),
                |input| insulated.advance(cx, input),
                || gate.is_requested(),
            )
            .unwrap();
        let curve_load = varying
            .point()
            .advance(
                shear(0.01),
                1.0,
                config(),
                |input| insulated.advance(cx, input),
                || gate.is_requested(),
            )
            .unwrap();
        let cold_load = original
            .advance(
                shear(0.01),
                1.0,
                config(),
                |input| cooled.advance(cx, input),
                || gate.is_requested(),
            )
            .unwrap();
        assert!(
            curve_load.next.caloric_state().temperature_k()
                < load.next.caloric_state().temperature_k()
        );
        assert!(curve_load.next.history().alpha < load.next.history().alpha);
        assert!(
            cold_load.next.caloric_state().temperature_k()
                < load.next.caloric_state().temperature_k()
        );
        assert!(cold_load.transport.external_heat_j < 0.0);
        let hot = environment(310.0, 5.0);
        let held = load
            .next
            .advance(
                shear(0.01),
                1.0,
                config(),
                |input| hot.advance(cx, input),
                || gate.is_requested(),
            )
            .unwrap();
        assert_eq!(held.mechanical_work_j, 0.0);
        assert!(held.physical_heat_j > 0.0);
        assert!(held.transport.external_heat_j > 0.0);
        assert!(held.next.history().alpha > load.next.history().alpha);
        assert!(held.mechanical.energy.stored_change_j_m3() < 0.0);
        assert!(held.combined_residual_j.abs() <= config().energy_tolerance_j);
        close(
            held.mechanical.energy.initial.total_j_m3(),
            load.mechanical.energy.final_state.total_j_m3(),
            1e-6,
        );
        let cool = environment(300.0, 5.0);
        let unload = load
            .next
            .advance(
                shear(0.009),
                1.0,
                config(),
                |input| cool.advance(cx, input),
                || gate.is_requested(),
            )
            .unwrap();
        assert_eq!(unload.next.history(), load.next.history());
        assert_eq!(unload.physical_heat_j, 0.0);
        assert!(unload.mechanical_work_j < 0.0);
        assert!(unload.transport.external_heat_j < 0.0);
        assert!(unload.combined_residual_j.abs() <= config().energy_tolerance_j);
        assert_eq!(unload.next.time_s(), 2.0);
    });
}

#[test]
fn removed_or_duplicated_heat_and_unbounded_report_allowances_refuse_atomically() {
    context(|gate, cx| {
        let fixture = Fixture::new(5.0, false);
        let foreign = Fixture::new(10.0, false);
        let point = fixture.point();
        let original = (
            *point.strain(),
            point.history().clone(),
            *point.caloric_state(),
            point.time_s(),
            point.resolved().clone(),
        );
        let ambient = environment(300.0, 0.0);
        for heat_factor in [0.0, 2.0] {
            let result = point.advance(
                shear(0.01),
                1.0,
                config(),
                |input| {
                    ambient.advance(cx, input).map(|mut report| {
                        // Keep the real owner's independent boundary transfer;
                        // corrupt only internal heat credit on the SAME chart.
                        report.state = input
                            .curve
                            .state_at_specific_enthalpy(
                                input.initial.specific_enthalpy_j_kg()
                                    + heat_factor * input.internal_heat_j / input.mass_kg,
                            )
                            .unwrap();
                        report
                    })
                },
                || gate.is_requested(),
            );
            assert!(matches!(result, Err(FeedbackError::Invalid(_))));
        }
        let oversized = point.advance(
            shear(0.01),
            1.0,
            config(),
            |input| {
                ambient.advance(cx, input).map(|mut report| {
                    report.energy_residual_tolerance_j = f64::MAX;
                    report
                })
            },
            || gate.is_requested(),
        );
        assert!(matches!(oversized, Err(FeedbackError::Invalid(_))));
        let foreign_state = point.advance(
            shear(0.01),
            1.0,
            config(),
            |input| {
                ambient.advance(cx, input).map(|mut report| {
                    report.state = foreign
                        .caloric
                        .curve()
                        .state_at_specific_enthalpy(0.0)
                        .unwrap();
                    report
                })
            },
            || gate.is_requested(),
        );
        assert!(matches!(foreign_state, Err(FeedbackError::Invalid(_))));
        let limited = point.advance(
            shear(0.01),
            1.0,
            FeedbackConfig {
                max_iterations: 1,
                ..config()
            },
            |input| ambient.advance(cx, input),
            || gate.is_requested(),
        );
        assert!(matches!(limited, Err(FeedbackError::NotConverged { .. })));
        let bad_biot = environment(300.0, 1e8);
        assert!(matches!(
            point.advance(
                shear(0.01),
                1.0,
                config(),
                |input| bad_biot.advance(cx, input),
                || gate.is_requested()
            ),
            Err(FeedbackError::Transport(_))
        ));
        let cancelled = point.advance(
            shear(0.01),
            1.0,
            config(),
            |input| {
                let result = ambient.advance(cx, input);
                gate.request();
                result
            },
            || gate.is_requested(),
        );
        assert!(matches!(cancelled, Err(FeedbackError::Cancelled)));
        assert_eq!(
            (
                *point.strain(),
                point.history().clone(),
                *point.caloric_state(),
                point.time_s(),
                point.resolved().clone()
            ),
            original
        );
    });
}

#[test]
fn thermal_evolution_keeps_initial_claims_and_refuses_disjoint_source_handoffs() {
    context(|gate, cx| {
        let mut data = properties(5.0, false);
        let mut competitor = scalar(YIELD_STRESS_PROPERTY, Pressure::DIMS, 5e5);
        competitor.validity = ValidityDomain::unconstrained().with("T", 303.0, 320.0);
        data.push(competitor.clone());
        let fixture = Fixture::from_properties(
            data,
            5.0,
            QueryPoint::new().with("T", 300.0).unwrap(),
            "T",
            SolidLiquidPhase::Solid,
        );
        let point = fixture.point();
        let original_claim = point
            .resolved()
            .material()
            .property(YIELD_STRESS_PROPERTY)
            .unwrap()
            .answer()
            .receipt
            .selected;
        let ambient = environment(300.0, 0.0);
        let step = point
            .advance(
                shear(0.01),
                1.0,
                config(),
                |input| ambient.advance(cx, input),
                || gate.is_requested(),
            )
            .unwrap();
        assert!(step.next.caloric_state().temperature_k() > 303.0);
        assert_eq!(
            step.next
                .resolved()
                .material()
                .property(YIELD_STRESS_PROPERTY)
                .unwrap()
                .answer()
                .receipt
                .selected,
            original_claim
        );
        let final_query = fixture
            .context
            .clone()
            .with("T", step.next.caloric_state().temperature_k())
            .unwrap();
        assert!(
            resolve_j2_plasticity(
                &fixture.card,
                &final_query,
                MaterialPropertySelection::SingleClaimOnly,
                BASIS,
                0.04
            )
            .is_err()
        );
        let mut short = properties(5.0, false);
        curve(
            &mut short[2],
            "T",
            vec![(280.0, 4e6), (300.0, 3e6), (302.0, 2.8e6)],
        );
        short.push(competitor);
        let unsupported = Fixture::from_properties(
            short,
            5.0,
            QueryPoint::new().with("T", 300.0).unwrap(),
            "T",
            SolidLiquidPhase::Solid,
        );
        assert!(matches!(
            CaloricPoint::try_new(unsupported.query(), 0.0, VOLUME, SURFACE),
            Err(FeedbackError::Support { .. })
        ));
    });
}

#[test]
fn nonaffine_moduli_variable_poisson_density_phase_and_foreign_context_refuse() {
    for (index, value) in [(0, 1.2e9), (1, 0.2), (3, 1e8), (5, 1000.0)] {
        let mut data = properties(5.0, false);
        // Initial and both span endpoints agree. Non-affine modulus slopes need
        // entropy jumps; fixed Poisson/density also cannot vary in the interior.
        curve(
            &mut data[index],
            "T",
            vec![
                (280.0, value),
                (300.0, value),
                (310.0, 1.1 * value),
                (320.0, value),
            ],
        );
        let fixture = Fixture::from_properties(
            data,
            5.0,
            QueryPoint::new().with("T", 300.0).unwrap(),
            "T",
            SolidLiquidPhase::Solid,
        );
        assert!(CaloricPoint::try_new(fixture.query(), 0.0, VOLUME, SURFACE).is_err());
    }
    let liquid = Fixture::from_properties(
        properties(5.0, false),
        5.0,
        QueryPoint::new().with("T", 300.0).unwrap(),
        "T",
        SolidLiquidPhase::Liquid,
    );
    assert!(CaloricPoint::try_new(liquid.query(), 0.0, VOLUME, SURFACE).is_err());
    let a = Fixture::new(5.0, false);
    let b = Fixture::new(10.0, false);
    let mut wrong_card = a.query();
    wrong_card.caloric = &b.caloric;
    assert!(CaloricPoint::try_new(wrong_card, 0.0, VOLUME, SURFACE).is_err());
    let pressure = Fixture::from_properties(
        properties(5.0, false),
        5.0,
        QueryPoint::new()
            .with("T", 300.0)
            .unwrap()
            .with("pressure", 1e5)
            .unwrap(),
        "T",
        SolidLiquidPhase::Solid,
    );
    let mut wrong_context = pressure.query();
    wrong_context.context = wrong_context.context.with("pressure", 2e5).unwrap();
    assert!(CaloricPoint::try_new(wrong_context, 0.0, VOLUME, SURFACE).is_err());
}

#[test]
fn no_source_role_may_freeze_an_evolving_plastic_strain_or_rate_coordinate() {
    for index in 0..6 {
        let mut data = properties(5.0, false);
        data[index].validity = data[index]
            .validity
            .clone()
            .with("plastic_strain", 0.0, 0.1);
        let fixture = Fixture::from_properties(
            data,
            5.0,
            QueryPoint::new()
                .with("T", 300.0)
                .unwrap()
                .with("plastic_strain", 0.0)
                .unwrap(),
            "T",
            SolidLiquidPhase::Solid,
        );
        assert!(matches!(
            CaloricPoint::try_new(fixture.query(), 0.0, VOLUME, SURFACE),
            Err(FeedbackError::Invalid(_))
        ));
    }
    let mut data = properties(5.0, false);
    data[4].validity = ValidityDomain::unconstrained().with("drive", 0.0, 10.0);
    let rate_context = QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with_quantity("drive", QuantitySpec::dimensional(Frequency::DIMS), 1.0)
        .unwrap();
    let fixture = Fixture::from_properties(data, 5.0, rate_context, "T", SolidLiquidPhase::Solid);
    assert!(matches!(
        CaloricPoint::try_new(fixture.query(), 0.0, VOLUME, SURFACE),
        Err(FeedbackError::Invalid(_))
    ));
}

#[test]
fn renamed_temperature_requires_semantics_and_retains_the_typed_source_context() {
    let axis = "specimen_kelvin";
    for quantity in [None, Some(QuantitySpec::dimensional(Temperature::DIMS))] {
        let mut data = properties(5.0, false);
        // All scalar claims make this exercise the explicit coordinate gate;
        // a temperature-curve classifier is not doing the refusal for it.
        data[2] = scalar(YIELD_STRESS_PROPERTY, Pressure::DIMS, 3e6);
        let point_context = match quantity {
            None => QueryPoint::new().with(axis, 300.0).unwrap(),
            Some(q) => QueryPoint::new().with_quantity(axis, q, 300.0).unwrap(),
        };
        let fixture =
            Fixture::from_properties(data, 5.0, point_context, axis, SolidLiquidPhase::Solid);
        assert!(matches!(
            CaloricPoint::try_new(fixture.query(), 0.0, VOLUME, SURFACE),
            Err(FeedbackError::Invalid(_))
        ));
    }
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let point_context = QueryPoint::new()
        .with_quantity(axis, absolute, 300.0)
        .unwrap();
    let mut data = properties(5.0, true);
    for property in &mut data {
        if let PropertyValue::Curve { abscissa, .. } = &mut property.value {
            *abscissa = axis.into();
            property.validity = ValidityDomain::unconstrained().with(axis, 280.0, 320.0);
        }
    }
    let fixture = Fixture::from_properties(data, 5.0, point_context, axis, SolidLiquidPhase::Solid);
    context(|gate, cx| {
        let ambient = environment(300.0, 0.0);
        let step = fixture
            .point()
            .advance(
                shear(0.01),
                1.0,
                config(),
                |input| ambient.advance(cx, input),
                || gate.is_requested(),
            )
            .unwrap();
        assert!(step.next.caloric_state().temperature_k() > 300.0);
        assert_eq!(
            step.next.resolved().material().axis_quantities().get(axis),
            Some(&absolute)
        );
        assert!(step.combined_residual_j.abs() <= config().energy_tolerance_j);
    });
}

fn pressure_surface(property: &mut PropertyClaim, values: [f64; 4]) {
    property.value = PropertyValue::Surface2D {
        axes: ["T".into(), "pressure".into()],
        axis_dims: [Temperature::DIMS, Pressure::DIMS],
        samples: [[280.0, 1e5], [320.0, 1e5], [280.0, 3e5], [320.0, 3e5]]
            .into_iter()
            .zip(values)
            .collect(),
        triangles: vec![[0, 1, 2], [1, 3, 2]],
        dims: property.key.quantity().dims(),
    };
    property.validity = ValidityDomain::unconstrained()
        .with("T", 280.0, 320.0)
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            1e5,
            3e5,
        );
    property.interpolation = InterpolationPolicy::LinearInside;
}

fn pressure_properties() -> Vec<PropertyClaim> {
    let mut data = properties(5.0, false);
    for (property, values) in data.iter_mut().zip([
        [1e9, 1e9, 1.4e9, 1.4e9],
        [0.125, 0.125, 0.375, 0.375],
        [5e6, 1e6, 5e6, 1e6],
        [0.5e8, 0.5e8, 1.5e8, 1.5e8],
        [5.0, 5.0, 15.0, 15.0],
        [800.0, 800.0, 1200.0, 1200.0],
    ]) {
        pressure_surface(property, values);
    }
    data
}

fn pressure_context(pressure_pa: f64) -> QueryPoint {
    QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            pressure_pa,
        )
        .unwrap()
}

#[test]
fn source_surfaces_admit_energetic_constants_over_the_complete_fixed_pressure_line() {
    let mut data = pressure_properties();
    let mut expansion = scalar(
        fs_material::state_point::LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        fs_material::state_point::INVERSE_TEMPERATURE_DIMS,
        2e-5,
    );
    pressure_surface(&mut expansion, [1e-5, 1e-5, 3e-5, 3e-5]);
    data.push(expansion);
    let fixture = Fixture::from_properties(
        data,
        10.0,
        pressure_context(2e5),
        "T",
        SolidLiquidPhase::Solid,
    );
    let mut query = fixture.query();
    query.thermal_expansion = Some(MaterialPropertySelection::SingleClaimOnly);
    let point = CaloricPoint::try_new(query, 0.0, VOLUME, SURFACE).unwrap();
    close(point.resolved().law().elastic.youngs, 1.2e9, 1e-6);
    close(point.resolved().law().elastic.poisson, 0.25, 1e-15);
    close(point.resolved().hardening_modulus_pa(), 1e8, 1e-7);
    close(point.mass_kg(), 0.001, 1e-15);
    assert_eq!(point.resolved().material().properties().len(), 7);
    for property in point.resolved().material().properties() {
        fixture
            .card
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
}

#[test]
fn every_surface_cell_crossing_checks_modulus_affinity_fixed_coefficients_and_positive_yield() {
    for (index, value) in [(0, 1.2e9), (1, 0.25), (3, 1e8), (5, 1000.0)] {
        let mut data = pressure_properties();
        // At this pressure the two endpoints agree, but the diagonal cut has
        // a non-affine modulus or varying fixed coefficient. Endpoint sampling is unsound.
        pressure_surface(&mut data[index], [value, 0.8 * value, 0.8 * value, value]);
        let fixture = Fixture::from_properties(
            data,
            10.0,
            pressure_context(2e5),
            "T",
            SolidLiquidPhase::Solid,
        );
        assert!(CaloricPoint::try_new(fixture.query(), 0.0, VOLUME, SURFACE).is_err());
    }
    let mut data = pressure_properties();
    pressure_surface(&mut data[2], [6e6, -1e6, -1e6, 6e6]);
    let fixture = Fixture::from_properties(
        data,
        7.5,
        pressure_context(1.5e5),
        "T",
        SolidLiquidPhase::Solid,
    );
    // Initial T=300 K and both thermal boundaries have positive strength;
    // the negative interior cut at T=310 K still forbids this source path.
    assert!(
        resolve_j2_plasticity(
            &fixture.card,
            &fixture.context,
            MaterialPropertySelection::SingleClaimOnly,
            BASIS,
            0.04
        )
        .is_ok()
    );
    assert!(CaloricPoint::try_new(fixture.query(), 0.0, VOLUME, SURFACE).is_err());
}

#[test]
fn surface_source_holes_and_second_temperature_bounds_cannot_be_frozen_away() {
    let mut data = pressure_properties();
    let right = f64::from_bits(310_f64.to_bits() + 1);
    data[2].value = PropertyValue::Surface2D {
        axes: ["T".into(), "pressure".into()],
        axis_dims: [Temperature::DIMS, Pressure::DIMS],
        samples: [
            [280.0, 1e5],
            [310.0, 1e5],
            [280.0, 3e5],
            [310.0, 3e5],
            [right, 1e5],
            [320.0, 1e5],
            [right, 3e5],
            [320.0, 3e5],
        ]
        .map(|at| (at, 3e6))
        .to_vec(),
        triangles: vec![[0, 1, 2], [1, 3, 2], [4, 5, 6], [5, 7, 6]],
        dims: Pressure::DIMS,
    };
    let gap = Fixture::from_properties(
        data,
        10.0,
        pressure_context(2e5),
        "T",
        SolidLiquidPhase::Solid,
    );
    assert!(
        resolve_j2_plasticity(
            &gap.card,
            &gap.context,
            MaterialPropertySelection::SingleClaimOnly,
            BASIS,
            0.04
        )
        .is_ok()
    );
    assert!(CaloricPoint::try_new(gap.query(), 0.0, VOLUME, SURFACE).is_err());

    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    for index in 0..6 {
        let mut data = properties(5.0, false);
        data[index].validity =
            data[index]
                .validity
                .clone()
                .with_quantity("other_temperature", absolute, 290.0, 310.0);
        let query = QueryPoint::new()
            .with("T", 300.0)
            .unwrap()
            .with_quantity("other_temperature", absolute, 300.0)
            .unwrap();
        let fixture = Fixture::from_properties(data, 5.0, query, "T", SolidLiquidPhase::Solid);
        assert!(CaloricPoint::try_new(fixture.query(), 0.0, VOLUME, SURFACE).is_err());
    }
}
