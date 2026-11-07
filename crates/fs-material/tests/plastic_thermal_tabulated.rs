//! Synthetic Cauchy-table, native thermal-history and reciprocal-energy regressions.
//! These manufactured coefficients are not an experimental alloy calibration.
#![deny(unsafe_code)]

use core::convert::Infallible;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::SmallStrainLaw;
use fs_material::phase::{
    HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase, UniformEnthalpyStep,
    UniformEnthalpyStepInput,
};
use fs_material::plastic::J2State;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::rate::{
    OVERSTRESS_EXPONENT_PROPERTY, OVERSTRESS_STRESS_SCALE_PROPERTY, OVERSTRESS_VISCOSITY_DIMS,
    OVERSTRESS_VISCOSITY_PROPERTY,
};
use fs_material::plastic::reduction::ReductionConfig;
use fs_material::plastic::section::thermal::{ThermoplasticSectionError, ThermoplasticSectionStep};
use fs_material::plastic::section::{PlaneStressLayer, SectionBudget, ThicknessQuadrature};
use fs_material::plastic::source::{HardeningPropertyBasis, KINEMATIC_HARDENING_MODULUS_PROPERTY};
use fs_material::plastic::tabulated::source::{
    ACCUMULATED_PLASTIC_STRAIN_AXIS as ALPHA, CAUCHY_YIELD_CURVE_PROPERTY,
    KIRCHHOFF_YIELD_CURVE_PROPERTY, YieldStressMeasure,
};
use fs_material::plastic::tabulated::thermal::AffineTabulatedHardening;
use fs_material::plastic::tabulated::{HardeningPoint, TabulatedHardening};
use fs_material::plastic::thermal::family::ThermoplasticSourceFamily;
use fs_material::plastic::thermal::tabulated::ThermalTabulatedSource;
use fs_material::plastic::thermal::{
    CaloricPoint, FeedbackConfig, FeedbackError, SourceQuery, ThermoplasticSource,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY,
    YOUNG_MODULUS_PROPERTY,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};

const E: f64 = 1.2e9;
const NU: f64 = 0.2;
const HK: f64 = 1.6e8;
const DE: f64 = -5e6;
const DHK: f64 = -1e6;
const CP: f64 = 500.0;
const RHO: f64 = 1000.0;
const CP_NAME: &str = "specific_heat_capacity";
const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);
const TEMPERATURE: &str = "specimen-temperature";
const BASIS: HardeningPropertyBasis = HardeningPropertyBasis::MixedPlasticStrainModuli;
const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const KNOTS: [(f64, f64); 5] = [
    (0.0, 3e6),
    (0.002, 4e6),
    (0.006, 4e6),
    (0.012, 7e6),
    (0.04, 8.4e6),
];
const SLOPES: [f64; 5] = [-10e3, -20e3, -20e3, -50e3, -80e3];
const SAMPLE_LIMIT: usize = 15;

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(TEMPERATURE, temperature_quantity(), t)
        .unwrap()
}
fn property(name: &str, dims: Dims, value: f64, slope: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: if slope == 0.0 {
            PropertyValue::Scalar { value, dims }
        } else {
            PropertyValue::Curve {
                abscissa: TEMPERATURE.into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, value - 20.0 * slope), (320.0, value + 20.0 * slope)],
                dims,
            }
        },
        validity: ValidityDomain::unconstrained().with_quantity(
            TEMPERATURE,
            temperature_quantity(),
            280.0,
            320.0,
        ),
        interpolation: if slope == 0.0 {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        },
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic thermal Cauchy-table experiment".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn surface() -> PropertyClaim {
    let mut claim = property(CAUCHY_YIELD_CURVE_PROPERTY, Pressure::DIMS, KNOTS[0].1, 0.0);
    claim.value = PropertyValue::Surface2D {
        axes: [ALPHA.into(), TEMPERATURE.into()],
        axis_dims: [Dims::NONE, Temperature::DIMS],
        samples: [280.0, 300.0, 320.0]
            .into_iter()
            .flat_map(|temperature| {
                KNOTS
                    .into_iter()
                    .zip(SLOPES)
                    .map(move |((alpha, y), slope)| {
                        ([alpha, temperature], slope.mul_add(temperature - 300.0, y))
                    })
            })
            .collect(),
        triangles: vec![],
        dims: Pressure::DIMS,
    };
    claim.validity =
        claim
            .validity
            .with_quantity(ALPHA, QuantitySpec::dimensional(Dims::NONE), 0.0, 0.04);
    claim.interpolation = InterpolationPolicy::TabulatedOnly;
    claim
}
fn properties(expansion: f64, variable_moduli: bool) -> Vec<PropertyClaim> {
    vec![
        property(
            YOUNG_MODULUS_PROPERTY,
            Pressure::DIMS,
            E,
            if variable_moduli { DE } else { 0.0 },
        ),
        property(POISSON_RATIO_PROPERTY, Dims::NONE, NU, 0.0),
        property(YIELD_STRESS_PROPERTY, Pressure::DIMS, KNOTS[0].1, SLOPES[0]),
        property(BASIS.property_name(), Pressure::DIMS, 0.0, 0.0),
        property(
            KINEMATIC_HARDENING_MODULUS_PROPERTY,
            Pressure::DIMS,
            HK,
            if variable_moduli { DHK } else { 0.0 },
        ),
        property(CP_NAME, CP_DIMS, CP, 0.0),
        property(DENSITY_PROPERTY, Density::DIMS, RHO, 0.0),
        property(
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            INVERSE_TEMPERATURE_DIMS,
            expansion,
            0.0,
        ),
        property(
            OVERSTRESS_VISCOSITY_PROPERTY,
            OVERSTRESS_VISCOSITY_DIMS,
            1e7,
            -50e3,
        ),
        property(OVERSTRESS_STRESS_SCALE_PROPERTY, Pressure::DIMS, 2e6, 0.0),
        property(OVERSTRESS_EXPONENT_PROPERTY, Dims::NONE, 2.0, 0.0),
        surface(),
    ]
}

struct Fixture {
    card: MaterialCard,
    caloric: ResolvedHeatCapacityCurve,
}
impl Fixture {
    fn new(properties: Vec<PropertyClaim>) -> Self {
        let mut claims = ClaimSet::new();
        for property in properties {
            claims.insert_claim(property).unwrap();
        }
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic thermal hardening specimen".into(),
                phase: "solid".into(),
                process: "numerical regression".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        let caloric = ResolvedHeatCapacityCurve::try_new(
            &card,
            HeatCapacitySourceQuery {
                heat_capacity: &PropertyKey::new(CP_NAME, CP_DIMS),
                density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
                temperature_axis: TEMPERATURE,
                lower: &point(280.0),
                upper: &point(320.0),
                selection: SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Solid,
            -20.0 * CP,
            1e-6,
        )
        .unwrap();
        Self { card, caloric }
    }
    fn query(&self) -> SourceQuery<'_> {
        SourceQuery {
            card: &self.card,
            caloric: &self.caloric,
            context: point(300.0),
            temperature_axis: TEMPERATURE.into(),
            selection: MaterialPropertySelection::SingleClaimOnly,
            hardening_basis: BASIS,
            thermal_expansion: Some(MaterialPropertySelection::SingleClaimOnly),
            strain_limit: 0.04,
        }
    }
    fn base(&self) -> ThermoplasticSource<'_> {
        ThermoplasticSource::try_new(self.query(), 0.0).unwrap()
    }
    fn source(&self) -> ThermalTabulatedSource<'_> {
        self.base()
            .with_tabulated_hardening(SINGLE, SAMPLE_LIMIT)
            .unwrap()
    }
}

fn h(t: f64) -> f64 {
    CP * (t - 300.0)
}
fn tolerance() -> EnergyTolerances {
    EnergyTolerances {
        energy_j_m3: 1e-5,
        yield_pa: 1e-3,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; allowance {tolerance:e}"
    );
}
fn inner(a: &[f64; 6], b: &[f64; 6]) -> f64 {
    (0..6)
        .map(|i| a[i] * b[i] * if i < 3 { 1.0 } else { 2.0 })
        .sum()
}
fn ordinates(t: f64) -> [f64; 5] {
    core::array::from_fn(|i| KNOTS[i].1 + SLOPES[i] * (t - 300.0))
}
fn primitive(alpha: f64, values: [f64; 5]) -> f64 {
    let mut integral = 0.0;
    for i in 0..KNOTS.len() - 1 {
        if alpha <= KNOTS[i].0 {
            break;
        }
        let dx = alpha.min(KNOTS[i + 1].0) - KNOTS[i].0;
        integral += dx
            * (values[i] - values[0]
                + 0.5 * dx * (values[i + 1] - values[i]) / (KNOTS[i + 1].0 - KNOTS[i].0));
    }
    integral
}
fn fields(t: f64, total: &[f64; 6], state: &J2State, expansion: f64) -> (f64, f64, [f64; 6]) {
    let youngs = E + DE * (t - 300.0);
    let mu = youngs / (2.0 * (1.0 + NU));
    let lambda = youngs * NU / ((1.0 + NU) * (1.0 - 2.0 * NU));
    let e: [f64; 6] = core::array::from_fn(|i| {
        total[i] - state.plastic_strain[i] - if i < 3 { expansion * (t - 300.0) } else { 0.0 }
    });
    let trace = e[0] + e[1] + e[2];
    let elastic = mu * inner(&e, &e) + 0.5 * lambda * trace * trace;
    let shape = inner(&state.plastic_strain, &state.plastic_strain) / 3.0;
    let stress: [f64; 6] =
        core::array::from_fn(|i| 2.0 * mu * e[i] + if i < 3 { lambda * trace } else { 0.0 });
    let energy = elastic + (HK + DHK * (t - 300.0)) * shape + primitive(state.alpha, ordinates(t));
    let entropy = expansion * (stress[0] + stress[1] + stress[2])
        - DE / youngs * elastic
        - DHK * shape
        - primitive(state.alpha, SLOPES);
    (energy, entropy, stress)
}

#[test]
fn heated_loaded_history_crosses_original_knot_and_retains_unloaded_shape_with_exact_entropy() {
    let expansion = 2e-5;
    let fixture = Fixture::new(properties(expansion, true));
    let source = fixture.source();
    let evidence = source.hardening_source().source();
    assert_eq!(evidence.card_identity(), fixture.card.content_hash());
    assert_eq!(evidence.stress_measure(), YieldStressMeasure::Cauchy);
    assert_eq!(evidence.claim(), &surface());
    assert_eq!(evidence.context(), &point(300.0));
    assert_eq!(evidence.observations().len(), SAMPLE_LIMIT);
    for (_, answer) in evidence.observations() {
        fixture
            .card
            .claims()
            .verify_receipt(&answer.receipt)
            .unwrap();
    }
    // A non-grid endpoint also checks consistent source-line versus affine-table
    // rounding, with the original receipts preserved rather than rewritten.
    source.resolve_enthalpy(h(310.3)).unwrap();
    let s = 0.00928;
    let strain = [s, -0.5 * s, -0.5 * s, 0.0, 0.0, 0.0];
    let initial = J2State {
        plastic_strain: [0.0; 6],
        alpha: 0.0,
    };
    let loaded = source
        .advance_with_energy(0.0, 0.0, &[0.0; 6], &initial, &strain, tolerance())
        .unwrap();
    let old = loaded.mechanical.state.clone();
    let c0 = 3.0 * E / (2.0 * (1.0 + NU));
    near(old.alpha, (c0 * s - 4e6) / (c0 + HK), 2e-14);
    assert!(old.alpha > 0.002 && old.alpha < 0.006);

    let (_, hot_law) = source.resolve_enthalpy(h(320.0)).unwrap();
    let hot_mechanical =
        core::array::from_fn(|i| strain[i] - if i < 3 { expansion * 20.0 } else { 0.0 });
    assert!(
        hot_law
            .tabulated_law()
            .advance_with_energy(&strain, &old, &hot_mechanical, tolerance())
            .is_err(),
        "public isothermal admission must still refuse an old state outside its own yield surface"
    );
    let heated = source
        .advance_with_energy(0.0, h(320.0), &strain, &old, &strain, tolerance())
        .unwrap();
    let next = &heated.mechanical.state;
    let c1 = 3.0 * (E + DE * 20.0) / (2.0 * (1.0 + NU));
    let hk1 = HK + DHK * 20.0;
    near(
        next.alpha,
        (c1 * s - 3.6e6 + 4e8 * 0.006) / (c1 + hk1 + 4e8),
        2e-14,
    );
    assert!(next.alpha > 0.006 && next.alpha < 0.012);
    let (psi0, entropy0, _) = fields(300.0, &strain, &old, expansion);
    let (psi1, entropy1, stress1) = fields(320.0, &strain, next, expansion);
    for (actual, expected) in heated.mechanical.stress_pa.iter().zip(stress1) {
        near(*actual, expected, 2e-6);
    }
    let fixed_old_mechanical =
        core::array::from_fn(|i| strain[i] + if i < 3 { expansion * 20.0 } else { 0.0 });
    let coefficient = fields(320.0, &fixed_old_mechanical, &old, expansion).0 - psi0;
    let entropy_storage = 320.0 * entropy1 - 300.0 * entropy0;
    near(heated.coefficient_storage_change_j_m3, coefficient, 1e-7);
    near(heated.stored_change_j_m3(), psi1 - psi0, 1e-7);
    near(
        heated.thermoelastic_storage_change_j_m3,
        entropy_storage,
        1e-7,
    );
    near(
        heated.reversible_heat_j_m3,
        (stress1[0] + stress1[1] + stress1[2]) * expansion * 20.0 - coefficient - entropy_storage,
        1e-7,
    );
    near(
        heated.mechanical.energy.physical_dissipation_j_m3,
        ordinates(320.0)[0] * (next.alpha - old.alpha),
        1e-7,
    );
    near(heated.mechanical_work_j_m3, 0.0, 0.0);
    near(
        heated.stored_change_j_m3()
            + heated.thermoelastic_storage_change_j_m3
            + heated.caloric_heat_j_m3()
            + heated.mechanical.energy.numerical_dissipation_j_m3,
        0.0,
        1e-5,
    );
    assert_eq!(old, loaded.mechanical.state);

    let residual_strain = next.plastic_strain;
    let unloaded = source
        .advance_with_energy(h(320.0), 0.0, &strain, next, &residual_strain, tolerance())
        .unwrap();
    near(unloaded.mechanical.state.alpha, next.alpha, 2e-14);
    for stress in unloaded.mechanical.stress_pa {
        near(stress, 0.0, 2e-6);
    }
    assert!(unloaded.mechanical.state.plastic_strain[0] > 0.006);
    near(
        unloaded.mechanical.energy.physical_dissipation_j_m3,
        0.0,
        1e-7,
    );
    source
        .validate_accepted(
            0.0,
            &residual_strain,
            &unloaded.mechanical.state,
            tolerance().yield_pa,
        )
        .unwrap();
}

#[test]
fn actual_duration_and_refreshed_kinetics_drive_viscous_holds_without_resetting_table_history() {
    let fixture = Fixture::new(properties(0.0, false));
    let strain = [0.012, -0.006, -0.006, 0.0, 0.0, 0.0];
    for power in [false, true] {
        let base = if power {
            fixture
                .base()
                .with_power_law_overstress(MaterialPropertySelection::SingleClaimOnly)
                .unwrap()
        } else {
            fixture
                .base()
                .with_linear_overstress(MaterialPropertySelection::SingleClaimOnly)
                .unwrap()
        };
        let source = base.with_tabulated_hardening(SINGLE, SAMPLE_LIMIT).unwrap();
        assert!(source.resolve_enthalpy(h(315.0)).is_err());
        assert!(source.resolve_enthalpy_for_step(h(315.0), 0.0).is_err());
        let initial = J2State {
            plastic_strain: [0.0; 6],
            alpha: 0.0,
        };
        let mut accumulated = Vec::new();
        for duration in [0.001, 0.1] {
            let (_, law) = source
                .resolve_enthalpy_for_step(h(315.0), duration)
                .unwrap();
            let rate = law.rate_law().unwrap();
            near(rate.viscosity_pa_s(), 9.25e6, 1e-7);
            assert_eq!(rate.duration_s(), duration);
            assert_eq!(rate.stress_scale_pa(), if power { Some(2e6) } else { None });
            let increment = source
                .advance_with_energy_for_step(
                    0.0,
                    h(315.0),
                    &[0.0; 6],
                    &initial,
                    &strain,
                    tolerance(),
                    duration,
                )
                .unwrap();
            let native = law.response(&strain, &initial).unwrap();
            assert_eq!(increment.mechanical.state, native.state);
            assert_eq!(increment.mechanical.stress_pa, native.stress_pa);
            let delta = native.state.alpha;
            let overstress = rate.endpoint_overstress_pa(delta).unwrap();
            near(
                law.tabulated_law()
                    .yield_function_with_state(&native.stress_pa, &native.state),
                overstress,
                1e-3,
            );
            assert!(overstress > 0.0);
            near(
                increment.mechanical.energy.physical_dissipation_j_m3,
                (ordinates(315.0)[0] + overstress) * delta,
                1e-6,
            );
            source
                .validate_accepted(h(315.0), &strain, &native.state, 1e-3)
                .unwrap();
            let held = source
                .advance_with_energy_for_step(
                    h(315.0),
                    h(315.0),
                    &strain,
                    &native.state,
                    &strain,
                    tolerance(),
                    duration,
                )
                .unwrap();
            assert!(held.mechanical.state.alpha > native.state.alpha);
            assert!(held.mechanical.stress_pa[0] < native.stress_pa[0]);
            assert!(held.mechanical.energy.physical_dissipation_j_m3 > 0.0);
            accumulated.push(delta);
        }
        assert!(accumulated[1] > accumulated[0]);
        assert_eq!(initial.alpha, 0.0);
    }
}

#[test]
fn existing_layered_plane_stress_owner_preserves_thermal_table_histories_and_heat_partition() {
    let fixture = Fixture::new(properties(0.0, false));
    let source = fixture.source();
    let layers = vec![
        PlaneStressLayer {
            bottom_z_m: -0.005,
            top_z_m: 0.0,
            subdivisions: 1,
            law: source.clone(),
        },
        PlaneStressLayer {
            bottom_z_m: 0.0,
            top_z_m: 0.005,
            subdivisions: 1,
            law: source,
        },
    ];
    let initial_h = [0.0; 2];
    let hot_h = [h(310.0), h(320.0)];
    let reduction = ReductionConfig {
        max_iterations: 64,
        stress_tolerance_pa: 1e-4,
        strain_limit: 0.04,
    };
    let budget = SectionBudget {
        maximum_layers: 2,
        maximum_quadrature_points: 6,
    };
    let loading = ThermoplasticSectionStep::new(
        &layers,
        &initial_h,
        &initial_h,
        0.1,
        ThicknessQuadrature::Gauss3,
        reduction,
        budget,
    )
    .unwrap();
    let virgin = loading.section().initial_state();
    let strain = [0.009, 0.0, 0.0006, 0.08, 0.0, 0.01];
    let loaded = loading.section().evaluate(strain, &virgin).unwrap();
    loading
        .audit_energy(&virgin, &loaded.state, tolerance(), 1e-6)
        .unwrap();
    let saved = loaded.state.clone();
    let heating = ThermoplasticSectionStep::new(
        &layers,
        &initial_h,
        &hot_h,
        0.2,
        ThicknessQuadrature::Gauss3,
        reduction,
        budget,
    )
    .unwrap();
    let heated = heating.section().evaluate(strain, &loaded.state).unwrap();
    let audit = heating
        .audit_energy(&loaded.state, &heated.state, tolerance(), 1e-6)
        .unwrap();
    assert_eq!(audit.total.resultants, heated.resultants);
    assert!(audit.total.physical_heat_j_m2 > 0.0);
    assert!(audit.total.coefficient_storage_change_j_m2 < 0.0);
    assert!(audit.total.thermoelastic_storage_change_j_m2.abs() > 1e-6);
    near(audit.total.generalized_work_j_m2, 0.0, 0.0);
    for (layer, energy) in audit.layers.iter().enumerate() {
        let yield_stress = heating.section().layers()[layer]
            .law
            .tabulated_law()
            .hardening()
            .initial_yield_stress_pa();
        let physical: f64 = heating
            .section()
            .points()
            .iter()
            .enumerate()
            .filter(|(_, point)| point.layer == layer)
            .map(|(index, point)| {
                point.weight_m
                    * yield_stress
                    * (heated.state.points[index].history.alpha
                        - loaded.state.points[index].history.alpha)
            })
            .sum();
        near(energy.physical_heat_j_m2, physical, 1e-8);
        near(energy.balance_residual_j_m2, 0.0, 1e-6);
    }
    assert_eq!(loaded.state, saved);
    assert_eq!(
        heated,
        heating.section().evaluate(strain, &loaded.state).unwrap()
    );
    let mut corrupt = heated.state.clone();
    corrupt.points[0].history.alpha += 1e-5;
    assert!(matches!(
        heating.audit_energy(&loaded.state, &corrupt, tolerance(), 1e-6),
        Err(ThermoplasticSectionError::State { .. })
    ));
}

fn feedback_config() -> FeedbackConfig {
    FeedbackConfig {
        max_iterations: 64,
        temperature_tolerance_k: 1e-9,
        energy_tolerance_j: 1e-8,
        yield_tolerance_pa: 1e-3,
        relaxation: 1.0,
    }
}
fn caloric_proposal(
    input: UniformEnthalpyStepInput<'_>,
    external_heat_j: f64,
) -> Result<UniformEnthalpyStep<()>, Infallible> {
    let h = input.initial.specific_enthalpy_j_kg()
        + (input.internal_heat_j + external_heat_j) / input.mass_kg;
    Ok(UniformEnthalpyStep {
        state: input.curve.state_at_specific_enthalpy(h).unwrap(),
        external_heat_j,
        energy_residual_tolerance_j: 1e-8,
        report: (),
    })
}

#[test]
fn original_uniform_feedback_keeps_loaded_heating_candidates_private_until_all_energy_gates_pass() {
    let fixture = Fixture::new(properties(0.0, false));
    let point = CaloricPoint::try_new(fixture.query(), 0.0, 1e-6, 1e-3)
        .unwrap()
        .with_tabulated_hardening(SINGLE, SAMPLE_LIMIT)
        .unwrap();
    let s = 0.00928;
    let strain = [s, -0.5 * s, -0.5 * s, 0.0, 0.0, 0.0];
    let loaded = point
        .advance(
            strain,
            1.0,
            feedback_config(),
            |input| caloric_proposal(input, 0.0),
            || false,
        )
        .unwrap();
    let saved = loaded.next.history().clone();
    let old_caloric = *loaded.next.caloric_state();
    assert!(saved.alpha > 0.002 && saved.alpha < 0.006);
    let heat = loaded.next.mass_kg() * CP * 10.0;
    let mut limited = feedback_config();
    limited.max_iterations = 1;
    assert!(matches!(
        loaded.next.advance(
            strain,
            1.0,
            limited,
            |input| caloric_proposal(input, heat),
            || false
        ),
        Err(FeedbackError::NotConverged { .. })
    ));
    assert_eq!(*loaded.next.history(), saved);
    assert_eq!(*loaded.next.caloric_state(), old_caloric);
    let heated = loaded
        .next
        .advance(
            strain,
            1.0,
            feedback_config(),
            |input| caloric_proposal(input, heat),
            || false,
        )
        .unwrap();
    assert!(heated.next.history().alpha > 0.006);
    assert!(heated.physical_heat_j > 0.0);
    assert!(heated.reversible_heat_j < 0.0);
    assert!(heated.iterations > 1);
    near(heated.thermal_residual_j, 0.0, 1e-8);
    near(heated.combined_residual_j, 0.0, 1e-8);
    assert!(std::ptr::eq(heated.next.caloric_source(), &fixture.caloric));
    assert_eq!(
        heated
            .next
            .tabulated_hardening_source()
            .source()
            .card_identity(),
        fixture.card.content_hash()
    );
    assert_eq!(*loaded.next.history(), saved);
}

#[test]
fn admission_refuses_competing_strength_source_and_stress_relabeling() {
    let fixture = Fixture::new(properties(0.0, false));
    assert!(
        fixture
            .base()
            .with_tabulated_hardening(SINGLE, SAMPLE_LIMIT - 1)
            .is_err()
    );
    assert!(fixture.source().resolve_enthalpy(h(320.001)).is_err());
    let mut nonzero_h = properties(0.0, false);
    nonzero_h[3] = property(BASIS.property_name(), Pressure::DIMS, 1e6, 0.0);
    assert!(
        Fixture::new(nonzero_h)
            .base()
            .with_tabulated_hardening(SINGLE, SAMPLE_LIMIT)
            .is_err()
    );
    let mut disagreement = properties(0.0, false);
    disagreement[2].value = PropertyValue::Curve {
        abscissa: TEMPERATURE.into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![
            (280.0, 3.2e6),
            (300.0, 3e6),
            (310.0, 2.95e6),
            (320.0, 2.8e6),
        ],
        dims: Pressure::DIMS,
    };
    assert!(
        Fixture::new(disagreement)
            .base()
            .with_tabulated_hardening(SINGLE, SAMPLE_LIMIT)
            .is_err(),
        "matching reference and endpoints cannot conceal an interior initial-yield discrepancy"
    );
    let mut kirchhoff = properties(0.0, false);
    kirchhoff.last_mut().unwrap().key =
        PropertyKey::new(KIRCHHOFF_YIELD_CURVE_PROPERTY, Pressure::DIMS);
    assert!(
        Fixture::new(kirchhoff)
            .base()
            .with_tabulated_hardening(SINGLE, SAMPLE_LIMIT)
            .is_err()
    );
}

#[test]
fn shared_affine_storage_retains_signed_interior_extrema_and_original_finite_api() {
    let reference = TabulatedHardening::new(vec![
        HardeningPoint {
            plastic_strain: 0.0,
            yield_stress_pa: 10.0,
        },
        HardeningPoint {
            plastic_strain: 1.0,
            yield_stress_pa: 11.0,
        },
        HardeningPoint {
            plastic_strain: 2.0,
            yield_stress_pa: 12.0,
        },
    ])
    .unwrap();
    let affine = AffineTabulatedHardening::new(reference.clone(), vec![0.0, 1.0, -1.0]).unwrap();
    near(affine.maximum_absolute_storage_slope_j_m3_k(), 0.75, 0.0);
    near(
        affine
            .hardening_energy_temperature_derivative_j_m3_k(1.5)
            .unwrap(),
        0.75,
        0.0,
    );
    near(
        affine
            .hardening_energy_temperature_derivative_j_m3_k(2.0)
            .unwrap(),
        0.5,
        0.0,
    );
    // Independent numeric fixture for the unchanged finite wrapper interface.
    // No Cauchy material/source claim is passed into a Kirchhoff binding.
    let finite = fs_material::plastic::finite::thermal::FiniteJ2ThermalTabulatedHardening::new(
        reference,
        vec![0.0, 1.0, -1.0],
    )
    .unwrap();
    near(
        finite
            .hardening_energy_temperature_derivative_j_m3_k(1.5)
            .unwrap(),
        0.75,
        0.0,
    );
    assert!(
        affine
            .hardening_energy_temperature_derivative_j_m3_k(2.001)
            .is_err()
    );
    assert!(
        finite
            .hardening_energy_temperature_derivative_j_m3_k(2.001)
            .is_err()
    );
}
