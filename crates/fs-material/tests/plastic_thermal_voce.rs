//! Synthetic material/source/energy regressions, not experimental calibration.

use core::convert::Infallible;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::SmallStrainLaw;
use fs_material::phase::{
    HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase, UniformEnthalpyStep,
    UniformEnthalpyStepInput,
};
use fs_material::plastic::J2State;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::rate::{OVERSTRESS_VISCOSITY_DIMS, OVERSTRESS_VISCOSITY_PROPERTY};
use fs_material::plastic::source::{HardeningPropertyBasis, KINEMATIC_HARDENING_MODULUS_PROPERTY};
use fs_material::plastic::thermal::{
    CaloricPoint, FeedbackConfig, FeedbackError, SourceQuery, ThermoplasticSource,
};
use fs_material::plastic::voce::source::{
    VOCE_SATURATION_INCREMENT_PROPERTY as Q, VOCE_SATURATION_RATE_PROPERTY as B,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY,
    YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Density, Dims, Pressure, Temperature};

const E: f64 = 1.2e9;
const NU: f64 = 0.2;
const H: f64 = 1e8;
const HK: f64 = 1.6e8;
const SY: f64 = 3e6;
const RATE: f64 = 300.0;
const Q0: f64 = 6e6;
const DQ: f64 = -2e5;
const CP: f64 = 500.0;
const RHO: f64 = 1000.0;
const BASIS: HardeningPropertyBasis = HardeningPropertyBasis::MixedPlasticStrainModuli;
const CP_NAME: &str = "specific_heat_capacity";
const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);

fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained(),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic thermal Voce source".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

fn curve(property: &mut PropertyClaim, samples: Vec<(f64, f64)>) {
    property.validity =
        ValidityDomain::unconstrained().with("T", samples[0].0, samples[samples.len() - 1].0);
    property.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: samples,
        dims: property.key.quantity().dims(),
    };
    property.interpolation = InterpolationPolicy::LinearInside;
}

fn properties(q: f64, slope: f64) -> Vec<PropertyClaim> {
    let mut data = vec![
        scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
        scalar(POISSON_RATIO_PROPERTY, Dims([0; 6]), NU),
        scalar(YIELD_STRESS_PROPERTY, Pressure::DIMS, SY),
        scalar(BASIS.property_name(), Pressure::DIMS, H),
        scalar(KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, HK),
        scalar(CP_NAME, CP_DIMS, CP),
        scalar(DENSITY_PROPERTY, Density::DIMS, RHO),
        scalar(Q, Pressure::DIMS, q),
        scalar(B, Dims([0; 6]), RATE),
        scalar(
            OVERSTRESS_VISCOSITY_PROPERTY,
            OVERSTRESS_VISCOSITY_DIMS,
            1e6,
        ),
    ];
    if slope != 0.0 {
        curve(
            &mut data[7],
            [280.0, 290.0, 300.0, 315.0, 320.0]
                .map(|t| (t, q + slope * (t - 300.0)))
                .to_vec(),
        );
    }
    data
}

struct Fixture {
    card: MaterialCard,
    caloric: ResolvedHeatCapacityCurve,
}

impl Fixture {
    fn new(data: Vec<PropertyClaim>) -> Self {
        let mut claims = ClaimSet::new();
        for property in data {
            claims.insert_claim(property).unwrap();
        }
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic saturating thermal solid".into(),
                phase: "solid".into(),
                process: "numerical fixture".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        let lower = QueryPoint::new().with("T", 280.0).unwrap();
        let upper = QueryPoint::new().with("T", 320.0).unwrap();
        let caloric = ResolvedHeatCapacityCurve::try_new(
            &card,
            HeatCapacitySourceQuery {
                heat_capacity: &PropertyKey::new(CP_NAME, CP_DIMS),
                density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
                temperature_axis: "T",
                lower: &lower,
                upper: &upper,
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
            context: QueryPoint::new().with("T", 300.0).unwrap(),
            temperature_axis: "T".into(),
            selection: MaterialPropertySelection::SingleClaimOnly,
            hardening_basis: BASIS,
            thermal_expansion: None,
            strain_limit: 0.04,
        }
    }

    fn source(&self) -> ThermoplasticSource<'_> {
        ThermoplasticSource::try_new(self.query(), 0.0).unwrap()
    }

    fn voce(&self) -> ThermoplasticSource<'_> {
        self.source()
            .with_voce_hardening(MaterialPropertySelection::SingleClaimOnly)
            .unwrap()
    }
}

fn tolerances() -> EnergyTolerances {
    EnergyTolerances {
        energy_j_m3: 1e-5,
        yield_pa: 1e-3,
    }
}
fn h(temperature: f64) -> f64 {
    CP * (temperature - 300.0)
}
fn inner(a: &[f64; 6], b: &[f64; 6]) -> f64 {
    (0..6)
        .map(|i| a[i] * b[i] * if i < 3 { 1.0 } else { 2.0 })
        .sum()
}
fn shape(alpha: f64) -> f64 {
    alpha + (-RATE * alpha).exp_m1() / RATE
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance={tolerance:e}"
    );
}
fn energy(strain: &[f64; 6], state: &J2State, saturation: f64) -> f64 {
    let e: [f64; 6] = core::array::from_fn(|i| strain[i] - state.plastic_strain[i]);
    let mu = E / (2.0 * (1.0 + NU));
    let lambda = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU));
    let trace = e[0] + e[1] + e[2];
    mu * inner(&e, &e)
        + 0.5 * lambda * trace * trace
        + 0.5 * H * state.alpha * state.alpha
        + HK * inner(&state.plastic_strain, &state.plastic_strain) / 3.0
        + saturation * shape(state.alpha)
}

#[test]
fn affine_saturation_preserves_sources_and_native_endpoint_stress_tangent_and_history() {
    let fixture = Fixture::new(properties(Q0, DQ));
    let source = fixture.voce();
    assert!(source.has_voce_hardening());
    let strain = [0.012, -0.006, -0.006, 0.002, 0.0, 0.0];
    let (_, cold) = source.resolve_expanded_enthalpy(0.0).unwrap();
    let old = cold.update_state(&strain, &cold.initial_state());
    assert!(old.alpha > 0.0);
    let (_, hot) = source.resolve_expanded_enthalpy(h(315.0)).unwrap();
    let native = hot.voce_law().unwrap();
    near(native.saturation_pa(), Q0 + DQ * 15.0, 0.0);
    assert_eq!(hot.resolved().material().properties().len(), 9);
    for property in hot.resolved().material().properties() {
        fixture
            .card
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
        assert_eq!(
            property.answer().receipt.selected,
            source
                .reference()
                .material()
                .property(property.requirement().name())
                .unwrap()
                .answer()
                .receipt
                .selected
        );
    }
    let stress = hot.stress(&strain, &old);
    let next = hot.update_state(&strain, &old);
    assert!(next.alpha > old.alpha);
    assert_eq!(stress, native.stress(&strain, &old));
    assert_eq!(next, native.update_state(&strain, &old));
    let tangent = hot.tangent(&strain, &old);
    assert_eq!(tangent, native.tangent(&strain, &old));
    for j in 0..6 {
        let mut a = strain;
        let mut b = strain;
        a[j] += 1e-8;
        b[j] -= 1e-8;
        let sa = hot.stress(&a, &old);
        let sb = hot.stress(&b, &old);
        for i in 0..6 {
            near(tangent[i][j], (sa[i] - sb[i]) / 2e-8, 2.0);
        }
    }
    source
        .validate_committed_history(hot.resolved(), &strain, &next, 1e-3)
        .unwrap();
    let (_, timed) = source
        .resolve_expanded_enthalpy_for_step(h(315.0), 7.0)
        .unwrap();
    assert_eq!(timed.stress(&strain, &old), stress);
    assert_eq!(timed.update_state(&strain, &old), next);
    assert_eq!(timed.tangent(&strain, &old), tangent);
}

#[test]
fn heating_a_loaded_point_yields_and_closes_independent_storage_entropy_and_work() {
    let fixture = Fixture::new(properties(Q0, DQ));
    let source = fixture.voce();
    let strain = [0.012, -0.006, -0.006, 0.002, 0.0, 0.0];
    let (_, cold) = source.resolve_expanded_enthalpy(0.0).unwrap();
    let old = cold.update_state(&strain, &cold.initial_state());
    let increment = source
        .advance_with_energy(0.0, h(315.0), &strain, &old, &strain, tolerances())
        .unwrap();
    let next = &increment.mechanical.state;
    let delta = next.alpha - old.alpha;
    assert!(delta > 0.0);
    let effective: [f64; 6] = core::array::from_fn(|i| {
        increment.mechanical.stress_pa[i] - (2.0 / 3.0) * HK * next.plastic_strain[i]
    });
    let radius = SY + H * next.alpha + (Q0 + 15.0 * DQ) * -(-RATE * next.alpha).exp_m1();
    near((1.5 * inner(&effective, &effective)).sqrt(), radius, 1e-3);
    let coefficient = 15.0 * DQ * shape(old.alpha);
    let entropy_storage = -DQ * (315.0 * shape(next.alpha) - 300.0 * shape(old.alpha));
    near(increment.coefficient_storage_change_j_m3, coefficient, 1e-7);
    near(
        increment.thermoelastic_storage_change_j_m3,
        entropy_storage,
        1e-7,
    );
    near(
        increment.stored_change_j_m3(),
        energy(&strain, next, Q0 + 15.0 * DQ) - energy(&strain, &old, Q0),
        1e-7,
    );
    near(
        increment.mechanical.energy.physical_dissipation_j_m3,
        SY * delta,
        1e-7,
    );
    near(
        increment.reversible_heat_j_m3,
        -coefficient - entropy_storage,
        1e-7,
    );
    near(increment.mechanical_work_j_m3, 0.0, 0.0);
    near(
        increment.stored_change_j_m3()
            + entropy_storage
            + increment.caloric_heat_j_m3()
            + increment.mechanical.energy.numerical_dissipation_j_m3,
        0.0,
        1e-5,
    );
    let timed = source
        .advance_with_energy_for_step(0.0, h(315.0), &strain, &old, &strain, tolerances(), 0.2)
        .unwrap();
    assert_eq!(timed.mechanical, increment.mechanical);
    assert_eq!(timed.reversible_heat_j_m3, increment.reversible_heat_j_m3);
    // A return at the accepted endpoint must not consume plastic history again.
    let repeated = source
        .advance_with_energy(h(315.0), h(315.0), &strain, next, &strain, tolerances())
        .unwrap();
    near(repeated.mechanical.state.alpha, next.alpha, 1e-16);
}

#[test]
fn frozen_history_heating_and_cooling_account_for_reversible_saturation_energy() {
    let fixture = Fixture::new(properties(Q0, DQ));
    let source = fixture.voce();
    let p = 0.001;
    let state = J2State {
        plastic_strain: [p, -p / 2.0, -p / 2.0, 0.0, 0.0, 0.0],
        alpha: 0.002,
    };
    let strain = [
        p + 0.0001,
        -p / 2.0 + 0.0001,
        -p / 2.0 + 0.0001,
        0.0,
        0.0,
        0.0,
    ];
    for (t0, t1) in [(290.0, 315.0), (315.0, 290.0)] {
        let next = source
            .advance_with_energy(h(t0), h(t1), &strain, &state, &strain, tolerances())
            .unwrap();
        assert_eq!(next.mechanical.state, state);
        let coefficient = DQ * (t1 - t0) * shape(state.alpha);
        near(next.coefficient_storage_change_j_m3, coefficient, 1e-7);
        near(next.stored_change_j_m3(), coefficient, 1e-7);
        near(next.thermoelastic_storage_change_j_m3, -coefficient, 1e-7);
        near(next.reversible_heat_j_m3, 0.0, 1e-7);
        near(next.mechanical.energy.physical_dissipation_j_m3, 0.0, 0.0);
        near(next.mechanical.energy.numerical_dissipation_j_m3, 0.0, 0.0);
    }
}

#[test]
fn saturation_combines_with_original_affine_moduli_and_integrated_expansion() {
    let (de, dh, dhk, expansion) = (-1e7, -1e6, -2e6, 1e-5);
    let mut data = properties(Q0, DQ);
    for (index, base, slope) in [(0, E, de), (3, H, dh), (4, HK, dhk)] {
        curve(
            &mut data[index],
            [280.0, 300.0, 320.0]
                .map(|t| (t, base + slope * (t - 300.0)))
                .to_vec(),
        );
    }
    data.push(scalar(
        LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        INVERSE_TEMPERATURE_DIMS,
        expansion,
    ));
    let fixture = Fixture::new(data);
    let mut query = fixture.query();
    query.thermal_expansion = Some(MaterialPropertySelection::SingleClaimOnly);
    let source = ThermoplasticSource::try_new(query, 0.0)
        .unwrap()
        .with_voce_hardening(MaterialPropertySelection::SingleClaimOnly)
        .unwrap();
    let p = 0.001;
    let state = J2State {
        plastic_strain: [p, -p / 2.0, -p / 2.0, 0.0, 0.0, 0.0],
        alpha: 0.002,
    };
    let elastic0 = [0.0004, 0.0003, 0.0002, 0.0001, 0.0, 0.0];
    let elastic1 = [0.0005, 0.0004, 0.0003, -0.00008, 0.0, 0.0];
    let (t0, t1) = (290.0, 315.0);
    let total = |t: f64, e: &[f64; 6]| {
        core::array::from_fn(|i| {
            e[i] + state.plastic_strain[i] + if i < 3 { expansion * (t - 300.0) } else { 0.0 }
        })
    };
    let old = total(t0, &elastic0);
    let new = total(t1, &elastic1);
    let account = source
        .advance_with_energy(h(t0), h(t1), &old, &state, &new, tolerances())
        .unwrap();
    assert_eq!(account.mechanical.state, state);
    let fields = |t: f64, e: &[f64; 6]| {
        let youngs = E + de * (t - 300.0);
        let mu = youngs / (2.0 * (1.0 + NU));
        let lambda = youngs * NU / ((1.0 + NU) * (1.0 - 2.0 * NU));
        let trace = e[0] + e[1] + e[2];
        let psi_elastic = mu * inner(e, e) + 0.5 * lambda * trace * trace;
        let iso_shape = 0.5 * state.alpha * state.alpha;
        let kin_shape = inner(&state.plastic_strain, &state.plastic_strain) / 3.0;
        let saturation_shape = shape(state.alpha);
        let trace_stress = youngs * trace / (1.0 - 2.0 * NU);
        let psi = psi_elastic
            + (H + dh * (t - 300.0)) * iso_shape
            + (HK + dhk * (t - 300.0)) * kin_shape
            + (Q0 + DQ * (t - 300.0)) * saturation_shape;
        let entropy = expansion * trace_stress
            - de / youngs * psi_elastic
            - dh * iso_shape
            - dhk * kin_shape
            - DQ * saturation_shape;
        (psi, entropy, trace_stress)
    };
    let (psi0, s0, _) = fields(t0, &elastic0);
    let (psi1, s1, trace_stress1) = fields(t1, &elastic1);
    let (psi_old_new_t, _, _) = fields(t1, &elastic0);
    let coefficient = psi_old_new_t - psi0;
    let entropy_storage = t1 * s1 - t0 * s0;
    near(account.stored_change_j_m3(), psi1 - psi0, 1e-7);
    near(account.coefficient_storage_change_j_m3, coefficient, 1e-7);
    near(
        account.thermoelastic_storage_change_j_m3,
        entropy_storage,
        1e-7,
    );
    near(
        account.reversible_heat_j_m3,
        trace_stress1 * expansion * (t1 - t0) - coefficient - entropy_storage,
        1e-7,
    );
    near(
        account.mechanical_work_j_m3
            - account.stored_change_j_m3()
            - entropy_storage
            - account.caloric_heat_j_m3()
            - account.mechanical.energy.numerical_dissipation_j_m3,
        0.0,
        1e-5,
    );
}

#[test]
fn a_zero_saturation_endpoint_keeps_nonzero_entropy_and_full_retained_history() {
    let fixture = Fixture::new(properties(Q0, -Q0 / 20.0));
    let source = fixture.voce();
    let p = 0.001;
    let state = J2State {
        plastic_strain: [p, -p / 2.0, -p / 2.0, 0.0, 0.0, 0.0],
        alpha: 0.002,
    };
    let strain = state.plastic_strain;
    let next = source
        .advance_with_energy(0.0, h(320.0), &strain, &state, &strain, tolerances())
        .unwrap();
    assert_eq!(next.mechanical.state, state);
    near(
        next.coefficient_storage_change_j_m3,
        -Q0 * shape(state.alpha),
        1e-7,
    );
    assert!(next.thermoelastic_storage_change_j_m3 > 0.0);
    near(
        next.thermoelastic_storage_change_j_m3,
        Q0 * shape(state.alpha),
        1e-7,
    );
    near(next.caloric_heat_j_m3(), 0.0, 1e-7);
}

#[test]
fn explicit_zero_saturation_is_bitwise_original_thermal_j2() {
    let fixture = Fixture::new(properties(0.0, 0.0));
    let original = fixture.source();
    let selected = fixture.voce();
    let strain = [0.012, -0.006, -0.006, 0.002, 0.0, 0.0];
    let old = original.reference().initial_state();
    let a = original
        .advance_with_energy(0.0, h(310.0), &[0.0; 6], &old, &strain, tolerances())
        .unwrap();
    let b = selected
        .advance_with_energy(0.0, h(310.0), &[0.0; 6], &old, &strain, tolerances())
        .unwrap();
    assert_eq!(a.mechanical, b.mechanical);
    assert_eq!(
        a.coefficient_storage_change_j_m3.to_bits(),
        b.coefficient_storage_change_j_m3.to_bits()
    );
    assert_eq!(
        a.thermoelastic_storage_change_j_m3.to_bits(),
        b.thermoelastic_storage_change_j_m3.to_bits()
    );
    assert_eq!(
        a.reversible_heat_j_m3.to_bits(),
        b.reversible_heat_j_m3.to_bits()
    );
    let (_, a) = original.resolve_expanded_enthalpy(h(310.0)).unwrap();
    let (_, b) = selected.resolve_expanded_enthalpy(h(310.0)).unwrap();
    assert_eq!(a.stress(&strain, &old), b.stress(&strain, &old));
    assert_eq!(a.tangent(&strain, &old), b.tangent(&strain, &old));
    assert_eq!(a.update_state(&strain, &old), b.update_state(&strain, &old));
}

#[test]
fn source_admission_checks_internal_breakpoints_and_refuses_unsupported_combinations() {
    let mut cases = Vec::new();
    let mut kinked = properties(Q0, 0.0);
    curve(
        &mut kinked[7],
        vec![(280.0, Q0), (290.0, Q0 / 2.0), (300.0, Q0), (320.0, Q0)],
    );
    cases.push(kinked);
    let mut changing_b = properties(Q0, 0.0);
    curve(
        &mut changing_b[8],
        vec![
            (280.0, RATE),
            (290.0, RATE + 1.0),
            (300.0, RATE),
            (320.0, RATE),
        ],
    );
    cases.push(changing_b);
    let mut negative = properties(Q0, 0.0);
    curve(
        &mut negative[7],
        vec![(280.0, -Q0), (300.0, Q0), (320.0, Q0)],
    );
    cases.push(negative);
    let mut short = properties(Q0, 0.0);
    curve(&mut short[7], vec![(290.0, Q0), (310.0, Q0)]);
    cases.push(short);
    let mut wrong_units = properties(Q0, 0.0);
    wrong_units[7] = scalar(Q, Density::DIMS, Q0);
    cases.push(wrong_units);
    let mut missing = properties(Q0, 0.0);
    missing.remove(8);
    cases.push(missing);
    for data in cases {
        let fixture = Fixture::new(data);
        assert!(
            fixture
                .source()
                .with_voce_hardening(MaterialPropertySelection::SingleClaimOnly)
                .is_err()
        );
    }
    let fixture = Fixture::new(properties(Q0, DQ));
    assert!(
        fixture
            .voce()
            .with_voce_hardening(MaterialPropertySelection::SingleClaimOnly)
            .is_err()
    );
    assert!(
        fixture
            .voce()
            .with_linear_overstress(MaterialPropertySelection::SingleClaimOnly)
            .is_err()
    );
    assert!(
        fixture
            .source()
            .with_linear_overstress(MaterialPropertySelection::SingleClaimOnly)
            .unwrap()
            .with_voce_hardening(MaterialPropertySelection::SingleClaimOnly)
            .is_err()
    );
}

fn config() -> FeedbackConfig {
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
    // Exact uniform, fixed-mass calorimetry. The original phase curve owns h->T;
    // the prescribed external heat is independent of mechanical residuals.
    let enthalpy = input.initial.specific_enthalpy_j_kg()
        + (input.internal_heat_j + external_heat_j) / input.mass_kg;
    Ok(UniformEnthalpyStep {
        state: input.curve.state_at_specific_enthalpy(enthalpy).unwrap(),
        external_heat_j,
        energy_residual_tolerance_j: 1e-8,
        report: (),
    })
}

#[test]
fn original_caloric_feedback_accepts_loaded_heating_once_and_keeps_failed_trials_private() {
    let fixture = Fixture::new(properties(Q0, DQ));
    let initial = CaloricPoint::try_new(fixture.query(), 0.0, 1e-6, 1e-3)
        .unwrap()
        .with_voce_hardening(MaterialPropertySelection::SingleClaimOnly)
        .unwrap();
    let strain = [0.012, -0.006, -0.006, 0.002, 0.0, 0.0];
    let loaded = initial
        .advance(
            strain,
            1.0,
            config(),
            |input| caloric_proposal(input, 0.0),
            || false,
        )
        .unwrap();
    assert!(loaded.next.history().alpha > 0.0);
    assert!(loaded.physical_heat_j > 0.0);
    // With this steep negative Q', reversible hardening entropy absorbs more
    // heat than the positive yield dissipation supplies. Preserve its sign.
    assert!(loaded.reversible_heat_j < -loaded.physical_heat_j);
    near(
        RHO * 1e-6 * CP * (loaded.next.caloric_state().temperature_k() - 300.0),
        loaded.physical_heat_j + loaded.reversible_heat_j,
        1e-8,
    );
    near(loaded.combined_residual_j, 0.0, 1e-8);
    let old_history = loaded.next.history().clone();
    let old_caloric = *loaded.next.caloric_state();
    let heat = loaded.next.mass_kg() * CP * 10.0;
    let mut limited = config();
    limited.max_iterations = 1;
    assert!(matches!(
        loaded.next.advance(
            strain,
            2.0,
            limited,
            |input| caloric_proposal(input, heat),
            || false
        ),
        Err(FeedbackError::NotConverged { .. })
    ));
    assert_eq!(*loaded.next.history(), old_history);
    assert_eq!(*loaded.next.caloric_state(), old_caloric);
    let heated = loaded
        .next
        .advance(
            strain,
            2.0,
            config(),
            |input| caloric_proposal(input, heat),
            || false,
        )
        .unwrap();
    assert!(heated.next.history().alpha > old_history.alpha);
    assert!(heated.next.caloric_state().temperature_k() > old_caloric.temperature_k() + 9.0);
    near(heated.combined_residual_j, 0.0, 1e-8);
    near(heated.thermal_residual_j, 0.0, 1e-8);
    let source = fixture.voce();
    let audited = source
        .advance_with_energy(
            old_caloric.specific_enthalpy_j_kg(),
            heated.next.caloric_state().specific_enthalpy_j_kg(),
            &strain,
            &old_history,
            &strain,
            tolerances(),
        )
        .unwrap();
    assert_eq!(*heated.next.history(), audited.mechanical.state);
    assert_eq!(*heated.next.stress_pa(), audited.mechanical.stress_pa);
    assert!(
        heated
            .next
            .clone()
            .with_voce_hardening(MaterialPropertySelection::SingleClaimOnly)
            .is_err()
    );
    assert!(matches!(
        heated.next.advance(
            strain,
            1.0,
            config(),
            |input| caloric_proposal(input, 0.0),
            || true
        ),
        Err(FeedbackError::Cancelled)
    ));
}
