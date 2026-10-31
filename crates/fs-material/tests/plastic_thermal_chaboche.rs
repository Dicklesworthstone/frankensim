//! Synthetic G0/G1/G3 thermal-history regressions, not measured calibration.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::SmallStrainLaw;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::plastic::chaboche::ChabocheState;
use fs_material::plastic::chaboche::source::BackstressProperties;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::rate::{
    OVERSTRESS_EXPONENT_PROPERTY, OVERSTRESS_STRESS_SCALE_PROPERTY, OVERSTRESS_VISCOSITY_DIMS,
    OVERSTRESS_VISCOSITY_PROPERTY,
};
use fs_material::plastic::source::{HardeningPropertyBasis, KINEMATIC_HARDENING_MODULUS_PROPERTY};
use fs_material::plastic::thermal::chaboche::ThermalChabocheSource;
use fs_material::plastic::thermal::{SourceQuery, ThermoplasticSource};
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
const Q0: f64 = 6e6;
const B0: f64 = 300.0;
const CP: f64 = 500.0;
const C: [f64; 2] = [4e8, 1e8];
const DC: [f64; 2] = [-8e6, 2e6];
const GAMMA: [f64; 2] = [120.0, 15.0];
const BASIS: HardeningPropertyBasis = HardeningPropertyBasis::MixedPlasticStrainModuli;
const CP_NAME: &str = "specific_heat_capacity";
const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);

fn schema() -> [BackstressProperties; 2] {
    [
        BackstressProperties::new("branch_fast_C", "branch_fast_gamma"),
        BackstressProperties::new("branch_slow_C", "branch_slow_gamma"),
    ]
}
fn scalar(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained(),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic thermal Chaboche source".into(),
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
fn affine(data: &mut [PropertyClaim], name: &str, reference: f64, slope: f64) {
    let property = data.iter_mut().find(|p| p.key.name() == name).unwrap();
    curve(
        property,
        [280.0, 290.0, 300.0, 315.0, 320.0]
            .map(|t| (t, reference + slope * (t - 300.0)))
            .to_vec(),
    );
}
fn properties(slopes: [f64; 2]) -> Vec<PropertyClaim> {
    let mut data = vec![
        scalar(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
        scalar(POISSON_RATIO_PROPERTY, Dims([0; 6]), NU),
        scalar(YIELD_STRESS_PROPERTY, Pressure::DIMS, SY),
        scalar(BASIS.property_name(), Pressure::DIMS, H),
        scalar(KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, HK),
        scalar(CP_NAME, CP_DIMS, CP),
        scalar(DENSITY_PROPERTY, Density::DIMS, 1000.0),
        scalar(Q, Pressure::DIMS, Q0),
        scalar(B, Dims([0; 6]), B0),
    ];
    for (i, names) in schema().iter().enumerate() {
        data.push(scalar(&names.modulus, Pressure::DIMS, C[i]));
        data.push(scalar(&names.recovery, Dims([0; 6]), GAMMA[i]));
        if slopes[i] != 0.0 {
            affine(&mut data, &names.modulus, C[i], slopes[i]);
        }
    }
    data
}
struct Fixture {
    card: MaterialCard,
    caloric: ResolvedHeatCapacityCurve,
    expansion: bool,
}
impl Fixture {
    fn new(data: Vec<PropertyClaim>) -> Self {
        let expansion = data
            .iter()
            .any(|p| p.key.name() == LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY);
        let mut claims = ClaimSet::new();
        for property in data {
            claims.insert_claim(property).unwrap();
        }
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic thermal cyclic solid".into(),
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
        Self {
            card,
            caloric,
            expansion,
        }
    }
    fn base(&self) -> ThermoplasticSource<'_> {
        ThermoplasticSource::try_new(
            SourceQuery {
                card: &self.card,
                caloric: &self.caloric,
                context: QueryPoint::new().with("T", 300.0).unwrap(),
                temperature_axis: "T".into(),
                selection: MaterialPropertySelection::SingleClaimOnly,
                hardening_basis: BASIS,
                thermal_expansion: self
                    .expansion
                    .then_some(MaterialPropertySelection::SingleClaimOnly),
                strain_limit: 0.04,
            },
            0.0,
        )
        .unwrap()
        .with_voce_hardening(MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
    }
    fn binding(&self) -> ThermalChabocheSource<'_, 2> {
        self.base()
            .with_chaboche_branches(MaterialPropertySelection::SingleClaimOnly, schema())
            .unwrap()
    }

    fn rate(&self, power_law: bool) -> ThermalChabocheSource<'_, 2> {
        let base = if power_law {
            self.base()
                .with_power_law_overstress(MaterialPropertySelection::SingleClaimOnly)
        } else {
            self.base()
                .with_linear_overstress(MaterialPropertySelection::SingleClaimOnly)
        }
        .unwrap();
        base.with_chaboche_branches(MaterialPropertySelection::SingleClaimOnly, schema())
            .unwrap()
    }
}

fn rate_properties() -> Vec<PropertyClaim> {
    let mut data = properties(DC);
    data.extend([
        scalar(
            OVERSTRESS_VISCOSITY_PROPERTY,
            OVERSTRESS_VISCOSITY_DIMS,
            2e6,
        ),
        scalar(OVERSTRESS_STRESS_SCALE_PROPERTY, Pressure::DIMS, 1e6),
        scalar(OVERSTRESS_EXPONENT_PROPERTY, Dims([0; 6]), 3.0),
    ]);
    for (name, reference, slope) in [
        (OVERSTRESS_VISCOSITY_PROPERTY, 2e6, -50_000.0),
        (OVERSTRESS_STRESS_SCALE_PROPERTY, 1e6, 10_000.0),
        (OVERSTRESS_EXPONENT_PROPERTY, 3.0, 0.01),
    ] {
        affine(&mut data, name, reference, slope);
    }
    data
}
fn tol() -> EnergyTolerances {
    EnergyTolerances {
        energy_j_m3: 2e-5,
        yield_pa: 1e-3,
    }
}
fn h(t: f64) -> f64 {
    CP * (t - 300.0)
}
fn inner(a: &[f64; 6], b: &[f64; 6]) -> f64 {
    (0..6)
        .map(|i| a[i] * b[i] * if i < 3 { 1.0 } else { 2.0 })
        .sum()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance={tolerance:e}"
    );
}
fn stored(strain: &[f64; 6], state: &ChabocheState<2>, t: f64, variable: bool) -> f64 {
    let dt = t - 300.0;
    let (youngs, hi, hk, q, expansion) = if variable {
        (
            E - 2e6 * dt,
            H - 1e6 * dt,
            HK - 2e6 * dt,
            Q0 - 1e5 * dt,
            1e-5 * dt,
        )
    } else {
        (E, H, HK, Q0, 0.0)
    };
    let elastic: [f64; 6] = core::array::from_fn(|i| {
        strain[i] - state.plastic.plastic_strain[i] - if i < 3 { expansion } else { 0.0 }
    });
    let mu = youngs / (2.0 * (1.0 + NU));
    let lambda = youngs * NU / ((1.0 + NU) * (1.0 - 2.0 * NU));
    let trace = elastic[0] + elastic[1] + elastic[2];
    let alpha = state.plastic.alpha;
    let mut energy = mu * inner(&elastic, &elastic)
        + 0.5 * lambda * trace * trace
        + 0.5 * hi * alpha * alpha
        + hk * inner(&state.plastic.plastic_strain, &state.plastic.plastic_strain) / 3.0
        + q * (alpha + (-B0 * alpha).exp_m1() / B0);
    for i in 0..2 {
        energy +=
            0.75 * inner(&state.backstress_pa[i], &state.backstress_pa[i]) / (C[i] + DC[i] * dt);
    }
    energy
}

#[test]
fn complete_source_pins_and_constant_temperature_native_reduction_are_exact() {
    let fixture = Fixture::new(properties(DC));
    let source = fixture.binding();
    let (_, law) = source.resolve_enthalpy(0.0, 0.0).unwrap();
    assert_eq!(source.properties(), &schema());
    assert_eq!(law.resolved().material().properties().len(), 13);
    for property in law.resolved().material().properties() {
        fixture
            .card
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
        assert_eq!(
            property.answer().receipt.selected,
            source
                .base()
                .reference()
                .material()
                .property(property.requirement().name())
                .unwrap()
                .answer()
                .receipt
                .selected
        );
    }
    let mut previous = [0.0; 6];
    let mut state = law.initial_state();
    for strain in [
        [0.01, -0.005, -0.005, 0.002, 0.0, 0.0],
        [-0.009, 0.0045, 0.0045, 0.003, 0.001, 0.0],
    ] {
        let native = law
            .law()
            .advance_with_energy(&previous, &state, &strain, tol())
            .unwrap();
        let next = source
            .advance_with_energy(0.0, 0.0, &previous, &state, &strain, tol())
            .unwrap();
        assert_eq!(native.state, next.state);
        assert_eq!(native.stress_pa, next.thermal.mechanical.stress_pa);
        assert_eq!(native.energy, next.thermal.mechanical.energy);
        assert_eq!(
            law.law().tangent(&strain, &state),
            law.tangent(&strain, &state)
        );
        previous = strain;
        state = next.state;
    }
}

#[test]
fn thermal_transport_preserves_strain_like_history_without_spurious_flow_or_heat() {
    let fixture = Fixture::new(properties(DC));
    let source = fixture.binding();
    let (_, law) = source.resolve_enthalpy(0.0, 0.0).unwrap();
    let loading = [0.008, -0.004, -0.004, 0.001, 0.0, 0.0];
    let loaded = source
        .advance_with_energy(0.0, 0.0, &[0.0; 6], &law.initial_state(), &loading, tol())
        .unwrap();
    let zero_stress_strain = loaded.state.plastic.plastic_strain;
    let unloaded = source
        .advance_with_energy(
            0.0,
            0.0,
            &loading,
            &loaded.state,
            &zero_stress_strain,
            tol(),
        )
        .unwrap();
    let old = unloaded.state.clone();
    let hot = source
        .advance_with_energy(
            0.0,
            h(315.0),
            &zero_stress_strain,
            &old,
            &zero_stress_strain,
            tol(),
        )
        .unwrap();
    assert_eq!(old, unloaded.state);
    assert_eq!(hot.state.plastic, old.plastic);
    let mut coefficient = 0.0;
    let mut entropy = 0.0;
    for i in 0..2 {
        let ratio = (C[i] + DC[i] * 15.0) / C[i];
        for j in 0..6 {
            near(
                hot.state.backstress_pa[i][j],
                old.backstress_pa[i][j] * ratio,
                2e-9,
            );
        }
        let shape = 0.75 * inner(&old.backstress_pa[i], &old.backstress_pa[i]) / (C[i] * C[i]);
        coefficient += DC[i] * 15.0 * shape;
        entropy -= DC[i] * shape;
    }
    near(
        hot.thermal.coefficient_storage_change_j_m3,
        coefficient,
        1e-10,
    );
    near(
        hot.thermal.thermoelastic_storage_change_j_m3,
        15.0 * entropy,
        1e-9,
    );
    near(
        hot.thermal.mechanical.energy.physical_dissipation_j_m3,
        0.0,
        0.0,
    );
    near(
        hot.thermal.mechanical.energy.numerical_dissipation_j_m3,
        0.0,
        0.0,
    );
    near(hot.thermal.caloric_heat_j_m3(), 0.0, 1e-9);
    let cold = source
        .advance_with_energy(
            h(315.0),
            0.0,
            &zero_stress_strain,
            &hot.state,
            &zero_stress_strain,
            tol(),
        )
        .unwrap();
    assert_eq!(cold.state.plastic, old.plastic);
    for i in 0..2 {
        for j in 0..6 {
            near(
                cold.state.backstress_pa[i][j],
                old.backstress_pa[i][j],
                2e-9,
            );
        }
    }
}

#[test]
fn hot_nonproportional_tangent_differentiates_transported_frozen_history() {
    let fixture = Fixture::new(properties(DC));
    let source = fixture.binding();
    let (_, cold) = source.resolve_enthalpy(0.0, 0.0).unwrap();
    let loaded = [0.011, -0.0055, -0.0055, 0.001, 0.0, 0.0];
    let state = cold.update_state(&loaded, &cold.initial_state());
    let (_, hot) = source.resolve_enthalpy(0.0, h(315.0)).unwrap();
    assert!(!hot.admissibility().tangent_symmetric);
    let strain = [-0.007, 0.002, 0.005, 0.004, -0.002, 0.001];
    let response = hot.response(&strain, &state).unwrap();
    assert_eq!(response.state, hot.update_state(&strain, &state));
    assert_eq!(response.stress_pa, hot.stress(&strain, &state));
    let step = 2e-8;
    for j in 0..6 {
        let mut plus = strain;
        let mut minus = strain;
        plus[j] += step;
        minus[j] -= step;
        let a = hot.stress(&plus, &state);
        let b = hot.stress(&minus, &state);
        for i in 0..6 {
            near(
                response.tangent_pa[i][j],
                (a[i] - b[i]) / (2.0 * step),
                15.0,
            );
        }
    }
    let next = source
        .advance_with_energy(0.0, h(315.0), &loaded, &state, &strain, tol())
        .unwrap();
    assert_eq!(next.state, response.state);
    assert_eq!(next.thermal.mechanical.stress_pa, response.stress_pa);
}

#[test]
fn simultaneous_affine_backbone_expansion_and_branch_energy_close_on_heated_reversal() {
    let mut data = properties(DC);
    for (name, value, slope) in [
        (YOUNG_MODULUS_PROPERTY, E, -2e6),
        (BASIS.property_name(), H, -1e6),
        (KINEMATIC_HARDENING_MODULUS_PROPERTY, HK, -2e6),
        (Q, Q0, -1e5),
    ] {
        affine(&mut data, name, value, slope);
    }
    data.push(scalar(
        LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
        INVERSE_TEMPERATURE_DIMS,
        1e-5,
    ));
    let fixture = Fixture::new(data);
    let source = fixture.binding();
    let (_, cold) = source.resolve_enthalpy(0.0, 0.0).unwrap();
    let old_strain = [0.012, -0.006, -0.006, 0.002, 0.0, 0.0];
    let state = cold.update_state(&old_strain, &cold.initial_state());
    let strain = [-0.008, 0.004, 0.004, -0.001, 0.003, 0.0];
    let next = source
        .advance_with_energy(0.0, h(315.0), &old_strain, &state, &strain, tol())
        .unwrap();
    let thermal = &next.thermal;
    near(
        thermal.stored_change_j_m3(),
        stored(&strain, &next.state, 315.0, true) - stored(&old_strain, &state, 300.0, true),
        2e-7,
    );
    let de: [f64; 6] = core::array::from_fn(|i| strain[i] - old_strain[i]);
    near(
        thermal.mechanical_work_j_m3,
        inner(&thermal.mechanical.stress_pa, &de),
        1e-9,
    );
    assert!(next.recovery_dissipation_j_m3 > 0.0);
    assert!(thermal.mechanical.energy.physical_dissipation_j_m3 > next.recovery_dissipation_j_m3);
    assert!(thermal.mechanical.energy.numerical_dissipation_j_m3 > 0.0);
    near(
        thermal.mechanical_work_j_m3
            - thermal.stored_change_j_m3()
            - thermal.thermoelastic_storage_change_j_m3
            - thermal.caloric_heat_j_m3()
            - thermal.mechanical.energy.numerical_dissipation_j_m3,
        0.0,
        2e-5,
    );
}

#[test]
fn heating_cooling_and_reversal_retain_each_history_and_change_the_cyclic_response() {
    let variable_fixture = Fixture::new(properties(DC));
    let fixed_fixture = Fixture::new(properties([0.0; 2]));
    let variable = variable_fixture.binding();
    let fixed = fixed_fixture.binding();
    let (_, initial) = variable.resolve_enthalpy(0.0, 0.0).unwrap();
    let mut state = initial.initial_state();
    let mut fixed_state = state.clone();
    let mut old_strain = [0.0; 6];
    let mut old_h = 0.0;
    let mut difference: f64 = 0.0;
    let mut plastic_steps = 0;
    for (t, axial, shear) in [
        (300.0, 0.011, 0.0),
        (315.0, 0.003, 0.002),
        (305.0, -0.009, 0.003),
        (285.0, 0.007, -0.002),
        (300.0, -0.007, 0.001),
    ] {
        let strain = [axial, -0.5 * axial, -0.5 * axial, shear, 0.0, 0.0];
        let next = variable
            .advance_with_energy(old_h, h(t), &old_strain, &state, &strain, tol())
            .unwrap();
        let frozen = fixed
            .advance_with_energy(old_h, h(t), &old_strain, &fixed_state, &strain, tol())
            .unwrap();
        assert!(next.state.plastic.alpha >= state.plastic.alpha);
        if next.state.plastic.alpha > state.plastic.alpha {
            plastic_steps += 1;
        }
        assert!(next.thermal.mechanical.energy.physical_dissipation_j_m3 >= 0.0);
        assert!(next.recovery_dissipation_j_m3 >= 0.0);
        difference = difference.max(
            (next.thermal.mechanical.stress_pa[0] - frozen.thermal.mechanical.stress_pa[0]).abs(),
        );
        state = next.state;
        fixed_state = frozen.state;
        old_strain = strain;
        old_h = h(t);
    }
    assert!(
        plastic_steps >= 3,
        "only {plastic_steps} plastic loading intervals"
    );
    eprintln!(
        "thermal cycle alpha={:.12e}, maximum stress effect={difference:.6e} Pa",
        state.plastic.alpha
    );
    assert!(
        difference > 20_000.0,
        "thermal branch effect was {difference} Pa"
    );
    assert_ne!(state.backstress_pa[0], state.backstress_pa[1]);
}

#[test]
fn complete_source_span_rejects_kinks_zero_crossing_variable_recovery_and_aliases() {
    for case in 0..4 {
        let mut data = properties(DC);
        let names = schema();
        let name = if case == 3 {
            &names[0].recovery
        } else {
            &names[0].modulus
        };
        let property = data.iter_mut().find(|p| p.key.name() == name).unwrap();
        curve(
            property,
            match case {
                0 => vec![(280.0, 4e8), (300.0, 2e8), (320.0, 4e8)],
                1 => vec![(280.0, 4e8), (320.0, 0.0)],
                2 => vec![(280.0, 4e8), (300.0, -1.0), (320.0, 4e8)],
                _ => vec![(280.0, 120.0), (300.0, 130.0), (320.0, 120.0)],
            },
        );
        let fixture = Fixture::new(data);
        assert!(
            fixture
                .base()
                .with_chaboche_branches(MaterialPropertySelection::SingleClaimOnly, schema())
                .is_err(),
            "admitted case {case}"
        );
    }
    let fixture = Fixture::new(properties(DC));
    for alias in [YOUNG_MODULUS_PROPERTY, CP_NAME, DENSITY_PROPERTY, Q] {
        let mut names = schema();
        names[0].modulus = alias.into();
        assert!(
            fixture
                .base()
                .with_chaboche_branches(MaterialPropertySelection::SingleClaimOnly, names)
                .is_err()
        );
    }
    let mut names = schema();
    names[1].modulus = names[0].modulus.clone();
    assert!(
        fixture
            .base()
            .with_chaboche_branches(MaterialPropertySelection::SingleClaimOnly, names)
            .is_err()
    );
    assert!(
        fixture
            .base()
            .with_chaboche_branches(MaterialPropertySelection::SingleClaimOnly, [])
            .is_err()
    );
    let mut data = properties(DC);
    data.push(scalar(
        OVERSTRESS_VISCOSITY_PROPERTY,
        OVERSTRESS_VISCOSITY_DIMS,
        1e6,
    ));
    let rate = Fixture::new(data);
    let rate_source = rate
        .base()
        .with_linear_overstress(MaterialPropertySelection::SingleClaimOnly)
        .unwrap();
    let rate_branches = rate_source
        .with_chaboche_branches(MaterialPropertySelection::SingleClaimOnly, schema())
        .unwrap();
    assert!(rate_branches.resolve_enthalpy(0.0, h(315.0)).is_err());
    assert!(
        rate_branches
            .resolve_enthalpy_for_step(0.0, h(315.0), 1e-3)
            .is_ok()
    );
}

#[test]
fn absent_branches_reduce_exactly_and_invalid_states_never_return_valid_raw_results() {
    let mut data = properties([0.0; 2]);
    for names in schema() {
        let p = data
            .iter_mut()
            .find(|p| p.key.name() == names.modulus)
            .unwrap();
        p.value = PropertyValue::Scalar {
            value: 0.0,
            dims: Pressure::DIMS,
        };
    }
    let fixture = Fixture::new(data);
    let source = fixture.binding();
    let (_, law) = source.resolve_enthalpy(0.0, h(315.0)).unwrap();
    let strain = [0.013, -0.0065, -0.0065, 0.002, 0.0, 0.0];
    let state = law.initial_state();
    let branch = source
        .advance_with_energy(0.0, h(315.0), &[0.0; 6], &state, &strain, tol())
        .unwrap();
    let base = source
        .base()
        .advance_with_energy(0.0, h(315.0), &[0.0; 6], &state.plastic, &strain, tol())
        .unwrap();
    assert_eq!(branch.thermal.mechanical, base.mechanical);
    assert_eq!(branch.state.plastic, base.mechanical.state);
    assert_eq!(branch.state.backstress_pa, [[0.0; 6]; 2]);
    assert_eq!(branch.thermal.caloric_heat_j_m3(), base.caloric_heat_j_m3());
    let mut bad = state;
    bad.backstress_pa[0][0] = 1.0;
    assert!(law.response(&strain, &bad).is_err());
    assert!(law.stress(&strain, &bad).iter().all(|v| v.is_nan()));
    assert!(law.update_state(&strain, &bad).plastic.alpha.is_nan());
    assert!(source.resolve_enthalpy(0.0, h(321.0)).is_err());
    assert!(law.response(&[0.1; 6], &law.initial_state()).is_err());
}

#[test]
fn joint_power_law_thermal_branches_use_one_clock_and_close_native_and_reversible_energy() {
    let fixture = Fixture::new(rate_properties());
    let source = fixture.rate(true);
    let (_, cold) = source.resolve_enthalpy_for_step(0.0, 0.0, 4e-4).unwrap();
    let old_strain = [0.011, -0.0055, -0.0055, 0.001, 0.0, 0.0];
    let loaded = source
        .advance_with_energy_for_step(
            0.0,
            0.0,
            &[0.0; 6],
            &cold.initial_state(),
            &old_strain,
            tol(),
            4e-4,
        )
        .unwrap();
    let state = loaded.state;
    let dt = 7e-4;
    let t = 315.0;
    let (_, hot) = source.resolve_enthalpy_for_step(0.0, h(t), dt).unwrap();
    let rate = hot.rate_law().unwrap();
    near(rate.duration_s(), dt, 0.0);
    near(rate.viscosity_pa_s(), 1.25e6, 0.0);
    near(rate.exponent(), 3.15, 1e-14);
    near(rate.stress_scale_pa().unwrap(), 1.15e6, 0.0);
    assert!(!hot.admissibility().tangent_symmetric);
    for property in hot.resolved().material().properties() {
        fixture
            .card
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
        assert_eq!(
            property.answer().receipt.selected,
            source
                .base()
                .reference()
                .material()
                .property(property.requirement().name())
                .unwrap()
                .answer()
                .receipt
                .selected
        );
    }
    let strain = [-0.009, 0.0045, 0.0045, 0.004, -0.002, 0.001];
    let response = hot.response(&strain, &state).unwrap();
    let frozen = state.clone();
    assert_eq!(response.state, hot.update_state(&strain, &state));
    assert_eq!(response.stress_pa, hot.stress(&strain, &state));
    assert_eq!(state, frozen);
    let step = 2e-8;
    for j in 0..6 {
        let mut plus = strain;
        let mut minus = strain;
        plus[j] += step;
        minus[j] -= step;
        let a = hot.stress(&plus, &state);
        let b = hot.stress(&minus, &state);
        for i in 0..6 {
            near(
                response.tangent_pa[i][j],
                (a[i] - b[i]) / (2.0 * step),
                20.0,
            );
        }
    }
    let next = source
        .advance_with_energy_for_step(0.0, h(t), &old_strain, &state, &strain, tol(), dt)
        .unwrap();
    assert_eq!(next.state, response.state);
    assert_eq!(next.thermal.mechanical.stress_pa, response.stress_pa);
    let delta = next.state.plastic.alpha - state.plastic.alpha;
    assert!(delta > 0.0);
    let endpoint_overstress = 1.15e6 * (1.25e6 * delta / (dt * 1.15e6)).powf(1.0 / 3.15);
    near(
        hot.law()
            .yield_function_with_state(&response.stress_pa, &next.state),
        endpoint_overstress,
        1e-3,
    );
    let mut recovery = 0.0;
    let mut coefficient = 0.0;
    let mut s0 = 0.0;
    let mut s1 = 0.0;
    for i in 0..2 {
        let c1 = C[i] + DC[i] * (t - 300.0);
        let a0 = inner(&state.backstress_pa[i], &state.backstress_pa[i]);
        let a1 = inner(&next.state.backstress_pa[i], &next.state.backstress_pa[i]);
        let shape0 = 0.75 * a0 / (C[i] * C[i]);
        let shape1 = 0.75 * a1 / (c1 * c1);
        recovery += 1.5 * GAMMA[i] * delta * a1 / c1;
        coefficient += (c1 - C[i]) * shape0;
        s0 -= DC[i] * shape0;
        s1 -= DC[i] * shape1;
    }
    near(next.recovery_dissipation_j_m3, recovery, 2e-8);
    near(
        next.thermal.mechanical.energy.physical_dissipation_j_m3,
        (SY + endpoint_overstress) * delta + recovery,
        2e-6,
    );
    near(
        next.thermal.coefficient_storage_change_j_m3,
        coefficient,
        2e-8,
    );
    near(
        next.thermal.thermoelastic_storage_change_j_m3,
        t * s1 - 300.0 * s0,
        2e-8,
    );
    near(
        next.thermal.stored_change_j_m3(),
        stored(&strain, &next.state, t, false) - stored(&old_strain, &state, 300.0, false),
        2e-7,
    );
    near(
        next.thermal.mechanical_work_j_m3
            - next.thermal.stored_change_j_m3()
            - next.thermal.thermoelastic_storage_change_j_m3
            - next.thermal.caloric_heat_j_m3()
            - next.thermal.mechanical.energy.numerical_dissipation_j_m3,
        0.0,
        tol().energy_j_m3,
    );
}

#[test]
fn joint_rate_relaxes_held_hot_strain_and_requires_the_actual_duration() {
    let fixture = Fixture::new(rate_properties());
    for power in [false, true] {
        let source = fixture.rate(power);
        let (_, cold) = source.resolve_enthalpy_for_step(0.0, 0.0, 2e-4).unwrap();
        let strain = [0.012, -0.006, -0.006, 0.0, 0.0, 0.0];
        let loaded = source
            .advance_with_energy_for_step(
                0.0,
                0.0,
                &[0.0; 6],
                &cold.initial_state(),
                &strain,
                tol(),
                2e-4,
            )
            .unwrap();
        let old = loaded.state.clone();
        let short = source
            .advance_with_energy_for_step(0.0, h(315.0), &strain, &old, &strain, tol(), 1e-5)
            .unwrap();
        let long = source
            .advance_with_energy_for_step(0.0, h(315.0), &strain, &old, &strain, tol(), 1e-2)
            .unwrap();
        assert!(short.state.plastic.alpha > old.plastic.alpha);
        assert!(long.state.plastic.alpha > short.state.plastic.alpha);
        assert!(long.thermal.mechanical.stress_pa[0] < short.thermal.mechanical.stress_pa[0]);
        assert!(long.thermal.mechanical.energy.physical_dissipation_j_m3 > 0.0);
        assert_eq!(loaded.state, old);
        let (_, trial) = source
            .resolve_enthalpy_for_step(0.0, h(315.0), 1e-5)
            .unwrap();
        for _ in 0..3 {
            assert_eq!(trial.update_state(&strain, &old), short.state);
            assert_eq!(
                trial.stress(&strain, &old),
                short.thermal.mechanical.stress_pa
            );
        }
        assert!(
            source
                .advance_with_energy(0.0, h(315.0), &strain, &old, &strain, tol())
                .is_err()
        );
        for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(source.resolve_enthalpy_for_step(0.0, h(315.0), dt).is_err());
            assert!(
                source
                    .advance_with_energy_for_step(0.0, h(315.0), &strain, &old, &strain, tol(), dt)
                    .is_err()
            );
        }
    }
}

#[test]
fn absent_rate_branches_recover_thermal_voce_and_power_one_recovers_linear_exactly() {
    for exponent in [1.0, 3.0] {
        let mut data = rate_properties();
        affine(&mut data, OVERSTRESS_EXPONENT_PROPERTY, exponent, 0.0);
        for names in schema() {
            affine(&mut data, &names.modulus, 0.0, 0.0);
        }
        let fixture = Fixture::new(data);
        for power in [false, true] {
            let source = fixture.rate(power);
            let (_, law) = source
                .resolve_enthalpy_for_step(0.0, h(315.0), 8e-4)
                .unwrap();
            let strain = [0.013, -0.0065, -0.0065, 0.002, 0.0, 0.0];
            let state = law.initial_state();
            let branch = source
                .advance_with_energy_for_step(
                    0.0,
                    h(315.0),
                    &[0.0; 6],
                    &state,
                    &strain,
                    tol(),
                    8e-4,
                )
                .unwrap();
            let base = source
                .base()
                .advance_with_energy_for_step(
                    0.0,
                    h(315.0),
                    &[0.0; 6],
                    &state.plastic,
                    &strain,
                    tol(),
                    8e-4,
                )
                .unwrap();
            assert_eq!(branch.thermal.mechanical, base.mechanical);
            assert_eq!(branch.state.plastic, base.mechanical.state);
            assert_eq!(branch.state.backstress_pa, [[0.0; 6]; 2]);
            assert_eq!(branch.thermal.caloric_heat_j_m3(), base.caloric_heat_j_m3());
            if exponent == 1.0 && power {
                let linear = fixture.rate(false);
                let next = linear
                    .advance_with_energy_for_step(
                        0.0,
                        h(315.0),
                        &[0.0; 6],
                        &state,
                        &strain,
                        tol(),
                        8e-4,
                    )
                    .unwrap();
                assert_eq!(next.state, branch.state);
                assert_eq!(next.thermal.mechanical, branch.thermal.mechanical);
            }
        }
    }
}
