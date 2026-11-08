//! G0/G1/G3: source-bound anisotropic heating on the original Hill and section owners.
//! The numerical coefficients below are authored fixtures, not measured sheet data.

#[path = "plastic_hill_thermal/mixed.rs"]
mod mixed;
#[path = "plastic_hill_thermal/voce.rs"]
mod voce;
#[path = "plastic_hill_thermal/rate.rs"]
mod rate;

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::hill::HillState;
use fs_material::plastic::hill::source::{
    HILL_COEFFICIENT_PROPERTIES, HILL_HARDENING_MODULUS_PROPERTY,
    HILL_REFERENCE_YIELD_STRESS_PROPERTY,
};
use fs_material::plastic::hill::thermal::{ThermalHillSource, ThermalHillSourceQuery};
use fs_material::plastic::reduction::ReductionConfig;
use fs_material::plastic::section::thermal::ThermoplasticSectionStep;
use fs_material::plastic::section::thermal::source::ThermalSectionSourceFamily;
use fs_material::plastic::section::{PlaneStressLayer, SectionBudget, ThicknessQuadrature};
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_material::{SmallStrainLaw, Voigt};
use fs_qty::{Density, Dims, Pressure, Temperature};

const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const QUARTER_TURN: [[f64; 3]; 3] = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
const E: f64 = 1.2e9;
const E_SLOPE: f64 = -1e7;
const NU: f64 = 0.2;
const H: f64 = 1e8;
const H_SLOPE: f64 = -1e6;
const Y: f64 = 3e6;
const Y_SLOPE: f64 = -1e5;
const CP: f64 = 500.0;
const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);

fn near(actual: f64, expected: f64, allowance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= allowance,
        "actual {actual:e}, expected {expected:e}, allowance {allowance:e}"
    );
}

fn property(name: &str, dims: Dims, value: f64, slope: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: if slope == 0.0 {
            PropertyValue::Scalar { value, dims }
        } else {
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, value - 20.0 * slope), (320.0, value + 20.0 * slope)],
                dims,
            }
        },
        validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
        interpolation: if slope == 0.0 {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        },
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "authored thermoplastic Hill sheet".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

fn properties(expansion: f64, varying_moduli: bool) -> Vec<PropertyClaim> {
    let mut values = vec![
        property(
            YOUNG_MODULUS_PROPERTY,
            Pressure::DIMS,
            E,
            if varying_moduli { E_SLOPE } else { 0.0 },
        ),
        property(POISSON_RATIO_PROPERTY, Dims([0; 6]), NU, 0.0),
        property(
            HILL_REFERENCE_YIELD_STRESS_PROPERTY,
            Pressure::DIMS,
            Y,
            Y_SLOPE,
        ),
        property(
            HILL_HARDENING_MODULUS_PROPERTY,
            Pressure::DIMS,
            H,
            if varying_moduli { H_SLOPE } else { 0.0 },
        ),
        property("specific_heat_capacity", CP_DIMS, CP, 0.0),
        property(DENSITY_PROPERTY, Density::DIMS, 1000.0, 0.0),
        property(
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            INVERSE_TEMPERATURE_DIMS,
            expansion,
            0.0,
        ),
    ];
    values.extend(
        HILL_COEFFICIENT_PROPERTIES
            .into_iter()
            .zip([0.35, 0.8, 0.45, 1.5, 1.1, 2.0])
            .map(|(name, value)| property(name, Dims([0; 6]), value, 0.0)),
    );
    values
}

struct Fixture {
    card: MaterialCard,
    caloric: ResolvedHeatCapacityCurve,
}

impl Fixture {
    fn new(properties: Vec<PropertyClaim>) -> Self {
        let mut claims = ClaimSet::new();
        for claim in properties {
            claims.insert_claim(claim).unwrap();
        }
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "authored rolled sheet".into(),
                phase: "solid".into(),
                process: "thermal anisotropy regression".into(),
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
                heat_capacity: &PropertyKey::new("specific_heat_capacity", CP_DIMS),
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

    fn query(&self, orientation: [[f64; 3]; 3]) -> ThermalHillSourceQuery<'_> {
        ThermalHillSourceQuery {
            card: &self.card,
            caloric: &self.caloric,
            context: QueryPoint::new().with("T", 300.0).unwrap(),
            temperature_axis: "T".into(),
            selection: MaterialPropertySelection::SingleClaimOnly,
            thermal_expansion: Some(MaterialPropertySelection::SingleClaimOnly),
            strain_limit: 0.04,
            material_to_spatial: orientation,
        }
    }

    fn source(&self, orientation: [[f64; 3]; 3]) -> ThermalHillSource<'_> {
        ThermalHillSource::try_new(self.query(orientation), 0.0).unwrap()
    }
}

fn tolerance() -> EnergyTolerances {
    EnergyTolerances {
        energy_j_m3: 1e-5,
        yield_pa: 1e-3,
    }
}

fn h(temperature_k: f64) -> f64 {
    CP * (temperature_k - 300.0)
}

fn dot(a: &Voigt, b: &Voigt) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum::<f64>() + 2.0 * (3..6).map(|i| a[i] * b[i]).sum::<f64>()
}

#[test]
fn constrained_heating_releases_energy_and_advances_the_original_directional_history() {
    let fixture = Fixture::new(properties(0.0, false));
    let mut query = fixture.query(IDENTITY);
    query.thermal_expansion = None;
    let source = ThermalHillSource::try_new(query, 0.0).unwrap();
    let strain = [0.008, -0.003, -0.005, 0.001, 0.0, 0.0];
    let loaded = source
        .section_increment(
            0.0,
            0.0,
            &[0.0; 6],
            &HillState::default(),
            &strain,
            tolerance(),
            0.1,
        )
        .unwrap();
    assert!(loaded.state.alpha > 0.0);
    source
        .validate_section_history(0.0, &strain, &loaded.state, 1e-3)
        .unwrap();
    let saved = loaded.state.clone();
    let endpoint = source.resolve_section_step(0.0, h(320.0), 0.1).unwrap();
    // This old state is physically valid at 300 K and outside the new yield
    // surface. The strict isothermal owner must still reject that rebinding.
    assert!(
        endpoint
            .law
            .resolved()
            .local_law()
            .advance_with_energy(&strain, &loaded.state, &strain, tolerance(),)
            .is_err()
    );
    let heated = source
        .section_increment(
            0.0,
            h(320.0),
            &strain,
            &loaded.state,
            &strain,
            tolerance(),
            0.1,
        )
        .unwrap();
    assert!(heated.state.alpha > loaded.state.alpha);
    assert!(heated.physical_heat_j_m3 > 0.0);
    assert!(heated.numerical_loss_j_m3 > 0.0);
    near(
        heated.physical_heat_j_m3,
        (Y + 20.0 * Y_SLOPE) * (heated.state.alpha - loaded.state.alpha),
        1e-6,
    );
    assert_eq!(heated.mechanical_work_j_m3, 0.0);
    assert_eq!(heated.coefficient_storage_change_j_m3, 0.0);
    assert_eq!(heated.reversible_heat_j_m3, 0.0);
    assert_eq!(heated.thermoelastic_storage_change_j_m3, 0.0);
    near(
        heated.stored_change_j_m3 + heated.physical_heat_j_m3 + heated.numerical_loss_j_m3,
        0.0,
        1e-5,
    );
    assert_eq!(
        heated.state,
        endpoint.law.update_state(&strain, &loaded.state)
    );
    assert_eq!(
        heated.stress_pa,
        endpoint.law.stress(&strain, &loaded.state)
    );
    source
        .validate_section_history(h(320.0), &strain, &heated.state, 1e-3)
        .unwrap();
    assert_eq!(loaded.state, saved);
    let replay = source
        .section_increment(
            0.0,
            h(320.0),
            &strain,
            &loaded.state,
            &strain,
            tolerance(),
            0.1,
        )
        .unwrap();
    assert_eq!(replay.state, heated.state);
    assert_eq!(replay.stress_pa, heated.stress_pa);
}

/// Independent tensor free energy and entropy with the authored affine moduli.
fn endpoint_oracle(
    temperature_k: f64,
    total: &Voigt,
    state: &HillState,
    rotation: &[[f64; 3]; 3],
    expansion: f64,
) -> (Voigt, Voigt, f64, f64) {
    let plastic = fs_material::tensor::rotate(&state.plastic_strain, rotation);
    let elastic: Voigt = core::array::from_fn(|i| {
        total[i]
            - plastic[i]
            - if i < 3 {
                expansion * (temperature_k - 300.0)
            } else {
                0.0
            }
    });
    let youngs = E + E_SLOPE * (temperature_k - 300.0);
    let hardening = H + H_SLOPE * (temperature_k - 300.0);
    let constitutive = |youngs: f64| -> Voigt {
        let lambda = youngs * NU / ((1.0 + NU) * (1.0 - 2.0 * NU));
        let twice_shear = youngs / (1.0 + NU);
        let trace = elastic[0] + elastic[1] + elastic[2];
        core::array::from_fn(|i| {
            twice_shear * elastic[i] + if i < 3 { lambda * trace } else { 0.0 }
        })
    };
    let stress = constitutive(youngs);
    let free_energy = 0.5 * dot(&elastic, &stress) + 0.5 * hardening * state.alpha.powi(2);
    let entropy = expansion * (stress[0] + stress[1] + stress[2])
        - 0.5 * dot(&elastic, &constitutive(E_SLOPE))
        - 0.5 * H_SLOPE * state.alpha.powi(2);
    (elastic, stress, free_energy, entropy)
}

#[test]
fn affine_moduli_expansion_and_rotated_plasticity_obey_the_tensor_first_law() {
    let expansion = 1e-5;
    let fixture = Fixture::new(properties(expansion, true));
    let rotation = fs_material::tensor::rotation([0.0, 0.0, 1.0], 0.37);
    let source = fixture.source(rotation);
    let old_strain = [0.008, -0.001, -0.002, 0.001, 0.0005, -0.0003];
    let new_strain = [0.009, -0.0008, -0.0025, 0.0018, 0.0006, -0.0001];
    let old = source
        .section_increment(
            0.0,
            0.0,
            &[0.0; 6],
            &HillState::default(),
            &old_strain,
            tolerance(),
            0.1,
        )
        .unwrap();
    let temperature = 311.0;
    let next = source
        .section_increment(
            0.0,
            h(temperature),
            &old_strain,
            &old.state,
            &new_strain,
            tolerance(),
            0.2,
        )
        .unwrap();
    let (e0, s0, psi0, entropy0) =
        endpoint_oracle(300.0, &old_strain, &old.state, &rotation, expansion);
    let (_, s1, psi1, entropy1) =
        endpoint_oracle(temperature, &new_strain, &next.state, &rotation, expansion);
    let delta_ts = temperature * entropy1 - 300.0 * entropy0;
    let coefficient = 0.5 * dot(&e0, &s0) * E_SLOPE * (temperature - 300.0) / E
        + 0.5 * H_SLOPE * (temperature - 300.0) * old.state.alpha.powi(2);
    let delta_strain: Voigt = core::array::from_fn(|i| new_strain[i] - old_strain[i]);
    let work = dot(&s1, &delta_strain);
    let reversible =
        expansion * (s1[0] + s1[1] + s1[2]) * (temperature - 300.0) - coefficient - delta_ts;
    for i in 0..6 {
        near(next.stress_pa[i], s1[i], 1e-3);
    }
    near(next.stored_change_j_m3, psi1 - psi0, 1e-5);
    near(next.coefficient_storage_change_j_m3, coefficient, 1e-5);
    near(next.thermoelastic_storage_change_j_m3, delta_ts, 1e-5);
    near(next.reversible_heat_j_m3, reversible, 1e-5);
    near(next.mechanical_work_j_m3, work, 1e-5);
    near(
        next.physical_heat_j_m3,
        (Y + Y_SLOPE * (temperature - 300.0)) * (next.state.alpha - old.state.alpha),
        1e-5,
    );
    near(
        work - (psi1 - psi0)
            - delta_ts
            - next.physical_heat_j_m3
            - reversible
            - next.numerical_loss_j_m3,
        0.0,
        1e-5,
    );
    let endpoint = source
        .resolve_section_step(0.0, h(temperature), 0.2)
        .unwrap();
    assert_eq!(
        endpoint.law.resolved().law().material_to_spatial(),
        &rotation
    );
    for property in endpoint.law.resolved().material().properties() {
        fixture
            .card
            .claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
    assert_eq!(
        endpoint.law.resolved().material().query().axes().get("T"),
        Some(&temperature)
    );
}

#[test]
fn rotated_layers_retain_unequal_histories_and_conservative_heat_during_a_thermal_gradient() {
    let fixture = Fixture::new(properties(1e-5, true));
    let make_layers = |top_orientation| {
        vec![
            PlaneStressLayer {
                bottom_z_m: -0.005,
                top_z_m: 0.0,
                subdivisions: 2,
                law: fixture.source(IDENTITY),
            },
            PlaneStressLayer {
                bottom_z_m: 0.0,
                top_z_m: 0.005,
                subdivisions: 2,
                law: fixture.source(top_orientation),
            },
        ]
    };
    let rotated = make_layers(QUARTER_TURN);
    let aligned = make_layers(IDENTITY);
    let cold = [0.0; 2];
    let hot = [h(305.0), h(315.0)];
    let reduction = ReductionConfig {
        max_iterations: 64,
        stress_tolerance_pa: 1e-4,
        strain_limit: 0.04,
    };
    let budget = SectionBudget {
        maximum_layers: 4,
        maximum_quadrature_points: 48,
    };
    let loading = ThermoplasticSectionStep::new(
        &rotated,
        &cold,
        &cold,
        0.1,
        ThicknessQuadrature::Gauss3,
        reduction,
        budget,
    )
    .unwrap();
    let reference = ThermoplasticSectionStep::new(
        &aligned,
        &cold,
        &cold,
        0.1,
        ThicknessQuadrature::Gauss3,
        reduction,
        budget,
    )
    .unwrap();
    let strain = [0.006, 0.0, 0.0, 0.0, 0.0, 0.0];
    let initial = loading.section().initial_state();
    let loaded = loading.section().evaluate(strain, &initial).unwrap();
    let uniform = reference
        .section()
        .evaluate(strain, &reference.section().initial_state())
        .unwrap();
    near(uniform.resultants[3], 0.0, 1e-8);
    assert!(
        loaded.resultants[3].abs() > 1.0,
        "rotating one layer must change its yielded membrane stress"
    );
    loading
        .audit_energy(&initial, &loaded.state, tolerance(), 1e-6)
        .unwrap();
    let saved = loaded.state.clone();
    let heating = ThermoplasticSectionStep::new(
        &rotated,
        &cold,
        &hot,
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
    assert_eq!(audit.total.generalized_work_j_m2, 0.0);
    assert!(audit.total.physical_heat_j_m2 > 0.0);
    for (layer, temperature) in [305.0, 315.0].into_iter().enumerate() {
        let expected: f64 = heating
            .section()
            .points()
            .iter()
            .enumerate()
            .filter(|(_, point)| point.layer == layer)
            .map(|(i, point)| {
                point.weight_m
                    * (Y + Y_SLOPE * (temperature - 300.0))
                    * (heated.state.points[i].history.alpha - loaded.state.points[i].history.alpha)
            })
            .sum();
        near(audit.layers[layer].physical_heat_j_m2, expected, 1e-7);
        near(audit.layers[layer].balance_residual_j_m2, 0.0, 1e-6);
    }
    near(
        audit.total.caloric_heat_j_m2(),
        audit
            .layers
            .iter()
            .map(|layer| layer.caloric_heat_j_m2())
            .sum(),
        1e-7,
    );
    assert_eq!(
        heated,
        heating.section().evaluate(strain, &loaded.state).unwrap()
    );
    assert_eq!(loaded.state, saved);
    let mut corrupted = heated.state.clone();
    corrupted.points[0].history.alpha += 1e-5;
    assert!(
        heating
            .audit_energy(&loaded.state, &corrupted, tolerance(), 1e-6)
            .is_err()
    );
}

#[test]
fn unowned_metric_changes_kinked_storage_and_incompatible_thermal_sources_refuse() {
    for (name, dims, values) in [
        (
            HILL_COEFFICIENT_PROPERTIES[0],
            Dims([0; 6]),
            [0.35, 0.4, 0.35],
        ),
        (
            YOUNG_MODULUS_PROPERTY,
            Pressure::DIMS,
            [1.4e9, 1.3e9, 1.0e9],
        ),
        (
            HILL_HARDENING_MODULUS_PROPERTY,
            Pressure::DIMS,
            [1.2e8, 1.1e8, 0.8e8],
        ),
        (POISSON_RATIO_PROPERTY, Dims([0; 6]), [0.2, 0.21, 0.2]),
        (DENSITY_PROPERTY, Density::DIMS, [1000.0, 1010.0, 1000.0]),
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            INVERSE_TEMPERATURE_DIMS,
            [1e-5, 2e-5, 1e-5],
        ),
    ] {
        let mut data = properties(1e-5, true);
        let claim = data
            .iter_mut()
            .find(|claim| claim.key.name() == name)
            .unwrap();
        claim.value = PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(280.0, values[0]), (300.0, values[1]), (320.0, values[2])],
            dims,
        };
        claim.interpolation = InterpolationPolicy::LinearInside;
        let fixture = Fixture::new(data);
        assert!(
            ThermalHillSource::try_new(fixture.query(IDENTITY), 0.0).is_err(),
            "must refuse {name}"
        );
    }
    let fixture = Fixture::new(properties(1e-5, true));
    let mut query = fixture.query(IDENTITY);
    query.context = query.context.with("T", 301.0).unwrap();
    assert!(ThermalHillSource::try_new(query, 0.0).is_err());
    let mut query = fixture.query(IDENTITY);
    query.context = query.context.with("pressure", 101325.0).unwrap();
    assert!(ThermalHillSource::try_new(query, 0.0).is_err());
    let other = Fixture::new(properties(0.0, false));
    let mut query = fixture.query(IDENTITY);
    query.caloric = &other.caloric;
    assert!(ThermalHillSource::try_new(query, 0.0).is_err());
    // The total free-expansion excursion stays inside the strain ceiling;
    // the coupled isochoric heat capacity is nevertheless not positive.
    let insufficient_capacity = Fixture::new(properties(4e-4, true));
    assert!(ThermalHillSource::try_new(insufficient_capacity.query(IDENTITY), 0.0).is_err());
}
