//! Synthetic thermal/source regressions, not measured material calibration.
use super::near;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::plastic::chaboche::source::BackstressProperties;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::rate::{OVERSTRESS_VISCOSITY_DIMS, OVERSTRESS_VISCOSITY_PROPERTY};
use fs_material::plastic::reduction::ReductionConfig;
use fs_material::plastic::section::thermal::{ThermoplasticSectionError, ThermoplasticSectionStep};
use fs_material::plastic::section::{PlaneStressLayer, SectionBudget, ThicknessQuadrature};
use fs_material::plastic::source::{HardeningPropertyBasis, KINEMATIC_HARDENING_MODULUS_PROPERTY};
use fs_material::plastic::thermal::{SourceQuery, ThermoplasticSource};
use fs_material::plastic::voce::source::{
    VOCE_SATURATION_INCREMENT_PROPERTY, VOCE_SATURATION_RATE_PROPERTY,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY,
    YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Density, Dims, Pressure, Temperature};

const E: f64 = 1.2e9;
const NU: f64 = 0.2;
const CP: f64 = 500.0;
const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);
const BASIS: HardeningPropertyBasis = HardeningPropertyBasis::MixedPlasticStrainModuli;

fn property(name: &str, dims: Dims, value: f64, slope: f64) -> PropertyClaim {
    let varying = slope != 0.0;
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: if varying {
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, value - 20.0 * slope), (320.0, value + 20.0 * slope)],
                dims,
            }
        } else {
            PropertyValue::Scalar { value, dims }
        },
        validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
        interpolation: if varying {
            InterpolationPolicy::LinearInside
        } else {
            InterpolationPolicy::ConstantWithinValidity
        },
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic layered thermal solid".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

struct Fixture {
    card: MaterialCard,
    caloric: ResolvedHeatCapacityCurve,
}
impl Fixture {
    fn new(expansion: f64) -> Self {
        let mut claims = ClaimSet::new();
        for claim in [
            property(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E, -1e7),
            property(POISSON_RATIO_PROPERTY, Dims([0; 6]), NU, 0.0),
            property(YIELD_STRESS_PROPERTY, Pressure::DIMS, 3e6, -1e5),
            property(BASIS.property_name(), Pressure::DIMS, 1e8, -1e6),
            property(
                KINEMATIC_HARDENING_MODULUS_PROPERTY,
                Pressure::DIMS,
                1.6e8,
                -5e5,
            ),
            property("specific_heat_capacity", CP_DIMS, CP, 0.0),
            property(DENSITY_PROPERTY, Density::DIMS, 1000.0, 0.0),
            property(
                LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
                INVERSE_TEMPERATURE_DIMS,
                expansion,
                0.0,
            ),
            property(
                VOCE_SATURATION_INCREMENT_PROPERTY,
                Pressure::DIMS,
                6e6,
                -1e5,
            ),
            property(VOCE_SATURATION_RATE_PROPERTY, Dims([0; 6]), 300.0, 0.0),
            property(
                OVERSTRESS_VISCOSITY_PROPERTY,
                OVERSTRESS_VISCOSITY_DIMS,
                1e7,
                0.0,
            ),
            property("fast_C", Pressure::DIMS, 4e8, -8e6),
            property("fast_gamma", Dims([0; 6]), 120.0, 0.0),
            property("slow_C", Pressure::DIMS, 1e8, 2e6),
            property("slow_gamma", Dims([0; 6]), 15.0, 0.0),
        ] {
            claims.insert_claim(claim).unwrap();
        }
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic thermal layered solid".into(),
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

    fn source(&self) -> ThermoplasticSource<'_> {
        ThermoplasticSource::try_new(
            SourceQuery {
                card: &self.card,
                caloric: &self.caloric,
                context: QueryPoint::new().with("T", 300.0).unwrap(),
                temperature_axis: "T".into(),
                selection: MaterialPropertySelection::SingleClaimOnly,
                hardening_basis: BASIS,
                thermal_expansion: Some(MaterialPropertySelection::SingleClaimOnly),
                strain_limit: 0.04,
            },
            0.0,
        )
        .unwrap()
    }
}

fn reduction() -> ReductionConfig {
    ReductionConfig {
        max_iterations: 64,
        stress_tolerance_pa: 1e-4,
        strain_limit: 0.04,
    }
}
fn budget() -> SectionBudget {
    SectionBudget {
        maximum_layers: 4,
        maximum_quadrature_points: 48,
    }
}
fn tolerance() -> EnergyTolerances {
    EnergyTolerances {
        energy_j_m3: 1e-5,
        yield_pa: 1e-3,
    }
}
fn h(t: f64) -> f64 {
    CP * (t - 300.0)
}
fn layers<S: Clone>(source: S) -> Vec<PlaneStressLayer<S>> {
    vec![
        PlaneStressLayer {
            bottom_z_m: -0.005,
            top_z_m: 0.0,
            subdivisions: 2,
            law: source.clone(),
        },
        PlaneStressLayer {
            bottom_z_m: 0.0,
            top_z_m: 0.005,
            subdivisions: 2,
            law: source,
        },
    ]
}

#[test]
fn thermal_gradient_recovers_plane_stress_force_moment_and_reciprocal_work() {
    let fixture = Fixture::new(1e-5);
    let layers = layers(fixture.source());
    let old_h = [0.0; 2];
    for (next_h, sign) in [([h(300.0), h(310.0)], 1.0), ([h(310.0), h(300.0)], -1.0)] {
        let step = ThermoplasticSectionStep::new(
            &layers,
            &old_h,
            &next_h,
            0.1,
            ThicknessQuadrature::Gauss3,
            reduction(),
            budget(),
        )
        .unwrap();
        let initial = step.section().initial_state();
        let next = step.section().evaluate([0.0; 6], &initial).unwrap();
        let audit = step
            .audit_energy(&initial, &next.state, tolerance(), 1e-6)
            .unwrap();
        let stress = -(E - 1e8) * 1e-5 * 10.0 / (1.0 - NU);
        for i in 0..2 {
            near(next.resultants[i], 0.005 * stress, 1e-8);
            near(
                next.resultants[i + 3],
                sign * 0.5 * 0.005_f64.powi(2) * stress,
                1e-10,
            );
        }
        assert_eq!(audit.total.resultants, next.resultants);
        assert_eq!(audit.total.physical_heat_j_m2, 0.0);
        assert!(audit.total.reversible_heat_j_m2 != 0.0);
        assert!(audit.total.thermoelastic_storage_change_j_m2 != 0.0);
        near(audit.total.balance_residual_j_m2, 0.0, 1e-6);
        assert_eq!(initial, step.section().initial_state());
    }
}

#[test]
fn heated_yielded_layers_preserve_local_heat_storage_and_exact_native_history() {
    let fixture = Fixture::new(1e-5);
    let layers = layers(fixture.source());
    let initial_h = [0.0; 2];
    let warm_h = [h(305.0), h(310.0)];
    let hot_h = [h(315.0), h(320.0)];
    let strain = [0.005, 0.0, 0.001, 0.2, 0.0, 0.03];
    let loading = ThermoplasticSectionStep::new(
        &layers,
        &initial_h,
        &warm_h,
        0.1,
        ThicknessQuadrature::Gauss3,
        reduction(),
        budget(),
    )
    .unwrap();
    let virgin = loading.section().initial_state();
    let loaded = loading.section().evaluate(strain, &virgin).unwrap();
    loading
        .audit_energy(&virgin, &loaded.state, tolerance(), 1e-6)
        .unwrap();
    let saved = loaded.state.clone();
    let heating = ThermoplasticSectionStep::new(
        &layers,
        &warm_h,
        &hot_h,
        0.2,
        ThicknessQuadrature::Gauss3,
        reduction(),
        budget(),
    )
    .unwrap();
    let heated = heating.section().evaluate(strain, &loaded.state).unwrap();
    let audit = heating
        .audit_energy(&loaded.state, &heated.state, tolerance(), 1e-6)
        .unwrap();
    assert_eq!(audit.total.resultants, heated.resultants);
    assert!(audit.total.coefficient_storage_change_j_m2.abs() > 1e-6);
    assert!(audit.total.thermoelastic_storage_change_j_m2.abs() > 1e-6);
    assert!(audit.total.physical_heat_j_m2 > 0.0);
    assert!(audit.total.numerical_loss_j_m2 > 0.0);
    assert_eq!(audit.total.generalized_work_j_m2, 0.0);
    for (layer, energy) in audit.layers.iter().enumerate() {
        let yield_stress = heating.section().layers()[layer]
            .law
            .resolved()
            .law()
            .yield_stress;
        let expected: f64 = heating
            .section()
            .points()
            .iter()
            .enumerate()
            .filter(|(_, point)| point.layer == layer)
            .map(|(i, point)| {
                point.weight_m
                    * yield_stress
                    * (heated.state.points[i].history.alpha - loaded.state.points[i].history.alpha)
            })
            .sum();
        near(energy.physical_heat_j_m2, expected, 1e-8);
        near(energy.balance_residual_j_m2, 0.0, 1e-6);
    }
    near(
        audit.total.caloric_heat_j_m2(),
        audit.layers.iter().map(|e| e.caloric_heat_j_m2()).sum(),
        1e-8,
    );
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
    let mut wrong_geometry = heated.state.clone();
    wrong_geometry.generalized_strain[0] += 1e-5;
    assert!(
        heating
            .audit_energy(&loaded.state, &wrong_geometry, tolerance(), 1e-6)
            .is_err()
    );
    assert_eq!(loaded.state, saved);
}

#[test]
fn timed_thermal_chaboche_retains_every_branch_and_one_physical_interval() {
    let fixture = Fixture::new(1e-5);
    let source = fixture
        .source()
        .with_voce_hardening(MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
        .with_linear_overstress(MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
        .with_chaboche_branches(
            MaterialPropertySelection::SingleClaimOnly,
            [
                BackstressProperties::new("fast_C", "fast_gamma"),
                BackstressProperties::new("slow_C", "slow_gamma"),
            ],
        )
        .unwrap();
    let layers = layers(source);
    let cold = [0.0; 2];
    let hot = [h(310.0); 2];
    let strain = [0.012, -0.002, 0.001, 0.0, 0.0, 0.0];
    let loading = ThermoplasticSectionStep::new(
        &layers,
        &cold,
        &hot,
        2e-4,
        ThicknessQuadrature::Gauss2,
        reduction(),
        budget(),
    )
    .unwrap();
    let virgin = loading.section().initial_state();
    let loaded = loading.section().evaluate(strain, &virgin).unwrap();
    loading
        .audit_energy(&virgin, &loaded.state, tolerance(), 1e-6)
        .unwrap();
    let saved = loaded.state.clone();
    let mut accumulated = vec![];
    for duration in [1e-4, 1e-2] {
        let hold = ThermoplasticSectionStep::new(
            &layers,
            &hot,
            &hot,
            duration,
            ThicknessQuadrature::Gauss2,
            reduction(),
            budget(),
        )
        .unwrap();
        assert_eq!(hold.duration_s(), duration);
        let held = hold.section().evaluate(strain, &loaded.state).unwrap();
        let audit = hold
            .audit_energy(&loaded.state, &held.state, tolerance(), 1e-6)
            .unwrap();
        accumulated.push(held.state.points[0].history.plastic.alpha);
        assert!(audit.total.physical_heat_j_m2 > 0.0);
        assert_eq!(audit.total.generalized_work_j_m2, 0.0);
        assert!(
            held.state.points[0]
                .history
                .backstress_pa
                .iter()
                .all(|branch| branch[0] > 0.0)
        );
        assert_ne!(
            held.state.points[0].history.backstress_pa[0],
            held.state.points[0].history.backstress_pa[1]
        );
        assert!(
            hold.audit_energy(&loaded.state, &loaded.state, tolerance(), 1e-6)
                .is_err()
        );
        assert_eq!(
            held,
            hold.section().evaluate(strain, &loaded.state).unwrap()
        );
    }
    assert!(accumulated[1] > accumulated[0]);
    assert_eq!(loaded.state, saved);
}

#[test]
fn interval_source_geometry_and_budget_mismatches_refuse_before_publication() {
    let fixture = Fixture::new(1e-5);
    let mut layers = layers(fixture.source());
    let old = [0.0; 2];
    let bad_h = [h(321.0); 2];
    for duration in [0.0, -1.0, f64::NAN] {
        assert!(
            ThermoplasticSectionStep::new(
                &layers,
                &old,
                &old,
                duration,
                ThicknessQuadrature::Gauss3,
                reduction(),
                budget()
            )
            .is_err()
        );
    }
    assert!(
        ThermoplasticSectionStep::new(
            &layers,
            &old[..1],
            &old,
            0.1,
            ThicknessQuadrature::Gauss3,
            reduction(),
            budget()
        )
        .is_err()
    );
    assert!(
        ThermoplasticSectionStep::new(
            &layers,
            &old,
            &bad_h,
            0.1,
            ThicknessQuadrature::Gauss3,
            reduction(),
            budget()
        )
        .is_err()
    );
    assert!(
        ThermoplasticSectionStep::new(
            &layers,
            &bad_h,
            &old,
            0.1,
            ThicknessQuadrature::Gauss3,
            reduction(),
            budget()
        )
        .is_err()
    );
    let tight = SectionBudget {
        maximum_layers: 2,
        maximum_quadrature_points: 11,
    };
    assert!(matches!(
        ThermoplasticSectionStep::new(
            &layers,
            &old,
            &old,
            0.1,
            ThicknessQuadrature::Gauss3,
            reduction(),
            tight
        ),
        Err(ThermoplasticSectionError::Section(_))
    ));
    layers[1].bottom_z_m = 1e-4;
    assert!(
        ThermoplasticSectionStep::new(
            &layers,
            &old,
            &old,
            0.1,
            ThicknessQuadrature::Gauss3,
            reduction(),
            budget()
        )
        .is_err()
    );
}
