use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyValue, Provenance,
    UncertaintyModel,
};

const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);

fn point(t: f64, phase: Option<f64>, pressure: f64) -> QueryPoint {
    let p = QueryPoint::new()
        .with_quantity("temperature", absolute_temperature(), t)
        .unwrap()
        .with_quantity("pressure", QuantitySpec::dimensional(PRESSURE), pressure)
        .unwrap();
    if let Some(phase) = phase {
        p.with_quantity(
            "source_phase_liquid",
            QuantitySpec::dimensional(Dims::NONE),
            phase,
        )
        .unwrap()
    } else {
        p
    }
}

fn card(latent: Option<f64>, liquid_density: bool, discrete: bool) -> MaterialCard {
    card_with_heat_capacities(
        latent,
        liquid_density,
        discrete,
        [1000.0, 2000.0],
        [2000.0, 4000.0],
    )
}

fn card_with_heat_capacities(
    latent: Option<f64>,
    liquid_density: bool,
    discrete: bool,
    solid_cp: [f64; 2],
    liquid_cp: [f64; 2],
) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let pressure = || {
        ValidityDomain::unconstrained().with_quantity(
            "pressure",
            QuantitySpec::dimensional(PRESSURE),
            100_000.0,
            200_000.0,
        )
    };
    for (phase, start, end, cp, rho) in [
        (0.0, 300.0, 400.0, solid_cp, 1000.0),
        (1.0, 400.0, 500.0, liquid_cp, 800.0),
    ] {
        for (name, dims, values) in [
            ("specific-heat-capacity", CP, cp),
            ("density", RHO, [rho, rho]),
        ] {
            if phase == 1.0 && name == "density" && !liquid_density {
                continue;
            }
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: PropertyValue::Curve {
                        abscissa: "temperature".into(),
                        abscissa_dims: fs_qty::Temperature::DIMS,
                        knots: vec![(start, values[0]), (end, values[1])],
                        dims,
                    },
                    validity: pressure()
                        .with_quantity("temperature", absolute_temperature(), start, end)
                        .with_quantity(
                            "source_phase_liquid",
                            QuantitySpec::dimensional(Dims::NONE),
                            phase,
                            phase,
                        ),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: if discrete {
                        InterpolationPolicy::TabulatedOnly
                    } else {
                        InterpolationPolicy::LinearInside
                    },
                    observations: vec![],
                    provenance: Provenance {
                        source: "synthetic melting source; not measured material data".into(),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
    }
    for (name, quantity, value) in [
        ("melting-point", absolute_temperature(), Some(400.0)),
        (
            "latent-heat-fusion",
            QuantitySpec::dimensional(SPECIFIC_ENERGY),
            latent,
        ),
    ] {
        if let Some(value) = value {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::with_quantity(name, quantity),
                    value: PropertyValue::Scalar {
                        value,
                        dims: quantity.dims(),
                    },
                    validity: pressure(),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    observations: vec![],
                    provenance: Provenance {
                        source: "synthetic equilibrium transition".into(),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic pure substance".into(),
            phase: "solid/liquid".into(),
            process: "analytic numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

struct Query {
    cp: PropertyKey,
    density: PropertyKey,
    tm: PropertyKey,
    latent: PropertyKey,
    points: [QueryPoint; 5],
}
impl Query {
    fn new() -> Self {
        Self {
            cp: PropertyKey::new("specific-heat-capacity", CP),
            density: PropertyKey::new("density", RHO),
            tm: PropertyKey::with_quantity("melting-point", absolute_temperature()),
            latent: PropertyKey::new("latent-heat-fusion", SPECIFIC_ENERGY),
            points: [
                point(300.0, Some(0.0), 100_000.0),
                point(400.0, Some(0.0), 100_000.0),
                point(400.0, Some(1.0), 100_000.0),
                point(500.0, Some(1.0), 100_000.0),
                point(400.0, None, 100_000.0),
            ],
        }
    }
    fn query(&self) -> MeltingSourceQuery<'_> {
        let branch = |index| HeatCapacitySourceQuery {
            heat_capacity: &self.cp,
            density: &self.density,
            temperature_axis: "temperature",
            lower: &self.points[index],
            upper: &self.points[index + 1],
            selection: SelectionPolicy::SingleClaimOnly,
        };
        MeltingSourceQuery {
            solid: branch(0),
            liquid: branch(2),
            melting_point: &self.tm,
            latent_heat_fusion: &self.latent,
            transition: &self.points[4],
            transition_selection: SelectionPolicy::SingleClaimOnly,
            phase_selector: Some("source_phase_liquid"),
        }
    }
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-9 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g1_source_cp_integrals_and_one_latent_jump_drive_original_phase_state() {
    let source = card(Some(200_000.0), true, false);
    let query = Query::new();
    let model = ResolvedMeltingCurve::try_new(&source, query.query(), 0.0, 0.1).unwrap();
    near(model.melting_temperature_k(), 400.0);
    near(model.latent_heat_fusion_j_kg(), 200_000.0);
    let [solid_h, liquid_h] = model.plateau_enthalpy_bounds_j_kg();
    near(solid_h, 150_000.0);
    near(liquid_h, 350_000.0);
    let curve = model.curve();
    near(
        curve.knots().last().unwrap().specific_enthalpy_j_kg,
        650_000.0,
    );
    let initial = curve.state_at_specific_enthalpy(0.0).unwrap();
    let halfway = curve.advance_specific_energy(initial, 250_000.0).unwrap();
    near(halfway.temperature_k(), 400.0);
    near(halfway.liquid_mass_fraction(), 0.5);
    // Additive specific volume at fixed mass; density is not averaged directly.
    near(
        halfway.bulk_density_kg_m3(),
        1.0 / (0.5 / 1000.0 + 0.5 / 800.0),
    );
    let mass = 2.0;
    near(
        mass / halfway.bulk_density_kg_m3(),
        0.5 * mass / 1000.0 + 0.5 * mass / 800.0,
    );
    assert_eq!(
        curve
            .temperature_derivative_at_specific_enthalpy(250_000.0)
            .unwrap(),
        0.0
    );
    near(
        curve
            .liquid_mass_fraction_derivative_at_specific_enthalpy(250_000.0)
            .unwrap(),
        1.0 / 200_000.0,
    );
    let hot = curve.advance_specific_energy(halfway, 400_000.0).unwrap();
    near(hot.temperature_k(), 500.0);
    near(hot.liquid_mass_fraction(), 1.0);
    let cooled = curve.advance_specific_energy(hot, -650_000.0).unwrap();
    assert_eq!(cooled, initial);
    assert_eq!(curve.material_card_identity(), source.content_hash());
    for receipt in model.receipts() {
        source.claims().verify_receipt(receipt).unwrap();
    }
}

#[test]
fn g3_reference_shift_and_refinement_preserve_energy_and_phase() {
    let source = card(Some(200_000.0), true, false);
    let query = Query::new();
    let original = ResolvedMeltingCurve::try_new(&source, query.query(), 0.0, 0.1).unwrap();
    let shifted =
        ResolvedMeltingCurve::try_new(&source, query.query(), -1_000_000.0, 0.01).unwrap();
    for h in [
        0.0, 62_500.0, 150_000.0, 250_000.0, 350_000.0, 475_000.0, 650_000.0,
    ] {
        let a = original.curve().state_at_specific_enthalpy(h).unwrap();
        let b = shifted
            .curve()
            .state_at_specific_enthalpy(h - 1_000_000.0)
            .unwrap();
        assert!((a.temperature_k() - b.temperature_k()).abs() < 0.0002);
        near(a.liquid_mass_fraction(), b.liquid_mass_fraction());
    }
    assert_ne!(original.curve().identity(), shifted.curve().identity());
    assert!(original.curve().state_at_specific_enthalpy(-1.0).is_err());
    assert!(
        original
            .curve()
            .state_at_specific_enthalpy(650_001.0)
            .is_err()
    );
}

#[test]
fn g0_missing_density_fusion_and_discrete_capacity_refuse_without_defaults() {
    let q = Query::new();
    for source in [
        card(None, true, false),
        card(Some(200_000.0), false, false),
        card(Some(200_000.0), true, true),
        card(Some(0.0), true, false),
        card(Some(-1.0), true, false),
    ] {
        assert!(ResolvedMeltingCurve::try_new(&source, q.query(), 0.0, 0.1).is_err());
    }
    let source = card(Some(200_000.0), true, false);
    assert!(ResolvedMeltingCurve::try_new(&source, q.query(), 1e300, 0.1).is_err());
    let mut wrong_units = Query::new();
    wrong_units.tm = PropertyKey::new("melting-point", fs_qty::Temperature::DIMS);
    assert!(ResolvedMeltingCurve::try_new(&source, wrong_units.query(), 0.0, 0.1).is_err());
}

#[test]
fn g0_finite_plateau_endpoints_cannot_hide_unrepresentable_latent_width() {
    let source =
        card_with_heat_capacities(Some(f64::MAX), true, false, [1e305, 1e305], [1e305, 1e305]);
    let query = Query::new();
    let solid_end_h = -4e307 + 1e305 * 100.0;
    let liquid_start_h = solid_end_h + f64::MAX;
    assert!(solid_end_h.is_finite() && liquid_start_h.is_finite());
    assert!(!(liquid_start_h - solid_end_h).is_finite());
    assert!(matches!(
        ResolvedMeltingCurve::try_new(&source, query.query(), -4e307, 0.1),
        Err(PhaseStateError::InvalidCurve {
            what: "fusion energy is not representable at the selected enthalpy datum"
        })
    ));
}

#[test]
fn g0_transition_gaps_pressure_changes_and_phase_relabeling_refuse() {
    let source = card(Some(200_000.0), true, false);
    for (index, replacement) in [
        (1, point(399.0, Some(0.0), 100_000.0)),
        (2, point(401.0, Some(1.0), 100_000.0)),
        (4, point(401.0, None, 100_000.0)),
        (4, point(400.0, None, 200_000.0)),
        (2, point(400.0, Some(0.0), 100_000.0)),
        (4, point(400.0, Some(0.0), 100_000.0)),
    ] {
        let mut query = Query::new();
        query.points[index] = replacement;
        assert!(ResolvedMeltingCurve::try_new(&source, query.query(), 0.0, 0.1).is_err());
    }
    let mut query = Query::new();
    for index in [2, 3] {
        query.points[index] = query.points[index]
            .clone()
            .with("pressure", 200_000.0)
            .unwrap();
    }
    // Every separate source query is admissible; their inconsistent pressure is not.
    assert!(ResolvedMeltingCurve::try_new(&source, query.query(), 0.0, 0.1).is_err());
}
