//! Evolving source viscosity reaches the original physical pour trajectory.

use super::*;
use fs_lbm::material::free_surface::rheology::{IsothermalRheologyLaw, RheologySolve};
use fs_material::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_vessel::MaterialRheologyPour;

fn solve() -> RheologySolve {
    RheologySolve {
        maximum_iterations: 128,
        relaxation: 0.5,
        absolute_relaxation_tolerance: 1e-13,
        relative_relaxation_tolerance: 1e-13,
    }
}

fn law<'a>(
    card: &'a MaterialCard,
    state: fs_material::state_point::ResolvedMaterialStatePoint,
) -> IsothermalRheologyLaw<'a> {
    IsothermalRheologyLaw::try_new(
        card, state, LiquidViscosityBasis::Dynamic, LiquidCapillarity::Sourced, "T", "rate",
    ).unwrap()
}

#[test]
fn g1_constant_source_reproduces_the_frozen_si_pour_and_its_render_field() {
    let (card, liquid) = source(LiquidViscosityBasis::Dynamic, Some(0.002));
    let query = point(300.0)
        .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 0.0)
        .unwrap();
    let binding = law(&card, liquid.material().at_query_point(&card, &query).unwrap());
    let evolving = MaterialRheologyPour::try_new(&binding, spec(), solve()).unwrap();
    let frozen = MaterialPour::try_new(&liquid, spec()).unwrap();
    assert_eq!(
        evolving.scaling().surface_tension_lattice(),
        frozen.scaling().surface_tension_lattice(),
    );
    // The existing frozen entry still refuses an evolving source marker.
    assert!(MaterialPour::try_new(evolving.scaling().material(), spec()).is_err());
    let gate = CancelGate::new_clock_free();
    for contact in [ContactModel::Neutral, ContactModel::Wetting] {
        let actual = evolving.run(contact, &gate).unwrap();
        let direct = frozen.run(contact, &gate).unwrap();
        assert_eq!(actual.pour().native().mass_field, direct.native().mass_field);
        assert_eq!(actual.pour().native().mass_drift, direct.native().mass_drift);
        assert_eq!(actual.pour().native().poured_mass, direct.native().poured_mass);
        assert_eq!(actual.pour().native().dribble_cells, direct.native().dribble_cells);
        assert_eq!(actual.pour().report(), direct.report());
        let receipt = actual.last_collision().material_at([2, 2]).unwrap();
        assert_eq!(receipt.card_identity(), card.content_hash());
        assert!(receipt.query_point().contains(&("surface_environment".into(), 7.0)));
        near(receipt.property(DYNAMIC_VISCOSITY_PROPERTY).unwrap().value_si(), 0.02);
    }
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        evolving.run(ContactModel::Neutral, &cancelled),
        Err(MaterialPourError::Cancelled),
    ));
}

fn rate_card(slope: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, initial) in [
        (DENSITY_PROPERTY, Density::DIMS, 1000.0),
        (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, 0.02),
        (SURFACE_TENSION_PROPERTY, SURFACE_TENSION_DIMS, 0.002),
    ] {
        let rate = name == DYNAMIC_VISCOSITY_PROPERTY;
        let mut validity = ValidityDomain::unconstrained()
            .with_quantity("T", *point(300.0).axis_quantities().get("T").unwrap(), 300.0, 400.0)
            .with("surface_environment", 7.0, 7.0);
        let value = if rate {
            validity = validity.with_quantity(
                "rate", QuantitySpec::dimensional(Frequency::DIMS), 0.0, 1000.0,
            );
            PropertyValue::Curve {
                abscissa: "rate".into(),
                abscissa_dims: Frequency::DIMS,
                knots: vec![(0.0, initial), (1000.0, initial + slope * 1000.0)],
                dims,
            }
        } else {
            PropertyValue::Scalar { value: initial, dims }
        };
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value,
            validity,
            uncertainty: UncertaintyModel::Unstated,
            interpolation: if rate {
                InterpolationPolicy::LinearInside
            } else {
                InterpolationPolicy::ConstantWithinValidity
            },
            observations: vec![],
            provenance: Provenance {
                source: "synthetic shear-response pour, not measurements".into(),
                license: "test-only".into(),
                artifact: None,
            },
        }).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic shear-dependent pour liquid".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    ).unwrap()
}

#[test]
fn g3_original_rate_curve_changes_the_native_pour_without_changing_its_physical_setup() {
    let sources = [rate_card(0.0), rate_card(0.01)];
    let query = point(300.0)
        .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 0.0)
        .unwrap();
    let requirements = [
        (DENSITY_PROPERTY, Density::DIMS),
        (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS),
        (SURFACE_TENSION_PROPERTY, SURFACE_TENSION_DIMS),
    ].map(|(name, dims)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive).unwrap()
    });
    let bindings: Vec<_> = sources.iter().map(|card| {
        law(card, resolve_material_state_point(
            card, &query, &requirements, MaterialPropertySelection::SingleClaimOnly,
        ).unwrap())
    }).collect();
    let gate = CancelGate::new_clock_free();
    let runs: Vec<_> = bindings.iter().map(|binding| {
        MaterialRheologyPour::try_new(binding, spec(), solve()).unwrap()
            .run(ContactModel::Neutral, &gate).unwrap()
    }).collect();
    assert_ne!(runs[0].pour().native().mass_field, runs[1].pour().native().mass_field);
    assert!((runs[0].pour().report().peak_speed_m_s - runs[1].pour().report().peak_speed_m_s).abs() > 1e-10);
    for (run, source) in runs.iter().zip(&sources) {
        let report = run.pour().report();
        assert!((report.final_mass_kg_per_m / report.initial_mass_kg_per_m - 1.0).abs() < 1e-10);
        let receipt = run.last_collision().material_at([2, 2]).unwrap();
        assert_eq!(receipt.card_identity(), source.content_hash());
        assert!(receipt.query_point().contains(&("T".into(), 300.0)));
    }
    let max_mu = (0..spec().cells[1]).flat_map(|y| (0..spec().cells[0]).map(move |x| [x, y]))
        .filter_map(|xy| runs[1].last_collision().material_at(xy))
        .map(|material| material.property(DYNAMIC_VISCOSITY_PROPERTY).unwrap().value_si())
        .fold(0.02_f64, f64::max);
    assert!(max_mu > 0.020001, "actual collision must query the evolving rate source");
}
