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

fn hydrostatic_law(card: &MaterialCard) -> IsothermalRheologyLaw<'_> {
    let query = point(300.0)
        .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 0.0)
        .unwrap();
    let reference = resolve_liquid_surface_state(
        card,
        &query,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    IsothermalRheologyLaw::try_new(
        card, reference.material().clone(), LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected, "T", "rate",
    ).unwrap()
}

#[test]
fn g1_hydrostatic_source_rheology_retains_rest_and_the_newtonian_limit() {
    let cards = [rate_card(0.0), rate_card(0.01)];
    let bindings: Vec<_> = cards.iter().map(hydrostatic_law).collect();
    let input = MaterialPourSpec {
        timestep_s: 0.00025,
        gravity_m_s2: 9.80665,
        tilt_final_radians: 0.0,
        maximum_speed_m_s: 1e-8,
        steps: 16,
        ..spec()
    };
    let query = point(300.0)
        .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 0.0)
        .unwrap();
    let fixed = resolve_liquid_surface_state(
        &cards[0], &query, LiquidViscosityBasis::Dynamic, LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let gate = CancelGate::new_clock_free();
    let expected = MaterialPour::try_new(&fixed, input).unwrap()
        .run_hydrostatic(100_000.0, &gate).unwrap();
    for (binding, card) in bindings.iter().zip(&cards) {
        let actual = MaterialRheologyPour::try_new(binding, input, solve()).unwrap()
            .run_hydrostatic(100_000.0, &gate).unwrap();
        assert!(actual.pour().report().peak_speed_m_s < 2e-11);
        near(actual.pour().report().initial_mass_kg_per_m, expected.report().initial_mass_kg_per_m);
        near(actual.pour().report().final_mass_kg_per_m, expected.report().final_mass_kg_per_m);
        for (&a, &b) in actual.pour().native().mass_field.iter().zip(&expected.native().mass_field) {
            if b == 0.0 {
                assert!(a.abs() < 2e-12);
            } else {
                near(a, b);
            }
        }
        let receipt = actual.last_collision().material_at([2, 2]).unwrap();
        let sample = actual.last_collision().rheology_at([2, 2]).unwrap();
        assert_eq!(receipt.card_identity(), card.content_hash());
        assert!(sample.shear_rate_per_s < 1e-7);
        assert!(receipt.query_point().contains(&("surface_environment".into(), 7.0)));
        assert_eq!(actual.pour().scaling().reference_pressure_pa(), Some(100_000.0));
    }
}

#[test]
fn g3_hydrostatic_pour_keeps_evolving_source_receipts_and_gauge_invariance() {
    let cards = [rate_card(0.0), rate_card(0.01)];
    let bindings: Vec<_> = cards.iter().map(hydrostatic_law).collect();
    let gate = CancelGate::new_clock_free();
    // Prescribed frame motion excites measurable shear after the initial
    // pressure head is balanced.
    let input = MaterialPourSpec { tilt_final_radians: 0.02, ..spec() };
    let runs: Vec<_> = bindings.iter().map(|binding| {
        MaterialRheologyPour::try_new(binding, input, solve()).unwrap()
            .run_hydrostatic(100_000.0, &gate).unwrap()
    }).collect();
    assert_ne!(runs[0].pour().native().mass_field, runs[1].pour().native().mass_field);
    let shifted = MaterialRheologyPour::try_new(&bindings[1], input, solve()).unwrap()
        .run_hydrostatic(104_096.0, &gate).unwrap();
    assert_eq!(shifted.pour().native().mass_field, runs[1].pour().native().mass_field);
    assert_eq!(shifted.pour().report(), runs[1].pour().report());
    assert_eq!(shifted.pour().scaling().reference_pressure_pa(), Some(104_096.0));
    let mut maximum_viscosity = 0.02_f64;
    for y in 0..spec().cells[1] {
        for x in 0..spec().cells[0] {
            let Some(receipt) = shifted.last_collision().material_at([x, y]) else { continue };
            let sample = shifted.last_collision().rheology_at([x, y]).unwrap();
            let viscosity = receipt.property(DYNAMIC_VISCOSITY_PROPERTY).unwrap().value_si();
            near(viscosity, 0.02 + 0.01 * sample.shear_rate_per_s);
            maximum_viscosity = maximum_viscosity.max(viscosity);
            assert_eq!(receipt.card_identity(), cards[1].content_hash());
            // The numerical pressure datum does not requery the fixed-pressure law.
            assert!(receipt.query_point().contains(&("pressure".into(), 100_000.0)));
        }
    }
    assert!(maximum_viscosity > 0.020001, "hydrostatic startup must retain actual rate queries");
    near(shifted.pour().report().final_mass_kg_per_m, shifted.pour().report().initial_mass_kg_per_m);
}

#[test]
fn g0_g4_hydrostatic_source_pour_retains_shared_refusals_and_retry() {
    let (card, liquid) = source(LiquidViscosityBasis::Dynamic, Some(0.002));
    let query = point(300.0)
        .with_quantity("rate", QuantitySpec::dimensional(Frequency::DIMS), 0.0)
        .unwrap();
    let capillary = law(&card, liquid.material().at_query_point(&card, &query).unwrap());
    let gate = CancelGate::new_clock_free();
    assert!(matches!(
        MaterialRheologyPour::try_new(&capillary, spec(), solve()).unwrap()
            .run_hydrostatic(100_000.0, &gate),
        Err(MaterialPourError::Invalid(_)),
    ));
    let binding = hydrostatic_law(&card);
    let prepared = MaterialRheologyPour::try_new(&binding, spec(), solve()).unwrap();
    for pressure in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(prepared.run_hydrostatic(pressure, &gate).is_err());
    }
    let unresolved = MaterialRheologyPour::try_new(
        &binding, MaterialPourSpec { maximum_density_excursion: 1e-10, ..spec() }, solve(),
    ).unwrap();
    assert!(unresolved.run_hydrostatic(100_000.0, &gate).is_err());
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        prepared.run_hydrostatic(100_000.0, &cancelled),
        Err(MaterialPourError::Cancelled),
    ));
    let a = prepared.run_hydrostatic(100_000.0, &gate).unwrap();
    let b = prepared.run_hydrostatic(100_000.0, &gate).unwrap();
    assert_eq!(a.pour().native().mass_field, b.pour().native().mass_field);
    assert_eq!(a.pour().report(), b.pour().report());
    assert_eq!(prepared.scaling().reference_pressure_pa(), None);
    assert!(core::ptr::eq(prepared.law(), &binding));
}
