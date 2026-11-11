//! Original-source heat and analytic SGS flux regressions.

use super::*;
use crate::conjugate::domain::Voxel;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    UncertaintyModel,
};
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, HeatCapacityKnot, SolidLiquidPhase,
};
use fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS as KD;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::{Dims, QuantitySpec, Temperature};

fn card(insulating: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new("k", KD),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: if insulating { vec![(300.0, 0.0), (500.0, 0.0)] }
                else { vec![(300.0, 1.0), (500.0, 3.0)] },
            dims: KD,
        },
        validity: ValidityDomain::unconstrained(),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured original SGS heat/Jacobian source".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }).unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured thermal source".into(),
            phase: "original liquid chart".into(),
            process: "SGS transport regression".into(),
            revision: 0,
        },
        claims,
        vec![],
    ).unwrap()
}

fn curve(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity_exact(
        card.content_hash(), SolidLiquidPhase::Liquid, 0.0,
        &[HeatCapacityKnot {
            temperature_k: 300.0, specific_heat_capacity_j_kg_k: 2.0, bulk_density_kg_m3: 900.0,
        }, HeatCapacityKnot {
            temperature_k: 400.0, specific_heat_capacity_j_kg_k: 4.0, bulk_density_kg_m3: 900.0,
        }, HeatCapacityKnot {
            temperature_k: 500.0, specific_heat_capacity_j_kg_k: 1.0, bulk_density_kg_m3: 900.0,
        }],
    ).unwrap()
}

fn prepare_source<'a>(
    card: &MaterialCard,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    bounds: [f64; 2],
) -> PreparedPhaseConductivity<'a> {
    let context = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0).unwrap()
        .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0).unwrap();
    PreparedPhaseConductivity::new(
        curve,
        ThermalMaterialPoint {
            material: card, context: &context, temperature_axis: "T", liquid_fraction_axis: Some("f"),
        },
        &PropertyKey::new("k", KD),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        bounds, || Ok(()),
    ).unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= 3e-12 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}");
}

fn geometry() -> ConductiveLink {
    ConductiveLink {
        link_id: 7, half_distances_m: [0.2, 0.7], contact_resistance_m2_k_w: 0.0,
    }
}

fn patch() -> ContactPatch {
    ContactPatch::symmetric(ThermalLink {
        id: 7, a: CellId(0), b: CellId(1), area_m2: 1.7, conductance_w_k: 0.0,
    })
}

/// Independent inverse of integral Cp(T)dT on both original sensible spans.
fn caloric(h: f64) -> (f64, f64) {
    if h <= 300.0 {
        let x = 2.0 * h / (2.0 + (4.0 + 0.04 * h).sqrt());
        (300.0 + x, 2.0 + 0.02 * x)
    } else {
        let dh = h - 300.0;
        let x = 2.0 * dh / (4.0 + (16.0 - 0.06 * dh).sqrt());
        (400.0 + x, 4.0 - 0.03 * x)
    }
}

#[test]
fn original_k_and_exact_cp_drive_total_resistance_and_complete_newton_tangents() {
    let material = card(false);
    let curve = curve(&material);
    let source = prepare_source(&material, &curve, [0.0, 550.0]);
    let law = SgsThermalLink { source: &source, density_diffusivity: [0.1, 0.5] };
    let h = [125.0, 462.5];
    let response = law.state_at_specific_enthalpies(h, geometry(), patch(), &mut || Ok(())).unwrap();
    // T=[350,450], Cp=[3,2.5], original k=[1.5,2.5]. The total
    // conductivities enter their common series resistance before forming a flux.
    let expected_g = 1.7 / (0.2 / 1.8 + 0.7 / 3.75);
    near(response.conductance_w_k, expected_g);
    let independent_flux = |h: [f64; 2]| {
        let state = h.map(caloric);
        let k = [
            1.0 + 0.01 * (state[0].0 - 300.0) + 0.1 * state[0].1,
            1.0 + 0.01 * (state[1].0 - 300.0) + 0.5 * state[1].1,
        ];
        1.7 * (state[0].0 - state[1].0) / (0.2 / k[0] + 0.7 / k[1])
    };
    for side in 0..2 {
        let mut lower = h;
        let mut upper = h;
        lower[side] -= 1e-4;
        upper[side] += 1e-4;
        let derivative = (independent_flux(upper) - independent_flux(lower)) / 2e-4;
        assert!((response.power_enthalpy_derivatives_w_kg_j[side] - derivative).abs() < 3e-7);
    }
    let gate = CancelGate::new_clock_free();
    near(law.conductance_w_k(h, geometry(), 1.7, &gate).unwrap(), expected_g);
    let samples = h.map(|h| source.state_at_specific_enthalpy(h).unwrap());
    near(law.conductance_from_samples_w_k(samples, geometry(), 1.7, &gate).unwrap(), expected_g);
    // A compiled source ending at the Cp kink selects its LEFT tangent.
    // A wider original source selects the RIGHT tangent at that same h.
    let truncated = prepare_source(&material, &curve, [0.0, 300.0]);
    let left = sensible_capacity(&truncated, truncated.state_at_specific_enthalpy(300.0).unwrap()).unwrap();
    let right = sensible_capacity(&source, source.state_at_specific_enthalpy(300.0).unwrap()).unwrap();
    near(left[0], 4.0);
    near(right[0], 4.0);
    near(left[1], 0.005);
    near(right[1], -0.0075);
}

#[test]
fn zero_molecular_k_can_mix_interior_fluid_but_solid_faces_remain_molecular() {
    let material = card(true);
    let curve = curve(&material);
    let source = prepare_source(&material, &curve, [0.0, 550.0]);
    let gate = CancelGate::new_clock_free();
    let mut domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    domain.set(2, 0, 0, Voxel::Solid(0));
    let mass = [2.0, 5.0, 0.0];
    let field = SgsThermalField {
        eddy_viscosity_m2_s: &[0.04, 0.08, 0.0], turbulent_prandtl_number: 0.8,
    };
    let mut events = Vec::new();
    for cell in 0..3 {
        events.push(TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
            id: CellId(cell), material_id: MaterialId(0), body_id: 0, segment_id: cell,
            mass_kg: 2.0, specific_h_j_kg: 125.0,
            geometry: CellGeometry {
                start_m: [cell as f64, 0.0, 0.0], end_m: [cell as f64 + 1.0, 0.0, 0.0],
                volume_m3: 1.0, surface_area_m2: 6.0,
            },
            boundary: Boundary { ambient_k: 350.0, htc_w_m2_k: 0.0, reservoir: None, source_w: 0.0 },
        })});
    }
    for id in 0..2 {
        events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
            id, a: CellId(id), b: CellId(id + 1), area_m2: 1.0, conductance_w_k: 0.0,
        })});
    }
    let caloric = [(MaterialId(0), &curve as &dyn ThermalLaw)];
    let initial = DepositionSimulation::new(events, &caloric, DepositionBudget {
        max_cells: 3, max_links: 2, max_events: 5, max_substeps: 1,
    }).unwrap();
    let (initial, _) = initial.advance_to(0.0, StepConfig { max_dt_s: 1.0, energy_tolerance_j: 0.0 }, || Ok(())).unwrap();
    let links = prepare_links(&domain, &source, &initial, &mass, field, 2, &gate).unwrap();
    assert_eq!(links.len(), 1);
    assert_eq!(links[0].0, 0); // fluid--solid link 1 has no SGS override.
    let g = ConductiveLink { link_id: 0, half_distances_m: [0.5; 2], contact_resistance_m2_k_w: 0.0 };
    let conductance = links[0].1.conductance_w_k([125.0, 462.5], g, 1.0, &gate).unwrap();
    near(conductance, 1.0 / (0.5 / 0.3 + 0.5 / 1.25));
    assert!(prepare_links(&domain, &source, &initial, &mass,
        SgsThermalField { eddy_viscosity_m2_s: &[0.0; 3], ..field }, 2, &gate).unwrap().is_empty());
    for pr in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(admit_prandtl_number(pr).is_err());
    }
    assert!(prepare_links(&domain, &source, &initial, &mass,
        SgsThermalField { eddy_viscosity_m2_s: &[0.04, 0.08, 0.01], ..field }, 2, &gate).is_err());
    assert!(prepare_links(&domain, &source, &initial, &mass, field, 1, &gate).is_err());
}

#[test]
fn conservative_flux_feedback_counts_opposing_links_and_refuses_bad_fields() {
    let material = card(false);
    let curve = curve(&material);
    let source = prepare_source(&material, &curve, [0.0, 550.0]);
    let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    let mass = [2.0; 3];
    let used = SgsThermalField {
        eddy_viscosity_m2_s: &[0.01; 3], turbulent_prandtl_number: 0.5,
    };
    let actual = [0.02; 3];
    let samples = [0.0, 300.0, 0.0].map(|h| Some(source.state_at_specific_enthalpy(h).unwrap()));
    let gate = CancelGate::new_clock_free();
    let run = |field, endpoint: &[Option<PhaseConductivitySample>], max_faces, gate| {
        heat_flux_mismatch_j(&domain, &source, &mass, used, field, endpoint, 0.2, max_faces, gate)
    };
    // Both physical face changes are equal/opposite around the warm middle
    // cell. A signed domain sum would hide them; 2|deltaQ| per face cannot.
    let expected: f64 = 4.0 * 0.2 * ((4.0 / 3.0) * (1.16 - 1.08)) * 100.0;
    let actual_bound = run(&actual, &samples, 2, &gate).unwrap();
    // Outward rounding bounds the native rounded conductances; the closed
    // form uses a different operation order, admitted within a few ulps here.
    assert!(actual_bound + 32.0 * f64::EPSILON * expected.abs() >= expected);
    near(actual_bound, expected);
    assert_eq!(run(used.eddy_viscosity_m2_s, &samples, 2, &gate).unwrap(), 0.0);
    let isothermal = [Some(source.state_at_specific_enthalpy(125.0).unwrap()); 3];
    assert_eq!(run(&actual, &isothermal, 2, &gate).unwrap(), 0.0);
    assert!(run(&[], &samples, 2, &gate).is_err());
    assert!(run(&[0.02, -0.01, 0.02], &samples, 2, &gate).is_err());
    assert!(run(&actual, &samples, 1, &gate).is_err());
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(run(&actual, &samples, 2, &cancelled), Err(ChtError::Cancelled)));
    assert_eq!(run(&actual, &samples, 2, &gate).unwrap(), actual_bound);

    // Latent capacity is not an inferred turbulent diffusivity. The exact
    // zero-scale branch can retain the original molecular plateau response.
    let latent_curve = EquilibriumEnthalpyPhaseCurve::try_new(material.content_hash(), vec![
        EnthalpyPhaseKnot { specific_enthalpy_j_kg: 0.0, temperature_k: 300.0, liquid_mass_fraction: 0.0, bulk_density_kg_m3: 900.0 },
        EnthalpyPhaseKnot { specific_enthalpy_j_kg: 100.0, temperature_k: 300.0, liquid_mass_fraction: 1.0, bulk_density_kg_m3: 900.0 },
        EnthalpyPhaseKnot { specific_enthalpy_j_kg: 200.0, temperature_k: 400.0, liquid_mass_fraction: 1.0, bulk_density_kg_m3: 900.0 },
    ]).unwrap();
    let latent = prepare_source(&material, &latent_curve, [0.0, 100.0]);
    let sample = latent.state_at_specific_enthalpy(100.0).unwrap();
    let zero = SgsThermalLink { source: &latent, density_diffusivity: [0.0; 2] };
    assert!(zero.response([sample; 2], geometry(), 1.7).is_ok());
    let active = SgsThermalLink { source: &latent, density_diffusivity: [0.1; 2] };
    assert!(active.response([sample; 2], geometry(), 1.7).is_err());
}
