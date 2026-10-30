//! G1/G3/G4: source-only saturation and temperature changes reach the actual
//! retained 3D field, native Joule heat and finite caloric cells together.
#[path = "nonlinear/heat.rs"]
mod heat;

use super::*;
use fs_em::spatial::{
    nonlinear::SpatialNonlinearMaterial,
    vector_potential::{
        induction::nonlinear::{SpatialNonlinearInduction, SpatialNonlinearInductionControls},
        nonlinear::{NonlinearVectorControls, SpatialNonlinearVectorMagnetostatics},
    },
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, PropertyClaim, PropertyValue, Provenance, QueryPoint,
    UncertaintyModel,
};
use fs_material::{
    field_storage::{
        MAGNETIC_PERMEABILITY_PROPERTY,
        nonlinear::{REVERSIBLE_FIELD_STRENGTH_PROPERTY, ReversibleMagneticMaterial},
        permanent::{FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS},
    },
    state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    },
};
use fs_qty::QuantitySpec;

fn context(t: f64) -> QueryPoint {
    source::absolute_point(t)
        .with_quantity("B", QuantitySpec::dimensional(FLUX_DENSITY_DIMS), 0.0)
        .unwrap()
}
fn h(b: f64, saturating: bool) -> f64 {
    if saturating && b > 1.0 {
        1.0 + 4.0 * (b - 1.0)
    } else {
        b
    }
}
fn energy(a: f64, saturating: bool) -> f64 {
    let b = 2.0 * a;
    if saturating && b > 1.0 {
        0.5 + (b - 1.0) + 2.0 * (b - 1.0).powi(2)
    } else {
        0.5 * b * b
    }
}
fn card(rho: f64, alpha: f64, basis: ElectricalPropertyBasis, saturating: bool) -> MaterialCard {
    let base = source::absolute_affine_magnetic_card(rho, alpha, 1.0, 0.0, basis);
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        if claim.key.name() != MAGNETIC_PERMEABILITY_PROPERTY {
            claims.insert_claim(claim.clone()).unwrap();
        }
    }
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(REVERSIBLE_FIELD_STRENGTH_PROPERTY, FIELD_STRENGTH_DIMS),
            value: PropertyValue::Curve {
                abscissa: "B".into(),
                abscissa_dims: FLUX_DENSITY_DIMS,
                knots: [0.0, 1.0, 4.0].map(|b| (b, h(b, saturating))).to_vec(),
                dims: FIELD_STRENGTH_DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "B",
                QuantitySpec::dimensional(FLUX_DENSITY_DIMS),
                0.0,
                4.0,
            ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured 3D saturation and thermal feedback".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}
fn conductor(
    card: &MaterialCard,
    t: f64,
    basis: ElectricalPropertyBasis,
) -> ResolvedConductorState {
    let material = resolve_material_state_point(
        card,
        &context(t),
        &[
            basis.requirement().unwrap(),
            ScalarPropertyRequirement::try_new(
                REVERSIBLE_FIELD_STRENGTH_PROPERTY,
                FIELD_STRENGTH_DIMS,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
            ScalarPropertyRequirement::try_new(
                "density",
                Density::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
            ScalarPropertyRequirement::try_new(
                "specific_heat_capacity",
                SPECIFIC_HEAT_CAPACITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    ResolvedConductorState::try_from_material_with_basis(material, basis).unwrap()
}
fn problem(
    card: &MaterialCard,
    basis: ElectricalPropertyBasis,
) -> SpatialNonlinearVectorMagnetostatics {
    // Reuse the original mesh fixture; only the material law is changed.
    let geometry = super::problem(
        &source::absolute_affine_magnetic_card(0.2, 0.0, 1.0, 0.0, basis),
        basis,
    );
    let law = ReversibleMagneticMaterial::resolve(
        card,
        &context(300.0),
        "B",
        4.0,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    SpatialNonlinearVectorMagnetostatics::new(
        geometry.positions_m().to_vec(),
        geometry.mesh().tets.clone(),
        vec![SpatialNonlinearMaterial::Reversible(law); 6],
        assembly(),
        &mut || false,
    )
    .unwrap()
}
fn field_controls() -> SpatialNonlinearInductionControls {
    SpatialNonlinearInductionControls {
        induction: super::field_controls(),
        maximum_newton_iterations: 40,
        maximum_backtracks: 30,
    }
}
fn initial<'m>(
    problem: &'m SpatialNonlinearVectorMagnetostatics,
    card: &MaterialCard,
    basis: ElectricalPropertyBasis,
    a: f64,
    saturating: bool,
) -> SpatialNonlinearInduction<'m> {
    let boundary = VectorMagneticBoundary {
        potential: problem
            .exterior_faces()
            .into_iter()
            .map(|face| VectorPotentialFace {
                face,
                edge_integrals_wb: [0.0; 3],
            })
            .collect(),
    };
    let mut guess = vec![0.0; problem.mesh().edges.len()];
    guess[diagonal(problem)] = a;
    let solve = CheckedCgConfig {
        tolerance: 1e-12,
        max_iterations: 1000,
        max_corrections: 2,
    };
    let field = problem
        .solve(
            &[[4.0 * h(2.0 * a, saturating); 3]; 6],
            &boundary,
            &guess,
            NonlinearVectorControls {
                field: VectorMagneticControls {
                    field: solve,
                    gauge: solve,
                    maximum_gauge_product_terms: 100000,
                    current_tolerance_a: 1e-9,
                    flux_tolerance_wb: 1e-10,
                    energy_tolerance_j: 1e-9,
                },
                relative_residual_tolerance: 1e-12,
                maximum_newton_iterations: 40,
                maximum_backtracks: 30,
            },
            &mut || false,
        )
        .unwrap();
    let state = SpatialNonlinearInduction::new(
        &field,
        vec![Some(conductor(card, 300.0, basis)); 6],
        SpatialInductionLimits {
            maximum_unknowns: 1000,
            maximum_conducting_cells: 100,
            maximum_assembly_terms: 10000,
        },
        1e-10,
        &mut || false,
    )
    .unwrap();
    near(state.potential_wb()[diagonal(problem)], a, 2e-10);
    state
}
fn diagonal(problem: &SpatialNonlinearVectorMagnetostatics) -> usize {
    problem.mesh().edges.binary_search(&[0, 7]).unwrap()
}
fn cells<'a>(
    card: &'a MaterialCard,
    curve: &'a ResolvedHeatCapacityCurve,
    query: &'a QueryPoint,
) -> [MaterialCell<'a>; 1] {
    [MaterialCell {
        cell: CellId(7),
        material_id: MaterialId(9),
        curve: curve.curve(),
        point: ThermalMaterialPoint {
            material: card,
            context: query,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }]
}

// Exact unit-cube Whitney integration gives the charge-eliminated conductor
// mass M=3*sigma/40 and magnetic effort d psi(a)/da=2 H(2a). Nested scalar
// bisection is independent of production matrices, constitutive and Newton code.
fn reference(
    a0: f64,
    t0: f64,
    dt: f64,
    rho: f64,
    alpha: f64,
    basis: ElectricalPropertyBasis,
    saturating: bool,
) -> (f64, f64, f64, f64) {
    let amplitude = |t| {
        bracket(0.0, a0, |a| {
            3.0 * sigma(t, rho, alpha, basis) / 40.0 * (a - a0) / dt + 2.0 * h(2.0 * a, saturating)
        })
    };
    let heat = |t: f64| 3.0 * sigma(t, rho, alpha, basis) / 40.0 * (amplitude(t) - a0).powi(2) / dt;
    let t = bracket(t0, t0 + energy(a0, saturating) / 4.0, |t| {
        4.0 * (t - t0) - heat(t)
    });
    let a = amplitude(t);
    let numerical =
        2.0 * h(2.0 * a, saturating) * (a - a0) - energy(a, saturating) + energy(a0, saturating);
    (a, t, heat(t), numerical)
}

#[test]
fn source_saturation_and_temperature_close_the_same_field_and_caloric_history() {
    let (rho, alpha, dt, a0) = (0.2, 0.009, 0.05, 1.2);
    for basis in [
        ElectricalPropertyBasis::Resistivity,
        ElectricalPropertyBasis::Conductivity,
    ] {
        let mut temperatures = Vec::new();
        for saturating in [false, true] {
            let card = card(rho, alpha, basis, saturating);
            let curve = curve(&card);
            let thermal = thermal(&curve, 0.0);
            let problem = problem(&card, basis);
            let state = initial(&problem, &card, basis, a0, saturating);
            let query = context(300.0);
            let cells = cells(&card, &curve, &query);
            let regions = regions();
            let network = SpatialMagneticThermalNetwork {
                cells: &cells,
                links: &[],
                magnetic_regions: &regions,
                volume_tolerance_m3: 1e-12,
            };
            let candidate = network
                .step_nonlinear_candidate(
                    &thermal,
                    dt,
                    &state,
                    &[[0.0; 3]; 6],
                    &[[0.0; 3]; 8],
                    thermal_controls(),
                    field_controls(),
                    &[0.0],
                    || false,
                )
                .unwrap();
            let (a, t, q, d) = reference(a0, 300.0, dt, rho, alpha, basis, saturating);
            near(
                candidate.field().potential_wb()[diagonal(&problem)],
                a,
                2e-9,
            );
            near(
                candidate.thermal().temperature_k(CellId(7)).unwrap(),
                t,
                2e-8,
            );
            let report = candidate.report();
            let native = &report.coupled.electrical;
            near(native.heat_j, q, 2e-8);
            near(native.numerical_dissipation_j, d, 2e-8);
            near(report.coupled.deposited_joule_j, q, 2e-8);
            near(report.coupled.thermal.stored_energy_change_j, q, 2e-8);
            near(native.stored_after_j, energy(a, saturating), 2e-8);
            near(report.total_energy_residual_j, 0.0, 2e-8);
            near(report.external_heat_j, 0.0, 2e-8);
            assert!(d > 1e-4 && q > 1e-4);
            for source in candidate.field().conductors().iter().flatten() {
                near(
                    source.resistivity_ohm_m(),
                    sigma(t, rho, alpha, basis).recip(),
                    2e-12,
                );
                near(source.material().query().axes()["T"], t, 1e-12);
            }
            assert!(std::ptr::eq(candidate.field().problem(), state.problem()));
            assert_eq!(state.time_s(), 0.0);
            near(thermal.temperature_k(CellId(7)).unwrap(), 300.0, 1e-12);
            let continued = network
                .step_nonlinear_candidate(
                    candidate.thermal(),
                    2.0 * dt,
                    candidate.field(),
                    &[[0.0; 3]; 6],
                    &[[0.0; 3]; 8],
                    thermal_controls(),
                    field_controls(),
                    &[0.0],
                    || false,
                )
                .unwrap();
            let (a2, t2, _, _) = reference(a, t, dt, rho, alpha, basis, saturating);
            near(
                continued.field().potential_wb()[diagonal(&problem)],
                a2,
                3e-9,
            );
            near(
                continued.thermal().temperature_k(CellId(7)).unwrap(),
                t2,
                3e-8,
            );
            near(continued.report().total_energy_residual_j, 0.0, 2e-8);
            temperatures.push(t);
        }
        assert!(temperatures[1] > temperatures[0] + 0.01);
    }
}

#[test]
fn actual_motional_current_heats_finite_cells_and_exports_lorentz_work_separately() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let (rho, alpha, a0, dt, dilation, external_w) = (0.2, 0.009, 1.2, 0.05, 0.6, 40.0);
    let card = card(rho, alpha, basis, true);
    let curve = curve(&card);
    let thermal = thermal(&curve, external_w);
    let problem = problem(&card, basis);
    let state = initial(&problem, &card, basis, a0, true);
    let query = context(300.0);
    let cells = cells(&card, &curve, &query);
    let regions = regions();
    let network = SpatialMagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let velocity: Vec<_> = problem
        .positions_m()
        .iter()
        .map(|p| p.map(|x| dilation * x))
        .collect();
    let candidate = network
        .step_nonlinear_candidate(
            &thermal,
            dt,
            &state,
            &[[0.0; 3]; 6],
            &velocity,
            thermal_controls(),
            field_controls(),
            &[0.0],
            || false,
        )
        .unwrap();
    // Exact charge elimination of v=dilation*x gives a conductor-rate
    // correction 2*dilation*a. Its Lorentz transfer is signed external work,
    // not another heat source; independent authored heat still enters once.
    let mass = |t| 3.0 * sigma(t, rho, alpha, basis) / 40.0;
    let amplitude = |t| {
        bracket(0.0, a0, |a| {
            mass(t) * ((a - a0) / dt + 2.0 * dilation * a) + 2.0 * h(2.0 * a, true)
        })
    };
    let heat = |t: f64| {
        let a = amplitude(t);
        dt * mass(t) * ((a - a0) / dt + 2.0 * dilation * a).powi(2)
    };
    let t = bracket(300.0, 305.0, |t| {
        4.0 * (t - 300.0) - heat(t) - dt * external_w
    });
    let a = amplitude(t);
    let work = -2.0 * dilation * a * mass(t) * ((a - a0) / dt + 2.0 * dilation * a) * dt;
    let numerical = 2.0 * h(2.0 * a, true) * (a - a0) - energy(a, true) + energy(a0, true);
    near(
        candidate.thermal().temperature_k(CellId(7)).unwrap(),
        t,
        2e-8,
    );
    near(
        candidate.field().potential_wb()[diagonal(&problem)],
        a,
        2e-9,
    );
    let report = candidate.report();
    near(report.coupled.electrical.heat_j, heat(t), 2e-8);
    near(report.coupled.electrical.mechanical_work_j, work, 2e-8);
    near(
        report.coupled.electrical.numerical_dissipation_j,
        numerical,
        2e-8,
    );
    near(report.external_heat_j, dt * external_w, 2e-8);
    near(
        report.coupled.thermal.stored_energy_change_j,
        heat(t) + dt * external_w,
        2e-8,
    );
    near(report.total_energy_residual_j, 0.0, 2e-8);
    assert!(work > 1e-3 && heat(t) > 1e-3 && numerical > 1e-3);
}

#[test]
fn lost_heat_source_substitution_budget_and_late_cancellation_leave_both_owners_unchanged() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let card = card(0.2, 0.009, basis, true);
    let curve = curve(&card);
    let thermal = thermal(&curve, 0.0);
    let problem = problem(&card, basis);
    let state = initial(&problem, &card, basis, 1.2, true);
    let query = context(300.0);
    let cells = cells(&card, &curve, &query);
    let regions = regions();
    let network = SpatialMagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let run = |network: &SpatialMagneticThermalNetwork<'_>, controls, cancel: fn() -> bool| {
        network.step_nonlinear_candidate(
            &thermal,
            0.05,
            &state,
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            controls,
            field_controls(),
            &[0.0],
            cancel,
        )
    };
    let mut budget = thermal_controls();
    budget.maximum_circuit_trials = 1;
    assert!(matches!(
        run(&network, budget, || false),
        Err(CircuitEnthalpyError::CircuitBudget)
    ));
    let before = state.potential_wb().to_vec();
    let mut calls = 0;
    assert!(matches!(
        network.step_nonlinear_candidate(
            &thermal,
            0.05,
            &state,
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            thermal_controls(),
            field_controls(),
            &[0.0],
            || {
                calls += 1;
                calls > 400
            }
        ),
        Err(CircuitEnthalpyError::Cancelled)
    ));
    let foreign = self::card(0.4, 0.009, basis, true);
    assert!(
        state
            .step_with_conductors_candidate(
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 8],
                0.05,
                vec![Some(conductor(&foreign, 300.0, basis)); 6],
                field_controls(),
                &mut || false
            )
            .is_err()
    );
    assert!(
        state
            .step_with_conductors_candidate(
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 8],
                0.05,
                vec![None; 6],
                field_controls(),
                &mut || false
            )
            .is_err()
    );
    let mut missing = regions.clone();
    missing[0] = SpatialMagneticRegion::Fixed { tetrahedron: 0 };
    let incomplete = SpatialMagneticThermalNetwork {
        magnetic_regions: &missing,
        ..network
    };
    assert!(run(&incomplete, thermal_controls(), || false).is_err());
    assert_eq!(state.potential_wb(), before);
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(thermal.time_s(), 0.0);
    near(thermal.temperature_k(CellId(7)).unwrap(), 300.0, 1e-12);
}
