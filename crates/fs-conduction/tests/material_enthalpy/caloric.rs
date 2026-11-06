//! Coarse public source boundaries must not replace the original Cp(T) law.

use super::*;
use fs_conduction::ThermalBc;
use fs_conduction::transient::enthalpy::material::{
    ExactMaterialEnthalpySource, MaterialEnthalpyError,
};

fn resolve_exact(
    card: &MaterialCard,
    reference: &QueryPoint,
    conductivity: &PropertyKey,
) -> Result<MaterialEnthalpyTransport, MaterialEnthalpyError> {
    MaterialEnthalpyTransport::from_isotropic_card_exact(
        card,
        ExactMaterialEnthalpySource {
            caloric: HeatCapacitySourceQuery {
                heat_capacity: &PropertyKey::new("cp", CP_DIMS),
                density: &PropertyKey::new("rho", Density::DIMS),
                temperature_axis: "temperature",
                lower: &point(300.0, 0.12),
                upper: &point(400.0, 0.12),
                selection: SelectionPolicy::SingleClaimOnly,
            },
            conductivity,
            reference,
            phase: SolidLiquidPhase::Solid,
            lower_specific_enthalpy_j_kg: 0.0,
        },
    )
}

fn exact_orthotropic(
    card: &MaterialCard,
    principal_axes: [[f64; 3]; 3],
) -> Result<MaterialEnthalpyTransport, MaterialEnthalpyError> {
    MaterialEnthalpyTransport::from_orthotropic_card_exact(
        card,
        ExactMaterialEnthalpySource {
            caloric: HeatCapacitySourceQuery {
                heat_capacity: &PropertyKey::new("cp", CP_DIMS),
                density: &PropertyKey::new("rho", Density::DIMS),
                temperature_axis: "temperature",
                lower: &point(300.0, 0.12),
                upper: &point(400.0, 0.12),
                selection: SelectionPolicy::SingleClaimOnly,
            },
            conductivity: [
                &PropertyKey::new("k_longitudinal", K_DIMS),
                &PropertyKey::new("k_radial", K_DIMS),
                &PropertyKey::new("k_tangential", K_DIMS),
            ],
            reference: &point(325.0, 0.12),
            phase: SolidLiquidPhase::Solid,
            lower_specific_enthalpy_j_kg: 0.0,
        },
        principal_axes,
    )
}

fn coarse_material(card: &MaterialCard) -> MaterialEnthalpyTransport {
    let source = MaterialEnthalpyTransport::from_isotropic_card(
        card,
        MaterialEnthalpySource {
            caloric: HeatCapacitySourceQuery {
                heat_capacity: &PropertyKey::new("cp", CP_DIMS),
                density: &PropertyKey::new("rho", Density::DIMS),
                temperature_axis: "temperature",
                lower: &point(300.0, 0.12),
                upper: &point(400.0, 0.12),
                selection: SelectionPolicy::SingleClaimOnly,
            },
            conductivity: &PropertyKey::new("k", K_DIMS),
            reference: &point(325.0, 0.12),
            phase: SolidLiquidPhase::Solid,
            lower_specific_enthalpy_j_kg: 0.0,
            maximum_interpolation_error_j_kg: 1e6,
        },
    )
    .unwrap();
    assert_eq!(source.caloric().curve().knots().len(), 2);
    assert_eq!(source.reference_density_kg_m3(), 8.75);
    source
}

#[test]
fn uniform_heating_uses_exact_source_capacity_with_two_public_caloric_knots() {
    let source = coarse_material(&card(1.0, true));
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let heat = ScalarField::Uniform(875.0);
    with_cx(|cx| {
        let stepper = source.prepare_uniform(cx, &mesh, budget()).unwrap();
        let old = [0.0; 4];
        let solved = stepper
            .advance(
                cx,
                ConductionProblem {
                    mesh: &mesh,
                    boundary: &boundary,
                    material: source.conductivity(),
                    element_materials: None,
                    source: &heat,
                },
                None,
                &old,
                10.0,
                config(),
            )
            .unwrap();
        // q*dt/rho0=1000 J/kg. Inverting the original integral
        // h=1000*dT+5*dT^2 gives this cancellation-free quadratic root.
        let expected_h = 875.0_f64 * 10.0 / 8.75;
        let expected_t =
            300.0 + 2.0 * expected_h / (1000.0 + (1_000_000.0 + 20.0 * expected_h).sqrt());
        for (&h, &t) in solved
            .specific_enthalpy_j_kg
            .iter()
            .zip(&solved.temperature)
        {
            assert!((h - expected_h).abs() < 1e-7);
            assert!((t - expected_t).abs() < 1e-9);
        }
        let source_heat = 875.0 * 10.0 / 6.0;
        assert!((solved.stored_energy_change_j - source_heat).abs() < 1e-7);
        assert!((solved.source_w * 10.0 - source_heat).abs() < 1e-7);
        assert!(solved.energy_residual_j.abs() < 1e-7);
        assert_eq!(solved.dirichlet_in_w, 0.0);
        for mass in stepper.reference_nodal_masses_kg() {
            assert!((mass - 8.75 / 24.0).abs() < 1e-15);
        }
        assert_eq!(old, [0.0; 4]);
    });
}

#[test]
fn prescribed_temperature_inverts_original_capacity_and_charges_boundary_reaction_once() {
    let source = coarse_material(&card(1.0, true));
    let mesh = mesh();
    let prescribed = 333.3;
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "prescribed source temperature",
            |_| true,
            ThermalBc::dirichlet(prescribed).unwrap(),
        )
        .unwrap()
        .finish()
        .unwrap();
    let heat = ScalarField::Uniform(0.0);
    let dt = 10.0;
    with_cx(|cx| {
        let stepper = source.prepare_uniform(cx, &mesh, budget()).unwrap();
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: source.conductivity(),
            element_materials: None,
            source: &heat,
        };
        let old = [0.0; 4];
        let solved = stepper
            .advance(cx, problem, None, &old, dt, config())
            .unwrap();
        let delta_temperature = prescribed - 300.0;
        let expected_h = 1000.0 * delta_temperature + 5.0 * delta_temperature * delta_temperature;
        let expected_energy = (8.75 / 6.0) * expected_h;
        for (&h, &t) in solved
            .specific_enthalpy_j_kg
            .iter()
            .zip(&solved.temperature)
        {
            assert!((h - expected_h).abs() < 1e-9);
            assert!((t - prescribed).abs() < 1e-11);
        }
        assert!((solved.stored_energy_change_j - expected_energy).abs() < 1e-7);
        assert!((solved.dirichlet_in_w * dt - expected_energy).abs() < 1e-7);
        assert_eq!(solved.source_w, 0.0);
        assert_eq!(solved.neumann_out_w, 0.0);
        assert_eq!(solved.robin_out_w, 0.0);
        assert!(solved.energy_residual_j.abs() < 1e-7);

        // Holding the same bath after the initial jump adds no further heat.
        let held = stepper
            .advance(
                cx,
                problem,
                None,
                &solved.specific_enthalpy_j_kg,
                dt,
                config(),
            )
            .unwrap();
        assert!(held.stored_energy_change_j.abs() < 1e-9);
        assert!(held.dirichlet_in_w.abs() < 1e-10);
        assert!(held.energy_residual_j.abs() < 1e-9);
        assert_eq!(old, [0.0; 4]);
    });
}

#[test]
fn compact_source_storage_heats_through_the_native_solver_without_subdivision() {
    let card = card(1.0, true);
    let source = resolve_exact(&card, &point(325.0, 0.12), &PropertyKey::new("k", K_DIMS)).unwrap();
    let bounded = material(&card);
    // The same strongly varying Cp needs thousands of boundary chords at the
    // existing 0.001 J/kg allowance. Analytical storage needs only its source.
    assert!(bounded.caloric().curve().knots().len() > 1000);
    assert_eq!(source.caloric().curve().knots().len(), 2);
    assert_eq!(source.caloric().receipts(), bounded.caloric().receipts());
    assert_eq!(source.reference_density(), bounded.reference_density());
    assert_eq!(source.card_identity(), card.content_hash());
    assert_eq!(
        source.reference_material().curve.identity(),
        source.caloric().curve().identity()
    );

    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let heat = ScalarField::Uniform(875.0);
    let dt = 1000.0;
    with_cx(|cx| {
        let stepper = source.prepare_uniform(cx, &mesh, budget()).unwrap();
        let old = [0.0; 4];
        let solved = stepper
            .advance(
                cx,
                ConductionProblem {
                    mesh: &mesh,
                    boundary: &boundary,
                    material: source.conductivity(),
                    element_materials: None,
                    source: &heat,
                },
                None,
                &old,
                dt,
                config(),
            )
            .unwrap();
        // m is frozen at rho0=8.75, while Cp grows by over 70% during this step.
        // Independent integral/inverse: h=1000*dT+5*dT^2, q*dt/rho0=100000.
        let expected_h = 875.0_f64 * dt / 8.75;
        let expected_t =
            300.0 + 2.0 * expected_h / (1000.0 + (1_000_000.0 + 20.0 * expected_h).sqrt());
        for (&h, &t) in solved
            .specific_enthalpy_j_kg
            .iter()
            .zip(&solved.temperature)
        {
            assert!((h - expected_h).abs() < 1e-7);
            assert!((t - expected_t).abs() < 1e-9);
            let expected_k = 2.0 + 0.02 * (t - 300.0);
            let tensor = source.conductivity().tensor_at(t).unwrap();
            assert!((tensor[0][0] - expected_k).abs() < 1e-13);
            assert_eq!(tensor, bounded.conductivity().tensor_at(t).unwrap());
        }
        let source_heat = 875.0 * dt / 6.0;
        assert!((solved.stored_energy_change_j - source_heat).abs() < 1e-7);
        assert!((solved.source_w * dt - source_heat).abs() < 1e-7);
        assert!(solved.energy_residual_j.abs() < 1e-7);
        assert_eq!(solved.dirichlet_in_w, 0.0);
        for mass in stepper.reference_nodal_masses_kg() {
            assert!((mass - 8.75 / 24.0).abs() < 1e-15);
        }
        assert_eq!(source.reference_density_kg_m3(), 8.75);
        assert_eq!(old, [0.0; 4]);
    });
    with_cancelled_cx(|cx| {
        assert!(source.prepare_uniform(cx, &mesh, budget()).is_err());
    });
}

#[test]
fn compact_directional_storage_preserves_nonlinear_spatial_equilibration() {
    let card = card(1.0, true);
    let compact = exact_orthotropic(&card, TURNED_AXES).unwrap();
    let bounded = orthotropic(&card, TURNED_AXES);
    assert_eq!(compact.caloric().curve().knots().len(), 2);
    assert_eq!(compact.caloric().receipts(), bounded.caloric().receipts());
    assert_eq!(compact.reference_density(), bounded.reference_density());
    let expected_tensor = [[3.0, 0.0, 0.0], [0.0, 30.0, 0.0], [0.0, 0.0, 1.5]];
    assert_eq!(
        compact.conductivity().tensor_at(350.0).unwrap(),
        expected_tensor
    );

    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let heat = ScalarField::Uniform(0.0);
    let old = [10_000.0, 20_000.0, 10_000.0, 10_000.0];
    with_cx(|cx| {
        let mut answers = Vec::new();
        for source in [&compact, &bounded] {
            let stepper = source.prepare_uniform(cx, &mesh, budget()).unwrap();
            let solved = stepper
                .advance(
                    cx,
                    ConductionProblem {
                        mesh: &mesh,
                        boundary: &boundary,
                        material: source.conductivity(),
                        element_materials: None,
                        source: &heat,
                    },
                    None,
                    &old,
                    10.0,
                    config(),
                )
                .unwrap();
            let energy_change: f64 = solved
                .specific_enthalpy_j_kg
                .iter()
                .zip(old)
                .map(|(now, before)| (8.75 / 24.0) * (now - before))
                .sum();
            assert!(energy_change.abs() < 1e-7);
            assert!(solved.stored_energy_change_j.abs() < 1e-7);
            assert!(solved.energy_residual_j.abs() < 1e-7);
            assert!(solved.specific_enthalpy_j_kg[1] < old[1]);
            answers.push(solved);
        }
        for (compact, bounded) in answers[0].temperature.iter().zip(&answers[1].temperature) {
            assert!((compact - bounded).abs() < 1e-9);
        }
        for (compact, bounded) in answers[0]
            .specific_enthalpy_j_kg
            .iter()
            .zip(&answers[1].specific_enthalpy_j_kg)
        {
            assert!((compact - bounded).abs() < 1e-6);
        }
        assert_eq!(old, [10_000.0, 20_000.0, 10_000.0, 10_000.0]);
    });
}

#[test]
fn compact_material_keeps_source_context_and_directional_admission() {
    let card = card(1.0, true);
    let conductivity = PropertyKey::new("k", K_DIMS);
    for reference in [point(325.0, 0.13), point(401.0, 0.12), point(299.0, 0.12)] {
        assert!(resolve_exact(&card, &reference, &conductivity).is_err());
    }
    assert!(
        resolve_exact(
            &card,
            &point(325.0, 0.12),
            &PropertyKey::new("rho", Density::DIMS),
        )
        .is_err()
    );
    assert!(resolve_exact(&super::card(1.0, false), &point(325.0, 0.12), &conductivity,).is_err());
    assert!(exact_orthotropic(&card, [[0.0; 3]; 3]).is_err());
}
