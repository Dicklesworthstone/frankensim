//! Coarse public source boundaries must not replace the original Cp(T) law.

use super::*;
use fs_conduction::ThermalBc;

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
