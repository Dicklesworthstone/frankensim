//! Graded reference inventory in the actual spatial heat/phase/adjoint owners.
use super::*;

const DENSITY: [f64; 4] = [3.0, 15.0, 3.0, 3.0];
const MASSES: [f64; 4] = [0.225, 0.325, 0.225, 0.225];

#[test]
fn graded_density_drives_the_native_nonuniform_sensible_heat_response() {
    let mesh = tetrahedron();
    let boundary = adiabatic(&mesh);
    let material = ConductivityModel::isotropic_declared(3.0).unwrap();
    let curve = sensible_curve();
    let old_h = [100.0; 4];
    let target_t = [305.0, 312.0, 307.0, 310.0];
    let target_h = target_t.map(|t| 5.0 * (t - 280.0));
    let dt = 0.2;
    // Independent affine-tet stiffness: k*V*(grad N_a dot grad N_b).
    let conduction = [
        0.5 * (3.0 * target_t[0] - target_t[1] - target_t[2] - target_t[3]),
        0.5 * (target_t[1] - target_t[0]),
        0.5 * (target_t[2] - target_t[0]),
        0.5 * (target_t[3] - target_t[0]),
    ];
    let load: [f64; 4] =
        std::array::from_fn(|a| MASSES[a] * (target_h[a] - old_h[a]) / dt + conduction[a]);
    // P1 volumetric-source matrix is V*(I+11^T)/20, independently inverted.
    let sum: f64 = load.iter().sum();
    let nodal_source = load.map(|value| (20.0 / VOLUME) * (value - sum / 5.0));
    let source =
        ScalarField::nodal("manufactured graded heating", 4, nodal_source.to_vec()).unwrap();
    with_cx(|cx| {
        let stepper = EnthalpyBackwardEuler::from_nodal_reference_density(
            cx,
            &mesh,
            &curve,
            &DENSITY,
            budget(),
        )
        .unwrap();
        for (actual, expected) in stepper.reference_nodal_masses_kg().iter().zip(MASSES) {
            assert!((actual - expected).abs() < 2e-15);
        }
        let p = problem(&mesh, &boundary, &material, &source);
        let got = stepper.advance(cx, p, None, &old_h, dt, config()).unwrap();
        for a in 0..4 {
            assert!((got.temperature[a] - target_t[a]).abs() < 2e-9);
        }
        let energy: f64 = MASSES
            .iter()
            .enumerate()
            .map(|(a, mass)| mass * (target_h[a] - old_h[a]))
            .sum();
        assert!((got.stored_energy_change_j - energy).abs() < 1e-9);
        assert!((got.stored_energy_change_j - dt * sum).abs() < 1e-9);
        assert!(got.energy_residual_j.abs() < 1e-9 && got.newton.converged);
        let averaged = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, 6.0, budget())
            .unwrap()
            .advance(cx, p, None, &old_h, dt, config())
            .unwrap();
        assert!((averaged.temperature[1] - got.temperature[1]).abs() > 1e-3);
    });
}

#[test]
fn graded_latent_storage_and_adjoint_use_the_same_immutable_reference_masses() {
    let mesh = tetrahedron();
    let boundary = adiabatic(&mesh);
    let material = ConductivityModel::isotropic_declared(3.0).unwrap();
    let curve = latent_curve();
    let source = ScalarField::uniform("latent deposition", 4.0).unwrap();
    let old = [150.0, 180.0, 210.0, 240.0];
    let dt = 3.0;
    let seed = [0.7, -0.4, 0.2, 1.1];
    with_cx(|cx| {
        let stepper = EnthalpyBackwardEuler::from_nodal_reference_density(
            cx,
            &mesh,
            &curve,
            &DENSITY,
            budget(),
        )
        .unwrap();
        let original_mass = bits(stepper.reference_nodal_masses_kg());
        let response = stepper
            .linearize_step(
                cx,
                problem(&mesh, &boundary, &material, &source),
                None,
                &old,
                dt,
                config(),
            )
            .unwrap();
        let got = response.primal();
        // The entire latent plateau has grad(T)=0 but positive material storage.
        for a in 0..4 {
            let h = old[a] + dt * (4.0 * VOLUME / 4.0) / MASSES[a];
            assert!((got.specific_enthalpy_j_kg[a] - h).abs() < 1e-9);
            assert!((got.temperature[a] - 350.0).abs() < 1e-12);
            assert!((got.liquid_mass_fraction[a] - (h - 100.0) / 200.0).abs() < 1e-10);
        }
        assert!((got.stored_energy_change_j - dt * 4.0 * VOLUME).abs() < 1e-9);
        assert_eq!(bits(stepper.reference_nodal_masses_kg()), original_mass);
        assert_eq!(
            response.temperature_pullback(cx, &seed).unwrap(),
            vec![0.0; 4]
        );
        let direction = [0.2, -0.3, 0.6, 1.0];
        let action = response.apply_jacobian(cx, &direction).unwrap();
        for a in 0..4 {
            assert!((action[a] - MASSES[a] * direction[a]).abs() < 1e-13);
        }
        let gradient = response
            .pullback(
                cx,
                &seed,
                LinearConfig {
                    tolerance: 1e-12,
                    max_iterations: 64,
                    restart: 4,
                },
            )
            .unwrap();
        let dual: [f64; 4] = std::array::from_fn(|a| seed[a] / MASSES[a]);
        let dual_sum: f64 = dual.iter().sum();
        for a in 0..4 {
            assert!((gradient.previous_specific_enthalpy[a] - seed[a]).abs() < 1e-12);
            let expected = dt * VOLUME / 20.0 * (dual_sum + dual[a]);
            assert!((gradient.source_density[a] - expected).abs() < 1e-12);
        }
    });
}

#[test]
fn prescribed_bath_charges_each_graded_storage_change_once() {
    let mesh = tetrahedron();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("bath", |_| true, ThermalBc::dirichlet(300.0).unwrap())
        .unwrap()
        .finish()
        .unwrap();
    let material = ConductivityModel::isotropic_declared(3.0).unwrap();
    let source = ScalarField::uniform("generation", 4.0).unwrap();
    let curve = sensible_curve();
    let old = [150.0, 200.0, 160.0, 180.0];
    let dt = 0.5;
    with_cx(|cx| {
        let stepper = EnthalpyBackwardEuler::from_nodal_reference_density(
            cx,
            &mesh,
            &curve,
            &DENSITY,
            budget(),
        )
        .unwrap();
        let p = problem(&mesh, &boundary, &material, &source);
        let got = stepper.advance(cx, p, None, &old, dt, config()).unwrap();
        assert_eq!(got.newton.iterations, 0);
        let energy: f64 = MASSES
            .iter()
            .enumerate()
            .map(|(a, mass)| mass * (100.0 - old[a]))
            .sum();
        assert!((got.stored_energy_change_j - energy).abs() < 1e-12);
        assert!((got.dirichlet_in_w - (energy / dt - 4.0 * VOLUME)).abs() < 1e-12);
        let next = stepper
            .advance(cx, p, None, &got.specific_enthalpy_j_kg, dt, config())
            .unwrap();
        assert_eq!(next.stored_energy_change_j, 0.0);
        assert!((next.dirichlet_in_w + 4.0 * VOLUME).abs() < 1e-12);
    });
}

#[test]
fn uniform_nodal_storage_replays_old_bits_and_original_admission() {
    let mesh = tetrahedron();
    let boundary = adiabatic(&mesh);
    let material = ConductivityModel::isotropic_declared(3.0).unwrap();
    let source = ScalarField::uniform("heating", 4.0).unwrap();
    let curve = sensible_curve();
    with_cx(|cx| {
        let uniform = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, 7.3, budget()).unwrap();
        let nodal = EnthalpyBackwardEuler::from_nodal_reference_density(
            cx,
            &mesh,
            &curve,
            &[7.3; 4],
            budget(),
        )
        .unwrap();
        assert_eq!(
            bits(uniform.reference_nodal_masses_kg()),
            bits(nodal.reference_nodal_masses_kg())
        );
        let p = problem(&mesh, &boundary, &material, &source);
        let old = [150.0, 200.0, 160.0, 180.0];
        let original = uniform.advance(cx, p, None, &old, 0.2, config()).unwrap();
        let replay = nodal.advance(cx, p, None, &old, 0.2, config()).unwrap();
        assert_eq!(
            bits(&original.specific_enthalpy_j_kg),
            bits(&replay.specific_enthalpy_j_kg)
        );
        for density in [
            &DENSITY[..3],
            &[3.0, 0.0, 3.0, 3.0],
            &[3.0, f64::NAN, 3.0, 3.0],
        ] {
            assert!(matches!(
                EnthalpyBackwardEuler::from_nodal_reference_density(
                    cx,
                    &mesh,
                    &curve,
                    density,
                    budget(),
                ),
                Err(EnthalpyError::InvalidInput(_))
            ));
        }
        assert!(matches!(
            EnthalpyBackwardEuler::from_nodal_reference_density(
                cx,
                &mesh,
                &curve,
                &DENSITY,
                EnthalpyBudget {
                    max_vertices: 3,
                    max_elements: 1
                },
            ),
            Err(EnthalpyError::Budget { .. })
        ));
    });
    with_cancelled_cx(|cx| {
        assert!(matches!(
            EnthalpyBackwardEuler::from_nodal_reference_density(
                cx,
                &mesh,
                &curve,
                &DENSITY,
                budget(),
            ),
            Err(EnthalpyError::Conduction(ConductionError::Cancelled { .. }))
        ));
    });
}
