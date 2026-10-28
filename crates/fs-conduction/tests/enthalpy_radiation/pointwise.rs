//! G1/G3 local surface integration, hot spots and enthalpy energy balances.
use super::*;
use fs_conduction::PointwiseAmbientRadiation;
use fs_conduction::transient::enthalpy::adjoint::EnthalpyAdjointError;

// Independent analytic simplex moments, not the production quadrature:
// integral lambda^alpha dA = 2 A product(alpha!)/(sum(alpha)+2)!.
// Expanding T^4 therefore gives integral lambda_i T^4 =
// A/105 sum_(a+b+c=4) (alpha_i+1) T0^a T1^b T2^c.
fn exact_face_load(t: [f64; 3], area: f64, ambient: f64) -> [f64; 3] {
    let mut out = [0.0; 3];
    for a in 0_i32..=4 {
        for b in 0_i32..=4 - a {
            let c = 4 - a - b;
            let term = t[0].powi(a) * t[1].powi(b) * t[2].powi(c);
            for (i, exponent) in [a, b, c].into_iter().enumerate() {
                out[i] += term * f64::from(exponent + 1) / 105.0;
            }
        }
    }
    out.map(|emission| {
        EPSILON * STEFAN_BOLTZMANN_W_M2_K4 * area * (emission - ambient.powi(4) / 3.0)
    })
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:.3e}"
    );
}

#[test]
fn local_surface_quartic_and_tangent_match_independent_simplex_moments() {
    let mesh = regular_tet();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .remainder("exterior", ThermalBc::adiabatic())
        .unwrap()
        .finish()
        .unwrap();
    let patches = [patch(300.0, 1300.0)];
    let t = [300.0, 400.0, 600.0, 1200.0];
    with_cx(|cx| {
        let radiation = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
        let got = radiation.assemble(cx, &t).unwrap();
        let mut expected = [0.0; 4];
        for face in mesh.boundary() {
            let local = exact_face_load(face.vertices.map(|v| t[v as usize]), face.area, 300.0);
            for (i, &v) in face.vertices.iter().enumerate() {
                expected[v as usize] += local[i];
            }
        }
        for (actual, reference) in got.outward_load_w.iter().zip(expected) {
            close(*actual, reference, 3e-11);
        }
        close(got.report.outward_heat_w, expected.into_iter().sum(), 8e-11);
        let report = &got.report.patches[0];
        let mean_model = patches[0]
            .heat_flux_w_m2(report.mean_surface_temperature_k)
            .unwrap()
            * report.area_m2;
        assert!(
            report.outward_heat_w > 1.4 * mean_model,
            "the resolved hot spot must radiate more than a uniform surface at the same mean"
        );

        let direction = [0.7, -0.9, 0.2, 1.3];
        let mut tangent = [0.0; 4];
        got.tangent_apply(cx, &direction, &mut tangent).unwrap();
        let delta = 1e-3;
        let plus = radiation
            .assemble(
                cx,
                &std::array::from_fn::<_, 4, _>(|i| t[i] + delta * direction[i]),
            )
            .unwrap();
        let minus = radiation
            .assemble(
                cx,
                &std::array::from_fn::<_, 4, _>(|i| t[i] - delta * direction[i]),
            )
            .unwrap();
        for i in 0..4 {
            close(
                tangent[i],
                (plus.outward_load_w[i] - minus.outward_load_w[i]) / (2.0 * delta),
                3e-8,
            );
        }
        assert!(
            direction
                .iter()
                .zip(tangent)
                .map(|(&a, b)| a * b)
                .sum::<f64>()
                > 0.0
        );
        let other = [-0.4, 0.8, 1.1, -0.2];
        let mut applied_other = [0.0; 4];
        got.tangent_apply(cx, &other, &mut applied_other).unwrap();
        close(
            other.iter().zip(tangent).map(|(&a, b)| a * b).sum(),
            direction
                .iter()
                .zip(applied_other)
                .map(|(&a, b)| a * b)
                .sum(),
            2e-13,
        );
    });
}

#[test]
fn pointwise_implicit_cooling_matches_uniform_scalar_balance_and_keeps_air_heat_separate() {
    let mesh = regular_tet();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .remainder("exterior", ThermalBc::robin(HTC, 300.0).unwrap())
        .unwrap()
        .finish()
        .unwrap();
    let material = ConductivityModel::isotropic_declared(K).unwrap();
    let source = ScalarField::Uniform(0.0);
    let curve = phase_curve(false);
    let patches = [patch(300.0, 600.0)];
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
        element_materials: None,
    };
    with_cx(|cx| {
        let radiation = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
        let stepper = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        let old = [CP * (450.0 - 200.0); 4];
        let dt = 0.2;
        let solution = stepper
            .advance_with_pointwise_ambient_radiation(
                cx,
                problem,
                None,
                &old,
                dt,
                config(),
                &radiation,
            )
            .unwrap();
        let reference = scalar_mean(450.0, 300.0, dt);
        for &temperature in &solution.temperature {
            close(temperature, reference, 2e-8);
        }
        close(
            solution.robin_out_w,
            HTC * 3.0_f64.sqrt() * (reference - 300.0),
            2e-8,
        );
        let radiant = solution.pointwise_radiation.as_ref().unwrap();
        close(
            solution.stored_energy_change_j + dt * (radiant.outward_heat_w + solution.robin_out_w),
            0.0,
            1e-8,
        );
        close(solution.energy_residual_j, 0.0, 1e-8);
        assert_eq!(old, [CP * (450.0 - 200.0); 4]);
    });
}

#[test]
fn radiation_changes_latent_enthalpy_without_creating_an_apparent_temperature_rise() {
    let mesh = regular_tet();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .remainder("exterior", ThermalBc::adiabatic())
        .unwrap()
        .finish()
        .unwrap();
    let material = ConductivityModel::isotropic_declared(K).unwrap();
    let source = ScalarField::Uniform(0.0);
    let curve = phase_curve(true);
    let patches = [patch(500.0, 550.0)];
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
        element_materials: None,
    };
    with_cx(|cx| {
        let radiation = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
        let stepper = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        let old = [1500.0; 4];
        let dt = 0.05;
        let solution = stepper
            .advance_with_pointwise_ambient_radiation(
                cx,
                problem,
                None,
                &old,
                dt,
                config(),
                &radiation,
            )
            .unwrap();
        let mass = RHO * 2.0_f64.sqrt() / 12.0;
        let heat_in = EPSILON
            * STEFAN_BOLTZMANN_W_M2_K4
            * 3.0_f64.sqrt()
            * (500.0_f64.powi(4) - 350.0_f64.powi(4));
        for i in 0..4 {
            close(solution.temperature[i], 350.0, 1e-12);
            close(
                solution.specific_enthalpy_j_kg[i],
                1500.0 + dt * heat_in / mass,
                1e-9,
            );
            close(
                solution.liquid_mass_fraction[i],
                (500.0 + dt * heat_in / mass) / 2000.0,
                1e-12,
            );
        }
        close(solution.stored_energy_change_j, dt * heat_in, 1e-8);
        close(solution.robin_out_w, 0.0, 0.0);
        assert!(
            solution
                .pointwise_radiation
                .as_ref()
                .unwrap()
                .outward_heat_w
                < 0.0
        );
    });
}

#[test]
fn nonuniform_radiating_wall_uses_local_quartic_and_accounts_for_prescribed_reaction() {
    let mesh = regular_tet();
    let wall = 400.0;
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "bath",
            |face| face.centroid[2] == 0.0,
            ThermalBc::dirichlet(wall).unwrap(),
        )
        .unwrap()
        .remainder("exterior", ThermalBc::robin(HTC, 300.0).unwrap())
        .unwrap()
        .finish()
        .unwrap();
    let material = ConductivityModel::isotropic_declared(K).unwrap();
    let source = ScalarField::Uniform(0.0);
    let curve = phase_curve(false);
    let patches = [patch(500.0, 600.0)];
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
        element_materials: None,
    };
    let dt = 0.2;
    let volume = 2.0_f64.sqrt() / 12.0;
    let face_area = 3.0_f64.sqrt() / 4.0;
    let (mut low, mut high) = (300.0, 500.0);
    for _ in 0..80 {
        let t = (low + high) / 2.0;
        let flux = 3.0 * exact_face_load([t, wall, wall], face_area, 500.0)[0];
        let convection = 3.0 * face_area * HTC * ((t + wall) / 6.0 - 300.0 / 3.0);
        let residual = RHO * CP * volume / 4.0 * (t - 300.0)
            + dt * (1.5 * K * volume * (t - wall) + convection + flux);
        if residual > 0.0 {
            high = t;
        } else {
            low = t;
        }
    }
    with_cx(|cx| {
        let radiation = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
        let stepper = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        let old = [1000.0; 4];
        let solution = stepper
            .advance_with_pointwise_ambient_radiation(
                cx,
                problem,
                None,
                &old,
                dt,
                config(),
                &radiation,
            )
            .unwrap();
        close(solution.temperature[3], (low + high) / 2.0, 2e-8);
        for &t in &solution.temperature[..3] {
            close(t, wall, 0.0);
        }
        let radiative_w = solution
            .pointwise_radiation
            .as_ref()
            .unwrap()
            .outward_heat_w;
        let independent_w: f64 = mesh
            .boundary()
            .iter()
            .filter(|f| f.centroid[2] != 0.0)
            .map(|f| {
                exact_face_load(
                    f.vertices.map(|v| solution.temperature[v as usize]),
                    f.area,
                    500.0,
                )
                .into_iter()
                .sum::<f64>()
            })
            .sum();
        close(radiative_w, independent_w, 1e-10);
        close(
            solution.stored_energy_change_j
                - dt * (solution.dirichlet_in_w - solution.robin_out_w - radiative_w),
            0.0,
            1e-8,
        );
    });
}

#[test]
fn unsupported_hot_spots_cancellation_and_radiation_blind_adjoints_refuse() {
    let mesh = regular_tet();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .remainder("exterior", ThermalBc::adiabatic())
        .unwrap()
        .finish()
        .unwrap();
    let material = ConductivityModel::isotropic_declared(K).unwrap();
    let source = ScalarField::Uniform(0.0);
    let curve = phase_curve(false);
    let patches = [patch(300.0, 600.0)];
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
        element_materials: None,
    };
    with_cx(|cx| {
        let radiation = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
        assert!(matches!(
            radiation.assemble(cx, &[300.0, 300.0, 300.0, 800.0]),
            Err(ConductionError::MaterialQuery { .. })
        ));
        with_cancelled_cx(|cancelled| {
            assert!(matches!(
                radiation.assemble(cancelled, &[300.0; 4]),
                Err(ConductionError::Cancelled { .. })
            ))
        });
        let stepper = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
        let old = [1000.0; 4];
        let equilibrium = stepper
            .advance_with_pointwise_ambient_radiation(
                cx,
                problem,
                None,
                &old,
                0.1,
                config(),
                &radiation,
            )
            .unwrap();
        close(
            equilibrium
                .pointwise_radiation
                .as_ref()
                .unwrap()
                .outward_heat_w,
            0.0,
            0.0,
        );
        assert!(matches!(
            stepper.linearize_accepted(cx, problem, None, &old, 0.1, config(), equilibrium),
            Err(EnthalpyAdjointError::InvalidInput(_))
        ));
    });
}
