//! G1/G3: independent directional phase laws on the existing enthalpy field.
use super::*;

const VALUES: [[f64; 3]; 3] = [[2.0, 8.0, 3.0], [7.0, 1.0, 9.0], [5.0, 4.0, 6.0]];
const AXES: [[f64; 3]; 3] = [[0.6, 0.8, 0.0], [-0.8, 0.6, 0.0], [0.0, 0.0, 1.0]];

fn sources<'a>(
    card: &MaterialCard,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
) -> [PreparedPhaseConductivity<'a>; 3] {
    ["k0", "k1", "k2"].map(|key| prepare_direction(card, curve, key, [-100.0, 500.0], 150_000.0))
}

// For h_h=-20+100*x on the unit right tetrahedron the cross-sectional
// probability density is 3*(1-x)^2. These are exact elementary hinge integrals.
fn sectional_means() -> [f64; 3] {
    VALUES.map(|[a, b, c]| {
        let s1 = (b - a) / 50.0;
        let s2 = (c - b) / 50.0;
        a + s1 * 80.0_f64.powi(4) / (4.0 * 100.0_f64.powi(3))
            + (s2 - s1) * 30.0_f64.powi(4) / (4.0 * 100.0_f64.powi(3))
    })
}

// Closed form for four distinct nodal enthalpies. This uses a truncated-power
// simplex integral, independently of the production B-spline/Gauss quadrature.
fn mean_for(values: [f64; 3], h: &[f64]) -> f64 {
    let hinge = |cut: f64| {
        (0..4)
            .map(|i| {
                let denominator = (0..4)
                    .filter(|&j| i != j)
                    .map(|j| h[i] - h[j])
                    .product::<f64>();
                assert_ne!(denominator, 0.0);
                (h[i] - cut).max(0.0).powi(4) / (4.0 * denominator)
            })
            .sum::<f64>()
    };
    let [a, b, c] = values;
    let s1 = (b - a) / 50.0;
    let s2 = (c - b) / 50.0;
    a + s1 * hinge(0.0) + (s2 - s1) * hinge(50.0) - s2 * hinge(100.0)
}

fn tensor_action(k: [f64; 3], gradient: [f64; 3]) -> [f64; 3] {
    // Analytic rotation about z by cos=3/5, sin=4/5; the third axis is unchanged.
    let xx = (9.0 * k[0] + 16.0 * k[1]) / 25.0;
    let yy = (16.0 * k[0] + 9.0 * k[1]) / 25.0;
    let xy = 12.0 * (k[0] - k[1]) / 25.0;
    [
        xx * gradient[0] + xy * gradient[1],
        xy * gradient[0] + yy * gradient[1],
        k[2] * gradient[2],
    ]
}

fn physical_residual(h: &[f64], old: &[f64]) -> [f64; 4] {
    let t = h
        .iter()
        .map(|&value| temperature(value))
        .collect::<Vec<_>>();
    let k = VALUES.map(|values| mean_for(values, h));
    let action = tensor_action(k, [t[1] - t[0], t[2] - t[0], t[3] - t[0]]);
    let divergence = [-action.iter().sum::<f64>(), action[0], action[1], action[2]];
    std::array::from_fn(|i| {
        MASS * (h[i] - old[i]) + DT * (divergence[i] / 6.0 - (Q.iter().sum::<f64>() + Q[i]) / 120.0)
    })
}

#[test]
fn principal_phase_sources_rotate_flux_and_preserve_the_nodal_heat_balance() {
    let mesh = mesh();
    let card = orthotropic_card(VALUES);
    let curve = caloric(&card);
    let sources = sources(&card, &curve);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let zero = ScalarField::Uniform(0.0);
    let h = [-20.0, 80.0, -20.0, -20.0];
    let expected = tensor_action(sectional_means(), [20.0, 0.0, 0.0]);
    assert!(expected[1].abs() > 1.0); // A scalar law cannot produce this transverse heat flux.
    with_cx(|cx| {
        let law =
            MaterialPhaseTransport::try_new_orthotropic(cx, sources.each_ref(), AXES, 2.0).unwrap();
        assert_eq!(law.principal_axes(), Some(AXES));
        assert_eq!(law.enthalpy_bounds_j_kg(), [-100.0, 500.0]);
        assert_eq!(law.conductivity_sources().len(), 3);
        for (retained, original) in law.conductivity_sources().iter().zip(&sources) {
            assert!(std::ptr::eq(*retained, original));
        }
        let solver = law.prepare_uniform(cx, &mesh, budget()).unwrap();
        let endpoint = solver
            .evaluate_endpoint(
                cx,
                law.problem(&mesh, &boundary, &zero),
                None,
                &h,
                &h,
                DT,
                config(),
                None,
            )
            .unwrap();
        near(
            endpoint.residual_j[0],
            -DT * expected.iter().sum::<f64>() / 6.0,
            2e-12,
        );
        for i in 0..3 {
            near(endpoint.residual_j[i + 1], DT * expected[i] / 6.0, 2e-12);
        }
        near(endpoint.energy_residual_j, 0.0, 2e-12);
        let flux = law.element_heat_flux(cx, &mesh, &h).unwrap()[0];
        for i in 0..3 {
            near(flux[i], -expected[i], 2e-11);
        }
        assert_eq!(solver.reference_nodal_masses_kg(), &[MASS; 4]);

        // Rotate geometry and every principal direction together by a cyclic
        // coordinate permutation. This is an exact proper rotation in binary64.
        let rotate = |v: [f64; 3]| [v[2], v[0], v[1]];
        let rotated_mesh = ConductionMesh::new(
            TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
            vec![
                [0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
            ],
        )
        .unwrap();
        let rotated = MaterialPhaseTransport::try_new_orthotropic(
            cx,
            sources.each_ref(),
            AXES.map(rotate),
            2.0,
        )
        .unwrap();
        let rotated_flux = rotated.element_heat_flux(cx, &rotated_mesh, &h).unwrap()[0];
        for (actual, expected) in rotated_flux.into_iter().zip(rotate(flux)) {
            near(actual, expected, 2e-11);
        }
        assert_ne!(law.identity(), rotated.identity());
    });
}

#[test]
fn latent_principal_derivatives_reach_the_native_newton_operator() {
    let mesh = mesh();
    let card = orthotropic_card(VALUES);
    let curve = caloric(&card);
    let sources = sources(&card, &curve);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let q = ScalarField::nodal("orthotropic phase heating", 4, Q.to_vec()).unwrap();
    with_cx(|cx| {
        let law =
            MaterialPhaseTransport::try_new_orthotropic(cx, sources.each_ref(), AXES, 2.0).unwrap();
        let solver = law.prepare_uniform(cx, &mesh, budget()).unwrap();
        let problem = law.problem(&mesh, &boundary, &q);
        let response = solver
            .linearize_step(cx, problem, None, &OLD, DT, config())
            .unwrap();
        let h = &response.primal().specific_enthalpy_j_kg;
        assert!(h[0] < 0.0 && h[3] > 100.0);
        assert!(h[1] > 50.0 && h[1] < 100.0 && h[2] > 0.0 && h[2] < 50.0);
        for residual in physical_residual(h, &OLD) {
            near(residual, 0.0, 3e-10);
        }
        near(
            response.primal().stored_energy_change_j,
            DT * Q.iter().sum::<f64>() / 24.0,
            1e-10,
        );
        for column in 0..4 {
            let mut direction = [0.0; 4];
            direction[column] = 1.0;
            let action = response.apply_jacobian(cx, &direction).unwrap();
            let mut plus = h.clone();
            let mut minus = h.clone();
            plus[column] += 0.01;
            minus[column] -= 0.01;
            let rp = solver
                .evaluate_endpoint(cx, problem, None, &OLD, &plus, DT, config(), None)
                .unwrap();
            let rm = solver
                .evaluate_endpoint(cx, problem, None, &OLD, &minus, DT, config(), None)
                .unwrap();
            let independent_plus = physical_residual(&plus, &OLD);
            let independent_minus = physical_residual(&minus, &OLD);
            for row in 0..4 {
                near(
                    action[row],
                    (rp.residual_j[row] - rm.residual_j[row]) / 0.02,
                    3e-9,
                );
                near(
                    action[row],
                    (independent_plus[row] - independent_minus[row]) / 0.02,
                    3e-9,
                );
            }
            if column == 1 || column == 2 {
                // These nodes stay at the melting temperature: only dK/dh can
                // produce a nonlocal Jacobian entry from their h perturbation.
                assert!(action[0].abs() > 1e-6);
            }
        }
    });
}

#[test]
fn equal_principal_sources_recover_isotropic_transport_in_a_rotated_frame() {
    let mesh = mesh();
    let card = card([2.0, 8.0, 3.0]);
    let curve = caloric(&card);
    let source = prepare(&card, &curve, [-100.0, 500.0]);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let q = ScalarField::nodal("isotropic comparison", 4, Q.to_vec()).unwrap();
    with_cx(|cx| {
        let isotropic = MaterialPhaseTransport::try_new_isotropic(cx, &source, 2.0).unwrap();
        let orthotropic =
            MaterialPhaseTransport::try_new_orthotropic(cx, [&source; 3], AXES, 2.0).unwrap();
        assert_eq!(isotropic.conductivity_sources().len(), 1);
        assert_eq!(isotropic.principal_axes(), None);
        let a = isotropic
            .prepare_uniform(cx, &mesh, budget())
            .unwrap()
            .advance(
                cx,
                isotropic.problem(&mesh, &boundary, &q),
                None,
                &OLD,
                DT,
                config(),
            )
            .unwrap();
        let b = orthotropic
            .prepare_uniform(cx, &mesh, budget())
            .unwrap()
            .advance(
                cx,
                orthotropic.problem(&mesh, &boundary, &q),
                None,
                &OLD,
                DT,
                config(),
            )
            .unwrap();
        for (&a, &b) in a
            .specific_enthalpy_j_kg
            .iter()
            .zip(&b.specific_enthalpy_j_kg)
        {
            near(a, b, 1e-10);
        }
        let a_flux = isotropic.element_heat_flux(cx, &mesh, &OLD).unwrap()[0];
        let b_flux = orthotropic.element_heat_flux(cx, &mesh, &OLD).unwrap()[0];
        for (a, b) in a_flux.into_iter().zip(b_flux) {
            near(a, b, 2e-11);
        }
    });
}

#[test]
fn directional_domains_intersect_and_incompatible_physical_sources_refuse() {
    let mesh = mesh();
    let card = orthotropic_card(VALUES);
    let curve = caloric(&card);
    let broad = prepare_direction(&card, &curve, "k0", [-100.0, 500.0], 150_000.0);
    let middle = prepare_direction(&card, &curve, "k1", [-20.0, 200.0], 150_000.0);
    let narrow = prepare_direction(&card, &curve, "k2", [-10.0, 100.0], 150_000.0);
    let wrong_pressure = prepare_direction(&card, &curve, "k1", [-20.0, 200.0], 160_000.0);
    let disjoint = prepare_direction(&card, &curve, "k2", [210.0, 300.0], 150_000.0);
    let other_curve = caloric(&card);
    let cloned_chart = prepare_direction(&card, &other_curve, "k2", [-10.0, 100.0], 150_000.0);
    // The original f coordinates differ at each clipped source's first point.
    // They still describe the same caloric phase path and fixed external state.
    let latent_sources = [
        prepare_direction(&card, &curve, "k0", [0.0, 90.0], 150_000.0),
        prepare_direction(&card, &curve, "k1", [10.0, 80.0], 150_000.0),
        prepare_direction(&card, &curve, "k2", [20.0, 70.0], 150_000.0),
    ];
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let zero = ScalarField::Uniform(0.0);
    with_cx(|cx| {
        let law =
            MaterialPhaseTransport::try_new_orthotropic(cx, [&broad, &middle, &narrow], AXES, 2.0)
                .unwrap();
        assert_eq!(law.enthalpy_bounds_j_kg(), [-10.0, 100.0]);
        assert_eq!(
            law.conductivity_sources()[0].enthalpy_bounds_j_kg(),
            [-100.0, 500.0]
        );
        let solver = law.prepare_uniform(cx, &mesh, budget()).unwrap();
        let problem = law.problem(&mesh, &boundary, &zero);
        let latent =
            MaterialPhaseTransport::try_new_orthotropic(cx, latent_sources.each_ref(), AXES, 2.0)
                .unwrap();
        assert_eq!(latent.enthalpy_bounds_j_kg(), [20.0, 70.0]);
        let equilibrium = [25.0, 30.0, 45.0, 65.0];
        assert_eq!(
            latent.element_heat_flux(cx, &mesh, &equilibrium).unwrap(),
            vec![[0.0; 3]]
        );
        let same = latent
            .prepare_uniform(cx, &mesh, budget())
            .unwrap()
            .advance(
                cx,
                latent.problem(&mesh, &boundary, &zero),
                None,
                &equilibrium,
                DT,
                config(),
            )
            .unwrap();
        assert_eq!(same.specific_enthalpy_j_kg, equilibrium);
        let outside = [-15.0; 4];
        assert!(matches!(
            solver.advance(cx, problem, None, &outside, DT, config()),
            Err(EnthalpyError::SourceEnthalpyOutsideDomain {
                vertex: 0,
                lower_j_kg: -10.0,
                ..
            })
        ));
        assert_eq!(outside, [-15.0; 4]);
        let upper = [100.0; 4];
        let accepted = solver
            .advance(cx, problem, None, &upper, DT, config())
            .unwrap();
        assert_eq!(accepted.specific_enthalpy_j_kg, upper);
        assert!(matches!(
            solver.linearize_accepted(cx, problem, None, &upper, DT, config(), accepted),
            Err(EnthalpyAdjointError::MaterialKink)
        ));
        for incompatible in [&wrong_pressure, &disjoint, &cloned_chart] {
            let rows = if std::ptr::eq(incompatible, &wrong_pressure) {
                [&broad, incompatible, &narrow]
            } else {
                [&broad, &middle, incompatible]
            };
            assert!(matches!(
                MaterialPhaseTransport::try_new_orthotropic(cx, rows, AXES, 2.0),
                Err(EnthalpyError::InvalidInput(_))
            ));
        }
        let mut invalid_axes = AXES;
        invalid_axes[1] = invalid_axes[0];
        assert!(
            MaterialPhaseTransport::try_new_orthotropic(
                cx,
                [&broad, &middle, &narrow],
                invalid_axes,
                2.0
            )
            .is_err()
        );
        assert!(
            MaterialPhaseTransport::try_new_orthotropic(cx, [&broad, &middle, &narrow], AXES, 0.0)
                .is_err()
        );
        with_cancelled_cx(|cancelled| {
            assert!(
                MaterialPhaseTransport::try_new_orthotropic(
                    cancelled,
                    [&broad, &middle, &narrow],
                    AXES,
                    2.0
                )
                .is_err()
            );
            assert!(law.element_heat_flux(cancelled, &mesh, &[20.0; 4]).is_err());
        });
    });
}

#[test]
fn a_uniform_corner_in_any_principal_source_has_no_classical_adjoint() {
    let mesh = mesh();
    // Only principal direction one has a kink at f=1/2.
    let card = orthotropic_card([[2.0; 3], [1.0, 5.0, 2.0], [3.0; 3]]);
    let curve = caloric(&card);
    let sources = sources(&card, &curve);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let zero = ScalarField::Uniform(0.0);
    with_cx(|cx| {
        let law =
            MaterialPhaseTransport::try_new_orthotropic(cx, sources.each_ref(), AXES, 2.0).unwrap();
        let solver = law.prepare_uniform(cx, &mesh, budget()).unwrap();
        let problem = law.problem(&mesh, &boundary, &zero);
        let h = [50.0; 4];
        let accepted = solver.advance(cx, problem, None, &h, DT, config()).unwrap();
        assert_eq!(accepted.specific_enthalpy_j_kg, h);
        assert!(matches!(
            solver.linearize_accepted(cx, problem, None, &h, DT, config(), accepted),
            Err(EnthalpyAdjointError::MaterialKink)
        ));
    });
    let zero_card = orthotropic_card([[2.0; 3], [0.0; 3], [3.0; 3]]);
    let zero_curve = caloric(&zero_card);
    let zero_sources = self::sources(&zero_card, &zero_curve);
    with_cx(|cx| {
        assert!(
            MaterialPhaseTransport::try_new_orthotropic(cx, zero_sources.each_ref(), AXES, 2.0)
                .is_err()
        );
    });
}
