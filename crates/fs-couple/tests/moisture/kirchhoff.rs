//! Source-integrated Fick fluxes: dry uptake, nonlinear steady flow and reservoirs.
//! Authored D(M) verifies the implementation; this is not measured wood evidence.
use super::*;

fn nonlinear() -> MoistureNonlinearConfig {
    MoistureNonlinearConfig {
        relative_tolerance: 1e-12,
        maximum_iterations: 128,
        relaxation: 1.0,
    }
}

fn state(initial: Vec<f64>, cx: &Cx<'_>) -> MoistureSlab {
    MoistureSlab::new_kirchhoff(
        card([0.0, 0.03], 0.3),
        point(0.1),
        "moisture_diffusivity",
        MaterialPropertySelection::SingleClaimOnly,
        geometry(1.0),
        initial,
        nonlinear(),
        cx,
    )
    .unwrap()
}

#[test]
fn g1_dry_degenerate_slab_wets_with_the_exact_nonlinear_face_potential() {
    with_cx(|cx| {
        let old = state(vec![0.0], cx);
        let (next, report) = old.advance(10.0, &faces(0.2), config(), cx).unwrap();
        // D(M)=0.1 M; Phi(M)=0.05 M². Two half-cell faces give
        // M=40*(Phi(.2)-Phi(M)), hence 2M²+M-.08=0.
        let expected = 0.16 / (1.0 + 1.64_f64.sqrt());
        near(next.moisture()[0], expected, 2e-12);
        assert!(expected > 0.0 && expected < 0.2);
        for transfer in report.boundary_water_kg {
            near(transfer, 0.25 * expected, 2e-12);
        }
        near(report.water_balance_residual_kg, 0.0, 2e-12);
        assert!(report.nonlinear.unwrap().relative_residual < nonlinear().relative_tolerance);
        assert!(report.kirchhoff_source.is_some());
        near(
            report.transport_start[0]
                .property("moisture_diffusivity")
                .unwrap()
                .value_si(),
            0.0,
            0.0,
        );
        near(
            report.transport_end[0]
                .property("moisture_diffusivity")
                .unwrap()
                .value_si(),
            0.1 * expected,
            2e-13,
        );
        assert_eq!(old.moisture(), &[0.0]);
        assert_eq!(old.time_s(), 0.0);
        assert_eq!(
            next.moisture(),
            old.advance(10.0, &faces(0.2), config(), cx)
                .unwrap()
                .0
                .moisture()
        );
    });
}

#[test]
fn g1_sqrt_profile_has_constant_physical_flux_and_remains_stationary() {
    with_cx(|cx| {
        let n = 16;
        let initial = (0..n)
            .map(|i| (0.04 * (i as f64 + 0.5) / n as f64).sqrt())
            .collect::<Vec<_>>();
        let old = state(initial.clone(), cx);
        let boundary = [
            MoistureFace::PrescribedDryBasis(0.0),
            MoistureFace::PrescribedDryBasis(0.2),
        ];
        let (next, report) = old.advance(2.0, &boundary, config(), cx).unwrap();
        // Phi is linear in x. Cell-centered differences and both physical
        // half-cell resistances therefore reproduce this nonlinear steady flow.
        for (&actual, &expected) in next.moisture().iter().zip(&initial) {
            near(actual, expected, 3e-12);
        }
        near(report.boundary_water_kg[0], -0.002, 3e-12);
        near(report.boundary_water_kg[1], 0.002, 3e-12);
        near(next.water_mass_kg(), old.water_mass_kg(), 3e-12);
        assert_eq!(report.nonlinear.unwrap().iterations, 1);
    });
}

#[test]
fn g1_finite_emission_solves_the_nonlinear_half_cell_during_uptake_and_desorption() {
    with_cx(|cx| {
        let dt = 10.0;
        let h = 0.004;
        let uptake = 0.16 / (11.0 + 21.64_f64.sqrt());
        let uptake_surface = 0.2 - uptake / 0.08;
        for (initial, expected, expected_surface, equilibrium) in [
            (0.0, uptake, uptake_surface, 0.2),
            // Manufactured endpoint: q/rho_dry = .1*(.11²-.12²)=-.00023.
            // Two faces over dt=10 remove .0046 kg water/kg dry solid.
            (0.1246, 0.12, 0.11, 0.0525),
        ] {
            let old = state(vec![initial], cx);
            let face = MoistureFace::SourceSorption {
                law: surface::law([h; 2]),
                relative_humidity: equilibrium / 0.25,
            };
            let (next, report) = old
                .advance(dt, &[face.clone(), face], config(), cx)
                .unwrap();
            near(next.moisture()[0], expected, 3e-12);
            for side in 0..2 {
                let actual_surface = report.surface_moisture[side].unwrap();
                near(actual_surface, expected_surface, 3e-12);
                let surface_transfer = old.dry_mass_kg() * h * dt * (equilibrium - actual_surface);
                let bulk_transfer = old.dry_mass_kg()
                    * 0.1
                    * dt
                    * (actual_surface * actual_surface - next.moisture()[0].powi(2));
                near(report.boundary_water_kg[side], surface_transfer, 3e-13);
                near(report.boundary_water_kg[side], bulk_transfer, 3e-13);
                assert!(report.surface_transfer[side].is_some());
            }
            near(report.water_balance_residual_kg, 0.0, 3e-12);
        }
    });
}

#[test]
fn g3_dry_uptake_depletes_the_original_finite_vapor_inventory() {
    with_cx(|cx| {
        let old = state(vec![0.0], cx);
        let air = atmosphere::atmosphere(0.8, 1.0);
        let law = atmosphere::sorption(0.25);
        let dt = 2.0;
        let candidate = air
            .advance(
                &old,
                dt,
                &[Some(law.clone()), Some(law)],
                config(),
                atmosphere::controls(),
                cx,
            )
            .unwrap();
        let capacity = air.vapor_mass_kg().unwrap() / air.relative_humidity().unwrap();
        let a = 0.25 * air.relative_humidity().unwrap();
        let b = 0.25 * old.dry_mass_kg() / capacity;
        // Eliminate RH through C*RH + m_dry*M = C*RH0, then use
        // M=2*dt*.1*((.25*RH)²-M²), independently of the slab solver.
        let qa = 2.0 * dt * 0.1 * (b * b - 1.0);
        let qb = -4.0 * dt * 0.1 * a * b - 1.0;
        let qc = 2.0 * dt * 0.1 * a * a;
        let expected = 2.0 * qc / (-qb + (qb * qb - 4.0 * qa * qc).sqrt());
        near(candidate.moisture().moisture()[0], expected, 2e-11);
        near(
            candidate.atmosphere().relative_humidity().unwrap(),
            0.8 - old.dry_mass_kg() * expected / capacity,
            2e-10,
        );
        near(
            candidate.atmosphere().vapor_mass_kg().unwrap() + candidate.moisture().water_mass_kg(),
            air.vapor_mass_kg().unwrap(),
            2e-13,
        );
        assert!(candidate.step().transport.nonlinear.is_some());
        near(
            atmosphere::carrier_mass(candidate.atmosphere()),
            atmosphere::carrier_mass(&air),
            2e-13,
        );
    });
}

#[test]
fn g0_source_path_and_nonlinear_work_refusals_publish_no_partial_wetting() {
    with_cx(|cx| {
        let old = state(vec![0.0], cx)
            .with_kirchhoff_diffusivity(MoistureNonlinearConfig {
                maximum_iterations: 1,
                ..nonlinear()
            })
            .unwrap();
        assert!(matches!(
            old.advance(10.0, &faces(0.2), config(), cx),
            Err(MoistureError::NonlinearNotConverged(_))
        ));
        assert_eq!(old.moisture(), &[0.0]);
        assert_eq!(old.time_s(), 0.0);
        let base = card([0.01; 2], 0.3);
        let mut claims = base.claims().clone();
        let mut claim = base.claims_for("moisture_diffusivity")[0].1.clone();
        claim.key = PropertyKey::new("unsupported_path_diffusivity", MOISTURE_DIFFUSIVITY_DIMS);
        claim.value = PropertyValue::Curve {
            abscissa: "moisture_content".into(),
            abscissa_dims: Dims::NONE,
            knots: vec![
                (0.0, 0.01),
                (0.12, 0.01),
                (0.15, -0.01),
                (0.18, 0.01),
                (0.3, 0.01),
            ],
            dims: MOISTURE_DIFFUSIVITY_DIMS,
        };
        claims.insert_claim(claim).unwrap();
        let material =
            MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap();
        let field = MoistureSlab::new_kirchhoff(
            material,
            point(0.1),
            "unsupported_path_diffusivity",
            MaterialPropertySelection::SingleClaimOnly,
            geometry(1.0),
            vec![0.05],
            nonlinear(),
            cx,
        )
        .unwrap();
        // Both sampled endpoints are positive; an intervening negative source
        // segment must still refuse before any physical state is published.
        assert!(matches!(
            field.advance(1.0, &faces(0.25), config(), cx),
            Err(MoistureError::Material(_))
        ));
        assert_eq!(field.moisture(), &[0.05]);
        assert_eq!(field.time_s(), 0.0);
    });
}
