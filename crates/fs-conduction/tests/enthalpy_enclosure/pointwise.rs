//! Independent simplex moments and finite-body balances for reflected local emission.
use super::*;
use fs_conduction::{
    ConductionError, PointwiseEnclosureRadiation, PointwiseEnclosureTolerance, ViewFactorEvidence,
    ViewFactorTolerance,
};

fn boundary(fixture: &Fixture, mixed: bool) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(&fixture.mesh)
        .region(
            "left",
            |f| f.vertices.iter().all(|v| *v < 4),
            if mixed {
                ThermalBc::neumann(-25.0).unwrap()
            } else {
                ThermalBc::adiabatic()
            },
        )
        .unwrap()
        .region(
            "right",
            |f| f.vertices.iter().all(|v| *v >= 4),
            if mixed {
                ThermalBc::robin(HTC, 300.0).unwrap()
            } else {
                ThermalBc::adiabatic()
            },
        )
        .unwrap()
        .finish()
        .unwrap()
}

fn problem<'a>(
    fixture: &'a Fixture,
    boundary: &'a ThermalBoundary,
    source: &'a ScalarField,
) -> ConductionProblem<'a> {
    ConductionProblem {
        boundary,
        ..fixture.problem(source)
    }
}

// Exact barycentric moments: int(N_i T_h^4) = A/105 sum_{|alpha|=4}
// (alpha_i+1) t^alpha. This does not use the production quadrature rule.
fn moments(area: f64, t: [f64; 3]) -> [f64; 3] {
    let mut out = [0.0; 3];
    for a in 0..=4 {
        for b in 0..=4 - a {
            let c = 4 - a - b;
            let powers = [a, b, c];
            let monomial = t[0].powi(a) * t[1].powi(b) * t[2].powi(c);
            for i in 0..3 {
                out[i] += area / 105.0 * f64::from(powers[i] + 1) * monomial;
            }
        }
    }
    out
}

// Explicit two-gray-plate radiosity solution, with exact local emission moments.
fn oracle(fixture: &Fixture, temperature: &[f64]) -> (Vec<f64>, [f64; 2], [f64; 2]) {
    let mut emission = vec![0.0; 8];
    let mut weights = vec![0.0; 8];
    let mut areas = [0.0; 2];
    for face in fixture.mesh.boundary() {
        let vertices = face.vertices.map(|v| v as usize);
        let body = vertices[0] / 4;
        areas[body] += face.area;
        let values = moments(face.area, vertices.map(|v| temperature[v]));
        for i in 0..3 {
            emission[vertices[i]] += STEFAN_BOLTZMANN_W_M2_K4 * values[i];
            weights[vertices[i]] += face.area / 3.0;
        }
    }
    let mean: [f64; 2] = std::array::from_fn(|body| {
        emission[4 * body..4 * body + 4].iter().sum::<f64>() / areas[body]
    });
    let jl = (EPS[0] * mean[0] + (1.0 - EPS[0]) * EPS[1] * mean[1])
        / (EPS[0] + EPS[1] - EPS[0] * EPS[1]);
    let jr = EPS[1] * mean[1] + (1.0 - EPS[1]) * jl;
    let irradiation = [jr, jl];
    let loads: Vec<_> = (0..8)
        .map(|v| EPS[v / 4] * (emission[v] - weights[v] * irradiation[v / 4]))
        .collect();
    let heat = std::array::from_fn(|body| loads[4 * body..4 * body + 4].iter().sum());
    (loads, heat, mean)
}

#[test]
fn heterogeneous_pointwise_enclosure_exposes_local_defects_when_total_energy_cancels() {
    use fs_conduction::transient::enthalpy::heterogeneous::{
        HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial,
    };
    let fixture = Fixture::new(true);
    let boundary = boundary(&fixture, false);
    let source = ScalarField::Uniform(0.0);
    let old = [997.0, 997.0, 997.0, 997.0, 501.0, 501.0, 501.0, 501.0];
    let candidate = [1000.0, 1000.0, 1000.0, 1000.0, 500.0, 500.0, 500.0, 500.0];
    with_cx(|cx| {
        let radiation = PointwiseEnclosureRadiation::new(
            cx, &fixture.mesh, &boundary, &fixture.enclosure,
            PointwiseEnclosureTolerance::default(),
        ).unwrap();
        let materials = [
            ReferenceEnthalpyMaterial { curve: &fixture.curve, reference_density_kg_m3: 8.0 },
            ReferenceEnthalpyMaterial { curve: &fixture.curve, reference_density_kg_m3: 24.0 },
        ];
        let engine = HeterogeneousEnthalpyBackwardEuler::new(
            cx, &fixture.mesh, &materials, &[0, 1],
            EnthalpyBudget { max_vertices: 8, max_elements: 2 },
        ).unwrap();
        let got = engine.evaluate_endpoint_with_pointwise_enclosure_radiation(
            cx, problem(&fixture, &boundary, &source), None,
            &old, &candidate, DT, step_config(), &radiation,
        ).unwrap();
        // Equal tetrahedral volumes and rho=[8,24], dh=[3,-1] give zero
        // total storage. Each hot node nevertheless has +6V joules, each
        // cold node -6V, and reflected radiation increases these local defects.
        let volume = 2.0_f64.sqrt() / 12.0;
        let heat = 3.0_f64.sqrt() * STEFAN_BOLTZMANN_W_M2_K4
            * (350.0_f64.powi(4) - 300.0_f64.powi(4))
            / (1.0 / EPS[0] + 1.0 / EPS[1] - 1.0);
        let nodal_defect = 6.0 * volume + DT * heat / 4.0;
        close(got.residual_norm_j, 8.0_f64.sqrt() * nodal_defect, 1e-10);
        close(got.stored_energy_change_j, 0.0, 1e-12);
        close(got.radiation_out_w, 0.0, 1e-9);
        close(got.energy_residual_j, 0.0, 1e-9);
        assert!(got.residual_norm_j > 1.0);
        assert_eq!(got.free_dofs, 8);
        assert!(got.convective_robin_fluxes.is_empty());
        for &temperature in &got.temperature[..4] { close(temperature, 350.0, 1e-12); }
        for &temperature in &got.temperature[4..] { close(temperature, 300.0, 1e-12); }
    });
    assert_eq!(old, [997.0, 997.0, 997.0, 997.0, 501.0, 501.0, 501.0, 501.0]);
    assert_eq!(candidate, [1000.0, 1000.0, 1000.0, 1000.0, 500.0, 500.0, 500.0, 500.0]);
}

#[test]
fn local_hot_spots_drive_exact_emission_reflection_and_a_nonsymmetric_tangent() {
    let fixture = Fixture::new(false);
    let boundary = boundary(&fixture, false);
    let t = [300.0, 310.0, 500.0, 550.0, 300.0, 320.0, 340.0, 360.0];
    with_cx(|cx| {
        let radiation = PointwiseEnclosureRadiation::new(
            cx,
            &fixture.mesh,
            &boundary,
            &fixture.enclosure,
            PointwiseEnclosureTolerance::default(),
        )
        .unwrap();
        let got = radiation.assemble(cx, &t).unwrap();
        let (loads, heat, mean) = oracle(&fixture, &t);
        for v in 0..8 {
            close(got.outward_load_w[v], loads[v], 2e-13);
        }
        for body in 0..2 {
            close(got.report.patches[body].outward_heat_w, heat[body], 2e-13);
            close(
                got.report.patches[body].mean_blackbody_emission_w_m2,
                mean[body],
                2e-14,
            );
        }
        close(heat[0], -heat[1], 2e-13);
        close(
            got.outward_load_w.iter().sum(),
            got.report.outward_heat_w,
            2e-12,
        );
        let patch_means: [f64; 2] =
            std::array::from_fn(|body| t[4 * body..4 * body + 4].iter().sum::<f64>() / 4.0);
        let mean_only = 3.0_f64.sqrt()
            * STEFAN_BOLTZMANN_W_M2_K4
            * (patch_means[0].powi(4) - patch_means[1].powi(4))
            / (1.0 / EPS[0] + 1.0 / EPS[1] - 1.0);
        assert!(
            heat[0] > 1.1 * mean_only,
            "local hot spots must survive patch exchange"
        );

        let direction = [0.2, -0.4, 0.6, 0.3, -0.1, 0.7, -0.2, 0.5];
        let mut applied = vec![0.0; 8];
        got.tangent_apply(cx, &direction, &mut applied).unwrap();
        let delta = 1e-3;
        let plus = std::array::from_fn::<_, 8, _>(|v| t[v] + delta * direction[v]);
        let minus = std::array::from_fn::<_, 8, _>(|v| t[v] - delta * direction[v]);
        let plus = oracle(&fixture, &plus).0;
        let minus = oracle(&fixture, &minus).0;
        for v in 0..8 {
            close(applied[v], (plus[v] - minus[v]) / (2.0 * delta), 2e-8);
        }
        close(applied.iter().sum(), 0.0, 1e-12);

        let mut a = [0.0; 8];
        a[0] = 1.0;
        let mut b = [0.0; 8];
        b[4] = 1.0;
        let mut ja = [0.0; 8];
        let mut jb = [0.0; 8];
        got.tangent_apply(cx, &a, &mut ja).unwrap();
        got.tangent_apply(cx, &b, &mut jb).unwrap();
        assert!(
            (ja[4] - jb[0]).abs() > 0.1,
            "reflected local T^3 columns are not symmetric"
        );
        for column in 0..8 {
            let mut unit = [0.0; 8];
            unit[column] = 1.0;
            let mut action = [0.0; 8];
            got.tangent_apply(cx, &unit, &mut action).unwrap();
            close(got.temperature_diagonal_w_k[column], action[column], 1e-13);
        }
    });
}

#[test]
fn implicit_nonuniform_bodies_close_each_nodal_balance_and_keep_authored_heat_separate() {
    let fixture = Fixture::new(false);
    let source = ScalarField::Uniform(0.0);
    let t = [300.0, 310.0, 500.0, 550.0, 300.0, 320.0, 340.0, 360.0];
    let old: Vec<_> = t.iter().map(|t| CP * (t - 200.0)).collect();
    let history = old.clone();
    let volume = 2.0_f64.sqrt() / 12.0;
    let area = 3.0_f64.sqrt();
    with_cx(|cx| {
        for mixed in [false, true] {
            let boundary = boundary(&fixture, mixed);
            let radiation = PointwiseEnclosureRadiation::new(
                cx,
                &fixture.mesh,
                &boundary,
                &fixture.enclosure,
                PointwiseEnclosureTolerance::default(),
            )
            .unwrap();
            let got = fixture
                .engine(cx, false)
                .advance_with_pointwise_enclosure_radiation(
                    cx,
                    problem(&fixture, &boundary, &source),
                    None,
                    &old,
                    DT,
                    step_config(),
                    &radiation,
                )
                .unwrap();
            let (radiative, _, _) = oracle(&fixture, &got.temperature);
            for body in 0..2 {
                let sum: f64 = got.temperature[4 * body..4 * body + 4].iter().sum();
                // Independent exact regular-tetrahedron P1 stiffness with k(mean T).
                let conductivity = 2.0 + 6.0 * (sum / 4.0 - 240.0) / 320.0;
                for v in 4 * body..4 * body + 4 {
                    let conduction = conductivity * volume * (2.0 * got.temperature[v] - 0.5 * sum);
                    let external = if !mixed {
                        0.0
                    } else if body == 0 {
                        -25.0 * area / 4.0
                    } else {
                        HTC * area / 24.0 * (2.0 * got.temperature[v] + sum - 1800.0)
                    };
                    let residual = RHO * volume / 4.0 * (got.specific_enthalpy_j_kg[v] - old[v])
                        + DT * (conduction + radiative[v] + external);
                    close(residual, 0.0, 3e-8);
                }
            }
            let prescribed = if mixed { -25.0 * area } else { 0.0 };
            let convection = if mixed {
                HTC * area * (got.temperature[4..].iter().sum::<f64>() / 4.0 - 300.0)
            } else {
                0.0
            };
            close(got.neumann_out_w, prescribed, 1e-13);
            close(got.robin_out_w, convection, 1e-12);
            close(
                got.stored_energy_change_j,
                -DT * (prescribed + convection),
                1e-8,
            );
            close(got.energy_residual_j, 0.0, 1e-8);
            close(
                got.pointwise_enclosure_radiation
                    .as_ref()
                    .unwrap()
                    .outward_heat_w,
                0.0,
                1e-9,
            );
            assert!(got.pointwise_radiation.is_none());
            assert_eq!(old, history);
        }
    });
}

#[test]
fn finite_donor_pays_for_receiver_latent_heat_at_constant_temperature() {
    let fixture = Fixture::new(true);
    let boundary = boundary(&fixture, false);
    let source = ScalarField::Uniform(0.0);
    let old = [
        3500.0, 3500.0, 3500.0, 3500.0, 1500.0, 1500.0, 1500.0, 1500.0,
    ];
    let mass = RHO * 2.0_f64.sqrt() / 12.0;
    let factor = 3.0_f64.sqrt() * STEFAN_BOLTZMANN_W_M2_K4 / (1.0 / EPS[0] + 1.0 / EPS[1] - 1.0);
    let (mut lo, mut hi) = (350.0_f64, 400.0_f64);
    for _ in 0..90 {
        let t = (lo + hi) / 2.0;
        if mass * CP * (t - 400.0) + DT * factor * (t.powi(4) - 350.0_f64.powi(4)) > 0.0 {
            hi = t;
        } else {
            lo = t;
        }
    }
    let donor = (lo + hi) / 2.0;
    let heat = factor * (donor.powi(4) - 350.0_f64.powi(4));
    with_cx(|cx| {
        let radiation = PointwiseEnclosureRadiation::new(
            cx,
            &fixture.mesh,
            &boundary,
            &fixture.enclosure,
            PointwiseEnclosureTolerance::default(),
        )
        .unwrap();
        let got = fixture
            .engine(cx, true)
            .advance_with_pointwise_enclosure_radiation(
                cx,
                problem(&fixture, &boundary, &source),
                None,
                &old,
                DT,
                step_config(),
                &radiation,
            )
            .unwrap();
        for v in 0..4 {
            close(got.temperature[v], donor, 2e-10);
            close(
                got.specific_enthalpy_j_kg[v],
                3500.0 - DT * heat / mass,
                2e-10,
            );
            assert_eq!(got.temperature[v + 4], 350.0);
            close(
                got.specific_enthalpy_j_kg[v + 4],
                1500.0 + DT * heat / mass,
                2e-10,
            );
            close(
                got.liquid_mass_fraction[v + 4],
                (500.0 + DT * heat / mass) / 2000.0,
                2e-11,
            );
        }
        close(got.stored_energy_change_j, 0.0, 1e-8);
    });
}

#[test]
fn equilibrium_stays_exact_but_local_source_support_and_adjoint_ownership_are_enforced() {
    let fixture = Fixture::new(false);
    let boundary = boundary(&fixture, false);
    let source = ScalarField::Uniform(0.0);
    with_cx(|cx| {
        let radiation = PointwiseEnclosureRadiation::new(
            cx,
            &fixture.mesh,
            &boundary,
            &fixture.enclosure,
            PointwiseEnclosureTolerance::default(),
        )
        .unwrap();
        let uniform = radiation.assemble(cx, &[350.0; 8]).unwrap();
        assert_eq!(uniform.outward_load_w, vec![0.0; 8]);
        assert_eq!(uniform.report.outward_heat_w, 0.0);
        let mut unsupported = [350.0; 8];
        unsupported[0] = 600.0;
        assert!(matches!(
            radiation.assemble(cx, &unsupported),
            Err(ConductionError::MaterialQuery { .. })
        ));
        with_cancelled_cx(|cancelled| {
            assert!(matches!(
                radiation.assemble(cancelled, &[350.0; 8]),
                Err(ConductionError::Cancelled { .. })
            ))
        });
        let old = [1500.0; 8];
        let engine = fixture.engine(cx, false);
        let accepted = engine
            .advance_with_pointwise_enclosure_radiation(
                cx,
                problem(&fixture, &boundary, &source),
                None,
                &old,
                DT,
                step_config(),
                &radiation,
            )
            .unwrap();
        assert_eq!(accepted.specific_enthalpy_j_kg, old);
        assert!(
            engine
                .linearize_accepted(
                    cx,
                    problem(&fixture, &boundary, &source),
                    None,
                    &old,
                    DT,
                    step_config(),
                    accepted,
                )
                .is_err()
        );
    });
}

#[test]
fn strong_self_view_keeps_negative_exact_diagonal_without_breaking_positive_preconditioning() {
    let mut fixture = Fixture::new(false);
    fixture.curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([0x71; 32]),
        SolidLiquidPhase::Solid,
        [200.0, 3500.0]
            .map(|t| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: CP * (t - 200.0),
                temperature_k: t,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: RHO,
            })
            .to_vec(),
    )
    .unwrap();
    fixture.material =
        ConductivityModel::constant_tensor([[0.05, 0.0, 0.0], [0.0, 0.05, 0.0], [0.0, 0.0, 0.05]])
            .unwrap();
    fixture.boundary = ThermalBoundaryBuilder::new(&fixture.mesh)
        .region(
            "left",
            |f| f.vertices.iter().all(|v| *v < 3),
            ThermalBc::adiabatic(),
        )
        .unwrap()
        .region(
            "right",
            |f| f.vertices.iter().all(|v| (4..7).contains(v)),
            ThermalBc::adiabatic(),
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let left = RadiationSurface::new(
        &fixture.mesh,
        "left",
        |f| f.vertices.iter().all(|v| *v < 3),
        emissivity_in("left", 1.0, 240.0, 3500.0),
    )
    .unwrap();
    let right = RadiationSurface::new(
        &fixture.mesh,
        "right",
        |f| f.vertices.iter().all(|v| (4..7).contains(v)),
        emissivity_in("right", 1.0, 240.0, 3500.0),
    )
    .unwrap();
    let area = left.area_m2();
    let factors = ViewFactorMatrix::admit(
        vec![area; 2],
        vec![vec![0.99, 0.01], vec![0.01, 0.99]],
        ViewFactorEvidence::Analytic {
            geometry: "prescribed reciprocal two-patch self-view".into(),
        },
        ViewFactorTolerance::default(),
    )
    .unwrap();
    fixture.enclosure = GrayDiffuseEnclosure::new(vec![left, right], factors).unwrap();
    let t = [300.0, 300.0, 3000.0, 300.0, 300.0, 300.0, 300.0, 300.0];
    let old: Vec<_> = t.iter().map(|t| CP * (t - 200.0)).collect();
    let source = ScalarField::Uniform(0.0);
    with_cx(|cx| {
        let radiation = PointwiseEnclosureRadiation::new(
            cx,
            &fixture.mesh,
            &fixture.boundary,
            &fixture.enclosure,
            PointwiseEnclosureTolerance::default(),
        )
        .unwrap();
        let assembly = radiation.assemble(cx, &t).unwrap();
        assert!(assembly.temperature_diagonal_w_k[0] < 0.0);
        // At dt=10 s the full reflected diagonal would make the old Picard
        // preconditioner negative; the true nonsymmetric action is still valid.
        assert!(
            RHO * 2.0_f64.sqrt() / 48.0
                + 10.0 * (0.05 * 2.0_f64.sqrt() / 8.0 + assembly.temperature_diagonal_w_k[0]) / CP
                < 0.0
        );
        let got = fixture
            .engine(cx, false)
            .advance_with_pointwise_enclosure_radiation(
                cx,
                fixture.problem(&source),
                None,
                &old,
                10.0,
                step_config(),
                &radiation,
            )
            .unwrap();
        assert!(got.temperature[2] < 3000.0);
        assert!(got.temperature[0] > 300.0);
        close(got.stored_energy_change_j, 0.0, 1e-8);
        close(got.energy_residual_j, 0.0, 1e-8);
    });
}


#[path = "pointwise/adjoint.rs"]
mod adjoint;
