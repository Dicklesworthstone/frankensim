//! The same physical enclosure tangent with independent convection controls.
//! The existing planar two-patch box fixture is an algebraic closure only.
use super::*;

fn physical_flux(
    cx: &fs_exec::Cx<'_>, c: Controls, all_flux: bool,
    fixed: Option<(usize, f64)>,
) -> Result<State, ConductionError> {
    let (complex, positions) = box_grid([3, 2, 2], [1.0; 3]);
    let mesh = ConductionMesh::new(complex, positions).unwrap();
    let mut boundary = ThermalBoundaryBuilder::new(&mesh);
    if let Some((axis, temperature)) = fixed {
        boundary = boundary.region("fixed", |face| on_box_face(face.centroid[axis], 0.0),
            ThermalBc::dirichlet(temperature).unwrap())?;
    }
    let boundary = boundary
        .region("lower", |f| on_box_face(f.centroid[0], 1.0) && f.centroid[1] < 0.5,
            ThermalBc::neumann(-12.0).unwrap())?
        .region("upper", |f| on_box_face(f.centroid[0], 1.0) && f.centroid[1] >= 0.5,
            if all_flux { ThermalBc::neumann(2.0).unwrap() }
            else { ThermalBc::robin(7.0, 310.0).unwrap() })?
        .adiabatic_remainder().finish()?;
    let lower = RadiationSurface::new(
        &mesh, "lower", |f| on_box_face(f.centroid[0], 1.0) && f.centroid[1] < 0.5,
        emissivity("lower", c.epsilon, "manufactured-lower"),
    ).unwrap();
    let upper = RadiationSurface::new(
        &mesh, "upper", |f| on_box_face(f.centroid[0], 1.0) && f.centroid[1] >= 0.5,
        emissivity("upper", c.second_epsilon, "manufactured-upper"),
    ).unwrap();
    let factors = ViewFactorMatrix::infinite_parallel_plates(lower.area_m2()).unwrap();
    let enclosure = GrayDiffuseEnclosure::new(vec![lower, upper], factors).unwrap();
    let material = ConductivityModel::isotropic(ConductivityTable::declared_curve(
        vec![(250.0, 6.0), (500.0, 16.0)],
    ).unwrap());
    let source = ScalarField::Nodal(mesh.positions().iter()
        .map(|p| c.density + 500.0 * p[1] + c.direction * direction(p)).collect());
    let mut cfg = config();
    cfg.stop.residual_rtol = 1e-12;
    cfg.stop.max_iterations = 80;
    let solved = solve_with_convective_gray_diffuse_enclosure(
        cx, ConductionProblem {
            mesh: &mesh, boundary: &boundary, material: &material,
            element_materials: None, source: &source,
        }, None, &enclosure, cfg, CoupledRadiationConfig {
            max_iterations: 256, surface_temperature_atol_k: 1e-10,
            surface_temperature_rtol: 0.0, relaxation: 0.5,
        },
    )?;
    Ok(State { mesh, boundary, material, source, enclosure, solved })
}

fn check_heat(cx: &fs_exec::Cx<'_>, state: &State, prescribed: f64) {
    assert!((state.solved.prescribed_neumann_out_w - prescribed).abs() < 1e-12);
    let means: Vec<_> = state.enclosure.surfaces().iter().map(|surface|
        surface.mean_temperature(&state.mesh, &state.solved.conduction.temperature).unwrap()).collect();
    let fresh = state.enclosure.solve(cx, &means).unwrap();
    let a = means[0];
    let b = means[1];
    let expected = STEFAN_BOLTZMANN_W_M2_K4 * 0.5 * (a - b) * (a + b) * (a*a + b*b)
        / (1.0 / nominal().epsilon + 1.0 / nominal().second_epsilon - 1.0);
    assert!(expected.abs() > 1e-3, "the two patches must exchange nonzero internal heat");
    assert!((fresh.net_outward_heat_w[0] - expected).abs() < 1e-8);
    assert!((fresh.net_outward_heat_w[1] + expected).abs() < 1e-8);
    let energy = &state.solved.conduction.report.energy;
    let applied = state.solved.radiation.radiosity.enclosure_energy_closure_w;
    let decomposition = energy.neumann_out_w + energy.robin_out_w
        - state.solved.prescribed_neumann_out_w - state.solved.convective_out_w - applied;
    assert!(decomposition.abs() < 1e-7, "applied boundary decomposition: {decomposition:e}");
    let physical_balance = energy.source_w + energy.dirichlet_in_w
        - state.solved.prescribed_neumann_out_w - state.solved.convective_out_w
        - fresh.enclosure_energy_closure_w;
    assert!(physical_balance.abs() < 1e-6, "original external heat balance: {physical_balance:e}");
}

#[test]
fn g1_mixed_neumann_robin_enclosure_dual_keeps_source_and_reflected_emissivity_controls() {
    with_cx(|cx| {
        let base = physical_flux(cx, nominal(), false, Some((0, 400.0))).unwrap();
        check_heat(cx, &base, -6.0);
        assert_eq!(base.solved.convective_robin_fluxes.len(), 1);
        assert_eq!(base.solved.convective_robin_fluxes[0].region, "upper");
        let before = base.solved.clone();
        let full = pullback_gray_diffuse_enclosure(
            cx, base.problem(), None, linear(), &base.solved.conduction.temperature,
            &["upper"], &base.enclosure, 1e-7, &base.weights(), 16 * base.mesh.vertex_count(),
        ).unwrap();
        assert_eq!(base.solved, before);
        assert_eq!(full.regions, ["upper"]);
        assert_eq!(full.radiation_regions, ["lower", "upper"]);
        assert_eq!(full.convection.log_htc.len(), 1);
        let analytic = [
            uniform_density_bar(&base, &full.convection.nodal_load),
            full.emissivities[0], full.emissivities[1],
        ];
        assert!(analytic[1].abs() > 1e-8 && analytic[2].abs() > 1e-8);
        for (index, step) in [0.25, 0.0001, 0.0001].into_iter().enumerate() {
            let probe = |sign: f64| {
                let mut controls = nominal();
                match index {
                    0 => controls.density += sign * step,
                    1 => controls.epsilon += sign * step,
                    _ => controls.second_epsilon += sign * step,
                }
                physical_flux(cx, controls, false, Some((0, 400.0))).unwrap().goal()
            };
            close(analytic[index], (probe(1.0) - probe(-1.0)) / (2.0 * step));
        }
        assert!(full.convection.relative_residual <= linear().tolerance);
    });
}

#[test]
fn g1_all_flux_enclosure_prescribed_lift_includes_fixed_radiating_vertices_without_robin_rows() {
    with_cx(|cx| {
        // The y=0 bath intersects the lower radiating trace. Its prescribed
        // vertices contribute to that patch mean and the complete J_FD block.
        let base = physical_flux(cx, nominal(), true, Some((1, 400.0))).unwrap();
        check_heat(cx, &base, -5.0);
        assert!(base.solved.convective_robin_fluxes.is_empty());
        assert_eq!(base.solved.convective_out_w, 0.0);
        let mut weights = base.weights();
        let fixed_goal = base.mesh.positions().iter().position(|p|
            on_box_face(p[0], 1.0) && on_box_face(p[1], 0.0)).unwrap();
        weights[fixed_goal] += 0.7;
        let cap = 16 * base.mesh.vertex_count();
        let full = pullback_gray_diffuse_enclosure(
            cx, base.problem(), None, linear(), &base.solved.conduction.temperature,
            &[], &base.enclosure, 1e-7, &weights, cap,
        ).unwrap();
        assert!(full.regions.is_empty());
        assert!(full.convection.references.is_empty());
        assert!(full.convection.log_htc.is_empty());
        assert_eq!(full.emissivities.len(), 2);
        let bound = EnclosureRadiationLinearization::bind(
            cx, base.problem(), linear(), &base.solved.conduction.temperature,
            &[], &[], &[], &base.enclosure, 1e-7, cap,
        ).unwrap();
        assert!(bound.flux_feedback().is_some());
        let combined = ConductionProblem { boundary: bound.boundary(), ..base.problem() };
        let lift = |flux| {
            let bars = RobinResponse::prescribed_temperature_pullback_with_flux_feedback(
                cx, combined, None, &base.solved.conduction.temperature, &weights,
                &full.convection.nodal_load, &[], bound.convective_mean_slopes(),
                bound.local_reference_slopes(), bound.reference_feedback(), flux, 65536,
            ).unwrap();
            base.boundary.dirichlet().iter().map(|&(vertex, _)| bars[vertex]).sum::<f64>()
        };
        let actual = lift(bound.flux_feedback());
        let step = 0.01;
        let probe = |temperature: f64| {
            let state = physical_flux(cx, nominal(), true, Some((1, temperature))).unwrap();
            weights.iter().zip(&state.solved.conduction.temperature)
                .map(|(&weight, &value)| weight * value).sum::<f64>()
        };
        close(actual, (probe(400.0 + step) - probe(400.0 - step)) / (2.0 * step));
        assert!((actual - lift(None)).abs() > 1e-7,
            "prescribed vertices must carry the Neumann radiation-mean feedback");
        assert!(actual > 0.7, "the direct fixed-node goal term is part of the lift");
        for &(vertex, _) in base.boundary.dirichlet() {
            assert_eq!(full.convection.nodal_load[vertex], 0.0);
        }
    });
}

#[test]
fn g0_closed_enclosure_fluxes_do_not_supply_a_steady_absolute_temperature_anchor() {
    with_cx(|cx| {
        assert!(physical_flux(cx, nominal(), true, None).is_err(),
            "an isolated all-flux steady enclosure must retain the anchor refusal");
    });
}

