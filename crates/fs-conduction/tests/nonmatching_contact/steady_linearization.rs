//! Independent steady heat-flow derivatives, including source feedback.
use super::*;
use super::steady::{contacts, coupling, reservoir_boundary, scalar_heat, solver_config};
use super::temperature::source_card;

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

fn contact_mean(mesh: &ConductionMesh, split: usize) -> Vec<f64> {
    let (a, _) = sides(mesh, split, 0.0);
    let area: f64 = a.iter().map(|&i| mesh.boundary()[i].area).sum();
    let mut weights = vec![0.0; mesh.vertex_count()];
    for face in a {
        let face = &mesh.boundary()[face];
        for vertex in face.vertices {
            weights[vertex as usize] += face.area / (3.0 * area);
        }
    }
    weights
}

#[test]
fn g3_steady_contact_adjoints_match_series_physics_and_fresh_coupled_solves() {
    with_cx(|cx| {
        for matching in [true, false] {
            let (mesh, split) = mesh(if matching { [1, 1, 1] } else { [1, 2, 3] }, 0.0);
            let bc = reservoir_boundary(&mesh, true);
            let card = source_card(1.0, false);
            let set = contacts(cx, &mesh, &bc, &card, split, matching);
            let material = ConductivityModel::isotropic_declared(10.0).unwrap();
            let source = ScalarField::Uniform(0.0);
            let problem = ConductionProblem {
                mesh: &mesh, boundary: &bc, material: &material,
                element_materials: None, source: &source,
            };
            let initial = vec![315.0; mesh.vertex_count()];
            let solved = set.solve_steady(cx, problem, &initial, solver_config(), coupling()).unwrap();
            let before = solved.conduction.temperature.clone();
            let response = set.linearize_steady_endpoint(
                cx, problem, &before, solver_config().linear,
                coupling().residual_tolerance_w, coupling().energy_tolerance_w,
                &["hot reservoir"], 10_000,
            ).unwrap();
            assert_eq!(response.temperature(), before);
            assert!(response.uses_nonlinear_jacobian());
            assert!(response.has_smooth_material_tangent());
            assert!(response.residual_norm_w() <= coupling().residual_tolerance_w);
            near(response.energy().closure_w, 0.0, coupling().energy_tolerance_w);

            // Independent uniform series equation, not a residual assembled
            // by the implementation under test. H is hot-to-cold heat in W.
            let heat = scalar_heat(true);
            let film = 0.1;
            let ta = 300.0 + 0.1 * heat;
            let tb = 330.0 - (0.1 + film) * heat;
            let r = 0.5 + (ta - 320.0) / 1024.0 + 3.0 * (tb - 320.0) / 1024.0;
            let denominator = 0.2 + film + r + heat * (-0.2 - 3.0 * film) / 1024.0;
            let dh_ds = -heat * r / denominator;
            let dh_dreference = (1.0 - 3.0 * heat / 1024.0) / denominator;
            let dh_dlogh = heat * film * (1.0 - 3.0 * heat / 1024.0) / denominator;
            let weights = contact_mean(&mesh, split);
            near(dot(&weights, response.temperature()), ta, 2e-7);
            near(response.dirichlet_reactions_w().iter().map(|(_, q)| q).sum(), -heat, 2e-7);

            let gradient = response.pullback_with_resistance_controls(
                cx, &weights, &[0.0], &[0.0], 10_000,
            ).unwrap();
            assert_eq!(gradient.resistance_multipliers.len(), 1);
            assert_eq!(gradient.resistance_multipliers[0].interface, "bond");
            near(gradient.resistance_multipliers[0].derivative, 0.1 * dh_ds, 2e-7);
            near(gradient.transport.references[0], 0.1 * dh_dreference, 2e-8);
            near(gradient.transport.log_htc[0], 0.1 * dh_dlogh, 2e-8);

            let mut direction = response.zero_direction();
            direction.references_k[0] = 1.0;
            let tangent = response.apply(cx, &direction).unwrap();
            near(dot(&weights, &tangent.temperature_k), 0.1 * dh_dreference, 2e-8);
            near(tangent.mean_wall_temperatures_k[0], 1.0 - film * dh_dreference, 2e-8);
            near(tangent.heat_rates_w[0], -dh_dreference, 2e-7);
            direction.references_k[0] = 0.0;
            direction.log_htc[0] = 1.0;
            let tangent = response.apply(cx, &direction).unwrap();
            near(dot(&weights, &tangent.temperature_k), 0.1 * dh_dlogh, 2e-8);
            near(tangent.mean_wall_temperatures_k[0], film * (heat - dh_dlogh), 2e-7);
            near(tangent.heat_rates_w[0], -dh_dlogh, 2e-7);

            // These are fresh source cards and nonlinear physical solves;
            // perturbing a frozen contact operator would miss the feedback.
            let objective = |scale| {
                let card = source_card(scale, false);
                let set = contacts(cx, &mesh, &bc, &card, split, matching);
                let next = set.solve_steady(cx, problem, &initial, solver_config(), coupling()).unwrap();
                dot(&weights, &next.conduction.temperature)
            };
            let h = 5e-4;
            let difference = (objective(1.0 + h) - objective(1.0 - h)) / (2.0 * h);
            near(gradient.resistance_multipliers[0].derivative, difference, 2e-5);

            // Mixed source/film/reference controls and heat observations also
            // exercise the true nonsymmetric transpose, including direct terms.
            direction.references_k[0] = 0.37;
            direction.log_htc[0] = -0.19;
            for (i, value) in direction.nodal_load_w.iter_mut().enumerate() {
                *value = (i % 7) as f64 / 11.0 - 0.2;
            }
            let tangent = response.apply(cx, &direction).unwrap();
            let gradient = response.pullback(cx, &weights, &[0.4], &[-0.1]).unwrap();
            let forward = dot(&weights, &tangent.temperature_k)
                + 0.4 * tangent.mean_wall_temperatures_k[0] - 0.1 * tangent.heat_rates_w[0];
            let reverse = dot(&gradient.nodal_load, &direction.nodal_load_w)
                + gradient.references[0] * direction.references_k[0]
                + gradient.log_htc[0] * direction.log_htc[0];
            near(forward, reverse, 2e-8);
            assert_eq!(solved.conduction.temperature, before);
        }
    });
}

#[test]
fn g0_steady_contact_derivatives_recheck_physical_endpoints_and_budgets() {
    with_cx(|cx| {
        let (mesh, split) = mesh([1, 2, 3], 0.0);
        let bc = boundary(&mesh);
        let card = source_card(1.0, false);
        let set = contacts(cx, &mesh, &bc, &card, split, false);
        let material = ConductivityModel::isotropic_declared(10.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem {
            mesh: &mesh, boundary: &bc, material: &material,
            element_materials: None, source: &source,
        };
        let initial = vec![315.0; mesh.vertex_count()];
        let solved = set.solve_steady(cx, problem, &initial, solver_config(), coupling()).unwrap();
        let good = solved.conduction.temperature;
        let dofs = DofMap::new(&bc, mesh.vertex_count()).unwrap();
        let linearize = |temperature: &[f64], residual, cap| {
            set.linearize_steady_endpoint(cx, problem, temperature, solver_config().linear,
                residual, coupling().energy_tolerance_w, &[], cap)
        };
        assert!(linearize(&good, 1e-8, 10_000).is_ok());
        let mut wrong_free = good.clone();
        wrong_free[dofs.free()[0]] += 1e-3;
        assert!(linearize(&wrong_free, 1e-8, 10_000).is_err());
        let mut wrong_prescribed = good.clone();
        let fixed = dofs.fixed()[0];
        wrong_prescribed[fixed] = wrong_prescribed[fixed].next_up();
        assert!(linearize(&wrong_prescribed, 1e-8, 10_000).is_err());
        assert!(linearize(&good, 0.0, 10_000).is_err());
        assert!(linearize(&good, 1e-8, 0).is_err());
        assert!(linearize(&vec![400.0; good.len()], 1e-8, 10_000).is_err());
        let checked = linearize(&good, 1e-8, 10_000).unwrap();
        assert!(checked.pullback_with_resistance_controls(
            cx, &vec![0.0; good.len()], &[], &[], 0,
        ).is_err());
        assert_eq!(checked.temperature(), good);
    });
}

#[test]
fn g0_insulated_steady_contact_equilibrium_cannot_bind_a_singular_response() {
    with_cx(|cx| {
        let (mesh, split) = mesh([1, 2, 3], 0.0);
        let bc = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let card = source_card(1.0, false);
        let set = contacts(cx, &mesh, &bc, &card, split, false);
        let material = ConductivityModel::isotropic_declared(10.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem {
            mesh: &mesh, boundary: &bc, material: &material,
            element_materials: None, source: &source,
        };
        // The primal residual vanishes for every constant temperature, yet
        // absolute steady temperatures and their derivatives are undetermined.
        assert!(matches!(
            set.linearize_steady_endpoint(cx, problem, &vec![320.0; mesh.vertex_count()],
                solver_config().linear, 1e-8, 1e-8, &[], 10_000),
            Err(ConductionError::SingularPureNeumann)
        ));
    });
}

#[test]
fn g3_steady_contact_and_robin_feedback_share_one_complete_dual() {
    with_cx(|cx| {
        for matching in [true, false] {
            let (mesh, split) = mesh(if matching { [1, 1, 1] } else { [1, 2, 3] }, 0.0);
            let bc = reservoir_boundary(&mesh, true);
            let card = source_card(1.0, false);
            let set = contacts(cx, &mesh, &bc, &card, split, matching);
            let material = ConductivityModel::isotropic_declared(10.0).unwrap();
            let source = ScalarField::Uniform(0.0);
            let problem = ConductionProblem { mesh: &mesh, boundary: &bc,
                material: &material, element_materials: None, source: &source };
            let solved = set.solve_steady(cx, problem, &vec![315.0; mesh.vertex_count()],
                solver_config(), coupling()).unwrap();
            let before = solved.conduction.temperature.clone();
            let response = set.linearize_steady_endpoint(cx, problem, &before,
                solver_config().linear, coupling().residual_tolerance_w,
                coupling().energy_tolerance_w, &["hot reservoir"], 10_000).unwrap();
            let weights = contact_mean(&mesh, split);
            let cap = 4 * mesh.vertex_count();
            let full = response.pullback_with_robin_feedback_and_resistance_controls(cx,
                &[0.05], &[0.1], &[0.2], &weights, cap, 10_000).unwrap();

            // Independent series physics for a physical reservoir with
            // h(W)=10+0.05(W-W0), ref(W)=330+0.3(W-W0). Splitting the
            // reference slope into local0.1 + additional0.2 must add once.
            // Eliminate the two bulk drops and the original linear R(TA,TB)
            // to obtain W(H,s), then solve only the scalar film balance.
            let wall0 = 330.0 - 0.1 * scalar_heat(true);
            let objective = |scale: f64, reference_offset: f64, log_h: f64| {
                let (mut lo, mut hi) = (0.0, 60.0);
                for _ in 0..80 {
                    let heat = 0.5 * (lo + hi);
                    let wall = (300.0 + 0.2 * heat
                        + heat * scale * (0.5 - 980.0 / 1024.0 - 0.2 * heat / 1024.0))
                        / (1.0 - 3.0 * heat * scale / 1024.0);
                    let h = log_h.exp() * (10.0 + 0.05 * (wall - wall0));
                    let reference = 330.0 + 0.3 * (wall - wall0) + reference_offset;
                    if heat > h * (reference - wall) { hi = heat; } else { lo = heat; }
                }
                300.0 + 0.1 * 0.5 * (lo + hi)
            };
            near(objective(1.0, 0.0, 0.0), dot(&weights, &before), 2e-7);
            let step = 5e-4;
            for (actual, high, low) in [
                (full.resistance_multipliers[0].derivative,
                    objective(1.0 + step, 0.0, 0.0), objective(1.0 - step, 0.0, 0.0)),
                (full.transport.references[0],
                    objective(1.0, step, 0.0), objective(1.0, -step, 0.0)),
                (full.transport.log_htc[0],
                    objective(1.0, 0.0, step), objective(1.0, 0.0, -step)),
            ] { near(actual, (high - low) / (2.0 * step), 2e-6); }
            let frozen = response.pullback_with_resistance_controls(cx, &weights, &[0.0], &[0.0], 10_000).unwrap();
            assert!((full.resistance_multipliers[0].derivative
                - frozen.resistance_multipliers[0].derivative).abs() > 0.05,
                "a contact-complete but boundary-frozen dual cannot pass this fixture");
            assert!(response.pullback_with_robin_feedback_and_resistance_controls(cx,
                &[0.05], &[0.1], &[0.2], &vec![0.0; weights.len()], cap - 1, 10_000).is_err());
            assert!(response.pullback_with_robin_feedback_and_resistance_controls(cx,
                &[0.05], &[0.1], &[0.2], &weights, cap, 0).is_err());
            assert_eq!(response.temperature(), before);
            assert_eq!(solved.conduction.temperature, before);
        }
    });
}
