//! G1: physical re-solves of live contact, k(T), and ambient radiation together.
use super::*;
use fs_conduction::radiation::AmbientRadiationLinearization;

// Whole contact-law scale, source scale, absolute emissivity, reservoir kelvin.
fn sample(cx: &Cx<'_>, convection: bool, controls: [f64; 4], adjoint: bool)
    -> (f64, Option<[f64; 4]>)
{
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("exterior", |f| f.centroid[2] < 0.0 && f.centroid[1] == 0.0,
            if convection { ThermalBc::robin(2.0, 297.0).unwrap() }
            else { ThermalBc::neumann(0.0).unwrap() })
        .unwrap().adiabatic_remainder().finish().unwrap();
    let original = card(false);
    let mut claim = original.claims().claims_ordered().next().unwrap().1.clone();
    let PropertyValue::Curve { knots, .. } = &mut claim.value else { panic!("curve fixture"); };
    for (_, value) in knots { *value *= controls[0]; }
    let mut claims = ClaimSet::new();
    claims.insert_claim(claim).unwrap();
    let card = InterfaceSystemCard::assemble(original.surface_a().clone(),
        original.surface_b().clone(), original.context().clone(), claims, vec![]).unwrap();
    let contacts = contacts(cx, &mesh, &boundary, &card);
    let material = material(true);
    let source = ScalarField::Nodal(vec![0., 0., 0., 0.,
        30. * controls[1], 30. * controls[1], 30. * controls[1], 30. * controls[1]]);
    let patch = AmbientRadiationPatch::new("exterior",
        patch_with_emissivity(controls[2]).emissivity().clone(), controls[3]).unwrap();
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
        material: &material, element_materials: None, source: &source };
    let mut config = config();
    config.residual_tolerance_w = 1e-9;
    config.energy_tolerance_w = 1e-9;
    let solved = contacts.solve_steady_newton_with_ambient_radiation(cx, problem,
        &[300.; 8], config, std::slice::from_ref(&patch), AmbientRadiationConfig {
            max_iterations: 256, temperature_tolerance_k: 1e-11,
            balance_tolerance_w: 1e-9, balance_relative_tolerance: 1e-12,
            relaxation: 0.5,
        }).unwrap();
    let temperature = &solved.contact.temperature;
    let goal = temperature[7];
    if !adjoint { return (goal, None); }
    let before: Vec<_> = temperature.iter().map(|t| t.to_bits()).collect();
    let linear = step_config().linear;
    let bound = AmbientRadiationLinearization::bind(cx, problem, temperature,
        &["exterior"], &[0.], Some(&[0.]), std::slice::from_ref(&patch), 4096).unwrap();
    let physical = ConductionProblem { boundary: bound.boundary(), ..problem };
    let response = contacts.linearize_steady_endpoint(cx, physical, temperature,
        linear, 1e-9, 1e-9, &["exterior"], 1).unwrap();
    let full = response.pullback_with_robin_feedback_and_resistance_controls(cx,
        bound.coefficient_mean_slopes(), bound.local_reference_slopes(),
        bound.reference_feedback().unwrap(), &[0., 0., 0., 0., 0., 0., 0., 1.],
        4096, 1).unwrap();
    let resistance = full.resistance_multipliers[0].derivative;
    let mut invalid = full.transport.clone();
    invalid.references.clear();
    assert!(bound.remap_gradient(cx, invalid, linear).is_err());
    let mut unresolved = full.transport.clone();
    unresolved.relative_residual = 2.0 * linear.tolerance;
    assert!(bound.remap_gradient(cx, unresolved, linear).is_err());
    let mut invalid_residual = full.transport.clone();
    invalid_residual.relative_residual = -1.0;
    assert!(bound.remap_gradient(cx, invalid_residual, linear).is_err());
    let gradient = bound.remap_gradient(cx, full.transport, linear).unwrap();
    assert_eq!(before, temperature.iter().map(|t| t.to_bits()).collect::<Vec<_>>(),
        "a derivative must retain the exact original primal field");
    if !convection {
        assert_eq!(gradient.convection.references, [0.0]);
        assert_eq!(gradient.convection.log_htc, [0.0]);
    }
    // Constant density 30 on the hot tet integrates to 5 W: each P1 vertex
    // receives 1.25 W per unit source multiplier, independent of contact/radiation.
    let power = 1.25 * gradient.convection.nodal_load[4..].iter().sum::<f64>();
    let mut stale = temperature.to_vec();
    stale[7] += 0.1;
    assert!(contacts.linearize_steady_endpoint(cx, physical, &stale,
        linear, 1e-9, 1e-9, &["exterior"], 1).is_err());
    assert!(AmbientRadiationLinearization::bind(cx, problem, temperature,
        &["exterior"], &[0.], Some(&[0.]), &[patch], 0).is_err());
    (goal, Some([resistance, power, gradient.emissivities[0],
        gradient.reservoir_temperatures[0]]))
}

#[test]
fn live_contact_ambient_dual_matches_physical_resolves_with_and_without_convection() {
    context(&CancelGate::new_clock_free(), |cx| {
        let nominal = [1.0, 1.0, 0.2, 300.0];
        for convection in [false, true] {
            let (goal, derivatives) = sample(cx, convection, nominal, true);
            near(goal, sample(cx, convection, nominal, false).0, 0.0);
            let derivatives = derivatives.unwrap();
            for (control, step) in [0.001, 0.001, 0.0002, 0.002].into_iter().enumerate() {
                let mut low = nominal;
                let mut high = nominal;
                low[control] -= step;
                high[control] += step;
                let expected = (sample(cx, convection, high, false).0
                    - sample(cx, convection, low, false).0) / (2.0 * step);
                assert!((derivatives[control] - expected).abs()
                    < 5e-5 * expected.abs().max(1e-3),
                    "convection={convection}, control={control}: {} vs {expected}",
                    derivatives[control]);
            }
        }
    });
}
