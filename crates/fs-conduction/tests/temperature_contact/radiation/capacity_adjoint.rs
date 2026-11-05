//! Full physical histories through live contact, capacity and ambient radiation.
use super::*;
use fs_conduction::radiation::AmbientRadiationLinearization;

fn scaled_card(scale: f64, surface: bool) -> InterfaceSystemCard {
    let original = if surface { surface_card() } else { card(false) };
    let mut claim = original.claims().claims_ordered().next().unwrap().1.clone();
    match &mut claim.value {
        PropertyValue::Curve { knots, .. } => {
            for (_, value) in knots { *value *= scale; }
        }
        PropertyValue::Surface2D { samples, .. } => {
            for (_, value) in samples { *value *= scale; }
        }
        _ => panic!("original temperature source"),
    }
    let mut claims = ClaimSet::new();
    claims.insert_claim(claim).unwrap();
    InterfaceSystemCard::assemble(original.surface_a().clone(), original.surface_b().clone(),
        original.context().clone(), claims, vec![]).unwrap()
}

// Whole contact-law scale, source scale, absolute emissivity, bath kelvin,
// total capacity scale and uniform initial-temperature offset in kelvin.
fn trajectory(cx: &Cx<'_>, surface: bool, selected: usize, control: [f64; 6], reverse: bool)
    -> (f64, Option<[f64; 6]>)
{
    let mesh = mesh();
    let boundary = capacity_boundary(&mesh, surface, None);
    let card = scaled_card(control[0], surface);
    let contacts = contacts(cx, &mesh, &boundary, &card);
    let material = material(true);
    let engine = BackwardEuler::per_element(cx, &mesh, &[10.0, 20.0]
        .map(|capacity| VolumetricHeatCapacity::declared(capacity * control[4]).unwrap())).unwrap();
    let patch = AmbientRadiationPatch::new("exterior",
        patch_with_emissivity(control[2]).emissivity().clone(), control[3]).unwrap();
    let patches = std::slice::from_ref(&patch);
    let config = StepConfig { energy_tolerance_j: 1e-9, ..step_config() };
    let coupling = ContactIterationConfig { residual_tolerance_j: 1e-9, ..coupling() };
    let radiation_config = AmbientRadiationConfig { max_iterations: 256,
        temperature_tolerance_k: 1e-12, balance_tolerance_w: 1e-9,
        ..AmbientRadiationConfig::default() };
    let steps = [(0.1, 1.0), (0.25, 0.3), (0.15, 0.0)];
    assert!(selected > 0 && selected <= steps.len());
    let source = |power: f64, scale: f64| ScalarField::Nodal(vec![0., 0., 0., 0.,
        30. * power * scale, 30. * power * scale, 30. * power * scale, 30. * power * scale]);
    let mut history = vec![old().into_iter().map(|value| value + control[5]).collect::<Vec<_>>()];
    let mut targets = Vec::new();
    for (dt, power) in steps {
        let source = source(power, control[1]);
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
            material: &material, element_materials: None, source: &source };
        let previous = history.last().unwrap();
        let before = previous.clone();
        let solved = if surface {
            contacts.advance_capacity_with_ambient_radiation(cx, &engine, problem,
                previous, dt, config, Some(nonlinear()), coupling, patches, radiation_config).unwrap()
        } else {
            contacts.advance_capacity_newton_with_ambient_radiation(cx, &engine, problem,
                previous, dt, config, nonlinear(), coupling, 1, patches, radiation_config).unwrap()
        };
        assert_eq!(&before, previous, "all trials must retain the original old field");
        assert!(solved.physical.residual_norm_j <= solved.physical_residual_tolerance_j);
        assert!(solved.physical.energy_residual_j.abs() <= config.energy_tolerance_j);
        targets.push(solved.physical_residual_tolerance_j);
        history.push(solved.conduction.temperature);
    }
    let goal = history[selected][7];
    if !reverse { return (goal, None); }
    let before = history.clone();
    let mut carry = vec![0.0; 8];
    carry[7] = 1.0;
    let mut derivatives = [0.0; 6];
    for index in (0..selected).rev() {
        let (dt, power) = steps[index];
        let source = source(power, control[1]);
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
            material: &material, element_materials: None, source: &source };
        let temperature = &history[index + 1];
        let radiation = AmbientRadiationLinearization::bind(cx, problem, temperature,
            &["exterior"], &[0.], Some(&[0.]), patches, 4096).unwrap();
        let physical = ConductionProblem { boundary: radiation.boundary(), ..problem };
        let response = engine.linearize_temperature_contact_endpoint(cx, physical,
            &contacts, &history[index], temperature, dt, config, targets[index], &["exterior"], 1).unwrap();
        assert_eq!(response.primal().temperature, *temperature);
        assert_eq!(response.primal().krylov_iterations, 0, "reverse must not replace the primal");
        let full = response.pullback_with_robin_and_flux_feedback_and_resistance_controls(cx,
            radiation.coefficient_mean_slopes(), radiation.local_reference_slopes(),
            radiation.reference_feedback().unwrap(), None, &carry, 4096, 1).unwrap();
        derivatives[0] += full.resistance_multipliers[0].derivative;
        let full = radiation.remap_gradient(cx, full.transport, config.linear).unwrap();
        if !surface {
            assert_eq!(full.convection.references, [0.0]);
            assert_eq!(full.convection.log_htc, [0.0]);
        }
        derivatives[1] += response.source_multiplier_pullback(cx, &full.convection.nodal_load).unwrap();
        derivatives[2] += full.emissivities[0];
        derivatives[3] += full.reservoir_temperatures[0];
        derivatives[4] += response.capacity_multiplier_pullback(cx, &full.convection.nodal_load).unwrap();
        carry = response.previous_temperature_pullback(cx, &full.convection.nodal_load).unwrap();
    }
    derivatives[5] = carry.iter().sum();
    assert_eq!(history, before, "the complete accepted history stays unchanged");
    (goal, Some(derivatives))
}

#[test]
fn capacity_ambient_contact_duals_match_full_and_selected_prefix_histories() {
    context(&CancelGate::new_clock_free(), |cx| {
        let nominal = [1.0, 1.0, 0.2, 295.0, 1.0, 0.0];
        // Vacuum curve/Newton uses the final field; Robin source-surface/Picard
        // uses an earlier retained state, as the sampled-peak reverse does.
        for (surface, selected) in [(false, 3), (true, 2)] {
            let (goal, derivatives) = trajectory(cx, surface, selected, nominal, true);
            near(goal, trajectory(cx, surface, selected, nominal, false).0, 0.0);
            let derivatives = derivatives.unwrap();
            for (coordinate, step) in [0.001, 0.001, 0.0002, 0.002, 0.001, 0.002].into_iter().enumerate() {
                let mut low = nominal;
                let mut high = nominal;
                low[coordinate] -= step;
                high[coordinate] += step;
                let expected = (trajectory(cx, surface, selected, high, false).0
                    - trajectory(cx, surface, selected, low, false).0) / (2.0 * step);
                assert!(derivatives[coordinate].abs() > 1e-7, "active physical control {coordinate}");
                assert!((derivatives[coordinate] - expected).abs() < 3e-4 * expected.abs().max(0.01),
                    "surface={surface}, state={selected}, control={coordinate}: complete history {} vs resolve {expected}",
                    derivatives[coordinate]);
            }
        }
    });
}
