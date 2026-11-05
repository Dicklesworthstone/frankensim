//! One complete contact/material/radiosity dual on mixed Robin/flux patches.
use super::*;
use fs_conduction::radiation::EnclosureRadiationLinearization;

fn scaled_card(scale: f64) -> InterfaceSystemCard {
    let original = card(false);
    let mut claim = original.claims().claims_ordered().next().unwrap().1.clone();
    let PropertyValue::Curve { knots, .. } = &mut claim.value else { panic!("curve fixture"); };
    for (_, resistance) in knots { *resistance *= scale; }
    let mut claims = ClaimSet::new();
    claims.insert_claim(claim).unwrap();
    InterfaceSystemCard::assemble(original.surface_a().clone(), original.surface_b().clone(),
        original.context().clone(), claims, vec![]).unwrap()
}

/// Return the same selected nodal goal after an ordinary physical solve.
/// Controls are whole contact-law scale, whole source scale, cold emissivity
/// and the external black-bath temperature, respectively.
fn sample(cx: &Cx<'_>, open: bool, controls: [f64; 4], adjoint: bool)
    -> (f64, Option<[f64; 4]>)
{
    let mesh = mesh();
    let boundary = boundary(&mesh, false);
    let material = material(true);
    let source = ScalarField::Nodal(vec![0., 0., 0., 0.,
        30. * controls[1], 30. * controls[1], 30. * controls[1], 30. * controls[1]]);
    let card = scaled_card(controls[0]);
    let contacts = contacts(cx, &mesh, &boundary, &card);
    let cold = RadiationSurface::new(&mesh, "cold",
        |f| f.centroid[2] < 0.0 && f.centroid[1] == 0.0,
        patch_with_emissivity(controls[2]).emissivity().clone()).unwrap();
    let hot = RadiationSurface::new(&mesh, "hot",
        |f| f.centroid[2] > 0.0 && f.centroid[1] == 0.0,
        patch().emissivity().clone()).unwrap();
    let factors = ViewFactorMatrix::infinite_parallel_plates(cold.area_m2()).unwrap();
    let closed = GrayDiffuseEnclosure::new(vec![cold, hot], factors).unwrap();
    let open_model = OpenDiffuseEnclosure::new(cx, closed.surfaces().to_vec(),
        vec![vec![0.0, 0.25], vec![0.25, 0.0]],
        ViewFactorEvidence::Analytic { geometry: "synthetic reciprocal open exchange".into() },
        ViewFactorTolerance::default()).unwrap();
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
        material: &material, element_materials: None, source: &source };
    let mut solve = config();
    solve.residual_tolerance_w = 1e-9;
    solve.energy_tolerance_w = 1e-9;
    let radiation = CoupledRadiationConfig { max_iterations: 256,
        surface_temperature_atol_k: 1e-11, ..super::controls() };
    let accepted = if open {
        contacts.solve_steady_newton_with_open_gray_diffuse_enclosure(cx, problem,
            &[300.; 8], solve, &open_model, controls[3], radiation).unwrap()
    } else {
        contacts.solve_steady_newton_with_gray_diffuse_enclosure(cx, problem,
            &[300.; 8], solve, &closed, radiation).unwrap()
    };
    let temperature = &accepted.contact.temperature;
    let goal = temperature[7];
    if !adjoint { return (goal, None); }
    let linear = step_config().linear;
    let radiation = if open {
        EnclosureRadiationLinearization::bind_open(cx, problem, linear, temperature,
            &["cold"], &[0.], &[0.], &open_model, controls[3], 1e-9, 4096).unwrap()
    } else {
        EnclosureRadiationLinearization::bind(cx, problem, linear, temperature,
            &["cold"], &[0.], &[0.], &closed, 1e-9, 4096).unwrap()
    };
    // The hot patch is truly prescribed flux, with no artificial Robin port.
    assert!(radiation.flux_feedback().is_some());
    let physical = ConductionProblem { boundary: radiation.boundary(), ..problem };
    let response = contacts.linearize_steady_endpoint(cx, physical, temperature,
        linear, 1e-9, 1e-9, &["cold"], 1).unwrap();
    let weights = [0., 0., 0., 0., 0., 0., 0., 1.];
    let full = response.pullback_with_robin_and_flux_feedback_and_resistance_controls(
        cx, radiation.convective_mean_slopes(), radiation.local_reference_slopes(),
        radiation.reference_feedback(), radiation.flux_feedback(), &weights, 4096, 1).unwrap();
    assert!(full.transport.relative_residual <= linear.tolerance);
    let resistance = full.resistance_multipliers[0].derivative;
    let full = radiation.remap_gradient(cx, full.transport, linear).unwrap();
    // A constant source on the hot tetrahedron integrates to 5 W, divided
    // equally among its four P1 nodes. This is an assembled load contraction,
    // not a mistaken source-density dot product.
    let source_scale = 1.25 * full.convection.nodal_load[4..].iter().sum::<f64>();
    assert_eq!(full.ambient_temperature.is_some(), open);
    let derivatives = [resistance, source_scale, full.emissivities[0],
        full.ambient_temperature.unwrap_or(0.0)];
    let mut stale = temperature.to_vec();
    stale[7] += 0.1;
    assert!(contacts.linearize_steady_endpoint(cx, physical, &stale,
        linear, 1e-9, 1e-9, &["cold"], 1).is_err(),
        "an unchecked radiative/contact field must not acquire a steady dual");
    (goal, Some(derivatives))
}

#[test]
fn complete_live_contact_enclosure_dual_matches_independent_physical_resolves() {
    context(&CancelGate::new_clock_free(), |cx| {
        let nominal = [1.0, 1.0, 0.2, 300.0];
        for open in [false, true] {
            let (_, derivatives) = sample(cx, open, nominal, true);
            let derivatives = derivatives.unwrap();
            for (control, step) in [0.001, 0.001, 0.0002, 0.002].into_iter().enumerate() {
                if control == 3 && !open { continue; }
                let mut low = nominal;
                let mut high = nominal;
                low[control] -= step;
                high[control] += step;
                let expected = (sample(cx, open, high, false).0
                    - sample(cx, open, low, false).0) / (2.0 * step);
                assert!((derivatives[control] - expected).abs()
                    < 5e-5 * expected.abs().max(1e-3),
                    "open={open}, control={control}: complete dual {} vs physical resolve {expected}",
                    derivatives[control]);
            }
        }
    });
}

/// Three original accepted capacity steps; controls are contact/source scales,
/// cold emissivity, black-bath kelvin, total capacity scale and initial kelvin.
fn capacity_sample(cx: &Cx<'_>, open: bool, control: [f64; 6], reverse: bool)
    -> (f64, Option<[f64; 6]>)
{
    use fs_conduction::transient::backward_euler::{ArealHeatCapacity, SurfaceHeatCapacity};
    let mesh = mesh();
    let boundary = boundary(&mesh, false);
    let material = material(true);
    let card = scaled_card(control[0]);
    let contacts = contacts(cx, &mesh, &boundary, &card);
    let cold = RadiationSurface::new(&mesh, "cold",
        |f| f.centroid[2] < 0.0 && f.centroid[1] == 0.0,
        patch_with_emissivity(control[2]).emissivity().clone()).unwrap();
    let hot = RadiationSurface::new(&mesh, "hot",
        |f| f.centroid[2] > 0.0 && f.centroid[1] == 0.0,
        patch().emissivity().clone()).unwrap();
    let coating = SurfaceHeatCapacity::nodal(cold.face_slots()[0],
        [1.0, 2.0, 4.0].map(|value| ArealHeatCapacity::declared(value * control[4]).unwrap()));
    let engine = BackwardEuler::per_element(cx, &mesh,
        &[10.0, 20.0].map(|value| VolumetricHeatCapacity::declared(value * control[4]).unwrap()))
        .unwrap().with_surface_heat_capacity(cx, &[coating]).unwrap();
    let factors = ViewFactorMatrix::infinite_parallel_plates(cold.area_m2()).unwrap();
    let closed = GrayDiffuseEnclosure::new(vec![cold, hot], factors).unwrap();
    let opened = OpenDiffuseEnclosure::new(cx, closed.surfaces().to_vec(),
        vec![vec![0.0, 0.25], vec![0.25, 0.0]],
        ViewFactorEvidence::Analytic { geometry: "synthetic reciprocal open exchange".into() },
        ViewFactorTolerance::default()).unwrap();
    let config = StepConfig { energy_tolerance_j: 1e-9, ..step_config() };
    let coupling = ContactIterationConfig { residual_tolerance_j: 1e-9, ..coupling() };
    let radiation_config = AmbientRadiationConfig { max_iterations: 256,
        temperature_tolerance_k: 1e-12, balance_tolerance_w: 1e-9,
        ..AmbientRadiationConfig::default() };
    let steps = [(0.1, 1.0), (0.25, 0.3), (0.15, 0.0)];
    let mut history = vec![old().into_iter().map(|t| t + control[5]).collect::<Vec<_>>()];
    let mut targets = Vec::new();
    for (dt, power) in steps {
        let source = ScalarField::Nodal(vec![0., 0., 0., 0.,
            30. * control[1] * power, 30. * control[1] * power,
            30. * control[1] * power, 30. * control[1] * power]);
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
            material: &material, element_materials: None, source: &source };
        let previous = history.last().unwrap();
        let solved = if open {
            contacts.advance_capacity_newton_with_open_enclosure_radiation(cx, &engine,
                problem, previous, dt, config, nonlinear(), coupling, 1,
                &opened, control[3], radiation_config).unwrap()
        } else {
            contacts.advance_capacity_newton_with_enclosure_radiation(cx, &engine,
                problem, previous, dt, config, nonlinear(), coupling, 1,
                &closed, radiation_config).unwrap()
        };
        targets.push(solved.physical_residual_tolerance_j);
        history.push(solved.conduction.temperature);
    }
    let goal = history.last().unwrap()[7];
    if !reverse { return (goal, None); }
    let mut carry = vec![0.0; 8];
    carry[7] = 1.0;
    let mut derivatives = [0.0; 6];
    for (index, (dt, power)) in steps.into_iter().enumerate().rev() {
        let source = ScalarField::Nodal(vec![0., 0., 0., 0.,
            30. * power, 30. * power, 30. * power, 30. * power]);
        let problem = ConductionProblem { mesh: &mesh, boundary: &boundary,
            material: &material, element_materials: None, source: &source };
        let temperature = &history[index + 1];
        let radiation = if open {
            EnclosureRadiationLinearization::bind_open(cx, problem, config.linear,
                temperature, &["cold"], &[0.], &[0.], &opened, control[3], 1e-9, 4096).unwrap()
        } else {
            EnclosureRadiationLinearization::bind(cx, problem, config.linear,
                temperature, &["cold"], &[0.], &[0.], &closed, 1e-9, 4096).unwrap()
        };
        assert!(radiation.flux_feedback().is_some(), "hot patch has no fictitious convection");
        let physical = ConductionProblem { boundary: radiation.boundary(), ..problem };
        let response = engine.linearize_temperature_contact_endpoint(cx, physical,
            &contacts, &history[index], temperature, dt, config, targets[index], &["cold"], 1).unwrap();
        assert_eq!(response.primal().temperature, *temperature);
        assert_eq!(response.primal().krylov_iterations, 0);
        let full = response.pullback_with_robin_and_flux_feedback_and_resistance_controls(cx,
            radiation.convective_mean_slopes(), radiation.local_reference_slopes(),
            radiation.reference_feedback(), radiation.flux_feedback(), &carry, 4096, 1).unwrap();
        derivatives[0] += full.resistance_multipliers[0].derivative;
        let full = radiation.remap_gradient(cx, full.transport, config.linear).unwrap();
        derivatives[1] += response.source_multiplier_pullback(cx, &full.convection.nodal_load).unwrap();
        derivatives[2] += full.emissivities[0];
        derivatives[3] += full.ambient_temperature.unwrap_or(0.0);
        derivatives[4] += response.capacity_multiplier_pullback(cx, &full.convection.nodal_load).unwrap();
        carry = response.previous_temperature_pullback(cx, &full.convection.nodal_load).unwrap();
        if index == 2 {
            let mut stale = temperature.clone();
            stale[7] += 0.1;
            assert!(engine.linearize_temperature_contact_endpoint(cx, physical, &contacts,
                &history[index], &stale, dt, config, targets[index], &["cold"], 1).is_err());
            assert!(response.pullback_with_robin_and_flux_feedback_and_resistance_controls(cx,
                radiation.convective_mean_slopes(), radiation.local_reference_slopes(),
                radiation.reference_feedback(), radiation.flux_feedback(), &[0.; 8], 0, 1).is_err());
        }
    }
    derivatives[5] = carry.iter().sum();
    (goal, Some(derivatives))
}

#[test]
fn capacity_contact_enclosure_reverse_keeps_coating_storage_and_complete_history() {
    context(&CancelGate::new_clock_free(), |cx| {
        let nominal = [1.0, 1.0, 0.2, 295.0, 1.0, 0.0];
        for open in [false, true] {
            let (goal, derivatives) = capacity_sample(cx, open, nominal, true);
            near(goal, capacity_sample(cx, open, nominal, false).0, 0.0);
            let derivatives = derivatives.unwrap();
            for (coordinate, step) in [0.001, 0.001, 0.0002, 0.002, 0.001, 0.002].into_iter().enumerate() {
                if coordinate == 3 && !open { continue; }
                let mut low = nominal;
                let mut high = nominal;
                low[coordinate] -= step;
                high[coordinate] += step;
                let expected = (capacity_sample(cx, open, high, false).0
                    - capacity_sample(cx, open, low, false).0) / (2.0 * step);
                assert!((derivatives[coordinate] - expected).abs()
                    < 3e-4 * expected.abs().max(0.01),
                    "open={open}, coordinate={coordinate}: complete history {} versus physical resolves {expected}",
                    derivatives[coordinate]);
            }
        }
    });
}
