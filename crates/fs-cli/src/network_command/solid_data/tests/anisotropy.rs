use super::*;

fn one_material(r: &mut Request, conductivity: &str) {
    let names = vec!["\"solid\""; r.mesh.element_count()].join(",");
    assign(r, &format!(r#"{{"materials":[{{"name":"solid",{conductivity},"source":"analytic fixture"}}],"element_materials":[{names}],"source_w_m3":0}}"#));
}

const ALIGNED: &str = r#""orthotropic":{"principal_axes":[[1,0,0],[0,1,0],[0,0,1]],"conductivity_w_m_k":[20,2,1]}"#;
const ROTATED: &str = r#""orthotropic":{"principal_axes":[[0.6,0.8,0],[-0.8,0.6,0],[0,0,1]],"conductivity_w_m_k":[20,2,1]}"#;

#[test]
fn isotropic_tensor_retains_the_legacy_field_and_implicit_derivatives() {
    let mut r = request();
    with_cx(|cx| {
        let flow = r.flow(cx).unwrap();
        let h = r.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
        let baseline = r.evaluate(cx, &flow, &h, true).unwrap();
        one_material(&mut r, r#""conductivity_tensor_w_m_k":[[10,0,0],[0,10,0],[0,0,10]]"#);
        // Deliberately poison the inactive compatibility coefficient with a
        // different POSITIVE value. ElementMaterials must remain authoritative.
        r.conductivity = 1.0e-6;
        let tensor = r.evaluate(cx, &flow, &h, true).unwrap();
        for (&a, &b) in baseline.temperatures.iter().zip(&tensor.temperatures) { close(a, b, 1e-7); }
        let a = baseline.gradient.unwrap();
        let b = tensor.gradient.unwrap();
        for (&a, &b) in a.log_htc.iter().zip(&b.log_htc) { close(a, b, 1e-7); }
        for (&a, &b) in a.inlets.iter().zip(&b.inlets) { close(a, b, 1e-7); }
    });
}

#[test]
fn turning_the_material_changes_the_heat_flow_by_the_physical_series_resistance() {
    let mut heats = Vec::new();
    for (declaration, axial_k) in [
        (ALIGNED, 20.0),
        (r#""orthotropic":{"principal_axes":[[0,1,0],[-1,0,0],[0,0,1]],"conductivity_w_m_k":[20,2,1]}"#, 2.0),
    ] {
        let mut r = request();
        one_material(&mut r, declaration);
        with_cx(|cx| {
            let flow = r.flow(cx).unwrap();
            let h = r.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
            let result = r.evaluate(cx, &flow, &h, false).unwrap();
            let capacity_first = 1.2 * 0.003 * 1007.0_f64;
            let capacity_mixed = 1.2 * 0.004 * 1007.0_f64;
            let film_first = capacity_first * -(-0.5 / capacity_first).exp_m1();
            let film_last = capacity_mixed * -(-0.8 / capacity_mixed).exp_m1();
            // Same mixed-inlet problem as the retained isotropic fixture,
            // but its solid resistance is L/(k_axial*A), not an average of k.
            let heat = 10.0 / (0.05 / (axial_k * 0.01)
                + 1.0 / film_first + 1.0 / film_last - 1.0 / capacity_mixed);
            let first_wall = 330.0 - heat / film_first;
            for (position, &temperature) in r.mesh.positions().iter().zip(&result.temperatures) {
                close(temperature, first_wall - heat * position[0] / (axial_k * 0.01), 1e-5);
            }
            close(result.coupled.solid[0].heat_rate_w, -heat, 1e-5);
            close(result.coupled.solid[1].heat_rate_w, heat, 1e-5);
            heats.push(result.coupled.solid[1].heat_rate_w);
        });
    }
    assert!(heats[0] - heats[1] > 0.1, "orientation must change physical heat flow");
}

#[test]
fn rotating_geometry_and_material_together_preserves_the_coupled_field() {
    let mut r = request();
    one_material(&mut r, ALIGNED);
    with_cx(|cx| {
        let flow = r.flow(cx).unwrap();
        let h = r.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
        let original = r.evaluate(cx, &flow, &h, true).unwrap();
        // R has columns equal to the declared principal-axis rows below.
        // Rotate the mesh itself, not just a material label or a test oracle.
        let positions = r.mesh.positions().iter().map(|p| {
            [0.6 * p[0] - 0.8 * p[1], 0.8 * p[0] + 0.6 * p[1], p[2]]
        }).collect();
        r.mesh = ConductionMesh::new(r.mesh.complex().clone(), positions).unwrap();
        one_material(&mut r, ROTATED);
        let rotated = r.evaluate(cx, &flow, &h, true).unwrap();
        for (&a, &b) in original.temperatures.iter().zip(&rotated.temperatures) { close(a, b, 1e-6); }
        for (a, b) in original.coupled.solid.iter().zip(&rotated.coupled.solid) {
            close(a.heat_rate_w, b.heat_rate_w, 1e-6);
        }
        for (&a, &b) in original.gradient.unwrap().log_htc.iter().zip(&rotated.gradient.unwrap().log_htc) {
            close(a, b, 1e-6);
        }
    });
}

#[test]
fn off_diagonal_tensor_adjoint_matches_full_perturbed_cooling_solves() {
    let mut r = request();
    one_material(&mut r, ROTATED);
    r.source = 1000.0;
    with_cx(|cx| {
        let flow = r.flow(cx).unwrap();
        let mut h: BTreeMap<_, _> = r.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
        let result = r.evaluate(cx, &flow, &h, true).unwrap();
        let delta = 1e-4_f64;
        h.insert("last-face".into(), 80.0 * delta.exp());
        let plus = r.evaluate(cx, &flow, &h, false).unwrap().objective;
        h.insert("last-face".into(), 80.0 * (-delta).exp());
        let minus = r.evaluate(cx, &flow, &h, false).unwrap().objective;
        close(result.gradient.unwrap().log_htc[1], (plus - minus) / (2.0 * delta), 5e-5);
    });
}
