//! Physical port reduction checks, not measured-instrument validation.
use fs_couple::bernoulli_aperture::plate::{PlateApertureOptions, PlateApertureReduction};
use fs_exec::CancelGate;
use fs_plate::{AssemblyOptions, EdgeSupport, PlateChart, PlateMesh, PlateSection, SliceOptions};

fn fixture(e: f64, density: f64) -> (PlateChart, PlateApertureOptions) {
    let mesh = PlateMesh::rectangle(0.025, 0.01, 6, 2);
    let root = mesh.nodes.iter().enumerate().filter(|(_, p)| p.0 == 0.0).map(|(i, _)| i).collect();
    let edges = mesh.boundary_edges().into_iter().filter(|&(a, b)|
        (mesh.nodes[a].0 - 0.025).abs() < 1e-12 && (mesh.nodes[b].0 - 0.025).abs() < 1e-12)
        .map(|(a, b)| [a, b]).collect();
    let section = PlateSection::isotropic(e, 0.3, 0.0003, density).unwrap();
    let chart = PlateChart::with_boundary_and_regions(mesh, section, root, vec![]).unwrap();
    (chart, PlateApertureOptions {
        assembly: AssemblyOptions { pretension: 0.0, support: EdgeSupport::Clamped },
        eigenvalue_window: (1.0, (core::f64::consts::TAU * 450.0).powi(2)), mode_index: 0,
        eigensolver: SliceOptions::default(), max_nodes: 100, max_triangles: 100,
        slit_edges: edges, rest_opening_m: 0.0002, damping_ratio: 0.02,
        max_slit_mode_variation: 0.25, max_slope: 0.1,
    })
}
fn reduction(e: f64, density: f64) -> PlateApertureReduction {
    let (c, o) = fixture(e, density);
    PlateApertureReduction::from_chart(c, o, &CancelGate::new()).unwrap()
}
fn near(actual: f64, expected: f64, tol: f64) {
    assert!((actual - expected).abs() <= tol * expected.abs().max(1e-30), "{actual:e} != {expected:e}");
}

#[test]
fn plate_mass_stiffness_and_pressure_work_share_the_reconstructed_opening_coordinate() {
    let r = reduction(4e9, 900.0);
    let model = r.chart().assemble(&[], &r.options().assembly).unwrap();
    let mut shape = vec![0.0; model.free];
    for (node, value) in r.shape_per_opening().iter().enumerate() {
        for c in 0..3 { if let Some(i) = model.dof_map[3*node+c] { shape[i] = value[c]; } }
    }
    let mut product = vec![0.0; shape.len()];
    model.m.spmv(&shape, &mut product);
    near(r.mass_kg(), shape.iter().zip(&product).map(|(a,b)|a*b).sum(), 1e-10);
    model.k.spmv(&shape, &mut product);
    near(r.stiffness_n_m(), shape.iter().zip(&product).map(|(a,b)|a*b).sum(), 1e-6);
    let mut area = 0.0;
    for &[a,b,c] in &r.chart().mesh.tris {
        let [p,q,t] = [a,b,c].map(|i| r.chart().mesh.nodes[i]);
        let da = 0.5*((q.0-p.0)*(t.1-p.1)-(q.1-p.1)*(t.0-p.0));
        area += da * (r.shape_per_opening()[a][0] + r.shape_per_opening()[b][0]
            + r.shape_per_opening()[c][0])/3.0;
    }
    near(r.pressure_area_m2(), area, 1e-12);
    near(r.width_m(), 0.01, 1e-12);
    let spec = r.dynamic_spec(1.2, 1e6, 1e-5, 1000);
    near(spec.stiffness_n_m * spec.aperture.rest_opening_m,
        spec.aperture.closing_pressure_pa * area, 1e-12);
    let (dp, velocity) = (700.0, -0.012);
    near((-dp * area) * velocity, dp * (-area * velocity), 1e-14);
    assert!(r.in_window_modes() >= 1);
    assert!(r.material_chart().is_none());
}

#[test]
fn changing_elasticity_and_density_changes_distinct_physical_coefficients_without_retuning() {
    let a = reduction(4e9, 900.0);
    let stiff = reduction(8e9, 900.0);
    let heavy = reduction(4e9, 1800.0);
    near(stiff.stiffness_n_m()/a.stiffness_n_m(), 2.0, 2e-5);
    near(stiff.mass_kg()/a.mass_kg(), 1.0, 2e-5);
    near(heavy.mass_kg()/a.mass_kg(), 2.0, 2e-5);
    near(heavy.stiffness_n_m()/a.stiffness_n_m(), 1.0, 2e-5);
    near(heavy.closing_pressure_pa()/a.closing_pressure_pa(), 1.0, 2e-5);
    near(stiff.pressure_area_m2()/a.pressure_area_m2(), 1.0, 2e-5);
    // Same static closure, different resonant mechanics; density is not a gain.
    assert!(heavy.stiffness_n_m()/heavy.mass_kg() < 0.51*a.stiffness_n_m()/a.mass_kg());
}

#[test]
fn regional_thickness_and_material_sections_are_used_in_the_actual_plate_pencil() {
    let baseline = reduction(4e9, 900.0);
    let (chart, options) = fixture(4e9, 900.0);
    let sections = chart.mesh.tris.iter().map(|tri| {
        let x = tri.iter().map(|&i| chart.mesh.nodes[i].0).sum::<f64>()/3.0;
        PlateSection::isotropic(if x < 0.0125 { 4e9 } else { 2e9 }, 0.3,
            if x < 0.0125 { 0.0003 } else { 0.0002 }, 900.0).unwrap()
    }).collect();
    let changed = PlateApertureReduction::from_chart(chart.with_element_sections(sections).unwrap(),
        options, &CancelGate::new()).unwrap();
    assert!((changed.stiffness_n_m()/baseline.stiffness_n_m()-1.0).abs() > 0.05);
    assert!((changed.mass_kg()/baseline.mass_kg()-1.0).abs() > 0.05);
    assert!((changed.pressure_area_m2()/baseline.pressure_area_m2()-1.0).abs() > 1e-4);
}

#[test]
fn retained_shape_obeys_the_declared_slope_and_slit_limits() {
    let r = reduction(4e9, 900.0);
    let rest = r.options().rest_opening_m;
    assert_eq!(r.max_slope_at(rest), 0.0);
    r.validate_opening(rest).unwrap();
    let at_unit = r.max_slope_at(rest + 1.0);
    r.validate_opening(rest + 0.5*r.options().max_slope/at_unit).unwrap();
    assert!(r.validate_opening(rest + 2.0*r.options().max_slope/at_unit).is_err());
    assert!(r.validate_opening(f64::NAN).is_err());
    let (c, mut o) = fixture(4e9, 900.0);
    // Add a stationary root edge to the moving slit: the normalized motion
    // cannot be uniform, so the explicit approximation bound must refuse.
    let (a,b) = c.mesh.boundary_edges().into_iter().find(|&(a,b)|
        c.mesh.nodes[a].0 == 0.0 && c.mesh.nodes[b].0 == 0.0).unwrap();
    o.slit_edges.push([a,b]);
    assert!(PlateApertureReduction::from_chart(c, o, &CancelGate::new()).is_err());
}

#[test]
fn malformed_geometry_budget_and_cancellation_refuse_without_substitute_modes() {
    let (c, o) = fixture(4e9, 900.0);
    let gate = CancelGate::new(); gate.request();
    assert!(PlateApertureReduction::from_chart(c.clone(), o.clone(), &gate).is_err());
    let mut bad = Vec::new();
    let mut a = o.clone(); a.max_nodes = 1; bad.push(a);
    let mut a = o.clone(); a.max_triangles = 1; bad.push(a);
    let mut a = o.clone(); a.slit_edges.push(a.slit_edges[0]); bad.push(a);
    let mut a = o.clone(); a.slit_edges = vec![[0,usize::MAX]]; bad.push(a);
    let mut a = o.clone(); a.slit_edges.clear(); bad.push(a);
    let mut a = o.clone(); a.mode_index = usize::MAX; bad.push(a);
    let mut a = o.clone(); a.damping_ratio = f64::NAN; bad.push(a);
    let mut a = o.clone(); a.rest_opening_m = 0.0; bad.push(a);
    for options in bad { assert!(PlateApertureReduction::from_chart(c.clone(), options, &CancelGate::new()).is_err()); }
}

#[path = "plate_aperture/runtime.rs"]
mod runtime;

#[path = "plate_aperture/closure.rs"]
mod closure;

#[path = "plate_aperture/relaxation.rs"]
mod relaxation;

#[path = "plate_aperture/performance.rs"]
mod performance;

#[path = "plate_aperture/file.rs"]
mod file;

#[path = "plate_aperture/viscothermal.rs"]
mod viscothermal;


#[path = "plate_aperture/regional_gas.rs"]
mod regional_gas;
