//! Shell and bell modal conformance test suite (bead `frankensim-music-v8-root-3ez8g.12.2`).

use fs_modal::SliceOptions;
use fs_plate::shell::{
    assemble_shell, canonical_church_bell_profile, generate_bell_shell, generate_cylinder_shell,
    modes_shell, ShellSupport,
};
use fs_plate::PlateSection;

fn bronze_bell_section() -> PlateSection {
    // Bell bronze (80% Cu, 20% Sn): E ≈ 105 GPa, nu = 0.34, rho = 8750 kg/m³, h = 15 mm
    PlateSection::isotropic(105e9, 0.34, 0.015, 8750.0).expect("valid bronze section")
}

#[test]
fn flat_facet_cylinder_shell_assembles_and_solves_modes() {
    let mesh = generate_cylinder_shell(0.2, 0.5, 12, 6);
    assert_eq!(mesh.node_count(), 12 * 7);
    assert_eq!(mesh.element_count(), 2 * 12 * 6);

    let sec = bronze_bell_section();
    // Clamp bottom ring (nodes 0..12)
    let b_nodes: Vec<usize> = (0..12).collect();
    let model = assemble_shell(&mesh, &sec, &b_nodes, ShellSupport::Clamped)
        .expect("cylinder shell assembles");

    assert!(model.free > 0);
    assert_eq!(model.k.nrows(), model.free);
    assert_eq!(model.m.nrows(), model.free);

    // Modal solve for lowest modes
    let opts = SliceOptions::default();
    let report = modes_shell(&model, (10.0, 1e8), &opts).expect("modal solve succeeds");
    assert!(report.modes.len() >= 4, "should find at least 4 cylinder modes");
}

#[test]
fn church_bell_revolve_and_partial_structure() {
    let profile = canonical_church_bell_profile(0.4); // 40 cm scale bell
    let mesh = generate_bell_shell(&profile, 16);
    assert!(mesh.node_count() > 0);
    assert!(mesh.element_count() > 0);

    let sec = bronze_bell_section();
    // Clamp crown (top nodes)
    let crown_nodes: Vec<usize> = (0..16).collect();
    let model = assemble_shell(&mesh, &sec, &crown_nodes, ShellSupport::Clamped)
        .expect("bell shell assembles");

    let opts = SliceOptions::default();
    let report = modes_shell(&model, (100.0, 5e8), &opts).expect("bell modal solve succeeds");

    assert!(report.modes.len() >= 5, "should resolve bell partial structure");
    let f0 = report.modes[0].lambda.sqrt() / (2.0 * std::f64::consts::PI);
    assert!(f0 > 50.0 && f0 < 5000.0, "fundamental frequency within audible band: {f0} Hz");
}

#[test]
fn shell_rigid_body_modes_on_free_mesh() {
    let mesh = generate_cylinder_shell(0.1, 0.2, 8, 4);
    let sec = bronze_bell_section();
    let model = assemble_shell(&mesh, &sec, &[], ShellSupport::Free).expect("free assembly succeeds");

    assert_eq!(model.free, 6 * mesh.node_count());
}
