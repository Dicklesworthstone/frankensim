use super::*;
use fs_material::thermovisco::{Limits, ThermalKnot};

fn law(branches: Vec<(f64, f64)>) -> ThermoMaxwell {
    ThermoMaxwell::new(1e8, branches, [300.0, 500.0].into_iter().map(|t| ThermalKnot {
        temperature_k: t, modulus_scale: 1.0, relaxation_scale: 1.0, free_strain: (t-300.0)*1e-5,
    }).collect(), Limits { max_strain: 0.03, max_strain_rate_s: 10.0 }).unwrap()
}
fn strip<'a>() -> LayeredStrip<'a> {
    LayeredStrip::new(0.01, Budget { max_layers: 32, max_branch_states: 1024 }).unwrap()
}
fn geometry(z: f64, thickness: f64) -> LayerGeometry {
    LayerGeometry { area_m2: 0.002*thickness, centroid_z_m: z, thickness_m: thickness }
}
fn config() -> Config {
    Config { force_tolerance_n: 1e-8, moment_tolerance_nm: 1e-10, energy_tolerance_j: 1e-9,
        max_strain: 0.03, max_rotation_rad: 0.1 }
}
fn close(a: f64, b: f64, tol: f64) { assert!((a-b).abs() < tol, "{a:e} != {b:e}"); }

#[test]
fn uniform_free_contraction_and_fixed_reaction_match_analytic_bar() {
    let law = law(vec![]);
    let strip = strip().deposit(1, geometry(0.0, 0.001), &law, 500.0).unwrap();
    let t = BTreeMap::from([(1, 300.0)]);
    let (free, r) = strip.advance_to(1.0, &t, Control::RELEASED, config(), || Ok(())).unwrap();
    close(free.axial_strain(), -0.002, 1e-14); close(free.curvature_per_m(), 0.0, 1e-12);
    close(free.deformed_position_m(0.01, 0.0).unwrap()[0], 0.00998, 1e-14);
    close(r.axial_force_n, 0.0, 1e-8);
    let (_, held) = strip.advance_to(1.0, &t, Control::HELD, config(), || Ok(())).unwrap();
    close(held.axial_reaction_n, 0.4, 1e-12); close(held.bending_reaction_nm, 0.0, 1e-12);
    close(held.energy_j.mechanical_work, 0.0, 1e-12);
}

#[test]
fn differential_layer_cooling_bends_after_prescribed_release() {
    let law = law(vec![]);
    let initial = strip().deposit(1, geometry(-0.0005, 0.001), &law, 500.0).unwrap()
        .deposit(2, geometry(0.0005, 0.001), &law, 500.0).unwrap();
    let t = BTreeMap::from([(1, 300.0), (2, 400.0)]);
    let (held, reactions) = initial.advance_to(1.0, &t, Control::HELD, config(), || Ok(())).unwrap();
    assert!(reactions.axial_reaction_n > 0.0 && reactions.bending_reaction_nm > 0.0);
    assert_eq!(held.deformed_position_m(0.01, 0.0).unwrap()[1], 0.0);
    let (released, result) = held.advance_to(1.1, &t, Control::RELEASED, config(), || Ok(())).unwrap();
    close(released.axial_strain(), -0.0015, 1e-14);
    close(released.curvature_per_m(), -0.75, 1e-12);
    close(released.deformed_position_m(0.01, 0.0).unwrap()[1], -3.75e-5, 1e-14);
    assert!(result.energy_j.stored_change < 0.0);
    assert!(result.force_residual_n.abs() < 1e-8 && result.moment_residual_nm.abs() < 1e-10);
    assert!(result.energy_j.residual.abs() < 1e-9);
}

#[test]
fn new_layer_uses_current_deformed_reference_without_resetting_old_stress() {
    let law = law(vec![(3e8, 2.0)]);
    let strip = strip().deposit(1, geometry(0.0, 0.001), &law, 400.0).unwrap();
    let control = Control { axial: AxialControl::Strain(0.001), bending: BendingControl::CurvaturePerM(0.2) };
    let (loaded, _) = strip.advance_to(1.0, &BTreeMap::from([(1, 350.0)]), control, config(), || Ok(())).unwrap();
    let born = loaded.deposit(2, geometry(0.001, 0.001), &law, 450.0).unwrap();
    assert_eq!(loaded.layers()[&1], born.layers()[&1]);
    for (point, z) in born.layers()[&2].material_points().iter().zip(born.layers()[&2].geometry.points()) {
        assert_eq!(point.stress_pa(), 0.0);
        close(point.birth_strain(), loaded.axial_strain()-z*loaded.curvature_per_m(), 1e-15);
    }
    assert_eq!(born.layers()[&2].birth_time_s, 1.0);
}

#[test]
fn rectangle_refinement_preserves_axial_bending_and_history_energy() {
    let law = law(vec![(3e8, 2.0)]);
    let a = strip().deposit(1, geometry(0.0, 0.002), &law, 400.0).unwrap();
    let b = strip().deposit(1, geometry(-0.0005, 0.001), &law, 400.0).unwrap()
        .deposit(2, geometry(0.0005, 0.001), &law, 400.0).unwrap();
    let load = Control { axial: AxialControl::Strain(0.001), bending: BendingControl::CurvaturePerM(0.2) };
    let (a, ar) = a.advance_to(1.0, &BTreeMap::from([(1, 350.0)]), load, config(), || Ok(())).unwrap();
    let (b, br) = b.advance_to(1.0, &BTreeMap::from([(1, 350.0), (2, 350.0)]), load, config(), || Ok(())).unwrap();
    close(ar.axial_force_n, br.axial_force_n, 1e-12); close(ar.bending_moment_nm, br.bending_moment_nm, 1e-12);
    close(a.stored_energy_j(), b.stored_energy_j(), 1e-12);
    close(ar.energy_j.dissipated, br.energy_j.dissipated, 1e-12);
}

#[test]
fn cold_relaxation_keeps_force_balance_and_irreversible_history_after_release() {
    let law = law(vec![(3e8, 2.0)]);
    let initial = strip().deposit(1, geometry(0.0, 0.001), &law, 500.0).unwrap();
    let t = BTreeMap::from([(1, 300.0)]);
    let (held, a) = initial.advance_to(1.0, &t, Control::HELD, config(), || Ok(())).unwrap();
    let (relaxed, b) = held.advance_to(3.0, &t, Control::HELD, config(), || Ok(())).unwrap();
    assert!(b.axial_reaction_n < a.axial_reaction_n && b.energy_j.dissipated > 0.0);
    let (free, r) = relaxed.advance_to(3.1, &t, Control::RELEASED, config(), || Ok(())).unwrap();
    assert!(free.axial_strain() < 0.0 && free.axial_strain() > -0.002);
    assert!(free.layers()[&1].material_points()[0].branch_strains()[0].abs() > 1e-6);
    close(r.axial_force_n, 0.0, 1e-8);
}

#[test]
fn cancellation_and_large_rotation_leave_original_checkpoint_unchanged() {
    let law = law(vec![(3e8, 2.0)]);
    let initial = strip().deposit(1, geometry(0.0, 0.001), &law, 400.0).unwrap();
    let old = initial.clone(); let t = BTreeMap::from([(1, 350.0)]);
    let mut checkpoints = 0;
    assert!(initial.advance_to(1.0, &t, Control::HELD, config(), || {
        checkpoints += 1; if checkpoints >= 3 { Err(Error::Cancelled) } else { Ok(()) }
    }).is_err());
    assert_eq!(initial, old);
    let excessive = Control { axial: AxialControl::Strain(0.0), bending: BendingControl::CurvaturePerM(20.0) };
    assert!(matches!(initial.advance_to(1.0, &t, excessive, config(), || Ok(())), Err(Error::Regime(_))));
    assert_eq!(initial, old);
}
