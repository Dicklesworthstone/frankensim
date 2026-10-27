//! Physical internal-coordinate transfer, not a normalized-amplitude reset.
use fs_duct::{Duct, Segment};
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, DgState, DgWallRelaxationTransfer};
use fs_material::gas::{GasSpec, GasState};
use fs_material::visco::GeneralizedMaxwell;
use fs_phs::WallPin;

fn pin() -> WallPin {
    WallPin { surface_density: 0.2, stiffness_per_area: 2e7, resistance: 2.0 }
}
fn spectrum() -> GeneralizedMaxwell {
    GeneralizedMaxwell::new(1e6, vec![(3e5, 0.0002)]).unwrap()
}
fn started() -> (DgDuct, DgState) {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let duct = DgDuct::new(
        &Duct { segments: vec![Segment::Cylinder { radius: 0.006, length: 0.3 }] },
        &gas, DgMesh { elements: 2, order: 2 }, DgBoundary::Closed, DgBoundary::Closed,
    ).unwrap().with_compliant_wall_relaxation(pin(), &spectrum()).unwrap();
    let propagation = duct.propagator(1e-6, 1e-13).unwrap();
    let mut state = duct.project_initial_state(|_| 1.0, |_| 0.0);
    for _ in 0..60 { propagation.step(&mut state, &[]).unwrap(); }
    propagation.sync_nodal(&mut state).unwrap();
    assert!(duct.compliant_wall_state(&state).unwrap().relaxation_stored_energy_j > 0.0);
    (duct, state)
}
fn mapping() -> [DgWallRelaxationTransfer; 2] {
    core::array::from_fn(|element| DgWallRelaxationTransfer {
        element, previous_arm: 0, replacement_arm: 0,
    })
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-11 * a.abs().max(b.abs()).max(1e-24), "{a:e} != {b:e}");
}

#[test]
fn identity_extension_maps_preserve_the_original_update_arithmetic() {
    let (mut original, mut a) = started();
    let mut mapped = original.clone();
    let mut b = a.clone();
    let next = WallPin { surface_density: 0.3, stiffness_per_area: 3e7, resistance: 3.0 };
    let law = GeneralizedMaxwell::new(1e6, vec![(5e5, 0.0003)]).unwrap();
    let report_a = original.rebind_compliant_wall_relaxation_regions(
        &mut a, &[next; 2], &[Some(&law); 2], &mapping(), 1.0, &mut || false,
    ).unwrap();
    let report_b = mapped.rebind_compliant_wall_relaxation_regions_with_extension_maps(
        &mut b, &[next; 2], &[Some(&law); 2], &mapping(), &[[1.0, 0.0]; 2], 1.0, &mut || false,
    ).unwrap();
    assert_eq!(report_a, report_b);
    assert_eq!(a, b);
    assert_eq!(original.operator_identity(), mapped.operator_identity());
}

#[test]
fn source_coordinate_maps_preserve_displacement_momentum_and_actual_storage() {
    let (mut duct, mut state) = started();
    let old_state = state.clone();
    let before = duct.compliant_wall_state(&state).unwrap();
    let next = WallPin { surface_density: 0.3, stiffness_per_area: 3e7, resistance: 3.0 };
    let law = GeneralizedMaxwell::new(1e6, vec![(5e5, 0.0003)]).unwrap();
    let transforms = [[-0.7, 0.25], [1.4, -0.1]];
    let report = duct.rebind_compliant_wall_relaxation_regions_with_extension_maps(
        &mut state, &[next; 2], &[Some(&law); 2], &mapping(), &transforms, 1.0, &mut || false,
    ).unwrap();
    let after = duct.compliant_wall_state(&state).unwrap();
    let old_k = pin().stiffness_per_area * spectrum().terms[0].0 / spectrum().e_inf;
    let new_k = next.stiffness_per_area * law.terms[0].0 / law.e_inf;
    for i in 0..before.displacement_m.len() {
        close(after.displacement_m[i], before.displacement_m[i]);
        close(next.surface_density * after.velocity_m_s[i], pin().surface_density * before.velocity_m_s[i]);
        let [a, b] = transforms[i / 3];
        close(after.relaxation_pressure_pa[0][i] / new_k,
            a * before.relaxation_pressure_pa[0][i] / old_k + b * before.displacement_m[i]);
    }
    assert_eq!(state.pressure(), old_state.pressure());
    assert_eq!(state.step_count(), old_state.step_count());
    close(report.parameter_work_j, report.wall_energy_after_j - report.wall_energy_before_j);
    close(report.energy_after_j - report.energy_before_j, report.parameter_work_j);
    // The replacement uses the same propagator and its physical dissipation.
    let propagation = duct.propagator(1e-6, 1e-13).unwrap();
    let initial = state.stored_energy_j();
    for _ in 0..50 { propagation.step(&mut state, &[]).unwrap(); }
    assert!(state.stored_energy_j() < initial);
}

#[test]
fn malformed_maps_budget_and_late_cancellation_leave_both_owners_unchanged() {
    let (original, state) = started();
    let law = spectrum();
    for transforms in [vec![], vec![[0.0, 1.0]; 2], vec![[f64::NAN, 0.0]; 2],
        vec![[1.0, f64::INFINITY]; 2], vec![[f64::MAX, f64::MAX]; 2]] {
        let (mut duct, mut candidate) = (original.clone(), state.clone());
        assert!(duct.rebind_compliant_wall_relaxation_regions_with_extension_maps(
            &mut candidate, &[pin(); 2], &[Some(&law); 2], &mapping(), &transforms, 1.0, &mut || false,
        ).is_err());
        assert_eq!(candidate, state);
        assert_eq!(duct.operator_identity(), original.operator_identity());
    }
    let transforms = [[2.0, 0.4]; 2];
    let (mut duct, mut candidate) = (original.clone(), state.clone());
    assert!(duct.rebind_compliant_wall_relaxation_regions_with_extension_maps(
        &mut candidate, &[pin(); 2], &[Some(&law); 2], &mapping(), &transforms, 0.0, &mut || false,
    ).is_err());
    assert_eq!(candidate, state);
    assert_eq!(duct.operator_identity(), original.operator_identity());
    let mut polls = 0;
    duct.rebind_compliant_wall_relaxation_regions_with_extension_maps(
        &mut candidate, &[pin(); 2], &[Some(&law); 2], &mapping(), &transforms, 1.0,
        &mut || { polls += 1; false },
    ).unwrap();
    for stop in [1, polls / 2, polls] {
        let (mut duct, mut candidate) = (original.clone(), state.clone());
        let mut count = 0;
        assert!(duct.rebind_compliant_wall_relaxation_regions_with_extension_maps(
            &mut candidate, &[pin(); 2], &[Some(&law); 2], &mapping(), &transforms, 1.0,
            &mut || { count += 1; count == stop },
        ).is_err());
        assert_eq!(candidate, state);
        assert_eq!(duct.operator_identity(), original.operator_identity());
    }
}
