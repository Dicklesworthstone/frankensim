//! Manufactured pressure-face geometry; no measured reed or threshold claim.
use super::*;
use fs_couple::bernoulli_aperture::modal::wetted::MaterialWettedSurface;
use fs_plate::shell::pressure::{OrientedShellFacet, ShellFace};

fn surface(valve: &ModalAperture, blade: usize, side: ShellFace, half: bool)
    -> MaterialWettedSurface
{
    let plate = &valve.material_plates()[blade];
    let mesh = &plate.source().chart().mesh;
    let facets = mesh.tris.iter().enumerate().filter(|(_, tri)| {
        !half || tri.iter().map(|&node| mesh.nodes[node].0).sum::<f64>() / 3.0 > 0.0125
    }).map(|(triangle, _)| OrientedShellFacet { triangle, side }).collect();
    MaterialWettedSurface::from_plate(plate, facets, "manufactured tagged pressure faces".into())
}
fn attached(valve: ModalAperture, blade: usize, side: ShellFace, half: bool) -> ModalAperture {
    let selected = surface(&valve, blade, side, half);
    valve.with_material_wetted_surface(blade, selected, 100000, &CancelGate::new_clock_free()).unwrap()
}

#[test]
fn exposed_original_faces_project_pressure_work_and_live_swept_volume_together() {
    let original = bound(vec![plate(900., 2e4), plate(1100., 2e4)], 1e6);
    let original_modes = original.spec().modes.clone();
    let mut valve = attached(original, 0, ShellFace::Positive, true);
    let plate = &valve.material_plates()[0];
    let port = plate.wetted_surface().unwrap();
    assert_eq!(port.surface().specimen_identity(), plate.source().specimen_identity());
    near(port.area_m2(), 0.025 * 0.01 / 2.0, 1e-14);
    let mesh = &plate.source().chart().mesh;
    for (local, index) in plate.mode_range().enumerate() {
        // Independent centroid quadrature is exact for uniform pressure times
        // this original linear displacement trace on each selected triangle.
        let reference: f64 = port.surface().upstream_facets().iter().map(|facet| {
            let tri = mesh.tris[facet.triangle];
            let [a,b,c] = tri.map(|node| mesh.nodes[node]);
            let area = ((b.0-a.0)*(c.1-a.1)-(b.1-a.1)*(c.0-a.0)) / 2.0;
            area * tri.iter().map(|&node| plate.displacement_shapes()[node][local]).sum::<f64>() / 3.0
        }).sum();
        assert!((reference-port.modal_area_m2()[local]).abs() < 3e-18);
        assert_eq!(valve.spec().modes[index].pressure_area_m2, port.modal_area_m2()[local]);
        let mut changed = valve.spec().modes[index];
        changed.pressure_area_m2 = original_modes[index].pressure_area_m2;
        assert_eq!(changed, original_modes[index]);
    }
    assert_eq!(valve.spec().modes[3..], original_modes[3..]);
    let dt = valve.spec().time_step_s;
    for _ in 0..32 {
        let before = valve.modal_states().to_vec();
        let frame = valve.step(drive(5.0), [0.0; 2]).unwrap();
        let displaced: f64 = valve.spec().modes.iter().zip(before.iter().zip(&frame.modes))
            .map(|(mode,(old,new))| -mode.pressure_area_m2 * (new.displacement_m-old.displacement_m)).sum();
        assert!((displaced-dt*frame.aperture.swept_flow_m3_s).abs() < 1e-24);
        let mechanical_pressure_work = frame.pressure_drop_pa * displaced;
        let received = frame.aperture.pressure_work_j
            - dt*frame.pressure_drop_pa*frame.aperture.jet_flow_m3_s;
        assert!((received-mechanical_pressure_work).abs() < 1e-22);
        let scale = frame.aperture.stored_energy_j + frame.aperture.dissipated_energy_j.abs()
            + frame.aperture.pressure_work_j.abs();
        assert!(frame.aperture.balance_residual_j().abs() < 1e-9*scale.max(1e-30));
    }
}

#[test]
fn swapping_upstream_physical_face_reverses_pressure_force_without_retuning_modes() {
    let new = || bound(vec![plate(900., 2e4)], 1e6);
    let mut front = attached(new(), 0, ShellFace::Positive, false);
    let mut back = attached(new(), 0, ShellFace::Negative, false);
    for (a,b) in front.spec().modes.iter().zip(&back.spec().modes) {
        assert_eq!(a.pressure_area_m2, -b.pressure_area_m2);
        assert_eq!(a.mass_kg, b.mass_kg);
        assert_eq!(a.stiffness_n_m, b.stiffness_n_m);
        assert_eq!(a.gap_projection, b.gap_projection);
        assert_eq!(a.force_projection, b.force_projection);
    }
    let a = front.step(drive(5.0), [0.0; 2]).unwrap();
    let b = back.step(drive(5.0), [0.0; 2]).unwrap();
    assert!(a.pressure_drop_pa > 0.0 && b.pressure_drop_pa > 0.0);
    for (a_mode,b_mode) in a.modes.iter().zip(&b.modes) {
        assert!((a_mode.displacement_m/a.pressure_drop_pa
            + b_mode.displacement_m/b.pressure_drop_pa).abs() < 1e-21);
    }
    assert_ne!(a.aperture.midpoint_opening_m, b.aperture.midpoint_opening_m);
}

#[test]
fn public_face_attachment_rebinds_existing_full_c_midpoint_and_exact_mobility() {
    use fs_couple::bernoulli_aperture::paired::BladePropagator;
    let gate = CancelGate::new_clock_free();
    // Manufactured positive definite, coupled mechanical viscosity. The fresh
    // raw-modal construction is an independent constructor-order reference for
    // the attachment's already-cached response, not another physical model.
    let matrix = vec![0.01,0.002,0.0, 0.002,0.02,0.003, 0.0,0.003,0.03];
    for propagator in [BladePropagator::Midpoint, BladePropagator::ExactZoh] {
        let material = ModalAperture::from_material_plates_with_propagator(
            vec![plate(900., 0.0)], 4e-4, super::super::network::gas(), 1e6,
            1.0/f64::from(MECHANICS_RATE), 4096, limits(), vec![BladeState::default(); 3],
            contact(), propagator, &gate).unwrap().with_viscosity(matrix.clone(), &gate).unwrap();
        let mut rebound = attached(material, 0, ShellFace::Negative, true);
        let mut fresh = ModalAperture::new_with_propagator(rebound.spec().clone(),
            rebound.modal_states().to_vec(), contact(), propagator).unwrap()
            .with_viscosity(matrix.clone(), &gate).unwrap();
        assert_eq!(rebound.viscosity_matrix_n_s_m(), fresh.viscosity_matrix_n_s_m());
        for _ in 0..32 {
            let a = rebound.step(drive(5.0), [0.0; 2]).unwrap();
            let b = fresh.step(drive(5.0), [0.0; 2]).unwrap();
            near(a.pressure_drop_pa, b.pressure_drop_pa, 1e-12);
            near(a.aperture.swept_flow_m3_s, b.aperture.swept_flow_m3_s, 1e-12);
            for (a,b) in a.modes.iter().zip(&b.modes) {
                near(a.displacement_m, b.displacement_m, 1e-12);
                near(a.velocity_m_s, b.velocity_m_s, 1e-12);
            }
        }
    }
}

#[test]
fn sheltered_blade_retains_mechanics_and_controls_gap_without_pressure_traction() {
    let valve = bound(vec![plate(900., 2e4)], 1e6);
    let selection = MaterialWettedSurface::from_plate(&valve.material_plates()[0], vec![],
        "manufactured sheltered face".into());
    let mut valve = valve.with_material_wetted_surface(0, selection, 100000,
        &CancelGate::new_clock_free()).unwrap();
    assert_eq!(valve.material_plates()[0].wetted_surface().unwrap().area_m2(), 0.0);
    assert!(valve.spec().modes.iter().all(|mode| mode.pressure_area_m2 == 0.0));
    let a = valve.step(drive(5.0), [0.0; 2]).unwrap();
    assert!(a.modes.iter().all(|state| *state == BladeState::default()));
    assert_eq!(a.aperture.swept_flow_m3_s, 0.0);
    assert!(a.aperture.jet_flow_m3_s > 0.0);
    let b = valve.step(drive(5.0), [1e-5, 0.0]).unwrap();
    assert!(b.modes.iter().any(|state| state.displacement_m != 0.0));
    assert_ne!(b.aperture.midpoint_opening_m, a.aperture.midpoint_opening_m);
    assert_eq!(b.aperture.swept_flow_m3_s, 0.0);
    assert!(b.aperture.mechanical_work_j > 0.0);
}

#[test]
fn foreign_duplicate_unknown_late_and_over_budget_faces_refuse() {
    let gate = CancelGate::new_clock_free();
    let new = || bound(vec![plate(900., 2e4)], 1e6);
    let other = bound(vec![plate(1100., 2e4)], 1e6);
    assert!(new().with_material_wetted_surface(0,
        surface(&other, 0, ShellFace::Positive, false), 100000, &gate).is_err());
    for facets in [vec![OrientedShellFacet { triangle: 24, side: ShellFace::Positive }],
        vec![OrientedShellFacet { triangle: 0, side: ShellFace::Positive },
             OrientedShellFacet { triangle: 0, side: ShellFace::Negative }]] {
        let valve = new();
        let selected = MaterialWettedSurface::from_plate(&valve.material_plates()[0], facets, "test".into());
        assert!(valve.with_material_wetted_surface(0, selected, 100000, &gate).is_err());
    }
    let valve = new();
    let selected = surface(&valve, 0, ShellFace::Positive, false);
    assert!(valve.with_material_wetted_surface(0, selected, 1, &gate).is_err());
    let valve = attached(new(), 0, ShellFace::Positive, false);
    let selected = surface(&valve, 0, ShellFace::Negative, false);
    assert!(valve.with_material_wetted_surface(0, selected, 100000, &gate).is_err());
    let mut valve = new();
    valve.step(drive(0.0), [0.0; 2]).unwrap();
    let selected = surface(&valve, 0, ShellFace::Positive, false);
    assert!(valve.with_material_wetted_surface(0, selected, 100000, &gate).is_err());
    let valve = new();
    let selected = surface(&valve, 0, ShellFace::Positive, false);
    gate.request();
    assert!(matches!(valve.with_material_wetted_surface(0, selected, 100000, &gate),
        Err(fs_couple::acoustic_realize::AcousticRealizeError::Cancelled)));
}
