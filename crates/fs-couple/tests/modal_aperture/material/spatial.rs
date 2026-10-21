use super::*;
use fs_couple::bernoulli_aperture::paired::BladePropagator;
use fs_couple::bernoulli_aperture::plate::closure::PlateClosureSpec;

fn spatial(propagator: BladePropagator, pair: bool, budget: usize) -> ModalAperture {
    let mut plates = vec![plate(900.0, 0.0)];
    if pair {
        plates.push(plate(1150.0, 0.0));
    }
    let mesh = &plates[0].0.chart().mesh;
    let closure = PlateClosureSpec {
        nodal_rest_gap_m: mesh
            .nodes
            .iter()
            .map(|p| 4e-6 + (p.1 - 0.005) * 0.001)
            .collect(),
        lay_triangles: (0..mesh.tris.len()).collect(),
        stiffness_pa_per_m_alpha: 1e11,
        alpha: 2.0,
        internal_loss_s_per_m: 2.0,
        provenance: "authored source-mesh modal contact regression".into(),
        max_penetration_m: 0.0005,
    };
    let n = plates.len() * 3;
    let mut limits = limits();
    limits.max_evaluations = budget;
    let mut viscosity = vec![0.0; n * n];
    for i in 0..n {
        viscosity[i * n + i] = 0.03;
        if i + 1 < n {
            viscosity[i * n + i + 1] = 0.004;
            viscosity[(i + 1) * n + i] = 0.004;
        }
    }
    let mut a = ModalAperture::from_material_plates_with_propagator(
        plates,
        4e-6,
        super::super::network::gas(),
        1e6,
        1.0 / f64::from(MECHANICS_RATE),
        512,
        limits,
        vec![BladeState::default(); n],
        contact(),
        propagator,
        &CancelGate::new_clock_free(),
    )
    .unwrap()
    .with_viscosity(viscosity, &CancelGate::new_clock_free())
    .unwrap()
    .with_plate_closure(closure)
    .unwrap();
    a.set_material_viscous_slit(0.001).unwrap();
    a
}

#[test]
fn original_material_columns_drive_distributed_contact_in_both_structural_images() {
    for propagator in [BladePropagator::Midpoint, BladePropagator::ExactZoh] {
        for pair in [false, true] {
            let mut a = spatial(propagator, pair, 512);
            let closure = a.plate_closure().unwrap();
            for plate in a.material_plates() {
                for (all, local) in closure
                    .displacement_shapes()
                    .iter()
                    .zip(plate.displacement_shapes())
                {
                    assert_eq!(&all[plate.mode_range()], local.as_slice());
                }
            }
            let before = closure.probe(&vec![0.0; a.spec().modes.len()]).unwrap();
            assert!(before.active_lay_points > 0);
            assert!(before.open_area_m2 > 0.0);
            let mut contacted = false;
            for sample in 0..64 {
                let f = a
                    .step(drive(if sample < 32 { 20.0 } else { -5.0 }), [0.0; 2])
                    .unwrap();
                contacted |= f.contact_force_n > 0.0;
                assert!(f.work_certificate.residual <= f.work_certificate.effective_budget);
                assert!(f.flow_certificate.residual <= f.flow_certificate.effective_budget);
                assert!(f.momentum_residual_n.iter().all(|f| f.abs() < 1e-7));
                assert!(f.aperture.dissipated_energy_j >= 0.0);
            }
            assert!(contacted);
            assert_eq!(a.accepted_steps(), 64);
            assert_eq!(a.modal_states().len(), if pair { 6 } else { 3 });
        }
    }
}

#[test]
fn nonlinear_budget_refusal_does_not_publish_modal_or_contact_state() {
    let mut a = spatial(BladePropagator::ExactZoh, true, 1);
    let states = a.modal_states().to_vec();
    let energy = a.stored_energy_j();
    let q: Vec<_> = states.iter().map(|s| s.displacement_m).collect();
    let contact = a.plate_closure().unwrap().probe(&q).unwrap();
    assert!(a.step(drive(20.0), [0.0; 2]).is_err());
    assert_eq!(a.accepted_steps(), 0);
    assert_eq!(a.modal_states(), states);
    assert_eq!(a.stored_energy_j(), energy);
    assert_eq!(a.plate_closure().unwrap().probe(&q).unwrap(), contact);
}
