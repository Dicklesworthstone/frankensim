//! A sourced hereditary solid remains in every original mode and acoustic port.
use super::*;
use crate::bernoulli_aperture::dynamic::ApertureDrive;
use crate::bernoulli_aperture::modal::relaxation::InitialModalMemory;
use fs_dcontact::Obstacle;

mod regional;
mod tensor;

const TENSOR: &str = "loss tensor eta11 eta22 eta66 eta12 eta16 eta26 1 1000000";
const MAXWELL: &str = "loss proportional-maxwell relaxed 1 boundary_modulus boundary_time";

fn source() -> String {
    input()
        .replace(TENSOR, MAXWELL)
        .replace("audio 768000 48000 512", "audio 768000 48000 256")
}
fn load(text: &str, block: usize, modulus: f64, time: f64) -> MaterialModalValvePerformance {
    MaterialModalValvePerformance::from_bytes_with_sources(
        text.as_bytes(),
        block,
        &CancelGate::new_clock_free(),
        None,
        |_| Ok(pack_with_relaxation(900.0, false, modulus, time)),
    )
    .unwrap()
}
fn bare(source: &ModalAperture, displacement_m: f64) -> ModalAperture {
    let spec = source.spec();
    let state = source
        .modal_states()
        .iter()
        .enumerate()
        .map(|(i, _)| BladeState {
            displacement_m: displacement_m / (i + 1) as f64,
            velocity_m_s: 0.0,
        })
        .collect();
    ModalAperture::from_material_plates_with_propagator(
        source
            .material_plates()
            .iter()
            .map(|p| (p.source().clone(), p.options().clone()))
            .collect(),
        spec.rest_opening_m,
        *source.gas().unwrap(),
        spec.impedance_pa_s_m3,
        spec.time_step_s,
        spec.max_steps,
        spec.limits,
        state,
        Obstacle::new(
            vec![-1.0],
            1,
            1,
            vec![0.0],
            vec![1.0],
            1e8,
            2.0,
            "manufactured multimode material-memory contact".into(),
        )
        .unwrap(),
        BladePropagator::ExactZoh,
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn properties() -> Vec<(String, String)> {
    vec![("boundary_modulus".into(), "boundary_time".into())]
}

#[test]
fn g1_independent_modal_memory_storage_loss_work_and_refused_step() {
    let text = source().replace(
        MAXWELL,
        "loss proportional-maxwell relaxed 2 boundary_modulus boundary_time fast_modulus fast_time",
    );
    let source = load(&text, 37, 1e8, 0.001);
    let mut selected = properties();
    selected.push(("fast_modulus".into(), "fast_time".into()));
    let mut reed = bare(valve(&source), 1e-7)
        .with_material_relaxation(
            &selected,
            InitialModalMemory::Unrelaxed,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let memory = reed.material_relaxation().unwrap();
    assert_eq!(memory.material_spectra().len(), 2);
    assert!(
        memory
            .material_spectra()
            .iter()
            .all(|regions| regions.len() == 2)
    );
    assert_eq!(memory.memory_sqrt_j(), &vec![vec![0.0; 2]; 6]);
    let initial: f64 = reed
        .spec()
        .modes
        .iter()
        .zip(reed.modal_states())
        .map(|(mode, state)| {
            0.5 * mode.stiffness_n_m * state.displacement_m.powi(2) * (1.0 + (1e8 + 3e7) / 4e9)
        })
        .sum();
    assert!((initial - reed.stored_energy_j()).abs() <= initial * 1e-14);
    let (mut total_loss, mut total_work, mut total_dissipation) = (0.0, 0.0, 0.0);
    for sample in 0..2048 {
        let forces = if sample < 256 {
            [1e-4, -7e-5]
        } else {
            [0.0; 2]
        };
        let before = reed.modal_states().to_vec();
        let frame = reed
            .step(
                ApertureDrive {
                    upstream_pressure_pa: if sample < 512 { 5.0 } else { 0.0 },
                    ..ApertureDrive::default()
                },
                forces,
            )
            .unwrap();
        let independent_work: f64 = reed
            .spec()
            .modes
            .iter()
            .zip(&before)
            .zip(&frame.modes)
            .map(|((mode, before), after)| {
                (mode.force_projection[0] * forces[0] + mode.force_projection[1] * forces[1])
                    * (after.displacement_m - before.displacement_m)
            })
            .sum();
        assert!((independent_work - frame.aperture.mechanical_work_j).abs() < 1e-22);
        assert!(frame.aperture.relaxation_loss_j >= -1e-30);
        assert!(
            frame.aperture.balance_residual_j().abs() <= frame.work_certificate.effective_budget
        );
        assert!(frame.momentum_residual_n.iter().all(|r| r.abs() < 1e-10));
        assert_eq!(reed.stored_energy_j(), frame.aperture.stored_energy_j);
        total_loss += frame.aperture.relaxation_loss_j;
        total_work += frame.aperture.pressure_work_j + independent_work;
        total_dissipation += frame.aperture.dissipated_energy_j;
    }
    assert!(total_loss > 0.0);
    assert!(
        reed.material_relaxation()
            .unwrap()
            .memory_sqrt_j()
            .iter()
            .flatten()
            .all(|z| *z != 0.0)
    );
    let balance = reed.stored_energy_j() - initial + total_dissipation - total_work;
    assert!(balance.abs() < 1e-11 * (initial + total_dissipation + total_work.abs()));
    let accepted = (
        reed.accepted_steps(),
        reed.modal_states().to_vec(),
        reed.material_relaxation().unwrap().memory_sqrt_j().to_vec(),
        reed.stored_energy_j(),
    );
    assert!(
        reed.step(
            ApertureDrive {
                upstream_pressure_pa: f64::NAN,
                ..ApertureDrive::default()
            },
            [0.0; 2]
        )
        .is_err()
    );
    assert_eq!(
        accepted,
        (
            reed.accepted_steps(),
            reed.modal_states().to_vec(),
            reed.material_relaxation().unwrap().memory_sqrt_j().to_vec(),
            reed.stored_energy_j()
        )
    );
}

#[test]
fn g3_material_spectrum_changes_returning_acoustics_and_microphone_without_changing_geometry() {
    let text = source();
    let original = load(&text, 37, 1e8, 0.001);
    let changed = load(&text, 37, 4e8, 0.0002);
    assert_eq!(valve(&original).spec(), valve(&changed).spec());
    assert_ne!(
        valve(&original).material_plates()[0]
            .source()
            .specimen_identity(),
        valve(&changed).material_plates()[0]
            .source()
            .specimen_identity()
    );
    let (pressure, source) = render(original);
    let (other, changed) = render(changed);
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(
        pressure
            .iter()
            .zip(&other)
            .any(|(a, b)| (a - b).abs() > 1e-15)
    );
    assert_ne!(
        source
            .source()
            .system()
            .node_frame(ValveSide::Downstream, 0),
        changed
            .source()
            .system()
            .node_frame(ValveSide::Downstream, 0)
    );
    let LoadedAperture::Modal(reed) = source.source().system().valve() else {
        panic!("original modes");
    };
    assert_eq!(reed.material_relaxation().unwrap().memory_sqrt_j().len(), 6);
    assert!(
        reed.material_relaxation()
            .unwrap()
            .memory_sqrt_j()
            .iter()
            .all(|m| m.len() == 1)
    );
    let (partitioned, _) = render(load(&text, 1, 1e8, 0.001));
    assert_eq!(pressure, partitioned);
}

#[test]
fn g1_zero_source_branch_reproduces_elastic_image_and_initial_memory_is_explicit() {
    let text = source();
    let elastic = load(&text.replace(MAXWELL, "loss elastic"), 37, 0.0, 0.001);
    let zero = load(&text, 37, 0.0, 0.001);
    assert_eq!(valve(&elastic).spec(), valve(&zero).spec());
    assert!(
        valve(&zero)
            .material_relaxation()
            .unwrap()
            .memory_sqrt_j()
            .iter()
            .all(Vec::is_empty)
    );
    assert_eq!(render(elastic).0, render(zero).0);
    let source = load(&text, 37, 1e8, 0.001);
    let relaxed = bare(valve(&source), 1e-7)
        .with_material_relaxation(
            &properties(),
            InitialModalMemory::Relaxed,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let unrelaxed = bare(valve(&source), 1e-7)
        .with_material_relaxation(
            &properties(),
            InitialModalMemory::Unrelaxed,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert!(unrelaxed.stored_energy_j() > relaxed.stored_energy_j());
    assert!(
        relaxed
            .material_relaxation()
            .unwrap()
            .memory_sqrt_j()
            .iter()
            .flatten()
            .all(|z| *z != 0.0)
    );
}

#[test]
fn g0_unsupported_material_laws_and_numerical_bands_refuse() {
    let text = source();
    for invalid in [
        text.replace("propagator exact-zoh", "propagator midpoint"),
        text.replacen("support clamped 0", "support clamped 0.00001", 1),
        text.replace(
            MAXWELL,
            "loss proportional-maxwell relaxed 1 missing_modulus boundary_time",
        ),
        text.replace(
            MAXWELL,
            "loss proportional-maxwell implicit 1 boundary_modulus boundary_time",
        ),
        text.replace(MAXWELL, "loss proportional-maxwell relaxed 0"),
    ] {
        assert!(
            MaterialModalValvePerformance::from_bytes_with_sources(
                invalid.as_bytes(),
                37,
                &CancelGate::new_clock_free(),
                None,
                |_| Ok(pack(900.0, false))
            )
            .is_err()
        );
    }
    for (modulus, time) in [(1e8, 1e-10), (1e14, 0.001)] {
        assert!(
            MaterialModalValvePerformance::from_bytes_with_sources(
                text.as_bytes(),
                37,
                &CancelGate::new_clock_free(),
                None,
                |_| Ok(pack_with_relaxation(900.0, false, modulus, time))
            )
            .is_err()
        );
    }
    // Only one original region changes: projecting its branch as a diagonal
    // modal decay would silently erase real cross-mode constitutive forces.
    let nonproportional = text.replacen("section blade.matpack", "section altered.matpack", 1);
    let error = MaterialModalValvePerformance::from_bytes_with_sources(
        nonproportional.as_bytes(),
        37,
        &CancelGate::new_clock_free(),
        None,
        |path| {
            Ok(pack_with_relaxation(
                900.0,
                false,
                if path == "altered.matpack" { 2e8 } else { 1e8 },
                0.001,
            ))
        },
    )
    .err()
    .unwrap();
    assert!(error.to_string().contains("nonproportional"));
}
