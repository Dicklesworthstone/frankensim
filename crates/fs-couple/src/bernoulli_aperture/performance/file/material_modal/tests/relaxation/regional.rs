//! Different original regions retain their complete constitutive cross terms.
use super::*;

const REGIONAL: &str = "loss regional-maxwell relaxed 1 boundary_modulus boundary_time";

fn regional_source() -> String {
    source().replace(MAXWELL, REGIONAL).replacen(
        "section blade.matpack",
        "section altered.matpack",
        1,
    )
}

fn load_regional(
    text: &str,
    block: usize,
    modulus: f64,
    time: f64,
) -> MaterialModalValvePerformance {
    MaterialModalValvePerformance::from_bytes_with_sources(
        text.as_bytes(),
        block,
        &CancelGate::new_clock_free(),
        None,
        |path| {
            let (modulus, time) = if path == "altered.matpack" {
                (modulus, time)
            } else {
                (1e8, 0.001)
            };
            Ok(pack_with_relaxation(900.0, false, modulus, time))
        },
    )
    .unwrap()
}

#[test]
fn regional_rows_keep_cross_mode_storage_work_loss_and_refused_history() {
    let source = load_regional(&regional_source(), 37, 4e8, 0.0002);
    let bare = bare(valve(&source), 1e-7);
    let equilibrium = bare.stored_energy_j();
    let mut reed = bare
        .with_regional_material_relaxation(
            &properties(),
            InitialModalMemory::Unrelaxed,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let memory = reed.material_relaxation().unwrap();
    assert!(memory.is_regional());
    assert_eq!(memory.memory_sqrt_j(), &vec![vec![0.0; 12]]);
    assert_eq!(memory.material_spectra()[0][0].terms, vec![(4e8, 0.0002)]);
    assert_eq!(memory.material_spectra()[0][1].terms, vec![(1e8, 0.001)]);
    assert_eq!(memory.material_spectra()[1][0].terms, vec![(1e8, 0.001)]);
    assert!(memory.regional_coupling_row(12).is_none());
    assert!(memory.regional_relaxation_time_s(12).is_none());
    let n = reed.modal_states().len();
    let mut stiffness = vec![0.0; n * n];
    for channel in 0..12 {
        let row = memory.regional_coupling_row(channel).unwrap();
        assert_eq!(row.len(), n);
        assert_eq!(
            memory.regional_relaxation_time_s(channel),
            Some(if channel < 3 { 0.0002 } else { 0.001 })
        );
        let other_blade = if channel < 6 { 3..6 } else { 0..3 };
        assert!(row[other_blade].iter().all(|v| *v == 0.0));
        for i in 0..n {
            for j in 0..n {
                stiffness[n * i + j] += row[i] * row[j];
            }
        }
    }
    assert!((0..3).any(|i| (i + 1..3).any(|j| stiffness[n * i + j].abs()
        > 1e-8 * (stiffness[n * i + i] * stiffness[n * j + j]).sqrt())));
    // A physical quadratic form includes BOTH signed cross terms; a sum of
    // fitted scalar modal losses cannot reproduce this initial stored energy.
    let mut additional = 0.0;
    for i in 0..n {
        for j in 0..n {
            additional += 0.5
                * reed.modal_states()[i].displacement_m
                * stiffness[n * i + j]
                * reed.modal_states()[j].displacement_m;
        }
    }
    let initial = reed.stored_energy_j();
    assert!((initial - equilibrium - additional).abs() <= 1e-12 * initial);
    let (mut loss, mut work, mut dissipation) = (0.0, 0.0, 0.0);
    for sample in 0..512 {
        let before = reed.modal_states().to_vec();
        let force = if sample < 128 {
            [1e-4, -7e-5]
        } else {
            [0.0; 2]
        };
        let frame = reed
            .step(
                ApertureDrive {
                    upstream_pressure_pa: if sample < 256 { 5.0 } else { 0.0 },
                    ..ApertureDrive::default()
                },
                force,
            )
            .unwrap();
        let mechanical: f64 = reed
            .spec()
            .modes
            .iter()
            .zip(&before)
            .zip(&frame.modes)
            .map(|((m, a), b)| {
                (m.force_projection[0] * force[0] + m.force_projection[1] * force[1])
                    * (b.displacement_m - a.displacement_m)
            })
            .sum();
        assert!((mechanical - frame.aperture.mechanical_work_j).abs() < 1e-22);
        assert!(frame.aperture.relaxation_loss_j >= -1e-30);
        assert!(frame.momentum_residual_n.iter().all(|v| v.abs() < 1e-10));
        assert!(
            frame.aperture.balance_residual_j().abs() <= frame.work_certificate.effective_budget
        );
        assert_eq!(reed.stored_energy_j(), frame.aperture.stored_energy_j);
        loss += frame.aperture.relaxation_loss_j;
        work += mechanical + frame.aperture.pressure_work_j;
        dissipation += frame.aperture.dissipated_energy_j;
    }
    assert!(loss > 0.0);
    assert!(
        (reed.stored_energy_j() - initial + dissipation - work).abs()
            < 1e-10 * (initial + dissipation + work.abs())
    );
    let accepted = (
        reed.accepted_steps(),
        reed.modal_states().to_vec(),
        reed.material_relaxation().unwrap().memory_sqrt_j().to_vec(),
        reed.stored_energy_j(),
    );
    assert!(accepted.2[0].iter().any(|z| *z != 0.0));
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
fn uniform_regional_spectra_recover_proportional_physics_without_refitting_modes() {
    let source = load(&source(), 37, 1e8, 0.001);
    let mut proportional = bare(valve(&source), 1e-7)
        .with_material_relaxation(
            &properties(),
            InitialModalMemory::Unrelaxed,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let mut regional = bare(valve(&source), 1e-7)
        .with_regional_material_relaxation(
            &properties(),
            InitialModalMemory::Unrelaxed,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(proportional.spec(), regional.spec());
    assert!(!proportional.material_relaxation().unwrap().is_regional());
    assert!(regional.material_relaxation().unwrap().is_regional());
    // Uniform excess stiffness is known independently of the new factor rows:
    // all physical bending energy scales by the original E_branch/E_equilibrium.
    let expected: f64 = regional
        .spec()
        .modes
        .iter()
        .zip(regional.modal_states())
        .map(|(m, s)| 0.5 * m.stiffness_n_m * s.displacement_m.powi(2) * (1.0 + 1e8 / 4e9))
        .sum();
    assert!((regional.stored_energy_j() - expected).abs() < 1e-7 * expected);
    let (mut a_loss, mut b_loss) = (0.0, 0.0);
    for sample in 0..256 {
        let drive = ApertureDrive {
            upstream_pressure_pa: if sample < 128 { 5.0 } else { 0.0 },
            ..ApertureDrive::default()
        };
        let force = if sample < 64 { [1e-4, -7e-5] } else { [0.0; 2] };
        let a = proportional.step(drive, force).unwrap();
        let b = regional.step(drive, force).unwrap();
        for (a, b) in a.modes.iter().zip(&b.modes) {
            assert!(
                (a.displacement_m - b.displacement_m).abs()
                    < 1e-7 * a.displacement_m.abs().max(b.displacement_m.abs()).max(1e-7)
            );
            assert!(
                (a.velocity_m_s - b.velocity_m_s).abs()
                    < 1e-7 * a.velocity_m_s.abs().max(b.velocity_m_s.abs()).max(1e-4)
            );
        }
        assert!(
            (a.pressure_drop_pa - b.pressure_drop_pa).abs()
                < 1e-7
                    * a.pressure_drop_pa
                        .abs()
                        .max(b.pressure_drop_pa.abs())
                        .max(1.0)
        );
        a_loss += a.aperture.relaxation_loss_j;
        b_loss += b.aperture.relaxation_loss_j;
    }
    assert!((a_loss - b_loss).abs() < 1e-7 * a_loss.abs().max(b_loss.abs()));
}

#[test]
fn one_region_material_swap_changes_returning_acoustics_and_block_invariant_pressure() {
    let text = regional_source();
    let original = load_regional(&text, 37, 1e8, 0.001);
    let changed = load_regional(&text, 37, 1e9, 0.0002);
    assert_eq!(valve(&original).spec(), valve(&changed).spec());
    for (a, b) in valve(&original)
        .material_plates()
        .iter()
        .zip(valve(&changed).material_plates())
    {
        assert_eq!(a.displacement_shapes(), b.displacement_shapes());
        assert_eq!(a.mode_range(), b.mode_range());
    }
    assert_ne!(
        valve(&original).material_plates()[0]
            .source()
            .specimen_identity(),
        valve(&changed).material_plates()[0]
            .source()
            .specimen_identity()
    );
    assert_eq!(
        valve(&original).material_plates()[1]
            .source()
            .specimen_identity(),
        valve(&changed).material_plates()[1]
            .source()
            .specimen_identity()
    );
    let (reference, unchanged_source) = render(original);
    let (pressure, changed) = render(changed);
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(
        reference
            .iter()
            .zip(&pressure)
            .any(|(a, b)| (a - b).abs() > 1e-15)
    );
    let (partitioned, partitioned_source) = render(load_regional(&text, 1, 1e9, 0.0002));
    assert_eq!(pressure, partitioned);
    assert_ne!(
        changed
            .source()
            .system()
            .node_frame(ValveSide::Downstream, 0),
        unchanged_source
            .source()
            .system()
            .node_frame(ValveSide::Downstream, 0)
    );
    assert_eq!(
        changed
            .source()
            .system()
            .node_frame(ValveSide::Downstream, 0),
        partitioned_source
            .source()
            .system()
            .node_frame(ValveSide::Downstream, 0)
    );
}

#[test]
fn localized_material_patch_reaches_modal_memory_and_the_actual_microphone() {
    // A patch adjacent to the clamped edge has only one unconstrained DKT
    // vertex (three physical coordinates). Its bending form is consequently
    // rank deficient on SIX retained modes. Such a local material is not an
    // invalid whole instrument, and its modal cross terms must remain active.
    let text = regional_source()
        .replace("blade 3 0.2 3000 0.04", "blade 6 0.2 12000 0.2")
        .replace("initial 3\nstate 0 0\nstate 0 0\nstate 0 0",
            "initial 6\nstate 0 0\nstate 0 0\nstate 0 0\nstate 0 0\nstate 0 0\nstate 0 0");
    let text = text.lines().map(|line| {
        if let Some(row) = line.strip_prefix("triangle ") {
            let fields: Vec<_> = row.split_whitespace().collect();
            let region = usize::from(fields[..3] != ["0", "8", "7"]);
            format!("triangle {} {} {} {region}", fields[0], fields[1], fields[2])
        } else { line.to_owned() }
    }).collect::<Vec<_>>().join("\n") + "\n";
    let reference = load_regional(&text, 37, 1e8, 0.001);
    let changed = load_regional(&text, 37, 6e8, 0.0002);
    let a = valve(&changed);
    assert_eq!(a.modal_states().len(), 12);
    let plate = &a.material_plates()[0];
    assert_eq!(plate.source().regions()[0].input.region.triangle_indices, vec![1]);
    let memory = a.material_relaxation().unwrap();
    assert!(memory.is_regional());
    assert_eq!(memory.material_spectra()[0][0].terms, vec![(6e8, 0.0002)]);
    let patch_channels: Vec<_> = (0..memory.memory_sqrt_j()[0].len())
        .filter(|&i| memory.regional_relaxation_time_s(i) == Some(0.0002)).collect();
    assert!((1..=3).contains(&patch_channels.len()),
        "one free physical vertex must not acquire six independent strain channels");
    assert!(patch_channels.iter().any(|&i|
        memory.regional_coupling_row(i).unwrap()[..6].iter().filter(|v| **v != 0.0).count() > 1));
    let (original_pressure, _) = render(reference);
    let (pressure, rendered) = render(changed);
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(pressure.iter().any(|p| p.abs() > 1e-15));
    assert!(pressure.iter().zip(original_pressure).any(|(a, b)| (a - b).abs() > 1e-15));
    let (partitioned, _) = render(load_regional(&text, 1, 6e8, 0.0002));
    assert_eq!(pressure, partitioned);
    let LoadedAperture::Modal(a) = rendered.source().system().valve() else { panic!("modal source") };
    assert!(a.material_relaxation().unwrap().memory_sqrt_j().iter().flatten().any(|z| *z != 0.0));
}

#[test]
fn regional_selection_refuses_missing_sources_unresolved_clocks_and_complete_history_overflow() {
    let text = regional_source();
    for invalid in [
        text.replace("propagator exact-zoh", "propagator midpoint"),
        text.replace(
            REGIONAL,
            "loss regional-maxwell implicit 1 boundary_modulus boundary_time",
        ),
        text.replace(REGIONAL, "loss regional-maxwell relaxed 0"),
        text.replace(
            REGIONAL,
            "loss regional-kelvin relaxed 1 boundary_modulus boundary_time",
        ),
        text.replace(
            "boundary_modulus boundary_time",
            "missing_modulus boundary_time",
        ),
        text.replace("boundary_modulus boundary_time", "density boundary_time"),
        text.replace("isotropic 293.15 100", "isotropic 330 100"),
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
    for (modulus, time) in [(-1.0, 0.001), (1e8, 1e-10), (1e14, 0.001)] {
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
    // Three finite-width regions per blade retain three bending channels each.
    // Eight physical terms need 144 histories: the full exact loss block admits
    // only 110 for these six retained mechanical modes. Tiny moduli keep this
    // a state-capacity refusal rather than an instantaneous-frequency refusal.
    let crowded = text.replace("sections 2","sections 3")
        .replace("section blade.matpack isotropic 293.15 100 0.00027 0",
            "section blade.matpack isotropic 293.15 100 0.00027 0\nsection blade.matpack isotropic 293.15 100 0.00027 0")
        .replace(REGIONAL,&format!("loss regional-maxwell relaxed 8{}",
            " boundary_modulus boundary_time".repeat(8)));
    let mut triangle = 0;
    let crowded = crowded
        .lines()
        .map(|line| {
            if let Some(prefix) = line.strip_prefix("triangle ") {
                let mut fields: Vec<_> = prefix.split_whitespace().collect();
                fields.pop();
                let result = format!("triangle {} {}", fields.join(" "), (triangle % 12) / 4);
                triangle += 1;
                result
            } else {
                line.to_string()
            }
        })
        .collect::<Vec<_>>()
        // lines() drops the final LF that canonical gesture bytes end with.
        .join("\n") + "\n";
    let error = MaterialModalValvePerformance::from_bytes_with_sources(
        crowded.as_bytes(),
        37,
        &CancelGate::new_clock_free(),
        None,
        |_| Ok(pack_with_relaxation(900.0, false, 1e3, 0.001)),
    )
    .err()
    .unwrap();
    assert!(error.to_string().contains("capacity"), "{error}");
}
