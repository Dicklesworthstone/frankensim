use super::*;
use fs_exec::CancelGate;
use fs_plate::shell::reduction::ReductionBudget;
use fs_plate::shell::{ShellMesh, ShellSupport, assemble_shell_sections};
use fs_plate::{ModePair, PlateSection};

fn shell() -> ShellReduction {
    let mesh = ShellMesh::new(
        vec![
            [0., 0., 0.],
            [0.1, 0., 0.01],
            [0., 0.1, 0.02],
            [0.08, 0.09, 0.035],
        ],
        vec![[0, 1, 3], [0, 3, 2]],
    )
    .unwrap();
    let sections = vec![
        PlateSection::isotropic(3.0e6, 0.3, 0.001, 1100.).unwrap(),
        PlateSection::isotropic(3.0e6, 0.3, 0.0007, 1100.).unwrap(),
    ];
    let model =
        assemble_shell_sections(&mesh, &sections, &[0, 1, 2], ShellSupport::Clamped).unwrap();
    let n = model.free;
    let mut k = vec![0.; n * n];
    let mut mass = k.clone();
    for i in 0..n {
        for j in 0..n {
            k[i * n + j] = model.k.get(i, j);
            mass[i * n + j] = model.m.get(i, j);
        }
    }
    let mut all = fs_modal::eigh_gen_dense(&k, &mass, n).unwrap();
    let mut normal_integral = [0.0; 3];
    for f in 0..mesh.tris.len() {
        let g = mesh.facet(f).unwrap();
        for a in 0..3 {
            normal_integral[a] += g.frame[2][a] * g.area_m2;
        }
    }
    let trace = |mode: &ModePair| -> f64 {
        (0..3)
            .map(|a| normal_integral[a] * model.dof_map[18 + a].map_or(0.0, |i| mode.phi[i]))
            .sum::<f64>()
            .abs()
    };
    // Select two of the actual numerical eigenvectors with resolved transverse
    // participation. The reducer still checks their mass/pencil normalization.
    all.sort_by(|a, b| trace(b).total_cmp(&trace(a)));
    assert!(trace(&all[1]) > 1e-10);
    let modes: Vec<ModePair> = all.into_iter().take(2).collect();
    ShellReduction::new(
        &mesh,
        &sections,
        &model,
        &modes,
        ReductionBudget {
            max_modes: 8,
            max_facet_modes: 100,
            relative_tolerance: 1e-6,
        },
    )
    .unwrap()
}
fn source_options() -> OpenScreenModesOptions {
    OpenScreenModesOptions {
        maximum_modes: 8,
        source_surface: OpenScreenRadiationBudget {
            maximum_vertices: 256,
            maximum_facets: 128,
            maximum_corner_weights: 768,
        },
        acoustic_refinements: 1,
        maximum_source_bytes: 2_000_000,
        maximum_setup_work: 10_000_000,
    }
}
fn source() -> OpenScreenLoadedModes {
    let damping = ModalPlateViscosity::from_psd_modal_matrix(
        vec![0.8, 0.3, 0.3, 0.7],
        2,
        1.,
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    OpenScreenLoadedModes::from_shell(
        &shell(),
        Medium {
            density: 1000.,
            sound_speed: 1500.,
        },
        Some(&damping),
        source_options(),
        || false,
    )
    .unwrap()
}
fn boundary_options() -> OpenOptions {
    OpenOptions {
        edge_level: 0,
        maximum_edge_level: 0,
        triangle_order: 4,
        angular_refinement: 0,
        maximum_unknowns: 256,
        maximum_cache_bytes: 4_000_000,
        maximum_work: 1_000_000_000,
        workers: 1,
    }
}
fn options() -> OpenScreenHarmonicOptions {
    OpenScreenHarmonicOptions {
        boundary: ScreenPortLimits {
            maximum_ports: 8,
            maximum_matrix_entries: 64,
            maximum_jump_entries: 2048,
            maximum_work: 1_000_000_000,
        },
        maximum_load_cases: 8,
        maximum_response_bytes: 100_000,
        maximum_work: 1_000_000,
        physical_residual_tolerance: 1e-10,
        power_relative_tolerance: 1e-10,
    }
}
fn close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()) + 1e-20,
        "{a:?} != {b:?}"
    );
}
fn matrix(source: &OpenScreenLoadedModes, omega: f64, z: &[C64]) -> [C64; 4] {
    core::array::from_fn(|at| {
        let mass = if at == 0 || at == 3 {
            omega * omega
        } else {
            0.
        };
        C64::new(source.stiffness[at] - mass, -omega * source.viscosity[at])
            + C64::new(0., -omega) * z[at]
    })
}
fn invert(a: [C64; 4], f: [C64; 2]) -> [C64; 2] {
    let determinant = a[0] * a[3] - a[1] * a[2];
    [
        (a[3] * f[0] - a[1] * f[1]) / determinant,
        (a[0] * f[1] - a[2] * f[0]) / determinant,
    ]
}

#[test]
fn g3_actual_shell_reaction_keeps_full_source_k_c_and_mutual_diffraction() {
    let original = shell();
    let source = source();
    let mut k = vec![0.; 4];
    original.hessian_matrix(&[0.; 2], &mut k, &mut [[0.; 3]; 2]);
    assert_eq!(source.stiffness_per_s2(), k);
    assert_eq!(source.viscosity_per_s(), [0.8, 0.3, 0.3, 0.7]);
    let prepared = source.prepare(800., boundary_options(), || false).unwrap();
    let omega = 450.;
    let force = [C64::new(0.5, 0.2), C64::new(-0.1, 0.3)];
    let batch = prepared
        .harmonic_response_batch(omega, &[&force], options(), || false)
        .unwrap();
    let r = &batch.responses()[0];
    let z = batch.radiation().impedance();
    assert!(z[1].abs() > 0. && z[2].abs() > 0.);
    let expected = invert(matrix(&source, omega, z), force);
    for i in 0..2 {
        close(r.modal_displacement[i], expected[i], 5e-12);
        close(
            r.modal_velocity[i],
            C64::new(0., -omega) * expected[i],
            5e-12,
        );
        close(
            r.fluid_reaction[i],
            -(z[2 * i] * r.modal_velocity[0] + z[2 * i + 1] * r.modal_velocity[1]),
            5e-12,
        );
    }
    let diagonal = [z[0], C64::ZERO, C64::ZERO, z[3]];
    let omitted = invert(matrix(&source, omega, &diagonal), force);
    assert!(
        (omitted[0] - expected[0]).abs() + (omitted[1] - expected[1]).abs()
            > 1e-8 * (expected[0].abs() + expected[1].abs()),
        "self-only radiation must not reproduce mutual loading"
    );
    let input =
        0.5 * (force[0] * r.modal_velocity[0].conj() + force[1] * r.modal_velocity[1].conj()).re;
    let screen = batch.radiation().input_power_w(&r.modal_velocity).unwrap();
    assert!((r.power.input_w - input).abs() < 1e-11 * input.abs() + 1e-20);
    assert!((r.power.radiated_w - screen).abs() < 1e-11 * screen.abs() + 1e-20);
    assert!(r.power.kelvin_w > 0. && screen >= 0.);
    assert!(r.power.relative_balance < 1e-10 && r.mechanical_residual < 1e-10);
    assert!(r.mean_structural_energy_j > 0.);
}

#[test]
fn g3_shared_batch_observes_actual_loaded_coherent_jumps_and_both_sides() {
    let source = source();
    let prepared = source.prepare(800., boundary_options(), || false).unwrap();
    let f = [C64::new(0.7, -0.2), C64::ZERO];
    let g = [C64::ZERO, C64::new(-0.1, 0.8)];
    let total = [f[0], g[1]];
    let zero = [C64::ZERO; 2];
    let batch = prepared
        .harmonic_response_batch(450., &[&f, &g, &total, &zero], options(), || false)
        .unwrap();
    for mode in 0..2 {
        close(
            batch.responses()[2].modal_velocity[mode],
            batch.responses()[0].modal_velocity[mode] + batch.responses()[1].modal_velocity[mode],
            1e-12,
        );
    }
    let jumps: Vec<_> = (0..4)
        .map(|i| batch.pressure_jump(i, || false).unwrap())
        .collect();
    for i in 0..jumps[0].len() {
        close(jumps[2][i], jumps[0][i] + jumps[1][i], 1e-12);
        assert_eq!(jumps[3][i], C64::ZERO);
        if source.surface().is_free_edge_vertex(i) {
            assert_eq!(jumps[2][i], C64::ZERO);
        }
    }
    for direction in [[0., 0., 1.], [0., 0., -1.]] {
        let expected = batch.far_field(0, direction, || false).unwrap()
            + batch.far_field(1, direction, || false).unwrap();
        close(
            batch.far_field(2, direction, || false).unwrap(),
            expected,
            2e-12,
        );
    }
    let points = [[0.04, 0.04, 0.3], [0.04, 0.04, -0.3]];
    let p = batch
        .pressure_at(2, &points, 0.01, near_field::Options::default(), || false)
        .unwrap();
    assert!(p.iter().all(|z| finite(*z) && z.abs() > 0.));
    let outgoing = batch.radiated_power_w(2, || false).unwrap();
    let work = batch.responses()[2].power.radiated_w;
    assert!((outgoing - work).abs() < 1e-7 * outgoing.abs() + 1e-20);
    assert_eq!(batch.responses()[3].power.input_w, 0.);
    assert_eq!(batch.responses()[3].power.relative_balance, 0.);
}

#[test]
fn g0_full_live_payload_and_original_boundary_limits_refuse_before_publication() {
    let source = source();
    let prepared = source.prepare(800., boundary_options(), || false).unwrap();
    let force = [C64::new(0.2, -0.4), C64::ONE];
    let baseline = prepared
        .harmonic_response_batch(450., &[&force], options(), || false)
        .unwrap();
    let bytes = (48 * 2 + 512) + 48 * 4 + 160 * 2 + 4096;
    let mut limited = options();
    limited.maximum_response_bytes = bytes - 1;
    assert!(
        prepared
            .harmonic_response_batch(450., &[&force], limited, || false)
            .is_err()
    );
    limited.maximum_response_bytes = bytes;
    assert!(
        prepared
            .harmonic_response_batch(450., &[&force], limited, || false)
            .is_ok()
    );
    let mut limited = options();
    limited.boundary.maximum_matrix_entries = 3;
    assert!(
        prepared
            .harmonic_response_batch(450., &[&force], limited, || false)
            .is_err()
    );
    assert!(
        prepared
            .harmonic_response_batch(801., &[&force], options(), || false)
            .is_err()
    );
    assert!(
        prepared
            .harmonic_response_batch(0., &[&force], options(), || false)
            .is_err()
    );
    assert!(matches!(
        prepared.harmonic_response_batch(450., &[&force], options(), || true),
        Err(BemModesError::Cancelled)
    ));
    let mut calls = 0;
    assert!(
        prepared
            .harmonic_response_batch(450., &[&force], options(), || {
                calls += 1;
                calls > 12
            })
            .is_err()
    );
    let replay = prepared
        .harmonic_response_batch(450., &[&force], options(), || false)
        .unwrap();
    assert_eq!(
        baseline.responses()[0].modal_displacement,
        replay.responses()[0].modal_displacement
    );
    assert_eq!(
        baseline.radiation().impedance(),
        replay.radiation().impedance()
    );
    assert!(
        baseline
            .pressure_at(1, &[], 0.01, near_field::Options::default(), || false)
            .is_err()
    );
}

#[test]
fn g0_source_capacity_and_original_open_topology_are_kept_separate() {
    let original = shell();
    let medium = Medium {
        density: 1.2,
        sound_speed: 343.,
    };
    let mut small = source_options();
    small.maximum_modes = 1;
    assert!(OpenScreenLoadedModes::from_shell(&original, medium, None, small, || false).is_err());
    let mut small = source_options();
    small.maximum_setup_work = 0;
    assert!(OpenScreenLoadedModes::from_shell(&original, medium, None, small, || false).is_err());
    assert!(matches!(
        OpenScreenLoadedModes::from_shell(&original, medium, None, source_options(), || true),
        Err(BemModesError::Cancelled)
    ));
    // The original two triangles have no interior pressure unknown. Preserve
    // this native refusal; acoustic subdivision must be authored explicitly.
    let mut unrefined = source_options();
    unrefined.acoustic_refinements = 0;
    let source =
        OpenScreenLoadedModes::from_shell(&original, medium, None, unrefined, || false).unwrap();
    assert!(source.prepare(800., boundary_options(), || false).is_err());
}
