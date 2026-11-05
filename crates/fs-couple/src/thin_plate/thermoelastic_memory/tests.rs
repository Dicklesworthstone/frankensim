//! G1/G3: tensor curvature, coupled dispersion, thermal histories and live ports.
use super::*;
use crate::thin_plate::{GasState, PlateBank};
use fs_math::c64::C64;
use fs_nlmodal::SosModalStorage;

fn law() -> SlabThermoelasticBending {
    SlabThermoelasticBending::new_coupled(
        [120.0, 0.0, 0.0, 0.0, 120.0, 0.0, 0.0, 0.0, 60.0],
        [0.1, 0.1, 0.0],
        1.0,
        20.0,
        4.0,
        1.0,
    )
    .unwrap()
}
fn mesh() -> PlateMesh {
    PlateMesh {
        nodes: vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)],
        tris: vec![[0, 1, 2]],
    }
}
fn patch(mesh: &PlateMesh, k: [f64; 3]) -> Vec<f64> {
    mesh.nodes
        .iter()
        .flat_map(|&(x, y)| {
            [
                0.5 * (k[0] * x * x + k[1] * y * y) + k[2] * x * y,
                k[0] * x + k[2] * y,
                k[1] * y + k[2] * x,
            ]
        })
        .collect()
}
fn spectrum(law: SlabThermoelasticBending) -> SlabPronySpectrum {
    law.lower_to_prony((0.0, 100.0), 0.01, 8).unwrap()
}
fn arms(
    fields: &[Vec<f64>],
    scales: &[f64],
    law: SlabThermoelasticBending,
    spectrum: &SlabPronySpectrum,
) -> Vec<fs_phs::RelaxationBranch> {
    project_thermal_relaxation(
        &mesh(),
        &(0..9).map(Some).collect::<Vec<_>>(),
        &fields.iter().map(Vec::as_slice).collect::<Vec<_>>(),
        scales,
        law,
        spectrum,
        64,
        1_000_000,
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn close(a: f64, b: f64, relative: f64, floor: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(floor),
        "{a:e} != {b:e}"
    );
}

#[test]
fn original_curvature_retains_full_cross_terms_and_null_direction() {
    let mesh = mesh();
    let law = law();
    let spectrum = spectrum(law);
    let fields = vec![
        patch(&mesh, [1.0, 0.0, 0.0]),
        patch(&mesh, [0.0, 1.0, 0.0]),
        patch(&mesh, [1.0, -1.0, 0.0]),
    ];
    let branches = arms(&fields, &[0.5, 2.0, 1.0], law, &spectrum);
    let mut actual = [0.0; 9];
    for arm in &branches {
        for i in 0..3 {
            for j in 0..3 {
                actual[3 * i + j] += arm.stiffness * arm.projection[2 * i] * arm.projection[2 * j];
            }
        }
    }
    // Independent constant-curvature slab integral: area=1/2, beta=(12,12,0),
    // C_epsilon=20-1*(.1*12+.1*12)=17.6, h=1, T=1.
    let relaxation = 0.5 * 12.0 * 12.0 / (12.0 * 17.6);
    let weight: f64 = spectrum.terms().iter().map(|t| t.weight()).sum();
    let drive = [0.5, 2.0, 0.0];
    for i in 0..3 {
        for j in 0..3 {
            let expected = relaxation * weight * drive[i] * drive[j];
            close(actual[3 * i + j], expected, 2e-13, 1e-16);
        }
    }
    assert!(actual[1] > 0.1, "cross-mode thermal coupling must survive");
    assert!(
        actual[8].abs() < 1e-28,
        "trace-free curvature must remain lossless"
    );
    assert!(branches.len() <= 3 * spectrum.terms().len());
}

#[test]
fn source_projection_and_qr_keep_weak_terms_without_a_rank_cutoff() {
    let small = 2.0_f64.powi(-55);
    let mut curvature = [[0.0; 9]; 3];
    curvature[0][0] = 1.0;
    curvature[0][1] = small;
    curvature[0][2] = -1.0;
    let dofs = core::array::from_fn(Some);
    assert_eq!(
        thermal_drive([1.0, 0.0, 0.0], &curvature, &dofs, &[1.0; 9], 1.0),
        small
    );
    curvature = [[0.0; 9]; 3];
    for (i, row) in curvature.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    assert_eq!(
        thermal_drive([1.0, small, -1.0], &curvature, &dofs, &[1.0; 9], 1.0),
        small
    );
    for weak in [1e-12, 1e-120, 1e-150] {
        let original = [1.0, weak, 0.0, 1.0, -weak, 0.0, 0.0, 0.0, 0.0];
        let mut rows = original;
        let factors = compress_rows(&mut rows, 3, 3, &CancelGate::new_clock_free()).unwrap();
        assert_eq!(factors.len(), 2);
        let strong: f64 = factors.iter().map(|r| r[0] * r[0]).sum();
        let small: f64 = factors.iter().map(|r| r[1] * r[1]).sum();
        close(strong, 2.0, 2e-14, 0.0);
        close(small / (weak * weak), 2.0, 2e-14, 0.0);
        assert!(factors.iter().all(|r| r[2] == 0.0));
    }
}

fn manufactured(omega: [f64; 2]) -> (VkBody, Vec<fs_phs::RelaxationBranch>, SlabPronySpectrum) {
    let mesh = mesh();
    let law = law();
    let spectrum = spectrum(law);
    let fields = vec![patch(&mesh, [1.0, 0.0, 0.0]), patch(&mesh, [0.0, 1.0, 0.0])];
    let branches = arms(&fields, &[1.0, 1.0], law, &spectrum);
    let body = VkBody::from_modal_storage_with_damping(
        SosModalStorage {
            omegas: omega.to_vec(),
            channels: vec![],
        },
        &[0.0; 4],
        &[1.0, 0.0],
        vec![0.4, -0.2],
    )
    .unwrap()
    .with_relaxation_branches(branches.clone(), 64)
    .unwrap();
    (body, branches, spectrum)
}

fn vector_field(body: &VkBody, x: &[f64], force: f64) -> Vec<f64> {
    let (j, r, g) = body.sys.structure();
    let n = x.len();
    let e = body.sys.effort(x);
    (0..n)
        .map(|i| {
            (0..n)
                .map(|k| (j[n * i + k] - r[n * i + k]) * e[k])
                .sum::<f64>()
                + g[i] * force
        })
        .collect()
}

#[test]
fn close_modes_follow_the_full_complex_thermal_operator() {
    let omega0 = [1.0, 1.02];
    let (body, branches, spectrum) = manufactured(omega0);
    let delta = 0.5 * 12.0 * 12.0 / (12.0 * 17.6);
    for w in [0.3, 0.95, 1.05, 2.0, 10.0] {
        // Analytic TWO-MODE harmonic resolvent. Off-diagonal terms alter both
        // modal responses; checking each diagonal loss ratio would miss this.
        let mut increment = C64::ZERO;
        for term in spectrum.terms() {
            let x = w * term.relaxation_time_s();
            increment = increment + C64::new(x * x, x).scale(delta * term.weight() / (1.0 + x * x));
        }
        let d0 = C64::from_re(omega0[0] * omega0[0] - w * w) + increment;
        let d1 = C64::from_re(omega0[1] * omega0[1] - w * w) + increment;
        let inverse = (d0 * d1 - increment * increment).recip();
        let q = [d1 * inverse, -increment * inverse];
        assert!(
            q[1].abs() > 1e-7,
            "the unforced coordinate is thermally coupled"
        );
        let mut state = vec![C64::ZERO; body.state().len()];
        for i in 0..2 {
            state[2 * i] = q[i];
            state[2 * i + 1] = C64::new(0.0, w) * q[i];
        }
        for (i, arm) in branches.iter().enumerate() {
            let strain = q[0].scale(arm.projection[0]) + q[1].scale(arm.projection[2]);
            state[4 + i] = strain.scale(det::sqrt(arm.stiffness))
                * C64::new(1.0, w * arm.relaxation_time_s).recip();
        }
        let real: Vec<_> = state.iter().map(|z| z.re).collect();
        let imag: Vec<_> = state.iter().map(|z| z.im).collect();
        let real_rhs = vector_field(&body, &real, 1.0);
        let imag_rhs = vector_field(&body, &imag, 0.0);
        let norm = state.iter().map(|z| z.abs()).fold(1.0_f64, f64::max) * w;
        for i in 0..state.len() {
            assert!((real_rhs[i] + w * imag[i]).abs() < 3e-11 * norm);
            assert!((imag_rhs[i] - w * real[i]).abs() < 3e-11 * norm);
        }
    }
}

#[test]
fn full_memory_dissipates_the_bright_motion_and_preserves_the_dark_motion() {
    let amplitude = 1e-3;
    let dt = 0.02;
    let steps = 50;
    for sign in [-1.0, 1.0] {
        let (mut body, _, _) = manufactured([1.0, 1.0]);
        body.x[0] = amplitude;
        body.x[2] = sign * amplitude;
        let energy0 = body.stored_energy_j();
        let mut body = body.into_prepared(dt).unwrap();
        let mut heat = 0.0;
        for _ in 0..steps {
            let step = body
                .drive_and_radiate_controlled(0.0, 1.2, 1.0, || false)
                .unwrap();
            heat += step.work.dissipated;
        }
        let energy1 = body.stored_energy_j();
        close(energy1 + heat, energy0, 2e-10, 1e-14);
        if sign < 0.0 {
            let t = dt * steps as f64;
            close(body.state()[0], amplitude * t.cos(), 2e-10, 1e-12);
            close(body.state()[1], -amplitude * t.sin(), 2e-10, 1e-12);
            assert!(
                heat < 1e-22,
                "a thermal null direction acquired heat: {heat:e}"
            );
        } else {
            assert!(
                heat > 1e-8,
                "the common curvature must produce physical heat"
            );
            assert!(energy1 < energy0);
        }
    }
}

fn plate() -> ThinPlate {
    ThinPlate {
        length_m: 0.20,
        width_m: 0.15,
        thickness_m: 0.002,
        density_kg_m3: 5000.0,
        e1_pa: 70e9,
        e2_pa: 40e9,
        nu12: 0.3,
        g12_pa: 20e9,
        material_angle_rad: 0.31,
        damping_ratio: 0.0,
        thermoelastic: None,
        kelvin_voigt_bending: None,
        n_modes: 2,
        geometric_nonlinearity: true,
        pretension_n_m: 1000.0,
        clamped: false,
    }
}
fn thermal() -> AnisotropicPlateThermal {
    AnisotropicPlateThermal {
        temperature_k: 300.0,
        engineering_expansion_per_k: [2e-5, 4e-5, 1e-5],
        specific_heat_j_kg_k: 600.0,
        thickness_conductivity_w_m_k: 100.0,
    }
}
fn options() -> ThermoelasticMemoryOptions {
    ThermoelasticMemoryOptions {
        band_rad_s: (0.0, 1e6),
        normalized_tolerance: 1e-3,
        max_thermal_terms: 8,
        max_memory_states: 64,
        max_projection_terms: 10_000_000,
    }
}

#[test]
fn original_nonlinear_plate_and_bank_image_retain_thermal_state_and_ports() {
    let gate = CancelGate::new_clock_free();
    let (body, spectrum) =
        VkBody::from_plate_with_thermoelastic_memory(plate(), thermal(), options(), &gate).unwrap();
    let extent = body.state().len();
    assert!(extent > 2 * body.areas.len());
    assert!(spectrum.normalized_error_bound() <= options().normalized_tolerance);
    let gas = GasState::try_new(
        &fs_material::gas::GasSpec::dry_air_ussa1976(),
        293.15,
        101_325.0,
    )
    .unwrap();
    let mut bank = PlateBank::default();
    bank.vk = Some(body);
    for _ in 0..8 {
        let step = bank
            .advance_controlled(1.0, 1e-6, gas.density, 1.0, || false)
            .unwrap();
        assert!(step.work.dissipated >= 0.0);
        let scale = step.work.delta_h.abs() + step.work.dissipated + step.work.supplied.abs();
        assert!(
            (step.work.delta_h + step.work.dissipated - step.work.supplied).abs()
                <= 1e-9 * scale.max(1e-24)
        );
    }
    let image = bank.physical_image(&gas, 48_000, &gate).unwrap().unwrap();
    assert_eq!(image.state_dim(), extent);
    assert!(
        image.state()[4..].iter().any(|v| *v != 0.0),
        "thermal histories must evolve"
    );
    close(
        image.system().hamiltonian(image.state()),
        bank.stored_energy_j(),
        1e-13,
        1e-24,
    );
    close(
        image.observation().volume_velocity(image.state()).unwrap(),
        bank.volume_velocity(),
        1e-13,
        1e-24,
    );
}

#[test]
fn missing_resources_duplicate_loss_and_zero_thermal_channels_are_explicit() {
    let mesh = mesh();
    let original = law();
    let source_spectrum = spectrum(original);
    let fields = [patch(&mesh, [1.0, 0.0, 0.0])];
    let project = |law, cap, work, gate: &CancelGate| {
        project_thermal_relaxation(
            &mesh,
            &(0..9).map(Some).collect::<Vec<_>>(),
            &[fields[0].as_slice()],
            &[1.0],
            law,
            &source_spectrum,
            cap,
            work,
            gate,
        )
    };
    let gate = CancelGate::new_clock_free();
    assert!(project(original, 0, 1_000_000, &gate).is_err());
    assert!(project(original, 64, 0, &gate).is_err());
    let stopped = CancelGate::new_clock_free();
    stopped.request();
    assert!(matches!(
        project(original, 64, 1_000_000, &stopped),
        Err(AcousticRealizeError::Cancelled)
    ));
    let mut duplicate = plate();
    duplicate.damping_ratio = 0.01;
    assert!(
        VkBody::from_plate_with_thermoelastic_memory(duplicate, thermal(), options(), &gate,)
            .is_err()
    );
    let zero = SlabThermoelasticBending::new_coupled(
        [120.0, 0.0, 0.0, 0.0, 120.0, 0.0, 0.0, 0.0, 60.0],
        [0.0; 3],
        1.0,
        20.0,
        4.0,
        1.0,
    )
    .unwrap();
    let zero_spectrum = spectrum(zero);
    let empty = project_thermal_relaxation(
        &mesh,
        &(0..9).map(Some).collect::<Vec<_>>(),
        &[fields[0].as_slice()],
        &[1.0],
        zero,
        &zero_spectrum,
        0,
        1_000_000,
        &gate,
    )
    .unwrap();
    assert!(empty.is_empty());
    assert!(
        admit_memory_band(
            &[1.0],
            &[fs_phs::RelaxationBranch {
                projection: vec![10.0, 0.0],
                stiffness: 1.0,
                relaxation_time_s: 1.0,
            }],
            (0.0, 2.0)
        )
        .is_err()
    );
}
