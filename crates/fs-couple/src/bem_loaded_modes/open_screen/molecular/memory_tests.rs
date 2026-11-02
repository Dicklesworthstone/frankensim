//! G1/G3: original structural memory in forced, molecular and illuminated shells.
use super::tests as fixture;
use super::*;
use crate::bem_loaded_modes::open_screen::tests as screen_fixture;
use fs_bem::helmholtz::bc::{IncidentField, PlaneWave};
use fs_bem::reciprocal_radiation::open_screen::MolecularScreenIlluminationLimits;
use fs_phs::RelaxationBranch;

fn limits() -> OpenScreenMemoryLimits {
    OpenScreenMemoryLimits {
        maximum_branches: 8,
        maximum_payload_bytes: 100_000,
        maximum_work: 1_000_000,
    }
}

fn arms() -> Vec<RelaxationBranch> {
    vec![
        RelaxationBranch {
            projection: vec![1.0, -0.7],
            stiffness: 50_000.0,
            relaxation_time_s: 0.01,
        },
        RelaxationBranch {
            projection: vec![0.3, 0.8],
            stiffness: 80_000.0,
            relaxation_time_s: 0.003,
        },
    ]
}

fn loaded(exterior: &MolecularExterior, branches: &[RelaxationBranch]) -> OpenScreenLoadedModes {
    fixture::source(exterior)
        .with_relaxation_branches(branches, limits(), || false)
        .unwrap()
}

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * actual.abs().max(expected.abs()) + 1e-22,
        "{actual:?} != {expected:?}"
    );
}

fn close_real(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * actual.abs().max(expected.abs()) + 1e-24,
        "{actual:e} != {expected:e}"
    );
}

/// Independent real-frequency constitutive law, eliminated ONLY in this oracle.
/// Production retains every internal coordinate and its original equation.
fn independent_displacement(
    source: &OpenScreenLoadedModes,
    branches: &[RelaxationBranch],
    omega: f64,
    impedance: &[C64],
    force: &[C64],
) -> [C64; 2] {
    let s = C64::new(0.0, -omega);
    let mut d: [C64; 4] = core::array::from_fn(|at| {
        C64::new(
            source.stiffness_per_s2()[at]
                - if at == 0 || at == 3 {
                    omega * omega
                } else {
                    0.0
                },
            -omega * source.viscosity_per_s()[at],
        ) + s * impedance[at]
    });
    for branch in branches {
        let st = s.scale(branch.relaxation_time_s);
        let factor = st / (C64::ONE + st);
        for i in 0..2 {
            for j in 0..2 {
                d[2 * i + j] = d[2 * i + j]
                    + factor.scale(branch.stiffness * branch.projection[i] * branch.projection[j]);
            }
        }
    }
    let determinant = d[0] * d[3] - d[1] * d[2];
    [
        (d[3] * force[0] - d[1] * force[1]) / determinant,
        (d[0] * force[1] - d[2] * force[0]) / determinant,
    ]
}

fn independent_memory(branches: &[RelaxationBranch], omega: f64, q: &[C64]) -> Vec<C64> {
    let s = C64::new(0.0, -omega);
    branches
        .iter()
        .map(|branch| {
            let strain = branch
                .projection
                .iter()
                .zip(q)
                .fold(C64::ZERO, |sum, (&b, &x)| sum + x.scale(b));
            let st = s.scale(branch.relaxation_time_s);
            (st / (C64::ONE + st)) * strain.scale(branch.stiffness.sqrt())
        })
        .collect()
}

fn heat(branches: &[RelaxationBranch], memory: &[C64]) -> f64 {
    branches
        .iter()
        .zip(memory)
        .map(|(branch, e)| 0.5 * e.norm_sq() / branch.relaxation_time_s)
        .sum()
}

fn structural_energy(source: &OpenScreenLoadedModes, q: &[C64], v: &[C64], e: &[C64]) -> f64 {
    let mut energy = 0.25 * v.iter().chain(e).map(|z| z.norm_sq()).sum::<f64>();
    for i in 0..2 {
        for j in 0..2 {
            energy += 0.25 * source.stiffness_per_s2()[2 * i + j] * (q[i].conj() * q[j]).re;
        }
    }
    energy
}

#[test]
fn g1_real_and_no_reservoir_shells_recover_full_maxwell_law_storage_and_heat() {
    let exterior = fixture::gas_exterior(false, 101325.0);
    let branches = arms();
    let source = loaded(&exterior, &branches);
    assert_eq!(source.memory_count(), branches.len());
    assert_eq!(source.relaxation_times_s(), [0.01, 0.003]);
    for (row, branch) in source.relaxation_rows_per_s().iter().zip(&branches) {
        for (&c, &b) in row.iter().zip(&branch.projection) {
            close_real(c, branch.stiffness.sqrt() * b, 1e-14);
        }
    }
    let prepared = source
        .prepare(800.0, screen_fixture::boundary_options(), || false)
        .unwrap();
    let omega = 450.0;
    let force = [C64::new(0.7, -0.2), C64::new(-0.1, 0.8)];
    let native = prepared
        .harmonic_response_batch(omega, &[&force], screen_fixture::options(), || false)
        .unwrap();
    let r = &native.responses()[0];
    let expected = independent_displacement(
        &source,
        &branches,
        omega,
        native.radiation().impedance(),
        &force,
    );
    let expected_memory = independent_memory(&branches, omega, &expected);
    for i in 0..2 {
        close(r.modal_displacement[i], expected[i], 2e-11);
        close(r.modal_velocity[i], C64::new(0.0, -omega) * expected[i], 2e-11);
        close(r.elastic_memory_sqrt_j[i], expected_memory[i], 2e-11);
    }
    close_real(r.power.relaxation_w, heat(&branches, &expected_memory), 3e-11);
    close_real(
        r.mean_structural_energy_j,
        structural_energy(&source, &expected, &r.modal_velocity, &expected_memory),
        3e-11,
    );
    assert!(r.power.relaxation_w > 0.0);
    assert!(r.memory_residual < 1e-10 && r.mechanical_residual < 1e-10);
    assert!(r.power.relative_balance < 1e-10);
    close_real(
        r.power.input_w,
        r.power.kelvin_w + r.power.relaxation_w + r.power.radiated_w,
        1e-10,
    );
    let molecular = prepared
        .harmonic_molecular_response_batch(&exterior, omega, &[&force], fixture::options(), || false)
        .unwrap();
    let converted = &molecular.responses()[0];
    for i in 0..2 {
        close(converted.modal_displacement[i], r.modal_displacement[i], 2e-11);
        close(converted.elastic_memory_sqrt_j[i], r.elastic_memory_sqrt_j[i], 2e-11);
        close(converted.fluid_reaction[i], r.fluid_reaction[i], 2e-11);
    }
    close_real(converted.power.relaxation_w, r.power.relaxation_w, 2e-11);
    close_real(converted.power.exterior_w, r.power.radiated_w, 2e-11);
    close_real(converted.mean_structural_energy_j, r.mean_structural_energy_j, 2e-11);
    assert!(converted.memory_residual < 1e-10);
    assert_eq!(converted.power.wall_w, 0.0);
}

#[test]
fn g3_molecular_gas_and_equal_time_material_arms_keep_coherent_motion_and_receivers() {
    let exterior = fixture::gas_exterior(true, 101325.0);
    let branches = vec![
        RelaxationBranch {
            projection: vec![1.0, -0.7],
            stiffness: 50_000.0,
            relaxation_time_s: 0.002,
        },
        RelaxationBranch {
            projection: vec![2.0, -1.4],
            stiffness: 50_000.0,
            relaxation_time_s: 0.002,
        },
    ];
    let source = loaded(&exterior, &branches);
    // Proportional, equal-time arms still own separate stored-energy coordinates.
    assert_eq!(source.memory_count(), 2);
    let prepared = source
        .prepare(800.0, screen_fixture::boundary_options(), || false)
        .unwrap();
    let f = [C64::new(0.7, -0.2), C64::ZERO];
    let g = [C64::ZERO, C64::new(-0.1, 0.8)];
    let total = [f[0], g[1]];
    let zero = [C64::ZERO; 2];
    let batch = prepared
        .harmonic_molecular_response_batch(
            &exterior,
            450.0,
            &[&f, &g, &total, &zero],
            fixture::options(),
            || false,
        )
        .unwrap();
    let r = &batch.responses()[2];
    let expected = independent_displacement(
        &source,
        &branches,
        450.0,
        batch.radiation().impedance(),
        &total,
    );
    let expected_memory = independent_memory(&branches, 450.0, &expected);
    for i in 0..2 {
        close(r.modal_displacement[i], expected[i], 2e-11);
        close(r.elastic_memory_sqrt_j[i], expected_memory[i], 2e-11);
        close(
            r.elastic_memory_sqrt_j[i],
            batch.responses()[0].elastic_memory_sqrt_j[i]
                + batch.responses()[1].elastic_memory_sqrt_j[i],
            2e-11,
        );
    }
    close(r.elastic_memory_sqrt_j[1], r.elastic_memory_sqrt_j[0].scale(2.0), 2e-11);
    close_real(r.power.relaxation_w, heat(&branches, &expected_memory), 3e-11);
    close_real(
        r.power.exterior_w,
        batch.radiation().boundary_work_w(&r.modal_velocity, || false).unwrap(),
        3e-11,
    );
    assert_eq!(batch.responses()[3].elastic_memory_sqrt_j, [C64::ZERO; 2]);
    assert_eq!(batch.responses()[3].power.relaxation_w, 0.0);
    assert_eq!(batch.responses()[3].power.relative_balance, 0.0);
    let points = [[0.04, 0.04, 0.12], [0.04, 0.04, -0.12]];
    let receiver = batch
        .prepare_receivers(
            &points,
            0.001,
            near_field::Options::default(),
            1_000_000,
            || false,
        )
        .unwrap();
    let count = receiver.kernel_evaluations();
    let fields: Vec<_> = (0..4)
        .map(|i| receiver.evaluate(i, || false).unwrap())
        .collect();
    for point in 0..2 {
        close(
            fields[2].scalar.pressure[point],
            fields[0].scalar.pressure[point] + fields[1].scalar.pressure[point],
            3e-11,
        );
        for axis in 0..3 {
            close(
                fields[2].particle_velocity_m_s[point][axis],
                fields[0].particle_velocity_m_s[point][axis]
                    + fields[1].particle_velocity_m_s[point][axis],
                3e-11,
            );
            assert_eq!(fields[3].particle_velocity_m_s[point][axis], C64::ZERO);
        }
        assert_eq!(fields[3].scalar.pressure[point], C64::ZERO);
    }
    assert_eq!(receiver.kernel_evaluations(), count);
    assert!(matches!(receiver.evaluate(2, || true), Err(BemModesError::Cancelled)));
    assert_eq!(
        fields[2].scalar.pressure,
        receiver.evaluate(2, || false).unwrap().scalar.pressure,
    );
}

fn incident_options() -> MolecularOpenScreenIncidentOptions {
    let original = fixture::options();
    MolecularOpenScreenIncidentOptions {
        boundary: MolecularScreenIlluminationLimits {
            boundary: original.boundary,
            maximum_incident_cases: 8,
            maximum_incident_jump_entries: 2048,
            maximum_source_point_products: 1_000_000,
            minimum_source_clearance_m: 0.001,
            maximum_live_bytes: 8_000_000,
        },
        maximum_load_cases: 8,
        maximum_response_bytes: 200_000,
        maximum_work: 2_000_000,
        physical_residual_tolerance: 1e-10,
        power_relative_tolerance: 1e-10,
    }
}

#[test]
fn g1_incident_sound_pays_structural_relaxation_heat_and_drives_the_observed_motion() {
    let exterior = fixture::gas_exterior(true, 101325.0);
    let branches = arms();
    let source = loaded(&exterior, &branches);
    let prepared = source
        .prepare(800.0, screen_fixture::boundary_options(), || false)
        .unwrap();
    let wave = IncidentField::PlaneWave(PlaneWave {
        direction: [0.0, 0.0, -1.0],
        amplitude_pa: C64::new(2.0, 0.4),
    });
    let zero = [C64::ZERO; 2];
    let external = [C64::new(0.02, 0.01), C64::new(-0.01, 0.03)];
    let forces: [&[C64]; 2] = [&zero, &external];
    let batch = prepared
        .harmonic_molecular_incident_response_batch(
            &exterior,
            450.0,
            &[wave, wave],
            &forces,
            incident_options(),
            || false,
        )
        .unwrap();
    for (case, force) in forces.iter().enumerate() {
        let r = &batch.responses()[case];
        let effective: Vec<_> = force
            .iter()
            .zip(&r.blocked_acoustic_force)
            .map(|(&a, &b)| a + b)
            .collect();
        let expected = independent_displacement(
            &source,
            &branches,
            450.0,
            batch.radiation().motion().impedance(),
            &effective,
        );
        let e = independent_memory(&branches, 450.0, &expected);
        for i in 0..2 {
            close(r.motion.modal_displacement[i], expected[i], 3e-11);
            close(r.motion.elastic_memory_sqrt_j[i], e[i], 3e-11);
        }
        close_real(r.power.structural_relaxation_w, heat(&branches, &e), 3e-11);
        close_real(r.motion.power.relaxation_w, r.power.structural_relaxation_w, 3e-11);
        assert!(r.motion.memory_residual < 1e-10);
        assert!(r.power.relative_balance < 1e-10);
        let residual = r.power.external_mechanical_input_w - r.power.structural_viscosity_w
            - r.power.structural_relaxation_w - r.power.total_exterior_w;
        assert!(residual.abs() < 1e-10 * r.power.apparent_power_scale_w);
    }
    let only = &batch.responses()[0].power;
    assert_eq!(only.external_mechanical_input_w, 0.0);
    assert!(only.structural_relaxation_w > 0.0 && only.total_exterior_w < 0.0);
    close_real(
        -only.total_exterior_w,
        only.structural_viscosity_w + only.structural_relaxation_w,
        1e-10,
    );
    let points = [[0.04, 0.04, 0.12], [0.04, 0.04, -0.12]];
    let rows = batch
        .radiation()
        .prepare_receivers(
            &points,
            0.001,
            near_field::Options::default(),
            1_000_000,
            || false,
        )
        .unwrap();
    let receiver = batch
        .prepare_receivers(
            &points,
            0.001,
            near_field::Options::default(),
            1_000_000,
            || false,
        )
        .unwrap();
    let independent = independent_displacement(
        &source,
        &branches,
        450.0,
        batch.radiation().motion().impedance(),
        &batch.responses()[0].blocked_acoustic_force,
    );
    let velocity = independent.map(|q| C64::new(0.0, -450.0) * q);
    let expected = rows.evaluate(0, &velocity, || false).unwrap();
    let actual = receiver.evaluate(0, || false).unwrap();
    for point in 0..2 {
        close(actual.scalar.pressure[point], expected.scalar.pressure[point], 3e-11);
        for axis in 0..3 {
            close(
                actual.particle_velocity_m_s[point][axis],
                expected.particle_velocity_m_s[point][axis],
                3e-11,
            );
        }
    }
}

#[test]
fn g0_memory_admission_empty_attachment_and_cancelled_batches_preserve_the_source() {
    let exterior = fixture::gas_exterior(false, 101325.0);
    let source = fixture::source(&exterior);
    let force = [C64::new(0.7, -0.2), C64::new(-0.1, 0.8)];
    let baseline = {
        let prepared = source
            .prepare(800.0, screen_fixture::boundary_options(), || false)
            .unwrap();
        let batch = prepared
            .harmonic_response_batch(450.0, &[&force], screen_fixture::options(), || false)
            .unwrap();
        batch.responses()[0].clone()
    };
    let source = source
        .with_relaxation_branches(&[], limits(), || false)
        .unwrap();
    let prepared = source
        .prepare(800.0, screen_fixture::boundary_options(), || false)
        .unwrap();
    let after = prepared
        .harmonic_response_batch(450.0, &[&force], screen_fixture::options(), || false)
        .unwrap();
    assert_eq!(baseline.modal_displacement, after.responses()[0].modal_displacement);
    assert_eq!(baseline.modal_velocity, after.responses()[0].modal_velocity);
    assert_eq!(
        baseline.mean_structural_energy_j.to_bits(),
        after.responses()[0].mean_structural_energy_j.to_bits(),
    );
    assert_eq!(after.responses()[0].elastic_memory_sqrt_j.len(), 0);
    let branches = arms();
    let mut tiny = limits();
    tiny.maximum_branches = 1;
    assert!(fixture::source(&exterior).with_relaxation_branches(&branches, tiny, || false).is_err());
    tiny = limits();
    tiny.maximum_payload_bytes = 1;
    assert!(fixture::source(&exterior).with_relaxation_branches(&branches, tiny, || false).is_err());
    tiny = limits();
    tiny.maximum_work = 1;
    assert!(fixture::source(&exterior).with_relaxation_branches(&branches, tiny, || false).is_err());
    assert!(matches!(
        fixture::source(&exterior).with_relaxation_branches(&branches, limits(), || true),
        Err(BemModesError::Cancelled)
    ));
    let source = loaded(&exterior, &branches);
    let prepared = source
        .prepare(800.0, screen_fixture::boundary_options(), || false)
        .unwrap();
    let mut budget = screen_fixture::options();
    // The former two-coordinate allowance cannot admit four physical states.
    budget.maximum_response_bytes = (16 * 2 + 32 * 2 + 640) + 96 * 4 + 320 * 2 + 4096;
    assert!(prepared.harmonic_response_batch(450.0, &[&force], budget, || false).is_err());
    let before = prepared
        .harmonic_response_batch(450.0, &[&force], screen_fixture::options(), || false)
        .unwrap();
    let mut calls = 0;
    assert!(
        prepared
            .harmonic_response_batch(450.0, &[&force], screen_fixture::options(), || {
                calls += 1;
                calls > 12
            })
            .is_err()
    );
    let replay = prepared
        .harmonic_response_batch(450.0, &[&force], screen_fixture::options(), || false)
        .unwrap();
    assert_eq!(before.responses()[0].modal_displacement, replay.responses()[0].modal_displacement);
    assert_eq!(before.responses()[0].elastic_memory_sqrt_j, replay.responses()[0].elastic_memory_sqrt_j);
}
