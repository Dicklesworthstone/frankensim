use super::*;
use crate::plate_viscosity::ModalPlateViscosity;
use fs_bem::observation_transfers::LocalWall;
use fs_exec::CancelGate;
use fs_material::{
    gas::{GasSpec, GasState},
    relaxation::{RelaxingMedium, input::read_relaxation},
};
use fs_phs::{RelaxationBranch, WallPin};
use fs_plate::shell::reduction::radiation::RadiationSurfaceBudget;

pub(super) fn exterior(relaxing: bool) -> MolecularExterior {
    let data = if relaxing {
        "declared 1 synthetic 1 vibration 1 100 1 power 101.325 293.15 0 200 500"
    } else {
        "declared 1 synthetic 0"
    };
    let law = read_relaxation(&mut data.split_ascii_whitespace()).unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let medium = RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    MolecularExterior::new(law.evaluate(&medium).unwrap(), [1.0, 1500.0]).unwrap()
}

fn settings() -> BemModesOptions {
    BemModesOptions {
        formulation: Formulation::PlainCbie,
        boundary: ComplexOptions {
            max_panels: 128,
            max_bytes: 64 * 128 * 128,
            max_work: 512 * 128 * 128,
        },
        panel_refinements: 0,
        reference_angular_frequency_rad_s: 700.0,
        physical_residual_tolerance: 1e-9,
        max_dimension: 256,
        max_dense_bytes: 20_000_000,
        max_operator_evaluations: 16,
        max_projection_work: 20_000_000,
    }
}
pub(super) fn options() -> BemHarmonicOptions {
    BemHarmonicOptions {
        max_load_cases: 8,
        max_response_bytes: 1_000_000,
        max_work: 1_000_000_000,
        power_relative_tolerance: 1e-9,
    }
}
pub(super) fn limits() -> MolecularResponseLimits {
    MolecularResponseLimits {
        maximum_retained_entries: 100_000,
        maximum_panel_products: 1_000_000,
    }
}
pub(super) fn receivers() -> Vec<ReceiverSpec> {
    vec![ReceiverSpec {
        id: "fixed physical receiver".into(),
        position_m: [0.1, 0.05, 0.5],
    }]
}
fn close(a: C64, b: C64, scale: f64) {
    assert!((a - b).abs() < 3e-8 * scale.max(1e-25), "{a:?} != {b:?}");
}

#[test]
fn actual_shell_source_matches_original_full_boundary_memory_solve() {
    // G3: use the existing FE-derived shell's ACTUAL generalized eigenvectors,
    // source Hessian and finite-thickness skin, not authored frequency/shape data.
    let shell = super::super::open_screen::tests::shell();
    let exterior = exterior(false);
    let damping = ModalPlateViscosity::from_psd_modal_matrix(
        vec![0.8, 0.3, 0.3, 0.7],
        2,
        1.0,
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let source = BemLoadedModes::from_shell(
        &shell,
        &[0.001; 4],
        RadiationSurfaceBudget {
            max_panels: 128,
            max_panel_modes: 256,
        },
        Medium::from(exterior.source().gas()),
        Some(&damping),
        settings(),
        || false,
    )
    .unwrap()
    .with_relaxation_branches(
        &[RelaxationBranch {
            projection: vec![0.5, -0.3],
            stiffness: 200.0,
            relaxation_time_s: 0.003,
        }],
        1,
        || false,
    )
    .unwrap();
    let omega = 700.0;
    let force = [C64::new(0.4, 0.1), C64::new(-0.2, 0.15)];
    let expected = source
        .harmonic_response_batch(omega, &[&force], &options(), || false)
        .unwrap()
        .remove(0);
    let prepared = source
        .prepare_molecular_harmonic(
            &exterior,
            omega,
            &receivers(),
            0.01,
            None,
            limits(),
            options(),
            || false,
        )
        .unwrap();
    assert!(std::ptr::eq(prepared.source(), &source));
    assert_eq!(prepared.radiation().source(), exterior.source());
    assert_eq!(prepared.radiation().sources().len(), shell.mode_count());
    let actual = prepared
        .harmonic_response_batch(&[&force], || false)
        .unwrap()
        .remove(0);
    let scale = expected
        .modal_displacement
        .iter()
        .map(|z| z.abs())
        .fold(0.0, f64::max);
    for (&a, &b) in actual
        .modal_displacement
        .iter()
        .zip(&expected.modal_displacement)
    {
        close(a, b, scale);
    }
    for (&a, &b) in actual.modal_velocity.iter().zip(&expected.modal_velocity) {
        close(a, b, scale * omega);
    }
    for (&a, &b) in actual
        .elastic_memory_sqrt_j
        .iter()
        .zip(&expected.elastic_memory_sqrt_j)
    {
        close(a, b, b.abs().max(scale));
    }
    assert!(actual.mechanical_residual < 1e-9 && actual.memory_residual < 1e-9);
    assert!(actual.power.relative_balance < 1e-9);
    assert!(
        (actual.power.exterior_w - expected.power.radiated_w).abs()
            < 3e-8 * expected.power.apparent_power_scale_w
    );
    assert_eq!(actual.power.wall_w, 0.0);
    // Observe a direct original-owner solve of the full loaded trace; this is
    // independent of the cached unit-source receiver superposition.
    let direct = exterior
        .solve_batch(
            source.surface(),
            omega,
            &[&expected.surface_normal_velocity_m_s],
            settings().formulation,
            &settings().boundary,
            || false,
        )
        .unwrap()
        .remove(0);
    let point = receivers()[0].position_m;
    let reference = direct
        .observe(source.surface(), &[point], 0.01, || false)
        .unwrap()
        .remove(0);
    close(
        actual.receivers[0].pressure_pa,
        reference.pressure_pa,
        reference.pressure_pa.abs(),
    );
    let momentum = omega * exterior.source().gas().density;
    for axis in 0..3 {
        close(
            actual.receivers[0].particle_velocity_m_s[axis],
            reference.gradient_pa_per_m[axis] * C64::new(0.0, -1.0 / momentum),
            reference.pressure_pa.abs()
                / (exterior.source().gas().density * exterior.source().gas().sound_speed),
        );
    }
}

pub(super) fn partly_driven(exterior: &MolecularExterior) -> (BemLoadedModes, LocalBoundary) {
    let surface = SpherePanels::icosphere(0.1, 0).unwrap();
    let medium = Medium::from(exterior.source().gas());
    let walls: Vec<_> = surface
        .centroids()
        .iter()
        .map(|p| usize::from(p[2] > 0.0))
        .collect();
    let first = surface
        .centroids()
        .iter()
        .map(|p| if p[2] <= 0.0 { 1.0 } else { 0.0 })
        .collect();
    let second = surface
        .centroids()
        .iter()
        .map(|p| if p[2] <= 0.0 { 0.5 + 3.0 * p[0] } else { 0.0 })
        .collect();
    let source = BemLoadedModes::build(
        surface,
        vec![first, second],
        vec![90_000.0, 5_000.0, 5_000.0, 150_000.0],
        vec![0.8, 0.2, 0.2, 0.5],
        medium,
        settings(),
        &mut || false,
    )
    .unwrap()
    .with_relaxation_branches(
        &[RelaxationBranch {
            projection: vec![0.8, -0.3],
            stiffness: 10_000.0,
            relaxation_time_s: 0.002,
        }],
        1,
        || false,
    )
    .unwrap();
    let boundary = LocalBoundary::try_new(
        medium,
        vec![
            LocalWall::Rigid,
            LocalWall::Mechanical(WallPin {
                surface_density: 0.02,
                stiffness_per_area: 50_000.0,
                resistance: 300.0,
            }),
        ],
        walls,
    )
    .unwrap();
    (source, boundary)
}

#[test]
fn local_wall_changes_actual_motion_while_all_material_work_and_coherent_sources_close() {
    let exterior = exterior(true);
    let (source, boundary) = partly_driven(&exterior);
    let prepared = source
        .prepare_molecular_harmonic(
            &exterior,
            700.0,
            &receivers(),
            0.01,
            Some(boundary),
            limits(),
            options(),
            || false,
        )
        .unwrap();
    let rigid = source
        .prepare_molecular_harmonic(
            &exterior,
            700.0,
            &receivers(),
            0.01,
            None,
            limits(),
            options(),
            || false,
        )
        .unwrap();
    let f = [C64::new(0.4, 0.1), C64::ZERO];
    let g = [C64::ZERO, C64::new(-0.2, 0.15)];
    let h = [f[0], g[1]];
    let zero = [C64::ZERO; 2];
    let batch = prepared
        .harmonic_response_batch(&[&f, &g, &h, &zero], || false)
        .unwrap();
    let no_wall = rigid
        .harmonic_response_batch(&[&h], || false)
        .unwrap()
        .remove(0);
    let total = &batch[2];
    let scale = total
        .modal_displacement
        .iter()
        .map(|z| z.abs())
        .fold(0.0, f64::max);
    assert!((total.modal_displacement[0] - no_wall.modal_displacement[0]).abs() > 1e-7 * scale);
    for mode in 0..2 {
        close(
            total.modal_displacement[mode],
            batch[0].modal_displacement[mode] + batch[1].modal_displacement[mode],
            scale,
        );
    }
    close(
        total.receivers[0].pressure_pa,
        batch[0].receivers[0].pressure_pa + batch[1].receivers[0].pressure_pa,
        total.receivers[0].pressure_pa.abs(),
    );
    let strain = source.relaxation_rows_per_s()[0]
        .iter()
        .zip(&total.modal_displacement)
        .fold(C64::ZERO, |sum, (&c, &x)| sum + x.scale(c));
    let s_tau = C64::new(0.0, -700.0 * source.relaxation_times_s()[0]);
    close(
        total.elastic_memory_sqrt_j[0],
        s_tau / (C64::ONE + s_tau) * strain,
        strain.abs(),
    );
    assert!(
        total.power.kelvin_w > 0.0 && total.power.relaxation_w > 0.0 && total.power.wall_w > 0.0
    );
    assert!(total.power.relative_balance < 1e-9);
    let consumed = total.power.kelvin_w
        + total.power.relaxation_w
        + total.power.exterior_w
        + total.power.wall_w;
    assert!((total.power.input_w - consumed).abs() < 1e-9 * total.power.apparent_power_scale_w);
    let mut reaction = [C64::ZERO; 2];
    prepared
        .radiation()
        .resisting_efforts(&total.modal_velocity, &mut reaction, || false)
        .unwrap();
    for (actual, expected) in total.fluid_reaction.iter().zip(reaction) {
        assert_eq!(*actual, -expected);
    }
    assert!(
        batch[3]
            .modal_displacement
            .iter()
            .chain(&batch[3].elastic_memory_sqrt_j)
            .all(|&v| v == C64::ZERO)
    );
    assert_eq!(batch[3].power.input_w, 0.0);
}

#[test]
fn complete_source_admission_and_late_cancellation_leave_preparation_reusable() {
    let exterior = exterior(true);
    let (source, wall) = partly_driven(&exterior);
    let mut tiny = options();
    tiny.max_work = 0;
    assert!(
        source
            .prepare_molecular_harmonic(
                &exterior,
                700.0,
                &receivers(),
                0.01,
                None,
                limits(),
                tiny,
                || false
            )
            .is_err()
    );
    assert!(
        source
            .prepare_molecular_harmonic(
                &exterior,
                700.0,
                &[],
                0.01,
                None,
                limits(),
                options(),
                || false
            )
            .is_err()
    );
    let source_wall = LocalBoundary::try_new(
        source.medium(),
        vec![LocalWall::Mechanical(WallPin {
            surface_density: 0.02,
            stiffness_per_area: 50_000.0,
            resistance: 300.0,
        })],
        vec![0; source.panel_count()],
    )
    .unwrap();
    assert!(
        source
            .prepare_molecular_harmonic(
                &exterior,
                700.0,
                &receivers(),
                0.01,
                Some(source_wall),
                limits(),
                options(),
                || false
            )
            .is_err()
    );
    let prepared = source
        .prepare_molecular_harmonic(
            &exterior,
            700.0,
            &receivers(),
            0.01,
            Some(wall),
            limits(),
            options(),
            || false,
        )
        .unwrap();
    let force = [C64::ONE, C64::new(0.1, 0.2)];
    let mut total_polls = 0;
    let before = prepared
        .harmonic_response_batch(&[&force], || {
            total_polls += 1;
            false
        })
        .unwrap()
        .remove(0);
    let mut polls = 0;
    assert!(
        prepared
            .harmonic_response_batch(&[&force], || {
                polls += 1;
                polls >= total_polls - 1
            })
            .is_err()
    );
    assert!(prepared.harmonic_response_batch(&[], || false).is_err());
    assert!(
        prepared
            .harmonic_response_batch(&[&force[..1]], || false)
            .is_err()
    );
    let after = prepared
        .harmonic_response_batch(&[&force], || false)
        .unwrap()
        .remove(0);
    assert_eq!(before.modal_displacement, after.modal_displacement);
    assert_eq!(before.receivers, after.receivers);
    assert_eq!(before.power.input_w, after.power.input_w);
}
