//! G1/G3: original two-mode mechanics, two wetted faces and one molecular gas.
use super::*;
use crate::bem_loaded_modes::open_screen::tests as fixture;
use fs_bem::reciprocal_radiation::open_screen::ScreenComplexValueLimits;
use fs_exec::CancelGate;
use fs_material::{
    gas::{GasSpec, GasState},
    relaxation::{RelaxingMedium, input::read_relaxation},
};

fn gas_exterior(relaxing: bool, pressure: f64) -> MolecularExterior {
    // Declared caloric fixture, not measured atmospheric validation.
    let data = if relaxing {
        "declared 1 synthetic 1 vibration 1 900 1 power 101.325 293.15 0 200 500"
    } else {
        "declared 1 synthetic 0"
    };
    let law = read_relaxation(&mut data.split_ascii_whitespace()).unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap();
    let medium = RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    MolecularExterior::new(law.evaluate(&medium).unwrap(), [100.0, 800.0]).unwrap()
}

fn source(exterior: &MolecularExterior) -> OpenScreenLoadedModes {
    let viscosity = ModalPlateViscosity::from_psd_modal_matrix(
        vec![0.8, 0.3, 0.3, 0.7],
        2,
        1.0,
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    OpenScreenLoadedModes::from_shell_molecular(
        &fixture::shell(),
        exterior,
        Some(&viscosity),
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
        },
        || false,
    )
    .unwrap()
}

fn options() -> MolecularOpenScreenHarmonicOptions {
    let native = fixture::options();
    MolecularOpenScreenHarmonicOptions {
        boundary: MolecularScreenPortLimits {
            ports: native.boundary,
            operator: ScreenComplexValueLimits {
                maximum_matrix_entries: 65536,
                maximum_workspace_bytes: 4_000_000,
                maximum_work: 1_000_000_000,
                maximum_growth_exponent: 8.0,
            },
            maximum_workspace_bytes: 4_000_000,
        },
        maximum_load_cases: native.maximum_load_cases,
        maximum_response_bytes: native.maximum_response_bytes,
        maximum_work: native.maximum_work,
        physical_residual_tolerance: native.physical_residual_tolerance,
        power_relative_tolerance: native.power_relative_tolerance,
    }
}

fn close(a: C64, b: C64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()) + 1e-22,
        "{a:?} != {b:?}"
    );
}

#[test]
fn no_reservoir_recovers_native_loaded_shell_and_original_two_sided_field() {
    let exterior = gas_exterior(false, 101325.0);
    let source = source(&exterior);
    let prepared = source
        .prepare(800.0, fixture::boundary_options(), || false)
        .unwrap();
    let force = [C64::new(0.5, 0.2), C64::new(-0.1, 0.3)];
    let native = prepared
        .harmonic_response_batch(450.0, &[&force], fixture::options(), || false)
        .unwrap();
    let actual = prepared
        .harmonic_molecular_response_batch(&exterior, 450.0, &[&force], options(), || false)
        .unwrap();
    let (a, b) = (&actual.responses()[0], &native.responses()[0]);
    for i in 0..2 {
        close(a.modal_displacement[i], b.modal_displacement[i], 5e-12);
        close(a.modal_velocity[i], b.modal_velocity[i], 5e-12);
        close(a.fluid_reaction[i], b.fluid_reaction[i], 5e-12);
    }
    assert!(
        (a.power.exterior_w - b.power.radiated_w).abs() < 5e-12 * b.power.radiated_w.abs() + 1e-24
    );
    assert_eq!(a.power.wall_w, 0.0);
    assert_eq!(a.power.relaxation_w, 0.0);
    let points = [[0.04, 0.04, 0.18], [0.04, 0.04, -0.18]];
    let receiver = actual
        .prepare_receivers(
            &points,
            0.001,
            near_field::Options::default(),
            1_000_000,
            || false,
        )
        .unwrap();
    let observed = receiver.evaluate(0, || false).unwrap();
    let expected = native
        .pressure_at(0, &points, 0.001, near_field::Options::default(), || false)
        .unwrap();
    for i in 0..points.len() {
        close(observed.scalar.pressure[i], expected[i], 5e-7);
    }
    assert_eq!(receiver.points_m(), points);
}

#[test]
fn molecular_loaded_cases_keep_full_mechanics_coherent_sound_and_signed_work() {
    let exterior = gas_exterior(true, 101325.0);
    let source = source(&exterior);
    let prepared = source
        .prepare(800.0, fixture::boundary_options(), || false)
        .unwrap();
    let omega = 450.0;
    let f = [C64::new(0.7, -0.2), C64::ZERO];
    let g = [C64::ZERO, C64::new(-0.1, 0.8)];
    let total = [f[0], g[1]];
    let zero = [C64::ZERO; 2];
    let forces: [&[C64]; 4] = [&f, &g, &total, &zero];
    let batch = prepared
        .harmonic_molecular_response_batch(&exterior, omega, &forces, options(), || false)
        .unwrap();
    let z = batch.radiation().impedance();
    assert!(z[1].abs() > 0.0 && z[2].abs() > 0.0);
    let d: [C64; 4] = core::array::from_fn(|at| {
        C64::new(
            source.stiffness[at]
                - if at == 0 || at == 3 {
                    omega * omega
                } else {
                    0.0
                },
            -omega * source.viscosity[at],
        ) + C64::new(0.0, -omega) * z[at]
    });
    let determinant = d[0] * d[3] - d[1] * d[2];
    for (case, force) in forces.iter().enumerate() {
        let r = &batch.responses()[case];
        let independent = [
            (d[3] * force[0] - d[1] * force[1]) / determinant,
            (d[0] * force[1] - d[2] * force[0]) / determinant,
        ];
        for i in 0..2 {
            close(r.modal_displacement[i], independent[i], 5e-12);
            close(
                r.fluid_reaction[i],
                -(z[2 * i] * r.modal_velocity[0] + z[2 * i + 1] * r.modal_velocity[1]),
                5e-12,
            );
        }
        let boundary_work = batch
            .radiation()
            .boundary_work_w(&r.modal_velocity, || false)
            .unwrap();
        assert!((r.power.exterior_w - boundary_work).abs() < 5e-12 * boundary_work.abs() + 1e-24);
        assert!(r.power.relative_balance < 1e-10 && r.mechanical_residual < 1e-10);
    }
    let points = [[0.04, 0.04, 0.08], [0.04, 0.04, -0.08]];
    let receiver = batch
        .prepare_receivers(
            &points,
            0.001,
            near_field::Options::default(),
            1_000_000,
            || false,
        )
        .unwrap();
    let evaluations = receiver.kernel_evaluations();
    let fields: Vec<_> = (0..4)
        .map(|case| receiver.evaluate(case, || false).unwrap())
        .collect();
    for r in 0..points.len() {
        close(
            fields[2].scalar.pressure[r],
            fields[0].scalar.pressure[r] + fields[1].scalar.pressure[r],
            5e-12,
        );
        assert_eq!(fields[3].scalar.pressure[r], C64::ZERO);
        for axis in 0..3 {
            close(
                fields[2].particle_velocity_m_s[r][axis],
                fields[0].particle_velocity_m_s[r][axis] + fields[1].particle_velocity_m_s[r][axis],
                5e-12,
            );
            assert_eq!(fields[3].particle_velocity_m_s[r][axis], C64::ZERO);
        }
    }
    assert_eq!(receiver.kernel_evaluations(), evaluations);
    assert!(receiver.evaluate(4, || false).is_err());
    assert!(matches!(
        receiver.evaluate(0, || true),
        Err(BemModesError::Cancelled)
    ));
    assert_eq!(
        fields[2].scalar.pressure,
        receiver.evaluate(2, || false).unwrap().scalar.pressure
    );
    let nonrelaxing = gas_exterior(false, 101325.0);
    let native = prepared
        .harmonic_molecular_response_batch(&nonrelaxing, omega, &[&total], options(), || false)
        .unwrap();
    assert!(
        batch
            .radiation()
            .impedance()
            .iter()
            .zip(native.radiation().impedance())
            .any(|(a, b)| (*a - *b).abs() > 1e-8 * a.abs())
    );
    let wrong_gas = gas_exterior(true, 90000.0);
    assert!(
        prepared
            .harmonic_molecular_response_batch(&wrong_gas, omega, &[&total], options(), || false)
            .is_err()
    );
    let mut limited = options();
    limited.maximum_response_bytes = 1;
    assert!(
        prepared
            .harmonic_molecular_response_batch(&exterior, omega, &[&total], limited, || false)
            .is_err()
    );
    assert!(matches!(
        prepared.harmonic_molecular_response_batch(&exterior, omega, &[&total], options(), || true),
        Err(BemModesError::Cancelled)
    ));
}
