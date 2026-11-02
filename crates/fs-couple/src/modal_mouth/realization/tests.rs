use super::*;
use crate::modal_mouth::ModalMouthBemBudget;
use fs_bem::panel3d::SpherePanels;
use fs_bem::radiation_bake::lathe_profile;
use fs_bem::reciprocal_radiation;
use fs_duct::modal::{ModalMouthImpedance, mm_input_impedance_matrix};
use fs_material::gas::GasSpec;
use fs_math::c64::C64;
use fs_vfit::causality::CausalityOptions;
use fs_vfit::passivity::CertificateClass;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}

fn geometry() -> (SpherePanels, CircularMouthPatch) {
    // Authored coarse closed cylinder, with two cap annuli. This exercises
    // the actual geometric operator; it is not a continuum-accuracy fixture.
    let radius = 0.02;
    let bell = lathe_profile(&[(0.0, radius), (0.01, radius), (0.02, radius)], 8).unwrap();
    let patch = CircularMouthPatch {
        centre_m: [0.02, 0.0, 0.0],
        outward_normal: [1.0, 0.0, 0.0],
        radius_m: radius,
        panels: (bell.mouth_start..bell.triangles.len()).collect(),
    };
    (SpherePanels::from_triangles(bell.triangles).unwrap(), patch)
}

fn budget() -> ModalMouthBemBudget {
    ModalMouthBemBudget {
        maximum_panels: 200,
        maximum_trace_entries: 400,
        maximum_matrix_entries: 2000,
        maximum_projection_terms: 1_000_000,
        maximum_basis_gram_defect: 0.45,
    }
}

fn route() -> ApertureRoute {
    ApertureRoute::Galerkin {
        budget: reciprocal_radiation::Budget {
            maximum_unknowns: 400,
            maximum_cache_bytes: 16 << 20,
            maximum_work: 4_000_000_000,
            angular_refinement: 0,
        },
    }
}

fn options() -> ModalMouthRealizationOptions {
    ModalMouthRealizationOptions {
        mechanics_rate_hz: 6000,
        source_band_hz: [1.0, 3000.0],
        operator_band_hz: [30.0, 3000.0],
        fit: BaffledMimoBakeOptions {
            frequency_samples: 32,
            orders: [2, 4, 6],
            fit_iterations: 12,
            maximum_relative_error: 0.03,
            max_fit_memory_bytes: 64 << 20,
            causality: CausalityOptions {
                subtractions: 2,
                tail_bound_ratio: 100.0,
                relative_data_tolerance: 1e-8,
            },
        },
        load: MimoRadiationLimits {
            max_states: 32,
            max_memory_bytes: 1 << 20,
            sampled_passivity_points: 257,
        },
        maximum_reciprocity_defect: 1e-10,
    }
}

fn relative_matrix(a: &[C64], b: &[C64]) -> f64 {
    let scale = b.iter().map(|z| z.norm_sq()).sum::<f64>().sqrt();
    a.iter()
        .zip(b)
        .map(|(x, y)| (*x - *y).norm_sq())
        .sum::<f64>()
        .sqrt()
        / scale
}

#[test]
fn geometric_signed_modal_memory_preserves_physical_samples_and_loaded_response() {
    let (surface, patch) = geometry();
    let gate = CancelGate::new();
    let source = BemModalMouth::new(&surface, patch, air(), ModeCount(2), budget(), &gate).unwrap();
    let options = options();
    let realized = source.realize(&route(), options, &gate).unwrap();
    let record = realized.record();
    let report = record.report();
    assert_eq!(record.modes(), ModeCount(2));
    assert_eq!(record.patch().centre_m, source.patch().centre_m);
    assert_eq!(record.exterior(), source.exterior());
    assert_eq!(record.options().source_band_hz, options.source_band_hz);
    assert_eq!(record.options().operator_band_hz, options.operator_band_hz);
    assert_eq!(record.training().loads().len(), 32);
    assert_eq!(record.held_out().loads().len(), 31);
    assert_eq!(record.causality().loads().len(), 32);
    assert_eq!(report.rungs.len(), 3);
    assert_eq!(
        report.rungs.iter().map(|r| r.order).collect::<Vec<_>>(),
        vec![2, 4, 6]
    );
    assert!(report.continuous.passive && report.sampled.passive);
    assert_eq!(report.sampled.class, CertificateClass::GridOnly);
    assert_eq!(realized.load().time_step_s(), 1.0 / 6000.0);
    assert_eq!(realized.load().ports(), 2);
    assert!(realized.load().source_model().order() > 0);
    assert!(
        realized
            .load()
            .pole_states()
            .iter()
            .all(|z| *z == C64::ZERO)
    );
    assert!(source.velocity_trace(1).unwrap().iter().any(|v| *v < 0.0));
    assert!(source.velocity_trace(1).unwrap().iter().any(|v| *v > 0.0));

    let duct = fs_duct::Duct {
        segments: vec![fs_duct::Segment::Cylinder {
            radius: source.patch().radius_m,
            length: 0.06,
        }],
    };
    let mut maximum_loaded_error = 0.0_f64;
    for actual in record.held_out().loads() {
        assert!(
            !record
                .training()
                .loads()
                .iter()
                .any(|row| row.omega_rad_s().to_bits() == actual.omega_rad_s().to_bits())
        );
        // The only convention change is complex conjugation, once.
        let matrix: Vec<_> = realized
            .load()
            .source_model()
            .eval_iw(actual.omega_rad_s())
            .iter()
            .map(|z| z.conj())
            .collect();
        assert!(relative_matrix(&matrix, actual.impedance()) <= options.fit.maximum_relative_error);
        let fitted = ModalMouthImpedance::new(
            actual.omega_rad_s(),
            source.patch().radius_m,
            *source.exterior(),
            ModeCount(2),
            matrix,
        )
        .unwrap();
        let raw_loop = mm_input_impedance_matrix(
            &duct,
            source.exterior(),
            actual.omega_rad_s(),
            fs_duct::LossModel::Lossless,
            actual,
            2,
            1,
        )
        .unwrap();
        let fitted_loop = mm_input_impedance_matrix(
            &duct,
            source.exterior(),
            actual.omega_rad_s(),
            fs_duct::LossModel::Lossless,
            &fitted,
            2,
            1,
        )
        .unwrap();
        maximum_loaded_error = maximum_loaded_error.max(relative_matrix(
            &fitted_loop.impedance_matrix,
            &raw_loop.impedance_matrix,
        ));
    }
    // A separate observed fixture check; local fit tolerance is no theorem
    // about an arbitrary lightly damped closed loop.
    assert!(
        maximum_loaded_error < 0.15,
        "loaded fixture error={maximum_loaded_error:e}"
    );
    let row = &record.held_out().source().rows[15];
    let gains = &record.held_out().source().volume_gains;
    assert!(gains[1].abs() < 0.25);
    assert!(row.impedance[1].abs() > 1e-5 * row.impedance[0].abs());
    let wrong: Vec<_> = row
        .impedance
        .iter()
        .enumerate()
        .map(|(i, z)| z.scale(1.0 / gains[i / 2] / gains[i % 2]))
        .collect();
    assert!(relative_matrix(&wrong, &row.impedance) > 1.0);

    let (mut load, retained) = realized.into_parts();
    let mut effort = [0.0; 2];
    load.step_mean_into(&[0.0, 1e-6], &mut effort).unwrap();
    assert!(
        effort[0].abs() > 1e-12,
        "higher-mode mutual reaction must be live"
    );
    assert!(load.pole_states().iter().any(|z| *z != C64::ZERO));
    load.step_mean_into(&[0.0, 0.0], &mut effort).unwrap();
    assert!(
        effort.iter().any(|p| p.abs() > 1e-12),
        "original radiation memory must continue"
    );
    assert_eq!(
        retained.report().selected_order,
        load.source_model().order()
    );
    eprintln!(
        "modal BEM memory: order={}, raw/held={:e}/{:e}, loaded={:e}, class={:?}",
        retained.report().selected_order,
        retained.report().maximum_fit_relative_error,
        retained.report().maximum_held_out_relative_error,
        maximum_loaded_error,
        retained.report().continuous.class
    );
}

#[test]
fn unsupported_clock_source_band_and_whole_sampling_budget_refuse_before_bem() {
    let (surface, patch) = geometry();
    let gate = CancelGate::new();
    let source = BemModalMouth::new(
        &surface,
        patch.clone(),
        air(),
        ModeCount(2),
        budget(),
        &gate,
    )
    .unwrap();
    let mut wrong = options();
    wrong.source_band_hz[1] = 2999.0;
    assert!(matches!(
        source.realize(&route(), wrong, &gate),
        Err(ModalMouthRealizationError::Invalid(_))
    ));
    wrong = options();
    wrong.source_band_hz[0] = 31.0;
    assert!(matches!(
        source.realize(&route(), wrong, &gate),
        Err(ModalMouthRealizationError::Invalid(_))
    ));
    wrong = options();
    wrong.operator_band_hz = [100.0, 3000.0];
    wrong.source_band_hz = wrong.operator_band_hz;
    assert!(matches!(
        source.realize(&route(), wrong, &gate),
        Err(ModalMouthRealizationError::Invalid(_))
    ));
    wrong = options();
    wrong.mechanics_rate_hz = 6001;
    assert!(matches!(
        source.realize(&route(), wrong, &gate),
        Err(ModalMouthRealizationError::Fit(_))
    ));
    wrong = options();
    wrong.maximum_reciprocity_defect = 0.1;
    assert!(matches!(
        source.realize(&route(), wrong, &gate),
        Err(ModalMouthRealizationError::Invalid(_))
    ));
    let mut small = budget();
    small.maximum_matrix_entries = 500; // Each sweep fits; their retained sum does not.
    let bounded = BemModalMouth::new(&surface, patch, air(), ModeCount(2), small, &gate).unwrap();
    assert!(matches!(
        bounded.realize(&route(), options(), &gate),
        Err(ModalMouthRealizationError::Invalid(_))
    ));
    gate.request();
    assert!(matches!(
        source.realize(&route(), options(), &gate),
        Err(ModalMouthRealizationError::Cancelled)
    ));
}
