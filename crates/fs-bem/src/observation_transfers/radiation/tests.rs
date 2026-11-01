use super::*;
use crate::helmholtz::{Formulation, IterativeOptions};
use crate::near_field::Options;
use core::f64::consts::PI;

const MEDIUM: Medium = Medium {
    density: 1.2,
    sound_speed: 340.0,
};
const BUDGET: HarmonicWorkBudget = HarmonicWorkBudget {
    maximum_products: MAX_PRODUCTS,
};
const DIRECTIONS: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, -1.0]];

fn options(refinements: usize) -> PreparedRadiationOptions {
    PreparedRadiationOptions {
        angular: helmholtz::DirectivityProjectionOptions::default(),
        initial_latitudes: 8,
        refinements,
        maximum_kernel_evaluations: usize::MAX,
        maximum_source_entries: MAX_OBSERVATION_ENTRIES,
    }
}
fn close(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * b.abs().max(1e-30),
        "{a:e} != {b:e}"
    );
}
fn analytic(
    sources: &[SourceId],
    refinements: usize,
    mut field: impl FnMut(usize, [f64; 3]) -> C64,
) -> PreparedRadiationObservation {
    PreparedRadiationObservation::prepare::<TransferError>(
        680.0,
        MEDIUM,
        sources,
        ContentHash([11; 32]),
        &DIRECTIONS,
        options(refinements),
        sources.len(),
        |directions| {
            let mut values = Vec::new();
            for source in 0..sources.len() {
                values.extend(directions.iter().map(|&direction| field(source, direction)));
            }
            Ok(values)
        },
        &|| false,
    )
    .unwrap()
}

/// Independent full-sphere moment identities: integral n_i n_j=4pi delta/3,
/// and fourth moments=4pi(delta_ij delta_kl+...)/15. Force is N, flow is m³/s.
#[test]
fn g1_coincident_flow_force_power_and_directional_cross_moments_are_physical() {
    let omega = 680.0;
    let k = omega / MEDIUM.sound_speed;
    let ids = [
        SourceId::PointVolumeFlow("Q".into()),
        SourceId::PointForce("Fz".into()),
    ];
    let bank = analytic(&ids, 4, |source, direction| {
        C64::new(
            0.0,
            if source == 0 {
                -omega * MEDIUM.density / (4.0 * PI)
            } else {
                -k * direction[2] / (4.0 * PI)
            },
        )
    });
    let q = C64::new(0.001, 0.002);
    let f = C64::new(0.8, -0.3);
    let a = C64::new(0.0, -omega * MEDIUM.density / (4.0 * PI)) * q;
    let b = C64::new(0.0, -k / (4.0 * PI)) * f;
    let z = MEDIUM.density * MEDIUM.sound_speed;
    let monopole = omega * omega * MEDIUM.density * q.norm_sq() / (8.0 * PI * MEDIUM.sound_speed);
    let dipole = k * k * f.norm_sq() / (24.0 * PI * z);
    let mut out = [C64::ZERO; 3];
    let power = bank
        .evaluate_coherent(&[q, f], BUDGET, &mut out, || false)
        .unwrap();
    close(power.total_power_w, monopole + dipole, 3e-13);
    close(
        power.first_directional_moment_w[2],
        4.0 * PI * (a * b.conj()).re / (3.0 * z),
        3e-13,
    );
    assert!(
        power.first_directional_moment_w[..2]
            .iter()
            .all(|x| x.abs() < power.total_power_w * 1e-14)
    );
    close(
        power.second_directional_moment_w[0][0],
        monopole / 3.0 + dipole / 5.0,
        3e-13,
    );
    close(
        power.second_directional_moment_w[2][2],
        monopole / 3.0 + 3.0 * dipole / 5.0,
        3e-13,
    );
    for (index, direction) in DIRECTIONS.iter().enumerate() {
        assert!((out[index] - a - b.scale(direction[2])).abs() < 1e-15);
    }
    let mut density = [0.0; 3];
    let partial = bank
        .evaluate_coherence_modes(&[q, f, q, -f], 2, BUDGET, &mut density, || false)
        .unwrap();
    close(partial.total_power_w, 2.0 * (monopole + dipole), 3e-13);
    assert!(
        partial
            .first_directional_moment_w
            .iter()
            .all(|x| x.abs() < partial.total_power_w * 1e-14)
    );
    // Unitary rotation of coherence factors: same source covariance, same field statistics.
    let diagonal = bank
        .evaluate_independent_sources(
            &[2.0 * q.norm_sq(), 2.0 * f.norm_sq()],
            BUDGET,
            &mut [0.0; 3],
            || false,
        )
        .unwrap();
    close(diagonal.total_power_w, partial.total_power_w, 3e-14);
    for i in 0..3 {
        close(
            diagonal.second_directional_moment_w[i][i],
            partial.second_directional_moment_w[i][i],
            3e-14,
        );
        close(
            density[i],
            (a.norm_sq() + b.norm_sq() * DIRECTIONS[i][2].powi(2)) / z,
            3e-14,
        );
    }
    // A phase rotation and a change of coordinate normalization leave power unchanged.
    let phase = C64::new(det::cos(0.7), det::sin(0.7));
    close(
        bank.evaluate_coherent(&[q * phase, f * phase], BUDGET, &mut out, || false)
            .unwrap()
            .total_power_w,
        power.total_power_w,
        3e-14,
    );
    let scaled = analytic(&ids, 4, |source, d| {
        C64::new(
            0.0,
            if source == 0 {
                -omega * MEDIUM.density / (4.0 * PI)
            } else {
                -1e6 * k * d[2] / (4.0 * PI)
            },
        )
    });
    close(
        scaled
            .evaluate_coherent(&[q, f.scale(1e-6)], BUDGET, &mut out, || false)
            .unwrap()
            .total_power_w,
        power.total_power_w,
        3e-14,
    );
}

#[test]
fn g1_quiet_source_cancellation_survives_before_positive_power_reduction() {
    let ids: Vec<_> = (0..3)
        .map(|i| SourceId::PointForce(format!("Fz-{i}")))
        .collect();
    let bank = analytic(&ids, 4, |_, d| C64::new(0.0, -2.0 * d[2] / (4.0 * PI)));
    let mut field = [C64::ZERO; 3];
    let power = bank
        .evaluate_coherent(
            &[C64::from_re(1e16), C64::ONE, C64::from_re(-1e16)],
            BUDGET,
            &mut field,
            || false,
        )
        .unwrap();
    close(
        power.total_power_w,
        4.0 / (24.0 * PI * MEDIUM.density * MEDIUM.sound_speed),
        3e-13,
    );
    assert_eq!(field[0], C64::ZERO);
    assert_eq!(field[1], C64::new(0.0, -2.0 / (4.0 * PI)));
    let null = bank
        .evaluate_coherent(
            &[C64::ONE, -C64::ONE, C64::ZERO],
            BUDGET,
            &mut field,
            || false,
        )
        .unwrap();
    assert_eq!(null.total_power_w, 0.0);
    assert_eq!(field, [C64::ZERO; 3]);

    // Equal opposite separated monopoles: P/P_single=2(1-sinc(kd)).
    // At kd=1e-8, a rounded source Gram matrix loses this entire physical power.
    for kd in [1e-8, 2.3] {
        let ids = [
            SourceId::PointVolumeFlow("left".into()),
            SourceId::PointVolumeFlow("right".into()),
        ];
        let pair = analytic(&ids, 4, |source, d| {
            let phase = if source == 0 { 0.0 } else { -kd * d[0] };
            C64::new(det::cos(phase), det::sin(phase))
        });
        let p = pair
            .evaluate_coherent(&[C64::ONE, -C64::ONE], BUDGET, &mut field, || false)
            .unwrap();
        let factor = if kd < 0.1 {
            kd * kd / 3.0 - kd.powi(4) / 60.0
        } else {
            2.0 * (1.0 - det::sin(kd) / kd)
        };
        let single = 4.0 * PI / (2.0 * MEDIUM.density * MEDIUM.sound_speed);
        assert!(p.total_power_w > 0.0);
        close(p.total_power_w, single * factor, 2e-10);
    }
}

fn legendre(degree: usize, x: f64) -> f64 {
    let (mut previous, mut current) = (1.0, x);
    for n in 2..=degree {
        let next = ((2 * n - 1) as f64 * x * current - (n - 1) as f64 * previous) / n as f64;
        previous = current;
        current = next;
    }
    if degree == 0 { previous } else { current }
}

/// High angular content is not certified by source self-power at a low grid.
/// Actual contractions must converge on two refinements; low SH degrees never
/// remove omitted physical multipoles from total power.
#[test]
fn g1_unresolved_high_multipole_refuses_then_refinement_recovers_its_full_power() {
    let ids = [SourceId::StructuralMode("l20".into())];
    let coarse = analytic(&ids, 3, |_, d| C64::from_re(legendre(20, d[2])));
    assert!(matches!(
        coarse.evaluate_coherent(&[C64::ONE], BUDGET, &mut [C64::ZERO; 3], || false),
        Err(TransferError::AngularPowerNotConverged { .. })
    ));
    // Requested rays are exact retained evaluations even when this bank's
    // angular integration is too coarse for a full-sphere power claim.
    let mut rays = [C64::ZERO; 3];
    coarse
        .observe_coherent(&[C64::ONE], BUDGET, &mut rays, || false)
        .unwrap();
    assert_eq!(rays[1], C64::ONE);
    let mut density = [0.0; 3];
    coarse
        .observe_coherence_modes(&[C64::ONE], 1, BUDGET, &mut density, || false)
        .unwrap();
    close(
        density[1],
        1.0 / (2.0 * MEDIUM.density * MEDIUM.sound_speed),
        3e-14,
    );
    let fine = analytic(&ids, 6, |_, d| C64::from_re(legendre(20, d[2])));
    let p = fine
        .evaluate_coherent(&[C64::ONE], BUDGET, &mut [C64::ZERO; 3], || false)
        .unwrap();
    close(
        p.total_power_w,
        4.0 * PI / (41.0 * 2.0 * MEDIUM.density * MEDIUM.sound_speed),
        3e-12,
    );
    assert!(p.evaluated_directions > coarse.cached_angular_directions());
}

#[test]
fn g3_retained_dense_and_iterative_surface_fields_share_power_and_source_units() {
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let area: f64 = surface.areas().iter().sum();
    let sources = [
        SourceId::ApertureVolumeFlow("unit-flow".into()),
        SourceId::StructuralMode("signed-x".into()),
    ];
    let traces = [
        vec![1.0 / area; surface.areas().len()],
        surface.normals().iter().map(|n| n[0]).collect::<Vec<_>>(),
    ];
    let receiver = [ReceiverSpec {
        id: "field".into(),
        position_m: [0.02, 0.03, 0.08],
    }];
    let owner = BemObservation::new(&surface, MEDIUM, &sources, &traces, &receiver, 1e-3).unwrap();
    let omega = 680.0;
    let solved = helmholtz::solve_radiation_batch(
        &surface,
        2.0,
        MEDIUM,
        &owner.fields.iter().map(Vec::as_slice).collect::<Vec<_>>(),
        Formulation::BurtonMiller,
    )
    .unwrap();
    let gate = CancelGate::new_clock_free();
    let dense = owner
        .prepare_radiation_from_solutions(omega, &solved, &DIRECTIONS, options(4), &gate)
        .unwrap();
    let iterative = IterativeBoundaryObservation::for_surface(
        &surface,
        MEDIUM,
        &sources,
        &traces,
        &receiver,
        LocalBoundary::try_new_iterative(
            MEDIUM,
            vec![LocalWall::Rigid],
            vec![0; surface.areas().len()],
        )
        .unwrap(),
        IterativeObservationOptions {
            highest_omega: omega,
            minimum_clearance_m: 1e-3,
            receiver_quadrature: Options::default(),
            solver: IterativeOptions {
                leaf_capacity: 64,
                fmm_order: 4,
                tol: 1e-11,
                restart: 40,
                max_iters: 200,
                ..IterativeOptions::default()
            },
            lower_window: None,
        },
        &gate,
    )
    .unwrap();
    let response = iterative
        .solve_frequency(omega, 4 * surface.areas().len(), &gate)
        .unwrap();
    let bank = response
        .prepare_radiation(&DIRECTIONS, options(4), &gate)
        .unwrap();
    let drives = [C64::new(0.002, -0.001), C64::new(0.4, 0.3)];
    let mut expected = [C64::ZERO; 3];
    let mut actual = expected;
    let dense_power = dense
        .evaluate_coherent(&drives, BUDGET, &mut expected, || false)
        .unwrap();
    let power = bank
        .evaluate_coherent(&drives, BUDGET, &mut actual, || false)
        .unwrap();
    close(power.total_power_w, dense_power.total_power_w, 2e-7);
    assert_eq!(bank.sources(), sources.as_slice());
    for i in 0..3 {
        assert!((actual[i] - expected[i]).abs() < 2e-7 * expected[i].abs());
    }
    close(
        (0..3)
            .map(|i| power.second_directional_moment_w[i][i])
            .sum(),
        power.total_power_w,
        3e-14,
    );
    let kernel_work = bank.kernel_evaluations();
    for _ in 0..3 {
        let repeated = bank
            .evaluate_coherent(&drives, BUDGET, &mut actual, || false)
            .unwrap();
        assert_eq!(power, repeated);
        assert_eq!(bank.kernel_evaluations(), kernel_work);
    }
    let mut wrong = solved.clone();
    wrong[0].velocity[0] = C64::ZERO;
    assert!(
        owner
            .prepare_radiation_from_solutions(omega, &wrong, &DIRECTIONS, options(4), &gate)
            .is_err()
    );
}

#[test]
fn g0_work_memory_shape_and_cancellation_refuse_without_partial_preparation() {
    let ids = [SourceId::PointForce("one".into())];
    for limited in [
        PreparedRadiationOptions {
            maximum_kernel_evaluations: 1,
            ..options(4)
        },
        PreparedRadiationOptions {
            maximum_source_entries: 1,
            ..options(4)
        },
    ] {
        let mut evaluated = false;
        let result = PreparedRadiationObservation::prepare::<TransferError>(
            680.0,
            MEDIUM,
            &ids,
            ContentHash([3; 32]),
            &DIRECTIONS,
            limited,
            1,
            |_| {
                evaluated = true;
                Ok(Vec::new())
            },
            &|| false,
        );
        assert!(result.is_err());
        assert!(!evaluated);
    }
    let bank = analytic(&ids, 4, |_, _| C64::ONE);
    let mut out = [C64::ZERO; 3];
    assert!(
        bank.evaluate_coherent(
            &[C64::ONE],
            HarmonicWorkBudget {
                maximum_products: 1
            },
            &mut out,
            || false
        )
        .is_err()
    );
    assert!(
        bank.evaluate_coherent(&[C64::new(f64::NAN, 0.0)], BUDGET, &mut out, || false)
            .is_err()
    );
    assert!(
        bank.evaluate_coherent(&[], BUDGET, &mut out, || false)
            .is_err()
    );
    assert!(
        bank.evaluate_coherent(&[C64::ONE], BUDGET, &mut out[..1], || false)
            .is_err()
    );
    assert!(matches!(
        bank.evaluate_coherent(&[C64::ONE], BUDGET, &mut out, || true),
        Err(TransferError::Cancelled)
    ));
    assert!(
        bank.evaluate_independent_sources(&[-1.0], BUDGET, &mut [0.0; 3], || false)
            .is_err()
    );
    let silence = bank
        .evaluate_coherence_modes(&[], 0, BUDGET, &mut [0.0; 3], || false)
        .unwrap();
    assert_eq!(silence.total_power_w, 0.0);
}
