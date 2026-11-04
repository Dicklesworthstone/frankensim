//! Physical local M/C/K tuning tangents of full complex mechanical modes.
//! Analytic references, original constraints, and independently refined nearby
//! physical pencils check the shared-border/gauge/scaling recovery together.
use fs_math::c64::C64;
use fs_modal::quadratic_modes::{
    PoleOptions, PoleResponseOptions, PreparedQuadraticModes, PreparedQuadraticSensitivities,
    QuadraticMode, QuadraticModeDerivative, QuadraticParameterDirection, QuadraticResponseError,
    QuadraticResponseOptions,
};
const I: C64 = C64::new(0.0, 1.0);
fn mode_options(rate: f64) -> PoleOptions {
    PoleOptions {
        spectral_scale: rate,
        step_tolerance: 2e-13,
        residual_tolerance: 2e-12,
        max_iterations: 40,
        max_dimension: 16,
        max_dense_bytes: 16 * 1024 * 1024,
        max_work: 100_000_000,
    }
}
fn response_options(rate: f64) -> QuadraticResponseOptions {
    QuadraticResponseOptions {
        numerical: PoleResponseOptions {
            spectral_scale: rate,
            pole_residual_tolerance: 3e-12,
            pole_step_tolerance: 3e-12,
            linear_residual_tolerance: 3e-12,
            max_dimension: 16,
            max_vectors: 8,
            max_dense_bytes: 16 * 1024 * 1024,
            max_work: 100_000_000,
        },
        max_observations: 4,
        max_physical_bytes: 16 * 1024 * 1024,
        max_physical_work: 100_000_000,
    }
}
fn direction<'a>(
    mass: &'a [f64],
    damping: &'a [f64],
    stiffness: &'a [f64],
) -> QuadraticParameterDirection<'a> {
    QuadraticParameterDirection {
        mass,
        damping,
        stiffness,
    }
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    let bound = if expected == C64::ZERO {
        tolerance
    } else {
        tolerance * expected.abs()
    };
    assert!(
        (actual - expected).abs() <= bound,
        "{actual:?} != {expected:?}; bound {bound}"
    );
}
fn check_receipt(owner: &PreparedQuadraticSensitivities<'_, '_>, d: &QuadraticModeDerivative) {
    assert!(d.linear_residual <= 3e-12);
    assert!(d.polynomial_relative_residual <= 3e-12);
    assert!(d.normalization_relative_residual <= 3e-12);
    assert_eq!(d.displacement[owner.phase_reference_coordinate()].im, 0.0);
    for i in 0..d.displacement.len() {
        close(
            d.velocity[i],
            owner.response().laplace_pole_per_s() * d.displacement[i]
                + d.laplace_frequency_per_s * owner.response().right_displacement()[i],
            2e-12,
        );
    }
}

#[test]
fn scalar_mass_damping_stiffness_and_common_scale_have_analytic_physical_tangents() {
    for c in [2.0, -2.0] {
        let (m, k) = (4.0, 36.0);
        let s = C64::new(
            -c / (2.0 * m),
            (k / m - (c / (2.0 * m)) * (c / (2.0 * m))).sqrt(),
        );
        let (mass, damping, stiffness) = ([m], [c], [k]);
        let model =
            PreparedQuadraticModes::new(&mass, &damping, &stiffness, 1, mode_options(3.0), || {
                false
            })
            .unwrap();
        let mode = model.refine(s, &[C64::ONE], || false).unwrap();
        let owner = model
            .prepare_sensitivities(&mode, response_options(3.0), || false)
            .unwrap();
        let (one, zero) = ([1.0], [0.0]);
        let dirs = [
            direction(&one, &zero, &zero),
            direction(&zero, &one, &zero),
            direction(&zero, &zero, &one),
            direction(&mass, &damping, &stiffness),
        ];
        let results = owner.parameter_derivatives(&dirs, || false).unwrap();
        let phi = C64::from_re(1.0 / m.sqrt());
        for (d, p) in results.iter().zip(&dirs) {
            let expected =
                -(s * s.scale(p.mass[0]) + s.scale(p.damping[0]) + C64::from_re(p.stiffness[0]))
                    / (s.scale(2.0 * m) + C64::from_re(c));
            close(d.laplace_frequency_per_s, expected, 2e-12);
            close(d.displacement[0], phi.scale(-0.5 * p.mass[0] / m), 2e-13);
            check_receipt(&owner, d);
        }
        // Scaling ALL physical coefficients changes the normalized shape even
        // though the pole and the pivot-one eigenvector do not move.
        close(results[3].laplace_frequency_per_s, C64::ZERO, 2e-12);
        close(results[3].displacement[0], phi.scale(-0.5), 2e-13);
        assert_eq!(owner.response().laplace_pole_per_s().re > 0.0, c < 0.0);
    }
}

#[test]
fn slow_overdamped_and_simple_zero_poles_keep_signed_parameter_motion() {
    for s in [-1.0, -2.0] {
        let model =
            PreparedQuadraticModes::new(&[1.0], &[3.0], &[2.0], 1, mode_options(2.0), || false)
                .unwrap();
        let mode = model
            .refine(C64::from_re(s), &[C64::ONE], || false)
            .unwrap();
        let owner = model
            .prepare_sensitivities(&mode, response_options(2.0), || false)
            .unwrap();
        let d = owner
            .parameter_derivatives(&[direction(&[0.0], &[0.0], &[1.0])], || false)
            .unwrap();
        close(
            d[0].laplace_frequency_per_s,
            C64::from_re(-1.0 / (2.0 * s + 3.0)),
            2e-13,
        );
        assert_eq!(d[0].displacement[0], C64::ZERO);
        check_receipt(&owner, &d[0]);
    }
    let model = PreparedQuadraticModes::new(&[1.0], &[2.0], &[0.0], 1, mode_options(2.0), || false)
        .unwrap();
    let mode = model.refine(C64::ZERO, &[C64::ONE], || false).unwrap();
    let owner = model
        .prepare_sensitivities(&mode, response_options(2.0), || false)
        .unwrap();
    let d = owner
        .parameter_derivatives(&[direction(&[0.0], &[0.0], &[1.0])], || false)
        .unwrap();
    close(d[0].laplace_frequency_per_s, C64::from_re(-0.5), 2e-13);
    close(d[0].velocity[0], C64::from_re(-0.5), 2e-13);
    check_receipt(&owner, &d[0]);
}

#[test]
fn nonnormal_direction_moves_the_pole_and_shape_despite_zero_right_quotient() {
    let m = [1.0, 0.0, 0.0, 1.0];
    let c = [0.4, 0.0, 0.0, 0.8];
    let k = [4.0, 7.0, 0.0, 9.0];
    let s = C64::new(-0.2, 3.96_f64.sqrt());
    let model = PreparedQuadraticModes::new(&m, &c, &k, 2, mode_options(3.0), || false).unwrap();
    let mode = model.refine(s, &[C64::ONE, C64::ZERO], || false).unwrap();
    let owner = model
        .prepare_sensitivities(&mode, response_options(3.0), || false)
        .unwrap();
    let zero = [0.0; 4];
    let dk = [0.0, 0.0, 1.0, 0.0];
    let negative = dk.map(|x| -x);
    let dirs = [
        direction(&zero, &zero, &dk),
        direction(&zero, &zero, &negative),
        direction(&zero, &zero, &zero),
    ];
    let forcing = [C64::ZERO, C64::ONE];
    let before = owner
        .response()
        .forcing_residues(&[&forcing], &[], || false)
        .unwrap();
    let d = owner.parameter_derivatives(&dirs, || false).unwrap();
    // p1(s)=s^2+.4s+4; p2(s)=s^2+.8s+9. dK21 produces
    // dphi2=-1/p2 and ds=7/(p1_s*p2), even though phi^H dK phi=0.
    let p2 = s * s + s.scale(0.8) + C64::from_re(9.0);
    let ds = C64::from_re(7.0) / ((s.scale(2.0) + C64::from_re(0.4)) * p2);
    close(d[0].laplace_frequency_per_s, ds, 2e-12);
    close(d[0].displacement[0], C64::ZERO, 2e-13);
    close(d[0].displacement[1], -C64::ONE / p2, 2e-12);
    assert!(d[0].laplace_frequency_per_s.abs() > 0.1);
    close(d[1].laplace_frequency_per_s, -ds, 2e-12);
    close(d[1].displacement[1], -d[0].displacement[1], 2e-12);
    assert_eq!(d[2].laplace_frequency_per_s, C64::ZERO);
    assert_eq!(d[2].displacement, vec![C64::ZERO; 2]);
    assert_eq!(d[2].velocity, vec![C64::ZERO; 2]);
    for result in &d {
        check_receipt(&owner, result);
    }
    assert_eq!(d, owner.parameter_derivatives(&dirs, || false).unwrap());
    let after = owner
        .response()
        .forcing_residues(&[&forcing], &[], || false)
        .unwrap();
    assert_eq!(before[0].displacement, after[0].displacement);
    close(before[0].displacement[0], -ds, 2e-12);
}

#[test]
fn signed_gyroscopic_coefficient_changes_frequency_without_fictitious_shape_or_decay() {
    let m = [2.0, 0.0, 0.0, 2.0];
    let k = [18.0, 0.0, 0.0, 18.0];
    let dc = [0.0, -1.0, 1.0, 0.0];
    let zero = [0.0; 4];
    for g in [-4.0_f64, 4.0] {
        let c = [0.0, -g, g, 0.0];
        let root = (g + (g * g + 144.0).sqrt()) / 4.0;
        let model =
            PreparedQuadraticModes::new(&m, &c, &k, 2, mode_options(4.0), || false).unwrap();
        let mode = model
            .refine(I.scale(root), &[C64::ONE, I], || false)
            .unwrap();
        let owner = model
            .prepare_sensitivities(&mode, response_options(4.0), || false)
            .unwrap();
        let d = owner
            .parameter_derivatives(&[direction(&zero, &dc, &zero)], || false)
            .unwrap();
        close(
            d[0].laplace_frequency_per_s,
            I.scale((1.0 + g / (g * g + 144.0).sqrt()) / 4.0),
            2e-12,
        );
        for value in &d[0].displacement {
            close(*value, C64::ZERO, 2e-12);
        }
        check_receipt(&owner, &d[0]);
    }
}

fn congruence(a: &[f64; 4]) -> [f64; 4] {
    // S=[[1,.5],[0,2]], q=S*z. Retain all original mutual mass.
    [
        a[0],
        0.5 * a[0] + 2.0 * a[1],
        0.5 * a[0] + 2.0 * a[2],
        0.25 * a[0] + a[1] + a[2] + 4.0 * a[3],
    ]
}
fn align(mode: &mut QuadraticMode, pivot: usize) {
    let z = mode.displacement[pivot];
    let phase = z.conj().scale(1.0 / z.abs());
    for q in &mut mode.displacement {
        *q = *q * phase;
    }
    mode.displacement[pivot].im = 0.0;
    for (v, &q) in mode.velocity.iter_mut().zip(&mode.displacement) {
        *v = mode.laplace_frequency_per_s * q;
    }
}

#[test]
fn full_mutual_mass_and_nonproportional_shape_match_independently_perturbed_pencils() {
    let base_m = [1.0, 0.0, 0.0, 1.0];
    let base_c = [1.0, 0.5, 0.5, 3.0];
    let base_k = [5.0, 2.5, 2.5, 5.0];
    let base_dm = [0.2, 0.1, 0.1, -0.3];
    let base_dc = [-0.2, 0.05, -0.15, 0.1];
    let base_dk = [0.4, -0.6, 0.3, 0.2];
    let (m, c, k) = (
        congruence(&base_m),
        congruence(&base_c),
        congruence(&base_k),
    );
    let (dm, dc, dk) = (
        congruence(&base_dm),
        congruence(&base_dc),
        congruence(&base_dk),
    );
    let s = C64::new(-1.0, 2.0);
    let model = PreparedQuadraticModes::new(&m, &c, &k, 2, mode_options(3.0), || false).unwrap();
    let mode = model
        .refine(s, &[C64::ONE - I.scale(0.25), I.scale(0.5)], || false)
        .unwrap();
    let owner = model
        .prepare_sensitivities(&mode, response_options(3.0), || false)
        .unwrap();
    let d = owner
        .parameter_derivatives(&[direction(&dm, &dc, &dk)], || false)
        .unwrap();
    let d = &d[0];
    let t: [C64; 4] = core::array::from_fn(|i| {
        s * s.scale(base_m[i]) + s.scale(base_c[i]) + C64::from_re(base_k[i])
    });
    let ts: [C64; 4] = core::array::from_fn(|i| s.scale(2.0 * base_m[i]) + C64::from_re(base_c[i]));
    let dt: [C64; 4] = core::array::from_fn(|i| {
        s * s.scale(base_dm[i]) + s.scale(base_dc[i]) + C64::from_re(base_dk[i])
    });
    let det_s = ts[0] * t[3] + t[0] * ts[3] - ts[1] * t[2] - t[1] * ts[2];
    let det_theta = dt[0] * t[3] + t[0] * dt[3] - dt[1] * t[2] - t[1] * dt[2];
    close(d.laplace_frequency_per_s, -det_theta / det_s, 3e-12);
    let step = 1e-5;
    let nearby = |theta: f64| {
        let mass = core::array::from_fn::<_, 4, _>(|i| m[i] + theta * dm[i]);
        let damping = core::array::from_fn::<_, 4, _>(|i| c[i] + theta * dc[i]);
        let stiffness = core::array::from_fn::<_, 4, _>(|i| k[i] + theta * dk[i]);
        let pencil =
            PreparedQuadraticModes::new(&mass, &damping, &stiffness, 2, mode_options(3.0), || {
                false
            })
            .unwrap();
        let mut mode = pencil
            .refine(
                s + d.laplace_frequency_per_s.scale(theta),
                owner.response().right_displacement(),
                || false,
            )
            .unwrap();
        align(&mut mode, owner.phase_reference_coordinate());
        mode
    };
    let lo = nearby(-step);
    let hi = nearby(step);
    for i in 0..2 {
        close(
            d.displacement[i],
            (hi.displacement[i] - lo.displacement[i]).scale(0.5 / step),
            5e-8,
        );
        close(
            d.velocity[i],
            (hi.velocity[i] - lo.velocity[i]).scale(0.5 / step),
            5e-8,
        );
    }
    let phi = owner.response().right_displacement();
    let mut normalization = C64::ZERO;
    for i in 0..2 {
        for j in 0..2 {
            normalization = normalization
                + (phi[i].conj() * d.displacement[j]).scale(2.0 * m[i * 2 + j])
                + (phi[i].conj() * phi[j]).scale(dm[i * 2 + j]);
        }
    }
    assert!(normalization.re.abs() < 2e-12);
    check_receipt(&owner, d);
}

#[test]
fn mixed_coordinate_rate_and_parameter_units_recover_original_physical_sensitivities() {
    let m = [1.0, 0.0, 0.0, 1.0];
    let c = [0.4, 0.0, 0.0, 0.8];
    let k = [4.0, 7.0, 0.0, 9.0];
    let dk = [0.0, 0.0, 1.0, 0.0];
    let zero = [0.0; 4];
    let s = C64::new(-0.2, 3.96_f64.sqrt());
    let p2 = s * s + s.scale(0.8) + C64::from_re(9.0);
    let ds = C64::from_re(7.0) / ((s.scale(2.0) + C64::from_re(0.4)) * p2);
    let dq = [C64::ZERO, -C64::ONE / p2];
    let dv = [ds, s * dq[1]];
    let coordinates = [2.0_f64.powi(300), -2.0_f64.powi(-300)];
    let map = |a: &[f64; 4]| {
        core::array::from_fn::<_, 4, _>(|i| a[i] * coordinates[i / 2] * coordinates[i % 2])
    };
    let (mm, cc, kk, dd) = (map(&m), map(&c), map(&k), map(&dk));
    let model = PreparedQuadraticModes::new(&mm, &cc, &kk, 2, mode_options(3.0), || false).unwrap();
    let mode = model.refine(s, &[C64::ONE, C64::ZERO], || false).unwrap();
    let owner = model
        .prepare_sensitivities(&mode, response_options(3.0), || false)
        .unwrap();
    let d = owner
        .parameter_derivatives(&[direction(&zero, &zero, &dd)], || false)
        .unwrap();
    close(d[0].laplace_frequency_per_s, ds, 2e-12);
    for i in 0..2 {
        close(d[0].displacement[i].scale(coordinates[i]), dq[i], 2e-12);
        close(d[0].velocity[i].scale(coordinates[i]), dv[i], 2e-12);
    }
    check_receipt(&owner, &d[0]);
    for exponent in [-500, 500] {
        let rate = 2.0_f64.powi(exponent);
        let mass_scale = 2.0_f64.powi(-exponent);
        let shape = 2.0_f64.powi(exponent / 2);
        let mm = m.map(|x| x * mass_scale);
        let kk = k.map(|x| x * rate);
        let dd = dk.map(|x| x * rate);
        let model =
            PreparedQuadraticModes::new(&mm, &c, &kk, 2, mode_options(3.0 * rate), || false)
                .unwrap();
        let mode = model
            .refine(s.scale(rate), &[C64::ONE, C64::ZERO], || false)
            .unwrap();
        let owner = model
            .prepare_sensitivities(&mode, response_options(3.0 * rate), || false)
            .unwrap();
        let d = owner
            .parameter_derivatives(&[direction(&zero, &zero, &dd)], || false)
            .unwrap();
        close(d[0].laplace_frequency_per_s.scale(1.0 / rate), ds, 3e-12);
        for i in 0..2 {
            close(d[0].displacement[i].scale(1.0 / shape), dq[i], 3e-12);
            close(d[0].velocity[i].scale(1.0 / shape / rate), dv[i], 3e-12);
        }
        check_receipt(&owner, &d[0]);
    }
    let model = PreparedQuadraticModes::new(&m, &c, &k, 2, mode_options(3.0), || false).unwrap();
    let mode = model.refine(s, &[C64::ONE, C64::ZERO], || false).unwrap();
    let owner = model
        .prepare_sensitivities(&mode, response_options(3.0), || false)
        .unwrap();
    for exponent in [-800, 800] {
        let parameter = 2.0_f64.powi(exponent);
        let dd = dk.map(|x| x * parameter);
        let d = owner
            .parameter_derivatives(&[direction(&zero, &zero, &dd)], || false)
            .unwrap();
        close(
            d[0].laplace_frequency_per_s.scale(1.0 / parameter),
            ds,
            3e-12,
        );
        close(d[0].displacement[1].scale(1.0 / parameter), dq[1], 3e-12);
        check_receipt(&owner, &d[0]);
    }
}

#[test]
fn direction_admission_resource_refusal_and_cancellation_do_not_mutate_the_factor() {
    let model = PreparedQuadraticModes::new(&[1.0], &[1.0], &[4.0], 1, mode_options(2.0), || false)
        .unwrap();
    let mode = model
        .refine(C64::new(-0.5, 3.75_f64.sqrt()), &[C64::ONE], || false)
        .unwrap();
    let owner = model
        .prepare_sensitivities(&mode, response_options(2.0), || false)
        .unwrap();
    let valid = direction(&[0.0], &[0.0], &[1.0]);
    let before = owner.parameter_derivatives(&[valid], || false).unwrap();
    assert!(owner.parameter_derivatives(&[], || false).is_err());
    assert!(owner.parameter_derivatives(&[valid; 9], || false).is_err());
    assert!(
        owner
            .parameter_derivatives(&[direction(&[], &[0.0], &[1.0])], || false)
            .is_err()
    );
    assert!(
        owner
            .parameter_derivatives(&[direction(&[0.0], &[f64::NAN], &[1.0])], || false)
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        owner.parameter_derivatives(&[valid, valid], || {
            polls += 1;
            polls > 12
        }),
        Err(QuadraticResponseError::Cancelled)
    ));
    assert_eq!(
        before,
        owner.parameter_derivatives(&[valid], || false).unwrap()
    );
    for kind in 0..4 {
        let mut options = response_options(2.0);
        match kind {
            0 => options.max_physical_bytes = 1,
            1 => options.max_physical_work = 1,
            2 => options.numerical.max_work = 1,
            _ => options.numerical.max_vectors = 0,
        };
        assert!(
            model
                .prepare_sensitivities(&mode, options, || false)
                .is_err()
        );
    }
    assert!(matches!(
        model.prepare_sensitivities(&mode, response_options(2.0), || true),
        Err(QuadraticResponseError::Cancelled)
    ));
    let m = [1.0, 0.0, 0.0, 1.0];
    let c = [0.4, 0.0, 0.0, 0.8];
    let k = [4.0, 7.0, 0.0, 9.0];
    let model = PreparedQuadraticModes::new(&m, &c, &k, 2, mode_options(3.0), || false).unwrap();
    let mode = model
        .refine(
            C64::new(-0.2, 3.96_f64.sqrt()),
            &[C64::ONE, C64::ZERO],
            || false,
        )
        .unwrap();
    let owner = model
        .prepare_sensitivities(&mode, response_options(3.0), || false)
        .unwrap();
    assert!(
        owner
            .parameter_derivatives(
                &[direction(&[0.0, 1.0, 0.0, 0.0], &[0.0; 4], &[0.0; 4])],
                || false
            )
            .is_err()
    );
    let critical =
        PreparedQuadraticModes::new(&[1.0], &[4.0], &[4.0], 1, mode_options(2.0), || false)
            .unwrap();
    let repeated = QuadraticMode {
        laplace_frequency_per_s: C64::from_re(-2.0),
        displacement: vec![C64::ONE],
        velocity: vec![C64::from_re(-2.0)],
        relative_residual: 0.0,
        solver_relative_residual: 0.0,
        relative_correction: 0.0,
        iterations: 0,
    };
    assert!(
        critical
            .prepare_sensitivities(&repeated, response_options(2.0), || false)
            .is_err()
    );
}
