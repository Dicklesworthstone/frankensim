//! G1/G3 substrate heating retains the original spatial coherence factors.
use super::*;

fn bank<'b, 's, 'm>(
    factors: &'b [PlanarAngularSpectrum<'s, 'm>],
) -> AngularSpectrumCoherence<'b, 's, 'm> {
    AngularSpectrumCoherence::try_new(&factors[0], factors, coherence_limits(), || false).unwrap()
}

// Independent TE interface transmission and outgoing plane-wave propagation.
// The production owner instead resolves a complete vector Fourier grid.
fn transmitted_te(epsilon: C64, kx: f64, depth: f64) -> C64 {
    let omega = 1e9;
    let incident_kz = det::sqrt(omega * omega * E * M - kx * kx);
    let kz = (epsilon.scale(omega * omega * M) - C64::from_re(kx * kx)).sqrt();
    assert!(kz.re > 0.0 && kz.im < 0.0);
    C64::from_re(2.0 * incident_kz) / (C64::from_re(incident_kz) + kz)
        * phase(-kz.re * depth).scale(det::exp(kz.im * depth))
}

#[test]
fn g1_substrate_coherence_preserves_lossy_fringes_and_integrated_watts() {
    let incident = material(E, M);
    let lossy = finite_layers::debye(2.0 * E, E, 1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&incident, DielectricResponse::Debye(&lossy), 1e9, &[])
            .unwrap();
    let g = grid();
    let kx = core::f64::consts::TAU / (g.width as f64 * g.spacing_m[0]);
    let epsilon = stack.exit_permittivity_f_m();
    // The helper's final arguments are area/gap; its pole has tau=1e-9 s.
    // At omega*tau=1 the substrate has substantial, explicitly known loss.
    complex_near(epsilon, C64::new(2.5 * E, -0.5 * E), 1e-25);
    let depth = 0.23;
    let u = transmitted_te(epsilon, 0.0, depth);
    let v = transmitted_te(epsilon, kx, depth);
    let mut previous_watts = None;
    for rho in [0.0, 0.6, 1.0] {
        let mut a = Vec::new();
        let mut b = Vec::new();
        for x in 0..g.width {
            let p = phase(-kx * x as f64 * g.spacing_m[0]);
            a.push([C64::ZERO, C64::ONE + p.scale(rho)]);
            b.push([C64::ZERO, p.scale(det::sqrt(1.0 - rho * rho))]);
        }
        let factors = [prepare(&stack, &a), prepare(&stack, &b)];
        let bank = bank(&factors);
        let plane = bank.transmitted_absorption_field(depth, || false).unwrap();
        for x in 0..g.width {
            let p = v * phase(-kx * x as f64 * g.spacing_m[0]);
            let expected = -0.5e9
                * epsilon.im
                * ((u + p.scale(rho)).norm_sq() + (1.0 - rho * rho) * p.norm_sq());
            close(plane.absorbed_power_density_w_m3()[x], expected);
        }
        let area = g.width as f64 * g.spacing_m[0] * g.spacing_m[1];
        close(
            plane.absorbed_power_per_depth_w_m(),
            -0.5e9 * epsilon.im * (u.norm_sq() + v.norm_sq()) * area,
        );
        let watts = bank
            .transmitted_absorbed_power_w(0.0, f64::INFINITY, || false)
            .unwrap();
        close(watts, bank.power().transmitted_w);
        if let Some(previous) = previous_watts {
            close(watts, previous);
        }
        previous_watts = Some(watts);
        assert!(std::ptr::eq(bank.source(), &stack));
    }
}

#[test]
fn g1_substrate_heat_distinguishes_coherent_cancellation_from_independent_factors() {
    let incident = material(E, M);
    let lossy = finite_layers::debye(2.0 * E, E, 1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&incident, DielectricResponse::Debye(&lossy), 1e9, &[])
            .unwrap();
    let input = vec![[C64::new(0.7, 0.2), C64::new(0.3, -0.6)]; grid().width];
    let opposite: Vec<_> = input.iter().map(|e| e.map(|v| -v)).collect();
    let cancelled: Vec<_> = input
        .iter()
        .zip(&opposite)
        .map(|(a, b)| std::array::from_fn(|c| a[c] + b[c]))
        .collect();
    let factors = [prepare(&stack, &input), prepare(&stack, &opposite)];
    let single = bank(&factors[..1]);
    let independent = bank(&factors);
    let zero = [prepare(&stack, &cancelled)];
    let dark = bank(&zero);
    let coherent = factors[0]
        .transmitted_absorption_field(0.2, limits(), || false)
        .unwrap();
    let rank_one = single.transmitted_absorption_field(0.2, || false).unwrap();
    assert_eq!(
        rank_one.absorbed_power_density_w_m3(),
        coherent.absorbed_power_density_w_m3()
    );
    let two = independent
        .transmitted_absorption_field(0.2, || false)
        .unwrap();
    for (&actual, &one) in two
        .absorbed_power_density_w_m3()
        .iter()
        .zip(rank_one.absorbed_power_density_w_m3())
    {
        assert!(actual > 0.0);
        close(actual, 2.0 * one);
    }
    for (begin, end) in [(0.0, 0.4), (0.2, f64::INFINITY)] {
        let one = factors[0]
            .transmitted_absorbed_power_w(begin, end, || false)
            .unwrap();
        assert_eq!(
            single
                .transmitted_absorbed_power_w(begin, end, || false)
                .unwrap(),
            one
        );
        close(
            independent
                .transmitted_absorbed_power_w(begin, end, || false)
                .unwrap(),
            2.0 * one,
        );
        assert_eq!(
            dark.transmitted_absorbed_power_w(begin, end, || false)
                .unwrap(),
            0.0
        );
    }
    assert!(
        dark.transmitted_absorption_field(0.2, || false)
            .unwrap()
            .absorbed_power_density_w_m3()
            .iter()
            .all(|&q| q == 0.0)
    );
}

#[test]
fn g3_complex_factor_mixing_preserves_substrate_heat_and_film_energy_balance() {
    let incident = material(E, M);
    let lossy = finite_layers::debye(2.0 * E, E, 1.0, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &incident,
        DielectricResponse::Debye(&lossy),
        1e9,
        &[DielectricLayer {
            response: DielectricResponse::Debye(&lossy),
            thickness_m: 0.23,
        }],
    )
    .unwrap();
    let (mut a, mut b, mut ua, mut ub) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let j = C64::new(0.0, 1.0);
    for x in 0..grid().width {
        let p = phase(-core::f64::consts::TAU * x as f64 / grid().width as f64);
        let av = [C64::new(0.8, 0.2) + p.scale(0.3), j];
        let bv = [p * C64::new(0.2, -0.6), p.scale(0.4) + C64::new(0.2, 0.1)];
        ua.push(std::array::from_fn(|c| {
            (av[c] + j * bv[c]).scale(core::f64::consts::FRAC_1_SQRT_2)
        }));
        ub.push(std::array::from_fn(|c| {
            (j * av[c] + bv[c]).scale(core::f64::consts::FRAC_1_SQRT_2)
        }));
        a.push(av);
        b.push(bv);
    }
    let factors = [prepare(&stack, &a), prepare(&stack, &b)];
    let rotated = [prepare(&stack, &ua), prepare(&stack, &ub)];
    let original = bank(&factors);
    let mixed = bank(&rotated);
    let plane = original
        .transmitted_absorption_field(0.19, || false)
        .unwrap();
    let changed = mixed.transmitted_absorption_field(0.19, || false).unwrap();
    same_field(plane.field(), changed.field());
    same_field(
        plane.field(),
        &original.field(Side::Transmitted, 0.19, || false).unwrap(),
    );
    for ((point, &heat), &other) in plane
        .field()
        .points()
        .iter()
        .zip(plane.absorbed_power_density_w_m3())
        .zip(changed.absorbed_power_density_w_m3())
    {
        close(heat, other);
        close(
            heat,
            -1e9 * stack.exit_permittivity_f_m().im
                * (0..3)
                    .map(|c| point.electric_coherency_v2_m2[c][c].re)
                    .sum::<f64>(),
        );
    }
    let begin = 0.11;
    let end = 0.73;
    let watts = original
        .transmitted_absorbed_power_w(begin, end, || false)
        .unwrap();
    // Independent TE/TM interface admittances and exponential-depth integrals
    // for the four authored (factor, transverse mode) channels above.
    close(watts, 0.224_091_889_160_567_9);
    close(
        watts,
        mixed
            .transmitted_absorbed_power_w(begin, end, || false)
            .unwrap(),
    );
    let entering = original.field(Side::Transmitted, begin, || false).unwrap();
    let leaving = original.field(Side::Transmitted, end, || false).unwrap();
    close(watts, entering.normal_power_w() - leaving.normal_power_w());
    let intervals = 64;
    let mut integral = 0.0;
    for i in 0..=intervals {
        let weight = if i == 0 || i == intervals {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        integral += weight
            * original
                .transmitted_absorption_field(
                    begin + (end - begin) * i as f64 / intervals as f64,
                    || false,
                )
                .unwrap()
                .absorbed_power_per_depth_w_m();
    }
    near(
        integral * (end - begin) / (3.0 * intervals as f64),
        watts,
        watts * 3e-8,
    );
    let tail = original
        .transmitted_absorbed_power_w(0.0, f64::INFINITY, || false)
        .unwrap();
    let power = original.power();
    assert!(power.film_absorbed_w > 0.0);
    close(tail, power.transmitted_w);
    close(
        power.reflected_w + power.film_absorbed_w + tail,
        power.incident_w,
    );
}

#[test]
fn g1_lossless_evanescent_coherence_deposits_no_substrate_heat() {
    let incident = material(4.0, 1.0);
    let exit = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&incident, DielectricResponse::Lossless(&exit), 1.0, &[])
            .unwrap();
    let grid = AngularSpectrumGrid {
        width: 2,
        height: 1,
        spacing_m: [2.0, 1.0],
    };
    let a = [[C64::ONE, C64::ONE], [-C64::ONE, -C64::ONE]];
    let b = a.map(|e| e.map(|v| v * C64::new(0.2, 0.3)));
    let factors = [
        PlanarAngularSpectrum::prepare(&stack, grid, &a, limits(), || false).unwrap(),
        PlanarAngularSpectrum::prepare(&stack, grid, &b, limits(), || false).unwrap(),
    ];
    let bank = bank(&factors);
    let surface = bank.transmitted_absorption_field(0.0, || false).unwrap();
    let deeper = bank.transmitted_absorption_field(1.0, || false).unwrap();
    assert!(
        deeper.field().points()[0].electric_coherency_v2_m2[0][0].re
            < surface.field().points()[0].electric_coherency_v2_m2[0][0].re
    );
    assert_eq!(surface.absorbed_power_density_w_m3(), &[0.0, 0.0]);
    assert_eq!(deeper.absorbed_power_per_depth_w_m(), 0.0);
    assert_eq!(
        bank.transmitted_absorbed_power_w(0.0, f64::INFINITY, || false)
            .unwrap(),
        0.0
    );
}

#[test]
fn g0_g4_substrate_queries_obey_bank_admission_darkness_and_transactional_refusal() {
    let incident = material(E, M);
    let lossy = finite_layers::debye(2.0 * E, E, 1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&incident, DielectricResponse::Debye(&lossy), 1e9, &[])
            .unwrap();
    let input = vec![[C64::ONE, C64::ONE]; grid().width];
    let factors = [prepare(&stack, &input), prepare(&stack, &input)];
    let bank = bank(&factors);
    let tight = AngularSpectrumCoherenceLimits {
        maximum_modes: factors.len(),
        maximum_point_mode_pairs: factors.len() * grid().width,
        maximum_workspace_bytes: bank.workspace_bytes(),
        maximum_fft_work: bank.field_fft_work(),
    };
    let admitted =
        AngularSpectrumCoherence::try_new(&factors[0], &factors, tight, || false).unwrap();
    assert!(admitted.transmitted_absorption_field(0.2, || false).is_ok());
    assert!(
        AngularSpectrumCoherence::try_new(
            &factors[0],
            &factors,
            AngularSpectrumCoherenceLimits {
                maximum_workspace_bytes: tight.maximum_workspace_bytes - 1,
                ..tight
            },
            || false,
        )
        .is_err()
    );
    let dark = AngularSpectrumCoherence::try_new(
        &factors[0],
        &[],
        AngularSpectrumCoherenceLimits {
            maximum_modes: 0,
            maximum_point_mode_pairs: 0,
            maximum_fft_work: 0,
            ..tight
        },
        || false,
    )
    .unwrap();
    let zero = dark.transmitted_absorption_field(0.2, || false).unwrap();
    assert!(zero.absorbed_power_density_w_m3().iter().all(|&q| q == 0.0));
    assert_eq!(zero.field().normal_power_w(), 0.0);
    assert_eq!(
        dark.transmitted_absorbed_power_w(0.0, f64::INFINITY, || false)
            .unwrap(),
        0.0
    );
    for source in [&bank, &dark] {
        for depth in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(
                source
                    .transmitted_absorption_field(depth, || false)
                    .is_err()
            );
        }
        for (begin, end) in [
            (f64::NAN, 1.0),
            (0.0, f64::NAN),
            (-1.0, 0.0),
            (1.0, 0.0),
            (f64::INFINITY, f64::INFINITY),
            (0.0, f64::NEG_INFINITY),
        ] {
            assert!(
                source
                    .transmitted_absorbed_power_w(begin, end, || false)
                    .is_err()
            );
        }
        assert_eq!(
            source
                .transmitted_absorbed_power_w(0.2, 0.2, || false)
                .unwrap(),
            0.0
        );
        assert!(matches!(
            source.transmitted_absorption_field(0.2, || true),
            Err(AngularSpectrumError::Cancelled)
        ));
        assert!(matches!(
            source.transmitted_absorbed_power_w(0.0, 1.0, || true),
            Err(AngularSpectrumError::Cancelled)
        ));
    }
    let before = bank.transmitted_absorption_field(0.2, || false).unwrap();
    let mut polls = 0;
    assert!(matches!(
        bank.transmitted_absorption_field(0.2, || {
            polls += 1;
            polls > 10
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    let replay = bank.transmitted_absorption_field(0.2, || false).unwrap();
    assert_eq!(
        before.absorbed_power_density_w_m3(),
        replay.absorbed_power_density_w_m3()
    );
    let mut polls = 0;
    assert!(matches!(
        bank.transmitted_absorbed_power_w(0.0, 1.0, || {
            polls += 1;
            polls > 4
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert!(
        bank.transmitted_absorbed_power_w(0.0, 1.0, || false)
            .unwrap()
            > 0.0
    );
}
