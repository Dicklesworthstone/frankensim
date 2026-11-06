use super::super::volume::EquilibriumVolumeEntropySolveV1;
use super::*;

fn controls() -> EquilibriumVolumeEntropySolveV1 {
    EquilibriumVolumeEntropySolveV1 {
        volume: EquilibriumVolumeSolveV1::default(),
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 2500.0,
        entropy_tolerance_j_k: 1e-7,
        maximum_temperature_iterations: 80,
    }
}

#[test]
fn g1_sv_inert_compression_matches_poisson_law_and_reversible_work() {
    let chemistry = owner(
        vec![model("A", 0.01, 3.5, 0.0, 0.0)],
        &["H"],
        vec![vec![1]],
    );
    let feed = [(id("A"), Amount::new(1.0))];
    let initial = chemistry
        .equilibrate_tv(
            Temperature::new(700.0),
            Volume::new(0.02),
            &feed,
            EquilibriumVolumeSolveV1::default(),
            || false,
        )
        .unwrap();
    let target = initial.thermodynamic().entropy_j_k();
    for ratio in [0.125_f64, 0.5, 1.0, 2.0, 8.0] {
        let state = chemistry
            .equilibrate_sv(
                Volume::new(0.02 * ratio),
                &feed,
                target,
                controls(),
                || false,
            )
            .unwrap();
        let expected_t = 700.0 * ratio.powf(-1.0 / 2.5);
        close(state.equilibrium().temperature().value(), expected_t, 2e-8);
        close(
            state.equilibrium().pressure().value(),
            R * expected_t / (0.02 * ratio),
            2e-8,
        );
        // The energy change is the reversible compression work. Its absolute
        // allowance also covers the zero-work return to the initial volume.
        let work = state.internal_energy().value() - initial.internal_energy().value();
        let expected_work = 2.5 * R * (expected_t - 700.0);
        assert!((work - expected_work).abs() <= 4e-7 * expected_t);
        assert!(
            state
                .thermodynamic()
                .entropy_residual_j_k()
                .unwrap()
                .abs()
                <= controls().entropy_tolerance_j_k
        );
        assert!(
            state.relative_volume_residual().abs() <= controls().volume.relative_volume_tolerance
        );
        assert!(state.target_internal_energy().is_none());
        assert!(state.equilibrium().species().iter().all(|s| s.source.is_some()));
    }
}

fn mixing_entropy(t: f64, q: f64) -> f64 {
    let b = 1.0 / (1.0 + (-q / (R * t)).exp());
    -b * b.ln() - (1.0 - b) * (1.0 - b).ln()
}

#[test]
fn g1_sv_reacting_compression_matches_independent_entropy_root() {
    let q = 20000.0;
    let chemistry = isomers(q);
    let feed = [(id("A"), Amount::new(1.0))];
    let initial = chemistry
        .equilibrate_tv(
            Temperature::new(1000.0),
            Volume::new(0.02),
            &feed,
            EquilibriumVolumeSolveV1::default(),
            || false,
        )
        .unwrap();
    let target = initial.thermodynamic().entropy_j_k();
    let ratio = 0.25_f64;
    // Equal source Cp and entropy offsets; reaction equilibrium changes only
    // the ideal-mixing entropy. The scalar oracle has no chemical Newton solve.
    let mut low = 200.0_f64;
    let mut high = 2500.0_f64;
    for _ in 0..90 {
        let t = low + 0.5 * (high - low);
        let ds = 3.0 * (t / 1000.0).ln() + ratio.ln() + mixing_entropy(t, q)
            - mixing_entropy(1000.0, q);
        if ds < 0.0 {
            low = t;
        } else {
            high = t;
        }
    }
    let expected = low + 0.5 * (high - low);
    let state = chemistry
        .equilibrate_sv(
            Volume::new(0.02 * ratio),
            &feed,
            target,
            controls(),
            || false,
        )
        .unwrap();
    close(state.equilibrium().temperature().value(), expected, 2e-8);
    close(
        species(state.equilibrium(), "B").amount.value(),
        1.0 / (1.0 + (-q / (R * expected)).exp()),
        2e-9,
    );
    assert!(expected < 1000.0 * ratio.powf(-1.0 / 3.0) - 20.0);
    let returned = chemistry
        .equilibrate_sv(
            initial.volume(),
            &feed,
            state.thermodynamic().entropy_j_k(),
            controls(),
            || false,
        )
        .unwrap();
    close(returned.equilibrium().temperature().value(), 1000.0, 3e-8);
    close(
        returned.equilibrium().pressure().value(),
        initial.equilibrium().pressure().value(),
        3e-8,
    );
}

#[test]
fn g3_sv_dissociation_matches_reacting_sound_speed_and_extensive_scaling() {
    let chemistry = dissociation();
    let feed = [(id("A2"), Amount::new(1.0))];
    let volume = 0.01;
    let initial = chemistry
        .equilibrate_tv(
            Temperature::new(1000.0),
            Volume::new(volume),
            &feed,
            EquilibriumVolumeSolveV1::default(),
            || false,
        )
        .unwrap();
    let entropy = initial.thermodynamic().entropy_j_k();
    let delta = 1e-3;
    let mut pressure = [0.0; 2];
    for (index, ratio) in [1.0 - delta, 1.0 + delta].into_iter().enumerate() {
        let state = chemistry
            .equilibrate_sv(
                Volume::new(volume * ratio),
                &feed,
                entropy,
                controls(),
                || false,
            )
            .unwrap();
        pressure[index] = state.equilibrium().pressure().value();
    }
    let dp_dv = (pressure[1] - pressure[0]) / (2.0 * delta * volume);
    let c_squared = -dp_dv * volume * volume / initial.equilibrium().mass().value();
    close(
        c_squared,
        initial.thermodynamic().response().sound_speed_m_s.powi(2),
        6e-6,
    );
    for scale in [1e-6, 1e6] {
        let mut limits = controls();
        limits.entropy_tolerance_j_k *= scale;
        let state = chemistry
            .equilibrate_sv(
                Volume::new(volume * scale),
                &[(id("A2"), Amount::new(scale))],
                entropy * scale,
                limits,
                || false,
            )
            .unwrap();
        close(state.equilibrium().temperature().value(), 1000.0, 2e-8);
        close(
            state.equilibrium().pressure().value(),
            initial.equilibrium().pressure().value(),
            2e-8,
        );
        close(
            state.internal_energy().value(),
            initial.internal_energy().value() * scale,
            3e-8,
        );
    }
}

#[test]
fn g4_sv_refuses_unclosed_entropy_invalid_sources_and_late_cancellation() {
    let chemistry = dissociation();
    let feed = [(id("A2"), Amount::new(1.0))];
    let initial = chemistry
        .equilibrate_tv(
            Temperature::new(1000.0),
            Volume::new(0.01),
            &feed,
            EquilibriumVolumeSolveV1::default(),
            || false,
        )
        .unwrap();
    let target = initial.thermodynamic().entropy_j_k();
    let mut limits = controls();
    limits.maximum_temperature_iterations = 1;
    assert!(matches!(
        chemistry.equilibrate_sv(initial.volume(), &feed, target, limits, || false),
        Err(EquilibriumVolumeErrorV1::Response(
            EquilibriumExpansionErrorV1::EntropyNoConvergence { .. }
        ))
    ));
    assert!(matches!(
        chemistry.equilibrate_sv(initial.volume(), &feed, target + 1e6, controls(), || false),
        Err(EquilibriumVolumeErrorV1::Response(
            EquilibriumExpansionErrorV1::EntropyNotBracketed { .. }
        ))
    ));
    limits = controls();
    limits.entropy_tolerance_j_k = f64::NAN;
    assert!(matches!(
        chemistry.equilibrate_sv(initial.volume(), &feed, target, limits, || false),
        Err(EquilibriumVolumeErrorV1::Invalid(_))
    ));
    assert!(matches!(
        chemistry.equilibrate_sv(initial.volume(), &feed, f64::NAN, controls(), || false),
        Err(EquilibriumVolumeErrorV1::Invalid(_))
    ));
    assert!(matches!(
        chemistry.equilibrate_sv(Volume::new(0.0), &feed, target, controls(), || false),
        Err(EquilibriumVolumeErrorV1::Invalid(_))
    ));
    limits = controls();
    limits.minimum_temperature_k = 100.0;
    assert!(matches!(
        chemistry.equilibrate_sv(initial.volume(), &feed, target, limits, || false),
        Err(EquilibriumVolumeErrorV1::Response(
            EquilibriumExpansionErrorV1::Invalid(_)
        ))
    ));
    let mut polls = 0;
    chemistry
        .equilibrate_sv(initial.volume(), &feed, target, controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut seen = 0;
    let cancelled = chemistry.equilibrate_sv(initial.volume(), &feed, target, controls(), || {
        seen += 1;
        seen >= polls
    });
    assert!(matches!(
        cancelled,
        Err(EquilibriumVolumeErrorV1::Equilibrium(
            CompositionEquilibriumErrorV1::Cancelled
        )) | Err(EquilibriumVolumeErrorV1::Response(
            EquilibriumExpansionErrorV1::Equilibrium(CompositionEquilibriumErrorV1::Cancelled)
        ))
    ));
    assert_eq!(feed[0].1.value(), 1.0);
}
