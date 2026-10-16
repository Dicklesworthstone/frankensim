use super::*;

fn scale() -> Vec<Course> { super::super::super::geometry::demonstration_scale().unwrap() }
fn mixed() -> Profile {
    Profile { bass: OctaveWeights::new([1.0, 2.0, 3.0, 0.5]).unwrap(),
        treble: OctaveWeights::new([3.0, 1.0, 0.2, 0.0]).unwrap(), max_partial_hz: None }
}
fn close(a: f64, b: f64) {
    assert!((a-b).abs() <= 2e-10 * a.abs().max(b.abs()).max(1.0), "{a} != {b}");
}

#[test]
fn two_one_profile_is_the_original_bit_exact_tuning_path() {
    for reference in [430.0, 440.0, 442.0, 450.0] {
        let mut original = scale();
        let mut profiled = original.clone();
        super::super::tune(&mut original, reference).unwrap();
        let report = tune(&mut profiled, reference, Profile::default()).unwrap();
        assert_eq!(original, profiled);
        assert_eq!(report.len(), 76);
        for row in report {
            assert!(row.signed_beats_hz[0].unwrap().abs() < 1e-8);
            assert!(row.signed_beats_hz[1..].iter().all(Option::is_none));
        }
    }
}

#[test]
fn each_single_partial_policy_matches_its_actual_pair_on_both_scales() {
    for mut courses in [scale(), super::super::super::steinway_scale::courses().unwrap()] {
        for index in 0..4 {
            let mut weights = [0.0; 4]; weights[index] = 1.0;
            let weights = OctaveWeights::new(weights).unwrap();
            let report = tune(&mut courses, 442.0,
                Profile { bass: weights, treble: weights, max_partial_hz: None }).unwrap();
            for row in report {
                let lower = &courses[usize::from(row.lower_midi - 21)];
                let upper = &courses[usize::from(row.lower_midi - 9)];
                let n = index + 1;
                close(lower.partial_hz(2*n, lower.tension_n), upper.partial_hz(n, upper.tension_n));
                assert!(row.signed_beats_hz[index].unwrap().abs() < 1e-7);
            }
        }
    }
}

#[test]
fn mixed_policy_minimizes_the_declared_objective_not_a_fitted_frequency() {
    let mut courses = scale();
    let profile = mixed();
    let report = tune(&mut courses, 440.0, profile).unwrap();
    for row in report {
        let low = courses[usize::from(row.lower_midi - 21)];
        let high = courses[usize::from(row.lower_midi - 9)];
        let upper = high.midi > 71;
        let (variable, neighbour, weights) = if upper { (high, low, profile.treble) }
            else { (low, high, profile.bass) };
        let cost = |tension: f64| -> f64 {
            weights.0.iter().enumerate().map(|(i, w)| {
                let n = i+1;
                let actual = variable.partial_hz(if upper { n } else { 2*n }, tension);
                let target = neighbour.partial_hz(if upper { 2*n } else { n }, neighbour.tension_n);
                w * ((actual-target) / n as f64).powi(2)
            }).sum::<f64>() / weights.0.iter().sum::<f64>()
        };
        let best = cost(variable.tension_n);
        close(best, row.weighted_rms_hz.powi(2));
        // Every individual exact partial match is a different admissible
        // tension candidate. A mixed optimum must be at least as good as all.
        for (i, w) in weights.0.iter().enumerate() {
            if *w == 0.0 { continue; }
            let n = i+1;
            let target = neighbour.partial_hz(if upper { 2*n } else { n }, neighbour.tension_n);
            if let Ok(single) = at_partial(variable, (if upper { n } else { 2*n }) as u8, target) {
                assert!(best <= cost(single.tension_n) + 1e-10);
            }
        }
        for factor in [0.9999, 1.0001] {
            assert!(best <= cost(variable.tension_n * factor) + 1e-10);
        }
        let stationary: f64 = weights.0.iter().enumerate().map(|(i, w)| {
            let n = i+1;
            let target = neighbour.partial_hz(if upper { 2*n } else { n }, neighbour.tension_n);
            w * (1.0 - target / variable.partial_hz(if upper { n } else { 2*n }, variable.tension_n))
        }).sum();
        assert!(stationary.abs() < 1e-10);
    }
}

#[test]
fn physical_fields_anchor_order_and_idempotence_are_preserved() {
    let before = scale(); let mut courses = before.clone();
    tune(&mut courses, 442.0, mixed()).unwrap();
    for (old, new) in before.iter().zip(&courses) {
        assert_eq!(*old, Course { tension_n: old.tension_n, ..*new });
        if (60..=71).contains(&new.midi) {
            close(new.partial_hz(1, new.tension_n), 442.0 * 2.0_f64.powf((f64::from(new.midi)-69.0)/12.0));
        }
    }
    let once = courses.clone();
    tune(&mut courses, 442.0, mixed()).unwrap(); assert_eq!(courses, once);
    let mut reverse: Vec<_> = before.iter().copied().rev().collect();
    tune(&mut reverse, 442.0, mixed()).unwrap();
    assert_eq!(reverse, courses.iter().copied().rev().collect::<Vec<_>>());
}

#[test]
fn target_band_exclusions_are_explicit_and_never_silently_fall_back() {
    let mut courses: Vec<_> = scale().into_iter().filter(|c| [69,81].contains(&c.midi)).collect();
    let all = OctaveWeights::new([1.0;4]).unwrap();
    let report = tune(&mut courses, 440.0, Profile { bass: all, treble: all,
        max_partial_hz: Some(1000.0) }).unwrap();
    assert_eq!(report.len(), 1);
    assert!(report[0].signed_beats_hz[0].is_some());
    assert!(report[0].signed_beats_hz[1..].iter().all(Option::is_none));
    for band in [Some(100.0), Some(0.0), Some(f64::NAN), Some(f64::INFINITY)] {
        let before = courses.clone();
        assert!(tune(&mut courses, 440.0, Profile { bass: all, treble: all,
            max_partial_hz: band }).is_err());
        assert_eq!(courses, before);
    }
}

#[test]
fn weights_are_scale_invariant_and_the_flexible_string_limit_is_equal_temperament() {
    assert_eq!(OctaveWeights::new([1.,2.,3.,4.]).unwrap(),
        OctaveWeights::new([1e200,2e200,3e200,4e200]).unwrap());
    let mut courses = scale();
    for c in &mut courses { c.flexural_rigidity_nm2 = 0.0; }
    tune(&mut courses, 442.0, mixed()).unwrap();
    for c in &courses {
        close(c.partial_hz(1,c.tension_n), 442.0 * 2.0_f64.powf((f64::from(c.midi)-69.0)/12.0));
    }
}

#[test]
fn bad_weights_and_missing_geometry_refuse_without_partial_tuning() {
    for text in ["", "1,2,3", "1,2,3,4,5", "0,0,0,0", "1,-1,0,0", "NaN,1,0,0", "inf,0,0,0"] {
        assert!(OctaveWeights::parse(text).is_err(), "{text}");
    }
    assert!(OctaveWeights::new([f64::MAX, f64::MIN_POSITIVE, 0.,0.]).is_err());
    for missing in [57, 69, 81] {
        let mut courses: Vec<_> = scale().into_iter().filter(|c| c.midi != missing).collect();
        let before = courses.clone();
        assert!(tune(&mut courses, 440.0, mixed()).is_err()); assert_eq!(courses, before);
    }
}

#[test]
fn a_bending_only_minimum_is_not_clamped_into_a_tensioned_string() {
    let mut variable = scale()[60]; // an otherwise valid treble course
    variable.flexural_rigidity_nm2 = 1e6;
    variable.validate().unwrap();
    let pairs = [Pair { partial:1, target_hz:1.0, weight:1.0 },
        Pair { partial:2, target_hz:2.0, weight:1.0 }, Pair::default(), Pair::default()];
    assert!(minimize(variable, &pairs).unwrap_err().contains("nonpositive tension"));
}

#[test]
fn excluding_the_dominant_weight_does_not_erase_the_remaining_objective() {
    let before: Vec<_> = scale().into_iter().filter(|c| [69,81].contains(&c.midi)).collect();
    let mut tiny = before.clone(); let mut reference = before;
    let profile = Profile { treble: OctaveWeights::new([1e-300,2e-300,1.,0.]).unwrap(),
        max_partial_hz: Some(2000.0), ..Profile::default() };
    let tiny_report = tune(&mut tiny, 440.0, profile).unwrap();
    let ordinary = Profile { treble: OctaveWeights::new([1.,2.,0.,0.]).unwrap(), ..profile };
    let reference_report = tune(&mut reference, 440.0, ordinary).unwrap();
    close(tiny[1].tension_n, reference[1].tension_n);
    close(tiny_report[0].weighted_rms_hz, reference_report[0].weighted_rms_hz);
    assert!(tiny_report[0].weighted_rms_hz > 0.0);
}
