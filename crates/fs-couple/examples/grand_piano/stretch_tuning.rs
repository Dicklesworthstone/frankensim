//! Opt-in 2:1 octave tuning of the actual pinned stiff-string scale.
//!
//! C4..B4 is an explicitly chosen first-partial equal-tempered temperament
//! octave, anchored at the requested A4. Outside it, each lower string's second
//! partial is matched to the upper string's first partial. This is a declared
//! tuning policy, not a measured Railsback curve or a perceptual calibration.
//! Only static tension changes; mass, length, EI, unison spread and all contact
//! geometry remain authoritative. Loaded-board poles and finite-rate numerical
//! dispersion are NOT fitted here. Do this before preparing any mechanics.

use super::geometry::Course;
use fs_math::det;
use std::f64::consts::{LN_2, PI, TAU};

const FIRST_KEY: u8 = 21;
const KEY_COUNT: usize = 88;
const TEMPERAMENT_FIRST: u8 = 60;
const TEMPERAMENT_LAST: u8 = 71;
const OCTAVE: usize = 12;

/// Retune a scale transactionally, preserving its order and every non-tension
/// field. A4 must be finite and in the front doors' existing 430..450 Hz range.
///
/// Partial scales are supported only when every supplied key has an unbroken
/// chain of supplied octave neighbours to C4..B4. Missing geometry is an error,
/// never a reason to synthesize a string or fall back to equal temperament.
/// No input course is changed if any validation or positive-tension solve fails.
///
/// # Errors
/// Returns an error for an invalid reference, invalid/duplicate course, missing
/// octave neighbour, or an unrepresentable positive-tension frequency match.
pub fn tune(courses: &mut [Course], concert_a_hz: f64) -> Result<(), String> {
    if !concert_a_hz.is_finite() || !(430.0..=450.0).contains(&concert_a_hz) {
        return Err("stretch tuning requires finite A4 in 430..450 Hz".into());
    }
    if courses.is_empty() {
        return Err("stretch tuning requires a nonempty scale".into());
    }

    // Course is Copy: all admission/solving is private until the final commit.
    let mut prepared = [None; KEY_COUNT];
    for course in courses.iter() {
        course.validate().map_err(|e| format!("stretch tuning key {}: {e}", course.midi))?;
        let index = usize::from(course.midi - FIRST_KEY);
        if prepared[index].is_some() {
            return Err(format!("stretch tuning has duplicate key {}", course.midi));
        }
        prepared[index] = Some(*course);
    }

    let first_anchor = usize::from(TEMPERAMENT_FIRST - FIRST_KEY);
    let after_anchor = usize::from(TEMPERAMENT_LAST - FIRST_KEY) + 1;
    for course in prepared[first_anchor..after_anchor].iter_mut().flatten() {
        let hz = concert_a_hz * det::exp((f64::from(course.midi) - 69.0) * LN_2 / 12.0);
        *course = at_partial(*course, 1, hz)?;
    }

    // Outward traversal guarantees that every neighbour has already been tuned.
    for index in after_anchor..KEY_COUNT {
        if let Some(course) = prepared[index] {
            let lower = prepared[index - OCTAVE].ok_or_else(|| format!(
                "stretch tuning key {} requires supplied octave neighbour {}",
                course.midi, course.midi - 12,
            ))?;
            prepared[index] = Some(at_partial(course, 1, lower.partial_hz(2, lower.tension_n))?);
        }
    }
    for index in (0..first_anchor).rev() {
        if let Some(course) = prepared[index] {
            let upper = prepared[index + OCTAVE].ok_or_else(|| format!(
                "stretch tuning key {} requires supplied octave neighbour {}",
                course.midi, course.midi + 12,
            ))?;
            prepared[index] = Some(at_partial(course, 2, upper.partial_hz(1, upper.tension_n))?);
        }
    }

    // This pass is infallible: all occupied input slots were admitted above.
    for course in courses {
        if let Some(tuned) = prepared[usize::from(course.midi - FIRST_KEY)] {
            course.tension_n = tuned.tension_n;
        }
    }
    Ok(())
}

/// Solve omega_n^2 = (T*k_n^2 + EI*k_n^4)/mu for T, not for an
/// independently assigned oscillator frequency. In particular EI is NOT scaled
/// along with tension. Positive tension is a physical admission, not a clamp.
fn at_partial(mut course: Course, partial: u8, target_hz: f64) -> Result<Course, String> {
    if !target_hz.is_finite() || target_hz <= 0.0 {
        return Err(format!("stretch tuning key {} has a nonfinite/invalid target", course.midi));
    }
    let k = f64::from(partial) * PI / course.length_m;
    let tension = (TAU * target_hz / k).powi(2) * course.linear_density_kg_m
        - course.flexural_rigidity_nm2 * k * k;
    if !tension.is_finite() || tension <= 0.0 {
        return Err(format!(
            "stretch tuning key {} partial {partial} at {target_hz} Hz requires nonpositive or nonfinite tension",
            course.midi,
        ));
    }
    course.tension_n = tension;
    course.validate().map_err(|e| format!("stretch tuning key {}: {e}", course.midi))?;
    // Guard catastrophic cancellation near the bending-only frequency floor.
    let actual_hz = course.partial_hz(usize::from(partial), tension);
    if !actual_hz.is_finite() || (actual_hz / target_hz - 1.0).abs() > 1e-10 {
        return Err(format!("stretch tuning key {} failed its partial-frequency residual", course.midi));
    }
    Ok(course)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scale() -> Vec<Course> {
        super::super::geometry::demonstration_scale().unwrap()
    }

    fn equal_hz(midi: u8, reference: f64) -> f64 {
        reference * 2.0_f64.powf((f64::from(midi) - 69.0) / 12.0)
    }

    fn fundamental(course: &Course) -> f64 {
        course.partial_hz(1, course.tension_n)
    }

    fn close(actual: f64, expected: f64) {
        assert!((actual / expected - 1.0).abs() < 1e-11, "{actual} != {expected}");
    }

    #[test]
    fn all_76_octaves_match_the_actual_second_and_first_partials() {
        let mut courses = scale();
        tune(&mut courses, 440.0).unwrap();
        for lower in &courses[..KEY_COUNT - OCTAVE] {
            let upper = &courses[usize::from(lower.midi - FIRST_KEY) + OCTAVE];
            close(lower.partial_hz(2, lower.tension_n), fundamental(upper));
        }
    }

    #[test]
    fn temperament_octave_and_concert_a_are_fixed_for_all_admitted_references() {
        for reference in [430.0, 440.0, 442.0, 450.0] {
            let mut courses = scale();
            tune(&mut courses, reference).unwrap();
            for course in courses.iter().filter(|c| (60..=71).contains(&c.midi)) {
                close(fundamental(course), equal_hz(course.midi, reference));
            }
            close(fundamental(&courses[48]), reference);
        }
    }

    #[test]
    fn zero_rigidity_recovers_equal_temperament_across_the_keyboard() {
        let mut courses = scale();
        for course in &mut courses {
            course.flexural_rigidity_nm2 = 0.0;
        }
        tune(&mut courses, 442.0).unwrap();
        for course in &courses {
            close(fundamental(course), equal_hz(course.midi, 442.0));
        }
    }

    #[test]
    fn positive_rigidity_flattens_bass_and_sharpens_treble() {
        let mut courses = scale();
        tune(&mut courses, 440.0).unwrap();
        for course in &courses {
            let hz = fundamental(course);
            if course.midi < TEMPERAMENT_FIRST {
                assert!(hz < equal_hz(course.midi, 440.0));
            } else if course.midi > TEMPERAMENT_LAST {
                assert!(hz > equal_hz(course.midi, 440.0));
            }
        }
    }

    #[test]
    fn only_tension_changes_and_repeated_preparation_does_not_accumulate_stretch() {
        let original = scale();
        let mut courses = original.clone();
        tune(&mut courses, 440.0).unwrap();
        for (old, new) in original.iter().zip(&courses) {
            assert_eq!(*old, Course { tension_n: old.tension_n, ..*new });
        }
        let once = courses.clone();
        tune(&mut courses, 440.0).unwrap();
        assert_eq!(courses, once);
    }

    #[test]
    fn unison_cents_still_retune_tension_without_scaling_bending() {
        let mut courses = scale();
        tune(&mut courses, 440.0).unwrap();
        for course in &courses {
            for cents in [-course.detune_cents, course.detune_cents] {
                let tension = course.tension_at_cents(cents).unwrap();
                let actual = 1200.0 * (course.partial_hz(1, tension) / fundamental(course)).log2();
                assert!((actual - cents).abs() < 1e-9);
            }
        }
    }

    #[test]
    fn arbitrary_input_order_is_preserved_and_does_not_affect_the_solve() {
        let mut forward = scale();
        let mut reversed = forward.iter().copied().rev().collect::<Vec<_>>();
        tune(&mut forward, 440.0).unwrap();
        tune(&mut reversed, 440.0).unwrap();
        for (a, b) in forward.iter().zip(reversed.iter().rev()) {
            assert_eq!(a, b);
        }
    }

    #[test]
    fn sparse_but_connected_octave_chains_are_supported() {
        let mut courses = scale().into_iter().filter(|c| [45, 57, 69, 81, 93].contains(&c.midi))
            .collect::<Vec<_>>();
        tune(&mut courses, 440.0).unwrap();
        for pair in courses.windows(2) {
            close(pair[0].partial_hz(2, pair[0].tension_n), fundamental(&pair[1]));
        }
    }

    #[test]
    fn missing_octave_geometry_refuses_transactionally_in_both_directions() {
        for missing in [57, 69, 81] {
            let mut courses = scale().into_iter().filter(|c| c.midi != missing).collect::<Vec<_>>();
            let before = courses.clone();
            let error = tune(&mut courses, 440.0).unwrap_err();
            assert!(error.contains(&format!("neighbour {missing}")), "{error}");
            assert_eq!(courses, before);
        }
    }

    #[test]
    fn duplicates_and_empty_scales_refuse() {
        let mut courses = scale();
        courses.push(courses[48]);
        let before = courses.clone();
        assert!(tune(&mut courses, 440.0).unwrap_err().contains("duplicate"));
        assert_eq!(courses, before);
        assert!(tune(&mut [], 440.0).is_err());
    }

    #[test]
    fn invalid_reference_never_changes_the_scale() {
        for reference in [f64::NAN, f64::INFINITY, 0.0, 429.9, 450.1] {
            let mut courses = scale();
            let before = courses.clone();
            assert!(tune(&mut courses, reference).is_err());
            assert_eq!(courses, before);
        }
    }

    #[test]
    fn invalid_key_and_physical_input_refuse_before_indexing_or_mutation() {
        for bad in [Course { midi: 0, ..scale()[0] }, Course { length_m: 0.0, ..scale()[0] }] {
            let mut courses = scale();
            courses[0] = bad;
            let before = courses.clone();
            assert!(tune(&mut courses, 440.0).is_err());
            assert_eq!(courses, before);
        }
    }

    #[test]
    fn impossible_bass_tension_rolls_back_even_after_the_treble_was_solved() {
        let mut courses = scale();
        courses[0].flexural_rigidity_nm2 *= 1e12;
        courses[0].validate().unwrap();
        let before = courses.clone();
        assert!(tune(&mut courses, 440.0).unwrap_err().contains("nonpositive"));
        assert_eq!(courses, before);
    }

    #[test]
    fn supplied_rigidity_drives_the_stretch_instead_of_a_key_number_curve() {
        let mut a = scale();
        let mut b = a.clone();
        b[36].flexural_rigidity_nm2 *= 3.0; // A3, below the temperament octave.
        tune(&mut a, 440.0).unwrap();
        tune(&mut b, 440.0).unwrap();
        assert!(fundamental(&b[36]) < fundamental(&a[36]));
        close(fundamental(&a[48]), fundamental(&b[48]));
        assert_eq!((a[36].flexural_rigidity_nm2 * 3.0).to_bits(),
            b[36].flexural_rigidity_nm2.to_bits());
    }
}
