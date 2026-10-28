//! Project the native solver's sampled transient peak into a probability study.
//! The producer includes the initial state and every accepted fine-grid state.
//! This observation carries no enclosure between the retained time samples.

use crate::json_read::JsonValue as J;

fn count(value: &J, name: &str) -> Result<usize, &'static str> {
    value.get(name).and_then(J::number_raw)
        .and_then(|raw| raw.parse::<usize>().ok())
        .ok_or("sampled peak has an invalid trajectory count")
}

/// Read only a complete observation bound to the same sealed run and region as
/// the final QoI. An old producer without peak coverage must refuse; replacing
/// the peak with the final value would change the probability event.
pub(super) fn temperature(
    receipt: &J, run: &str, region: &str, final_k: f64,
) -> Result<f64, &'static str> {
    if receipt.str_field("run") != Some(run) || !final_k.is_finite() || final_k < 0.0 {
        return Err("sampled peak belongs to another run or has an invalid final temperature");
    }
    let time = receipt.get("transient")
        .ok_or("temperature-peak requires retained native transient evidence")?;
    let peak = time.get("sampled_peak")
        .ok_or("native trajectory has no sampled peak including the initial state")?;
    let steps = count(time, "fine_steps")?;
    let samples = count(peak, "samples")?;
    let state = count(peak, "state")?;
    let horizon = time.f64_field("final_time_s")
        .filter(|v| v.is_finite() && *v > 0.0)
        .ok_or("sampled peak has an invalid trajectory horizon")?;
    let at = peak.f64_field("time_s")
        .filter(|v| v.is_finite() && *v >= 0.0 && *v <= horizon)
        .ok_or("sampled peak time is outside its completed trajectory")?;
    let value = peak.f64_field("value_k")
        .filter(|v| v.is_finite() && *v >= 0.0)
        .ok_or("sampled peak is not finite nonnegative Kelvin")?;
    if time.str_field("status") != Some("completed")
        || peak.str_field("region") != Some(region)
        || peak.str_field("time_scope") != Some("accepted-fine-grid")
        || peak.str_field("authority") != Some("Estimated")
        || peak.get("includes_initial") != Some(&J::Bool(true))
        || peak.get("continuous_time_bound") != Some(&J::Bool(false))
        || steps == 0 || steps.checked_add(1) != Some(samples) || state >= samples
        || (state == 0) != (at == 0.0) || (state == steps) != (at == horizon)
        || value < final_k || (state == steps && value != final_k)
    {
        return Err("sampled peak does not cover this run's full fine trajectory and selected region");
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    const RECEIPT: &str = r#"{"run":"child-1","transient":{
        "status":"completed","final_time_s":2,"fine_steps":8,
        "sampled_peak":{"region":"solid","value_k":330,"time_s":0,
            "state":0,"samples":9,"includes_initial":true,
            "time_scope":"accepted-fine-grid","continuous_time_bound":false,
            "authority":"Estimated"}}}"#;

    fn read(source: &str, final_k: f64) -> Result<f64, &'static str> {
        temperature(&J::parse(source).unwrap(), "child-1", "solid", final_k)
    }

    #[test]
    fn g1_initial_overheating_is_not_lost_when_the_endpoint_has_cooled() {
        assert_eq!(read(RECEIPT, 300.0), Ok(330.0));
        let interior = RECEIPT.replace("\"time_s\":0", "\"time_s\":0.75")
            .replace("\"state\":0", "\"state\":3");
        assert_eq!(read(&interior, 300.0), Ok(330.0));
        let final_peak = RECEIPT.replace("\"time_s\":0", "\"time_s\":2")
            .replace("\"state\":0", "\"state\":8");
        assert_eq!(read(&final_peak, 330.0), Ok(330.0));
        assert!(read(&final_peak, 300.0).is_err());
    }

    #[test]
    fn g3_peak_projection_refuses_changed_runs_regions_and_incomplete_coverage() {
        for (from, to) in [
            ("child-1", "child-2"), ("solid", "other"),
            ("completed", "paused"), ("\"samples\":9", "\"samples\":8"),
            ("\"samples\":9", "\"samples\":9.0"),
            ("\"fine_steps\":8", "\"fine_steps\":0"),
            ("\"fine_steps\":8", "\"fine_steps\":18446744073709551615"),
            ("\"state\":0", "\"state\":9"),
            ("\"includes_initial\":true", "\"includes_initial\":false"),
            ("accepted-fine-grid", "coarse-grid"), ("Estimated", "Verified"),
            ("\"continuous_time_bound\":false", "\"continuous_time_bound\":true"),
            ("\"value_k\":330", "\"value_k\":299"),
        ] {
            assert!(read(&RECEIPT.replace(from, to), 300.0).is_err(), "accepted {to}");
        }
        assert!(read(r#"{"run":"child-1"}"#, 300.0).is_err());
        assert!(read(&RECEIPT.replace("sampled_peak", "old_peak"), 300.0).is_err());
    }

    #[test]
    fn g0_peak_state_and_time_must_identify_the_same_accepted_endpoint() {
        for (from, to) in [
            ("\"time_s\":0", "\"time_s\":-0.1"),
            ("\"time_s\":0", "\"time_s\":0.25"),
            ("\"state\":0", "\"state\":1"),
            ("\"time_s\":0", "\"time_s\":2.1"),
            ("\"final_time_s\":2", "\"final_time_s\":0"),
        ] {
            assert!(read(&RECEIPT.replace(from, to), 300.0).is_err());
        }
        let wrong_end = RECEIPT.replace("\"time_s\":0", "\"time_s\":2")
            .replace("\"state\":0", "\"state\":7");
        assert!(read(&wrong_end, 300.0).is_err());
        for value in [f64::NAN, f64::INFINITY, -1.0, 331.0] {
            assert!(read(RECEIPT, value).is_err());
        }
    }
}
