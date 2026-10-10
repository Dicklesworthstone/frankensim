//! Joint temperature/enthalpy discrepancy on the actual coupled endpoint map.
//! Only the accepted pair carries h; the coarse result and rejected pairs never
//! contribute physical history, phase summaries or energy accounting.
use super::*;
use super::super::adaptive as time;

#[derive(Debug, Clone, Copy)]
pub(super) struct Config {
    pub(super) time: time::Config,
    pub(super) absolute_specific_enthalpy_tolerance_j_kg: f64,
}

impl Config {
    pub(super) fn parse(value: &J, max_step_s: f64) -> Result<Self> {
        Ok(Self {
            time: time::Config::parse_for_storage(value, max_step_s, true)?,
            absolute_specific_enthalpy_tolerance_j_kg: positive(
                get(value, "absolute_specific_enthalpy_tolerance_j_kg")?,
                "adaptive.absolute_specific_enthalpy_tolerance_j_kg",
            )?,
        })
    }
}

pub(super) type Stats = time::Stats;

#[allow(clippy::too_many_arguments)]
pub(super) fn step(
    cx: &Cx<'_>,
    old_h: &[f64],
    old_temperature: &[f64],
    current_time: f64,
    interval_end: f64,
    suggested_s: f64,
    max_step_s: f64,
    config: Config,
    stats: &mut Stats,
    advance: impl FnMut(&[f64], f64) -> Result<Endpoint>,
) -> Result<time::Accepted<Endpoint>> {
    time::step_with(
        cx, old_h, current_time, interval_end, suggested_s, max_step_s,
        config.time, stats, advance,
        |endpoint: &Endpoint| endpoint.solid.specific_enthalpy_j_kg.as_slice(),
        |old, coarse: &Endpoint, fine: &Endpoint| discrepancy(
            cx, old, old_temperature,
            &coarse.solid.specific_enthalpy_j_kg, &coarse.solid.temperature,
            &fine.solid.specific_enthalpy_j_kg, &fine.solid.temperature, config,
        ),
    )
}

#[allow(clippy::too_many_arguments)]
fn discrepancy(
    cx: &Cx<'_>, old_h: &[f64], old_temperature: &[f64],
    coarse_h: &[f64], coarse_temperature: &[f64],
    fine_h: &[f64], fine_temperature: &[f64], config: Config,
) -> Result<f64> {
    if old_h.len() != old_temperature.len() {
        return Err(bad("adaptive enthalpy requires paired nodal h and temperature fields"));
    }
    let temperature_ratio = time::field_error_ratio(
        cx, old_temperature, coarse_temperature, fine_temperature,
        config.time.absolute_tolerance_k, config.time.relative_tolerance,
    )?;
    let enthalpy_ratio = time::field_error_ratio(
        cx, old_h, coarse_h, fine_h,
        config.absolute_specific_enthalpy_tolerance_j_kg,
        config.time.relative_tolerance,
    )?;
    Ok(temperature_ratio.max(enthalpy_ratio))
}

pub(super) fn render(stats: &Stats, config: Option<Config>) -> Result<String> {
    let Some(config) = config else { return Ok("null".into()); };
    Ok(format!(
        "{{\"method\":\"backward-euler-enthalpy-step-doubling\",\"trials\":{},\"rejected_trials\":{},\"accepted_trials\":{},\"largest_accepted_error_ratio\":{},\"absolute_tolerance_k\":{},\"absolute_specific_enthalpy_tolerance_j_kg\":{},\"relative_temperature_change_tolerance\":{},\"relative_specific_enthalpy_change_tolerance\":{},\"minimum_trial_step_s\":{},\"max_trials\":{},\"scope\":\"maximum normalized nodal discrepancy in both temperature and specific enthalpy; relative scales use changes from the old state, not reference offsets; accepts two complete coupled half-steps without extrapolation; only accepted h is history and only accepted heat enters energy accounting; all attempted solid work is charged; no global, midpoint, inter-step-peak or continuum error bound; adaptive-grid derivatives are not supplied\"}}",
        stats.trials, stats.rejected, stats.trials - stats.rejected,
        num(stats.largest_accepted_ratio)?,
        num(config.time.absolute_tolerance_k)?,
        num(config.absolute_specific_enthalpy_tolerance_j_kg)?,
        num(config.time.relative_tolerance)?,
        num(config.time.relative_tolerance)?,
        num(config.time.minimum_trial_step_s)?,
        config.time.max_trials,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::super::super::tests::with_cx;

    const POLICY: &str = r#"{"absolute_tolerance_k":0.1,"absolute_specific_enthalpy_tolerance_j_kg":0.25,"relative_tolerance":0.01,"minimum_trial_step_s":0.00001,"max_trials":100}"#;

    #[test]
    fn adaptive_enthalpy_requires_its_explicit_positive_tolerance() {
        let value = J::parse(POLICY).unwrap();
        assert!(Config::parse(&value, 1.0).is_ok());
        assert!(time::Config::parse(&value, 1.0).is_err());
        let missing = POLICY.replace("\"absolute_specific_enthalpy_tolerance_j_kg\":0.25,", "");
        assert!(Config::parse(&J::parse(&missing).unwrap(), 1.0).is_err());
        assert!(time::Config::parse(&J::parse(&missing).unwrap(), 1.0).is_ok());
        for replacement in ["0", "-1"] {
            let value = POLICY.replace(
                "\"absolute_specific_enthalpy_tolerance_j_kg\":0.25",
                &format!("\"absolute_specific_enthalpy_tolerance_j_kg\":{replacement}"),
            );
            assert!(Config::parse(&J::parse(&value).unwrap(), 1.0).is_err());
        }
    }

    #[test]
    fn latent_state_error_and_reference_offsets_use_the_full_field() {
        with_cx(|cx| {
            let config = Config::parse(&J::parse(POLICY).unwrap(), 1.0).unwrap();
            let t = [350.0, 350.0];
            let old = [2000.0, 2050.0];
            let coarse = [2020.0, 2070.0];
            let fine = [2021.0, 2072.0];
            let ratio = discrepancy(cx, &old, &t, &coarse, &t, &fine, &t, config).unwrap();
            assert!(ratio > 1.0, "flat temperature must not hide h discrepancy");
            let shift = |field: [f64; 2]| field.map(|h| h + 1_000_000.0);
            let shifted = discrepancy(
                cx, &shift(old), &t, &shift(coarse), &t, &shift(fine), &t, config,
            ).unwrap();
            assert_eq!(ratio, shifted, "enthalpy reference zero must not affect refinement");
            let mut relaxed = config;
            relaxed.absolute_specific_enthalpy_tolerance_j_kg = 10.0;
            assert!(discrepancy(cx, &old, &t, &coarse, &t, &fine, &t, relaxed).unwrap() < 1.0);
            assert!(discrepancy(
                cx, &old, &t, &coarse, &t, &fine, &[350.0, 351.0], relaxed,
            ).unwrap() > 1.0, "temperature discrepancy must also pass");
            assert!(discrepancy(cx, &old, &t[..1], &coarse, &t, &fine, &t, config).is_err());
        });
    }
}
