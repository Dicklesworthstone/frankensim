//! Explicit causal observation policy; the original physical clock is data.

use super::{BandlimitedObservationPolicy, ResolvedDgInputError, Rows, bad};

pub(super) fn read(
    rows: &mut Rows<'_>,
) -> Result<Option<BandlimitedObservationPolicy>, ResolvedDgInputError> {
    if !rows.next_is("observation_filter") {
        return Ok(None);
    }
    let mut row = rows.row("observation_filter")?;
    if row.word()? != "bandlimited" {
        return Err(bad(row.line, "observation_filter requires bandlimited"));
    }
    let policy = BandlimitedObservationPolicy {
        passband_hz: row.scalar()?,
        stopband_hz: row.scalar()?,
        phase_count: row.number()?,
        max_preparation_terms: row.number()?,
        max_history_values: row.number()?,
        max_output_terms: row.number()?,
    };
    row.end()?;
    // The original physical dt, both Nyquist limits and the complete channel
    // layout are admitted together by the shared cold resampler constructor.
    Ok(Some(policy))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn optional_policy_retains_exact_authored_parameters_and_record_order() {
        let text = "observation_filter bandlimited 6000 8000 1024 2000000 65536 65536\ngas dry-air 300 100000\n";
        let mut rows = Rows {
            lines: text.lines(),
            line: 3,
        };
        let policy = read(&mut rows).unwrap().unwrap();
        assert_eq!(policy.passband_hz, 6000.0);
        assert_eq!(policy.stopband_hz, 8000.0);
        assert_eq!(policy.phase_count, 1024);
        assert_eq!(policy.max_preparation_terms, 2_000_000);
        assert_eq!(policy.max_history_values, 65_536);
        assert_eq!(policy.max_output_terms, 65_536);
        assert!(rows.next_is("gas"));
        assert!(read(&mut rows).unwrap().is_none());
        assert_eq!(rows.line, 4);
        assert!(rows.next_is("gas"));
    }

    #[test]
    fn malformed_policy_never_falls_back_to_unfiltered_observation() {
        for fields in [
            "unknown 6000 8000 1024 2000000 65536 65536",
            "bandlimited 6000 8000 1024 2000000 65536",
            "bandlimited 6000 8000 1024 2000000 65536 65536 trailing",
            "bandlimited NaN 8000 1024 2000000 65536 65536",
            "bandlimited 6000 inf 1024 2000000 65536 65536",
            "bandlimited 6000 8000 1.5 2000000 65536 65536",
            "bandlimited 6000 8000 1024 -1 65536 65536",
        ] {
            let text = format!("observation_filter {fields}\n");
            let mut rows = Rows {
                lines: text.lines(),
                line: 3,
            };
            assert!(read(&mut rows).is_err(), "{fields}");
        }
    }
}
