//! Opt-in local coating checks on the original mechanical playback clock.
use super::Result;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::MultiBowPerformance;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::coating::ShellCoatingMonitorOptions;
use std::{ffi::OsString, io::Write};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct Options {
    pub maximum_bytes: usize,
    pub maximum_work: usize,
}

impl Options {
    /// Remove this variable-length option before the existing path-pair parser.
    /// The same syntax works alongside a harmonic report and a WAV render.
    pub(super) fn extract(args: &[OsString]) -> Result<(Vec<OsString>, Option<Self>)> {
        let mut retained = args.iter().take(2).cloned().collect::<Vec<_>>();
        let mut selected = None;
        let mut index = 2;
        while index < args.len() {
            if args[index] != "--coating-strain-guard" {
                retained.push(args[index].clone());
                index += 1;
                continue;
            }
            if selected.is_some() {
                return Err("select --coating-strain-guard only once".into());
            }
            let values = args
                .get(index + 1..index + 3)
                .ok_or("--coating-strain-guard needs MAX_BYTES MAX_WORK")?;
            let positive = |value: &OsString| -> Result<usize> {
                let value: usize = value
                    .to_str()
                    .ok_or("coating strain budgets must be UTF-8 integers")?
                    .parse()?;
                if value == 0 {
                    return Err(
                        "coating strain requires positive explicit byte and work allowances".into(),
                    );
                }
                Ok(value)
            };
            selected = Some(Self {
                maximum_bytes: positive(&values[0])?,
                maximum_work: positive(&values[1])?,
            });
            index += 3;
        }
        Ok((retained, selected))
    }

    pub(super) fn validate_context(self, analysis_only: bool, shell: bool) -> Result<()> {
        if analysis_only {
            return Err("--coating-strain-guard requires a WAV performance; harmonic response uses --body-response-strain".into());
        }
        if !shell {
            return Err(
                "--coating-strain-guard requires a curved or assembled shell body profile".into(),
            );
        }
        Ok(())
    }

    pub(super) fn budget(self) -> ShellCoatingMonitorOptions {
        ShellCoatingMonitorOptions {
            maximum_bytes: self.maximum_bytes,
            maximum_work: self.maximum_work,
        }
    }
}

/// Borrow the monitor's aggregate; neither the report nor I/O advances physics.
/// Also called after a failed render so an envelope refusal retains its location.
pub(super) fn describe(performance: &MultiBowPerformance, output: &mut impl Write) -> Result<()> {
    let Some(report) = performance.coating_monitor_report() else {
        return Ok(());
    };
    writeln!(
        output,
        "coating_strain_guard={:?}; resources={:?}; initial_absolute_sample={}; checked_endpoints={}; last_checked_sample={}; failure={:?}",
        report.options,
        report.resources,
        report.initial_absolute_sample,
        report.checked_endpoints,
        report.last_checked_sample,
        report.failure
    )?;
    writeln!(
        output,
        "coating_strain_scope=actual retained body displacement at initial and accepted mechanical endpoints; full coated facets and film thickness; deterministic floating-point Bernstein bounds; component_order=[xx,yy,zz,xy,yz,zx]; authored constitutive envelope, not a material failure threshold; no between-step or omitted prestress/cure/continuum certificate"
    )?;
    if let Some(memory) = report.relaxation {
        writeln!(
            output,
            "coating_memory_branches={}; coating_memory_stored_energy_j={:.17e}; coating_memory_instantaneous_dissipation_w={:.17e}; maximum_coating_memory_stored_energy_j={:.17e}; maximum_coating_memory_stored_energy_sample={}; maximum_coating_memory_instantaneous_dissipation_w={:.17e}; maximum_coating_memory_instantaneous_dissipation_sample={}",
            memory.current.branch_count,
            memory.current.stored_energy_j,
            memory.current.dissipated_power_w,
            memory.maximum_stored_energy_j,
            memory.maximum_stored_energy_sample,
            memory.maximum_dissipated_power_w,
            memory.maximum_dissipated_power_sample
        )?;
        writeln!(
            output,
            "coating_memory_scope=original accepted material histories at the same checked endpoints; relaxing-arm storage only, equilibrium film stiffness remains in shell energy; instantaneous W is not integrated heat; equal-time factors may mix films; no thermal feedback or temperature rise is inferred"
        )?;
    } else {
        writeln!(
            output,
            "coating_memory_scope=no hereditary coating image selected; no coating-memory energy or power result"
        )?;
    }
    for (index, film) in report.films().enumerate() {
        let original = film.source.source;
        let input = original.input();
        let bound = film.maximum;
        writeln!(
            output,
            "coating_strain_film={index}; component_index={:?}; component_name={:?}; source_film={}; region={:?}; source={}; faces_m={:?}; authored_frobenius_limit={:.17e}; covered_elements={}; maximum_frobenius_bound={:.17e}; maximum_limit_ratio={:.17e}; maximum_sample={}; component_abs_upper_bounds={:?}; ruling_element={}; ruling_control={}; ruling_z_face={}",
            film.source.component_index,
            film.source.component_name,
            film.source.film_index,
            input.region,
            input.material.source().identity(),
            original.faces_m(),
            input.material.strain_limit(),
            bound.covered_elements,
            bound.frobenius_upper_bound,
            bound.limit_ratio,
            bound.sample,
            bound.component_abs_upper_bounds,
            bound.ruling_element,
            bound.ruling_control_point,
            bound.ruling_z_face
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(values: &[&str]) -> Vec<OsString> {
        values.iter().map(OsString::from).collect()
    }

    #[test]
    fn absent_guard_preserves_the_existing_command() {
        let input = args(&[
            "material_strings",
            "source",
            "new.wav",
            "--body-profile",
            "body",
        ]);
        let (retained, guard) = Options::extract(&input).unwrap();
        assert_eq!(retained, input);
        assert_eq!(guard, None);
    }

    #[test]
    fn guard_keeps_response_trace_and_observation_options_in_original_order() {
        let input = args(&[
            "material_strings",
            "source",
            "new.wav",
            "--body-profile",
            "body",
            "--coating-strain-guard",
            "1048576",
            "1000000000",
            "--bow-trace",
            "new.tsv",
            "--decimate",
        ]);
        let (retained, guard) = Options::extract(&input).unwrap();
        assert_eq!(
            retained,
            args(&[
                "material_strings",
                "source",
                "new.wav",
                "--body-profile",
                "body",
                "--bow-trace",
                "new.tsv",
                "--decimate"
            ])
        );
        let options = guard.unwrap();
        assert_eq!(options.maximum_bytes, 1048576);
        assert_eq!(options.maximum_work, 1000000000);
        options.validate_context(false, true).unwrap();
    }

    #[test]
    fn malformed_missing_or_duplicate_allowances_refuse() {
        for suffix in [
            vec![],
            vec!["10"],
            vec!["0", "1"],
            vec!["1", "0"],
            vec!["-1", "1"],
            vec!["1", "1.5"],
            vec!["1", "--decimate"],
            vec!["1", "2", "--coating-strain-guard", "1", "2"],
        ] {
            let mut input = args(&[
                "material_strings",
                "source",
                "new.wav",
                "--coating-strain-guard",
            ]);
            input.extend(args(&suffix));
            assert!(Options::extract(&input).is_err(), "accepted {input:?}");
        }
    }

    #[test]
    fn a_guard_requires_playback_and_a_physical_shell() {
        let options = Options {
            maximum_bytes: 1,
            maximum_work: 1,
        };
        assert!(options.validate_context(true, true).is_err());
        assert!(options.validate_context(false, false).is_err());
    }
}
