//! Read-only structural response of the same source-bound shell image.
//! This linear transfer is normalized to a unit peak bridge force; WAV rendering
//! continues through the complete string/body/cavity/bow owner.
use super::{Result, input, poll};
use fs_couple::shear_shell::coating::harmonic::{
    ShellCoatingHarmonicOptions, ShellCoatingHarmonicResponse,
};
use fs_couple::shear_shell::{MaterialShellImage, assembly::MaterialShellAssemblyImage};
use fs_math::c64::C64;
use std::{ffi::OsString, fmt::Write as _, fs::OpenOptions, io::Write as _, path::PathBuf};

#[path = "body_response/loaded.rs"]
mod loaded;

pub(super) struct Options {
    analysis_only: bool,
    loaded: bool,
    output: PathBuf,
    minimum_hz: f64,
    maximum_hz: f64,
    points: usize,
    /// None requests every unit-force column of the loaded bridge matrix.
    bridge: Option<usize>,
}

impl Options {
    pub(super) fn extract(args: &[OsString]) -> Result<(Vec<OsString>, Option<Self>)> {
        let mut retained = args.iter().take(2).cloned().collect::<Vec<_>>();
        let mut selected = None;
        let mut index = 2;
        while index < args.len() {
            if args[index] != "--body-response"
                && args[index] != "--body-response-only"
                && args[index] != "--loaded-body-response"
                && args[index] != "--loaded-body-response-only"
            {
                retained.push(args[index].clone());
                index += 1;
                continue;
            }
            if selected.is_some() {
                return Err("select exactly one structural or loaded body-response option".into());
            }
            let values = args
                .get(index + 1..index + 6)
                .ok_or("body-response needs NEW.tsv MIN_HZ MAX_HZ POINTS BRIDGE_INDEX")?;
            let number = |value: &OsString| -> Result<f64> {
                Ok(value
                    .to_str()
                    .ok_or("body-response numbers must be UTF-8")?
                    .parse()?)
            };
            let minimum_hz = number(&values[1])?;
            let maximum_hz = number(&values[2])?;
            let points: usize = values[3]
                .to_str()
                .ok_or("invalid body-response point count")?
                .parse()?;
            let loaded = args[index] == "--loaded-body-response"
                || args[index] == "--loaded-body-response-only";
            let bridge = if values[4] == "all" {
                if !loaded {
                    return Err("bridge selection 'all' requires a loaded body-response command".into());
                }
                None
            } else {
                Some(values[4].to_str().ok_or("invalid body-response bridge index")?.parse::<usize>()?)
            };
            if !minimum_hz.is_finite()
                || !maximum_hz.is_finite()
                || minimum_hz <= 0.0
                || maximum_hz <= minimum_hz
                || !(2..=1024).contains(&points)
                || bridge.is_some_and(|bridge| bridge >= 8)
            {
                return Err("body-response requires finite 0 < MIN_HZ < MAX_HZ, 2..=1024 points and a bridge index in 0..8 (or 'all' for a loaded response)".into());
            }
            selected = Some(Self {
                analysis_only: args[index] == "--body-response-only"
                    || args[index] == "--loaded-body-response-only",
                loaded,
                output: PathBuf::from(&values[0]),
                minimum_hz,
                maximum_hz,
                points,
                bridge,
            });
            index += 6;
        }
        Ok((retained, selected))
    }

    pub(super) const fn analysis_only(&self) -> bool {
        self.analysis_only
    }
    pub(super) const fn loaded(&self) -> bool {
        self.loaded
    }

    fn frequency_hz(&self, point: usize) -> f64 {
        if point == 0 {
            self.minimum_hz
        } else if point + 1 == self.points {
            self.maximum_hz
        } else {
            let lo = fs_math::det::ln(self.minimum_hz);
            let span = fs_math::det::ln(self.maximum_hz) - lo;
            fs_math::det::exp(lo + span * point as f64 / (self.points - 1) as f64)
        }
    }

    /// Finish every source-band and solve check before creating a new report.
    pub(super) fn write(
        &self,
        shell: Option<&MaterialShellImage>,
        assembly: Option<&MaterialShellAssemblyImage>,
        gate: &fs_exec::CancelGate,
    ) -> Result<()> {
        if self.loaded {
            return Err("loaded body-response requires the actual coupled performer".into());
        }
        let bridge = self.bridge.ok_or("the unloaded structural response requires one bridge")?;
        let (ports, memory) = if let Some(image) = shell {
            (image.force_ports(), image.coating_relaxation())
        } else if let Some(image) = assembly {
            (image.force_ports(), image.coating_relaxation())
        } else {
            return Err("body-response requires a curved shell or assembled shell body".into());
        };
        let force: Vec<_> = ports
            .get(bridge)
            .ok_or("body-response bridge index is absent from this body")?
            .iter()
            .map(|&value| C64::from_re(value))
            .collect();
        let films = memory.map_or(0, |image| image.films().len());
        let columns = 2 * ports.len() + films + 12;
        let bytes = self
            .points
            .checked_mul(columns)
            .and_then(|n| n.checked_mul(32))
            .and_then(|n| n.checked_add(films.checked_mul(256)?));
        if bytes.is_none_or(|bytes| bytes > 8 * 1_048_576) {
            return Err("body-response exceeds 8 MiB; reduce its point count".into());
        }
        let budget = ShellCoatingHarmonicOptions {
            maximum_work: input::MAX_WORK as usize / self.points,
            solve_relative_tolerance: 1e-10,
        };
        let response = |omega| -> Result<ShellCoatingHarmonicResponse> {
            Ok(if let Some(image) = shell {
                image.coating_harmonic_response(omega, &force, budget, gate)?
            } else {
                assembly
                    .expect("checked shell image")
                    .coating_harmonic_response(omega, &force, budget, gate)?
            })
        };
        let mut text = String::new();
        writeln!(
            text,
            "# frankensim-body-response-v1\n# scope: unloaded structural body\n# excluded: strings, bow/contact, cavity pressure feedback, exterior radiation reaction and configured microphone filters"
        )?;
        writeln!(
            text,
            "# peak phasors: exp(-i omega t); linear transfer normalized to a 1 N peak force at bridge {}\n# amplitude range: no per-element strain certificate\n# film law: exact source response at each admitted frequency; wood/glue: fixed reduced viscosity",
            bridge
        )?;
        writeln!(
            text,
            "# modes: {}; total work allowance: {}; source support checked at every frequency",
            force.len(),
            input::MAX_WORK
        )?;
        if let Some(memory) = memory {
            let band = memory.source_band_rad_s();
            writeln!(
                text,
                "# source-fit band Hz: {} {}",
                band.0 / core::f64::consts::TAU,
                band.1 / core::f64::consts::TAU
            )?;
            for (index, film) in memory.films().iter().enumerate() {
                let input = film.source().input();
                writeln!(
                    text,
                    "# film {index}: region={} source={} thickness_m={}",
                    input.region,
                    input.material.source().identity(),
                    input.thickness_m
                )?;
            }
        }
        if text.len() > 8 * 1_048_576 {
            return Err("body-response source metadata exceeds 8 MiB".into());
        }
        write!(text, "frequency_hz")?;
        for port in 0..ports.len() {
            write!(
                text,
                "\tbridge{port}_mobility_re_m_s_per_N\tbridge{port}_mobility_im_m_s_per_N"
            )?;
        }
        write!(
            text,
            "\texterior_volume_re_m3_s_per_N\texterior_volume_im_m3_s_per_N\tcavity_volume_re_m3_s_per_N\tcavity_volume_im_m3_s_per_N\texternal_power_W_per_N2\twood_glue_loss_W_per_N2\tinterface_loss_W_per_N2"
        )?;
        for film in 0..films {
            write!(text, "\tfilm{film}_loss_W_per_N2")?;
        }
        writeln!(
            text,
            "\ttotal_loss_W_per_N2\tpower_residual_W_per_N2\tsolve_relative_residual"
        )?;
        for point in 0..self.points {
            poll(gate)?;
            let frequency = self.frequency_hz(point);
            let value = response(core::f64::consts::TAU * frequency)?;
            write!(text, "{frequency:.17e}")?;
            for port in &value.force_port_velocity_m_s {
                write!(text, "\t{:.17e}\t{:.17e}", port.re, port.im)?;
            }
            for volume in [
                value.exterior_volume_velocity_m3_s,
                value.cavity_volume_velocity_m3_s,
            ] {
                write!(text, "\t{:.17e}\t{:.17e}", volume.re, volume.im)?;
            }
            write!(
                text,
                "\t{:.17e}\t{:.17e}\t{:.17e}",
                value.external_power_w,
                value.background_dissipation_w,
                value.interface_dissipation_w
            )?;
            for loss in &value.film_dissipation_w {
                write!(text, "\t{loss:.17e}")?;
            }
            writeln!(
                text,
                "\t{:.17e}\t{:.17e}\t{:.17e}",
                value.total_dissipation_w(),
                value.power_residual_w,
                value.solve_relative_residual
            )?;
            if text.len() > 8 * 1_048_576 {
                return Err("body-response exceeds 8 MiB; reduce its point count".into());
            }
        }
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&self.output)?
            .write_all(text.as_bytes())?;
        eprintln!(
            "body_response={:?}; frequencies={}; drive_bridge={}; normalization_peak_force_N=1; scope=linear_unloaded_structural_response",
            self.output, self.points, bridge
        );
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn words(values: &[&str]) -> Vec<OsString> {
        values.iter().map(OsString::from).collect()
    }

    #[test]
    fn response_only_needs_no_wav_and_preserves_source_options() {
        let args = words(&[
            "renderer",
            "performance",
            "--body-response-only",
            "response.tsv",
            "100",
            "2000",
            "32",
            "0",
            "--body-profile",
            "body.profile",
            "--self-contained",
        ]);
        let (retained, option) = Options::extract(&args).unwrap();
        assert!(option.unwrap().analysis_only());
        assert_eq!(
            retained,
            words(&[
                "renderer",
                "performance",
                "--body-profile",
                "body.profile",
                "--self-contained"
            ])
        );
        let args = words(&[
            "renderer",
            "performance",
            "render.wav",
            "--body-response",
            "response.tsv",
            "100",
            "2000",
            "32",
            "0",
        ]);
        let (retained, option) = Options::extract(&args).unwrap();
        assert!(!option.unwrap().analysis_only());
        assert_eq!(retained, words(&["renderer", "performance", "render.wav"]));
        let args = words(&[
            "renderer",
            "performance",
            "render.wav",
            "--loaded-body-response",
            "loaded.tsv",
            "100",
            "2000",
            "32",
            "0",
        ]);
        let (retained, option) = Options::extract(&args).unwrap();
        assert!(option.as_ref().unwrap().loaded());
        assert!(!option.unwrap().analysis_only());
        assert_eq!(retained, words(&["renderer", "performance", "render.wav"]));
        let args = words(&[
            "renderer", "performance", "--loaded-body-response-only",
            "loaded.tsv", "101", "2000", "32", "0",
            "--body-profile", "body.profile", "--self-contained",
        ]);
        let (retained, option) = Options::extract(&args).unwrap();
        let option = option.unwrap();
        assert!(option.loaded() && option.analysis_only());
        assert_eq!(retained, words(&[
            "renderer", "performance", "--body-profile", "body.profile", "--self-contained",
        ]));
    }

    #[test]
    fn response_modes_are_exclusive_and_require_complete_arguments() {
        let first = [
            "renderer",
            "performance",
            "--body-response-only",
            "response.tsv",
            "100",
            "2000",
            "32",
            "0",
        ];
        let mut repeated = words(&first);
        repeated.extend(words(&[
            "--body-response",
            "other.tsv",
            "100",
            "2000",
            "32",
            "0",
        ]));
        assert!(Options::extract(&repeated).is_err());
        repeated.truncate(first.len());
        repeated.extend(words(&[
            "--loaded-body-response-only", "other.tsv", "101", "2000", "32", "0",
        ]));
        assert!(Options::extract(&repeated).is_err());
        assert!(Options::extract(&words(&first[..7])).is_err());
    }

    #[test]
    fn all_bridges_selects_a_loaded_matrix_without_changing_source_arguments() {
        for command in ["--loaded-body-response", "--loaded-body-response-only"] {
            let args = words(&[
                "renderer", "performance", command, "matrix.tsv", "101", "2000", "32", "all",
                "--body-profile", "body.profile",
            ]);
            let (retained, option) = Options::extract(&args).unwrap();
            let option = option.unwrap();
            assert!(option.loaded());
            assert_eq!(option.bridge, None);
            assert_eq!(retained, words(&["renderer", "performance", "--body-profile", "body.profile"]));
        }
        for command in ["--body-response", "--body-response-only"] {
            assert!(Options::extract(&words(&[
                "renderer", "performance", command, "matrix.tsv", "101", "2000", "32", "all",
            ])).is_err());
        }
        for bridge in ["-1", "8", "ALL", "0,1", ""] {
            assert!(Options::extract(&words(&[
                "renderer", "performance", "--loaded-body-response-only", "matrix.tsv",
                "101", "2000", "32", bridge,
            ])).is_err());
        }
    }
}
