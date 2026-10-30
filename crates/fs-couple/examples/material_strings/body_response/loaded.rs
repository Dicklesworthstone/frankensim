//! The actual performer's continuous released response, through its original
//! harmonic owner. This adapter adds only a bounded report, not a physical law.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::MultiBowStringBank;

impl Options {
    /// One complete sweep shares the performance's explicitly declared work
    /// allowance. Every sample is validated before a new report is created.
    pub(crate) fn write_loaded(
        &self,
        bank: &MultiBowStringBank,
        maximum_work: u64,
        gate: &fs_exec::CancelGate,
    ) -> Result<()> {
        if !self.loaded {
            return Err("loaded export requires --loaded-body-response".into());
        }
        let source = bank.source().bank();
        let ports = if let Some(image) = source.material_shell_image() {
            image.force_ports()
        } else if let Some(image) = source.material_shell_assembly_image() {
            image.force_ports()
        } else {
            return Err("loaded body-response requires a curved or assembled shell".into());
        };
        if self.bridge >= ports.len() {
            return Err("body-response bridge index is absent from this body".into());
        }
        let per_frequency = usize::try_from(maximum_work / self.points as u64)
            .map_err(|_| "loaded response work allowance exceeds this host's address size")?;
        if per_frequency == 0 {
            return Err("loaded response work allowance is empty".into());
        }
        let mut text = String::new();
        writeln!(text, "# frankensim-loaded-body-response-v1")?;
        writeln!(
            text,
            "# scope: continuous linear response of the actual coupled performer with bow/contact traction released"
        )?;
        writeln!(
            text,
            "# peak phasors: exp(-i omega t); force: 1 N at physical body bridge {}",
            self.bridge
        )?;
        writeln!(
            text,
            "# loading: retained strings, bridge connections, selected cavity, material memory and admitted exterior reaction"
        )?;
        writeln!(
            text,
            "# receiver: original coupled time-owner pressure; no additional microphone filter; no finite-amplitude strain certificate"
        )?;
        writeln!(
            text,
            "# material law: actual retained constitutive realization; source bands enforced; no state or history is changed"
        )?;
        writeln!(
            text,
            "# total work allowance: {maximum_work}; per-frequency allowance: {per_frequency}"
        )?;
        writeln!(
            text,
            "# exterior gas density_kg_m3={} sound_speed_m_s={}",
            source.gas().density,
            source.gas().sound_speed
        )?;
        if let Some(memory) = source.material_shell_coating_relaxation() {
            let band = memory.source_band_rad_s();
            writeln!(
                text,
                "# coating source-fit band Hz: {} {}",
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
        write!(text, "frequency_hz")?;
        for port in 0..ports.len() {
            write!(
                text,
                "\tbridge{port}_mobility_re_m_s_per_N\tbridge{port}_mobility_im_m_s_per_N"
            )?;
        }
        writeln!(
            text,
            "\tshell_exterior_volume_re_m3_s_per_N\tshell_exterior_volume_im_m3_s_per_N\tshell_cavity_volume_re_m3_s_per_N\tshell_cavity_volume_im_m3_s_per_N\treceiver_pressure_re_Pa_per_N\treceiver_pressure_im_Pa_per_N\texternal_power_W_per_N2\tcomponent_loss_W_per_N2\tconnection_loss_W_per_N2\tcoating_loss_W_per_N2\tother_material_loss_W_per_N2\tradiation_loss_W_per_N2\ttotal_loss_W_per_N2\tmean_stored_energy_J_per_N2\tcoating_mean_stored_energy_J_per_N2\tpower_residual_W_per_N2\tsolve_relative_residual"
        )?;
        if text.len() > 8 * 1_048_576 {
            return Err("loaded response source metadata exceeds 8 MiB".into());
        }
        for point in 0..self.points {
            poll(gate)?;
            let frequency = self.frequency_hz(point);
            let value = bank.material_shell_bridge_harmonic_response(
                core::f64::consts::TAU * frequency,
                self.bridge,
                C64::from_re(1.0),
                per_frequency,
                gate,
            )?;
            let dynamics = &value.released.dynamics;
            write!(text, "{frequency:.17e}")?;
            for port in &value.bridge_velocity_m_s {
                write!(text, "\t{:.17e}\t{:.17e}", port.re, port.im)?;
            }
            for phasor in [
                value.exterior_volume_velocity_m3_s,
                value.cavity_volume_velocity_m3_s,
                value.released.observer_pressure_pa,
            ] {
                write!(text, "\t{:.17e}\t{:.17e}", phasor.re, phasor.im)?;
            }
            for scalar in [
                dynamics.external_power_w,
                dynamics.component_dissipation_w,
                dynamics.connection_dissipation_w,
                value.coating_dissipation_w,
                value.other_material_dissipation_w,
                dynamics.radiation_dissipation_w,
                dynamics.total_dissipation_w(),
                dynamics.mean_stored_energy_j,
                value.coating_mean_stored_energy_j,
                dynamics.power_residual_w,
                dynamics.solve_relative_residual,
            ] {
                write!(text, "\t{scalar:.17e}")?;
            }
            writeln!(text)?;
            if text.len() > 8 * 1_048_576 {
                return Err("loaded response exceeds 8 MiB; reduce its point count".into());
            }
        }
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&self.output)?
            .write_all(text.as_bytes())?;
        eprintln!(
            "loaded_body_response={:?}; frequencies={}; drive_bridge={}; normalization_peak_force_N=1; scope=linear_released_coupled_response",
            self.output, self.points, self.bridge
        );
        Ok(())
    }
}
