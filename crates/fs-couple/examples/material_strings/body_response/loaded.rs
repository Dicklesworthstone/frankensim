//! The actual performer's continuous released response, through its original
//! harmonic owner. This adapter adds only a bounded report, not a physical law.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::harmonic::ShellBridgeHarmonicResponse;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::MultiBowStringBank;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::observation::harmonic::{
    ShellHarmonicReceiver, ShellHarmonicReceiverSelection,
};
use fs_couple::shear_shell::coating::strain::ShellCoatingStrainReport;

impl Options {
    /// Observe the same microphone declaration selected for playback, using
    /// its original frequency producer instead of a sampled FIR realization.
    pub(crate) fn write_loaded_for_spec(
        &self, bank: &MultiBowStringBank, spec: &input::Spec,
        profile: Option<&super::super::body::Profile>, maximum_work: u64,
        gate: &fs_exec::CancelGate,
    ) -> Result<()> {
        if spec.plate_observer.is_some() {
            return Err("loaded shell response requires a shell or acoustic-skin receiver".into());
        }
        match (&spec.shell_observer, &spec.one_way_observer) {
            (Some(options), None) => self.write_loaded(bank,
                Some(ShellHarmonicReceiverSelection::Loaded {
                    options, quadrature: spec.shell_receiver_quadrature,
                }), maximum_work, gate),
            (None, Some(options)) => {
                let skin = profile.and_then(|body| body.acoustic_skin(spec))
                    .ok_or("harmonic skin receiver requires the original profile's acoustic skin")?;
                self.write_loaded(bank, Some(ShellHarmonicReceiverSelection::Skin {
                    options, skin: &skin,
                }), maximum_work, gate)
            }
            (None, None) => self.write_loaded(bank, None, maximum_work, gate),
            _ => Err("loaded response requires one physical receiver declaration".into()),
        }
    }

    /// One complete sweep shares the performance's explicitly declared work
    /// allowance. Every sample is validated before a new report is created.
    pub(crate) fn write_loaded(
        &self,
        bank: &MultiBowStringBank,
        receiver: Option<ShellHarmonicReceiverSelection<'_>>,
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
        if ports.is_empty() || self.bridge.is_some_and(|bridge| bridge >= ports.len()) {
            return Err("body-response bridge index is absent from this body".into());
        }
        let drives = self.bridge.map_or(ports.len(), |_| 1);
        let evaluations = self.points.checked_mul(drives)
            .ok_or("coating strain evaluation count overflow")?;
        let mut strain = self.strain.map(|options| options.prepare(
            source.material_shell_image(), source.material_shell_assembly_image(), evaluations, gate))
            .transpose()?;
        let strain_columns = strain.as_ref().map(|observer| strain::additional_columns(observer.film_count()))
            .transpose()?.unwrap_or(0);
        let columns = (2 * ports.len() + 23 + if receiver.is_some() { 6 } else { 0 })
            .checked_add(strain_columns).ok_or("loaded response column count overflow")?;
        if self.points.checked_mul(drives).and_then(|n| n.checked_mul(columns))
            .and_then(|n| n.checked_mul(32)).is_none_or(|bytes| bytes > 8 * 1_048_576) {
            return Err("loaded response exceeds 8 MiB; reduce its point count or select one bridge".into());
        }
        let per_frequency = usize::try_from(maximum_work / self.points as u64)
            .map_err(|_| "loaded response work allowance exceeds this host's address size")?;
        if per_frequency == 0 {
            return Err("loaded response work allowance is empty".into());
        }
        let receiver = receiver.map(|selection| {
            ShellHarmonicReceiver::new(bank, selection, self.points, gate)
        }).transpose()?;
        let prepared = receiver.as_ref().map(|receiver| receiver.prepare(gate)).transpose()?;
        let mut text = String::new();
        writeln!(text, "# frankensim-loaded-body-response-v{}", if strain.is_some() { 3 } else { 2 })?;
        writeln!(
            text,
            "# scope: continuous linear response of the actual coupled performer with bow/contact traction released"
        )?;
        writeln!(
            text,
            "# peak phasors: exp(-i omega t); each row applies 1 N only at its drive_bridge; all bridge velocities retain their physical order"
        )?;
        writeln!(
            text,
            "# loading: retained strings, bridge connections, selected cavity, material memory and admitted exterior reaction"
        )?;
        if let Some(receiver) = &prepared {
            writeln!(text, "# receiver: original exterior source at the forcing frequency; raw physical phase, no FIR latency or output normalization")?;
            writeln!(text, "# receiver_binding: {:?}", receiver.report())?;
        } else {
            writeln!(text, "# receiver: original coupled time-owner pressure; no additional exterior receiver")?;
        }
        writeln!(text, "# native_pressure retains the original modal/compact readout; no finite-amplitude strain or exterior convergence certificate")?;
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
        if let Some(observer) = &strain { strain::metadata(&mut text, observer)?; }
        write!(text, "frequency_hz\tdrive_bridge")?;
        for port in 0..ports.len() {
            write!(
                text,
                "\tbridge{port}_mobility_re_m_s_per_N\tbridge{port}_mobility_im_m_s_per_N"
            )?;
        }
        write!(
            text,
            "\tshell_exterior_volume_re_m3_s_per_N\tshell_exterior_volume_im_m3_s_per_N\tshell_cavity_volume_re_m3_s_per_N\tshell_cavity_volume_im_m3_s_per_N\tnative_pressure_re_Pa_per_N\tnative_pressure_im_Pa_per_N\treceiver_pressure_re_Pa_per_N\treceiver_pressure_im_Pa_per_N\texternal_power_W_per_N2\tcomponent_loss_W_per_N2\tconnection_loss_W_per_N2\tcoating_loss_W_per_N2\tother_material_loss_W_per_N2\tradiation_loss_W_per_N2\ttotal_loss_W_per_N2\tmean_stored_energy_J_per_N2\tcoating_mean_stored_energy_J_per_N2\tpower_residual_W_per_N2\tsolve_relative_residual"
        )?;
        if prepared.is_some() {
            for axis in ["x", "y", "z"] {
                write!(text, "\treceiver_gradient_{axis}_re_Pa_per_m_per_N\treceiver_gradient_{axis}_im_Pa_per_m_per_N")?;
            }
        }
        if let Some(observer) = &strain { strain::header(&mut text, observer.film_count())?; }
        writeln!(text)?;
        if text.len() > 8 * 1_048_576 {
            return Err("loaded response source metadata exceeds 8 MiB".into());
        }
        for point in 0..self.points {
            poll(gate)?;
            let frequency = self.frequency_hz(point);
            let omega = core::f64::consts::TAU * frequency;
            if let Some(receiver) = &prepared {
                let values = if let Some(bridge) = self.bridge {
                    vec![receiver.bridge_harmonic_response(omega, bridge, C64::from_re(1.0), per_frequency, gate)?]
                } else {
                    receiver.bridge_harmonic_matrix(omega, per_frequency, gate)?
                };
                for (column, value) in values.iter().enumerate() {
                    let strain_report = if let Some(observer) = &mut strain {
                        Some(observer.evaluate(value.body.body_displacement_m_sqrt_kg()?, gate)?)
                    } else { None };
                    append_row(&mut text, frequency, self.bridge.unwrap_or(column), &value.body,
                        value.receiver_pressure_pa, Some(&value.receiver_gradient_pa_per_m), strain_report.as_ref())?;
                }
            } else {
                let values = if let Some(bridge) = self.bridge {
                    vec![bank.material_shell_bridge_harmonic_response(omega, bridge, C64::from_re(1.0), per_frequency, gate)?]
                } else {
                    bank.material_shell_bridge_harmonic_matrix(omega, per_frequency, gate)?
                };
                for (column, value) in values.iter().enumerate() {
                    let strain_report = if let Some(observer) = &mut strain {
                        Some(observer.evaluate(value.body_displacement_m_sqrt_kg()?, gate)?)
                    } else { None };
                    append_row(&mut text, frequency, self.bridge.unwrap_or(column), value,
                        value.released.observer_pressure_pa, None, strain_report.as_ref())?;
                }
            }
        }
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&self.output)?
            .write_all(text.as_bytes())?;
        eprintln!(
            "loaded_body_response={:?}; frequencies={}; drive_columns={drives}; exterior_receiver={}; normalization_peak_force_N=1; scope=linear_released_coupled_response",
            self.output, self.points, prepared.is_some()
        );
        Ok(())
    }
}

fn append_row(text: &mut String, frequency: f64, bridge: usize,
    value: &ShellBridgeHarmonicResponse, receiver_pressure: C64, gradient: Option<&[C64; 3]>,
    strain_report: Option<&ShellCoatingStrainReport<'_>>) -> Result<()> {
    let dynamics = &value.released.dynamics;
    write!(text, "{frequency:.17e}\t{bridge}")?;
    for port in &value.bridge_velocity_m_s {
        write!(text, "\t{:.17e}\t{:.17e}", port.re, port.im)?;
    }
    for phasor in [value.exterior_volume_velocity_m3_s, value.cavity_volume_velocity_m3_s,
        value.released.observer_pressure_pa, receiver_pressure] {
        write!(text, "\t{:.17e}\t{:.17e}", phasor.re, phasor.im)?;
    }
    for scalar in [dynamics.external_power_w, dynamics.component_dissipation_w,
        dynamics.connection_dissipation_w, value.coating_dissipation_w,
        value.other_material_dissipation_w, dynamics.radiation_dissipation_w,
        dynamics.total_dissipation_w(), dynamics.mean_stored_energy_j,
        value.coating_mean_stored_energy_j, dynamics.power_residual_w,
        dynamics.solve_relative_residual] {
        write!(text, "\t{scalar:.17e}")?;
    }
    if let Some(gradient) = gradient {
        for component in gradient {
            write!(text, "\t{:.17e}\t{:.17e}", component.re, component.im)?;
        }
    }
    if let Some(report) = strain_report { strain::row(text, report)?; }
    writeln!(text)?;
    if text.len() > 8 * 1_048_576 {
        return Err("loaded response exceeds 8 MiB; reduce its point count".into());
    }
    Ok(())
}
