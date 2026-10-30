//! The actual performer's continuous released response, through its original
//! harmonic owner. This adapter adds only a bounded report, not a physical law.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::harmonic::ShellBridgeHarmonicResponse;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::MultiBowStringBank;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::observation::harmonic::{
    ShellHarmonicReceiver, ShellHarmonicReceiverSelection,
};
use fs_couple::shear_shell::coating::strain::ShellCoatingStrainReport;

/// Column identity for normalized unit drives, or actual SI values for one
/// complete coherent physical load. Numeric response transport is shared.
#[derive(Clone, Copy)]
enum RowKind {
    Unit(usize),
    Physical,
}

fn coherent_forces(ports: usize, authored: &[(usize, C64)]) -> Result<Vec<C64>> {
    if ports == 0 || authored.is_empty() || authored.len() > MAX_COHERENT_FORCE_RECORDS {
        return Err("coherent response requires an original bridge and 1..=64 physical force records".into());
    }
    let mut forces = vec![C64::ZERO; ports];
    for &(bridge, force) in authored {
        let total = forces.get_mut(bridge)
            .ok_or("coherent force names a bridge absent from this body")?;
        if !force.re.is_finite() || !force.im.is_finite() {
            return Err("coherent bridge forces must have finite real and imaginary newtons".into());
        }
        // Authored order is the deterministic summation order, including
        // duplicate actuators. Never normalize or discard a cancelling port.
        total.re += force.re;
        total.im += force.im;
        if !total.re.is_finite() || !total.im.is_finite() {
            return Err("coherent bridge force accumulation overflow".into());
        }
    }
    Ok(forces)
}


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
        let unit_bridge = match &self.excitation {
            Excitation::Unit(bridge) => Some(*bridge),
            _ => None,
        };
        if ports.is_empty() || unit_bridge.is_some_and(|bridge| bridge >= ports.len()) {
            return Err("body-response bridge index is absent from this body".into());
        }
        // Validate every physical input and complete its source-bound vector
        // before preparing receivers, strain maps or an output file.
        let coherent = match &self.excitation {
            Excitation::Coherent(authored) => Some(coherent_forces(ports.len(), authored)?),
            _ => None,
        };
        let units = if coherent.is_some() { strain::Units::Physical } else { strain::Units::PerNewton };
        let drives = if matches!(&self.excitation, Excitation::All) { ports.len() } else { 1 };
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
        writeln!(text, "# frankensim-loaded-body-response-v{}", if coherent.is_some() { 4 } else if strain.is_some() { 3 } else { 2 })?;
        writeln!(
            text,
            "# scope: continuous linear response of the actual coupled performer with bow/contact traction released"
        )?;
        if let Some(forces) = &coherent {
            writeln!(text, "# peak phasors: exp(-i omega t); each row solves the complete coherent physical force vector; actual SI values without force normalization")?;
            writeln!(text, "# duplicate actuators sum in authored order before solving; bridge velocities and summed forces retain original source order")?;
            for (bridge, force) in forces.iter().enumerate() {
                writeln!(text, "# physical_peak_force_N bridge={bridge} real_N={:.17e} imaginary_N={:.17e}", force.re, force.im)?;
            }
        } else {
            writeln!(
                text,
                "# peak phasors: exp(-i omega t); each row applies 1 N only at its drive_bridge; all bridge velocities retain their physical order"
            )?;
        }
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
        if let Some(observer) = &strain { strain::metadata(&mut text, observer, units)?; }
        append_header(&mut text, ports.len(), prepared.is_some(), units)?;
        if let Some(observer) = &strain { strain::header(&mut text, observer.film_count(), units)?; }
        writeln!(text)?;
        if text.len() > 8 * 1_048_576 {
            return Err("loaded response source metadata exceeds 8 MiB".into());
        }
        for point in 0..self.points {
            poll(gate)?;
            let frequency = self.frequency_hz(point);
            let omega = core::f64::consts::TAU * frequency;
            if let Some(receiver) = &prepared {
                let values = if let Some(forces) = &coherent {
                    vec![receiver.bridge_forces_harmonic_response(omega, forces, per_frequency, gate)?]
                } else if let Some(bridge) = unit_bridge {
                    vec![receiver.bridge_harmonic_response(omega, bridge, C64::from_re(1.0), per_frequency, gate)?]
                } else {
                    receiver.bridge_harmonic_matrix(omega, per_frequency, gate)?
                };
                for (column, value) in values.iter().enumerate() {
                    let strain_report = if let Some(observer) = &mut strain {
                        Some(observer.evaluate(value.body.body_displacement_m_sqrt_kg()?, gate)?)
                    } else { None };
                    let kind = if coherent.is_some() { RowKind::Physical }
                        else { RowKind::Unit(unit_bridge.unwrap_or(column)) };
                    append_row(&mut text, frequency, kind, &value.body,
                        value.receiver_pressure_pa, Some(&value.receiver_gradient_pa_per_m), strain_report.as_ref())?;
                }
            } else {
                let values = if let Some(forces) = &coherent {
                    vec![bank.material_shell_bridge_forces_harmonic_response(omega, forces, per_frequency, gate)?]
                } else if let Some(bridge) = unit_bridge {
                    vec![bank.material_shell_bridge_harmonic_response(omega, bridge, C64::from_re(1.0), per_frequency, gate)?]
                } else {
                    bank.material_shell_bridge_harmonic_matrix(omega, per_frequency, gate)?
                };
                for (column, value) in values.iter().enumerate() {
                    let strain_report = if let Some(observer) = &mut strain {
                        Some(observer.evaluate(value.body_displacement_m_sqrt_kg()?, gate)?)
                    } else { None };
                    let kind = if coherent.is_some() { RowKind::Physical }
                        else { RowKind::Unit(unit_bridge.unwrap_or(column)) };
                    append_row(&mut text, frequency, kind, value,
                        value.released.observer_pressure_pa, None, strain_report.as_ref())?;
                }
            }
        }
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&self.output)?
            .write_all(text.as_bytes())?;
        if let Some(forces) = &coherent {
            eprintln!("loaded_body_response={:?}; frequencies={}; coherent_peak_forces_N={forces:?}; exterior_receiver={}; scope=linear_released_coupled_physical_response",
                self.output, self.points, prepared.is_some());
        } else {
            eprintln!(
                "loaded_body_response={:?}; frequencies={}; drive_columns={drives}; exterior_receiver={}; normalization_peak_force_N=1; scope=linear_released_coupled_response",
                self.output, self.points, prepared.is_some()
            );
        }
        Ok(())
    }
}

fn append_header(text: &mut String, ports: usize, receiver: bool, units: strain::Units) -> Result<()> {
    let physical = units == strain::Units::Physical;
    let linear = if physical { "" } else { "_per_N" };
    let quadratic = if physical { "" } else { "_per_N2" };
    write!(text, "frequency_hz")?;
    if !physical { write!(text, "\tdrive_bridge")?; }
    let motion = if physical { "velocity" } else { "mobility" };
    for port in 0..ports {
        write!(text, "\tbridge{port}_{motion}_re_m_s{linear}\tbridge{port}_{motion}_im_m_s{linear}")?;
    }
    for (name, unit) in [("shell_exterior_volume", "m3_s"), ("shell_cavity_volume", "m3_s"),
        ("native_pressure", "Pa"), ("receiver_pressure", "Pa")] {
        write!(text, "\t{name}_re_{unit}{linear}\t{name}_im_{unit}{linear}")?;
    }
    for (name, unit) in [("external_power", "W"), ("component_loss", "W"),
        ("connection_loss", "W"), ("coating_loss", "W"), ("other_material_loss", "W"),
        ("radiation_loss", "W"), ("total_loss", "W"), ("mean_stored_energy", "J"),
        ("coating_mean_stored_energy", "J"), ("power_residual", "W")] {
        write!(text, "\t{name}_{unit}{quadratic}")?;
    }
    write!(text, "\tsolve_relative_residual")?;
    if receiver {
        for axis in ["x", "y", "z"] {
            write!(text, "\treceiver_gradient_{axis}_re_Pa_per_m{linear}\treceiver_gradient_{axis}_im_Pa_per_m{linear}")?;
        }
    }
    Ok(())
}

fn append_row(text: &mut String, frequency: f64, kind: RowKind,
    value: &ShellBridgeHarmonicResponse, receiver_pressure: C64, gradient: Option<&[C64; 3]>,
    strain_report: Option<&ShellCoatingStrainReport<'_>>) -> Result<()> {
    let dynamics = &value.released.dynamics;
    write!(text, "{frequency:.17e}")?;
    if let RowKind::Unit(bridge) = kind { write!(text, "\t{bridge}")?; }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coherent_inputs_sum_complex_actuators_in_authored_source_order() {
        let forces = coherent_forces(3, &[
            (2, C64::new(1.0, 2.0)), (0, C64::new(-3.0, 4.0)),
            (2, C64::new(-0.5, -2.0)),
        ]).unwrap();
        assert_eq!(forces, vec![C64::new(-3.0, 4.0), C64::ZERO, C64::new(0.5, 0.0)]);
        let first = coherent_forces(1, &[
            (0, C64::from_re(1e16)), (0, C64::from_re(-1e16)), (0, C64::from_re(1.0)),
        ]).unwrap();
        let second = coherent_forces(1, &[
            (0, C64::from_re(1e16)), (0, C64::from_re(1.0)), (0, C64::from_re(-1e16)),
        ]).unwrap();
        assert_eq!(first[0], C64::from_re(1.0));
        assert_eq!(second[0], C64::ZERO);
        assert_eq!(coherent_forces(2, &[
            (1, C64::new(2.0, -3.0)), (1, C64::new(-2.0, 3.0)),
        ]).unwrap(), vec![C64::ZERO; 2]);
    }

    #[test]
    fn coherent_inputs_refuse_missing_source_ports_nonfinite_values_and_sum_overflow() {
        assert!(coherent_forces(0, &[(0, C64::ZERO)]).is_err());
        assert!(coherent_forces(2, &[]).is_err());
        assert!(coherent_forces(2, &[(2, C64::ZERO)]).is_err());
        assert!(coherent_forces(2, &vec![(0, C64::ZERO); MAX_COHERENT_FORCE_RECORDS + 1]).is_err());
        for force in [C64::new(f64::NAN, 0.0), C64::new(0.0, f64::INFINITY)] {
            assert!(coherent_forces(2, &[(0, force)]).is_err());
        }
        for force in [C64::new(f64::MAX, 0.0), C64::new(0.0, f64::MAX)] {
            assert!(coherent_forces(2, &[(0, force), (0, force)]).is_err());
        }
    }

    #[test]
    fn physical_and_unit_force_headers_distinguish_actual_motion_power_and_strain_scale() {
        let mut unit = String::new();
        append_header(&mut unit, 1, true, strain::Units::PerNewton).unwrap();
        strain::header(&mut unit, 1, strain::Units::PerNewton).unwrap();
        assert_eq!(unit.split('\t').count(), 44);
        assert!(unit.starts_with("frequency_hz\tdrive_bridge\tbridge0_mobility_re_m_s_per_N"));
        assert!(unit.contains("\treceiver_pressure_re_Pa_per_N\t"));
        assert!(unit.contains("\tcoating_loss_W_per_N2\t"));
        assert!(unit.contains("\treceiver_gradient_x_re_Pa_per_m_per_N\t"));
        assert!(unit.contains("\tcoating_strain_sufficient_peak_force_N\t"));

        let mut physical = String::new();
        append_header(&mut physical, 1, true, strain::Units::Physical).unwrap();
        strain::header(&mut physical, 1, strain::Units::Physical).unwrap();
        assert_eq!(physical.split('\t').count(), 43);
        assert!(physical.starts_with("frequency_hz\tbridge0_velocity_re_m_s\t"));
        assert!(physical.contains("\treceiver_pressure_re_Pa\t"));
        assert!(physical.contains("\tcoating_loss_W\t"));
        assert!(physical.contains("\tcoating_mean_stored_energy_J\t"));
        assert!(physical.contains("\treceiver_gradient_x_re_Pa_per_m\t"));
        assert!(physical.contains("\tcoating_strain_sufficient_common_amplitude_scale\t"));
        assert!(physical.contains("\tfilm0_strain_sufficient_amplitude_scale\t"));
        for normalized in ["drive_bridge", "mobility", "_per_N", "_at_1N", "peak_force_N"] {
            assert!(!physical.contains(normalized), "{normalized}");
        }
    }
}
