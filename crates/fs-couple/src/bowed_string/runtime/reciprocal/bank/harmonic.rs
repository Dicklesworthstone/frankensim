//! Bridge-to-sound response of the original source-bound performer.
//! No body, sympathetic string, cavity, Maxwell history or selected air pole is
//! removed. Contacts are explicitly released; this is not bowed linearization.

use super::multiple::MultiBowStringBank;
use super::*;
use crate::render::schedule::force::coupled::contact::friction::harmonic::ReleasedHarmonicResponse;
use crate::shear_shell::ResolvedShearShell;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_math::c64::C64;

/// Peak phasors in exp(-i omega t), for the caller's physical bridge force.
#[derive(Clone, Debug, PartialEq)]
pub struct ShellBridgeHarmonicResponse {
    // Original body coordinates inside the already retained coupled response.
    // Keep this private so an observer never has to reconstruct string offsets.
    body_mode_range: core::ops::Range<usize>,
    /// All original coordinates, losses, material/air histories and pressure.
    pub released: ReleasedHarmonicResponse,
    /// Physical velocities in the original string/bridge declaration order.
    pub bridge_velocity_m_s: Vec<C64>,
    /// Signed structural exterior displacement-volume rate; necks are separate.
    pub exterior_volume_velocity_m3_s: C64,
    /// Signed structural outward-from-cavity displacement-volume rate.
    pub cavity_volume_velocity_m3_s: C64,
    /// Heat of the actually bound positive Prony coating branches [W].
    /// Equal-time factors mix films, so no per-film allocation is invented.
    pub coating_dissipation_w: f64,
    /// Mean AC storage of those same coating branches [J].
    pub coating_mean_stored_energy_j: f64,
    /// All remaining material-branch heat (strings, gas thermal memory, etc.).
    pub other_material_dissipation_w: f64,
    /// Intersection of retained source/fitted-operator frequency supports.
    /// Discretization and finite-amplitude validity are separate questions.
    pub source_band_rad_s: (f64, f64),
}

impl ShellBridgeHarmonicResponse {
    /// Borrow the original retained shell displacement without allocating or
    /// reconstructing it from bridge velocities. This is the body basis used
    /// by the source image's coating-strain observer, with every connected
    /// string/cavity/load contribution already present in the solved motion.
    /// The checked slice also refuses a caller-truncated public response.
    pub fn body_displacement_m_sqrt_kg(&self) -> Result<&[C64], ReciprocalStringError> {
        self.released.dynamics.displacement_m_sqrt_kg
            .get(self.body_mode_range.clone())
            .ok_or_else(|| invalid("harmonic response lost its original body coordinates"))
    }
}

#[derive(Clone, Copy)]
struct Trace<'a> {
    ports: &'a [Vec<f64>],
    exterior: &'a [f64],
    cavity: &'a [f64],
    source_band: (f64, f64),
    coating_histories: usize,
}

fn intersect(band: &mut (f64, f64), next: (f64, f64)) -> Result<(), ReciprocalStringError> {
    if !next.0.is_finite() || next.0 < 0.0 || next.1.is_nan() || next.1 < next.0 {
        return Err(invalid(
            "harmonic source has an invalid angular-frequency support",
        ));
    }
    band.0 = band.0.max(next.0);
    band.1 = band.1.min(next.1);
    if band.1 < band.0 {
        return Err(invalid(
            "harmonic sources have no common angular-frequency support",
        ));
    }
    Ok(())
}

fn material_band(
    source: &ResolvedMaterialStatePoint,
    band: &mut (f64, f64),
) -> Result<(), ReciprocalStringError> {
    for property in source.properties() {
        if let Some(next) = property.angular_frequency_band_rad_s().map_err(|_| {
            invalid("harmonic material property has invalid frequency-axis semantics")
        })? {
            intersect(band, next)?;
        }
    }
    Ok(())
}

fn string_band(
    source: &MaterialStringSource,
    band: &mut (f64, f64),
) -> Result<(), ReciprocalStringError> {
    match source {
        MaterialStringSource::WithPolarizations(p) => string_band(p.source(), band)?,
        MaterialStringSource::WithTorsion(t) => {
            string_band(t.bending(), band)?;
            intersect(band, t.section().omega_band_rad_s())?;
            if let Some(law) = t.relaxation() {
                intersect(band, law.omega_band_rad_s())?;
            }
        }
        MaterialStringSource::Homogeneous(s) => material_band(s.material(), band)?,
        MaterialStringSource::Wound(s) => intersect(band, s.section().omega_band_rad_s())?,
        MaterialStringSource::Concentric { section, .. } => {
            intersect(band, section.omega_band_rad_s())?;
            for layer in section.elastic().layers() {
                material_band(&layer.material, band)?;
            }
        }
        MaterialStringSource::ConcentricProny { section, .. } => {
            intersect(band, section.omega_band_rad_s())?;
            for layer in section.elastic().layers() {
                material_band(&layer.material, band)?;
            }
        }
    }
    let string = source.resolved_transverse_string()?;
    if let Some(law) = string.kelvin_voigt_bending {
        intersect(band, law.omega_band_rad_s)?;
    }
    if let Some(law) = string.relaxation_bending {
        intersect(band, law.omega_band_rad_s)?;
    }
    Ok(())
}

fn shell_band(
    shell: &ResolvedShearShell,
    band: &mut (f64, f64),
) -> Result<(), ReciprocalStringError> {
    if !shell.coatings().is_empty() && shell.coating_relaxation_options().is_none() {
        return Err(invalid(
            "loaded coating response requires the explicit hereditary shell option",
        ));
    }
    for region in shell.regions() {
        intersect(band, region.omega_band_rad_s())?;
    }
    for interface in shell.interfaces() {
        intersect(band, interface.omega_band_rad_s())?;
    }
    Ok(())
}

impl MaterialStringBankSource {
    fn harmonic_trace(
        &self,
        omega: f64,
        gate: &CancelGate,
    ) -> Result<Trace<'_>, ReciprocalStringError> {
        poll(gate)?;
        if !omega.is_finite() || omega <= 0.0 {
            return Err(invalid(
                "loaded harmonic response requires positive finite angular frequency",
            ));
        }
        if self.residual_flexibility {
            return Err(invalid(
                "loaded bridge response requires explicit residual-flexibility force columns",
            ));
        }
        self.gas
            .validate_coherence()
            .map_err(|_| invalid("loaded harmonic response requires the original coherent gas"))?;
        let mut band = (0.0, f64::INFINITY);
        let (ports, exterior, cavity, projected_memory) =
            if let Some(image) = &self.material_shell_image {
                shell_band(image.source(), &mut band)?;
                (
                    image.force_ports(),
                    image.volume_shapes_m2_per_sqrt_kg(),
                    image.cavity_volume_shapes_m2_per_sqrt_kg(),
                    image.coating_relaxation(),
                )
            } else if let Some(image) = &self.material_shell_assembly_image {
                for part in &image.source().components {
                    poll(gate)?;
                    shell_band(&part.specimen, &mut band)?;
                }
                (
                    image.force_ports(),
                    image.volume_shapes_m2_per_sqrt_kg(),
                    image.cavity_volume_shapes_m2_per_sqrt_kg(),
                    image.coating_relaxation(),
                )
            } else {
                return Err(invalid(
                    "loaded bridge response requires its original material shell image",
                ));
            };
        if projected_memory.is_some() != self.material_shell_coating_relaxation.is_some() {
            return Err(invalid(
                "loaded response requires coating memory to be bound before analysis",
            ));
        }
        let coating_histories = if let Some(memory) = &self.material_shell_coating_relaxation {
            intersect(&mut band, memory.source_band_rad_s())?;
            memory.relaxation_times_s().len()
        } else {
            0
        };
        for member in &self.strings {
            poll(gate)?;
            string_band(&member.source, &mut band)?;
        }
        for binding in &self.material_bridge_interfaces {
            intersect(&mut band, binding.interface.omega_band_rad_s())?;
        }
        if let Some(memory) = &self.string_bending_relaxation {
            for &support in memory.omega_bands_rad_s() {
                intersect(&mut band, support)?;
            }
        }
        if let Some(fluid) = &self.fluid_loading {
            for string in fluid.strings() {
                poll(gate)?;
                string
                    .flow()
                    .impedance_per_length(omega)
                    .map_err(ReciprocalStringError::Material)?;
            }
        }
        if let Some(air) = &self.body_radiation {
            let support = if let Some(report) = air.material_shell_bem() {
                (report.options.minimum_hz, report.options.maximum_hz)
            } else if let Some(report) = air.material_plate_rayleigh() {
                (report.options.minimum_hz, report.options.maximum_hz)
            } else {
                return Err(invalid(
                    "loaded source response requires the exterior model's retained source-fit band",
                ));
            };
            intersect(
                &mut band,
                (
                    support.0 * core::f64::consts::TAU,
                    support.1 * core::f64::consts::TAU,
                ),
            )?;
        }
        if omega < band.0 || omega > band.1 {
            return Err(invalid(
                "loaded harmonic drive exceeds retained string, wood, film, interface or exterior support",
            ));
        }
        Ok(Trace {
            ports,
            exterior,
            cavity,
            source_band: band,
            coating_histories,
        })
    }
}

fn prepare_force(
    source: &MaterialStringBankSource,
    components: &[ModalAcousticTimeModel],
    modes: usize,
    trace: &Trace<'_>,
    bridge: usize,
    peak_force_n: C64,
) -> Result<(usize, Vec<C64>), ReciprocalStringError> {
    if !peak_force_n.re.is_finite() || !peak_force_n.im.is_finite() {
        return Err(invalid("loaded harmonic bridge force must be finite"));
    }
    let body = source.strings.len();
    let component = components
        .get(body)
        .ok_or_else(|| invalid("loaded harmonic response lost the original body component"))?;
    let row = trace
        .ports
        .get(bridge)
        .ok_or_else(|| invalid("loaded harmonic bridge index is absent from this body"))?;
    let n = row.len();
    if n == 0
        || n > component.modes().len()
        || trace.exterior.len() != n
        || trace.cavity.len() != n
        || trace.ports.iter().any(|port| port.len() != n)
    {
        return Err(invalid(
            "loaded harmonic response lost complete shell force and volume maps",
        ));
    }
    let start: usize = components[..body]
        .iter()
        .map(|component| component.modes().len())
        .sum();
    if start.checked_add(n).is_none_or(|end| end > modes) {
        return Err(invalid(
            "loaded harmonic body extent exceeds original mechanical coordinates",
        ));
    }
    let mut force = vec![C64::ZERO; modes];
    for (out, &shape) in force[start..start + n].iter_mut().zip(row) {
        *out = peak_force_n.scale(shape);
    }
    Ok((start, force))
}

fn prepare_bridge_forces(
    source: &MaterialStringBankSource, components: &[ModalAcousticTimeModel], modes: usize,
    trace: &Trace<'_>, peak_forces_n: &[C64],
) -> Result<(usize, Vec<C64>), ReciprocalStringError> {
    if peak_forces_n.is_empty() || peak_forces_n.len() != trace.ports.len()
        || peak_forces_n.iter().any(|f| !f.re.is_finite() || !f.im.is_finite())
    {
        return Err(invalid("coherent shell forcing needs one finite complex force per original bridge"));
    }
    let (start, mut force) = prepare_force(source, components, modes, trace, 0, peak_forces_n[0])?;
    let n = trace.exterior.len();
    for (&drive, row) in peak_forces_n[1..].iter().zip(&trace.ports[1..]) {
        for (out, &shape) in force[start..start + n].iter_mut().zip(row) {
            *out = *out + drive.scale(shape);
        }
    }
    if force.iter().any(|f| !f.re.is_finite() || !f.im.is_finite()) {
        return Err(invalid("coherent shell force superposition is not representable"));
    }
    Ok((start, force))
}

fn prepare_force_matrix(
    source: &MaterialStringBankSource, components: &[ModalAcousticTimeModel], modes: usize,
    trace: &Trace<'_>, gate: &CancelGate,
) -> Result<(usize, Vec<Vec<C64>>), ReciprocalStringError> {
    use crate::render::schedule::force::coupled::harmonic::{
        MAX_HARMONIC_RIGHT_HAND_SIDES, MAX_HARMONIC_UNKNOWNS,
    };
    if trace.ports.is_empty() || trace.ports.len() > MAX_HARMONIC_RIGHT_HAND_SIDES
        || modes > MAX_HARMONIC_UNKNOWNS
    {
        return Err(invalid("complete bridge matrix exceeds the original harmonic source capacity"));
    }
    let mut forces = Vec::with_capacity(trace.ports.len());
    let mut start = 0;
    for bridge in 0..trace.ports.len() {
        poll(gate)?;
        let (offset, force) = prepare_force(source, components, modes, trace, bridge, C64::ONE)?;
        start = offset;
        forces.push(force);
    }
    Ok((start, forces))
}

fn finish(
    source: &MaterialStringBankSource,
    cards: &[BowedStringCard],
    trace: Trace<'_>,
    start: usize,
    released: ReleasedHarmonicResponse,
    gate: &CancelGate,
) -> Result<ShellBridgeHarmonicResponse, ReciprocalStringError> {
    poll(gate)?;
    let dynamics = &released.dynamics;
    let velocity = &dynamics.velocity_m_sqrt_kg_per_s[start..start + trace.exterior.len()];
    let project = |row: &[f64]| {
        row.iter()
            .zip(velocity)
            .fold(C64::ZERO, |sum, (&shape, &v)| sum + v.scale(shape))
    };
    let range = source.shell_coating_memory_range(cards, gate)?.unwrap_or(0..0);
    if range.len() != trace.coating_histories {
        return Err(invalid("coating harmonic observation lost its source history range"));
    }
    let (begin, end) = (range.start, range.end);
    let heat = dynamics
        .material_branch_dissipation_w
        .get(begin..end)
        .ok_or_else(|| invalid("coating harmonic heat lost its original branch order"))?;
    let energy = dynamics
        .material_branch_mean_storage_j
        .get(begin..end)
        .ok_or_else(|| invalid("coating harmonic storage lost its original branch order"))?;
    let coating_dissipation_w: f64 = heat.iter().sum();
    // Sum the complement directly, avoiding cancellation of two nearly equal
    // positive aggregate losses when varnish dominates the material heat.
    let other_material_dissipation_w = dynamics.material_branch_dissipation_w[..begin]
        .iter()
        .chain(&dynamics.material_branch_dissipation_w[end..])
        .sum();
    let bridge_velocity_m_s: Vec<_> = trace.ports.iter().map(|row| project(row)).collect();
    let exterior_volume_velocity_m3_s = project(trace.exterior);
    let cavity_volume_velocity_m3_s = project(trace.cavity);
    let coating_mean_stored_energy_j = energy.iter().sum();
    if bridge_velocity_m_s
        .iter()
        .copied()
        .chain([exterior_volume_velocity_m3_s, cavity_volume_velocity_m3_s])
        .any(|value| !value.re.is_finite() || !value.im.is_finite())
    {
        return Err(invalid(
            "loaded harmonic physical trace is not representable",
        ));
    }
    Ok(ShellBridgeHarmonicResponse {
        body_mode_range: start..start + trace.exterior.len(),
        released,
        bridge_velocity_m_s,
        exterior_volume_velocity_m3_s,
        cavity_volume_velocity_m3_s,
        coating_dissipation_w,
        coating_mean_stored_energy_j,
        other_material_dissipation_w,
        source_band_rad_s: trace.source_band,
    })
}

impl ReciprocalStringBank {
    /// Apply a peak physical force at one retained shell bridge and inspect the
    /// entire released performer. The actual positive Prony film law, all
    /// strings, cavity and selected exterior reaction are retained. Pressure is
    /// the existing bank observation; no finite-distance receiver is inferred.
    /// Source-band admission does not certify mesh/modal or amplitude accuracy.
    /// A residual-flexibility extension refuses until its complete physical
    /// bridge-force columns can be supplied; retained-mode forcing alone would
    /// silently omit that extension's directly driven displacement.
    pub fn material_shell_bridge_harmonic_response(
        &self,
        omega_rad_s: f64,
        bridge_index: usize,
        peak_force_n: C64,
        maximum_work: usize,
        gate: &CancelGate,
    ) -> Result<ShellBridgeHarmonicResponse, ReciprocalStringError> {
        let trace = self.source.harmonic_trace(omega_rad_s, gate)?;
        let (start, force) = prepare_force(
            &self.source,
            self.system.components(),
            self.system.mode_count(),
            &trace,
            bridge_index,
            peak_force_n,
        )?;
        let response = self
            .system
            .released_harmonic_response(omega_rad_s, &force, maximum_work, gate)
            .map_err(ReciprocalStringError::Friction)?;
        finish(&self.source, &self.cards, trace, start, response, gate)
    }

    /// Coherent physical peak forces in original bridge order. Their signed
    /// complex columns are combined BEFORE solving and accounting, preserving
    /// interference in pressure, stored energy and every heat channel.
    pub fn material_shell_bridge_forces_harmonic_response(
        &self, omega_rad_s: f64, peak_forces_n: &[C64], maximum_work: usize,
        gate: &CancelGate,
    ) -> Result<ShellBridgeHarmonicResponse, ReciprocalStringError> {
        let trace = self.source.harmonic_trace(omega_rad_s, gate)?;
        let (start, force) = prepare_bridge_forces(&self.source, self.system.components(),
            self.system.mode_count(), &trace, peak_forces_n)?;
        let response = self.system.released_harmonic_response(omega_rad_s, &force, maximum_work, gate)
            .map_err(ReciprocalStringError::Friction)?;
        finish(&self.source, &self.cards, trace, start, response, gate)
    }

    /// Unit-peak-force columns of the full physical bridge response matrix.
    /// Column j drives original bridge j by 1 N and observes EVERY bridge plus
    /// pressure and losses. All columns share one operator assembly/LU; source
    /// admission, contact refusals and accounting are the single-drive rules.
    /// Losses belong to individual unit load cases and cannot be superposed.
    pub fn material_shell_bridge_harmonic_matrix(
        &self, omega_rad_s: f64, maximum_work: usize, gate: &CancelGate,
    ) -> Result<Vec<ShellBridgeHarmonicResponse>, ReciprocalStringError> {
        let trace = self.source.harmonic_trace(omega_rad_s, gate)?;
        let (start, forces) = prepare_force_matrix(&self.source, self.system.components(),
            self.system.mode_count(), &trace, gate)?;
        let rows: Vec<_> = forces.iter().map(Vec::as_slice).collect();
        self.system.released_harmonic_responses(omega_rad_s, &rows, maximum_work, gate)
            .map_err(ReciprocalStringError::Friction)?.into_iter()
            .map(|response| finish(&self.source, &self.cards, trace, start, response, gate)).collect()
    }
}

impl MultiBowStringBank {
    /// Same loaded response on the actual multiple-bow performer, with every
    /// contact explicitly released and every original mechanical coordinate
    /// preserved. Bound nonlinear/contact storage refuses rather than vanish.
    pub fn material_shell_bridge_harmonic_response(
        &self,
        omega_rad_s: f64,
        bridge_index: usize,
        peak_force_n: C64,
        maximum_work: usize,
        gate: &CancelGate,
    ) -> Result<ShellBridgeHarmonicResponse, ReciprocalStringError> {
        let source = self.source().bank();
        let trace = source.harmonic_trace(omega_rad_s, gate)?;
        let (start, force) = prepare_force(
            source,
            self.system().components(),
            self.system().mode_count(),
            &trace,
            bridge_index,
            peak_force_n,
        )?;
        let response = self
            .system()
            .released_harmonic_response(omega_rad_s, &force, maximum_work, gate)
            .map_err(ReciprocalStringError::Friction)?;
        finish(source, self.cards(), trace, start, response, gate)
    }

    /// Coherent drive of every original bridge on the actual multiple-bow bank.
    /// Combined phasors are solved once, so cross terms in loss and power remain
    /// physical. Contacts are released with the same admission as a single drive.
    pub fn material_shell_bridge_forces_harmonic_response(
        &self, omega_rad_s: f64, peak_forces_n: &[C64], maximum_work: usize,
        gate: &CancelGate,
    ) -> Result<ShellBridgeHarmonicResponse, ReciprocalStringError> {
        let source = self.source().bank();
        let trace = source.harmonic_trace(omega_rad_s, gate)?;
        let (start, force) = prepare_bridge_forces(source, self.system().components(),
            self.system().mode_count(), &trace, peak_forces_n)?;
        let response = self.system().released_harmonic_response(omega_rad_s, &force, maximum_work, gate)
            .map_err(ReciprocalStringError::Friction)?;
        finish(source, self.cards(), trace, start, response, gate)
    }

    /// Full bridge mobility/pressure matrix on the complete multiple-bow bank.
    /// Returned unit-drive columns retain declaration order and share one LU.
    /// The allowance covers the entire frequency's batch; no mode or port is
    /// omitted to meet it. Each column retains its own source and power checks.
    pub fn material_shell_bridge_harmonic_matrix(
        &self, omega_rad_s: f64, maximum_work: usize, gate: &CancelGate,
    ) -> Result<Vec<ShellBridgeHarmonicResponse>, ReciprocalStringError> {
        let source = self.source().bank();
        let trace = source.harmonic_trace(omega_rad_s, gate)?;
        let (start, forces) = prepare_force_matrix(source, self.system().components(),
            self.system().mode_count(), &trace, gate)?;
        let rows: Vec<_> = forces.iter().map(Vec::as_slice).collect();
        self.system().released_harmonic_responses(omega_rad_s, &rows, maximum_work, gate)
            .map_err(ReciprocalStringError::Friction)?.into_iter()
            .map(|response| finish(source, self.cards(), trace, start, response, gate)).collect()
    }
}
