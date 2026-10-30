//! Material-derived planar exterior loading of the original shared body.
//!
//! The retained DKT displacement trace supplies the SAME signed, mass-normalized
//! velocity fields used by the bridges and material loss. The existing Rayleigh
//! owner supplies their full reciprocal impedance; the existing PSD-residue
//! fitter supplies a causal passive approximation. No new eigenproblem, modal
//! damping ratio, acoustic gain, or time integrator is introduced.
use super::super::{ModalPreparation, ReciprocalStringBank, ReciprocalStringError};
use crate::render::plate::impact::radiation::{Model, fit};
use crate::thin_plate::RetainedPlateBasis;
use fs_bem::helmholtz::{HelmholtzError, Medium, rayleigh_modal::{AreaSample, PreparedSurface}};
use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_plate::PlateMesh;

/// Explicit cold geometry, frequency and realization work limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlateRayleighOptions {
    /// One or two independent infinite baffled half-spaces. With a separately
    /// modeled interior cavity, select one exterior half-space.
    pub sides: u8,
    /// Positive lower edge of the sampled impedance band [Hz].
    pub minimum_hz: f64,
    /// Upper edge [Hz], covering every retained body reference oscillator.
    pub maximum_hz: f64,
    /// Odd 33..=257 samples; even rows train and odd rows independently check.
    pub frequency_samples: usize,
    /// All triangles contribute three area points; no cell is discarded.
    pub maximum_surface_cells: usize,
    /// Ceiling on frequencies * (cells^2*modes + cells*modes^2).
    pub maximum_sweep_terms: u64,
    /// Complete fitted poles, before the shared exact owner's stricter capacity.
    pub maximum_poles: usize,
}

impl PlateRayleighOptions {
    /// Admission before source I/O or spectral preparation. The 8x off-band
    /// storage poles used by the shared fitter must resolve on the SAME clock.
    pub fn validate(&self, sample_rate_hz: u32) -> Result<(), PlateRayleighError> {
        if !matches!(self.sides, 1 | 2)
            || !self.minimum_hz.is_finite() || self.minimum_hz <= 0.0
            || !self.maximum_hz.is_finite() || self.maximum_hz <= self.minimum_hz
            || !(33..=257).contains(&self.frequency_samples) || self.frequency_samples % 2 == 0
            || !(3..=4096).contains(&self.maximum_surface_cells)
            || !(1..=1_000_000_000).contains(&self.maximum_sweep_terms)
            || !(1..=256).contains(&self.maximum_poles)
            || sample_rate_hz == 0
            || 8.0 * core::f64::consts::TAU * self.maximum_hz
                >= 0.9 * core::f64::consts::PI * f64::from(sample_rate_hz)
        {
            return Err(PlateRayleighError::Invalid(
                "planar radiation needs 1 or 2 sides, a positive resolved band, odd 33..257 frequencies and finite cell/work/pole ceilings"));
        }
        Ok(())
    }
}

/// Original material/geometry, coherent gas, and sampled fit evidence.
/// These are numerical estimates, not a continuum or measured-instrument claim.
#[derive(Clone, Debug)]
pub struct PlateRayleighReport {
    /// Exact regional source used for the retained mechanical eigenbasis.
    pub specimen_identity: ContentHash,
    /// Same atmosphere as strings, body cavity and the existing pressure observer.
    pub gas: GasState,
    /// Original complete numerical selection, without relaxed error tolerances.
    pub options: PlateRayleighOptions,
    /// All retained structural body modes, including signed cancelling fields.
    pub body_modes: usize,
    /// Actual quadrature cells in the existing plate mesh.
    pub surface_cells: usize,
    /// Realized acoustic oscillator count; each contributes one q/p pair.
    pub acoustic_poles: usize,
    /// Complete estimated scalar work of the frequency sweep.
    pub sweep_terms: u64,
    /// Existing Rayleigh area-resolution diagnostic at the upper band edge.
    pub minimum_equivalent_disc_ppw: f64,
    /// Dictionary spacing set by the body's largest cell separation and gas c.
    pub resonator_spacing_rad_s: f64,
    /// Maximum whole-matrix complex error, normalized by the sampled load.
    pub peak_error: f64,
    /// Whole-matrix RMS complex error under the original fitter's limit.
    pub rms_error: f64,
    /// Independent resistance error; weak radiation loss cannot hide in inertia.
    pub resistance_peak_error: f64,
    /// Original fitter's RMS resistance error.
    pub resistance_rms_error: f64,
}

/// Preparation/attachment refusals retain the original physical owner.
#[derive(Debug)]
pub enum PlateRayleighError {
    /// Incomplete geometry, unsupported source, clock or work budget.
    Invalid(&'static str),
    /// Original Rayleigh geometry/frequency owner refused.
    Surface(HelmholtzError),
    /// Original positive-real fit refused its unchanged error gates.
    Fit(String),
    /// The complete load exceeded the mechanical owner's admission.
    Bank(ReciprocalStringError),
    /// Cancellation before publication of a complete loaded bank.
    Cancelled,
}
impl core::fmt::Display for PlateRayleighError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material plate radiation: {self:?}")
    }
}
impl std::error::Error for PlateRayleighError {}

fn poll(gate: &CancelGate) -> Result<(), PlateRayleighError> {
    if gate.is_requested() { Err(PlateRayleighError::Cancelled) } else { Ok(()) }
}

impl ReciprocalStringBank {
    /// Compute and cold-attach a passive exterior reaction from THIS regional
    /// material plate, its original retained modes and the selected atmosphere.
    ///
    /// All signed cross-mode loads participate in the existing simultaneous bow,
    /// string, bridge, cavity and material-memory solve. The sampled band must
    /// cover all retained free-body frequencies; this is not a bound on every
    /// coupled instrument pole or its out-of-band forcing. Geometry/gas/material
    /// changes require a fresh realization and automatically rebuild the load.
    ///
    /// Scope: piecewise-linear transverse traces of a flat DKT body in one or
    /// two independent infinite baffles. The optional material plate observer
    /// consumes this same signed surface basis through the shared Rayleigh
    /// receiver. Without it, output remains compact loaded acceleration. This
    /// is not a curved/unbaffled violin exterior or a
    /// quadrature-error certificate. Shells, prepared-only propagation, and a
    /// second half-space over a modeled cavity refuse explicitly.
    pub fn with_material_plate_rayleigh(
        self, options: PlateRayleighOptions, gate: &CancelGate,
    ) -> Result<(Self, PlateRayleighReport), PlateRayleighError> {
        let (model, report) = self.prepare_material_plate_rayleigh(options, gate)?;
        let mut bank = self.with_body_radiation_load(model, options.maximum_poles, gate)
            .map_err(PlateRayleighError::Bank)?;
        bank.source.body_radiation.as_mut().ok_or(PlateRayleighError::Invalid(
            "computed planar radiation lost its original source binding"))?.plate_rayleigh = Some(report.clone());
        poll(gate)?;
        Ok((bank, report))
    }

    /// Prepare the complete positive-real image without changing the bank.
    /// The returned model uses every mode of its original body component, in
    /// exactly the same order. It must not be reused on a different material
    /// realization; `with_material_plate_rayleigh` keeps that binding automatic.
    pub fn prepare_material_plate_rayleigh(
        &self, options: PlateRayleighOptions, gate: &CancelGate,
    ) -> Result<(Model, PlateRayleighReport), PlateRayleighError> {
        poll(gate)?;
        options.validate(self.source.options.sample_rate_hz)?;
        if self.system.samples_rendered() != 0
            || matches!(self.source.prepared_linear, ModalPreparation::Prepared(_)) {
            return Err(PlateRayleighError::Invalid("computed radiation requires an unstarted ordinary or harmonic-only modal bank"));
        }
        if options.sides == 2 && self.source.material_cavity.is_some() {
            return Err(PlateRayleighError::Invalid("a modeled interior cavity permits only one exterior half-space"));
        }
        let body = self.source.material_chart.as_ref().ok_or(PlateRayleighError::Invalid(
            "computed planar radiation requires an original regional material plate chart"))?;
        let basis = self.source.material_chart_basis.as_ref().ok_or(PlateRayleighError::Invalid(
            "computed radiation requires the original retained body eigenbasis"))?;
        let modes = self.system.components().get(self.source.strings.len())
            .ok_or(PlateRayleighError::Invalid("missing original structural body component"))?.modes();
        if basis.source_identity != body.specimen.specimen_identity()
            || basis.mass_normalized_modes.len() != modes.len()
            || modes.iter().any(|m| {
                let hz = m.angular_frequency_rad_s / core::f64::consts::TAU;
                hz < options.minimum_hz || hz > options.maximum_hz
            })
        {
            return Err(PlateRayleighError::Invalid("radiation must cover the same complete material body basis and its reference frequencies"));
        }
        let mesh = &body.specimen.chart().mesh;
        let cells = mesh.tris.len().checked_mul(3)
            .ok_or(PlateRayleighError::Invalid("radiation cell count overflow"))?;
        let n = modes.len() as u64;
        let c = cells as u64;
        let sweep_terms = c.checked_mul(n).and_then(|v| c.checked_add(n).and_then(|s| v.checked_mul(s)))
            .and_then(|v| v.checked_mul(options.frequency_samples as u64))
            .ok_or(PlateRayleighError::Invalid("radiation sweep work overflow"))?;
        if cells > options.maximum_surface_cells || sweep_terms > options.maximum_sweep_terms {
            return Err(PlateRayleighError::Invalid("complete plate radiation exceeds its cell or sweep-work ceiling"));
        }
        let samples = surface_samples(mesh, basis, gate)?;
        let surface = PreparedSurface::new(&samples).map_err(PlateRayleighError::Surface)?;
        let gas = self.source.gas;
        gas.validate_coherence().map_err(|_| PlateRayleighError::Invalid("incoherent exterior gas"))?;
        let medium = Medium { density: gas.density, sound_speed: gas.sound_speed };
        let low = core::f64::consts::TAU * options.minimum_hz;
        let high = core::f64::consts::TAU * options.maximum_hz;
        let ppw = surface.equivalent_disc_ppw(high, medium).map_err(PlateRayleighError::Surface)?;
        // Linear samples resolve physical travel-time detail at every frequency.
        let mut omega: Vec<_> = (0..options.frequency_samples).map(|i|
            low + (high-low) * i as f64 / (options.frequency_samples-1) as f64).collect();
        omega[0] = low;
        omega[options.frequency_samples-1] = high;
        let extent = surface.maximum_distance_m();
        let spacing = if extent > 0.0 {
            (core::f64::consts::TAU * gas.sound_speed / (6.0 * extent)).min(high-low)
        } else { high-low };
        let mut impedance = Vec::with_capacity(omega.len());
        for &w in &omega {
            poll(gate)?;
            let mut row = surface.power_consistent_impedance(w, medium)
                .map_err(PlateRayleighError::Surface)?;
            for value in &mut row { *value = value.scale(f64::from(options.sides)); }
            impedance.push(row);
        }
        poll(gate)?;
        let fitted = fit::fit_resolved(&omega, &impedance, modes.len(), spacing)
            .map_err(PlateRayleighError::Fit)?;
        poll(gate)?;
        if fitted.model.poles.len() > options.maximum_poles {
            return Err(PlateRayleighError::Invalid("complete passive radiation fit exceeds the declared pole ceiling"));
        }
        let report = PlateRayleighReport {
            specimen_identity: basis.source_identity, gas, options,
            body_modes: modes.len(), surface_cells: cells, acoustic_poles: fitted.model.poles.len(),
            sweep_terms, minimum_equivalent_disc_ppw: ppw, resonator_spacing_rad_s: spacing,
            peak_error: fitted.peak_error, rms_error: fitted.rms_error,
            resistance_peak_error: fitted.resistance_peak_error,
            resistance_rms_error: fitted.resistance_rms_error,
        };
        Ok((fitted.model, report))
    }
}

// Three-point degree-two triangle quadrature. It integrates the P1 trace's
// displacement/pressure products consistently with the original plate ports;
// no normalization by signed radiating area can erase a cancelling mode.
pub(in crate::bowed_string::runtime::reciprocal::bank) fn surface_samples(
    mesh: &PlateMesh, basis: &RetainedPlateBasis, gate: &CancelGate)
    -> Result<Vec<AreaSample>, PlateRayleighError>
{
    if basis.dof_map.len() != 3 * mesh.nodes.len() || basis.mass_normalized_modes.is_empty() {
        return Err(PlateRayleighError::Invalid("radiation trace requires the original complete nodal basis"));
    }
    let mut samples = Vec::with_capacity(3 * mesh.tris.len());
    for tri in &mesh.tris {
        poll(gate)?;
        if tri.iter().any(|&i| i >= mesh.nodes.len()) {
            return Err(PlateRayleighError::Invalid("radiation triangle names an absent original node"));
        }
        let [a,b,c] = tri.map(|i| mesh.nodes[i]);
        let area = 0.5 * ((b.0-a.0)*(c.1-a.1) - (b.1-a.1)*(c.0-a.0));
        if !area.is_finite() || area <= 0.0 {
            return Err(PlateRayleighError::Invalid("radiation requires positive oriented plate triangles"));
        }
        for weights in [[2.0/3.0,1.0/6.0,1.0/6.0], [1.0/6.0,2.0/3.0,1.0/6.0], [1.0/6.0,1.0/6.0,2.0/3.0]] {
            let mut position_m = [0.0;2];
            for (&node, weight) in tri.iter().zip(weights) {
                position_m[0] += weight * mesh.nodes[node].0;
                position_m[1] += weight * mesh.nodes[node].1;
            }
            let mode_shape = basis.mass_normalized_modes.iter().map(|mode| {
                tri.iter().zip(weights).try_fold(0.0, |sum, (&node,weight)| {
                    let value = match basis.dof_map[3*node] {
                        Some(dof) => *mode.get(dof).ok_or(PlateRayleighError::Invalid("incomplete retained radiation eigenvector"))?,
                        None => 0.0,
                    };
                    Ok(sum + weight * value)
                })
            }).collect::<Result<Vec<_>, PlateRayleighError>>()?;
            samples.push(AreaSample { position_m, area_m2: area/3.0, mode_shape });
        }
    }
    Ok(samples)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn original_signed_trace_preserves_volume_and_consistent_work() {
        let mesh = PlateMesh::rectangle(0.4, 0.3, 2, 2);
        let basis = RetainedPlateBasis {
            source_identity: ContentHash([1;32]),
            dof_map: (0..3*mesh.nodes.len()).map(Some).collect(),
            mass_normalized_modes: vec![mesh.nodes.iter().flat_map(|&(x,_)| [x-0.2,0.0,0.0]).collect()],
        };
        let cells = surface_samples(&mesh, &basis, &CancelGate::new_clock_free()).unwrap();
        assert_eq!(cells.len(), 3*mesh.tris.len());
        let area: f64 = cells.iter().map(|c|c.area_m2).sum();
        let volume: f64 = cells.iter().map(|c|c.area_m2*c.mode_shape[0]).sum();
        let work: f64 = cells.iter().map(|c|c.area_m2*c.mode_shape[0].powi(2)).sum();
        assert!((area-0.12).abs() < 1e-15);
        assert!(volume.abs() < 1e-16);
        assert!((work-0.3*0.4_f64.powi(3)/12.0).abs() < 1e-16);
        let surface = PreparedSurface::new(&cells).unwrap();
        let air = Medium { density: 1.2, sound_speed: 343.0 };
        let z = surface.power_consistent_impedance(1000.0, air).unwrap()[0];
        assert!(z.re > 0.0, "a cancelling surface field still radiates");
        let dense = Medium { density: 2.4, ..air };
        let twice = surface.power_consistent_impedance(1000.0, dense).unwrap()[0];
        assert!((twice.re-2.0*z.re).abs() < 1e-13);
        assert!((twice.im-2.0*z.im).abs() < 1e-13);
    }
}
