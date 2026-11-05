//! Source-bound thermoelastic memory on the original shear-shell modal basis.
//!
//! Homogeneous symmetric regions retain the insulated slab's reversible
//! membrane increment and full bending diffusion memory. The existing MITC
//! owner supplies physical strain factors; the existing QR owner compresses
//! their rows, and the coupled held-force owner advances every state and heat.
//! No additional eigenbasis, scalar thermal damping ratio or time engine is used.

use super::{MaterialShellImage, invalid, owner, poll};
use crate::acoustic_realize::AcousticRealizeError;
use crate::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use crate::render::schedule::force::coupled::{
    CoupledModalSystem, MAX_MODAL_CONNECTIONS, ModalColumnConnection, ModalCouplingConfig,
    ModalPreparation, PreparedModalBudget,
};
use crate::thin_plate::thermoelastic_memory::compress_rows;
use fs_exec::CancelGate;
use fs_material::visco::slab_ted::SlabPronySpectrum;
use fs_math::{c64::C64, det};
use fs_phs::RelaxationBranch;
use fs_plate::shell::mitc::thermoelastic::{
    ShellThermoelasticFactors, project_thermoelastic_factors, projection_terms,
};

mod source;
pub use source::{ShellThermalExpansionProperties, ShellThermalProperties};

/// Explicit constitutive accuracy and cold preparation resources.
#[derive(Clone, Copy, Debug)]
pub struct ShellThermoelasticOptions {
    /// Requested angular-frequency band [rad/s]. Each original source must
    /// cover the nominal modes and represented instantaneous stiffness bound.
    pub band_rad_s: (f64, f64),
    /// Absolute normalized constitutive-response tolerance per regional law.
    /// This does not bound relative error in an almost-null loss direction.
    pub normalized_tolerance: f64,
    /// Through-thickness diffusion terms per distinct regional law, 1..=64.
    pub max_thermal_terms: usize,
    /// Complete bending-memory state allowance, 0..=8192.
    pub max_memory_states: usize,
    /// Original MITC reads, row compression and memory expansion allowance.
    pub max_projection_terms: usize,
    /// Conservative live f64 factor-storage allowance, including original
    /// quadrature, QR workspace, retained membrane columns and memory maps.
    /// Excludes the original immutable mechanical image and the element
    /// owner's fixed-size local workspace.
    pub max_factor_entries: usize,
    /// Positive spring scale [N/m] for the membrane coordinate factorization.
    /// Changing it leaves physical stiffness unchanged. Reported thermal
    /// connection forces are virtual coordinates, not attachment reactions;
    /// the native connection-force ceiling applies at this declared scale.
    pub reference_stiffness_n_m: f64,
}

/// Original mechanical source and its full linear thermal realization.
///
/// The isothermal pencil, mass normalization, mode order, force footprints
/// and both volume maps remain unchanged. Independent Kelvin and interface
/// losses are retained once; they must exclude this explicitly selected
/// thermal contribution. Membrane storage has no dissipative thermal pole.
///
/// This local flat-slab constitutive approximation assumes a homogeneous
/// symmetric plane-stress substrate, a held source state, insulated faces and
/// through-thickness Fourier conduction. It does not include lateral heat
/// conduction, coatings, evolving temperature, finite strain or a thermal
/// initial field. Initial thermal and mechanical perturbations are zero.
/// The original elastic director condensation remains the spatial trial space.
///
/// Each spectrum's bound covers its normalized harmonic constitutive law on
/// its reported source-intersected band. It does not certify spatial/modal
/// truncation, projection roundoff, an arbitrary forced trajectory or the
/// frequencies produced after adding another physical coupling.
#[derive(Clone)]
pub struct ThermoelasticShellImage {
    original: MaterialShellImage,
    properties: Vec<ShellThermalProperties>,
    spectra: Vec<SlabPronySpectrum>,
    membrane: Vec<Vec<f64>>,
    branches: Vec<RelaxationBranch>,
    reference_stiffness_n_m: f64,
}

impl MaterialShellImage {
    /// Attach explicit thermal roles from these SAME original regional sources.
    ///
    /// All selected physical modes and off-diagonal thermal terms are retained.
    /// Regions with identical finite spectra share one row compression, even
    /// when their expansion, orientation or elastic stiffness differ.
    /// Distinct thermal rates never share a memory coordinate.
    ///
    /// # Errors
    /// Missing/unheld source data, incompatible substrate/loss policy, invalid
    /// resources, exhausted spectral/projection/state allowance, cancellation,
    /// or represented instantaneous stiffness outside original source support.
    pub fn with_thermoelastic_memory(
        &self,
        properties: &[ShellThermalProperties],
        options: ShellThermoelasticOptions,
        gate: &CancelGate,
    ) -> Result<ThermoelasticShellImage, AcousticRealizeError> {
        poll(gate)?;
        let n = self.physical_modes().len();
        if n == 0 || n > super::MAX_SHELL_MODES
            || options.max_memory_states > 8192 || options.max_projection_terms == 0
            || options.max_factor_entries == 0
            || !options.reference_stiffness_n_m.is_finite()
            || options.reference_stiffness_n_m <= 0.0
        {
            return Err(invalid("shell thermal memory requires finite bounded original coordinates and resources"));
        }
        let bound = source::bind_regions(self, properties, options.band_rad_s, gate)?;
        let mut spectra = Vec::with_capacity(bound.len());
        for region in &bound {
            poll(gate)?;
            spectra.push(region.law.lower_to_prony(
                region.source_band_rad_s, options.normalized_tolerance, options.max_thermal_terms,
            ).map_err(owner)?);
        }
        // Only identical represented kernels may mix spatial thermal rows.
        // Keep every regional source and its own applicability/error report.
        let mut group_spectrum: Vec<usize> = Vec::new();
        let mut region_group = Vec::with_capacity(spectra.len());
        let mut kernels = std::collections::BTreeMap::new();
        for (region, spectrum) in spectra.iter().enumerate() {
            poll(gate)?;
            let key: Vec<_> = spectrum.terms().iter().map(|term|
                (term.weight().to_bits(), term.relaxation_time_s().to_bits())).collect();
            let group = match kernels.get(&key) {
                Some(&group) => group,
                None => {
                    let group = group_spectrum.len();
                    kernels.insert(key, group);
                    group_spectrum.push(region);
                    group
                }
            };
            region_group.push(group);
        }
        let material = self.source();
        let elements = material.mesh.tris.len();
        let count = elements.checked_mul(7)
            .ok_or_else(|| invalid("shell thermal quadrature extent overflows"))?.max(n);
        let entries = count.checked_mul(n)
            .ok_or_else(|| invalid("shell thermal factor extent overflows"))?;
        let scratch_entries = entries.checked_mul(4)
            .and_then(|v| n.checked_mul(n)?.checked_mul(4)?.checked_add(v))
            .and_then(|v| n.checked_mul(20)?.checked_add(v))
            .and_then(|v| elements.checked_mul(6)?.checked_add(v))
            .ok_or_else(|| invalid("shell thermal workspace extent overflows"))?;
        if scratch_entries > options.max_factor_entries {
            return Err(invalid("shell thermal quadrature and QR exceed their storage allowance"));
        }
        let mut group_facets = vec![0usize; group_spectrum.len()];
        for &region in &material.element_regions {
            group_facets[region_group[region]] += 1;
        }
        let qr_work = |rows: usize| -> Option<usize> {
            rows.checked_mul(n)?.checked_mul(n.checked_mul(8)?.checked_add(64)?)
        };
        let mut work = projection_terms(material.mesh.nodes.len(), elements, elements, n)
            .and_then(|v| qr_work(count)?.checked_add(v))
            .ok_or_else(|| invalid("shell thermal preparation work overflows"))?;
        for facets in group_facets {
            let rows = facets.checked_mul(7)
                .ok_or_else(|| invalid("regional thermal quadrature extent overflows"))?.max(n);
            work = work.checked_add(qr_work(rows)
                .ok_or_else(|| invalid("regional thermal QR work overflows"))?)
                .ok_or_else(|| invalid("shell thermal preparation work overflows"))?;
        }
        if work > options.max_projection_terms {
            return Err(invalid("complete shell thermal projection exceeds its work allowance"));
        }
        let mut factors = Vec::with_capacity(elements);
        for (element, &region) in material.element_regions.iter().enumerate() {
            poll(gate)?;
            let input = material.regions[region].input();
            let angle = fs_plate::shell::survey::material_axis_angle(
                &material.mesh.facet(element).map_err(owner)?, input.axis1_global,
            ).map_err(owner)? + input.region.material_angle_rad;
            let bending = rotate_factor(bound[region].law.relaxation_factor(), angle)?;
            let ratio = det::sqrt(12.0) / material.regions[region].thickness_m();
            let membrane = bending.map(|value| ratio * value);
            if !ratio.is_finite() || membrane.iter().any(|v| !v.is_finite())
                || membrane.iter().zip(bending).any(|(a, d)| *a == 0.0 && d != 0.0)
            {
                return Err(invalid("shell thermal membrane factor scaling is unrepresentable"));
            }
            factors.push(Some(ShellThermoelasticFactors { membrane, bending }));
        }
        let mut projected = project_thermoelastic_factors(
            &material.mesh, self.directors(), &material.sections, &factors,
            self.physical_modes(), material.budget, options.max_projection_terms,
            options.max_factor_entries, || gate.is_requested(),
        ).map_err(owner)?;
        let mut membrane = compress_rows(&mut projected.membrane, projected.row_count, n, gate)?;
        // Drop the original membrane quadrature before regional QR allocation.
        drop(projected.membrane);
        let scale = det::sqrt(options.reference_stiffness_n_m);
        for column in &mut membrane {
            for value in column {
                let original = *value;
                *value /= scale;
                if !value.is_finite() || (*value == 0.0 && original != 0.0) {
                    return Err(invalid("shell thermal membrane coordinate is unrepresentable"));
                }
            }
        }
        let mut group_rows = vec![Vec::new(); group_spectrum.len()];
        for (row, element) in projected.row_elements.iter().enumerate() {
            if let Some(element) = element {
                group_rows[region_group[material.element_regions[*element]]].push(row);
            }
        }
        let mut branches = Vec::new();
        for (group, indices) in group_rows.iter().enumerate() {
            poll(gate)?;
            let count = indices.len().max(n);
            let mut rows = vec![0.0; count * n];
            for (to, &from) in indices.iter().enumerate() {
                rows[to * n..(to + 1) * n]
                    .copy_from_slice(&projected.bending[from * n..(from + 1) * n]);
            }
            let columns = compress_rows(&mut rows, count, n, gate)?;
            let spectrum = &spectra[group_spectrum[group]];
            let extra = columns.len().checked_mul(spectrum.terms().len())
                .ok_or_else(|| invalid("shell thermal memory extent overflows"))?;
            let states = branches.len().checked_add(extra)
                .ok_or_else(|| invalid("shell thermal memory extent overflows"))?;
            let stored = states.checked_mul(2 * n)
                .and_then(|v| v.checked_add(scratch_entries));
            work = extra.checked_mul(2 * n).and_then(|v| v.checked_add(work))
                .ok_or_else(|| invalid("shell thermal memory expansion work overflows"))?;
            if states > options.max_memory_states
                || stored.is_none_or(|v| v > options.max_factor_entries)
                || work > options.max_projection_terms
            {
                return Err(invalid("complete shell thermal memory exceeds its state, storage or work allowance"));
            }
            for term in spectrum.terms() {
                for column in &columns {
                    poll(gate)?;
                    let mut projection = vec![0.0; 2 * n];
                    for (i, &value) in column.iter().enumerate() {
                        projection[2 * i] = value;
                    }
                    branches.push(RelaxationBranch {
                        projection, stiffness: term.weight(),
                        relaxation_time_s: term.relaxation_time_s(),
                    });
                }
            }
        }
        let omega: Vec<_> = self.radiators().iter().map(|r| r.omega).collect();
        let mut band = options.band_rad_s;
        for support in bound.iter().map(|r| r.source_band_rad_s)
            .chain(material.interfaces.iter().map(|r| r.omega_band_rad_s()))
        {
            band = (band.0.max(support.0), band.1.min(support.1));
        }
        admit_band(&omega, &membrane, options.reference_stiffness_n_m, &branches, band)?;
        poll(gate)?;
        Ok(ThermoelasticShellImage {
            original: self.clone(), properties: properties.to_vec(), spectra,
            membrane, branches, reference_stiffness_n_m: options.reference_stiffness_n_m,
        })
    }
}

impl ThermoelasticShellImage {
    /// Original isothermal image, including unchanged force/volume/source maps.
    #[must_use]
    pub const fn original(&self) -> &MaterialShellImage { &self.original }

    /// Explicit thermal property selections, in original material-region order.
    #[must_use]
    pub fn properties(&self) -> &[ShellThermalProperties] { &self.properties }

    /// Regional finite spectra with their own source bands and error bounds.
    #[must_use]
    pub fn spectra(&self) -> &[SlabPronySpectrum] { &self.spectra }

    /// Complete passive bending arms, interleaved [q0,0,q1,0,...].
    #[must_use]
    pub fn relaxation_branches(&self) -> &[RelaxationBranch] { &self.branches }

    /// Conservative membrane columns in the unchanged mass-normalized basis.
    /// Their stiffness is reference_stiffness_n_m times the sum of outer products.
    #[must_use]
    pub fn membrane_columns(&self) -> &[Vec<f64>] { &self.membrane }

    #[must_use]
    pub const fn reference_stiffness_n_m(&self) -> f64 { self.reference_stiffness_n_m }

    /// Fill generalized force from the original physical port forces, without
    /// allocation. Every signed port contribution is retained. Invalid size,
    /// nonfinite input or an unrepresentable contraction leaves output unchanged.
    pub fn fill_force(
        &self, forces_n: &[f64], output: &mut [f64],
    ) -> Result<(), AcousticRealizeError> {
        let ports = self.original.force_ports();
        let n = self.original.physical_modes().len();
        if forces_n.len() != ports.len() || output.len() != n
            || forces_n.iter().any(|v| !v.is_finite())
        {
            return Err(invalid("shell thermal drive requires every original finite port force"));
        }
        let force = |i| compensated_sum(ports.iter().zip(forces_n).map(|(p, f)| p[i] * f));
        if (0..n).any(|i| !force(i).is_finite()) {
            return Err(invalid("shell thermal drive projection is unrepresentable"));
        }
        for (i, out) in output.iter_mut().enumerate() { *out = force(i); }
        Ok(())
    }

    /// Prepare the existing complete held-force owner at zero perturbation.
    ///
    /// Original material/interface Kelvin loss, reversible thermal membrane
    /// storage and every bending memory coordinate share one coupled operator.
    /// All original setup, instantaneous-stiffness Nyquist, work and heat gates
    /// remain active. Thermal diffusion poles impose no artificial dt/tau limit.
    /// Pressure observation is zero until an explicit acoustic owner is bound.
    /// The original image exposes both reciprocal volume maps for that purpose.
    pub fn prepare_coupled(
        &self,
        sample_rate_hz: u32,
        modal_budget: ModalAcousticTimeBudget,
        coupling: ModalCouplingConfig,
        preparation: PreparedModalBudget,
        gate: &CancelGate,
    ) -> Result<CoupledModalSystem, AcousticRealizeError> {
        poll(gate)?;
        let n = self.original.physical_modes().len();
        let material = self.original.source();
        let interfaces = material.interfaces.iter()
            .filter(|source| source.law().damping_n_s_m != 0.0).count();
        let channels = self.original.damping().channel_count()
            .checked_add(interfaces).and_then(|v| v.checked_add(self.membrane.len()))
            .ok_or_else(|| invalid("shell thermal connection extent overflows"))?;
        let terms = channels.checked_add(1)
            .and_then(|v| n.checked_mul(v)?.checked_mul(v));
        preparation.validate().map_err(owner)?;
        // Admit the FINAL augmented extent before constructing a memoryless
        // exponential that would certainly be discarded and then refused.
        let dimension = n.checked_mul(3)
            .and_then(|v| v.checked_add(self.branches.len()))
            .and_then(|v| v.checked_mul(2));
        if dimension.is_none_or(|v| v > fs_phs::linear_exact::MAX_PREPARED_EXPONENTIAL_DIM) {
            return Err(invalid("shell thermal memory exceeds the original prepared exponential capacity"));
        }
        if n > coupling.max_modes || n > preparation.max_modes
            || coupling.max_modes > fs_modal::time::MAX_TIME_DOMAIN_ACOUSTIC_MODES
            || channels > coupling.max_connections || coupling.max_connections > MAX_MODAL_CONNECTIONS
            || terms.is_none_or(|v| v > coupling.max_setup_terms)
        {
            return Err(invalid("shell thermal runtime exceeds original mode, connection or work limits"));
        }
        let model = ModalAcousticTimeModel::try_new(
            sample_rate_hz, self.original.radiators().iter().map(|r| ModalAcousticMode {
                angular_frequency_rad_s: r.omega,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::ZERO,
            }).collect(), modal_budget,
        ).map_err(owner)?;
        let mut connections = Vec::with_capacity(channels);
        for column in self.original.damping().channel_shapes() {
            connections.push(ModalColumnConnection {
                column: column.clone(), stiffness_n_m: 0.0,
                damping_n_s_m: self.original.damping().reference_resistance_n_s_m(),
                rest_extension_m: 0.0,
            });
        }
        for (source, column) in material.interfaces.iter().zip(self.original.interface_columns()) {
            if source.law().damping_n_s_m != 0.0 {
                connections.push(ModalColumnConnection {
                    column: column.clone(), stiffness_n_m: 0.0,
                    damping_n_s_m: source.law().damping_n_s_m, rest_extension_m: 0.0,
                });
            }
        }
        for column in &self.membrane {
            connections.push(ModalColumnConnection {
                column: column.clone(), stiffness_n_m: self.reference_stiffness_n_m,
                damping_n_s_m: 0.0, rest_extension_m: 0.0,
            });
        }
        poll(gate)?;
        let system = CoupledModalSystem::new_with_columns_preparation(
            vec![model], connections, &[], coupling, ModalPreparation::Prepared(preparation), gate,
        ).map_err(owner)?;
        if self.branches.is_empty() { return Ok(system); }
        system.with_material_relaxation(self.branches.clone(), None, gate).map_err(owner)
    }
}

// Stress-like covector transformation R^T f, matching the original section's
// engineering-strain rotation exactly. The source's selected axes are not
// inferred again from the geometry and the low-level projector does not rotate.
fn rotate_factor(factor: [f64; 3], angle: f64) -> Result<[f64; 3], AcousticRealizeError> {
    let (s, c) = (det::sin(angle), det::cos(angle));
    let t = [
        [c * c, s * s, c * s],
        [s * s, c * c, -c * s],
        [-2.0 * c * s, 2.0 * c * s, c * c - s * s],
    ];
    let rotated = core::array::from_fn(|i|
        compensated_sum((0..3).map(|j| t[j][i] * factor[j])));
    if rotated.iter().any(|v| !v.is_finite()) {
        return Err(invalid("shell thermal material-frame factor is unrepresentable"));
    }
    Ok(rotated)
}

fn compensated_sum(values: impl Iterator<Item = f64>) -> f64 {
    let mut sum = 0.0_f64;
    let mut correction = 0.0;
    for value in values {
        let next = sum + value;
        correction += if sum.abs() >= value.abs() {
            (sum - next) + value
        } else {
            (value - next) + sum
        };
        sum = next;
    }
    sum + correction
}

// lambda_max(diag(omega²) + sum k b b^T) <= max(omega²) + sum k ||b||².
// Positive operations are rounded outward. This bounds only represented
// instantaneous stiffness; it does not turn arbitrary forcing into a bandlimit.
fn admit_band(
    omega: &[f64], membrane: &[Vec<f64>], stiffness: f64,
    branches: &[RelaxationBranch], band: (f64, f64),
) -> Result<(), AcousticRealizeError> {
    if band.0 >= band.1
        || omega.iter().any(|w| !w.is_finite() || *w < band.0 || *w > band.1)
    {
        return Err(invalid("shell thermal source band excludes an original retained mode"));
    }
    let mut upper = omega.iter().map(|w| (w * w).next_up()).fold(0.0, f64::max);
    let mut add = |k: f64, values: &[f64]| {
        let norm = values.iter().filter(|v| **v != 0.0)
            .fold(0.0_f64, |sum, v| (sum + (v * v).next_up()).next_up());
        if norm != 0.0 { upper = (upper + (k * norm).next_up()).next_up(); }
    };
    for column in membrane { add(stiffness, column); }
    for branch in branches { add(branch.stiffness, &branch.projection); }
    if !upper.is_finite() || upper > (band.1 * band.1).next_down() {
        return Err(invalid("shell thermal stiffening exceeds the original constitutive band"));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
