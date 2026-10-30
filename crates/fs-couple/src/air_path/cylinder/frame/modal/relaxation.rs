//! Original layer E(t) on the SAME complete frame/fluid mechanical timeline.
//!
//! This adapter uses the existing source spectrum resolver and native Maxwell
//! internal-strain owner. Three geometry-owned strain factors integrate each
//! longitudinal modulus increment: delta_E*A/L, 3*delta_E*I/L^3, delta_E*I/L^3.
//! They share the original loaded-inertia congruence with physical forcing.
//! Neither a fitted modal damping ratio nor a fictitious internal mass is used.

/// Source-pinned temperature changes of relaxation rates at retained history.
pub mod temperature;

use super::{FrameFluidCoordinates, FrameFluidPlacement, MaterialFluidFrame};
use crate::air_path::cylinder::frame::FrameFluidError;
use crate::render::schedule::force::coupled::{CoupledModalSystem, ModalCouplingError};
use fs_exec::CancelGate;
use fs_la::factor::cholesky;
use fs_material::section::relaxation::{ConcentricPronySection, LayerRelaxationProperties};
use fs_material::state_point::MaterialStatePointError;
use fs_phs::RelaxationBranch;
use core::ops::Range;
use std::collections::BTreeSet;

/// One scalar longitudinal E(t) spectrum per original layer, used for BOTH
/// axial strain and bending of this prismatic isotropic-section member.
/// This is not an anisotropic bending-only loss law or a torsional spectrum.
#[derive(Clone, Debug)]
pub struct FrameLongitudinalRelaxation {
    pub element: usize,
    pub layers: Vec<LayerRelaxationProperties>,
}

/// Additional cold work bounds. Native COMPLETE assembly/memory/energy and
/// instantaneous-stiffness Nyquist admissions remain independently binding.
#[derive(Clone, Copy, Debug)]
pub struct FrameRelaxationBudget {
    pub max_members: usize,
    /// All declared source pairs, including explicit zero modulus increments.
    pub max_source_arms: usize,
    /// Three strain coordinates per positive source arm. Admission conservatively
    /// counts three per declared arm BEFORE evaluating or allocating source data.
    pub max_memory_coordinates: usize,
    /// Algebraic work scale, including original-state conversion, equilibrium
    /// port admission, one Cholesky chart and complete native projections.
    /// Native solver preparation and source-law evaluation retain their OWN caps.
    pub max_setup_terms: usize,
}

#[derive(Debug)]
pub enum FrameRelaxationError {
    Invalid(&'static str),
    Frame(FrameFluidError),
    Source(MaterialStatePointError),
    Native(ModalCouplingError),
    Cancelled,
}
impl core::fmt::Display for FrameRelaxationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source frame relaxation: {self:?}")
    }
}
impl std::error::Error for FrameRelaxationError {}
impl From<FrameFluidError> for FrameRelaxationError {
    fn from(error: FrameFluidError) -> Self {
        match error { FrameFluidError::Cancelled => Self::Cancelled, other => Self::Frame(other) }
    }
}
impl From<ModalCouplingError> for FrameRelaxationError {
    fn from(error: ModalCouplingError) -> Self {
        match error { ModalCouplingError::Cancelled => Self::Cancelled, other => Self::Native(other) }
    }
}
impl From<MaterialStatePointError> for FrameRelaxationError {
    fn from(error: MaterialStatePointError) -> Self { Self::Source(error) }
}
type Result<T> = core::result::Result<T, FrameRelaxationError>;
fn invalid(message: &'static str) -> FrameRelaxationError { FrameRelaxationError::Invalid(message) }
fn poll(gate: &CancelGate) -> Result<()> {
    if gate.is_requested() { Err(FrameRelaxationError::Cancelled) } else { Ok(()) }
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() { Ok(value) } else { Err(invalid("frame relaxation arithmetic is unrepresentable")) }
}
fn positive(value: f64) -> Result<f64> {
    finite(value)?;
    if value > 0.0 { Ok(value) } else { Err(invalid("positive frame relaxation coefficient required")) }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FrameStrainChannel { Axial, SymmetricBending, AntisymmetricBending }

/// Identity of one actual native constitutive coordinate; source-arm indices
/// include zero increments even though an exact zero creates no dynamic state.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FrameRelaxationCoordinate {
    pub element: usize,
    pub layer: usize,
    pub source_arm: usize,
    pub channel: FrameStrainChannel,
    pub native_memory: usize,
}

/// Same-frame source material and its exact equilibrium section, no anonymous
/// effective material. Its layer selectors, receipts and original tau survive.
#[derive(Clone, Debug)]
pub struct FrameRelaxationMember {
    pub element: usize,
    pub source: ConcentricPronySection,
}

/// Source identity accompanies native evolving state; this object does NOT
/// integrate or own a second copy of constitutive history.
pub struct FrameRelaxationBinding<'m, 'f, F> {
    source: &'m MaterialFluidFrame<'f, F>,
    members: Vec<FrameRelaxationMember>,
    coordinates: Vec<FrameRelaxationCoordinate>,
    branches: Vec<RelaxationBranch>,
    memory: Range<usize>,
    native_modes: usize,
}

#[derive(Clone, Copy, Debug)]
pub struct FrameRelaxationState {
    pub coordinate: FrameRelaxationCoordinate,
    /// Original integrated strain coordinate (metres), not a point fibre strain.
    pub extension_m: f64,
    /// Viscous extension z/sqrt(k) of the SAME native Maxwell arm [m].
    pub viscous_extension_m: f64,
    pub stored_energy_j: f64,
    /// Instantaneous physical loss rate, NOT the integrated interval heat.
    pub dissipated_power_w: f64,
}

impl<'m, 'f, F> FrameFluidPlacement<'m, 'f, F> {
    /// Attach all selected source-member longitudinal spectra before playback.
    /// The original frame supplies EQUILIBRIUM stiffness; arms add only the
    /// relaxing increments. No second equilibrium spring or bulk viscosity is
    /// inferred. Nonzero Kelvin-Voigt loss on a selected member refuses, rather
    /// than silently counting an unspecified second material loss mechanism.
    ///
    /// `None` explicitly starts these arms fully relaxed at current motion.
    /// `Some(z)` supplies energy-normalized internal strains in returned
    /// coordinate order. Zero z on a displaced body is an unrelaxed loading.
    /// Internal strains outside the represented axial/cubic-bending subspace
    /// are not inferred. Zero-projection positive arms are KEPT, including any
    /// initial stored energy. Only exactly zero source increments have no state.
    ///
    /// All existing fluid momenta and foreign material memories remain intact.
    /// Declare this frame's selected members together: an already participating
    /// Maxwell projection refuses a second frame binding. Distinct frames may
    /// append their own spectra through the SAME native augmented solver.
    ///
    /// Fixed temperature, reference density/geometry, perfect layer bonding and
    /// small longitudinal strain apply. This is not a live thermoelastic Prony
    /// update or automatic heat deposition. Use the returned native owner for
    /// BOTH transient and harmonic response, not the source's frozen pHS image.
    /// The consuming native constructor returns no partial owner on refusal.
    pub fn with_longitudinal_relaxation(
        &self,
        mechanical: CoupledModalSystem,
        selections: &[FrameLongitudinalRelaxation],
        initial_memory: Option<Vec<f64>>,
        budget: FrameRelaxationBudget,
        gate: &CancelGate,
    ) -> Result<(FrameRelaxationBinding<'m, 'f, F>, CoupledModalSystem)> {
        poll(gate)?;
        let source = self.source();
        let m = source.basis().modes().len();
        let total = mechanical.mode_count();
        if mechanical.samples_rendered() != 0 || selections.is_empty()
            || selections.len() > budget.max_members
            || source.basis().reference_displacement().iter().any(|&u| u != 0.0) {
            return Err(invalid("source relaxation requires bounded cold member bindings and zero frozen preload"));
        }
        let mut arms = 0usize;
        let mut unique = BTreeSet::new();
        for selection in selections {
            poll(gate)?;
            if !unique.insert(selection.element) { return Err(invalid("member relaxation was declared twice")); }
            for layer in &selection.layers {
                arms = arms.checked_add(layer.branches.len())
                    .ok_or_else(|| invalid("source relaxation count overflow"))?;
            }
        }
        let memory_cap = arms.checked_mul(3).ok_or_else(|| invalid("strain memory count overflow"))?;
        let extent = total.checked_mul(2).ok_or_else(|| invalid("native projection extent overflow"))?;
        let work = m.checked_mul(m).and_then(|square| {
            square.checked_mul(3).and_then(|n| n.checked_add(m.checked_mul(18)?))
                .and_then(|n| n.checked_mul(selections.len()))
                .and_then(|n| n.checked_add(square.checked_mul(m)?))
                .and_then(|n| n.checked_add(extent.checked_mul(memory_cap)?))
                .and_then(|n| n.checked_add(source.system().state_dim()
                    .checked_mul(source.system().state_dim())?.checked_mul(2)?))
                .and_then(|n| n.checked_add(m.checked_mul(total)?))
                .and_then(|n| n.checked_add(m.checked_mul(mechanical.material_memory().len())?))
        });
        if arms > budget.max_source_arms || memory_cap > budget.max_memory_coordinates
            || work.is_none_or(|n| n > budget.max_setup_terms) {
            return Err(invalid("complete source relaxation exceeds declared memory/setup budget"));
        }
        // Check the retained physical chart before using native component indices.
        self.physical_state(mechanical.components(), gate)?;
        let component = self.structural_component();
        let offset: usize = mechanical.components()[..component].iter().map(|c| c.modes().len()).sum();
        if mechanical.material_relaxation_branches().iter().any(|branch|
            branch.projection.iter().step_by(2).skip(offset).take(m).any(|&b| b != 0.0)) {
            return Err(invalid("declare this frame's material memories together before any participating relaxation binding"));
        }
        let factor = cholesky(source.modal_mass_matrix(), m)
            .map_err(|error| FrameFluidError::Factor(error))?;
        let chart = FrameFluidCoordinates { source, factor };
        admit_equilibrium(self, &chart, &mechanical, offset, gate)?;
        let mut members = Vec::with_capacity(selections.len());
        let mut branches = Vec::with_capacity(memory_cap);
        let mut coordinates = Vec::with_capacity(memory_cap);
        let start = mechanical.material_memory().len();
        let nyquist = positive(core::f64::consts::PI / mechanical.sample_period_s())?;
        for selection in selections {
            poll(gate)?;
            let frame = source.frame();
            if frame.material_viscous_section(selection.element).is_some_and(|section|
                section.axial().viscous_rigidity_n_s() != 0.0
                    || section.bending().viscous_stiffness_n_m2_s() != 0.0) {
                return Err(invalid("selected Prony members must not already carry nonzero Kelvin-Voigt material loss"));
            }
            let section = frame.material_section(selection.element)
                .ok_or_else(|| invalid("relaxation requires the ORIGINAL source-bound member section"))?
                .clone().with_prony_bending(&selection.layers)?;
            let band = section.omega_band_rad_s();
            if band.0 != 0.0 || band.1 < nyquist {
                return Err(invalid("longitudinal source spectrum must cover DC through the actual native Nyquist band"));
            }
            let (dofs, weights, length) = frame.reference_strain_traces(selection.element)
                .map_err(FrameFluidError::from)?;
            // Compute each physical map ONCE, not once per layer or source arm.
            let mut projections = Vec::with_capacity(3);
            for column in weights {
                poll(gate)?;
                let dry = source.basis().modes().iter().map(|mode|
                    column.iter().zip(dofs).try_fold(0.0, |sum, (&b, dof)|
                        finite(sum + b * mode.displacement[dof]))).collect::<Result<Vec<_>>>()?;
                let transformed = chart.project_column(&dry)?;
                let mut full = vec![0.0; extent];
                for (i, &b) in transformed.iter().enumerate() { full[2 * (offset + i)] = b; }
                projections.push(full);
            }
            for (layer, (geometry, law)) in section.elastic().layer_integrals().iter()
                .zip(section.spectra()).enumerate() {
                for (arm, &(delta, tau)) in law.terms.iter().enumerate() {
                    poll(gate)?;
                    if delta == 0.0 { continue; }
                    let axial = positive(delta * geometry.area_m2() / length)?;
                    let bending = positive(delta * geometry.second_moment_m4() / length / length / length)?;
                    for (channel, stiffness, projection) in [
                        (FrameStrainChannel::Axial, axial, &projections[0]),
                        (FrameStrainChannel::SymmetricBending, positive(3.0 * bending)?, &projections[1]),
                        (FrameStrainChannel::AntisymmetricBending, bending, &projections[2]),
                    ] {
                        coordinates.push(FrameRelaxationCoordinate { element: selection.element,
                            layer, source_arm: arm, channel, native_memory: start.checked_add(branches.len())
                                .ok_or_else(|| invalid("native memory address overflow"))? });
                        branches.push(RelaxationBranch { projection: projection.to_vec(), stiffness,
                            relaxation_time_s: tau });
                    }
                }
            }
            members.push(FrameRelaxationMember { element: selection.element, source: section });
        }
        let end = start.checked_add(branches.len()).ok_or_else(|| invalid("native memory range overflow"))?;
        poll(gate)?;
        // Existing native owner controls all causal evolution, instantaneous
        // stiffness admission, work/storage/loss and prepared exponential reuse.
        let mechanical = mechanical.with_additional_material_relaxation(branches.clone(), initial_memory, gate)?;
        poll(gate)?;
        Ok((FrameRelaxationBinding { source, members, coordinates, branches,
            memory: start..end, native_modes: total }, mechanical))
    }
}

// Do not let source E_infinity be counted twice or replace the declared dry
// equilibrium by arbitrary modified springs. Foreign assembly links remain valid.
fn admit_equilibrium<F>(placement: &FrameFluidPlacement<'_, '_, F>,
    chart: &FrameFluidCoordinates<'_, '_, F>, mechanical: &CoupledModalSystem,
    offset: usize, gate: &CancelGate) -> Result<()> {
    let source = placement.source();
    let count = source.basis().modes().len();
    let range = &placement.connections().elastic;
    if range.len() != count || mechanical.constant_loss_connections().iter()
        .any(|fit| range.contains(&fit.connection)) {
        return Err(invalid("source frame requires its original equilibrium elastic ports"));
    }
    let mut unit = vec![0.0; count];
    for (i, mode) in source.basis().modes().iter().enumerate() {
        poll(gate)?;
        let index = range.start + i;
        let law = mechanical.connection_coefficients(index)
            .ok_or_else(|| invalid("source equilibrium spring is absent"))?;
        let column = mechanical.connection_column(index)
            .ok_or_else(|| invalid("source equilibrium column is absent"))?;
        unit[i] = 1.0;
        let expected = chart.project_column(&unit)?;
        unit[i] = 0.0;
        if law.stiffness_n_m.to_bits() != mode.eigenvalue_per_s2.to_bits()
            || law.damping_n_s_m != 0.0 || law.rest_extension_m != 0.0
            || column.len() != mechanical.mode_count()
            || column.iter().enumerate().any(|(j, &b)| b != if (offset..offset + count).contains(&j) {
                expected[j - offset]
            } else { 0.0 }) {
            return Err(invalid("original source equilibrium stiffness or physical coordinate map changed"));
        }
    }
    Ok(())
}

impl<'m, 'f, F> FrameRelaxationBinding<'m, 'f, F> {
    pub const fn source(&self) -> &'m MaterialFluidFrame<'f, F> { self.source }
    pub fn members(&self) -> &[FrameRelaxationMember] { &self.members }
    pub fn coordinates(&self) -> &[FrameRelaxationCoordinate] { &self.coordinates }
    pub fn memory_range(&self) -> Range<usize> { self.memory.clone() }

    /// Read the ACTUAL accepted native histories without advancing mechanics,
    /// recomputing source coefficients or changing the physical clock. Frozen
    /// source identity and branch definitions are checked; a same-sized foreign
    /// branch block is not accepted merely because its dimensions agree.
    pub fn states(&self, mechanical: &CoupledModalSystem, gate: &CancelGate)
        -> Result<Vec<FrameRelaxationState>> {
        poll(gate)?;
        if mechanical.mode_count() != self.native_modes {
            return Err(invalid("complete mechanical layout changed after source memory binding"));
        }
        let native = mechanical.material_relaxation_branches().get(self.memory.clone())
            .ok_or_else(|| invalid("source material branch block is absent"))?;
        let histories = mechanical.material_memory().get(self.memory.clone())
            .ok_or_else(|| invalid("source internal strains are absent"))?;
        let mut states = Vec::with_capacity(self.coordinates.len());
        for (((coordinate, expected), branch), &z) in self.coordinates.iter().zip(&self.branches)
            .zip(native).zip(histories) {
            poll(gate)?;
            if branch.stiffness.to_bits() != expected.stiffness.to_bits()
                || branch.relaxation_time_s.to_bits() != expected.relaxation_time_s.to_bits()
                || branch.projection != expected.projection {
                return Err(invalid("source branch law or strain projection changed"));
            }
            let extension = branch.projection.iter().step_by(2)
                .zip(mechanical.components().iter().flat_map(|c| c.states()))
                .try_fold(0.0, |sum, (&b, q)| finite(sum + b * q.displacement_m_sqrt_kg))?;
            let root = positive(fs_math::det::sqrt(branch.stiffness))?;
            let mismatch = finite(root * extension - z)?;
            let energy = finite(0.5 * mismatch * mismatch)?;
            states.push(FrameRelaxationState { coordinate: *coordinate, extension_m: extension,
                viscous_extension_m: finite(z / root)?, stored_energy_j: energy,
                dissipated_power_w: finite(mismatch * mismatch / branch.relaxation_time_s)? });
        }
        poll(gate)?;
        Ok(states)
    }
}

#[cfg(test)]
mod tests;
