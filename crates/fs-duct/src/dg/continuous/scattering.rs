//! Incident and reflected pressure waves on the complete original DG state.
//!
//! Source loading precedes LU: (s I - L + B Z^-1 B^T)x = 2 B Z^-1 a.
//! No bare duct impedance is formed. Original pressure/flow equations and
//! material losses are observed afterwards, without changing the source owner.

use super::{
    DgContinuous, DgHarmonicResponse, admit_tolerance, finite, invalid, lu_error, poll,
};
use crate::DuctError;
use crate::dg::{DgDuct, DgMedium, MAX_DG_PROPAGATOR_DOFS, PortId};
use fs_blake3::ContentHash;
use fs_la::eigen_complex::{LuComplex, lu_complex_with_cancel};
use fs_material::gas::GasState;
use fs_math::{c64::C64, det};

/// Peak pressure waves at one original port, in the declared reference order.
#[derive(Clone, Debug, PartialEq)]
pub struct DgScatteringPort {
    pub id: PortId,
    /// Positive real pressure/volume-flow reference [Pa s/m3].
    pub reference_impedance_pa_s_m3: f64,
    pub incident_pressure_pa: C64,
    pub reflected_pressure_pa: C64,
    /// Independent squared incident and reflected wave amplitudes [W].
    pub incident_power_w: f64,
    pub reflected_power_w: f64,
    /// Signed original Re(p conj(Q))/2 [W], never a difference of wave powers.
    pub net_input_power_w: f64,
    /// Incident minus reflected minus physical input power [W].
    /// This is a rounding diagnostic, not an estimate of weak physical loss.
    pub power_residual_w: f64,
}

/// Coherent steady-state response with all original material coordinates.
#[derive(Clone, Debug, PartialEq)]
pub struct DgScatteringResponse {
    /// Original DG fields, pressure/flow ports, independent loss and storage.
    pub physical: DgHarmonicResponse,
    /// Selected ports in the order supplied when the factor was prepared.
    pub waves: Vec<DgScatteringPort>,
    /// Native loss minus the sum of incident-minus-reflected power [W].
    /// Diagnostic only; physical.power_residual_w uses pressure/flow work and
    /// the independently evaluated original material/boundary loss owners.
    pub wave_power_residual_w: f64,
}

/// One loaded full-state LU, borrowed from one immutable original DG owner.
///
/// Repeated responses share this factor. No history, material, port basis or
/// physical coordinate is reduced or modified. The numerical residual remains
/// a backward-error test, not a conditioning or forward-error certificate.
pub struct DgScatteringFactor<'a, M: DgMedium = GasState> {
    original: DgContinuous<'a, M>,
    omega_rad_s: f64,
    references: Vec<(PortId, f64)>,
    indices: Vec<usize>,
    roots: Vec<f64>,
    factor: LuComplex,
    scale: f64,
    maximum_work: usize,
    relative_tolerance: f64,
}

fn admit_work<M: DgMedium>(
    view: DgContinuous<'_, M>,
    selected: usize,
    maximum: usize,
) -> Result<(), DuctError> {
    let n = view.dimension();
    let terms = view.ports().len().checked_add(selected)
        .and_then(|v| v.checked_add(n))
        .and_then(|v| v.checked_add(16));
    let work = n.checked_mul(n)
        .and_then(|v| terms.and_then(|terms| v.checked_mul(terms)));
    if n == 0 || n > MAX_DG_PROPAGATOR_DOFS || work.is_none_or(|v| v > maximum) {
        return Err(invalid("complete DG scattering operation exceeds its explicit work budget"));
    }
    Ok(())
}

impl<M: DgMedium> DgDuct<M> {
    /// Prepare coherent incident-wave excitation at selected original ports.
    ///
    /// Each reference is positive, finite and real in Pa s/m3. For peak
    /// phasors under exp(-i omega t), Q=(2a-p)/Z flows INTO the original duct
    /// and b=p-a leaves it. Every unselected original port has exactly zero
    /// flow. Omitted incident amplitudes leave their selected reference loads
    /// connected, with a=0.
    ///
    /// The full generator receives B Z^-1 B^T before the one scaled LU.
    /// This avoids isolated impedance poles; genuinely singular or
    /// uncontrollable systems still refuse. All original material frequency
    /// bands and the complete-state cap apply. Preparation charges
    /// n^2*(n+original_ports+selected_ports+16) before allocation; each
    /// response uses the same allowance, not a cumulative lifetime budget.
    /// Cancellation publishes no partial factor or response.
    pub fn scattering_factor(
        &self,
        omega_rad_s: f64,
        references: &[(PortId, f64)],
        maximum_work: usize,
        relative_tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgScatteringFactor<'_, M>, DuctError> {
        poll(cancelled)?;
        let original = self.continuous_operator();
        original.admit_frequency(omega_rad_s)?;
        admit_tolerance(relative_tolerance)?;
        if references.is_empty() || references.len() > original.ports().len() {
            return Err(invalid("DG scattering needs a nonempty selection of original ports"));
        }
        admit_work(original, references.len(), maximum_work)?;
        let mut indices = Vec::with_capacity(references.len());
        let mut roots = Vec::with_capacity(references.len());
        for (i, &(id, reference)) in references.iter().enumerate() {
            poll(cancelled)?;
            if !reference.is_finite() || reference <= 0.0
                || references[..i].iter().any(|(old, _)| *old == id)
            {
                return Err(invalid("DG scattering references need distinct ports and positive finite impedances"));
            }
            let index = original.ports().iter().position(|(port, _)| *port == id)
                .ok_or_else(|| invalid("DG scattering reference names an undeclared port"))?;
            indices.push(index);
            roots.push(det::sqrt(reference));
        }
        let n = original.dimension();
        let mut matrix = vec![C64::ZERO; n * n];
        for i in 0..n {
            poll(cancelled)?;
            for j in 0..n {
                let value = original.generator_per_s()[i * n + j];
                if !value.is_finite() {
                    return Err(invalid("DG scattering original generator is not finite"));
                }
                matrix[i * n + j] = C64::new(
                    -value, if i == j { -omega_rad_s } else { 0.0 },
                );
            }
        }
        // B/sqrt(Z) avoids first forming 1/Z or B_i*B_j, either of which
        // can overflow even when the admitted rank-one load is representable.
        let mut column = vec![0.0; n];
        for (&index, &root) in indices.iter().zip(&roots) {
            poll(cancelled)?;
            for (value, &b) in column.iter_mut().zip(&original.ports()[index].1) {
                *value = b / root;
                if !value.is_finite() {
                    return Err(invalid("DG scattering loading column is not representable"));
                }
            }
            for i in 0..n {
                poll(cancelled)?;
                if column[i] == 0.0 { continue; }
                for j in 0..n {
                    let value = matrix[i * n + j].re + column[i] * column[j];
                    if !value.is_finite() {
                        return Err(invalid("DG scattering loaded matrix is not representable"));
                    }
                    matrix[i * n + j].re = value;
                }
            }
        }
        let scale = matrix.iter().fold(omega_rad_s, |scale, value|
            scale.max(value.re.abs()).max(value.im.abs()));
        for row in matrix.chunks_mut(n) {
            poll(cancelled)?;
            for value in row {
                value.re /= scale;
                value.im /= scale;
            }
        }
        let factor = lu_complex_with_cancel(&matrix, n, &mut *cancelled).map_err(lu_error)?;
        poll(cancelled)?;
        Ok(DgScatteringFactor {
            original, omega_rad_s, references: references.to_vec(), indices, roots,
            factor, scale, maximum_work, relative_tolerance,
        })
    }
}

impl<M: DgMedium> DgScatteringFactor<'_, M> {
    /// Exact declared order and original physical units; no implicit sorting.
    #[must_use]
    pub fn reference_impedances(&self) -> &[(PortId, f64)] {
        &self.references
    }

    /// Identity of the original duct, separately from its declared wave references.
    #[must_use]
    pub fn operator_identity(&self) -> ContentHash {
        self.original.operator_identity()
    }

    /// Reuse the full loaded LU for one coherent incident-pressure vector [Pa].
    ///
    /// Inputs may arrive in any order. Missing selected waves are zero;
    /// unknown, unselected, repeated or non-finite inputs refuse. The original
    /// observer receives the recovered physical flows and checks the original
    /// equations, storage and native losses without interpreting source
    /// reference loading as material heat.
    pub fn response(
        &self,
        incident: &[(PortId, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgScatteringResponse, DuctError> {
        poll(cancelled)?;
        if incident.len() > self.references.len() {
            return Err(invalid("too many DG scattering incident inputs"));
        }
        let mut amplitudes = vec![C64::ZERO; self.references.len()];
        for (i, &(id, value)) in incident.iter().enumerate() {
            poll(cancelled)?;
            if !finite(value) || incident[..i].iter().any(|(old, _)| *old == id) {
                return Err(invalid("DG incident waves need distinct ports and finite amplitudes"));
            }
            let index = self.references.iter().position(|(port, _)| *port == id)
                .ok_or_else(|| invalid("DG incident wave names an unselected port"))?;
            amplitudes[index] = value;
        }
        let normalized: Vec<_> = amplitudes.iter().zip(&self.roots)
            .map(|(a, &root)| C64::new(a.re / root, a.im / root)).collect();
        if normalized.iter().any(|&value| !finite(value)) {
            return Err(invalid("DG incident wave is not representable in its reference"));
        }
        let n = self.original.dimension();
        let mut state = vec![C64::ZERO; n];
        for i in 0..n {
            poll(cancelled)?;
            for (k, &index) in self.indices.iter().enumerate() {
                let b = self.original.ports()[index].1[i] / self.roots[k];
                state[i] = state[i] + normalized[k].scale((2.0 * b) / self.scale);
            }
            if !finite(state[i]) {
                return Err(invalid("DG scattering forcing is not representable"));
            }
        }
        self.factor.solve_with_cancel(&mut state, &mut *cancelled).map_err(lu_error)?;
        if state.iter().any(|&value| !finite(value)) {
            return Err(invalid("DG scattering solution is not finite"));
        }
        let mut flows = Vec::with_capacity(self.indices.len());
        for (k, &index) in self.indices.iter().enumerate() {
            poll(cancelled)?;
            let pressure = self.original.ports()[index].1.iter().zip(&state)
                .fold(C64::ZERO, |sum, (&b, &x)| sum + x.scale(b));
            let root = self.roots[k];
            let scaled_flow = normalized[k].scale(2.0)
                - C64::new(pressure.re / root, pressure.im / root);
            let flow = C64::new(scaled_flow.re / root, scaled_flow.im / root);
            if !finite(pressure) || !finite(flow) {
                return Err(invalid("DG scattering pressure or flow is not representable"));
            }
            flows.push((self.references[k].0, flow));
        }
        let physical = self.original.observe_harmonic(
            self.omega_rad_s, &state, &flows, self.maximum_work,
            self.relative_tolerance, cancelled,
        )?;
        let mut waves = Vec::with_capacity(self.references.len());
        let mut wave_supply = 0.0;
        for (k, &(id, reference)) in self.references.iter().enumerate() {
            poll(cancelled)?;
            let original = &physical.ports[self.indices[k]];
            let reflected = original.pressure_pa - amplitudes[k];
            let root = self.roots[k];
            let reflected_normalized = C64::new(reflected.re / root, reflected.im / root);
            let incident_power_w = normalized[k].scale(core::f64::consts::FRAC_1_SQRT_2).norm_sq();
            let reflected_power_w = reflected_normalized
                .scale(core::f64::consts::FRAC_1_SQRT_2).norm_sq();
            let net_input_power_w = original.mean_input_power_w;
            let power_residual_w = incident_power_w - reflected_power_w - net_input_power_w;
            if !finite(reflected) || [
                incident_power_w, reflected_power_w, net_input_power_w, power_residual_w,
            ].iter().any(|value| !value.is_finite()) {
                return Err(invalid("DG scattering wave power is not representable"));
            }
            wave_supply += incident_power_w - reflected_power_w;
            waves.push(DgScatteringPort {
                id, reference_impedance_pa_s_m3: reference,
                incident_pressure_pa: amplitudes[k], reflected_pressure_pa: reflected,
                incident_power_w, reflected_power_w, net_input_power_w, power_residual_w,
            });
        }
        let native_loss = physical.wall_dissipation_w + physical.molecular_dissipation_w
            + physical.matched_absorption_w;
        let wave_power_residual_w = native_loss - wave_supply;
        if !wave_power_residual_w.is_finite() {
            return Err(invalid("DG scattering wave balance is not representable"));
        }
        poll(cancelled)?;
        Ok(DgScatteringResponse { physical, waves, wave_power_residual_w })
    }
}
