//! Finite reflected waves without an intermediate input-impedance inverse.

use super::*;
use crate::network::{PortExcitation, coherent_sum};

/// Original inlet fields and source work from one coherent modal boundary solve.
///
/// All vectors retain the original radial mode order. Peak phasors use
/// `exp(-i omega t)` and positive inlet volume flow enters the duct. Source
/// impedance dissipation is separate from the gas, wall and mouth losses.
#[derive(Debug, Clone)]
pub struct ModalSourceResponse {
    /// Actual complex inlet modal pressures [Pa].
    pub inlet_pressure: Vec<C64>,
    /// Actual complex inlet modal volume flows [m³/s], into the duct.
    pub inlet_flow: Vec<C64>,
    /// Signed cycle-mean generator work [W] for each mode. A pressure source
    /// uses `Re(p_open conj(U))/2`; a flow source uses its actual inlet work.
    /// A driven generator may absorb work from the other coherent sources.
    pub generator_powers_w: Vec<f64>,
    /// Nonnegative cycle-mean loss [W] in each source's declared resistance,
    /// `Re(Zs) |U|²/2`. Ideal flow sources have no internal impedance loss.
    pub source_dissipations_w: Vec<f64>,
    /// Signed total power entering the duct, `sum Re(P conj(U))/2` [W].
    pub inlet_power_w: f64,
    /// Signed total generator power [W], before internal source dissipation.
    pub generator_power_w: f64,
    /// Total loss [W] in the source resistances. Generator power equals this
    /// loss plus inlet power, up to the boundary solve's numerical error.
    pub source_dissipation_w: f64,
}

/// Staircase input reflection in the inlet station's native modal waves.
///
/// Pressure and volume flow satisfy `P = (I + R) a` and
/// `U = diag(1/Zc) (I - R) a`. The response is retained before the impedance
/// conversion, so an inlet flow node (`R` with eigenvalue +1) need not fail a
/// reflection query. Native evanescent waves are not power-normalized; use
/// [`Self::at_references`] for a real-reference scattering port.
#[derive(Debug, Clone)]
pub struct ModalReflectionResponse {
    /// Angular frequency [rad/s].
    pub omega: f64,
    /// Retained radial mode count, including the plane mode.
    pub n_modes: usize,
    /// Native pressure-wave reflection matrix, row-major N x N (`b = R a`).
    pub reflection_matrix: Vec<C64>,
    /// Native per-mode characteristic impedances [Pa s/m^3].
    pub characteristic_impedances: Vec<C64>,
    /// Smallest plane-mode shear number along the chain.
    pub min_shear_number: f64,
    /// Mouth `ka`, using the exterior gas for a full flanged radiation load.
    pub mouth_ka: f64,
    /// Cutoff frequencies of modes 1..N-1 at the narrowest station [Hz].
    pub cutoffs_at_throat_hz: Vec<f64>,
    /// Number of staircase stations used in the recurrence.
    pub steps_used: usize,
}

impl ModalReflectionResponse {
    /// Reflection against positive real per-mode pressure-wave references.
    ///
    /// For `D = diag(references)`, form `A = (P + D U)/2` and
    /// `B = (P - D U)/2` from the same native waves and solve `R_new A = B`.
    /// This right solve preserves modal coupling and never inverts `I - R`.
    /// The returned waves carry signed power
    /// `sum((|a_i|^2 - |b_i|^2) / references[i]) / 2`; unequal references
    /// require that metric when testing passivity or reciprocity.
    ///
    /// # Errors
    /// Invalid shapes/references or a singular wave conversion. Modal cutoff,
    /// load and staircase validity limits are unchanged.
    pub fn at_references(&self, references: &[f64]) -> Result<Vec<C64>, DuctError> {
        let (a, mut reflected, _) = self.reference_maps(references)?;
        let n = self.n_modes;
        let lu = lu_complex(&transpose(&a, n), n).map_err(|_| DuctError::Singular)?;
        for row in reflected.chunks_exact_mut(n) {
            lu.solve(row);
        }
        finite_waves(&reflected)?;
        Ok(reflected)
    }

    // Return twice the incident/reflected real-reference waves, together with
    // the native incident-basis scale used in BOTH maps. Observers reuse this
    // same boundary solve rather than reconstructing a native wave through Z.
    fn reference_maps(&self, references: &[f64])
        -> Result<(Vec<C64>, Vec<C64>, Vec<f64>), DuctError>
    {
        let (p, q) = self.pressure_and_scaled_flow(references)?;
        let n = self.n_modes;
        // The factor 1/2 cancels between B and A. Solving A^T x = B_row^T
        // avoids an explicit inverse and the incorrect noncommuting A^-1 B.
        let mut a: Vec<_> = p.iter().zip(&q).map(|(p, q)| *p + *q).collect();
        let mut reflected: Vec<_> = p.iter().zip(&q).map(|(p, q)| *p - *q).collect();
        finite_waves(&a)?;
        finite_waves(&reflected)?;
        // Re-scale the native incident basis in both maps. This leaves B A^-1
        // unchanged and keeps quasi-static evanescent ports from setting the
        // entire factorization's scale through their very small native Zc.
        let mut scales = Vec::with_capacity(n);
        for j in 0..n {
            let scale = (0..n).map(|i| a[i * n + j].abs()).fold(0.0, f64::max);
            if scale == 0.0 {
                return Err(DuctError::Singular);
            }
            for i in 0..n {
                a[i * n + j] = a[i * n + j] / C64::from_re(scale);
                reflected[i * n + j] = reflected[i * n + j] / C64::from_re(scale);
            }
            scales.push(scale);
        }
        Ok((a, reflected, scales))
    }

    // One coherent incident vector: retain the solved native amplitudes for
    // mouth observation while returning the reflected real-reference waves.
    pub(super) fn excitation(&self, references: &[f64], incident: &[C64])
        -> Result<(Vec<C64>, Vec<C64>), DuctError>
    {
        let (a, b, scales) = self.reference_maps(references)?;
        let n = scales.len();
        if incident.len() != n {
            return Err(DuctError::BadParameter { what: "modal incident pressure must cover every retained mode" });
        }
        finite_waves(incident)?;
        let lu = lu_complex(&a, n).map_err(|_| DuctError::Singular)?;
        let mut native: Vec<_> = incident.iter().map(|value| value.scale(2.0)).collect();
        finite_waves(&native)?;
        lu.solve(&mut native);
        let reflected: Vec<_> = (0..n).map(|i| (0..n)
            .fold(C64::ZERO, |sum, j| sum + b[i*n+j]*native[j]).scale(0.5)).collect();
        for (value, scale) in native.iter_mut().zip(scales) {
            *value = *value / C64::from_re(scale);
        }
        finite_waves(&native)?;
        finite_waves(&reflected)?;
        Ok((reflected, native))
    }

    /// Plane-source reflection with zero higher-mode inlet volume flow.
    ///
    /// This is the closure associated with the impedance matrix's `Z00`:
    /// `a0 = (P0 + reference U0)/2`, `b0 = (P0 - reference U0)/2`, and
    /// `U_i = 0` for i > 0. It differs from the (0,0) entry of a full
    /// scattering matrix, which instead matches the other incident ports.
    /// Solving those boundary equations directly retains finite plane
    /// reflection at a pressure antinode without forming `Z00` or the full Z.
    /// Exactly disconnected higher-mode blocks have zero source-driven
    /// amplitude, including at their free resonances. No nonzero coupling is
    /// discarded, and a singular source-connected block still refuses.
    ///
    /// # Errors
    /// Invalid response/reference or a singular source-boundary solve.
    pub fn plane_reflection(&self, reference: f64) -> Result<C64, DuctError> {
        self.plane_excitation(reference).map(|(reflected, _)| reflected)
    }

    pub(super) fn plane_excitation(&self, reference: f64)
        -> Result<(C64, Vec<C64>), DuctError>
    {
        let n = self.n_modes;
        let references = [reference; MAX_MODES];
        let references = references.get(..n).ok_or(DuctError::BadParameter {
            what: "modal reflection mode count exceeds MAX_MODES",
        })?;
        let (p, q) = self.pressure_and_scaled_flow(references)?;
        // A closed higher mode may have R=+1 while the driven plane response
        // remains finite. Its homogeneous equation then reads 0=0, not a
        // singular source response. Retain exactly the connected modal block
        // and select zero free amplitude outside it. Test BOTH directions so
        // nonreciprocal user-authored matrices cannot hide a coupling. Exact
        // zeros alone separate blocks; arbitrarily weak coupling is retained.
        let mut connected = [false; MAX_MODES];
        connected[0] = true;
        let (modes, count) = self.connected_modes(connected);
        let modes = &modes[..count];
        let mut boundary = vec![C64::ZERO; count * count];
        for (row, &i) in modes.iter().enumerate() {
            for (column, &j) in modes.iter().enumerate() {
                boundary[row * count + column] = if i == 0 {
                    p[j] + q[j]
                } else {
                    let value = -self.reflection_matrix[i * n + j];
                    if i == j { value + C64::ONE } else { value }
                };
            }
        }
        let lu = lu_complex(&boundary, count).map_err(|_| DuctError::Singular)?;
        let mut source = [C64::ZERO; MAX_MODES];
        source[0] = C64::from_re(2.0);
        lu.solve(&mut source[..count]);
        let mut outgoing = vec![C64::ZERO; n];
        for (&mode, &amplitude) in modes.iter().zip(&source) {
            outgoing[mode] = amplitude;
        }
        let reflected = (0..n)
            .fold(C64::ZERO, |sum, j| sum + (p[j] - q[j]) * outgoing[j])
            .scale(0.5);
        finite_waves(&[reflected])?;
        finite_waves(&outgoing)?;
        Ok((reflected, outgoing))
    }

    /// Drive every retained inlet mode with the existing physical source law.
    ///
    /// Sources are in original modal order. `VolumeFlow(0)` closes that mode;
    /// a `Pressure` source enforces `P + Zs U = p_open`, including its finite
    /// passive complex impedance. Zero `p_open` retains that load, while zero
    /// `Zs` is an ideal pressure source. All nonzero drives share this
    /// response's frequency and phase clock. Coherent higher-mode drives and
    /// mixed pressure/flow conditions retain every original modal coupling.
    ///
    /// The solve uses `(I+R)a` and `diag(1/Zc)(I-R)a` directly, without an
    /// input-impedance inverse. Exactly disconnected unforced components have
    /// zero amplitude, including free resonances; a singular driven component
    /// still refuses. This is one harmonic sample, not a causal model for an
    /// arbitrary frequency-dependent source impedance.
    ///
    /// # Errors
    /// Invalid response/source dimensions, nonfinite drive or impedance,
    /// negative source resistance, nonrepresentable fields/work, or a singular
    /// source-connected boundary solve. Original modal admissions still apply.
    pub fn drive_sources(&self, sources: &[PortExcitation])
        -> Result<ModalSourceResponse, DuctError>
    {
        self.source_excitation(sources).map(|(response, _)| response)
    }

    pub(super) fn source_excitation(&self, sources: &[PortExcitation])
        -> Result<(ModalSourceResponse, Vec<C64>), DuctError>
    {
        self.admit_native()?;
        let n = self.n_modes;
        if sources.len() != n {
            return Err(DuctError::BadParameter {
                what: "modal sources must prescribe every retained inlet mode",
            });
        }
        let mut connected = [false; MAX_MODES];
        let mut p_weight = [C64::ZERO; MAX_MODES];
        let mut u_weight = [C64::ZERO; MAX_MODES];
        let mut drive = [C64::ZERO; MAX_MODES];
        for (i, &source) in sources.iter().enumerate() {
            let zc = self.characteristic_impedances[i];
            match source {
                PortExcitation::VolumeFlow(flow) => {
                    finite_waves(&[flow])?;
                    connected[i] = flow != C64::ZERO;
                    u_weight[i] = C64::ONE;
                    drive[i] = zc * flow;
                }
                PortExcitation::Pressure { open_circuit_pressure_pa: pressure,
                    impedance_pa_s_m3: impedance } => {
                    finite_waves(&[pressure, impedance])?;
                    if impedance.re < 0.0 {
                        return Err(DuctError::BadParameter {
                            what: "modal source impedance must have nonnegative resistance",
                        });
                    }
                    connected[i] = pressure != C64::ZERO;
                    // Multiply P+Zs U=p_open by Zc/scale. Bounded impedance
                    // weights avoid forming an overflowing Zs/Zc ratio and
                    // keep evanescent native impedances in their own phase.
                    let scale = zc.re.abs().max(zc.im.abs())
                        .max(impedance.re.abs()).max(impedance.im.abs());
                    p_weight[i] = zc / C64::from_re(scale);
                    u_weight[i] = impedance / C64::from_re(scale);
                    if p_weight[i] == C64::ZERO {
                        return Err(DuctError::Singular);
                    }
                    drive[i] = p_weight[i] * pressure;
                }
            }
        }
        finite_waves(&drive[..n])?;
        let (modes, count) = self.connected_modes(connected);
        let modes = &modes[..count];
        let mut native = vec![C64::ZERO; n];
        if count != 0 {
            let mut boundary = [C64::ZERO; MAX_MODES * MAX_MODES];
            let mut rhs = [C64::ZERO; MAX_MODES];
            for (row, &i) in modes.iter().enumerate() {
                for (column, &j) in modes.iter().enumerate() {
                    let r = self.reflection_matrix[i*n+j];
                    let delta = if i == j { C64::ONE } else { C64::ZERO };
                    boundary[row*count+column] = p_weight[i]*(delta+r) + u_weight[i]*(delta-r);
                }
                let values = &mut boundary[row*count..(row+1)*count];
                finite_waves(values)?;
                let scale = values.iter().fold(0.0_f64,
                    |scale, z| scale.max(z.re.abs()).max(z.im.abs()));
                if scale == 0.0 { return Err(DuctError::Singular); }
                for value in values { *value = *value / C64::from_re(scale); }
                rhs[row] = drive[i] / C64::from_re(scale);
            }
            finite_waves(&rhs[..count])?;
            let lu = lu_complex(&boundary[..count*count], count).map_err(|_| DuctError::Singular)?;
            lu.solve(&mut rhs[..count]);
            for (&i, &amplitude) in modes.iter().zip(&rhs) { native[i] = amplitude; }
            finite_waves(&native)?;
        }
        let mut inlet_pressure = Vec::with_capacity(n);
        let mut inlet_flow = Vec::with_capacity(n);
        for i in 0..n {
            inlet_pressure.push(coherent_sum((0..n).map(|j| {
                let delta = if i == j { C64::ONE } else { C64::ZERO };
                (delta+self.reflection_matrix[i*n+j])*native[j]
            })));
            inlet_flow.push(coherent_sum((0..n).map(|j| {
                let delta = if i == j { C64::ONE } else { C64::ZERO };
                (delta-self.reflection_matrix[i*n+j])*native[j]
            })) / self.characteristic_impedances[i]);
        }
        finite_waves(&inlet_pressure)?;
        finite_waves(&inlet_flow)?;
        let mut generator_powers_w = Vec::with_capacity(n);
        let mut source_dissipations_w = Vec::with_capacity(n);
        for (i, &source) in sources.iter().enumerate() {
            let (pressure, resistance) = match source {
                PortExcitation::VolumeFlow(_) => (inlet_pressure[i], 0.0),
                PortExcitation::Pressure { open_circuit_pressure_pa, impedance_pa_s_m3 } =>
                    (open_circuit_pressure_pa, impedance_pa_s_m3.re),
            };
            let flow = inlet_flow[i];
            generator_powers_w.push((pressure*flow.conj()).re*0.5);
            // Form the dissipative pressure first, so a tiny flow's squared
            // magnitude does not underflow before a large resistance acts.
            source_dissipations_w.push((flow.scale(resistance)*flow.conj()).re*0.5);
        }
        let generator_power_w = coherent_sum(generator_powers_w.iter().copied().map(C64::from_re)).re;
        let source_dissipation_w = coherent_sum(source_dissipations_w.iter().copied().map(C64::from_re)).re;
        let inlet_power_w = coherent_sum(inlet_pressure.iter().zip(&inlet_flow)
            .map(|(&pressure, &flow)| (pressure*flow.conj()).scale(0.5))).re;
        if generator_powers_w.iter().chain(&source_dissipations_w)
            .chain(&[generator_power_w, source_dissipation_w, inlet_power_w])
            .any(|value| !value.is_finite())
        {
            return Err(DuctError::Singular);
        }
        Ok((ModalSourceResponse { inlet_pressure, inlet_flow, generator_powers_w,
            source_dissipations_w, inlet_power_w, generator_power_w, source_dissipation_w }, native))
    }

    // The exact undirected coupling graph includes every nonzero entry in
    // either direction. Callers have already admitted the native dimensions.
    fn connected_modes(&self, mut connected: [bool; MAX_MODES]) -> ([usize; MAX_MODES], usize) {
        let n = self.n_modes;
        loop {
            let mut changed = false;
            for i in 0..n {
                if connected[i] {
                    for j in 0..n {
                        if !connected[j]
                            && (self.reflection_matrix[i * n + j] != C64::ZERO
                                || self.reflection_matrix[j * n + i] != C64::ZERO)
                        {
                            connected[j] = true;
                            changed = true;
                        }
                    }
                }
            }
            if !changed { break; }
        }
        let mut modes = [0; MAX_MODES];
        let mut count = 0;
        for (i, &included) in connected[..n].iter().enumerate() {
            if included {
                modes[count] = i;
                count += 1;
            }
        }
        (modes, count)
    }

    fn admit_native(&self) -> Result<(), DuctError> {
        let n = self.n_modes;
        if !(1..=MAX_MODES).contains(&n)
            || self.reflection_matrix.len() != n * n
            || self.characteristic_impedances.len() != n
            || self.characteristic_impedances.iter().any(|z| {
                !z.re.is_finite() || !z.im.is_finite() || z.abs() == 0.0
            })
        {
            return Err(DuctError::BadParameter {
                what: "modal reflection needs an N x N matrix and finite nonzero native impedances",
            });
        }
        finite_waves(&self.reflection_matrix)?;
        Ok(())
    }

    fn pressure_and_scaled_flow(&self, references: &[f64]) -> Result<(Vec<C64>, Vec<C64>), DuctError> {
        self.admit_native()?;
        let n = self.n_modes;
        if references.len() != n || references.iter().any(|z| !z.is_finite() || *z <= 0.0) {
            return Err(DuctError::BadParameter { what: "modal reflection needs positive real references for every mode" });
        }
        let mut p = self.reflection_matrix.clone();
        let mut q = vec![C64::ZERO; n * n];
        for i in 0..n {
            let ratio = C64::from_re(references[i]) / self.characteristic_impedances[i];
            for j in 0..n {
                let delta = if i == j { C64::ONE } else { C64::ZERO };
                p[i * n + j] = delta + p[i * n + j];
                q[i * n + j] = ratio * (delta - self.reflection_matrix[i * n + j]);
            }
        }
        finite_waves(&p)?;
        finite_waves(&q)?;
        Ok((p, q))
    }
}

fn finite_waves(waves: &[C64]) -> Result<(), DuctError> {
    if waves.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
        return Err(DuctError::BadParameter { what: "modal reflection waves must be finite" });
    }
    Ok(())
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod source_tests;
