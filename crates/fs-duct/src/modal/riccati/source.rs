//! Physical inlet sources in the original continuous pressure/flow chart.
//!
//! Preparation retains the original Magnus station maps once. A source query
//! solves `zref X a + Zs Y a = p_open` or `Y a = U_prescribed`, then walks the
//! same reference-wave maps toward the mouth. Neither an inlet impedance nor
//! a native characteristic impedance at transverse cutoff is required.

use super::*;
use crate::modal::ModalSourceResponse;
use crate::network::{PortExcitation, coherent_sum};

/// One coherent drive and its complete continuous modal pressure/flow field.
#[derive(Debug, Clone)]
pub struct ModalSourceField {
    /// Physical inlet fields, signed source work and power entering the duct.
    pub inlet: ModalSourceResponse,
    /// Every original integration station, including the mouth and inlet.
    pub field: ModalField,
}

/// One frozen smooth profile, material, frequency and terminal boundary.
///
/// Source changes reuse all original propagation and backward station maps.
/// Retained storage is O(steps N²); each drive takes one at-most-N³ boundary
/// solve and O(steps N²) field reconstruction, with no matrix exponentials.
/// `steps` and `ModeCount` retain their existing finite admission caps and
/// independent refinement meaning. These are harmonic, axisymmetric fields;
/// no causal source-impedance realization or modal truncation bound is inferred.
pub struct PreparedModalProfile {
    record: Stations,
    axis_weight: Vec<f64>,
    gas: GasState,
    omega: f64,
    zref: f64,
    length: f64,
    min_shear: f64,
    mouth_ka: f64,
}

impl PreparedModalProfile {
    /// Prepare scalar/ideal/tabulated/matched termination with the original
    /// continuous Magnus propagator. Cancellation is checked at every step.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        profile: &RadiusProfile,
        state: &GasState,
        omega: f64,
        loss: LossModel,
        load: MouthLoad<'_>,
        modes: ModeCount,
        steps: RiccatiSteps,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, DuctError> {
        Self::prepare(
            profile,
            state,
            omega,
            loss,
            load.into(),
            modes,
            steps,
            cancelled,
        )
    }

    /// Prepare the complete flanged radiation matrix in its original exterior
    /// gas. All signed mutual mouth loads act on the same continuous field.
    #[allow(clippy::too_many_arguments)]
    pub fn with_flanged_mouth(
        profile: &RadiusProfile,
        state: &GasState,
        exterior: &GasState,
        omega: f64,
        loss: LossModel,
        modes: ModeCount,
        steps: RiccatiSteps,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, DuctError> {
        Self::prepare(
            profile,
            state,
            omega,
            loss,
            ProfileLoad::Flanged(exterior),
            modes,
            steps,
            cancelled,
        )
    }

    /// Prepare an original full physical mouth matrix at its declared radius,
    /// frequency, gas and retained mode count. No diagonal approximation occurs.
    #[allow(clippy::too_many_arguments)]
    pub fn with_mouth_matrix(
        profile: &RadiusProfile,
        state: &GasState,
        omega: f64,
        loss: LossModel,
        load: &ModalMouthImpedance,
        modes: ModeCount,
        steps: RiccatiSteps,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, DuctError> {
        poll(cancelled, 0)?;
        load.admit(state, omega, profile.outlet_radius(), modes.0)?;
        Self::prepare(
            profile,
            state,
            omega,
            loss,
            ProfileLoad::Matrix(load),
            modes,
            steps,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        profile: &RadiusProfile,
        state: &GasState,
        omega: f64,
        loss: LossModel,
        load: ProfileLoad<'_>,
        modes: ModeCount,
        steps: RiccatiSteps,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, DuctError> {
        poll(cancelled, 0)?;
        modes.validate()?;
        let run = propagate_controlled(
            profile, state, omega, loss, load, modes.0, steps, true, cancelled,
        )?;
        Ok(Self {
            record: run.stations.ok_or(DuctError::Singular)?,
            axis_weight: axis_weights(modes.0),
            gas: *state,
            omega,
            zref: run.zref,
            length: profile.length(),
            min_shear: run.min_shear,
            mouth_ka: run.mouth_ka,
        })
    }

    /// Original material snapshot, never an independently substituted rho/c pair.
    pub const fn gas_state(&self) -> &GasState {
        &self.gas
    }
    /// Frozen positive angular frequency [rad/s].
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega
    }
    /// Retained radial modes, including the plane mode.
    pub fn modes(&self) -> ModeCount {
        ModeCount(self.axis_weight.len())
    }
    /// Original number of Magnus steps; no source query refines or discards them.
    pub fn steps_used(&self) -> usize {
        self.record.wave_back.len()
    }
    /// Positive real computational pressure/volume-flow reference [Pa s/m³].
    pub const fn reference_impedance_pa_s_m3(&self) -> f64 {
        self.zref
    }
    /// Minimum local plane-mode shear number retained from propagation.
    pub const fn minimum_shear_number(&self) -> f64 {
        self.min_shear
    }
    /// Original mouth acoustic radius, including its explicit exterior gas.
    pub const fn mouth_ka(&self) -> f64 {
        self.mouth_ka
    }

    /// Apply one coherent physical source boundary in original radial order.
    ///
    /// A zero volume-flow drive closes a mode. A zero open-circuit pressure
    /// retains its finite passive source impedance. Every nonzero source and
    /// every cross-mode coupling participates in one solve. Exactly disconnected
    /// unforced components take zero forced amplitude, including free resonances;
    /// no nonzero coupling is pruned. A singular driven component still refuses.
    ///
    /// The complete field uses the same source solution, including at pressure
    /// or flow nodes and transverse cutoffs. Source-resistance heat, signed
    /// generator work and power entering the duct remain separate. Cancellation
    /// at reconstruction stations publishes no partial field or mutable state.
    pub fn drive_sources(
        &self,
        sources: &[PortExcitation],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ModalSourceField, DuctError> {
        poll(cancelled, 0)?;
        let wave = self.source_wave(sources)?;
        let field = reconstruct_field(
            &self.record,
            self.zref,
            self.length,
            &self.axis_weight,
            wave,
            cancelled,
        )?;
        let inlet = ModalSourceResponse::from_fields(
            sources,
            field.pressure[0].clone(),
            field.flow[0].clone(),
        )?;
        Ok(ModalSourceField { inlet, field })
    }

    fn source_wave(&self, sources: &[PortExcitation]) -> Result<Vec<C64>, DuctError> {
        let n = self.modes().0;
        if sources.len() != n {
            return Err(bad(
                "continuous modal sources must cover every retained mode",
            ));
        }
        let mut connected = vec![false; n];
        let mut p_weight = vec![C64::ZERO; n];
        let mut u_weight = vec![C64::ZERO; n];
        let mut drive = vec![C64::ZERO; n];
        for (i, source) in sources.iter().enumerate() {
            match *source {
                PortExcitation::VolumeFlow(flow) => {
                    if !finite(flow) {
                        return Err(bad("continuous modal source flow must be finite"));
                    }
                    u_weight[i] = C64::ONE;
                    drive[i] = flow;
                }
                PortExcitation::Pressure {
                    open_circuit_pressure_pa: pressure,
                    impedance_pa_s_m3: z,
                } => {
                    if !finite(pressure) || !finite(z) || z.re < 0.0 {
                        return Err(bad(
                            "continuous modal source needs finite pressure and passive impedance",
                        ));
                    }
                    let scale = self.zref.max(z.re.abs()).max(z.im.abs());
                    p_weight[i] = C64::from_re(self.zref / scale);
                    u_weight[i] = z / C64::from_re(scale);
                    drive[i] = pressure / C64::from_re(scale);
                    if p_weight[i] == C64::ZERO || (pressure != C64::ZERO && drive[i] == C64::ZERO)
                    {
                        return Err(DuctError::Singular);
                    }
                }
            }
            if !finite(drive[i]) {
                return Err(DuctError::Singular);
            }
            connected[i] = drive[i] != C64::ZERO;
        }
        let p = self.record.pressure.last().ok_or(DuctError::Singular)?;
        let u = self.record.flow.last().ok_or(DuctError::Singular)?;
        loop {
            let mut changed = false;
            for i in 0..n {
                if connected[i] {
                    for j in 0..n {
                        if !connected[j]
                            && [p[i * n + j], p[j * n + i], u[i * n + j], u[j * n + i]]
                                .iter()
                                .any(|z| *z != C64::ZERO)
                        {
                            connected[j] = true;
                            changed = true;
                        }
                    }
                }
            }
            if !changed {
                break;
            }
        }
        let modes: Vec<_> = (0..n).filter(|&i| connected[i]).collect();
        let m = modes.len();
        let mut wave = vec![C64::ZERO; n];
        if m == 0 {
            return Ok(wave);
        }
        let mut matrix = vec![C64::ZERO; m * m];
        let mut rhs = vec![C64::ZERO; m];
        for (row, &i) in modes.iter().enumerate() {
            let entries = &mut matrix[row * m..(row + 1) * m];
            for (entry, &j) in entries.iter_mut().zip(&modes) {
                *entry = p_weight[i] * p[i * n + j] + u_weight[i] * u[i * n + j];
                if !finite(*entry) {
                    return Err(DuctError::Singular);
                }
            }
            let scale = entries
                .iter()
                .fold(0.0_f64, |s, z| s.max(z.re.abs()).max(z.im.abs()));
            if scale == 0.0 {
                return Err(DuctError::Singular);
            }
            for entry in entries {
                *entry = *entry / C64::from_re(scale);
            }
            rhs[row] = drive[i] / C64::from_re(scale);
            if !finite(rhs[row]) {
                return Err(DuctError::Singular);
            }
        }
        let lu = fs_la::eigen_complex::lu_complex(&matrix, m).map_err(|_| DuctError::Singular)?;
        lu.solve(&mut rhs);
        for (&i, amplitude) in modes.iter().zip(rhs) {
            if !finite(amplitude) {
                return Err(DuctError::Singular);
            }
            wave[i] = amplitude;
        }
        Ok(wave)
    }

    /// Complete reflected-pressure matrix at positive real modal references.
    ///
    /// `b = R a`, `a_i = (P_i + reference_i U_i)/2`. This uses one factorization
    /// of the retained pressure/flow pair, never an input-impedance inverse.
    /// Power uses `sum((|a_i|²-|b_i|²)/reference_i)/2`; native evanescent modes
    /// and cutoff fields remain in the original underlying station maps.
    pub fn at_references(&self, references: &[f64]) -> Result<Vec<C64>, DuctError> {
        let n = self.modes().0;
        if references.len() != n || references.iter().any(|r| !r.is_finite() || *r <= 0.0) {
            return Err(bad(
                "continuous modal reflection needs positive real modal references",
            ));
        }
        let p = self.record.pressure.last().ok_or(DuctError::Singular)?;
        let u = self.record.flow.last().ok_or(DuctError::Singular)?;
        let scales: Vec<_> = references.iter().map(|r| self.zref.max(*r)).collect();
        let mut incoming = vec![C64::ZERO; n * n];
        let mut outgoing = vec![C64::ZERO; n * n];
        let mut row_scales = vec![0.0_f64; n];
        for i in 0..n {
            if self.zref / scales[i] == 0.0 || references[i] / scales[i] == 0.0 {
                return Err(DuctError::Singular);
            }
            for j in 0..n {
                incoming[i * n + j] = p[i * n + j].scale(self.zref / scales[i])
                    + u[i * n + j].scale(references[i] / scales[i]);
                outgoing[i * n + j] = p[i * n + j].scale(self.zref / scales[i])
                    - u[i * n + j].scale(references[i] / scales[i]);
            }
            let row = &mut incoming[i * n..(i + 1) * n];
            let scale = row
                .iter()
                .fold(0.0_f64, |s, z| s.max(z.re.abs()).max(z.im.abs()));
            if scale == 0.0 || !scale.is_finite() {
                return Err(DuctError::Singular);
            }
            for value in row {
                *value = *value / C64::from_re(scale);
            }
            row_scales[i] = scale;
        }
        let lu = fs_la::eigen_complex::lu_complex(&incoming, n).map_err(|_| DuctError::Singular)?;
        let mut result = vec![C64::ZERO; n * n];
        let mut wave = vec![C64::ZERO; n];
        for j in 0..n {
            wave.fill(C64::ZERO);
            wave[j] = C64::from_re((2.0 / scales[j]) / row_scales[j]);
            if !finite(wave[j]) {
                return Err(DuctError::Singular);
            }
            lu.solve(&mut wave);
            for i in 0..n {
                result[i * n + j] = coherent_sum((0..n).map(|k| outgoing[i * n + k] * wave[k]))
                    .scale(0.5 * scales[i]);
                if !finite(result[i * n + j]) {
                    return Err(DuctError::Singular);
                }
            }
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
