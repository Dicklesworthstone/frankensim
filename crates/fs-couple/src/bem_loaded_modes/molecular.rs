//! Original source mechanics loaded by the same molecular field they radiate.
//!
//! At one positive frequency the boundary's complete unit-source solutions
//! give its exact discrete Schur complement, not a fitted or truncated load.
//! Every original mechanical and Maxwell-memory coordinate remains. Receiver
//! fields and wall/exterior work use the same molecular source-solution bank.
use super::*;
use fs_bem::{
    helmholtz::molecular::{MolecularError, MolecularExterior, MolecularNearFieldOptions},
    observation_transfers::{
        HarmonicField, LocalBoundary, MolecularBemObservation, MolecularObservationError,
        MolecularObservationFormulation, MolecularReceiverError, MolecularResponseLimits,
        PreparedMolecularResponse, ReceiverSpec, SourceId, TransferError,
    },
};
use fs_la::eigen_complex::LuComplex;

mod incident;
pub use incident::{MolecularIncidentLimits, MolecularIncidentPower, MolecularIncidentResponse};

/// Original mechanics/solve refusal or original material/boundary observation
/// refusal. No source error is replaced with a generic successful response.
#[derive(Debug)]
pub enum MolecularModesError {
    Modes(BemModesError),
    Observation(MolecularObservationError),
}
impl From<BemModesError> for MolecularModesError {
    fn from(value: BemModesError) -> Self {
        Self::Modes(value)
    }
}
impl From<MolecularObservationError> for MolecularModesError {
    fn from(value: MolecularObservationError) -> Self {
        Self::Observation(value)
    }
}
impl From<MolecularError> for MolecularModesError {
    fn from(value: MolecularError) -> Self {
        Self::Observation(value.into())
    }
}
impl From<TransferError> for MolecularModesError {
    fn from(value: TransferError) -> Self {
        Self::Observation(value.into())
    }
}
impl core::fmt::Display for MolecularModesError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Modes(e) => e.fmt(f),
            Self::Observation(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for MolecularModesError {}

/// Cycle-average powers for peak complex amplitudes, in watts. Exterior work
/// includes outgoing flux and eventual gas relaxation heat. It is not the
/// distant flux alone; no duplicate downstream absorption is applied.
#[derive(Debug, Clone, Copy)]
pub struct MolecularModesPower {
    pub input_w: f64,
    pub kelvin_w: f64,
    pub relaxation_w: f64,
    pub exterior_w: f64,
    pub wall_w: f64,
    pub balance_w: f64,
    pub apparent_power_scale_w: f64,
    pub relative_balance: f64,
}

/// Actual force-scaled motion and matched receivers in original source order.
/// Modal displacement/velocity use the original unit-mass coordinates, and
/// fluid reaction is force ON that source: -Z times its solved modal velocity.
#[derive(Debug, Clone)]
pub struct MolecularLoadedResponse {
    pub angular_frequency_rad_s: f64,
    pub modal_displacement: Vec<C64>,
    pub modal_velocity: Vec<C64>,
    pub elastic_memory_sqrt_j: Vec<C64>,
    pub fluid_reaction: Vec<C64>,
    /// Same order and physical locations as `radiation().harmonic().receivers()`.
    pub receivers: Vec<HarmonicField>,
    /// Absolute adaptive-quadrature estimates for these actual coherent fields.
    /// They do not bound boundary discretization, source error, or relative
    /// error near cancellation. None denotes the original centroid receiver.
    pub receiver_quadrature_error_estimates: Option<Vec<MolecularReceiverError>>,
    pub mechanical_residual: f64,
    pub memory_residual: f64,
    pub power: MolecularModesPower,
}

/// One molecular boundary solve bank and one cached mechanical/memory LU.
/// Borrowing the original source prevents swapping its mode basis, geometry or
/// constitutive memory after preparation. This is a real-frequency response;
/// no holomorphic pole continuation or time-domain realization is asserted.
pub struct PreparedMolecularLoadedModes<'s> {
    source: &'s BemLoadedModes,
    exterior: MolecularExterior,
    boundary: Option<LocalBoundary>,
    minimum_clearance_m: f64,
    near_field_options: Option<MolecularNearFieldOptions>,
    radiation: PreparedMolecularResponse,
    raw: PhysicalPencil,
    lu: LuComplex,
    rows: Vec<f64>,
    columns: Vec<f64>,
    options: BemHarmonicOptions,
}

impl BemLoadedModes {
    /// Couple the ORIGINAL source skin, complete mechanics and Maxwell arms to
    /// one declared relaxing exterior and optional original local wall laws.
    /// Source-supported panels remain rigid under the wall owner's admission;
    /// other passive patches may absorb, store or reradiate their pressure work.
    /// The supplied gas must retain the model's original equilibrium rho/c.
    ///
    /// D = K-omega² I-i*omega*(C+Z), with every original Maxwell coordinate
    /// appended by its unchanged owner. Z is assembled from this exact source's
    /// normal-velocity weights, never an unrelated externally supplied matrix.
    /// Boundary elimination is exact for the selected discrete equations; no
    /// source mode, material arm or self/mutual acoustic entry is dropped.
    ///
    /// Receivers must be nonempty and satisfy the original field-clearance
    /// admission; no dummy receiver or default position is manufactured.
    ///
    /// The original `max_dense_bytes` bounds mechanical preparation/solve scratch's
    /// 96*d²+256*d+4096 logical bytes, with d=modes+memories. Boundary/source-bank
    /// storage uses its separate `response_limits` and original ComplexOptions.
    /// `solve_options` bounds every force batch and is retained unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_molecular_harmonic(
        &self,
        exterior: &MolecularExterior,
        omega: f64,
        receivers: &[ReceiverSpec],
        minimum_clearance_m: f64,
        boundary: Option<LocalBoundary>,
        response_limits: MolecularResponseLimits,
        solve_options: BemHarmonicOptions,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMolecularLoadedModes<'_>, MolecularModesError> {
        self.prepare_molecular_harmonic_with_receivers(
            exterior,
            omega,
            receivers,
            minimum_clearance_m,
            boundary,
            response_limits,
            solve_options,
            None,
            cancelled,
        )
    }

    /// The same original loaded mechanics with adaptive retained-triangle
    /// pressure/particle-velocity receivers. One frequency/receiver row bank
    /// serves every source column; later force cases perform no new quadrature.
    /// Selecting these receivers does not change the boundary load, mechanical
    /// operator, material memory, or work balance. Quadrature estimates remain
    /// separate from the unchanged boundary-discretization limitations.
    ///
    /// Receiver rows use the supplied independent storage/work/depth limits.
    /// Each force response additionally reserves 32 bytes per receiver for its
    /// four absolute quadrature estimates in max_response_bytes.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_molecular_harmonic_adaptive(
        &self,
        exterior: &MolecularExterior,
        omega: f64,
        receivers: &[ReceiverSpec],
        minimum_clearance_m: f64,
        boundary: Option<LocalBoundary>,
        response_limits: MolecularResponseLimits,
        solve_options: BemHarmonicOptions,
        near_field_options: MolecularNearFieldOptions,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMolecularLoadedModes<'_>, MolecularModesError> {
        self.prepare_molecular_harmonic_with_receivers(
            exterior,
            omega,
            receivers,
            minimum_clearance_m,
            boundary,
            response_limits,
            solve_options,
            Some(near_field_options),
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare_molecular_harmonic_with_receivers(
        &self,
        exterior: &MolecularExterior,
        omega: f64,
        receivers: &[ReceiverSpec],
        minimum_clearance_m: f64,
        boundary: Option<LocalBoundary>,
        response_limits: MolecularResponseLimits,
        solve_options: BemHarmonicOptions,
        near_field_options: Option<MolecularNearFieldOptions>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMolecularLoadedModes<'_>, MolecularModesError> {
        poll(&mut cancelled)?;
        let m = self.mode_count();
        let n = m + self.memory_count();
        let bytes = 96 * (n as u128).pow(2) + 256 * n as u128 + 4096;
        let work = 32 * (n as u128).pow(3) + 128 * (n as u128).pow(2);
        if !omega.is_finite()
            || omega <= 0.0
            || !(omega * omega).is_finite()
            || omega * omega == 0.0
            || self.medium != Medium::from(exterior.source().gas())
            || solve_options.max_load_cases == 0
            || solve_options.max_load_cases > 256
            || !solve_options.power_relative_tolerance.is_finite()
            || solve_options.power_relative_tolerance <= 0.0
            || solve_options.power_relative_tolerance >= 1.0
            || bytes > self.settings.max_dense_bytes as u128
            || work > solve_options.max_work as u128
        {
            return Err(bad("molecular loading needs the original gas, positive frequency and admitted preparation/solve limits").into());
        }
        let formulation = match self.settings.formulation {
            Formulation::PlainCbie => MolecularObservationFormulation::TriangleCbie,
            Formulation::BurtonMiller => MolecularObservationFormulation::CentroidBurtonMiller,
            Formulation::BurtonMillerWrongAlphaSign => {
                return Err(bad(
                    "molecular structural loading refuses the wrong-sign mutation formulation",
                )
                .into());
            }
        };
        let sources: Vec<_> = (0..m)
            .map(|i| SourceId::StructuralMode(format!("original-structural-mode-{i}")))
            .collect();
        let mut observer = MolecularBemObservation::new_with_formulation(
            &self.surface,
            exterior.source().clone(),
            &sources,
            &self.weights,
            receivers,
            minimum_clearance_m,
            exterior.band_rad_s(),
            formulation,
            self.settings.boundary,
        )?;
        if let Some(boundary) = &boundary {
            observer = observer.with_local_boundary(boundary.clone())?;
        }
        if let Some(options) = near_field_options {
            observer = observer.with_near_field(options)?;
        }
        let radiation = observer.prepare_response(omega, response_limits, &mut cancelled)?;
        let z = radiation.impedance();
        let mut matrix = vec![C64::ZERO; n * n];
        let mut envelope = vec![0.0; n * n];
        let mut derivative_scratch = vec![C64::ZERO; n * n];
        for i in 0..m {
            poll(&mut cancelled)?;
            for j in 0..m {
                let index = i * n + j;
                let (k, c) = (self.stiffness[i * m + j], self.viscosity[i * m + j]);
                matrix[index] =
                    C64::from_re(k) + C64::new(0.0, -omega) * (C64::from_re(c) + z[i * m + j]);
                envelope[index] = k.abs() + omega * (c.abs() + z[i * m + j].abs());
                if i == j {
                    matrix[index] = matrix[index] - C64::from_re(omega * omega);
                    envelope[index] += omega * omega;
                }
            }
        }
        self.memory.add_to_pencil(
            C64::from_re(omega),
            m,
            0,
            &mut matrix,
            &mut derivative_scratch,
            &mut envelope,
            &mut cancelled,
        )?;
        if matrix.iter().any(|&z| !finite(z)) || envelope.iter().any(|v| !v.is_finite()) {
            return Err(bad("molecular loaded mechanical operator is unrepresentable").into());
        }
        // This carrier is used ONLY for fixed balancing and original-row
        // residuals. A frequency derivative of the eliminated load is neither
        // computed nor exposed, and this object cannot enter the pole API.
        let raw = PhysicalPencil {
            value: PoleEvaluation {
                matrix,
                derivative: Vec::new(),
                residual_scale: 1.0,
            },
            envelope,
        };
        let (rows, columns) = raw.fixed_scales(n, &mut cancelled)?;
        let mut balanced = raw.value.matrix.clone();
        for i in 0..n {
            poll(&mut cancelled)?;
            for j in 0..n {
                balanced[i * n + j] = balanced[i * n + j].scale(columns[j] / rows[i]);
                if !finite(balanced[i * n + j]) {
                    return Err(
                        bad("balanced molecular mechanical operator is unrepresentable").into(),
                    );
                }
            }
        }
        let lu = lu_complex_with_cancel(&balanced, n, &mut cancelled)
            .map_err(BemModesError::ResponseSolve)?;
        poll(&mut cancelled)?;
        Ok(PreparedMolecularLoadedModes {
            source: self,
            exterior: exterior.clone(),
            boundary,
            minimum_clearance_m,
            near_field_options,
            radiation,
            raw,
            lu,
            rows,
            columns,
            options: solve_options,
        })
    }
}

impl PreparedMolecularLoadedModes<'_> {
    pub fn source(&self) -> &BemLoadedModes {
        self.source
    }
    pub fn radiation(&self) -> &PreparedMolecularResponse {
        &self.radiation
    }
    /// Explicit adaptive receiver choice retained for incident and motion
    /// fields alike. None keeps the original centroid receiver path.
    pub fn near_field_options(&self) -> Option<MolecularNearFieldOptions> {
        self.near_field_options
    }

    fn admit_force_batch(&self, forces: &[&[C64]]) -> Result<(), MolecularModesError> {
        let (m, n) = (self.source.mode_count(), self.rows.len());
        let r = self.radiation.harmonic().receivers().len();
        let cases = forces.len();
        let receiver_bytes = if self.near_field_options.is_some() {
            96
        } else {
            64
        };
        let bytes =
            cases as u128 * (32 * n as u128 + 64 * m as u128 + receiver_bytes * r as u128 + 512);
        let work = 32 * (n as u128).pow(3)
            + 128 * cases as u128 * (n as u128).pow(2)
            + 64 * cases as u128 * m as u128 * (r as u128 + 1);
        if cases == 0
            || cases > self.options.max_load_cases
            || forces
                .iter()
                .any(|f| f.len() != m || f.iter().any(|&z| !finite(z)))
            || bytes > self.options.max_response_bytes as u128
            || work > self.options.max_work as u128
        {
            return Err(bad("molecular loaded response needs a complete finite force batch within output/work limits").into());
        }
        Ok(())
    }

    /// Solve complete coherent force columns against the cached operator, then
    /// use each ACTUAL loaded velocity for the matched receiver and power bank.
    /// All Maxwell arms survive, including mechanically cancelling histories.
    /// Failed batches publish no partial responses; preparation is immutable.
    ///
    /// Work admission is 32*d³+128*cases*d²+64*cases*modes*(receivers+1),
    /// including matched receiver evaluation and retaining the conservative
    /// factorization allowance already paid during preparation. Response bytes
    /// additionally include source reactions and all receivers; original model
    /// payload admission bounds the temporary solve/audit buffers separately.
    /// Adaptive receiver estimates add 32*cases*receivers bytes while remaining
    /// within the existing conservative matched-observation work reservation.
    pub fn harmonic_response_batch(
        &self,
        forces: &[&[C64]],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<MolecularLoadedResponse>, MolecularModesError> {
        poll(&mut cancelled)?;
        self.admit_force_batch(forces)?;
        let (m, n) = (self.source.mode_count(), self.rows.len());
        let r = self.radiation.harmonic().receivers().len();
        let cases = forces.len();
        let omega = self.radiation.harmonic().angular_frequency_rad_s();
        let mut output = Vec::with_capacity(cases);
        for &force in forces {
            poll(&mut cancelled)?;
            let mut state = vec![C64::ZERO; n];
            for i in 0..m {
                state[i] = divide(force[i], self.rows[i]);
            }
            self.lu
                .solve_with_cancel(&mut state, &mut cancelled)
                .map_err(BemModesError::ResponseSolve)?;
            for (value, &scale) in state.iter_mut().zip(&self.columns) {
                *value = value.scale(scale);
            }
            if state.iter().any(|&z| !finite(z)) {
                return Err(bad("molecular loaded state overflow").into());
            }
            let residuals = self.raw.physical_residuals(
                m,
                0,
                &state,
                Some(force),
                self.source.settings.physical_residual_tolerance,
                &mut cancelled,
            )?;
            let velocity: Vec<_> = state[..m]
                .iter()
                .map(|&x| C64::new(0.0, -omega) * x)
                .collect();
            let mut reaction = vec![C64::ZERO; m];
            self.radiation
                .resisting_efforts(&velocity, &mut reaction, &mut cancelled)?;
            for value in &mut reaction {
                *value = -*value;
            }
            let mut receivers = vec![HarmonicField::default(); r];
            self.radiation.harmonic().evaluate_coherent(
                &velocity,
                &mut receivers,
                &mut cancelled,
            )?;
            let receiver_quadrature_error_estimates = if self.near_field_options.is_some() {
                self.radiation
                    .receiver_quadrature_for_rates(&velocity, &mut cancelled)?
            } else {
                None
            };
            let power = self.power(force, &velocity, &state[m..], &mut cancelled)?;
            let memory = state.split_off(m);
            output.push(MolecularLoadedResponse {
                angular_frequency_rad_s: omega,
                modal_displacement: state,
                modal_velocity: velocity,
                elastic_memory_sqrt_j: memory,
                fluid_reaction: reaction,
                receivers,
                receiver_quadrature_error_estimates,
                mechanical_residual: residuals[0],
                memory_residual: residuals[2],
                power,
            });
        }
        poll(&mut cancelled)?;
        Ok(output)
    }

    fn power(
        &self,
        force: &[C64],
        velocity: &[C64],
        memory: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MolecularModesPower, MolecularModesError> {
        let amplitude = force.iter().map(|z| z.abs()).fold(0.0, f64::max);
        let amplitude = if amplitude == 0.0 { 1.0 } else { amplitude };
        let normalize = |values: &[C64]| -> Result<Vec<C64>, MolecularModesError> {
            values
                .iter()
                .map(|&value| {
                    let v = divide(value, amplitude);
                    if !finite(v)
                        || (value.re != 0.0 && v.re == 0.0)
                        || (value.im != 0.0 && v.im == 0.0)
                    {
                        Err(bad("molecular work normalization is unrepresentable").into())
                    } else {
                        Ok(v)
                    }
                })
                .collect()
        };
        let f = normalize(force)?;
        let v = normalize(velocity)?;
        let e = normalize(memory)?;
        let boundary = self.radiation.power_w(&v, &mut *cancelled)?;
        let m = v.len();
        let mut sums = [0.0; 4];
        let mut corrections = [0.0; 4];
        let mut add = |channel: usize, value: f64| {
            let adjusted = value - corrections[channel];
            let next = sums[channel] + adjusted;
            corrections[channel] = (next - sums[channel]) - adjusted;
            sums[channel] = next;
        };
        for i in 0..m {
            poll(cancelled)?;
            add(0, 0.5 * (f[i] * v[i].conj()).re);
            add(3, 0.5 * f[i].abs() * v[i].abs());
            for j in 0..m {
                let c = self.source.viscosity[i * m + j];
                add(1, 0.5 * (v[i].conj() * v[j]).re * c);
                add(
                    3,
                    0.5 * v[i].abs()
                        * v[j].abs()
                        * (c.abs() + self.radiation.impedance()[i * m + j].abs()),
                );
            }
        }
        for (&value, &tau) in e.iter().zip(&self.source.memory.times) {
            poll(cancelled)?;
            let heat = 0.5 * value.norm_sq() / tau;
            add(2, heat);
            add(3, heat);
        }
        let balance = sums[0]
            - sums[1]
            - sums[2]
            - boundary.exterior_power_w
            - boundary.absorbed_wall_power_w;
        let relative = if sums[3] == 0.0 {
            balance.abs()
        } else {
            balance.abs() / sums[3]
        };
        if sums.iter().any(|v| !v.is_finite()) || !relative.is_finite() {
            return Err(bad("molecular loaded power overflow").into());
        }
        if relative > self.options.power_relative_tolerance {
            return Err(BemModesError::PowerBalance {
                relative,
                tolerance: self.options.power_relative_tolerance,
            }
            .into());
        }
        let restore = |value: f64| -> Result<f64, MolecularModesError> {
            let scaled = (value * amplitude) * amplitude;
            if !scaled.is_finite() || (value != 0.0 && scaled == 0.0) {
                Err(bad("physical molecular power is unrepresentable").into())
            } else {
                Ok(scaled)
            }
        };
        Ok(MolecularModesPower {
            input_w: restore(sums[0])?,
            kelvin_w: restore(sums[1])?,
            relaxation_w: restore(sums[2])?,
            exterior_w: restore(boundary.exterior_power_w)?,
            wall_w: restore(boundary.absorbed_wall_power_w)?,
            balance_w: restore(balance)?,
            apparent_power_scale_w: restore(sums[3])?,
            relative_balance: relative,
        })
    }
}

#[cfg(test)]
mod adaptive_tests;
#[cfg(test)]
mod sphere_tests;
#[cfg(test)]
mod tests;
