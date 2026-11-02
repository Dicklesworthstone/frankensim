//! Loading and sound from the same original signed modal velocity sources.
//!
//! The native Galerkin owner retains its P1 pressure/P0 auxiliary densities.
//! Its adaptive receiver integrates those densities directly; panel pressure
//! means are used only for conjugate mouth work. All sources share one boundary
//! factor and one receiver row bank at a frequency. Coherent source cases then
//! reuse the cached physical fields, including higher modes with zero net flow.

use super::{BemModalMouth, ModalMouthImpedance};
use fs_bem::helmholtz::{HelmholtzError, Medium};
use fs_bem::near_field::{Observation, Options, VectorObservation};
use fs_bem::reciprocal_radiation::{Budget, Prepared, ReceiverGeometry, VelocityPortSolution};
use fs_duct::DuctError;
use fs_exec::CancelGate;
use fs_math::{c64::C64, dd::Dd};

/// Separate limits for retained source/field payload and receiver contractions.
/// Native Galerkin storage/work and adaptive quadrature keep their own budgets.
#[derive(Clone, Copy, Debug)]
pub struct ModalMouthObservationBudget {
    pub maximum_receivers: usize,
    /// Conservative source solutions, receiver rows, field bank and output
    /// envelope. Excludes borrowed geometry and the separately bounded native
    /// Galerkin geometry/factor workspace. No receiver or mode is discarded.
    pub maximum_observation_bytes: usize,
    /// Conservative scalar terms for projecting all source fields and one
    /// coherent observation, excluding native solve/quadrature work.
    pub maximum_contraction_terms: u64,
}

#[derive(Debug)]
pub enum ModalMouthObservationError {
    Invalid(&'static str),
    Source(HelmholtzError),
    Load(DuctError),
    Cancelled,
}
impl core::fmt::Display for ModalMouthObservationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "modal mouth field: {self:?}")
    }
}
impl std::error::Error for ModalMouthObservationError {}
type Result<T> = core::result::Result<T, ModalMouthObservationError>;
fn invalid(message: &'static str) -> ModalMouthObservationError {
    ModalMouthObservationError::Invalid(message)
}
fn source(error: HelmholtzError) -> ModalMouthObservationError {
    match error {
        HelmholtzError::Cancelled => ModalMouthObservationError::Cancelled,
        other => ModalMouthObservationError::Source(other),
    }
}
fn poll(gate: &CancelGate) -> Result<()> {
    if gate.is_requested() {
        Err(ModalMouthObservationError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

/// Original modal mouth and one reusable native Galerkin geometry preparation.
/// The native closed-surface resonance-exclusion and mesh/angular limits apply.
pub struct PreparedModalMouthField<'mouth, 'surface> {
    mouth: &'mouth BemModalMouth<'surface>,
    boundary: Prepared<'surface>,
    highest_omega_rad_s: f64,
    medium: Medium,
}

impl<'surface> BemModalMouth<'surface> {
    /// Prepare matched loading/observation over an explicitly bounded band.
    /// This selects the native Galerkin route explicitly; no route fallback,
    /// remeshing, trace normalization or continuum-accuracy claim is made.
    pub fn prepare_galerkin_field<'mouth>(
        &'mouth self,
        highest_omega_rad_s: f64,
        budget: Budget,
        gate: &CancelGate,
    ) -> Result<PreparedModalMouthField<'mouth, 'surface>> {
        poll(gate)?;
        let n = self.basis.modes().0 as u128;
        if !highest_omega_rad_s.is_finite()
            || highest_omega_rad_s <= 0.0
            || 2 * n * n > self.budget.maximum_matrix_entries as u128
            || 2 * self.surface.areas().len() as u128 * n * n
                > self.budget.maximum_projection_terms as u128
        {
            return Err(invalid(
                "complete modal matrix and positive source band required",
            ));
        }
        let medium = Medium {
            density: self.exterior.density,
            sound_speed: self.exterior.sound_speed,
        };
        let boundary = Prepared::new(
            self.surface,
            highest_omega_rad_s / medium.sound_speed,
            budget,
            || gate.is_requested(),
        )
        .map_err(source)?;
        // Admit the complete physical source batch before any frequency solve.
        if boundary.velocity_port_work(n as usize).map_err(source)? > budget.maximum_work {
            return Err(invalid(
                "complete modal source batch exceeds native work budget",
            ));
        }
        poll(gate)?;
        Ok(PreparedModalMouthField {
            mouth: self,
            boundary,
            highest_omega_rad_s,
            medium,
        })
    }
}

/// Original exterior receivers, reusable over the prepared source band.
pub struct ModalMouthReceivers<'prepared, 'mouth, 'surface> {
    source: &'prepared PreparedModalMouthField<'mouth, 'surface>,
    geometry: ReceiverGeometry<'prepared, 'surface>,
    positions_m: Vec<[f64; 3]>,
}
impl<'mouth, 'surface> PreparedModalMouthField<'mouth, 'surface> {
    pub fn mouth(&self) -> &BemModalMouth<'surface> {
        self.mouth
    }
    pub fn highest_omega_rad_s(&self) -> f64 {
        self.highest_omega_rad_s
    }
    pub fn medium(&self) -> Medium {
        self.medium
    }

    /// Admit all original receiver positions and their exact surface clearance.
    /// Payload/work refusal precedes geometry and receiver-row allocation.
    pub fn receivers<'prepared>(
        &'prepared self,
        positions_m: &[[f64; 3]],
        minimum_clearance_m: f64,
        budget: ModalMouthObservationBudget,
        gate: &CancelGate,
    ) -> Result<ModalMouthReceivers<'prepared, 'mouth, 'surface>> {
        poll(gate)?;
        let r = positions_m.len() as u128;
        let n = self.mouth.basis.modes().0 as u128;
        let panels = self.mouth.surface.areas().len() as u128;
        // Includes native adaptive P1 row banks and componentwise discrepancies,
        // original velocity/pressure/density solutions, and complete outputs.
        let bytes = 2048 * r * panels
            + 256 * (n * r + n + r)
            + 96 * n * self.boundary.unknowns() as u128
            + 32 * n * n;
        let terms = 256 * r * n * panels + 64 * r * n + 16 * n * n;
        if r == 0
            || r > budget.maximum_receivers as u128
            || bytes > budget.maximum_observation_bytes as u128
            || terms > budget.maximum_contraction_terms as u128
        {
            return Err(invalid(
                "complete source/receiver field exceeds its declared limits",
            ));
        }
        let geometry = self
            .boundary
            .receivers(positions_m, minimum_clearance_m)
            .map_err(source)?;
        poll(gate)?;
        Ok(ModalMouthReceivers {
            source: self,
            geometry,
            positions_m: positions_m.to_vec(),
        })
    }
}

/// Complete work-conjugate mouth matrix and matched receiver transfer bank.
/// Original source densities and solve/flux diagnostics remain accessible.
pub struct ModalMouthFieldResponse<'prepared, 'surface> {
    load: ModalMouthImpedance,
    source: VelocityPortSolution<'prepared, 'surface>,
    positions_m: Vec<[f64; 3]>,
    unit_fields: Vec<VectorObservation>,
    kernel_evaluations: usize,
}
impl<'prepared, 'mouth, 'surface> ModalMouthReceivers<'prepared, 'mouth, 'surface> {
    pub fn positions_m(&self) -> &[[f64; 3]] {
        &self.positions_m
    }
    pub fn clearances_m(&self) -> &[f64] {
        self.geometry.clearances_m()
    }

    /// One shared source factor and one adaptive receiver-row preparation.
    /// Pressure AND particle velocity use the original P1/P0 Green field; each
    /// derivative retains the native componentwise quadrature admission.
    pub fn solve(
        &self,
        omega_rad_s: f64,
        quadrature: Options,
        gate: &CancelGate,
    ) -> Result<ModalMouthFieldResponse<'prepared, 'surface>> {
        poll(gate)?;
        quadrature.validate().map_err(source)?;
        if !omega_rad_s.is_finite()
            || omega_rad_s <= 0.0
            || omega_rad_s > self.source.highest_omega_rad_s
        {
            return Err(invalid(
                "frequency lies outside the original prepared source band",
            ));
        }
        let native: &'prepared Prepared<'surface> = &self.source.boundary;
        let mouth = self.source.mouth;
        let shapes: Vec<_> = mouth.apertures.iter().map(|a| a.shape.as_slice()).collect();
        let k = omega_rad_s / self.source.medium.sound_speed;
        let original = native
            .solve_velocity_ports(k, self.source.medium, &shapes, || gate.is_requested())
            .map_err(source)?;
        let load = ModalMouthImpedance::new(
            omega_rad_s,
            mouth.patch.radius_m,
            mouth.exterior,
            mouth.basis.modes(),
            original.impedance().to_vec(),
        )
        .map_err(ModalMouthObservationError::Load)?;
        let rows = self
            .geometry
            .prepare(k, self.source.medium, quadrature, true, || {
                gate.is_requested()
            })
            .map_err(source)?;
        let mut unit_fields = Vec::with_capacity(shapes.len());
        for solution in original.solutions() {
            poll(gate)?;
            unit_fields.push(rows.evaluate_velocity(solution).map_err(source)?);
        }
        poll(gate)?;
        Ok(ModalMouthFieldResponse {
            load,
            source: original,
            positions_m: self.positions_m.clone(),
            unit_fields,
            kernel_evaluations: rows.kernel_evaluations(),
        })
    }
}

/// One coherent modal excitation at the response's frequency and reference plane.
pub struct ModalMouthObservation {
    /// Physical p/u phasors under exp(-i omega t), with propagated absolute
    /// quadrature discrepancy estimates. These omit mesh/solve/model error.
    pub field: VectorObservation,
    /// Complete P=Z U at the physical mouth, in the original mode ordering.
    pub modal_pressure_pa: Vec<C64>,
    /// (1/2) sum P_n conj(U_n). Real part is outward mean work [W]. The signed
    /// imaginary part retains the exp(-i omega t) reactive-power convention.
    pub complex_mouth_power_w: C64,
    /// (1/2) Re(p conj(u)) in the original Cartesian frame [W/m²].
    pub mean_intensity_w_m2: Vec<[f64; 3]>,
}

fn sum_complex(values: impl Iterator<Item = C64>) -> C64 {
    let (mut re, mut im) = (Dd::ZERO, Dd::ZERO);
    for value in values {
        re = re + Dd::from_f64(value.re);
        im = im + Dd::from_f64(value.im);
    }
    C64::new(re.to_f64(), im.to_f64())
}
impl<'prepared, 'surface> ModalMouthFieldResponse<'prepared, 'surface> {
    pub fn load(&self) -> &ModalMouthImpedance {
        &self.load
    }
    pub fn source(&self) -> &VelocityPortSolution<'prepared, 'surface> {
        &self.source
    }
    pub fn positions_m(&self) -> &[[f64; 3]] {
        &self.positions_m
    }
    pub fn unit_fields(&self) -> &[VectorObservation] {
        &self.unit_fields
    }
    /// Cold quadrature work. Repeated coherent observations add no quadrature.
    pub fn kernel_evaluations(&self) -> usize {
        self.kernel_evaluations
    }

    /// Apply complete complex modal volume flows [m³/s], retaining interference.
    /// Use the same load in the interior solve that supplied these mouth flows.
    /// No BEM solve, factorization or receiver quadrature occurs here.
    pub fn observe(&self, flows: &[C64], gate: &CancelGate) -> Result<ModalMouthObservation> {
        poll(gate)?;
        let n = self.unit_fields.len();
        if flows.len() != n || flows.iter().any(|&v| !finite(v)) {
            return Err(invalid(
                "one finite complex flow per original mouth mode required",
            ));
        }
        let modal_pressure_pa: Vec<_> = self
            .load
            .impedance()
            .chunks_exact(n)
            .map(|row| sum_complex(row.iter().zip(flows).map(|(&z, &u)| z * u)))
            .collect();
        let complex_mouth_power_w = sum_complex(
            modal_pressure_pa
                .iter()
                .zip(flows)
                .map(|(&p, &u)| (p * u.conj()).scale(0.5)),
        );
        let r = self.positions_m.len();
        let mut pressure = Vec::with_capacity(r);
        let mut pressure_errors = Vec::with_capacity(r);
        let mut velocity = Vec::with_capacity(r);
        let mut velocity_errors = Vec::with_capacity(r);
        let mut mean_intensity_w_m2 = Vec::with_capacity(r);
        for receiver in 0..r {
            poll(gate)?;
            let p = sum_complex(
                self.unit_fields
                    .iter()
                    .zip(flows)
                    .map(|(f, &u)| f.scalar.pressure[receiver] * u),
            );
            let v = std::array::from_fn(|c| {
                sum_complex(
                    self.unit_fields
                        .iter()
                        .zip(flows)
                        .map(|(f, &u)| f.particle_velocity_m_s[receiver][c] * u),
                )
            });
            let pe: f64 = self
                .unit_fields
                .iter()
                .zip(flows)
                .map(|(f, u)| f.scalar.quadrature_error_estimate_pa[receiver] * u.abs())
                .sum();
            let ve: [f64; 3] = std::array::from_fn(|c| {
                self.unit_fields
                    .iter()
                    .zip(flows)
                    .map(|(f, u)| f.quadrature_error_estimate_m_s[receiver][c] * u.abs())
                    .sum()
            });
            let intensity: [f64; 3] = v.map(|u: C64| 0.5 * (p * u.conj()).re);
            if !finite(p)
                || v.iter().any(|&u| !finite(u))
                || !pe.is_finite()
                || ve.iter().chain(&intensity).any(|v| !v.is_finite())
            {
                return Err(invalid(
                    "physical receiver field, discrepancy or intensity overflow",
                ));
            }
            pressure.push(p);
            pressure_errors.push(pe);
            velocity.push(v);
            velocity_errors.push(ve);
            mean_intensity_w_m2.push(intensity);
        }
        if !finite(complex_mouth_power_w) || modal_pressure_pa.iter().any(|&p| !finite(p)) {
            return Err(invalid("physical modal mouth work overflow"));
        }
        poll(gate)?;
        Ok(ModalMouthObservation {
            field: VectorObservation {
                scalar: Observation {
                    pressure,
                    quadrature_error_estimate_pa: pressure_errors,
                },
                particle_velocity_m_s: velocity,
                quadrature_error_estimate_m_s: velocity_errors,
            },
            modal_pressure_pa,
            complex_mouth_power_w,
            mean_intensity_w_m2,
        })
    }
}

#[cfg(test)]
mod tests;
