//! Driven end loads for the original uniform-flow acoustic cylinder.
//!
//! The conjugate variables are stagnation-enthalpy perturbation `H'` [J/kg]
//! and OUTWARD perturbation mass flow `m_out'` [kg/s]. At each end the declared
//! single-frequency load obeys `H' = Z m_out' + H_drive'`. For `Re(Z) >= 0`,
//! its absorbed disturbance power is `Re(Z) |m_out'|²/2`; the source supplies
//! `-Re(H_drive' conj(m_out'))/2` to the duct. Both are evaluated from the
//! solved physical fields, independently of the multiple-reflection solve.
//!
//! The characteristic impedance in these variables is `R = c/(rho A)`.
//! Writing the incoming/outgoing enthalpy characteristics as `a` and `b`
//! gives `H'=a+b`, `m_out'=(b-a)/R`, and
//! `a = (Z-R)/(Z+R) b + R/(Z+R) H_drive'` at EITHER end. This common basis
//! permits a fixed-size, power-consistent feedback solve despite different
//! upstream/downstream pressure-wave amplitudes and propagation delays.
//!
//! These are explicitly prescribed perturbation loads in a maintained mean
//! through-flow. They do not model a sealed physical cap, exit vortex shedding,
//! a stationary pressure/volume-flow radiation impedance, or an unsteady
//! mean-flow supply. A passive sample at one frequency is not a causal
//! broadband realization. Source amplitudes must remain in the original
//! small-signal regime. See the parent module's material and transverse-mode
//! admission, and Rienstra's Myers energy pairing cited there.

use super::{UniformFlowCylinder, UniformFlowExcitation, UniformFlowField, bad, finite, positive};
use crate::DuctError;
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_math::{c64::C64, det};

/// A passive harmonic stagnation-enthalpy / outward-mass-flow impedance.
///
/// Units are `(J/kg)/(kg/s) = J s/kg²`. This is NOT the stationary acoustic
/// impedance `p/Q` [Pa s/m³]. At zero mean flow only, `Z_Hm = Z_pQ/rho²`.
/// The constructor admits a single frequency sample, not its causal extension.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FlowImpedance {
    value: C64,
}

impl FlowImpedance {
    /// Declare a finite load with nonnegative real part. Zero is a prescribed
    /// enthalpy boundary. Infinite impedance is not encoded as a large number.
    pub fn try_new(impedance_j_s_kg2: C64) -> Result<Self, DuctError> {
        if !finite(impedance_j_s_kg2) || impedance_j_s_kg2.re < 0.0 {
            return Err(bad(
                "flow end impedance must be finite and have nonnegative resistance",
            ));
        }
        Ok(Self {
            value: impedance_j_s_kg2,
        })
    }

    /// Original authored impedance [J s/kg²].
    pub const fn value_j_s_kg2(self) -> C64 {
        self.value
    }
}

/// One solved end section, using the normal OUT OF the cylinder.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FlowTerminalObservation {
    /// Total stagnation-enthalpy perturbation H' [J/kg].
    pub enthalpy_j_kg: C64,
    /// Total perturbation mass flow leaving the cylinder [kg/s].
    pub outward_mass_flow_kg_s: C64,
    /// Authored in-series enthalpy source H_drive' [J/kg].
    pub source_enthalpy_j_kg: C64,
    /// Independent residual H' - Z m_out' - H_drive' [J/kg].
    pub boundary_residual_j_kg: C64,
    /// Nonnegative real-load absorption Re(Z) |m_out'|²/2 [W].
    pub absorbed_disturbance_power_w: f64,
    /// Signed work from this source INTO the cylinder [W]. A coherent source
    /// may absorb work from the other source, so this value may be negative.
    pub source_work_into_duct_w: f64,
    /// Direct physical flux Re(H' conj(m_out'))/2 OUT OF the cylinder [W].
    pub outward_disturbance_power_w: f64,
}

/// Independent terminal work account for a driven, lossless cylinder.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FlowTerminationPower {
    /// Sum of both sources' signed work into the cylinder [W].
    pub source_work_into_duct_w: f64,
    /// Sum of the two original loads' nonnegative absorption [W].
    pub absorbed_disturbance_power_w: f64,
    /// Source work minus load absorption [W]. This is a floating-point balance
    /// residual, not a heat source or a bound on unmodeled flow physics.
    pub balance_residual_w: f64,
    /// Sum of the independently observed end fluxes [W]; exactly zero for
    /// the lossless interior. It is not inferred from the load/source account.
    pub interior_flux_residual_w: f64,
}

/// Fixed-size solve for both end loads, including all multiple reflections.
///
/// Every drive reuses four complex source-to-pressure coefficients. There is
/// no heap allocation, material re-query, iteration, or wave-series truncation
/// in preparation or evaluation. The original cylinder retains all geometry,
/// material, flow, frequency, and plane-wave assumptions.
#[derive(Debug)]
pub struct PreparedFlowTerminations<'duct, 'source, M: FluidAcoustics + ?Sized> {
    duct: &'duct UniformFlowCylinder<'source, M>,
    impedances: [FlowImpedance; 2],
    reference_impedance: f64,
    reflection: [C64; 2],
    source_to_pressure: [C64; 4],
    condition_estimate: f64,
}

/// Solved excitation and independently evaluated physical end observations.
/// Its original cylinder remains borrowed even if the prepared load solve is
/// dropped. Observation does not call the material or solve the feedback again.
#[derive(Debug)]
pub struct TerminatedFlowExcitation<'duct, 'source, M: FluidAcoustics + ?Sized> {
    field: UniformFlowExcitation<'duct, 'source, M>,
    terminals: [FlowTerminalObservation; 2],
    power: FlowTerminationPower,
}

impl<'source, M: FluidAcoustics + ?Sized> UniformFlowCylinder<'source, M> {
    /// Matched enthalpy/mass-flow impedance c/(rho A) [J s/kg²].
    /// It accepts a pure outgoing characteristic at either end, including
    /// mean flow. An unrepresentable reference impedance refuses explicitly.
    pub fn characteristic_flow_impedance(&self) -> Result<FlowImpedance, DuctError> {
        let reference = (self.state.sound_speed_m_s() / self.state.density_kg_m3()) / self.area_m2;
        if !positive(reference) {
            return Err(bad("flow end characteristic impedance is unrepresentable"));
        }
        FlowImpedance::try_new(C64::from_re(reference))
    }

    /// Prepare left and right conjugate end loads for repeated coherent drives.
    ///
    /// `maximum_condition >= 1` bounds cancellation in `1-r_left*r_right*t+*t-`.
    /// It is a local feedback condition estimate, not a rigorous error bound.
    /// Undamped resonances and excessive cancellation return `DuctError::Singular`;
    /// no artificial damping or denominator regularization is introduced.
    pub fn prepare_terminations(
        &self,
        impedances: [FlowImpedance; 2],
        maximum_condition: f64,
    ) -> Result<PreparedFlowTerminations<'_, 'source, M>, DuctError> {
        if !maximum_condition.is_finite() || maximum_condition < 1.0 {
            return Err(bad(
                "flow end condition bound must be finite and at least one",
            ));
        }
        let reference = self.characteristic_flow_impedance()?.value.re;
        let mut reflection = [C64::ZERO; 2];
        let mut source_gain = [C64::ZERO; 2];
        for i in 0..2 {
            let z = impedances[i].value;
            // Scale BEFORE adding Z and R; a valid large load must not overflow
            // just because an unscaled intermediate sum would do so.
            let scale = z.re.max(z.im.abs()).max(reference);
            let r = reference / scale;
            let z = C64::new(z.re / scale, z.im / scale);
            let denominator = z + C64::from_re(r);
            reflection[i] = (z - C64::from_re(r)) / denominator;
            source_gain[i] = C64::from_re(r) / denominator;
            if !positive(r)
                || !finite(reflection[i])
                || !finite(source_gain[i])
                || source_gain[i].abs() == 0.0
            {
                return Err(bad("flow end reflection or source gain is unrepresentable"));
            }
        }
        let left_return = reflection[0] * self.transmission[1];
        let right_return = reflection[1] * self.transmission[0];
        let (denominator, condition) = crate::compose::feedback(left_return, right_return);
        if !condition.is_finite() || condition > maximum_condition {
            return Err(DuctError::Singular);
        }
        let gains = source_gain.map(|gain| gain / denominator);
        let pressure_scale = self
            .mach_factors
            .map(|factor| self.state.density_kg_m3() / factor);
        let source_to_pressure = [
            gains[0].scale(pressure_scale[0]),
            (left_return * gains[1]).scale(pressure_scale[0]),
            (right_return * gains[0]).scale(pressure_scale[1]),
            gains[1].scale(pressure_scale[1]),
        ];
        if !pressure_scale.into_iter().all(positive)
            || !source_to_pressure.into_iter().all(finite)
            || source_to_pressure[0].abs() == 0.0
            || source_to_pressure[3].abs() == 0.0
        {
            return Err(bad(
                "flow end source-to-pressure response is unrepresentable",
            ));
        }
        Ok(PreparedFlowTerminations {
            duct: self,
            impedances,
            reference_impedance: reference,
            reflection,
            source_to_pressure,
            condition_estimate: condition,
        })
    }
}

impl<'duct, 'source, M: FluidAcoustics + ?Sized> PreparedFlowTerminations<'duct, 'source, M> {
    /// Original cylinder, retaining its original material owner.
    pub const fn duct(&self) -> &'duct UniformFlowCylinder<'source, M> {
        self.duct
    }

    /// Original loads in [left, right] order.
    pub const fn impedances(&self) -> &[FlowImpedance; 2] {
        &self.impedances
    }

    /// Characteristic H'/m' reference c/(rho A) [J s/kg²].
    pub const fn reference_impedance_j_s_kg2(&self) -> f64 {
        self.reference_impedance
    }

    /// Reflected/incident ENTHALPY characteristics at [left, right]. Their
    /// squared magnitudes are disturbance-power reflection fractions. Pressure
    /// reflection adds the direction-dependent (1±M) conversion factors.
    pub const fn enthalpy_reflection(&self) -> [C64; 2] {
        self.reflection
    }

    /// Local multiple-reflection cancellation estimate from preparation.
    pub const fn condition_estimate(&self) -> f64 {
        self.condition_estimate
    }

    /// Solve coherent in-series enthalpy drives [J/kg], ordered [left, right].
    /// Both loads use the normal out of the cylinder. Zero drives produce
    /// silence; finite unrepresentable fields/work refuse without publication.
    pub fn solve(
        &self,
        source_enthalpy_j_kg: [C64; 2],
    ) -> Result<TerminatedFlowExcitation<'duct, 'source, M>, DuctError> {
        if !source_enthalpy_j_kg.into_iter().all(finite) {
            return Err(bad("flow end source enthalpies must be finite"));
        }
        let h = source_enthalpy_j_kg;
        let t = self.source_to_pressure;
        let incoming_pressure = [t[0] * h[0] + t[1] * h[1], t[2] * h[0] + t[3] * h[1]];
        let field = self.duct.prepare_excitation(incoming_pressure)?;
        let sections = [
            field.at_position(0.0)?,
            field.at_position(self.duct.length_m)?,
        ];
        let terminals = core::array::from_fn(|i| {
            let enthalpy = sections[i].stagnation_enthalpy_perturbation_j_kg;
            let mass_flow =
                sections[i]
                    .mass_flow_perturbation_kg_s
                    .scale(if i == 0 { -1.0 } else { 1.0 });
            let load = self.impedances[i].value;
            FlowTerminalObservation {
                enthalpy_j_kg: enthalpy,
                outward_mass_flow_kg_s: mass_flow,
                source_enthalpy_j_kg: h[i],
                boundary_residual_j_kg: enthalpy - load * mass_flow - h[i],
                // Scale before squaring: a small mass flow through a large
                // real impedance can dissipate representable physical power.
                absorbed_disturbance_power_w: mass_flow
                    .scale(det::sqrt(load.re) * core::f64::consts::FRAC_1_SQRT_2)
                    .norm_sq(),
                source_work_into_duct_w: -0.5 * (h[i] * mass_flow.conj()).re,
                outward_disturbance_power_w: 0.5 * (enthalpy * mass_flow.conj()).re,
            }
        });
        let source_work =
            terminals[0].source_work_into_duct_w + terminals[1].source_work_into_duct_w;
        let absorbed =
            terminals[0].absorbed_disturbance_power_w + terminals[1].absorbed_disturbance_power_w;
        let power = FlowTerminationPower {
            source_work_into_duct_w: source_work,
            absorbed_disturbance_power_w: absorbed,
            balance_residual_w: source_work - absorbed,
            interior_flux_residual_w: terminals[0].outward_disturbance_power_w
                + terminals[1].outward_disturbance_power_w,
        };
        if terminals.iter().any(|end| {
            !finite(end.boundary_residual_j_kg)
                || ![
                    end.absorbed_disturbance_power_w,
                    end.source_work_into_duct_w,
                    end.outward_disturbance_power_w,
                ]
                .into_iter()
                .all(f64::is_finite)
        }) || ![
            power.source_work_into_duct_w,
            power.absorbed_disturbance_power_w,
            power.balance_residual_w,
            power.interior_flux_residual_w,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return Err(bad("flow end boundary work or residual is unrepresentable"));
        }
        Ok(TerminatedFlowExcitation {
            field,
            terminals,
            power,
        })
    }
}

impl<'duct, 'source, M: FluidAcoustics + ?Sized> TerminatedFlowExcitation<'duct, 'source, M> {
    /// Solved original characteristic field, including actual incident and
    /// outgoing disturbance powers rather than the externally authored drives.
    pub const fn field(&self) -> &UniformFlowExcitation<'duct, 'source, M> {
        &self.field
    }

    /// Physical observations and boundary-law residuals, [left, right].
    pub const fn terminals(&self) -> &[FlowTerminalObservation; 2] {
        &self.terminals
    }

    /// Independent total source work, load absorption and conservation residuals.
    pub const fn power(&self) -> &FlowTerminationPower {
        &self.power
    }

    /// Observe the original pressure, velocity and disturbance energy at x [m].
    pub fn at_position(&self, x_m: f64) -> Result<UniformFlowField, DuctError> {
        self.field.at_position(x_m)
    }
}

#[cfg(test)]
mod tests;
