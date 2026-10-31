//! Nonlinear monopole interactions between restrained spherical bubbles.
//!
//! Every member retains its original Rayleigh--Plesset owner. For
//! S_i = d(R_i² R_i')/dt, the joint incompressible acceleration solve is
//!
//! S_i/R_i + sum_(j!=i) S_j/d_ij = R_i a_isolated,i + 2 R_i'^2.
//!
//! The radii, gas pressures, viscosity and capillarity are fully nonlinear.
//! This is the leading separated-sphere potential-flow model, not a model of
//! translating bubbles, jets, close contact, cloud shocks or acoustic radiation
//! loss. Centers are externally restrained; their velocities and constraint
//! work are zero. No translational force or momentum-conservation claim is made.
//!
//! Source: Harkin, Kaper & Nadim (2001), DOI 10.1017/S0022112001005857;
//! the monopole pressure rho*d(R²R')/dt/d is also Eq. 1 of Yuan, Sankin & Zhong
//! (2011), https://pmc.ncbi.nlm.nih.gov/articles/PMC3248064/.
//! No finite-compressibility pressure-derivative terms are mixed into this law.
//!
//! The SAME symmetric inertia gives the liquid kinetic energy
//! T = rho/(8pi) [sum Q_i²/R_i + 2 sum_(i<j) Q_i Q_j/d_ij], Q_i=V_i'.
//! Original gas volume work, surface energy, and (p_ref-p_v)V complete the
//! mechanical storage. Authored pressure work and native viscous loss are
//! integrated on the same joint DP54 trajectory and checked before publication.

use super::{
    BubbleNumerics, BubbleObservation, BubbleState, CavitationError, RayleighPlessetBubble, finite,
    poll, positive, stepper,
};
use fs_material::cavitation::CavitationState;

#[path = "cluster/dynamics.rs"]
mod dynamics;
use dynamics::Workspace;

/// Fixed bound on both the joint nonlinear solve and fixed-stack integration.
pub const MAX_INTERACTING_BUBBLES: usize = 16;
const WORK: usize = 2 * MAX_INTERACTING_BUBBLES;
const LOSS: usize = WORK + 1;
const DIMENSION: usize = LOSS + 1;

/// An original physical bubble at a fixed center in an unbounded liquid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClusterBubble {
    pub bubble: RayleighPlessetBubble,
    pub center_m: [f64; 3],
}

/// Explicit asymptotic domain and bounded system size.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BubbleClusterLimits {
    /// Admitted number of original bubbles, in 1..=16.
    pub maximum_bubbles: usize,
    /// max(R_max,i,R_max,j)/d_ij ceiling, in (0,0.1].
    pub maximum_radius_to_separation: f64,
    /// Row sum of sqrt(R_max,i R_max,j)/d_ij, in (0,0.5].
    /// The normalized symmetric inertia then has eigenvalues in [1-s,1+s].
    pub maximum_inertia_row_sum: f64,
    /// Solved RP wall-pressure jump divided by rho*c², in (0,0.01].
    /// Low wall Mach alone does not admit a large initial pressure imbalance.
    pub maximum_liquid_pressure_ratio: f64,
    /// CALLER-DECLARED bandwidth of the actual nonlinear volume sources [rad/s],
    /// including transients/harmonics, not merely the imposed-pressure drive.
    /// The solver cannot certify that a waveform satisfies this assumption.
    pub maximum_source_angular_frequency_rad_s: f64,
    /// omega_source * max(cluster diameter, R_max)/c ceiling, in (0,0.1].
    /// It admits negligible propagation delay within the retained approximation.
    pub maximum_retardation_phase_rad: f64,
}

/// Original radial accuracy/work controls plus the mechanical energy scale.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BubbleClusterNumerics {
    pub radial: BubbleNumerics,
    /// Positive tolerance for each integrated work/loss coordinate and the
    /// endpoint energy balance. Relative tolerance is the original radial value.
    pub absolute_energy_tolerance_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum BubbleClusterError {
    InvalidInput(&'static str),
    Native(CavitationError),
    Member {
        index: usize,
        cause: CavitationError,
    },
    MaterialMismatch {
        index: usize,
    },
    Geometry {
        first: usize,
        second: usize,
        reason: &'static str,
    },
    OutsideDomain(&'static str),
    LinearSolve,
    EnergyBalance {
        residual_j: f64,
        tolerance_j: f64,
    },
    Cancelled,
}
impl From<CavitationError> for BubbleClusterError {
    fn from(error: CavitationError) -> Self {
        if error == CavitationError::Cancelled {
            Self::Cancelled
        } else {
            Self::Native(error)
        }
    }
}
impl core::fmt::Display for BubbleClusterError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "interacting bubble cluster refused: {self:?}")
    }
}
impl std::error::Error for BubbleClusterError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Native(error) | Self::Member { cause: error, .. } => Some(error),
            _ => None,
        }
    }
}

/// Coupled sources and mechanical storage at one common accepted time.
#[derive(Clone, Debug, PartialEq)]
pub struct BubbleClusterObservation {
    pub time_s: f64,
    pub external_absolute_pressure_pa: f64,
    /// Each observation uses its SOLVED acceleration and the external pressure
    /// plus neighbor monopole pressure. Its kinetic-energy field remains that
    /// bubble's self contribution; the aggregate below includes cross terms.
    pub bubbles: Vec<BubbleObservation>,
    pub total_volume_m3: f64,
    pub total_volume_rate_m3_s: f64,
    pub total_volume_acceleration_m3_s2: f64,
    pub liquid_kinetic_energy_j: f64,
    pub bubble_potential_j: f64,
    pub total_mechanical_energy_j: f64,
    pub viscous_dissipation_w: f64,
    /// -(p_external-p_reference) sum V_i' [W], signed work into this storage.
    pub drive_power_w: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BubbleClusterStep {
    pub observation: BubbleClusterObservation,
    pub pressure_work_j: f64,
    pub viscous_loss_j: f64,
    pub energy_balance_residual_j: f64,
    pub accepted_steps: usize,
    pub rejected_steps: usize,
    pub rhs_evaluations: usize,
    pub maximum_normalized_local_error: f64,
}

/// Fixed original members and an atomically published joint physical state.
#[derive(Debug)]
pub struct BubbleCluster {
    members: Vec<ClusterBubble>,
    states: Vec<BubbleState>,
    material: CavitationState,
    time_s: f64,
    reference_pressure_pa: f64,
    limits: BubbleClusterLimits,
    inverse_distances: Vec<f64>,
    maximum_inertia_row_sum: f64,
    maximum_retardation_phase_rad: f64,
}

impl BubbleCluster {
    /// Initial states must have the exact supplied common clock. No equilibrium
    /// radius, gas inventory, pressure or velocity is manufactured or projected.
    pub fn try_new(
        members: Vec<ClusterBubble>,
        states: Vec<BubbleState>,
        time_s: f64,
        reference_pressure_pa: f64,
        limits: BubbleClusterLimits,
        initial_absolute_pressure_pa: f64,
    ) -> Result<Self, BubbleClusterError> {
        let count = members.len();
        if !(1..=MAX_INTERACTING_BUBBLES).contains(&limits.maximum_bubbles)
            || count == 0
            || count > limits.maximum_bubbles
            || states.len() != count
            || !time_s.is_finite()
            || !reference_pressure_pa.is_finite()
            || !positive(limits.maximum_radius_to_separation)
            || limits.maximum_radius_to_separation > 0.1
            || !positive(limits.maximum_inertia_row_sum)
            || limits.maximum_inertia_row_sum > 0.5
            || !positive(limits.maximum_liquid_pressure_ratio)
            || limits.maximum_liquid_pressure_ratio > 0.01
            || !positive(limits.maximum_source_angular_frequency_rad_s)
            || !positive(limits.maximum_retardation_phase_rad)
            || limits.maximum_retardation_phase_rad > 0.1
        {
            return Err(BubbleClusterError::InvalidInput(
                "cluster size, clock or domain controls",
            ));
        }
        let material = members[0].bubble.liquid();
        let mut extent: f64 = 0.0;
        for (index, (member, state)) in members.iter().zip(&states).enumerate() {
            if member.bubble.liquid() != material {
                return Err(BubbleClusterError::MaterialMismatch { index });
            }
            if member
                .center_m
                .into_iter()
                .any(|coordinate| !coordinate.is_finite())
                || state.time_s() != time_s
            {
                return Err(BubbleClusterError::InvalidInput(
                    "finite centers and exact shared clock",
                ));
            }
            extent = extent.max(member.bubble.limits().maximum_radius_m);
        }
        let mut inverse_distances = vec![0.0; count * count];
        let mut row_sums = [0.0; MAX_INTERACTING_BUBBLES];
        for i in 0..count {
            for j in i + 1..count {
                let a = members[i].center_m;
                let b = members[j].center_m;
                let distance = (a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2]);
                let ri = members[i].bubble.limits().maximum_radius_m;
                let rj = members[j].bubble.limits().maximum_radius_m;
                if !positive(distance)
                    || !positive(1.0 / distance)
                    || ri + rj >= distance
                    || ri.max(rj) / distance > limits.maximum_radius_to_separation
                {
                    return Err(BubbleClusterError::Geometry {
                        first: i,
                        second: j,
                        reason: "maximum-radius spheres are not well separated",
                    });
                }
                inverse_distances[i * count + j] = 1.0 / distance;
                inverse_distances[j * count + i] = 1.0 / distance;
                let coupling = (ri.sqrt() * rj.sqrt()) / distance;
                if !positive(coupling) {
                    return Err(BubbleClusterError::Geometry {
                        first: i,
                        second: j,
                        reason: "interaction is not representable",
                    });
                }
                row_sums[i] += coupling;
                row_sums[j] += coupling;
                extent = extent.max(distance);
            }
        }
        let maximum_inertia_row_sum = row_sums[..count].iter().copied().fold(0.0, f64::max);
        if !maximum_inertia_row_sum.is_finite()
            || maximum_inertia_row_sum > limits.maximum_inertia_row_sum
        {
            return Err(BubbleClusterError::OutsideDomain(
                "collective inertia coupling exceeds bound",
            ));
        }
        let maximum_retardation_phase_rad = (extent / material.liquid_sound_speed_m_s())
            * limits.maximum_source_angular_frequency_rad_s;
        if !maximum_retardation_phase_rad.is_finite()
            || maximum_retardation_phase_rad > limits.maximum_retardation_phase_rad
        {
            return Err(BubbleClusterError::OutsideDomain(
                "declared source bandwidth is not long-wave",
            ));
        }
        let owner = Self {
            members,
            states,
            material,
            time_s,
            reference_pressure_pa,
            limits,
            inverse_distances,
            maximum_inertia_row_sum,
            maximum_retardation_phase_rad,
        };
        owner.observe(initial_absolute_pressure_pa)?;
        Ok(owner)
    }

    pub fn members(&self) -> &[ClusterBubble] {
        &self.members
    }
    pub fn states(&self) -> &[BubbleState] {
        &self.states
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn material(&self) -> CavitationState {
        self.material
    }
    pub const fn reference_pressure_pa(&self) -> f64 {
        self.reference_pressure_pa
    }
    pub const fn limits(&self) -> BubbleClusterLimits {
        self.limits
    }
    pub const fn maximum_inertia_row_sum(&self) -> f64 {
        self.maximum_inertia_row_sum
    }
    pub const fn maximum_retardation_phase_rad(&self) -> f64 {
        self.maximum_retardation_phase_rad
    }

    pub fn observe(
        &self,
        pressure_pa: f64,
    ) -> Result<BubbleClusterObservation, BubbleClusterError> {
        self.observe_values(
            self.time_s,
            &self.values(),
            pressure_pa,
            &mut Workspace::new(self.members.len())?,
        )
    }

    pub fn advance<Pressure: Fn(f64) -> f64>(
        &mut self,
        end_time_s: f64,
        numerics: &BubbleClusterNumerics,
        pressure: &Pressure,
    ) -> Result<BubbleClusterStep, BubbleClusterError> {
        self.advance_controlled(end_time_s, numerics, pressure, &mut || false)
    }

    /// One jointly accepted trajectory; every error and cancellation leaves all
    /// member states and the common clock unchanged. The deterministic pressure
    /// callback may revisit times after rejected trials. The native maximum step
    /// must resolve the drive; it does not certify the declared SOURCE bandwidth.
    pub fn advance_controlled<Pressure, Cancel>(
        &mut self,
        end_time_s: f64,
        numerics: &BubbleClusterNumerics,
        pressure: &Pressure,
        cancelled: &mut Cancel,
    ) -> Result<BubbleClusterStep, BubbleClusterError>
    where
        Pressure: Fn(f64) -> f64,
        Cancel: FnMut() -> bool,
    {
        poll(cancelled)?;
        numerics.radial.validate()?;
        if !positive(numerics.absolute_energy_tolerance_j) {
            return Err(BubbleClusterError::InvalidInput(
                "positive mechanical energy tolerance",
            ));
        }
        let initial = self.values();
        let mut workspace = Workspace::new(self.members.len())?;
        let before =
            self.observe_values(self.time_s, &initial, pressure(self.time_s), &mut workspace)?;
        let mut absolute = [1.0; DIMENSION];
        for i in 0..self.members.len() {
            absolute[2 * i] = numerics.radial.radius_absolute_tolerance_m;
            absolute[2 * i + 1] = numerics.radial.velocity_absolute_tolerance_m_s;
        }
        absolute[WORK] = numerics.absolute_energy_tolerance_j;
        absolute[LOSS] = numerics.absolute_energy_tolerance_j;
        let mut rhs = |time, values: &[f64; DIMENSION]| {
            let frame = self.evaluate(time, values, pressure(time), &mut workspace)?;
            let mut result = [0.0; DIMENSION];
            for i in 0..self.members.len() {
                result[2 * i] = values[2 * i + 1];
                result[2 * i + 1] = frame.accelerations[i];
            }
            result[WORK] = frame.drive_power_w;
            result[LOSS] = frame.viscous_dissipation_w;
            Ok::<_, BubbleClusterError>(result)
        };
        let result = stepper::advance_array(
            self.time_s,
            initial,
            end_time_s,
            &numerics.radial,
            absolute,
            &mut rhs,
            cancelled,
            &refinable,
        )?;
        poll(cancelled)?;
        let observation = self.observe_values(
            result.time_s,
            &result.values,
            pressure(result.time_s),
            &mut workspace,
        )?;
        let pressure_work_j = result.values[WORK];
        let viscous_loss_j = result.values[LOSS];
        let residual = observation.total_mechanical_energy_j
            - before.total_mechanical_energy_j
            - pressure_work_j
            + viscous_loss_j;
        let tolerance = numerics.absolute_energy_tolerance_j
            + numerics.radial.relative_tolerance
                * (observation.total_mechanical_energy_j.abs()
                    + before.total_mechanical_energy_j.abs()
                    + pressure_work_j.abs()
                    + viscous_loss_j.abs());
        if !residual.is_finite()
            || !tolerance.is_finite()
            || !viscous_loss_j.is_finite()
            || viscous_loss_j < 0.0
            || residual.abs() > tolerance
        {
            return Err(BubbleClusterError::EnergyBalance {
                residual_j: residual,
                tolerance_j: tolerance,
            });
        }
        poll(cancelled)?;
        for (state, current) in self.states.iter_mut().zip(&observation.bubbles) {
            *state = current.state;
        }
        self.time_s = result.time_s;
        Ok(BubbleClusterStep {
            observation,
            pressure_work_j,
            viscous_loss_j,
            energy_balance_residual_j: residual,
            accepted_steps: result.accepted_steps,
            rejected_steps: result.rejected_steps,
            rhs_evaluations: result.rhs_evaluations,
            maximum_normalized_local_error: result.maximum_normalized_local_error,
        })
    }

    fn values(&self) -> [f64; DIMENSION] {
        let mut result = [0.0; DIMENSION];
        for (i, state) in self.states.iter().enumerate() {
            result[2 * i] = state.radius_m();
            result[2 * i + 1] = state.wall_velocity_m_s();
        }
        result
    }
}

fn member(index: usize, cause: CavitationError) -> BubbleClusterError {
    BubbleClusterError::Member { index, cause }
}
fn refinable(error: &BubbleClusterError) -> bool {
    match error {
        BubbleClusterError::Native(cause) | BubbleClusterError::Member { cause, .. } => {
            stepper::refinable(cause)
        }
        _ => false,
    }
}

#[cfg(test)]
#[path = "cluster/tests.rs"]
mod tests;
