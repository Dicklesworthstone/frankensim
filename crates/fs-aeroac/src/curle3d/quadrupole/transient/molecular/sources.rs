//! Causal, coherent molecular propagation of original volume, force and stress.
//!
//! With G(s,r)=exp(-kappa*r)/(4*pi*r), the volume and force contributions are
//! p=rho*s*Q*G - F.grad(G), u=-Q*grad(G) + Hessian(G)*F/(rho*s).
//! This is the original compact-source Green operator, with
//! kappa=s*sqrt(rho/K(s)) from the SAME material as the stress propagator.
//! See de Hoop, *Handbook of Radiation and Scattering of Waves*, chapter 5:
//! <https://www.atdehoop.com/Books/DeHoopAT_book_2008/DeHoopAT1995Handbook_chap05.pdf>.
//!
//! Extracting r/c_inf leaves only the existing F_0, F_1 and F_2 radial symbols
//! for volume and force; stress also uses F_3. All mechanisms at one path are
//! combined before those filters, then sources/images combine before p*u.
//! There is one material continuation, contour, kernel bank and convolution
//! engine. No plane-wave velocity or pressure-only attenuation is substituted.

use super::super::{PreparedTransientQuadrupoleObserver, Transfer, TransientQuadrupoleSource};
use super::{PreparedMolecularQuadrupoleObserver, invalid, nonfinite, observe, poll};
use crate::AeroacError;
use crate::curle3d::moving::CompactSourceStrength;
use crate::material::MaterialAeroacError;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::relaxation::RelaxationState;
use fs_math::dd::Dd;

pub use super::super::RetardedQuadrupoleField as RetardedSourceField;
pub use super::{
    MolecularQuadrupoleRecord as MolecularSourceRecord, MolecularTransientBudget,
    MolecularTransientConfig,
};
pub use crate::curle3d::quadrupole::recorded::QuadrupoleSampleClock as MolecularSampleClock;

/// Original accepted samples with explicit finite C2 history reconstruction.
pub mod history;
/// Complete causal returns between two opposing ideal acoustic planes.
pub mod slab;

const SCOPE: &str = "Initially resting prescribed stationary compact 3D volume, force and \
stress sources in one homogeneous molecularly relaxing gas. Complete coherent pressure \
and vector velocity retain exact ideal-image paths and explicit frozen-speed fronts. \
Finite-prefix BDF2 convolution quadrature has separate timestep, contour and original \
history errors. No classical viscothermal loss, motion, finite-impedance scattering, \
radiation reaction, source/bandwidth qualification or measured-SPL authority.";

/// All physical mechanisms at one source position and SAME emission event.
/// Separate positions remain separate source indices; set absent mechanisms
/// to zero. Values, derivatives and integrals must come from one consistent
/// original history in the medium's fixed Cartesian frame.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TransientAcousticSource {
    /// Original volume injection [m³/s], its rate [m³/s²], force ON the fluid
    /// [N] and its fixed-frame rate [N/s]. No extra density multiplies force.
    pub compact: CompactSourceStrength,
    /// Original causal force integral [N s], retaining its integration constant
    /// across every call/block. This supplies the transverse dipole near field.
    pub force_impulse_n_s: [f64; 3],
    /// Original symmetric stress moment, rate, acceleration and integral.
    /// Tensor order is xx, yy, zz, xy, xz, yz, with physical shear entries.
    pub stress: TransientQuadrupoleSource,
}

fn finite_compact(source: &TransientAcousticSource) -> bool {
    [
        source.compact.volume_flow_m3_s,
        source.compact.volume_flow_rate_m3_s2,
    ]
    .iter()
    .chain(&source.compact.force_on_fluid_n)
    .chain(&source.compact.force_rate_on_fluid_n_s)
    .chain(&source.force_impulse_n_s)
    .all(|v| v.is_finite())
}

fn initially_resting(source: &TransientAcousticSource) -> bool {
    finite_compact(source)
        && source.compact.volume_flow_m3_s == 0.0
        && source
            .compact
            .force_on_fluid_n
            .iter()
            .chain(&source.force_impulse_n_s)
            .all(|&v| v == 0.0)
        && observe::compatible_initial(&source.stress)
}

/// The six distinct volume/force scales, including the original image sign.
/// These complement, rather than duplicate, the stress owner's path geometry.
#[derive(Clone, Copy)]
struct CompactTransfer {
    volume_pressure: f64,
    near: f64,
    far: f64,
    velocity: [f64; 3],
}

impl CompactTransfer {
    fn new(path: &Transfer, medium: FluidAcousticState) -> Result<Self, AeroacError> {
        let c = medium.sound_speed_m_s();
        let density = medium.density_kg_m3();
        let range = path.delay_s * c;
        // Native direct/image coefficients already carry the pressure character.
        let sign = path.pressure_scales[0].signum();
        let green = ((0.25 / core::f64::consts::PI) / range) * sign;
        let volume_pressure = density * green;
        let near = green / range;
        let far = green / c;
        let velocity = [
            (near / range) / density,
            (near / c) / density,
            (far / c) / density,
        ];
        if [volume_pressure, near, far]
            .iter()
            .chain(&velocity)
            .any(|v| !v.is_finite() || *v == 0.0)
        {
            return Err(invalid(
                "molecular volume/force coefficients are unrepresentable",
            ));
        }
        Ok(Self {
            volume_pressure,
            near,
            far,
            velocity,
        })
    }

    fn inputs(
        &self,
        path: &Transfer,
        source: TransientAcousticSource,
    ) -> Result<[[f64; 4]; 4], AeroacError> {
        if source.compact == CompactSourceStrength::default()
            && source.force_impulse_n_s == [0.0; 3]
        {
            return observe::inputs(path, source.stress);
        }
        if !finite_compact(&source) {
            return Err(nonfinite("molecular volume/force history"));
        }
        let n = path.direction.map(Dd::from_f64);
        let force = source.compact.force_on_fluid_n.map(Dd::from_f64);
        let rate = source.compact.force_rate_on_fluid_n_s.map(Dd::from_f64);
        let impulse = source.force_impulse_n_s.map(Dd::from_f64);
        let radial =
            |vector: [Dd; 3]| (0..3).fold(Dd::ZERO, |sum, axis| sum + n[axis] * vector[axis]);
        let rf = radial(force);
        let rd = radial(rate);
        let ri = radial(impulse);
        let q = Dd::from_f64(source.compact.volume_flow_m3_s);
        let qdot = Dd::from_f64(source.compact.volume_flow_rate_m3_s2);
        let near = Dd::from_f64(self.near);
        let far = Dd::from_f64(self.far);
        let velocity = self.velocity.map(Dd::from_f64);
        let three = Dd::from_f64(3.0);
        let mut rows = [[Dd::ZERO; 4]; 4];
        rows[0][0] = qdot * Dd::from_f64(self.volume_pressure) + rf * near;
        rows[1][0] = rd * far;
        for axis in 0..3 {
            rows[0][axis + 1] =
                n[axis] * q * near + (three * n[axis] * ri - impulse[axis]) * velocity[0];
            rows[1][axis + 1] =
                n[axis] * qdot * far + (three * n[axis] * rf - force[axis]) * velocity[1];
            rows[2][axis + 1] = n[axis] * rd * velocity[2];
        }
        observe::rotate_velocity(path, &mut rows);
        if source.stress != TransientQuadrupoleSource::default() {
            let stress = observe::stress_inputs(path, source.stress)?;
            for order in 0..4 {
                for field in 0..4 {
                    rows[order][field] = rows[order][field] + stress[order][field];
                }
            }
        }
        observe::finish_inputs(rows)
    }
}

/// Complete source composition over the ORIGINAL stress propagation cache.
/// One kernel bank and one evaluation cover co-located volume, force and stress;
/// independent source positions retain their original ordering and delays.
pub struct PreparedMolecularSourceObserver {
    propagation: PreparedMolecularQuadrupoleObserver,
    compact: Vec<CompactTransfer>,
    workspace_bytes: usize,
}

impl PreparedMolecularSourceObserver {
    /// Prepare direct paths in the original material's frozen-speed medium.
    pub fn new(
        molecular: RelaxationState,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        config: MolecularTransientConfig,
        budget: MolecularTransientBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MaterialAeroacError> {
        poll(&mut cancelled)?;
        let frozen = FluidAcousticState::try_new(
            molecular.gas().density,
            molecular.frozen_sound_speed_m_s(),
        )
        .map_err(MaterialAeroacError::Material)?;
        let paths = PreparedTransientQuadrupoleObserver::new_cancellable(
            frozen,
            observer_m,
            source_positions_m,
            budget.maximum_paths,
            &mut cancelled,
        )?;
        Self::from_paths(molecular, paths, config, budget, cancelled)
    }

    /// Retain the complete original direct/ideal-wall/wedge path owner.
    /// Its density and speed must match this material's frozen-speed medium.
    /// All six added scales per path count against the SAME numerical byte cap.
    pub fn from_paths(
        molecular: RelaxationState,
        paths: PreparedTransientQuadrupoleObserver,
        config: MolecularTransientConfig,
        mut budget: MolecularTransientBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MaterialAeroacError> {
        poll(&mut cancelled)?;
        let projection_bytes = paths
            .path_count()
            .checked_mul(core::mem::size_of::<CompactTransfer>())
            .ok_or_else(|| invalid("molecular source projection extent overflow"))?;
        budget.maximum_workspace_bytes = budget
            .maximum_workspace_bytes
            .checked_sub(projection_bytes)
            .ok_or_else(|| invalid("molecular source projection byte budget"))?;
        let propagation = PreparedMolecularQuadrupoleObserver::from_paths(
            molecular,
            paths,
            config,
            budget,
            &mut cancelled,
        )?;
        let mut compact = Vec::new();
        compact
            .try_reserve_exact(propagation.geometry.path_count())
            .map_err(|_| invalid("molecular source projections exceed available memory"))?;
        for path in &propagation.geometry.transfers {
            poll(&mut cancelled)?;
            compact.push(CompactTransfer::new(path, propagation.geometry.acoustic())?);
        }
        let workspace_bytes = propagation.workspace_bytes + projection_bytes;
        poll(&mut cancelled)?;
        Ok(Self {
            propagation,
            compact,
            workspace_bytes,
        })
    }

    /// Observe the coherent complete source assembly on its original clock.
    ///
    /// Calling declares zero prehistory and unexcited exterior/material state.
    /// At the epoch, volume flow, force, force impulse, stress, stress rate and
    /// stress integral MUST be zero. Finite volume-flow rate, force rate and
    /// stress acceleration jumps are admitted; their onset is not assigned a
    /// smooth-start second-order accuracy claim. Missing delta impulses are
    /// never fabricated. Every derivative/integral must be consistent with
    /// its original source history; this adapter neither differentiates nor
    /// integrates arbitrary samples.
    ///
    /// Queries may be out of order across paths and always lie in the finite
    /// configured prefix. Each path's frozen-speed delay is applied once using
    /// the native DD clock. Later source/received tails are unspecified.
    pub fn observe_initially_resting(
        &self,
        clock: MolecularSampleClock,
        source_at: impl FnMut(usize, Dd) -> Result<TransientAcousticSource, AeroacError>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularSourceRecord, AeroacError> {
        self.propagation.observe_projected(
            clock,
            source_at,
            initially_resting,
            |index, path, source| self.compact[index].inputs(path, source),
            SCOPE,
            cancelled,
        )
    }

    /// The same prepared cache remains available to the original stress-only
    /// API, including its native recorded-history binding, without rebuilding.
    pub const fn stress_observer(&self) -> &PreparedMolecularQuadrupoleObserver {
        &self.propagation
    }
    pub const fn molecular_state(&self) -> &RelaxationState {
        self.propagation.molecular_state()
    }
    pub const fn paths(&self) -> &PreparedTransientQuadrupoleObserver {
        self.propagation.paths()
    }
    pub const fn config(&self) -> MolecularTransientConfig {
        self.propagation.config()
    }
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    pub const fn scope(&self) -> &'static str {
        SCOPE
    }
}

#[cfg(test)]
mod tests;
