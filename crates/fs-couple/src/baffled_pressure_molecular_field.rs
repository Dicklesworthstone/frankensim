//! Complete causal fields of the original Rayleigh surface in a relaxing gas.
//!
//! A surface sample has normal velocity v_n = sum(phi_m * v_m). Its source
//! for the existing full-space molecular Green operator is Q = 2 A v_n and
//! Qdot = 2 A sum(phi_m * a_m). The factor two is the rigid-baffle image, not
//! an extra source or an intensity multiplier. Every quadrature point retains
//! its own position, propagation direction and frozen-speed arrival front.
//!
//! The original molecular owner evaluates p = rho s Q G and u = -Q grad(G),
//! with G = exp(-kappa(s) R)/(4 pi R). Thus absorption and dispersion enter
//! the velocity gradient as well as pressure; the reactive near field is not
//! reconstructed from p/(rho c). Signed instantaneous intensity follows the
//! coherent field sum. Rayleigh's half-space source normalization is given in
//! <https://euphonics.org/4-3-2-the-rayleigh-integral-and-the-baffled-piston/>.
//!
//! This is an offline finite-prefix observation of prescribed, initially
//! resting motion in a stationary z=0 infinite rigid baffle. It composes the
//! existing BDF2 convolution quadrature and the ORIGINAL caloric-reservoir
//! Laplace response. No classical viscothermal loss, thermal/shear boundary
//! layers, finite-baffle scattering, source motion or backreaction is added.
//! Unlike the effective harmonic/gas pressure path, this does not add classical
//! attenuation to the reservoir law. Source, spatial, timestep and contour
//! refinement remain separate; no automatic bandwidth or measured-SPL claim.

use super::{MAX_RAYLEIGH_INPUTS, RayleighMedium, SurfaceSample};
use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::moving::CompactSourceStrength;
use fs_aeroac::curle3d::quadrupole::transient::molecular::sources::{
    MolecularSourceRecord, PreparedMolecularSourceObserver, TransientAcousticSource,
};
use fs_aeroac::material::MaterialAeroacError;
use fs_bem::{helmholtz::Medium, near_field::FirstOrder};
use fs_material::relaxation::RelaxationState;
use fs_math::{c64::C64, dd::Dd};

pub use fs_aeroac::curle3d::quadrupole::transient::molecular::sources::{
    MolecularSampleClock, MolecularTransientBudget, MolecularTransientConfig, RetardedSourceField,
};

const SCOPE: &str = "Initially resting prescribed stationary flat Rayleigh surface in one \
homogeneous molecularly relaxing gas. Original signed area/basis coefficients, individual \
frozen-speed fronts, complete coherent pressure/vector velocity and signed instantaneous \
intensity. Finite-prefix BDF2 convolution quadrature; no classical viscothermal loss, \
finite-baffle scattering, moving geometry, backreaction or measured-SPL qualification.";

/// Original coordinate motion at ONE physical emission time.
///
/// These are the first and second time derivatives of the same input
/// coordinate, in the units selected by `SurfaceSample::mode_shape`.
/// `mode_shape * velocity` must be normal velocity [m/s]; for a shape of
/// 1/area the input is volume flow [m³/s] and acceleration is its rate [m³/s²].
/// The observer never differences an arbitrary velocity record or invents
/// an integration constant. Supply consistent derivatives from its owner.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct InputMotion {
    pub velocity: f64,
    pub acceleration: f64,
}

/// Surface admission and the additional work of projecting its input basis.
#[derive(Clone, Copy, Debug)]
pub struct Options {
    /// Positive receiver height above z=0 [m], as in the harmonic observer.
    /// Smaller clearances require finer surface quadrature.
    pub minimum_clearance_m: f64,
    /// Bound on DD products for coefficient preparation plus ONE observation:
    /// 2 * samples_on_surface * inputs * (received_samples + 2).
    /// This includes the native per-source onset queries. Callback work is
    /// caller-owned; the separate native transform budget remains unchanged.
    pub maximum_projection_terms: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            minimum_clearance_m: 0.05,
            maximum_projection_terms: 1 << 26,
        }
    }
}

/// One received finite prefix, bound to its ORIGINAL gas microphone reference.
/// Fields retain the native propagation owner's scope alongside the narrower
/// surface scope available through `scope()`. Later received tails are unknown.
#[derive(Clone, Debug)]
pub struct MolecularBaffledRecord {
    original: MolecularSourceRecord,
    medium: RayleighMedium,
}

impl MolecularBaffledRecord {
    pub const fn clock(&self) -> MolecularSampleClock {
        self.original.clock()
    }

    pub fn fields(&self) -> &[RetardedSourceField] {
        self.original.fields()
    }

    /// Original integer clock in DD; no accumulated timestep or new epoch.
    pub fn time_at_clock_dd(&self, index: usize) -> Option<Dd> {
        self.original.time_at_clock_dd(index)
    }

    /// Equilibrium rho and c from the same immutable gas that propagated sound.
    /// The wavefront instead uses that gas's calorically frozen sound speed.
    pub const fn microphone_medium(&self) -> RayleighMedium {
        self.medium
    }

    /// Ideal pressure/velocity projection at one received sample. Pattern and
    /// orientation may change without rerunning or retiming the source history.
    /// Its real equilibrium rho*c sensitivity is a microphone calibration, not
    /// a replacement for the dispersive field impedance or a plane-wave claim.
    pub fn microphone_pressure_pa(
        &self,
        index: usize,
        microphone: FirstOrder,
    ) -> Result<f64, String> {
        let field = self
            .fields()
            .get(index)
            .ok_or("molecular Rayleigh sample index")?;
        microphone
            .observe(
                C64::from_re(field.pressure_pa),
                field.particle_velocity_m_s.map(C64::from_re),
                Medium {
                    density: self.medium.density,
                    sound_speed: self.medium.sound_speed,
                },
            )
            .map(|value| value.re)
            .map_err(|error| error.to_string())
    }

    pub const fn scope(&self) -> &'static str {
        SCOPE
    }
}

/// One prepared material Green cache and the ORIGINAL signed surface basis.
/// Each source point is a separate propagation path, even when its delay equals
/// another point's. Modal products remain DD until the complete local volume
/// source is formed. The native observer then sums coherent Cartesian fields.
pub struct PreparedMolecularBaffledField {
    propagation: PreparedMolecularSourceObserver,
    volume_coefficients: Vec<Dd>,
    inputs: usize,
    workspace_bytes: usize,
    projection_terms: usize,
}

fn invalid(what: &'static str) -> AeroacError {
    AeroacError::InvalidParameter { what }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AeroacError> {
    if cancelled() {
        Err(AeroacError::Cancelled {
            what: "molecular Rayleigh field",
        })
    } else {
        Ok(())
    }
}

fn finite(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite() && value.to_f64().is_finite()
}

impl PreparedMolecularBaffledField {
    /// Prepare the original positive-area z=0 quadrature and complete input
    /// basis for one receiver in z>0. No independent density, sound speed or
    /// path attenuation is accepted: all derive from `molecular` exactly once.
    ///
    /// The native byte cap also includes retained DD area/basis coefficients
    /// and temporary source positions. The original surface, caller histories,
    /// material metadata and allocator overhead remain caller-owned. The native
    /// path, sample, transform and work limits apply before CQ preparation;
    /// `Options` additionally bounds all surface/input projection work.
    /// Cancellation, invalid geometry, unrepresentable coefficients or any
    /// exhausted budget returns no partially prepared observer.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        surface: &[SurfaceSample],
        inputs: usize,
        position_m: [f64; 3],
        molecular: RelaxationState,
        config: MolecularTransientConfig,
        mut budget: MolecularTransientBudget,
        options: Options,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MaterialAeroacError> {
        poll(&mut cancelled)?;
        if !(1..=MAX_RAYLEIGH_INPUTS).contains(&inputs)
            || surface.is_empty()
            || surface.len() > 120_000
            || surface.len() > budget.maximum_paths
            || !options.minimum_clearance_m.is_finite()
            || options.minimum_clearance_m <= 0.0
            || position_m.iter().any(|v| !v.is_finite())
            || position_m[2] < options.minimum_clearance_m
        {
            return Err(
                invalid("molecular Rayleigh surface, input, path count or clearance").into(),
            );
        }
        let coefficients = surface
            .len()
            .checked_mul(inputs)
            .ok_or_else(|| invalid("molecular Rayleigh basis extent overflow"))?;
        let projection_terms = config
            .samples
            .checked_add(2)
            .and_then(|n| n.checked_mul(coefficients))
            .and_then(|n| n.checked_mul(2))
            .filter(|&n| n <= options.maximum_projection_terms)
            .ok_or_else(|| invalid("molecular Rayleigh complete projection work budget"))?;
        let mapping_bytes = coefficients
            .checked_mul(core::mem::size_of::<Dd>())
            .and_then(|n| {
                n.checked_add(
                    surface
                        .len()
                        .checked_mul(core::mem::size_of::<[f64; 3]>())?,
                )
            })
            .ok_or_else(|| invalid("molecular Rayleigh mapping extent overflow"))?;
        budget.maximum_workspace_bytes = budget
            .maximum_workspace_bytes
            .checked_sub(mapping_bytes)
            .ok_or_else(|| invalid("molecular Rayleigh mapping byte budget"))?;
        let gas = molecular.gas();
        let impedance = gas.density * gas.sound_speed;
        if !impedance.is_finite() || impedance <= 0.0 {
            return Err(invalid("molecular Rayleigh equilibrium microphone impedance").into());
        }
        let mut volume_coefficients = Vec::new();
        let mut positions = Vec::new();
        volume_coefficients
            .try_reserve_exact(coefficients)
            .map_err(|_| invalid("molecular Rayleigh basis allocation"))?;
        positions
            .try_reserve_exact(surface.len())
            .map_err(|_| invalid("molecular Rayleigh position allocation"))?;
        for (index, sample) in surface.iter().enumerate() {
            if index.is_multiple_of(64) {
                poll(&mut cancelled)?;
            }
            if !sample.area_m2.is_finite()
                || sample.area_m2 <= 0.0
                || sample.position_m.iter().any(|v| !v.is_finite())
                || sample.position_m[2] != 0.0
                || sample.mode_shape.len() != inputs
                || sample.mode_shape.iter().any(|v| !v.is_finite())
            {
                return Err(invalid(
                    "molecular Rayleigh requires finite flat areas and complete shapes",
                )
                .into());
            }
            positions.push(sample.position_m);
            for &shape in &sample.mode_shape {
                let coefficient =
                    Dd::from_f64(sample.area_m2) * Dd::from_f64(shape) * Dd::from_f64(2.0);
                if !finite(coefficient) {
                    return Err(invalid("molecular Rayleigh volume coefficient overflow").into());
                }
                volume_coefficients.push(coefficient);
            }
        }
        let propagation = PreparedMolecularSourceObserver::new(
            molecular,
            position_m,
            &positions,
            config,
            budget,
            &mut cancelled,
        )?;
        let workspace_bytes = propagation
            .workspace_bytes()
            .checked_add(mapping_bytes)
            .ok_or_else(|| invalid("molecular Rayleigh total workspace overflow"))?;
        poll(&mut cancelled)?;
        Ok(Self {
            propagation,
            volume_coefficients,
            inputs,
            workspace_bytes,
            projection_terms,
        })
    }

    /// Observe ONE consistent original motion history on its physical clock.
    /// The callback takes the original input index and DD emission time, which
    /// can be out of order between surface paths. All queries lie in the native
    /// finite source prefix; no received or source tail is extrapolated.
    ///
    /// Calling declares zero velocity before the epoch and an unexcited gas.
    /// Every input velocity at the epoch MUST be zero, including coordinates
    /// whose surface contributions cancel. A finite acceleration jump is
    /// allowed without a smooth-onset second-order accuracy claim. Nonfinite
    /// values and missing physical onset impulses are refused, not repaired.
    ///
    /// Velocity and acceleration must be deterministic, consistent derivatives
    /// of the same coordinate; no differentiation, integration or antialiasing
    /// is performed here. Refine source resolution and CQ timestep separately.
    /// Cancellation returns no partial record and leaves this cache reusable.
    pub fn observe_initially_resting(
        &self,
        clock: MolecularSampleClock,
        mut input_at: impl FnMut(usize, Dd) -> Result<InputMotion, AeroacError>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularBaffledRecord, AeroacError> {
        let epoch = Dd::from_f64(clock.epoch_s());
        let original = self.propagation.observe_initially_resting(
            clock,
            |source, time| {
                let row =
                    &self.volume_coefficients[source * self.inputs..(source + 1) * self.inputs];
                let mut volume = Dd::ZERO;
                let mut rate = Dd::ZERO;
                for (input, &coefficient) in row.iter().enumerate() {
                    let motion = input_at(input, time)?;
                    if !motion.velocity.is_finite() || !motion.acceleration.is_finite() {
                        return Err(AeroacError::NonFinite {
                            what: "molecular Rayleigh input motion",
                        });
                    }
                    if time == epoch && motion.velocity != 0.0 {
                        return Err(invalid(
                            "molecular Rayleigh requires every input initially at rest",
                        ));
                    }
                    volume = volume + coefficient * Dd::from_f64(motion.velocity);
                    rate = rate + coefficient * Dd::from_f64(motion.acceleration);
                }
                if !finite(volume) || !finite(rate) {
                    return Err(AeroacError::NonFinite {
                        what: "molecular Rayleigh projected volume source",
                    });
                }
                Ok(TransientAcousticSource {
                    compact: CompactSourceStrength {
                        volume_flow_m3_s: volume.to_f64(),
                        volume_flow_rate_m3_s2: rate.to_f64(),
                        ..CompactSourceStrength::default()
                    },
                    ..TransientAcousticSource::default()
                })
            },
            cancelled,
        )?;
        Ok(MolecularBaffledRecord {
            original,
            medium: RayleighMedium::from(self.molecular_state().gas()),
        })
    }

    pub const fn input_count(&self) -> usize {
        self.inputs
    }
    pub fn surface_sample_count(&self) -> usize {
        self.propagation.paths().source_count()
    }
    pub fn source_position_m(&self, index: usize) -> Option<[f64; 3]> {
        self.propagation.paths().source_position_m(index)
    }
    pub fn source_delay_s(&self, index: usize) -> Option<f64> {
        self.propagation.paths().source_delay_s(index)
    }
    pub const fn molecular_state(&self) -> &RelaxationState {
        self.propagation.molecular_state()
    }
    pub const fn config(&self) -> MolecularTransientConfig {
        self.propagation.config()
    }
    /// Peak numerical payload bound, including the surface mapping.
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    /// Admitted DD projection product bound for preparation plus one record.
    pub const fn projection_terms(&self) -> usize {
        self.projection_terms
    }
    pub const fn scope(&self) -> &'static str {
        SCOPE
    }
}

#[cfg(test)]
#[path = "baffled_pressure_molecular_field_tests.rs"]
mod tests;
