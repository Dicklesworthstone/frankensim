//! Complete moving volume, force and stress fields on one original causal event.
//!
//! This composes the native compact-source and moving-stress Green derivatives.
//! It neither differentiates a new source record nor replaces either physical
//! past-worldline integral with a scalar impulse or delayed stationary field.

use super::{MovingSourceArrayField, PreparedMovingSourceArray, ReflectedFieldSum, Source};
use crate::AeroacError;
use crate::curle3d::moving::trajectory::{invalid, poll};
use crate::curle3d::moving::{CompactSourceStrength, MovingCompactObserver};
use crate::curle3d::quadrupole::transient::TransientQuadrupoleSource;
use fs_math::dd::Dd;

/// Original physical source mechanisms carried by one accepted moving path.
///
/// All values and derivatives belong to the SAME returned emission time and
/// fixed fluid-frame axes. Explicit zero fields omit a mechanism. A producer
/// must keep its volume/traction/stress sources disjoint; a permeable-surface
/// description must not count a volume source a second time.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MovingMultipoleSource {
    /// Volume injection [m³/s], force ON fluid [N], and their first derivatives.
    pub compact: CompactSourceStrength,
    /// Physical symmetric stress moment [N m] and its first two derivatives.
    /// Its scalar tensor integral is UNUSED here: motion needs the separate
    /// spatially weighted past-worldline integral in `MovingMultipoleHistory`.
    pub stress: TransientQuadrupoleSource,
}

/// Two distinct past-source integrals at one current virtual microphone.
///
/// Integrate on the original source clock through the returned emission time,
/// using this EXACT receiver position for the entire past source worldline.
/// `MovingCompactObserver::{force_history_integrand,quadrupole_history_integrand}`
/// supply the original kernels; neither includes image parity or rotation.
/// Continuous history clearance, integration accuracy and prior forcing remain
/// caller-owned. Zero explicitly declares no prior contribution; it does not
/// reconstruct missing history or an unknown homogeneous initial acoustic field.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MovingMultipoleHistory {
    /// Integral of [3 n(n.F)-F]/R³ [N s/m³].
    pub force_n_s_m3: [Dd; 3],
    /// Integral of [15 n(n.Q.n)-6 Q n-3 tr(Q)n]/R⁴ [N s/m³].
    pub stress_n_s_m3: [Dd; 3],
}

impl PreparedMovingSourceArray<'_> {
    /// Coherent complete Q/F/stress fields at one original microphone event.
    ///
    /// Each original source/image solves its causal root ONCE. The same root
    /// and analytic position/velocity/acceleration/jerk feed the native compact
    /// and stress owners; both callbacks run once on that emission clock.
    /// The history callback receives source and image indices plus the actual
    /// virtual observer at this receiver event. Direct sound has image zero.
    /// Supply untransformed physical histories: the native image owner applies
    /// pressure parity and the full polar-vector map once to the complete field.
    ///
    /// All source mechanisms and paths accumulate in DD before output rounding
    /// and signed instantaneous p*u. Microphone velocity is not added to acoustic
    /// velocity or used as a pressure multiplier. Uniform-flow preparations keep
    /// their original fluid-frame axes; p*u excludes mean-flow energy advection.
    /// Each source may move independently and span different emission times.
    /// No relative phase, initial field, missing tail or source independence is
    /// inferred. The native source-model/compactness and lossless-fluid limits
    /// remain; this does not provide finite-body scattering or source feedback.
    ///
    /// Queries allocate nothing. Original aggregate source/path/root limits
    /// bound the complete operation; the added local field contractions have
    /// fixed work per path. Callback integration work retains its own limits.
    /// Cancellation, missing/nonfinite histories or any native refusal returns
    /// no partial field. Source/path ordering is the constructor's exact order.
    pub fn multipole_field_at(
        &self,
        receiver_time_s: f64,
        source_at: &mut impl FnMut(usize, f64) -> Result<MovingMultipoleSource, AeroacError>,
        history_at: &mut impl FnMut(
            usize,
            usize,
            &MovingCompactObserver,
            f64,
        ) -> Result<MovingMultipoleHistory, AeroacError>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MovingSourceArrayField, AeroacError> {
        poll(cancelled)?;
        let [first, last] = self.receiver_time_interval_s;
        if !receiver_time_s.is_finite() || receiver_time_s < first || receiver_time_s > last {
            return Err(invalid(
                "receiver time lacks complete moving multipole source history",
            ));
        }
        let mut sum = ReflectedFieldSum::default();
        for index in 0..self.source_count() {
            let source = self.sources.at(index);
            for path in 0..source.path_count() {
                poll(cancelled)?;
                let pair = match source {
                    Source::Direct(pair) => pair,
                    Source::Reflected(reflected) => reflected.path(path)?,
                };
                let original = pair.prepare_receiver(receiver_time_s, cancelled)?.pressure;
                let time = original.source_time_s();
                let residual = original.retarded_time_residual_s();
                let observer = pair.observer_at(receiver_time_s)?;
                let (motion, jerk) = pair.source().kinematics_with_jerk_at(time)?;
                let compact = observer.prepare_force_sample(original, motion)?;
                let stress = observer.prepare_quadrupole_field_sample(original, motion, jerk)?;
                poll(cancelled)?;
                let strength = source_at(index, time)?;
                poll(cancelled)?;
                let history = history_at(index, path, &observer, time)?;
                let (_, p_compact, v_compact) =
                    compact.evaluate_dd(strength.compact, history.force_n_s_m3)?;
                let (p_stress, v_stress) =
                    stress.evaluate_dd(strength.stress, history.stress_n_s_m3)?;
                let pressure = p_compact + p_stress;
                let velocity = core::array::from_fn(|i| v_compact[i] + v_stress[i]);
                match source {
                    Source::Direct(_) => sum.add(pressure, velocity, residual),
                    Source::Reflected(reflected) => {
                        reflected.add_field_path(path, pressure, velocity, residual, &mut sum)?;
                    }
                }
            }
        }
        poll(cancelled)?;
        let field = sum.finish(receiver_time_s)?;
        Ok(MovingSourceArrayField {
            receiver_time_s,
            source_count: self.source_count(),
            source_path_count: self.source_path_count,
            pressure_pa: field.pressure_pa,
            particle_velocity_m_s: field.particle_velocity_m_s,
            instantaneous_intensity_w_m2: field.instantaneous_intensity_w_m2,
            maximum_time_residual_s: field.maximum_time_residual_s,
        })
    }
}
