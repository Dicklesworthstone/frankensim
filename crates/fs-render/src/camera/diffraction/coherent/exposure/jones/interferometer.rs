//! Full polarized spectral interference through the original finite shutter.
//!
//! The native Jones source owner supplies one two-component field per shared
//! frequency/covariance column. Uniform arm optics act before each frequency's
//! balanced combiner; its output fields then use the original temporal
//! coherency contraction. Independent source columns add only afterward.

use super::super::super::{
    CoherentDiffractionInterferometer, CoherentJonesArm,
    interferometer::visibility,
    polarization::{PolarizedDiffractionCoherency, interferometer::product},
};
use super::super::interferometer::{CoherentInterferometerExposureLimits, energy as scalar_energy};
use super::coherency::{energy, physical};
use super::{
    CoherentDiffractionError, Dd, JonesFields, MAX_COHERENT_DIFFRACTION_TONES,
    PreparedCoherentJonesDiffractionExposure, admit_count, limit, poll, rounded,
};
use fs_math::c64::C64;

/// Full polarized, shutter-averaged output densities of a passive two-arm receiver.
///
/// Each coherency retains both diagonal populations, its complex cross moment,
/// and native Stokes. Units are per tangent-coordinate area; a collected pixel
/// requires spatial integration. Visibility denotes an achromatic combiner-phase
/// scan, not the contrast of a broadband optical-path scan.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarizedCoherentInterferometerExposureSample {
    /// Arm powers after the passive split and uniform Jones arm elements.
    pub arm_power_density_w: [f64; 2],
    /// Full time/ensemble-averaged coherency at the plus and minus ports.
    pub output_coherency_density_w: [PolarizedDiffractionCoherency; 2],
    /// Trace of <E0(t) E1(t)†>, after the complete relative arm phase.
    pub mutual_coherence_density_w: C64,
    pub fringe_visibility: f64,
    pub arm_energy_density_j: [f64; 2],
    pub output_coherency_energy_density_j: [PolarizedDiffractionCoherency; 2],
    /// Original pupil transfers, including every retained patterned-pupil mode.
    pub pupil_evaluations: usize,
    /// Real frequency contractions, charged before any pupil transfer.
    pub pair_contractions: usize,
}

#[derive(Clone, Copy)]
struct InterferometerMoments {
    arms: [Dd; 2],
    outputs: [[Dd; 4]; 2],
    mutual: [Dd; 2],
}

impl InterferometerMoments {
    const ZERO: Self = Self {
        arms: [Dd::ZERO; 2],
        outputs: [[Dd::ZERO; 4]; 2],
        mutual: [Dd::ZERO; 2],
    };
}

impl PreparedCoherentJonesDiffractionExposure<'_, '_> {
    /// Recombine polarized spectral fields before integrating their finite shutter.
    ///
    /// The original pupil and its optional patterned Jones response remain
    /// frequency-specific. Each supplied `CoherentJonesArm` is an explicitly
    /// frequency-independent, uniform Jones sequence; additional dispersion is
    /// not inferred. The shared receiver's nondispersive optical path and sensor
    /// plane use the same authored-Hz phase owner as scalar interferometry.
    ///
    /// Covariance-column identities span all frequencies, polarizations and
    /// arms. Equal-frequency contributions combine before any outer product.
    /// Orthogonal arm fields may carry polarization interference even when
    /// their trace mutual coherence is zero. Passive loss is never normalized.
    ///
    /// With N distinct frequencies and B source banks, two active arms cost
    /// B*2*N*(6*N+1) real contractions; one active arm costs B*N*(3*N+1).
    /// Closed/empty receivers do no pupil or pair work. Explicit observation
    /// limits admit this query independently of preparation's original
    /// single-detector budget. Scratch is fixed and observation allocates nothing.
    /// Cancellation publishes no partial result and leaves preparation reusable.
    pub fn observe_interferometer(
        &self,
        receiver: CoherentDiffractionInterferometer,
        arms: [CoherentJonesArm; 2],
        limits: CoherentInterferometerExposureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PolarizedCoherentInterferometerExposureSample, CoherentDiffractionError> {
        poll(&mut cancelled)?;
        let (pupil_evaluations, pair_contractions) =
            self.jones_interferometer_work(receiver, limits)?;
        let moments =
            self.jones_interferometer_moments(receiver, &arms, [0.0; 2], &mut cancelled)?;
        let duration = Dd::from_f64(self.shutter().duration_s());
        let result = PolarizedCoherentInterferometerExposureSample {
            arm_power_density_w: [rounded(moments.arms[0])?, rounded(moments.arms[1])?],
            output_coherency_density_w: [
                physical(moments.outputs[0])?,
                physical(moments.outputs[1])?,
            ],
            mutual_coherence_density_w: C64::new(
                rounded(moments.mutual[0])?,
                rounded(moments.mutual[1])?,
            ),
            fringe_visibility: visibility(moments.arms, moments.mutual)?,
            arm_energy_density_j: [
                scalar_energy(moments.arms[0], duration)?,
                scalar_energy(moments.arms[1], duration)?,
            ],
            output_coherency_energy_density_j: [
                energy(moments.outputs[0], duration)?,
                energy(moments.outputs[1], duration)?,
            ],
            pupil_evaluations,
            pair_contractions,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    #[allow(clippy::float_cmp)] // The prepared sensor plane is an authored value.
    fn jones_interferometer_work(
        &self,
        receiver: CoherentDiffractionInterferometer,
        limits: CoherentInterferometerExposureLimits,
    ) -> Result<(usize, usize), CoherentDiffractionError> {
        if receiver.sensor_distance_m() != self.optics.sensor_distance_m {
            return Err(CoherentDiffractionError::InvalidSensorPlane);
        }
        let evaluations =
            receiver.pupil_evaluations(self.pupil_evaluations, limits.maximum_pupil_evaluations)?;
        let active = receiver
            .arm_amplitudes()
            .iter()
            .filter(|&&a| a != Dd::ZERO)
            .count();
        let quantity = "Jones interferometer exposure pair contractions";
        let contractions = if evaluations == 0 {
            0
        } else {
            let n = self.frequency_count();
            // Arm traces: 2*N*(N+1). Output coherencies: 2*N*(3*N+1).
            // Trace mutual coherence: 2*N*(2*N-1). A single active arm
            // shares one coherency contraction between its two half-power ports.
            let per_bank = if active == 1 {
                n * (3 * n + 1)
            } else {
                2 * n * (6 * n + 1)
            };
            self.factors
                .checked_add(1)
                .and_then(|banks| banks.checked_mul(per_bank))
                .ok_or_else(|| limit(quantity, limits.maximum_pair_contractions))?
        };
        admit_count(quantity, contractions, limits.maximum_pair_contractions)?;
        Ok((evaluations, contractions))
    }

    fn jones_interferometer_moments(
        &self,
        receiver: CoherentDiffractionInterferometer,
        arms: &[CoherentJonesArm; 2],
        relay_offset: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<InterferometerMoments, CoherentDiffractionError> {
        let mut result = InterferometerMoments::ZERO;
        let amplitudes = receiver.arm_amplitudes();
        let active = amplitudes.iter().filter(|&&a| a != Dd::ZERO).count();
        if active == 0 || self.pupil_evaluations == 0 {
            return Ok(result);
        }
        let mut coordinates = receiver.coordinates();
        for at in &mut coordinates {
            for axis in 0..2 {
                at[axis] = rounded(Dd::from_f64(at[axis]) + Dd::from_f64(relay_offset[axis]))?;
            }
        }
        let rotations =
            self.optics
                .interferometer_rotations_offset(receiver, relay_offset, cancelled)?;
        let identity = [[Dd::ONE, Dd::ZERO]; MAX_COHERENT_DIFFRACTION_TONES];
        let combiner = [Dd::from_f64(0.5).sqrt(), Dd::ZERO];
        for bank in 0..=self.factors {
            poll(cancelled)?;
            let mut fields: [JonesFields; 2] =
                [[[[Dd::ZERO; 2]; 2]; MAX_COHERENT_DIFFRACTION_TONES]; 2];
            for arm in 0..2 {
                if amplitudes[arm] == Dd::ZERO {
                    continue;
                }
                fields[arm] = self.group_fields(bank, coordinates[arm], cancelled)?;
                for i in 0..self.frequency_count() {
                    poll(cancelled)?;
                    let split = [amplitudes[arm], Dd::ZERO];
                    let field = arms[arm].apply([
                        product(fields[arm][i][0], split)?,
                        product(fields[arm][i][1], split)?,
                    ])?;
                    fields[arm][i] = [
                        product(field[0], rotations[arm][i])?,
                        product(field[1], rotations[arm][i])?,
                    ];
                }
            }
            if active == 1 {
                let arm = usize::from(amplitudes[0] == Dd::ZERO);
                let moment = self.contract_coherency(&mut fields[arm], &identity, cancelled)?;
                result.arms[arm] = result.arms[arm] + moment[0] + moment[1];
                let port = moment.map(|value| value * Dd::from_f64(0.5));
                for axis in 0..2 {
                    if moment[axis] != Dd::ZERO && port[axis] == Dd::ZERO {
                        return Err(CoherentDiffractionError::Unrepresentable);
                    }
                }
                for output in &mut result.outputs {
                    for i in 0..4 {
                        output[i] = output[i] + port[i];
                    }
                }
                continue;
            }
            let mut bank_arms = [Dd::ZERO; 2];
            for arm in 0..2 {
                for axis in 0..2 {
                    poll(cancelled)?;
                    let mut component = fields[arm].map(|field| field[axis]);
                    bank_arms[arm] = bank_arms[arm]
                        + self
                            .optics
                            .contract_component(&mut component, &identity, cancelled)?;
                }
                result.arms[arm] = result.arms[arm] + bank_arms[arm];
            }
            for (port, sign) in [Dd::ONE, -Dd::ONE].into_iter().enumerate() {
                let mut output = [[[Dd::ZERO; 2]; 2]; MAX_COHERENT_DIFFRACTION_TONES];
                for i in 0..self.frequency_count() {
                    poll(cancelled)?;
                    for axis in 0..2 {
                        output[i][axis] = product(
                            std::array::from_fn(|part| {
                                fields[0][i][axis][part] + sign * fields[1][i][axis][part]
                            }),
                            combiner,
                        )?;
                    }
                }
                let moment = self.contract_coherency(&mut output, &identity, cancelled)?;
                for i in 0..4 {
                    result.outputs[port][i] = result.outputs[port][i] + moment[i];
                }
            }
            let mut mutual = [Dd::ZERO; 2];
            for axis in 0..2 {
                let cross = self.optics.interferometer_cross(
                    &fields[0].map(|field| field[axis]),
                    &fields[1].map(|field| field[axis]),
                    cancelled,
                )?;
                for part in 0..2 {
                    mutual[part] = mutual[part] + cross[part];
                }
            }
            // A failed covariance bank cannot be hidden by unrelated sources.
            visibility(bank_arms, mutual)?;
            for part in 0..2 {
                result.mutual[part] = result.mutual[part] + mutual[part];
            }
        }
        Ok(result)
    }
}
