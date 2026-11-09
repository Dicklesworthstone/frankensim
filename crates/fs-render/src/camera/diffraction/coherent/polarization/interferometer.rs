//! Polarization-resolved two-arm detection on the original camera pupil.
//!
//! Every source bank retains both transverse fields through independent passive
//! arm optics. The balanced output fields form full positive coherencies before
//! independent covariance banks add. Orthogonal arm polarizations have zero
//! trace interference but can still carry phase-dependent output polarization.

use super::super::{
    CoherentDiffractionError, CoherentDiffractionInterferometer, CoherentJonesPupilLayer,
    admit_count,
    exposure::admit_sensor,
    interferometer::{InterferometerTemporalCoherence, TemporalWeights, mul, visibility},
    poll, rounded,
};
use super::{
    JonesField, Moments, PolarizedDiffractionCoherency, PreparedPolarizedCoherentDiffraction,
    coherency, is_zero, outer,
};
use fs_math::{c64::C64, dd::Dd};

type Matrix = [[[Dd; 2]; 2]; 2];
const ZERO_MATRIX: Matrix = [[[Dd::ZERO; 2]; 2]; 2];

/// Maximum uniform, passive Jones elements retained in one interferometer arm.
pub const MAX_COHERENT_JONES_ARM_LAYERS: usize = 8;

/// Ordered uniform optical elements acting on one arm's two Jones components.
///
/// Every element uses the original admitted `CoherentJonesPupilLayer`: principal
/// power transmissions lie in [0,1], with authored phases and transverse axes.
/// Uniform polarizers, retarders and their ordered compositions are supported.
/// The retained layer sequence preserves aligned weak principal populations.
/// No relative accuracy is certified for arbitrary rotations or interference nulls.
/// These optics act after the pupil-to-sensor propagation; spatial
/// pupil gradients and a change of ray geometry are outside this arm model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentJonesArm {
    layers: [Matrix; MAX_COHERENT_JONES_ARM_LAYERS],
    count: usize,
}

impl CoherentJonesArm {
    /// A clear arm with no polarization transformation or additional loss.
    pub const IDENTITY: Self = Self {
        layers: [ZERO_MATRIX; MAX_COHERENT_JONES_ARM_LAYERS],
        count: 0,
    };

    /// Cache up to eight native uniform Jones layers in optical traversal order.
    ///
    /// Principal passivity is enforced by the original layer constructor.
    /// Nonzero spatial gradients are refused: an image arm has one uniform
    /// matrix per element, rather than the pupil's spatial Fourier expansion.
    pub fn try_new(layers: &[CoherentJonesPupilLayer]) -> Result<Self, CoherentDiffractionError> {
        admit_count(
            "Jones arm layers",
            layers.len(),
            MAX_COHERENT_JONES_ARM_LAYERS,
        )?;
        let mut result = Self::IDENTITY;
        for (index, layer) in layers.iter().enumerate() {
            result.layers[index] = layer.uniform_matrix()?;
        }
        result.count = layers.len();
        Ok(result)
    }

    pub(in super::super) fn apply(
        &self,
        mut field: JonesField,
    ) -> Result<JonesField, CoherentDiffractionError> {
        for matrix in &self.layers[..self.count] {
            let mut transformed = [[Dd::ZERO; 2]; 2];
            for row in 0..2 {
                for column in 0..2 {
                    let term = product(matrix[row][column], field[column])?;
                    for axis in 0..2 {
                        transformed[row][axis] = transformed[row][axis] + term[axis];
                    }
                }
                for value in transformed[row] {
                    rounded(value)?;
                }
            }
            field = transformed;
        }
        Ok(field)
    }
}

/// Ensemble point detection per tangent-coordinate area in the common frame.
///
/// Both ports retain separate diagonal populations, complex cross moments and
/// native Stokes. No rounded Stokes vector is converted back into Jones fields.
/// Equal-area unit-magnification relays have the shared scalar receiver's scope.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarizedCoherentInterferometerSample {
    /// Arm powers after the passive split and all authored arm optics.
    pub arm_power_density_w: [f64; 2],
    /// Total physical coherencies at the plus and minus combiner outputs.
    pub output_coherency_density_w: [PolarizedDiffractionCoherency; 2],
    /// Trace mutual coherence sum_p <E0_p conj(E1_p)> after all arm phases.
    pub mutual_coherence_density_w: C64,
    /// Total-power fringe visibility 2*|trace(W01)|/(P0+P1), or zero when dark.
    /// Polarization orthogonality, unequal powers and partial coherence reduce it.
    pub fringe_visibility: f64,
    /// Original pupil transfers actually evaluated, shared by both components.
    pub pupil_evaluations: usize,
}

impl PreparedPolarizedCoherentDiffraction<'_, '_> {
    /// Recombine original Jones fields after independent passive arm optics.
    ///
    /// Source banks remain independent; matching wavelengths do not establish
    /// coherence. Optical path difference and physical sensor chirp use the
    /// shared scalar receiver's phase owner and original medium wavelength.
    /// Output fields retain DD differences before their outer products, and
    /// the combiner's amplitude normalization precedes squaring.
    ///
    /// Work is preflighted, allocation-free and bounded by the original wave
    /// count plus at most eight uniform layers per active arm and source bank.
    /// Cancellation publishes no partial observation and leaves the query reusable.
    pub fn observe_interferometer(
        &self,
        receiver: CoherentDiffractionInterferometer,
        arms: [CoherentJonesArm; 2],
        maximum_pupil_evaluations: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PolarizedCoherentInterferometerSample, CoherentDiffractionError> {
        self.observe_interferometer_with_coherence(
            receiver, arms, None, maximum_pupil_evaluations, cancelled,
        )
    }

    /// Stationary finite-linewidth interference through the original Jones arms.
    ///
    /// The temporal correlation attenuates the complete cross-arm coherency,
    /// including phase-dependent output polarization when the arm polarizations
    /// are orthogonal. It leaves each arm's polarization and power unchanged.
    /// Positive coherent/independent-arm outer products preserve weak populations
    /// and the original source-bank independence. Jones optics and pupil transfer
    /// remain frozen at the central wavelength under the explicit narrowband law.
    pub fn observe_interferometer_with_linewidth(
        &self,
        receiver: CoherentDiffractionInterferometer,
        arms: [CoherentJonesArm; 2],
        coherence: InterferometerTemporalCoherence,
        maximum_pupil_evaluations: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PolarizedCoherentInterferometerSample, CoherentDiffractionError> {
        self.observe_interferometer_with_coherence(
            receiver, arms, Some(coherence), maximum_pupil_evaluations, cancelled,
        )
    }

    fn observe_interferometer_with_coherence(
        &self,
        receiver: CoherentDiffractionInterferometer,
        arms: [CoherentJonesArm; 2],
        coherence: Option<InterferometerTemporalCoherence>,
        maximum_pupil_evaluations: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PolarizedCoherentInterferometerSample, CoherentDiffractionError> {
        poll(&mut cancelled)?;
        let temporal = if let Some(coherence) = coherence {
            coherence.admit(self.wavelength_nm)?;
            coherence.weights([0.0; 2])?
        } else {
            TemporalWeights::COHERENT
        };
        admit_sensor(self.camera, receiver.sensor_distance_m())?;
        let pupil_evaluations =
            receiver.pupil_evaluations(self.waves, maximum_pupil_evaluations)?;
        let coordinates = receiver.coordinates();
        let amplitudes = receiver.arm_amplitudes();
        let rotation = receiver.relative_rotation(
            self.wavelength_nm,
            self.camera.pupil_wavelength_nm(self.wavelength_nm)?,
            coordinates,
        )?;
        let combiner = [Dd::from_f64(0.5).sqrt(), Dd::ZERO];
        let mut power = [Dd::ZERO; 2];
        let mut output: [Moments; 2] = [[Dd::ZERO; 4]; 2];
        let mut mutual = [Dd::ZERO; 2];
        for bank in core::iter::once(self.illumination.mean)
            .chain(self.illumination.factors.iter().map(|factor| factor.waves))
        {
            poll(&mut cancelled)?;
            let mut fields = [[[Dd::ZERO; 2]; 2]; 2];
            for arm in 0..2 {
                if !is_zero(amplitudes[arm]) {
                    let field = self.field(bank, coordinates[arm], &mut cancelled)?;
                    let split = [amplitudes[arm], Dd::ZERO];
                    fields[arm] =
                        arms[arm].apply([product(field[0], split)?, product(field[1], split)?])?;
                }
            }
            for component in &mut fields[1] {
                *component = product(*component, rotation)?;
            }
            let mut independent = [Dd::ZERO; 4];
            for arm in 0..2 {
                let moment = outer(fields[arm])?;
                power[arm] = power[arm] + (moment[0] + moment[1]);
                for index in 0..4 {
                    independent[index] = independent[index]
                        + moment[index] * (temporal.incoherent * Dd::from_f64(0.5));
                }
            }
            for component in 0..2 {
                let cross = product(
                    fields[0][component],
                    [fields[1][component][0], -fields[1][component][1]],
                )?;
                for axis in 0..2 {
                    mutual[axis] = mutual[axis] + cross[axis] * temporal.coherent;
                }
            }
            for (port, sign) in [Dd::ONE, -Dd::ONE].into_iter().enumerate() {
                let mut field = [[Dd::ZERO; 2]; 2];
                for component in 0..2 {
                    field[component] = product(
                        std::array::from_fn(|axis| {
                            (fields[0][component][axis] + sign * fields[1][component][axis])
                                * temporal.field_scale
                        }),
                        combiner,
                    )?;
                }
                let moment = outer(field)?;
                for index in 0..4 {
                    output[port][index] = output[port][index]
                        + (moment[index] + independent[index]);
                }
            }
        }
        let result = PolarizedCoherentInterferometerSample {
            arm_power_density_w: [rounded(power[0])?, rounded(power[1])?],
            output_coherency_density_w: [coherency(output[0])?, coherency(output[1])?],
            mutual_coherence_density_w: C64::new(rounded(mutual[0])?, rounded(mutual[1])?),
            fringe_visibility: visibility(power, mutual)?,
            pupil_evaluations,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }
}

pub(in super::super) fn product(
    a: [Dd; 2],
    b: [Dd; 2],
) -> Result<[Dd; 2], CoherentDiffractionError> {
    let result = mul(a, b);
    for value in result {
        rounded(value)?;
    }
    if result.into_iter().all(is_zero)
        && a.into_iter().any(|value| !is_zero(value))
        && b.into_iter().any(|value| !is_zero(value))
    {
        return Err(CoherentDiffractionError::Unrepresentable);
    }
    Ok(result)
}

