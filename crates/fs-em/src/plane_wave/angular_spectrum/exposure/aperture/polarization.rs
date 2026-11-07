//! Real finite-shutter correlations and real or complex field analyzers.
//!
//! These are uncentered moments of the REAL instantaneous fields, retaining
//! both carrier sum and difference terms where the original source law requires
//! them. They are not an analytic-signal coherency matrix: equal-time real
//! moments alone do not determine circular handedness or a Stokes vector.
//! Direct complex analyzers act on the original phasors before real products,
//! retaining the phase information needed to select circular/elliptical fields.

use super::*;

mod spectral;
pub use spectral::AngularSpectrumExposureApertureJonesSpectrum;

/// Complete real second moments of E and H over the authored area and shutter.
///
/// EE and HH are symmetric; EH[i][j] is <E_i H_j>, whose antisymmetric part
/// determines the mean Poynting vector. Means have not been subtracted, so
/// these are second moments rather than centered source covariance. Components
/// retain the original Cartesian basis, including longitudinal/evanescent fields.
/// No impedance, transverse-plane or optical-cycle approximation is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumExposureApertureMoments {
    parts: [Dd; 21],
    reading: AngularSpectrumExposureApertureReading,
}

impl AngularSpectrumExposureApertureMoments {
    pub fn electric_second_moment_v2_m2(&self) -> [[f64; 3]; 3] {
        self.block(0, 0)
    }
    pub fn magnetic_second_moment_a2_m2(&self) -> [[f64; 3]; 3] {
        self.block(3, 3)
    }
    pub fn electric_magnetic_second_moment_w_m2(&self) -> [[f64; 3]; 3] {
        self.block(0, 3)
    }
    /// Original diagonal field moments, signed vector flux and normal watts.
    pub fn reading(&self) -> AngularSpectrumExposureApertureReading {
        self.reading
    }
    /// Contract retained DD moments with real, dimensionless analyzer weights.
    /// For very deep analyzer nulls, the source owner's direct projection query
    /// projects fields before products and avoids subtracting large moments.
    /// Neither route gives a general relative-error guarantee near every null.
    pub fn project(
        &self,
        projection: AngularSpectrumExposureApertureProjection,
    ) -> Result<AngularSpectrumExposureApertureProjectionReading, AngularSpectrumError> {
        let mut sums = [Dd::ZERO; 3];
        for i in 0..3 {
            for j in 0..3 {
                let (e, h) = (
                    Dd::from_f64(projection.electric[i]),
                    Dd::from_f64(projection.magnetic[i]),
                );
                sums[0] =
                    sums[0] + e * self.parts[index(i, j)] * Dd::from_f64(projection.electric[j]);
                sums[1] = sums[1]
                    + h * self.parts[index(i + 3, j + 3)] * Dd::from_f64(projection.magnetic[j]);
                sums[2] = sums[2]
                    + e * self.parts[index(i, j + 3)] * Dd::from_f64(projection.magnetic[j]);
            }
        }
        projection_reading(sums)
    }

    fn block(&self, row: usize, column: usize) -> [[f64; 3]; 3] {
        core::array::from_fn(|i| {
            core::array::from_fn(|j| self.parts[index(row + i, column + j)].to_f64())
        })
    }
}

/// Two authored REAL linear field observations: e = a.E and h = b.H.
///
/// Weights are finite and dimensionless, with no implicit normalization. A unit
/// vector selects a Cartesian field component along that axis; zero disables a
/// channel. Scaling a weight scales its signal and squares its mean square.
/// This describes field projection, not a material polarizer, complex Jones
/// filter, detector responsivity, absorption or output optical power.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumExposureApertureProjection {
    electric: [f64; 3],
    magnetic: [f64; 3],
}
impl AngularSpectrumExposureApertureProjection {
    pub fn try_new(electric: [f64; 3], magnetic: [f64; 3]) -> Result<Self, AngularSpectrumError> {
        if electric.iter().chain(&magnetic).any(|v| !v.is_finite()) {
            return Err(bad("aperture analyzer weights must be finite"));
        }
        Ok(Self { electric, magnetic })
    }
    pub fn electric_weights(self) -> [f64; 3] {
        self.electric
    }
    pub fn magnetic_weights(self) -> [f64; 3] {
        self.magnetic
    }
}

/// Complex field analyzers applied to each original positive-frequency phasor.
///
/// The two real observed signals are `Re(sum_j a† E_j exp(i omega_j t))`
/// and `Re(sum_j b† H_j exp(i omega_j t))`. The dagger CONJUGATES the
/// supplied weights, matching the monochromatic interference owner. A
/// quadrature weight therefore distinguishes opposite circular or elliptical
/// fields that have identical real, cycle-averaged second moments.
///
/// Weights are finite and dimensionless, with no implicit normalization. They
/// describe an explicitly authored ideal, frequency-independent observation
/// over the supplied spectral lines, not a causal broadband filter, material
/// polarizer, Stokes reconstruction, absorbed detector power or back-action.
/// The original real/proper source statistics and finite shutter remain exact;
/// carrier sums are not silently discarded for deterministic/real sources.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumExposureApertureJonesProjection {
    electric: [C64; 3],
    magnetic: [C64; 3],
}

impl AngularSpectrumExposureApertureJonesProjection {
    pub fn try_new(electric: [C64; 3], magnetic: [C64; 3]) -> Result<Self, AngularSpectrumError> {
        if electric.iter().chain(&magnetic).any(|&v| !finite(v)) {
            return Err(bad("aperture Jones analyzer weights must be finite"));
        }
        Ok(Self { electric, magnetic })
    }

    pub fn electric_weights(self) -> [C64; 3] {
        self.electric
    }

    pub fn magnetic_weights(self) -> [C64; 3] {
        self.magnetic
    }
}

impl From<AngularSpectrumExposureApertureProjection>
    for AngularSpectrumExposureApertureJonesProjection
{
    fn from(value: AngularSpectrumExposureApertureProjection) -> Self {
        Self {
            electric: value.electric.map(C64::from_re),
            magnetic: value.magnetic.map(C64::from_re),
        }
    }
}

/// Area/shutter averages of e², h² and e*h for the two observed real signals.
/// The mixed field moment retains its sign; it is not a Poynting component.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumExposureApertureProjectionReading {
    pub electric_mean_squared_v2_m2: f64,
    pub magnetic_mean_squared_a2_m2: f64,
    pub electric_magnetic_mean_product_w_m2: f64,
}

impl AngularSpectrumExposure<'_, '_, '_> {
    /// Observe all 21 independent real second moments with the original
    /// finite-area/shutter kernel. Each ordered modal/frequency pair needs 21
    /// scalar products; the ordinary aperture reading retains its cheaper path.
    pub fn measure_aperture_moments(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureMoments, AngularSpectrumError> {
        measure(
            self.anchor,
            || core::iter::once((self, false)),
            ExposureSide::OneSided(side),
            distance_from_face_m,
            aperture,
            window,
            limits,
            Full([Dd::ZERO; 21]),
            &mut cancelled,
        )
    }

    /// Full real correlations of incoming, outgoing or total fields at either
    /// original face, preserving the source's actual two-port phase relation.
    pub fn measure_illumination_aperture_moments(
        &self,
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureMoments, AngularSpectrumError> {
        measure(
            self.anchor,
            || core::iter::once((self, false)),
            ExposureSide::TwoSided(side),
            distance_from_face_m,
            aperture,
            window,
            limits,
            Full([Dd::ZERO; 21]),
            &mut cancelled,
        )
    }

    /// Project each original grouped DD field BEFORE quadratic products.
    /// This retains weak analyzer residuals that can disappear when contracting
    /// an already rounded correlation matrix. The two real signals require
    /// three products per modal/frequency pair and no new source solve. Accepts
    /// real weights or complex Jones weights with the authored a†E convention.
    #[allow(clippy::too_many_arguments)]
    pub fn measure_aperture_projection(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        projection: impl Into<AngularSpectrumExposureApertureJonesProjection>,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureProjectionReading, AngularSpectrumError> {
        measure(
            self.anchor,
            || core::iter::once((self, false)),
            ExposureSide::OneSided(side),
            distance_from_face_m,
            aperture,
            window,
            limits,
            Projected::new(projection),
            &mut cancelled,
        )
    }

    /// Direct real or complex analyzers at either original physical face.
    /// They read the native fields without replacing a material boundary.
    #[allow(clippy::too_many_arguments)]
    pub fn measure_illumination_aperture_projection(
        &self,
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        projection: impl Into<AngularSpectrumExposureApertureJonesProjection>,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureProjectionReading, AngularSpectrumError> {
        measure(
            self.anchor,
            || core::iter::once((self, false)),
            ExposureSide::TwoSided(side),
            distance_from_face_m,
            aperture,
            window,
            limits,
            Projected::new(projection),
            &mut cancelled,
        )
    }
}

impl AngularSpectrumExposureEnsemble<'_, '_, '_, '_> {
    /// Full real correlations of the optional mean and authored real/proper
    /// factors; within-factor interference remains and factor moments add.
    pub fn measure_aperture_moments(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureMoments, AngularSpectrumError> {
        measure(
            self.anchor,
            || self.banks(),
            ExposureSide::OneSided(side),
            distance_from_face_m,
            aperture,
            window,
            limits,
            Full([Dd::ZERO; 21]),
            &mut cancelled,
        )
    }
    /// Full real correlations for a selected incoming/outgoing/total face.
    pub fn measure_illumination_aperture_moments(
        &self,
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureMoments, AngularSpectrumError> {
        measure(
            self.anchor,
            || self.banks(),
            ExposureSide::TwoSided(side),
            distance_from_face_m,
            aperture,
            window,
            limits,
            Full([Dd::ZERO; 21]),
            &mut cancelled,
        )
    }
    /// Expected real-signal moments for authored real or complex analyzers.
    /// Every original factor is projected before products; real factors retain
    /// sum terms and proper factors retain beats.
    #[allow(clippy::too_many_arguments)]
    pub fn measure_aperture_projection(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        projection: impl Into<AngularSpectrumExposureApertureJonesProjection>,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureProjectionReading, AngularSpectrumError> {
        measure(
            self.anchor,
            || self.banks(),
            ExposureSide::OneSided(side),
            distance_from_face_m,
            aperture,
            window,
            limits,
            Projected::new(projection),
            &mut cancelled,
        )
    }
    /// Expected real-signal moments at either original physical face, with
    /// authored real or complex analyzer weights.
    #[allow(clippy::too_many_arguments)]
    pub fn measure_illumination_aperture_projection(
        &self,
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        projection: impl Into<AngularSpectrumExposureApertureJonesProjection>,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureProjectionReading, AngularSpectrumError> {
        measure(
            self.anchor,
            || self.banks(),
            ExposureSide::TwoSided(side),
            distance_from_face_m,
            aperture,
            window,
            limits,
            Projected::new(projection),
            &mut cancelled,
        )
    }
}

fn index(a: usize, b: usize) -> usize {
    let (row, column) = (a.min(b), a.max(b));
    row * (13 - row) / 2 + column - row
}

pub(super) struct Full(pub(super) [Dd; 21]);
impl Accumulator for Full {
    type Output = AngularSpectrumExposureApertureMoments;
    fn add<K: Kernel>(&mut self, kernel: K, a: &FieldParts, b: &FieldParts) {
        for i in 0..6 {
            for j in i..6 {
                self.0[index(i, j)] = self.0[index(i, j)] + kernel.product(a[i], b[j]);
            }
        }
    }
    fn finish(
        self,
        aperture: AngularSpectrumAperture,
    ) -> Result<Self::Output, AngularSpectrumError> {
        if self.0.iter().any(|&v| !finite_dd(v)) {
            return Err(bad("aperture field correlation is unrepresentable"));
        }
        let mut basic = [Dd::ZERO; 9];
        for c in 0..3 {
            basic[c] = self.0[index(c, c)];
            basic[c + 3] = self.0[index(c + 3, c + 3)];
            let (u, v) = ((c + 1) % 3, (c + 2) % 3);
            basic[c + 6] = self.0[index(u, v + 3)] - self.0[index(v, u + 3)];
        }
        Ok(AngularSpectrumExposureApertureMoments {
            parts: self.0,
            reading: basic.finish(aperture)?,
        })
    }
}

trait ProjectionLaw {
    fn admit_frequencies(
        &self,
        _angular_frequencies_rad_s: &[f64],
        _cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AngularSpectrumError> {
        Ok(())
    }
    fn at(
        &self,
        angular_frequency_rad_s: f64,
    ) -> Result<(AngularSpectrumExposureApertureJonesProjection, bool), AngularSpectrumError>;
}

struct ConstantProjection {
    projection: AngularSpectrumExposureApertureJonesProjection,
    real: bool,
}
impl ProjectionLaw for ConstantProjection {
    fn at(
        &self,
        _angular_frequency_rad_s: f64,
    ) -> Result<(AngularSpectrumExposureApertureJonesProjection, bool), AngularSpectrumError> {
        Ok((self.projection, self.real))
    }
}

fn real_weights(projection: AngularSpectrumExposureApertureJonesProjection) -> bool {
    projection
        .electric
        .iter()
        .chain(&projection.magnetic)
        .all(|z| z.im == 0.0)
}

struct Projected<P = ConstantProjection> {
    projection: P,
    sums: [Dd; 3],
}
impl Projected {
    fn new(projection: impl Into<AngularSpectrumExposureApertureJonesProjection>) -> Self {
        let projection = projection.into();
        Self {
            projection: ConstantProjection {
                projection,
                real: real_weights(projection),
            },
            sums: [Dd::ZERO; 3],
        }
    }
}
impl<P: ProjectionLaw> Accumulator for Projected<P> {
    type Output = AngularSpectrumExposureApertureProjectionReading;
    fn admit_frequencies(
        &self,
        angular_frequencies_rad_s: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AngularSpectrumError> {
        self.projection
            .admit_frequencies(angular_frequencies_rad_s, cancelled)
    }
    fn prepare_field(
        &self,
        angular_frequency_rad_s: f64,
        field: &mut [FieldParts],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), AngularSpectrumError> {
        let (projection, real) = self.projection.at(angular_frequency_rad_s)?;
        for (i, value) in field.iter_mut().enumerate() {
            if i % 32 == 0 {
                poll(cancelled)?;
            }
            let (mut e, mut h) = ([Dd::ZERO; 2], [Dd::ZERO; 2]);
            for component in 0..3 {
                for part in 0..2 {
                    e[part] = e[part]
                        + Dd::from_f64(projection.electric[component].re) * value[component][part];
                    h[part] = h[part]
                        + Dd::from_f64(projection.magnetic[component].re)
                            * value[component + 3][part];
                }
                if !real {
                    // Conjugate analyzer weights before the quadratic kernel.
                    // Keep DD low parts through crossed-circular dark ports.
                    let ei = Dd::from_f64(projection.electric[component].im);
                    let hi = Dd::from_f64(projection.magnetic[component].im);
                    e[0] = e[0] + ei * value[component][1];
                    e[1] = e[1] - ei * value[component][0];
                    h[0] = h[0] + hi * value[component + 3][1];
                    h[1] = h[1] - hi * value[component + 3][0];
                }
            }
            if e.iter().chain(&h).any(|&v| !finite_dd(v)) {
                return Err(bad("aperture analyzer field is unrepresentable"));
            }
            *value = [
                e,
                [Dd::ZERO; 2],
                [Dd::ZERO; 2],
                h,
                [Dd::ZERO; 2],
                [Dd::ZERO; 2],
            ];
        }
        Ok(())
    }
    fn add<K: Kernel>(&mut self, kernel: K, a: &FieldParts, b: &FieldParts) {
        self.sums[0] = self.sums[0] + kernel.product(a[0], b[0]);
        self.sums[1] = self.sums[1] + kernel.product(a[3], b[3]);
        self.sums[2] = self.sums[2] + kernel.product(a[0], b[3]);
    }
    fn finish(
        self,
        _aperture: AngularSpectrumAperture,
    ) -> Result<Self::Output, AngularSpectrumError> {
        projection_reading(self.sums)
    }
}

fn projection_reading(
    sums: [Dd; 3],
) -> Result<AngularSpectrumExposureApertureProjectionReading, AngularSpectrumError> {
    if sums.iter().any(|&v| !finite_dd(v)) || sums[..2].iter().any(|v| v.to_f64() < 0.0) {
        return Err(bad("aperture analyzer moment is unrepresentable"));
    }
    Ok(AngularSpectrumExposureApertureProjectionReading {
        electric_mean_squared_v2_m2: sums[0].to_f64(),
        magnetic_mean_squared_a2_m2: sums[1].to_f64(),
        electric_magnetic_mean_product_w_m2: sums[2].to_f64(),
    })
}

#[cfg(test)]
mod tests;
