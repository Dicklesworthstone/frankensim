//! Finite-area stress and energy from the original real Maxwell fields.
//!
//! The shared aperture engine integrates continuous modal products and the
//! actual shutter before observation. Only the twelve constitutive dyads and
//! six mixed products needed here are formed; no full correlation matrix,
//! spatial FFT, constitutive fit or point-sampled stress grid is introduced.

use super::super::stress::{admit_normal, finite_vector, stress_from_dyad, symmetric};
use super::*;

/// Area/shutter-averaged Maxwell stress, field energy density and energy flux.
///
/// The original observed half-space supplies positive lossless nondispersive
/// epsilon and mu. Real means and real factors retain carrier-sum and beat
/// terms; proper-complex factors retain their original second-moment law.
/// Incoming and outgoing fields are combined BEFORE the quadratic contraction.
/// Longitudinal, evanescent and off-diagonal polarization terms are retained.
///
/// This is an electromagnetic surface observation, not a mechanical body
/// force. A complete force balance needs the closed-surface stress and the
/// field-momentum rate, plus the material momentum law where applicable.
/// Opposite finite patches alone do not close a surface: their side walls
/// cannot be discarded as they can for a complete periodic cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumExposureApertureStress {
    aperture: AngularSpectrumAperture,
    stress: [Dd; 6],
    energy: Dd,
    poynting: [Dd; 3],
}

impl AngularSpectrumExposureApertureStress {
    pub const fn aperture(&self) -> AngularSpectrumAperture {
        self.aperture
    }

    /// Original Cartesian tensor: epsilon <EE> + mu <HH> - I <u> [Pa].
    pub fn mean_maxwell_stress_pa(&self) -> [[f64; 3]; 3] {
        symmetric(self.stress.map(Dd::to_f64))
    }

    /// <u> = (epsilon <E.E> + mu <H.H>)/2 [J/m³]. No cycle half is
    /// applied here: the original temporal/statistical products already own it.
    pub fn mean_field_energy_density_j_m3(&self) -> f64 {
        self.energy.to_f64()
    }

    /// Signed <E x H> [W/m²], with the actual near-field direction retained.
    pub fn mean_poynting_w_m2(&self) -> [f64; 3] {
        self.poynting.map(Dd::to_f64)
    }

    /// Continuous integral of the stress tensor over the declared rectangle
    /// [N]. Area multiplication retains DD until the public conversion.
    pub fn integrated_stress_n(&self) -> Result<[[f64; 3]; 3], AngularSpectrumError> {
        let area = Dd::from_f64(self.aperture.area_m2());
        let tensor = symmetric(self.stress.map(|v| (v * area).to_f64()));
        for row in tensor {
            finite_vector(row)?;
        }
        Ok(tensor)
    }

    /// Average traction T*n [Pa] for a finite unit normal. This contraction
    /// does not rotate the original transverse observation plane.
    pub fn mean_traction_pa(&self, normal: [f64; 3]) -> Result<[f64; 3], AngularSpectrumError> {
        self.traction(normal, Dd::ONE)
    }

    /// Integral of T*n over the rectangle [N]. For its actual z-plane surface
    /// choose n = +/-ez. Reversing the normal reverses the result.
    pub fn integrated_traction_n(
        &self,
        normal: [f64; 3],
    ) -> Result<[f64; 3], AngularSpectrumError> {
        self.traction(normal, Dd::from_f64(self.aperture.area_m2()))
    }

    /// Signed mean z-directed power through the rectangle [W].
    pub fn normal_power_w(&self) -> Result<f64, AngularSpectrumError> {
        let power = (self.poynting[2] * Dd::from_f64(self.aperture.area_m2())).to_f64();
        if !power.is_finite() {
            return Err(bad("aperture energy flux is unrepresentable"));
        }
        Ok(power)
    }

    fn traction(&self, normal: [f64; 3], area: Dd) -> Result<[f64; 3], AngularSpectrumError> {
        admit_normal(normal)?;
        finite_vector(symmetric(self.stress).map(|row| {
            (row.into_iter()
                .zip(normal)
                .fold(Dd::ZERO, |s, (t, n)| s + t * Dd::from_f64(n))
                * area)
                .to_f64()
        }))
    }
}

impl AngularSpectrumExposure<'_, '_, '_> {
    /// Finite-area, finite-shutter stress and energy in the original observed
    /// lossless nondispersive half-space. A lossless left query is valid even
    /// when the exit is lossy; a query inside that lossy exit refuses.
    pub fn measure_aperture_stress(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureStress, AngularSpectrumError> {
        let side = ExposureSide::OneSided(side);
        let sum = accumulator(self.anchor, side, &mut cancelled)?;
        measure(
            self.anchor,
            || core::iter::once((self, false)),
            side,
            distance_from_face_m,
            aperture,
            window,
            limits,
            sum,
            &mut cancelled,
        )
    }

    /// Two-face counterpart: the total field includes both propagation
    /// directions before products. Separate incoming/outgoing stresses omit
    /// their interference and cannot replace this observation.
    pub fn measure_illumination_aperture_stress(
        &self,
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureStress, AngularSpectrumError> {
        let side = ExposureSide::TwoSided(side);
        let sum = accumulator(self.anchor, side, &mut cancelled)?;
        measure(
            self.anchor,
            || core::iter::once((self, false)),
            side,
            distance_from_face_m,
            aperture,
            window,
            limits,
            sum,
            &mut cancelled,
        )
    }
}

impl AngularSpectrumExposureEnsemble<'_, '_, '_, '_> {
    /// Expected stress, energy and flux under the original mean and complete
    /// real/proper factor bank. Cross terms within each bank are retained;
    /// uncorrelated banks are never coherently added to each other.
    pub fn measure_aperture_stress(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureStress, AngularSpectrumError> {
        let side = ExposureSide::OneSided(side);
        let sum = accumulator(self.anchor, side, &mut cancelled)?;
        measure(
            self.anchor,
            || self.banks(),
            side,
            distance_from_face_m,
            aperture,
            window,
            limits,
            sum,
            &mut cancelled,
        )
    }

    /// Expected two-face finite-patch observation of the same original ensemble.
    pub fn measure_illumination_aperture_stress(
        &self,
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureStress, AngularSpectrumError> {
        let side = ExposureSide::TwoSided(side);
        let sum = accumulator(self.anchor, side, &mut cancelled)?;
        measure(
            self.anchor,
            || self.banks(),
            side,
            distance_from_face_m,
            aperture,
            window,
            limits,
            sum,
            &mut cancelled,
        )
    }
}

fn accumulator(
    anchor: &PlanarAngularSpectrum<'_, '_>,
    side: ExposureSide,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Stress, AngularSpectrumError> {
    poll(cancelled)?;
    let medium = side.momentum_material(anchor.source())?;
    let epsilon = linear_electric_permittivity(medium).map_err(ReflectionError::from)?;
    let mu = linear_magnetic_permeability(medium).map_err(ReflectionError::from)?;
    Ok(Stress::new(epsilon, mu))
}

struct Stress {
    epsilon: Dd,
    mu: Dd,
    dyad: [Dd; 6],
    poynting: [Dd; 3],
}
impl Stress {
    fn new(epsilon: f64, mu: f64) -> Self {
        Self {
            epsilon: Dd::from_f64(epsilon),
            mu: Dd::from_f64(mu),
            dyad: [Dd::ZERO; 6],
            poynting: [Dd::ZERO; 3],
        }
    }
}
impl Accumulator for Stress {
    type Output = AngularSpectrumExposureApertureStress;

    fn add<K: Kernel>(&mut self, kernel: K, a: &FieldParts, b: &FieldParts) {
        for (c, (u, v)) in [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
            .into_iter()
            .enumerate()
        {
            self.dyad[c] = self.dyad[c]
                + self.epsilon * kernel.product(a[u], b[v])
                + self.mu * kernel.product(a[u + 3], b[v + 3]);
        }
        for c in 0..3 {
            let (u, v) = ((c + 1) % 3, (c + 2) % 3);
            self.poynting[c] =
                self.poynting[c] + kernel.product(a[u], b[v + 3]) - kernel.product(a[v], b[u + 3]);
        }
    }

    fn finish(
        self,
        aperture: AngularSpectrumAperture,
    ) -> Result<Self::Output, AngularSpectrumError> {
        if self
            .dyad
            .iter()
            .chain(&self.poynting)
            .any(|&v| !finite_dd(v) || !v.to_f64().is_finite())
            || self.dyad[..3]
                .iter()
                .any(|v| v.hi < 0.0 || (v.hi == 0.0 && v.lo < 0.0))
        {
            return Err(bad("aperture stress or field energy is unrepresentable"));
        }
        let energy = (self.dyad[0] + self.dyad[1] + self.dyad[2]) * Dd::from_f64(0.5);
        let stress = stress_from_dyad(self.dyad);
        if !finite_dd(energy)
            || !energy.to_f64().is_finite()
            || stress
                .iter()
                .any(|&v| !finite_dd(v) || !v.to_f64().is_finite())
        {
            return Err(bad("aperture stress or field energy is unrepresentable"));
        }
        Ok(AngularSpectrumExposureApertureStress {
            aperture,
            stress,
            energy,
            poynting: self.poynting,
        })
    }
}

#[cfg(test)]
#[path = "stress/tests.rs"]
mod tests;
