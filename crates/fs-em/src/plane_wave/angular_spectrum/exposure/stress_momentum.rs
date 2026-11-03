//! Joint exterior stress and field-momentum observation.
//!
//! A transient momentum balance needs both the Maxwell stress and the field
//! storage term. These observations share one original frequency bank,
//! temporal-kernel cache and FFT synthesis; their physical contractions and
//! constitutive scaling remain owned by the individual stress/momentum paths.
//! The result is not a mechanical force: a consumer still owns its enclosing
//! surface/volume integration and any material momentum or striction law.

use super::momentum::MomentumKernel;
use super::*;

/// Complete joint-query numerical storage. The original reconstruction/FFT
/// workspace is counted once, together with both compensated accumulation
/// planes, both outputs, the momentum frequency coefficients and modal scratch.
/// Borrowed owners, earlier outputs and allocator overhead are excluded.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumExposureStressMomentumLimits {
    pub maximum_workspace_bytes: usize,
}

/// Stress, Minkowski field momentum and its weighted instantaneous rate from
/// the same original source, exterior plane and normalized observation window.
/// Both outputs retain original point samples and continuous periodic-cell
/// integrals, including signed Nyquist modes without sampled-product aliasing.
#[derive(Debug)]
pub struct AngularSpectrumExposureStressMomentumField {
    stress: AngularSpectrumExposureStressField,
    momentum: AngularSpectrumExposureMomentumField,
}

impl AngularSpectrumExposureStressMomentumField {
    pub fn grid(&self) -> AngularSpectrumGrid {
        self.stress.grid()
    }

    /// Maxwell stress/traction in Pa and continuous transverse integrals in N.
    pub fn stress(&self) -> &AngularSpectrumExposureStressField {
        &self.stress
    }

    /// Mean D cross B and mean d(D cross B)/dt, in N s/m^3 and N/m^3,
    /// with continuous transverse integrals in N s/m and N/m. In vacuum
    /// the momentum is S/c^2. A material mechanical momentum is not inferred.
    pub fn momentum(&self) -> &AngularSpectrumExposureMomentumField {
        &self.momentum
    }

    /// Transfer both original result owners without cloning their point planes.
    pub fn into_fields(
        self,
    ) -> (
        AngularSpectrumExposureStressField,
        AngularSpectrumExposureMomentumField,
    ) {
        (self.stress, self.momentum)
    }
}

impl AngularSpectrumExposure<'_, '_, '_> {
    pub fn stress_momentum_workspace_bytes(&self) -> Result<usize, AngularSpectrumError> {
        stress::stress_workspace(self.anchor.modes.len(), self.momentum_workspace_bytes()?)
    }

    /// Observe Maxwell stress, field momentum and its rate together from the
    /// original incoming/reflected/total/transmitted field. The observed side
    /// supplies its own original scalar, nondispersive lossless material.
    ///
    /// One frequency reconstruction, one temporal-kernel cache and one FFT
    /// synthesis per frequency serve all observables. The original N*F^2 work
    /// admission remains authoritative: each real pair uses thirty local
    /// products and at most sixty modal products; a proper pair uses at most
    /// twenty-four local and forty-eight modal products. Work remains
    /// O(L*N + F*N*log N + F^2*N), with no new film or constitutive solve.
    pub fn stress_momentum_field(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureStressMomentumLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureStressMomentumField, AngularSpectrumError> {
        stress_momentum_field(
            self.anchor,
            core::iter::once((self, false)),
            self.stress_momentum_workspace_bytes()?,
            ExposureSide::OneSided(side),
            distance_from_face_m,
            window,
            limits,
            &mut cancelled,
        )
    }

    /// Observe one incoming, outgoing or TOTAL field at either physical face.
    /// Both propagation directions interfere before any stress or momentum
    /// contraction; separate incoming/outgoing observables cannot replace it.
    pub fn illumination_stress_momentum_field(
        &self,
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureStressMomentumLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureStressMomentumField, AngularSpectrumError> {
        stress_momentum_field(
            self.anchor,
            core::iter::once((self, false)),
            self.stress_momentum_workspace_bytes()?,
            ExposureSide::TwoSided(side),
            distance_from_face_m,
            window,
            limits,
            &mut cancelled,
        )
    }
}

impl AngularSpectrumExposureEnsemble<'_, '_, '_, '_> {
    pub fn stress_momentum_workspace_bytes(&self) -> Result<usize, AngularSpectrumError> {
        stress::stress_workspace(self.anchor.modes.len(), self.momentum_workspace_bytes()?)
    }

    /// Observe the authored mean and mutually uncorrelated real/proper
    /// covariance factors. One original bank is retained at a time and serves
    /// all three observables. Only the authored proper law removes carrier
    /// sum terms; no stationary or independent-face source is inferred.
    pub fn stress_momentum_field(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureStressMomentumLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureStressMomentumField, AngularSpectrumError> {
        stress_momentum_field(
            self.anchor,
            self.banks(),
            self.stress_momentum_workspace_bytes()?,
            ExposureSide::OneSided(side),
            distance_from_face_m,
            window,
            limits,
            &mut cancelled,
        )
    }

    pub fn illumination_stress_momentum_field(
        &self,
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureStressMomentumLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureStressMomentumField, AngularSpectrumError> {
        stress_momentum_field(
            self.anchor,
            self.banks(),
            self.stress_momentum_workspace_bytes()?,
            ExposureSide::TwoSided(side),
            distance_from_face_m,
            window,
            limits,
            &mut cancelled,
        )
    }
}

fn stress_momentum_field<'e, 'b, 's, 'm>(
    anchor: &PlanarAngularSpectrum<'s, 'm>,
    banks: impl IntoIterator<Item = (&'e AngularSpectrumExposure<'b, 's, 'm>, bool)>,
    workspace: usize,
    side: ExposureSide,
    distance: f64,
    window: HarmonicWindow,
    limits: AngularSpectrumExposureStressMomentumLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<AngularSpectrumExposureStressMomentumField, AngularSpectrumError>
where
    'b: 'e,
    's: 'b,
    'm: 's,
{
    poll(cancelled)?;
    if workspace > limits.maximum_workspace_bytes {
        return Err(bad(
            "exposure stress and momentum exceed their workspace budget",
        ));
    }
    if !distance.is_finite() || distance < 0.0 {
        return Err(bad(
            "stress and momentum field distance must be finite and nonnegative",
        ));
    }
    let medium = side.momentum_material(anchor.source())?;
    let epsilon = linear_electric_permittivity(medium).map_err(ReflectionError::from)?;
    let mu = linear_magnetic_permeability(medium).map_err(ReflectionError::from)?;
    let n = anchor.modes.len();
    let mut stress_sums = Vec::new();
    let mut momentum_sums = Vec::new();
    stress_sums
        .try_reserve_exact(n)
        .map_err(|_| bad("exposure joint stress accumulation allocation failed"))?;
    momentum_sums
        .try_reserve_exact(n)
        .map_err(|_| bad("exposure joint momentum accumulation allocation failed"))?;
    for i in 0..n {
        if i % 64 == 0 {
            poll(cancelled)?;
        }
        stress_sums.push([Dd::ZERO; 6]);
        momentum_sums.push([Dd::ZERO; 6]);
    }
    let mut stress_integrated = [Dd::ZERO; 6];
    let mut momentum_integrated = [Dd::ZERO; 6];
    for (bank, proper) in banks {
        poll(cancelled)?;
        let (stress, momentum) = if proper {
            accumulate_joint_bank(
                bank,
                side,
                distance,
                epsilon,
                mu,
                true,
                |a, b| window.proper_product_kernel(a, b),
                proper_product,
                ProperHarmonicProduct::mean_real_product_parts_dd,
                &mut stress_sums,
                &mut momentum_sums,
                cancelled,
            )?
        } else {
            accumulate_joint_bank(
                bank,
                side,
                distance,
                epsilon,
                mu,
                false,
                |a, b| window.product_kernel(a, b),
                product,
                HarmonicProduct::mean_real_product_parts_dd,
                &mut stress_sums,
                &mut momentum_sums,
                cancelled,
            )?
        };
        for c in 0..6 {
            stress_integrated[c] = stress_integrated[c] + stress[c];
            momentum_integrated[c] = momentum_integrated[c] + momentum[c];
        }
    }
    let stress =
        stress::finish_stress_field(anchor.grid(), stress_sums, stress_integrated, cancelled)?;
    let momentum = momentum::finish_momentum_field(
        anchor.grid(),
        momentum_sums,
        momentum_integrated,
        epsilon,
        mu,
        cancelled,
    )?;
    poll(cancelled)?;
    Ok(AngularSpectrumExposureStressMomentumField { stress, momentum })
}

fn accumulate_joint_bank<P: Copy>(
    bank: &AngularSpectrumExposure<'_, '_, '_>,
    side: ExposureSide,
    distance: f64,
    epsilon: f64,
    mu: f64,
    proper: bool,
    kernel: impl Fn(f64, f64) -> Result<P, fs_fft::harmonic_window::WindowError>,
    product: impl Fn(P, C64, C64) -> Result<Dd, MomentError>,
    product_parts: impl Fn(P, [Dd; 2], [Dd; 2]) -> Result<Dd, fs_fft::harmonic_window::WindowError>,
    stress_sums: &mut [[Dd; 6]],
    momentum_sums: &mut [[Dd; 6]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<([Dd; 6], [Dd; 6]), AngularSpectrumError> {
    let f = bank.frequency_count();
    if f == 0 {
        return Ok(([Dd::ZERO; 6], [Dd::ZERO; 6]));
    }
    let mut observation = ExposureObservationBank::prepare(
        bank,
        side,
        distance,
        |a, b| kernel(a, b).map(|product| MomentumKernel::new(a, b, product)),
        cancelled,
    )?;
    // Stress projects the original product out of the shared momentum kernel.
    // No second cache or phase-kernel evaluation is needed.
    let stress_cell = stress::integrate_cell_dyad(
        bank.anchor.grid(),
        &observation.fields,
        &observation.kernels,
        epsilon,
        mu,
        |kernel, a, b| product_parts(kernel.product, a, b),
        cancelled,
    )?;
    let momentum_cell = momentum::integrate_cell_momentum(
        bank.anchor.grid(),
        &observation.fields,
        &observation.kernels,
        proper,
        &product_parts,
        cancelled,
    )?;
    observation.synthesize(bank, cancelled)?;
    for (point, (stress_sum, momentum_sum)) in stress_sums
        .iter_mut()
        .zip(momentum_sums.iter_mut())
        .enumerate()
    {
        poll(cancelled)?;
        let stress = stress::integrate_dyad(
            f,
            &observation.kernels,
            |i| observation.fields[i][point],
            epsilon,
            mu,
            |kernel, a, b| product(kernel.product, a, b),
            cancelled,
        )
        .map_err(stress::dyad_error)?;
        let momentum = momentum::integrate_momentum(
            f,
            &observation.kernels,
            |i| observation.fields[i][point],
            proper,
            momentum::rotate_i,
            &product,
            cancelled,
        )
        .map_err(momentum::momentum_error)?;
        for c in 0..6 {
            stress_sum[c] = stress_sum[c] + stress[c];
            momentum_sum[c] = momentum_sum[c] + momentum[c];
        }
    }
    Ok((stress_cell, momentum_cell))
}
