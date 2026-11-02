//! Finite exposures of an authored proper-complex spectral covariance.
//!
//! A source is its deterministic mean plus sum_r z_r B_r(omega, x), with
//! E[z_r]=0, E[z_r z_s*]=delta_rs and E[z_r z_s]=0. These second moments are
//! explicit inputs through the factor fields, not sampled or inferred phases.
//! Properness removes the covariance's sum-frequency terms even for a short
//! shutter; correlations between distinct frequencies still produce beats.
//! Properness alone implies neither Gaussianity nor stationarity. See Neeser
//! and Massey, IEEE Trans. Information Theory 39 (1993), equations (2)-(5):
//! <https://www.isiweb.ee.ethz.ch/papers/arch/fneese-mass-inspec-1993-1.pdf>.

use super::*;

/// Whole-ensemble work limits, including the optional deterministic mean.
/// Workspace conservatively includes the largest original exposure query,
/// the ensemble's compensated accumulation and the final output. Other
/// borrowed banks/beam owners, earlier outputs and allocator overhead are
/// outside this logical numerical-storage bound.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumExposureEnsembleLimits {
    pub maximum_factors: usize,
    pub maximum_point_line_pairs: usize,
    pub maximum_point_frequency_pairs: usize,
    pub maximum_workspace_bytes: usize,
    pub maximum_fft_work: u64,
}

/// Mean illumination plus mutually uncorrelated proper-complex factors.
///
/// Each borrowed exposure bank is one factor across its authored frequencies.
/// Its fields include their square-root statistical weights. Equal-frequency
/// lines interfere within the factor through the original Maxwell operator;
/// different factors add second moments and need not be orthogonal. The
/// resulting covariance is positive semidefinite by construction. In symbols,
/// W(omega,x;omega',x') = sum_r B_r(omega,x) B_r(omega',x')*.
/// Only the stated joint second moments are required; statistical independence
/// and a particular probability distribution are not assumed.
///
/// An optional deterministic mean retains the existing exact real-product
/// exposure, including carrier sum terms. Its moments add to the zero-mean
/// covariance without a mean/fluctuation cross term. Cross-frequency source
/// correlations can be nonstationary; this API does not infer a stationary
/// spectrum, source bandwidth, pseudocovariance or an individual realization.
pub struct AngularSpectrumExposureEnsemble<'e, 'b, 's, 'm> {
    pub(super) anchor: &'b PlanarAngularSpectrum<'s, 'm>,
    pub(super) mean: Option<&'e AngularSpectrumExposure<'b, 's, 'm>>,
    pub(super) factors: &'e [AngularSpectrumExposure<'b, 's, 'm>],
    line_count: usize,
    point_line_pairs: usize,
    point_frequency_pairs: usize,
    workspace_bytes: usize,
    field_fft_work: u64,
}

impl<'e, 'b, 's, 'm> AngularSpectrumExposureEnsemble<'e, 'b, 's, 'm> {
    /// The anchor fixes the original material owners, films and grid without
    /// adding illumination. Empty factors and no mean represent darkness.
    /// Preparation is O(R) for R already prepared covariance-factor banks;
    /// no source fields, frequency groups or covariance matrices are copied.
    pub fn try_new(
        anchor: &'b PlanarAngularSpectrum<'s, 'm>,
        mean: Option<&'e AngularSpectrumExposure<'b, 's, 'm>>,
        factors: &'e [AngularSpectrumExposure<'b, 's, 'm>],
        limits: AngularSpectrumExposureEnsembleLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        poll(&mut cancelled)?;
        if factors.len() > 4096 || factors.len() > limits.maximum_factors {
            return Err(bad("exposure covariance factor count exceeds its budget"));
        }
        let mut line_count = 0_usize;
        let mut point_frequency_pairs = 0_usize;
        let mut field_fft_work = 0_u64;
        let mut bank_workspace = anchor.workspace_bytes();
        for (bank, proper) in mean
            .into_iter()
            .map(|bank| (bank, false))
            .chain(factors.iter().map(|bank| (bank, true)))
        {
            poll(&mut cancelled)?;
            if bank.anchor.grid() != anchor.grid()
                || !bank
                    .anchor
                    .source()
                    .same_sources_and_geometry(anchor.source())
            {
                return Err(bad(
                    "exposure ensemble requires the same original materials, layers and grid",
                ));
            }
            line_count = line_count
                .checked_add(bank.line_count())
                .ok_or_else(|| bad("exposure ensemble line work overflow"))?;
            point_frequency_pairs = point_frequency_pairs
                .checked_add(bank.point_frequency_pairs())
                .ok_or_else(|| bad("exposure ensemble pair work overflow"))?;
            field_fft_work = field_fft_work
                .checked_add(bank.field_fft_work())
                .ok_or_else(|| bad("exposure ensemble FFT work overflow"))?;
            // Existing bank budgets contain HarmonicProduct storage. Account
            // for any larger proper kernel without assuming its representation.
            let extra_kernel_bytes = if proper {
                let extra = core::mem::size_of::<ProperHarmonicProduct>()
                    .saturating_sub(core::mem::size_of::<HarmonicProduct>());
                bank.frequency_count()
                    .checked_mul(bank.frequency_count())
                    .and_then(|v| v.checked_mul(extra))
                    .ok_or_else(|| bad("exposure ensemble workspace overflow"))?
            } else {
                0
            };
            bank_workspace = bank_workspace.max(
                bank.workspace_bytes()
                    .checked_add(extra_kernel_bytes)
                    .ok_or_else(|| bad("exposure ensemble workspace overflow"))?,
            );
        }
        let n = anchor.modes.len();
        let point_line_pairs = n
            .checked_mul(line_count)
            .ok_or_else(|| bad("exposure ensemble line work overflow"))?;
        let workspace_bytes = n
            .checked_mul(
                core::mem::size_of::<[Dd; 9]>()
                    + core::mem::size_of::<AngularSpectrumExposurePoint>(),
            )
            .and_then(|v| v.checked_add(bank_workspace))
            .ok_or_else(|| bad("exposure ensemble workspace overflow"))?;
        if point_line_pairs > limits.maximum_point_line_pairs
            || point_frequency_pairs > limits.maximum_point_frequency_pairs
            || field_fft_work > limits.maximum_fft_work
            || workspace_bytes > limits.maximum_workspace_bytes
        {
            return Err(bad(
                "exposure ensemble exceeds its work or workspace budget",
            ));
        }
        Ok(Self {
            anchor,
            mean,
            factors,
            line_count,
            point_line_pairs,
            point_frequency_pairs,
            workspace_bytes,
            field_fft_work,
        })
    }

    pub fn grid(&self) -> AngularSpectrumGrid {
        self.anchor.grid()
    }
    pub fn factor_count(&self) -> usize {
        self.factors.len()
    }
    /// Total authored lines in the covariance factors and optional mean.
    pub fn line_count(&self) -> usize {
        self.line_count
    }
    pub fn point_line_pairs(&self) -> usize {
        self.point_line_pairs
    }
    pub fn point_frequency_pairs(&self) -> usize {
        self.point_frequency_pairs
    }
    pub fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    pub fn field_fft_work(&self) -> u64 {
        self.field_fft_work
    }

    /// Ensemble and finite-time averages through the original film responses.
    /// Stream one factor at a time, retaining O(N*max F_r + max F_r² + N)
    /// query storage. Work is O(N + R + sum_r(L_r*N + F_r*N*log N + F_r²*N)),
    /// including the optional mean as another bank. Temporal kernels are
    /// cached once per factor/frequency pair, never per spatial point.
    pub fn field(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        window: HarmonicWindow,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureField, AngularSpectrumError> {
        self.field_on(
            ExposureSide::OneSided(side),
            distance_from_face_m,
            window,
            &mut cancelled,
        )
    }

    /// Joint left/right source covariance through the original two-port fields.
    /// Mean and factors may each be one- or two-sided; all retain the same
    /// physical materials, grid and common source clock. Faces are correlated
    /// WITHIN each factor, never treated as independent by their location.
    pub fn illumination_field(
        &self,
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        window: HarmonicWindow,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureField, AngularSpectrumError> {
        self.field_on(
            ExposureSide::TwoSided(side),
            distance_from_face_m,
            window,
            &mut cancelled,
        )
    }

    fn field_on(
        &self,
        side: ExposureSide,
        distance_from_face_m: f64,
        window: HarmonicWindow,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureField, AngularSpectrumError> {
        poll(&mut cancelled)?;
        if !distance_from_face_m.is_finite() || distance_from_face_m < 0.0 {
            return Err(bad("field distance must be finite and nonnegative"));
        }
        let n = self.anchor.modes.len();
        let mut sums = Vec::new();
        sums.try_reserve_exact(n)
            .map_err(|_| bad("exposure ensemble accumulation allocation failed"))?;
        for i in 0..n {
            if i % 64 == 0 {
                poll(&mut cancelled)?;
            }
            sums.push([Dd::ZERO; 9]);
        }
        if let Some(mean) = self.mean {
            accumulate_bank(
                mean,
                side,
                distance_from_face_m,
                |a, b| window.product_kernel(a, b),
                product,
                &mut sums,
                &mut cancelled,
            )?;
        }
        for factor in self.factors {
            accumulate_bank(
                factor,
                side,
                distance_from_face_m,
                |a, b| window.proper_product_kernel(a, b),
                proper_product,
                &mut sums,
                &mut cancelled,
            )?;
        }
        let mut points = Vec::new();
        points
            .try_reserve_exact(n)
            .map_err(|_| bad("exposure ensemble result allocation failed"))?;
        let mut total = Dd::ZERO;
        for (i, sum) in sums.into_iter().enumerate() {
            if i % 64 == 0 {
                poll(&mut cancelled)?;
            }
            let value = sum.map(Dd::to_f64);
            if value.iter().any(|x| !x.is_finite()) || value[..6].iter().any(|x| *x < 0.0) {
                return Err(bad("exposure ensemble second moment is unrepresentable"));
            }
            total = total + sum[8];
            points.push(AngularSpectrumExposurePoint {
                electric_mean_squared_v2_m2: [value[0], value[1], value[2]],
                magnetic_mean_squared_a2_m2: [value[3], value[4], value[5]],
                poynting_w_m2: [value[6], value[7], value[8]],
            });
        }
        let grid = self.grid();
        let normal_power_w = total.to_f64() * grid.spacing_m[0] * grid.spacing_m[1];
        if !normal_power_w.is_finite() {
            return Err(bad("integrated exposure ensemble power is unrepresentable"));
        }
        poll(&mut cancelled)?;
        Ok(AngularSpectrumExposureField {
            grid,
            points,
            normal_power_w,
        })
    }
}

// Monomorphized for deterministic/proper products: no statistical-mode branch
// in the spatial/frequency contraction, and no second Maxwell implementation.
fn accumulate_bank<P: Copy>(
    bank: &AngularSpectrumExposure<'_, '_, '_>,
    side: ExposureSide,
    distance: f64,
    kernel: impl Fn(f64, f64) -> Result<P, fs_fft::harmonic_window::WindowError>,
    product: impl Fn(P, C64, C64) -> Result<Dd, MomentError>,
    sums: &mut [[Dd; 9]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AngularSpectrumError> {
    poll(cancelled)?;
    if bank.line_count() == 0 {
        return Ok(());
    }
    let f = bank.frequency_count();
    let mut kernels = Vec::new();
    kernels
        .try_reserve_exact(f * f)
        .map_err(|_| bad("exposure ensemble kernel allocation failed"))?;
    for &a in &bank.frequencies {
        for &b in &bank.frequencies {
            poll(cancelled)?;
            kernels.push(
                kernel(a, b)
                    .map_err(|_| bad("exposure ensemble temporal kernel is unrepresentable"))?,
            );
        }
    }
    let mut fields = Vec::new();
    fields
        .try_reserve_exact(f)
        .map_err(|_| bad("exposure ensemble field allocation failed"))?;
    for group in 0..f {
        fields.push(bank.frequency_field(group, side, distance, cancelled)?);
    }
    for (point, sum) in sums.iter_mut().enumerate() {
        poll(cancelled)?;
        let value = integrate_point_dd(f, &kernels, |i| fields[i][point], &product, cancelled)
            .map_err(|error| match error {
                MomentError::Cancelled => AngularSpectrumError::Cancelled,
                MomentError::Product(_) => {
                    bad("exposure ensemble field product is unrepresentable")
                }
                MomentError::Unrepresentable => {
                    bad("exposure ensemble second moment is unrepresentable")
                }
            })?;
        for (sum, value) in sum.iter_mut().zip(value) {
            *sum = *sum + value;
        }
    }
    Ok(())
}
