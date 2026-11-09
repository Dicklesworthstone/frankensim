//! Continuous first spatial moments of finite-aperture Maxwell traction.
//!
//! A fringe's position contributes to angular momentum. Multiplying the
//! aperture-mean traction by a centre lever arm misses that contribution.
//! This observer integrates the first moments of the SAME modal products,
//! with exact rectangle kernels and the original shutter/statistical law.

use super::*;

/// Integrated Maxwell traction and its moment about a finite aperture centre.
/// Quantities use the original positive-z surface normal and Cartesian axes.
/// This surface observation does not itself establish a mechanical body wrench.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumExposureApertureWrench {
    aperture: AngularSpectrumAperture,
    force: [Dd; 3],
    moment: [Dd; 3],
}

impl AngularSpectrumExposureApertureWrench {
    /// The original rectangle and its physical coordinate system.
    pub const fn aperture(&self) -> AngularSpectrumAperture {
        self.aperture
    }

    /// Integral of T*ez over the original physical rectangle [N].
    /// This is surface traction with the POSITIVE z normal.
    pub fn integrated_traction_n(&self) -> [f64; 3] {
        self.force.map(Dd::to_f64)
    }

    /// Integral of (x-centre) cross (T*ez) [N m], about the rectangle centre.
    /// Both in-plane fringe lever arms are integrated before contraction.
    pub fn centered_traction_moment_nm(&self) -> [f64; 3] {
        self.moment.map(Dd::to_f64)
    }

    /// Translate the moment origin in the original Cartesian axes.
    /// offset_m is the vector FROM the aperture centre TO the new origin [m].
    /// The exact input offset is contracted in DD before public rounding.
    /// Negating BOTH returned force and moment reverses the surface normal.
    pub fn traction_moment_about_nm(
        &self,
        offset_m: [f64; 3],
    ) -> Result<[f64; 3], AngularSpectrumError> {
        finite_vector(offset_m)?;
        let value = core::array::from_fn(|axis| {
            let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
            self.moment[axis] - Dd::from_f64(offset_m[u]) * self.force[v]
                + Dd::from_f64(offset_m[v]) * self.force[u]
        });
        finish_vector(value)
    }
}

impl AngularSpectrumExposure<'_, '_, '_> {
    /// Exact continuous aperture/shutter traction and its first spatial moment.
    ///
    /// The original observed half-space must be lossless and nondispersive.
    /// Includes longitudinal fields, incoming/outgoing interference and real
    /// carrier-sum terms. The result is surface transport, NOT a body torque:
    /// a closed-surface balance and field-angular-momentum rate remain required.
    pub fn measure_aperture_wrench(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureWrench, AngularSpectrumError> {
        measure_moments(
            self.anchor,
            || core::iter::once((self, false)),
            ExposureSide::OneSided(side),
            distance_from_face_m,
            aperture,
            window,
            limits,
            &mut cancelled,
        )
    }

    /// Two-sided version; complete incoming and outgoing fields combine first.
    pub fn measure_illumination_aperture_wrench(
        &self,
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureWrench, AngularSpectrumError> {
        measure_moments(
            self.anchor,
            || core::iter::once((self, false)),
            ExposureSide::TwoSided(side),
            distance_from_face_m,
            aperture,
            window,
            limits,
            &mut cancelled,
        )
    }
}

impl AngularSpectrumExposureEnsemble<'_, '_, '_, '_> {
    /// Expected traction and moment from the original mean and real/proper bank.
    /// Distinct random factors remain separate; no source statistics are inferred.
    pub fn measure_aperture_wrench(
        &self,
        side: AngularSpectrumSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureWrench, AngularSpectrumError> {
        measure_moments(
            self.anchor,
            || self.banks(),
            ExposureSide::OneSided(side),
            distance_from_face_m,
            aperture,
            window,
            limits,
            &mut cancelled,
        )
    }

    /// Two-sided expected first moment, retaining the original coherence law.
    pub fn measure_illumination_aperture_wrench(
        &self,
        side: AngularSpectrumIlluminationSide,
        distance_from_face_m: f64,
        aperture: AngularSpectrumAperture,
        window: HarmonicWindow,
        limits: AngularSpectrumExposureApertureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularSpectrumExposureApertureWrench, AngularSpectrumError> {
        measure_moments(
            self.anchor,
            || self.banks(),
            ExposureSide::TwoSided(side),
            distance_from_face_m,
            aperture,
            window,
            limits,
            &mut cancelled,
        )
    }
}

#[allow(clippy::too_many_arguments)]
fn measure_moments<'e, 'b: 'e, 's: 'b, 'm: 's, I>(
    anchor: &PlanarAngularSpectrum<'s, 'm>,
    banks: impl Fn() -> I,
    side: ExposureSide,
    distance: f64,
    aperture: AngularSpectrumAperture,
    window: HarmonicWindow,
    limits: AngularSpectrumExposureApertureLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<AngularSpectrumExposureApertureWrench, AngularSpectrumError>
where
    I: Iterator<Item = (&'e AngularSpectrumExposure<'b, 's, 'm>, bool)>,
{
    poll(cancelled)?;
    if !distance.is_finite() || distance < 0.0 {
        return Err(bad(
            "aperture field distance must be finite and nonnegative",
        ));
    }
    let medium = accumulator(anchor, side, cancelled)?;
    let mut sum = Wrench::new(medium.epsilon, medium.mu);
    let mut frequency_pairs = 0usize;
    let mut bank_workspace = anchor.workspace_bytes();
    let (mut real, mut proper) = (false, false);
    for (bank, is_proper) in banks() {
        poll(cancelled)?;
        let f = bank.frequency_count();
        frequency_pairs = count(frequency_pairs.checked_add(count(f.checked_mul(f))?))?;
        let extra = count(anchor.modes.len().checked_mul(f).and_then(|n| {
            n.checked_mul(core::mem::size_of::<FieldParts>() - core::mem::size_of::<[C64; 6]>())
        }))?;
        bank_workspace = bank_workspace.max(count(bank.workspace_bytes().checked_add(extra))?);
        real |= f != 0 && !is_proper;
        proper |= f != 0 && is_proper;
    }
    let grid = anchor.grid();
    let extra = moment_cache_bytes(grid, real, proper)?;
    let (pair_work, workspace) = admit::<Wrench>(
        grid,
        frequency_pairs,
        count(bank_workspace.checked_add(extra))?,
        real,
        proper,
        limits,
    )?;
    let mut projection_headroom = limits.maximum_mode_frequency_pairs - pair_work;
    let workspace_headroom = limits.maximum_workspace_bytes - workspace;
    let real_cache = real
        .then(|| MomentSpatial::<RealKernel>::prepare(anchor, aperture, cancelled))
        .transpose()?;
    let proper_cache = proper
        .then(|| MomentSpatial::<ProperKernel>::prepare(anchor, aperture, cancelled))
        .transpose()?;
    for (bank, is_proper) in banks() {
        poll(cancelled)?;
        if bank.frequency_count() == 0 {
            continue;
        }
        let projected_modes = temporal::admit(
            window,
            &bank.frequencies,
            anchor.modes.len(),
            is_proper,
            &mut projection_headroom,
            workspace_headroom,
        )?;
        if is_proper {
            observe_bank(
                bank,
                proper_cache
                    .as_ref()
                    .ok_or_else(|| bad("missing proper first-moment cache"))?,
                side,
                distance,
                window,
                projected_modes,
                &mut sum,
                cancelled,
            )?;
        } else {
            observe_bank(
                bank,
                real_cache
                    .as_ref()
                    .ok_or_else(|| bad("missing real first-moment cache"))?,
                side,
                distance,
                window,
                projected_modes,
                &mut sum,
                cancelled,
            )?;
        }
    }
    let result = sum.finish(aperture)?;
    poll(cancelled)?;
    Ok(result)
}

fn moment_cache_bytes(
    grid: AngularSpectrumGrid,
    real: bool,
    proper: bool,
) -> Result<usize, AngularSpectrumError> {
    let pairs = count(grid.width.checked_mul(grid.width).and_then(|x| {
        grid.height
            .checked_mul(grid.height)
            .and_then(|y| x.checked_add(y))
    }))?;
    let kernel_bytes = usize::from(real) * core::mem::size_of::<OddKernel>()
        + usize::from(proper) * core::mem::size_of::<ProperKernel>();
    // Include headers conservatively; the base admission separately counts its
    // original spatial cache headers, grouped fields and temporal kernels.
    count(
        pairs
            .checked_mul(kernel_bytes)
            .and_then(|n| n.checked_add(core::mem::size_of::<MomentSpatial<RealKernel>>()))
            .and_then(|n| n.checked_add(core::mem::size_of::<MomentSpatial<ProperKernel>>())),
    )
}

trait FirstKernel<K: Kernel>: Copy {
    fn compose_first(self, other: K) -> Self;
    fn product_first(self, a: [Dd; 2], b: [Dd; 2]) -> Dd;
}

trait MomentKernel: Kernel {
    type First: FirstKernel<Self>;
    fn prepare_first(
        window: HarmonicWindow,
        a: Dd,
        b: Dd,
    ) -> Result<Self::First, AngularSpectrumError>;
}

// Odd weighted Gram entries <xi*cos(a xi)*sin(b xi)> and its transpose.
// The even/even and odd/odd entries vanish on the centered interval.
#[derive(Clone, Copy)]
struct OddKernel([Dd; 2]);

impl FirstKernel<RealKernel> for OddKernel {
    fn compose_first(self, other: RealKernel) -> Self {
        let [p, q] = self.0;
        let [c, s] = other.0;
        Self([p * c - q * s, q * c - p * s])
    }
    fn product_first(self, a: [Dd; 2], b: [Dd; 2]) -> Dd {
        a[0] * b[1] * self.0[0] + a[1] * b[0] * self.0[1]
    }
}

impl MomentKernel for RealKernel {
    type First = OddKernel;
    fn prepare_first(
        window: HarmonicWindow,
        a: Dd,
        b: Dd,
    ) -> Result<OddKernel, AngularSpectrumError> {
        let difference = centered_sine_moment(a - b, window.duration_s())?;
        let sum = centered_sine_moment(a + b, window.duration_s())?;
        Ok(OddKernel([
            (sum - difference) * Dd::from_f64(0.5),
            (sum + difference) * Dd::from_f64(0.5),
        ]))
    }
}

impl FirstKernel<ProperKernel> for ProperKernel {
    fn compose_first(self, other: Self) -> Self {
        self.compose(other)
    }
    fn product_first(self, a: [Dd; 2], b: [Dd; 2]) -> Dd {
        self.product(a, b)
    }
}

impl MomentKernel for ProperKernel {
    type First = ProperKernel;
    fn prepare_first(window: HarmonicWindow, a: Dd, b: Dd) -> Result<Self, AngularSpectrumError> {
        let difference = a - b;
        let moment = centered_sine_moment(difference, window.duration_s())?;
        let start = difference * window.start_clock_dd();
        let half_width = half_product(difference, window.duration_s());
        if !finite_dd(start) || !finite_dd(half_width) {
            return Err(bad("first-moment phase is unrepresentable"));
        }
        // Keep a tiny aperture phase out of a distant DD coordinate, just as
        // the original real-field center rotation does.
        let (s, c) = fs_math::det::sin_cos_dd(start);
        let (t, d) = fs_math::det::sin_cos_dd(half_width);
        let (s, c) = (s * d + c * t, c * d - s * t);
        // <xi exp(-i q x)> = -i <xi sin(q xi)> exp(-i q centre).
        Ok(Self([-moment * s, -moment * c]))
    }
}

fn half_product(a: Dd, b: f64) -> Dd {
    if a.hi.abs() >= b.abs() {
        (a * Dd::from_f64(0.5)) * Dd::from_f64(b)
    } else {
        a * Dd::from_f64(0.5 * b)
    }
}

// <xi sin(k xi)> = L/2 * (sin(z)-z*cos(z))/z^2, z=k*L/2.
// A DD series avoids the cancellation at small width; split divisions avoid
// forming z^2 for large z. No finite-difference derivative or spatial sampling.
fn centered_sine_moment(k: Dd, width: f64) -> Result<Dd, AngularSpectrumError> {
    if !finite_dd(k) || !width.is_finite() || width <= 0.0 {
        return Err(bad("invalid aperture first-moment interval"));
    }
    if k.hi == 0.0 && k.lo == 0.0 {
        return Ok(Dd::ZERO);
    }
    let z = half_product(k, width);
    if !finite_dd(z) || (z.hi == 0.0 && z.lo == 0.0) {
        return Err(bad("aperture first-moment phase is unrepresentable"));
    }
    let ratio = if z.hi.abs() <= 0.5 {
        let zz = z * z;
        let mut term = z / Dd::from_f64(3.0);
        let mut sum = term;
        for n in 0..13 {
            term = -term * zz * Dd::from_f64((n + 2) as f64)
                / Dd::from_f64(((n + 1) * (2 * n + 4) * (2 * n + 5)) as f64);
            sum = sum + term;
        }
        sum
    } else {
        let (s, c) = fs_math::det::sin_cos_dd(z);
        (s / z - c) / z
    };
    let value = half_product(ratio, width);
    if !finite_dd(value) || !value.to_f64().is_finite() || (value.hi == 0.0 && value.lo == 0.0) {
        return Err(bad("aperture first-moment weight is unrepresentable"));
    }
    Ok(value)
}

struct MomentSpatial<K: MomentKernel> {
    base: Spatial<K>,
    x: Vec<K::First>,
    y: Vec<K::First>,
}

impl<K: MomentKernel> MomentSpatial<K> {
    fn prepare(
        anchor: &PlanarAngularSpectrum<'_, '_>,
        aperture: AngularSpectrumAperture,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        Self::from_axes(
            anchor.grid(),
            aperture,
            |i| transverse_wavenumber_dd(anchor, 0, i),
            |i| transverse_wavenumber_dd(anchor, 1, i),
            cancelled,
        )
    }

    fn from_axes(
        grid: AngularSpectrumGrid,
        aperture: AngularSpectrumAperture,
        kx: impl Fn(usize) -> Dd,
        ky: impl Fn(usize) -> Dd,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AngularSpectrumError> {
        let base = Spatial::from_axes_dd(grid, aperture, &kx, &ky, cancelled)?;
        let axis = |count: usize,
                    which: usize,
                    frequency: &dyn Fn(usize) -> Dd,
                    cancelled: &mut dyn FnMut() -> bool|
         -> Result<Vec<K::First>, AngularSpectrumError> {
            let window =
                HarmonicWindow::new(aperture.origin_m()[which], aperture.extent_m()[which])
                    .map_err(|_| bad("aperture first-moment interval is unrepresentable"))?;
            let length = self::count(count.checked_mul(count))?;
            let mut values = Vec::new();
            values
                .try_reserve_exact(length)
                .map_err(|_| bad("first-moment allocation failed"))?;
            for i in 0..count {
                for j in 0..count {
                    if j % 32 == 0 && cancelled() {
                        return Err(AngularSpectrumError::Cancelled);
                    }
                    values.push(K::prepare_first(window, frequency(i), frequency(j))?);
                }
            }
            Ok(values)
        };
        let x = axis(grid.width, 0, &kx, cancelled)?;
        let y = axis(grid.height, 1, &ky, cancelled)?;
        Ok(Self { base, x, y })
    }
}

#[allow(clippy::too_many_arguments)]
fn observe_bank<K: MomentKernel>(
    bank: &AngularSpectrumExposure<'_, '_, '_>,
    spatial: &MomentSpatial<K>,
    side: ExposureSide,
    distance: f64,
    window: HarmonicWindow,
    projected_modes: Option<usize>,
    sum: &mut Wrench,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AngularSpectrumError> {
    let f = bank.frequency_count();
    let mut fields = Vec::new();
    fields
        .try_reserve_exact(f)
        .map_err(|_| bad("first-moment frequency field allocation failed"))?;
    for group in 0..f {
        let mut field = bank.frequency_modal_parts(group, side, distance, cancelled)?;
        if K::CENTER_FIELDS {
            center_fields_dd(
                &mut field,
                bank.frequencies[group],
                |i| transverse_wavenumber_dd(bank.anchor, 0, i % bank.anchor.grid.width),
                |i| transverse_wavenumber_dd(bank.anchor, 1, i / bank.anchor.grid.width),
                spatial.base.aperture,
                window,
                cancelled,
            )?;
        }
        fields.push(field);
    }
    if let Some(modes) = projected_modes {
        let points = spatial.base.grid.width * spatial.base.grid.height;
        let fields =
            temporal::project::<K>(window, &bank.frequencies, &fields, points, modes, cancelled)?;
        integrate_projected_moments(spatial, &fields, modes, sum, cancelled)
    } else {
        let temporal = pair_kernels(f, window, |i| -bank.frequencies[i], cancelled)?;
        integrate_moments(spatial, &temporal, &fields, sum, cancelled)
    }
}

fn integrate_projected_moments<K: MomentKernel>(
    spatial: &MomentSpatial<K>,
    fields: &[FieldParts],
    modes: usize,
    sum: &mut Wrench,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AngularSpectrumError> {
    let grid = spatial.base.grid;
    let points = grid.width * grid.height;
    for i in 0..points {
        poll(cancelled)?;
        for j in 0..points {
            let x = (i % grid.width) * grid.width + j % grid.width;
            let y = (i / grid.width) * grid.height + j / grid.width;
            let normal = spatial.base.x[x].compose(spatial.base.y[y]);
            let first_x = spatial.x[x].compose_first(spatial.base.y[y]);
            let first_y = spatial.y[y].compose_first(spatial.base.x[x]);
            for mode in 0..modes {
                poll(cancelled)?;
                let time = K::temporal_mode(mode);
                let (a, b) = (&fields[mode * points + i], &fields[mode * points + j]);
                sum.add(normal.compose(time), a, b);
                sum.add_first::<K>(0, first_x.compose_first(time), a, b);
                sum.add_first::<K>(1, first_y.compose_first(time), a, b);
            }
        }
    }
    Ok(())
}

fn integrate_moments<K: MomentKernel>(
    spatial: &MomentSpatial<K>,
    temporal: &[K],
    fields: &[Vec<FieldParts>],
    sum: &mut Wrench,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), AngularSpectrumError> {
    let grid = spatial.base.grid;
    let n = grid.width * grid.height;
    let f = fields.len();
    for i in 0..n {
        poll(cancelled)?;
        for j in 0..n {
            let x = (i % grid.width) * grid.width + j % grid.width;
            let y = (i / grid.width) * grid.height + j / grid.width;
            let normal = spatial.base.x[x].compose(spatial.base.y[y]);
            let first_x = spatial.x[x].compose_first(spatial.base.y[y]);
            let first_y = spatial.y[y].compose_first(spatial.base.x[x]);
            for (a, left) in fields.iter().enumerate() {
                for (b, right) in fields.iter().enumerate() {
                    if b % 16 == 0 {
                        poll(cancelled)?;
                    }
                    let time = temporal[a * f + b];
                    sum.add(normal.compose(time), &left[i], &right[j]);
                    sum.add_first::<K>(0, first_x.compose_first(time), &left[i], &right[j]);
                    sum.add_first::<K>(1, first_y.compose_first(time), &left[i], &right[j]);
                }
            }
        }
    }
    Ok(())
}

struct Wrench {
    epsilon: Dd,
    mu: Dd,
    force: [Dd; 3],
    first: [[Dd; 3]; 2],
}

impl Wrench {
    fn new(epsilon: Dd, mu: Dd) -> Self {
        Self {
            epsilon,
            mu,
            force: [Dd::ZERO; 3],
            first: [[Dd::ZERO; 3]; 2],
        }
    }

    fn traction(
        &self,
        product: impl Fn([Dd; 2], [Dd; 2]) -> Dd,
        a: &FieldParts,
        b: &FieldParts,
    ) -> [Dd; 3] {
        let dyad = |u: usize, v: usize| {
            self.epsilon * product(a[u], b[v]) + self.mu * product(a[u + 3], b[v + 3])
        };
        [
            dyad(0, 2),
            dyad(1, 2),
            (dyad(2, 2) - dyad(0, 0) - dyad(1, 1)) * Dd::from_f64(0.5),
        ]
    }

    fn add_first<K: Kernel>(
        &mut self,
        axis: usize,
        kernel: impl FirstKernel<K>,
        a: &FieldParts,
        b: &FieldParts,
    ) {
        let value = self.traction(|a, b| kernel.product_first(a, b), a, b);
        for (total, v) in self.first[axis].iter_mut().zip(value) {
            *total = *total + v;
        }
    }
}

impl Accumulator for Wrench {
    type Output = AngularSpectrumExposureApertureWrench;
    fn add<K: Kernel>(&mut self, kernel: K, a: &FieldParts, b: &FieldParts) {
        let value = self.traction(|a, b| kernel.product(a, b), a, b);
        for (total, v) in self.force.iter_mut().zip(value) {
            *total = *total + v;
        }
    }
    fn finish(
        self,
        aperture: AngularSpectrumAperture,
    ) -> Result<Self::Output, AngularSpectrumError> {
        let area = Dd::from_f64(aperture.area_m2());
        let force = self.force.map(|v| v * area);
        let moment = [
            self.first[1][2] * area,
            -self.first[0][2] * area,
            (self.first[0][1] - self.first[1][0]) * area,
        ];
        finish_vector(force)?;
        finish_vector(moment)?;
        Ok(AngularSpectrumExposureApertureWrench {
            aperture,
            force,
            moment,
        })
    }
}

fn count(value: Option<usize>) -> Result<usize, AngularSpectrumError> {
    value.ok_or_else(|| bad("aperture first-moment work/storage overflow"))
}

fn finish_vector(value: [Dd; 3]) -> Result<[f64; 3], AngularSpectrumError> {
    if value.iter().any(|&v| !finite_dd(v)) {
        return Err(bad("aperture traction moment is unrepresentable"));
    }
    finite_vector(value.map(Dd::to_f64))
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod temporal_tests;
