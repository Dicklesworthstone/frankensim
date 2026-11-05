//! Complete modal loading of a planar radiator in an infinite rigid baffle.
//!
//! This projects the same half-space Rayleigh kernel and equivalent-disc self
//! integral as `baffled_piston_impedance`. It does not solve a closed-body BEM
//! problem. Positive area cells carry signed, work-conjugate velocity shapes;
//! every retained field contributes to the full reciprocal impedance matrix.
//! Fields use `exp(-i omega t)`: a mass-like load has negative imaginary part.
//!
//! `power_consistent_impedance` instead uses the regular resistance kernel's
//! exact diagonal limit at the SAME quadrature points as outgoing pressure.
//! Its resistance is a hemispherical far-field Gram matrix (PSD in exact
//! arithmetic), including signed cancelling fields. The singular reactance
//! retains the equivalent-disc self integral; this is NOT a complete causal
//! passivity certificate or a finite-cell Galerkin quadrature-error bound.
use super::{HelmholtzError, Medium, MIN_PANELS_PER_WAVELENGTH, green};
use fs_math::{c64::C64, dd::Dd, det};
use core::f64::consts::PI;
use std::sync::OnceLock;

#[path = "rayleigh_modal_resistance.rs"]
mod resistance;

/// Complete structural velocity fields retained by the loading integral.
/// Downstream passive realizations may project these fields into fewer acoustic
/// ports; that reduction must not remove coordinates from this physical matrix.
pub const MAX_MODAL_RAYLEIGH_FIELDS: usize = 128;

/// Area cells retained by one prepared surface. No dense factorization happens on
/// this path (the dense Helmholtz LU keeps its own `MAX_DENSE_PANELS`): storage is
/// the strict pair-distance triangle, about 1.07 GB at this cap, plus one
/// `cells * modes` complex pressure field per frequency being evaluated.
pub const MAX_MODAL_RAYLEIGH_CELLS: usize = 16_384;

/// Bound on `cells^2 * modes + cells * modes^2` per frequency, plus the bounded
/// source-moment/Gram work when low-frequency resistance needs it: a cost
/// ceiling, not an accuracy gate. It admits the complete 128-field image up to about
/// 12,000 cells (a board resolved to 2.4 kHz needs about 9,400); larger surfaces
/// retain fewer fields. A caller preparing many frequencies must additionally
/// bound its complete sweep. No cells or fields are silently discarded.
pub const MAX_MODAL_RAYLEIGH_WORK: u64 = 20_000_000_000;

/// Highest polar order admitted by direct outgoing-power integration. The
/// azimuth rule has twice this order. An unresolved angular field refuses;
/// it is never replaced by a negative or clipped matrix quadratic form.
pub const MAX_RAYLEIGH_POWER_POLAR_ORDER: usize = 128;

/// Mean power of the coherent point-area field into one radiating half-space.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentRadiatedPower {
    /// Mean outgoing power [W] for peak modal velocity phasors.
    pub power_w: f64,
    /// Analytic angular truncation bound [W], in exact arithmetic. This does
    /// not enclose floating-point error, finite-cell quadrature error, or
    /// uncertainty in geometry, source phases or material properties.
    pub angular_truncation_bound_w: f64,
    /// Number of distinct directions evaluated after planar symmetries.
    pub angular_directions: usize,
}

/// Admit a declared number of per-frequency Rayleigh scalar-work units.
/// This is also the allocation-free boundary used by sweep planners.
pub fn admit_modal_work(work: u128) -> Result<(), HelmholtzError> {
    if work > u128::from(MAX_MODAL_RAYLEIGH_WORK) {
        return Err(invalid("Rayleigh loading exceeds the per-frequency modal work limit"));
    }
    Ok(())
}

/// One positive-area quadrature cell on the radiator's plane, in SI units.
#[derive(Clone, Debug)]
pub struct AreaSample {
    /// Coordinates in the common baffle plane [m].
    pub position_m: [f64; 2],
    /// Physical integration weight [m^2], included once in each work integral.
    pub area_m2: f64,
    /// Normal velocity per generalized modal velocity; all cells cover all modes.
    pub mode_shape: Vec<f64>,
}

/// Geometry prepared once for a bounded frequency sweep; no spatial truncation.
pub struct PreparedSurface {
    count: usize,
    modes: usize,
    areas: Vec<f64>,
    radii: Vec<f64>,
    // Keep the area product's low component until physical source coordinates
    // have been combined. Rounding each nearly identical mode first would
    // manufacture radiation when those modes interfere destructively.
    weighted_shape: Vec<Dd>,
    positions: Vec<[f64; 2]>,
    /// Signed source moments are projected once into direction/degree
    /// coefficients for reuse across low-frequency sweeps. Every modal
    /// coordinate is retained; this only changes evaluation of the sinc kernel.
    resistance_coefficients: OnceLock<Result<Vec<Dd>, HelmholtzError>>,
    /// Strict upper triangle, in row order: about 1.07 GB at the cell cap.
    distances: Vec<f64>,
    maximum_radius_m: f64,
    maximum_distance_m: f64,
}

/// One coherent motion of the retained surface, reusable across angular grids.
/// Modal velocities have already been combined at each physical area cell;
/// subsequent observation costs O(cells * directions), independent of mode count.
pub struct PreparedCoherentRadiation<'a> {
    surface: &'a PreparedSurface,
    flow: Vec<[Dd; 2]>,
    total_flow: [Dd; 2],
    k: f64,
    pressure_scale: f64,
}

impl PreparedCoherentRadiation<'_> {
    /// Complex far-field amplitude F [Pa m], with p(r d) -> F(d) exp(i k r)/r
    /// on the same exp(-i omega t) clock as the modal loading. All directions
    /// must be unit vectors. The infinite rigid baffle radiates into z >= 0;
    /// z < 0 returns zero, so this evaluator also fits the existing
    /// `directivity_sh_table_from_evaluator_with_options` full-sphere adapter.
    /// Its angular convergence gates still apply to the baffle-plane jump.
    ///
    /// Relative source phases are evaluated about a retained cell, then the
    /// common origin phase is restored. Cosine-minus-one is formed without
    /// subtracting near-unit values, retaining compact cancelling source fields.
    /// No modal power, absolute value or far-receiver distance approximation is
    /// substituted for a coherent pressure. This is the SAME point-area model
    /// as `power_consistent_impedance`, not a finite-cell integration bound.
    ///
    /// # Errors
    /// Invalid directions, explicit complete-request kernel budget, the shared
    /// scalar-work ceiling, allocation, cancellation or nonfinite arithmetic.
    /// Cancellation is polled every 64 directions during validation before
    /// allocation, then at each direction and each 64-cell tile. A failed
    /// request publishes no partial image. Empty direction lists are admitted.
    pub fn far_field(
        &self,
        directions: &[[f64; 3]],
        maximum_kernel_evaluations: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<C64>, HelmholtzError> {
        if cancelled() { return Err(HelmholtzError::Cancelled); }
        let work = directions.len().checked_mul(self.surface.count);
        if work.is_none_or(|n| n > maximum_kernel_evaluations) {
            return Err(HelmholtzError::FarFieldBudget {
                required_evaluations: work.unwrap_or(usize::MAX),
                maximum_evaluations: maximum_kernel_evaluations,
            });
        }
        admit_modal_work(8 * work.unwrap_or(0) as u128)?;
        for (index, d) in directions.iter().enumerate() {
            if index % 64 == 0 && cancelled() { return Err(HelmholtzError::Cancelled); }
            if d.iter().any(|v| !v.is_finite())
                || (det::hypot(det::hypot(d[0], d[1]), d[2]) - 1.0).abs() > 1e-10 {
                return Err(invalid("Rayleigh far-field directions must be finite unit vectors"));
            }
        }
        let mut output = reserved(directions.len(), "coherent Rayleigh far-field image")?;
        let origin = self.surface.positions[0];
        let k = Dd::from_f64(-self.k);
        for d in directions {
            if cancelled() { return Err(HelmholtzError::Cancelled); }
            if d[2] < 0.0 {
                output.push(C64::ZERO);
                continue;
            }
            let direction = [Dd::from_f64(d[0]), Dd::from_f64(d[1])];
            let reference_phase = (k * (direction[0] * Dd::from_f64(origin[0])
                + direction[1] * Dd::from_f64(origin[1]))).to_f64();
            if !reference_phase.is_finite() {
                return Err(invalid("Rayleigh far-field phase exceeds finite numerical range"));
            }
            let mut sum = self.total_flow;
            for (cell, (point, flow)) in self.surface.positions.iter().zip(&self.flow).enumerate() {
                if cell % 64 == 0 && cancelled() { return Err(HelmholtzError::Cancelled); }
                let phase = (k * (direction[0] * (Dd::from_f64(point[0]) - Dd::from_f64(origin[0]))
                    + direction[1] * (Dd::from_f64(point[1]) - Dd::from_f64(origin[1])))).to_f64();
                if !phase.is_finite() {
                    return Err(invalid("Rayleigh far-field phase exceeds finite numerical range"));
                }
                let sine = Dd::from_f64(det::sin(phase));
                let half_sine = Dd::from_f64(det::sin(0.5 * phase));
                let cosine_minus_one = Dd::from_f64(-2.0) * half_sine * half_sine;
                sum[0] = sum[0] + flow[0] * cosine_minus_one - flow[1] * sine;
                sum[1] = sum[1] + flow[1] * cosine_minus_one + flow[0] * sine;
            }
            let (sine, cosine) = det::sin_cos(reference_phase);
            let (sine, cosine) = (Dd::from_f64(sine), Dd::from_f64(cosine));
            let scale = Dd::from_f64(self.pressure_scale);
            let value = C64::new(
                ((sum[0] * sine + sum[1] * cosine) * scale).to_f64(),
                (-(sum[0] * cosine - sum[1] * sine) * scale).to_f64(),
            );
            if !value.re.is_finite() || !value.im.is_finite() {
                return Err(invalid("Rayleigh far-field pressure exceeds finite numerical range"));
            }
            output.push(value);
        }
        if cancelled() { return Err(HelmholtzError::Cancelled); }
        Ok(output)
    }
}

fn invalid(what: &'static str) -> HelmholtzError {
    HelmholtzError::BadParameter { what }
}

fn reserved<T>(capacity: usize, what: &'static str) -> Result<Vec<T>, HelmholtzError> {
    let mut values = Vec::new();
    values.try_reserve_exact(capacity)
        .map_err(|_| HelmholtzError::AllocationFailed { what })?;
    Ok(values)
}

impl PreparedSurface {
    /// Retain every supplied cell and mode. Coincident cells refuse: their
    /// distinct physical self integrals cannot be inferred from area alone.
    /// The complete geometry and one impedance evaluation use about 1.1 GB at
    /// the cell cap; every reservation is fallible.
    /// The base modal work cap is checked before geometry storage is allocated;
    /// low-frequency moment work is also admitted before its first evaluation.
    pub fn new(samples: &[AreaSample]) -> Result<Self, HelmholtzError> {
        let count = samples.len();
        if count == 0 || count > MAX_MODAL_RAYLEIGH_CELLS {
            return Err(HelmholtzError::WorkCap { panels: count });
        }
        let modes = samples[0].mode_shape.len();
        if !(1..=MAX_MODAL_RAYLEIGH_FIELDS).contains(&modes) {
            return Err(invalid("Rayleigh loading requires 1..128 complete modal fields"));
        }
        let cells = count as u64;
        let fields = modes as u64;
        admit_modal_work(u128::from(cells) * u128::from(fields) * u128::from(cells + fields))?;
        let mut areas = reserved(count, "complete Rayleigh cell areas")?;
        let mut radii = reserved(count, "complete Rayleigh cell radii")?;
        let mut weighted_shape = reserved(count * modes, "complete Rayleigh weighted shapes")?;
        let mut positions = reserved(count, "complete Rayleigh cell positions")?;
        let mut maximum_radius_m = 0.0_f64;
        for cell in samples {
            if !cell.area_m2.is_finite() || cell.area_m2 <= 0.0
                || cell.position_m.iter().any(|v| !v.is_finite())
                || cell.mode_shape.len() != modes
                || cell.mode_shape.iter().any(|v| !v.is_finite()) {
                return Err(invalid("Rayleigh cells need finite coordinates, positive areas and complete finite shapes"));
            }
            let radius = det::sqrt(cell.area_m2 / PI);
            if !radius.is_finite() || radius <= 0.0 {
                return Err(invalid("Rayleigh cell radius exceeds finite numerical range"));
            }
            areas.push(cell.area_m2);
            radii.push(radius);
            positions.push(cell.position_m);
            maximum_radius_m = maximum_radius_m.max(radius);
            for &shape in &cell.mode_shape {
                let weighted = Dd::from_f64(cell.area_m2) * Dd::from_f64(shape);
                if !weighted.hi.is_finite() || !weighted.lo.is_finite()
                    || (shape != 0.0 && weighted.to_f64() == 0.0) {
                    return Err(invalid("Rayleigh area-weighted shape exceeds finite numerical range"));
                }
                weighted_shape.push(weighted);
            }
        }
        let pairs = count * (count - 1) / 2;
        let mut distances = reserved(pairs, "complete Rayleigh pair distances")?;
        let mut maximum_distance_m = 0.0_f64;
        for i in 0..count {
            for j in i + 1..count {
                let a = samples[i].position_m;
                let b = samples[j].position_m;
                let r = (a[0] - b[0]).hypot(a[1] - b[1]);
                if !r.is_finite() || r <= 0.0 {
                    return Err(invalid("Rayleigh area cells must have distinct, finitely separated centres"));
                }
                maximum_distance_m = maximum_distance_m.max(r);
                distances.push(r);
            }
        }
        Ok(Self { count, modes, areas, radii, weighted_shape, positions,
            resistance_coefficients: OnceLock::new(), distances, maximum_radius_m, maximum_distance_m })
    }

    /// Largest distance between two cell centres [m]. The kernel's longest travel
    /// time is this over the sound speed, which sets how finely its impedance varies
    /// in frequency (detail on a scale of c / distance, at every frequency).
    pub fn maximum_distance_m(&self) -> f64 { self.maximum_distance_m }

    /// Wavelength divided by the largest equivalent-disc diameter. This is an
    /// area-resolution diagnostic, not a bound on shape, aspect-ratio or
    /// quadrature error. The same six-cells-per-wavelength floor as the BEM
    /// owner is enforced before a sampled impedance is formed.
    pub fn equivalent_disc_ppw(&self, omega: f64, medium: Medium) -> Result<f64, HelmholtzError> {
        if !omega.is_finite() || omega <= 0.0 || !medium.density.is_finite()
            || medium.density <= 0.0 || !medium.sound_speed.is_finite() || medium.sound_speed <= 0.0 {
            return Err(invalid("Rayleigh frequency, density and sound speed must be positive and finite"));
        }
        let k = omega / medium.sound_speed;
        let diameter_in_wavelengths = k * self.maximum_radius_m / PI;
        let ppw = 1.0 / diameter_in_wavelengths;
        if !k.is_finite() || k <= 0.0 || !ppw.is_finite() || ppw <= 0.0 {
            return Err(invalid("Rayleigh wavelength exceeds finite numerical range"));
        }
        Ok(ppw)
    }

    /// Generalized force per generalized velocity for ONE radiating half-space.
    /// This is `Phi^T A K A Phi`, with the Rayleigh `2G` kernel off the diagonal
    /// and the existing equivalent-disc source integral on each diagonal.
    /// Both off-diagonal directions share one kernel evaluation. Projecting
    /// the pressure fields costs O(cells^2*modes + cells*modes^2), not one
    /// spatial double integral per matrix entry. Each reciprocal entry is
    /// evaluated once; data are never averaged or clipped to appear passive.
    pub fn impedance(&self, omega: f64, medium: Medium) -> Result<Vec<C64>, HelmholtzError> {
        self.impedance_image(omega, medium, false)
    }

    /// Loading consistent with the power of the outgoing area-point field.
    /// For `b_i = area_i * shape_i`, the resistance kernel is
    /// `rho*omega^2/(2*pi*c) * sinc(k*distance_ij)`, including `sinc(0)=1`.
    /// For peak complex velocity amplitudes, half the Hermitian quadratic
    /// form of its real part equals mean outgoing hemisphere pressure power.
    /// The sinc matrix is the angular Gram matrix of
    /// `exp(-i*k*direction.position_i)` on this planar surface.
    ///
    /// Integrating the self-cell resistance over a disc but keeping all other
    /// resistance entries at points subtracts a positive diagonal from that
    /// Gram matrix. Cancelling fields can then radiate NEGATIVE power, even
    /// on wavelength-resolved meshes. This image integrates the smooth real
    /// kernel consistently; it does not clip eigenvalues or add empirical loss.
    /// The original [`Self::impedance`] image remains available for comparison.
    ///
    /// At `k * maximum_distance <= 0.5`, a compensated source-moment expansion
    /// evaluates this same angular Gram matrix before squaring. This preserves
    /// weak multipole radiation that a difference of near-unit sinc entries
    /// loses to roundoff. The degree-16 exponential remainder is below 2.2e-20
    /// times the sum of absolute area-weighted shapes per angular amplitude;
    /// the positive angular rule integrates the polynomial products exactly.
    /// This is an absolute truncation estimate, not relative accuracy for an
    /// arbitrarily dark mode or a remedy for ill-conditioned quadratic forms
    /// subsequently assembled from the returned matrix.
    ///
    /// Singular self reactance, work bounds, all cells/fields and reciprocal
    /// projection are unchanged. Positive resistance alone does not establish
    /// a causal positive-real time model: consumers must still pass their
    /// existing passive-realization and frequency-fit gates. Floating-point
    /// nullspace roundoff is not removed or presented as extra physical loss.
    pub fn power_consistent_impedance(&self, omega: f64, medium: Medium)
        -> Result<Vec<C64>, HelmholtzError> {
        self.impedance_image(omega, medium, true)
    }

    /// Evaluate coherent outgoing power without first rounding a resistance
    /// matrix. Every modal velocity is a peak phasor on the same exp(-i omega t)
    /// clock. Mode amplitudes are combined BEFORE squaring; independent modal
    /// powers must not be added for a coherent excitation.
    ///
    /// Compact sources reuse the cached degree-16 directional coefficients,
    /// contracting modes at each degree before frequency evaluation. Larger
    /// sources use positive angular quadrature of the combined cell field;
    /// its order is chosen from the exponential Taylor remainder, up to
    /// [`MAX_RAYLEIGH_POWER_POLAR_ORDER`]. Both use compensated arithmetic and
    /// retain signed cross terms without a loss floor. The returned analytic
    /// bound is absolute, not a relative guarantee for arbitrarily dark fields.
    ///
    /// This evaluates the SAME area-point radiation as
    /// [`Self::power_consistent_impedance`], without assembling its reactance
    /// or solving a boundary problem. It neither changes loading nor certifies
    /// a causal time model. Work is admitted before allocation. Cancellation
    /// is checked at entry and each angular direction; initial compact-source
    /// cache construction has the existing bounded, noninterruptible projection.
    pub fn coherent_radiated_power(
        &self,
        omega: f64,
        medium: Medium,
        velocity: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentRadiatedPower, HelmholtzError> {
        resistance::coherent_power(self, omega, medium, velocity, &mut cancelled)
    }

    /// Prepare a complex far-field observer of this same loaded surface motion.
    /// Peak modal velocity phasors are combined into physical cell volume flows
    /// with retained double-double area products before any propagation or
    /// angular observation. All signed modes and relative phases are retained.
    /// This supplies directional pressure where `coherent_radiated_power`
    /// supplies only an angular integral; it performs no boundary solve.
    ///
    /// Source contraction costs O(cells * modes) once and stores two
    /// double-doubles per cell. Reuse the result for directivity refinement.
    /// The six-cell-per-wavelength admission is the existing area-resolution
    /// diagnostic, not a spatial quadrature or finite-baffle accuracy promise.
    /// # Errors
    /// Invalid medium/frequency, incomplete or nonfinite modal motion, spatial
    /// resolution, work, allocation, cancellation or nonfinite cell volume flow.
    pub fn prepare_coherent_radiation(
        &self,
        omega: f64,
        medium: Medium,
        velocity: &[C64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedCoherentRadiation<'_>, HelmholtzError> {
        if cancelled() { return Err(HelmholtzError::Cancelled); }
        let ppw = self.equivalent_disc_ppw(omega, medium)?;
        if ppw < MIN_PANELS_PER_WAVELENGTH {
            return Err(HelmholtzError::TooCoarse { panels_per_wavelength: ppw });
        }
        if velocity.len() != self.modes
            || velocity.iter().any(|v| !v.re.is_finite() || !v.im.is_finite()) {
            return Err(invalid("Rayleigh far-field requires every finite modal velocity"));
        }
        admit_modal_work(self.count as u128 * self.modes as u128)?;
        let pressure_scale = medium.density * (omega / (2.0 * PI));
        if !pressure_scale.is_finite() || pressure_scale <= 0.0 {
            return Err(invalid("Rayleigh far-field scale exceeds finite numerical range"));
        }
        let mut flow = reserved(self.count, "coherent Rayleigh cell volume velocities")?;
        let mut total_flow = [Dd::ZERO; 2];
        for cell in 0..self.count {
            if cell % 64 == 0 && cancelled() { return Err(HelmholtzError::Cancelled); }
            let mut value = [Dd::ZERO; 2];
            for (mode, v) in velocity.iter().enumerate() {
                let weight = self.weighted_shape[cell * self.modes + mode];
                value[0] = value[0] + weight * Dd::from_f64(v.re);
                value[1] = value[1] + weight * Dd::from_f64(v.im);
            }
            if value.iter().any(|x| !x.hi.is_finite() || !x.lo.is_finite()) {
                return Err(invalid("Rayleigh cell volume flow exceeds finite numerical range"));
            }
            for axis in 0..2 { total_flow[axis] = total_flow[axis] + value[axis]; }
            flow.push(value);
        }
        if total_flow.iter().any(|x| !x.hi.is_finite() || !x.lo.is_finite()) {
            return Err(invalid("Rayleigh total volume flow exceeds finite numerical range"));
        }
        if cancelled() { return Err(HelmholtzError::Cancelled); }
        Ok(PreparedCoherentRadiation { surface: self, flow, total_flow,
            k: omega / medium.sound_speed, pressure_scale })
    }

    fn impedance_image(&self, omega: f64, medium: Medium, power_consistent: bool)
        -> Result<Vec<C64>, HelmholtzError> {
        let ppw = self.equivalent_disc_ppw(omega, medium)?;
        if ppw < MIN_PANELS_PER_WAVELENGTH {
            return Err(HelmholtzError::TooCoarse { panels_per_wavelength: ppw });
        }
        let k = omega / medium.sound_speed;
        let prefactor = C64::new(0.0, -2.0 * omega * medium.density);
        if !prefactor.im.is_finite() {
            return Err(invalid("Rayleigh pressure scale exceeds finite numerical range"));
        }
        let resistance_scale = if power_consistent {
            let value = (-prefactor.im * (0.25 / PI)) * k;
            if !value.is_finite() || value <= 0.0 {
                return Err(invalid("Rayleigh regular resistance exceeds finite numerical range"));
            }
            value
        } else { 0.0 };
        let low_frequency_resistance = if power_consistent && self.count > 1
            && k * self.maximum_distance_m <= 0.5 {
            Some(resistance::matrix(self, k, resistance_scale)?)
        } else { None };
        let mut pressure = reserved(self.count * self.modes, "complete Rayleigh pressure fields")?;
        pressure.resize(self.count * self.modes, C64::ZERO);
        for i in 0..self.count {
            let phase = k * self.radii[i];
            // rho*c*(1-exp(i*k*a))/A is the existing disc integral after
            // multiplying by -i*omega*rho/(2*pi) and dividing by source area.
            // 1-cos(x)=2*sin(x/2)^2 retains weak low-frequency resistance.
            let scale = medium.density * medium.sound_speed / self.areas[i];
            let mut self_kernel = C64::new(2.0 * det::sin(0.5 * phase).powi(2), -det::sin(phase)).scale(scale);
            if power_consistent { self_kernel.re = resistance_scale; }
            if low_frequency_resistance.is_some() { self_kernel.re = 0.0; }
            for mode in 0..self.modes {
                pressure[i * self.modes + mode] = self_kernel.scale(self.weighted_shape[i * self.modes + mode].to_f64());
            }
        }
        let mut pair = 0;
        for i in 0..self.count {
            for j in i + 1..self.count {
                let r = self.distances[pair];
                pair += 1;
                if !(k * r).is_finite() {
                    return Err(invalid("Rayleigh travel phase exceeds finite numerical range"));
                }
                let mut kernel = green(k, r) * prefactor;
                // The shared off-diagonal kernel already has real part
                // rho*omega*k/(2*pi)*sinc(k*r). Reuse its sine evaluation.
                // Retain the regular limit if only the phase underflows.
                if power_consistent && k * r == 0.0 { kernel.re = resistance_scale; }
                if low_frequency_resistance.is_some() { kernel.re = 0.0; }
                for mode in 0..self.modes {
                    let a = i * self.modes + mode;
                    let b = j * self.modes + mode;
                    pressure[a] = pressure[a] + kernel.scale(self.weighted_shape[b].to_f64());
                    pressure[b] = pressure[b] + kernel.scale(self.weighted_shape[a].to_f64());
                }
            }
        }
        let mut result = reserved(self.modes * self.modes, "complete Rayleigh modal impedance")?;
        result.resize(self.modes * self.modes, C64::ZERO);
        for row in 0..self.modes {
            for column in row..self.modes {
                let mut value = C64::ZERO;
                for cell in 0..self.count {
                    value = value + pressure[cell * self.modes + column].scale(self.weighted_shape[cell * self.modes + row].to_f64());
                }
                if let Some(resistance) = &low_frequency_resistance {
                    value.re = resistance[row * self.modes + column];
                }
                if !value.re.is_finite() || !value.im.is_finite() {
                    return Err(invalid("Rayleigh modal impedance exceeds finite numerical range"));
                }
                result[row * self.modes + column] = value;
                result[column * self.modes + row] = value;
            }
        }
        Ok(result)
    }
}

/// Exact static pressure/volume-flow inertance of a uniform rectangular opening.
/// The same double Rayleigh integral supplies the independent shape check for
/// shared surface quadrature. The actual rectangle shape is retained without
/// an equivalent circular end correction. The stable slender-aspect arithmetic is unchanged
/// from the original compact-mouth owner.
pub fn rectangular_static_inertance(size: [f64; 2], density: f64) -> Result<f64, HelmholtzError> {
    if size.iter().any(|value| !value.is_finite() || *value <= 0.0)
        || !density.is_finite() || density <= 0.0 {
        return Err(invalid("rectangular inertia requires positive finite geometry and density"));
    }
    let length = size[0].max(size[1]);
    let t = size[0].min(size[1]) / length;
    if !t.is_finite() || t <= 0.0 {
        return Err(invalid("rectangular aspect ratio is unrepresentable"));
    }
    let root = fs_math::det::sqrt(1.0 + t * t);
    // asinh(t)/t near zero; first omitted term <3.1e-18 on this interval.
    let asinh_ratio = if t < 0.01 {
        let t2 = t * t;
        1.0 + t2 * (-1.0 / 6.0 + t2 * (3.0 / 40.0 - t2 * (5.0 / 112.0)))
    } else {
        fs_math::det::ln(t + root) / t
    };
    let inverse_asinh = fs_math::det::ln(1.0 + root) - fs_math::det::ln(t);
    let shape = asinh_ratio + inverse_asinh + (t - (root * root + root + 1.0) / (root + 1.0)) / 3.0;
    let inertia = (density / length) * (shape / core::f64::consts::PI);
    if !inertia.is_finite() || inertia <= 0.0 {
        return Err(invalid("rectangular inertia is unrepresentable"));
    }
    Ok(inertia)
}

#[cfg(test)]
#[path = "rayleigh_modal_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "rayleigh_modal_power_tests.rs"]
mod power_tests;
