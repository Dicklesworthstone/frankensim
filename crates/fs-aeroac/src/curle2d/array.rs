//! Coherent superposition with the full near-field pressure gradient.
//!
//! Each source is a compact line dipole, but their separation need not be
//! acoustically compact. Apply the outgoing Green function at each actual
//! position BEFORE summing. Collapsing the forces to one point erases path
//! phase; adding source intensities erases interference. Statistical observations
//! accept explicit independent coherence factors, or independent source vectors;
//! unknown relative phase is never silently replaced with full coherence.
//!
//! Peak phasors use exp(-i omega t) and ONE caller-defined time origin. From
//! linear momentum, u = grad(p)/(i rho omega), so complex intensity is
//! p conj(u)/2. In particular, |p|²/(2 rho c) is not the intensity of a
//! standing wave or a dipole near field. The standard intensity definition is
//! also documented at <https://doc.comsol.com/6.3/doc/com.comsol.help.aco/aco_ug_pressure.05.140.html>
//! (with the opposite harmonic convention there).
//!
//! This is homogeneous, stationary, lossless free-space propagation of
//! prescribed sources, without scattering or source-to-source fluid loading.
//! It does not extend the original Curle source model's calibration authority.

use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, dd::Dd, det};

use crate::bessel::{hankel0_outgoing, hankel1_outgoing};
use crate::{AeroacError, SCOPE_STATEMENT};

/// Total field at one observer, retaining interference and reactive motion.
#[derive(Clone, Debug)]
pub struct CoherentDipoleField {
    /// Complex peak pressure [Pa], exp(-i omega t).
    pub pressure_pa: C64,
    /// Complex peak particle velocity [m/s] in the two spatial directions.
    pub particle_velocity_m_s: [C64; 2],
    /// Time-averaged energy flux [W/m²], Re(p conj(u))/2. Components may be
    /// negative; no plane-wave approximation or positivity clamp is applied.
    pub active_intensity_w_m2: [f64; 2],
    /// Im(p conj(u))/2 [W/m²], in the declared negative-time convention.
    pub reactive_intensity_w_m2: [f64; 2],
    /// Existing two-dimensional source/absolute-level authority boundary.
    pub scope: &'static str,
}

/// Work limits for a statistical observation, checked before propagation.
#[derive(Clone, Copy, Debug)]
pub struct CoherenceBudget {
    /// Maximum number of independent coherence factors, including empty ones.
    pub max_modes: usize,
    /// Maximum product of factor count and prepared source count.
    pub max_source_mode_pairs: usize,
}

/// Second moments at one frequency, averaged over time and source statistics.
///
/// These are spectral-line powers, not power spectral densities per hertz.
/// The input factors use peak-force units, so every quadratic quantity includes
/// the harmonic time-average factor 1/2. No random realization, probability law,
/// confidence interval or source-coherence estimate is inferred.
#[derive(Clone, Debug)]
pub struct DipoleFieldStatistics {
    /// E[|p|²]/2 [Pa²]. Its square root is pressure RMS, not energy flux.
    pub mean_square_pressure_pa2: f64,
    /// E[u_i conj(u_j)]/2 [m²/s²], a Hermitian spectral covariance.
    /// The real diagonal entries are mean-square component velocities.
    pub velocity_covariance_m2_s2: [[C64; 2]; 2],
    /// Re(E[p conj(u)])/2 [W/m²]. Signed components retain interference.
    pub active_intensity_w_m2: [f64; 2],
    /// Im(E[p conj(u)])/2 [W/m²], in exp(-i omega t) convention.
    pub reactive_intensity_w_m2: [f64; 2],
    /// Existing two-dimensional source/absolute-level authority boundary.
    pub scope: &'static str,
}

#[derive(Clone, Debug)]
struct Transfer {
    pressure: [C64; 2],
    velocity: [[C64; 2]; 2],
}

/// Reusable source-to-observer operator at one frequency and material state.
///
/// Preparation evaluates the Hankel kernels and stores six complex coefficients
/// per source. Applying new force phasors is allocation-free and O(source count),
/// with no special-function evaluations. Input order fixes reduction order;
/// double-double accumulation preserves small residuals between strong sources.
#[derive(Clone, Debug)]
pub struct PreparedDipoleObserver {
    acoustic: FluidAcousticState,
    frequency_hz: f64,
    transfers: Vec<Transfer>,
}

impl PreparedDipoleObserver {
    /// Prepare for source positions [m] in their permanent force-array order.
    /// `max_sources` is the caller's allocation/work budget. An empty source
    /// list represents silence. The acoustic snapshot must be that of the
    /// original source medium; it supplies both propagation and fluid inertia.
    ///
    /// # Errors
    /// Refuses a budget overrun before allocation, nonpositive/nonfinite
    /// frequency, invalid geometry, a source at the observer, or coefficients
    /// outside the representable numerical range.
    pub fn new(
        acoustic: FluidAcousticState,
        frequency_hz: f64,
        observer_m: [f64; 2],
        source_positions_m: &[[f64; 2]],
        max_sources: usize,
    ) -> Result<Self, AeroacError> {
        if source_positions_m.len() > max_sources {
            return Err(invalid(
                "dipole source count exceeds the preparation budget",
            ));
        }
        let omega = core::f64::consts::TAU * frequency_hz;
        let k = omega / acoustic.sound_speed_m_s();
        let rho_omega = acoustic.density_kg_m3() * omega;
        if [frequency_hz, omega, k, rho_omega]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(invalid(
                "dipole frequency and material wave scales must be positive and finite",
            ));
        }
        if observer_m.iter().any(|v| !v.is_finite()) {
            return Err(nonfinite("dipole observer position"));
        }
        let mut transfers = Vec::new();
        transfers
            .try_reserve_exact(source_positions_m.len())
            .map_err(|_| invalid("dipole transfer allocation exceeds available memory"))?;
        for source in source_positions_m {
            if source.iter().any(|v| !v.is_finite()) {
                return Err(nonfinite("dipole source position"));
            }
            let dx = observer_m[0] - source[0];
            let dy = observer_m[1] - source[1];
            let r = det::hypot(dx, dy);
            let kr = k * r;
            if r == 0.0 {
                return Err(invalid("observer coincides with a dipole source"));
            }
            if !r.is_finite() || !kr.is_finite() || kr <= 0.0 {
                return Err(nonfinite("dipole source-observer wave scale"));
            }
            let direction = [dx / r, dy / r];
            let h1 = hankel1_outgoing(kr);
            let ik4 = C64::new(0.0, 0.25 * k);
            let radial = ik4 * h1;
            // dH1(z)/dz = H0(z) - H1(z)/z. Differentiate both
            // the radial Green kernel AND the dipole direction cosine.
            let radial_derivative = ik4 * (hankel0_outgoing(kr).scale(k) - h1.scale(1.0 / r));
            let tangential_derivative = radial.scale(1.0 / r);
            let mut transfer = Transfer {
                pressure: direction.map(|d| radial.scale(d)),
                velocity: [[C64::ZERO; 2]; 2],
            };
            for (axis, row) in transfer.velocity.iter_mut().enumerate() {
                for (force_axis, coefficient) in row.iter_mut().enumerate() {
                    let nn = direction[axis] * direction[force_axis];
                    let identity = if axis == force_axis { 1.0 } else { 0.0 };
                    let gradient =
                        radial_derivative.scale(nn) + tangential_derivative.scale(identity - nn);
                    *coefficient = C64::new(gradient.im / rho_omega, -gradient.re / rho_omega);
                }
            }
            if transfer
                .pressure
                .iter()
                .chain(transfer.velocity.iter().flatten())
                .any(|v| !finite(*v))
            {
                return Err(nonfinite("dipole field transfer"));
            }
            transfers.push(transfer);
        }
        Ok(Self {
            acoustic,
            frequency_hz,
            transfers,
        })
    }

    /// Number of source force vectors required on every application.
    #[must_use]
    pub fn source_count(&self) -> usize {
        self.transfers.len()
    }

    /// The shared physical frequency [Hz]. Distinct frequencies are separate
    /// observations; they cannot be coherently combined as one spectral line.
    #[must_use]
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }

    /// Exact homogeneous material snapshot used to prepare every source path.
    #[must_use]
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }

    /// Apply complex peak force per span [N/m] to each source. Phases must
    /// reference the SAME physical clock and frequency. For example, delaying
    /// one history by tau multiplies its force by exp(+i omega tau).
    /// No relative phases are synthesized and no pressure normalization occurs.
    ///
    /// # Errors
    /// Refuses mismatched source counts and nonfinite or overflowing fields.
    pub fn evaluate(&self, forces_n_m: &[[C64; 2]]) -> Result<CoherentDipoleField, AeroacError> {
        let [pressure, velocity_x, velocity_y] = self.amplitudes(forces_n_m)?;
        let velocity = [velocity_x, velocity_y];
        // Form quadratic observables only AFTER coherent linear superposition.
        let intensity = velocity.map(|v| pressure.scale(0.5) * v.conj());
        if intensity.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("coherent dipole field or intensity"));
        }
        Ok(CoherentDipoleField {
            pressure_pa: pressure,
            particle_velocity_m_s: velocity,
            active_intensity_w_m2: intensity.map(|v| v.re),
            reactive_intensity_w_m2: intensity.map(|v| v.im),
            scope: SCOPE_STATEMENT,
        })
    }

    /// Propagate an explicit factorization of source coherence at this frequency.
    ///
    /// For factor `b_r[j][axis]` [N/m], the source force is
    /// `F[j][axis] = sum_r b_r[j][axis] xi_r`, with
    /// `E[xi_r conj(xi_s)] = delta_rs`. Thus the source second-moment matrix is
    /// `B B*`, positive semidefinite by construction. Each factor retains the
    /// relative phases of ALL sources and force directions on the same clock.
    /// Propagate and coherently sum each factor first, then add its quadratic
    /// observables to those of the independent factors. One factor reproduces a
    /// deterministic coherent field's power; zero factors represent silence.
    ///
    /// Factors may come from a caller's measured or modeled cross-spectrum, but
    /// this method does not estimate them or invent their correlation. To use a
    /// one-sided force PSD over a frequency bin of width df, supply factors of
    /// `2 df S_F` in (N/m)²; the factor two converts to this peak convention.
    /// Discrete tones already supplied as peak phasors need no bandwidth factor.
    /// A sum of bin powers assumes a stationary, spectrally resolved observation;
    /// finite-time cross-frequency interference is not modeled here.
    ///
    /// Application is allocation-free, O(source count × factor count), reusing
    /// the cached spatial transfers without forming a dense covariance matrix.
    /// Fixed input order and double-double reductions apply within and between
    /// factors. For independent sources, [`Self::evaluate_independent_sources`]
    /// avoids the quadratic storage/work of one sparse factor per source.
    ///
    /// # Errors
    /// Refuses either work-budget overrun (including usize product overflow),
    /// a factor with the wrong source count, and nonfinite or overflowing values.
    pub fn evaluate_coherence_modes(
        &self,
        force_modes_n_m: &[&[[C64; 2]]],
        budget: CoherenceBudget,
    ) -> Result<DipoleFieldStatistics, AeroacError> {
        let pairs = force_modes_n_m.len().checked_mul(self.source_count());
        if force_modes_n_m.len() > budget.max_modes
            || pairs.is_none_or(|count| count > budget.max_source_mode_pairs)
        {
            return Err(invalid("dipole coherence factors exceed the work budget"));
        }
        if force_modes_n_m
            .iter()
            .any(|mode| mode.len() != self.source_count())
        {
            return Err(invalid(
                "each coherence factor requires one vector per dipole source",
            ));
        }
        let mut statistics = StatisticsSum::default();
        for mode in force_modes_n_m {
            statistics.add(self.amplitudes(mode)?)?;
        }
        statistics.finish()
    }

    /// Propagate sources with mutually uncorrelated unit-power amplitudes.
    ///
    /// Each vector [N/m] fixes that source's peak amplitude, polarization and
    /// relative x/y phase. The two components of one source remain coherent;
    /// distinct sources have zero cross second moment. This is an explicit
    /// statistical assumption, not the treatment of merely unknown phases.
    ///
    /// This allocation-free O(source count) path applies each cached transfer
    /// once and adds source second moments. It is equivalent to one factor per
    /// source with zeros at every other source, without constructing those
    /// O(source count²) arrays. The preparation's source budget bounds the work.
    ///
    /// # Errors
    /// Refuses a mismatched source count, nonfinite forces or moment overflow.
    pub fn evaluate_independent_sources(
        &self,
        forces_n_m: &[[C64; 2]],
    ) -> Result<DipoleFieldStatistics, AeroacError> {
        if forces_n_m.len() != self.source_count() {
            return Err(invalid(
                "one force vector is required per prepared dipole source",
            ));
        }
        let mut statistics = StatisticsSum::default();
        for (transfer, force) in self.transfers.iter().zip(forces_n_m) {
            let mut sums = [ComplexSum::default(); 3];
            transfer.accumulate(*force, &mut sums)?;
            statistics.add(sums.map(ComplexSum::value))?;
        }
        statistics.finish()
    }

    fn amplitudes(&self, forces_n_m: &[[C64; 2]]) -> Result<[C64; 3], AeroacError> {
        if forces_n_m.len() != self.source_count() {
            return Err(invalid(
                "one force vector is required per prepared dipole source",
            ));
        }
        let mut sums = [ComplexSum::default(); 3];
        for (transfer, force) in self.transfers.iter().zip(forces_n_m) {
            transfer.accumulate(*force, &mut sums)?;
        }
        let values = sums.map(ComplexSum::value);
        if values.iter().any(|value| !finite(*value)) {
            return Err(nonfinite("coherent dipole field"));
        }
        Ok(values)
    }
}

impl Transfer {
    fn accumulate(&self, force: [C64; 2], sums: &mut [ComplexSum; 3]) -> Result<(), AeroacError> {
        for (axis, value) in force.into_iter().enumerate() {
            if !finite(value) {
                return Err(nonfinite("dipole source force"));
            }
            sums[0].add(self.pressure[axis] * value)?;
            for (velocity_axis, row) in self.velocity.iter().enumerate() {
                sums[velocity_axis + 1].add(row[axis] * value)?;
            }
        }
        Ok(())
    }
}

#[derive(Default)]
struct StatisticsSum {
    // Diagonal of the pressure/velocity covariance, always real.
    mean_squares: [Dd; 3],
    // Its upper triangle: p conj(ux), p conj(uy), ux conj(uy).
    cross: [ComplexSum; 3],
}

impl StatisticsSum {
    fn add(&mut self, field: [C64; 3]) -> Result<(), AeroacError> {
        for (sum, value) in self.mean_squares.iter_mut().zip(field) {
            // Scale before multiplication so a representable half-square is
            // not refused merely because the unhalved square would overflow.
            let square = (0.5 * value.re) * value.re + (0.5 * value.im) * value.im;
            if !square.is_finite() {
                return Err(nonfinite("dipole field second moment"));
            }
            *sum = *sum + Dd::from_f64(square);
        }
        self.cross[0].add(field[0].scale(0.5) * field[1].conj())?;
        self.cross[1].add(field[0].scale(0.5) * field[2].conj())?;
        self.cross[2].add(field[1].scale(0.5) * field[2].conj())?;
        Ok(())
    }

    fn finish(self) -> Result<DipoleFieldStatistics, AeroacError> {
        let squares = self.mean_squares.map(Dd::to_f64);
        let cross = self.cross.map(ComplexSum::value);
        if squares.iter().any(|value| !value.is_finite())
            || cross.iter().any(|value| !finite(*value))
        {
            return Err(nonfinite("accumulated dipole field statistics"));
        }
        Ok(DipoleFieldStatistics {
            mean_square_pressure_pa2: squares[0],
            velocity_covariance_m2_s2: [
                [C64::from_re(squares[1]), cross[2]],
                [cross[2].conj(), C64::from_re(squares[2])],
            ],
            active_intensity_w_m2: [cross[0].re, cross[1].re],
            reactive_intensity_w_m2: [cross[0].im, cross[1].im],
            scope: SCOPE_STATEMENT,
        })
    }
}

#[derive(Clone, Copy, Default)]
struct ComplexSum {
    re: Dd,
    im: Dd,
}

impl ComplexSum {
    fn add(&mut self, value: C64) -> Result<(), AeroacError> {
        if !finite(value) {
            return Err(nonfinite("dipole source contribution"));
        }
        self.re = self.re + Dd::from_f64(value.re);
        self.im = self.im + Dd::from_f64(value.im);
        Ok(())
    }
    fn value(self) -> C64 {
        C64::new(self.re.to_f64(), self.im.to_f64())
    }
}

fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}
fn invalid(what: &'static str) -> AeroacError {
    AeroacError::InvalidParameter { what }
}
fn nonfinite(what: &'static str) -> AeroacError {
    AeroacError::NonFinite { what }
}
