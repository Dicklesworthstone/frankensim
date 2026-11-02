//! Coherent 3D fields from prescribed compact volume and force sources.
//!
//! With peak phasors exp(-i omega t), the outgoing Green function is
//! `G = exp(ikr)/(4 pi r)`. A positive volume injection Q [m³/s] and force
//! **on the fluid** F [N] give `p = -i omega rho G Q - F.grad(G)`.
//! Away from sources, momentum gives `u = grad(p)/(i omega rho)`.
//! Differentiating the dipole direction as well as its radial kernel retains
//! tangential near-field motion, including at a dipole's pressure null.
//!
//! All source positions and phases enter the linear field before forming
//! `p conj(u)/2`. Neither summing source powers nor replacing intensity by
//! `|p|²/(2 rho c)` is valid for general interference and near fields.
//! Independent or partially coherent sources require explicit statistical
//! input; an unknown phase is never silently replaced by a common phase.
//!
//! The ordinary constructor uses a homogeneous lossless fluid. The molecular
//! constructor uses the original material's relaxing bulk modulus, with
//! `k² = rho omega²/K(omega)`, in the SAME Green function and its derivatives.
//! Momentum remains `grad(p) = i omega rho u`; inserting a lossy wave number
//! into `rho c`, or filtering only pressure, would change the energy transport.
//! This molecular path omits shear viscosity, thermal conduction and wall
//! layers, as does the molecular exterior BEM owner. No extra attenuation
//! coefficient or fitted absorption filter is introduced.
//!
//! Both paths are homogeneous, stationary, free-space point-source models.
//! Source size, motion, scattering and mutual fluid loading are not resolved.
//! In particular, an integrated force from a periodic flow cell is not
//! automatically an isolated radiator. These SI transfers do not turn 2D
//! force-per-span data into 3D sources or validate a physical source model.

use super::statistics::{ComplexSum, StatisticsSum, finite, invalid, nonfinite};
use super::{CoherenceBudget, CompactFieldStatistics};
use crate::AeroacError;
use crate::material::MaterialAeroacError;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::relaxation::RelaxationState;
use fs_math::{c64::C64, det};

/// Applicability of the prescribed three-dimensional source transfer.
pub const COMPACT_FIELD_SCOPE: &str = "Prescribed compact 3D volume flow [m³/s] and force on fluid [N]; \
homogeneous stationary lossless free-space observation, with no source-size validation, \
moving sources, scattering, absorption, mutual loading or measured-SPL qualification.";

/// Applicability of propagation through the original molecular bulk response.
pub const MOLECULAR_COMPACT_FIELD_SCOPE: &str = "Prescribed compact 3D volume flow [m³/s] and force on fluid [N]; \
homogeneous stationary free-space observation with material-owned molecular bulk relaxation and inviscid momentum; \
no shear/conduction or wall losses, source-size validation, moving sources, scattering, mutual loading or measured-SPL qualification.";

/// One compact source's complex peak amplitudes at a common source time origin.
///
/// The four components are coherent within this source, including monopole–
/// dipole cross terms. Distinct physical mechanisms should only share a vector
/// when their location and relative phase are actually known. Zero is silence.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CompactSourcePhasor {
    /// Outward volume injection Q [m³/s], not volume acceleration or velocity.
    pub volume_flow_m3_s: C64,
    /// Actual force ON the fluid in Cartesian coordinates [N], never N/m.
    pub force_on_fluid_n: [C64; 3],
}

impl CompactSourcePhasor {
    /// A pure physical volume source, without an inferred force contribution.
    #[must_use]
    pub const fn monopole(volume_flow_m3_s: C64) -> Self {
        Self {
            volume_flow_m3_s,
            force_on_fluid_n: [C64::ZERO; 3],
        }
    }

    /// A pure physical force source, without an inferred volume contribution.
    #[must_use]
    pub const fn force_dipole(force_on_fluid_n: [C64; 3]) -> Self {
        Self {
            volume_flow_m3_s: C64::ZERO,
            force_on_fluid_n,
        }
    }

    fn components(self) -> [C64; 4] {
        [
            self.volume_flow_m3_s,
            self.force_on_fluid_n[0],
            self.force_on_fluid_n[1],
            self.force_on_fluid_n[2],
        ]
    }
}

/// The total harmonic field, including signed energy transport and reactive motion.
#[derive(Clone, Debug)]
pub struct CoherentCompactField {
    /// Complex peak pressure [Pa], exp(-i omega t).
    pub pressure_pa: C64,
    /// Complex peak particle velocity [m/s] along all three Cartesian axes.
    pub particle_velocity_m_s: [C64; 3],
    /// Re(p conj(u))/2 [W/m²]. No plane-wave approximation or sign clamp.
    pub active_intensity_w_m2: [f64; 3],
    /// Im(p conj(u))/2 [W/m²], in the negative-time harmonic convention.
    pub reactive_intensity_w_m2: [f64; 3],
    /// Explicit prescribed-source model and observation boundaries.
    pub scope: &'static str,
}

#[derive(Clone, Debug)]
struct Transfer {
    source_m: [f64; 3],
    delay_s: f64,
    phase_delay_s: f64,
    // Rows: pressure, ux, uy, uz. Columns: Q, Fx, Fy, Fz, with SI units above.
    field: [[C64; 4]; 4],
}

/// Cached field transfer from multiple compact 3D sources to one observer.
///
/// Source separation may span many wavelengths; only each constituent source
/// must admit a compact model. Preparation checks an explicit source budget and
/// caches sixteen complex coefficients per source. Application is O(N), has no
/// allocation or transcendental evaluation, and uses fixed input order with
/// double-double reductions to preserve residuals between strong sources.
#[derive(Clone, Debug)]
pub struct PreparedCompactSourceObserver {
    acoustic: FluidAcousticState,
    molecular: Option<RelaxationState>,
    frequency_hz: f64,
    wavenumber_per_m: C64,
    observer_m: [f64; 3],
    transfers: Vec<Transfer>,
}

impl PreparedCompactSourceObserver {
    /// Prepare all actual source paths in one material snapshot and coordinate frame.
    ///
    /// # Errors
    /// Refuses a source count over `max_sources`, nonpositive frequency, invalid
    /// coordinates, coincident source/observer positions, failed allocation,
    /// and wave scales or transfer coefficients outside finite representation.
    pub fn new(
        acoustic: FluidAcousticState,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        max_sources: usize,
    ) -> Result<Self, AeroacError> {
        let k = C64::from_re(core::f64::consts::TAU * frequency_hz / acoustic.sound_speed_m_s());
        Self::prepare(
            acoustic,
            frequency_hz,
            observer_m,
            source_positions_m,
            max_sources,
            k,
            None,
        )
    }

    /// Prepare molecular absorption AND dispersion from one original gas state.
    ///
    /// The material's causal reservoir law supplies `K(omega)` and
    /// `k = omega sqrt(rho/K)`. Density, reservoir strengths and relaxation times
    /// are retained together; no air law or extra Stokes–Kirchhoff attenuation
    /// is added. This solves relaxing compression with inviscid momentum, not
    /// the full thermoviscous equations. Classical shear/conduction losses,
    /// spatially varying states and source feedback remain outside this model.
    ///
    /// A single spectral line has a phase delay, not a broadband propagation
    /// time: use [`Self::source_phase_delay_s`] for the actual retained phase.
    /// All original source/observer/work refusals still apply, and material
    /// refusals remain distinguishable as [`MaterialAeroacError::Material`].
    pub fn new_molecular(
        source: RelaxationState,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        max_sources: usize,
    ) -> Result<Self, MaterialAeroacError> {
        let acoustic = source
            .gas()
            .acoustic_state()
            .map_err(MaterialAeroacError::Material)?;
        let k = source
            .relaxing_wavenumber(core::f64::consts::TAU * frequency_hz)
            .map_err(MaterialAeroacError::Material)?;
        Ok(Self::prepare(
            acoustic,
            frequency_hz,
            observer_m,
            source_positions_m,
            max_sources,
            k,
            Some(source),
        )?)
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare(
        acoustic: FluidAcousticState,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        max_sources: usize,
        k: C64,
        molecular: Option<RelaxationState>,
    ) -> Result<Self, AeroacError> {
        if source_positions_m.len() > max_sources {
            return Err(invalid("compact source count exceeds preparation budget"));
        }
        let omega = core::f64::consts::TAU * frequency_hz;
        let rho_omega = acoustic.density_kg_m3() * omega;
        if [frequency_hz, omega, k.re, rho_omega]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
            || !k.im.is_finite()
            || k.im < 0.0
        {
            return Err(invalid(
                "compact-source frequency and wave scales must be positive finite",
            ));
        }
        if observer_m.iter().any(|v| !v.is_finite()) {
            return Err(nonfinite("compact observer position"));
        }
        let mut transfers = Vec::new();
        transfers
            .try_reserve_exact(source_positions_m.len())
            .map_err(|_| invalid("compact transfer allocation exceeds available memory"))?;
        for source in source_positions_m {
            if source.iter().any(|v| !v.is_finite()) {
                return Err(nonfinite("compact source position"));
            }
            let dx = core::array::from_fn::<_, 3, _>(|i| observer_m[i] - source[i]);
            let r = det::hypot(det::hypot(dx[0], dx[1]), dx[2]);
            let kr = k.scale(r);
            let inv_r = 1.0 / r;
            let delay_s = r / acoustic.sound_speed_m_s();
            let phase_delay_s = if molecular.is_some() {
                r * (k.re / omega)
            } else {
                delay_s
            };
            if [r, kr.re, inv_r, delay_s, phase_delay_s]
                .iter()
                .any(|v| !v.is_finite() || *v <= 0.0)
                || !kr.im.is_finite()
            {
                return Err(invalid(
                    "compact source-observer separation is not representable",
                ));
            }
            let n = dx.map(|v| v / r);
            let (sin, cos) = det::sin_cos(kr.re);
            let decay = if k.im == 0.0 { 1.0 } else { det::exp(-kr.im) };
            let green = C64::new(cos, sin).scale(((0.25 / core::f64::consts::PI) / r) * decay);
            // a = -G'; pressure of a force dipole is a (n.F).
            let radial = green * C64::new(inv_r + k.im, -k.re);
            // a' = G (k² + 2 i k/r - 2/r²), including the derivative
            // of exp(-Im(k) r). The angular derivative remains a/r.
            let radial_derivative = green
                * C64::new(
                    k.re * k.re - k.im * k.im - 2.0 * inv_r * inv_r - 2.0 * k.im * inv_r,
                    2.0 * k.re * (k.im + inv_r),
                );
            let tangential_derivative = radial.scale(inv_r);
            let mut field = [[C64::ZERO; 4]; 4];
            field[0][0] = green * C64::new(0.0, -rho_omega);
            if radial == C64::ZERO || field[0][0] == C64::ZERO {
                return Err(invalid("compact source transfer underflows representation"));
            }
            for axis in 0..3 {
                field[0][axis + 1] = radial.scale(n[axis]);
                field[axis + 1][0] = radial.scale(n[axis]);
                for force_axis in 0..3 {
                    let nn = n[axis] * n[force_axis];
                    let identity = if axis == force_axis { 1.0 } else { 0.0 };
                    let gradient =
                        radial_derivative.scale(nn) + tangential_derivative.scale(identity - nn);
                    field[axis + 1][force_axis + 1] =
                        C64::new(gradient.im / rho_omega, -gradient.re / rho_omega);
                }
            }
            if field.iter().flatten().any(|value| !finite(*value)) {
                return Err(nonfinite("compact source field transfer"));
            }
            transfers.push(Transfer {
                source_m: *source,
                delay_s,
                phase_delay_s,
                field,
            });
        }
        Ok(Self {
            acoustic,
            molecular,
            frequency_hz,
            wavenumber_per_m: k,
            observer_m,
            transfers,
        })
    }

    /// Exact checked equilibrium material used for every path. In the molecular
    /// path its sound speed is the equilibrium reference, not the phase speed.
    #[must_use]
    pub const fn acoustic(&self) -> FluidAcousticState {
        self.acoustic
    }

    /// Original evaluated molecular reservoirs, or None for lossless propagation.
    #[must_use]
    pub const fn molecular_state(&self) -> Option<&RelaxationState> {
        self.molecular.as_ref()
    }

    /// Actual outgoing complex wave number [1/m], in exp(-i omega t) convention.
    #[must_use]
    pub const fn wavenumber_per_m(&self) -> C64 {
        self.wavenumber_per_m
    }

    /// Physical propagation model and its explicit applicability boundary.
    #[must_use]
    pub const fn scope(&self) -> &'static str {
        if self.molecular.is_some() {
            MOLECULAR_COMPACT_FIELD_SCOPE
        } else {
            COMPACT_FIELD_SCOPE
        }
    }

    /// Frequency shared by every source [Hz]; different bins are separate observations.
    #[must_use]
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }

    /// Observer position [m] in the caller's original coordinate frame.
    #[must_use]
    pub const fn observer_position_m(&self) -> [f64; 3] {
        self.observer_m
    }

    /// Number of source vectors required on every application.
    #[must_use]
    pub fn source_count(&self) -> usize {
        self.transfers.len()
    }

    /// Original source position [m], or None for an index outside the array.
    #[must_use]
    pub fn source_position_m(&self, index: usize) -> Option<[f64; 3]> {
        self.transfers.get(index).map(|transfer| transfer.source_m)
    }

    /// Equilibrium reference time r/c [s]. This is the path travel time for
    /// lossless propagation. With molecular dispersion it is only a reference;
    /// the transfer uses [`Self::source_phase_delay_s`] and is not a pure delay.
    #[must_use]
    pub fn source_delay_s(&self, index: usize) -> Option<f64> {
        self.transfers.get(index).map(|transfer| transfer.delay_s)
    }

    /// Phase delay r Re(k)/omega [s] at this frequency, already included in the
    /// outgoing transfer. It equals r/c for the lossless constructor. In a
    /// dispersive gas it is neither the group delay nor the signal-front time;
    /// broadband propagation requires the frequency-dependent transfer itself.
    #[must_use]
    pub fn source_phase_delay_s(&self, index: usize) -> Option<f64> {
        self.transfers
            .get(index)
            .map(|transfer| transfer.phase_delay_s)
    }

    /// Coherently observe source peak phasors on one physical clock.
    ///
    /// Source phasors refer to the emission clock, not already-retarded observer
    /// time. The prepared outgoing kernel applies the path phase exactly once.
    /// No spectral normalization, source correlation or extra gas density is
    /// inserted into a force dipole. Refuses wrong counts or nonfinite results.
    pub fn evaluate(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<CoherentCompactField, AeroacError> {
        let [pressure, ux, uy, uz] = self.amplitudes(sources)?;
        let velocity = [ux, uy, uz];
        let intensity = velocity.map(|v| pressure.scale(0.5) * v.conj());
        if intensity.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("compact field intensity"));
        }
        Ok(CoherentCompactField {
            pressure_pa: pressure,
            particle_velocity_m_s: velocity,
            active_intensity_w_m2: intensity.map(|v| v.re),
            reactive_intensity_w_m2: intensity.map(|v| v.im),
            scope: self.scope(),
        })
    }

    /// Observe this frequency for finite-duration, potentially beating fields.
    ///
    /// Feed lines at the same receiver and on one physical source clock to
    /// [`crate::observation::PreparedHarmonicField`]. All propagation phases are
    /// already present; no second delay or amplitude normalization is required.
    /// The prescribed compact-source model's applicability remains unchanged.
    pub fn evaluate_harmonic_line(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<crate::observation::HarmonicFieldLine<3>, AeroacError> {
        let [pressure, ux, uy, uz] = self.amplitudes(sources)?;
        Ok(crate::observation::HarmonicFieldLine {
            frequency_hz: self.frequency_hz,
            pressure_pa: pressure,
            particle_velocity_m_s: [ux, uy, uz],
        })
    }

    /// Observe explicit source coherence factors `S = B xi`, `E[xi xi*] = I`.
    ///
    /// Each mode contains one source vector per position, with physically known
    /// cross-source and volume/force phases. Each is coherently propagated first;
    /// quadratic moments then add across independent factors. This represents a
    /// positive-semidefinite source second moment without a dense source matrix
    /// or positivity repair. One mode reproduces coherent-field powers; zero
    /// modes give silence. No coherence is estimated or invented.
    ///
    /// Factors use peak Q/F units. For one-sided source PSD in a bin of width df,
    /// factor `2 df S`, respecting the different Q and F units in each row.
    /// These outputs are spectral-line powers, not PSD/Hz. Adding distinct-bin
    /// powers assumes a stationary spectrally resolved observation, not general
    /// finite-time cross-frequency interference. This allocation-free O(N r)
    /// path checks factor count and source-factor work before propagation.
    pub fn evaluate_coherence_modes(
        &self,
        modes: &[&[CompactSourcePhasor]],
        budget: CoherenceBudget,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        let pairs = modes.len().checked_mul(self.source_count());
        if modes.len() > budget.max_modes
            || pairs.is_none_or(|count| count > budget.max_source_mode_pairs)
        {
            return Err(invalid("compact source coherence exceeds work budget"));
        }
        if modes.iter().any(|mode| mode.len() != self.source_count()) {
            return Err(invalid(
                "each coherence factor requires one vector per compact source",
            ));
        }
        let mut statistics = StatisticsSum::default();
        for mode in modes {
            statistics.add(self.amplitudes(mode)?)?;
        }
        statistics.finish(self.scope())
    }

    /// Observe mutually independent sources in O(N), without N sparse N-vector factors.
    ///
    /// Q and all force components of each individual source remain coherent;
    /// different source vectors have zero cross second moment by the caller's
    /// explicit assumption. Unknown relative phase alone does not justify this
    /// assumption. The preparation's source budget bounds allocation-free work.
    pub fn evaluate_independent_sources(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<CompactFieldStatistics, AeroacError> {
        self.check_count(sources)?;
        let mut statistics = StatisticsSum::default();
        for (transfer, source) in self.transfers.iter().zip(sources) {
            let mut sums = [ComplexSum::default(); 4];
            transfer.accumulate(*source, &mut sums)?;
            statistics.add(sums.map(ComplexSum::value))?;
        }
        statistics.finish(self.scope())
    }

    fn check_count(&self, sources: &[CompactSourcePhasor]) -> Result<(), AeroacError> {
        if sources.len() != self.source_count() {
            return Err(invalid(
                "one vector is required per prepared compact source",
            ));
        }
        Ok(())
    }

    fn amplitudes(&self, sources: &[CompactSourcePhasor]) -> Result<[C64; 4], AeroacError> {
        self.check_count(sources)?;
        let mut sums = [ComplexSum::default(); 4];
        for (transfer, source) in self.transfers.iter().zip(sources) {
            transfer.accumulate(*source, &mut sums)?;
        }
        let values = sums.map(ComplexSum::value);
        if values.iter().any(|value| !finite(*value)) {
            return Err(nonfinite("accumulated compact source field"));
        }
        Ok(values)
    }
}

impl Transfer {
    fn accumulate(
        &self,
        source: CompactSourcePhasor,
        sums: &mut [ComplexSum; 4],
    ) -> Result<(), AeroacError> {
        for (component, value) in source.components().into_iter().enumerate() {
            if !finite(value) {
                return Err(nonfinite("compact source amplitude"));
            }
            // Pure monopoles and dipoles need not multiply the absent mechanism.
            if value != C64::ZERO {
                for (row, sum) in self.field.iter().zip(sums.iter_mut()) {
                    sum.add(row[component] * value)?;
                }
            }
        }
        Ok(())
    }
}
