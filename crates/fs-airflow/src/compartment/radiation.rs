//! Reciprocal gray, nonscattering radiation through the original two gas layers.
//!
//! The opaque boundary geometry and emissivities remain owned by
//! `fs-conduction`. Authored reciprocal chords subdivide its admitted view
//! factors and carry ordered path lengths through the lower/upper layers.
//! No geometry, gas opacity, soot yield or spectral average is inferred.
//!
//! On each isothermal segment the formal transfer solution is
//! `I_out = exp(-kappa L) I_in + (1-exp(-kappa L)) sigma T^4`, with diffuse
//! surface reflection `J = epsilon sigma T^4 + (1-epsilon) H`.
//! Here `I` is pi times bolometric radiance, so its blackbody value is
//! `sigma T^4`; the reciprocal area weight already contains the angular measure.
//! References: NIST TN 1889v1 (2015), sections 5.1 and 5.1.4,
//! <https://doi.org/10.6028/NIST.TN.1889v1>, equations 5.1 and 5.19.
//! The reciprocal exchange coefficients below follow by expanding that
//! piecewise-constant formal solution in both directions along each chord.
//!
//! This is a finite optical-path quadrature and gray-gas model, not a spectral
//! gas database, RADCAL implementation, scattering model, or geometry generator.
//! The caller must rebuild paths when the layer interface or geometry changes.
//! View-factor closure and conservation do not establish geometric or physical
//! accuracy. Direct flame radiation remains a separate, explicitly booked source.

use fs_conduction::radiation::{STEFAN_BOLTZMANN_W_M2_K4, SurfaceEmissivity, ViewFactorMatrix};
use fs_gas::multicomponent::compartment::{Layer, LayerRates, TwoZone, ZoneError};
use fs_la::prepared_lu::{LuWorkspace, LuWorkspaceError};
use fs_qty::{Area, Length, Power, Qty, Temperature};

/// Source-explicit Rayleigh soot absorption and its Planck-weighted gray mean.
pub mod soot;

/// Wavelength-resolved band/soot transfer through these same reciprocal paths.
pub mod spectral;

/// Source-bounded molecular Planck means from actual gas partial pressures.
pub mod molecular;

/// Participating radiation shared by finite walls and receding reactive slabs.
pub mod reactive;

/// Participating radiation coupled to original finite-wall endpoints.
pub mod walls;

/// Maximum surface axis for the bounded dense reflection solve.
pub const MAX_RADIATION_SURFACES: usize = 64;
/// Maximum authored reciprocal chords per enclosure.
pub const MAX_OPTICAL_PATHS: usize = 16_384;
/// Gray absorption coefficient [1/m], excluding scattering.
pub type AbsorptionCoefficient = Qty<-1, 0, 0, 0, 0, 0>;

/// Explicit absorption source at the actual evaluated layer state.
/// Zero declares a transparent layer; it never means missing material data.
#[derive(Clone, Debug, PartialEq)]
pub struct LayerAbsorption {
    coefficient: AbsorptionCoefficient,
    source: String,
}
impl LayerAbsorption {
    /// Bind a finite nonnegative coefficient to its stated source/model.
    /// The owner must resolve its temperature/composition validity before calling.
    pub fn new(
        coefficient: AbsorptionCoefficient,
        source: impl Into<String>,
    ) -> Result<Self, RadiationError> {
        let source = source.into();
        if !nonnegative(coefficient.value()) || source.trim().is_empty() {
            return Err(RadiationError::Invalid(
                "finite nonnegative sourced gas absorption required",
            ));
        }
        Ok(Self {
            coefficient,
            source,
        })
    }
    /// Original supplied absorption coefficient.
    pub const fn coefficient(&self) -> AbsorptionCoefficient {
        self.coefficient
    }
    /// Original source/model identity, without an inferred uncertainty claim.
    pub fn source(&self) -> &str {
        &self.source
    }
}

/// One ordered, uniform portion of a reciprocal optical chord.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OpticalSegment {
    /// Original gas-layer axis.
    pub layer: Layer,
    /// Nonnegative distance traversed inside this layer.
    pub length: Length,
}

/// A pair of reversed diffuse optical directions sharing one exchange area.
/// For different endpoints, the sum of these weights is `A_i F_ij = A_j F_ji`.
/// A self-view chord contributes TWICE its weight to `A_i F_ii`, because both
/// endpoint directions belong to that same patch. No factor is renormalized.
#[derive(Clone, Debug, PartialEq)]
pub struct OpticalPath {
    /// Ordered boundary endpoints, in the original view-factor axis.
    pub surfaces: [usize; 2],
    /// Reciprocal area weight [m²], strictly positive.
    pub exchange_area: Area,
    /// Traversal from the first surface to the second, at most two segments.
    /// Empty means vacuum/transparent direct exchange. Two equal layer IDs
    /// may explicitly subdivide one isothermal optical segment.
    pub segments: Vec<OpticalSegment>,
}

/// Original admitted boundary geometry plus its explicitly supplied chord model.
#[derive(Clone, Debug, PartialEq)]
pub struct ParticipatingEnclosure {
    geometry: ViewFactorMatrix,
    paths: Vec<OpticalPath>,
    path_source: String,
    maximum_geometry_relative_residual: f64,
}

/// Explicit algebraic acceptance gates, separate from physical model error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RadiationSolve {
    /// Largest diffuse-surface equation defect [W/m²].
    pub surface_residual_w_m2: f64,
    /// Largest absolute sum of all gas and boundary powers [W].
    pub energy_residual_w: f64,
}
impl Default for RadiationSolve {
    fn default() -> Self {
        Self {
            surface_residual_w_m2: 1e-8,
            energy_residual_w: 1e-7,
        }
    }
}

/// Reusable numerical storage. No physical state is mutated by a radiation solve.
pub struct RadiationWorkspace {
    n: usize,
    exchange: Vec<f64>,
    matrix: Vec<f64>,
    rhs: Vec<f64>,
    emission: Vec<f64>,
    radiosity: Vec<f64>,
    powers: Vec<f64>,
    lu: LuWorkspace,
}

/// Signed heat exchange at one frozen set of gas and boundary temperatures.
#[derive(Clone, Debug, PartialEq)]
pub struct RadiationReport {
    /// Watts INTO boundary patches, in the original view-factor order.
    pub surface_heat_power: Vec<Power>,
    /// Watts INTO the original `[lower, upper]` gas layers.
    pub gas_heat_power: [Power; 2],
    /// Outgoing surface radiosities, including reflected irradiation [W/m²].
    pub surface_radiosity_w_m2: Vec<f64>,
    /// Original surface source objects used in this evaluation.
    pub surface_emissivities: Vec<SurfaceEmissivity>,
    /// Actual supplied gas optical coefficients and sources.
    pub absorption: [LayerAbsorption; 2],
    /// Largest independent surface reflection-equation defect [W/m²].
    pub maximum_surface_residual_w_m2: f64,
    /// Signed sum of boundary and gas powers; retained, never redistributed [W].
    pub energy_residual_w: f64,
}
impl RadiationReport {
    /// Add the identical gas exchanges to the existing source owner, atomically.
    /// The caller must also apply `surface_heat_power` to the corresponding
    /// walls exactly once. These are internal transfers, not fire heat release.
    pub fn add_gas_heat(&self, rates: &mut [LayerRates; 2]) -> Result<(), RadiationError> {
        let mut candidate = *rates;
        for (rate, &heat) in candidate.iter_mut().zip(&self.gas_heat_power) {
            rate.add_heat(heat)?;
        }
        *rates = candidate;
        Ok(())
    }
}

/// Input, source-domain, numerical or cancellation refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum RadiationError {
    /// Invalid source, geometry, extent, temperature or arithmetic.
    Invalid(&'static str),
    /// A participating molecular gas lies outside the original fit domain.
    MolecularTemperature {
        /// Actual blackbody-weighting/local gas temperature [K].
        temperature_k: f64,
    },
    /// Original dense solver refusal.
    Linear(LuWorkspaceError),
    /// Original gas source addition refusal.
    Zone(ZoneError),
    /// The actual equations exceed the supplied numerical budget.
    Residual {
        /// Largest surface equation defect [W/m²].
        surface_w_m2: f64,
        /// Net enclosure energy defect [W].
        energy_w: f64,
    },
    /// The declared spectral work budget cannot resolve the requested outputs.
    SpectralBudget {
        /// Actual wavelength evaluations already performed.
        evaluations: usize,
        /// Caller-declared maximum number of wavelength evaluations.
        maximum_evaluations: usize,
    },
    /// The selected wavelength window omits too much blackbody emission.
    SpectralTruncation {
        /// Largest omitted fraction among the actual surface/layer spectra.
        maximum_fraction: f64,
        /// Caller-declared maximum omitted blackbody fraction.
        allowed_fraction: f64,
    },
    /// No physical state or report was published.
    Cancelled,
}
impl From<LuWorkspaceError> for RadiationError {
    fn from(e: LuWorkspaceError) -> Self {
        Self::Linear(e)
    }
}
impl From<ZoneError> for RadiationError {
    fn from(e: ZoneError) -> Self {
        Self::Zone(e)
    }
}
impl core::fmt::Display for RadiationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "participating compartment radiation: {self:?}")
    }
}
impl std::error::Error for RadiationError {}
fn nonnegative(x: f64) -> bool {
    x.is_finite() && x >= 0.0
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn finite(x: f64) -> Result<f64, RadiationError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(RadiationError::Invalid("non-finite radiation arithmetic"))
    }
}
fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), RadiationError> {
    if cancelled() {
        Err(RadiationError::Cancelled)
    } else {
        Ok(())
    }
}
fn emission_difference(t: f64, reference: f64) -> Result<f64, RadiationError> {
    finite(
        STEFAN_BOLTZMANN_W_M2_K4
            * (t - reference)
            * (t + reference)
            * (t * t + reference * reference),
    )
}

impl ParticipatingEnclosure {
    /// Bind chords to the original admitted opaque view-factor geometry.
    ///
    /// Each ordered entry must match within `geometry_relative_tolerance`
    /// times that source surface's area. Row closure is checked independently.
    /// This allowance is numerical; it is not a ray-quadrature error estimate.
    /// Tolerances larger than 1e-8 refuse rather than hiding missing ray power.
    pub fn new(
        geometry: ViewFactorMatrix,
        paths: Vec<OpticalPath>,
        path_source: impl Into<String>,
        geometry_relative_tolerance: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RadiationError> {
        checkpoint(&mut cancelled)?;
        let n = geometry.len();
        let path_source = path_source.into();
        if !(2..=MAX_RADIATION_SURFACES).contains(&n)
            || paths.is_empty()
            || paths.len() > MAX_OPTICAL_PATHS
            || path_source.trim().is_empty()
            || !nonnegative(geometry_relative_tolerance)
            || geometry_relative_tolerance > 1e-8
        {
            return Err(RadiationError::Invalid(
                "bounded sourced optical paths and geometry tolerance required",
            ));
        }
        let mut coverage = vec![0.0; n * n];
        for path in &paths {
            checkpoint(&mut cancelled)?;
            let [a, b] = path.surfaces;
            if a >= n
                || b >= n
                || !positive(path.exchange_area.value())
                || path.segments.len() > 2
                || path.segments.iter().any(|s| !nonnegative(s.length.value()))
            {
                return Err(RadiationError::Invalid(
                    "invalid optical chord endpoints, weight or segments",
                ));
            }
            coverage[a * n + b] = finite(coverage[a * n + b] + path.exchange_area.value())?;
            coverage[b * n + a] = finite(coverage[b * n + a] + path.exchange_area.value())?;
        }
        let mut maximum = 0.0_f64;
        for i in 0..n {
            checkpoint(&mut cancelled)?;
            let area = geometry.areas_m2()[i];
            let mut row = 0.0;
            for j in 0..n {
                row = finite(row + coverage[i * n + j])?;
                let expected = finite(area * geometry.factors()[i][j])?;
                maximum = maximum.max((coverage[i * n + j] - expected).abs() / area);
            }
            maximum = maximum.max((row - area).abs() / area);
        }
        if maximum > geometry_relative_tolerance {
            return Err(RadiationError::Invalid(
                "optical paths do not close the original reciprocal view factors",
            ));
        }
        Ok(Self {
            geometry,
            paths,
            path_source,
            maximum_geometry_relative_residual: maximum,
        })
    }

    /// Original admitted boundary areas, factors and their source identity.
    pub const fn geometry(&self) -> &ViewFactorMatrix {
        &self.geometry
    }
    /// Retained optical-path generation source; no geometric accuracy is inferred.
    pub fn path_source(&self) -> &str {
        &self.path_source
    }
    /// Original reciprocal chords, including their declared traversal order.
    pub fn paths(&self) -> &[OpticalPath] {
        &self.paths
    }
    /// Measured numerical chord/geometry closure discrepancy.
    pub const fn maximum_geometry_relative_residual(&self) -> f64 {
        self.maximum_geometry_relative_residual
    }
    /// Prepare numerical buffers once for this boundary count.
    pub fn workspace(&self) -> Result<RadiationWorkspace, RadiationError> {
        let n = self.geometry.len();
        Ok(RadiationWorkspace {
            n,
            exchange: vec![0.0; (n + 2) * (n + 2)],
            matrix: vec![0.0; n * n],
            rhs: vec![0.0; n],
            emission: vec![0.0; n + 2],
            radiosity: vec![0.0; n],
            powers: vec![0.0; n + 2],
            lu: LuWorkspace::new(n)?,
        })
    }

    /// Use the actual original layer temperatures and refuse nonzero path length
    /// through an absent layer. Geometry must describe this exact interface.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_room(
        &self,
        room: &TwoZone,
        surface_temperatures: &[Temperature],
        emissivities: &[SurfaceEmissivity],
        absorption: &[LayerAbsorption; 2],
        controls: RadiationSolve,
        work: &mut RadiationWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RadiationReport, RadiationError> {
        for path in &self.paths {
            checkpoint(&mut cancelled)?;
            for segment in &path.segments {
                if segment.length.value() > 0.0 && room.layer(segment.layer).volume().value() == 0.0
                {
                    return Err(RadiationError::Invalid(
                        "optical path crosses an absent gas layer",
                    ));
                }
            }
        }
        self.evaluate(
            surface_temperatures,
            emissivities,
            [
                room.layer(Layer::Lower).temperature(),
                room.layer(Layer::Upper).temperature(),
            ],
            absorption,
            controls,
            work,
            cancelled,
        )
    }

    /// Solve gray absorption/emission and all diffuse surface reflections.
    ///
    /// Coefficients and sources remain frozen for this evaluation. The routine
    /// publishes no gas/wall update and establishes no time-integration order.
    /// Dense and chord buffers are reused; the returned source-bearing report
    /// allocates only in proportion to the boundary count.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate(
        &self,
        surface_temperatures: &[Temperature],
        emissivities: &[SurfaceEmissivity],
        gas_temperatures: [Temperature; 2],
        absorption: &[LayerAbsorption; 2],
        controls: RadiationSolve,
        work: &mut RadiationWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RadiationReport, RadiationError> {
        self.admit_state(
            surface_temperatures,
            emissivities,
            gas_temperatures,
            controls,
            work,
            &mut cancelled,
        )?;
        let n = self.geometry.len();
        let m = n + 2;
        self.prepare_optics(
            emissivities,
            absorption.each_ref().map(|a| a.coefficient.value()),
            work,
            &mut cancelled,
        )?;
        // Subtract one common blackbody background before solving. A transparent,
        // disconnected gas layer must not obscure adjacent surface temperatures.
        let mut reference = surface_temperatures
            .iter()
            .map(|t| t.value())
            .fold(f64::INFINITY, f64::min);
        for (layer, t) in gas_temperatures.iter().enumerate() {
            if work.exchange[(n + layer) * m..(n + layer + 1) * m]
                .iter()
                .any(|&g| g > 0.0)
            {
                reference = reference.min(t.value());
            }
        }
        let reference_emission = finite(STEFAN_BOLTZMANN_W_M2_K4 * reference.powi(4))?;
        for (i, t) in surface_temperatures
            .iter()
            .chain(gas_temperatures.iter())
            .enumerate()
        {
            work.emission[i] = emission_difference(t.value(), reference)?;
        }
        let (maximum, balance) =
            self.solve_emission(emissivities, controls, work, &mut cancelled)?;
        let mut radiosity = Vec::with_capacity(n);
        for &delta in &work.radiosity {
            let value = finite(reference_emission + delta)?;
            if value < 0.0 {
                return Err(RadiationError::Invalid("negative surface radiosity"));
            }
            radiosity.push(value);
        }
        checkpoint(&mut cancelled)?;
        Ok(RadiationReport {
            surface_heat_power: work.powers[..n].iter().map(|&q| Power::new(q)).collect(),
            gas_heat_power: [Power::new(work.powers[n]), Power::new(work.powers[n + 1])],
            surface_radiosity_w_m2: radiosity,
            surface_emissivities: emissivities.to_vec(),
            absorption: (*absorption).clone(),
            maximum_surface_residual_w_m2: maximum,
            energy_residual_w: balance,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn admit_state(
        &self,
        surface_temperatures: &[Temperature],
        emissivities: &[SurfaceEmissivity],
        gas_temperatures: [Temperature; 2],
        controls: RadiationSolve,
        work: &RadiationWorkspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RadiationError> {
        checkpoint(cancelled)?;
        let n = self.geometry.len();
        if work.n != n
            || surface_temperatures.len() != n
            || emissivities.len() != n
            || !nonnegative(controls.surface_residual_w_m2)
            || !nonnegative(controls.energy_residual_w)
            || surface_temperatures
                .iter()
                .chain(gas_temperatures.iter())
                .any(|t| !positive(t.value()))
        {
            return Err(RadiationError::Invalid(
                "radiation temperatures, workspace, sources or controls",
            ));
        }
        for i in 0..n {
            if let Some((lo, hi)) = emissivities[i].temperature_validity_k()
                && !(lo..=hi).contains(&surface_temperatures[i].value())
            {
                return Err(RadiationError::Invalid(
                    "surface temperature is outside original emissivity validity",
                ));
            }
        }
        Ok(())
    }

    // Both gray and spectral radiation use this exact optical operator. The
    // retained factor can serve many emission spectra when band optics are fixed.
    fn prepare_optics(
        &self,
        emissivities: &[SurfaceEmissivity],
        absorption: [f64; 2],
        work: &mut RadiationWorkspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RadiationError> {
        let n = self.geometry.len();
        let m = n + 2;
        if absorption.iter().any(|&k| !nonnegative(k)) {
            return Err(RadiationError::Invalid("invalid spectral absorption"));
        }
        work.exchange.fill(0.0);
        for path in &self.paths {
            checkpoint(cancelled)?;
            let count = path.segments.len();
            let mut nodes = [0_usize; 4];
            let mut absorptance = [1.0; 4];
            let mut transmission = [0.0; 4];
            nodes[0] = path.surfaces[0];
            nodes[count + 1] = path.surfaces[1];
            for (k, segment) in path.segments.iter().enumerate() {
                let layer = segment.layer.index();
                let optical_depth = finite(absorption[layer] * segment.length.value())?;
                nodes[k + 1] = n + layer;
                transmission[k + 1] = (-optical_depth).exp();
                absorptance[k + 1] = -(-optical_depth).exp_m1();
            }
            // Pair exchange is area * both absorptances * all intervening
            // transmittances. Opaque endpoints have absorptance one; their
            // reflected fraction is handled by the surface radiosity equation.
            for a in 0..count + 2 {
                let mut through = 1.0;
                for b in a + 1..count + 2 {
                    let i = nodes[a];
                    let j = nodes[b];
                    if i != j {
                        let g = finite(
                            path.exchange_area.value() * absorptance[a] * through * absorptance[b],
                        )?;
                        let next = finite(work.exchange[i * m + j] + g)?;
                        work.exchange[i * m + j] = next;
                        work.exchange[j * m + i] = next;
                    }
                    through *= transmission[b];
                }
            }
        }
        work.matrix.fill(0.0);
        for i in 0..n {
            checkpoint(cancelled)?;
            let e = emissivities[i].value();
            let reflected_per_area = (1.0 - e) / self.geometry.areas_m2()[i];
            let mut diagonal = e;
            for j in 0..m {
                let coefficient = finite(reflected_per_area * work.exchange[i * m + j])?;
                diagonal = finite(diagonal + coefficient)?;
                if j < n {
                    work.matrix[i * n + j] = -coefficient;
                }
            }
            work.matrix[i * n + i] = diagonal;
        }
        work.lu.factor(&work.matrix)?;
        checkpoint(cancelled)
    }

    // `emission` contains the blackbody input minus one common background.
    // Pair powers are independent of that background at every wavelength.
    fn solve_emission(
        &self,
        emissivities: &[SurfaceEmissivity],
        controls: RadiationSolve,
        work: &mut RadiationWorkspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), RadiationError> {
        let n = self.geometry.len();
        let m = n + 2;
        for i in 0..n {
            checkpoint(cancelled)?;
            let e = emissivities[i].value();
            let reflected_per_area = (1.0 - e) / self.geometry.areas_m2()[i];
            let mut rhs = e * work.emission[i];
            for j in n..m {
                let coefficient = finite(reflected_per_area * work.exchange[i * m + j])?;
                rhs = finite(coefficient.mul_add(work.emission[j], rhs))?;
            }
            work.rhs[i] = finite(rhs)?;
        }
        work.lu
            .solve_factored_into(&work.rhs, &mut work.radiosity)?;
        checkpoint(cancelled)?;
        work.powers.fill(0.0);
        for i in 0..m {
            checkpoint(cancelled)?;
            let a = if i < n {
                work.radiosity[i]
            } else {
                work.emission[i]
            };
            for j in i + 1..m {
                let b = if j < n {
                    work.radiosity[j]
                } else {
                    work.emission[j]
                };
                let into_i = finite(work.exchange[i * m + j] * (b - a))?;
                work.powers[i] = finite(work.powers[i] + into_i)?;
                work.powers[j] = finite(work.powers[j] - into_i)?;
            }
        }
        let mut maximum = 0.0_f64;
        for i in 0..n {
            let e = emissivities[i].value();
            let residual = finite(
                e * (work.radiosity[i] - work.emission[i])
                    - (1.0 - e) * work.powers[i] / self.geometry.areas_m2()[i],
            )?;
            maximum = maximum.max(residual.abs());
        }
        let balance = finite(work.powers.iter().copied().sum())?;
        if maximum > controls.surface_residual_w_m2 || balance.abs() > controls.energy_residual_w {
            return Err(RadiationError::Residual {
                surface_w_m2: maximum,
                energy_w: balance,
            });
        }
        checkpoint(cancelled)?;
        Ok((maximum, balance))
    }
}

#[cfg(test)]
mod tests;
