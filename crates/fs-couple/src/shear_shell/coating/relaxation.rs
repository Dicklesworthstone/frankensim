//! Passive film memory in the original equilibrium MITC shell basis.
//!
//! The measured scalar spectrum and the authored constant real Poisson ratio
//! make every film tensor proportional to E. The original MITC viscosity
//! projector therefore supplies its complete unit-E strain form, including
//! membrane/bending coupling, transverse shear and the equilibrium bubble
//! transformation. Positive Prony increments multiply this same form.
//!
//! Equal relaxation times are combined BEFORE positive factorization across
//! all regions and components. This is exact state-space aggregation for
//! relaxed or uniformly unrelaxed initial material, not diagonal modal loss.
//! A transient selects the globally passive Prony continuation explicitly;
//! fidelity to the measured law is claimed only on the retained fit band.

use super::{ResolvedShearShell, ResolvedShellCoating, apply_coatings, compile_material_shear_shell,
    invalid, owner, poll};
use crate::acoustic_realize::AcousticRealizeError;
use fs_exec::CancelGate;
use fs_material::visco::LoweredModel;
use fs_plate::shell::mitc::{MitcBudget, NodalDirectors, viscosity::project_viscosity};
use std::collections::BTreeMap;

/// Explicit source-fit accuracy and complete retained material-state limits.
/// Selection means a passive Prony continuation outside the measured band;
/// it does not widen any source claim or make a transient band-limited.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoatingRelaxationOptions {
    /// Candidate relaxation times in the original positive lowering (1..=32).
    pub terms: usize,
    /// Independent sampled relative storage-error ceiling, in (0,1).
    pub storage_relative_tolerance: f64,
    /// Independent sampled relative loss-error ceiling, in (0,1).
    pub loss_relative_tolerance: f64,
    /// Complete scalar Maxwell coordinates after equal-time aggregation.
    pub max_memory_coordinates: usize,
}

impl CoatingRelaxationOptions {
    pub(crate) fn validate(self) -> Result<(), AcousticRealizeError> {
        if !(1..=32).contains(&self.terms)
            || !(1..=8192).contains(&self.max_memory_coordinates)
            || !self.storage_relative_tolerance.is_finite()
            || !self.loss_relative_tolerance.is_finite()
            || self.storage_relative_tolerance <= 0.0 || self.storage_relative_tolerance >= 1.0
            || self.loss_relative_tolerance <= 0.0 || self.loss_relative_tolerance >= 1.0 {
            return Err(invalid("coating memory requires 1..=32 terms, explicit relative fit errors and 1..=8192 states"));
        }
        Ok(())
    }

    pub(crate) fn projection_terms(self, n: usize, films: usize) -> Option<usize> {
        // Conservative before identical-time and identical-law aggregation.
        let squares = n.checked_mul(n)?;
        let branches = films.checked_mul(self.terms)?;
        let accumulation = squares.checked_mul(branches)?;
        let factor = squares.checked_mul(n)?.checked_mul(branches)?.checked_mul(4)?;
        accumulation.checked_add(factor)
    }
}

/// One original film and its full positive strain form in the shared body basis.
#[derive(Clone, Debug)]
pub struct ProjectedShellCoating {
    source: ResolvedShellCoating,
    unit_stiffness: Vec<f64>,
}
impl ProjectedShellCoating {
    /// Original material, thickness, signed faces and equilibrium basis storage.
    pub const fn source(&self) -> &ResolvedShellCoating { &self.source }
    /// Row-major H [1/(Pa s²)] such that the film contributes E H to modal K.
    pub fn unit_stiffness_per_pa_s2(&self) -> &[f64] { &self.unit_stiffness }
}

/// Complete source-bound memory on the existing mass-normalized shell modes.
/// Each row c stores `(c.q-z)^2/2` and dissipates `(c.q-z)^2/tau`.
#[derive(Clone, Debug)]
pub struct ShellCoatingRelaxation {
    films: Vec<ProjectedShellCoating>,
    fits: Vec<LoweredModel>,
    rows: Vec<Vec<f64>>,
    times: Vec<f64>,
    modes: usize,
    band: (f64, f64),
    options: CoatingRelaxationOptions,
}
impl ShellCoatingRelaxation {
    pub const fn mode_count(&self) -> usize { self.modes }
    pub fn coupling_rows_per_s(&self) -> &[Vec<f64>] { &self.rows }
    pub fn relaxation_times_s(&self) -> &[f64] { &self.times }
    /// Intersected measured fit/substrate support, not the continuation domain.
    pub const fn source_band_rad_s(&self) -> (f64, f64) { self.band }
    /// Unique original positive lowering results, including separate sampled
    /// storage/loss errors. Their grid checks are not continuous error bounds.
    pub fn fit_certificates(&self) -> &[LoweredModel] { &self.fits }
    pub fn films(&self) -> &[ProjectedShellCoating] { &self.films }
    pub const fn options(&self) -> CoatingRelaxationOptions { self.options }

    pub(crate) fn from_projected(films: Vec<ProjectedShellCoating>, omega: &[f64],
        options: CoatingRelaxationOptions, max_work: usize, gate: &CancelGate)
        -> Result<Option<Self>, AcousticRealizeError>
    {
        poll(gate)?;
        options.validate()?;
        if films.is_empty() { return Ok(None); }
        let n = omega.len();
        if n == 0 || n > crate::shear_shell::MAX_SHELL_MODES
            || options.projection_terms(n, films.len()).is_none_or(|w| w > max_work)
            || films.iter().any(|f| f.unit_stiffness.len() != n*n) {
            return Err(invalid("coating memory requires a complete bounded common modal basis and factor work"));
        }
        let mut band = (0.0_f64, f64::MAX);
        let mut fits = Vec::<LoweredModel>::new();
        let mut indices = BTreeMap::<[u64; 7], usize>::new();
        let mut forms = BTreeMap::<u64, Vec<f64>>::new();
        for film in &films {
            poll(gate)?;
            let source = &film.source;
            band.0 = band.0.max(source.omega_band.0);
            band.1 = band.1.min(source.omega_band.1);
            if omega.iter().any(|w| !w.is_finite() || *w < band.0 || *w > band.1) {
                return Err(invalid("retained coating modes exceed the original measured fit support"));
            }
            let law = source.input.material.law();
            let source_band = law.band_rad_s();
            let key = [law.equilibrium_modulus_pa(), law.fit_modulus_pa(), law.alpha(),
                law.beta(), law.relaxation_time_s(), source_band.0, source_band.1].map(f64::to_bits);
            let index = if let Some(&index) = indices.get(&key) { index } else {
                let fit = law.lower_to_prony(options.terms, options.storage_relative_tolerance,
                    options.loss_relative_tolerance).map_err(owner)?;
                poll(gate)?;
                let index = fits.len();
                fits.push(fit); indices.insert(key, index); index
            };
            let fit = &fits[index];
            if fit.model.e_inf != source.reference_storage_pa {
                return Err(invalid("coating memory must use its original equilibrium shell stiffness exactly once"));
            }
            for &(increment, tau) in &fit.model.terms {
                if increment == 0.0 { continue; }
                poll(gate)?;
                let form = forms.entry(tau.to_bits()).or_insert_with(|| vec![0.0; n*n]);
                for (out, h) in form.iter_mut().zip(&film.unit_stiffness) { *out += increment*h; }
                if form.iter().any(|v| !v.is_finite()) {
                    return Err(invalid("coating memory form accumulation overflow"));
                }
            }
        }
        let mut rows = Vec::new();
        let mut times = Vec::new();
        for (tau, form) in forms {
            poll(gate)?;
            let factors = crate::plate_viscosity::relaxation::factor(&form, n, gate)?;
            if rows.len().checked_add(factors.len()).is_none_or(|count|
                count > options.max_memory_coordinates) {
                return Err(invalid("complete coating memory exceeds its declared state allowance"));
            }
            times.extend(std::iter::repeat_n(f64::from_bits(tau), factors.len()));
            rows.extend(factors);
        }
        poll(gate)?;
        Ok(Some(Self { films, fits, rows, times, modes: n, band, options }))
    }
}

impl ResolvedShearShell {
    /// The selected continuation policy, if any. Absence retains frequency-slice
    /// storage and modal film viscosity; the two paths never run together.
    pub const fn coating_relaxation_options(&self) -> Option<CoatingRelaxationOptions> {
        self.coating_relaxation
    }

    /// Select passive hereditary films and rebuild their equilibrium A/B/D at
    /// the same signed offsets, with unchanged physical mass and pressure skin.
    /// The finite source band stays authoritative for fit comparisons; choosing
    /// this method explicitly authors the Prony continuation for transients.
    pub fn with_coating_relaxation(&self, options: CoatingRelaxationOptions,
        gate: &CancelGate) -> Result<Self, AcousticRealizeError>
    {
        poll(gate)?; options.validate()?;
        let substrate = compile_material_shear_shell(self.mesh.clone(),
            self.regions.iter().map(|r| r.input.clone()).collect(), self.budget, gate)?;
        let mut rebuilt = apply_coatings(substrate,
            self.coatings.iter().map(|film| film.input.clone()).collect(), Some(options), gate)?
            .with_interfaces(self.interfaces.iter().map(|i| i.input().clone()).collect(),
                self.interface_budget, gate)?;
        rebuilt.pressure_surfaces = self.pressure_surfaces.clone();
        poll(gate)?;
        Ok(rebuilt)
    }

    /// Project source memory through the same physical modes and original MITC
    /// bubble transformation used for elastic, force, pressure and wood loss.
    pub fn project_coating_relaxation(&self, directors: &NodalDirectors,
        physical_modes: &[Vec<f64>], omega: &[f64], budget: MitcBudget,
        max_projection_terms: usize, gate: &CancelGate)
        -> Result<Option<ShellCoatingRelaxation>, AcousticRealizeError>
    {
        let Some(options) = self.coating_relaxation else { poll(gate)?; return Ok(None); };
        let forms = self.project_coating_forms(directors, physical_modes, omega,
            budget, max_projection_terms, gate)?;
        ShellCoatingRelaxation::from_projected(forms, omega, options, max_projection_terms, gate)
    }

    pub(crate) fn project_coating_forms(&self, directors: &NodalDirectors,
        physical_modes: &[Vec<f64>], omega: &[f64], budget: MitcBudget,
        max_projection_terms: usize, gate: &CancelGate)
        -> Result<Vec<ProjectedShellCoating>, AcousticRealizeError>
    {
        poll(gate)?;
        if self.coating_relaxation.is_none() { return Ok(Vec::new()); }
        let n = physical_modes.len();
        if n == 0 || n > crate::shear_shell::MAX_SHELL_MODES || omega.len() != n
            || self.loss_projection_terms(n).is_none_or(|work| work > max_projection_terms) {
            return Err(invalid("coating projection needs complete bounded modes and work"));
        }
        self.validate_coating_frequencies(omega, gate)?;
        let mut projected = Vec::with_capacity(self.coatings.len());
        for (source, viscosity) in self.coatings.iter().zip(&self.coating_viscosity) {
            poll(gate)?;
            let scale = source.reference_viscosity_pa_s;
            if !scale.is_finite() || scale <= 0.0 {
                return Err(invalid("coating unit strain form needs a positive representable source loss slice"));
            }
            let mut matrix = project_viscosity(&self.mesh, directors, &self.sections,
                viscosity, physical_modes, budget, max_projection_terms,
                || gate.is_requested()).map_err(owner)?;
            for value in &mut matrix { *value /= scale; }
            if matrix.iter().any(|v| !v.is_finite()) {
                return Err(invalid("coating unit modulus projection is unrepresentable"));
            }
            projected.push(ProjectedShellCoating { source: source.clone(), unit_stiffness: matrix });
        }
        poll(gate)?;
        Ok(projected)
    }
}
