//! Exterior air reaction of the SAME material shell and explicit cavity mouths.
//!
//! The original signed structural/flow traces drive the shared symmetric P1/P0
//! Galerkin owner. Its transpose work projection produces every body/mouth cross
//! term. The unchanged strict reciprocity/passivity gate precedes the existing
//! causal passive fit and complete mechanical propagator. No collocation result
//! is symmetrized, no resistance is added, and no mode is discarded.
//!
//! This is the original owner's bounded, below-resonance chart. Its geometric
//! slab exclusion and angular-resolution guard refuse wider unresolved bands;
//! there is no alternate-solver fallback. Spatial/angular refinement and the
//! unchanged finite-band realization gates remain necessary. The same P1/P0
//! solution supplies the optional receiver, not a pressure-mean conversion.
//! These discrete physical identities do not qualify a measured instrument.

use super::super::{ReciprocalStringBank, ReciprocalStringError};
use crate::render::plate::impact::radiation::{Model, fit};
use crate::shear_shell::ResolvedShearShell;
use crate::cavity_modal::CavityNeckLoss;
use super::super::multiple::performance::observation::{
    ShellObserverError, one_way::{ShellAcousticSkin, skin_identity, skin_traces},
};
use fs_bem::helmholtz::{HelmholtzError, Medium};
use fs_bem::modal_loading::{
    ModalLoadingError, ReciprocalLoadingDiagnostics, modal_radiation_loading_reciprocal,
};
use fs_bem::reciprocal_radiation::{
    Budget as ReciprocalBudget, MAX_ANGULAR_REFINEMENT, MAX_RECIPROCAL_PANELS, MAX_RECIPROCAL_UNKNOWNS,
};
#[cfg(test)]
use fs_bem::{helmholtz::{IterativeOptions,PanelRule},modal_loading::{
    LoadingRoute,ModalLoadingBake,SurfaceModeTraces,modal_radiation_loading_routed,
}};
use fs_bem::panel3d::SpherePanels;
use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_plate::shell::pressure::ShellFace;

/// Independent geometry, frequency, work and causal-realization limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShellBemOptions {
    /// Positive lower frequency [Hz].
    pub minimum_hz: f64,
    /// Upper frequency [Hz], covering every retained body reference mode.
    pub maximum_hz: f64,
    /// Odd 33..=257 frequencies: alternating training and withheld rows.
    pub frequency_samples: usize,
    /// All selected exterior triangles must fit; none are discarded.
    pub maximum_surface_panels: usize,
    /// Ceiling on the original Galerkin velocity_port_work times every sampled
    /// frequency, including singular moments, angular Gram assembly and full LU.
    pub maximum_sweep_terms: u64,
    /// Complete acoustic pole count, before the exact owner's stricter limits.
    pub maximum_poles: usize,
}
impl ShellBemOptions {
    /// The declared panel ceiling bounds the original P1/P0 unknowns and cached
    /// singular moments. The original owner still checks its exact byte/work
    /// requirements; these are conservative storage allowances, not estimates
    /// of accuracy or measured execution time.
    pub(in crate::bowed_string::runtime::reciprocal::bank) fn reciprocal_budget(
        &self, frequencies:usize, angular_refinement:u32,
    ) -> Result<ReciprocalBudget,ShellBemError> {
        if frequencies==0 || angular_refinement>MAX_ANGULAR_REFINEMENT {
            return Err(ShellBemError::Invalid("reciprocal source requires a nonempty grid and an admitted angular refinement"));
        }
        let panels=self.maximum_surface_panels;
        let cache=panels.checked_mul(panels).and_then(|n|n.checked_mul(128))
            .and_then(|n|panels.checked_mul(4096).and_then(|m|n.checked_add(m)))
            .ok_or(ShellBemError::Invalid("reciprocal cache allowance overflow"))?;
        Ok(ReciprocalBudget {
            maximum_unknowns:panels.saturating_mul(2).min(MAX_RECIPROCAL_UNKNOWNS),
            maximum_cache_bytes:cache,
            maximum_work:self.maximum_sweep_terms/(frequencies as u64),
            angular_refinement,
        })
    }

    /// Admission before geometry allocation or a Helmholtz factorization.
    pub fn validate(&self, rate: u32) -> Result<(), ShellBemError> {
        if !self.minimum_hz.is_finite()
            || self.minimum_hz <= 0.0
            || !self.maximum_hz.is_finite()
            || self.maximum_hz <= self.minimum_hz
            || !(33..=257).contains(&self.frequency_samples)
            || self.frequency_samples % 2 == 0
            || !(4..=MAX_RECIPROCAL_PANELS).contains(&self.maximum_surface_panels)
            || !(1..=100_000_000_000_000).contains(&self.maximum_sweep_terms)
            || !(1..=crate::render::plate::impact::radiation::MAX_POLES).contains(&self.maximum_poles)
            || rate == 0
            || 8.0 * core::f64::consts::TAU * self.maximum_hz
                >= 0.9 * core::f64::consts::PI * f64::from(rate)
        {
            return Err(ShellBemError::Invalid(
                "shell radiation needs a positive resolved band, odd 33..257 frequencies and bounded panel/work/pole allowances",
            ));
        }
        Ok(())
    }
}

/// Source binding and numerical diagnostics of the accepted physical load.
#[derive(Clone, Debug)]
pub struct ShellBemReport {
    /// Same signed material/geometry/modal basis pinned by the loaded bank.
    pub basis_identity: ContentHash,
    /// Shared loading declaration binding surface, traces, gas, grid and panel rule.
    pub loading_identity: ContentHash,
    /// Complete atmosphere shared by strings, body, cavity and pressure output.
    pub gas: GasState,
    pub options: ShellBemOptions,
    pub body_modes: usize,
    /// Actual compact aperture-flow ports after the complete structural basis.
    pub opening_ports: usize,
    pub surface_panels: usize,
    /// Actual original variational source diagnostics, including refinement/work.
    pub reciprocal: ReciprocalLoadingDiagnostics,
    pub acoustic_poles: usize,
    pub sweep_terms: u64,
    /// Lowest panel-resolution diagnostic on the complete sampled band.
    pub minimum_panels_per_wavelength: f64,
    /// Universal condition-number lower bound 1. The variational owner reports
    /// its actual backward residual instead of a condition probe or enclosure.
    pub maximum_condition_lower_bound: f64,
    /// Largest relative reciprocity defect admitted by the shared loading owner.
    pub reciprocity_relative_residual: f64,
    /// Withheld whole-matrix complex fit error under the existing fitter's gate.
    pub peak_error: f64,
    pub rms_error: f64,
    /// Independently gated resistance error; weak loss cannot hide in inertia.
    pub resistance_peak_error: f64,
    pub resistance_rms_error: f64,
    /// The declared acoustic skin the load was computed on; none on the shell path.
    pub skin_identity: Option<ContentHash>,
    /// Body modes above the loaded band, after the `body_modes` it loads. They
    /// carry no exterior load; the shell path loads every mode, so zero there.
    pub one_way_body_modes: usize,
}

#[derive(Debug)]
pub enum ShellBemError {
    Invalid(&'static str),
    Geometry(fs_bem::BemError),
    Solve(HelmholtzError),
    /// Shared full-matrix admission refused the sampled loading.
    Loading(ModalLoadingError),
    Fit(String),
    Bank(ReciprocalStringError),
    Cancelled,
}
impl core::fmt::Display for ShellBemError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material shell radiation: {self:?}")
    }
}
impl std::error::Error for ShellBemError {}
fn poll(gate: &CancelGate) -> Result<(), ShellBemError> {
    if gate.is_requested() {
        Err(ShellBemError::Cancelled)
    } else {
        Ok(())
    }
}

impl ReciprocalStringBank {
    /// Compute a closed shell's air load and attach it to the existing reciprocal
    /// mechanics. String/body/material changes rebuild the original modes and
    /// their radiation together. The same gas supplies c, rho and string loss.
    ///
    /// Body pressure reacts on the bridges and every string inside the existing
    /// exact propagator, friction solve, work balance and scheduled performance.
    /// No independent body oscillator, gain normalization or new stepper is used.
    pub fn with_material_shell_bem(
        self,
        options: ShellBemOptions,
        gate: &CancelGate,
    ) -> Result<(Self, ShellBemReport), ShellBemError> {
        self.with_material_shell_bem_with_angular_refinement(options,0,gate)
    }

    /// Independently refine the original outgoing angular Gram quadrature
    /// (0..=2), retaining the same declared geometry, band and work ceilings.
    pub fn with_material_shell_bem_with_angular_refinement(
        self, options:ShellBemOptions, angular_refinement:u32, gate:&CancelGate,
    ) -> Result<(Self,ShellBemReport),ShellBemError> {
        let (model, report) = self.prepare_material_shell_bem_with_angular_refinement(
            options,angular_refinement,gate)?;
        let mut bank = self
            .with_body_radiation_load(model, options.maximum_poles, gate)
            .map_err(ShellBemError::Bank)?;
        // Only this combined source-derived path can bind the preparation report
        // to the load that is now present in the original mechanical owner.
        let binding = bank.source.body_radiation.as_mut().ok_or(
            ShellBemError::Invalid("computed shell load lost its source binding"),
        )?;
        binding.shell_bem = Some(report.clone());
        poll(gate)?;
        Ok((bank, report))
    }

    /// Cold preparation without changing any mechanical state. Every original
    /// body mode participates, including those with zero net displaced volume.
    /// The result belongs to this exact material basis and gas; prefer the
    /// combined constructor above to preserve that binding automatically.
    ///
    /// A closed outward skin is required. Separate assembly parts may meet at
    /// identical coordinates; no tolerance-based welding or hole filling occurs.
    /// Exact edge closure is checked, not self-intersection or manifold quality.
    /// BEM is synchronous within one bounded frequency solve; cancellation is
    /// polled between frequencies, material parts and complete source batches.
    pub fn prepare_material_shell_bem(
        &self,
        options: ShellBemOptions,
        gate: &CancelGate,
    ) -> Result<(Model, ShellBemReport), ShellBemError> {
        self.prepare_material_shell_bem_with_angular_refinement(options,0,gate)
    }

    /// Preparation with explicitly selected angular refinement of the SAME
    /// variational source. Unresolved bands/budgets refuse without changing
    /// owner, weakening admission or returning a partial fitted model.
    pub fn prepare_material_shell_bem_with_angular_refinement(
        &self, options:ShellBemOptions, angular_refinement:u32, gate:&CancelGate,
    ) -> Result<(Model,ShellBemReport),ShellBemError> {
        poll(gate)?;
        options.validate(self.source.options.sample_rate_hz)?;
        // A prepared-linear bank carries the load as its air bank beside the transition.
        if self.system.samples_rendered() != 0 || self.source.body_radiation.is_some() {
            return Err(ShellBemError::Invalid(
                "computed shell radiation requires an unstarted bank without a prior exterior load",
            ));
        }
        if self.source.material_cavity.as_ref().is_some_and(|c| c.cavity.loss() == CavityNeckLoss::ExteriorMouth) {
            return Err(ShellBemError::Invalid(
                "an exterior-owned mouth radiates through the caps of a declared skin",
            ));
        }
        let (basis_identity, _) = self
            .source
            .body_radiation_basis(&self.system, gate)
            .map_err(ShellBemError::Bank)?;
        let modes = self
            .system
            .components()
            .get(self.source.strings.len())
            .ok_or(ShellBemError::Invalid(
                "missing original material body component",
            ))?
            .modes();
        if modes.is_empty()
            || modes.len() > crate::render::plate::impact::radiation::MAX_PORTS
            || modes.iter().any(|m| {
                let hz = m.angular_frequency_rad_s / core::f64::consts::TAU;
                hz < options.minimum_hz || hz > options.maximum_hz
            })
        {
            return Err(ShellBemError::Invalid(
                "shell radiation must cover every original body mode within the shared passive port capacity",
            ));
        }
        let opening_ports = self.source.spatial_cavity.as_ref().map_or(0, |c|c.opening_count());
        if modes.len()+opening_ports > crate::render::plate::impact::radiation::MAX_PORTS {
            return Err(ShellBemError::Invalid("complete body/mouth basis exceeds passive port capacity"));
        }
        if let Some(cavity) = self.source.spatial_cavity.as_ref().filter(|c|c.opening_count() != 0) {
            if cavity.options().maximum_hz > options.maximum_hz {
                return Err(ShellBemError::Invalid("joint exterior band must cover the complete interior pressure window"));
            }
            for index in 0..cavity.opening_count() {
                let reference = cavity.coupling().neck_fixed_wall_angular_frequency_rad_s(index)
                    .map_err(|_|ShellBemError::Invalid("missing original compact neck frequency"))?;
                let extent = cavity.opening_extent_m(index)
                    .ok_or(ShellBemError::Invalid("missing original mouth extent"))?;
                if reference/core::f64::consts::TAU > options.maximum_hz
                    || core::f64::consts::TAU*options.maximum_hz*extent/self.source.gas.sound_speed > 0.3
                {
                    return Err(ShellBemError::Invalid("exterior band must cover neck dynamics inside the compact mouth chart"));
                }
            }
        }
        let (surface, fields, _) = self.source.shell_surface_traces(
            modes.len(), options.maximum_surface_panels, options.maximum_sweep_terms,
            options.frequency_samples, gate,
        )?;
        self.fit_exterior_load(options, angular_refinement, basis_identity, modes.len(),
            opening_ports, &surface, &fields, gate)
    }

    /// The shared Galerkin loading, passive fit and report of one complete
    /// exterior source on `surface`: `body_modes` structural traces first,
    /// then `opening_ports` mouth traces, in the order of the bank's ports.
    #[allow(clippy::too_many_arguments)]
    fn fit_exterior_load(
        &self, options:ShellBemOptions, angular_refinement:u32, basis_identity:ContentHash,
        body_modes:usize, opening_ports:usize, surface:&SpherePanels, fields:&[Vec<f64>],
        gate:&CancelGate,
    ) -> Result<(Model,ShellBemReport),ShellBemError> {
        let panels = surface.areas().len();
        let gas = self.source.gas;
        gas.validate_coherence()
            .map_err(|_| ShellBemError::Invalid("incoherent shell exterior gas"))?;
        let medium = Medium {
            density: gas.density,
            sound_speed: gas.sound_speed,
        };
        let low = core::f64::consts::TAU * options.minimum_hz;
        let high = core::f64::consts::TAU * options.maximum_hz;
        let mut omega: Vec<_> = (0..options.frequency_samples)
            .map(|i| low + (high - low) * i as f64 / (options.frequency_samples - 1) as f64)
            .collect();
        omega[0] = low;
        omega[options.frequency_samples - 1] = high;
        let entries=fields.len().checked_mul(fields.len()).and_then(|n|n.checked_mul(omega.len()))
            .ok_or(ShellBemError::Invalid("complete reciprocal matrix extent overflow"))?;
        let source=modal_radiation_loading_reciprocal(surface,fields,&omega,medium,
            options.reciprocal_budget(omega.len(),angular_refinement)?,basis_identity,entries,
            ||gate.is_requested()).map_err(|error|match error {
                ModalLoadingError::Helmholtz(HelmholtzError::Cancelled)=>ShellBemError::Cancelled,
                ModalLoadingError::Helmholtz(error)=>ShellBemError::Solve(error),
                other=>ShellBemError::Loading(other),
            })?;
        let loading=&source.loading;
        let reciprocal=source.diagnostics;
        let mut impedance = Vec::with_capacity(omega.len());
        let mut reciprocity = 0.0_f64;
        for row in loading.rows() {
            poll(gate)?;
            impedance.push(row.impedance().to_vec());
            reciprocity = reciprocity.max(row.reciprocity_relative_defect());
        }
        let mut extent = 0.0_f64;
        for a in surface.centroids() {
            poll(gate)?;
            for b in surface.centroids() {
                extent = extent.max(fs_math::det::hypot(
                    fs_math::det::hypot(a[0]-b[0],a[1]-b[1]),a[2]-b[2]));
            }
        }
        let spacing = if extent > 0.0 {
            (core::f64::consts::TAU * gas.sound_speed / (6.0 * extent)).min(high - low)
        } else {
            high - low
        };
        let fitted = fit::fit_resolved(&omega, &impedance, fields.len(), spacing)
            .map_err(ShellBemError::Fit)?;
        poll(gate)?;
        if fitted.model.poles.len() > options.maximum_poles {
            return Err(ShellBemError::Invalid(
                "complete shell acoustic fit exceeds its pole allowance",
            ));
        }
        let report = ShellBemReport {
            basis_identity,
            loading_identity: loading.declaration_hash(),
            gas,
            options,
            body_modes,
            opening_ports,
            surface_panels: panels,
            reciprocal,
            acoustic_poles: fitted.model.poles.len(),
            sweep_terms: reciprocal.sweep_work,
            minimum_panels_per_wavelength: reciprocal.minimum_panels_per_wavelength,
            maximum_condition_lower_bound: 1.0,
            reciprocity_relative_residual: reciprocity,
            peak_error: fitted.peak_error,
            rms_error: fitted.rms_error,
            resistance_peak_error: fitted.resistance_peak_error,
            resistance_rms_error: fitted.resistance_rms_error,
            skin_identity: None,
            one_way_body_modes: 0,
        };
        Ok((fitted.model, report))
    }

    /// The same Galerkin air load on a declared closed acoustic skin of this
    /// body, the skin its one-way observer reads, attached to the original
    /// mechanics as by [`Self::with_material_shell_bem`].
    pub fn with_material_skin_bem(
        self, options:ShellBemOptions, skin:&ShellAcousticSkin, maximum_anchor_distance_m:f64,
        angular_refinement:u32, gate:&CancelGate,
    ) -> Result<(Self,ShellBemReport),ShellBemError> {
        let (model, report) = self.prepare_material_skin_bem(
            options, skin, maximum_anchor_distance_m, angular_refinement, gate)?;
        // Body coordinates relative to the body's first mode: the loaded prefix,
        // then the neck, the first coordinate of the component after the body.
        let structural = report.body_modes + report.one_way_body_modes;
        let mut coordinates: Vec<_> = (0..report.body_modes).collect();
        if report.opening_ports == 1 {
            coordinates.push(structural);
        }
        let mut bank = self
            .attach_body_radiation(model, options.maximum_poles, Some(coordinates), gate)
            .map_err(ShellBemError::Bank)?;
        let binding = bank.source.body_radiation.as_mut().ok_or(
            ShellBemError::Invalid("computed skin load lost its source binding"),
        )?;
        binding.shell_bem = Some(report.clone());
        poll(gate)?;
        Ok((bank, report))
    }

    /// Cold preparation on `skin` without changing any mechanical state.
    /// Every skin vertex follows its anchored points on the original parts,
    /// as in the observer. A lumped cavity's mouths are the skin's aperture
    /// caps, carrying its neck volume flow in the declared shares: one port on
    /// the one neck coordinate. Its outer end must be exterior-owned
    /// ([`CavityNeckLoss::ExteriorMouth`]); a compact outer end already holds
    /// that exterior mass and resistance, and refuses here.
    pub fn prepare_material_skin_bem(
        &self, options:ShellBemOptions, skin:&ShellAcousticSkin, maximum_anchor_distance_m:f64,
        angular_refinement:u32, gate:&CancelGate,
    ) -> Result<(Model,ShellBemReport),ShellBemError> {
        poll(gate)?;
        options.validate(self.source.options.sample_rate_hz)?;
        if self.system.samples_rendered() != 0 || self.source.body_radiation.is_some() {
            return Err(ShellBemError::Invalid(
                "computed skin radiation requires an unstarted bank without a prior exterior load",
            ));
        }
        if self.source.spatial_cavity.is_some() {
            return Err(ShellBemError::Invalid(
                "a spatial cavity's mouths are ports of the shell path, not skin apertures",
            ));
        }
        let (basis_identity, _) = self.source.body_radiation_basis(&self.system, gate)
            .map_err(ShellBemError::Bank)?;
        let body = self.source.strings.len();
        let modes = self.system.components().get(body)
            .ok_or(ShellBemError::Invalid("missing original material body component"))?
            .modes();
        // The load covers the ascending prefix of body modes inside the band. The
        // owner's uniqueness chart caps the band; modes above it stay one-way.
        let hz = |m: &crate::modal_acoustic_time::ModalAcousticMode| {
            m.angular_frequency_rad_s / core::f64::consts::TAU
        };
        let loaded = modes.iter().take_while(|m| hz(m) <= options.maximum_hz).count();
        if loaded == 0
            || loaded > crate::render::plate::impact::radiation::MAX_PORTS
            || modes[..loaded].iter().any(|m| hz(m) < options.minimum_hz)
            || modes[loaded..].iter().any(|m| hz(m) <= options.maximum_hz)
        {
            return Err(ShellBemError::Invalid(
                "skin radiation loads an ascending prefix of the body modes from the band's lower edge",
            ));
        }
        let mut parts = Vec::new();
        if let Some(image) = &self.source.material_shell_image {
            parts.push((image.source().mesh(), image.physical_modes()));
        } else if let Some(image) = &self.source.material_shell_assembly_image {
            for (part, rows) in image.source().components.iter().zip(image.physical_modes()) {
                parts.push((part.specimen.mesh(), rows.as_slice()));
            }
        } else {
            return Err(ShellBemError::Invalid(
                "computed skin radiation requires an original material shell basis",
            ));
        }
        let traces = skin_traces(&parts, modes.len(), skin, options.maximum_surface_panels,
            maximum_anchor_distance_m, gate).map_err(|error| match error {
                ShellObserverError::Cancelled => ShellBemError::Cancelled,
                ShellObserverError::Invalid(what) => ShellBemError::Invalid(what),
                _ => ShellBemError::Invalid("acoustic skin traces refused"),
            })?;
        let mut fields = traces.modes;
        fields.truncate(loaded);
        match (self.source.material_cavity.as_ref(), skin.apertures.is_empty()) {
            (None, true) => {}
            (Some(cavity), false) if cavity.cavity.loss() == CavityNeckLoss::ExteriorMouth => {
                let neck = self.system.components().get(body + 1)
                    .ok_or(ShellBemError::Invalid("skin load lost the cavity neck component"))?;
                if neck.modes().len() != 1 {
                    return Err(ShellBemError::Invalid("skin load needs the one physical neck coordinate"));
                }
                // The band must reach the neck's own resonance with its declared compact ends.
                let geometry = cavity.cavity.geometry();
                let gas = cavity.cavity.gas();
                let neck_hz = fs_phs::helmholtz_natural_frequency(geometry.volume_m3,
                    geometry.neck_radius_m, geometry.neck_length_m, gas.density, gas.sound_speed,
                    crate::cavity_modal::neck_ends(geometry))
                    .map_err(|_| ShellBemError::Invalid("missing original compact neck frequency"))?
                    / core::f64::consts::TAU;
                if options.minimum_hz > neck_hz {
                    return Err(ShellBemError::Invalid("exterior band must reach the cavity neck resonance"));
                }
                let shape = cavity.cavity.neck_volume_shape_m2_per_sqrt_kg();
                let mut field = vec![0.0; traces.surface.areas().len()];
                for (aperture, trace) in skin.apertures.iter().zip(&traces.apertures) {
                    for (value, cap) in field.iter_mut().zip(trace) {
                        *value += aperture.flow_share * shape * cap;
                    }
                }
                fields.push(field);
            }
            (Some(_), false) => return Err(ShellBemError::Invalid(
                "a compact outer neck end already carries the mouth's exterior load; declare it exterior-owned",
            )),
            (Some(_), true) => return Err(ShellBemError::Invalid(
                "a lumped cavity radiates through declared skin apertures",
            )),
            (None, false) => return Err(ShellBemError::Invalid(
                "skin apertures need a lumped cavity neck",
            )),
        }
        if fields.len() > crate::render::plate::impact::radiation::MAX_PORTS {
            return Err(ShellBemError::Invalid("complete body/mouth skin source exceeds passive port capacity"));
        }
        let opening_ports = fields.len() - loaded;
        let work = sweep_work(traces.surface.areas().len(), fields.len(), options.frequency_samples)
            .ok_or(ShellBemError::Invalid("complete BEM sweep extent overflow"))?;
        if work > options.maximum_sweep_terms {
            return Err(ShellBemError::Invalid("complete skin exterior exceeds its sweep-work allowance"));
        }
        let (model, mut report) = self.fit_exterior_load(options, angular_refinement,
            basis_identity, loaded, opening_ports, &traces.surface, &fields, gate)?;
        report.skin_identity = Some(skin_identity(skin));
        report.one_way_body_modes = modes.len() - loaded;
        Ok((model, report))
    }
}

impl super::super::MaterialStringBankSource {
    /// One original geometry/normal-trace extraction shared by the air reaction
    /// and its exterior observation. The complete source basis is retained.
    pub(in crate::bowed_string::runtime::reciprocal::bank) fn shell_surface_traces(
        &self,
        modes: usize,
        maximum_surface_panels: usize,
        maximum_sweep_terms: u64,
        frequency_samples: usize,
        gate: &CancelGate,
    ) -> Result<(SpherePanels, Vec<Vec<f64>>, u64), ShellBemError> {
        if modes == 0 || modes > crate::render::plate::impact::radiation::MAX_PORTS
            || frequency_samples == 0 || !(4..=MAX_RECIPROCAL_PANELS).contains(&maximum_surface_panels)
        {
            return Err(ShellBemError::Invalid("invalid complete shell trace/work dimensions"));
        }
        let mut parts = Vec::new();
        if let Some(image) = &self.material_shell_image {
            parts.push((image.source(), image.physical_modes()));
        } else if let Some(image) = &self.material_shell_assembly_image {
            for (part, rows) in image.source().components.iter().zip(image.physical_modes()) {
                parts.push((&part.specimen, rows.as_slice()));
            }
        } else {
            return Err(ShellBemError::Invalid(
                "computed closed-body radiation requires an original material shell basis",
            ));
        }
        let body_panels = parts
            .iter()
            .try_fold(0usize, |count, (shell, _)| {
                count.checked_add(
                    shell
                        .pressure_surfaces()
                        .map_or(shell.mesh().tris.len(), |s| s.exterior.len()),
                )
            })
            .ok_or(ShellBemError::Invalid(
                "shell exterior panel extent overflow",
            ))?;
        let cavity = self.spatial_cavity.as_ref();
        let opening_ports = cavity.map_or(0, |c|c.opening_count());
        let opening_panels = if let Some(cavity) = cavity {
            (0..opening_ports).try_fold(0usize,|sum,index| {
                cavity.opening_triangles(index).and_then(|faces|sum.checked_add(faces.len()))
            }).ok_or(ShellBemError::Invalid("mouth exterior panel extent overflow"))?
        } else { 0 };
        let panels = body_panels.checked_add(opening_panels)
            .ok_or(ShellBemError::Invalid("complete exterior panel extent overflow"))?;
        let ports = modes.checked_add(opening_ports)
            .ok_or(ShellBemError::Invalid("complete exterior port extent overflow"))?;
        if ports > crate::render::plate::impact::radiation::MAX_PORTS {
            return Err(ShellBemError::Invalid("complete body/mouth source exceeds passive port capacity"));
        }
        let work = sweep_work(panels, ports, frequency_samples)
            .ok_or(ShellBemError::Invalid("complete BEM sweep extent overflow"))?;
        if panels < 4
            || panels > maximum_surface_panels
            || work > maximum_sweep_terms
        {
            return Err(ShellBemError::Invalid(
                "complete shell exterior exceeds its panel or sweep-work allowance",
            ));
        }
        let (surface, fields) = surface_trace(&parts, modes, cavity, gate)?;
        Ok((surface, fields, work))
    }
}

// Dense sweep extent of one exterior source: p^2 (p + n) + p n^2 per frequency.
fn sweep_work(panels: usize, ports: usize, frequency_samples: usize) -> Option<u64> {
    let p = panels as u64;
    let n = ports as u64;
    p.checked_mul(p)
        .and_then(|p2| p.checked_add(n).and_then(|sum| p2.checked_mul(sum)))
        .and_then(|v| {
            n.checked_mul(n)
                .and_then(|n2| p.checked_mul(n2))
                .and_then(|w| v.checked_add(w))
        })
        .and_then(|v| v.checked_mul(frequency_samples as u64))
}

// BEM velocity and pressure use the SAME constant-panel normal trace. Rotations
// affect the elastic modes, but have no extra midsurface translation to add.
fn surface_trace(
    parts: &[(&ResolvedShearShell, &[Vec<f64>])],
    modes: usize,
    cavity: Option<&super::super::cavity::MaterialShellCavity>,
    gate: &CancelGate,
) -> Result<(SpherePanels, Vec<Vec<f64>>), ShellBemError> {
    let mut triangles = Vec::new();
    let mut node_rows = Vec::new();
    for &(shell, rows) in parts {
        poll(gate)?;
        let mesh = shell.mesh();
        if rows.len() != modes
            || rows
                .iter()
                .any(|r| r.len() != 6 * mesh.nodes.len() || r.iter().any(|v| !v.is_finite()))
        {
            return Err(ShellBemError::Invalid(
                "exterior trace must retain the complete original shell basis",
            ));
        }
        let facets: Vec<_> = shell.pressure_surfaces().map_or_else(
            || {
                (0..mesh.tris.len())
                    .map(|i| (i, ShellFace::Positive))
                    .collect()
            },
            |s| s.exterior.iter().map(|f| (f.triangle, f.side)).collect(),
        );
        for (index, side) in facets {
            poll(gate)?;
            let mut nodes = mesh.tris[index];
            if side == ShellFace::Negative {
                nodes.swap(1, 2);
            }
            triangles.push(nodes.map(|n| mesh.nodes[n]));
            let mut velocities = Vec::with_capacity(modes);
            for row in rows {
                velocities.push(core::array::from_fn::<_, 3, _>(|axis| {
                    (row[6 * nodes[0] + axis] + row[6 * nodes[1] + axis] + row[6 * nodes[2] + axis])
                        / 3.0
                }));
            }
            node_rows.push(velocities);
        }
    }
    let body_panels = triangles.len();
    let opening_ports = cavity.map_or(0, |c|c.opening_count());
    if let Some(cavity) = cavity {
        for index in 0..opening_ports {
            poll(gate)?;
            triangles.extend_from_slice(cavity.opening_triangles(index)
                .ok_or(ShellBemError::Invalid("lost original mouth boundary triangles"))?);
        }
    }
    let surface = SpherePanels::from_triangles(triangles).map_err(ShellBemError::Geometry)?;
    let mut fields = vec![vec![0.0;surface.areas().len()]; modes+opening_ports];
    for (panel, (normal, velocities)) in surface.normals().iter().zip(node_rows).enumerate() {
        poll(gate)?;
        for (field, velocity) in fields[..modes].iter_mut().zip(velocities) {
            let value = normal[0] * velocity[0] + normal[1] * velocity[1] + normal[2] * velocity[2];
            if !value.is_finite() {
                return Err(ShellBemError::Invalid(
                    "nonfinite material shell normal trace",
                ));
            }
            field[panel] = value;
        }
    }
    if let Some(cavity) = cavity {
        let mut first = body_panels;
        for index in 0..opening_ports {
            poll(gate)?;
            let count = cavity.opening_triangles(index)
                .ok_or(ShellBemError::Invalid("lost original mouth boundary triangles"))?.len();
            let port = cavity.coupling().neck_radiation_port(index)
                .map_err(|_|ShellBemError::Invalid("lost original mouth volume-flow port"))?;
            let value = port.volume_weight_m2_per_sqrt_kg/port.area_m2;
            if !value.is_finite() || value <= 0.0 {
                return Err(ShellBemError::Invalid("mouth normal velocity scaling is unrepresentable"));
            }
            fields[modes+index][first..first+count].fill(value);
            first += count;
        }
    }
    Ok((surface, fields))
}

// Retain the instrument's physical trace and work ceiling while delegating
// geometry admission, formulation, projection and matrix assessment together.
#[cfg(test)]
fn surface_loading(
    surface: &SpherePanels,
    fields: &[Vec<f64>],
    omega: &[f64],
    medium: Medium,
    basis_identity: ContentHash,
    gate: &CancelGate,
) -> Result<ModalLoadingBake, ShellBemError> {
    let entries = fields
        .len()
        .checked_mul(fields.len())
        .and_then(|n| n.checked_mul(omega.len()))
        .ok_or(ShellBemError::Invalid("shell modal matrix extent overflow"))?;
    // The existing 2048-panel ceiling stays below the shared dense cap. Only
    // this panel rule reaches the dense owner; there is no iterative fallback.
    let numerical = IterativeOptions {
        panel_rule: PanelRule::Triangle {
            max_pairs: surface.areas().len().saturating_pow(2),
        },
        ..IterativeOptions::default()
    };
    modal_radiation_loading_routed(
        SurfaceModeTraces::ClosedBody {
            surface,
            shapes: fields,
        },
        LoadingRoute::ClosedBodyGeometryPolicy,
        omega,
        medium,
        basis_identity,
        entries,
        Some(&numerical),
        || gate.is_requested(),
    )
    .map_err(|error| match error {
        ModalLoadingError::Helmholtz(HelmholtzError::Cancelled) => ShellBemError::Cancelled,
        ModalLoadingError::Helmholtz(error) => ShellBemError::Solve(error),
        other => ShellBemError::Loading(other),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_bem::helmholtz::Formulation;
    use fs_material::gas::GasSpec;

    fn source() -> ContentHash {
        fs_blake3::hash_domain("material-shell-bem-test", b"signed physical surface modes")
    }

    #[test]
    fn signed_zero_volume_fields_keep_mutual_air_reaction_and_the_selected_medium() {
        let gate = CancelGate::new_clock_free();
        let surface = SpherePanels::icosphere(0.1, 1).unwrap();
        // A translating sphere is a dipole: zero total displaced volume does
        // not mean zero exterior pressure, inertia, or radiation resistance.
        let field: Vec<_> = surface.normals().iter().map(|n| n[0]).collect();
        let volume: f64 = field.iter().zip(surface.areas()).map(|(v, a)| v * a).sum();
        assert!(volume.abs() < 1e-15);
        let opposite: Vec<_> = field.iter().map(|v| -v).collect();
        let fields = vec![field, opposite];
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
        let medium = Medium::from(&gas);
        let omega = medium.sound_speed / 0.1;
        let loading =
            surface_loading(&surface, &fields, &[omega], medium, source(), &gate).unwrap();
        let row = &loading.rows()[0];
        let diagnostics = row.bem_diagnostics().unwrap();
        assert_eq!(diagnostics.formulation, Formulation::PlainCbie);
        assert!(diagnostics.minimum_panels_per_wavelength >= 6.0);
        assert!(diagnostics.maximum_condition_lower_bound >= 1.0);
        assert!(row.reciprocity_relative_defect() <= 1e-12);
        let z = row.impedance();
        assert!(z[0].re > 0.0 && z[0].abs() > 1e-6);
        for (a, b) in [(z[0], z[3]), (z[0], z[1].scale(-1.0)), (z[1], z[2])] {
            assert!((a - b).abs() < 1e-12 * z[0].abs());
        }
        let denser = surface_loading(
            &surface,
            &fields,
            &[omega],
            Medium {
                density: 2.0 * medium.density,
                ..medium
            },
            source(),
            &gate,
        )
        .unwrap();
        assert_ne!(loading.declaration_hash(), denser.declaration_hash());
        for (a, b) in z.iter().zip(denser.rows()[0].impedance()) {
            assert!((*b - a.scale(2.0)).abs() < 1e-12 * a.abs());
        }
        let faster = surface_loading(
            &surface,
            &fields,
            &[omega],
            Medium {
                sound_speed: 2.0 * medium.sound_speed,
                ..medium
            },
            source(),
            &gate,
        )
        .unwrap();
        let faster = faster.rows()[0].impedance();
        assert!((faster[0] - z[0]).abs() > 0.01 * z[0].abs());
        // This is a phase/impedance change, not a scalar gain adjustment.
        assert!((faster[0].re / faster[0].im - z[0].re / z[0].im).abs() > 0.01);
    }

    #[test]
    fn every_surface_edge_and_outward_component_is_required() {
        let gate = CancelGate::new_clock_free();
        let medium = Medium::air();
        let admit = |triangles: Vec<[[f64; 3]; 3]>| -> Result<ModalLoadingBake, ShellBemError> {
            let surface =
                SpherePanels::from_triangles(triangles).map_err(ShellBemError::Geometry)?;
            let fields = vec![vec![1.0; surface.areas().len()]];
            surface_loading(
                &surface,
                &fields,
                &[medium.sound_speed / 0.1],
                medium,
                source(),
                &gate,
            )
        };
        let surface = SpherePanels::icosphere(0.1, 0).unwrap();
        let triangles = surface.triangles().unwrap();
        admit(triangles.to_vec()).unwrap();
        assert!(admit(triangles[1..].to_vec()).is_err());
        let mut inside_out = triangles.to_vec();
        for tri in &mut inside_out {
            tri.swap(1, 2);
        }
        assert!(admit(inside_out).is_err());
        let mut bad_seam = triangles.to_vec();
        bad_seam[0][0][0] += 1e-12;
        assert!(admit(bad_seam).is_err());
        let mut duplicate = triangles.to_vec();
        duplicate.push(triangles[0]);
        assert!(admit(duplicate).is_err());
        assert!(admit(vec![]).is_err());
    }
}
