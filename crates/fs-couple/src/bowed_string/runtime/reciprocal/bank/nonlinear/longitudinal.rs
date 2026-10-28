//! Material-derived axial inertia in the same bowed-string/body/fluid network.
//!
//! Reuses fs-nlmodal's positive finite-strain energy and the existing exact
//! linear/friction/DG owners. Retained longitudinal sine modes carry the SOLID
//! mass; transverse slope maps carry their already admitted gas-loaded mass.
//! The fixed-end elastic axial extension is conservative. Original bending,
//! torsion, body, bridge and fluid losses remain; no axial viscosity is inferred
//! from a bending-loss property. Omitted axial channels equilibrate statically.
use super::super::multiple::MultiBowStringBank;
use super::super::{MaterialStringBankSource, MaterialStringSource, ReciprocalStringBank};
use super::{ImpactError, invalid};
use crate::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeModel};
use crate::render::schedule::force::coupled::ModalCouplingError;
use crate::render::schedule::force::coupled::contact::friction::nonlinear::{
    MappedModalPotentialBinding, ModalCoordinate, ModalPotential, NonlinearModalConfig,
};
use fs_exec::CancelGate;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::section::axial::ConcentricAxialViscosity;
use fs_math::{c64::C64, det};
use fs_nlmodal::{KcStringParams, longitudinal::LongitudinalString};
use fs_phs::Storage;

mod loss;
pub use loss::LongitudinalMaterialLoss;
pub mod bearings;
pub use bearings::{AxialBearingGeometry, MaterialAxialBearing};

/// Explicit retained physics and constitutive domain, independent of instrument.
#[derive(Clone, Debug, PartialEq)]
pub struct LongitudinalStringResolution {
    /// Strictly increasing one-based longitudinal sine indices, at most 32.
    /// A selected index is never dropped to fit the clock or a material band.
    pub axial_modes: Vec<usize>,
    /// Conservative continuous-span transverse slope bound in (0, 0.3].
    pub maximum_slope: f64,
    /// Maximum absolute dynamic axial strain in (0, 0.1]. Also require positive
    /// total local tension throughout the continuous span.
    pub maximum_strain: f64,
}
impl LongitudinalStringResolution {
    /// Validate before allocation or source I/O.
    pub fn validate(&self) -> Result<(), ImpactError> {
        if self.axial_modes.is_empty()
            || self.axial_modes.len() > fs_nlmodal::longitudinal::MAX_LONGITUDINAL_MODES
            || self.axial_modes[0] == 0
            || self
                .axial_modes
                .last()
                .is_some_and(|&n| n > fs_nlmodal::longitudinal::MAX_LONGITUDINAL_INDEX)
            || self.axial_modes.windows(2).any(|w| w[0] >= w[1])
            || !self.maximum_slope.is_finite()
            || !(0.0 < self.maximum_slope && self.maximum_slope <= 0.3)
            || !self.maximum_strain.is_finite()
            || !(0.0 < self.maximum_strain && self.maximum_strain <= 0.1)
        {
            return Err(invalid(
                "longitudinal waves require ordered explicit modes and finite moderate-slope/strain bounds",
            ));
        }
        Ok(())
    }
}

/// One immutable constitutive map; original material receipts stay on its member.
#[derive(Clone, Debug)]
pub struct MaterialLongitudinalString {
    model: LongitudinalString,
    resolution: LongitudinalStringResolution,
    transverse_modes: usize,
    planes: usize,
    slope_scales: Vec<f64>,
    solid_linear_density_kg_m: f64,
    transverse_linear_density_kg_m: f64,
    axial_rigidity_n: f64,
    prestress_n: f64,
    axial_component: usize,
    elastic_frequency_band_rad_s: Option<(f64, f64)>,
    axial_material_loss: Option<ConcentricAxialViscosity>,
}
impl MaterialLongitudinalString {
    /// Existing finite-strain owner, including the material axial frequencies.
    pub const fn model(&self) -> &LongitudinalString {
        &self.model
    }
    /// Actual solid mass per length, never borrowed from transverse gas loading.
    pub const fn solid_linear_density_kg_m(&self) -> f64 {
        self.solid_linear_density_kg_m
    }
    /// Actual transverse mass per length used by the original bank coordinates.
    pub const fn transverse_linear_density_kg_m(&self) -> f64 {
        self.transverse_linear_density_kg_m
    }
    /// EA from the same original material/section as the transverse source.
    pub const fn axial_rigidity_n(&self) -> f64 {
        self.axial_rigidity_n
    }
    /// New axial component follows every original string, body and fluid state.
    pub const fn axial_component(&self) -> usize {
        self.axial_component
    }
    /// Intersection of explicit elasticity/density frequency applicability.
    pub const fn elastic_frequency_band_rad_s(&self) -> Option<(f64, f64)> {
        self.elastic_frequency_band_rad_s
    }
    /// Original per-layer axial loss, absent only for an explicit elastic selection.
    /// The retained law supplies linear u_xt loss, not finite-slope strain-rate loss.
    pub const fn axial_material_loss(&self) -> Option<&ConcentricAxialViscosity> {
        self.axial_material_loss.as_ref()
    }
    /// Declared modal selection and finite-motion envelope.
    pub const fn resolution(&self) -> &LongitudinalStringResolution {
        &self.resolution
    }
    fn transverse_count(&self) -> usize {
        self.transverse_modes * self.planes
    }
    fn slopes(&self, q: &[f64]) -> Vec<f64> {
        let mut slopes = vec![0.0; self.model.slope_count()];
        for plane in 0..self.planes {
            for mode in 0..self.transverse_modes {
                slopes[plane * (self.transverse_modes + 1) + mode + 1] =
                    self.slope_scales[mode] * q[plane * self.transverse_modes + mode];
            }
        }
        slopes
    }
    fn observe(&self, q: &[f64]) -> Result<LongitudinalStringObservation, ModalCouplingError> {
        self.validate_coordinates(q)?;
        let slopes = self.slopes(q);
        let axial = &q[self.transverse_count()..];
        let (mean, variation) = self.model.strain_bounds(&slopes, axial).map_err(owner)?;
        Ok(LongitudinalStringObservation {
            mean_dynamic_strain: mean,
            maximum_abs_dynamic_strain: mean.abs() + variation,
            minimum_tension_n: self.prestress_n + self.axial_rigidity_n * (mean - variation),
            maximum_tension_n: self.prestress_n + self.axial_rigidity_n * (mean + variation),
            strain_energy_j: self.model.energy(&slopes, axial).map_err(owner)?,
        })
    }
}
/// Conservative continuous-span bounds, not samples at a few string stations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LongitudinalStringObservation {
    pub mean_dynamic_strain: f64,
    pub maximum_abs_dynamic_strain: f64,
    pub minimum_tension_n: f64,
    pub maximum_tension_n: f64,
    /// Full positive strain energy, including retained axial elastic storage.
    pub strain_energy_j: f64,
}
#[derive(Clone, Debug)]
pub struct StringBankLongitudinal {
    strings: Vec<MaterialLongitudinalString>,
    config: NonlinearModalConfig,
    bearings: Vec<MaterialAxialBearing>,
}
impl StringBankLongitudinal {
    /// Original material-string declaration order.
    pub fn strings(&self) -> &[MaterialLongitudinalString] {
        &self.strings
    }
    pub const fn config(&self) -> NonlinearModalConfig {
        self.config
    }
    /// Original axial bearing geometry and complete resolved surface laws.
    pub fn bearings(&self) -> &[MaterialAxialBearing] {
        &self.bearings
    }
}
impl MaterialStringBankSource {
    pub const fn longitudinal_waves(&self) -> Option<&StringBankLongitudinal> {
        self.longitudinal_waves.as_ref()
    }
}
fn owner(_: fs_nlmodal::NlModalError) -> ModalCouplingError {
    ModalCouplingError::Invalid {
        what: "material longitudinal strain is outside its finite constitutive domain",
    }
}
impl Storage for MaterialLongitudinalString {
    fn hamiltonian(&self, q: &[f64]) -> f64 {
        if q.len() != self.coordinate_count() {
            return f64::NAN;
        }
        self.model
            .energy(&self.slopes(q), &q[self.transverse_count()..])
            .unwrap_or(f64::NAN)
    }
    fn gradient(&self, q: &[f64], out: &mut [f64]) {
        if q.len() != self.coordinate_count() || out.len() != q.len() {
            out.fill(f64::NAN);
            return;
        }
        let slopes = self.slopes(q);
        let axial = &q[self.transverse_count()..];
        let mut gradient = vec![0.0; slopes.len()];
        let (transverse_gradient, axial_gradient) = out.split_at_mut(self.transverse_count());
        if self
            .model
            .discrete_gradient(
                &slopes,
                &slopes,
                axial,
                axial,
                &mut gradient,
                axial_gradient,
            )
            .is_err()
        {
            out.fill(f64::NAN);
            return;
        }
        for plane in 0..self.planes {
            for mode in 0..self.transverse_modes {
                transverse_gradient[plane * self.transverse_modes + mode] = self.slope_scales[mode]
                    * gradient[plane * (self.transverse_modes + 1) + mode + 1];
            }
        }
    }
}
impl ModalPotential for MaterialLongitudinalString {
    fn coordinate_count(&self) -> usize {
        self.transverse_count() + self.model.axial_modes().len()
    }
    fn validate_coordinates(&self, q: &[f64]) -> Result<(), ModalCouplingError> {
        if q.len() != self.coordinate_count() || q.iter().any(|v| !v.is_finite()) {
            return Err(owner(fs_nlmodal::NlModalError::Parameter {
                what: "finite complete coordinates",
            }));
        }
        let slopes = self.slopes(q);
        let mut slope_bound = 0.0;
        for mode in 1..=self.transverse_modes {
            let norm_squared: f64 = (0..self.planes)
                .map(|plane| slopes[plane * (self.transverse_modes + 1) + mode].powi(2))
                .sum();
            slope_bound += det::sqrt(norm_squared);
        }
        let (mean, variation) = self
            .model
            .strain_bounds(&slopes, &q[self.transverse_count()..])
            .map_err(owner)?;
        if !slope_bound.is_finite()
            || slope_bound > self.resolution.maximum_slope
            || mean.abs() + variation > self.resolution.maximum_strain
            || self.prestress_n + self.axial_rigidity_n * (mean - variation) <= 0.0
        {
            return Err(ModalCouplingError::Invalid {
                what: "longitudinal string exceeds its slope/strain envelope or loses local tension",
            });
        }
        Ok(())
    }
}

fn elastic_source(
    state: &ResolvedMaterialStatePoint,
    band: &mut Option<(f64, f64)>,
) -> Result<(), ImpactError> {
    for key in ["density", "young_modulus"] {
        let coefficient = state.property(key).ok_or_else(|| {
            invalid("longitudinal elasticity requires original density and Young modulus receipts")
        })?;
        if !coefficient.is_constant_at_fixed_temperature() {
            return Err(invalid(
                "longitudinal elastic coefficients must be constant at the frozen material temperature",
            ));
        }
        if let Some(next) = coefficient
            .angular_frequency_band_rad_s()
            .map_err(|e| ImpactError::Owner(e.to_string()))?
        {
            *band = Some(band.map_or(next, |old| (old.0.max(next.0), old.1.min(next.1))));
        }
    }
    Ok(())
}
impl MaterialStringSource {
    pub(super) fn longitudinal_material(&self, band: &mut Option<(f64, f64)>) -> Result<f64, ImpactError> {
        match self {
            Self::WithPolarizations(p) => p.source().longitudinal_material(band),
            Self::WithTorsion(t) => t.bending().longitudinal_material(band),
            Self::Wound(_) => Err(invalid("resolved winding axial waves are outside the quasi-static spring reduction")),
            Self::Homogeneous(s) => {
                elastic_source(s.material(), band)?;
                Ok(s.string().lin_density_kg_m)
            }
            Self::Concentric { section, .. } => {
                for layer in section.elastic().layers() {
                    elastic_source(&layer.material, band)?;
                }
                Ok(section.elastic().line_density_kg_m())
            }
            Self::ConcentricProny { section, .. } => {
                for layer in section.elastic().layers() {
                    elastic_source(&layer.material, band)?;
                }
                Ok(section.elastic().line_density_kg_m())
            }
        }
    }
}
impl ReciprocalStringBank {
    /// Retain explicit longitudinal waves for every original material string.
    ///
    /// Geometry and the original homogeneous/bonded material determine EA,
    /// solid axial mass and every axial frequency. Both transverse directions
    /// contribute to the SAME strain. Original body, twist, gas memory, material
    /// interfaces, pressure observers and nonzero initial states survive.
    /// Select connections/gas/prepared propagation first. Subsequent bows,
    /// compliant stops, schedules and PCM use the unchanged coupled owner.
    ///
    /// This replaces instantaneous geometric stretching, so the two selections
    /// are mutually exclusive. The retained axial modes are elastic/lossless;
    /// original transverse/torsion/bridge/body/fluid loss still dissipates energy.
    /// No fitted axial damping, automatic retuning, winding slip or moving
    /// endpoints are supplied. Refine mode selection and clock for convergence;
    /// source bands, mode ceilings, strain limits and nonlinear iteration can
    /// refuse rather than drop modes, clip strain or invent a missing property.
    pub fn with_longitudinal_waves(
        self,
        resolutions: Vec<LongitudinalStringResolution>,
        config: NonlinearModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, ImpactError> {
        self.with_longitudinal_waves_and_loss(resolutions, LongitudinalMaterialLoss::Elastic, config, gate)
    }

    /// Retain the same nonlinear axial/transverse coupling with explicit axial loss.
    /// Every solid layer supplies its own uniaxial viscosity; missing data refuses.
    /// The original exact linear owner integrates the derived modal damping and
    /// includes its dissipation in the same friction/body/fluid energy account.
    /// No output filtering, second damping step or extra state is introduced.
    ///
    /// Kelvin--Voigt here is the small-strain linear u_xt operator. The nonlinear
    /// elastic strain law is unchanged, but finite-slope viscous y_x*y_xt terms,
    /// creep, heat evolution and interlayer slip are not modeled. Source bands
    /// screen retained linear axial frequencies, not all nonlinear harmonics.
    pub fn with_longitudinal_waves_and_loss(
        self,
        resolutions: Vec<LongitudinalStringResolution>,
        loss: LongitudinalMaterialLoss,
        config: NonlinearModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, ImpactError> {
        self.with_longitudinal_waves_and_bearings(resolutions, loss, Vec::new(), config, gate)
    }

    /// Add source-bound INTERIOR axial bearings while preparing the original
    /// longitudinal modes. The full supported string may include afterlength;
    /// its fixed endpoints do not move. Every body row is derived from the
    /// actual material shell node, global string axis and contact arm:
    /// t dot (U + Theta cross r). The same signed row applies the reciprocal
    /// force, so axial strain energy can excite the original radiating body.
    ///
    /// Every interface must match the real outer string and all incident body
    /// surface states, frozen temperatures and complete coupled frequency band.
    /// Constant-loss interfaces use the existing prepared propagation owner.
    /// All original modes, source identities, material memory, gas, radiation
    /// and nonlinear coordinate indices survive. No new propagator, output
    /// gain, endpoint mass approximation or inferred joint coefficient exists.
    pub fn with_longitudinal_waves_and_bearings(
        mut self,
        resolutions: Vec<LongitudinalStringResolution>,
        loss: LongitudinalMaterialLoss,
        bearings: Vec<MaterialAxialBearing>,
        config: NonlinearModalConfig,
        gate: &CancelGate,
    ) -> Result<Self, ImpactError> {
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        config
            .validate()
            .map_err(|e| ImpactError::Owner(e.to_string()))?;
        if self.system.samples_rendered() != 0
            || self.system.has_nonlinear_potentials()
            || resolutions.len() != self.cards.len()
        {
            return Err(invalid(
                "longitudinal waves require one cold selection per original string before nonlinear binding",
            ));
        }
        let mut strings = Vec::with_capacity(resolutions.len());
        let mut components = Vec::with_capacity(resolutions.len());
        let mut bindings = Vec::with_capacity(resolutions.len());
        for (index, ((member, card), resolution)) in self
            .source
            .strings
            .iter()
            .zip(&self.cards)
            .zip(resolutions)
            .enumerate()
        {
            if gate.is_requested() {
                return Err(ImpactError::Cancelled);
            }
            resolution.validate()?;
            let mut band = None;
            let solid_mass = member.source.longitudinal_material(&mut band)?;
            let ea = member.source.axial_rigidity_n()?;
            let axial_material_loss = loss.bind(&member.source)?;
            let planes = if matches!(&member.source, MaterialStringSource::WithPolarizations(_)) {
                2
            } else {
                1
            };
            let model = LongitudinalString::new(
                KcStringParams {
                    length: card.length_m,
                    lin_density: solid_mass,
                    tension: card.tension_n,
                    ea,
                },
                card.mode_count,
                planes,
                &resolution.axial_modes,
            )
            .map_err(|e| ImpactError::Owner(e.to_string()))?;
            if model
                .axial_omegas()
                .iter()
                .any(|&omega| band.is_some_and(|(lo, hi)| omega < lo || omega > hi))
            {
                return Err(invalid(
                    "retained longitudinal frequency exceeds original material elasticity/density applicability",
                ));
            }
            let modes = model
                .axial_omegas()
                .iter()
                .map(|&omega| -> Result<ModalAcousticMode, ImpactError> {
                    let damping_ratio = match &axial_material_loss {
                        Some(law) => law.longitudinal_damping_ratio(omega)
                            .map_err(|e| ImpactError::Owner(e.to_string()))?,
                        None => 0.0, // Explicit Elastic selection; never missing data.
                    };
                    Ok(ModalAcousticMode {
                        angular_frequency_rad_s: omega,
                        damping_ratio,
                        pressure_per_modal_velocity: C64::ZERO,
                    })
                })
                .collect::<Result<Vec<_>, _>>()?;
            components.push(
                ModalAcousticTimeModel::try_new(
                    self.source.options.sample_rate_hz,
                    modes,
                    self.source.options.string_budget,
                )
                .map_err(|e| ImpactError::Owner(format!("longitudinal modal admission: {e:?}")))?,
            );
            let law = MaterialLongitudinalString {
                model,
                resolution,
                transverse_modes: card.mode_count,
                planes,
                slope_scales: (1..=card.mode_count)
                    .map(|n| {
                        det::sqrt(2.0 / (card.linear_density_kg_m * card.length_m))
                            * n as f64
                            * core::f64::consts::PI
                            / card.length_m
                    })
                    .collect(),
                solid_linear_density_kg_m: solid_mass,
                transverse_linear_density_kg_m: card.linear_density_kg_m,
                axial_rigidity_n: ea,
                prestress_n: card.tension_n,
                axial_component: self.system.components().len() + index,
                elastic_frequency_band_rad_s: band,
                axial_material_loss,
            };
            let transverse = law.transverse_count();
            if self.system.components()[index].modes().len() != member.source.modes()
                || transverse > member.source.modes()
            {
                return Err(invalid(
                    "longitudinal storage lost the original complete string layout",
                ));
            }
            let coordinates = (0..transverse)
                .map(|mode| ModalCoordinate {
                    component: index,
                    mode,
                })
                .chain(
                    (0..law.model.axial_modes().len()).map(|mode| ModalCoordinate {
                        component: law.axial_component,
                        mode,
                    }),
                )
                .collect();
            bindings.push(MappedModalPotentialBinding {
                coordinates,
                replaces_linear: (transverse..law.coordinate_count()).collect(),
                potential: Box::new(law.clone()),
            });
            strings.push(law);
        }
        let (connections, anchors) = bearings::lower(
            &self.source, &strings, self.system.components(), &bearings, gate)?;
        self.system = self
            .system
            .with_internal_components_and_connections(components, connections, &anchors, gate)
            .and_then(|system| system.with_mapped_nonlinear_potentials(bindings, config, gate))
            .map_err(|e| ImpactError::Owner(e.to_string()))?;
        bearings::validate_bands(&self.source, &strings, &bearings, &self.system, gate)?;
        self.source.validate_body_memory_frequencies(&self.system, gate)
            .map_err(|e| ImpactError::Owner(e.to_string()))?;
        if !bearings.is_empty() {
            self.source.validate_material_bridge_frequencies(&self.system, gate)
                .map_err(|e| ImpactError::Owner(e.to_string()))?;
        }
        self.source.longitudinal_waves = Some(StringBankLongitudinal { strings, config, bearings });
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        Ok(self)
    }
    /// Original strain and conservative tension bounds at the accepted state.
    pub fn longitudinal_observation(
        &self,
        string: usize,
    ) -> Result<LongitudinalStringObservation, ImpactError> {
        observe(&self.source, self.system.components(), string)
    }
}
impl MultiBowStringBank {
    pub fn longitudinal_observation(
        &self,
        string: usize,
    ) -> Result<LongitudinalStringObservation, ImpactError> {
        observe(self.source().bank(), self.system().components(), string)
    }
}
fn observe(
    source: &MaterialStringBankSource,
    components: &[ModalAcousticTimeModel],
    string: usize,
) -> Result<LongitudinalStringObservation, ImpactError> {
    let law = source
        .longitudinal_waves
        .as_ref()
        .and_then(|bank| bank.strings.get(string))
        .ok_or_else(|| invalid("string has no bound longitudinal-wave law"))?;
    let q: Vec<_> = components[string].states()[..law.transverse_count()]
        .iter()
        .chain(components[law.axial_component].states())
        .map(|s| s.displacement_m_sqrt_kg)
        .collect();
    law.observe(&q)
        .map_err(|e| ImpactError::Owner(e.to_string()))
}
