//! Source-bound coherent absorption in perfectly bonded viscoelastic layers.
//!
//! Each physical solid keeps its original bulk/shear Maxwell laws, density,
//! thickness, material identity and strain/frequency support. A normal wave
//! crosses each bonded interface with continuous traction and normal velocity.
//! Positive bulk and shear strain-work integrals locate the heat in its actual
//! material, including weak loss and opaque shielding. No averaged loss law or
//! effective homogeneous specimen replaces the heterogeneous stack.
//!
//! Preparation and a complete drive take O(n) work; preparation retains O(n)
//! coefficients. Repeated drives, per-layer heat queries and local field queries
//! allocate nothing and never query the original material or exterior sources.
//! These remain frozen, small-strain, isotropic, planar, steady harmonic models
//! at normal incidence: no oblique P/S conversion, interface slip/compliance,
//! finite aperture, poroelastic frame or thermal evolution is inferred.

use super::{
    MAX_VISCOELASTIC_SLAB_BRANCHES, ViscoelasticSlabError, ViscoelasticSlabField,
    ViscoelasticSlabSource, kernel, poll, source_moduli,
};
use crate::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_blake3::ContentHash;
use fs_math::{c64::C64, det};

/// Hard bound on the number of original solids in one prepared stack.
pub const MAX_VISCOELASTIC_LAYERS: usize = 256;

/// One physical layer in left-to-right order, with no inferred interface law.
#[derive(Clone, Copy, Debug)]
pub struct ViscoelasticLayer<'a> {
    /// Original solid laws, density, identity and applicability.
    pub source: ViscoelasticSlabSource<'a>,
    /// Nonnegative physical thickness [m]. A zero layer has no interior or heat.
    pub thickness_m: f64,
}

/// Explicit preparation work and storage allowance; no hidden growth.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ViscoelasticLayerBudget {
    /// At most 256 original layers, with at least one required in the stack.
    pub maximum_layers: usize,
    /// Total bulk plus shear branches across the stack. Every individual layer
    /// also retains the original slab's hard limit of 64 branches.
    pub maximum_total_branches: usize,
}

/// Independent physical powers over the common authored solid-stack area.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ViscoelasticLayersPower {
    /// Incoming cycle-mean watts at the physical left and right exterior faces.
    pub incident_power_w: [f64; 2],
    /// Outgoing cycle-mean watts at those same two exterior faces.
    pub outgoing_power_w: [f64; 2],
    /// Sum of positive bulk strain-work integrals in the original solids [W].
    pub bulk_absorbed_power_w: f64,
    /// Sum of positive shear strain-work integrals in the original solids [W].
    pub shear_absorbed_power_w: f64,
    /// Incoming minus outgoing minus independently integrated local heat [W].
    pub balance_residual_w: f64,
}

/// Local material contribution, retaining the actual source that heats.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ViscoelasticLayerAbsorption {
    /// Original authored material-state identity of this physical layer.
    pub source_identity: ContentHash,
    /// Positive bulk strain-work integral over this layer [W].
    pub bulk_absorbed_power_w: f64,
    /// Positive shear strain-work integral over this layer [W].
    pub shear_absorbed_power_w: f64,
    /// Conservative all-depth peak-strain bound admitted by this layer's law.
    pub maximum_absolute_strain_bound: f64,
}

/// Prepared normal-incidence stack, retaining all original material owners.
pub struct ViscoelasticLayers<'a> {
    layers: Vec<ViscoelasticLayer<'a>>,
    exterior_sources: [&'a dyn FluidAcoustics; 2],
    exterior_states: [FluidAcousticState; 2],
    angular_frequency_rad_s: f64,
    wave: kernel::layers::Stack,
}

/// Immutable coherent drive of one prepared stack, with no per-drive allocation.
pub struct DrivenViscoelasticLayers<'s, 'a> {
    stack: &'s ViscoelasticLayers<'a>,
    input: [C64; 2],
    area_m2: f64,
    power: kernel::layers::Power,
}

fn error(error: kernel::layers::Error) -> ViscoelasticSlabError {
    match error {
        kernel::layers::Error::Invalid(what) => ViscoelasticSlabError::Invalid(what),
        kernel::layers::Error::Cancelled => ViscoelasticSlabError::Cancelled,
        kernel::layers::Error::OutsideSourceSupport => ViscoelasticSlabError::OutsideSourceSupport,
    }
}

impl<'a> ViscoelasticLayers<'a> {
    /// Prepare every original source once and compose its coherent scattering.
    /// The two exterior sources are queried once each, or once in total if the
    /// references are identical. Zero-thickness layers still validate their
    /// original source. Preparation is atomic on any refusal or cancellation.
    pub fn try_new(
        layers: &[ViscoelasticLayer<'a>],
        exterior_sources: [&'a dyn FluidAcoustics; 2],
        angular_frequency_rad_s: f64,
        budget: ViscoelasticLayerBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ViscoelasticSlabError> {
        poll(cancelled)?;
        if budget.maximum_layers > MAX_VISCOELASTIC_LAYERS
            || budget.maximum_total_branches
                > MAX_VISCOELASTIC_LAYERS * MAX_VISCOELASTIC_SLAB_BRANCHES
            || layers.len() > budget.maximum_layers
        {
            return Err(ViscoelasticSlabError::Budget);
        }
        if layers.is_empty() {
            return Err(ViscoelasticSlabError::Invalid(
                "a solid stack needs at least one layer",
            ));
        }
        let mut total_branches = 0_usize;
        for layer in layers {
            poll(cancelled)?;
            total_branches = total_branches
                .checked_add(layer.source.bulk.terms.len())
                .and_then(|n| n.checked_add(layer.source.shear.terms.len()))
                .filter(|n| *n <= budget.maximum_total_branches)
                .ok_or(ViscoelasticSlabError::Budget)?;
        }
        let mut moduli = Vec::with_capacity(layers.len());
        let mut reference = Vec::with_capacity(layers.len());
        for layer in layers {
            poll(cancelled)?;
            let values = source_moduli(
                layer.source,
                angular_frequency_rad_s,
                MAX_VISCOELASTIC_SLAB_BRANCHES,
                cancelled,
            )?;
            let modulus = values[0] + values[1].scale(4.0 / 3.0);
            let impedance = modulus
                .sqrt()
                .scale(det::sqrt(layer.source.density_kg_m3))
                .abs();
            if !impedance.is_finite() || impedance <= 0.0 {
                return Err(ViscoelasticSlabError::Invalid(
                    "solid reference impedance is unrepresentable",
                ));
            }
            moduli.push(values);
            reference.push(impedance);
        }
        let left = exterior_sources[0].acoustic_state()?;
        let right = if core::ptr::eq(exterior_sources[0], exterior_sources[1]) {
            left
        } else {
            exterior_sources[1].acoustic_state()?
        };
        let mut waves = Vec::with_capacity(layers.len());
        let mut previous = left.characteristic_impedance_pa_s_m();
        for (i, layer) in layers.iter().enumerate() {
            poll(cancelled)?;
            // Balanced real wave coordinates at the exact interface. These
            // zero-distance ports are mathematical, not extra fluid layers.
            let next = if i + 1 == layers.len() {
                right.characteristic_impedance_pa_s_m()
            } else {
                det::sqrt(reference[i]) * det::sqrt(reference[i + 1])
            };
            let wave = kernel::Slab::new(
                layer.source.density_kg_m3,
                moduli[i],
                angular_frequency_rad_s,
                layer.thickness_m,
                [previous, next],
            )
            .map_err(ViscoelasticSlabError::Invalid)?;
            waves.push(kernel::layers::Layer {
                wave,
                maximum_strain: layer.source.maximum_absolute_strain,
            });
            previous = next;
        }
        let wave = kernel::layers::Stack::new(waves, cancelled).map_err(error)?;
        poll(cancelled)?;
        Ok(Self {
            layers: layers.to_vec(),
            exterior_sources,
            exterior_states: [left, right],
            angular_frequency_rad_s,
            wave,
        })
    }

    /// Number of original physical layers.
    pub fn layer_count(&self) -> usize {
        self.layers.len()
    }
    /// Original source references, density, support, identity and thickness.
    pub fn layer(&self, index: usize) -> Option<ViscoelasticLayer<'a>> {
        self.layers.get(index).copied()
    }
    /// Original prepared bulk and shear complex moduli for a physical layer.
    pub fn layer_moduli_pa(&self, index: usize) -> Option<[C64; 2]> {
        self.wave.layer(index).map(|wave| wave.moduli)
    }
    /// Original exterior references, without re-evaluation.
    pub const fn exterior_sources(&self) -> [&'a dyn FluidAcoustics; 2] {
        self.exterior_sources
    }
    /// Original frozen exterior states used at the physical outer boundaries.
    pub const fn exterior_states(&self) -> [FluidAcousticState; 2] {
        self.exterior_states
    }
    /// Common positive angular frequency [rad/s].
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.angular_frequency_rad_s
    }
    /// Row-major pressure scattering: outgoing [left,right] from incoming
    /// [left,right] at the two physical exterior faces.
    pub const fn pressure_scattering(&self) -> [C64; 4] {
        self.wave.scattering
    }

    /// Illuminate both physical faces coherently with peak pressure phasors.
    /// Both inputs cover the same positive area. Every original layer's strain
    /// ceiling and independent local energy balance must admit the full field.
    pub fn drive<'s>(
        &'s self,
        incident_peak_pressure_pa: [C64; 2],
        area_m2: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DrivenViscoelasticLayers<'s, 'a>, ViscoelasticSlabError> {
        let power = self
            .wave
            .drive(incident_peak_pressure_pa, area_m2, cancelled)
            .map_err(error)?;
        Ok(DrivenViscoelasticLayers {
            stack: self,
            input: incident_peak_pressure_pa,
            area_m2,
            power,
        })
    }
}

impl DrivenViscoelasticLayers<'_, '_> {
    /// Physical exterior powers and sums of original positive local work.
    pub const fn power(&self) -> ViscoelasticLayersPower {
        ViscoelasticLayersPower {
            incident_power_w: self.power.incident,
            outgoing_power_w: self.power.outgoing,
            bulk_absorbed_power_w: self.power.heat[0],
            shear_absorbed_power_w: self.power.heat[1],
            balance_residual_w: self.power.residual,
        }
    }
    /// Outgoing peak pressures at the original physical exterior reference faces.
    pub const fn outgoing_pressure_pa(&self) -> [C64; 2] {
        self.power.outgoing_pressure
    }
    /// Positive heat in one physical layer, with its original source identity.
    pub fn layer_absorption(
        &self,
        index: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ViscoelasticLayerAbsorption, ViscoelasticSlabError> {
        poll(cancelled)?;
        let local = self
            .stack
            .wave
            .layer_drive(index, self.input, self.area_m2)
            .map_err(error)?;
        poll(cancelled)?;
        Ok(ViscoelasticLayerAbsorption {
            source_identity: self.stack.layers[index].source.source_identity,
            bulk_absorbed_power_w: local.heat[0],
            shear_absorbed_power_w: local.heat[1],
            maximum_absolute_strain_bound: local.strain_bound,
        })
    }

    /// Positive bulk/shear heat over one physical depth interval [m].
    ///
    /// The exact pressure Gram integral is restricted inside the SAME original
    /// layer; this inserts no interface and uses no difference of cumulative
    /// heats. The full [0, thickness] interval retains the whole-layer result.
    /// A zero-width interval has zero heat. The returned strain bound still
    /// conservatively covers the entire physical layer.
    pub fn absorption_between(
        &self,
        index: usize,
        depths_m: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ViscoelasticLayerAbsorption, ViscoelasticSlabError> {
        poll(cancelled)?;
        let local = self
            .stack
            .wave
            .layer_drive(index, self.input, self.area_m2)
            .map_err(error)?;
        let wave = self.stack.wave.layer(index).expect("validated solid layer");
        let heat = wave
            .interval_heat(local.inward, local.parity, depths_m, self.area_m2)
            .map_err(ViscoelasticSlabError::Invalid)?;
        poll(cancelled)?;
        Ok(ViscoelasticLayerAbsorption {
            source_identity: self.stack.layers[index].source.source_identity,
            bulk_absorbed_power_w: heat[0],
            shear_absorbed_power_w: heat[1],
            maximum_absolute_strain_bound: local.strain_bound,
        })
    }

    /// Local depth is measured from this physical layer's left face, including
    /// both faces. Zero-thickness layers have no interior and refuse a query.
    /// Interface traction/velocity agree; constitutive strain and heat may jump.
    pub fn field(
        &self,
        index: usize,
        depth_m: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ViscoelasticSlabField, ViscoelasticSlabError> {
        poll(cancelled)?;
        let field = self
            .stack
            .wave
            .field(index, self.input, depth_m)
            .map_err(error)?;
        poll(cancelled)?;
        Ok(ViscoelasticSlabField {
            compressive_traction_pa: field.pressure,
            normal_velocity_m_s: field.velocity,
            longitudinal_strain: field.strain,
            bulk_heat_density_w_m3: field.heat_density[0],
            shear_heat_density_w_m3: field.heat_density[1],
        })
    }
}
