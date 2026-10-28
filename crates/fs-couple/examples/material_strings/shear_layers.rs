//! Source-pack geometry for the original coupled material shear-layer owner.
//! This adapter supplies data only. Constitutive integration, accepted history,
//! physical loss and pressure remain in the shared instrument/contact solve.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::MultiBowStringBank;
use fs_material::state_point::{
    PLASTIC_SHEAR_VISCOSITY_PROPERTY, ResolvedShearLayer, SHEAR_LAYER_MODULUS_PROPERTY,
    SHEAR_LAYER_YIELD_PROPERTY,
};

/// Area is the complete nominal layer area on one authored bow track. A uniform
/// ribbon partitions this area once across its cells, retaining the thickness.
#[derive(Clone, Debug)]
pub struct Declaration {
    pub track: String,
    pub pack: String,
    pub area_m2: f64,
    pub thickness_m: f64,
    pub maximum_elastic_shear_strain: f64,
}
#[derive(Clone, Debug)]
pub struct Inputs {
    declarations: Vec<Declaration>,
}
impl Inputs {
    pub fn new(declarations: Vec<Declaration>) -> Result<Self> {
        if declarations.is_empty() || declarations.len() > 16 {
            return Err(
                "material shear layers require one declaration per bow track, at most 16".into(),
            );
        }
        for (i, d) in declarations.iter().enumerate() {
            if d.track.is_empty()
                || d.track.len() > 128
                || declarations[..i].iter().any(|other| other.track == d.track)
                || d.pack.is_empty()
                || d.pack.len() > 512
                || !Path::new(&d.pack)
                    .components()
                    .all(|part| matches!(part, std::path::Component::Normal(_)))
                || !d.area_m2.is_finite()
                || d.area_m2 <= 0.0
                || !d.thickness_m.is_finite()
                || d.thickness_m <= 0.0
                || !d.maximum_elastic_shear_strain.is_finite()
                || d.maximum_elastic_shear_strain <= 0.0
                || d.maximum_elastic_shear_strain > 0.1
            {
                return Err("shear layers need distinct tracks, relative pack paths, positive area/thickness and an elastic shear limit in (0,0.1]".into());
            }
        }
        Ok(Self { declarations })
    }
    pub fn validate_tracks(&self, schedule: &GestureSchedule) -> Result<()> {
        if self.declarations.len() != schedule.tracks().len()
            || schedule
                .tracks()
                .iter()
                .any(|track| !self.declarations.iter().any(|d| d.track == track.id))
        {
            return Err(
                "material shear layers require exactly one named layer for every bow track".into(),
            );
        }
        Ok(())
    }
    pub fn preflight(&self, spec: &input::Spec) -> Result<()> {
        self.validate_tracks(&spec.schedule)
    }
    /// Reuse the original bulk-pack catalog and shared byte budget. Resolve all
    /// declarations before attaching them in one complete immutable binding.
    /// Missing data never leaves an inert construction law active at playback.
    pub fn bind(
        &self,
        spec: &input::Spec,
        bank: MultiBowStringBank,
        catalog: &mut Catalog,
        loader: &mut impl FnMut(&str) -> Result<MaterialPack>,
        gate: &CancelGate,
    ) -> Result<MultiBowStringBank> {
        poll(gate)?;
        self.preflight(spec)?;
        let mut requirements = vec![
            Requirement::try_new(
                SHEAR_LAYER_MODULUS_PROPERTY,
                Pressure::DIMS,
                Admit::StrictlyPositive,
            )?,
            Requirement::try_new(
                SHEAR_LAYER_YIELD_PROPERTY,
                Pressure::DIMS,
                Admit::NonNegative,
            )?,
            Requirement::try_new(
                PLASTIC_SHEAR_VISCOSITY_PROPERTY,
                DynViscosity::DIMS,
                Admit::StrictlyPositive,
            )?,
        ];
        if spec.contact_thermal {
            requirements.push(Requirement::try_new("density", Density::DIMS, Admit::StrictlyPositive)?);
        }
        let mut bindings = Vec::with_capacity(spec.bow_contact_count());
        for (port, track) in spec.schedule.tracks().iter().enumerate() {
            poll(gate)?;
            let declaration = self
                .declarations
                .iter()
                .find(|d| d.track == track.id)
                .ok_or("missing material shear-layer declaration")?;
            let state =
                resolve_properties(&declaration.pack, catalog, loader, spec, &requirements)?;
            let cells = spec.bow_ribbon.map_or(1, |(_, cells)| cells);
            let source = ResolvedShearLayer::new(
                state,
                declaration.area_m2 / cells as f64,
                declaration.thickness_m,
                declaration.maximum_elastic_shear_strain,
            )?;
            source.require_angular_frequency_rad_s(core::f64::consts::PI * f64::from(spec.rate))?;
            if spec.bow_ribbon.is_some() {
                for cell in 0..cells {
                    bindings.push((cell, source.clone()));
                }
            } else {
                bindings.push((port, source));
            }
        }
        poll(gate)?;
        Ok(bank.with_bow_shear_layers(bindings, gate)?)
    }
}
