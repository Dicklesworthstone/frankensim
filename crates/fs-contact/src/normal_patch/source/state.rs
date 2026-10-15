//! Re-resolve the existing point-contact ladder at one complete physical state.
//! This composes original bulk, ordered interface, model and geometry owners;
//! it contains no constitutive formulas or material-name defaults.
use fs_matdb::{InterfaceSystemCard,MaterialCard,QueryPoint};
use fs_material::state_point::{MaterialPropertySelection,MaterialStatePointError,
    ScalarAdmissibility,ScalarPropertyRequirement,resolve_interface_state_point,
    resolve_isotropic_solid_state_point};
use fs_tribo::InputAuthority;
use crate::interface_binding::{InterfaceBindingError,NormalContactModelSelection,
    ADHESION_ENERGY_PROPERTY,ADHESION_ENERGY_DIMS,bind_dry_interface_system_card,
    bind_isotropic_elastic_interface,bind_normal_interface_state,bind_normal_contact_model};
use super::{SourcePointContact,PointContactGeometry,PointContactEnvelope,NormalPatchError};

/// Immutable data and explicit source selections, reusable across state points.
/// Both bulk materials and the interface use the same declared physical point,
/// as required by the existing isotropic contact binding. This does not average
/// different temperatures or infer missing pressure, frequency or history axes.
pub struct PointContactMaterialCards<'a> {
    pub surface_a:&'a MaterialCard,
    pub surface_b:&'a MaterialCard,
    pub interface:&'a InterfaceSystemCard,
    pub bulk_selection:[MaterialPropertySelection;2],
    pub interface_selection:MaterialPropertySelection,
    pub model_selection:NormalContactModelSelection,
    pub authority:InputAuthority,
}

/// Preserve which original owner refused a complete contact state.
#[derive(Debug)]
pub enum ContactMaterialResolveError {
    Material(MaterialStatePointError),
    Interface(InterfaceBindingError),
    Contact(NormalPatchError),
}
impl core::fmt::Display for ContactMaterialResolveError {
    fn fmt(&self,f:&mut core::fmt::Formatter<'_>)->core::fmt::Result {
        match self {
            Self::Material(e)=>write!(f,"contact material: {e}"),
            Self::Interface(e)=>write!(f,"contact interface: {e}"),
            Self::Contact(e)=>write!(f,"contact geometry/envelope: {e}"),
        }
    }
}
impl std::error::Error for ContactMaterialResolveError {}

impl PointContactMaterialCards<'_> {
    /// Resolve complete sourced elasticity, yield, interface state and the
    /// selected normal model, then compile through the original Hertz evaluator.
    /// Geometry and operating bounds remain explicit. Each returned source owns
    /// its original usage receipts and model; no previous source is mutated.
    pub fn resolve(&self,point:&QueryPoint,geometry:PointContactGeometry,envelope:PointContactEnvelope)
        ->Result<SourcePointContact,ContactMaterialResolveError> {
        use ContactMaterialResolveError::{Material,Interface,Contact};
        let a=resolve_isotropic_solid_state_point(self.surface_a,point,self.bulk_selection[0].clone()).map_err(Material)?;
        let b=resolve_isotropic_solid_state_point(self.surface_b,point,self.bulk_selection[1].clone()).map_err(Material)?;
        let dry=bind_dry_interface_system_card(self.interface,self.authority).map_err(Interface)?;
        let elastic=bind_isotropic_elastic_interface(&dry,&a,&b).map_err(Interface)?;
        let requirement=ScalarPropertyRequirement::try_new(ADHESION_ENERGY_PROPERTY,
            ADHESION_ENERGY_DIMS,ScalarAdmissibility::NonNegative).map_err(Material)?;
        let interface=resolve_interface_state_point(self.interface,point,&[requirement],
            self.interface_selection.clone()).map_err(Material)?;
        let normal=bind_normal_interface_state(&elastic,&interface).map_err(Interface)?;
        let model=bind_normal_contact_model(&normal,self.model_selection).map_err(Interface)?;
        SourcePointContact::new(model,geometry,envelope).map_err(Contact)
    }
}
