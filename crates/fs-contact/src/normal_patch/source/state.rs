//! Re-resolve the existing point-contact ladder at complete physical states.
//! This composes original bulk, ordered interface, model and geometry owners;
//! it contains no constitutive formulas or material-name defaults.
use fs_matdb::{InterfaceSystemCard,MaterialCard,QueryPoint};
use fs_material::state_point::{IsotropicSolidStatePoint,MaterialPropertySelection,MaterialStatePointError,
    ScalarAdmissibility,ScalarPropertyRequirement,resolve_interface_state_point,
    resolve_isotropic_solid_state_point};
use fs_tribo::InputAuthority;
use crate::interface_binding::{InterfaceBindingError,NormalContactModelSelection,
    ADHESION_ENERGY_PROPERTY,ADHESION_ENERGY_DIMS,bind_dry_interface_system_card,
    bind_normal_interface_state,bind_normal_contact_model,
    bind_isotropic_elastic_interface_at_fixed_environments,
    temperature::bind_isotropic_elastic_interface_at_temperatures};
use super::{SourcePointContact,PointContactGeometry,PointContactEnvelope,NormalPatchError};

/// Immutable data and explicit source selections, reusable across state points.
/// `resolve` uses one common point. `resolve_at_temperatures` retains independent
/// bulk temperatures and a separately supplied interface-law temperature. Neither
/// averages bodies or infers missing pressure, frequency or history coordinates.
/// `resolve_at_fixed_environments` additionally admits independently held,
/// source-qualified bulk and interface conditions without copying their states.
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
    /// selected normal model at a common point through the original evaluator.
    pub fn resolve(&self,point:&QueryPoint,geometry:PointContactGeometry,envelope:PointContactEnvelope)
        ->Result<SourcePointContact,ContactMaterialResolveError> {
        self.resolve_at_temperatures([point,point],point,geometry,envelope)
    }

    /// Resolve each bulk at its own temperature, retaining the declared state
    /// for interface properties/model selection. Only T may differ: all other
    /// coordinates and axis schemas must match. Every original property and
    /// model validity gate remains active, including both bulk temperatures.
    /// No effective-temperature averaging or thermal-stress correction is used.
    pub fn resolve_at_temperatures(&self,bulk_points:[&QueryPoint;2],interface_point:&QueryPoint,
        geometry:PointContactGeometry,envelope:PointContactEnvelope)
        ->Result<SourcePointContact,ContactMaterialResolveError> {
        self.resolve_points(bulk_points,interface_point,geometry,envelope,false)
    }

    /// Resolve original bulk and interface properties at their independently
    /// held environments. Only source-qualified fixed external coordinates are
    /// admitted, with one common optional angular-frequency applicability point.
    /// Geometry, original model domains and normal-contact limits are unchanged.
    pub fn resolve_at_fixed_environments(&self,bulk_points:[&QueryPoint;2],interface_point:&QueryPoint,
        geometry:PointContactGeometry,envelope:PointContactEnvelope)
        ->Result<SourcePointContact,ContactMaterialResolveError> {
        self.resolve_points(bulk_points,interface_point,geometry,envelope,true)
    }

    /// Bind complete, already selected solid states, including their explicit
    /// elastic parameterizations and yield sources. Both must belong to these
    /// exact immutable bulk cards. Only temperature may differ between bodies;
    /// no bulk property is reselected or removed from either source bundle.
    pub fn resolve_from_bulk_states_at_temperatures(&self,bulk:[&IsotropicSolidStatePoint;2],
        interface_point:&QueryPoint,geometry:PointContactGeometry,envelope:PointContactEnvelope)
        ->Result<SourcePointContact,ContactMaterialResolveError> {
        self.resolve_states(bulk,interface_point,geometry,envelope,false)
    }

    /// Bind the original complete solids at independently held source-qualified
    /// environments through the same native contact owner. Original model,
    /// interface, geometry and yield admission remain mandatory.
    pub fn resolve_from_bulk_states_at_fixed_environments(&self,bulk:[&IsotropicSolidStatePoint;2],
        interface_point:&QueryPoint,geometry:PointContactGeometry,envelope:PointContactEnvelope)
        ->Result<SourcePointContact,ContactMaterialResolveError> {
        self.resolve_states(bulk,interface_point,geometry,envelope,true)
    }

    fn resolve_points(&self,bulk_points:[&QueryPoint;2],interface_point:&QueryPoint,
        geometry:PointContactGeometry,envelope:PointContactEnvelope,independent_environment:bool)
        ->Result<SourcePointContact,ContactMaterialResolveError> {
        use ContactMaterialResolveError::Material;
        let a=resolve_isotropic_solid_state_point(self.surface_a,bulk_points[0],self.bulk_selection[0].clone()).map_err(Material)?;
        let b=resolve_isotropic_solid_state_point(self.surface_b,bulk_points[1],self.bulk_selection[1].clone()).map_err(Material)?;
        self.resolve_states([&a,&b],interface_point,geometry,envelope,independent_environment)
    }

    fn resolve_states(&self,bulk:[&IsotropicSolidStatePoint;2],interface_point:&QueryPoint,
        geometry:PointContactGeometry,envelope:PointContactEnvelope,independent_environment:bool)
        ->Result<SourcePointContact,ContactMaterialResolveError> {
        use ContactMaterialResolveError::{Material,Interface,Contact};
        for (state,card) in bulk.into_iter().zip([self.surface_a,self.surface_b]) {
            if state.resolved().card_identity()!=card.content_hash() {
                return Err(Material(MaterialStatePointError::InvalidDerived {
                    quantity:"contact solid requires its exact configured immutable bulk card",
                }));
            }
        }
        let [a,b]=bulk;
        let dry=bind_dry_interface_system_card(self.interface,self.authority).map_err(Interface)?;
        let elastic=if independent_environment {
            bind_isotropic_elastic_interface_at_fixed_environments(&dry,a,b,interface_point)
        } else {
            bind_isotropic_elastic_interface_at_temperatures(&dry,a,b,interface_point)
        }.map_err(Interface)?;
        let requirement=ScalarPropertyRequirement::try_new(ADHESION_ENERGY_PROPERTY,
            ADHESION_ENERGY_DIMS,ScalarAdmissibility::NonNegative).map_err(Material)?;
        let interface=resolve_interface_state_point(self.interface,interface_point,&[requirement],
            self.interface_selection.clone()).map_err(Material)?;
        if independent_environment && interface.properties().iter().any(|property|
            !property.is_constant_at_fixed_environment()) {
            return Err(Interface(InterfaceBindingError::InvalidDerived {
                quantity: "fixed_environment_interface",
            }));
        }
        let normal=bind_normal_interface_state(&elastic,&interface).map_err(Interface)?;
        let model=bind_normal_contact_model(&normal,self.model_selection).map_err(Interface)?;
        SourcePointContact::new(model,geometry,envelope).map_err(Contact)
    }
}
