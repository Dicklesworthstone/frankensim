//! Authored demonstration cards, resolved through the shared material owners.
//! These values are synthetic; they do not describe measured named substances.

use super::Result;
use fs_evidence::ValidityDomain;
use fs_ga::{Se3, So3, Vec3};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid::{ResolvedLiquidState, resolve_liquid_state};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, DynViscosity, QuantitySpec, Temperature};

pub(super) const SIDES_M: Vec3 = Vec3 {
    x: 2.0,
    y: 2.0,
    z: 1.0,
};

fn card(phase: &str, properties: Vec<(PropertyKey, f64)>) -> Result<MaterialCard> {
    let mut claims = ClaimSet::new();
    for (key, value) in properties {
        let dims = key.dims();
        claims.insert_claim(PropertyClaim {
            key,
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(200.0, value), (400.0, value)],
                dims,
            },
            validity: ValidityDomain::unconstrained().with("T", 200.0, 400.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "authored heterogeneous flotation demonstration".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })?;
    }
    Ok(MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-flotation".into(),
            phase: phase.into(),
            process: "authored data; not measured".into(),
            revision: 0,
        },
        claims,
        vec![],
    )?)
}

pub(super) fn resolve(
    left_density_kg_m3: f64,
    right_density_kg_m3: f64,
    liquid_density_kg_m3: f64,
) -> Result<(ResolvedBody, ResolvedLiquidState)> {
    let left = card(
        "solid",
        vec![(
            PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            left_density_kg_m3,
        )],
    )?;
    let right = card(
        "solid",
        vec![(
            PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            right_density_kg_m3,
        )],
    )?;
    let point = QueryPoint::new().with("T", 300.0)?;
    let half = VolumeMoments::cuboid(Vec3::new(1.0, 2.0, 1.0))?;
    // These disjoint regions exactly tile the displacement envelope. Density
    // changes move the real material COM and change the original inertia.
    let body = resolve_body(
        &[
            RegionBinding {
                key: "left-half",
                geometry: half
                    .transformed(Se3::from_parts(So3::identity(), Vec3::new(-0.5, 0.0, 0.0))?)?,
                material: &left,
                point: &point,
                selection: MaterialPropertySelection::SingleClaimOnly,
            },
            RegionBinding {
                key: "right-half",
                geometry: half
                    .transformed(Se3::from_parts(So3::identity(), Vec3::new(0.5, 0.0, 0.0))?)?,
                material: &right,
                point: &point,
                selection: MaterialPropertySelection::SingleClaimOnly,
            },
        ],
        2,
        || false,
    )?;

    // Heat/viscosity are explicitly authored requirements of the common liquid
    // transport snapshot; the hydrostatic owner consumes only its density.
    let conductivity = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::ThermalConductivity,
        ValueForm::Static,
    ));
    let liquid = card(
        "liquid",
        vec![
            (
                PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
                liquid_density_kg_m3,
            ),
            (
                PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS),
                4000.0,
            ),
            (
                PropertyKey::with_quantity(THERMAL_CONDUCTIVITY_PROPERTY, conductivity),
                0.5,
            ),
            (
                PropertyKey::new("dynamic_viscosity", DynViscosity::DIMS),
                0.001,
            ),
        ],
    )?;
    Ok((
        body,
        resolve_liquid_state(&liquid, &point, MaterialPropertySelection::SingleClaimOnly)?,
    ))
}
