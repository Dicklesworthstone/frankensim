//! Uniform, independently queried solids in the original isotropic contact law.
//! Only the temperature coordinates may differ. No representative temperature,
//! material mixture, thermal stress or contact heat-transfer law is invented.
use super::*;
use fs_matdb::QueryPoint;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

/// Bind the two actual bulk states and an explicit interface-law state point.
///
/// Each solid supplies its own E(T), nu(T) and yield limit through the original
/// state-point resolver. The existing reduced-modulus owner combines them in
/// their ordered roles. `interface_point` defines the state for interface claims
/// and the normal model; it is NOT substituted for either bulk temperature.
///
/// All three complete coordinate sets and quantity schemas must match, and
/// every non-T coordinate must be bit-identical. T must be positive kelvin,
/// either explicitly dimensional/absolute or consistently legacy-untyped.
/// The common-temperature limit retains the original binding identity exactly.
/// Subsequent normal-model admission checks the model domain at all three points.
///
/// This is frozen, uniform-per-solid elasticity for a declared current geometry.
/// It does not derive thermal expansion/prestress, flash temperatures, a thermal
/// boundary law or continuous loaded thermomechanical work. The caller supplies
/// and justifies the interface-law temperature; there is no default average.
pub fn bind_isotropic_elastic_interface_at_temperatures(
    interface: &BoundDryInterfaceSystem,
    surface_a: &IsotropicSolidStatePoint,
    surface_b: &IsotropicSolidStatePoint,
    interface_point: &QueryPoint,
) -> Result<BoundIsotropicElasticInterface, InterfaceBindingError> {
    check_surface_roles(interface, surface_a, surface_b)?;
    let point: Vec<_> = interface_point.axes().iter().map(|(k,v)|(k.clone(),*v)).collect();
    let quantities = interface_point.axis_quantities();
    if quantities.get("T").is_some_and(|q| {
        *q != QuantitySpec::dimensional(Temperature::DIMS)
            && *q != QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature, ValueForm::Static))
    }) {
        return Err(InterfaceBindingError::StatePointMismatch);
    }
    for bulk in [surface_a, surface_b] {
        let p = bulk.resolved().query_point();
        if bulk.resolved().axis_quantities() != quantities || p.len() != point.len()
            || !point.iter().any(|(k,t)| k == "T" && t.is_finite() && *t > 0.0)
            || p.iter().zip(&point).any(|((k,a),(name,b))| {
                k != name || if k == "T" { !a.is_finite() || *a <= 0.0 }
                    else { a.to_bits() != b.to_bits() }
            })
        {
            return Err(InterfaceBindingError::StatePointMismatch);
        }
    }
    let mut bound = bind_elastic_properties(interface, surface_a, surface_b)?;
    if state_points_exact_eq(surface_a.resolved().query_point(), &point)
        && state_points_exact_eq(surface_b.resolved().query_point(), &point)
    {
        return Ok(bound);
    }
    // The original property-pair identity already contains both exact queried
    // states and schemas. Bind the independently selected interface coordinates
    // too; changing only interface T must never reuse an old model identity.
    let mut identity = DomainHasher::new("org.frankensim.fs-contact.temperature-interface.v1");
    identity.update(bound.identity.as_bytes());
    identity.update(&(point.len() as u64).to_le_bytes());
    for (axis, value) in &point {
        identity.update(&(axis.len() as u64).to_le_bytes());
        identity.update(axis.as_bytes());
        identity.update(&value.to_bits().to_le_bytes());
    }
    bound.state_point = point;
    bound.identity = identity.finalize();
    Ok(bound)
}
