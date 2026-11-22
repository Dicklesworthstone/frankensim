//! Authored lossless N-port relations, executed by the native component.

use fs_ports::{
    CoordinateBinding, JunctionLaw, JunctionSchema, PortOrientation, PortTimestamp, StableId,
};
use fs_system::elements::DiracJunction;

use super::{Domain, J, Refusal, coefficient_matrix, fields};

/// Bound dense scene preparation and interface solve work, as for PHS ports.
const MAX_PORTS: usize = 64;

/// All endpoints share the scene's scalar SI coordinates. The owner admits
/// exact skew symmetry and retains port ordering and the complete relation
/// in its snapshot identity. No matrix is repaired or inverted here.
pub(super) fn junction(definition: &J, name: &str, what: &str) -> Result<DiracJunction, Refusal> {
    fields(definition, &["kind", "domain", "ports", "incidence"], what)?;
    let domain = definition
        .str_field("domain")
        .and_then(Domain::parse)
        .ok_or_else(|| {
            format!("{what}: domain is one of translation, rotation, electrical, hydraulic")
        })?;
    let declared = definition
        .get("ports")
        .and_then(J::as_array)
        .filter(|ports| (2..=MAX_PORTS).contains(&ports.len()))
        .ok_or_else(|| format!("{what}.ports: name 2 to {MAX_PORTS} distinct scalar ports"))?;
    let identity = |value: &str| StableId::new(value).map_err(|error| format!("{what}: {error:?}"));
    let coordinates = CoordinateBinding::new(
        identity("system-scalar-si")?,
        identity("system-scene")?,
        PortOrientation::OutwardFromOwner,
    );
    let timestamp = PortTimestamp::new(identity("system-scene-clock")?, 0);
    let ports = declared
        .iter()
        .enumerate()
        .map(|(index, port)| {
            let label = port
                .as_str()
                .ok_or_else(|| format!("{what}.ports[{index}]: expected a port name"))?;
            // The scene addresses ports as element.port, splitting at the
            // last dot. A dotted port name could never be connected.
            if label.contains('.') {
                return Err(format!(
                    "{what}.ports[{index}]: port names cannot contain '.'"
                ));
            }
            let id = StableId::new(label).map_err(|_| {
                format!("{what}.ports[{index}]: {label:?} is not a port identifier")
            })?;
            domain
                .port_kind()
                .scalar_seed_schema(id, coordinates.clone(), timestamp.clone())
                .map_err(|error| format!("{what}.ports[{index}]: {error:?}"))
        })
        .collect::<Result<Vec<_>, String>>()?;
    let incidence = coefficient_matrix(
        definition.get("incidence"),
        ports.len(),
        ports.len(),
        &format!("{what}.incidence"),
    )?;
    let law = JunctionLaw::Dirac {
        rows: ports.len(),
        incidence,
    };
    let schema = JunctionSchema::try_new(identity(name)?, ports, law)
        .map_err(|error| format!("{what}: {error:?}"))?;
    DiracJunction::new(schema).map_err(|error| format!("{what}: {error}").into())
}
