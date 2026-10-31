//! Original interface-source response at each contact's solved operating point.
//! The native series-contact owner still computes the only physical drop/heat.

use super::{ConductionCoordinate, MaterialCircuit, MaterialCircuitError, invalid, poll};
use crate::material_circuit::contacts::{ContactHeat, SeriesElectricalContact};
use fs_matdb::InterfaceSystemCard;
use fs_qty::{Dims, QuantitySpec};

/// Electrical orientation of the interface card's ordered first/second surfaces.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ContactDirection {
    /// Positive branch current enters the card's first surface.
    AlongBranch,
    /// Positive branch current enters the card's second surface.
    AgainstBranch,
}
impl ContactDirection {
    fn sign(self) -> f64 {
        match self { Self::AlongBranch => 1.0, Self::AgainstBranch => -1.0 }
    }
}

/// A coordinate owned by a zero-thickness electrical interface.
/// Signed quantities use the explicitly bound ordered-surface direction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ContactCoordinate {
    /// Individual contact voltage jump [V], excluding bulk and other contacts.
    Voltage,
    /// Magnitude of the individual contact voltage jump [V].
    VoltageMagnitude,
    /// Series contact current [A].
    Current,
    /// Magnitude of contact current [A].
    CurrentMagnitude,
    /// Current divided by this contact's physical conducting area [A/m2].
    CurrentDensity,
    /// Magnitude of contact current density [A/m2].
    CurrentDensityMagnitude,
}
impl ContactCoordinate {
    pub(crate) fn dimensions(self) -> Dims {
        match self {
            Self::Voltage | Self::VoltageMagnitude => Dims([2, 1, -3, 0, -1, 0]),
            Self::Current | Self::CurrentMagnitude => Dims([0, 0, 0, 0, 1, 0]),
            Self::CurrentDensity | Self::CurrentDensityMagnitude => Dims([-2, 0, 0, 0, 1, 0]),
        }
    }
    fn value(self, point: ContactConductionPoint) -> f64 {
        match self {
            Self::Voltage => point.voltage_v,
            Self::VoltageMagnitude => point.voltage_v.abs(),
            Self::Current => point.current_a,
            Self::CurrentMagnitude => point.current_a.abs(),
            Self::CurrentDensity => point.current_density_a_m2,
            Self::CurrentDensityMagnitude => point.current_density_a_m2.abs(),
        }
    }
}

/// Bind an existing, typed source coordinate; names carry no implicit meaning.
#[derive(Clone, Copy, Debug)]
pub struct ContactAxis<'a> {
    pub name: &'a str,
    pub coordinate: ContactCoordinate,
}

/// Original complete ordered-interface card and its electrical coordinates.
pub struct ContactResponse<'a> {
    /// Stable contact ID, not vector position or containing branch index.
    pub id: u64,
    pub card: &'a InterfaceSystemCard,
    pub direction: ContactDirection,
    /// At most one each of voltage, current and current density. Empty is valid
    /// only if every property and applicability coordinate is electrically fixed.
    pub axes: &'a [ContactAxis<'a>],
}

/// Actual midpoint contact state, oriented from source surface A to B.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactConductionPoint {
    pub id: u64,
    pub branch: usize,
    pub direction: ContactDirection,
    pub voltage_v: f64,
    pub current_a: f64,
    pub current_density_a_m2: f64,
}

/// Applied-versus-evaluated interface residual; it is never an extra heat source.
#[derive(Clone, Debug, PartialEq)]
pub struct ContactConductivitySample {
    pub point: ContactConductionPoint,
    /// Complete original interface bundle evaluated at the actual operating point.
    pub contact: SeriesElectricalContact,
    pub resistance_defect_ohm: f64,
    pub voltage_defect_v: f64,
    pub power_defect_w: f64,
}

pub(super) fn electrical(dims: Dims) -> bool {
    [ContactCoordinate::Voltage.dimensions(), ContactCoordinate::Current.dimensions(),
        ContactCoordinate::CurrentDensity.dimensions(), ConductionCoordinate::ElectricField.dimensions()]
        .contains(&dims)
}

pub(super) fn admit<'a, 'b>(
    circuit: &MaterialCircuit,
    contacts: &[SeriesElectricalContact],
    sources: Option<&'a [ContactResponse<'b>]>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<&'a ContactResponse<'b>>, MaterialCircuitError> {
    let Some(sources) = sources else { return Ok(Vec::new()); };
    if sources.len() != contacts.len() {
        return Err(invalid("every contact requires one original interface response"));
    }
    let mut ordered: Vec<_> = sources.iter().collect();
    ordered.sort_by_key(|binding| binding.id);
    for (binding, seed) in ordered.iter().zip(contacts) {
        poll(cancelled)?;
        let state = seed.source().resolved();
        if binding.id != seed.id() || binding.card.content_hash() != state.card_identity()
            || binding.axes.len() > 3 {
            return Err(invalid("contact response must retain its identity and original card"));
        }
        let original = circuit.contacts().iter().find(|old| old.id() == binding.id)
            .ok_or_else(|| invalid("contact response lost its original interface"))?;
        let replay = original.source().at_query_point(binding.card, &state.query())?;
        if replay.resolved().properties().iter().zip(state.properties()).any(|(a, b)|
            a.requirement() != b.requirement()
                || a.answer().receipt.selected != b.answer().receipt.selected
                || a.value_si().to_bits() != b.value_si().to_bits()) {
            return Err(invalid("contact seed is not an original-source refresh"));
        }
        for (index, axis) in binding.axes.iter().enumerate() {
            if state.state_coordinate(axis.name).is_none()
                || state.axis_quantities().get(axis.name)
                    != Some(&QuantitySpec::dimensional(axis.coordinate.dimensions()))
                || binding.axes[..index].iter().any(|other|
                    other.name == axis.name || other.coordinate.dimensions() == axis.coordinate.dimensions()) {
                return Err(invalid("contact response axes require unique existing SI voltage/current/density coordinates"));
            }
        }
        let covered = |name: &str, dims: Dims| !electrical(dims)
            || binding.axes.iter().any(|axis| axis.name == name && axis.coordinate.dimensions() == dims);
        if state.axis_quantities().iter().any(|(name, q)| !covered(name, q.dims())) {
            return Err(invalid("contact electrical coordinates cannot remain frozen or imply an electric field"));
        }
        for property in state.properties() {
            let claim = binding.card.claims().claim(property.answer().receipt.selected)
                .ok_or_else(|| invalid("contact response lost its selected source"))?;
            if property.source_axes().iter().any(|(name, dims)| !covered(name, *dims))
                || claim.validity.axis_quantities().iter().any(|(name, q)| !covered(name, q.dims())) {
                return Err(invalid("every electrical source and applicability axis needs the contact closure"));
            }
        }
    }
    Ok(ordered)
}

pub(super) fn evaluate(
    original: &SeriesElectricalContact,
    binding: &ContactResponse<'_>,
    heat: &ContactHeat,
) -> Result<ContactConductivitySample, MaterialCircuitError> {
    if heat.contact.id() != original.id() || heat.contact.branch() != original.branch()
        || heat.contact.area_m2().to_bits() != original.area_m2().to_bits() {
        return Err(invalid("native contact response changed identity or geometry"));
    }
    let sign = binding.direction.sign();
    let point = ContactConductionPoint {
        id: original.id(), branch: original.branch(), direction: binding.direction,
        voltage_v: sign * heat.voltage_v, current_a: sign * heat.current_a,
        current_density_a_m2: sign * heat.current_a / original.area_m2(),
    };
    if !point.current_density_a_m2.is_finite()
        || point.current_a != 0.0 && point.current_density_a_m2 == 0.0 {
        return Err(invalid("contact current density is unrepresentable"));
    }
    let mut query = original.source().resolved().query();
    for axis in binding.axes {
        query = query.with(axis.name, axis.coordinate.value(point))
            .map_err(|_| invalid("contact electrical coordinate is unrepresentable"))?;
    }
    let contact = SeriesElectricalContact::new(original.id(), original.branch(),
        original.source().at_query_point(binding.card, &query)?, original.area_m2())?;
    let evaluated_voltage = point.current_a * contact.resistance_ohm();
    let evaluated_power = evaluated_voltage * point.current_a;
    let resistance_defect_ohm = contact.resistance_ohm() - heat.contact.resistance_ohm();
    let voltage_defect_v = point.voltage_v - evaluated_voltage;
    let power_defect_w = heat.power_w - evaluated_power;
    if ![evaluated_voltage, evaluated_power, resistance_defect_ohm, voltage_defect_v, power_defect_w]
        .iter().all(|v| v.is_finite())
        || point.current_a != 0.0 && contact.resistance_ohm() != 0.0
            && (evaluated_voltage == 0.0 || evaluated_power == 0.0) {
        return Err(invalid("contact constitutive residual is unrepresentable"));
    }
    Ok(ContactConductivitySample { point, contact, resistance_defect_ohm, voltage_defect_v, power_defect_w })
}
