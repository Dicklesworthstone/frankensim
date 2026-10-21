//! Ordered material data and original caloric/history ownership for dry sliding.
use super::*;
use fs_contact::interface_binding::{BoundDryInterfaceSystem, bind_dry_interface_system_card};
use fs_matdb::QueryPoint;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedInterfaceStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_interface_state_point,
};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec};
use fs_thermal::enthalpy::CellId;
use fs_tribo::{FrictionLaw, HeatPartition, InputAuthority, WorkLedger};

pub const STATIC_FRICTION_PROPERTY: &str = "static-friction-coefficient";
pub const KINETIC_FRICTION_PROPERTY: &str = "kinetic-friction-coefficient";
const SPEED_DIMS: Dims = Dims([1, 0, -1, 0, 0, 0]);
const LENGTH_DIMS: Dims = Dims([1, 0, 0, 0, 0, 0]);

/// A declared finite, dry nominal patch with independent surface temperatures.
/// The heat split is caller data, not an inferred effusivity or flash model.
#[derive(Clone)]
pub struct SlidingInterface<'a> {
    pub point: ThermalInterfacePoint<'a>,
    pub caloric: [&'a EquilibriumEnthalpyPhaseCurve; 2],
    /// Distinct receiving thermal cells in the ORIGINAL surface A/B order.
    pub cells: [CellId; 2],
    pub nominal_area_m2: f64,
    pub pressure_axis: &'a str,
    pub speed_axis: &'a str,
    /// Optional actual nonnegative sliding arc length since initialization [m].
    /// This is distinct from the signed vector slip retained by fs-tribo.
    pub accumulated_distance_axis: Option<&'a str>,
    /// Fixed caller tangent frame for both relative displacement components.
    pub tangent_frame_id: &'a str,
    pub selection: MaterialPropertySelection,
    pub partition: HeatPartition,
    /// Required when the partition exports work outside these finite bodies.
    /// The destination's temperature/dynamics is not modeled by this reduction.
    pub outside_destination: Option<&'a str>,
}

/// One immutable source binding and accepted original dry-law history.
#[derive(Clone)]
pub struct MaterialSlidingContact<'a> {
    pub(super) data: SlidingInterface<'a>,
    bound: BoundDryInterfaceSystem,
    pub(super) history: DryTangentialState,
    pub(super) work: WorkLedger,
    pub(super) distance_m: f64,
    pub(super) time_s: f64,
}

#[derive(Clone, Debug)]
pub struct SlidingMaterialSample {
    pub resolved: ResolvedInterfaceStatePoint,
    pub phases: [EquilibriumPhaseState; 2],
    pub law: FrictionLaw,
    pub pressure_pa: f64,
    pub speed_mps: f64,
    pub accumulated_distance_m: f64,
}

impl<'a> MaterialSlidingContact<'a> {
    /// Select friction sources at the actual initial temperatures and declared
    /// initial load/speed, then pin those exact claims for every later trial.
    /// Source-domain failures never select a replacement. The initial slip and
    /// dissipative work are zero; existing runs continue using the returned
    /// contact rather than reconstructing its material/history binding.
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        data: SlidingInterface<'a>,
        initial_normal_load_n: f64,
        initial_speed_mps: f64,
        authority: InputAuthority,
    ) -> Result<Self, SlidingError> {
        validate(&data)?;
        let bound = bind_dry_interface_system_card(data.point.interface, authority)
            .map_err(SlidingError::Binding)?;
        let initial = sample(
            &data,
            thermal,
            initial_normal_load_n,
            initial_speed_mps,
            0.0,
        )?;
        let pins = MaterialPropertySelection::PinnedByProperty(
            initial
                .resolved
                .properties()
                .iter()
                .map(|p| {
                    (
                        p.requirement().name().to_owned(),
                        p.answer().receipt.selected,
                    )
                })
                .collect(),
        );
        let history = initial.law.initial_state();
        Ok(Self {
            data: SlidingInterface {
                selection: pins,
                ..data
            },
            bound,
            history,
            work: WorkLedger::default(),
            distance_m: 0.0,
            time_s: thermal.time_s(),
        })
    }
    pub const fn binding(&self) -> &SlidingInterface<'a> {
        &self.data
    }
    pub const fn interface(&self) -> &BoundDryInterfaceSystem {
        &self.bound
    }
    pub const fn state(&self) -> &DryTangentialState {
        &self.history
    }
    pub const fn work(&self) -> &WorkLedger {
        &self.work
    }
    pub const fn accumulated_distance_m(&self) -> f64 {
        self.distance_m
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }

    pub(super) fn sample<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
        normal_n: f64,
        speed: f64,
        distance: f64,
    ) -> Result<SlidingMaterialSample, SlidingError<E>> {
        sample(&self.data, thermal, normal_n, speed, distance)
    }
}

fn axis<E>(point: &QueryPoint, name: &str, dims: Dims) -> Result<(), SlidingError<E>> {
    let pressure =
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
    if !point.axes().contains_key(name)
        || point.axis_quantities().get(name).is_some_and(|q| {
            *q != QuantitySpec::dimensional(dims) && !(dims == Pressure::DIMS && *q == pressure)
        })
    {
        return Err(SlidingError::Invalid(
            "sliding coordinates require existing, dimensionally correct source axes",
        ));
    }
    Ok(())
}

fn validate<E>(data: &SlidingInterface<'_>) -> Result<(), SlidingError<E>> {
    if !data.nominal_area_m2.is_finite()
        || data.nominal_area_m2 <= 0.0
        || data.cells[0] == data.cells[1]
        || data.tangent_frame_id.trim().is_empty()
        || data
            .outside_destination
            .is_some_and(|s| s.trim().is_empty())
        || data.partition.other() > 0.0 && data.outside_destination.is_none()
    {
        return Err(SlidingError::Invalid(
            "sliding needs positive area, distinct cells/frame and a named outside destination",
        ));
    }
    let coordinates: Vec<_> = data
        .point
        .temperature_axes
        .into_iter()
        .chain(data.point.liquid_fraction_axes.into_iter().flatten())
        .chain(data.point.mean_temperature_axis)
        .chain([data.pressure_axis, data.speed_axis])
        .chain(data.accumulated_distance_axis)
        .collect();
    for (i, name) in coordinates.iter().enumerate() {
        if name.trim().is_empty() || coordinates[..i].contains(name) {
            return Err(SlidingError::Invalid(
                "sliding, pressure and thermal coordinates must be distinct",
            ));
        }
    }
    axis(data.point.context, data.pressure_axis, Pressure::DIMS)?;
    axis(data.point.context, data.speed_axis, SPEED_DIMS)?;
    if let Some(name) = data.accumulated_distance_axis {
        axis(data.point.context, name, LENGTH_DIMS)?;
    }
    Ok(())
}

fn sample<E>(
    data: &SlidingInterface<'_>,
    thermal: &DepositionSimulation<'_>,
    normal_n: f64,
    speed: f64,
    distance: f64,
) -> Result<SlidingMaterialSample, SlidingError<E>> {
    if !normal_n.is_finite()
        || normal_n <= 0.0
        || !speed.is_finite()
        || speed < 0.0
        || !distance.is_finite()
        || distance < 0.0
    {
        return Err(SlidingError::Invalid(
            "dry sliding requires a positive held load and finite nonnegative speed/distance",
        ));
    }
    let phases = data
        .cells
        .into_iter()
        .zip(data.caloric)
        .map(|(id, curve)| {
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(SlidingError::Invalid("inactive sliding thermal cell"))?;
            if !thermal.uses_material_law(cell.material_id, curve) {
                return Err(SlidingError::Invalid(
                    "sliding must retain each original caloric object",
                ));
            }
            let phase = curve.state_at_specific_enthalpy(cell.specific_h_j_kg)?;
            if phase.liquid_mass_fraction() != 0.0 {
                return Err(SlidingError::Invalid(
                    "dry sliding does not admit liquid or partially melted surfaces",
                ));
            }
            Ok(phase)
        })
        .collect::<Result<Vec<_>, SlidingError<E>>>()?;
    let phases = [phases[0], phases[1]];
    let pressure = number(normal_n / data.nominal_area_m2)?;
    if pressure <= 0.0 {
        return Err(SlidingError::Invalid("positive nominal pressure underflow"));
    }
    let mut point = data.point.at(phases)?;
    for (name, value) in [(data.pressure_axis, pressure), (data.speed_axis, speed)]
        .into_iter()
        .chain(data.accumulated_distance_axis.map(|name| (name, distance)))
    {
        point = point
            .with(name, value)
            .map_err(|e| ThermalMaterialError::Coordinate(e.to_string()))?;
    }
    let requirements = [STATIC_FRICTION_PROPERTY, KINETIC_FRICTION_PROPERTY]
        .into_iter()
        .map(|name| {
            ScalarPropertyRequirement::try_new(name, Dims::NONE, ScalarAdmissibility::NonNegative)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let resolved = resolve_interface_state_point(
        data.point.interface,
        &point,
        &requirements,
        data.selection.clone(),
    )?;
    for property in resolved.properties() {
        let claim = data
            .point
            .interface
            .claims()
            .claim(property.answer().receipt.selected)
            .ok_or(SlidingError::Invalid(
                "selected friction claim is absent from its original interface",
            ))?;
        if !claim.value.axes().all(|(axis, dims)| {
            data.point.is_thermal_axis(axis, dims)
                || axis == data.pressure_axis && dims == Pressure::DIMS
                || axis == data.speed_axis && dims == SPEED_DIMS
                || Some(axis) == data.accumulated_distance_axis && dims == LENGTH_DIMS
        }) {
            return Err(SlidingError::Invalid(
                "friction source must follow declared physical state coordinates",
            ));
        }
        if property
            .angular_frequency_band_rad_s()?
            .is_some_and(|(minimum, _)| minimum > 0.0)
        {
            return Err(SlidingError::Invalid(
                "sliding coefficient source does not admit a static force response",
            ));
        }
    }
    let law = FrictionLaw::Coulomb {
        static_mu: resolved
            .property(STATIC_FRICTION_PROPERTY)
            .expect("resolved requirement")
            .value_si(),
        kinetic_mu: resolved
            .property(KINETIC_FRICTION_PROPERTY)
            .expect("resolved requirement")
            .value_si(),
    };
    law.validate()?;
    Ok(SlidingMaterialSample {
        resolved,
        phases,
        law,
        pressure_pa: pressure,
        speed_mps: speed,
        accumulated_distance_m: distance,
    })
}
