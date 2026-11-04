//! Shared original modal point rows for native rigid and articulated receivers.

use super::{
    CoupledModalSystem, ModalAcousticTimeModel, ModalCouplingError, checked, invalid, poll,
};
use fs_contact::frictional_multicontact::{FrictionalContact, NormalVelocityLaw};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_mbd::Vec3;

/// One physical material point in a component's complete original modal basis.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalImpactPoint {
    pub component: usize,
    /// World position at zero modal displacement [m].
    pub reference_position_world_m: Vec3,
    /// World displacement per mass-normalized coordinate [1/sqrt(kg)].
    /// Every original mode is present, including explicit zero traces.
    pub shapes_world: Vec<Vec3>,
}

pub(super) struct ModalContactInput<'a> {
    pub id: u64,
    pub point: &'a ModalImpactPoint,
    pub normal_world: Vec3,
    pub tangent_world: Vec3,
    pub activation_tolerance_m: f64,
    pub friction_coefficient: f64,
    pub normal_law: NormalVelocityLaw,
}

/// The same signed modal rows form geometry, flow, mobility and impulse maps.
pub(super) struct ModalContactRows {
    pub axes: Vec<[Vec3; 3]>,
    pub columns: Vec<Vec<f64>>,
    pub gram: Vec<Dd>,
    pub layout: Vec<usize>,
    pub sources: Vec<FrictionalContact>,
}

impl ModalContactRows {
    /// Inputs are already in canonical contact-ID order. The receiving owner
    /// bounds preparation/storage and validates its mechanical endpoint first.
    pub fn new(
        system: &CoupledModalSystem,
        contacts: &[ModalContactInput<'_>],
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        let dimension = 3 * contacts.len();
        let mut axes = Vec::with_capacity(contacts.len());
        let mut columns = Vec::with_capacity(dimension);
        let mut sources = Vec::with_capacity(contacts.len());
        for (index, contact) in contacts.iter().enumerate() {
            poll(Some(gate))?;
            if (index > 0 && contacts[index - 1].id >= contact.id)
                || !contact.point.reference_position_world_m.is_finite()
                || !contact.activation_tolerance_m.is_finite()
                || contact.activation_tolerance_m < 0.0
                || !contact.friction_coefficient.is_finite()
                || contact.friction_coefficient < 0.0
                || !match contact.normal_law {
                    NormalVelocityLaw::NewtonRestitution { coefficient } => {
                        coefficient.is_finite() && (0.0..=1.0).contains(&coefficient)
                    }
                    NormalVelocityLaw::Target { speed } => speed.is_finite() && speed >= 0.0,
                }
            {
                return Err(invalid(
                    "modal contact identity, witness, friction or normal law",
                ));
            }
            let model = system
                .models
                .get(contact.point.component)
                .ok_or_else(|| invalid("modal material point is outside the original network"))?;
            if contact.point.shapes_world.len() != model.modes().len()
                || contact
                    .point
                    .shapes_world
                    .iter()
                    .any(|shape| !shape.is_finite())
            {
                return Err(invalid(
                    "world point traces must cover the complete original modal basis",
                ));
            }
            let normal = unit(contact.normal_world)?;
            let tangent = unit(
                contact
                    .tangent_world
                    .sub(normal.scale(checked(dot(normal, contact.tangent_world))?)),
            )?;
            let frame = [normal, tangent, normal.cross(tangent)];
            for axis in frame {
                let mut column = vec![0.0; system.mode_count()];
                for (mode, &shape) in contact.point.shapes_world.iter().enumerate() {
                    if mode % 64 == 0 {
                        poll(Some(gate))?;
                    }
                    column[system.offsets[contact.point.component] + mode] =
                        -checked(dot(axis, shape))?;
                }
                columns.push(column);
            }
            axes.push(frame);
            sources.push(FrictionalContact {
                id: contact.id,
                gap: 0.0,
                activation_tolerance: contact.activation_tolerance_m,
                relative_velocity: [0.0; 3],
                friction_coefficient: contact.friction_coefficient,
                normal_law: contact.normal_law,
            });
        }
        let mut gram = vec![Dd::ZERO; dimension * dimension];
        for row in 0..dimension {
            poll(Some(gate))?;
            for col in row..dimension {
                let mut value = Dd::ZERO;
                for (mode, (&a, &b)) in columns[row].iter().zip(&columns[col]).enumerate() {
                    if mode % 64 == 0 {
                        poll(Some(gate))?;
                    }
                    value = value + Dd::from_f64(a) * Dd::from_f64(b);
                }
                checked(value)?;
                gram[row * dimension + col] = value;
                gram[col * dimension + row] = value;
            }
        }
        poll(Some(gate))?;
        Ok(Self {
            axes,
            columns,
            gram,
            sources,
            layout: system
                .models
                .iter()
                .map(|model| model.modes().len())
                .collect(),
        })
    }

    pub fn matches(&self, system: &CoupledModalSystem) -> bool {
        self.layout.len() == system.models.len()
            && system
                .models
                .iter()
                .zip(&self.layout)
                .all(|(model, &count)| model.modes().len() == count)
    }

    pub fn trace(
        &self,
        index: usize,
        models: &[ModalAcousticTimeModel],
        velocity: bool,
        gate: &CancelGate,
    ) -> Result<[Dd; 3], ModalCouplingError> {
        let mut result = [Dd::ZERO; 3];
        for (mode, state) in models.iter().flat_map(|model| model.states()).enumerate() {
            if mode % 64 == 0 {
                poll(Some(gate))?;
            }
            let value = Dd::from_f64(if velocity {
                state.velocity_m_sqrt_kg_per_s
            } else {
                state.displacement_m_sqrt_kg
            });
            for (axis, sum) in result.iter_mut().enumerate() {
                *sum = *sum + Dd::from_f64(self.columns[3 * index + axis][mode]) * value;
            }
        }
        Ok(result)
    }
}

pub(super) fn components(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}

pub(super) fn dot(a: Vec3, b: Vec3) -> Dd {
    components(a)
        .into_iter()
        .zip(components(b))
        .fold(Dd::ZERO, |sum, (a, b)| {
            sum + Dd::from_f64(a) * Dd::from_f64(b)
        })
}

pub(super) fn product(a: f64, b: f64) -> Result<Dd, ModalCouplingError> {
    let value = Dd::from_f64(a) * Dd::from_f64(b);
    checked(value)?;
    if a != 0.0 && b != 0.0 && value.hi == 0.0 && value.lo == 0.0 {
        return Err(invalid("contact projection loses a nonzero momentum term"));
    }
    Ok(value)
}

fn unit(value: Vec3) -> Result<Vec3, ModalCouplingError> {
    if !value.is_finite() {
        return Err(invalid("finite contact frame required"));
    }
    let scale = value.x.abs().max(value.y.abs()).max(value.z.abs());
    if scale == 0.0 {
        return Err(invalid(
            "nonzero independent contact frame directions required",
        ));
    }
    let value = Vec3::new(value.x / scale, value.y / scale, value.z / scale);
    Ok(value.scale(value.x.hypot(value.y).hypot(value.z).recip()))
}
