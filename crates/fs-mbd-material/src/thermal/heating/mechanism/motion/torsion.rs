//! Circular material shafts on the existing articulated and finite-heat timeline.
//!
//! Saint-Venant torsion and the original Kelvin–Voigt section supply
//! `k=sum(G_i J_i)/L`, `c=sum(eta_i J_i)/L`, torque `-k theta-c omega`,
//! and potential `k theta²/2`. Each layer receives its own `eta_i J_i` share
//! of integrated shear work. Reference geometry and coefficients freeze over
//! one existing thermal interval; endpoint stiffness changes are parameter work.
//! This reduction adds no distributed shaft inertia or second motion integrator.
use super::{Error, finite, member::MemberInventory, nonnegative, positive};
use crate::thermal::{ExpandedArticulated, heating::ThermalRegion};
use fs_matdb::{ClaimSelection, InterpolationPolicy, PropertyValue};
use fs_material::phase::EquilibriumPhaseState;
use fs_material::section::torsion::ConcentricTorsionalSection;
use fs_material::section_temperature::ConcentricSectionTemperatureLaw;
use fs_material::state_point::DENSITY_PROPERTY;
use fs_mbd::articulated::JointKind;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::collections::BTreeMap;

/// Lumped circular shaft compliance/loss on one actual revolute coordinate.
/// Its layer masses already belong to the original child/parent rigid links.
/// Axial/hoop stress, noncircular warping, slip and distributed torsional waves
/// require other owners. Thermal expansion of carrier inertia remains governed
/// by the original mechanism; reference section geometry stays fixed here.
#[derive(Clone, Debug)]
pub struct TorsionalJointShaft<'a> {
    pub link: usize,
    pub cells: Vec<CellId>,
    pub law: ConcentricSectionTemperatureLaw<'a>,
    /// Explicit solid shear modulus and Kelvin–Voigt shear viscosity roles.
    pub shear_modulus_property: String,
    pub shear_viscosity_property: String,
    pub reference_length_m: f64,
    pub reference_joint_angle_rad: f64,
    /// Outer-fibre engineering shear bound in (0, 0.1].
    pub maximum_absolute_shear_strain: f64,
    /// Declared complete represented/excited band [0, maximum] in rad/s.
    /// Every original coefficient must support it. This caller obligation is
    /// not inferred from a sampled trajectory or one natural frequency.
    pub maximum_angular_frequency_rad_s: f64,
    /// Accuracy scale of this shaft's integrated physical loss coordinate [J].
    pub absolute_loss_tolerance_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct TorsionalShaftState {
    pub link: usize,
    pub cells: Vec<CellId>,
    pub phases: Vec<EquilibriumPhaseState>,
    pub section: ConcentricTorsionalSection,
    pub twist_rad: f64,
    pub stiffness_nm_rad: f64,
    pub viscous_resistance_nm_s: f64,
    pub generalized_torque_nm: f64,
    pub elastic_energy_j: f64,
    pub dissipated_power_w: f64,
}

#[derive(Clone)]
pub(super) struct BoundTorsion<'a> {
    pub definition: TorsionalJointShaft<'a>,
    pub dof: usize,
    pub state: TorsionalShaftState,
    inventory: MemberInventory<'a>,
}
impl<'a> BoundTorsion<'a> {
    pub fn new(
        definition: TorsionalJointShaft<'a>,
        configuration: &ExpandedArticulated<'_>,
        regions: &[Vec<ThermalRegion<'a>>],
        thermal: &DepositionSimulation<'_>,
        position: &[f64],
        velocity: &[f64],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, Error> {
        if !definition.reference_joint_angle_rad.is_finite()
            || !positive(definition.maximum_absolute_shear_strain)
            || definition.maximum_absolute_shear_strain > 0.1
            || !positive(definition.maximum_angular_frequency_rad_s)
            || !positive(definition.absolute_loss_tolerance_j)
        {
            return Err(Error::Invalid(
                "shaft needs finite reference angle, small shear, frequency and loss bounds",
            ));
        }
        let inventory = MemberInventory::new(
            definition.link,
            JointKind::Revolute,
            &definition.cells,
            definition.law.reference(),
            definition.reference_length_m,
            configuration,
            regions,
            thermal,
            checkpoint,
        )?;
        let phases = inventory.phases(thermal, checkpoint)?;
        let section = resolve(&definition, &inventory, &phases, checkpoint)?;
        let state = TorsionalShaftState {
            link: definition.link,
            cells: definition.cells.clone(),
            phases,
            stiffness_nm_rad: finite(section.stiffness_n_m2() / definition.reference_length_m)?,
            viscous_resistance_nm_s: finite(
                section.viscous_stiffness_n_m2_s() / definition.reference_length_m,
            )?,
            section,
            twist_rad: 0.0,
            generalized_torque_nm: 0.0,
            elastic_energy_j: 0.0,
            dissipated_power_w: 0.0,
        };
        let mut result = Self {
            dof: inventory.dof,
            definition,
            state,
            inventory,
        };
        result.update(position[result.dof], velocity[result.dof])?;
        checkpoint()?;
        Ok(result)
    }

    /// Allocation-free effort, physical dissipation and the same potential.
    pub fn response(&self, q: f64, rate: f64) -> Result<(f64, f64, f64), Error> {
        let twist = finite(q - self.definition.reference_joint_angle_rad)?;
        let curvature = finite(twist / self.definition.reference_length_m)?;
        let shear = finite(self.state.section.shear().elastic().outer_radius_m() * curvature)?;
        if shear.abs() > self.definition.maximum_absolute_shear_strain
            || twist != 0.0 && (curvature == 0.0 || shear == 0.0)
            || !positive(self.state.stiffness_nm_rad)
            || !nonnegative(self.state.viscous_resistance_nm_s)
        {
            return Err(Error::Invalid(
                "shaft left its small-strain or representable constitutive domain",
            ));
        }
        let elastic = finite(self.state.stiffness_nm_rad * twist)?;
        let viscous = finite(self.state.viscous_resistance_nm_s * rate)?;
        let energy = finite(0.5 * elastic * twist)?;
        let power = finite(viscous * rate)?;
        if twist != 0.0 && energy == 0.0 {
            return Err(Error::Invalid("shaft elastic energy underflow"));
        }
        Ok((finite(-elastic - viscous)?, power, energy))
    }

    pub fn deposit_loss(
        &self,
        loss_j: f64,
        dt_s: f64,
        heat_w: &mut BTreeMap<CellId, f64>,
    ) -> Result<(), Error> {
        if !nonnegative(loss_j) || !positive(dt_s) {
            return Err(Error::Budget(
                "shaft integrated loss must be nonnegative; refine RK tolerances",
            ));
        }
        let total = self.state.section.viscous_stiffness_n_m2_s();
        if total == 0.0 {
            return if loss_j == 0.0 {
                Ok(())
            } else {
                Err(Error::Invalid("lossless shaft cannot generate heat"))
            };
        }
        let mut deposited = 0.0;
        for (&cell, &coefficient) in self
            .definition
            .cells
            .iter()
            .zip(self.state.section.layer_viscous_stiffness_n_m2_s())
        {
            let heat = finite(loss_j * (coefficient / total))?;
            deposited = finite(deposited + heat)?;
            let watts = heat_w.entry(cell).or_insert(0.0);
            *watts = finite(*watts + heat / dt_s)?;
        }
        if (deposited - loss_j).abs() > 128.0 * f64::EPSILON * loss_j {
            return Err(Error::Invalid(
                "shaft layer heat does not close its original shear work",
            ));
        }
        Ok(())
    }

    pub fn refresh(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        position: &[f64],
        velocity: &[f64],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Error> {
        let phases = self.inventory.phases(thermal, checkpoint)?;
        let section = resolve(&self.definition, &self.inventory, &phases, checkpoint)?;
        self.state.stiffness_nm_rad =
            finite(section.stiffness_n_m2() / self.definition.reference_length_m)?;
        self.state.viscous_resistance_nm_s =
            finite(section.viscous_stiffness_n_m2_s() / self.definition.reference_length_m)?;
        self.state.phases = phases;
        self.state.section = section;
        self.update(position[self.dof], velocity[self.dof])?;
        checkpoint()?;
        Ok(())
    }
    fn update(&mut self, q: f64, rate: f64) -> Result<(), Error> {
        let (torque, power, energy) = self.response(q, rate)?;
        self.state.twist_rad = q - self.definition.reference_joint_angle_rad;
        self.state.generalized_torque_nm = torque;
        self.state.dissipated_power_w = power;
        self.state.elastic_energy_j = energy;
        Ok(())
    }
}

fn resolve(
    definition: &TorsionalJointShaft<'_>,
    inventory: &MemberInventory<'_>,
    phases: &[EquilibriumPhaseState],
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<ConcentricTorsionalSection, Error> {
    let temperatures = phases.iter().map(|p| p.temperature_k()).collect::<Vec<_>>();
    let section = definition
        .law
        .at_temperatures_k(&temperatures)?
        .with_shear_modulus(&definition.shear_modulus_property, 1.0)?
        .with_torsional_viscosity(&definition.shear_viscosity_property)?;
    let band = section.omega_band_rad_s();
    if band.0 != 0.0 || band.1 < definition.maximum_angular_frequency_rad_s {
        return Err(Error::Invalid(
            "shaft sources must support DC through its complete declared frequency band",
        ));
    }
    for (i, ((current, reference), region)) in section
        .shear()
        .elastic()
        .layers()
        .iter()
        .zip(definition.law.reference().layers())
        .zip(&inventory.regions)
        .enumerate()
    {
        checkpoint()?;
        let original = reference.material.query();
        let actual = current.material.query();
        if original
            .clone()
            .with("T", temperatures[i])
            .map_err(|_| Error::Invalid("shaft requires canonical absolute T"))?
            != actual
        {
            return Err(Error::Invalid(
                "shaft must evolve only its original absolute T coordinate",
            ));
        }
        let reference_temperature = original
            .axes()
            .get("T")
            .ok_or(Error::Invalid("shaft requires canonical absolute T"))?;
        let (lower, upper) = if temperatures[i] < *reference_temperature {
            (&actual, &original)
        } else {
            (&original, &actual)
        };
        let density = current
            .material
            .property(DENSITY_PROPERTY)
            .ok_or(Error::Invalid("shaft lost its reference density"))?;
        let claim = region
            .card
            .claims()
            .claim(density.answer().receipt.selected)
            .ok_or(Error::Invalid("shaft reference density source is absent"))?;
        if !matches!(claim.value, PropertyValue::Scalar { .. })
            || claim.interpolation != InterpolationPolicy::ConstantWithinValidity
        {
            return Err(Error::Invalid(
                "shaft requires an explicit scalar reference density",
            ));
        }
        for name in [
            &definition.shear_modulus_property,
            &definition.shear_viscosity_property,
        ] {
            let property = current
                .material
                .property(name)
                .ok_or(Error::Invalid("shaft lost an original shear coefficient"))?;
            let claim = region
                .card
                .claims()
                .claim(property.answer().receipt.selected)
                .ok_or(Error::Invalid("shaft shear source is absent"))?;
            let line = region
                .card
                .claims()
                .query_line_typed(
                    &claim.key,
                    "T",
                    lower,
                    upper,
                    ClaimSelection::Pinned(property.answer().receipt.selected),
                )
                .map_err(|_| {
                    Error::Invalid("shaft shear source does not support its complete thermal path")
                })?;
            if line.knots.iter().any(|(_, answer)| {
                let value = answer.evidence.value.value;
                if name == &definition.shear_modulus_property {
                    !positive(value)
                } else {
                    !nonnegative(value)
                }
            }) {
                return Err(Error::Invalid(
                    "shaft shear coefficients leave their physical domain along the thermal path",
                ));
            }
        }
    }
    Ok(section)
}
