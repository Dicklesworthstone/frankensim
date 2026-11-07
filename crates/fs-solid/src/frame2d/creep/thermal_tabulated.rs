//! Original thermal hardening on the existing timed fiber-frame owner.
//!
//! Preparation and constitutive work stay in fs-material. This adapter exposes
//! reference-volume integration measures and converts the native thermal account
//! without another return map, source query, heat law or equilibrium solve.

use super::{
    Frame2d, FrameCreepError, FrameInelasticEnergy, SectionModel, poll, positive, require,
};
use crate::fiber::{FiberLaw, FiberState};
use crate::fiber::chaboche::thermal::{ThermalChabocheFiberLaw, ThermalChabocheFiberState};
use fs_blake3::ContentHash;
use fs_material::creep::CreepTolerances;
use fs_material::plastic::energy::{EnergyError, EnergyTolerances};
use fs_material::plastic::thermal::tabulated::axial::{
    UniaxialThermalTabulated, UniaxialThermalTabulatedState,
};

/// One original integration point, observed without consuming a thermal interval.
/// Records have the same member/section/fiber order as `with_creep_laws` and the
/// native work audit. The mechanical member mass is independently inspectable
/// so a caloric composition can check it against the sourced fiber inventories.
#[derive(Clone, Copy, Debug)]
pub struct FrameThermalTabulatedFiber<'a> {
    pub element: usize,
    pub section: usize,
    pub fiber: usize,
    pub reference_volume_m3: f64,
    pub member_reference_mass_kg: f64,
    pub axial_strain: f64,
    pub law: &'a UniaxialThermalTabulated,
    pub state: &'a UniaxialThermalTabulatedState,
}

/// A retained source-bound thermal law and its original accepted history.
/// Heterogeneous members keep each native family and exact branch count.
#[derive(Clone, Copy, Debug)]
pub enum FrameThermalFiberLaw<'a> {
    Tabulated {
        law: &'a UniaxialThermalTabulated,
        state: &'a UniaxialThermalTabulatedState,
    },
    Chaboche {
        law: &'a ThermalChabocheFiberLaw,
        state: &'a ThermalChabocheFiberState,
    },
}

impl FrameThermalFiberLaw<'_> {
    pub fn temperature_k(self) -> f64 {
        match self {
            Self::Tabulated { state, .. } => state.temperature_k(),
            Self::Chaboche { state, .. } => state.temperature_k(),
        }
    }
    pub fn old_temperature_k(self) -> f64 {
        match self {
            Self::Tabulated { law, .. } => law.old_temperature_k(),
            Self::Chaboche { law, .. } => law.old_temperature_k(),
        }
    }
    pub fn law_temperature_k(self) -> f64 {
        match self {
            Self::Tabulated { law, .. } => law.temperature_k(),
            Self::Chaboche { law, .. } => law.temperature_k(),
        }
    }
    pub fn density_kg_m3(self) -> f64 {
        match self {
            Self::Tabulated { law, .. } => law.density_kg_m3(),
            Self::Chaboche { law, .. } => law.density_kg_m3(),
        }
    }
    pub fn caloric_curve_identity(self) -> ContentHash {
        match self {
            Self::Tabulated { law, .. } => law.caloric_curve_identity(),
            Self::Chaboche { law, .. } => law.caloric_curve_identity(),
        }
    }
    /// Admit the next original-source interval without evolving any history.
    pub fn admit_successor(self, next: &FiberLaw) -> Result<(), FrameCreepError> {
        let next_old_temperature = match (self, next) {
            (Self::Tabulated { law, .. }, FiberLaw::ThermalTabulated(next)) => {
                law.admit_successor(next).map_err(FrameCreepError::Thermal)?;
                next.old_temperature_k()
            }
            (Self::Chaboche { law, .. }, FiberLaw::ThermalChaboche(next)) => {
                law.admit_successor(next).map_err(FrameCreepError::ThermalChaboche)?;
                next.old_temperature_k()
            }
            _ => return Err(FrameCreepError::Invalid(
                "thermal frame history cannot change its native material family",
            )),
        };
        require(
            next_old_temperature == self.temperature_k(),
            "thermal frame interval must start at the original accepted temperature",
        )
    }

    fn admit_strain(self, strain: f64) -> Result<(), FrameCreepError> {
        match self {
            Self::Tabulated { law, state } => law.tensor_strain(strain, state)
                .map(|_| ()).map_err(FrameCreepError::Thermal),
            Self::Chaboche { law, state } => law.tensor_strain(strain, state)
                .map(|_| ()).map_err(FrameCreepError::ThermalChaboche),
        }
    }
}

/// Original geometric and inertial inventory, shared by all admitted thermal
/// inelastic families. Reading this record consumes no physical interval.
#[derive(Clone, Copy, Debug)]
pub struct FrameThermalFiber<'a> {
    pub element: usize,
    pub section: usize,
    pub fiber: usize,
    pub reference_volume_m3: f64,
    pub member_reference_mass_kg: f64,
    pub axial_strain: f64,
    pub law: FrameThermalFiberLaw<'a>,
}

impl Frame2d {
    /// Read one accepted thermal plastic history and its prepared source interval.
    /// This never evaluates another physical interval to observe a loaded state.
    pub fn thermal_tabulated_fiber(
        &self,
        element: usize,
        section: usize,
        fiber: usize,
    ) -> Option<(&UniaxialThermalTabulated, &UniaxialThermalTabulatedState)> {
        let SectionModel::Fiber(section) = self.elements.get(element)?.sections.get(section)?
        else {
            return None;
        };
        match (&section.fibers.get(fiber)?.law, section.states.get(fiber)?) {
            (FiberLaw::ThermalTabulated(law), FiberState::ThermalTabulated(state)) => {
                Some((law, state))
            }
            _ => None,
        }
    }

    /// Inspect a complete thermal-table frame through the same original
    /// inventory used for heterogeneous thermal materials.
    pub fn thermal_tabulated_fibers(
        &self,
        maximum_fiber_points: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<FrameThermalTabulatedFiber<'_>>, FrameCreepError> {
        self.thermal_inelastic_fibers(maximum_fiber_points, &mut cancelled)?
            .into_iter()
            .map(|point| {
                poll(&mut cancelled)?;
                let FrameThermalFiberLaw::Tabulated { law, state } = point.law else {
                    return Err(FrameCreepError::Invalid(
                        "thermal table inventory cannot substitute another native fiber family",
                    ));
                };
                Ok(FrameThermalTabulatedFiber {
                    element: point.element,
                    section: point.section,
                    fiber: point.fiber,
                    reference_volume_m3: point.reference_volume_m3,
                    member_reference_mass_kg: point.member_reference_mass_kg,
                    axial_strain: point.axial_strain,
                    law,
                    state,
                })
            })
            .collect()
    }

    /// Inspect every original thermal plastic fiber before caloric composition.
    /// Mixed tabulated and cyclic materials retain their exact native histories.
    /// Missing families refuse; no physical volume or inertia is omitted.
    /// Work and allocation are bounded by the original fiber-point budget.
    pub fn thermal_inelastic_fibers(
        &self,
        maximum_fiber_points: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<FrameThermalFiber<'_>>, FrameCreepError> {
        poll(&mut cancelled)?;
        require(
            maximum_fiber_points > 0 && !self.elements.is_empty(),
            "thermal frame inventory needs nonempty members and a positive point budget",
        )?;
        let mut inventory = Vec::new();
        for (element, member) in self.elements.iter().enumerate() {
            let member_reference_mass_kg = member.mass_per_length * member.l0;
            require(
                !member.sections.is_empty()
                    && member.sections.len() == member.w.len()
                    && member.sections.len() == member.commit.d.len()
                    && positive(member.l0)
                    && member_reference_mass_kg.is_finite()
                    && member_reference_mass_kg >= 0.0,
                "invalid thermal member reference integration or mass",
            )?;
            for (section_index, section) in member.sections.iter().enumerate() {
                let SectionModel::Fiber(section) = section else {
                    return Err(FrameCreepError::Invalid(
                        "thermal frame inventory requires every section to retain native thermal fibers",
                    ));
                };
                require(
                    !section.fibers.is_empty() && section.fibers.len() == section.states.len(),
                    "thermal frame inventory requires matching nonempty fiber histories",
                )?;
                for (fiber_index, (fiber, state)) in
                    section.fibers.iter().zip(&section.states).enumerate()
                {
                    poll(&mut cancelled)?;
                    require(
                        inventory.len() < maximum_fiber_points,
                        "thermal frame inventory exceeds its fiber point budget",
                    )?;
                    let law = match (&fiber.law, state) {
                        (FiberLaw::ThermalTabulated(law), FiberState::ThermalTabulated(state)) => {
                            FrameThermalFiberLaw::Tabulated { law, state }
                        }
                        (FiberLaw::ThermalChaboche(law), FiberState::ThermalChaboche(state)) => {
                            FrameThermalFiberLaw::Chaboche { law, state }
                        }
                        _ => {
                            return Err(FrameCreepError::Invalid(
                                "thermal frame inventory cannot omit or substitute a native fiber family",
                            ));
                        }
                    };
                    let reference_volume_m3 = fiber.area * member.w[section_index] * member.l0;
                    let (extension, curvature) = member.commit.d[section_index];
                    let axial_strain = extension - fiber.y * curvature;
                    require(
                        positive(fiber.area)
                            && fiber.y.is_finite()
                            && positive(reference_volume_m3)
                            && axial_strain.is_finite(),
                        "thermal frame inventory requires finite accepted strain and positive reference volume",
                    )?;
                    law.admit_strain(axial_strain)?;
                    inventory.push(FrameThermalFiber {
                        element,
                        section: section_index,
                        fiber: fiber_index,
                        reference_volume_m3,
                        member_reference_mass_kg,
                        axial_strain,
                        law,
                    });
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(inventory)
    }
}

pub(super) fn audit(
    law: &UniaxialThermalTabulated,
    old: &UniaxialThermalTabulatedState,
    accepted: &UniaxialThermalTabulatedState,
    old_strain: f64,
    new_strain: f64,
    tolerances: CreepTolerances,
) -> Result<FrameInelasticEnergy, FrameCreepError> {
    let increment = law
        .advance_with_energy(
            old_strain,
            old,
            new_strain,
            EnergyTolerances {
                yield_pa: tolerances.stress_pa,
                energy_j_m3: tolerances.energy_j_m3,
            },
        )
        .map_err(FrameCreepError::Thermal)?;
    require(
        increment.state == *accepted,
        "accepted thermal fiber history disagrees with its native interval work account",
    )?;
    // Read old storage at its own coefficients. Subtracting a large thermal
    // coefficient change from new-modulus storage could erase the old inventory.
    let initial_stored_j_m3 = law
        .stored_energy(old_strain, old)
        .map_err(FrameCreepError::Thermal)?
        .total_j_m3();
    let thermal = increment.thermal;
    let native = thermal.mechanical.energy;
    let delta = accepted.mechanical.equivalent_plastic_strain()
        - old.mechanical.equivalent_plastic_strain();
    let plastic_dissipation_j_m3 = law.tensor_law().tabulated_law().backbone().yield_stress * delta;
    let overstress_dissipation_j_m3 = match law.tensor_law().rate_law() {
        Some(rate) => {
            rate.endpoint_overstress_pa(delta).map_err(|error| {
                FrameCreepError::Mechanical(super::SolidError::MaterialRefused {
                    what: error.to_string(),
                })
            })? * delta
        }
        None => 0.0,
    };
    let mut energy = FrameInelasticEnergy {
        initial_stored_j_m3,
        final_stored_j_m3: native.final_state.total_j_m3(),
        mechanical_work_j_m3: thermal.mechanical_work_j_m3,
        physical_dissipation_j_m3: native.physical_dissipation_j_m3,
        creep_dissipation_j_m3: 0.0,
        plastic_dissipation_j_m3,
        overstress_dissipation_j_m3,
        recovery_dissipation_j_m3: 0.0,
        numerical_dissipation_j_m3: native.numerical_dissipation_j_m3,
        thermoelastic_storage_change_j_m3: thermal.thermoelastic_storage_change_j_m3,
        reversible_heat_j_m3: thermal.reversible_heat_j_m3,
        residual_j_m3: 0.0,
    };
    let partition_residual = energy.physical_dissipation_j_m3
        - energy.plastic_dissipation_j_m3
        - energy.overstress_dissipation_j_m3;
    energy.residual_j_m3 = energy.mechanical_work_j_m3
        - (energy.final_stored_j_m3 - energy.initial_stored_j_m3)
        - energy.thermoelastic_storage_change_j_m3
        - energy.physical_dissipation_j_m3
        - energy.reversible_heat_j_m3
        - energy.numerical_dissipation_j_m3;
    for residual_j_m3 in [partition_residual, energy.residual_j_m3] {
        if !residual_j_m3.is_finite() || residual_j_m3.abs() > tolerances.energy_j_m3 {
            return Err(FrameCreepError::Plastic(EnergyError::EnergyResidual {
                residual_j_m3,
                tolerance_j_m3: tolerances.energy_j_m3,
            }));
        }
    }
    Ok(energy)
}
