//! Permanent creep of native force-based fiber members over one physical interval.
//!
//! The existing section inversion, member compatibility, global Newton solve and
//! backtracking remain the mechanical owners. Every trial uses accepted fiber
//! history. Work is integrated with the original fiber areas, Lobatto weights
//! and reference lengths. Physical creep/plastic loss, thermoelastic storage,
//! and signed reciprocal heat remain distinct.
//!
//! This checked path admits linear geometry, stationary supports and Norton,
//! joint Norton/J2, or their sourced thermal-expansion wrappers. Temperature
//! is prescribed through original source-bound laws between intervals; this
//! owner reports reciprocal heat but does not itself solve a caloric field.
//! Nodal rotations and the sum of absolute quadrature curvature are bounded;
//! unsampled spatial extrema still require section/member refinement.

use super::{Frame2d, Geometry, SectionModel, SolidError, StaticStep};
use crate::fiber::{FiberLaw, FiberState};
use fs_material::creep::plastic::uniaxial::{UniaxialNortonJ2, UniaxialNortonJ2State};
use fs_material::creep::thermal::uniaxial::{
    UniaxialThermalNorton, UniaxialThermalNortonJ2, UniaxialThermalNortonJ2State,
    UniaxialThermalNortonState,
};
use fs_material::creep::uniaxial::{UniaxialNorton, UniaxialNortonState};
use fs_material::creep::{CreepEnergy, CreepError, CreepTolerances, NortonCreep};
use fs_material::state_point::ResolvedMaterialStatePoint;
use std::cell::Cell;

/// Independent physical-unit gates and work bound for one native frame step.
#[derive(Clone, Copy, Debug)]
pub struct FrameCreepTolerances {
    pub material: CreepTolerances,
    pub energy_j: f64,
    pub force_n: f64,
    pub moment_n_m: f64,
    /// Total fiber points over every member and original Lobatto section.
    pub maximum_fiber_points: usize,
}

#[derive(Debug)]
pub enum FrameCreepError {
    Invalid(&'static str),
    Mechanical(SolidError),
    Material(CreepError),
    Cancelled,
    Equilibrium {
        dof: usize,
        residual: f64,
        tolerance: f64,
    },
    Energy {
        residual_j: f64,
        tolerance_j: f64,
    },
}
impl core::fmt::Display for FrameCreepError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native frame creep refused: {self:?}")
    }
}
impl std::error::Error for FrameCreepError {}

/// One material-point account and its original reference integration measure.
#[derive(Clone, Debug)]
pub struct FrameFiberCreepEnergy {
    pub element: usize,
    pub section: usize,
    pub fiber: usize,
    pub reference_volume_m3: f64,
    pub energy: CreepEnergy,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct FrameCreepEnergy {
    pub initial_stored_j: f64,
    pub final_stored_j: f64,
    pub applied_work_j: f64,
    pub support_work_j: f64,
    pub creep_dissipation_j: f64,
    pub plastic_dissipation_j: f64,
    pub physical_dissipation_j: f64,
    /// Change in T times mechanical entropy, separate from free-energy storage.
    pub thermoelastic_storage_change_j: f64,
    /// Signed reciprocal transfer into a caloric owner; may be negative.
    pub reversible_heat_j: f64,
    pub numerical_dissipation_j: f64,
    pub equilibrium_work_residual_j: f64,
    pub balance_residual_j: f64,
}

impl FrameCreepEnergy {
    /// Irreversible material heat plus signed reciprocal thermal exchange.
    pub fn caloric_heat_j(self) -> f64 {
        self.physical_dissipation_j + self.reversible_heat_j
    }
}

/// Immutable proposal; accepting `frame` advances history exactly once.
#[derive(Clone, Debug)]
pub struct FrameCreepStep {
    pub frame: Frame2d,
    pub duration_s: f64,
    pub equilibrium: StaticStep,
    pub fibers: Vec<FrameFiberCreepEnergy>,
    pub energy: FrameCreepEnergy,
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn require(ok: bool, reason: &'static str) -> Result<(), FrameCreepError> {
    if ok {
        Ok(())
    } else {
        Err(FrameCreepError::Invalid(reason))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), FrameCreepError> {
    if cancelled() {
        Err(FrameCreepError::Cancelled)
    } else {
        Ok(())
    }
}
fn energy_gate(residual: f64, tolerance: f64) -> Result<(), FrameCreepError> {
    if residual.is_finite() && residual.abs() <= tolerance {
        Ok(())
    } else {
        Err(FrameCreepError::Energy {
            residual_j: residual,
            tolerance_j: tolerance,
        })
    }
}

/// Borrow the original law and history; this dispatch adds no constitutive solve.
#[derive(Clone, Copy)]
enum CreepPoint<'a> {
    Norton(&'a UniaxialNorton, &'a UniaxialNortonState),
    NortonJ2(&'a UniaxialNortonJ2, &'a UniaxialNortonJ2State),
    ThermalNorton(&'a UniaxialThermalNorton, &'a UniaxialThermalNortonState),
    ThermalNortonJ2(&'a UniaxialThermalNortonJ2, &'a UniaxialThermalNortonJ2State),
}
impl<'a> CreepPoint<'a> {
    fn new(law: &'a FiberLaw, state: &'a FiberState) -> Result<Self, FrameCreepError> {
        match (law, state) {
            (FiberLaw::Norton(law), FiberState::Norton(state)) => Ok(Self::Norton(law, state)),
            (FiberLaw::NortonJ2(law), FiberState::NortonJ2(state)) => {
                Ok(Self::NortonJ2(law, state))
            }
            (FiberLaw::ThermalNorton(law), FiberState::ThermalNorton(state)) => {
                Ok(Self::ThermalNorton(law, state))
            }
            (FiberLaw::ThermalNortonJ2(law), FiberState::ThermalNortonJ2(state)) => {
                Ok(Self::ThermalNortonJ2(law, state))
            }
            _ => Err(FrameCreepError::Invalid(
                "checked creep work requires matching timed fiber histories",
            )),
        }
    }
    fn kinetics(self) -> &'a NortonCreep {
        match self {
            Self::Norton(law, _) => law.backbone(),
            Self::NortonJ2(law, _) => law.backbone().creep(),
            Self::ThermalNorton(law, _) => law.kinetics(),
            Self::ThermalNortonJ2(law, _) => law.kinetics(),
        }
    }
    fn material(self) -> Option<&'a ResolvedMaterialStatePoint> {
        match self {
            Self::Norton(law, _) => law.material(),
            Self::NortonJ2(law, _) => law.material(),
            Self::ThermalNorton(law, _) => law.material(),
            Self::ThermalNortonJ2(law, _) => law.material(),
        }
    }
    fn expansion_material(self) -> Option<&'a ResolvedMaterialStatePoint> {
        match self {
            Self::ThermalNorton(law, _) => Some(law.expansion_material()),
            Self::ThermalNortonJ2(law, _) => Some(law.expansion_material()),
            Self::Norton(..) | Self::NortonJ2(..) => None,
        }
    }
    fn hardening(self) -> Option<f64> {
        match self {
            Self::Norton(..) | Self::ThermalNorton(..) => None,
            Self::NortonJ2(law, _) => Some(law.backbone().plastic().hardening),
            Self::ThermalNortonJ2(law, _) => Some(law.hardening()),
        }
    }
    fn yield_stress(self) -> Option<f64> {
        match self {
            Self::Norton(..) | Self::ThermalNorton(..) => None,
            Self::NortonJ2(law, _) => Some(law.backbone().plastic().yield_stress),
            Self::ThermalNortonJ2(law, _) => Some(law.yield_stress()),
        }
    }
    fn tensor_strain(self, strain: f64) -> Result<(), FrameCreepError> {
        match self {
            Self::Norton(law, state) => law.tensor_strain(strain, state),
            Self::NortonJ2(law, state) => law.tensor_strain(strain, state),
            Self::ThermalNorton(law, state) => law.tensor_strain(strain, state),
            Self::ThermalNortonJ2(law, state) => law.tensor_strain(strain, state),
        }
        .map(|_| ())
        .map_err(FrameCreepError::Material)
    }
    fn audit(
        self,
        committed: Self,
        old_strain: f64,
        new_strain: f64,
        tolerances: CreepTolerances,
    ) -> Result<CreepEnergy, FrameCreepError> {
        let modulus = self.kinetics().elastic().youngs;
        let close = |a: &[f64; 6], b: &[f64; 6], pa: f64, pb: f64| {
            a.iter()
                .zip(b)
                .all(|(a, b)| (a - b).abs() * modulus <= tolerances.stress_pa)
                && (pa - pb).abs() * modulus <= tolerances.stress_pa
        };
        let (energy, agrees) = match (self, committed) {
            (Self::Norton(law, state), Self::Norton(_, accepted)) => {
                let increment = law
                    .advance_with_energy(old_strain, state, new_strain, tolerances)
                    .map_err(FrameCreepError::Material)?;
                let history = accepted.tensor_history();
                let agrees = close(
                    &history.creep_strain,
                    &increment.state.creep_strain,
                    history.equivalent_creep_strain,
                    increment.state.equivalent_creep_strain,
                );
                (increment.energy, agrees)
            }
            (Self::NortonJ2(law, state), Self::NortonJ2(_, accepted)) => {
                let increment = law
                    .advance_with_energy(old_strain, state, new_strain, tolerances)
                    .map_err(FrameCreepError::Material)?;
                let history = accepted.tensor_history();
                let agrees = close(
                    &history.creep.creep_strain,
                    &increment.state.creep.creep_strain,
                    history.creep.equivalent_creep_strain,
                    increment.state.creep.equivalent_creep_strain,
                ) && close(
                    &history.plastic.plastic_strain,
                    &increment.state.plastic.plastic_strain,
                    history.plastic.alpha,
                    increment.state.plastic.alpha,
                );
                (increment.energy, agrees)
            }
            (Self::ThermalNorton(law, state), Self::ThermalNorton(_, accepted)) => {
                let increment = law
                    .advance_with_energy(old_strain, state, new_strain, tolerances)
                    .map_err(FrameCreepError::Material)?;
                (increment.energy, increment.state == *accepted)
            }
            (Self::ThermalNortonJ2(law, state), Self::ThermalNortonJ2(_, accepted)) => {
                let increment = law
                    .advance_with_energy(old_strain, state, new_strain, tolerances)
                    .map_err(FrameCreepError::Material)?;
                (increment.energy, increment.state == *accepted)
            }
            _ => {
                return Err(FrameCreepError::Invalid(
                    "timed fiber family changed during native solve",
                ));
            }
        };
        require(
            agrees,
            "published native history disagrees with the checked constitutive step",
        )?;
        Ok(energy)
    }
}

impl Frame2d {
    /// Existing static, modal and inertial routes do not own a creep clock.
    /// A duration-bearing fiber trial must stay inside the checked physical step.
    pub(super) fn require_untimed_materials(&self) -> Result<(), SolidError> {
        let has_creep = self.elements.iter().any(|element| {
            element.sections.iter().any(|section| match section {
                SectionModel::Fiber(section) => section
                    .fibers
                    .iter()
                    .any(|fiber| matches!(
                        &fiber.law,
                        FiberLaw::Norton(_)
                            | FiberLaw::NortonJ2(_)
                            | FiberLaw::ThermalNorton(_)
                            | FiberLaw::ThermalNortonJ2(_)
                    )),
                SectionModel::Elastic { .. } => false,
            })
        });
        if has_creep {
            Err(SolidError::MaterialRefused {
                what: "timed creep fibers require solve_creep_step with an explicit physical interval; untimed static/modal and inertial frame paths do not own creep history".into(),
            })
        } else {
            Ok(())
        }
    }

    /// Read accepted material history and retained source receipts without
    /// evaluating a duration-bearing constitutive trial.
    pub fn creep_fiber(
        &self,
        element: usize,
        section: usize,
        fiber: usize,
    ) -> Option<(&UniaxialNorton, &UniaxialNortonState)> {
        let SectionModel::Fiber(s) = self.elements.get(element)?.sections.get(section)? else {
            return None;
        };
        match (&s.fibers.get(fiber)?.law, s.states.get(fiber)?) {
            (FiberLaw::Norton(law), FiberState::Norton(state)) => Some((law, state)),
            _ => None,
        }
    }

    /// Read both accepted memories and source receipts without another trial.
    pub fn creep_plastic_fiber(
        &self,
        element: usize,
        section: usize,
        fiber: usize,
    ) -> Option<(&UniaxialNortonJ2, &UniaxialNortonJ2State)> {
        let SectionModel::Fiber(s) = self.elements.get(element)?.sections.get(section)? else {
            return None;
        };
        match (&s.fibers.get(fiber)?.law, s.states.get(fiber)?) {
            (FiberLaw::NortonJ2(law), FiberState::NortonJ2(state)) => Some((law, state)),
            _ => None,
        }
    }

    /// Read accepted sourced thermocreep history without another physical trial.
    pub fn thermal_creep_fiber(
        &self,
        element: usize,
        section: usize,
        fiber: usize,
    ) -> Option<(&UniaxialThermalNorton, &UniaxialThermalNortonState)> {
        let SectionModel::Fiber(s) = self.elements.get(element)?.sections.get(section)? else {
            return None;
        };
        match (&s.fibers.get(fiber)?.law, s.states.get(fiber)?) {
            (FiberLaw::ThermalNorton(law), FiberState::ThermalNorton(state)) => {
                Some((law, state))
            }
            _ => None,
        }
    }

    /// Read accepted thermal creep/plastic histories and temperature without aging them.
    pub fn thermal_creep_plastic_fiber(
        &self,
        element: usize,
        section: usize,
        fiber: usize,
    ) -> Option<(&UniaxialThermalNortonJ2, &UniaxialThermalNortonJ2State)> {
        let SectionModel::Fiber(s) = self.elements.get(element)?.sections.get(section)? else {
            return None;
        };
        match (&s.fibers.get(fiber)?.law, s.states.get(fiber)?) {
            (FiberLaw::ThermalNortonJ2(law), FiberState::ThermalNortonJ2(state)) => {
                Some((law, state))
            }
            _ => None,
        }
    }

    /// Rebind the next explicit interval/temperature without evolving history.
    /// Laws are ordered by member, original Lobatto section, then fiber. Elastic
    /// and hardening coefficients, kinetic parameters, card identities, selected claims and
    /// nonthermal source coordinates must remain identical. `temperature_axis`
    /// is the original explicit coordinate used by each sourced family; authored
    /// laws have no source coordinate. Reference geometry and inertia stay fixed.
    /// Norton and joint Norton/J2 laws retain their original history family.
    /// Sourced yield stress may follow temperature; authored yield stays fixed.
    pub fn with_creep_laws<L: Clone + Into<FiberLaw>>(
        &self,
        laws: &[L],
        temperature_axis: &str,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, FrameCreepError> {
        poll(&mut cancelled)?;
        let mut next = self.clone();
        let mut index = 0;
        for e in &mut next.elements {
            for (p, section) in e.sections.iter_mut().enumerate() {
                let SectionModel::Fiber(section) = section else {
                    return Err(FrameCreepError::Invalid(
                        "creep rebinding requires native timed fibers",
                    ));
                };
                require(
                    section.fibers.len() == section.states.len(),
                    "fiber/history count mismatch",
                )?;
                for (fiber, state) in section.fibers.iter_mut().zip(&section.states) {
                    poll(&mut cancelled)?;
                    let old = CreepPoint::new(&fiber.law, state)?;
                    let next_law: FiberLaw = laws
                        .get(index)
                        .ok_or(FrameCreepError::Invalid("missing next fiber law"))?
                        .clone()
                        .into();
                    let law = CreepPoint::new(&next_law, state)?;
                    require(
                        old.kinetics().elastic() == law.kinetics().elastic()
                            && old.kinetics().parameters() == law.kinetics().parameters()
                            && old.hardening() == law.hardening(),
                        "creep history cannot silently change elastic/hardening storage or kinetic parameters",
                    )?;
                    match (old.material(), law.material()) {
                        (None, None) => require(
                            old.yield_stress() == law.yield_stress(),
                            "changing authored yield requires an explicit source evolution owner",
                        )?,
                        (Some(a), Some(b)) => {
                            require(
                                a.card_identity() == b.card_identity()
                                    && a.axis_quantities() == b.axis_quantities()
                                    && a.properties().len() == b.properties().len()
                                    && a.properties().iter().zip(b.properties()).all(|(a, b)| {
                                        a.requirement() == b.requirement()
                                            && a.answer().receipt.selected
                                                == b.answer().receipt.selected
                                    })
                                    && a.query_point().len() == b.query_point().len()
                                    && a.query_point().iter().zip(b.query_point()).all(
                                        |((ka, va), (kb, vb))| {
                                            ka == kb
                                                && if ka == temperature_axis {
                                                    *va == old.kinetics().temperature_k()
                                                        && *vb == law.kinetics().temperature_k()
                                                } else {
                                                    va.to_bits() == vb.to_bits()
                                                }
                                        },
                                    )
                                    && a.state_coordinate(temperature_axis).is_some(),
                                "creep refresh must retain the original card, claims and fixed source context",
                            )?;
                        }
                        _ => {
                            return Err(FrameCreepError::Invalid(
                                "authored and sourced histories cannot be interchanged",
                            ));
                        }
                    }
                    match (old.expansion_material(), law.expansion_material()) {
                        (None, None) => {}
                        (Some(a), Some(b)) => {
                            require(
                                a.card_identity() == b.card_identity()
                                    && a.axis_quantities() == b.axis_quantities()
                                    && a.properties().len() == b.properties().len()
                                    && a.properties().iter().zip(b.properties()).all(|(a, b)| {
                                        a.requirement() == b.requirement()
                                            && a.answer().receipt.selected
                                                == b.answer().receipt.selected
                                    }),
                                "thermal creep refresh must retain its original expansion source",
                            )?;
                        }
                        _ => {
                            return Err(FrameCreepError::Invalid(
                                "thermal and zero-expansion creep histories cannot be interchanged",
                            ));
                        }
                    }
                    law.tensor_strain(e.commit.d[p].0 - fiber.y * e.commit.d[p].1)?;
                    fiber.law = next_law;
                    index += 1;
                }
            }
        }
        require(
            index == laws.len() && index > 0,
            "one next law per original timed fiber required",
        )?;
        next.reset_rayleigh_reference();
        poll(&mut cancelled)?;
        Ok(next)
    }

    /// Solve one absolute load endpoint over the common duration carried by all
    /// fiber laws. The existing force-based solver sees the same old history in
    /// section Newton iterations, compatibility substeps and backtracking.
    /// Independent SI force/moment and whole-frame work gates precede return.
    /// Failure/cancellation leaves `self` unchanged. Use the returned `frame`
    /// for continuation; `element_forces` reads its stored balanced endpoint.
    pub fn solve_creep_step(
        &self,
        load: &[f64],
        tolerances: FrameCreepTolerances,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FrameCreepStep, FrameCreepError> {
        poll(&mut cancelled)?;
        require(
            self.geometry == Geometry::Linear
                && !self.elements.is_empty()
                && self.vel.iter().chain(&self.acc).all(|v| *v == 0.0),
            "checked frame creep requires nonempty linear geometry and stationary inertia",
        )?;
        require(
            load.len() == self.ndof()
                && load.iter().all(|x| x.is_finite())
                && self.u.iter().all(|x| x.is_finite())
                && [
                    tolerances.energy_j,
                    tolerances.force_n,
                    tolerances.moment_n_m,
                    tolerances.material.stress_pa,
                    tolerances.material.energy_j_m3,
                ]
                .into_iter()
                .all(positive)
                && tolerances.maximum_fiber_points > 0,
            "invalid creep loads, controls or accepted displacement",
        )?;
        let mut count = 0usize;
        let mut duration = None;
        for e in &self.elements {
            require(
                e.viscous_section.is_none(),
                "frame viscosity needs a separate timed work owner",
            )?;
            for (p, section) in e.sections.iter().enumerate() {
                let SectionModel::Fiber(section) = section else {
                    return Err(FrameCreepError::Invalid(
                        "checked creep work requires timed fiber sections",
                    ));
                };
                require(
                    !section.fibers.is_empty() && section.fibers.len() == section.states.len(),
                    "nonempty matching fiber/history arrays required",
                )?;
                for (fiber, state) in section.fibers.iter().zip(&section.states) {
                    poll(&mut cancelled)?;
                    count = count
                        .checked_add(1)
                        .ok_or(FrameCreepError::Invalid("fiber count overflow"))?;
                    require(
                        count <= tolerances.maximum_fiber_points,
                        "fiber work budget exceeded",
                    )?;
                    let law = CreepPoint::new(&fiber.law, state)?;
                    require(
                        positive(fiber.area) && fiber.y.is_finite(),
                        "invalid fiber area or position",
                    )?;
                    let dt = law.kinetics().duration_s();
                    if let Some(previous) = duration {
                        require(
                            dt == previous,
                            "all fibers need one common physical duration",
                        )?;
                    } else {
                        duration = Some(dt);
                    }
                    let strain = e.commit.d[p].0 - fiber.y * e.commit.d[p].1;
                    law.tensor_strain(strain)?;
                    self.admit_creep_geometry(e, law.kinetics().elastic().strain_limit)?;
                }
            }
        }
        let mut next = self.clone();
        let cancelled_seen = Cell::new(false);
        let mut checkpoint = || {
            if cancelled_seen.get() || cancelled() {
                cancelled_seen.set(true);
                Err(SolidError::InvalidInput {
                    what: "frame creep cancelled".into(),
                })
            } else {
                Ok(())
            }
        };
        let map_native = |e| {
            if cancelled_seen.get() {
                FrameCreepError::Cancelled
            } else {
                FrameCreepError::Mechanical(e)
            }
        };
        let mut steps = next
            .static_load_controlled(load, 1, &mut checkpoint)
            .map_err(map_native)?;
        // Reconstruct endpoint forces with OLD history, never age the newly
        // committed state just to observe its endpoint stress or support work.
        let (forces, _, _) = self
            .assemble_controlled(&next.u, &mut checkpoint)
            .map_err(map_native)?;
        drop(checkpoint);
        let mut energy = FrameCreepEnergy::default();
        for (dof, (&force, &applied)) in forces.iter().zip(load).enumerate() {
            poll(&mut cancelled)?;
            let residual = force - applied;
            let tolerance = if dof % 3 == 2 {
                tolerances.moment_n_m
            } else {
                tolerances.force_n
            };
            if !residual.is_finite() || (!self.fixed[dof] && residual.abs() > tolerance) {
                return Err(FrameCreepError::Equilibrium {
                    dof,
                    residual,
                    tolerance,
                });
            }
            let delta = next.u[dof] - self.u[dof];
            energy.applied_work_j += applied * delta;
            if self.fixed[dof] {
                energy.support_work_j += residual * delta;
            }
        }
        let mut fibers = Vec::with_capacity(count);
        let mut material_work = 0.0;
        for (element, (old, new)) in self.elements.iter().zip(&next.elements).enumerate() {
            for (p, section) in old.sections.iter().enumerate() {
                let SectionModel::Fiber(section) = section else {
                    unreachable!()
                };
                for (f, (fiber, state)) in section.fibers.iter().zip(&section.states).enumerate() {
                    poll(&mut cancelled)?;
                    let law = CreepPoint::new(&fiber.law, state)?;
                    next.admit_creep_geometry(new, law.kinetics().elastic().strain_limit)?;
                    let old_strain = old.commit.d[p].0 - fiber.y * old.commit.d[p].1;
                    let new_strain = new.commit.d[p].0 - fiber.y * new.commit.d[p].1;
                    let SectionModel::Fiber(committed) = &new.sections[p] else {
                        unreachable!()
                    };
                    let committed =
                        CreepPoint::new(&committed.fibers[f].law, &committed.states[f])?;
                    let account =
                        law.audit(committed, old_strain, new_strain, tolerances.material)?;
                    let volume = fiber.area * old.w[p] * old.l0;
                    require(
                        positive(volume),
                        "fiber reference volume is not representable",
                    )?;
                    energy.initial_stored_j += volume * account.initial_stored_j_m3;
                    energy.final_stored_j += volume * account.final_stored_j_m3;
                    energy.creep_dissipation_j += volume * account.creep_dissipation_j_m3;
                    energy.plastic_dissipation_j += volume * account.plastic_dissipation_j_m3;
                    energy.physical_dissipation_j += volume * account.physical_dissipation_j_m3;
                    energy.thermoelastic_storage_change_j +=
                        volume * account.thermoelastic_storage_change_j_m3;
                    energy.reversible_heat_j += volume * account.reversible_heat_j_m3;
                    energy.numerical_dissipation_j += volume * account.numerical_dissipation_j_m3;
                    material_work += volume * account.mechanical_work_j_m3;
                    fibers.push(FrameFiberCreepEnergy {
                        element,
                        section: p,
                        fiber: f,
                        reference_volume_m3: volume,
                        energy: account,
                    });
                }
            }
        }
        let work = energy.applied_work_j + energy.support_work_j;
        energy.equilibrium_work_residual_j = work - material_work;
        energy.balance_residual_j = work
            - (energy.final_stored_j - energy.initial_stored_j)
            - energy.thermoelastic_storage_change_j
            - energy.physical_dissipation_j
            - energy.reversible_heat_j
            - energy.numerical_dissipation_j;
        energy_gate(energy.equilibrium_work_residual_j, tolerances.energy_j)?;
        energy_gate(energy.balance_residual_j, tolerances.energy_j)?;
        poll(&mut cancelled)?;
        Ok(FrameCreepStep {
            frame: next,
            duration_s: duration.expect("admitted nonempty fibers"),
            equilibrium: steps.remove(0),
            fibers,
            energy,
        })
    }

    fn admit_creep_geometry(
        &self,
        e: &super::FrameElement,
        limit: f64,
    ) -> Result<(), FrameCreepError> {
        let dx = (self.u[3 * e.j] - self.u[3 * e.i]) / e.l0;
        let dy = (self.u[3 * e.j + 1] - self.u[3 * e.i + 1]) / e.l0;
        // Endpoint angles alone can hide large opposing internal curvatures.
        let rotation_variation = e
            .commit
            .d
            .iter()
            .zip(&e.w)
            .map(|((_, curvature), weight)| e.l0 * weight * curvature.abs())
            .sum::<f64>();
        require(
            dx.is_finite()
                && dy.is_finite()
                && dx.hypot(dy) <= limit
                && rotation_variation.is_finite()
                && self.u[3 * e.i + 2].abs() + rotation_variation <= limit
                && self.u[3 * e.j + 2].abs() <= limit,
            "member displacements or rotations exceed the admitted small-geometry envelope",
        )
    }
}
