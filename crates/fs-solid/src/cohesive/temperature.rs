//! Temperature evolution of the original deformable cohesive interfaces.
//!
//! This is the fixed-jump coefficient-change and irreversible-relaxation model
//! owned by `fs_material::cohesive::temperature`, applied at every original
//! quadrature site. It preserves geometry, elastic operator, displacement and
//! accumulated fracture history. A subsequent native equilibrium or inertial
//! step responds to the changed bond tractions; transfer itself performs no
//! displacement work and does not assert equilibrium at the new temperature.
//!
//! Parameter work and fixed-jump damage release are separate surface-energy
//! accounts. Neither is automatically heat. Bulk thermal expansion, changes of
//! bulk elasticity, moving frames, remeshing and thermal transport require
//! their own operators and work accounts; they are not inferred here.

use fs_exec::Cx;
use fs_material::cohesive::source::ResolvedBilinearCohesive;

use super::{CohesivePatch, CohesiveState, CohesiveSystem, Error, finite, positive};

/// Explicit energy acceptance gates for this zero-duration source transition.
#[derive(Debug, Clone, Copy)]
pub struct TemperatureConfig {
    /// Pointwise constitutive work/storage/loss gate per original bonded area.
    pub local_energy_tolerance_j_m2: f64,
    /// Gate on the area-integrated interface energy residual.
    pub total_energy_tolerance_j: f64,
}

/// Area-integrated constitutive accounts, in joules, on original bonded areas.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TemperatureEnergy {
    /// Signed work of changing coefficients at retained damage and fixed jump.
    pub parameter_work_j: f64,
    /// Irreversible loss from any additional damage at that same fixed jump.
    pub damage_dissipated_j: f64,
    /// Complete recoverable interface-energy change across both stages.
    pub stored_change_j: f64,
    /// Parameter work minus storage change minus damage dissipation.
    pub defect_j: f64,
}

/// Complete immutable candidate. The old accepted state remains retryable.
pub struct TemperatureStep<'m, 'a> {
    /// Original displacement and irreversible history at the new source state.
    pub state: CohesiveState<'m, 'a>,
    /// Integrated fixed-jump coefficient-work and damage-release accounts.
    pub energy: TemperatureEnergy,
}

impl<'a> CohesiveSystem<'a> {
    /// Prepare new temperature sources on exactly the same reference patches,
    /// trace maps, areas, DOF ordering, elastic operator and admission budgets.
    /// One successor is required per original patch, in its existing order.
    /// Each successor must retain its original card, selected claims, ordered
    /// surfaces and all nonthermal conditions. This prepares coefficients only;
    /// accepted history must use `CohesiveState::transfer_temperature` below.
    pub fn with_temperature_sources<'b>(
        &self,
        successors: &[&'b ResolvedBilinearCohesive],
        temperature_axis: &str,
        cx: &Cx<'_>,
    ) -> Result<CohesiveSystem<'b>, Error>
    where
        'a: 'b,
    {
        self.with_temperature_sources_controlled(successors, temperature_axis, &mut || {
            cx.checkpoint().is_err()
        })
    }

    /// Predicate variant with the same source checks and bounded per-site polls.
    pub fn with_temperature_sources_controlled<'b, C: FnMut() -> bool>(
        &self,
        successors: &[&'b ResolvedBilinearCohesive],
        temperature_axis: &str,
        cancelled: &mut C,
    ) -> Result<CohesiveSystem<'b>, Error>
    where
        'a: 'b,
    {
        poll(cancelled)?;
        if successors.len() != self.patches.len() {
            return Err(Error::Invalid(
                "one temperature successor per original patch required",
            ));
        }
        let mut patches = Vec::with_capacity(self.patches.len());
        for (old, &next) in self.patches.iter().zip(successors) {
            poll(cancelled)?;
            require_source(old, next, temperature_axis)?;
            patches.push(CohesivePatch {
                id: old.id,
                area_m2: old.area_m2,
                law: next.law(),
                source: Some(next),
                terms: old.terms.clone(),
            });
        }
        poll(cancelled)?;
        Ok(CohesiveSystem {
            elastic: self.elastic,
            patches,
            n: self.n,
            budget: self.budget,
        })
    }
}

impl CohesiveState<'_, '_> {
    /// Transfer the original accepted history to source-defined temperatures.
    /// The target may be prepared with `with_temperature_sources`; source and
    /// mechanics identity are checked again against THIS accepted state.
    /// All points and integrated accounts publish together, after cancellation
    /// and local/global energy checks. No time advances and no bond is reset.
    pub fn transfer_temperature<'n, 'b>(
        &self,
        target: &'n CohesiveSystem<'b>,
        temperature_axis: &str,
        config: TemperatureConfig,
        cx: &Cx<'_>,
    ) -> Result<TemperatureStep<'n, 'b>, Error> {
        self.transfer_temperature_controlled(target, temperature_axis, config, &mut || {
            cx.checkpoint().is_err()
        })
    }

    /// Predicate variant. Original numerical work bounds already limit the
    /// number of sites, DOFs and trace columns. Polls every site and at most
    /// 256 copied displacement entries, including before final publication.
    pub fn transfer_temperature_controlled<'n, 'b, C: FnMut() -> bool>(
        &self,
        target: &'n CohesiveSystem<'b>,
        temperature_axis: &str,
        config: TemperatureConfig,
        cancelled: &mut C,
    ) -> Result<TemperatureStep<'n, 'b>, Error> {
        poll(cancelled)?;
        positive(config.local_energy_tolerance_j_m2)?;
        positive(config.total_energy_tolerance_j)?;
        if !core::ptr::eq(self.system.elastic, target.elastic)
            || self.system.n != target.n
            || self.system.patches.len() != target.patches.len()
            || target.elastic.n() != target.n
        {
            return Err(Error::Invalid(
                "temperature transfer requires the original elastic operator and DOFs",
            ));
        }
        let mut energy = TemperatureEnergy {
            parameter_work_j: 0.0,
            damage_dissipated_j: 0.0,
            stored_change_j: 0.0,
            defect_j: 0.0,
        };
        let mut points = Vec::with_capacity(self.points.len());
        for ((old, next), point) in self
            .system
            .patches
            .iter()
            .zip(&target.patches)
            .zip(&self.points)
        {
            poll(cancelled)?;
            if old.id != next.id || old.area_m2 != next.area_m2 || old.terms != next.terms {
                return Err(Error::Invalid(
                    "temperature transfer cannot change patch identity, area or trace map",
                ));
            }
            let source = next.source.ok_or(Error::Invalid(
                "temperature transfer requires sourced interface patches",
            ))?;
            require_source(old, source, temperature_axis)?;
            let step = point
                .transfer_temperature(next.law, config.local_energy_tolerance_j_m2)
                .map_err(|error| Error::Material {
                    patch: old.id,
                    error,
                })?;
            energy.parameter_work_j =
                finite(energy.parameter_work_j + old.area_m2 * step.parameter_work_j_m2)?;
            energy.damage_dissipated_j =
                finite(energy.damage_dissipated_j + old.area_m2 * step.damage_dissipated_j_m2)?;
            energy.stored_change_j =
                finite(energy.stored_change_j + old.area_m2 * step.stored_change_j_m2)?;
            points.push(step.point);
        }
        energy.defect_j =
            finite(energy.parameter_work_j - energy.stored_change_j - energy.damage_dissipated_j)?;
        if energy.defect_j.abs() > config.total_energy_tolerance_j {
            return Err(Error::CohesiveEnergyDefect(energy.defect_j));
        }
        let mut displacements = Vec::with_capacity(self.displacements.len());
        for chunk in self.displacements.chunks(256) {
            poll(cancelled)?;
            displacements.extend_from_slice(chunk);
        }
        poll(cancelled)?;
        Ok(TemperatureStep {
            state: CohesiveState {
                system: target,
                displacements,
                points,
            },
            energy,
        })
    }
}

fn require_source(
    original: &CohesivePatch<'_>,
    successor: &ResolvedBilinearCohesive,
    temperature_axis: &str,
) -> Result<(), Error> {
    original
        .source
        .ok_or(Error::Invalid(
            "temperature transfer requires sourced interface patches",
        ))?
        .require_temperature_successor(successor, temperature_axis)
        .map_err(|error| Error::Source {
            patch: original.id,
            error,
        })
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
