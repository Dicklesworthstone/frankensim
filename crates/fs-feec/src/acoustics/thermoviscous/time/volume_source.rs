//! Causal body-force and heat excitation of the original complete fluid state.
//!
//! The harmonic source owners provide the original weak spatial loads. Held
//! real force amplitudes [N] enter the vector momentum equations, including
//! the constrained moving-wall rows, before the SAME full mixed mass inverse.
//! Heating [W] enters the translational thermal rows. Every pressure, velocity,
//! temperature, molecular, finite-wall and structural coordinate is retained.
//!
//! Source work is integrated along that same prepared trajectory:
//! integral(f . v dV dt) and integral(q' T'/T0 dV dt), with no harmonic
//! factor of one-half. Signed work can remove energy. Geometry, equilibrium
//! and the force/heat mechanisms remain caller-owned; no mean flow, ALE,
//! nonlinear momentum law, source feedback or convergence claim is inferred.

use super::super::body_force::ThermoviscousBodyForceSource;
use super::super::heat_source::ThermoviscousHeatSource;
use super::*;

mod interface;
pub use interface::{
    PreparedConformingThermoviscousVolumeTimePair, PreparedThermoviscousVolumeTimePair,
    ThermoviscousTimeRegion,
};

/// Immutable fixed-clock evolution with original distributed volume inputs.
///
/// The underlying field cache stays private: its extra force and heat ports
/// cannot enter the legacy pressure-only interface APIs. Typed pair bindings
/// carry those source ports through the original joint trajectory. Preparation changes
/// input columns and their integral observations, not the physical state.
pub struct PreparedThermoviscousVolumeTime<'a, 'c> {
    fields: PreparedThermoviscousTime<'a, 'c>,
    force_count: usize,
    heat_count: usize,
}

/// One complete candidate and its original spatial source observations/work.
#[derive(Clone, Debug)]
pub struct ThermoviscousVolumeTimeStep {
    /// Original endpoint, mechanical flow means and complete energy account.
    /// Its work includes both sums of the source work reported below.
    pub field: ThermoviscousTimeStep,
    /// Time mean of integral(shape_s . v dV) [m/s], in force-source order.
    /// Shapes retain their actual sign and magnitude, including zero resultant.
    /// The full physical velocity includes the original moving-wall lifting.
    pub mean_body_force_velocity_m_s: Vec<f64>,
    /// Signed integral(force_amplitude_s * shape_s . v dV dt) [J] per source.
    pub body_force_work_by_source_j: Vec<f64>,
    /// Time mean of integral(shape_s * T'/T0 dV), in heat-source order.
    /// This is dimensionless and is not divided by the shape's integral.
    pub mean_heat_source_temperature_ratio: Vec<f64>,
    /// Signed integral(heat_amplitude_s * shape_s * T'/T0 dV dt) [J].
    pub heat_source_work_by_source_j: Vec<f64>,
}

impl<'c> ThermoviscousFe<'c> {
    /// Prepare held original body forces and heating with every resolved field.
    ///
    /// At least one volume source is required. Sources must belong to this
    /// exact fluid/structure owner; their order, sign and normalization survive.
    /// Original callbacks are never invoked during preparation or warm steps.
    /// Force-source count and the complete prepared port/state/work budgets
    /// are admitted before dense preparation; no source or field is truncated.
    ///
    /// Frozen Cp plus every molecular capacity must equal the medium's Cp.
    /// Finite conducting solids require their original depth discretization;
    /// optional tangential coefficients select the existing spatial wall owner.
    /// Pass None for both wall arguments when there are no conducting regions.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_time_with_volume_sources<'a>(
        &'a self,
        duration_s: f64,
        frozen_cp_j_kg_k: f64,
        reservoirs: &[ThermoviscousCaloricReservoir],
        walls: Option<ThermoviscousWallDiscretization<'_>>,
        tangential: Option<&[TangentialWallConductivity]>,
        force_sources: &[&ThermoviscousBodyForceSource<'_, '_>],
        heat_sources: &[&ThermoviscousHeatSource<'_, '_>],
        limits: ThermoviscousTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedThermoviscousVolumeTime<'a, 'c>, ThermoviscousTimeError> {
        poll(cancelled)?;
        if (force_sources.is_empty() && heat_sources.is_empty())
            || (tangential.is_some() && walls.is_none())
        {
            return Err(invalid(
                "original volume sources and explicit spatial-wall resolution required",
            ));
        }
        let fields = PreparedThermoviscousTime::new_with_volume_sources(
            self,
            duration_s,
            frozen_cp_j_kg_k,
            reservoirs,
            walls,
            tangential,
            force_sources,
            heat_sources,
            limits,
            cancelled,
        )?;
        poll(cancelled)?;
        Ok(PreparedThermoviscousVolumeTime {
            fields,
            force_count: force_sources.len(),
            heat_count: heat_sources.len(),
        })
    }
}

impl<'c> PreparedThermoviscousVolumeTime<'_, 'c> {
    /// Original resolved spatial owner.
    #[must_use]
    pub fn inner(&self) -> &ThermoviscousFe<'c> {
        self.fields.inner()
    }

    /// Complete original pHS system. Port order is pressure, structural force,
    /// original body-force profiles, then original heat profiles.
    #[must_use]
    pub fn system(&self) -> &PortHamiltonian {
        self.fields.system()
    }

    #[must_use]
    pub const fn duration_s(&self) -> f64 {
        self.fields.duration_s
    }

    #[must_use]
    pub const fn body_force_source_count(&self) -> usize {
        self.force_count
    }

    #[must_use]
    pub const fn heat_source_count(&self) -> usize {
        self.heat_count
    }

    /// Complete all-rest state, including molecular, solid and structural fields.
    #[must_use]
    pub fn zero_state(&self) -> ThermoviscousTimeState {
        self.fields.zero_state()
    }

    /// Stored acoustic energy from the original complete mixed mass [J].
    pub fn stored_energy_j(
        &self,
        state: &ThermoviscousTimeState,
    ) -> Result<f64, ThermoviscousTimeError> {
        self.fields.stored_energy_j(state)
    }

    /// Original frozen translational heat capacity [J/(kg K)].
    #[must_use]
    pub const fn frozen_cp_j_kg_k(&self) -> f64 {
        self.fields.frozen_cp_j_kg_k
    }

    #[must_use]
    pub fn caloric_reservoirs(&self) -> &[ThermoviscousCaloricReservoir] {
        self.fields.caloric_reservoirs()
    }

    /// Original full surface/depth ordering of each finite conducting solid.
    #[must_use]
    pub fn conductive_wall_layouts(&self) -> &[ConductiveWallTimeLayout] {
        self.fields.conductive_wall_layouts()
    }

    /// Upper frequency retained by the original finite-wall depth model.
    #[must_use]
    pub fn maximum_wall_frequency_hz(&self) -> Option<f64> {
        self.fields.maximum_wall_frequency_hz()
    }

    /// Original in-plane solid conductivity, when spatial walls were requested.
    #[must_use]
    pub fn tangential_wall_conductivity(&self) -> Option<&[TangentialWallConductivity]> {
        self.fields.tangential_wall_conductivity()
    }

    /// Exact interval-mean pressure-port relation with all original volume drives.
    /// The selected traction slot is ignored; every force, heat input and
    /// complete initial field contributes to the free inward flow. Preparation
    /// already computed these maps, so this performs no trial propagation.
    #[allow(clippy::too_many_arguments)]
    pub fn average_relation(
        &self,
        initial: &ThermoviscousTimeState,
        traction_pa: &[f64],
        wall_forces_n_per_sqrt_kg: &[f64],
        force_drive_n: &[f64],
        heat_drive_w: &[f64],
        selected_port: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), ThermoviscousTimeError> {
        let (free, mobility) = self.average_relations(
            initial,
            traction_pa,
            wall_forces_n_per_sqrt_kg,
            force_drive_n,
            heat_drive_w,
            &[selected_port],
            cancelled,
        )?;
        Ok((free[0], mobility[0]))
    }

    /// Complete `Q = free + mobility * pressure` relation in selected-port order.
    /// Every cross-port term and every original force/heat column is retained.
    /// Selected traction slots are ignored (and may be NaN); all other drives
    /// must be finite and complete. No source is resampled or state changed.
    #[allow(clippy::too_many_arguments)]
    pub fn average_relations(
        &self,
        initial: &ThermoviscousTimeState,
        traction_pa: &[f64],
        wall_forces_n_per_sqrt_kg: &[f64],
        force_drive_n: &[f64],
        heat_drive_w: &[f64],
        selected_ports: &[usize],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, Vec<f64>), ThermoviscousTimeError> {
        poll(cancelled)?;
        if traction_pa.len() != self.fields.model.loads.len()
            || wall_forces_n_per_sqrt_kg.len() != self.fields.wall_scale.len()
            || force_drive_n.len() != self.force_count
            || heat_drive_w.len() != self.heat_count
        {
            return Err(invalid(
                "complete original pressure, structural, body-force and heat drives required",
            ));
        }
        self.fields.average_relations_from_inputs(
            initial,
            &[
                traction_pa,
                wall_forces_n_per_sqrt_kg,
                force_drive_n,
                heat_drive_w,
            ],
            selected_ports,
            cancelled,
        )
    }

    /// One complete candidate under real held mechanical and thermal inputs.
    ///
    /// Body forces [N] and heating [W] multiply their original spatial profiles.
    /// Source work uses exact interval integrals from the cached complete flow,
    /// not endpoint samples. State and clock remain caller-owned and unchanged
    /// on every refusal, including cancellation after propagation.
    pub fn step(
        &self,
        initial: &ThermoviscousTimeState,
        traction_pa: &[f64],
        wall_forces_n_per_sqrt_kg: &[f64],
        force_drive_n: &[f64],
        heat_drive_w: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermoviscousVolumeTimeStep, ThermoviscousTimeError> {
        poll(cancelled)?;
        let prepared = &self.fields;
        let pressure_ports = prepared.model.loads.len();
        let force_start = pressure_ports + prepared.wall_scale.len();
        let heat_start = force_start + self.force_count;
        if traction_pa.len() != pressure_ports
            || wall_forces_n_per_sqrt_kg.len() != prepared.wall_scale.len()
            || force_drive_n.len() != self.force_count
            || heat_drive_w.len() != self.heat_count
            || traction_pa
                .iter()
                .chain(wall_forces_n_per_sqrt_kg)
                .chain(force_drive_n)
                .chain(heat_drive_w)
                .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "complete finite held pressure, structural, body-force and heat drives required",
            ));
        }
        let input: Vec<_> = traction_pa
            .iter()
            .chain(wall_forces_n_per_sqrt_kg)
            .chain(force_drive_n)
            .chain(heat_drive_w)
            .copied()
            .collect();
        let x = prepared.encode(initial)?;
        let candidate = prepared.flow.evaluate(&x, &input)?;
        poll(cancelled)?;
        let body_force_work_by_source_j = source_work(
            force_drive_n,
            &candidate.integrals[HEAT_CHANNELS + force_start..HEAT_CHANNELS + heat_start],
        )?;
        let heat_source_work_by_source_j = source_work(
            heat_drive_w,
            &candidate.integrals[HEAT_CHANNELS + heat_start..],
        )?;
        let work = ThermoviscousTimeWork {
            initial_energy_j: prepared.system.hamiltonian(&x),
            final_energy_j: prepared.system.hamiltonian(&candidate.state),
            delta_energy_j: energy_change_j(&prepared.system, &x, &candidate.state),
            pressure_work_j: traction_pa
                .iter()
                .enumerate()
                .map(|(port, &drive)| drive * candidate.integrals[HEAT_CHANNELS + port])
                .sum(),
            structural_work_j: wall_forces_n_per_sqrt_kg
                .iter()
                .enumerate()
                .map(|(mode, &drive)| {
                    drive * candidate.integrals[HEAT_CHANNELS + pressure_ports + mode]
                })
                .sum(),
            body_force_work_j: sum_work(&body_force_work_by_source_j),
            heat_source_work_j: sum_work(&heat_source_work_by_source_j),
            conductive_boundary_work_j: 0.0,
            fluid_viscous_heat_j: candidate.integrals[0],
            fluid_conductive_heat_j: candidate.integrals[1],
            structural_damping_heat_j: candidate.integrals[2],
            molecular_relaxation_heat_j: candidate.integrals[3],
            solid_conductive_heat_j: candidate.integrals[4],
        };
        admit_work(
            work,
            &candidate.integrals,
            prepared.maximum_relative_energy_residual,
        )?;
        let means: Vec<_> = candidate.integrals[HEAT_CHANNELS..]
            .iter()
            .map(|integral| integral / prepared.duration_s)
            .collect();
        if means.iter().any(|mean| !mean.is_finite()) {
            return Err(invalid("nonfinite original transient source/flow means"));
        }
        let state = prepared.decode(&candidate.state);
        prepared.encode(&state)?;
        poll(cancelled)?;
        Ok(ThermoviscousVolumeTimeStep {
            field: ThermoviscousTimeStep {
                state,
                mean_inward_flow_m3_s: means[..pressure_ports].to_vec(),
                mean_wall_velocity_m_sqrt_kg_per_s: means[pressure_ports..force_start].to_vec(),
                work,
            },
            mean_body_force_velocity_m_s: means[force_start..heat_start].to_vec(),
            body_force_work_by_source_j,
            mean_heat_source_temperature_ratio: means[heat_start..].to_vec(),
            heat_source_work_by_source_j,
        })
    }
}

fn source_work(drives: &[f64], integrals: &[f64]) -> Result<Vec<f64>, ThermoviscousTimeError> {
    let work: Vec<_> = drives
        .iter()
        .zip(integrals)
        .map(|(&drive, &integral)| (Dd::from_f64(drive) * Dd::from_f64(integral)).to_f64())
        .collect();
    if work.iter().any(|value| !value.is_finite()) {
        return Err(invalid("nonfinite original volumetric source work"));
    }
    Ok(work)
}

fn sum_work(work: &[f64]) -> f64 {
    work.iter()
        .fold(Dd::ZERO, |sum, &value| sum + Dd::from_f64(value))
        .to_f64()
}
