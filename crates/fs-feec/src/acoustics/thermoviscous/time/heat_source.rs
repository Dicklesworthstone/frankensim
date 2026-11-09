//! Causal volume/surface heating on the original complete thermoviscous state.
//!
//! Spatial heat loads are the SAME original P2 volume/surface integrals as the harmonic
//! heat-source owner. Real held amplitudes [W] enter the translational thermal
//! equation before the full mixed mass inverse is applied. No fitted source
//! response, reduced state, endpoint loss or prescribed pressure increment is
//! substituted for the fluid, molecular, solid and structural trajectory.

use super::super::heat_source::ThermoviscousHeatSource;
use super::*;

/// Immutable fixed-clock heat-enabled evolution of all original fields.
///
/// A separate type keeps the additional thermal ports out of the legacy
/// pressure-only interface APIs. Preparation retains the full original state;
/// sources add input columns and their exact integral observations only.
pub struct PreparedThermoviscousHeatTime<'a, 'c> {
    fields: PreparedThermoviscousTime<'a, 'c>,
    source_count: usize,
}

/// One admitted complete field trajectory and its spatial thermal-source work.
#[derive(Clone, Debug)]
pub struct ThermoviscousHeatTimeStep {
    /// Original endpoint, mechanical flow means and complete energy account.
    /// Its work includes the sum of the volumetric heat-source work below.
    pub field: ThermoviscousTimeStep,
    /// Time mean of integral(shape_s * T'/T0 dV), in source order.
    /// This dimensionless work-conjugate observation retains the actual signed
    /// shape; it is not divided by its possibly zero spatial integral.
    pub mean_heat_source_temperature_ratio: Vec<f64>,
    /// Signed integral(drive_s * shape_s * T'/T0 dV dt) [J] per source.
    /// A source can remove acoustic energy. No sign clipping is performed.
    pub heat_source_work_by_source_j: Vec<f64>,
}

impl<'c> ThermoviscousFe<'c> {
    /// Prepare real held volume or surface heating alongside every original field.
    ///
    /// Sources must be nonempty and belong to this exact FE owner. Their order,
    /// sign and normalization are retained. The original weak loads are copied
    /// into the physical input map once; the source callbacks are never invoked
    /// during preparation or warm stepping.
    ///
    /// Frozen Cp plus all molecular capacities must equal the medium's Cp.
    /// Finite solid fields require their original depth discretization; optional
    /// tangential coefficients select the existing spatial wall owner. Pass
    /// None for both wall arguments when no finite conductive regions exist.
    /// The preparation budget includes every appended heat port; exceeding it
    /// is a refusal and does not truncate any field or source.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_time_with_heat_sources<'a>(
        &'a self,
        duration_s: f64,
        frozen_cp_j_kg_k: f64,
        reservoirs: &[ThermoviscousCaloricReservoir],
        walls: Option<ThermoviscousWallDiscretization<'_>>,
        tangential: Option<&[TangentialWallConductivity]>,
        heat_sources: &[&ThermoviscousHeatSource<'_, '_>],
        limits: ThermoviscousTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedThermoviscousHeatTime<'a, 'c>, ThermoviscousTimeError> {
        poll(cancelled)?;
        if heat_sources.is_empty() || (tangential.is_some() && walls.is_none()) {
            return Err(invalid(
                "nonempty original heat sources and explicit spatial-wall resolution required",
            ));
        }
        let fields = PreparedThermoviscousTime::new_with_volume_sources(
            self,
            duration_s,
            frozen_cp_j_kg_k,
            reservoirs,
            walls,
            tangential,
            &[],
            heat_sources,
            limits,
            cancelled,
        )?;
        poll(cancelled)?;
        Ok(PreparedThermoviscousHeatTime {
            fields,
            source_count: heat_sources.len(),
        })
    }
}

impl<'c> PreparedThermoviscousHeatTime<'_, 'c> {
    /// Original resolved spatial owner.
    #[must_use]
    pub fn inner(&self) -> &ThermoviscousFe<'c> {
        self.fields.inner()
    }

    /// Complete immutable pHS system. Pressure ports precede structural-force
    /// ports, which precede the original volume heat sources in supplied order.
    #[must_use]
    pub fn system(&self) -> &PortHamiltonian {
        self.fields.system()
    }

    /// Fixed interval of the cached complete trajectory [s].
    #[must_use]
    pub const fn duration_s(&self) -> f64 {
        self.fields.duration_s
    }

    /// Number of independently held spatial heat drives.
    #[must_use]
    pub const fn source_count(&self) -> usize {
        self.source_count
    }

    /// Complete all-rest field, including molecular, solid and structural state.
    #[must_use]
    pub fn zero_state(&self) -> ThermoviscousTimeState {
        self.fields.zero_state()
    }

    /// Stored acoustic energy from the original full mixed mass [J].
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

    /// Original molecular capacities and times, in supplied order.
    #[must_use]
    pub fn caloric_reservoirs(&self) -> &[ThermoviscousCaloricReservoir] {
        self.fields.caloric_reservoirs()
    }

    /// Original full surface/depth ordering of each finite conducting solid.
    #[must_use]
    pub fn conductive_wall_layouts(&self) -> &[ConductiveWallTimeLayout] {
        self.fields.conductive_wall_layouts()
    }

    /// One complete candidate under real held pressure, force and heat drives.
    ///
    /// Heat amplitudes [W] multiply original volume [m^-3] or surface [m^-2] shapes.
    /// Source work uses the exact interval integral of the resolved temperature
    /// observation, with NO harmonic factor of one-half. It is acoustic
    /// available-energy input, not the first-order enthalpy supplied to the
    /// fixed equilibrium. State and accepted clock remain caller-owned and
    /// unchanged on every refusal, including cancellation after propagation.
    pub fn step(
        &self,
        initial: &ThermoviscousTimeState,
        traction_pa: &[f64],
        wall_forces_n_per_sqrt_kg: &[f64],
        heat_drive_w: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermoviscousHeatTimeStep, ThermoviscousTimeError> {
        poll(cancelled)?;
        let prepared = &self.fields;
        let pressure_ports = prepared.model.loads.len();
        let mechanical_ports = pressure_ports + prepared.wall_scale.len();
        if traction_pa.len() != pressure_ports
            || wall_forces_n_per_sqrt_kg.len() != prepared.wall_scale.len()
            || heat_drive_w.len() != self.source_count
            || traction_pa
                .iter()
                .chain(wall_forces_n_per_sqrt_kg)
                .chain(heat_drive_w)
                .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "complete finite held pressure, structural and heat drives required",
            ));
        }
        let input: Vec<_> = traction_pa
            .iter()
            .chain(wall_forces_n_per_sqrt_kg)
            .chain(heat_drive_w)
            .copied()
            .collect();
        let x = prepared.encode(initial)?;
        let candidate = prepared.flow.evaluate(&x, &input)?;
        poll(cancelled)?;
        let heat_source_work_by_source_j: Vec<_> = heat_drive_w
            .iter()
            .enumerate()
            .map(|(source, &drive)| {
                (Dd::from_f64(drive)
                    * Dd::from_f64(candidate.integrals[HEAT_CHANNELS + mechanical_ports + source]))
                .to_f64()
            })
            .collect();
        if heat_source_work_by_source_j
            .iter()
            .any(|work| !work.is_finite())
        {
            return Err(invalid("nonfinite original heat-source work"));
        }
        let heat_source_work_j = heat_source_work_by_source_j
            .iter()
            .fold(Dd::ZERO, |sum, &work| sum + Dd::from_f64(work))
            .to_f64();
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
            body_force_work_j: 0.0,
            heat_source_work_j,
            conductive_boundary_work_j: 0.0,
            fluid_viscous_heat_j: candidate.integrals[0],
            fluid_conductive_heat_j: candidate.integrals[1],
            structural_damping_heat_j: candidate.integrals[2],
            molecular_relaxation_heat_j: candidate.integrals[3],
            solid_conductive_heat_j: candidate.integrals[4],
            structural_relaxation_heat_j: candidate.integrals[5],
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
        Ok(ThermoviscousHeatTimeStep {
            field: ThermoviscousTimeStep {
                state,
                mean_inward_flow_m3_s: means[..pressure_ports].to_vec(),
                mean_wall_velocity_m_sqrt_kg_per_s: means[pressure_ports..mechanical_ports]
                    .to_vec(),
                work,
            },
            mean_heat_source_temperature_ratio: means[mechanical_ports..].to_vec(),
            heat_source_work_by_source_j,
        })
    }
}
