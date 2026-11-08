//! Physical conducting electrodes on the original settled AC thermal residual.

use super::{
    CircuitEnthalpySolve, Cx, DepositionSimulation, HeatAssembly, HeatError,
    SpatialConductorMaterial, SpatialHarmonicInductionControls, SpatialHarmonicThermalCandidate,
    SpatialInduction, SpatialMagneticMaterial, SpatialMagneticRegion, WeightedAssemblyLimits,
    heat_assembly_with_drive,
};
use fs_em::spatial::{
    SpatialMagneticError,
    vector_potential::{
        SpatialVectorMagnetostatics,
        induction::{
            harmonic::SpatialHarmonicInductionResponse,
            terminals::{SpatialTerminalInduction, current::SpatialCurrentInduction},
        },
    },
};
use fs_math::c64::C64;

// Only this choice of native physical excitation differs. No electrical space,
// field solve, thermal residual or Joule integration is implemented here.
pub(super) enum HarmonicHeatingDrive<'s, 'm, C: SpatialConductorMaterial> {
    Isolated(&'s SpatialInduction<'m, C>),
    Voltage {
        source: &'s SpatialTerminalInduction<'m, C>,
        voltage_v: &'s [C64],
    },
    Current {
        source: &'s SpatialCurrentInduction<'m, C>,
        incoming_current_a: &'s [C64],
    },
}

impl<'m, C: SpatialConductorMaterial> HarmonicHeatingDrive<'_, 'm, C> {
    pub(super) fn problem(&self) -> &SpatialVectorMagnetostatics {
        match self {
            Self::Isolated(source) => source.problem(),
            Self::Voltage { source, .. } => source.problem(),
            Self::Current { source, .. } => source.problem(),
        }
    }
    pub(super) fn conductors(&self) -> &[Option<C>] {
        match self {
            Self::Isolated(source) => source.conductors(),
            Self::Voltage { source, .. } => source.conductors(),
            Self::Current { source, .. } => source.conductors(),
        }
    }

    // Observe while the actual source-refreshed owner is alive. The shared
    // caloric path retains only its owned heat, terminal and diagnostic receipt.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn with_materials<R>(
        &self,
        materials: Vec<SpatialMagneticMaterial>,
        conductors: Vec<Option<C>>,
        impressed_current_a_m2: &[[C64; 3]],
        velocity_m_s: &[[f64; 3]],
        omega_rad_s: f64,
        controls: SpatialHarmonicInductionControls,
        assembly: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
        observe: impl FnOnce(&SpatialHarmonicInductionResponse<'_, 'm, C>) -> R,
    ) -> Result<R, SpatialMagneticError> {
        match self {
            Self::Isolated(source) => {
                let field = source
                    .prepare_materials_trial(materials, conductors, assembly, cancelled)?
                    .into_harmonic(cancelled)?;
                let response = field.solve_harmonic(
                    impressed_current_a_m2,
                    velocity_m_s,
                    omega_rad_s,
                    controls,
                    cancelled,
                )?;
                Ok(observe(&response))
            }
            Self::Voltage { source, voltage_v } => {
                let field = source
                    .prepare_harmonic_materials_trial(materials, conductors, assembly, cancelled)?;
                let response = field.solve_harmonic_voltage(
                    voltage_v,
                    impressed_current_a_m2,
                    velocity_m_s,
                    omega_rad_s,
                    controls,
                    cancelled,
                )?;
                Ok(observe(&response))
            }
            Self::Current {
                source,
                incoming_current_a,
            } => {
                let field = source
                    .prepare_harmonic_materials_trial(materials, conductors, assembly, cancelled)?;
                let response = field.solve_harmonic_current(
                    incoming_current_a,
                    impressed_current_a_m2,
                    velocity_m_s,
                    omega_rad_s,
                    controls,
                    cancelled,
                )?;
                Ok(observe(&response))
            }
        }
    }
}

/// Heat finite bodies with source-temperature feedback under held AC voltages.
///
/// Every proposed temperature refreshes the same physical electrode owner's
/// conductivity and solves its full sparse A–phi response. The imposed peak
/// voltages are in original electrode order; large common-mode offsets are
/// removed by the original per-conductor reference. Impressed current and
/// electrode excitation act coherently in the SAME solve. Their separate signed
/// supplies and the actual electrode reactions remain in the accepted report.
///
/// Source identity, fixed permeability, stationary geometry, complete material
/// heat transport and the settled-carrier/slow-thermal approximation are exactly
/// those of [`super::heat_assembly_candidate`]. The native transient owner stays
/// unchanged. This adds no electrode contact impedance, lead field, circuit
/// graph or external source loss; only resolved conductor heat is deposited.
#[allow(clippy::too_many_arguments)]
pub fn voltage_heat_assembly_candidate<
    't,
    A: HeatAssembly + ?Sized,
    C: SpatialConductorMaterial,
>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    template: &SpatialTerminalInduction<'_, C>,
    voltage_v: &[C64],
    impressed_current_a_m2: &[[C64; 3]],
    omega_rad_s: f64,
    regions: &[SpatialMagneticRegion],
    volume_tolerance_m3: f64,
    controls: CircuitEnthalpySolve,
    field_controls: SpatialHarmonicInductionControls,
    assembly: WeightedAssemblyLimits,
    initial_energy_increments_j: &[f64],
) -> Result<SpatialHarmonicThermalCandidate<'t, A::Sample>, HeatError<A::Error>> {
    heat_assembly_with_drive(
        provider,
        cx,
        thermal,
        target_s,
        HarmonicHeatingDrive::Voltage {
            source: template,
            voltage_v,
        },
        impressed_current_a_m2,
        omega_rad_s,
        regions,
        volume_tolerance_m3,
        controls,
        field_controls,
        assembly,
        initial_energy_increments_j,
    )
}

/// Heat finite bodies with source-temperature feedback under total AC currents.
///
/// Peak incoming phasors retain original electrode order and must balance on
/// each conducting component. The original grouped-potential owner solves the
/// voltage and spatial current distribution at every temperature; neither a
/// uniform current profile nor a frozen terminal impedance is imposed. Actual
/// reactions and separate impressed/terminal supplies accompany cellwise heat.
/// All material, thermal and settled-carrier boundaries of
/// [`voltage_heat_assembly_candidate`] apply, including atomic refusal.
#[allow(clippy::too_many_arguments)]
pub fn current_heat_assembly_candidate<
    't,
    A: HeatAssembly + ?Sized,
    C: SpatialConductorMaterial,
>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    template: &SpatialCurrentInduction<'_, C>,
    incoming_current_a: &[C64],
    impressed_current_a_m2: &[[C64; 3]],
    omega_rad_s: f64,
    regions: &[SpatialMagneticRegion],
    volume_tolerance_m3: f64,
    controls: CircuitEnthalpySolve,
    field_controls: SpatialHarmonicInductionControls,
    assembly: WeightedAssemblyLimits,
    initial_energy_increments_j: &[f64],
) -> Result<SpatialHarmonicThermalCandidate<'t, A::Sample>, HeatError<A::Error>> {
    heat_assembly_with_drive(
        provider,
        cx,
        thermal,
        target_s,
        HarmonicHeatingDrive::Current {
            source: template,
            incoming_current_a,
        },
        impressed_current_a_m2,
        omega_rad_s,
        regions,
        volume_tolerance_m3,
        controls,
        field_controls,
        assembly,
        initial_energy_increments_j,
    )
}
