//! Finite, source-bound gas thermosiphons on the original inertial loop.
//!
//! The registered GasEnthalpyLaw owns NASA-9 h(T), T(h), caloric support and the
//! formation-energy datum. Native cells store M*h and native upwind links carry
//! donor enthalpy, never Cp(T)*T or averaged temperature. The same immutable law
//! supplies local rho, mu, k and Cp at its original thermodynamic pressure,
//! composition, material sources and nonthermal query conditions.
//!
//! This is a fixed-reference-mass, low-Mach Boussinesq approximation for a closed
//! gas circuit. Only gravity uses local density. Geometry, inertia and continuity
//! retain rho0; there is no acoustic/compressible inventory, pressure evolution,
//! reaction, condensation or inferred wall-film closure. The common source,
//! density, Reynolds, work, momentum and thermal acceptance gates remain active.

use super::*;
use fs_conduction::transient::deposition::backward_euler::BoundaryStepCoefficients;
use fs_material::thermochemical_gas::{
    ResolvedThermochemicalGasState, energy::enthalpy::GasEnthalpyLaw,
};

/// Explicit gas regime admission beside the common finite-inventory limits.
#[derive(Clone, Copy, Debug)]
pub struct GasLoopLimits {
    pub inventory: FiniteLoopLimits,
    /// Maximum speed/sound speed at the largest absolute old/new/average flow.
    /// Must be finite in (0, 0.3]. It bounds the admitted discrete interval's
    /// endpoint-held source approximation, not unresolved continuous extrema.
    pub maximum_mach: f64,
}

/// Source-backed gas pipe storage and ONE native retained circulation momentum.
///
/// Register this exact source object as each fluid cell's native ThermalLaw.
/// Initial cell mass is rho_reference*V. Warming leaves these reference masses
/// unchanged; it does not silently replace them with local EOS density*V.
pub struct FiniteGasLoop<'a> {
    source: &'a GasEnthalpyLaw,
    reference: ResolvedThermochemicalGasState,
    maximum_mach: f64,
    core: FiniteLoopCore<'a>,
}

fn resolve<E>(
    source: &GasEnthalpyLaw,
    temperature: Temperature,
) -> Result<ResolvedThermochemicalGasState, FiniteLoopError<E>> {
    // Same operation as the native isobaric source owner: replace only the
    // temperature coordinate, retaining pressure, composition and every pin.
    let point = source
        .context()
        .clone()
        .with("temperature", temperature.value())
        .map_err(|_| FiniteLoopError::Invalid("invalid original gas temperature coordinate"))?;
    source
        .model()
        .evaluate(&point)
        .map_err(FiniteLoopError::Transport)
}

fn admit_mach<E>(
    cx: &Cx<'_>,
    sources: &[ResolvedThermochemicalGasState],
    losses: &[ResistanceEvaluation],
    maximum: f64,
) -> Result<(), FiniteLoopError<E>> {
    if sources.len() != losses.len() {
        return Err(FiniteLoopError::Invalid(
            "one original gas state per hydraulic section required",
        ));
    }
    for (section, (source, loss)) in sources.iter().zip(losses).enumerate() {
        poll(cx)?;
        let sound_speed = source.state().sound_speed;
        let speed = loss.volume_flow_m3_s.abs() / loss.law.flow_area().value();
        let mach = speed / sound_speed;
        if !positive(sound_speed) || !nonnegative(mach) || mach > maximum {
            return Err(FiniteLoopError::GasMach {
                section,
                value: mach,
                maximum,
            });
        }
    }
    Ok(())
}

impl<'a> FiniteGasLoop<'a> {
    /// Bind source pressure/composition and actual finite-cell caloric identity.
    /// Reference density is resolved from the law's unchanged original context.
    /// Unsupported reference/initial conditions refuse before any owner exists.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        cx: &Cx<'_>,
        source: &'a GasEnthalpyLaw,
        thermal: &DepositionSimulation<'_>,
        sections: &[PipeCell],
        gravity: SeriesGravity,
        initial_volume_flow_m3_s: f64,
        limits: GasLoopLimits,
    ) -> Result<Self, FiniteLoopError> {
        poll(cx)?;
        if !positive(limits.maximum_mach) || limits.maximum_mach > 0.3 {
            return Err(FiniteLoopError::Invalid(
                "gas circulation requires an explicit Mach limit in (0,0.3]",
            ));
        }
        let reference = source
            .model()
            .evaluate(source.context())
            .map_err(FiniteLoopError::Transport)?;
        // A context outside this law's admitted caloric bracket is not a valid
        // inertia/continuity reference, even if the underlying source is broader.
        source
            .enthalpy_at_temperature(reference.state().temperature)
            .map_err(FiniteLoopError::GasCaloric)?;
        let transport = reference
            .transport_state()
            .map_err(FiniteLoopError::Transport)?;
        let core = FiniteLoopCore::try_new(
            cx,
            source,
            transport,
            thermal,
            sections,
            gravity,
            initial_volume_flow_m3_s,
            limits.inventory,
        )?;
        let (sources, transports) =
            core.resolve(cx, thermal, &|_, t| resolve::<Infallible>(source, t))?;
        let losses = core.losses::<Infallible>(cx, &transports, initial_volume_flow_m3_s)?;
        admit_mach::<Infallible>(cx, &sources, &losses, limits.maximum_mach)?;
        Ok(Self {
            source,
            reference,
            maximum_mach: limits.maximum_mach,
            core,
        })
    }

    pub const fn source(&self) -> &'a GasEnthalpyLaw {
        self.source
    }
    /// Original source-bearing state used to fix inertia and reference mass.
    pub const fn reference(&self) -> &ResolvedThermochemicalGasState {
        &self.reference
    }
    pub const fn reference_density_kg_m3(&self) -> f64 {
        self.core.reference_density
    }
    pub const fn inertance_pa_s2_m3(&self) -> f64 {
        self.core.inertance
    }
    pub const fn accepted_time_s(&self) -> f64 {
        self.core.accepted_time_s
    }
    pub fn volume_flow_m3_s(&self) -> f64 {
        self.core.volume_flow_m3_s()
    }
    pub fn kinetic_energy_j(&self) -> f64 {
        self.core.kinetic_energy_j()
    }

    /// Advance finite gas enthalpy and signed inertial circulation atomically.
    ///
    /// Every trial resolves the same NASA/transport owner at native T(h).
    /// Returned snapshots retain original pressure, species and source receipts.
    /// Heat callbacks see actual signed-average flow and may be repeated by the
    /// common bounded nonlinear solve. They must be pure physical proposals.
    ///
    /// At rest the original heat network still evolves pipe inventory. During
    /// reversal native signed-average donor transport remains unchanged, while
    /// Mach/Reynolds and omitted-work admission use the endpoint-speed envelope.
    /// Source, domain, regime, budget or cancellation failure publishes neither
    /// pipe momentum nor any native thermal cell/history.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'m, T, E>(
        &mut self,
        cx: &Cx<'_>,
        thermal: &mut DepositionSimulation<'m>,
        target_s: f64,
        driving_pressure_pa: f64,
        heat: impl Fn(
            &DepositionSimulation<'m>,
            &FiniteLoopHydraulics<ResolvedThermochemicalGasState>,
        ) -> Result<BoundaryStepCoefficients<T>, E>,
        controls: FiniteLoopControls,
        energy_guess_j: &[f64],
    ) -> Result<FiniteLoopStep<T, ResolvedThermochemicalGasState>, FiniteLoopError<E>> {
        let source = self.source;
        let maximum = self.maximum_mach;
        self.core.advance(
            cx,
            thermal,
            target_s,
            driving_pressure_pa,
            heat,
            |_, t| resolve(source, t),
            |hydraulic| admit_mach(cx, &hydraulic.sources, &hydraulic.envelope_losses, maximum),
            controls,
            energy_guess_j,
        )
    }
}
