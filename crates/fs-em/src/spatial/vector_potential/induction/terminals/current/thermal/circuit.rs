//! Loaded AC circuits on an admitted source-temperature electrode trial.

use super::SpatialHarmonicCurrentMaterialTrial;
use crate::spatial::vector_potential::induction::harmonic::terminals::{
    HarmonicCircuitControls, SpatialHarmonicTerminalCircuitResponse,
    SpatialHarmonicTerminalControls, SpatialHarmonicTerminalError,
    SpatialHarmonicTerminalFrequency, SpatialHarmonicTerminalPort, SpatialHarmonicTerminalResponse,
};
use fs_material::conductor::{ResolvedConductorState, spatial::SpatialConductorMaterial};
use fs_math::c64::C64;
use fs_phs::circuit::{CircuitGraph, MutualInductance};
use fs_solver::SolveReport;

type Result<T> = core::result::Result<T, SpatialHarmonicTerminalError>;

/// Settled AC port responses borrowing one original material-temperature trial.
///
/// The private native frequency owner retains the complete coupled impedance,
/// independent impressed excitation and cached cochains. Only harmonic responses
/// are exposed: this view cannot yield a refreshed transient owner or advance
/// magnetic history without its parameter-work account.
pub struct SpatialHarmonicCurrentMaterialFrequency<
    's,
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    frequency: SpatialHarmonicTerminalFrequency<'s, 'm, C>,
}

impl<'m, C: SpatialConductorMaterial> SpatialHarmonicCurrentMaterialTrial<'m, C> {
    /// Prepare the original physical electrode impedance at this temperature.
    /// Material operators have already been admitted and are not prepared again.
    /// Native unit-current probes retain zero impressed excitation; a nonzero
    /// independent source adds one open-circuit solve and its affine voltage.
    /// Every later current or circuit response recovers one actual native field.
    ///
    /// A caller bounding total field work must reserve one trial per independent
    /// electrode port, one for a nonzero impressed source, and one per requested
    /// actual response before starting preparation. Native field/interface
    /// storage and solver allowances remain in `controls`.
    pub fn prepare_harmonic_frequency(
        &self,
        impressed_current_a_m2: &[[C64; 3]],
        omega_rad_s: f64,
        velocity_m_s: &[[f64; 3]],
        controls: SpatialHarmonicTerminalControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicCurrentMaterialFrequency<'_, 'm, C>> {
        let frequency = SpatialHarmonicTerminalFrequency::with_impressed_current(
            &self.source,
            impressed_current_a_m2,
            omega_rad_s,
            velocity_m_s,
            controls,
            cancelled,
        )?;
        Ok(SpatialHarmonicCurrentMaterialFrequency { frequency })
    }
}

impl<'s, 'm, C: SpatialConductorMaterial> SpatialHarmonicCurrentMaterialFrequency<'s, 'm, C> {
    /// Original non-reference electrodes and their physical return electrodes.
    pub fn ports(&self) -> &[SpatialHarmonicTerminalPort] {
        self.frequency.ports()
    }

    /// Complete row-major impedance at this admitted source temperature [ohm].
    pub fn impedance_ohm(&self) -> &[C64] {
        self.frequency.impedance_ohm()
    }

    /// Original induced open-circuit voltage in physical port order [V].
    pub fn open_circuit_voltages_v(&self) -> &[C64] {
        self.frequency.open_circuit_voltages_v()
    }

    pub fn unit_reports(&self) -> &[SolveReport] {
        self.frequency.unit_reports()
    }

    pub fn impressed_report(&self) -> Option<&SolveReport> {
        self.frequency.impressed_report()
    }

    /// Actual native field solves used by this completed preparation.
    /// Each subsequent current/circuit response adds one native field solve,
    /// including when its cached seed needs no further Krylov iterations.
    pub fn preparation_field_trials(&self) -> usize {
        self.unit_reports().len() + usize::from(self.impressed_report().is_some())
    }

    /// Recover the actual combined impressed and port-current field once.
    /// Original material, charge, Ampere and complex-power checks remain native.
    pub fn current_response(
        &self,
        currents_a: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicTerminalResponse<'s, 'm, C>> {
        self.frequency.current_response(currents_a, cancelled)
    }

    /// Close the original network and recover its actual loaded field once.
    /// Numerical reference resistors cancel; graph resistor heat, conductor
    /// heat, impressed supply and graph supply retain their original accounts.
    #[allow(clippy::too_many_arguments)]
    pub fn circuit_response(
        &self,
        graph: &CircuitGraph,
        terminal_branches: &[usize],
        external_mutual: &[MutualInductance],
        source_phasors: &[C64],
        controls: HarmonicCircuitControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHarmonicTerminalCircuitResponse<'s, 'm, C>> {
        self.frequency.circuit_response(
            graph,
            terminal_branches,
            external_mutual,
            source_phasors,
            controls,
            cancelled,
        )
    }
}
