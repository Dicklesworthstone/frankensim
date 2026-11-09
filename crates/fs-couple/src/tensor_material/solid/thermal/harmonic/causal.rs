//! Atomic original material-memory and caloric candidates for harmonic forcing.

use super::{
    C64, CircuitEnthalpyError, CircuitEnthalpySolve, Cx, Dd, DepositionSimulation,
    ElementThermalCell, Error, HarmonicSolidHeatCandidate, HarmonicSolidHeatReport,
    HeatAssembly, ViscoHarmonicMaterialSpectrumField, candidate_from_observation, finite, poll,
};
use fs_solid::linear3::visco::{
    ViscoHarmonicMaterialCausalObserver, ViscoHarmonicMaterialHistory,
};

/// Both unpublished endpoints. Accepting them together preserves the original
/// constitutive memory, physical clock and thermal heat handoff exactly once.
pub struct CausalHarmonicSolidHeatCandidate<'s, 't, S> {
    history: ViscoHarmonicMaterialHistory<'s>,
    thermal: HarmonicSolidHeatCandidate<'t, S>,
}

impl<'s, 't, S> CausalHarmonicSolidHeatCandidate<'s, 't, S> {
    pub const fn history(&self) -> &ViscoHarmonicMaterialHistory<'s> {
        &self.history
    }
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        self.thermal.thermal()
    }
    /// Original causal constitutive storage/work/heat and native caloric audit.
    pub fn report(&self) -> &HarmonicSolidHeatReport<S> {
        self.thermal.report()
    }
    pub fn into_parts(self) -> (
        ViscoHarmonicMaterialHistory<'s>,
        DepositionSimulation<'t>,
        HarmonicSolidHeatReport<S>,
    ) {
        let (thermal, report) = self.thermal.into_parts();
        (self.history, thermal, report)
    }
}

/// Advance retained Maxwell memory and deposit its actual finite-interval heat.
///
/// The displacement banks and static bias use ALL original node-major xyz rows,
/// including prescribed supports, as required by the causal solid observer.
/// Initial nodal motion and material strain must join the prescribed waveform
/// continuously. Original histories are borrowed and neither accepted state is
/// changed by preparation, a failed thermal trial, or cancellation.
///
/// The full material source, temperature independence, reference mass/volume,
/// receiving-cell ownership, nonlinear transport and energy gates are identical
/// to the settled harmonic candidate. The causal material observer runs once
/// before thermal iteration; its actual Kelvin/Maxwell heat is deposited once.
/// Retained spring energy can relax into heat even with no AC drive. Stress work,
/// direct storage change and caloric storage remain independently accounted.
///
/// The original thermal owner has an f64 clock. Its start must equal the full
/// DD material clock, and the f64 duration must reach the target exactly in DD.
/// This admits ordinary adjacent thermal steps while explicitly refusing an
/// interval whose endpoint would silently lose material time. It does not round
/// or reset a retained clock. Standalone causal material evolution can continue
/// with finer DD intervals. Temperature feedback into material coefficients,
/// held thermal eigenstrain, mechanical force balance and startup-constrained
/// plane stress remain outside this physical handoff.
#[allow(clippy::too_many_arguments)]
pub fn causal_heat_assembly_candidate<'s, 't, A: HeatAssembly + ?Sized>(
    provider: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    observer: &ViscoHarmonicMaterialCausalObserver<'s>,
    accepted: &ViscoHarmonicMaterialHistory<'s>,
    full_displacements: &[&[C64]],
    full_bias_m: &[f64],
    cells: &[ElementThermalCell<'_>],
    controls: CircuitEnthalpySolve,
    initial_energy_increments_j: &[f64],
) -> Result<CausalHarmonicSolidHeatCandidate<'s, 't, A::Sample>, Error<A::Error>> {
    poll(cx)?;
    let duration = finite(target_s - thermal.time_s())?;
    if duration <= 0.0
        || !target_s.is_finite()
        || accepted.time_s_dd() != Dd::from_f64(thermal.time_s())
        || accepted.time_s_dd() + Dd::from_f64(duration) != Dd::from_f64(target_s)
    {
        return Err(CircuitEnthalpyError::Invalid(
            "causal material and thermal endpoints must share the complete physical clock",
        ));
    }
    let mut next_history = None;
    let thermal = candidate_from_observation(
        provider,
        cx,
        thermal,
        target_s,
        observer.source(),
        cells,
        controls,
        initial_energy_increments_j,
        |window| {
            let step = observer.advance_window(
                accepted, full_displacements, full_bias_m, window.duration_s(), cx,
            )?;
            let field = step.field;
            next_history = Some(step.history);
            Ok(ViscoHarmonicMaterialSpectrumField {
                window,
                angular_frequencies_rad_s: field.angular_frequencies_rad_s,
                points: field.points,
                start_stored_energy_j: field.start_stored_energy_j,
                end_stored_energy_j: field.end_stored_energy_j,
                stored_energy_change_j: field.stored_energy_change_j,
                kelvin_heat_j: field.kelvin_heat_j,
                maxwell_heat_j: field.maxwell_heat_j,
                material_heat_j: field.material_heat_j,
                stress_work_j: field.stress_work_j,
                energy_balance_residual_j: field.energy_balance_residual_j,
            })
        },
    )?;
    let history = next_history.ok_or(CircuitEnthalpyError::Invalid(
        "causal material observation did not produce an endpoint",
    ))?;
    poll(cx)?;
    Ok(CausalHarmonicSolidHeatCandidate { history, thermal })
}
