//! A nonlinear instantaneous field law and the original passive Debye arms.
//!
//! q_inf=A(D-sum P_k), q_k=A P_k make the physical fixed-geometry potential
//! H=A*d*[w_inf(D-sum P_k)+sum P_k^2/(2 chi_k)] separable in the SAME native
//! branch charges. The existing RC arms give tau_k Pdot_k+P_k=chi_k E and
//! positive heat A*d*sum(tau_k/chi_k)*Pdot_k^2. No relaxation fit, history
//! reset, extra integrator or frequency-dependent energy formula is used.

use super::*;

fn admit_source(
    source: &ResolvedDebyeCapacitor,
    instantaneous: &MaterialDielectricCapacitor,
) -> Result<(), MaterialCircuitError> {
    let original = source.material();
    let field = instantaneous.material().material();
    let axis = instantaneous.material().source_path().axis.as_str();
    if original.card_identity() != field.card_identity()
        || source.area_m2().to_bits() != instantaneous.area_m2().to_bits()
        || source.gap_m().to_bits() != instantaneous.gap_m().to_bits()
        || original
            .query_point()
            .iter()
            .filter(|(name, _)| name != axis)
            .ne(field.query_point().iter().filter(|(name, _)| name != axis))
        || original
            .axis_quantities()
            .iter()
            .filter(|(name, _)| name.as_str() != axis)
            .ne(field
                .axis_quantities()
                .iter()
                .filter(|(name, _)| name.as_str() != axis))
    {
        return Err(invalid(
            "instantaneous and Debye laws require the same original card, geometry and non-field context",
        ));
    }
    // The bridge already admits field-only instantaneous storage. Constant
    // source increments keep this closed history's polarization potential free
    // of missing temperature/other parameter work; tau may be a fixed-context
    // kinetic snapshot, since this owner never changes that context.
    for name in std::iter::once(source.high_frequency_property()).chain(
        source
            .pole_properties()
            .iter()
            .map(|p| p.permittivity_increment.as_str()),
    ) {
        let property = original
            .property(name)
            .ok_or_else(|| invalid("original Debye energy source absent"))?;
        if !property.source_axes().is_empty()
            || property.answer().receipt.decision
                != fs_matdb::EvaluationDecision::ConstantWithinValidity
        {
            return Err(invalid(
                "closed nonlinear Debye storage requires constant original background and pole increments",
            ));
        }
    }
    let declared = original
        .property(source.high_frequency_property())
        .ok_or_else(|| invalid("original instantaneous permittivity absent"))?
        .value_si();
    let actual = instantaneous
        .material()
        .evaluate(0.0)?
        .differential_permittivity_f_m;
    if (actual - declared).abs() > 32.0 * f64::EPSILON * actual.max(declared) {
        return Err(invalid(
            "nonlinear instantaneous zero-field tangent contradicts the original Debye background",
        ));
    }
    Ok(())
}

impl DebyeCircuit {
    /// Assemble a NEW nonlinear instantaneous dielectric with original Debye
    /// polarization histories. The supplied D(E) law explicitly represents the
    /// instantaneous response; feeding an equilibrium total D(E) would count
    /// the retained polarization again and is not the declared physical model.
    ///
    /// Both source owners must retain the same card, non-E context and geometry,
    /// and the zero-field tangent must match the declared epsilon_infinity to
    /// arithmetic roundoff. Each original charge/flux and polarization charge
    /// enters unchanged. The active entry in `charges_c` is q_inf, not total
    /// terminal charge; `instantaneous.branch()` selects the original branch.
    ///
    /// Geometry is fixed. Source changes require the explicit conservative
    /// reversible-source rebind, which retains every charge and returns its
    /// parameter work. The finite thermal owner admits temperature-independent
    /// field/polarization potentials with supported tau(T); actual temperature-
    /// dependent potentials need their own field entropy and remain refused.
    /// Native capacitor efforts close each original store's energy change;
    /// nonlinear work is not redistributed into another changing capacitor.
    /// Whole-circuit and dielectric terminal-work audits remain strict, and a
    /// failed account never becomes extra polarization heat. Resolve relaxation
    /// and source knots and check time-step refinement; discrete energy closure
    /// is not a temporal error bound.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_reversible_instantaneous(
        graph: CircuitGraph,
        source: ResolvedDebyeCapacitor,
        instantaneous: MaterialDielectricCapacitor,
        fluxes_wb: Vec<f64>,
        charges_c: Vec<f64>,
        polarization_charges_c: Vec<f64>,
        initial_inputs: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MaterialCircuitError> {
        poll(&mut cancelled)?;
        admit_source(&source, &instantaneous)?;
        Self::assemble_with_instantaneous(
            graph,
            instantaneous.branch(),
            source,
            Some(instantaneous),
            fluxes_wb,
            charges_c,
            polarization_charges_c,
            initial_inputs,
            cancelled,
        )
    }

    /// Source-backed instantaneous law, present only for the nonlinear profile.
    /// Ordinary Debye constructors keep their original linear background.
    pub fn reversible_instantaneous(&self) -> Option<&MaterialDielectricCapacitor> {
        self.prepared.reversible_instantaneous.as_ref()
    }

    /// Total free charge at the dielectric's original terminal orientation [C].
    /// Internal polarization charges are counted once alongside q_inf.
    pub fn terminal_charge_c(&self) -> Result<f64, MaterialCircuitError> {
        let (mut sum, mut correction) = (self.instantaneous_charge_c(), 0.0);
        for &slot in &self.prepared.polarization_charge_indices {
            accumulate(&mut sum, &mut correction, self.state[slot])?;
        }
        finite(sum + correction)
    }
}
