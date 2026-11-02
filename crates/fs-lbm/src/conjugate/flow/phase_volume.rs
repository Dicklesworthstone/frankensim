//! Isothermal phase-volume evolution on the original conservative flux owner.

use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState};

use super::density_transport::{admit, advect_projected_material};
use super::{
    AdvectedMaterial, CancelGate, ChtError, FlowFace, FlowField, MaterialAdvectionControls,
    MaterialAdvectionInput, MaterialAdvectionStep, MaterialVolumeSource, VoxelDomain, finite,
    finite_positive, poll, poll_every,
};

/// Original material inventories and externally determined caloric increments.
#[derive(Clone, Copy, Debug)]
pub struct CoexistenceVolumeInput<'a> {
    /// Authoritative cell mass [kg], positive in fluids and zero in solids.
    pub mass_kg: &'a [f64],
    /// Original chart-resolved states; solid grid cells must contain None.
    pub phase_states: &'a [Option<EquilibriumPhaseState>],
    /// Prescribed NET material enthalpy increment [J] before transport.
    /// This is not a raw heat input: the caller owns its thermal/work closure.
    /// Solid grid cells require exactly zero.
    pub net_enthalpy_increment_j: &'a [f64],
    /// Incoming/backflow states on x-, x+, y-, y+, z-, z+ boundaries.
    /// Their density is obtained from the same original coexistence chart.
    pub boundary_inflow: [Option<EquilibriumPhaseState>; 6],
}

/// Existing pressure/transport controls plus actual material-volume admission.
#[derive(Clone, Copy, Debug)]
pub struct CoexistenceVolumeControls {
    /// Native pressure and conservative mass/enthalpy transport controls.
    /// Its enthalpy budget also includes caloric-increment representation error.
    pub transport: MaterialAdvectionControls,
    /// Maximum absolute M/rho_chart(h) - cell volume at either endpoint [m3].
    /// This independently checks the actual returned inventory and phase state.
    pub maximum_cell_volume_defect_m3: f64,
}

/// Caloric, occupied-volume and separately required mechanical source work.
#[derive(Clone, Debug, PartialEq)]
pub struct CoexistenceVolumeReport {
    /// Sum of the prescribed net material enthalpy increments [J].
    pub requested_enthalpy_increment_j: f64,
    /// Sum of M*(h_star-h_old), with the actual representable h_star [J].
    pub represented_enthalpy_increment_j: f64,
    /// Sum of absolute requested-versus-represented increment defects [J].
    pub enthalpy_increment_defect_j: f64,
    /// Actual total material enthalpy change, including carried enthalpy [J].
    pub material_enthalpy_change_j: f64,
    /// Increment representation error plus the native transport defect [J].
    pub absolute_enthalpy_defect_j: f64,
    /// Largest absolute original mass/chart occupied-volume defect [m3].
    pub maximum_initial_volume_defect_m3: f64,
    /// Largest absolute returned mass/chart occupied-volume defect [m3].
    pub maximum_final_volume_defect_m3: f64,
    /// Additional dt*sum(p_c*S_c) required by the pressure projection [J].
    /// This is distinct from the prescribed net material enthalpy increment;
    /// the external energy owner must supply/account for it. Numerical
    /// projection loss is retained separately in the native projection report.
    /// Uses correction pressure in the native gauge; fixed-boundary pressure
    /// work remains in transport.projection.report and must be included there.
    pub required_source_pressure_work_j: f64,
}

/// One atomically proposed phase, material and corrected-flow endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct CoexistenceVolumeStep {
    /// Actual corrected flow and conservative material inventory.
    pub transport: MaterialAdvectionStep,
    /// Original-chart phase states evaluated at the returned specific h.
    pub phase_states: Vec<Option<EquilibriumPhaseState>>,
    /// Net outward volume rate supplied to the native pressure constraint.
    /// These thermodynamic volume sources contain NO injected material mass.
    pub volume_source_m3_s: Vec<f64>,
    /// Independently admitted caloric, material-volume and source-work account.
    pub report: CoexistenceVolumeReport,
}

impl FlowField {
    /// Evolve a common-velocity isothermal coexistence mixture on one original
    /// adjacent knot interval, retaining its mass, enthalpy and occupied volume.
    ///
    /// The selected knots must have exactly equal temperature and increasing
    /// liquid mass fraction. The native chart then has v(h)=a+b*h in this
    /// interval, including charts with a different sensible-density rule.
    /// Each original M,h must satisfy M*v(h)=V within the volume budget.
    /// A prescribed NET caloric increment gives h_star=h+delta_H/M and
    /// S=M*(v(h_star)-v(h))/dt. S enters pressure continuity with ZERO local
    /// material mass source. The existing implicit donor owner transports M
    /// and M*h together on those actual corrected face fluxes.
    ///
    /// At an exact solve its shared operator A obeys
    /// A*(a*rho+b*rho*h)=V+dt*S=A*1, including incoming boundary states from
    /// the same interval. Thus the endpoint retains rho*v(h)=1. Actual returned
    /// masses and chart states independently pass the volume gate; no density,
    /// mass or phase fraction is repaired to force this identity.
    /// A tolerated initial occupied-volume residual is propagated and checked
    /// again, never converted into an artificial physical volume source. For
    /// exact initial closure, the source is equivalently (M*v(h_star)-V)/dt.
    ///
    /// This is a kinematic caloric/continuity substep for an explicitly mobile
    /// single-velocity mixture. The caller supplies the momentum predictor and
    /// its admissible phase mobility, and accounts for the reported additional
    /// source pressure work. The net enthalpy increments are not raw heat.
    /// The original chart's pressure, composition and source applicability stay
    /// fixed; no EOS pressure feedback, solid stress, mushy drag, nucleation,
    /// thermal solve or free-surface geometry is inferred. Incoming enthalpy
    /// already follows its declared pressure-flow convention. Projection loss
    /// remains numerical and is never deposited as latent or sensible heat.
    ///
    /// All cells and incoming states must remain in this ONE interval. Sensible
    /// branch crossing, another chart, incompatible sealed expansion, a failed
    /// native solve, exceeded budgets or cancellation returns no candidate.
    /// The original flow, inventories and material states remain unchanged.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance_coexistence_volume(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        curve: &EquilibriumEnthalpyPhaseCurve,
        lower_knot: usize,
        input: CoexistenceVolumeInput<'_>,
        controls: CoexistenceVolumeControls,
        gate: &CancelGate,
    ) -> Result<CoexistenceVolumeStep, ChtError> {
        poll(gate)?;
        let n = domain.cell_count();
        if n > controls.transport.max_cells
            || input.mass_kg.len() != n
            || input.phase_states.len() != n
            || input.net_enthalpy_increment_j.len() != n
        {
            return Err(refuse("phase inventory shape or cell budget is invalid"));
        }
        let dt = controls.transport.projection.duration_s;
        finite_positive("flow.coexistence_duration", dt)?;
        let volume_budget = controls.maximum_cell_volume_defect_m3;
        if !volume_budget.is_finite() || volume_budget < 0.0 {
            return Err(refuse("invalid cell-volume defect budget"));
        }
        let upper_knot = lower_knot
            .checked_add(1)
            .ok_or_else(|| refuse("coexistence interval index overflow"))?;
        let (lower, upper) = curve
            .knots()
            .get(lower_knot)
            .zip(curve.knots().get(upper_knot))
            .ok_or_else(|| refuse("one adjacent original knot interval required"))?;
        if lower.temperature_k != upper.temperature_k
            || lower.liquid_mass_fraction >= upper.liquid_mass_fraction
        {
            return Err(refuse("selected interval must be isothermal coexistence"));
        }
        let bounds = [lower.specific_enthalpy_j_kg, upper.specific_enthalpy_j_kg];
        finite_positive("flow.coexistence_enthalpy_span", bounds[1] - bounds[0])?;
        let mut boundary_inflow = [None; 6];
        for (output, supplied) in boundary_inflow.iter_mut().zip(input.boundary_inflow) {
            if let Some(state) = supplied {
                admit_phase(curve, bounds, state)?;
                *output = Some(AdvectedMaterial {
                    density_kg_m3: state.bulk_density_kg_m3(),
                    specific_enthalpy_j_kg: state.specific_enthalpy_j_kg(),
                });
            }
        }
        let sources = vec![MaterialVolumeSource::None; n];
        let mut h_star = vec![0.0; n];
        let mut volume_source = vec![0.0; n];
        let mut requested = Sum::default();
        let mut represented = Sum::default();
        let mut increment_defect = Sum::default();
        let mut initial_defect = 0.0_f64;
        for c in 0..n {
            poll_every(gate, c)?;
            let delta_h_total = input.net_enthalpy_increment_j[c];
            finite("flow.net_enthalpy_increment_j", delta_h_total)?;
            if !domain.is_fluid(c) {
                if input.phase_states[c].is_some()
                    || delta_h_total != 0.0
                    || input.mass_kg[c] != 0.0
                {
                    return Err(refuse(
                        "solid grid cells cannot contain phase inventory or caloric increments",
                    ));
                }
                continue;
            }
            let old = input.phase_states[c]
                .ok_or_else(|| refuse("every fluid cell requires an original phase state"))?;
            admit_phase(curve, bounds, old)?;
            let mass = input.mass_kg[c];
            finite_positive("flow.coexistence_mass", mass)?;
            let old_volume = mass / old.bulk_density_kg_m3();
            finite_positive("flow.coexistence_initial_volume", old_volume)?;
            let volume_error = old_volume - domain.volume(c);
            finite("flow.coexistence_initial_volume_defect", volume_error)?;
            initial_defect = initial_defect.max(volume_error.abs());
            if volume_error.abs() > volume_budget {
                return Err(refuse(
                    "initial inventory and phase density do not fill the original cell",
                ));
            }
            let dh = delta_h_total / mass;
            finite("flow.coexistence_specific_increment", dh)?;
            if delta_h_total != 0.0 && dh == 0.0 {
                return Err(refuse("net specific-enthalpy increment underflowed"));
            }
            let next = curve
                .advance_specific_energy(old, dh)
                .map_err(|error| refuse(error.to_string()))?;
            admit_phase(curve, bounds, next)?;
            let accepted_dh = next.specific_enthalpy_j_kg() - old.specific_enthalpy_j_kg();
            if delta_h_total != 0.0 && accepted_dh == 0.0 {
                return Err(refuse(
                    "net enthalpy increment is lost at the caloric datum",
                ));
            }
            let accepted_delta_h = mass * accepted_dh;
            finite("flow.coexistence_represented_increment", accepted_delta_h)?;
            requested.add(delta_h_total);
            represented.add(accepted_delta_h);
            increment_defect.add((delta_h_total - accepted_delta_h).abs());
            h_star[c] = next.specific_enthalpy_j_kg();
            let next_volume = mass / next.bulk_density_kg_m3();
            finite_positive("flow.coexistence_source_volume", next_volume)?;
            // Stable M*(1/rho_star-1/rho_old), retaining small phase changes.
            let density_change = old.bulk_density_kg_m3() - next.bulk_density_kg_m3();
            let relative = density_change / old.bulk_density_kg_m3().max(next.bulk_density_kg_m3());
            let delta_volume = old_volume.max(next_volume) * relative;
            finite("flow.coexistence_volume_increment", delta_volume)?;
            if density_change != 0.0 && delta_volume == 0.0 {
                return Err(refuse("phase-volume increment underflowed"));
            }
            // An admitted initial EOS roundoff defect is not physical phase
            // expansion; carry it to the independent final volume check.
            volume_source[c] = delta_volume / dt;
            finite("flow.coexistence_volume_source", volume_source[c])?;
        }
        for value in [requested.value, represented.value, increment_defect.value] {
            finite("flow.coexistence_caloric_source", value)?;
        }
        if increment_defect.value > controls.transport.maximum_enthalpy_defect_j {
            return Err(refuse(
                "represented net enthalpy increments exceed the energy-defect budget",
            ));
        }
        let transport_input = MaterialAdvectionInput {
            mass_kg: input.mass_kg,
            specific_enthalpy_j_kg: &h_star,
            volume_sources: &sources,
            boundary_inflow,
        };
        let density = admit(domain, transport_input, controls.transport, gate)?;
        let projected = self.project_density_volume_source(
            domain,
            faces,
            &density,
            &volume_source,
            controls.transport.projection,
            gate,
        )?;
        let transport = advect_projected_material(
            domain,
            projected,
            transport_input,
            controls.transport,
            gate,
        )?;
        let mut phase_states = vec![None; n];
        let mut final_defect = 0.0_f64;
        for c in 0..n {
            poll_every(gate, c)?;
            if !domain.is_fluid(c) {
                continue;
            }
            let state = curve
                .state_at_specific_enthalpy(transport.specific_enthalpy_j_kg[c])
                .map_err(|error| refuse(error.to_string()))?;
            admit_phase(curve, bounds, state)?;
            let volume = transport.mass_kg[c] / state.bulk_density_kg_m3();
            finite_positive("flow.coexistence_final_volume", volume)?;
            let defect = volume - domain.volume(c);
            finite("flow.coexistence_final_volume_defect", defect)?;
            final_defect = final_defect.max(defect.abs());
            if defect.abs() > volume_budget {
                return Err(refuse(
                    "transported mass and phase density exceed the cell-volume defect budget",
                ));
            }
            phase_states[c] = Some(state);
        }
        let absolute_enthalpy_defect_j =
            increment_defect.value + transport.report.absolute_enthalpy_defect_j;
        finite(
            "flow.coexistence_total_enthalpy_defect",
            absolute_enthalpy_defect_j,
        )?;
        if absolute_enthalpy_defect_j > controls.transport.maximum_enthalpy_defect_j {
            return Err(refuse(
                "caloric source and transport exceed the combined enthalpy-defect budget",
            ));
        }
        let material_enthalpy_change_j = represented.value + transport.report.enthalpy_change_j;
        finite(
            "flow.coexistence_enthalpy_change",
            material_enthalpy_change_j,
        )?;
        let report = CoexistenceVolumeReport {
            requested_enthalpy_increment_j: requested.value,
            represented_enthalpy_increment_j: represented.value,
            enthalpy_increment_defect_j: increment_defect.value,
            material_enthalpy_change_j,
            absolute_enthalpy_defect_j,
            maximum_initial_volume_defect_m3: initial_defect,
            maximum_final_volume_defect_m3: final_defect,
            required_source_pressure_work_j: transport.projection.report.pressure_source_work_j,
        };
        poll(gate)?;
        Ok(CoexistenceVolumeStep {
            transport,
            phase_states,
            volume_source_m3_s: volume_source,
            report,
        })
    }
}

fn admit_phase(
    curve: &EquilibriumEnthalpyPhaseCurve,
    bounds: [f64; 2],
    state: EquilibriumPhaseState,
) -> Result<(), ChtError> {
    let h = state.specific_enthalpy_j_kg();
    if state.phase_curve_identity() != curve.identity()
        || state.material_card_identity() != curve.material_card_identity()
        || h < bounds[0]
        || h > bounds[1]
    {
        return Err(refuse(
            "phase state must retain the original chart and selected coexistence interval",
        ));
    }
    Ok(())
}

#[derive(Default)]
struct Sum {
    value: f64,
    compensation: f64,
}

impl Sum {
    fn add(&mut self, value: f64) {
        let corrected = value - self.compensation;
        let next = self.value + corrected;
        self.compensation = (next - self.value) - corrected;
        self.value = next;
    }
}

fn refuse(reason: impl Into<String>) -> ChtError {
    ChtError::InvalidInput {
        field: "flow.coexistence_volume",
        reason: reason.into(),
    }
}

#[cfg(test)]
mod tests;
