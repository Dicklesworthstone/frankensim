//! Original-chart thermal volume evolution on the conservative flux owner.

mod conduction;
pub use conduction::{
    ConductiveThermalVolumeControls, ConductiveThermalVolumeMomentumControls,
    ConductiveThermalVolumeMomentumStep, ConductiveThermalVolumeReport,
    ConductiveThermalVolumeStep, ConductiveViscousThermalVolumeMomentumStep,
};
mod momentum;
pub use momentum::{
    ResistiveThermalVolumeMomentumControls, ResistiveThermalVolumeMomentumStep,
    ResolvedViscousThermalVolumeMomentumStep, ThermalResistanceCouplingReport,
    ThermalVolumeMomentumControls, ThermalVolumeMomentumInput, ThermalVolumeMomentumStep,
    ViscousThermalVolumeMomentumStep,
};

use fs_material::phase::{
    EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, PhaseIntervalDensityLaw, PhaseStateError,
};

use super::density_transport::{admit, advect_projected_material};
use super::{
    AdvectedMaterial, CancelGate, ChtError, DensityProjection, FlowFace, FlowField,
    MaterialAdvectionControls, MaterialAdvectionInput, MaterialAdvectionStep, MaterialVolumeSource,
    VoxelDomain, finite, finite_positive, poll, poll_every,
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
    /// Their density is obtained from the same admitted original chart range.
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

/// Budgets for sensible or latent expansion on one original chart interval.
#[derive(Clone, Copy, Debug)]
pub struct ThermalVolumeControls {
    /// Existing pressure, mass, enthalpy and actual occupied-volume budgets.
    pub volume: CoexistenceVolumeControls,
    /// Maximum complete pressure/mass/enthalpy trials. All trials start from
    /// the same original inventory; the caloric increment is applied once.
    pub max_volume_iterations: usize,
    /// Relaxation in (0, 1] for the nonlinear mixing-volume source update.
    /// This changes convergence only, not the accepted EOS or energy budget.
    pub source_relaxation: f64,
}

/// Full original-chart thermal evolution, including sensible/latent crossings.
#[derive(Clone, Copy, Debug)]
pub struct PiecewiseThermalVolumeControls {
    /// Native transport, occupied-volume and nonlinear source-trial budgets.
    pub thermal: ThermalVolumeControls,
    /// Maximum TOTAL original intervals traversed by mixing-volume queries
    /// across all original/incoming donors and all source trials combined.
    /// Each query is preflighted against the remaining budget. Zero permits
    /// only candidates needing no unequal-enthalpy mixing-source evaluation.
    pub max_mixing_intervals: usize,
}

/// Original-chart thermal volume result with its nonlinear closure evidence.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalVolumeStep {
    /// Accepted native flow, material inventories, chart states and work.
    pub volume: CoexistenceVolumeStep,
    /// Number of complete source trials, including the accepted trial.
    pub source_iterations: usize,
    /// Additional source actually applied for conservative enthalpy mixing
    /// [m3/s]. Zero when mixing stays in one affine-specific-volume interval;
    /// crossings of slope-changing knots can have either sign. This is not heat.
    pub mixing_volume_source_m3_s: Vec<f64>,
}

#[derive(Clone, Copy)]
enum ChartVolumeRange {
    Interval {
        lower_knot: usize,
        coexistence_only: bool,
    },
    Full {
        max_mixing_intervals: usize,
    },
}

impl FlowField {
    /// Evolve mobile material across any original sensible or latent intervals.
    /// Uses the original chart and the same conservative pressure/mass/enthalpy
    /// owner as [`Self::advance_thermal_volume`], without fitting a new density
    /// law or splitting/reapplying the prescribed net caloric increment.
    ///
    /// For each original or incoming donor mass m, the source iteration adds
    /// `m * (v(hr) + v'(hr)*(hd-hr) - v(hd))`. The receiver tangent cancels
    /// under the native enthalpy balance. The chart evaluates this remainder
    /// stably on its actual crossed intervals, including signed slope changes
    /// at phase boundaries. Exact accepted donor masses supply the weights.
    /// Both nonlinear trials and total crossed-interval work are bounded; a
    /// failed EOS, conservation, work, domain or cancellation gate returns no
    /// candidate. All trials retain the original mass, heat increment and flow.
    ///
    /// The full original h domain bounds every initial, calorically advanced,
    /// boundary and transported state. Its pressure/composition/source claims
    /// remain fixed. The caller owns thermal/work closure and phase mobility;
    /// this adds no thermal PDE, pressure-dependent EOS or solid stress law.
    pub fn advance_piecewise_thermal_volume(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        curve: &EquilibriumEnthalpyPhaseCurve,
        input: CoexistenceVolumeInput<'_>,
        controls: PiecewiseThermalVolumeControls,
        gate: &CancelGate,
    ) -> Result<ThermalVolumeStep, ChtError> {
        self.advance_chart_volume(
            domain,
            faces,
            curve,
            ChartVolumeRange::Full {
                max_mixing_intervals: controls.max_mixing_intervals,
            },
            input,
            controls.thermal,
            gate,
        )
    }

    /// Expand or contract an explicitly mobile material on one original
    /// caloric interval, including sourced sensible rho(T) and latent volume.
    /// Uses the same inventory carrier as [`Self::advance_coexistence_volume`].
    ///
    /// The chart owns its density interpolation. An affine specific-volume
    /// interval closes in one native transport solve. When rho(h) is affine,
    /// conservative enthalpy mixing has an additional nonlinear volume change.
    /// For receiver density rho_r and donor density rho_d, each original or
    /// incoming mass m contributes the exact concave mixing remainder
    /// `-m/rho_d * ((rho_d-rho_r)/rho_r)^2`. It follows from mass and enthalpy
    /// balance, not from a density fit, volume renormalization or heat source.
    ///
    /// Iteration updates ONLY the continuity source, reusing the original
    /// flow, mass, prescribed net caloric increment and boundary states. Each
    /// trial uses the original paired donor transfers; no trial is a physical
    /// timestep. The actual returned M/rho(h), M and M*h independently pass
    /// their budgets. Only the accepted trial's pressure work is reported.
    /// An admitted initial volume residual never becomes a physical source.
    ///
    /// All material states must remain on the selected adjacent interval.
    /// Nonlinear iteration is bounded and may refuse; it does not infer EOS
    /// pressure feedback to make a sealed incompatible component admissible.
    /// Mixing follows first-order donor transport and requires refinement.
    /// The caller still owns thermal/work closure, momentum prediction and
    /// admitted mobility. This method supplies no solid stress or free surface.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_thermal_volume(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        curve: &EquilibriumEnthalpyPhaseCurve,
        lower_knot: usize,
        input: CoexistenceVolumeInput<'_>,
        controls: ThermalVolumeControls,
        gate: &CancelGate,
    ) -> Result<ThermalVolumeStep, ChtError> {
        self.advance_chart_volume(
            domain,
            faces,
            curve,
            ChartVolumeRange::Interval {
                lower_knot,
                coexistence_only: false,
            },
            input,
            controls,
            gate,
        )
    }

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
        self.advance_chart_volume(
            domain,
            faces,
            curve,
            ChartVolumeRange::Interval {
                lower_knot,
                coexistence_only: true,
            },
            input,
            ThermalVolumeControls {
                volume: controls,
                max_volume_iterations: 1,
                source_relaxation: 1.0,
            },
            gate,
        )
        .map(|step| step.volume)
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_chart_volume(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        curve: &EquilibriumEnthalpyPhaseCurve,
        range: ChartVolumeRange,
        input: CoexistenceVolumeInput<'_>,
        thermal_controls: ThermalVolumeControls,
        gate: &CancelGate,
    ) -> Result<ThermalVolumeStep, ChtError> {
        Self::advance_chart_volume_with_projection(
            domain,
            curve,
            range,
            input,
            thermal_controls,
            gate,
            |density, source| {
                self.project_density_volume_source(
                    domain,
                    faces,
                    density,
                    source,
                    thermal_controls.volume.transport.projection,
                    gate,
                )
                .map(|projection| (projection, ()))
            },
        )
        .map(|(thermal, ())| thermal)
    }

    /// The original chart/transport loop retains only the accepted projection's
    /// private readback. A mechanics consumer can audit that same pressure
    /// impulse after acceptance, without another solve or a projection clone.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn advance_chart_volume_with_projection<T>(
        domain: &VoxelDomain,
        curve: &EquilibriumEnthalpyPhaseCurve,
        range: ChartVolumeRange,
        input: CoexistenceVolumeInput<'_>,
        thermal_controls: ThermalVolumeControls,
        gate: &CancelGate,
        mut project: impl FnMut(&[f64], &[f64]) -> Result<(DensityProjection, T), ChtError>,
    ) -> Result<(ThermalVolumeStep, T), ChtError> {
        poll(gate)?;
        if thermal_controls.max_volume_iterations == 0
            || !thermal_controls.source_relaxation.is_finite()
            || thermal_controls.source_relaxation <= 0.0
            || thermal_controls.source_relaxation > 1.0
        {
            return Err(refuse(
                "invalid thermal-volume iteration budget or relaxation",
            ));
        }
        let controls = thermal_controls.volume;
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
        let (bounds, density_law, mut remaining_mixing_intervals) = match range {
            ChartVolumeRange::Interval {
                lower_knot,
                coexistence_only,
            } => {
                let upper_knot = lower_knot
                    .checked_add(1)
                    .ok_or_else(|| refuse("coexistence interval index overflow"))?;
                let (lower, upper) = curve
                    .knots()
                    .get(lower_knot)
                    .zip(curve.knots().get(upper_knot))
                    .ok_or_else(|| refuse("one adjacent original knot interval required"))?;
                if coexistence_only
                    && (lower.temperature_k != upper.temperature_k
                        || lower.liquid_mass_fraction >= upper.liquid_mass_fraction)
                {
                    return Err(refuse("selected interval must be isothermal coexistence"));
                }
                let law = curve
                    .density_law_on_interval(lower_knot)
                    .map_err(|error| refuse(error.to_string()))?;
                let bounds = [lower.specific_enthalpy_j_kg, upper.specific_enthalpy_j_kg];
                finite_positive("flow.coexistence_enthalpy_span", bounds[1] - bounds[0])?;
                (bounds, Some(law), 0)
            }
            ChartVolumeRange::Full {
                max_mixing_intervals,
            } => (
                [
                    curve.knots()[0].specific_enthalpy_j_kg,
                    curve.knots()[curve.knots().len() - 1].specific_enthalpy_j_kg,
                ],
                None,
                max_mixing_intervals,
            ),
        };
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
        let mut rho_star = vec![0.0; n];
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
            rho_star[c] = next.bulk_density_kg_m3();
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
        let caloric_source = volume_source.clone();
        let mut mixing_source = vec![0.0; n];
        let mut endpoint_density = vec![0.0; n];
        for source_iterations in 1..=thermal_controls.max_volume_iterations {
            poll(gate)?;
            let (projected, readback) = project(&density, &volume_source)?;
            let transport = advect_projected_material(
                domain,
                projected,
                transport_input,
                controls.transport,
                gate,
            )?;
            let mut final_defect = 0.0_f64;
            for c in 0..n {
                poll_every(gate, c)?;
                if !domain.is_fluid(c) {
                    continue;
                }
                let h = transport.specific_enthalpy_j_kg[c];
                if h < bounds[0] || h > bounds[1] {
                    return Err(refuse(
                        "transported enthalpy left the admitted original chart range",
                    ));
                }
                // Trial EOS queries allocate and hash nothing. Construct the
                // identity-bearing phase states only for the accepted trial.
                let value = curve
                    .values_at_specific_enthalpy(h)
                    .map_err(|error| refuse(error.to_string()))?;
                endpoint_density[c] = value.bulk_density_kg_m3;
                let volume = transport.mass_kg[c] / value.bulk_density_kg_m3;
                finite_positive("flow.coexistence_final_volume", volume)?;
                let defect = volume - domain.volume(c);
                finite("flow.coexistence_final_volume_defect", defect)?;
                final_defect = final_defect.max(defect.abs());
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
            if final_defect <= volume_budget {
                let mut phase_states = vec![None; n];
                for c in 0..n {
                    poll_every(gate, c)?;
                    if domain.is_fluid(c) {
                        phase_states[c] = Some(
                            curve
                                .state_at_specific_enthalpy(transport.specific_enthalpy_j_kg[c])
                                .map_err(|error| refuse(error.to_string()))?,
                        );
                    }
                }
                let material_enthalpy_change_j =
                    represented.value + transport.report.enthalpy_change_j;
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
                    required_source_pressure_work_j: transport
                        .projection
                        .report
                        .pressure_source_work_j,
                };
                poll(gate)?;
                return Ok((
                    ThermalVolumeStep {
                        volume: CoexistenceVolumeStep {
                            transport,
                            phase_states,
                            volume_source_m3_s: volume_source,
                            report,
                        },
                        source_iterations,
                        mixing_volume_source_m3_s: mixing_source,
                    },
                    readback,
                ));
            }
            if density_law == Some(PhaseIntervalDensityLaw::AffineSpecificVolume) {
                return Err(refuse(
                    "transported mass and phase density exceed the cell-volume defect budget",
                ));
            }
            if source_iterations == thermal_controls.max_volume_iterations {
                return Err(refuse(
                    "nonlinear thermal-volume iteration budget exhausted",
                ));
            }

            // Original inventory plus incoming exact native donor masses are
            // the convex enthalpy-mixing weights. Outgoing material has the
            // receiver's own endpoint h and therefore zero Taylor remainder.
            let stencil = transport.accepted_mass_transport(gate)?;
            let mut mixed_volume = vec![Sum::default(); n];
            for (row, &c) in stencil.cells.iter().enumerate() {
                poll_every(gate, row)?;
                mixed_volume[c].add(if density_law.is_some() {
                    mixing_volume(
                        stencil.initial_mass_kg[row],
                        rho_star[c],
                        endpoint_density[c],
                    )?
                } else {
                    chart_mixing_volume(
                        curve,
                        stencil.initial_mass_kg[row],
                        transport.specific_enthalpy_j_kg[c],
                        h_star[c],
                        &mut remaining_mixing_intervals,
                        gate,
                    )?
                });
            }
            for (index, route) in stencil.routes.iter().enumerate() {
                poll_every(gate, index)?;
                if let Some(to) = route.to {
                    let c = stencil.cells[to];
                    mixed_volume[c].add(if density_law.is_some() {
                        let donor_density = route.from.map_or_else(
                            || route.inlet.expect("native admitted inlet").density_kg_m3,
                            |from| endpoint_density[stencil.cells[from]],
                        );
                        mixing_volume(route.mass, donor_density, endpoint_density[c])?
                    } else {
                        let donor_h = route.from.map_or_else(
                            || {
                                route
                                    .inlet
                                    .expect("native admitted inlet")
                                    .specific_enthalpy_j_kg
                            },
                            |from| transport.specific_enthalpy_j_kg[stencil.cells[from]],
                        );
                        chart_mixing_volume(
                            curve,
                            route.mass,
                            transport.specific_enthalpy_j_kg[c],
                            donor_h,
                            &mut remaining_mixing_intervals,
                            gate,
                        )?
                    });
                }
            }
            for c in 0..n {
                poll_every(gate, c)?;
                let target = mixed_volume[c].value / dt;
                finite("flow.thermal_mixing_volume_source", target)?;
                mixing_source[c] +=
                    thermal_controls.source_relaxation * (target - mixing_source[c]);
                volume_source[c] = caloric_source[c] + mixing_source[c];
                finite("flow.thermal_total_volume_source", volume_source[c])?;
            }
        }
        unreachable!("positive bounded source loop returns its final trial or refusal")
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
            "phase state must retain the original chart and selected interval",
        ));
    }
    Ok(())
}

/// Negative first-order Taylor remainder of v(h)=1/(a+b*h), evaluated
/// without subtracting nearly equal enthalpies or reciprocal densities.
fn mixing_volume(mass: f64, donor_density: f64, receiver_density: f64) -> Result<f64, ChtError> {
    let relative = (donor_density - receiver_density) / receiver_density;
    let volume = -((mass / donor_density) * relative) * relative;
    finite("flow.thermal_mixing_volume", volume)?;
    if donor_density != receiver_density && volume == 0.0 {
        return Err(refuse("nonlinear mixing-volume increment underflowed"));
    }
    Ok(volume)
}

fn chart_mixing_volume(
    curve: &EquilibriumEnthalpyPhaseCurve,
    mass: f64,
    receiver_h: f64,
    donor_h: f64,
    remaining_intervals: &mut usize,
    gate: &CancelGate,
) -> Result<f64, ChtError> {
    let (specific_volume, traversed) = curve
        .mixing_specific_volume_m3_kg(receiver_h, donor_h, *remaining_intervals, || {
            gate.is_requested()
        })
        .map_err(|error| match error {
            PhaseStateError::MixingVolumeCancelled => ChtError::Cancelled,
            other => refuse(other.to_string()),
        })?;
    *remaining_intervals -= traversed;
    let volume = mass * specific_volume;
    finite("flow.chart_mixing_volume", volume)?;
    if specific_volume != 0.0 && mass != 0.0 && volume == 0.0 {
        return Err(refuse("original-chart mixing-volume increment underflowed"));
    }
    Ok(volume)
}

#[derive(Clone, Copy, Default)]
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
