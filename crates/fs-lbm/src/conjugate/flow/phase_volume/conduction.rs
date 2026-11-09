//! Original source conductivity drives the existing caloric/volume transport.

use fs_blake3::ContentHash;
use fs_material::thermal::conduction::phase::{PhaseConductivitySample, PreparedPhaseConductivity};

use super::{
    CancelGate, ChtError, CoexistenceVolumeInput, FlowFace, FlowField,
    PiecewiseThermalVolumeControls, Sum, ThermalVolumeStep, VoxelDomain, finite, finite_positive,
    poll, poll_every,
};

mod boundary;
#[cfg(feature = "implicit-enthalpy-coupling")]
mod implicit;
mod momentum;
pub use boundary::ThermalVolumeBathFace;
#[cfg(feature = "implicit-enthalpy-coupling")]
pub use implicit::{
    ImplicitConductiveResistiveThermalVolumeMomentumStep,
    ImplicitConductiveThermalVolumeMomentumStep, ImplicitConductiveThermalVolumeReport,
    ImplicitConductiveThermalVolumeStep, ImplicitConductiveViscousThermalVolumeMomentumStep,
    ImplicitConjugateResistiveThermalVolumeMomentumStep,
    ImplicitConjugateThermalVolumeMomentumStep, ImplicitConjugateThermalVolumeReport,
    ImplicitConjugateThermalVolumeStep, ImplicitConjugateViscousResistiveThermalVolumeMomentumStep,
    ImplicitConjugateViscousThermalVolumeMomentumStep, ImplicitResistanceHeatWork,
    ImplicitThermalVolumeLinearControls, ImplicitThermalVolumeSolveControls,
    ImplicitThermalVolumeWorkspace, ImplicitThermorheologicalThermalVolumeMomentumStep,
    ImplicitThermorheologicalResistiveThermalVolumeMomentumStep,
    ImplicitThermoviscousResistiveThermalVolumeMomentumStep,
    ImplicitThermoviscousThermalVolumeMomentumStep, ThermalViscosityFeedbackControls,
    ThermalViscosityFeedbackReport, ThermalVolumeBathExchange, ThermalVolumeBoundaries,
    ThermalVolumeBoundaryFace, ThermalVolumeContact, ThermalVolumeSolidCell,
    ThermalVolumeSolidEndpoint,
};
pub use momentum::{
    ConductiveResistiveThermalVolumeMomentumControls, ConductiveResistiveThermalVolumeMomentumStep,
    ConductiveThermalVolumeMomentumControls, ConductiveThermalVolumeMomentumStep,
    ConductiveViscousThermalVolumeMomentumStep,
};

/// The existing full-chart budgets and an explicit bound on conductive work.
#[derive(Clone, Copy, Debug)]
pub struct ConductiveThermalVolumeControls {
    /// Original pressure, mass, enthalpy, occupied-volume and chart-work gates.
    /// Its enthalpy-defect budget covers heat accumulation, adding authored
    /// increments, caloric representation and native transport TOGETHER.
    pub volume: PiecewiseThermalVolumeControls,
    /// Maximum internal fluid/fluid plus explicitly supplied bath faces inspected,
    /// including insulating and isothermal faces. Each face is visited once,
    /// independent of source trials.
    pub max_conductive_faces: usize,
}

/// Actual original-source diffusion and its separate numerical heat defects.
#[derive(Clone, Debug, PartialEq)]
pub struct ConductiveThermalVolumeReport {
    /// Identity already retained by the prepared original k(T,f) source.
    pub conductivity_source_identity: ContentHash,
    pub conductive_faces: usize,
    /// Number of explicitly supplied fluid-boundary bath faces, including zero heat.
    pub boundary_faces: usize,
    /// Signed bath heat by face orientation, ordered x-, x+, y-, y+, z-, z+ [J].
    /// This is external heat, separate from native advective enthalpy transport.
    pub boundary_heat_j: [f64; 6],
    /// Sum of absolute individual bath-face transfers [J]. Opposed heating and
    /// cooling remain visible even when their signed supply cancels.
    pub absolute_boundary_heat_j: f64,
    /// Sum of |dt G (T_b-T_a)|, counting each internal face once [J].
    pub exchanged_heat_j: f64,
    /// Largest sum of internal enthalpy-secant weights and bath conductances
    /// weighted by the original chart's maximum dT/dh. At most one; this is a
    /// frozen-coefficient explicit-step bound, not a fitted Cp model.
    pub maximum_enthalpy_mixing_weight: f64,
    /// Absolute net cell heat minus independently summed bath supply [J].
    pub conduction_balance_defect_j: f64,
    /// Sum of absolute floating-point addition residuals when accumulating
    /// face heat and combining it with the caller's authored net increments.
    /// Includes heat lost beside a much larger authored increment [J].
    pub increment_assembly_defect_j: f64,
    /// Assembly plus conductive balance plus original caloric/transport
    /// defects [J], including composed molecular-heat assembly or coupled
    /// resistance heat/constitutive defects,
    /// all checked against the SAME native enthalpy budget.
    pub absolute_enthalpy_defect_j: f64,
}

/// An atomic original-source conductive and conservative volume-flow endpoint.
#[derive(Clone, Debug, PartialEq)]
pub struct ConductiveThermalVolumeStep {
    /// Unmodified original chart, mass/enthalpy, pressure and volume result.
    pub thermal: ThermalVolumeStep,
    /// Internal and bath heat received by each cell [J], before authored increments.
    /// Solid entries are zero. Numerical pressure loss is never placed here.
    pub conductive_enthalpy_increment_j: Vec<f64>,
    pub report: ConductiveThermalVolumeReport,
}

// All public consumers retain the same original heat preparation and final
// source/enthalpy gates. The evolution value owns its one accepted transport.
struct ConductiveAdvance<T> {
    evolution: T,
    increments: Vec<f64>,
    report: ConductiveThermalVolumeReport,
}

impl FlowField {
    /// Conduct original material heat and evolve its resulting density/volume.
    /// Calling this method explicitly declares the supplied scalar source
    /// isotropic. Its SAME borrowed caloric chart owns T(h), phase and rho(h).
    /// No density, heat capacity, phase mixing rule or new source is fitted.
    ///
    /// At the ORIGINAL cell states, each internal fluid/fluid face supplies
    /// `G = area / (half_width_a/k_a + half_width_b/k_b)` and paired heat
    /// `Q = dt*G*(T_b-T_a)`. One source sample per fluid cell serves every face;
    /// zero k insulates. Solid interfaces and exterior thermal faces are
    /// insulated. Graded Cartesian half widths and face areas remain physical.
    ///
    /// Before native transport, each cell must satisfy
    /// `sum_j Q_ij / (M_i*(h_j-h_i)) <= 1`, ignoring zero-temperature-difference
    /// faces. The nonnegative secants make diffusion a convex combination in
    /// h, including true latent plateaus with dT/dh=0. An excessive requested
    /// step refuses: there is no enthalpy clipping or hidden thermal substep.
    /// This explicit diffusion and the subsequent native implicit donor
    /// transport are first-order operator splitting and require refinement.
    ///
    /// Conductive heat is added to the caller's NET caloric increments for the
    /// same interval. Conduction uses the original state, not an intermediate
    /// state after those sources. The original source must cover initial,
    /// incoming, calorically advanced and accepted enthalpies. Addition
    /// residuals, including a small Q lost beside a large authored increment,
    /// share the native enthalpy-defect budget with the complete volume step.
    ///
    /// The inherited caloric/reference-pressure convention still applies.
    /// The caller owns phase mobility, the momentum predictor and thermal/work
    /// closure. Native required source-pressure work remains separate; this
    /// adds no pressure-dependent EOS, closed thermomechanical feedback,
    /// fresh momentum integration, solid mechanics or free-surface geometry.
    /// Any domain, source, budget, convergence or cancellation refusal returns
    /// no endpoint and changes neither the original field nor its inventories.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_isotropic_conductive_thermal_volume(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        conductivity: &PreparedPhaseConductivity<'_>,
        input: CoexistenceVolumeInput<'_>,
        controls: ConductiveThermalVolumeControls,
        gate: &CancelGate,
    ) -> Result<ConductiveThermalVolumeStep, ChtError> {
        self.advance_isotropic_conductive_thermal_volume_with_baths(
            domain,
            faces,
            conductivity,
            &[],
            input,
            controls,
            gate,
        )
    }

    /// Advance the same source-owned conduction and volume transport with
    /// explicit thermal baths on selected fluid-domain boundary faces.
    ///
    /// Each bath adds `dt*A*(T_bath-T_cell)/(d/k_cell+R)`, where the original
    /// graded grid supplies the physical area A and cell-centre-to-face distance
    /// d (half the local width). R is the authored area-specific external
    /// resistance; zero gives a prescribed surface temperature. Conductivity
    /// is frozen at the SAME original cell state as internal conduction.
    /// Exterior faces and faces adjoining fixed solids admit a bath; all
    /// unlisted faces stay insulated. A fixed-solid neighbor supplies geometry
    /// only. A bath supplies no material mass or extra caloric storage.
    ///
    /// Bath conductances share the internal explicit-step bound, using the
    /// maximum dT/dh of the original chart over its admitted conductivity range.
    /// True latent plateaus retain their original energy and zero thermal slope.
    /// Bath temperature is an external condition, not a new material query;
    /// actual heated and transported cells must still stay in source support.
    /// This is frozen-coefficient, first-order heat/volume splitting. It does
    /// not solve a finite wall, a nonlinear through-wall profile or a bath EOS.
    ///
    /// Signed bath heat enters native enthalpy once and is reported separately
    /// from advective inlet/outlet enthalpy and required pressure-source work.
    /// Input/source, topology, work, energy, stability and cancellation refusals
    /// publish no endpoint and leave every supplied inventory unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_isotropic_conductive_thermal_volume_with_baths(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: &[ThermalVolumeBathFace],
        input: CoexistenceVolumeInput<'_>,
        controls: ConductiveThermalVolumeControls,
        gate: &CancelGate,
    ) -> Result<ConductiveThermalVolumeStep, ChtError> {
        let accepted = self.conduct_then_advance(
            domain,
            conductivity,
            baths,
            input,
            controls,
            gate,
            |heated| {
                self.advance_piecewise_thermal_volume(
                    domain,
                    faces,
                    conductivity.caloric_curve(),
                    heated,
                    controls.volume,
                    gate,
                )
            },
            |thermal| (thermal, 0.0),
        )?;
        Ok(ConductiveThermalVolumeStep {
            thermal: accepted.evolution,
            conductive_enthalpy_increment_j: accepted.increments,
            report: accepted.report,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn conduct_then_advance<T>(
        &self,
        domain: &VoxelDomain,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: &[ThermalVolumeBathFace],
        input: CoexistenceVolumeInput<'_>,
        controls: ConductiveThermalVolumeControls,
        gate: &CancelGate,
        advance: impl FnOnce(CoexistenceVolumeInput<'_>) -> Result<T, ChtError>,
        thermal: impl Fn(&T) -> (&ThermalVolumeStep, f64),
    ) -> Result<ConductiveAdvance<T>, ChtError> {
        let n = domain.cell_count();
        let transport = controls.volume.thermal.volume.transport;
        let dt = transport.projection.duration_s;
        let energy_budget = transport.maximum_enthalpy_defect_j;
        let bounds = conductivity.enthalpy_bounds_j_kg();
        let PreparedInput { samples, baths } =
            prepare_input(domain, conductivity, baths, input, controls, gate)?;

        let mut increments = vec![0.0; n];
        let mut mixing_weights = vec![Sum::default(); n];
        let mut assembly_defect = Sum::default();
        let mut exchanged = Sum::default();
        let mut conductive_faces = 0usize;
        for c in 0..n {
            poll_every(gate, c)?;
            let Some(a) = samples[c] else { continue };
            for face in [1usize, 3, 5] {
                let Some(neighbor) = domain.neighbor(c, face) else {
                    continue;
                };
                let Some(b) = samples[neighbor] else { continue };
                conductive_faces = conductive_faces
                    .checked_add(1)
                    .filter(|&count| count <= controls.max_conductive_faces)
                    .ok_or_else(|| refuse("conductive face-work budget exhausted"))?;
                let ka = a.conductivity_w_m_k();
                let kb = b.conductivity_w_m_k();
                let delta_t = b.temperature_k() - a.temperature_k();
                if ka == 0.0 || kb == 0.0 || delta_t == 0.0 {
                    continue;
                }
                let axis = face / 2;
                let ra = (0.5 * domain.widths(c)[axis]) / ka;
                let rb = (0.5 * domain.widths(neighbor)[axis]) / kb;
                finite_positive("flow.conductive_half_resistance", ra)?;
                finite_positive("flow.conductive_half_resistance", rb)?;
                // Scale the harmonic denominator so adding two individually
                // representable resistances cannot overflow needlessly.
                let largest = ra.max(rb);
                let g = (domain.face_area(c, axis) / largest) / (1.0 + ra.min(rb) / largest);
                finite_positive("flow.conductive_face_conductance", g)?;
                let q = dt * (g * delta_t);
                finite("flow.conductive_face_heat", q)?;
                if q == 0.0 {
                    return Err(refuse("nonzero conductive face heat underflowed"));
                }
                let delta_h = b.specific_enthalpy_j_kg() - a.specific_enthalpy_j_kg();
                finite("flow.conductive_face_enthalpy_difference", delta_h)?;
                if delta_h == 0.0 || delta_h.is_sign_positive() != delta_t.is_sign_positive() {
                    return Err(refuse(
                        "conductive source is not monotone in original enthalpy",
                    ));
                }
                let weight_mass = q / delta_h;
                finite_positive("flow.conductive_secant_mass", weight_mass)?;
                for cell in [c, neighbor] {
                    let weight = weight_mass / input.mass_kg[cell];
                    finite_positive("flow.conductive_enthalpy_weight", weight)?;
                    mixing_weights[cell].add(weight);
                    if !mixing_weights[cell].value.is_finite() || mixing_weights[cell].value > 1.0 {
                        return Err(refuse(
                            "explicit conductive enthalpy mixing weight exceeds one",
                        ));
                    }
                }
                for (cell, heat) in [(c, q), (neighbor, -q)] {
                    let (sum, residual) = accounted_add(increments[cell], heat)?;
                    increments[cell] = sum;
                    assembly_defect.add(residual.abs());
                }
                exchanged.add(q.abs());
                finite("flow.conductive_exchanged_heat", exchanged.value)?;
            }
        }

        let boundary = boundary::apply(
            domain,
            conductivity,
            &baths,
            &samples,
            input.mass_kg,
            dt,
            conductive_faces,
            controls.max_conductive_faces,
            &mut increments,
            &mut mixing_weights,
            &mut assembly_defect,
            gate,
        )?;

        let mut combined = Vec::with_capacity(n);
        let mut conductive_balance = Sum::default();
        let mut maximum_weight = 0.0_f64;
        for c in 0..n {
            poll_every(gate, c)?;
            conductive_balance.add(increments[c]);
            maximum_weight = maximum_weight.max(mixing_weights[c].value);
            let (sum, residual) = accounted_add(input.net_enthalpy_increment_j[c], increments[c])?;
            assembly_defect.add(residual.abs());
            combined.push(sum);
            if let Some(sample) = samples[c] {
                let h = sample.specific_enthalpy_j_kg() + sum / input.mass_kg[c];
                if !(bounds[0]..=bounds[1]).contains(&h) {
                    return Err(refuse(
                        "calorically advanced state leaves the original conductivity source",
                    ));
                }
            }
        }
        for heat in boundary.heat_j {
            if heat != 0.0 {
                conductive_balance.add(-heat);
            }
        }
        let balance_defect = conductive_balance.value.abs();
        let extra_defect = assembly_defect.value + balance_defect;
        finite("flow.conductive_heat_assembly_defect", extra_defect)?;
        if extra_defect > energy_budget {
            return Err(refuse(
                "conductive heat assembly exceeds the shared enthalpy-defect budget",
            ));
        }
        poll(gate)?;
        let evolution = advance(CoexistenceVolumeInput {
            net_enthalpy_increment_j: &combined,
            ..input
        })?;
        let (thermal, evolution_heat_defect) = thermal(&evolution);
        if !evolution_heat_defect.is_finite() || evolution_heat_defect < 0.0 {
            return Err(refuse("invalid composed thermal heat-assembly defect"));
        }
        let absolute_enthalpy_defect_j =
            extra_defect + evolution_heat_defect + thermal.volume.report.absolute_enthalpy_defect_j;
        finite(
            "flow.conductive_total_enthalpy_defect",
            absolute_enthalpy_defect_j,
        )?;
        if absolute_enthalpy_defect_j > energy_budget {
            return Err(refuse(
                "conduction and native volume transport exceed the shared enthalpy-defect budget",
            ));
        }
        for (c, state) in thermal.volume.phase_states.iter().enumerate() {
            poll_every(gate, c)?;
            if let Some(state) = state {
                admit_state(conductivity, *state)?;
            }
        }
        poll(gate)?;
        Ok(ConductiveAdvance {
            evolution,
            increments,
            report: ConductiveThermalVolumeReport {
                conductivity_source_identity: conductivity.identity(),
                conductive_faces,
                boundary_faces: baths.len(),
                boundary_heat_j: boundary.heat_j,
                absolute_boundary_heat_j: boundary.absolute_heat_j,
                exchanged_heat_j: exchanged.value,
                maximum_enthalpy_mixing_weight: maximum_weight,
                conduction_balance_defect_j: balance_defect,
                increment_assembly_defect_j: assembly_defect.value,
                absolute_enthalpy_defect_j,
            },
        })
    }
}

struct PreparedInput<'a> {
    samples: Vec<Option<PhaseConductivitySample>>,
    baths: std::collections::BTreeMap<(usize, usize), &'a ThermalVolumeBathFace>,
}

fn admit_state(
    conductivity: &PreparedPhaseConductivity<'_>,
    state: fs_material::phase::EquilibriumPhaseState,
) -> Result<(), ChtError> {
    let curve = conductivity.caloric_curve();
    let bounds = conductivity.enthalpy_bounds_j_kg();
    if state.phase_curve_identity() != curve.identity()
        || state.material_card_identity() != curve.material_card_identity()
        || !(bounds[0]..=bounds[1]).contains(&state.specific_enthalpy_j_kg())
    {
        Err(refuse(
            "phase state is outside the original conductivity source/chart",
        ))
    } else {
        Ok(())
    }
}

fn prepare_input<'a>(
    domain: &VoxelDomain,
    conductivity: &PreparedPhaseConductivity<'_>,
    baths: &'a [ThermalVolumeBathFace],
    input: CoexistenceVolumeInput<'_>,
    controls: ConductiveThermalVolumeControls,
    gate: &CancelGate,
) -> Result<PreparedInput<'a>, ChtError> {
    poll(gate)?;
    let n = domain.cell_count();
    let transport = controls.volume.thermal.volume.transport;
    if n > transport.max_cells
        || input.mass_kg.len() != n
        || input.phase_states.len() != n
        || input.net_enthalpy_increment_j.len() != n
    {
        return Err(refuse(
            "conductive inventory shape or cell budget is invalid",
        ));
    }
    finite_positive("flow.conductive_duration", transport.projection.duration_s)?;
    if !transport.maximum_enthalpy_defect_j.is_finite() || transport.maximum_enthalpy_defect_j < 0.0
    {
        return Err(refuse("invalid combined conductive enthalpy-defect budget"));
    }
    let baths = boundary::prepare(domain, baths, controls.max_conductive_faces, gate)?;
    for state in input.boundary_inflow.into_iter().flatten() {
        admit_state(conductivity, state)?;
    }
    let mut samples = Vec::with_capacity(n);
    for c in 0..n {
        poll_every(gate, c)?;
        finite(
            "flow.authored_net_enthalpy_increment",
            input.net_enthalpy_increment_j[c],
        )?;
        if !domain.is_fluid(c) {
            if input.mass_kg[c] != 0.0
                || input.phase_states[c].is_some()
                || input.net_enthalpy_increment_j[c] != 0.0
            {
                return Err(refuse(
                    "solid grid cells cannot contain conductive phase inventory",
                ));
            }
            samples.push(None);
            continue;
        }
        finite_positive("flow.conductive_cell_mass", input.mass_kg[c])?;
        let state = input.phase_states[c]
            .ok_or_else(|| refuse("each fluid cell needs an original conductive phase state"))?;
        admit_state(conductivity, state)?;
        samples.push(Some(
            conductivity
                .state_at_specific_enthalpy(state.specific_enthalpy_j_kg())
                .map_err(|error| refuse(error.to_string()))?,
        ));
    }
    Ok(PreparedInput { samples, baths })
}

/// Error-free addition residual for finite nonoverflowing binary64 operands.
/// The heat owner gates its absolute value rather than hiding lost small heat
/// inside a large caller-authored net increment or a cancelling neighbor sum.
fn accounted_add(a: f64, b: f64) -> Result<(f64, f64), ChtError> {
    let sum = a + b;
    finite("flow.conductive_increment_sum", sum)?;
    let b_virtual = sum - a;
    let residual = (a - (sum - b_virtual)) + (b - b_virtual);
    finite("flow.conductive_increment_addition_residual", residual)?;
    Ok((sum, residual))
}

fn refuse(reason: impl Into<String>) -> ChtError {
    ChtError::InvalidInput {
        field: "flow.conductive_thermal_volume",
        reason: reason.into(),
    }
}

#[cfg(test)]
mod tests;
