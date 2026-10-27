//! Source-backed contact-area exchange between two canonical thermal inputs.
//! The original finite-body solver owns the single signed transfer Q. No heat
//! capacity, conduction law, latent-heat rule or nonlinear solver is duplicated.
use super::{InterfaceResistance,SeriesThermalResistance};
use crate::ConductionError;
use crate::lumped::{LumpedContactStep,LumpedContactStepConfig,LumpedEnthalpyBody,
    LumpedThermalTransport,solve_lumped_contact_step};
use fs_exec::Cx;
use fs_material::phase::{UniformEnthalpyStep,UniformEnthalpyStepInput};

/// Exact interface source and finite-body solve retained by both proposals.
#[derive(Clone, Debug, PartialEq)]
pub struct ContactPairReport {
    pub interface: InterfaceResistance,
    /// Physical area held for this thermal step [m2]; zero is separation.
    pub area_m2: f64,
    /// Original source-bearing R''/area series network, absent when separated.
    pub resistance: Option<SeriesThermalResistance>,
    /// Original implicit finite-body result, absent for exactly zero exchange area.
    pub contact: Option<LumpedContactStep>,
}

/// Propose both finite reservoirs through their original material charts.
///
/// `inputs` supply both old states, invariant masses/geometries, internal heat
/// and one common positive duration. Both source terms enter the SAME implicit
/// finite-body solve; no intermediate heat-only state must lie in a chart.
/// The interface is an already resolved, source-backed area-specific resistance.
/// Its exact receipt is retained; R''/area is built by the original interface
/// owner and held for this step. Callers own selection at the appropriate source
/// coordinates; no temperature/pressure dependence or perfect contact is inferred.
///
/// At zero physical area the link is exactly open, not a very large finite R.
/// Each reservoir still receives its internal heat, with source/energy checks.
/// A returning contact resumes from these actual thermal endpoints. One signed
/// Q is -Q in side A's boundary heat and +Q in B's, including latent storage.
/// Publish BOTH proposals together; this routine mutates neither input.
///
/// Reference volume/whole-surface area supplies the original reduced Biot gate
/// with the total contact conductance spread over that surface. This is NOT a
/// qualification of local spreading or flash temperature beneath a small patch.
/// Fixed thermal geometry and negligible body pV work remain assumptions.
/// No ambient/radiation exchange or thermal feedback to mechanics is included.
pub fn advance_contact_pair(
    cx:&Cx<'_>, inputs:[UniformEnthalpyStepInput<'_>;2],
    transports:[&LumpedThermalTransport;2], interface:&InterfaceResistance,
    area_m2:f64, energy_tolerance_j:f64,
) -> Result<[UniformEnthalpyStep<ContactPairReport>;2],ConductionError> {
    poll(cx)?;
    if !area_m2.is_finite() || area_m2<0.0 || !positive(energy_tolerance_j)
        || !positive(inputs[0].duration_s) || inputs[0].duration_s!=inputs[1].duration_s {
        return Err(invalid("contact area, thermal allowance or common clock is invalid"));
    }
    for i in &inputs {
        if ![i.mass_kg,i.volume_m3,i.surface_area_m2].into_iter().all(positive)
            || !i.internal_heat_j.is_finite() || area_m2>i.surface_area_m2
            || i.curve.state_at_specific_enthalpy(i.initial.specific_enthalpy_j_kg())
                .map_err(|e|invalid(&e.to_string()))?!=i.initial {
            return Err(invalid("thermal pair needs valid geometry/heat and each original chart state"));
        }
    }
    // Source identity and continuous transport support remain load-bearing even
    // at separation. Do not skip admission merely because Q or supplied heat is zero.
    let a=body("contact side A",inputs[0],transports[0])?;
    let b=body("contact side B",inputs[1],transports[1])?;
    let (states,q,resistance,contact)=if area_m2==0.0 {
        let state=|i:UniformEnthalpyStepInput<'_>|i.curve.state_at_specific_enthalpy(
            i.initial.specific_enthalpy_j_kg()+i.internal_heat_j/i.mass_kg)
                .map_err(|e|invalid(&e.to_string()));
        ([state(inputs[0])?,state(inputs[1])?],0.0,None,None)
    } else {
        let resistance=SeriesThermalResistance::new(vec![interface.term_for_area("physical contact",area_m2)?])?;
        let contact=solve_lumped_contact_step(cx,&a,&b,&resistance,LumpedContactStepConfig {
            specific_enthalpy_a_j_kg:inputs[0].initial.specific_enthalpy_j_kg(),
            specific_enthalpy_b_j_kg:inputs[1].initial.specific_enthalpy_j_kg(),
            supplied_heat_a_j:inputs[0].internal_heat_j,supplied_heat_b_j:inputs[1].internal_heat_j,
            duration_s:inputs[0].duration_s,energy_tolerance_j,
        })?;
        ([contact.state_a,contact.state_b],contact.heat_a_to_b_j,Some(resistance),Some(contact))
    };
    let boundary=[-q,q];let mut total_residual=0.0;
    for side in 0..2 {
        let residual=inputs[side].mass_kg*(states[side].specific_enthalpy_j_kg()
            -inputs[side].initial.specific_enthalpy_j_kg())-inputs[side].internal_heat_j-boundary[side];
        if !residual.is_finite() || residual.abs()>energy_tolerance_j {
            return Err(invalid("contact pair exceeds a reservoir's energy allowance"));
        }
        total_residual+=residual;
    }
    if !total_residual.is_finite() || total_residual.abs()>energy_tolerance_j {
        return Err(invalid("contact pair exceeds its combined energy allowance"));
    }
    poll(cx)?;
    let report=ContactPairReport {interface:interface.clone(),area_m2,resistance,contact};
    Ok([0,1].map(|side|UniformEnthalpyStep {state:states[side],external_heat_j:boundary[side],
        energy_residual_tolerance_j:energy_tolerance_j,report:report.clone()}))
}
fn body<'a>(name:&str,input:UniformEnthalpyStepInput<'a>,transport:&LumpedThermalTransport)
    ->Result<LumpedEnthalpyBody<'a>,ConductionError> {
    LumpedEnthalpyBody::try_new_with_transport(name,input.mass_kg,input.surface_area_m2,0.0,
        input.volume_m3/input.surface_area_m2,transport.clone(),input.curve)
}
fn positive(v:f64)->bool {v.is_finite() && v>0.0}
fn poll(cx:&Cx<'_>)->Result<(),ConductionError> {
    cx.checkpoint().map_err(|_|ConductionError::Cancelled {stage:"lumped-contact-pair",at:0})
}
fn invalid(what:&str)->ConductionError {
    ConductionError::Interface {interface:"lumped-contact-pair".into(),what:what.into(),
        fix:"supply both original finite thermal states and one admitted source-backed contact area".into()}
}

use super::temperature::ContactIterationConfig;
use crate::lumped::LUMPED_BIOT_CEILING;
use fs_material::phase::EquilibriumPhaseState;
use fs_material::state_point::{MaterialPropertySelection, ResolvedInterfaceStatePoint};
use fs_material::thermal::interface::ThermalInterfacePoint;

/// One ordered material contact law with its original selected source pinned.
///
/// The declared side temperatures and liquid mass fractions follow the actual
/// caloric states. Other supplied source coordinates remain fixed. The binding
/// owns no thermal history and never switches claims at a support boundary.
#[derive(Clone, Debug)]
pub struct MaterialContactLaw<'a> {
    binding: ThermalInterfacePoint<'a>,
    reference: ResolvedInterfaceStatePoint,
}

impl<'a> MaterialContactLaw<'a> {
    /// Select and admit the static interface law at two original material states.
    /// Bulk identities, ordered surfaces and coordinate conventions are checked
    /// by the shared thermal-interface owner.
    pub fn new(
        binding: ThermalInterfacePoint<'a>,
        states: [EquilibriumPhaseState; 2],
        selection: MaterialPropertySelection,
    ) -> Result<Self, ConductionError> {
        let reference = binding
            .resolve_resistance(states, selection)
            .map_err(|error| invalid(&error.to_string()))?;
        Ok(Self { binding, reference })
    }

    /// Original source selection and complete reference conditions.
    #[must_use]
    pub fn reference(&self) -> &ResolvedInterfaceStatePoint {
        &self.reference
    }

    /// Propose both reservoirs using resistance at their actual endpoint states.
    ///
    /// Each bounded Picard trial delegates to the original conservative
    /// finite-body solver with the original states and supplied heats. Only
    /// trial enthalpies are relaxed; resistance values and accepted histories
    /// are never blended. Acceptance checks the physical heat-law residual
    /// against a fresh endpoint source answer and rechecks the reduced Biot
    /// bound there. Both proposals must be published together by their owner.
    ///
    /// Zero area leaves the contact open while preserving internal heating and
    /// source admission. A failed query, cancellation or exhausted iteration
    /// budget returns no pair. Source checks cover the visited points, not
    /// continuous paths through sparse multiaxis support. This fixed-geometry
    /// reduction supplies no uniqueness guarantee, local spreading/flash
    /// temperature model, evolving pressure or ambient/radiative exchange.
    pub fn advance(
        &self,
        cx: &Cx<'_>,
        inputs: [UniformEnthalpyStepInput<'_>; 2],
        transports: [&LumpedThermalTransport; 2],
        area_m2: f64,
        energy_tolerance_j: f64,
        iteration: ContactIterationConfig,
    ) -> Result<[UniformEnthalpyStep<MaterialContactPairReport>; 2], ConductionError> {
        poll(cx)?;
        if iteration.max_iterations == 0
            || !positive(iteration.relaxation)
            || iteration.relaxation > 1.0
            || !positive(iteration.residual_tolerance_j)
        {
            return Err(invalid(
                "material contact needs a positive iteration budget, relaxation in (0,1] and a finite positive physical residual allowance",
            ));
        }
        let selected = self.reference.properties()[0].answer().receipt.selected;
        let mut trial = [inputs[0].initial, inputs[1].initial];
        for index in 0..iteration.max_iterations {
            poll(cx)?;
            let point = self
                .binding
                .at(trial)
                .map_err(|error| invalid(&error.to_string()))?;
            // Retain the complete coordinate schema and exact selected source.
            self.reference
                .at_query_point(self.binding.interface, &point)
                .map_err(|error| invalid(&error.to_string()))?;
            let interface = InterfaceResistance::from_card_pinned(
                "material contact",
                self.binding.interface,
                &point,
                selected,
            )?;
            let pair = advance_contact_pair(
                cx,
                inputs,
                transports,
                &interface,
                area_m2,
                energy_tolerance_j,
            )?;
            let states = [pair[0].state, pair[1].state];
            let endpoint_point = self
                .binding
                .at(states)
                .map_err(|error| invalid(&error.to_string()))?;
            let endpoint = self
                .reference
                .at_query_point(self.binding.interface, &endpoint_point)
                .map_err(|error| invalid(&error.to_string()))?;
            let endpoint_resistance = endpoint.properties()[0].value_si();
            let heat_a_to_b_j = pair[1].external_heat_j;
            let (physical_contact_residual_j, endpoint_maximum_biot) =
                if let Some(contact) = &pair[0].report.contact {
                    let resistance_k_per_w = endpoint_resistance / area_m2;
                    if !positive(resistance_k_per_w) {
                        return Err(invalid("endpoint contact resistance is not representable"));
                    }
                    let residual = heat_a_to_b_j
                        - inputs[0].duration_s
                            * ((states[0].temperature_k() - states[1].temperature_k())
                                / resistance_k_per_w);
                    // The inner solver uses each body's whole-chart minimum
                    // conductivity. Its conservative bound scales with 1/R''.
                    let biot =
                        contact.maximum_biot * (interface.value_m2_k_per_w() / endpoint_resistance);
                    if !residual.is_finite() || !biot.is_finite() || biot > LUMPED_BIOT_CEILING {
                        return Err(invalid(
                            "endpoint heat law or lumped contact Biot bound is inadmissible",
                        ));
                    }
                    (residual, biot)
                } else {
                    (heat_a_to_b_j, 0.0)
                };
            if physical_contact_residual_j.abs() <= iteration.residual_tolerance_j {
                poll(cx)?;
                let report = MaterialContactPairReport {
                    applied: pair[0].report.clone(),
                    endpoint,
                    physical_contact_residual_j,
                    endpoint_maximum_biot,
                    iterations: index + 1,
                };
                return Ok(pair.map(|proposal| UniformEnthalpyStep {
                    state: proposal.state,
                    external_heat_j: proposal.external_heat_j,
                    energy_residual_tolerance_j: proposal.energy_residual_tolerance_j,
                    report: report.clone(),
                }));
            }
            for side in 0..2 {
                let enthalpy = (1.0 - iteration.relaxation) * trial[side].specific_enthalpy_j_kg()
                    + iteration.relaxation * states[side].specific_enthalpy_j_kg();
                trial[side] = inputs[side]
                    .curve
                    .state_at_specific_enthalpy(enthalpy)
                    .map_err(|error| invalid(&error.to_string()))?;
            }
        }
        Err(invalid(
            "material contact exhausted its bounded iterations before closing the physical heat law",
        ))
    }
}

/// Applied conservative inner solve and the actual endpoint material answer.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialContactPairReport {
    /// Original source-bearing fixed-resistance inner solve, without alteration.
    pub applied: ContactPairReport,
    /// The same selected interface claim at both actual endpoint states.
    pub endpoint: ResolvedInterfaceStatePoint,
    /// Q - dt (T_A - T_B) / (R''_endpoint / area), in joules.
    pub physical_contact_residual_j: f64,
    /// Conservative contact Biot bound using the endpoint resistance.
    pub endpoint_maximum_biot: f64,
    /// Number of inner finite-body solves used by the accepted proposal.
    pub iterations: usize,
}
