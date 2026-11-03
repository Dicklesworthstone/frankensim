//! Native steady projects consume the existing source-backed contact owner.

use fs_conduction::interface::temperature::steady::SteadyContactConfig;
use fs_conduction::interface::temperature::{TemperatureContact, TemperatureContacts};
use fs_conduction::{ConductionProblem, ConductionSolver, DofMap, SolveConfig};
use fs_exec::Cx;
use fs_matdb::{ClaimId, PropertyValue, QueryPoint};
use fs_project::{BindingRequirements, BindingTarget, CardLibrary, ProjectSpec};

use super::{SolveRefusal, conduction_error, json_string, parse_claim_id};
use crate::cards::CardPackSet;

const MAX_RESPONSES: usize = 64;
pub(super) const NO_CLAIM: &str = "Source-backed matching P1 contact uses the arithmetic mean of the two independent face temperature means as the declared T coordinate. Each local response retains its original card, selected claim and manufactured joint state. The existing steady spatial solver retains heterogeneous k(T), prescribed boundaries and natural/transported-air consumers. The returned field passes freshly queried physical watt residual and heat-balance gates. This is a bounded contact fixed point, not an arbitrary-law convergence or uniqueness proof. Radiation, transient storage, adaptive-goal refinement, adjoints and inverse/roundoff certificates for this nonlinear contact path are not supplied. Uniform spatial ladder differences retain their existing Estimated meaning. No evolving pressure, geometry or new material authority is inferred.";

struct Source {
    name: String,
    card: String,
    point: QueryPoint,
    selected: ClaimId,
}

pub(super) struct Prepared {
    library: CardLibrary,
    sources: Vec<Source>,
    temperature_axis: String,
}

/// Counted contact work and the last physically accepted field on one mesh.
#[derive(Default)]
pub(super) struct Work {
    responses: usize,
    nonlinear: usize,
    krylov: usize,
    last_temperature: Option<Vec<f64>>,
}

impl Work {
    fn record(
        &mut self,
        responses: usize,
        nonlinear: usize,
        krylov: usize,
    ) -> Result<(), SolveRefusal> {
        let add = |old: usize, new: usize| {
            old.checked_add(new)
                .ok_or_else(|| bad("counted contact work overflowed"))
        };
        self.responses = add(self.responses, responses)?;
        self.nonlinear = add(self.nonlinear, nonlinear)?;
        self.krylov = add(self.krylov, krylov)?;
        Ok(())
    }
}

fn bad(message: impl Into<String>) -> SolveRefusal {
    conduction_error(
        "cli-solve-conduction-contact",
        message,
        "inspect the contact source support, named matching surfaces and steady solver controls",
    )
}

fn lower(error: fs_conduction::ConductionError) -> SolveRefusal {
    match error {
        fs_conduction::ConductionError::Cancelled { .. } => conduction_error(
            "cli-solve-cancelled",
            "source-backed thermal contact was cancelled",
            "resume the accepted pipeline prefix",
        ),
        other => bad(format!("source-backed thermal contact refused: {other}")),
    }
}

fn number(value: f64) -> Result<String, SolveRefusal> {
    super::canonical_f64(value).ok_or_else(|| bad("nonfinite contact evidence"))
}

fn poll(cx: &Cx<'_>) -> Result<(), SolveRefusal> {
    cx.checkpoint().map_err(|_| {
        conduction_error(
            "cli-solve-cancelled",
            "source-backed contact work was cancelled",
            "resume the accepted pipeline prefix",
        )
    })
}

/// Resolve once before numerical work. A curve's equal endpoints do not erase
/// interior variation, and a nonconstant two-coordinate law naming T is never frozen.
pub(super) fn prepare(
    spec: &ProjectSpec,
    cards: &CardPackSet,
) -> Result<Option<Prepared>, SolveRefusal> {
    if spec.interface_cards.as_deref().unwrap_or(&[]).is_empty() {
        return Ok(None);
    }
    let library = cards.library();
    let requirements = BindingRequirements::thermal_steady_v1();
    let resolution = fs_project::resolve_bindings(spec, &library, &requirements);
    if let Some(error) = resolution.violations.first() {
        return Err(conduction_error(
            error.code,
            error.what.clone(),
            error.fix.clone(),
        ));
    }
    let mut sources = Vec::new();
    let mut variable = false;
    for binding in resolution.bindings {
        let BindingTarget::Interface(name) = &binding.target else {
            continue;
        };
        let property = binding
            .properties
            .iter()
            .find(|row| row.property == fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY)
            .ok_or_else(|| bad("resolved contact has no resistance property"))?;
        let card = library
            .interface(&binding.card)
            .ok_or_else(|| bad("resolved contact lost its original source card"))?;
        let selected = parse_claim_id(&property.selected_claim, name)?;
        let claim = card
            .claims()
            .claim(selected)
            .ok_or_else(|| bad("resolved contact lost its selected source claim"))?;
        variable |= property.value_lo.to_bits() != property.value_hi.to_bits()
            || match &claim.value {
                PropertyValue::Curve {
                    abscissa, knots, ..
                } if abscissa == &requirements.temperature_axis => {
                    knots.iter().any(|(t, value)| {
                        (binding.range_lo..=binding.range_hi).contains(t)
                            && value.to_bits() != property.value_lo.to_bits()
                    })
                }
                PropertyValue::Surface2D { axes, samples, .. } => {
                    // A globally flat source table has no temperature feedback.
                    // Every other T surface stays conservative: equal queries
                    // alone cannot prove a constant path through its cells.
                    axes.contains(&requirements.temperature_axis)
                        && samples.first().is_none_or(|(_, reference)| {
                            samples
                                .iter()
                                .any(|(_, value)| value.to_bits() != reference.to_bits())
                        })
                }
                _ => false,
            };
        let declared = spec
            .interface_cards
            .as_deref()
            .unwrap_or(&[])
            .iter()
            .find(|row| &row.interface == name)
            .ok_or_else(|| bad("resolved contact lost its manufactured joint state"))?;
        let point = fs_project::interface_state::query_point(
            &declared.state,
            &requirements.temperature_axis,
            binding.range_lo,
            fs_project::interface_state::StatePoint::Nominal,
        )
        .map_err(|error| bad(format!("contact query point refused: {}", error.detail)))?;
        sources.push(Source {
            name: name.clone(),
            card: binding.card,
            point,
            selected,
        });
    }
    if variable {
        Ok(Some(Prepared {
            library,
            sources,
            temperature_axis: requirements.temperature_axis,
        }))
    } else {
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cards::{CardPackKind, RawCardPack};
    use fs_matdb::{ClaimSet, InterpolationPolicy, NormalizedInterfacePack, NormalizedPack};

    #[test]
    fn constant_surface_contact_keeps_legacy_path_after_binding() {
        let original = NormalizedInterfacePack::from_bytes(include_bytes!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../examples/contact-pair/cold-hot.fsintpk"
        )))
        .unwrap();
        for off_point_resistance in [0.1, 0.2] {
            let mut claims = ClaimSet::new();
            for id in original.card().claims().observation_ids() {
                claims
                    .register_observation(original.card().claims().observation(id).unwrap().clone())
                    .unwrap();
            }
            let mut claim = original
                .card()
                .claims()
                .claims_ordered()
                .next()
                .unwrap()
                .1
                .clone();
            claim.value = PropertyValue::Surface2D {
                axes: ["T".into(), "normal_pressure".into()],
                axis_dims: [
                    fs_project::spec::dims::TEMPERATURE,
                    fs_project::spec::dims::PRESSURE,
                ],
                samples: vec![
                    ([280.0, 0.9e6], off_point_resistance),
                    ([320.0, 0.9e6], 0.1),
                    ([300.0, 1.0e6], 0.1),
                ],
                triangles: vec![[0, 1, 2]],
                dims: fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
            claims.insert_claim(claim).unwrap();
            let pack = NormalizedPack::new(
                "synthetic-contact-surface",
                "unit-test",
                original.source_artifact(),
                original.redistribution_terms(),
                claims,
                Vec::new(),
                Vec::new(),
            )
            .unwrap();
            let pack = NormalizedInterfacePack::new(
                original.card().surface_a().clone(),
                original.card().surface_b().clone(),
                original.card().context().clone(),
                pack,
            )
            .unwrap();
            let cards = CardPackSet::admit(vec![
                RawCardPack {
                    kind: CardPackKind::Interface,
                    source: "surface-test".into(),
                    bytes: pack.to_bytes(),
                    expect: None,
                },
                RawCardPack {
                    kind: CardPackKind::Material,
                    source: "solid-test".into(),
                    bytes: include_bytes!(concat!(
                        env!("CARGO_MANIFEST_DIR"),
                        "/../../data/reference-project/aa6061.fsmcdpk"
                    ))
                    .to_vec(),
                    expect: None,
                },
            ])
            .unwrap();
            let mut spec = fs_project::parse_sexpr_migrating(include_str!(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../../examples/contact-pair/contact-pair.fsim"
            )))
            .unwrap()
            .decoded
            .spec;
            // Native box-support admission currently accepts a Surface2D at
            // a point; source-cell corners alone do not prove an interval.
            for material in spec.materials.as_mut().unwrap() {
                material.temp_lo.value = 300.0;
                material.temp_hi.value = 300.0;
            }
            let joint = &mut spec.interface_cards.as_mut().unwrap()[0];
            joint.card = cards.interfaces()[0].card().to_hex();
            if let fs_project::InterfaceState::DryContact {
                pressure_half_width,
                ..
            } = &mut joint.state
            {
                pressure_half_width.value = 0.0;
            }
            let resolution = fs_project::resolve_bindings(
                &spec,
                &cards.library(),
                &BindingRequirements::thermal_steady_v1(),
            );
            assert!(resolution.admissible(), "{:?}", resolution.violations);
            let binding = resolution
                .bindings
                .iter()
                .find(|row| matches!(row.target, BindingTarget::Interface(_)))
                .unwrap();
            assert_eq!(binding.properties[0].value_lo.to_bits(), 0.1_f64.to_bits());
            assert_eq!(
                binding.properties[0].value_lo,
                binding.properties[0].value_hi
            );
            assert_eq!(
                prepare(&spec, &cards).unwrap().is_some(),
                off_point_resistance != 0.1,
                "globally flat tables keep constant-contact compositions; other T tables stay conservative"
            );
        }
    }
}

pub(super) fn admit(
    spec: &ProjectSpec,
    adjoint: bool,
    mean_bound: bool,
) -> Result<(), SolveRefusal> {
    let setup = spec
        .cooling
        .as_ref()
        .and_then(|row| row.conduction.as_ref())
        .ok_or_else(|| bad("temperature-varying contacts need an explicit conduction setup"))?;
    if super::transient::requested(spec)
        || setup.radiation.is_some()
        || setup.enclosure_radiation.is_some()
        || adjoint
        || mean_bound
        || spec
            .solver
            .as_ref()
            .is_some_and(|solver| solver.fidelity == super::SOLVER_FIDELITY_ADAPTIVE)
    {
        return Err(conduction_error(
            "cli-solve-conduction-contact-composition",
            "temperature-varying native contact currently admits steady base/uniform-ladder solves with prescribed, natural or airflow convection; this requested composition needs the complete contact feedback",
            "use the steady base or ladder path without radiation or derivative/bound outputs; never replace the contact law by a frozen resistance",
        ));
    }
    Ok(())
}

impl Prepared {
    pub(super) fn solve(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        resolution: &fs_project::ConductionInterfaceResolution,
        mut config: SolveConfig,
        max_pairs: usize,
        work: &mut Work,
    ) -> Result<(super::radiation::SolidSolution, String), SolveRefusal> {
        poll(cx)?;
        let mut laws = Vec::with_capacity(self.sources.len());
        for source in &self.sources {
            poll(cx)?;
            let pairs = resolution
                .pairs
                .iter()
                .filter(|pair| pair.interface == source.name)
                .map(fs_project::ResolvedConductionInterfacePair::face_pair)
                .collect::<Vec<_>>();
            if pairs.is_empty() {
                return Err(bad(format!(
                    "contact {} has no matching faces",
                    source.name
                )));
            }
            let card = self
                .library
                .interface(&source.card)
                .ok_or_else(|| bad("contact lost its admitted source card"))?;
            laws.push(
                TemperatureContact::new_pinned(
                    &source.name,
                    pairs,
                    card,
                    &source.point,
                    &self.temperature_axis,
                    source.selected,
                )
                .map_err(lower)?,
            );
        }
        let contacts = TemperatureContacts::new(
            cx,
            problem.mesh,
            problem.boundary,
            Vec::new(),
            laws,
            max_pairs,
        )
        .map_err(lower)?;
        // Use the spatial owner's actual initial-field and material-span policy.
        // This prepares state only; it executes no uncounted thermal response.
        let initial = ConductionSolver::new(problem, config.clone()).map_err(lower)?;
        let dofs = DofMap::new(problem.boundary, problem.mesh.vertex_count()).map_err(lower)?;
        let initial = dofs.scatter(&initial.state().free_temperature);
        let mut power_scale = 0.0;
        for (element, tet) in problem.mesh.complex().tets.iter().enumerate() {
            if element % 256 == 0 {
                poll(cx)?;
            }
            let nodal = tet
                .iter()
                .map(|&vertex| problem.source.at(vertex as usize).abs())
                .sum::<f64>();
            power_scale += problem.mesh.element_volume(element) * nodal / 4.0;
        }
        let tolerance =
            (config.stop.residual_rtol * power_scale.max(1.0)).max(config.stop.residual_atol);
        if !tolerance.is_finite() || tolerance <= 0.0 {
            return Err(bad("contact watt allowance is not positive and finite"));
        }
        // Inner step-size convergence must not replace the original physical
        // contact residual gate. The source-backed owner checks it independently.
        config.stop.step_atol = 0.0;
        config.stop.residual_rtol = config.stop.residual_rtol.min(1e-13);
        config.stop.residual_atol = tolerance * 0.01;
        config.linear.tolerance = config.linear.tolerance.min(1e-13);
        let coupling = SteadyContactConfig {
            max_iterations: MAX_RESPONSES,
            relaxation: 0.5,
            residual_tolerance_w: tolerance,
            energy_tolerance_w: tolerance,
        };
        let accepted = contacts
            // The next natural/airflow callback starts from its last accepted
            // thermal field. The owner restores current prescribed values and
            // rechecks material and contact support, residual and heat balance.
            .solve_steady(
                cx,
                problem,
                work.last_temperature.as_deref().unwrap_or(&initial),
                config,
                coupling,
            )
            .map_err(lower)?;
        work.record(
            accepted.iterations,
            accepted.nonlinear_iterations,
            accepted.krylov_iterations,
        )?;
        let mut states = Vec::with_capacity(accepted.physical_contacts.states().len());
        for (index, state) in accepted.physical_contacts.states().iter().enumerate() {
            if index % 256 == 0 {
                poll(cx)?;
            }
            states.push(format!(concat!("{{\"interface\":{},\"side_a_face\":{},\"side_b_face\":{},",
                "\"side_a_temperature_k\":{},\"side_b_temperature_k\":{},\"query_temperature_k\":{},",
                "\"resistance_m2_k_w\":{},\"selected_claim\":{},\"usage_receipt\":{}}}"),
                json_string(&state.interface), state.pair.side_a, state.pair.side_b,
                number(state.side_a_temperature_k)?, number(state.side_b_temperature_k)?,
                number(state.interface_temperature_k)?, number(state.resistance.value_m2_k_per_w())?,
                json_string(&state.resistance.receipt().selected.0.to_hex()),
                json_string(&state.resistance.receipt().content_hash().to_hex())));
        }
        let receipt = format!(
            concat!(
                "{{\"model\":\"source-backed-face-mean-temperature\",",
                "\"thermal_responses\":{},\"nonlinear_iterations\":{},\"krylov_iterations\":{},",
                "\"total_thermal_responses\":{},\"total_nonlinear_iterations\":{},\"total_krylov_iterations\":{},",
                "\"residual_norm_w\":{},\"residual_tolerance_w\":{},\"energy_closure_w\":{},",
                "\"energy_tolerance_w\":{},\"face_states\":[{}],\"scope\":{}}}"
            ),
            accepted.iterations,
            accepted.nonlinear_iterations,
            accepted.krylov_iterations,
            work.responses,
            work.nonlinear,
            work.krylov,
            number(accepted.residual_norm_w)?,
            number(tolerance)?,
            number(accepted.energy.closure_w)?,
            number(tolerance)?,
            states.join(","),
            json_string(NO_CLAIM)
        );
        let mut conduction = accepted.conduction;
        // The published physical rows describe this unchanged field with fresh
        // contacts. Inner numerical history remains the counted last response;
        // all outer work is retained separately in the contact receipt above.
        conduction.report.final_residual = accepted.residual_norm_w;
        conduction.report.residual_threshold = tolerance;
        conduction.report.energy = accepted.energy;
        conduction.report.robin_fluxes = accepted.robin_fluxes;
        conduction.report.interface_fluxes = accepted
            .physical_contacts
            .frozen_interfaces()
            .fluxes(&conduction.temperature)
            .map_err(lower)?;
        poll(cx)?;
        work.last_temperature = Some(conduction.temperature.clone());
        Ok((
            super::radiation::SolidSolution::from_conduction(conduction),
            receipt,
        ))
    }
}
