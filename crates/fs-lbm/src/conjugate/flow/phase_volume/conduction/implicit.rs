//! One native implicit enthalpy interval before the original volume transport.

use super::*;
use fs_conduction::enthalpy::conduction::implicit::source::ImplicitHeatSourceLaw;
use fs_conduction::enthalpy::conduction::implicit::{
    ImplicitConductionControls, ImplicitConductionGraph, ImplicitConductionLaw,
    ImplicitConductionReport,
};
use fs_conduction::enthalpy::conduction::{ConductionObserver, ConductiveLink, ConductivityLaw};
use fs_conduction::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};
use fs_math::dd::Dd;
use std::collections::BTreeMap;

pub use fs_conduction::enthalpy::conduction::implicit::{
    ImplicitConductionLinearControls as ImplicitThermalVolumeLinearControls,
    ImplicitConductionWorkspace as ImplicitThermalVolumeWorkspace,
};

mod momentum;
pub use momentum::ImplicitConductiveThermalVolumeMomentumStep;
mod contact;
pub use contact::ThermalVolumeContact;
mod surface;
pub use surface::{ThermalVolumeBathExchange, ThermalVolumeBoundaries, ThermalVolumeBoundaryFace};
mod viscous;
pub use viscous::{
    ImplicitConductiveViscousThermalVolumeMomentumStep,
    ImplicitConjugateViscousThermalVolumeMomentumStep,
    ImplicitThermorheologicalThermalVolumeMomentumStep,
    ImplicitThermoviscousThermalVolumeMomentumStep, ThermalViscosityFeedbackControls,
    ThermalViscosityFeedbackReport,
};
mod solid;
pub use solid::{
    ImplicitConjugateThermalVolumeMomentumStep, ImplicitConjugateThermalVolumeReport,
    ImplicitConjugateThermalVolumeStep, ThermalVolumeSolidCell, ThermalVolumeSolidEndpoint,
};
mod resistance;
pub use resistance::{
    ImplicitConductiveResistiveThermalVolumeMomentumStep,
    ImplicitConjugateResistiveThermalVolumeMomentumStep,
    ImplicitConjugateViscousResistiveThermalVolumeMomentumStep, ImplicitResistanceHeatWork,
    ImplicitThermoviscousResistiveThermalVolumeMomentumStep,
};

/// Bounded native analytic Newton/Krylov work for ONE backward-Euler interval.
/// The existing flow duration and shared enthalpy-defect budget own the clock
/// and total energy gate. There are no hidden thermal substeps or time retries.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitThermalVolumeSolveControls {
    pub linear: ImplicitThermalVolumeLinearControls,
    pub maximum_newton_iterations: usize,
    pub maximum_backtracks: usize,
    pub maximum_residual_evaluations: usize,
    pub cell_tolerance_j: f64,
}

impl Default for ImplicitThermalVolumeSolveControls {
    fn default() -> Self {
        let native = ImplicitConductionControls::default();
        Self {
            linear: ImplicitThermalVolumeLinearControls::default(),
            maximum_newton_iterations: native.maximum_newton_iterations,
            maximum_backtracks: native.maximum_backtracks,
            maximum_residual_evaluations: native.maximum_residual_evaluations,
            cell_tolerance_j: native.cell_tolerance_j,
        }
    }
}

/// Endpoint-evaluated physical heat and its original native nonlinear solve.
#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitConductiveThermalVolumeReport {
    pub conductivity_source_identity: ContentHash,
    pub conductive_faces: usize,
    pub boundary_faces: usize,
    /// Signed inward bath heat by x-, x+, y-, y+, z-, z+ face normal [J].
    pub boundary_heat_j: [f64; 6],
    /// Actual prescribed bath exchanges in cell/face/channel order. Opposing
    /// reservoirs on one shared surface remain separate, with original Kelvin.
    pub bath_exchanges: Vec<ThermalVolumeBathExchange>,
    /// Sum of absolute channel heats, including opposing baths on one face.
    pub absolute_boundary_heat_j: f64,
    pub exchanged_heat_j: f64,
    pub conduction_balance_defect_j: f64,
    pub increment_assembly_defect_j: f64,
    /// Sum of absolute cell mismatches between the native implicit endpoint's
    /// stored change and independently reconstructed physical heat [J].
    pub endpoint_heat_defect_j: f64,
    /// Assembly, endpoint, balance and original volume-transport defects [J].
    pub absolute_enthalpy_defect_j: f64,
    /// Exact native nonlinear work and energy account. Distinct bath faces are
    /// attached as local feedback powers, so their entropy is marked unclassified
    /// by this native report; it is not a complete fluid-plus-bath entropy audit.
    pub solve: ImplicitConductionReport,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ImplicitConductiveThermalVolumeStep {
    pub thermal: ThermalVolumeStep,
    /// Accepted endpoint Fourier heat, before authored increments [J].
    pub conductive_enthalpy_increment_j: Vec<f64>,
    pub report: ImplicitConductiveThermalVolumeReport,
}

struct ImplicitAdvance<T> {
    evolution: T,
    increments: Vec<f64>,
    report: ImplicitConductiveThermalVolumeReport,
    solid_endpoints: Vec<ThermalVolumeSolidEndpoint>,
    solid_authored_heat_j: f64,
    solid_enthalpy_change_j: f64,
    solid_to_fluid_heat_j: f64,
    solid_boundary_heat_j: f64,
    absolute_solid_boundary_heat_j: f64,
}

/// Original-state heat interval, before any pressure/material-flow trial.
struct PreparedImplicitHeat {
    /// Exact native pretransport endpoint, including bound finite solids.
    endpoint_samples: Vec<Option<PhaseConductivitySample>>,
    combined: Vec<f64>,
    increments: Vec<f64>,
    report: ImplicitConductiveThermalVolumeReport,
    solid_endpoints: Vec<ThermalVolumeSolidEndpoint>,
    solid_authored_heat_j: f64,
    solid_enthalpy_change_j: f64,
    solid_to_fluid_heat_j: f64,
    solid_boundary_heat_j: f64,
    absolute_solid_boundary_heat_j: f64,
}

impl FlowField {
    /// Solve source-owned k(T,f), original latent storage and all prescribed
    /// baths at one backward-Euler endpoint, then evolve native thermal volume.
    ///
    /// Actual mass and graded geometry are held during the heat solve. The SAME
    /// borrowed caloric curve and conductivity supply analytic h tangents to the
    /// existing sparse implicit graph. Accepted endpoint face heat joins authored
    /// heat once before the original density/pressure/material transport. Fluid
    /// half distances, separate bath temperatures/resistances and embedded-wall
    /// ownership match [`Self::advance_isotropic_conductive_thermal_volume_with_baths`].
    /// No explicit diffusion CFL is imposed; bounded nonlinear convergence and
    /// source/energy admission remain required. This first-order conduction/flow
    /// split still needs time refinement. It does not infer phase mobility,
    /// pressure work, finite wall storage or simultaneous advection-diffusion.
    /// `workspace` must have one entry per fluid cell; fixed solids add none.
    ///
    /// # Errors
    /// Any original source, geometry, work, nonlinear, heat, volume or cancellation
    /// refusal returns no endpoint and changes none of the supplied inventories.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_isotropic_implicit_conductive_thermal_volume_with_baths(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: &[ThermalVolumeBathFace],
        input: CoexistenceVolumeInput<'_>,
        controls: ConductiveThermalVolumeControls,
        solve: ImplicitThermalVolumeSolveControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
    ) -> Result<ImplicitConductiveThermalVolumeStep, ChtError> {
        let accepted = self.implicit_conduct_then_advance(
            domain,
            conductivity,
            baths,
            input,
            controls,
            solve,
            workspace,
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
            |thermal| thermal,
        )?;
        Ok(ImplicitConductiveThermalVolumeStep {
            thermal: accepted.evolution,
            conductive_enthalpy_increment_j: accepted.increments,
            report: accepted.report,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn implicit_conduct_then_advance<T>(
        &self,
        domain: &VoxelDomain,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: &[ThermalVolumeBathFace],
        input: CoexistenceVolumeInput<'_>,
        controls: ConductiveThermalVolumeControls,
        solve: ImplicitThermalVolumeSolveControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
        advance: impl FnOnce(CoexistenceVolumeInput<'_>) -> Result<T, ChtError>,
        thermal: impl Fn(&T) -> &ThermalVolumeStep,
    ) -> Result<ImplicitAdvance<T>, ChtError> {
        self.implicit_conjugate_then_advance(
            domain,
            conductivity,
            baths.into(),
            &[],
            &[],
            input,
            controls,
            solve,
            workspace,
            gate,
            advance,
            thermal,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn implicit_conjugate_then_advance<T>(
        &self,
        domain: &VoxelDomain,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: ThermalVolumeBoundaries<'_>,
        solids: &[ThermalVolumeSolidCell<'_, '_>],
        contacts: &[ThermalVolumeContact<'_>],
        input: CoexistenceVolumeInput<'_>,
        controls: ConductiveThermalVolumeControls,
        solve: ImplicitThermalVolumeSolveControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
        advance: impl FnOnce(CoexistenceVolumeInput<'_>) -> Result<T, ChtError>,
        thermal: impl Fn(&T) -> &ThermalVolumeStep,
    ) -> Result<ImplicitAdvance<T>, ChtError> {
        let prepared = self.prepare_implicit_conjugate_heat(
            domain,
            conductivity,
            baths,
            solids,
            contacts,
            input,
            controls,
            solve,
            workspace,
            gate,
        )?;
        self.finish_implicit_conjugate_heat(
            prepared,
            conductivity,
            input,
            controls,
            gate,
            advance,
            thermal,
        )
    }

    /// Finish one already accepted native heat interval without replaying it.
    #[allow(clippy::too_many_arguments)]
    fn finish_implicit_conjugate_heat<T>(
        &self,
        prepared: PreparedImplicitHeat,
        conductivity: &PreparedPhaseConductivity<'_>,
        input: CoexistenceVolumeInput<'_>,
        controls: ConductiveThermalVolumeControls,
        gate: &CancelGate,
        advance: impl FnOnce(CoexistenceVolumeInput<'_>) -> Result<T, ChtError>,
        thermal: impl Fn(&T) -> &ThermalVolumeStep,
    ) -> Result<ImplicitAdvance<T>, ChtError> {
        let PreparedImplicitHeat {
            endpoint_samples,
            combined,
            increments,
            mut report,
            solid_endpoints,
            solid_authored_heat_j,
            solid_enthalpy_change_j,
            solid_to_fluid_heat_j,
            solid_boundary_heat_j,
            absolute_solid_boundary_heat_j,
        } = prepared;
        drop(endpoint_samples);
        poll(gate)?;
        let evolution = advance(CoexistenceVolumeInput {
            net_enthalpy_increment_j: &combined,
            ..input
        })?;
        let native = thermal(&evolution);
        let total =
            report.absolute_enthalpy_defect_j + native.volume.report.absolute_enthalpy_defect_j;
        finite("flow.implicit_total_enthalpy_defect", total)?;
        if total
            > controls
                .volume
                .thermal
                .volume
                .transport
                .maximum_enthalpy_defect_j
        {
            return Err(refuse(
                "implicit conduction and volume transport exceed shared enthalpy budget",
            ));
        }
        for (cell, state) in native.volume.phase_states.iter().enumerate() {
            poll_every(gate, cell)?;
            if let Some(state) = state {
                admit_state(conductivity, *state)?;
            }
        }
        poll(gate)?;
        report.absolute_enthalpy_defect_j = total;
        Ok(ImplicitAdvance {
            evolution,
            increments,
            report,
            solid_endpoints,
            solid_authored_heat_j,
            solid_enthalpy_change_j,
            solid_to_fluid_heat_j,
            solid_boundary_heat_j,
            absolute_solid_boundary_heat_j,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn prepare_implicit_conjugate_heat(
        &self,
        domain: &VoxelDomain,
        conductivity: &PreparedPhaseConductivity<'_>,
        baths: ThermalVolumeBoundaries<'_>,
        solids: &[ThermalVolumeSolidCell<'_, '_>],
        contacts: &[ThermalVolumeContact<'_>],
        input: CoexistenceVolumeInput<'_>,
        controls: ConductiveThermalVolumeControls,
        solve: ImplicitThermalVolumeSolveControls,
        workspace: &mut ImplicitThermalVolumeWorkspace,
        gate: &CancelGate,
    ) -> Result<PreparedImplicitHeat, ChtError> {
        let boundary_faces = baths
            .baths
            .len()
            .checked_add(baths.sources.len())
            .filter(|&count| count <= controls.max_conductive_faces)
            .ok_or_else(|| refuse("conductive face-work budget exhausted by thermal boundaries"))?;
        // Retain the original fluid-only admission order and arithmetic when
        // no wall inventory is bound. Conjugate baths need the complete active
        // set so that neither an inactive receiver nor an internal face passes.
        let PreparedInput {
            mut samples,
            baths: fluid_baths,
        } = prepare_input(
            domain,
            conductivity,
            if solids.is_empty() { baths.baths } else { &[] },
            input,
            controls,
            gate,
        )?;
        let wall_at = solid::prepare(domain, solids, &mut samples, gate)?;
        let contact = contact::PreparedContacts::new(
            domain,
            contacts,
            &samples,
            &wall_at,
            controls.max_conductive_faces,
            gate,
        )?;
        let source_faces = baths.sources;
        let baths = if solids.is_empty() {
            fluid_baths
        } else {
            boundary::prepare_for_active_cells(
                domain,
                baths.baths,
                controls.max_conductive_faces,
                |cell| samples[cell].is_some(),
                gate,
            )?
        };
        let wall = |cell: usize| wall_at.get(cell).copied().flatten().map(|i| &solids[i]);
        let source = |cell| wall(cell).map_or(conductivity, |s| s.conductivity);
        let mass = |cell| wall(cell).map_or(input.mass_kg[cell], |s| s.mass_kg);
        let authored = |cell| {
            wall(cell).map_or(input.net_enthalpy_increment_j[cell], |s| {
                s.net_enthalpy_increment_j
            })
        };
        let bath_sources = surface::prepare(
            domain,
            baths,
            source_faces,
            &samples,
            &source,
            controls.max_conductive_faces,
            gate,
        )?;
        // Pool a source only when it retains the SAME original caloric object:
        // native feedback admission checks that borrowed law by address. Equal
        // content in separately allocated curves must keep separate bindings.
        // Pointers are lookup keys only; material IDs follow input order and
        // neither source identities nor returned diagnostics contain addresses.
        let material_key = |law: &PreparedPhaseConductivity<'_>| {
            (law.identity(), core::ptr::from_ref(law.caloric_curve()))
        };
        let mut sources = Vec::new();
        let mut material_ids = BTreeMap::new();
        if !solids.is_empty() {
            sources.push(conductivity);
            material_ids.insert(material_key(conductivity), MaterialId(0));
        }
        for (index, solid) in solids.iter().enumerate() {
            poll_every(gate, index)?;
            let key = material_key(solid.conductivity);
            if !material_ids.contains_key(&key) {
                material_ids.insert(key, MaterialId(sources.len() as u64));
                sources.push(solid.conductivity);
            }
        }
        let transport = controls.volume.thermal.volume.transport;
        let dt = transport.projection.duration_s;
        let budget = transport.maximum_enthalpy_defect_j;
        let n = domain.cell_count();
        let mut events = Vec::new();
        let mut links = Vec::new();
        for (c, sample) in samples.iter().enumerate() {
            poll_every(gate, c)?;
            let Some(sample) = sample else { continue };
            let at = domain.coords(c);
            let center = domain.center(at[0], at[1], at[2]);
            let mut start = center;
            let mut end = center;
            start[0] -= 0.5 * domain.widths(c)[0];
            end[0] += 0.5 * domain.widths(c)[0];
            let power = authored(c) / dt;
            finite("flow.implicit_authored_power", power)?;
            if power == 0.0 && authored(c) != 0.0 {
                return Err(refuse("authored implicit heat power underflowed"));
            }
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(c as u64),
                    material_id: if solids.is_empty() {
                        MaterialId(0)
                    } else {
                        material_ids[&material_key(source(c))]
                    },
                    body_id: 0,
                    segment_id: c as u64,
                    mass_kg: mass(c),
                    specific_h_j_kg: sample.specific_enthalpy_j_kg(),
                    geometry: CellGeometry {
                        start_m: start,
                        end_m: end,
                        volume_m3: domain.volume(c),
                        // Filled from actual link coverage and remaining faces
                        // below, in the native owner's accumulation order.
                        surface_area_m2: 0.0,
                    },
                    boundary: Boundary {
                        ambient_k: sample.temperature_k(),
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: power,
                    },
                }),
            });
        }
        let active = events.len();
        if active == 0 || workspace.cells() != active {
            return Err(refuse(
                "implicit workspace needs one entry per active fluid or finite-solid cell",
            ));
        }
        let mut covered = vec![0.0; n];
        for c in 0..n {
            poll_every(gate, c)?;
            if samples[c].is_none() {
                continue;
            }
            for side in [1, 3, 5] {
                let Some(other) = domain.neighbor(c, side) else {
                    continue;
                };
                if samples[other].is_none() {
                    continue;
                }
                if links.len() >= controls.max_conductive_faces.saturating_sub(boundary_faces) {
                    return Err(refuse("implicit conductive face-work budget exhausted"));
                }
                let id = links.len() as u64;
                let area = domain.face_area(c, side / 2);
                covered[c] += area;
                covered[other] += area;
                events.push(TimedEvent {
                    time_s: 0.0,
                    event: Event::Contact(ThermalLink {
                        id,
                        a: CellId(c as u64),
                        b: CellId(other as u64),
                        area_m2: area,
                        conductance_w_k: 0.0,
                    }),
                });
                links.push(ConductiveLink {
                    link_id: id,
                    half_distances_m: [
                        0.5 * domain.widths(c)[side / 2],
                        0.5 * domain.widths(other)[side / 2],
                    ],
                    contact_resistance_m2_k_w: contact.fixed_resistance(c, other),
                });
            }
        }
        for (index, event) in events.iter_mut().take(active).enumerate() {
            poll_every(gate, index)?;
            if let Event::Birth(cell) = &mut event.event {
                let c = cell.id.0 as usize;
                let mut area = covered[c];
                for side in 0..6 {
                    if domain
                        .neighbor(c, side)
                        .is_none_or(|other| samples[other].is_none())
                    {
                        area += domain.face_area(c, side / 2);
                    }
                }
                finite_positive("flow.implicit_cell_surface", area)?;
                cell.geometry.surface_area_m2 = area;
            }
        }
        let event_count = events.len();
        let mut caloric_laws = Vec::with_capacity(sources.len());
        let mut conductive_laws = Vec::with_capacity(sources.len());
        let mut implicit_laws = Vec::with_capacity(sources.len());
        for (index, &law) in sources.iter().enumerate() {
            poll_every(gate, index)?;
            let id = MaterialId(index as u64);
            caloric_laws.push((id, law.caloric_curve() as &dyn ThermalLaw));
            conductive_laws.push((id, law as &dyn ConductivityLaw));
            implicit_laws.push((id, law as &dyn ImplicitConductionLaw));
        }
        let fluid_caloric = [(
            MaterialId(0),
            conductivity.caloric_curve() as &dyn ThermalLaw,
        )];
        let fluid_conductivity = [(MaterialId(0), conductivity as &dyn ConductivityLaw)];
        let fluid_implicit = [(MaterialId(0), conductivity as &dyn ImplicitConductionLaw)];
        let initial = DepositionSimulation::new(
            events,
            if sources.is_empty() {
                &fluid_caloric
            } else {
                &caloric_laws
            },
            DepositionBudget {
                max_cells: active,
                max_links: links.len(),
                max_events: event_count,
                max_substeps: 1,
            },
        )
        .map_err(|e| native_error(gate, e))?;
        let (initial, _) = initial
            .advance_to(
                0.0,
                StepConfig {
                    max_dt_s: dt,
                    energy_tolerance_j: 0.0,
                },
                || checkpoint(gate),
            )
            .map_err(|e| native_error(gate, e))?;
        let observer = ConductionObserver::new(
            &initial,
            if sources.is_empty() {
                &fluid_conductivity
            } else {
                &conductive_laws
            },
            &links,
            &[],
            || checkpoint(gate),
        )
        .map_err(|e| native_error(gate, e))?;
        let mut graph = ImplicitConductionGraph::new(
            &initial,
            &observer,
            if sources.is_empty() {
                &fluid_implicit
            } else {
                &implicit_laws
            },
            solve.linear,
            || checkpoint(gate),
        )
        .map_err(|e| native_error(gate, e))?;
        let contact_bindings = contact.bindings(&initial, gate)?;
        if !contact_bindings.is_empty() {
            graph = graph
                .with_contacts(&initial, &contact_bindings, || checkpoint(gate))
                .map_err(|e| native_error(gate, e))?;
        }
        if !bath_sources.is_empty() {
            let mut bindings = Vec::with_capacity(bath_sources.len());
            for (index, (&cell, law)) in bath_sources.iter().enumerate() {
                poll_every(gate, index)?;
                bindings.push((CellId(cell as u64), law as &dyn ImplicitHeatSourceLaw));
            }
            graph = graph
                .with_heat_sources(&initial, &bindings, || checkpoint(gate))
                .map_err(|e| native_error(gate, e))?;
        }
        let (endpoint, solve) = graph
            .advance_to(
                &initial,
                dt,
                ImplicitConductionControls {
                    maximum_dt_s: dt,
                    minimum_dt_s: dt,
                    maximum_steps: 1,
                    maximum_halvings: 0,
                    maximum_newton_iterations: solve.maximum_newton_iterations,
                    maximum_backtracks: solve.maximum_backtracks,
                    maximum_residual_evaluations: solve.maximum_residual_evaluations,
                    cell_tolerance_j: solve.cell_tolerance_j,
                    total_tolerance_j: budget,
                },
                workspace,
                || checkpoint(gate),
            )
            .map_err(|e| native_error(gate, e))?;
        let mut final_samples = vec![None; n];
        for (&CellId(cell), state) in endpoint.cells() {
            poll(gate)?;
            final_samples[cell as usize] = Some(
                source(cell as usize)
                    .state_at_specific_enthalpy(state.specific_h_j_kg)
                    .map_err(|e| refuse(e.to_string()))?,
            );
        }
        let mut increments = vec![0.0; n];
        let mut assembly = Sum::default();
        let mut exchanged = Sum::default();
        let mut solid_to_fluid = Sum::default();
        for (index, (&_, link)) in initial.links().iter().enumerate() {
            poll_every(gate, index)?;
            let a = link.a.0 as usize;
            let b = link.b.0 as usize;
            let ka = final_samples[a].unwrap().conductivity_w_m_k();
            let kb = final_samples[b].unwrap().conductivity_w_m_k();
            let rc = if contact.has_sources() {
                contact.resistance_at(
                    a,
                    b,
                    [
                        final_samples[a].unwrap().specific_enthalpy_j_kg(),
                        final_samples[b].unwrap().specific_enthalpy_j_kg(),
                    ],
                    gate,
                )?
            } else {
                links[index].contact_resistance_m2_k_w
            };
            let q = if ka == 0.0 || kb == 0.0 {
                0.0
            } else {
                let geometry = links[index];
                let ra = geometry.half_distances_m[0] / ka;
                let rb = geometry.half_distances_m[1] / kb;
                finite_positive("flow.implicit_half_resistance", ra)?;
                finite_positive("flow.implicit_half_resistance", rb)?;
                let g = if rc == 0.0 {
                    let largest = ra.max(rb);
                    (link.area_m2 / largest) / (1.0 + ra.min(rb) / largest)
                } else {
                    let largest = ra.max(rb).max(rc);
                    (link.area_m2 / largest) / (ra / largest + rb / largest + rc / largest)
                };
                finite_positive("flow.implicit_face_conductance", g)?;
                let delta = final_samples[b].unwrap().temperature_k()
                    - final_samples[a].unwrap().temperature_k();
                let heat = dt * (g * delta);
                if delta != 0.0 && heat == 0.0 {
                    return Err(refuse("nonzero implicit conductive heat underflowed"));
                }
                heat
            };
            finite("flow.implicit_face_heat", q)?;
            for (cell, heat) in [(a, q), (b, -q)] {
                let (value, residual) = accounted_add(increments[cell], heat)?;
                increments[cell] = value;
                assembly.add(residual.abs());
            }
            exchanged.add(q.abs());
            if domain.is_fluid(a) != domain.is_fluid(b) {
                solid_to_fluid.add(if domain.is_fluid(a) { q } else { -q });
            }
        }
        let mut boundary_heat = [Sum::default(); 6];
        let mut bath_exchanges = Vec::with_capacity(boundary_faces);
        let mut absolute_boundary_heat = Sum::default();
        let mut solid_boundary_heat = Sum::default();
        let mut absolute_solid_boundary_heat = Sum::default();
        for (&cell, law) in &bath_sources {
            poll(gate)?;
            let sample = final_samples[cell].unwrap();
            for face in &law.faces {
                let response = face
                    .response(sample, &mut || checkpoint(gate))
                    .map_err(|e| native_error(gate, e))?;
                for (channel, bath) in response.baths[..response.bath_count].iter().enumerate() {
                    let power = -bath.outward_power_w;
                    let heat = dt * power;
                    finite("flow.implicit_bath_heat", heat)?;
                    if power != 0.0 && heat == 0.0 {
                        return Err(refuse("nonzero implicit bath heat underflowed"));
                    }
                    let (value, residual) = accounted_add(increments[cell], heat)?;
                    increments[cell] = value;
                    assembly.add(residual.abs());
                    boundary_heat[face.side].add(heat);
                    absolute_boundary_heat.add(heat.abs());
                    if wall(cell).is_some() {
                        solid_boundary_heat.add(heat);
                        absolute_solid_boundary_heat.add(heat.abs());
                    }
                    bath_exchanges.push(ThermalVolumeBathExchange {
                        cell,
                        face: face.face,
                        channel,
                        temperature_k: bath.temperature_k,
                        heat_j: heat,
                    });
                }
            }
        }
        let mut combined = vec![0.0; n];
        let mut balance = Sum::default();
        let mut endpoint_defect = Sum::default();
        let mut solid_endpoints = Vec::with_capacity(solids.len());
        let mut solid_authored = Sum::default();
        let mut solid_change = Sum::default();
        for c in 0..n {
            poll_every(gate, c)?;
            balance.add(increments[c]);
            let (heat, residual) = accounted_add(authored(c), increments[c])?;
            assembly.add(residual.abs());
            combined[c] = heat;
            if let Some(old) = samples[c] {
                let law = source(c);
                let bounds = law.enthalpy_bounds_j_kg();
                let original = wall(c).map_or_else(|| input.phase_states[c].unwrap(), |s| s.phase);
                let advanced = law
                    .caloric_curve()
                    .advance_specific_energy(original, heat / mass(c))
                    .map_err(|e| refuse(e.to_string()))?;
                if !(bounds[0]..=bounds[1]).contains(&advanced.specific_enthalpy_j_kg()) {
                    return Err(refuse(
                        "implicit heated state leaves original conductivity source",
                    ));
                }
                let change = Dd::from_f64(mass(c))
                    * (Dd::from_f64(final_samples[c].unwrap().specific_enthalpy_j_kg())
                        - Dd::from_f64(old.specific_enthalpy_j_kg()));
                endpoint_defect.add((change - Dd::from_f64(heat)).to_f64().abs());
                if let Some(solid) = wall(c) {
                    solid_authored.add(solid.net_enthalpy_increment_j);
                    solid_change.add(change.to_f64());
                    solid_endpoints.push(ThermalVolumeSolidEndpoint {
                        cell: c,
                        mass_kg: solid.mass_kg,
                        phase: law
                            .caloric_curve()
                            .state_at_specific_enthalpy(
                                final_samples[c].unwrap().specific_enthalpy_j_kg(),
                            )
                            .map_err(|e| refuse(e.to_string()))?,
                        conductivity_source_identity: law.identity(),
                        conductive_enthalpy_increment_j: increments[c],
                    });
                    // Wall storage stays wholly outside the fluid's original
                    // inventory arrays and its conservative material routes.
                    combined[c] = 0.0;
                    increments[c] = 0.0;
                }
            }
        }
        let boundary_heat_j = boundary_heat.map(|sum| sum.value);
        for heat in boundary_heat_j {
            if heat != 0.0 {
                balance.add(-heat);
            }
        }
        let extra_defect = assembly.value + balance.value.abs() + endpoint_defect.value;
        finite("flow.implicit_heat_defect", extra_defect)?;
        finite("flow.implicit_internal_exchange", exchanged.value)?;
        for value in [
            solid_authored.value,
            solid_change.value,
            solid_to_fluid.value,
            solid_boundary_heat.value,
            absolute_solid_boundary_heat.value,
        ] {
            finite("flow.implicit_finite_solid_energy", value)?;
        }
        finite(
            "flow.implicit_absolute_bath_heat",
            absolute_boundary_heat.value,
        )?;
        if extra_defect > budget {
            return Err(refuse(
                "implicit endpoint and heat assembly exceed shared enthalpy budget",
            ));
        }
        poll(gate)?;
        Ok(PreparedImplicitHeat {
            endpoint_samples: final_samples,
            combined,
            increments,
            solid_endpoints,
            solid_authored_heat_j: solid_authored.value,
            solid_enthalpy_change_j: solid_change.value,
            solid_to_fluid_heat_j: solid_to_fluid.value,
            solid_boundary_heat_j: solid_boundary_heat.value,
            absolute_solid_boundary_heat_j: absolute_solid_boundary_heat.value,
            report: ImplicitConductiveThermalVolumeReport {
                conductivity_source_identity: conductivity.identity(),
                conductive_faces: links.len(),
                boundary_faces,
                boundary_heat_j,
                bath_exchanges,
                absolute_boundary_heat_j: absolute_boundary_heat.value,
                exchanged_heat_j: exchanged.value,
                conduction_balance_defect_j: balance.value.abs(),
                increment_assembly_defect_j: assembly.value,
                endpoint_heat_defect_j: endpoint_defect.value,
                absolute_enthalpy_defect_j: extra_defect,
                solve,
            },
        })
    }
}

fn checkpoint(gate: &CancelGate) -> Result<(), DepositionError> {
    if gate.is_requested() {
        Err(DepositionError::Cancelled)
    } else {
        Ok(())
    }
}

fn native_error(gate: &CancelGate, error: impl core::fmt::Display) -> ChtError {
    if gate.is_requested() {
        ChtError::Cancelled
    } else {
        refuse(error.to_string())
    }
}
