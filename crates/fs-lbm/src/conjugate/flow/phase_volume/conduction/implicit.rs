//! One native implicit enthalpy interval before the original volume transport.

use super::*;
use fs_conduction::enthalpy::conduction::implicit::source::{
    ImplicitHeatSourceLaw, ImplicitHeatSourceState,
};
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
        let PreparedInput { samples, baths } =
            prepare_input(domain, conductivity, baths, input, controls, gate)?;
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
            let power = input.net_enthalpy_increment_j[c] / dt;
            finite("flow.implicit_authored_power", power)?;
            if power == 0.0 && input.net_enthalpy_increment_j[c] != 0.0 {
                return Err(refuse("authored implicit heat power underflowed"));
            }
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: CellId(c as u64),
                    material_id: MaterialId(0),
                    body_id: 0,
                    segment_id: c as u64,
                    mass_kg: input.mass_kg[c],
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
                "implicit workspace needs one entry per active fluid cell",
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
                if links.len() >= controls.max_conductive_faces.saturating_sub(baths.len()) {
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
                    contact_resistance_m2_k_w: 0.0,
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
        let initial = DepositionSimulation::new(
            events,
            &[(MaterialId(0), conductivity.caloric_curve())],
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
            &[(MaterialId(0), conductivity as &dyn ConductivityLaw)],
            &links,
            &[],
            || checkpoint(gate),
        )
        .map_err(|e| native_error(gate, e))?;
        let mut bath_sources = BTreeMap::<usize, BathSource<'_, '_>>::new();
        for (&(cell, side), bath) in &baths {
            poll(gate)?;
            bath_sources
                .entry(cell)
                .or_insert_with(|| BathSource {
                    conductivity,
                    faces: Vec::new(),
                })
                .faces
                .push(BoundBath {
                    side,
                    temperature: bath.temperature_k,
                    area: domain.face_area(cell, side / 2),
                    distance: 0.5 * domain.widths(cell)[side / 2],
                    resistance: bath.resistance_m2_k_w,
                });
        }
        let mut graph = ImplicitConductionGraph::new(
            &initial,
            &observer,
            &[(MaterialId(0), conductivity as &dyn ImplicitConductionLaw)],
            solve.linear,
            || checkpoint(gate),
        )
        .map_err(|e| native_error(gate, e))?;
        if !bath_sources.is_empty() {
            let bindings: Vec<_> = bath_sources
                .iter()
                .map(|(&cell, law)| (CellId(cell as u64), law as &dyn ImplicitHeatSourceLaw))
                .collect();
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
                conductivity
                    .state_at_specific_enthalpy(state.specific_h_j_kg)
                    .map_err(|e| refuse(e.to_string()))?,
            );
        }
        let mut increments = vec![0.0; n];
        let mut assembly = Sum::default();
        let mut exchanged = Sum::default();
        for (index, (&_, link)) in initial.links().iter().enumerate() {
            poll_every(gate, index)?;
            let a = link.a.0 as usize;
            let b = link.b.0 as usize;
            let ka = final_samples[a].unwrap().conductivity_w_m_k();
            let kb = final_samples[b].unwrap().conductivity_w_m_k();
            let q = if ka == 0.0 || kb == 0.0 {
                0.0
            } else {
                let geometry = links[index];
                let ra = geometry.half_distances_m[0] / ka;
                let rb = geometry.half_distances_m[1] / kb;
                finite_positive("flow.implicit_half_resistance", ra)?;
                finite_positive("flow.implicit_half_resistance", rb)?;
                let largest = ra.max(rb);
                let g = (link.area_m2 / largest) / (1.0 + ra.min(rb) / largest);
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
        }
        let mut boundary_heat = [Sum::default(); 6];
        let mut absolute_boundary_heat = Sum::default();
        for (&cell, law) in &bath_sources {
            poll(gate)?;
            let sample = final_samples[cell].unwrap();
            for face in &law.faces {
                let power = face
                    .response(sample)
                    .map_err(|e| native_error(gate, e))?
                    .power_w;
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
            }
        }
        let mut combined = vec![0.0; n];
        let mut balance = Sum::default();
        let mut endpoint_defect = Sum::default();
        let bounds = conductivity.enthalpy_bounds_j_kg();
        for c in 0..n {
            poll_every(gate, c)?;
            balance.add(increments[c]);
            let (heat, residual) = accounted_add(input.net_enthalpy_increment_j[c], increments[c])?;
            assembly.add(residual.abs());
            combined[c] = heat;
            if let Some(old) = samples[c] {
                let advanced = conductivity
                    .caloric_curve()
                    .advance_specific_energy(
                        input.phase_states[c].unwrap(),
                        heat / input.mass_kg[c],
                    )
                    .map_err(|e| refuse(e.to_string()))?;
                if !(bounds[0]..=bounds[1]).contains(&advanced.specific_enthalpy_j_kg()) {
                    return Err(refuse(
                        "implicit heated state leaves original conductivity source",
                    ));
                }
                let change = Dd::from_f64(input.mass_kg[c])
                    * (Dd::from_f64(final_samples[c].unwrap().specific_enthalpy_j_kg())
                        - Dd::from_f64(old.specific_enthalpy_j_kg()));
                endpoint_defect.add((change - Dd::from_f64(heat)).to_f64().abs());
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
        let evolution = advance(CoexistenceVolumeInput {
            net_enthalpy_increment_j: &combined,
            ..input
        })?;
        let native = thermal(&evolution);
        let total = extra_defect + native.volume.report.absolute_enthalpy_defect_j;
        finite("flow.implicit_total_enthalpy_defect", total)?;
        if total > budget {
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
        Ok(ImplicitAdvance {
            evolution,
            increments,
            report: ImplicitConductiveThermalVolumeReport {
                conductivity_source_identity: conductivity.identity(),
                conductive_faces: links.len(),
                boundary_faces: baths.len(),
                boundary_heat_j,
                absolute_boundary_heat_j: absolute_boundary_heat.value,
                exchanged_heat_j: exchanged.value,
                conduction_balance_defect_j: balance.value.abs(),
                increment_assembly_defect_j: assembly.value,
                endpoint_heat_defect_j: endpoint_defect.value,
                absolute_enthalpy_defect_j: total,
                solve,
            },
        })
    }
}

struct BoundBath {
    side: usize,
    temperature: f64,
    area: f64,
    distance: f64,
    resistance: f64,
}

impl BoundBath {
    fn response(
        &self,
        sample: PhaseConductivitySample,
    ) -> Result<ImplicitHeatSourceState, DepositionError> {
        let k = sample.conductivity_w_m_k();
        let dk = sample.conductivity_enthalpy_derivative_kg_w_m_k_j();
        if k == 0.0 {
            // The original prepared source refuses zero k with a nonzero selected
            // tangent. Do not silently erase the finite one-sided dG/dh limit.
            if dk != 0.0 {
                return Err(DepositionError::MaterialDomain);
            }
            return Ok(ImplicitHeatSourceState {
                power_w: 0.0,
                power_enthalpy_derivative_w_kg_j: 0.0,
            });
        }
        let bulk = self.distance / k;
        if !bulk.is_finite() || bulk <= 0.0 {
            return Err(DepositionError::Arithmetic);
        }
        let largest = bulk.max(self.resistance);
        let scaled_total = 1.0 + bulk.min(self.resistance) / largest;
        let g = (self.area / largest) / scaled_total;
        let delta = self.temperature - sample.temperature_k();
        let dg_term = if dk == 0.0 || delta == 0.0 {
            0.0
        } else {
            g * ((bulk / largest) / scaled_total) * (dk / k) * delta
        };
        let power = g * delta;
        let derivative = dg_term - g * sample.temperature_enthalpy_derivative_kg_k_j();
        if !g.is_finite()
            || g <= 0.0
            || !power.is_finite()
            || !derivative.is_finite()
            || (delta != 0.0 && power == 0.0)
        {
            return Err(DepositionError::Arithmetic);
        }
        Ok(ImplicitHeatSourceState {
            power_w: power,
            power_enthalpy_derivative_w_kg_j: derivative,
        })
    }
}

struct BathSource<'a, 'm> {
    conductivity: &'a PreparedPhaseConductivity<'m>,
    faces: Vec<BoundBath>,
}

impl ImplicitHeatSourceLaw for BathSource<'_, '_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.conductivity.caloric_curve()
    }
    fn state_at_specific_enthalpy(
        &self,
        h: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitHeatSourceState, DepositionError> {
        checkpoint()?;
        let sample = self
            .conductivity
            .state_at_specific_enthalpy(h)
            .map_err(|_| DepositionError::MaterialDomain)?;
        let mut power = Sum::default();
        let mut derivative = Sum::default();
        for face in &self.faces {
            let state = face.response(sample)?;
            power.add(state.power_w);
            derivative.add(state.power_enthalpy_derivative_w_kg_j);
        }
        if !power.value.is_finite() || !derivative.value.is_finite() {
            return Err(DepositionError::Arithmetic);
        }
        Ok(ImplicitHeatSourceState {
            power_w: power.value,
            power_enthalpy_derivative_w_kg_j: derivative.value,
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
