//! Temperature-dependent interface sources through the original nonlinear wall.
//! Only contact query temperatures are iterated; bulk integration and shooting
//! remain owned by the original source-free and internally generated wall solvers.

use fs_exec::Cx;
use fs_matdb::{ClaimSelection, InterfaceSystemCard, QueryPoint};

use super::{
    LayeredWall, NonlinearGeneratedLayeredWallSolution, NonlinearGenerationOptions,
    NonlinearLayeredWallSolution, NonlinearWallOptions, checkpoint, finite, invalid, positive,
};
use crate::ConductionError;
use crate::interface::InterfaceResistance;
use crate::interface::temperature::TemperatureContactCoordinates;

/// A selected source law at one explicitly indexed interface of a planar stack.
/// Side A is the physical left layer and side B the right layer, even when heat
/// flows right-to-left. All non-temperature source coordinates remain fixed.
#[derive(Debug, Clone)]
pub struct LayeredTemperatureContact<'a> {
    interface_index: usize,
    name: String,
    card: &'a InterfaceSystemCard,
    reference_point: QueryPoint,
    coordinates: TemperatureContactCoordinates,
    reference: InterfaceResistance,
}

impl<'a> LayeredTemperatureContact<'a> {
    /// Bind the interface after layer `interface_index`. The selected claim is
    /// pinned for every subsequent query; source support is never extended.
    /// The original wall validates the index when this binding is used.
    pub fn new(
        interface_index: usize,
        name: impl Into<String>,
        card: &'a InterfaceSystemCard,
        reference_point: &QueryPoint,
        coordinates: TemperatureContactCoordinates,
        selection: ClaimSelection,
    ) -> Result<Self, ConductionError> {
        let name = name.into();
        if name.trim().is_empty() || name.len() > 128 {
            return Err(invalid("contact name must contain 1..=128 UTF-8 bytes"));
        }
        coordinates.admit_reference(&name, reference_point)?;
        let reference = match selection {
            ClaimSelection::Policy(policy) => {
                InterfaceResistance::from_card(&name, card, reference_point, policy)?
            }
            ClaimSelection::Pinned(selected) => {
                InterfaceResistance::from_card_pinned(&name, card, reference_point, selected)?
            }
        };
        Ok(Self {
            interface_index,
            name,
            card,
            reference_point: reference_point.clone(),
            coordinates,
            reference,
        })
    }

    /// Original selected interface answer and complete reference conditions.
    #[must_use]
    pub const fn reference_resistance(&self) -> &InterfaceResistance {
        &self.reference
    }

    fn at(&self, point: &QueryPoint) -> Result<InterfaceResistance, ConductionError> {
        InterfaceResistance::from_card_pinned(
            &self.name,
            self.card,
            point,
            self.reference.receipt().selected,
        )
    }
}

/// Work limits for source-contact feedback around the original wall solver.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TemperatureContactWallOptions {
    /// Existing bulk shooting and physical flux/temperature tolerances. The
    /// same flux tolerance also gates freshly evaluated contact laws.
    pub solver: NonlinearWallOptions,
    /// Maximum complete wall responses, in 1..=256, including discarded ones.
    pub max_contact_iterations: usize,
    /// Convex update of source temperature coordinates, in (0,1]. Resistance
    /// values, returned temperatures and source receipts are never blended.
    pub relaxation: f64,
}

/// Applied and actual source states at one physical interface.
#[derive(Debug, Clone, PartialEq)]
pub struct LayeredContactState {
    /// Interface after this layer, in original left-to-right stack order.
    pub interface_index: usize,
    /// Original caller-supplied interface name.
    pub name: String,
    /// Explicit source-temperature convention, without sorting hot/cold sides.
    pub coordinates: TemperatureContactCoordinates,
    /// Actual left and right interface temperatures from the returned solve, K.
    pub temperatures_k: [f64; 2],
    /// Original source answer used by the last bulk response.
    pub applied: InterfaceResistance,
    /// Same pinned source, freshly queried at the returned interface temperatures.
    pub physical: InterfaceResistance,
    /// (T_left - T_right) / R''_physical, W/m2.
    pub physical_flux_w_per_m2: f64,
    /// Original local outgoing flux of the layer to the left, W/m2. It is
    /// common across a source-free stack; generation changes it between layers.
    pub solved_flux_w_per_m2: f64,
    /// Physical contact flux minus the local solved interface flux, W/m2.
    pub flux_residual_w_per_m2: f64,
    /// Left outgoing minus right incoming flux at this interface, W/m2.
    pub continuity_residual_w_per_m2: f64,
}

/// One conservative wall response with independently checked contact physics.
#[derive(Debug, Clone, PartialEq)]
pub struct TemperatureContactWallSolution {
    /// Original bulk result, including all material sources. Its retained wall
    /// contains the applied contact coefficients; actual answers are below.
    pub conduction: NonlinearLayeredWallSolution,
    /// Fresh source answers and actual interface temperatures, in binding order.
    pub contacts: Vec<LayeredContactState>,
    /// Largest absolute bulk/film/contact flux residual, W/m2, using the actual
    /// contact sources at the unchanged returned temperatures.
    pub maximum_flux_residual_w_per_m2: f64,
    /// Count of complete native wall responses, including rejected candidates.
    pub contact_iterations: usize,
    /// Sum of shooting trials in every native response.
    pub total_shooting_iterations: usize,
    /// Original caller-declared work, relaxation and residual limits.
    pub options: TemperatureContactWallOptions,
}

/// Original generated-bulk limits and bounded contact-source feedback.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeneratedTemperatureContactWallOptions {
    /// Original material envelope, shooting budget and physical residual limits.
    pub generation: NonlinearGenerationOptions,
    /// Maximum complete generated-wall responses, in 1..=256.
    pub max_contact_iterations: usize,
    /// Convex update of source temperature coordinates, in (0,1].
    pub relaxation: f64,
}

/// Internally generated wall with original source and contact heat accounts.
#[derive(Debug, Clone, PartialEq)]
pub struct GeneratedTemperatureContactWallSolution {
    /// Original bulk response, including local fluxes, interior extrema, both
    /// reservoir heat rates, and signed generation counted exactly once.
    pub conduction: NonlinearGeneratedLayeredWallSolution,
    /// Applied and actual source answers in the caller's binding order.
    pub contacts: Vec<LayeredContactState>,
    /// Largest absolute bulk/film/local-contact flux residual, W/m2.
    pub maximum_flux_residual_w_per_m2: f64,
    /// All native wall responses, including discarded contact-source trials.
    pub contact_iterations: usize,
    /// Sum of shooting trials in every native response.
    pub total_shooting_iterations: usize,
    /// Unchanged caller-declared generation and contact controls.
    pub options: GeneratedTemperatureContactWallOptions,
}

// The two native bulk owners retain their equations, source admission, boundary
// gates and profiles. Only this small view is needed by the shared contact loop.
trait ContactWallResponse {
    fn shooting_iterations(&self) -> usize;
    fn maximum_flux_residual(&self) -> f64;
    fn flux_scale(&self) -> f64;
    fn interface_state(&self, index: usize) -> ([f64; 2], [f64; 2]);
}

impl ContactWallResponse for NonlinearLayeredWallSolution {
    fn shooting_iterations(&self) -> usize {
        self.iterations
    }
    fn maximum_flux_residual(&self) -> f64 {
        self.maximum_flux_residual_w_per_m2
    }
    fn flux_scale(&self) -> f64 {
        self.heat_flux_w_per_m2.abs()
    }
    fn interface_state(&self, index: usize) -> ([f64; 2], [f64; 2]) {
        (
            [
                self.layer_temperatures[index].right_k,
                self.layer_temperatures[index + 1].left_k,
            ],
            [self.heat_flux_w_per_m2; 2],
        )
    }
}

impl ContactWallResponse for NonlinearGeneratedLayeredWallSolution {
    fn shooting_iterations(&self) -> usize {
        self.iterations
    }
    fn maximum_flux_residual(&self) -> f64 {
        self.maximum_flux_residual_w_per_m2
    }
    fn flux_scale(&self) -> f64 {
        self.layer_temperatures
            .iter()
            .flat_map(|row| row.face_fluxes_w_per_m2)
            .map(f64::abs)
            .fold(0.0, f64::max)
    }
    fn interface_state(&self, index: usize) -> ([f64; 2], [f64; 2]) {
        let left = &self.layer_temperatures[index];
        let right = &self.layer_temperatures[index + 1];
        (
            [left.face_temperatures_k[1], right.face_temperatures_k[0]],
            [left.face_fluxes_w_per_m2[1], right.face_fluxes_w_per_m2[0]],
        )
    }
}

struct ContactFeedback<S> {
    conduction: S,
    contacts: Vec<LayeredContactState>,
    maximum_flux_residual_w_per_m2: f64,
    contact_iterations: usize,
    total_shooting_iterations: usize,
}

impl LayeredWall {
    /// Solve nonlinear bulk conduction and source-dependent contact resistance
    /// at their own distinct, actual interface temperatures.
    ///
    /// Listed bindings replace only their explicitly indexed contact values;
    /// other contacts and both surface films retain the original wall inputs.
    /// Every trial delegates bulk k(T), source-domain admission and common-flux
    /// shooting to the existing nonlinear wall owner. The working stack and
    /// vector storage are reused; only contact coefficients and temperature
    /// coordinates change between responses. The original inputs are immutable.
    ///
    /// The reference queries initialize the iteration. Only declared source
    /// temperature coordinates are relaxed. Acceptance re-queries every pinned
    /// contact law at the unchanged solved temperatures and checks its physical
    /// flux against the common stack flux. Equal reservoir temperatures query
    /// the actual equilibrium temperature before solving, so zero heat cannot
    /// conceal an incorrect local resistance or transmittance.
    ///
    /// This is a bounded Picard composition, not a uniqueness or convergence
    /// guarantee for arbitrary source feedback. Out-of-domain queries refuse;
    /// there is no extrapolation, claim switching or fallback to perfect contact.
    /// Source admission covers visited points, not continuous paths through
    /// sparse multiaxis support. Geometry, pressure, phase and other nonthermal
    /// context stay fixed. Generation, radial geometry, transient storage and
    /// uncertainty propagation are outside this entry point.
    pub fn solve_temperature_dependent_with_contacts(
        &self,
        cx: &Cx<'_>,
        reservoir_temperatures_k: [f64; 2],
        area_m2: f64,
        contacts: &[LayeredTemperatureContact<'_>],
        options: TemperatureContactWallOptions,
    ) -> Result<TemperatureContactWallSolution, ConductionError> {
        let result = self.solve_contact_feedback(
            cx,
            reservoir_temperatures_k,
            area_m2,
            contacts,
            options,
            true,
            |wall| {
                wall.solve_temperature_dependent(
                    cx,
                    reservoir_temperatures_k,
                    area_m2,
                    options.solver,
                )
            },
        )?;
        Ok(TemperatureContactWallSolution {
            conduction: result.conduction,
            contacts: result.contacts,
            maximum_flux_residual_w_per_m2: result.maximum_flux_residual_w_per_m2,
            contact_iterations: result.contact_iterations,
            total_shooting_iterations: result.total_shooting_iterations,
            options,
        })
    }

    /// Close original temperature-dependent contacts around the existing
    /// internally generated nonlinear wall solve. Signed uniform layer sources
    /// can include independently supplied volumetric Joule loss; no electrical
    /// or temperature-dependent generation law is inferred here.
    ///
    /// Each interface uses its own outgoing/incoming flux from the generated
    /// bulk response. Both native sides must agree with the actual contact law
    /// under their local flux tolerances. Equal boundary temperatures, or zero
    /// NET generation, do not imply an isothermal or zero-flux interior. Only a
    /// complete all-zero source vector enables equilibrium initialization.
    ///
    /// The unchanged bulk owner admits interior extrema and both boundaries.
    /// Its original reservoir-minus-generation energy residual must also meet
    /// the flux tolerance at the largest native face flux, multiplied by area.
    /// Generation is passed once to that owner and never deposited again at a
    /// contact. All source receipts, work limits and cancellation remain active.
    /// As with the source-free entry point, no global convergence, extrapolation,
    /// radial geometry or transient contact feedback is implied.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_temperature_dependent_with_generation_and_contacts(
        &self,
        cx: &Cx<'_>,
        reservoir_temperatures_k: [f64; 2],
        area_m2: f64,
        volumetric_generation_w_m3: &[f64],
        contacts: &[LayeredTemperatureContact<'_>],
        options: GeneratedTemperatureContactWallOptions,
    ) -> Result<GeneratedTemperatureContactWallSolution, ConductionError> {
        let source_free = volumetric_generation_w_m3.len() == self.layers.len()
            && volumetric_generation_w_m3
                .iter()
                .all(|&source| source == 0.0);
        let feedback = TemperatureContactWallOptions {
            solver: options.generation.solver,
            max_contact_iterations: options.max_contact_iterations,
            relaxation: options.relaxation,
        };
        let result = self.solve_contact_feedback(
            cx,
            reservoir_temperatures_k,
            area_m2,
            contacts,
            feedback,
            source_free,
            |wall| {
                wall.solve_temperature_dependent_with_generation(
                    cx,
                    reservoir_temperatures_k,
                    area_m2,
                    volumetric_generation_w_m3,
                    options.generation,
                )
            },
        )?;
        let energy_tolerance_w = finite(
            "generated wall energy tolerance",
            area_m2
                * options
                    .generation
                    .solver
                    .flux_tolerance(result.conduction.flux_scale())?,
        )?;
        if result.conduction.energy_residual_w.abs() > energy_tolerance_w {
            return Err(invalid(
                "generated contact wall reservoir and source heat do not close",
            ));
        }
        checkpoint(cx, result.contact_iterations)?;
        Ok(GeneratedTemperatureContactWallSolution {
            conduction: result.conduction,
            contacts: result.contacts,
            maximum_flux_residual_w_per_m2: result.maximum_flux_residual_w_per_m2,
            contact_iterations: result.contact_iterations,
            total_shooting_iterations: result.total_shooting_iterations,
            options,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_contact_feedback<S: ContactWallResponse>(
        &self,
        cx: &Cx<'_>,
        reservoir_temperatures_k: [f64; 2],
        area_m2: f64,
        contacts: &[LayeredTemperatureContact<'_>],
        options: TemperatureContactWallOptions,
        source_free: bool,
        solve_bulk: impl Fn(&LayeredWall) -> Result<S, ConductionError>,
    ) -> Result<ContactFeedback<S>, ConductionError> {
        checkpoint(cx, 0)?;
        options.solver.validate()?;
        positive("wall area", area_m2)?;
        for temperature in reservoir_temperatures_k {
            positive("reservoir absolute temperature", temperature)?;
        }
        if options.max_contact_iterations == 0
            || options.max_contact_iterations > 256
            || !options.relaxation.is_finite()
            || options.relaxation <= 0.0
            || options.relaxation > 1.0
        {
            return Err(invalid(
                "contact iterations must lie in 1..=256 and relaxation in (0,1]",
            ));
        }
        if contacts.len() > self.contact_resistances_m2_k_per_w.len() {
            return Err(invalid(
                "source contact bindings exceed the number of wall interfaces",
            ));
        }
        let mut indices = std::collections::BTreeSet::new();
        let mut queries = Vec::with_capacity(contacts.len());
        for contact in contacts {
            if contact.interface_index >= self.contact_resistances_m2_k_per_w.len()
                || !indices.insert(contact.interface_index)
            {
                return Err(invalid(
                    "source contact indices must be unique and name existing interfaces",
                ));
            }
            queries.push(
                if source_free && reservoir_temperatures_k[0] == reservoir_temperatures_k[1] {
                    contact.coordinates.query_point_at(
                        &contact.name,
                        &contact.reference_point,
                        reservoir_temperatures_k,
                    )?
                } else {
                    contact.reference_point.clone()
                },
            );
        }
        let mut working = self.clone();
        let mut applied = Vec::with_capacity(contacts.len());
        let mut states = Vec::with_capacity(contacts.len());
        let mut physical_points = Vec::with_capacity(contacts.len());
        let mut total_shooting_iterations = 0;
        let mut last_residual = 0.0;
        for iteration in 1..=options.max_contact_iterations {
            checkpoint(cx, iteration)?;
            for (index, (contact, point)) in contacts.iter().zip(&queries).enumerate() {
                checkpoint(cx, index)?;
                let source = contact.at(point)?;
                working.contact_resistances_m2_k_per_w[contact.interface_index] =
                    source.value_m2_k_per_w();
                applied.push(source);
            }
            let conduction = solve_bulk(&working)?;
            total_shooting_iterations += conduction.shooting_iterations();
            let mut maximum_residual = conduction.maximum_flux_residual();
            let mut converged =
                maximum_residual <= options.solver.flux_tolerance(conduction.flux_scale())?;
            states.clear();
            physical_points.clear();
            for (index, (contact, applied)) in contacts.iter().zip(applied.drain(..)).enumerate() {
                checkpoint(cx, index)?;
                let i = contact.interface_index;
                let (temperatures_k, interface_fluxes) = conduction.interface_state(i);
                let point = contact.coordinates.query_point_at(
                    &contact.name,
                    &contact.reference_point,
                    temperatures_k,
                )?;
                let physical = contact.at(&point)?;
                let physical_flux_w_per_m2 = finite(
                    "physical contact flux",
                    (temperatures_k[0] - temperatures_k[1]) / physical.value_m2_k_per_w(),
                )?;
                let flux_residual_w_per_m2 = finite(
                    "physical contact flux residual",
                    physical_flux_w_per_m2 - interface_fluxes[0],
                )?;
                let right_residual = finite(
                    "physical contact downstream flux residual",
                    physical_flux_w_per_m2 - interface_fluxes[1],
                )?;
                let continuity_residual_w_per_m2 = finite(
                    "original local contact continuity",
                    interface_fluxes[0] - interface_fluxes[1],
                )?;
                let tolerance = [
                    options.solver.flux_tolerance(interface_fluxes[0])?,
                    options.solver.flux_tolerance(interface_fluxes[1])?,
                ];
                converged &= flux_residual_w_per_m2.abs() <= tolerance[0]
                    && right_residual.abs() <= tolerance[1]
                    && continuity_residual_w_per_m2.abs() <= tolerance[0].max(tolerance[1]);
                maximum_residual = maximum_residual
                    .max(flux_residual_w_per_m2.abs())
                    .max(right_residual.abs())
                    .max(continuity_residual_w_per_m2.abs());
                states.push(LayeredContactState {
                    interface_index: i,
                    name: contact.name.clone(),
                    coordinates: contact.coordinates.clone(),
                    temperatures_k,
                    applied,
                    physical,
                    physical_flux_w_per_m2,
                    solved_flux_w_per_m2: interface_fluxes[0],
                    flux_residual_w_per_m2,
                    continuity_residual_w_per_m2,
                });
                physical_points.push(point);
            }
            if converged {
                checkpoint(cx, iteration)?;
                return Ok(ContactFeedback {
                    conduction,
                    contacts: states,
                    maximum_flux_residual_w_per_m2: maximum_residual,
                    contact_iterations: iteration,
                    total_shooting_iterations,
                });
            }
            last_residual = maximum_residual;
            if iteration == options.max_contact_iterations {
                break;
            }
            for ((contact, query), physical) in
                contacts.iter().zip(&mut queries).zip(&physical_points)
            {
                let mut next = query.clone();
                for axis in contact.coordinates.axes() {
                    let old = query.axes()[axis];
                    let target = physical.axes()[axis];
                    let temperature = options.relaxation.mul_add(target - old, old);
                    next = next.with(axis, temperature).map_err(|error| {
                        invalid(&format!("contact temperature update refused: {error}"))
                    })?;
                }
                *query = next;
            }
        }
        Err(invalid(&format!(
            "source-contact wall exhausted {} responses; physical flux residual {last_residual:e} W/m2",
            options.max_contact_iterations,
        )))
    }
}
