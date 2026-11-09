//! Curved, polarized refraction between original nonuniform P1 thermal cells.
//!
//! Native clearance-controlled GRIN integration owns every incoming and outgoing
//! arc, transverse frame, Beer integral, phase length and group length. The two
//! original optical laws meet at a numerically localized contact: no temperature
//! average, position projection, unintegrated offset or permanent face skip is
//! used. A ray can bend back toward the face it just left.
//!
//! The selected material's original affine P1 temperature polynomial supports
//! only the admitted numerical departure collar outside its element. This is
//! bounded by both a thin metric slab and a ball excluding every other face of
//! both cells. It is a local numerical interface approximation, not a physical
//! exterior temperature field. `contact_resolution_m` reports its localization
//! scale; neither that number nor the ODE estimates certify global path error.

use super::*;
use crate::gradient_index::event::{
    BoundaryOracle, BoundaryQuery, BoundarySide, BoundaryWitness, DepartureNeighborhood,
    DepartureRequest, GradientBoundaryPropagation,
};
use fs_ivl::Interval;

/// Arc, local numerical accuracy and whole-request work limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct P1CurvedContactOptions {
    /// Attempted ODE steps and boundary queries are shared across the incoming
    /// leg, collar preparation, and both nonzero outgoing legs. Rejected steps
    /// count. A chained request does not charge an already retained history again.
    pub event: GradientEventOptions,
    pub max_incident_distance_m: f64,
    pub max_outgoing_distance_m: f64,
}

impl Default for P1CurvedContactOptions {
    fn default() -> Self {
        Self {
            event: GradientEventOptions::default(),
            max_incident_distance_m: 1.0,
            max_outgoing_distance_m: 1.0,
        }
    }
}

impl P1CurvedContactOptions {
    fn validate(self) -> Result<(), P1ContactError> {
        self.event.validate()?;
        if [self.max_incident_distance_m, self.max_outgoing_distance_m]
            .into_iter()
            .any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(GradientIndexError::InvalidOptions.into());
        }
        Ok(())
    }
}

/// Original paired faces, with no isothermal restriction on either element.
/// Preparation validates the same local geometry and source ownership as the
/// exact owner. It does not certify the global embedding of the thermal mesh.
pub struct P1CurvedContactInterface<'a> {
    contact: P1ContactInterface<'a>,
}

/// Both nonzero power branches at one numerical material discontinuity.
#[derive(Clone, Copy)]
pub struct P1CurvedContactSplit<'a> {
    pub incident: GradientIndexPropagation,
    pub incident_stokes: Stokes,
    pub incident_frame: PolarizationFrame,
    pub fresnel: FresnelDielectric,
    /// Positive film work at THIS interface relative to the original launch.
    pub coating_absorbed_power_fraction: f64,
    pub reflected: Option<P1CurvedContactBranch<'a>>,
    pub transmitted: Option<P1CurvedContactBranch<'a>>,
    /// Unique work in this request, including rejected trials on either branch.
    pub attempted_steps: u32,
    /// Unique geometry work, including one paired-collar preparation query.
    pub boundary_evaluations: u32,
    /// Largest native incidence-scaled numerical contact scale on these paths.
    /// This is not an interval enclosure of the true interface or global error.
    pub contact_resolution_m: f64,
}

/// Original-material branch with a private native numerical arrival witness.
#[derive(Clone, Copy)]
pub struct P1CurvedContactBranch<'a> {
    propagation: GradientIndexPropagation,
    stokes: Stokes,
    frame: PolarizationFrame,
    binding: &'a Binding<'a>,
    wavelength_nm: f64,
    termination: BoundaryTermination,
    contact_resolution_m: f64,
}

impl P1CurvedContactBranch<'_> {
    pub const fn propagation(&self) -> &GradientIndexPropagation {
        &self.propagation
    }
    pub const fn power_fraction(&self) -> f64 {
        self.stokes.i
    }
    pub const fn stokes(&self) -> Stokes {
        self.stokes
    }
    pub const fn frame(&self) -> PolarizationFrame {
        self.frame
    }
    pub fn element_index(&self) -> usize {
        self.binding.temperature_source().element_index()
    }
    pub const fn termination(&self) -> BoundaryTermination {
        self.termination
    }
    pub const fn arrival_face(&self) -> Option<usize> {
        match self.termination {
            BoundaryTermination::Surface(witness) => Some(witness.feature_id),
            BoundaryTermination::DistanceLimit => None,
        }
    }
    pub const fn contact_resolution_m(&self) -> f64 {
        self.contact_resolution_m
    }
}

impl<'a> P1CurvedContactInterface<'a> {
    pub fn try_new(
        cx: &Cx<'_>,
        bindings: [&'a Binding<'a>; 2],
        local_faces: [usize; 2],
    ) -> Result<Self, P1ContactError> {
        Ok(Self {
            contact: P1ContactInterface::prepare(cx, bindings, local_faces, false)?,
        })
    }

    /// Preserve the original films and side-0-to-side-1 physical ordering.
    #[must_use]
    pub fn with_coating_stack(mut self, coating: DielectricCoatingStack<'a>) -> Self {
        self.contact = self.contact.with_coating_stack(coating);
        self
    }

    pub const fn coating_stack(&self) -> Option<DielectricCoatingStack<'a>> {
        self.contact.coating_stack()
    }

    #[allow(clippy::too_many_arguments)]
    pub fn split(
        &self,
        cx: &Cx<'_>,
        incident_side: usize,
        ray: Ray,
        transverse_axis: Vec3,
        wavelength_nm: f64,
        options: P1CurvedContactOptions,
    ) -> Result<P1CurvedContactSplit<'a>, P1ContactError> {
        self.split_polarized(
            cx,
            incident_side,
            ray,
            transverse_axis,
            Stokes::unpolarized(1.0),
            wavelength_nm,
            options,
        )
    }

    /// Launch strictly inside the selected original element. Stokes I must be
    /// one; every returned component is relative to this original launch.
    #[allow(clippy::too_many_arguments)]
    pub fn split_polarized(
        &self,
        cx: &Cx<'_>,
        incident_side: usize,
        ray: Ray,
        transverse_axis: Vec3,
        polarization: Stokes,
        wavelength_nm: f64,
        options: P1CurvedContactOptions,
    ) -> Result<P1CurvedContactSplit<'a>, P1ContactError> {
        poll(cx)?;
        options.validate()?;
        if polarization.i != 1.0
            || Stokes::try_new(
                polarization.i,
                polarization.q,
                polarization.u,
                polarization.v,
            )
            .is_none()
        {
            return Err(BeamError::InvalidStokes.into());
        }
        let side = self
            .contact
            .sides
            .get(incident_side)
            .ok_or(P1ContactError::InvalidLaunch)?;
        if side
            .binding
            .temperature_source()
            .barycentric(ray.origin)
            .map_err(|_| P1ContactError::InvalidLaunch)?
            .iter()
            .any(|&value| value <= 0.0)
        {
            return Err(P1ContactError::InvalidLaunch);
        }
        let field = side.binding.prepare(wavelength_nm)?;
        let mut boundary = side.binding.temperature_source().boundary_oracle(
            side.binding.temperature_source().element_index(),
            options.event.surface_tolerance_m,
        )?;
        let event = propagate_to_boundary_with_control(
            &field,
            &mut boundary,
            ray,
            transverse_axis,
            wavelength_nm,
            options.max_incident_distance_m,
            None,
            options.event,
            || cx.checkpoint().is_ok(),
        )?;
        let mut work = Work::default();
        work.record(&event)?;
        let witness = match event.termination {
            BoundaryTermination::Surface(witness) => witness,
            BoundaryTermination::DistanceLimit => {
                return Err(P1ContactError::IncidentDistanceLimit);
            }
        };
        if witness.feature_id != side.face {
            return Err(P1ContactError::ContactNotFirstFace {
                local_face: witness.feature_id,
            });
        }
        self.split_hit(
            cx,
            incident_side,
            event.propagation,
            scale_stokes(polarization, event.propagation.power_transmittance()),
            witness,
            wavelength_nm,
            event.contact_resolution_m,
            work,
            options,
        )
    }

    /// Continue only a retained arrival on this original face/material/snapshot.
    /// A numerical boundary witness is never promoted to an exact mesh hit.
    pub fn split_arrival(
        &self,
        cx: &Cx<'_>,
        arrival: &P1CurvedContactBranch<'a>,
        options: P1CurvedContactOptions,
    ) -> Result<P1CurvedContactSplit<'a>, P1ContactError> {
        poll(cx)?;
        options.validate()?;
        let BoundaryTermination::Surface(witness) = arrival.termination else {
            return Err(P1ContactError::ForeignArrival);
        };
        let side = self
            .contact
            .sides
            .iter()
            .position(|side| {
                side.face == witness.feature_id && same_region(side.binding, arrival.binding)
            })
            .ok_or(P1ContactError::ForeignArrival)?;
        self.split_hit(
            cx,
            side,
            arrival.propagation,
            arrival.stokes,
            witness,
            arrival.wavelength_nm,
            arrival.contact_resolution_m,
            Work::default(),
            options,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn split_hit(
        &self,
        cx: &Cx<'_>,
        incident_side: usize,
        incident: GradientIndexPropagation,
        incident_stokes: Stokes,
        witness: BoundaryWitness,
        wavelength_nm: f64,
        previous_resolution_m: f64,
        mut work: Work,
        options: P1CurvedContactOptions,
    ) -> Result<P1CurvedContactSplit<'a>, P1ContactError> {
        poll(cx)?;
        work.charge_geometry(options.event)?;
        let collar = ContactCollar::prepare(
            &self.contact,
            incident_side,
            incident.ray.origin,
            witness,
            options.event,
        )?;
        let sources = self.contact.sides.each_ref().map(|side| CollarTemperature {
            element: *side.binding.temperature_source(),
            face: side.face,
            collar,
        });
        let a = &self.contact.sides[incident_side];
        let b = &self.contact.sides[1 - incident_side];
        let source = a.binding.with_temperature_source(&sources[incident_side]);
        let target = b.binding.with_temperature_source(&sources[1 - incident_side]);
        // Prepare each wavelength law once for both its interface sample and
        // outgoing ODE stages. Temperature and its gradient remain spatially live.
        let source = source.prepare(wavelength_nm)?;
        let target = target.prepare(wavelength_nm)?;
        let sample_i = source.optical_sample(incident.ray.origin)?;
        let sample_t = target.optical_sample(incident.ray.origin)?;
        if sample_i != incident.final_sample {
            return Err(P1ContactError::ForeignArrival);
        }
        let channels = contact_channels(
            self.contact.coating,
            incident_side,
            incident,
            incident_stokes,
            witness.normal,
            sample_i.phase_index(),
            sample_t.phase_index(),
            wavelength_nm,
        )?;
        let mut resolution = previous_resolution_m;
        let reflected = if let Some((power, frame)) = channels.reflected {
            let branch = self.branch(
                cx,
                incident_side,
                incident,
                power,
                frame,
                wavelength_nm,
                collar,
                &source,
                BoundarySide::Negative,
                resolution,
                &mut work,
                options,
            )?;
            resolution = resolution.max(branch.contact_resolution_m);
            Some(branch)
        } else {
            None
        };
        let transmitted = if let Some((power, frame)) = channels.transmitted {
            let branch = self.branch(
                cx,
                1 - incident_side,
                incident,
                power,
                frame,
                wavelength_nm,
                collar,
                &target,
                BoundarySide::Positive,
                previous_resolution_m,
                &mut work,
                options,
            )?;
            resolution = resolution.max(branch.contact_resolution_m);
            Some(branch)
        } else {
            None
        };
        Ok(P1CurvedContactSplit {
            incident,
            incident_stokes,
            incident_frame: channels.incident_frame,
            fresnel: channels.fresnel,
            coating_absorbed_power_fraction: channels.coating_absorbed_power_fraction,
            reflected,
            transmitted,
            attempted_steps: work.steps,
            boundary_evaluations: work.queries,
            contact_resolution_m: resolution,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn branch<F: GradientIndexField>(
        &self,
        cx: &Cx<'_>,
        side_index: usize,
        previous: GradientIndexPropagation,
        power: Stokes,
        frame: PolarizationFrame,
        wavelength_nm: f64,
        collar: ContactCollar,
        field: &F,
        selected_side: BoundarySide,
        previous_resolution_m: f64,
        work: &mut Work,
        options: P1CurvedContactOptions,
    ) -> Result<P1CurvedContactBranch<'a>, P1ContactError> {
        poll(cx)?;
        let side = &self.contact.sides[side_index];
        let departure = GradientDeparture {
            witness: collar.witness,
            side: selected_side,
        };
        let mut boundary = ContactBoundary {
            current: side.binding.temperature_source().boundary_oracle(
                side.binding.temperature_source().element_index(),
                options.event.surface_tolerance_m,
            )?,
            collar,
            departure,
        };
        let event = propagate_to_boundary_with_control(
            field,
            &mut boundary,
            Ray {
                origin: previous.ray.origin,
                dir: frame.direction(),
            },
            frame.horizontal(),
            wavelength_nm,
            options.max_outgoing_distance_m,
            Some(departure),
            work.remaining(options.event)?,
            || cx.checkpoint().is_ok(),
        )?;
        work.record(&event)?;
        let final_frame = PolarizationFrame::try_new(
            event.propagation.ray.dir,
            event.propagation.transverse_axis,
        )?;
        // Q/U/V stay in the integrator's parallel-transported frame. Rotating
        // between the nonparallel endpoint planes would invent a bulk analyzer.
        let stokes = physical_roundoff(
            scale_stokes(power, event.propagation.power_transmittance()),
            power.i,
        )?;
        Ok(P1CurvedContactBranch {
            propagation: append(previous, event.propagation)?,
            stokes,
            frame: final_frame,
            binding: side.binding,
            wavelength_nm,
            termination: event.termination,
            contact_resolution_m: previous_resolution_m.max(event.contact_resolution_m),
        })
    }
}

#[derive(Default)]
struct Work {
    steps: u32,
    queries: u32,
}

impl Work {
    fn remaining(
        &self,
        mut options: GradientEventOptions,
    ) -> Result<GradientEventOptions, P1ContactError> {
        options.integration.max_steps = options
            .integration
            .max_steps
            .checked_sub(self.steps)
            .filter(|&v| v > 0)
            .ok_or(GradientIndexError::StepBudgetExceeded)?;
        options.max_boundary_queries = options
            .max_boundary_queries
            .checked_sub(self.queries)
            .filter(|&v| v > 0)
            .ok_or(GradientIndexError::BoundaryQueryBudgetExceeded)?;
        Ok(options)
    }
    fn record(&mut self, event: &GradientBoundaryPropagation) -> Result<(), P1ContactError> {
        self.steps = self
            .steps
            .checked_add(event.propagation.accepted_steps)
            .and_then(|v| v.checked_add(event.propagation.rejected_steps))
            .ok_or(GradientIndexError::StepBudgetExceeded)?;
        self.queries = self
            .queries
            .checked_add(event.boundary_evaluations)
            .ok_or(GradientIndexError::BoundaryQueryBudgetExceeded)?;
        Ok(())
    }
    fn charge_geometry(&mut self, options: GradientEventOptions) -> Result<(), P1ContactError> {
        if self.queries >= options.max_boundary_queries {
            return Err(GradientIndexError::BoundaryQueryBudgetExceeded.into());
        }
        self.queries += 1;
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct ContactCollar {
    origin: Point3,
    witness: BoundaryWitness,
    radius_m: f64,
    slab_half_width_m: f64,
    surface_tolerance_m: f64,
}

impl ContactCollar {
    fn prepare(
        contact: &P1ContactInterface<'_>,
        incident_side: usize,
        origin: Point3,
        witness: BoundaryWitness,
        options: GradientEventOptions,
    ) -> Result<Self, P1ContactError> {
        let source = &contact.sides[incident_side];
        if witness.boundary_id != source.binding.temperature_source().element_index()
            || witness.feature_id != source.face
            || witness.side != BoundarySide::Negative
            || witness.distance_upper_m > options.surface_tolerance_m
        {
            return Err(P1ContactError::ForeignArrival);
        }
        let (_, original) = triangle_proximity(
            origin,
            source.binding.temperature_source().faces_m()[source.face],
            source.face as u32,
        )
        .map_err(|_| P1ContactError::UnresolvedFace)?;
        if original.surface_point != witness.surface_point
            || original.normal != Some(witness.normal)
            || original.side != ProximitySide::Negative
            || original.distance_upper_m > witness.distance_upper_m
        {
            return Err(P1ContactError::ForeignArrival);
        }
        let destination = &contact.sides[1 - incident_side];
        let mut clearance = f64::INFINITY;
        for (face, vertices) in destination
            .binding
            .temperature_source()
            .faces_m()
            .into_iter()
            .enumerate()
        {
            if face != destination.face {
                clearance = clearance.min(
                    triangle_proximity(witness.surface_point, vertices, face as u32)
                        .map_err(|_| P1ContactError::UnresolvedFace)?
                        .0,
                );
            }
        }
        let radius_m = ((Interval::point(clearance)
            - Interval::point(original.surface_point_error_m))
            * Interval::point(0.5))
        .lo()
        .min(
            witness
                .departure_radius_m
                .ok_or(GradientIndexError::UnsupportedBoundaryDeparture)?,
        );
        let slab_half_width_m = (Interval::point(options.surface_tolerance_m)
            + Interval::point(options.integration.position_tolerance_m)
            + Interval::point(2.0) * Interval::point(original.surface_point_error_m))
        .hi();
        if !radius_m.is_finite()
            || radius_m <= 4.0 * options.surface_tolerance_m
            || !slab_half_width_m.is_finite()
            || slab_half_width_m <= 0.0
        {
            return Err(GradientIndexError::BoundaryResolutionLimit.into());
        }
        Ok(Self {
            origin,
            witness,
            radius_m,
            slab_half_width_m,
            surface_tolerance_m: options.surface_tolerance_m,
        })
    }

    fn contains(self, position: Point3) -> bool {
        if ![position.x, position.y, position.z]
            .into_iter()
            .all(f64::is_finite)
        {
            return false;
        }
        let p = [position.x, position.y, position.z].map(Interval::point);
        let center = [
            self.witness.surface_point.x,
            self.witness.surface_point.y,
            self.witness.surface_point.z,
        ]
        .map(Interval::point);
        let d = core::array::from_fn::<_, 3, _>(|i| p[i] - center[i]);
        let a = d.map(Interval::abs);
        let radius = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt().hi();
        let n = [
            self.witness.normal.x,
            self.witness.normal.y,
            self.witness.normal.z,
        ]
        .map(Interval::point);
        let slab = (d[0] * n[0] + d[1] * n[1] + d[2] * n[2]).abs().hi();
        radius <= self.radius_m && slab <= self.slab_half_width_m
    }
}

// The public P1 element keeps its strict original domain. Only this privately
// admitted collar evaluates its one-sided affine continuation, never an average
// with the adjacent material and never a newly fitted temperature law.
struct CollarTemperature<'a> {
    element: P1ConductionElement<'a>,
    face: usize,
    collar: ContactCollar,
}

impl CollarTemperature<'_> {
    fn sample(&self, position: Point3) -> Result<TemperatureGradientSample, GradientIndexError> {
        match self.element.sample_temperature_jet(position) {
            Ok(jet) => return Ok(jet.sample),
            Err(TransientTemperatureFieldError::OutsideElement) => {}
            Err(error) => return Err(error.into_gradient_error()),
        }
        if !self.collar.contains(position) {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        let displacement = position.delta_from(self.element.positions[0]);
        let mut lambda = [0.0; 4];
        for (entry, gradient) in lambda.iter_mut().zip(self.element.gradients).skip(1) {
            *entry = gradient.dot(displacement);
        }
        lambda[0] = ((1.0 - lambda[1]) - lambda[2]) - lambda[3];
        if lambda
            .iter()
            .enumerate()
            .any(|(i, &v)| !v.is_finite() || (i != self.face && v < 0.0))
        {
            return Err(GradientIndexError::OutsideFieldDomain);
        }
        let values = self.element.values;
        let temperature_k = lambda[1].mul_add(
            values[1] - values[0],
            lambda[2].mul_add(
                values[2] - values[0],
                lambda[3].mul_add(values[3] - values[0], values[0]),
            ),
        );
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(GradientIndexError::InvalidSample);
        }
        Ok(TemperatureGradientSample {
            temperature_k,
            gradient_k_per_m: self.element.gradient,
        })
    }
}
impl TemperatureField for CollarTemperature<'_> {
    fn temperature_k(&self, position: Point3) -> Result<f64, GradientIndexError> {
        Ok(self.sample(position)?.temperature_k)
    }
}
impl TemperatureGradientField for CollarTemperature<'_> {
    fn temperature_gradient(
        &self,
        position: Point3,
    ) -> Result<TemperatureGradientSample, GradientIndexError> {
        self.sample(position)
    }
}

struct ContactBoundary {
    current: ConductionElementBoundary,
    collar: ContactCollar,
    departure: GradientDeparture,
}
impl BoundaryOracle for ContactBoundary {
    fn query(&mut self, position: Point3) -> Result<BoundaryQuery, GradientIndexError> {
        // Includes the departure face again once the native collar is cleared.
        self.current.query(position)
    }
    fn prepare_departure(
        &mut self,
        request: DepartureRequest,
    ) -> Result<DepartureNeighborhood, GradientIndexError> {
        if request.origin != self.collar.origin
            || request.contact != self.departure
            || request.surface_tolerance_m != self.collar.surface_tolerance_m
        {
            return Err(GradientIndexError::UnsupportedBoundaryDeparture);
        }
        Ok(DepartureNeighborhood::Planar {
            radius_m: self.collar.radius_m,
        })
    }
}

#[cfg(test)]
mod tests;
