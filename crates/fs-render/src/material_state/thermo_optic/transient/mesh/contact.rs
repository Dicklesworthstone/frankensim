//! Refraction across an explicitly paired contact of two isothermal elements.
//!
//! Each side retains its original P1 snapshot, reference glass and thermo-optic
//! law. Uniformity is proved from all four original nodal temperatures; a
//! gradient is never replaced by an element average. Native watertight triangle
//! hits locate the contact and the next face exactly for these straight rays.
//! There is no GRIN arrival collar, thermal extrapolation or ray-position offset.
//! The native dielectric owner supplies Fresnel Mueller matrices, Snell
//! directions and total-internal-reflection phase. Actual s/p frames and the
//! complete Stokes state survive every contact. Thermal contact resistance supplies
//! no optical absorption: bulk Beer loss remains with each original glass.
//! An explicitly bound coherent coating uses its original spectral layer laws,
//! ordered from side 0 to side 1. Its positive dissipative work is reported
//! separately from either bulk path, including losses below the R+T residual.
//!
//! Preparation stores eight one-triangle native charts. Queries allocate nothing
//! and make at most twelve triangle intersection queries. Returned branches can
//! be passed to another explicitly prepared contact on the same original
//! thermal snapshot. This is geometrical, incoherent, smooth-interface power
//! transport; a coating is a local planar scattering boundary, without a ray
//! displacement or film transit delay. It adds no coherence between distinct
//! paths, roughness, gradient-element continuation or bundle derivatives.
//! [`curved`] adds a separate bounded owner for nonuniform P1 elements.

use super::*;
use crate::charts::proximity::segment_clearance;
use crate::charts::{MeshHit, TriMesh};
use crate::dielectric::{
    DielectricCoatingStack, DielectricError, FresnelDielectric, reflect, refract,
};
use crate::gradient_index::{
    GradientIndexField, GradientIndexOptions, GradientIndexSample, propagate_with_control,
};
use crate::polarization::dielectric::fresnel_dielectric_polarized;
use crate::polarization::dielectric::transport::{BeamError, PolarizationFrame, physical_roundoff};
use crate::polarization::{DiagonalJonesError, Mueller, Stokes};

/// Curved propagation and polarized refraction through nonuniform P1 contacts.
pub mod curved;

type Binding<'a> = MaterialThermoOpticField<'a, P1ConductionElement<'a>>;

struct ContactSide<'a> {
    binding: &'a Binding<'a>,
    face: usize,
    faces: [TriMesh; 4],
}

/// Two original coincident triangular faces and their distinct optical owners.
/// Faces may share vertex IDs or use duplicated thermal-contact nodes. Their
/// geometric vertices must match exactly, and the two interiors must lie on
/// opposite sides. Pairing establishes local contact, not global mesh validity.
pub struct P1ContactInterface<'a> {
    sides: [ContactSide<'a>; 2],
    coating: Option<DielectricCoatingStack<'a>>,
}

/// Complete deterministic power split at one actual contact.
#[derive(Clone, Copy)]
pub struct P1ContactSplit<'a> {
    /// Complete path to this interface, before applying its Fresnel powers.
    pub incident: GradientIndexPropagation,
    /// Incident Stokes power, including the complete preceding path.
    pub incident_stokes: Stokes,
    /// Actual transverse reference of `incident_stokes`.
    pub incident_frame: PolarizationFrame,
    /// Native unpolarized diagnostic coefficients. Actual branch powers use
    /// the retained Stokes state, not repeated products of these averages.
    pub fresnel: FresnelDielectric,
    /// Power absorbed in the finite coating at THIS contact, relative to the
    /// original unit-intensity launch. Uses positive film work weighted by the
    /// actual incident s/p powers, not the rounded `1 - R - T` residual. Bare
    /// and lossless contacts return exact zero. Bulk attenuation, earlier
    /// contact losses and any later branch propagation are excluded.
    pub coating_absorbed_power_fraction: f64,
    /// Absent for an exactly null reflected channel.
    pub reflected: Option<P1ContactBranch<'a>>,
    /// Absent for total internal reflection or zero transmitted power.
    pub transmitted: Option<P1ContactBranch<'a>>,
}

/// One immutable outgoing branch, including every preceding path segment.
/// The private original-owner/face witness makes subsequent contact composition
/// explicit; caller-authored rays cannot masquerade as a retained arrival.
#[derive(Clone, Copy)]
pub struct P1ContactBranch<'a> {
    propagation: GradientIndexPropagation,
    stokes: Stokes,
    frame: PolarizationFrame,
    binding: &'a Binding<'a>,
    wavelength_nm: f64,
    arrival: Option<(usize, MeshHit)>,
}

impl P1ContactBranch<'_> {
    pub const fn propagation(&self) -> &GradientIndexPropagation {
        &self.propagation
    }
    /// Power relative to the original incident beam, including all interface
    /// fractions and Beer attenuation. No radiance eta-squared factor belongs
    /// in this normal-power fraction.
    pub const fn power_fraction(&self) -> f64 {
        self.stokes.i
    }
    /// Full Stokes power relative to unit input intensity, expressed in `frame()`.
    pub const fn stokes(&self) -> Stokes {
        self.stokes
    }
    pub const fn frame(&self) -> PolarizationFrame {
        self.frame
    }
    pub fn element_index(&self) -> usize {
        self.binding.temperature_source().element_index()
    }
    /// Original local face at termination, or None for the distance limit.
    pub fn arrival_face(&self) -> Option<usize> {
        self.arrival.map(|(face, _)| face)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum P1ContactError {
    Mesh(P1MeshError),
    Dielectric(DielectricError),
    Polarization(BeamError),
    InvalidPair,
    NonuniformElement {
        element: usize,
    },
    InvalidLaunch,
    ContactNotFirstFace {
        local_face: usize,
    },
    /// The incident arc budget ended before the selected contact was reached.
    IncidentDistanceLimit,
    ForeignArrival,
    UnresolvedFace,
}

impl From<GradientIndexError> for P1ContactError {
    fn from(error: GradientIndexError) -> Self {
        Self::Mesh(error.into())
    }
}
impl From<ThermoOpticFieldError> for P1ContactError {
    fn from(error: ThermoOpticFieldError) -> Self {
        Self::Mesh(error.into())
    }
}
impl From<DielectricError> for P1ContactError {
    fn from(error: DielectricError) -> Self {
        Self::Dielectric(error)
    }
}
impl From<BeamError> for P1ContactError {
    fn from(error: BeamError) -> Self {
        Self::Polarization(error)
    }
}
impl From<DiagonalJonesError> for P1ContactError {
    fn from(error: DiagonalJonesError) -> Self {
        Self::Polarization(error.into())
    }
}
impl core::fmt::Display for P1ContactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "P1 contact refraction: {self:?}")
    }
}
impl core::error::Error for P1ContactError {}

impl<'a> P1ContactInterface<'a> {
    pub fn try_new(
        cx: &Cx<'_>,
        bindings: [&'a Binding<'a>; 2],
        local_faces: [usize; 2],
    ) -> Result<Self, P1ContactError> {
        Self::prepare(cx, bindings, local_faces, true)
    }

    fn prepare(
        cx: &Cx<'_>,
        bindings: [&'a Binding<'a>; 2],
        local_faces: [usize; 2],
        require_uniform: bool,
    ) -> Result<Self, P1ContactError> {
        poll(cx)?;
        let elements = bindings.map(|b| b.temperature_source());
        if local_faces.iter().any(|&f| f >= 4)
            || !same_snapshot(bindings[0], bindings[1])
            || elements[0].element == elements[1].element
        {
            return Err(P1ContactError::InvalidPair);
        }
        for element in elements {
            if require_uniform && element.values.iter().any(|&t| t != element.values[0]) {
                return Err(P1ContactError::NonuniformElement {
                    element: element.element,
                });
            }
        }
        let faces = core::array::from_fn::<_, 2, _>(|i| elements[i].faces_m());
        let a = faces[0][local_faces[0]];
        let b = faces[1][local_faces[1]];
        if !a.iter().all(|vertex| b.contains(vertex))
            || triangle_proximity(elements[1].positions[local_faces[1]], a, 0)
                .map_err(|_| P1ContactError::InvalidPair)?
                .1
                .side
                != ProximitySide::Positive
            || triangle_proximity(elements[0].positions[local_faces[0]], b, 0)
                .map_err(|_| P1ContactError::InvalidPair)?
                .1
                .side
                != ProximitySide::Positive
        {
            return Err(P1ContactError::InvalidPair);
        }
        let sides = core::array::from_fn(|i| ContactSide {
            binding: bindings[i],
            face: local_faces[i],
            faces: faces[i].map(|vertices| TriMesh::new(vertices.to_vec(), vec![[0, 1, 2]])),
        });
        Ok(Self {
            sides,
            coating: None,
        })
    }

    /// Bind original coherent films in physical order from side 0 to side 1.
    /// Incidence from side 1 reverses their traversal automatically. The stack
    /// retains its declared wavelength band, including zero-thickness films;
    /// each adjacent index still comes from that side's accepted thermal state.
    ///
    /// Films are modeled as a local planar boundary: their coherent phases
    /// affect R/T and polarization, but no geometric thickness, group delay or
    /// thermal-state update is added to the bulk path. Returned absorbed power
    /// is an optical result, not an applied conduction source.
    #[must_use]
    pub fn with_coating_stack(mut self, coating: DielectricCoatingStack<'a>) -> Self {
        self.coating = Some(coating);
        self
    }

    /// Original borrowed film laws, with side-0-to-side-1 ordering.
    pub const fn coating_stack(&self) -> Option<DielectricCoatingStack<'a>> {
        self.coating
    }

    /// Trace an initially unpolarized interior launch to this contact, then propagate both nonzero
    /// power branches to their next original face or the supplied OUTGOING arc
    /// limit. Each side's original thermo-optic law is evaluated at the actual
    /// vacuum wavelength and its own accepted uniform temperature.
    #[allow(clippy::too_many_arguments)]
    pub fn split(
        &self,
        cx: &Cx<'_>,
        incident_side: usize,
        ray: Ray,
        transverse_axis: Vec3,
        wavelength_nm: f64,
        max_outgoing_distance_m: f64,
    ) -> Result<P1ContactSplit<'a>, P1ContactError> {
        self.split_polarized(
            cx,
            incident_side,
            ray,
            transverse_axis,
            Stokes::unpolarized(1.0),
            wavelength_nm,
            max_outgoing_distance_m,
        )
    }

    /// The same exact path for an explicitly polarized unit-intensity launch.
    /// `(Q,U,V)` refer to the input ray's transverse axis. Requiring I=1 keeps
    /// every returned Stokes component and power fraction relative to the
    /// original launch; the caller scales by physical incident spectral power.
    #[allow(clippy::too_many_arguments)]
    pub fn split_polarized(
        &self,
        cx: &Cx<'_>,
        incident_side: usize,
        ray: Ray,
        transverse_axis: Vec3,
        polarization: Stokes,
        wavelength_nm: f64,
        max_outgoing_distance_m: f64,
    ) -> Result<P1ContactSplit<'a>, P1ContactError> {
        poll(cx)?;
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
            .sides
            .get(incident_side)
            .ok_or(P1ContactError::InvalidLaunch)?;
        let barycentric = side
            .binding
            .temperature_source()
            .barycentric(ray.origin)
            .map_err(|_| P1ContactError::InvalidLaunch)?;
        if barycentric.iter().any(|&b| b <= 0.0) {
            return Err(P1ContactError::InvalidLaunch);
        }
        let sample = optical_sample(side.binding, wavelength_nm)?;
        let launch = advance(cx, sample, ray, transverse_axis, wavelength_nm, 0.0)?;
        let (face, hit) = side.next_face(cx, &launch.ray, None, f64::INFINITY)?;
        if face != side.face {
            return Err(P1ContactError::ContactNotFirstFace { local_face: face });
        }
        let mut incident = advance(
            cx,
            sample,
            launch.ray,
            launch.transverse_axis,
            wavelength_nm,
            hit.hit.t,
        )?;
        // This is the same analytic straight-ray endpoint, retained in the
        // native triangle owner's arithmetic; no unintegrated offset is added.
        incident.ray.origin = hit.hit.point;
        self.split_hit(
            cx,
            incident_side,
            incident,
            scale_stokes(polarization, incident.power_transmittance()),
            hit,
            wavelength_nm,
            max_outgoing_distance_m,
        )
    }

    /// Continue a caller-selected branch that reached this interface's original
    /// face. Snapshot, material ownership, wavelength and accumulated integrals
    /// remain those of the original path. No geometry search, epsilon offset,
    /// resampling of the previous event or source-state refresh is performed.
    pub fn split_arrival(
        &self,
        cx: &Cx<'_>,
        arrival: &P1ContactBranch<'a>,
        max_outgoing_distance_m: f64,
    ) -> Result<P1ContactSplit<'a>, P1ContactError> {
        poll(cx)?;
        let (face, hit) = arrival.arrival.ok_or(P1ContactError::ForeignArrival)?;
        let incident_side = self
            .sides
            .iter()
            .position(|s| s.face == face && same_region(s.binding, arrival.binding))
            .ok_or(P1ContactError::ForeignArrival)?;
        self.split_hit(
            cx,
            incident_side,
            arrival.propagation,
            arrival.stokes,
            hit,
            arrival.wavelength_nm,
            max_outgoing_distance_m,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn split_hit(
        &self,
        cx: &Cx<'_>,
        incident_side: usize,
        incident: GradientIndexPropagation,
        incident_stokes: Stokes,
        hit: MeshHit,
        wavelength_nm: f64,
        max_distance_m: f64,
    ) -> Result<P1ContactSplit<'a>, P1ContactError> {
        if !max_distance_m.is_finite() || max_distance_m < 0.0 {
            return Err(GradientIndexError::InvalidOptions.into());
        }
        let sample_i = optical_sample(self.sides[incident_side].binding, wavelength_nm)?;
        let sample_t = optical_sample(self.sides[1 - incident_side].binding, wavelength_nm)?;
        if incident.final_sample != sample_i {
            return Err(P1ContactError::ForeignArrival);
        }
        let outward = hit.hit.normal.ok_or(P1ContactError::UnresolvedFace)?;
        let channels = contact_channels(
            self.coating,
            incident_side,
            incident,
            incident_stokes,
            outward,
            sample_i.phase_index(),
            sample_t.phase_index(),
            wavelength_nm,
        )?;
        let reflected = if let Some((power, frame)) = channels.reflected {
            Some(self.branch(
                cx,
                incident_side,
                incident,
                power,
                frame,
                sample_i,
                wavelength_nm,
                max_distance_m,
            )?)
        } else {
            None
        };
        let transmitted = if let Some((power, frame)) = channels.transmitted {
            Some(self.branch(
                cx,
                1 - incident_side,
                incident,
                power,
                frame,
                sample_t,
                wavelength_nm,
                max_distance_m,
            )?)
        } else {
            None
        };
        Ok(P1ContactSplit {
            incident,
            incident_stokes,
            incident_frame: channels.incident_frame,
            fresnel: channels.fresnel,
            coating_absorbed_power_fraction: channels.coating_absorbed_power_fraction,
            reflected,
            transmitted,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn branch(
        &self,
        cx: &Cx<'_>,
        side: usize,
        previous: GradientIndexPropagation,
        power: Stokes,
        frame: PolarizationFrame,
        sample: GradientIndexSample,
        wavelength_nm: f64,
        max_distance_m: f64,
    ) -> Result<P1ContactBranch<'a>, P1ContactError> {
        let ray = Ray {
            origin: previous.ray.origin,
            dir: frame.direction(),
        };
        let launch = advance(cx, sample, ray, frame.horizontal(), wavelength_nm, 0.0)?;
        let side = &self.sides[side];
        // A straight ray selected into this convex element cannot return to
        // its departure plane. All other original faces remain in the search.
        let (face, hit) = side.next_face(cx, &launch.ray, Some(side.face), max_distance_m)?;
        let reached = hit.hit.t <= max_distance_m;
        let mut segment = advance(
            cx,
            sample,
            launch.ray,
            launch.transverse_axis,
            wavelength_nm,
            hit.hit.t.min(max_distance_m),
        )?;
        if reached {
            segment.ray.origin = hit.hit.point;
        }
        let final_frame = PolarizationFrame::try_new(segment.ray.dir, segment.transverse_axis)?;
        let stokes = physical_roundoff(
            frame.express_in(
                scale_stokes(power, segment.power_transmittance()),
                final_frame,
            ),
            power.i,
        )?;
        let propagation = append(previous, segment)?;
        Ok(P1ContactBranch {
            propagation,
            stokes,
            frame: final_frame,
            binding: side.binding,
            wavelength_nm,
            arrival: reached.then_some((face, hit)),
        })
    }
}

// Both exact straight and numerical curved contacts use this same optical
// boundary response. Geometry and bulk propagation remain with their owners.
struct ContactChannels {
    incident_frame: PolarizationFrame,
    fresnel: FresnelDielectric,
    coating_absorbed_power_fraction: f64,
    reflected: Option<(Stokes, PolarizationFrame)>,
    transmitted: Option<(Stokes, PolarizationFrame)>,
}

#[allow(clippy::too_many_arguments)]
fn contact_channels(
    coating: Option<DielectricCoatingStack<'_>>,
    incident_side: usize,
    incident: GradientIndexPropagation,
    incident_stokes: Stokes,
    outward: Vec3,
    ni: f64,
    nt: f64,
    wavelength_nm: f64,
) -> Result<ContactChannels, P1ContactError> {
    let cosine = incident.ray.dir.dot(outward);
    if !cosine.is_finite() || cosine <= 0.0 || cosine > 1.0 + 2.0e-10 {
        return Err(P1ContactError::InvalidLaunch);
    }
    // One native solve supplies branch matrices and positive film work.
    let balance = coating
        .map(|stack| {
            stack.layer_power_balance(wavelength_nm, incident_side == 1, cosine.min(1.0), ni, nt)
        })
        .transpose()?;
    let response = match &balance {
        Some(balance) => balance.fresnel,
        None => fresnel_dielectric_polarized(cosine.min(1.0), ni, nt)?,
    };
    let fresnel = response.scalar;
    let normal = outward.scale(-1.0);
    let wo = incident.ray.dir.scale(-1.0);
    let incident_frame = PolarizationFrame::try_new(incident.ray.dir, incident.transverse_axis)?;
    let perpendicular = cross(incident.ray.dir, normal);
    let s = if perpendicular == Vec3::new(0.0, 0.0, 0.0) {
        incident_frame.horizontal()
    } else {
        perpendicular
    };
    let interface_frame = PolarizationFrame::try_new(incident.ray.dir, s)?;
    let input = physical_roundoff(
        incident_frame.express_in(incident_stokes, interface_frame),
        incident_stokes.i,
    )?;
    let reflected_power = apply_power(response.reflection, input)?;
    let transmitted_power = apply_power(response.transmission, input)?;
    let coating_absorbed_power_fraction = balance.map_or(0.0, |balance| {
        let s_power = (0.5 * input.i + 0.5 * input.q).max(0.0);
        let p_power = (0.5 * input.i - 0.5 * input.q).max(0.0);
        s_power * balance.absorbed[0] + p_power * balance.absorbed[1]
    });
    let reflected = if reflected_power.i > 0.0 {
        Some((
            reflected_power,
            PolarizationFrame::try_new(reflect(wo, normal), interface_frame.horizontal())?,
        ))
    } else {
        None
    };
    let transmitted = if transmitted_power.i > 0.0 {
        Some((
            transmitted_power,
            PolarizationFrame::try_new(
                refract(wo, normal, ni, nt, fresnel.transmitted_cosine)?,
                interface_frame.horizontal(),
            )?,
        ))
    } else {
        None
    };
    Ok(ContactChannels {
        incident_frame,
        fresnel,
        coating_absorbed_power_fraction,
        reflected,
        transmitted,
    })
}

impl ContactSide<'_> {
    fn next_face(
        &self,
        cx: &Cx<'_>,
        ray: &Ray,
        skip: Option<usize>,
        limit_m: f64,
    ) -> Result<(usize, MeshHit), P1ContactError> {
        let mut best: Option<(usize, MeshHit)> = None;
        for (face, mesh) in self.faces.iter().enumerate() {
            if skip == Some(face) {
                continue;
            }
            let hit = mesh
                .intersect_surface_with_cx(cx, ray)
                .map_err(|_| GradientIndexError::Cancelled)?;
            if let Some(hit) = hit {
                if hit.hit.t <= limit_m
                    && best.is_some_and(|(_, previous)| previous.hit.t == hit.hit.t)
                {
                    return Err(P1ContactError::UnresolvedFace);
                }
                if best.is_none_or(|(_, previous)| hit.hit.t < previous.hit.t) {
                    best = Some((face, hit));
                }
            }
        }
        let (face, mut hit) = best.ok_or(P1ContactError::UnresolvedFace)?;
        if hit.hit.t > limit_m {
            return Ok((face, hit));
        }
        let vertices = self.binding.temperature_source().faces_m()[face];
        let (_, witness) = triangle_proximity(hit.hit.point, vertices, face as u32)
            .map_err(|_| P1ContactError::UnresolvedFace)?;
        for edge in 0..3 {
            if segment_clearance(hit.hit.point, vertices[edge], vertices[(edge + 1) % 3])
                .map_err(|_| P1ContactError::UnresolvedFace)?
                <= witness.surface_point_error_m
            {
                return Err(P1ContactError::UnresolvedFace);
            }
        }
        hit.hit.normal = witness.normal;
        Ok((face, hit))
    }
}

fn poll(cx: &Cx<'_>) -> Result<(), P1ContactError> {
    cx.checkpoint()
        .map_err(|_| GradientIndexError::Cancelled.into())
}
fn same_snapshot(a: &Binding<'_>, b: &Binding<'_>) -> bool {
    let a_source = a.temperature_source().snapshot();
    let b_source = b.temperature_source().snapshot();
    core::ptr::eq(a_source.mesh(), b_source.mesh())
        && core::ptr::eq(
            a_source.nodal_temperature_k(),
            b_source.nodal_temperature_k(),
        )
        && a.snapshot_identity() == b.snapshot_identity()
}
fn same_region(a: &Binding<'_>, b: &Binding<'_>) -> bool {
    same_snapshot(a, b)
        && a.temperature_source().element == b.temperature_source().element
        && core::ptr::eq(a.reference_glass(), b.reference_glass())
        && core::ptr::eq(a.law(), b.law())
        && match (a.absorption_source(), b.absorption_source()) {
            (None, None) => true,
            (Some(a), Some(b)) => core::ptr::eq(a, b),
            _ => false,
        }
}
fn optical_sample(
    binding: &Binding<'_>,
    wavelength_nm: f64,
) -> Result<GradientIndexSample, P1ContactError> {
    Ok(binding.optical_sample(binding.temperature_source().positions[0], wavelength_nm)?)
}

// A private constant coefficient law used only after native convex-cell
// geometry bounds the exact straight segment; it is not a public unbounded
// extension of the original temperature provider.
struct UniformSample(GradientIndexSample);
impl GradientIndexField for UniformSample {
    fn sample(&self, _: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
        Ok(self.0)
    }
    fn is_spatially_uniform(&self) -> bool {
        true
    }
}
fn advance(
    cx: &Cx<'_>,
    sample: GradientIndexSample,
    ray: Ray,
    axis: Vec3,
    wavelength_nm: f64,
    distance_m: f64,
) -> Result<GradientIndexPropagation, P1ContactError> {
    Ok(propagate_with_control(
        &UniformSample(sample),
        ray,
        axis,
        wavelength_nm,
        distance_m,
        GradientIndexOptions::default(),
        || cx.checkpoint().is_ok(),
    )?)
}
fn scale_stokes(input: Stokes, factor: f64) -> Stokes {
    Stokes {
        i: input.i * factor,
        q: input.q * factor,
        u: input.u * factor,
        v: input.v * factor,
    }
}
fn apply_power(matrix: Mueller, input: Stokes) -> Result<Stokes, P1ContactError> {
    let magnitude = matrix
        .0
        .iter()
        .flatten()
        .fold(0.0_f64, |m, v| m.max(v.abs()));
    // Every native bare/coated interface is diagonal Jones in this s/p frame.
    // Preserve its weak polarization channel near Brewster incidence instead
    // of subtracting two rounded, nearly equal Mueller intensity entries.
    Ok(physical_roundoff(
        matrix.apply_diagonal_jones(input)?,
        input.i * magnitude,
    )?)
}
fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}
fn append(
    previous: GradientIndexPropagation,
    mut next: GradientIndexPropagation,
) -> Result<GradientIndexPropagation, P1ContactError> {
    let add = |a, b| (Dd::from_f64(a) + Dd::from_f64(b)).to_f64();
    next.distance_m = add(previous.distance_m, next.distance_m);
    next.optical_path_length_m = add(previous.optical_path_length_m, next.optical_path_length_m);
    next.optical_depth = add(previous.optical_depth, next.optical_depth);
    next.group_optical_path_length_m = previous
        .group_optical_path_length_m
        .zip(next.group_optical_path_length_m)
        .map(|(a, b)| add(a, b));
    if !next.distance_m.is_finite()
        || !next.optical_path_length_m.is_finite()
        || !next.optical_depth.is_finite()
        || next
            .group_optical_path_length_m
            .is_some_and(|v| !v.is_finite())
    {
        return Err(GradientIndexError::NumericalRange.into());
    }
    next.initial_sample = previous.initial_sample;
    next.accepted_steps = previous
        .accepted_steps
        .checked_add(next.accepted_steps)
        .ok_or(GradientIndexError::StepBudgetExceeded)?;
    next.rejected_steps = previous
        .rejected_steps
        .checked_add(next.rejected_steps)
        .ok_or(GradientIndexError::StepBudgetExceeded)?;
    next.field_evaluations = previous
        .field_evaluations
        .checked_add(next.field_evaluations)
        .ok_or(GradientIndexError::StepBudgetExceeded)?;
    next.max_normalized_error_estimate = previous
        .max_normalized_error_estimate
        .max(next.max_normalized_error_estimate);
    Ok(next)
}

#[cfg(test)]
mod tests;
