//! Source-backed contact resistance at explicitly declared trace coordinates.
//!
//! Only the declared temperature coordinates change. Pressure, finish, medium,
//! source card and selected claim remain fixed. The two bulk traces stay distinct;
//! an arithmetic mean is used only by an explicitly selected mean-temperature law.
//! Frozen snapshots are NOT complete nonlinear-contact tangents or adjoints.

pub mod capacity;
/// Steady source-backed contact coupling with physical watt residuals.
pub mod steady;
/// Primary-enthalpy contact evolution, including latent phase states.
pub mod enthalpy;
/// Complete source-backed contact temperature Jacobian and its transpose.
pub mod tangent;

use fs_exec::Cx;
use fs_matdb::{InterfaceSystemCard, QueryPoint, SelectionPolicy};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

use super::{InterfaceFacePair, InterfaceResistance, InterfaceSurface, NonmatchingOptions,
    NonmatchingSurface, ThermalInterfaces};
use crate::{ConductionError, ConductionMesh, ThermalBoundary};

/// Thermal coordinates owned by a contact source law. Matching contacts use
/// face means; nonmatching contacts use both traces at each common quadrature point.
///
/// Face-pair A/B order is the source card's declared surface order. A reversal
/// of the physical interface must reverse both the face assignment and the
/// corresponding source coordinates; the solver never sorts the two sides.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TemperatureContactCoordinates {
    /// One source axis follows (trace A + trace B)/2 at the declared samples.
    Mean(String),
    /// Separate source axes follow the independent traces in A/B order.
    TwoSided([String; 2]),
}

#[derive(Debug, Clone)]
enum ContactSurface {
    Matching(InterfaceSurface),
    Nonmatching(NonmatchingSurface),
}

enum ContactFaces {
    Matching(Vec<InterfaceFacePair>),
    Nonmatching { side_a:Vec<usize>, side_b:Vec<usize>, options:NonmatchingOptions },
}

/// A named surface and one source law, selected once at a reference point.
#[derive(Debug, Clone)]
pub struct TemperatureContact<'a> {
    card: &'a InterfaceSystemCard,
    reference_point: QueryPoint,
    coordinates: TemperatureContactCoordinates,
    reference: ContactSurface,
    resistance: InterfaceResistance,
}

impl<'a> TemperatureContact<'a> {
    /// Select a source claim at the explicit reference temperature, then pin it.
    /// Later evaluations cannot switch claims when the original support ends.
    /// Typed axes must describe absolute temperature (or dimensional kelvin);
    /// only the historical untyped axis `T` is also accepted.
    pub fn new(
        name: impl Into<String>,
        face_pairs: Vec<InterfaceFacePair>,
        card: &'a InterfaceSystemCard,
        point: &QueryPoint,
        temperature_axis: &str,
        selection: SelectionPolicy,
    ) -> Result<Self, ConductionError> {
        Self::bind(name.into(), ContactFaces::Matching(face_pairs), card, point,
            TemperatureContactCoordinates::Mean(temperature_axis.into()), selection)
    }

    /// Bind R(T_A,T_B) to the independently resolved P1 trace means.
    ///
    /// The two named axes must be distinct, present in the original query,
    /// and describe positive absolute temperatures. A source may depend on
    /// either axis or both; its original interpolation cells and validity
    /// domain still govern every query. Pressure, medium, finish and history
    /// remain fixed at the explicitly supplied source conditions.
    pub fn new_two_sided(
        name: impl Into<String>,
        face_pairs: Vec<InterfaceFacePair>,
        card: &'a InterfaceSystemCard,
        point: &QueryPoint,
        temperature_axes: [&str; 2],
        selection: SelectionPolicy,
    ) -> Result<Self, ConductionError> {
        Self::bind(name.into(), ContactFaces::Matching(face_pairs), card, point,
            TemperatureContactCoordinates::TwoSided(temperature_axes.map(String::from)), selection)
    }

    /// Bind a local source law to two independently triangulated planar traces.
    /// Each retained three-point, degree-two overlap rule queries the original
    /// source at both P1 temperatures at the SAME physical point. The nonlinear
    /// reciprocal law is quadrature-approximated, not integrated exactly or
    /// replaced by a face-mean source. Geometry/work policies remain explicit.
    #[allow(clippy::too_many_arguments)]
    pub fn new_nonmatching(
        name:impl Into<String>, side_a:Vec<usize>, side_b:Vec<usize>,
        card:&'a InterfaceSystemCard, point:&QueryPoint,
        coordinates:TemperatureContactCoordinates, selection:SelectionPolicy,
        options:NonmatchingOptions,
    )->Result<Self,ConductionError> {
        Self::bind(name.into(),ContactFaces::Nonmatching{side_a,side_b,options},
            card,point,coordinates,selection)
    }

    fn bind(
        name: String,
        faces: ContactFaces,
        card: &'a InterfaceSystemCard,
        point: &QueryPoint,
        coordinates: TemperatureContactCoordinates,
        selection: SelectionPolicy,
    ) -> Result<Self, ConductionError> {
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature, ValueForm::Static,
        ));
        let axes = match &coordinates {
            TemperatureContactCoordinates::Mean(axis) => std::slice::from_ref(axis),
            TemperatureContactCoordinates::TwoSided(axes) => axes.as_slice(),
        };
        for (i, axis) in axes.iter().enumerate() {
            if axis.trim().is_empty() || axes[..i].contains(axis) {
                return Err(invalid(&name, "contact temperature coordinates must be nonblank and distinct"));
            }
            let descriptor_ok = match point.axis_quantities().get(axis) {
                Some(q) => *q == absolute || *q == QuantitySpec::dimensional(Temperature::DIMS),
                None => axis == "T",
            };
            if !descriptor_ok || !point.axes().get(axis)
                .is_some_and(|t| t.is_finite() && *t > 0.0)
            {
                return Err(invalid(&name, "each contact temperature coordinate needs a positive absolute-temperature value"));
            }
        }
        let resistance = InterfaceResistance::from_card(&name, card, point, selection)?;
        let reference = match faces {
            ContactFaces::Matching(pairs)=>ContactSurface::Matching(InterfaceSurface::new(name,pairs,resistance.clone())?),
            ContactFaces::Nonmatching{side_a,side_b,options}=>ContactSurface::Nonmatching(
                NonmatchingSurface::new(name,side_a,side_b,resistance.clone(),options)?),
        };
        Ok(Self {
            card, reference_point: point.clone(), coordinates,
            reference, resistance,
        })
    }

    /// Stable interface name, retaining the ordered card's surface convention.
    #[must_use]
    pub fn name(&self) -> &str {
        match &self.reference {ContactSurface::Matching(surface)=>surface.name(),
            ContactSurface::Nonmatching(surface)=>&surface.name}
    }
    /// Original query including fixed pressure and all other source conditions.
    #[must_use]
    pub fn reference_point(&self) -> &QueryPoint { &self.reference_point }
    /// Explicit source-coordinate convention, including ordered axis names.
    #[must_use]
    pub const fn coordinates(&self) -> &TemperatureContactCoordinates { &self.coordinates }
    /// Reference source answer. Its selected claim stays pinned at every trial.
    #[must_use]
    pub fn reference_resistance(&self) -> &InterfaceResistance { &self.resistance }

    fn query_point_at(&self, temperatures_k: [f64; 2]) -> Result<QueryPoint, ConductionError> {
        let [a, b] = [positive(temperatures_k[0])?, positive(temperatures_k[1])?];
        let mut point = self.reference_point.clone();
        match &self.coordinates {
            TemperatureContactCoordinates::Mean(axis) => {
                point = point.with(axis, positive(0.5 * a + 0.5 * b)?)
                    .map_err(|error| invalid(self.name(), format!("temperature query refused: {error}")))?;
            }
            TemperatureContactCoordinates::TwoSided(axes) => {
                for (axis, temperature) in axes.iter().zip([a, b]) {
                    point = point.with(axis, temperature)
                        .map_err(|error| invalid(self.name(), format!("temperature query refused: {error}")))?;
                }
            }
        }
        Ok(point)
    }

    fn at_traces(&self, temperatures_k: [f64; 2]) -> Result<InterfaceResistance, ConductionError> {
        let point = self.query_point_at(temperatures_k)?;
        InterfaceResistance::from_card_pinned(
            self.name(), self.card, &point, self.resistance.receipt().selected,
        )
    }
}

/// Actual local state and source answer for one explicitly paired P1 triangle.
#[derive(Debug, Clone, PartialEq)]
pub struct ContactFaceState {
    /// Named physical interface.
    pub interface: String,
    /// Original oriented face slots; no faces are paired by nearest distance.
    pub pair: InterfaceFacePair,
    /// P1 area mean of the independent A temperature trace, K.
    pub side_a_temperature_k: f64,
    /// P1 area mean of the independent B temperature trace, K.
    pub side_b_temperature_k: f64,
    /// Arithmetic mean of the two trace means, K. This is the query coordinate
    /// only for a mean-temperature law; two-sided laws use the separate traces.
    pub interface_temperature_k: f64,
    /// Fresh positive source answer and its exact local query receipt.
    pub resistance: InterfaceResistance,
}

/// Contact operator frozen at a stated field, with the complete local source data.
#[derive(Debug, Clone)]
pub struct TemperatureContactEvaluation {
    interfaces: ThermalInterfaces,
    states: Vec<ContactFaceState>,
    quadrature_states:Vec<ContactQuadratureState>,
}

/// One physical common-refinement quadrature query, with its original A/B roles.
#[derive(Debug,Clone,PartialEq)]
pub struct ContactQuadratureState {
    /// Named physical interface.
    pub interface:String,
    /// Ordered original boundary-face slots containing this common point.
    pub pair:InterfaceFacePair,
    /// Canonical triangle/point index within this named interface.
    pub quadrature_index:usize,
    /// Positive physical area weight of this integration point, m².
    pub weight_m2:f64,
    /// Interpolated P1 surface temperature on the original A trace, K.
    pub side_a_temperature_k:f64,
    /// Interpolated P1 surface temperature on the original B trace, K.
    pub side_b_temperature_k:f64,
    /// Original selected source evaluated at these local coordinates.
    pub resistance:InterfaceResistance,
}

impl TemperatureContactEvaluation {
    /// Frozen operator for a primal trial or residual evaluation only.
    /// Its fixed-resistance adjoint omits dR/dT and is NOT a total derivative.
    #[must_use]
    pub fn frozen_interfaces(&self) -> &ThermalInterfaces { &self.interfaces }
    /// Source-backed state for each temperature-dependent face pair.
    #[must_use]
    pub fn states(&self) -> &[ContactFaceState] { &self.states }
    /// Pointwise nonmatching source states, including exactly coincident pieces.
    /// These are actual source queries; a surface aggregate uses only a
    /// representative receipt and does not invent an averaged material claim.
    #[must_use]
    pub fn quadrature_states(&self)->&[ContactQuadratureState] {&self.quadrature_states}
}

/// Complete contact set, including any unchanged constant interfaces.
/// Geometry is checked once. Frozen trials reuse retained common-refinement
/// stencils, original cards and ordered traces; they never re-intersect meshes.
#[derive(Debug)]
pub struct TemperatureContacts<'m, 'c> {
    mesh: &'m ConductionMesh,
    fixed: Vec<InterfaceSurface>,
    variable: Vec<TemperatureContact<'c>>,
    variable_pairs: usize,
    reference:ThermalInterfaces,
    nonmatching_triangles:usize,
    variable_quadrature_points:usize,
}

impl<'m, 'c> TemperatureContacts<'m, 'c> {
    /// Bind all contact faces, with an explicit total face-pair work/storage cap.
    /// At most 64 named interfaces are admitted. These count limits do not
    /// measure source-card storage, allocator overhead or peak process memory.
    pub fn new(
        cx: &Cx<'_>,
        mesh: &'m ConductionMesh,
        boundary: &ThermalBoundary,
        fixed: Vec<InterfaceSurface>,
        variable: Vec<TemperatureContact<'c>>,
        max_face_pairs: usize,
    ) -> Result<Self, ConductionError> {
        if variable.iter().any(|law|matches!(law.reference,ContactSurface::Nonmatching(_))) {
            return Err(invalid("<contact-set>","nonmatching laws require an explicit common-quadrature allowance"));
        }
        Self::with_nonmatching(cx,mesh,boundary,fixed,Vec::new(),variable,max_face_pairs,0)
    }

    /// Bind matching and nonmatching laws with independent trace-work limits.
    /// max_face_pairs bounds matching pairs plus retained overlap triangles;
    /// max_quadrature_points bounds all nonmatching samples, including fixed
    /// patches. Each planar surface also retains its own geometry-work policy.
    /// Exact pieces of variable nonmatching laws stay on the pointwise rule;
    /// they are admitted by the original matching owner but assembled once.
    #[allow(clippy::too_many_arguments)]
    pub fn with_nonmatching(
        cx:&Cx<'_>,mesh:&'m ConductionMesh,boundary:&ThermalBoundary,
        fixed:Vec<InterfaceSurface>,mut fixed_nonmatching:Vec<NonmatchingSurface>,
        mut variable:Vec<TemperatureContact<'c>>,max_face_pairs:usize,max_quadrature_points:usize,
    )->Result<Self,ConductionError> {
        poll(cx, 0)?;
        if variable.is_empty() || fixed.len().checked_add(variable.len())
            .and_then(|n|n.checked_add(fixed_nonmatching.len())).is_none_or(|n| n > 64) {
            return Err(invalid("<contact-set>", "one to 64 interfaces and at least one temperature law are required"));
        }
        let mut count = 0usize;
        for surface in fixed.iter().chain(variable.iter().filter_map(|law|match &law.reference {
            ContactSurface::Matching(surface)=>Some(surface),ContactSurface::Nonmatching(_)=>None})) {
            poll(cx, count)?;
            count = count.checked_add(surface.face_pairs().len())
                .ok_or_else(|| invalid("<contact-set>", "contact face count overflows"))?;
            if count > max_face_pairs {
                return Err(invalid(surface.name(), "contact faces exceed the declared pair allowance"));
            }
        }
        variable.sort_by(|a, b| a.name().cmp(b.name()));
        let variable_pairs = variable.iter().filter_map(|law|match &law.reference {
            ContactSurface::Matching(surface)=>Some(surface.face_pairs().len()),ContactSurface::Nonmatching(_)=>None}).sum();
        let mut all = fixed.clone();
        for law in &variable {match &law.reference {
            ContactSurface::Matching(surface)=>all.push(surface.clone()),
            ContactSurface::Nonmatching(surface)=>fixed_nonmatching.push(surface.clone()),
        }}
        // The original owner rejects missing, duplicated, shared-node,
        // noncoincident and externally occupied contact faces.
        let reference=if fixed_nonmatching.is_empty() {ThermalInterfaces::new(mesh,boundary,all)?}
            else {ThermalInterfaces::with_nonmatching(cx,mesh,boundary,all,fixed_nonmatching)?};
        let nonmatching_triangles=reference.nonmatching.iter().try_fold(0usize,|n,s|
            n.checked_add(s.triangle_count())).ok_or_else(||invalid("<contact-set>","overlap-triangle count overflows"))?;
        if count.checked_add(nonmatching_triangles).is_none_or(|n|n>max_face_pairs)
            || nonmatching_triangles.checked_mul(3).is_none_or(|n|n>max_quadrature_points
                || n.checked_mul(std::mem::size_of::<ContactQuadratureState>()).is_none_or(|bytes|bytes>isize::MAX as usize)) {
            return Err(invalid("<contact-set>","contact triangles or quadrature states exceed their explicit allowances"));
        }
        let variable_quadrature_points=reference.nonmatching.iter().filter(|surface|
            variable.iter().any(|law|law.name()==surface.name)).map(|surface|surface.triangle_count()*3).sum();
        poll(cx, count)?;
        Ok(Self { mesh, fixed, variable, variable_pairs, reference, nonmatching_triangles,variable_quadrature_points })
    }

    /// Query every declared local thermal coordinate and build one complete snapshot.
    /// All nonthermal coordinates and selected source claims stay unchanged.
    /// A single failed local query refuses the whole snapshot.
    pub fn evaluate(
        &self, cx: &Cx<'_>, boundary: &ThermalBoundary, temperatures: &[f64],
    ) -> Result<TemperatureContactEvaluation, ConductionError> {
        poll(cx, 0)?;
        if temperatures.len() != self.mesh.vertex_count() {
            return Err(invalid("<contact-set>", "temperature field and contact mesh differ"));
        }
        for (i, &t) in temperatures.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            positive(t)?;
        }
        self.reference.validate_for(self.mesh,boundary)?;
        let mut surfaces = Vec::new();
        let mut states = Vec::new();
        states.try_reserve_exact(self.variable_pairs)
            .map_err(|_| invalid("<contact-set>", "local source-state allocation refused"))?;
        for law in &self.variable {
            let ContactSurface::Matching(reference)=&law.reference else {continue;};
            let mut resistances = Vec::new();
            resistances.try_reserve_exact(reference.face_pairs().len())
                .map_err(|_| invalid(law.name(), "local resistance allocation refused"))?;
            for &pair in reference.face_pairs() {
                poll(cx, states.len())?;
                let mean = |slot: usize| -> Result<f64, ConductionError> {
                    let face = &self.mesh.boundary()[slot];
                    let mut total = 0.0;
                    for &vertex in &face.vertices {
                        total = finite(total + temperatures[vertex as usize] / 3.0)?;
                    }
                    positive(total)
                };
                let a = mean(pair.side_a)?;
                let b = mean(pair.side_b)?;
                let temperature = positive(0.5 * a + 0.5 * b)?;
                let resistance = law.at_traces([a, b])?;
                resistances.push(resistance.clone());
                states.push(ContactFaceState {
                    interface: law.name().into(), pair, side_a_temperature_k: a,
                    side_b_temperature_k: b, interface_temperature_k: temperature, resistance,
                });
            }
            surfaces.push(InterfaceSurface::new_mapped(
                law.name(), reference.face_pairs().to_vec(), resistances,
            )?);
        }
        let mut interfaces=self.reference.clone();
        let mut quadrature_states=Vec::new();
        quadrature_states.try_reserve_exact(self.variable_quadrature_points)
            .map_err(|_|invalid("<contact-set>","local quadrature state allocation refused"))?;
        let mut quadrature_owned=std::collections::BTreeSet::new();
        for surface in &mut interfaces.nonmatching {
            let Some(law)=self.variable.iter().find(|law|law.name()==surface.name) else {continue;};
            *surface=surface.at_quadrature(cx,temperatures,|pair,index,weight,traces| {
                let resistance=law.at_traces(traces)?;
                quadrature_states.push(ContactQuadratureState {interface:law.name().into(),pair,
                    quadrature_index:index,weight_m2:weight,side_a_temperature_k:traces[0],
                    side_b_temperature_k:traces[1],resistance:resistance.clone()});
                Ok(resistance)
            })?;
            quadrature_owned.insert(law.name().to_string());
        }
        interfaces.refresh_matching(&surfaces,&quadrature_owned)?;
        poll(cx, states.len())?;
        Ok(TemperatureContactEvaluation { interfaces, states,quadrature_states })
    }

    pub(crate) fn check_mesh(&self, mesh: &ConductionMesh) -> Result<(), ConductionError> {
        if std::ptr::eq(self.mesh, mesh) { Ok(()) }
        else { Err(invalid("<contact-set>", "contact and spatial problem must use the identical mesh")) }
    }
}

/// Fixed physical residual target and bounded Picard contact iterations.
/// Relaxation applies only to the next query field, never to accepted history.
#[derive(Debug, Clone, Copy)]
pub struct ContactIterationConfig {
    /// Maximum counted thermal responses. Nonzero; inner work remains separate.
    pub max_iterations: usize,
    /// Query-field relaxation in (0,1]. This does not regularize the source law.
    pub relaxation: f64,
    /// Absolute full free-node residual target, J, unchanged throughout the step.
    pub residual_tolerance_j: f64,
}

impl ContactIterationConfig {
    fn admit(self, inner_krylov_limit: usize) -> Result<usize, ConductionError> {
        if self.max_iterations == 0 || inner_krylov_limit == 0
            || !(self.relaxation.is_finite() && self.relaxation > 0.0 && self.relaxation <= 1.0)
            || !(self.residual_tolerance_j.is_finite() && self.residual_tolerance_j > 0.0)
        {
            return Err(invalid("<contact-step>", "positive contact work/residual budgets and relaxation in (0,1] are required"));
        }
        self.max_iterations.checked_mul(inner_krylov_limit)
            .ok_or_else(|| invalid("<contact-step>", "combined contact and Krylov work allowance overflows"))
    }

    fn relax(self, cx: &Cx<'_>, guess: &mut [f64], candidate: &[f64]) -> Result<(), ConductionError> {
        if guess.len() != candidate.len() {
            return Err(invalid("<contact-step>", "query and candidate field layouts differ"));
        }
        for (i, (old, &new)) in guess.iter_mut().zip(candidate).enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            *old = positive((1.0 - self.relaxation) * *old + self.relaxation * new)?;
        }
        poll(cx, guess.len())
    }
}

fn invalid(name: &str, what: impl Into<String>) -> ConductionError {
    ConductionError::Interface {
        interface: name.into(), what: what.into(),
        fix: "supply complete ordered contact traces and in-domain source laws; inspect the unchanged physical residual and work limits".into(),
    }
}
fn finite(value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() { Ok(value) }
    else { Err(invalid("<contact-step>", "nonfinite contact arithmetic")) }
}
fn positive(value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() && value > 0.0 { Ok(value) }
    else { Err(invalid("<contact-step>", "absolute temperatures must be finite and positive")) }
}
fn poll(cx: &Cx<'_>, at: usize) -> Result<(), ConductionError> {
    cx.checkpoint().map_err(|_| ConductionError::Cancelled { stage: "temperature-contact", at })
}
