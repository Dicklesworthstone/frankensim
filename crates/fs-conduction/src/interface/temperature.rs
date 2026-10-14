//! Source-backed contact resistance evaluated at each pair's mean trace temperature.
//!
//! Only the explicit temperature coordinate changes. Pressure, finish, medium,
//! source card and selected claim remain fixed. The two bulk traces stay distinct;
//! their arithmetic mean is a declared interface closure, not thermal equilibrium.
//! Frozen snapshots are NOT complete nonlinear-contact tangents or adjoints.

pub mod capacity;
/// Primary-enthalpy contact evolution, including latent phase states.
pub mod enthalpy;

use fs_exec::Cx;
use fs_matdb::{InterfaceSystemCard, QueryPoint, SelectionPolicy};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

use super::{InterfaceFacePair, InterfaceResistance, InterfaceSurface, ThermalInterfaces};
use crate::{ConductionError, ConductionMesh, ThermalBoundary};

/// A named matching surface and one source law, selected once at a reference point.
#[derive(Debug, Clone)]
pub struct TemperatureContact<'a> {
    card: &'a InterfaceSystemCard,
    reference_point: QueryPoint,
    temperature_axis: String,
    reference: InterfaceSurface,
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
        let name = name.into();
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature, ValueForm::Static,
        ));
        let descriptor_ok = match point.axis_quantities().get(temperature_axis) {
            Some(q) => *q == absolute || *q == QuantitySpec::dimensional(Temperature::DIMS),
            None => temperature_axis == "T",
        };
        if !descriptor_ok || !point.axes().get(temperature_axis)
            .is_some_and(|t| t.is_finite() && *t > 0.0)
        {
            return Err(invalid(&name, "an explicit positive absolute-temperature coordinate is required"));
        }
        let resistance = InterfaceResistance::from_card(&name, card, point, selection)?;
        let reference = InterfaceSurface::new(name, face_pairs, resistance.clone())?;
        Ok(Self {
            card, reference_point: point.clone(), temperature_axis: temperature_axis.into(),
            reference, resistance,
        })
    }

    /// Stable interface name, retaining the ordered card's surface convention.
    #[must_use]
    pub fn name(&self) -> &str { self.reference.name() }
    /// Original query including fixed pressure and all other source conditions.
    #[must_use]
    pub fn reference_point(&self) -> &QueryPoint { &self.reference_point }
    /// Reference source answer. Its selected claim stays pinned at every trial.
    #[must_use]
    pub fn reference_resistance(&self) -> &InterfaceResistance { &self.resistance }

    fn at(&self, temperature_k: f64) -> Result<InterfaceResistance, ConductionError> {
        positive(temperature_k)?;
        let point = self.reference_point.clone().with(&self.temperature_axis, temperature_k)
            .map_err(|error| invalid(self.name(), format!("temperature query refused: {error}")))?;
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
    /// Declared closure: half of each trace mean, K.
    pub interface_temperature_k: f64,
    /// Fresh positive source answer and its exact local query receipt.
    pub resistance: InterfaceResistance,
}

/// Contact operator frozen at a stated field, with the complete local source data.
#[derive(Debug, Clone)]
pub struct TemperatureContactEvaluation {
    interfaces: ThermalInterfaces,
    states: Vec<ContactFaceState>,
}

impl TemperatureContactEvaluation {
    /// Frozen operator for a primal trial or residual evaluation only.
    /// Its fixed-resistance adjoint omits dR/dT and is NOT a total derivative.
    #[must_use]
    pub fn frozen_interfaces(&self) -> &ThermalInterfaces { &self.interfaces }
    /// Source-backed state for each temperature-dependent face pair.
    #[must_use]
    pub fn states(&self) -> &[ContactFaceState] { &self.states }
}

/// Complete matching-contact set, including any unchanged constant interfaces.
/// Geometry is checked at construction. Rebinding reuses the original mesh,
/// source cards and paired faces; it never substitutes perfect contact.
#[derive(Debug)]
pub struct TemperatureContacts<'m, 'c> {
    mesh: &'m ConductionMesh,
    fixed: Vec<InterfaceSurface>,
    variable: Vec<TemperatureContact<'c>>,
    variable_pairs: usize,
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
        mut variable: Vec<TemperatureContact<'c>>,
        max_face_pairs: usize,
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        if variable.is_empty() || fixed.len().checked_add(variable.len()).is_none_or(|n| n > 64) {
            return Err(invalid("<contact-set>", "one to 64 interfaces and at least one temperature law are required"));
        }
        let mut count = 0usize;
        for surface in fixed.iter().chain(variable.iter().map(|law| &law.reference)) {
            poll(cx, count)?;
            count = count.checked_add(surface.face_pairs().len())
                .ok_or_else(|| invalid("<contact-set>", "contact face count overflows"))?;
            if count > max_face_pairs {
                return Err(invalid(surface.name(), "contact faces exceed the declared pair allowance"));
            }
        }
        variable.sort_by(|a, b| a.name().cmp(b.name()));
        let variable_pairs = variable.iter().map(|law| law.reference.face_pairs().len()).sum();
        let mut all = fixed.clone();
        all.extend(variable.iter().map(|law| law.reference.clone()));
        // The original owner rejects missing, duplicated, shared-node,
        // noncoincident and externally occupied contact faces.
        ThermalInterfaces::new(mesh, boundary, all)?;
        poll(cx, count)?;
        Ok(Self { mesh, fixed, variable, variable_pairs })
    }

    /// Query every local interface temperature and build one complete snapshot.
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
        let mut surfaces = self.fixed.clone();
        let mut states = Vec::new();
        states.try_reserve_exact(self.variable_pairs)
            .map_err(|_| invalid("<contact-set>", "local source-state allocation refused"))?;
        for law in &self.variable {
            let mut resistances = Vec::new();
            resistances.try_reserve_exact(law.reference.face_pairs().len())
                .map_err(|_| invalid(law.name(), "local resistance allocation refused"))?;
            for &pair in law.reference.face_pairs() {
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
                let resistance = law.at(temperature)?;
                resistances.push(resistance.clone());
                states.push(ContactFaceState {
                    interface: law.name().into(), pair, side_a_temperature_k: a,
                    side_b_temperature_k: b, interface_temperature_k: temperature, resistance,
                });
            }
            surfaces.push(InterfaceSurface::new_mapped(
                law.name(), law.reference.face_pairs().to_vec(), resistances,
            )?);
        }
        let interfaces = ThermalInterfaces::new(self.mesh, boundary, surfaces)?;
        poll(cx, states.len())?;
        Ok(TemperatureContactEvaluation { interfaces, states })
    }

    fn check_mesh(&self, mesh: &ConductionMesh) -> Result<(), ConductionError> {
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
        fix: "supply complete matching contacts and in-domain source laws; inspect the unchanged physical residual and work limits".into(),
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
