//! Source-resolved electrical contacts on matching original exterior traces.
//!
//! The surface term integrates (V_a-V_b)(test_a-test_b)/R'' with the same P1
//! trace basis as the bulk FEEC operator. Duplicate, coincident vertices retain
//! the physical voltage jump; a small fictitious conducting layer is unnecessary.

use super::contact_quadrature::{
    ConductionContactQuadraturePoint, ContactCoefficient, ContactQuadrature, QuadratureError,
    edge_point_count, edge_rule,
};
use super::*;
use crate::planar::{MagneticError, PlanarMagnetostatics};
use fs_material::electrical_interface::{
    ELECTRICAL_CONTACT_RESISTANCE_PROPERTY, ResolvedElectricalInterface,
};
use fs_math::dd::Dd;
use fs_qty::{Dims, QuantitySpec};
use fs_sparse::Coo;

/// One ordered, zero-thickness material interface. The two original exterior
/// edges must have distinct vertices with matching endpoint coordinates, in
/// this exact pairing order, and their adjacent cells must lie on opposite sides.
/// Source surface A belongs to `first_vertices`; B belongs to `second_vertices`.
#[derive(Clone, Debug, PartialEq)]
pub struct ConductionContact {
    pub first_vertices: [usize; 2],
    pub second_vertices: [usize; 2],
    pub source: ResolvedElectricalInterface,
}

/// Physical contact response, oriented from the first surface to the second.
/// Voltage jump is P1 along the paired edge. Uniform coefficients give P1
/// current density and exact surface integrals. Source-sampled coefficients
/// retain endpoint densities and their admitted quadrature observations.
#[derive(Clone, Debug, PartialEq)]
pub struct ConductionContactResponse {
    pub voltage_jump_v: [f64; 2],
    pub current_density_a_m2: [f64; 2],
    pub current_a: f64,
    pub power_w: f64,
    /// Actual source-sampled integration points, empty for the exact uniform law.
    pub quadrature: Vec<ConductionContactQuadraturePoint<2>>,
}

// Cache the exact two-node surface moments once with the native stiffness.
// The endpoint weights are 1,1 for a planar edge and 2*r_i/(r_0+r_1)
// for a surface of revolution. With g=A/R'', the mass entries are
// g/6 * [[1+w_0, 1], [1, 1+w_1]]. Row sums integrate current;
// the positive-square quadratic form integrates irreversible heat.
#[derive(Clone)]
pub(super) struct UniformContactConductance {
    total_s: f64,
    endpoint_weights: [f64; 2],
}

impl UniformContactConductance {
    fn matrix_entry(&self, i: usize, j: usize) -> f64 {
        (self.total_s / 6.0)
            * if i == j {
                1.0 + self.endpoint_weights[i]
            } else {
                1.0
            }
    }

    fn current_a(&self, jump: [f64; 2]) -> f64 {
        self.total_s
            * (((2.0 + self.endpoint_weights[0]) / 6.0) * jump[0]
                + ((2.0 + self.endpoint_weights[1]) / 6.0) * jump[1])
    }

    fn power_w(&self, jump: [f64; 2]) -> f64 {
        let [d0, d1] = jump;
        let sum = d0 + d1;
        (self.total_s / 6.0)
            * (self.endpoint_weights[0] * d0 * d0 + self.endpoint_weights[1] * d1 * d1 + sum * sum)
    }
}

#[derive(Clone)]
pub(super) enum ContactConductance {
    Uniform(UniformContactConductance),
    Sampled(ContactQuadrature<2>),
}

impl ContactConductance {
    fn matrix_entry(&self, i: usize, j: usize) -> f64 {
        match self {
            Self::Uniform(exact) => exact.matrix_entry(i, j),
            Self::Sampled(quadrature) => quadrature.matrix_entry(i, j),
        }
    }
}

fn quadrature_error(error: QuadratureError) -> ConductionError {
    match error {
        QuadratureError::Invalid(reason) => ConductionError::Invalid(reason),
        QuadratureError::Cancelled => ConductionError::Cancelled,
    }
}

// The trace has positive swept area and nonnegative radii. Scaling before
// summing keeps physical moments usable when r_0+r_1 would overflow.
fn radial_contact_weights(radii: [f64; 2]) -> [f64; 2] {
    let scale = radii[0].max(radii[1]);
    let scaled = radii.map(|r| r / scale);
    let sum = scaled[0] + scaled[1];
    scaled.map(|r| 2.0 * (r / sum))
}

fn edge_key(mut vertices: [usize; 2]) -> [usize; 2] {
    vertices.sort_unstable();
    vertices
}

// Planar and tetrahedral traces share the same original interface-source law.
// Only their geometry and trace quadrature differ.
pub(crate) fn admit_contact_source(
    source: &ResolvedElectricalInterface,
    axes: &[Option<String>; 2],
) -> Result<f64, ConductionError> {
    let mut owned: Vec<(&str, QuantitySpec)> = Vec::new();
    for axis in axes.iter().flatten() {
        if axis.trim().is_empty()
            || owned.iter().any(|(previous, _)| *previous == axis.as_str())
            || source.resolved().axis_quantities().get(axis)
                != Some(&QuantitySpec::dimensional(Dims::NONE))
            || !source
                .resolved()
                .state_coordinate(axis)
                .is_some_and(|value| value.is_finite() && (0.0..=1.0).contains(&value))
        {
            return Err(ConductionError::Invalid(
                "owned contact phase needs a unique existing dimensional mass fraction in [0,1]",
            ));
        }
        owned.push((axis.as_str(), QuantitySpec::dimensional(Dims::NONE)));
    }
    if !source
        .resolved()
        .property(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY)
        .is_some_and(|property| {
            if owned.is_empty() {
                property.is_static_at_fixed_external_conditions()
            } else {
                property.has_only_owned_coordinates_and_fixed_external_conditions(&owned)
            }
        })
    {
        return Err(ConductionError::Invalid(
            "contact resistance requires static source support at fixed external conditions",
        ));
    }
    let resistance = source.resistance_ohm_m2();
    if !positive(resistance) {
        return Err(ConductionError::Invalid(
            "positive specific contact resistance required; perfect contact needs shared mesh vertices",
        ));
    }
    Ok(resistance)
}

impl<C: PlanarConductorMaterial> PlanarConduction<C> {
    /// Assemble bulk and source-backed surface conduction in one native field.
    /// The specific contact resistance R'' [ohm m2] is divided by actual area,
    /// with the exact P1 surface mass matrix, not a lumped endpoint resistor.
    /// Material, pressure or temperature changes require resolving the original
    /// interface source again and preparing a new operator.
    ///
    /// Finite positive, static R'' is required. Explicit zero resistance belongs
    /// on a conforming shared-node mesh; this owner does not invent a penalty.
    /// No surface capacitance, moving contact, or thermal heat partition is inferred.
    pub fn new_with_contacts(
        mesh: TriComplex2,
        conductors: Vec<C>,
        contacts: Vec<ConductionContact>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ConductionError> {
        poll(cancelled)?;
        if contacts.len() > limits.max_triplets / 16 {
            return Err(ConductionError::Invalid(
                "electrical contacts exceed the declared triplet capacity",
            ));
        }
        let axes = vec![[None, None]; contacts.len()];
        Self::new_with_phase_contacts(mesh, conductors, contacts, axes, limits, cancelled)
    }

    /// The SAME native surface operator with explicit externally owned liquid
    /// mass fractions on either ordered side. Each declared coordinate must be
    /// present as a dimensional scalar in [0,1]. Every other source response
    /// and applicability coordinate retains the fixed-external admission rule.
    ///
    /// A thermal/phase owner must resolve these exact values from its accepted
    /// caloric state before assembly, then refresh the same sources on every
    /// proposed change. This method neither determines phase from temperature
    /// nor adds phase mixing, motion, storage or an electrical source policy.
    /// Axis rows correspond to INPUT contacts and are sorted with their traces;
    /// returned axes use the same canonical order as `contacts()`.
    pub fn new_with_phase_contacts(
        mesh: TriComplex2,
        conductors: Vec<C>,
        contacts: Vec<ConductionContact>,
        liquid_mass_fraction_axes: Vec<[Option<String>; 2]>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ConductionError> {
        Self::assemble_with_phase_contacts(
            mesh,
            conductors,
            contacts,
            liquid_mass_fraction_axes,
            limits,
            cancelled,
            |source| source.has_admitted_static_sources(),
        )
    }

    // Constitutive response owners admit their own bulk coordinates; geometry,
    // interface support, surface quadrature and current/heat stay native here.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn assemble_with_phase_contacts(
        mesh: TriComplex2,
        conductors: Vec<C>,
        contacts: Vec<ConductionContact>,
        liquid_mass_fraction_axes: Vec<[Option<String>; 2]>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
        admit: impl Fn(&C) -> bool,
    ) -> Result<Self, ConductionError> {
        Self::assemble_with_admitted_contacts(
            mesh,
            conductors,
            contacts,
            liquid_mass_fraction_axes,
            limits,
            cancelled,
            admit,
            |_, source, axes| admit_contact_source(source, axes),
        )
    }

    // The response owner may admit its explicitly owned local electrical axes.
    // All geometry, passivity, assembly and heat remain in this SAME owner.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn assemble_with_admitted_contacts(
        mesh: TriComplex2,
        conductors: Vec<C>,
        contacts: Vec<ConductionContact>,
        liquid_mass_fraction_axes: Vec<[Option<String>; 2]>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
        admit: impl Fn(&C) -> bool,
        admit_contact: impl Fn(
            usize,
            &ResolvedElectricalInterface,
            &[Option<String>; 2],
        ) -> Result<f64, ConductionError>,
    ) -> Result<Self, ConductionError> {
        poll(cancelled)?;
        if !matches!(mesh.metric(), Metric2::Planar { .. }) {
            return Err(ConductionError::Invalid(
                "this conductor owner requires a planar metric",
            ));
        }
        Self::assemble_metric_with_admitted_contacts(
            mesh,
            conductors,
            contacts,
            liquid_mass_fraction_axes,
            limits,
            cancelled,
            admit,
            admit_contact,
        )
    }

    // Metric-specific facades admit their geometry before entering this shared
    // source, topology, surface-moment and field-recovery owner.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn assemble_metric_with_admitted_contacts(
        mesh: TriComplex2,
        conductors: Vec<C>,
        contacts: Vec<ConductionContact>,
        liquid_mass_fraction_axes: Vec<[Option<String>; 2]>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
        admit: impl Fn(&C) -> bool,
        admit_contact: impl Fn(
            usize,
            &ResolvedElectricalInterface,
            &[Option<String>; 2],
        ) -> Result<f64, ConductionError>,
    ) -> Result<Self, ConductionError> {
        Self::assemble_metric_with_contact_coefficients(
            mesh,
            conductors,
            contacts,
            liquid_mass_fraction_axes,
            limits,
            cancelled,
            admit,
            |index, source, axes| {
                admit_contact(index, source, axes).map(ContactCoefficient::Uniform)
            },
        )
    }

    // A nonlinear owner admits each original source sample before returning
    // its coefficient. Geometry, positive quadrature, assembly and heat remain
    // native; the exact fixed-coefficient constructor keeps its original path.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn assemble_metric_with_contact_coefficients(
        mesh: TriComplex2,
        conductors: Vec<C>,
        contacts: Vec<ConductionContact>,
        liquid_mass_fraction_axes: Vec<[Option<String>; 2]>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
        admit: impl Fn(&C) -> bool,
        admit_contact: impl Fn(
            usize,
            &ResolvedElectricalInterface,
            &[Option<String>; 2],
        ) -> Result<ContactCoefficient, ConductionError>,
    ) -> Result<Self, ConductionError> {
        poll(cancelled)?;
        if liquid_mass_fraction_axes.len() != contacts.len() {
            return Err(ConductionError::Invalid(
                "one phase-coordinate row per electrical contact required",
            ));
        }
        let triplets = mesh.faces().len().checked_mul(9).and_then(|bulk| {
            contacts
                .len()
                .checked_mul(16)
                .and_then(|surface| bulk.checked_add(surface))
        });
        if triplets.is_none_or(|count| count > limits.max_triplets) {
            return Err(ConductionError::Invalid(
                "bulk and contact assembly exceed the declared triplet capacity",
            ));
        }
        let mut problem = Self::assemble_metric(mesh, conductors, limits, cancelled, admit)?;
        if contacts.is_empty() {
            return Ok(problem);
        }
        let trace = problem.mesh.boundary_trace().map_err(|_| {
            ConductionError::Invalid("electrical contacts require original exterior traces")
        })?;
        let exterior: BTreeMap<_, _> = trace
            .edges
            .iter()
            .map(|edge| (edge_key(edge.oriented_vertices.map(|v| v as usize)), edge))
            .collect();
        let mut bound: Vec<_> = contacts
            .into_iter()
            .zip(liquid_mass_fraction_axes)
            .collect();
        bound.sort_by_key(|(contact, _)| {
            (
                edge_key(contact.first_vertices),
                edge_key(contact.second_vertices),
            )
        });
        let (contacts, liquid_mass_fraction_axes): (Vec<_>, Vec<_>) = bound.into_iter().unzip();
        let n = problem.mesh.vertices().len();
        let mut matrix = Coo::new(n, n);
        for row in 0..n {
            poll(cancelled)?;
            let (columns, values) = problem.stiffness.row(row);
            for (&column, &value) in columns.iter().zip(values) {
                matrix.push(row, column, value);
            }
        }
        for (contact_index, (contact, axes)) in
            contacts.iter().zip(&liquid_mass_fraction_axes).enumerate()
        {
            poll(cancelled)?;
            let first = edge_key(contact.first_vertices);
            let second = edge_key(contact.second_vertices);
            let Some(first_edge) = exterior.get(&first) else {
                return Err(ConductionError::Invalid(
                    "first contact trace is not an exterior edge",
                ));
            };
            let Some(second_edge) = exterior.get(&second) else {
                return Err(ConductionError::Invalid(
                    "second contact trace is not an exterior edge",
                ));
            };
            if first.iter().any(|vertex| second.contains(vertex))
                || !problem.contact_edges.insert(first)
                || !problem.contact_edges.insert(second)
            {
                return Err(ConductionError::Invalid(
                    "contact traces require distinct original vertices and one interface per edge",
                ));
            }
            for endpoint in 0..2 {
                if problem.mesh.vertices()[contact.first_vertices[endpoint]]
                    != problem.mesh.vertices()[contact.second_vertices[endpoint]]
                {
                    return Err(ConductionError::Invalid(
                        "contact endpoint pairs must have identical physical coordinates",
                    ));
                }
            }
            let [a, b] = contact.first_vertices.map(|v| problem.mesh.vertices()[v]);
            let side = |face: usize| -> Result<f64, ConductionError> {
                let vertex = problem.mesh.faces()[face]
                    .iter()
                    .copied()
                    .find(|&v| {
                        let p = problem.mesh.vertices()[v as usize];
                        p != a && p != b
                    })
                    .ok_or(ConductionError::Invalid(
                        "contact has no adjacent cell interior",
                    ))?;
                let p = problem.mesh.vertices()[vertex as usize];
                finite((b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]))
            };
            let first_side = side(first_edge.source_face)?;
            let second_side = side(second_edge.source_face)?;
            if first_side == 0.0
                || second_side == 0.0
                || first_side.is_sign_positive() == second_side.is_sign_positive()
            {
                return Err(ConductionError::Invalid(
                    "contact cells must lie on opposite sides of the matched trace",
                ));
            }
            contact
                .source
                .require_surface_materials(
                    problem.conductors[first_edge.source_face]
                        .material()
                        .material(),
                    problem.conductors[second_edge.source_face]
                        .material()
                        .material(),
                )
                .map_err(|_| {
                    ConductionError::Invalid(
                        "contact source must match both ordered adjacent material states",
                    )
                })?;
            let coefficient = admit_contact(contact_index, &contact.source, axes)?;
            let area = problem.mesh.edge_measure(first_edge.global_edge).ok_or(
                ConductionError::Invalid("contact trace lacks physical area"),
            )?;
            let endpoint_weights = match problem.mesh.metric() {
                Metric2::Planar { .. } => [1.0; 2],
                Metric2::Axisymmetric { .. } => radial_contact_weights([a[0], b[0]]),
            };
            let moments = match coefficient {
                ContactCoefficient::Uniform(resistance) => {
                    let conductance = finite(area / resistance)?;
                    if !positive(conductance / 6.0) {
                        return Err(ConductionError::Invalid(
                            "contact conductance is not representable",
                        ));
                    }
                    ContactConductance::Uniform(UniformContactConductance {
                        total_s: conductance,
                        endpoint_weights,
                    })
                }
                ContactCoefficient::Sampled {
                    subdivisions,
                    resistance_ohm_m2,
                } => {
                    if resistance_ohm_m2.len()
                        != edge_point_count(subdivisions).map_err(quadrature_error)?
                    {
                        return Err(ConductionError::Invalid(
                            "complete sampled edge contact coefficients required",
                        ));
                    }
                    ContactConductance::Sampled(
                        ContactQuadrature::new(
                            edge_rule(subdivisions, cancelled).map_err(quadrature_error)?,
                            area,
                            endpoint_weights,
                            resistance_ohm_m2,
                            cancelled,
                        )
                        .map_err(quadrature_error)?,
                    )
                }
            };
            for i in 0..2 {
                for j in 0..2 {
                    let value = finite(moments.matrix_entry(i, j))?;
                    let ai = contact.first_vertices[i];
                    let aj = contact.first_vertices[j];
                    let bi = contact.second_vertices[i];
                    let bj = contact.second_vertices[j];
                    matrix.push(ai, aj, value);
                    matrix.push(bi, bj, value);
                    matrix.push(ai, bj, -value);
                    matrix.push(bi, aj, -value);
                }
            }
            problem.contact_conductance_s.push(moments);
        }
        poll(cancelled)?;
        problem.stiffness = matrix.assemble();
        for row in 0..n {
            poll(cancelled)?;
            if problem
                .stiffness
                .row(row)
                .1
                .iter()
                .any(|value| !value.is_finite())
            {
                return Err(ConductionError::Invalid(
                    "unrepresentable bulk/contact stiffness",
                ));
            }
        }
        problem.contacts = contacts;
        problem.contact_liquid_mass_fraction_axes = liquid_mass_fraction_axes;
        poll(cancelled)?;
        Ok(problem)
    }
}

impl<C: PlanarOhmicResponse> PlanarConduction<C> {
    /// Original ordered source bindings in canonical first/second edge order.
    pub fn contacts(&self) -> &[ConductionContact] {
        &self.contacts
    }

    /// Explicit phase-owner coordinates in the same canonical order as contacts.
    /// The fixed-external constructor records two `None` entries per contact.
    pub fn contact_liquid_mass_fraction_axes(&self) -> &[[Option<String>; 2]] {
        &self.contact_liquid_mass_fraction_axes
    }

    pub(super) fn is_contact_edge(&self, vertices: [usize; 2]) -> bool {
        self.contact_edges.contains(&edge_key(vertices))
    }

    pub(super) fn prescribed_boundary(
        &self,
        boundary: &[(usize, f64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<Option<f64>>, ConductionError> {
        if self.contacts.is_empty() {
            return Ok(scalar_potential::prescribed_boundary(
                &self.mesh, boundary, cancelled,
            )?);
        }
        let prescribed = PlanarMagnetostatics::admit_mesh_traces(&self.mesh, boundary, cancelled)
            .map_err(|error| match error {
            MagneticError::Cancelled => ConductionError::Cancelled,
            MagneticError::Invalid(reason) => ConductionError::Invalid(reason),
            _ => ConductionError::Invalid("invalid conductor boundary trace"),
        })?;
        for &[a, b] in &self.contact_edges {
            poll(cancelled)?;
            if prescribed[a].is_some() && prescribed[b].is_some() {
                return Err(ConductionError::Invalid(
                    "internal contact traces cannot also be voltage electrodes",
                ));
            }
        }
        let n = prescribed.len();
        let mut neighbors = vec![Vec::new(); n];
        for &[a, b] in self.mesh.edges() {
            poll(cancelled)?;
            neighbors[a as usize].push(b as usize);
            neighbors[b as usize].push(a as usize);
        }
        for contact in &self.contacts {
            poll(cancelled)?;
            for endpoint in 0..2 {
                let a = contact.first_vertices[endpoint];
                let b = contact.second_vertices[endpoint];
                neighbors[a].push(b);
                neighbors[b].push(a);
            }
        }
        let mut visited = vec![false; n];
        let mut stack = Vec::new();
        for start in 0..n {
            if visited[start] {
                continue;
            }
            visited[start] = true;
            stack.push(start);
            let mut anchored = false;
            while let Some(vertex) = stack.pop() {
                poll(cancelled)?;
                anchored |= prescribed[vertex].is_some();
                for &next in &neighbors[vertex] {
                    if !visited[next] {
                        visited[next] = true;
                        stack.push(next);
                    }
                }
            }
            if !anchored {
                return Err(ConductionError::Invalid(
                    "every contact-connected component needs a prescribed exterior edge",
                ));
            }
        }
        Ok(prescribed)
    }

    pub(super) fn contact_response(
        &self,
        voltage: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<ConductionContactResponse>, ConductionError> {
        let mut response = Vec::with_capacity(self.contacts.len());
        for (contact, moments) in self.contacts.iter().zip(&self.contact_conductance_s) {
            poll(cancelled)?;
            let [a0, a1] = contact.first_vertices;
            let [b0, b1] = contact.second_vertices;
            let d0 = finite(voltage[a0] - voltage[b0])?;
            let d1 = finite(voltage[a1] - voltage[b1])?;
            let (density, current, power, quadrature) = match moments {
                ContactConductance::Uniform(exact) => {
                    let resistance = contact.source.resistance_ohm_m2();
                    (
                        [finite(d0 / resistance)?, finite(d1 / resistance)?],
                        finite(exact.current_a([d0, d1]))?,
                        finite(exact.power_w([d0, d1]))?,
                        Vec::new(),
                    )
                }
                ContactConductance::Sampled(rule) => {
                    let points = rule
                        .response([d0, d1], cancelled)
                        .map_err(quadrature_error)?;
                    let first = points
                        .first()
                        .ok_or(ConductionError::Invalid("empty edge quadrature"))?;
                    let last = points
                        .last()
                        .ok_or(ConductionError::Invalid("empty edge quadrature"))?;
                    if first.basis != [1.0, 0.0] || last.basis != [0.0, 1.0] {
                        return Err(ConductionError::Invalid(
                            "edge quadrature lost original endpoints",
                        ));
                    }
                    let density = [first.current_density_a_m2, last.current_density_a_m2];
                    let (mut current, mut power) = (Dd::from_f64(0.0), Dd::from_f64(0.0));
                    for point in &points {
                        poll(cancelled)?;
                        current = current
                            + Dd::from_f64(point.area_weight_m2)
                                * Dd::from_f64(point.current_density_a_m2);
                        power = power + Dd::from_f64(point.power_w);
                    }
                    (
                        density,
                        finite(current.to_f64())?,
                        finite(power.to_f64())?,
                        points,
                    )
                }
            };
            if power < 0.0 || ((d0 != 0.0 || d1 != 0.0) && power == 0.0) {
                return Err(ConductionError::Invalid(
                    "contact Joule heat is not representably positive",
                ));
            }
            response.push(ConductionContactResponse {
                voltage_jump_v: [d0, d1],
                current_density_a_m2: density,
                current_a: current,
                power_w: power,
                quadrature,
            });
        }
        Ok(response)
    }
}
