//! Original material phase states on the retained ordered contact traces.
//!
//! The declared discretization interpolates each side's nodal temperature and
//! liquid MASS fraction separately. It never averages h across materials or
//! fabricates an equilibrium state from interpolated trace coordinates.

use super::*;
use crate::transient::enthalpy::{
    EnthalpyError, heterogeneous::HeterogeneousEnthalpyBackwardEuler,
};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState};
use fs_material::thermal::interface::ThermalInterfacePoint;

/// Explicit original A/B material and caloric owners for a contact law.
/// The exact curves must also be assigned to these vertices by the h solver.
#[derive(Debug, Clone, Copy)]
pub struct MaterialContactState<'a> {
    pub point: ThermalInterfacePoint<'a>,
    pub curves: [&'a EquilibriumEnthalpyPhaseCurve; 2],
}

/// One original caloric state, before P1 trace interpolation.
#[derive(Debug, Clone, PartialEq)]
pub struct ContactMaterialNodeState {
    pub vertex: usize,
    pub phase: EquilibriumPhaseState,
}

/// Explicit P1 fraction coordinates used by one original interface query.
/// `None` identifies the declared face-mean law; `Some` is the retained common
/// point index. Temperatures and the complete receipt remain in the matching
/// `ContactFaceState` or `ContactQuadratureState`.
#[derive(Debug, Clone, PartialEq)]
pub struct ContactMaterialTraceState {
    pub interface: String,
    pub pair: InterfaceFacePair,
    pub quadrature_index: Option<usize>,
    pub liquid_mass_fractions: [f64; 2],
}

#[derive(Debug, Clone, Copy)]
pub(super) struct BoundPhaseCurve<'a> {
    pub curve: &'a EquilibriumEnthalpyPhaseCurve,
    pub fraction_coordinate: bool,
}

pub(super) struct PhaseField {
    pub fractions: Vec<f64>,
    pub temperature_slopes: Vec<f64>,
    pub fraction_slopes: Vec<f64>,
    pub states: Vec<ContactMaterialNodeState>,
}

impl<'a> TemperatureContact<'a> {
    /// Follow the original material T(h) and liquid-fraction f(h) at every trial.
    /// The existing temperature convention must identify the corresponding
    /// ordered axes (or the explicitly declared mean). All other coordinates,
    /// selected claims, geometry and original mapped A-face contexts stay fixed.
    /// Temperature-only solves refuse this binding because T cannot identify h
    /// on a latent plateau. Use `advance_enthalpy_newton` and its accepted adjoint.
    pub fn with_material_states(
        mut self,
        material: MaterialContactState<'a>,
    ) -> Result<Self, ConductionError> {
        if self.material_states.is_some()
            || self.card.content_hash() != material.point.interface.content_hash()
            || self.reference_point != *material.point.context
        {
            return Err(invalid(
                self.name(),
                "material contact requires its original interface and complete reference context exactly once",
            ));
        }
        let matches = match &self.coordinates {
            TemperatureContactCoordinates::Mean(axis) => {
                material.point.mean_temperature_axis == Some(axis.as_str())
            }
            TemperatureContactCoordinates::TwoSided(axes) => axes
                .iter()
                .zip(material.point.temperature_axes)
                .all(|(a, b)| a == b),
        };
        if !matches {
            return Err(invalid(
                self.name(),
                "material contact temperature axes differ from the original ordered source binding",
            ));
        }
        let states = [0, 1].map(|side| {
            material.curves[side]
                .state_at_specific_enthalpy(material.curves[side].knots()[0].specific_enthalpy_j_kg)
        });
        let [a, b] = states;
        let states = [
            a.map_err(|e| invalid(self.name(), e.to_string()))?,
            b.map_err(|e| invalid(self.name(), e.to_string()))?,
        ];
        for (context, resistance) in std::iter::once((&self.reference_point, &self.resistance))
            .chain(
                self.face_references
                    .iter()
                    .map(|entry| (&entry.reference_point, &entry.resistance)),
            )
        {
            let binding = ThermalInterfacePoint {
                context,
                ..material.point
            };
            // This original owner checks ordered surface IDs, card identities,
            // distinct coordinates and exact quantity descriptors. Its result
            // is deliberately not queried: a chart endpoint need not lie in R's
            // source support when the actual reference query does.
            binding
                .at(states)
                .map_err(|e| invalid(self.name(), e.to_string()))?;
            let claim = self
                .card
                .claims()
                .claim(resistance.receipt().selected)
                .ok_or_else(|| {
                    invalid(
                        self.name(),
                        "material contact lost its original selected claim",
                    )
                })?;
            if !claim
                .value
                .axes()
                .all(|(axis, dims)| binding.is_response_axis(axis, dims))
                || context.axes().iter().any(|(axis, value)| {
                    (axis == "omega"
                        || context
                            .axis_quantities()
                            .get(axis)
                            .and_then(|q| q.semantic_type())
                            .is_some_and(|s| matches!(s.kind(), QuantityKind::Frequency(_))))
                        && *value != 0.0
                })
            {
                return Err(invalid(
                    self.name(),
                    "material contact needs static response in declared thermal/phase or typed fixed external coordinates",
                ));
            }
        }
        self.material_states = Some(material);
        Ok(self)
    }

    pub(super) fn query_point_at_state(
        &self,
        side_a: Option<usize>,
        temperatures: [f64; 2],
        fractions: Option<[f64; 2]>,
    ) -> Result<QueryPoint, ConductionError> {
        let Some(material) = self.material_states else {
            return self.query_point_at(side_a, temperatures);
        };
        let fractions = fractions.ok_or_else(|| {
            invalid(
                self.name(),
                "material contact requires original enthalpy states, not a temperature-only query",
            )
        })?;
        let mut point = self.source_context(side_a)?.0.clone();
        for side in 0..2 {
            let t = positive(temperatures[side])?;
            if !fractions[side].is_finite() || !(0.0..=1.0).contains(&fractions[side]) {
                return Err(invalid(
                    self.name(),
                    "P1 liquid mass fraction must lie in [0,1]",
                ));
            }
            point = point
                .with(material.point.temperature_axes[side], t)
                .map_err(|e| invalid(self.name(), e.to_string()))?;
            if let Some(axis) = material.point.liquid_fraction_axes[side] {
                point = point
                    .with(axis, fractions[side])
                    .map_err(|e| invalid(self.name(), e.to_string()))?;
            }
        }
        if let Some(axis) = material.point.mean_temperature_axis {
            point = point
                .with(
                    axis,
                    positive(0.5 * temperatures[0] + 0.5 * temperatures[1])?,
                )
                .map_err(|e| invalid(self.name(), e.to_string()))?;
        }
        Ok(point)
    }
}

pub(super) fn prepare_curves<'a>(
    cx: &Cx<'_>,
    mesh: &ConductionMesh,
    laws: &[TemperatureContact<'a>],
) -> Result<Vec<Option<BoundPhaseCurve<'a>>>, ConductionError> {
    if laws.iter().all(|law| law.material_states.is_none()) {
        return Ok(Vec::new());
    }
    let mut curves: Vec<Option<BoundPhaseCurve<'a>>> = vec![None; mesh.vertex_count()];
    for (index, law) in laws.iter().enumerate() {
        poll(cx, index)?;
        let Some(material) = law.material_states else {
            continue;
        };
        let mut bind = |slot: usize, side: usize| -> Result<(), ConductionError> {
            poll(cx, slot)?;
            let face = mesh.boundary().get(slot).ok_or_else(|| {
                invalid(
                    law.name(),
                    "material contact face is outside the original mesh",
                )
            })?;
            for vertex in face.vertices {
                let vertex = vertex as usize;
                let entry = &mut curves[vertex];
                if entry.is_some_and(|old| !std::ptr::eq(old.curve, material.curves[side])) {
                    return Err(invalid(
                        law.name(),
                        "a shared contact node requires the same original caloric object",
                    ));
                }
                let fraction_coordinate = material.point.liquid_fraction_axes[side].is_some()
                    || entry.is_some_and(|old| old.fraction_coordinate);
                *entry = Some(BoundPhaseCurve {
                    curve: material.curves[side],
                    fraction_coordinate,
                });
            }
            Ok(())
        };
        match &law.reference {
            ContactSurface::Matching(surface) | ContactSurface::PointwiseMatching(surface) => {
                for pair in surface.face_pairs() {
                    bind(pair.side_a, 0)?;
                    bind(pair.side_b, 1)?;
                }
            }
            ContactSurface::Nonmatching(surface) => {
                for &slot in &surface.side_a {
                    bind(slot, 0)?;
                }
                for &slot in &surface.side_b {
                    bind(slot, 1)?;
                }
            }
        }
    }
    Ok(curves)
}

impl TemperatureContacts<'_, '_> {
    pub(crate) fn check_material_curves<'c>(
        &self,
        cx: &Cx<'_>,
        mut curve_at: impl FnMut(usize) -> Option<&'c EquilibriumEnthalpyPhaseCurve>,
    ) -> Result<(), ConductionError> {
        for (vertex, expected) in self.phase_curves.iter().enumerate() {
            if vertex % 256 == 0 {
                poll(cx, vertex)?;
            }
            if let Some(expected) = expected
                && !curve_at(vertex).is_some_and(|actual| std::ptr::eq(actual, expected.curve))
            {
                return Err(invalid(
                    "<contact-set>",
                    "material contact must retain the native vertex's exact caloric object",
                ));
            }
        }
        Ok(())
    }

    pub(super) fn phase_field(
        &self,
        cx: &Cx<'_>,
        h: &[f64],
        temperature: &[f64],
        derivatives: bool,
        strict: bool,
    ) -> Result<Option<PhaseField>, ConductionError> {
        if self.phase_curves.is_empty() {
            return Ok(None);
        }
        let n = self.mesh.vertex_count();
        if h.len() != n || temperature.len() != n {
            return Err(invalid(
                "<contact-set>",
                "material contact requires complete original h and temperature fields",
            ));
        }
        let mut result = PhaseField {
            fractions: vec![0.0; n],
            temperature_slopes: vec![0.0; n],
            fraction_slopes: vec![0.0; n],
            states: Vec::new(),
        };
        for (vertex, binding) in self.phase_curves.iter().enumerate() {
            if vertex % 256 == 0 {
                poll(cx, vertex)?;
            }
            let Some(binding) = binding else {
                continue;
            };
            let curve = binding.curve;
            let failure = |error: fs_material::phase::PhaseStateError| {
                invalid("<contact-set>", error.to_string())
            };
            let phase = curve
                .state_at_specific_enthalpy(h[vertex])
                .map_err(failure)?;
            if phase.temperature_k().to_bits() != temperature[vertex].to_bits() {
                return Err(invalid(
                    "<contact-set>",
                    "contact temperature differs from the original nodal caloric state",
                ));
            }
            result.fractions[vertex] = phase.liquid_mass_fraction();
            if derivatives {
                if strict {
                    admit_chart_derivative(curve, h[vertex], binding.fraction_coordinate)?;
                }
                result.temperature_slopes[vertex] = curve
                    .temperature_derivative_at_specific_enthalpy(h[vertex])
                    .map_err(failure)?;
                result.fraction_slopes[vertex] = curve
                    .liquid_mass_fraction_derivative_at_specific_enthalpy(h[vertex])
                    .map_err(failure)?;
            }
            result
                .states
                .push(ContactMaterialNodeState { vertex, phase });
        }
        Ok(Some(result))
    }

    /// Re-evaluate original material phase contacts without advancing history.
    /// The two independent P1 T/f traces use the same quadrature as their heat.
    pub fn evaluate_enthalpy(
        &self,
        cx: &Cx<'_>,
        boundary: &ThermalBoundary,
        engine: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
        h: &[f64],
    ) -> Result<TemperatureContactEvaluation, EnthalpyError> {
        let (temperature, _) = self.enthalpy_fields(cx, engine, h, false)?;
        Ok(self.evaluate_enthalpy_with_query_refusal(cx, boundary, h, &temperature, &mut |_| {})?)
    }

    pub(crate) fn evaluate_enthalpy_with_query_refusal(
        &self,
        cx: &Cx<'_>,
        boundary: &ThermalBoundary,
        h: &[f64],
        temperature: &[f64],
        on_refusal: &mut dyn FnMut(&fs_matdb::MatDbError),
    ) -> Result<TemperatureContactEvaluation, ConductionError> {
        let phase = self.phase_field(cx, h, temperature, false, false)?;
        self.evaluate_state(cx, boundary, temperature, phase.as_ref(), on_refusal)
    }

    pub(super) fn enthalpy_fields(
        &self,
        cx: &Cx<'_>,
        engine: &HeterogeneousEnthalpyBackwardEuler<'_, '_>,
        h: &[f64],
        slopes: bool,
    ) -> Result<(Vec<f64>, Vec<f64>), EnthalpyError> {
        self.check_mesh(engine.mesh())?;
        self.check_material_curves(cx, |vertex| engine.phase_curve_at(vertex))?;
        if h.len() != self.mesh.vertex_count() {
            return Err(EnthalpyError::InvalidInput(
                "contact enthalpy must match its original mesh",
            ));
        }
        let mut temperature = Vec::with_capacity(h.len());
        let mut derivative = if slopes {
            Vec::with_capacity(h.len())
        } else {
            Vec::new()
        };
        for (vertex, &value) in h.iter().enumerate() {
            if vertex % 256 == 0 {
                poll(cx, vertex)?;
            }
            let curve = engine
                .phase_curve_at(vertex)
                .expect("mesh and native chart layout checked");
            temperature.push(
                curve
                    .temperature_at_specific_enthalpy(value)
                    .map_err(|source| EnthalpyError::Phase { vertex, source })?,
            );
            if slopes {
                derivative.push(
                    curve
                        .temperature_derivative_at_specific_enthalpy(value)
                        .map_err(|source| EnthalpyError::Phase { vertex, source })?,
                );
            }
        }
        Ok((temperature, derivative))
    }
}

pub(super) fn admit_chart_derivative(
    curve: &EquilibriumEnthalpyPhaseCurve,
    h: f64,
    fraction: bool,
) -> Result<(), ConductionError> {
    let knots = curve.knots();
    if h <= knots[0].specific_enthalpy_j_kg || h >= knots[knots.len() - 1].specific_enthalpy_j_kg {
        return Err(invalid(
            "<contact-adjoint>",
            "material contact needs an interior two-sided caloric derivative",
        ));
    }
    let upper = knots
        .partition_point(|knot| knot.specific_enthalpy_j_kg <= h)
        .min(knots.len() - 1);
    if upper > 1 && h == knots[upper - 1].specific_enthalpy_j_kg {
        let slope = |a: usize, b: usize| {
            let dh = knots[b].specific_enthalpy_j_kg - knots[a].specific_enthalpy_j_kg;
            [
                (knots[b].temperature_k - knots[a].temperature_k) / dh,
                (knots[b].liquid_mass_fraction - knots[a].liquid_mass_fraction) / dh,
            ]
        };
        let left = slope(upper - 2, upper - 1);
        let right = slope(upper - 1, upper);
        if left[0] != right[0] || fraction && left[1] != right[1] {
            return Err(invalid(
                "<contact-adjoint>",
                "material contact has a caloric or liquid-fraction slope corner",
            ));
        }
    }
    Ok(())
}
