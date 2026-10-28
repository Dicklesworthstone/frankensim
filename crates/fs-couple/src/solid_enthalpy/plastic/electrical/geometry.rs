//! A material current path follows the original native P1 displacement field.
//!
//! For a reference direction n, F = I + grad(u), J = det(F), and lambda = |F n|.
//! The uniform scalar conductor has L = lambda L0 and A = J A0 / lambda ONLY
//! when n is a principal stretch direction: the deformed electrodes must stay
//! perpendicular to the deformed current path. A sheared nonprincipal path
//! requires a spatial current solve and is refused here. This is a bounded
//! small-strain geometry correction, not finite-strain solid mechanics.
//!
//! Reference inventory and thermal geometry do not change. The original
//! conductor evaluates R from the same endpoint material bundle; the circuit
//! deposits I^2 R once. A resistor has no stored field energy, so changing its
//! geometry creates neither parameter work nor an additional heat source.

use super::{ElectricalSolidError, ElectricalSolidSourceError, ThermoplasticSolid};
use crate::material_circuit::{MAX_MATERIAL_CIRCUIT_BRANCHES, MaterialResistor};
use fs_exec::Cx;
use fs_material::conductor::ResolvedConductorState;
use fs_math::det::sqrt;
use fs_thermal::enthalpy::CellId;
use std::collections::{BTreeMap, BTreeSet};

/// One whole-cell uniform bulk resistor, addressed by its original circuit branch.
/// Its thermal cell and source geometry come from the prepared heat network.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DeformingConductorPath {
    pub branch: usize,
    /// Nonzero direction in the undeformed global material frame. Normalized
    /// once at binding; its sign has no effect on the scalar resistance.
    pub reference_direction: [f64; 3],
}

/// Explicit admission bounds for the uniform-path geometry reduction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductorGeometrySolve {
    /// Maximum Frobenius norm of grad(u), in (0,1). The native material's own
    /// small-strain/rotation limit is still enforced independently.
    pub maximum_displacement_gradient: f64,
    /// Relative residual |F^T F n - lambda^2 n| / lambda^2, in (0,1).
    /// A small numerical tolerance admits roundoff in a principal direction;
    /// this is not an accuracy certificate for a spatial current field.
    pub orthogonality_tolerance: f64,
}

#[derive(Clone, Debug)]
pub(in crate::solid_enthalpy::plastic) struct ConductorGeometryBinding {
    paths: Vec<DeformingConductorPath>,
    controls: ConductorGeometrySolve,
}

impl ThermoplasticSolid<'_> {
    /// Follow explicitly selected bulk current paths through the SAME mechanical
    /// endpoint as conductivity, Joule heat and enthalpy. No separate element or
    /// material mapping is accepted: the original heat route selects the tet.
    /// Reference L0/A0 come from each fresh prepared-network source evaluation,
    /// never from a previously stretched conductor. This prevents accumulated
    /// stretching when continuing an accepted trajectory.
    ///
    /// Only the named bulk paths deform. Contact areas, thermal transport,
    /// capacitances, inductances and topology retain their declared geometry.
    /// Isotropic DC response, fully solid P1 kinematics and principal paths only;
    /// source strain dependence, motional EMF and electromagnetic forces need
    /// their own constitutive/field owners. Empty paths retain fixed geometry.
    pub fn with_deforming_conductors(
        mut self,
        mut paths: Vec<DeformingConductorPath>,
        controls: ConductorGeometrySolve,
        cx: &Cx<'_>,
    ) -> Result<Self, ElectricalSolidError> {
        cx.checkpoint()
            .map_err(|_| ElectricalSolidError::Cancelled)?;
        if paths.len() > MAX_MATERIAL_CIRCUIT_BRANCHES
            || !controls.maximum_displacement_gradient.is_finite()
            || !(0.0..1.0).contains(&controls.maximum_displacement_gradient)
            || controls.maximum_displacement_gradient == 0.0
            || !controls.orthogonality_tolerance.is_finite()
            || !(0.0..1.0).contains(&controls.orthogonality_tolerance)
            || controls.orthogonality_tolerance == 0.0
        {
            return Err(ElectricalSolidError::Invalid(
                "deforming conductors need bounded positive geometry controls and branch count",
            ));
        }
        let mut seen = BTreeSet::new();
        for path in &mut paths {
            cx.checkpoint()
                .map_err(|_| ElectricalSolidError::Cancelled)?;
            let norm2 = dot(path.reference_direction, path.reference_direction);
            if !seen.insert(path.branch) || !norm2.is_finite() || norm2 <= 0.0 {
                return Err(ElectricalSolidError::Invalid(
                    "deforming conductors need distinct branches and finite nonzero reference directions",
                ));
            }
            let norm = sqrt(norm2);
            path.reference_direction = path.reference_direction.map(|v| v / norm);
        }
        self.conductor_geometry = if paths.is_empty() {
            None
        } else {
            Some(ConductorGeometryBinding { paths, controls })
        };
        Ok(self)
    }
}

pub(super) struct PreparedConductorGeometry {
    branch: usize,
    element: usize,
    direction: [f64; 3],
    length_m: f64,
    area_m2: f64,
    controls: ConductorGeometrySolve,
}

fn check(cx: &Cx<'_>) -> Result<(), ElectricalSolidSourceError> {
    cx.checkpoint()
        .map_err(|_| ElectricalSolidSourceError::Mechanical(super::ThermoplasticError::Cancelled))
}

pub(super) fn prepare(
    solid: &ThermoplasticSolid<'_>,
    resistors: &[MaterialResistor],
    routes: &BTreeMap<usize, CellId>,
    cx: &Cx<'_>,
) -> Result<Vec<PreparedConductorGeometry>, ElectricalSolidSourceError> {
    let Some(binding) = &solid.conductor_geometry else {
        return Ok(Vec::new());
    };
    let mut prepared = Vec::with_capacity(binding.paths.len());
    for path in &binding.paths {
        check(cx)?;
        let original = resistors.iter().find(|r| r.branch == path.branch).ok_or(
            ElectricalSolidSourceError::Invalid("deforming path has no original material resistor"),
        )?;
        let cell = routes
            .get(&path.branch)
            .ok_or(ElectricalSolidSourceError::Invalid(
                "deforming path has no original caloric destination",
            ))?;
        let element = solid
            .cells
            .iter()
            .position(|b| b.accepted.id == *cell)
            .ok_or(ElectricalSolidSourceError::Invalid(
                "deforming path must belong to a native solid element",
            ))?;
        let conductor = &original.conductor;
        if !crate::solid_enthalpy::same_measure(
            conductor.length_m() * conductor.area_m2(),
            solid.cells[element].volume_m3,
        ) {
            return Err(ElectricalSolidSourceError::Invalid(
                "deforming conductor must retain its whole-cell reference volume",
            ));
        }
        prepared.push(PreparedConductorGeometry {
            branch: path.branch,
            element,
            direction: path.reference_direction,
            length_m: conductor.length_m(),
            area_m2: conductor.area_m2(),
            controls: binding.controls,
        });
    }
    Ok(prepared)
}

pub(super) fn refresh(
    paths: &[PreparedConductorGeometry],
    gradients: &[[[f64; 3]; 3]],
    resistors: &mut [MaterialResistor],
    cx: &Cx<'_>,
) -> Result<(), ElectricalSolidSourceError> {
    for path in paths {
        check(cx)?;
        let gradient = gradients
            .get(path.element)
            .ok_or(ElectricalSolidSourceError::Invalid(
                "deforming conductor lost its native element gradient",
            ))?;
        let norm2: f64 = gradient.iter().flatten().map(|g| g * g).sum();
        if !norm2.is_finite()
            || norm2
                > path.controls.maximum_displacement_gradient
                    * path.controls.maximum_displacement_gradient
        {
            return Err(ElectricalSolidSourceError::Invalid(
                "conductor displacement gradient exceeds its admitted affine geometry bound",
            ));
        }
        let mut f = *gradient;
        for (i, row) in f.iter_mut().enumerate() {
            row[i] += 1.0;
        }
        let j = f[0][0] * (f[1][1] * f[2][2] - f[1][2] * f[2][1])
            - f[0][1] * (f[1][0] * f[2][2] - f[1][2] * f[2][0])
            + f[0][2] * (f[1][0] * f[2][1] - f[1][1] * f[2][0]);
        let tangent = f.map(|row| dot(row, path.direction));
        let stretch2 = dot(tangent, tangent);
        let cn: [f64; 3] = std::array::from_fn(|a| {
            f[0][a] * tangent[0] + f[1][a] * tangent[1] + f[2][a] * tangent[2]
        });
        let residual: [f64; 3] = std::array::from_fn(|a| cn[a] - stretch2 * path.direction[a]);
        let defect = sqrt(dot(residual, residual));
        if !j.is_finite()
            || j <= 0.0
            || !stretch2.is_finite()
            || stretch2 <= 0.0
            || !defect.is_finite()
            || defect > path.controls.orthogonality_tolerance * stretch2
        {
            return Err(ElectricalSolidSourceError::Invalid(
                "deformed uniform conductor requires positive volume and a principal current direction",
            ));
        }
        let stretch = sqrt(stretch2);
        let resistor = resistors
            .iter_mut()
            .find(|r| r.branch == path.branch)
            .ok_or(ElectricalSolidSourceError::Invalid(
                "deforming conductor disappeared from the endpoint material state",
            ))?;
        // Adopt the actual same-point source bundle; do not resolve or fit a
        // separate material coefficient, and do not change any circuit history.
        resistor.conductor = ResolvedConductorState::try_from_material_with_property(
            resistor.conductor.material().clone(),
            resistor.conductor.source_basis(),
            resistor.conductor.source_key(),
        )
        .and_then(|state| {
            state.into_uniform_conductor(path.length_m * stretch, path.area_m2 * j / stretch)
        })
        .map_err(|e| ElectricalSolidSourceError::Circuit(e.into()))?;
    }
    check(cx)
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
