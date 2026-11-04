//! Component positions and explicitly selected directional pressure drag.
//!
//! OrthogonalFaces uses each face's NORMAL relative speed, rather than the
//! aircraft's total speed. The caller supplies the projected areas and their
//! orientation; the original ledger owns Cd and its uncertainty/provenance.
//! This is a reduced face-pressure law, not resolved bluff-body interference,
//! skin friction, or a calibration of the registered Flyer in crossflow.

use super::{DragLedger, DragWrench, Refusal, hash_domain};

/// Directional area model at one physical component's application point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ComponentDragShape {
    /// The original scalar area and Cd, with force opposite the full local flow.
    Isotropic,
    /// Three orthogonal, two-sided pressure faces. For each normal n,
    /// F = -rho*Cd*A*|v dot n|*(v dot n)*n/2. Thus a zero-thickness face
    /// resists normal flow and generates no pressure force in tangential flow.
    /// The faces share the ledger's declared Cd; their areas must be supplied
    /// for that coefficient convention, not inferred from the scalar area.
    OrthogonalFaces {
        /// Orthonormal unit face normals, expressed in body FRD axes.
        normals_body: [[f64; 3]; 3],
        /// Nonnegative projected face areas [m²], at least one positive.
        areas_m2: [f64; 3],
    },
}

/// A component's original ledger identity and its body-fixed load geometry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ComponentDragGeometry {
    /// Must match the corresponding ledger item, in original ledger order.
    pub component_id: &'static str,
    /// Application point in the aircraft's body FRD geometry frame [m].
    pub point_body_m: [f64; 3],
    /// Explicit area/orientation model; no empirical values are manufactured.
    pub shape: ComponentDragShape,
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn invalid(message: &'static str) -> Refusal {
    Refusal {
        code: "drag-geometry-invalid",
        message: message.into(),
        ranked_repairs: vec![
            "supply one finite application point and admitted area model per ledger component, in ledger order".into(),
        ],
    }
}

impl DragLedger {
    /// Admit a complete component layout before any incident-field sampling.
    /// No component can be omitted, duplicated, or rebound by a permutation.
    pub fn admit_geometry(&self, geometry: &[ComponentDragGeometry]) -> Result<(), Refusal> {
        self.admit()?;
        if geometry.len() != self.components.len() {
            return Err(invalid(
                "component drag geometry must cover the complete ledger",
            ));
        }
        for (component, geometry) in self.components.iter().zip(geometry) {
            if component.component_id != geometry.component_id
                || !geometry.point_body_m.iter().all(|v| v.is_finite())
            {
                return Err(invalid(
                    "component drag identity, order or application point is invalid",
                ));
            }
            if let ComponentDragShape::OrthogonalFaces {
                normals_body,
                areas_m2,
            } = geometry.shape
            {
                if !areas_m2.iter().all(|a| a.is_finite() && *a >= 0.0)
                    || !areas_m2.iter().any(|a| *a > 0.0)
                    || !normals_body.iter().flatten().all(|v| v.is_finite())
                {
                    return Err(invalid(
                        "directional drag needs finite nonnegative areas and finite normals",
                    ));
                }
                for i in 0..3 {
                    if (dot(normals_body[i], normals_body[i]) - 1.0).abs() > 1e-12
                        || (0..i).any(|j| dot(normals_body[i], normals_body[j]).abs() > 1e-12)
                    {
                        return Err(invalid("directional drag face normals must be orthonormal"));
                    }
                }
            }
        }
        Ok(())
    }

    /// Sum actual component loads about origin_body_m, without allocating
    /// reports or hashing in the success path. relative_velocity_at supplies
    /// the full body-point velocity relative to local air, INCLUDING rigid
    /// angular motion exactly once. Sampling refusals propagate unchanged.
    ///
    /// The explicit unresolved-interference fraction scales each component's
    /// load and loss at its own point. It remains an estimated proportional
    /// correction, not resolved flow between components. drag_n is the sum
    /// of local wind-axis drag projections; it need not equal the norm of the
    /// resultant force when component flows differ or pressure creates sideforce.
    ///
    /// # Errors
    /// Ledger/geometry admission and sample refusals; air-state-invalid for
    /// nonfinite flow/origin or nonpositive density; drag-wrench-overflow for
    /// nonrepresentable force, moment or relative power.
    pub fn wrench_with_geometry(
        &self,
        rho_kg_m3: f64,
        geometry: &[ComponentDragGeometry],
        origin_body_m: [f64; 3],
        relative_velocity_at: &dyn Fn([f64; 3]) -> Result<[f64; 3], Refusal>,
    ) -> Result<DragWrench, Refusal> {
        self.admit_geometry(geometry)?;
        if !rho_kg_m3.is_finite()
            || rho_kg_m3 <= 0.0
            || !origin_body_m.iter().all(|v| v.is_finite())
        {
            return Err(air_invalid());
        }
        let mut result = DragWrench::default();
        for (component, geometry) in self.components.iter().zip(geometry) {
            let velocity = relative_velocity_at(geometry.point_body_m)?;
            if !velocity.iter().all(|v| v.is_finite()) {
                return Err(air_invalid());
            }
            let speed = velocity[0].hypot(velocity[1]).hypot(velocity[2]);
            if !speed.is_finite() {
                return Err(overflow());
            }
            if speed == 0.0 {
                continue;
            }
            let (force, drag, power) = match geometry.shape {
                ComponentDragShape::Isotropic => {
                    let q = 0.5 * rho_kg_m3 * speed * speed;
                    let drag = q * component.area_m2 * component.cd;
                    (velocity.map(|v| -drag * (v / speed)), drag, drag * speed)
                }
                ComponentDragShape::OrthogonalFaces {
                    normals_body,
                    areas_m2,
                } => {
                    let mut force = [0.0; 3];
                    let mut power = 0.0;
                    for (normal, area) in normals_body.into_iter().zip(areas_m2) {
                        if area == 0.0 {
                            continue;
                        }
                        let vn = dot(velocity, normal);
                        let q = 0.5 * rho_kg_m3 * vn.abs() * vn.abs();
                        let magnitude = q * area * component.cd;
                        let signed_force = -magnitude.copysign(vn);
                        for i in 0..3 {
                            force[i] += signed_force * normal[i];
                        }
                        power += magnitude * vn.abs();
                    }
                    (force, power / speed, power)
                }
            };
            let arm =
                std::array::from_fn::<_, 3, _>(|i| geometry.point_body_m[i] - origin_body_m[i]);
            let moment = [
                arm[1] * force[2] - arm[2] * force[1],
                arm[2] * force[0] - arm[0] * force[2],
                arm[0] * force[1] - arm[1] * force[0],
            ];
            for i in 0..3 {
                result.force_n[i] += force[i];
                result.moment_nm[i] += moment[i];
            }
            result.drag_n += drag;
            result.dissipated_power_w += power;
        }
        // Apply the same explicit fractional allowance to the resolved loads
        // and work. This is equivalent to scaling every component at its point.
        let fraction = self.unresolved_interference_frac;
        for value in result.force_n.iter_mut().chain(result.moment_nm.iter_mut()) {
            *value += *value * fraction;
        }
        result.drag_n += result.drag_n * fraction;
        result.dissipated_power_w += result.dissipated_power_w * fraction;
        if !result
            .force_n
            .iter()
            .chain(result.moment_nm.iter())
            .all(|v| v.is_finite())
            || !result.drag_n.is_finite()
            || !result.dissipated_power_w.is_finite()
        {
            return Err(overflow());
        }
        Ok(result)
    }

    /// Separate spatial model identity. The original scalar ledger digest is
    /// unchanged; callers selecting component geometry must include this one.
    #[must_use]
    pub fn spatial_digest(&self, geometry: &[ComponentDragGeometry]) -> String {
        let mut bytes = self.digest().into_bytes();
        bytes.extend_from_slice(&(geometry.len() as u64).to_le_bytes());
        for geometry in geometry {
            bytes.extend_from_slice(&(geometry.component_id.len() as u64).to_le_bytes());
            bytes.extend_from_slice(geometry.component_id.as_bytes());
            for value in geometry.point_body_m {
                bytes.extend_from_slice(&value.to_bits().to_le_bytes());
            }
            match geometry.shape {
                ComponentDragShape::Isotropic => bytes.push(0),
                ComponentDragShape::OrthogonalFaces {
                    normals_body,
                    areas_m2,
                } => {
                    bytes.push(1);
                    for value in normals_body.into_iter().flatten().chain(areas_m2) {
                        bytes.extend_from_slice(&value.to_bits().to_le_bytes());
                    }
                }
            }
        }
        hash_domain("org.frankensim.fs-flyer.component-drag.v1", &bytes).to_hex()
    }
}

fn air_invalid() -> Refusal {
    Refusal {
        code: "air-state-invalid",
        message:
            "component drag needs positive finite density and finite local relative motion/origin"
                .into(),
        ranked_repairs: vec!["supply each material-point velocity relative to local air".into()],
    }
}

fn overflow() -> Refusal {
    Refusal {
        code: "drag-wrench-overflow",
        message: "component drag force, moment or power is not representable".into(),
        ranked_repairs: vec![
            "use representable geometry, density and local relative motion".into(),
        ],
    }
}
