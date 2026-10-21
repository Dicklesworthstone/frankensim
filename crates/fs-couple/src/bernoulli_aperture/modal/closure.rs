//! Physical nodal lay/blade contact and slit flow on the retained DKT basis.
//!
//! Every original displacement column participates in local gap and contact
//! work, including modes with zero mean slit motion. The same interface owns
//! one contact potential and one set of slit edges for one or two blades.

use super::material::MaterialModalPlate;
use super::{AcousticRealizeError, BladeState, MAX_MODES, invalid, nonlinear};
use crate::bernoulli_aperture::plate::closure::{
    PlateClosureProbe, PlateClosureSpec, contact_quadrature, positive_linear_mean,
};
use crate::bernoulli_aperture::viscous::{ViscousProfileFlow, ViscousSlit};
use fs_dcontact::{ContactStorage, Obstacle, OpeningContactStep};
use fs_phs::Storage;

/// Physical point reactions pulled back into the original modal coordinates.
#[derive(Clone, Copy, Debug)]
pub(crate) struct ClosureResponse {
    pub forces_n: [f64; MAX_MODES],
    pub normal_force_n: f64,
    pub dissipated_power_w: f64,
}

/// One matching planar interface attached to the original material modes.
/// Contact uses triangle area quadrature; flow integrates the positive P1 gap
/// along each slit edge. Pressure faces and swept-volume projections remain
/// those of the source plates. This is a fixed small-slope basis, not changing
/// pressure footprints, nonmatching mesh contact or resolved moving-wall flow.
pub struct ModalPlateClosure {
    spec: PlateClosureSpec,
    shapes: Vec<Vec<f64>>,
    edges: Vec<([usize; 2], f64)>,
    nodes: Vec<usize>,
    points: Vec<Obstacle>,
    storage: ContactStorage,
    modes: usize,
    width_m: f64,
}
pub(super) struct ModalContactStep<'a> {
    closure: &'a ModalPlateClosure,
    points: Vec<OpeningContactStep>,
}

struct ZeroStorage;
impl Storage for ZeroStorage {
    fn hamiltonian(&self, _: &[f64]) -> f64 {
        0.0
    }
    fn gradient(&self, _: &[f64], out: &mut [f64]) {
        out.fill(0.0);
    }
}
impl ModalPlateClosure {
    pub(super) fn from_plates(
        plates: &[MaterialModalPlate],
        spec: PlateClosureSpec,
    ) -> Result<Self, AcousticRealizeError> {
        let first = plates
            .first()
            .ok_or_else(|| invalid("modal closure needs original material plates"))?;
        let mesh = &first.source().chart().mesh;
        let edges = |p: &MaterialModalPlate| {
            p.options()
                .slit_edges
                .iter()
                .map(|&[a, b]| (a.min(b), a.max(b)))
                .collect::<std::collections::BTreeSet<_>>()
        };
        if plates.len() > 2
            || plates.iter().any(|p| {
                p.source().chart().mesh.nodes != mesh.nodes
                    || p.source().chart().mesh.tris != mesh.tris
                    || edges(p) != edges(first)
            })
        {
            return Err(invalid(
                "modal contact requires matching original planar meshes and slit edges",
            ));
        }
        let modes = plates.last().expect("nonempty plates").mode_range().end;
        if modes == 0 || modes > MAX_MODES {
            return Err(invalid(
                "modal contact exceeds its original coordinate budget",
            ));
        }
        let mut shapes = vec![vec![0.0; modes]; mesh.node_count()];
        for plate in plates {
            for (row, shape) in shapes.iter_mut().zip(plate.displacement_shapes()) {
                row[plate.mode_range()].copy_from_slice(shape);
            }
        }
        Self::build(mesh, first.options().slit_edges.as_slice(), shapes, spec)
    }

    fn build(
        mesh: &fs_plate::PlateMesh,
        slit: &[[usize; 2]],
        shapes: Vec<Vec<f64>>,
        spec: PlateClosureSpec,
    ) -> Result<Self, AcousticRealizeError> {
        let modes = shapes.first().map_or(0, Vec::len);
        if modes == 0
            || modes > MAX_MODES
            || shapes.len() != mesh.node_count()
            || shapes
                .iter()
                .any(|s| s.len() != modes || s.iter().any(|x| !x.is_finite()))
        {
            return Err(invalid(
                "modal contact requires complete finite original displacement columns",
            ));
        }
        let (weights, nodes) = contact_quadrature(mesh, &spec)?;
        let mut rows = Vec::with_capacity(nodes.len() * modes);
        let mut gaps = Vec::with_capacity(nodes.len());
        let mut areas = Vec::with_capacity(nodes.len());
        let mut points = Vec::with_capacity(nodes.len());
        for &node in &nodes {
            let gap = spec.nodal_rest_gap_m[node];
            let law = Obstacle::new(
                vec![-1.0],
                1,
                1,
                vec![gap],
                vec![weights[node]],
                spec.stiffness_pa_per_m_alpha,
                spec.alpha,
                spec.provenance.clone(),
            )
            .and_then(|p| p.with_internal_loss(spec.internal_loss_s_per_m))
            .map_err(nonlinear)?;
            points.push(law);
            rows.extend(shapes[node].iter().map(|x| -x));
            gaps.push(gap);
            areas.push(weights[node]);
        }
        let law = Obstacle::new(
            rows,
            nodes.len(),
            modes,
            gaps,
            areas,
            spec.stiffness_pa_per_m_alpha,
            spec.alpha,
            spec.provenance.clone(),
        )
        .and_then(|p| p.with_internal_loss(spec.internal_loss_s_per_m))
        .map_err(nonlinear)?;
        let storage =
            ContactStorage::new(Box::new(ZeroStorage), modes, vec![law]).map_err(nonlinear)?;
        let edges: Vec<_> = slit
            .iter()
            .map(|&edge| {
                let [a, b] = edge.map(|node| mesh.nodes[node]);
                (edge, (a.0 - b.0).hypot(a.1 - b.1))
            })
            .collect();
        let width_m: f64 = edges.iter().map(|(_, length)| length).sum();
        if !width_m.is_finite() || width_m <= 0.0 {
            return Err(invalid(
                "modal closure needs a finite positive physical slit width",
            ));
        }
        Ok(Self {
            spec,
            shapes,
            edges,
            nodes,
            points,
            storage,
            modes,
            width_m,
        })
    }

    /// Original total rest clearances, physical contact patches and pair law.
    pub const fn spec(&self) -> &PlateClosureSpec {
        &self.spec
    }
    /// Original interface node for each physical contact quadrature point.
    pub fn contact_nodes(&self) -> &[usize] {
        &self.nodes
    }
    /// Original signed nodal displacement columns in complete modal order.
    pub fn displacement_shapes(&self) -> &[Vec<f64>] {
        &self.shapes
    }
    fn validate_coordinates(&self, q: &[f64]) -> Result<(), AcousticRealizeError> {
        if q.len() != self.modes || q.iter().any(|x| !x.is_finite()) {
            return Err(invalid(
                "modal closure needs all finite retained coordinates",
            ));
        }
        Ok(())
    }
    fn gap(&self, node: usize, q: &[f64]) -> f64 {
        self.spec.nodal_rest_gap_m[node]
            + self.shapes[node]
                .iter()
                .zip(q)
                .map(|(s, x)| s * x)
                .sum::<f64>()
    }
    /// Local gap includes every mode on every supplied blade, once.
    pub fn nodal_gap_m(&self, node: usize, q: &[f64]) -> Result<f64, AcousticRealizeError> {
        self.validate_coordinates(q)?;
        if node >= self.shapes.len() {
            return Err(invalid("modal closure node is outside its source mesh"));
        }
        let gap = self.gap(node, q);
        if !gap.is_finite() {
            return Err(invalid("modal local gap overflowed"));
        }
        Ok(gap)
    }
    /// Integrate max(gap(s),0), including the triangular open part of a partially
    /// closed edge. Clipping the mean or endpoint values is not this integral.
    pub fn open_area_m2(&self, q: &[f64]) -> Result<f64, AcousticRealizeError> {
        self.validate_coordinates(q)?;
        let area: f64 = self
            .edges
            .iter()
            .map(|&(edge, length)| {
                length * positive_linear_mean(self.gap(edge[0], q), self.gap(edge[1], q))
            })
            .sum();
        if !area.is_finite() || area < 0.0 {
            return Err(invalid("modal slit area overflowed"));
        }
        Ok(area)
    }
    /// Width-mean signed clearance, distinct from the positive flow-equivalent gap.
    pub fn mean_gap_m(&self, q: &[f64]) -> Result<f64, AcousticRealizeError> {
        self.validate_coordinates(q)?;
        let gap: f64 = self
            .edges
            .iter()
            .map(|&(edge, length)| {
                length * f64::midpoint(self.gap(edge[0], q), self.gap(edge[1], q))
            })
            .sum::<f64>()
            / self.width_m;
        if !gap.is_finite() {
            return Err(invalid("modal mean clearance overflowed"));
        }
        Ok(gap)
    }
    /// Same analytic laminar/Bernoulli channel integral as the scalar and paired
    /// plate owners, evaluated on the actual P1 gap. Material gas supplies mu/rho.
    pub fn viscous_flow(
        &self,
        law: ViscousSlit,
        q: &[f64],
        density_kg_m3: f64,
        pressure_drop_pa: f64,
    ) -> Result<ViscousProfileFlow, AcousticRealizeError> {
        self.validate_coordinates(q)?;
        let mut total = ViscousProfileFlow::default();
        for &(edge, length) in &self.edges {
            let local = law.linear_gap_flow(
                length,
                [self.gap(edge[0], q), self.gap(edge[1], q)],
                density_kg_m3,
                pressure_drop_pa,
            )?;
            total.flow_m3_s += local.flow_m3_s;
            total.flow_per_pressure_m3_s_pa += local.flow_per_pressure_m3_s_pa;
        }
        if !total.flow_m3_s.is_finite() || !total.flow_per_pressure_m3_s_pa.is_finite() {
            return Err(invalid("modal spatial slit flow overflowed"));
        }
        Ok(total)
    }
    /// Original fs-dcontact potential on every physical contact point, once.
    pub fn stored_energy_j(&self, q: &[f64]) -> Result<f64, AcousticRealizeError> {
        self.validate_coordinates(q)?;
        let mut state = [0.0; 2 * MAX_MODES];
        for (i, x) in q.iter().enumerate() {
            state[2 * i] = *x;
        }
        let energy = self.storage.hamiltonian(&state[..2 * self.modes]);
        if !energy.is_finite() || energy < 0.0 {
            return Err(invalid("modal contact storage overflowed"));
        }
        Ok(energy)
    }
    /// Physical area, contact activity and penetration from the same source map.
    pub fn probe(&self, q: &[f64]) -> Result<PlateClosureProbe, AcousticRealizeError> {
        let mut probe = PlateClosureProbe {
            open_area_m2: self.open_area_m2(q)?,
            active_lay_points: 0,
            active_lay_area_m2: 0.0,
            max_penetration_m: 0.0,
        };
        for (&node, point) in self.nodes.iter().zip(&self.points) {
            let gap = self.gap(node, q);
            if !gap.is_finite() {
                return Err(invalid("modal contact gap overflowed"));
            }
            if gap < 0.0 {
                probe.active_lay_points += 1;
                probe.active_lay_area_m2 += point.weights()[0];
                probe.max_penetration_m = probe.max_penetration_m.max(-gap);
            }
        }
        Ok(probe)
    }
    pub(super) fn validate_state(&self, states: &[BladeState]) -> Result<(), AcousticRealizeError> {
        let q = displacements(states);
        if self.probe(&q[..self.modes])?.max_penetration_m > self.spec.max_penetration_m {
            return Err(invalid(
                "modal interface exceeds its declared local penetration allowance",
            ));
        }
        self.stored_energy_j(&q[..self.modes])?;
        Ok(())
    }
    pub(super) fn energy(&self, states: &[BladeState]) -> f64 {
        let q = displacements(states);
        self.stored_energy_j(&q[..self.modes]).unwrap_or(f64::NAN)
    }
    pub(super) fn state_gap(&self, states: &[BladeState]) -> f64 {
        let q = displacements(states);
        self.mean_gap_m(&q[..self.modes]).unwrap_or(f64::NAN)
    }
    pub(super) fn contact_step(
        &self,
        before: &[f64],
    ) -> Result<ModalContactStep<'_>, AcousticRealizeError> {
        self.validate_coordinates(before)?;
        let points = self
            .points
            .iter()
            .zip(&self.nodes)
            .map(|(law, &node)| {
                let opening = self.shapes[node]
                    .iter()
                    .zip(before)
                    .map(|(s, x)| s * x)
                    .sum();
                OpeningContactStep::new(law, opening).map_err(nonlinear)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(ModalContactStep {
            closure: self,
            points,
        })
    }
    pub(super) fn response(
        &self,
        before: &[f64],
        after: &[f64],
        velocity: &[f64],
    ) -> Result<ClosureResponse, AcousticRealizeError> {
        self.contact_step(before)?.response(after, velocity)
    }
}
impl ModalContactStep<'_> {
    // The existing allocation-free scalar power-law secant is cached once per
    // point and mechanical sample. Newton probes retain that same old state;
    // they cannot advance a contact or Maxwell history.
    pub(super) fn response(
        &self,
        after: &[f64],
        velocity: &[f64],
    ) -> Result<ClosureResponse, AcousticRealizeError> {
        self.closure.validate_coordinates(after)?;
        self.closure.validate_coordinates(velocity)?;
        let mut result = ClosureResponse {
            forces_n: [0.0; MAX_MODES],
            normal_force_n: 0.0,
            dissipated_power_w: 0.0,
        };
        for (step, &node) in self.points.iter().zip(&self.closure.nodes) {
            let shape = &self.closure.shapes[node];
            let dot = |q: &[f64]| shape.iter().zip(q).map(|(s, x)| s * x).sum::<f64>();
            let (elastic, damping) = step.coefficients(dot(after)).map_err(nonlinear)?;
            let local_velocity = dot(velocity);
            let force = (elastic - damping * local_velocity).max(0.0);
            let power = (elastic - force) * local_velocity;
            if !force.is_finite() || !power.is_finite() || power < 0.0 {
                return Err(invalid("modal local contact power is unresolved"));
            }
            for (generalized, s) in result.forces_n.iter_mut().zip(shape) {
                *generalized += s * force;
            }
            result.normal_force_n += force;
            result.dissipated_power_w += power;
        }
        if result
            .forces_n
            .iter()
            .chain([&result.normal_force_n, &result.dissipated_power_w])
            .any(|x| !x.is_finite())
        {
            return Err(invalid("modal pointwise contact work is unresolved"));
        }
        Ok(result)
    }
}

pub(super) fn displacements(states: &[BladeState]) -> [f64; MAX_MODES] {
    let mut q = [0.0; MAX_MODES];
    for (out, state) in q.iter_mut().zip(states) {
        *out = state.displacement_m;
    }
    q
}

#[cfg(test)]
mod tests {
    use super::*;
    fn fixture() -> ModalPlateClosure {
        let mesh = fs_plate::PlateMesh::from_unstructured(
            vec![(0.0, 0.0), (0.01, 0.0), (0.0, 0.02)],
            vec![[0, 1, 2]],
        )
        .unwrap();
        ModalPlateClosure::build(
            &mesh,
            &[[0, 1]],
            vec![
                vec![1.0, 1.0, 0.3],
                vec![1.0, -1.0, 1.2],
                vec![0.2, -0.4, 0.7],
            ],
            PlateClosureSpec {
                nodal_rest_gap_m: vec![0.0001; 3],
                lay_triangles: vec![0],
                stiffness_pa_per_m_alpha: 3e9,
                alpha: 2.0,
                internal_loss_s_per_m: 2.0,
                provenance: "explicit numerical multimode interface fixture".into(),
                max_penetration_m: 0.01,
            },
        )
        .unwrap()
    }
    #[test]
    fn zero_mean_mode_opens_part_of_an_otherwise_closed_slit() {
        let c = fixture();
        let q = [-0.0001, 0.0002, 0.0];
        assert!(c.mean_gap_m(&q).unwrap().abs() < 1e-20);
        assert!((c.open_area_m2(&q).unwrap() - 5e-7).abs() < 1e-20);
        let response = c.response(&q, &q, &[0.0; 3]).unwrap();
        assert!(response.normal_force_n > 0.0);
        assert!(response.forces_n[1] < 0.0); // zero-mean mode still touches the lay
        let law = ViscousSlit::new(1.81e-5, 0.001).unwrap();
        assert!(c.viscous_flow(law, &q, 1.2, 100.0).unwrap().flow_m3_s > 0.0);
    }
    #[test]
    fn physical_contact_potential_and_pointwise_secants_close_all_modal_work() {
        let c = fixture();
        let dt = 0.002;
        for a in [[0.0; 3], [-0.001, 0.0002, 0.0004], [0.002, -0.003, 0.001]] {
            for b in [[0.0; 3], [-0.002, 0.001, -0.0003], [0.003, -0.002, 0.0002]] {
                let v = std::array::from_fn::<_, 3, _>(|i| (b[i] - a[i]) / dt);
                let f = c.response(&a, &b, &v).unwrap();
                let work: f64 = (0..3).map(|i| f.forces_n[i] * (b[i] - a[i])).sum();
                let e0 = c.stored_energy_j(&a).unwrap();
                let e1 = c.stored_energy_j(&b).unwrap();
                let loss = dt * f.dissipated_power_w;
                assert!(loss >= 0.0);
                assert!(
                    (e1 - e0 + work + loss).abs()
                        < 3e-12 * (e0 + e1 + work.abs() + loss).max(1e-30)
                );
                // Independent physical potential, using the triangle area / 3.
                let expected: f64 = (0..3)
                    .map(|node| {
                        1e-4 / 3.0 * 3e9 * (-c.nodal_gap_m(node, &b).unwrap()).max(0.0).powi(3)
                            / 3.0
                    })
                    .sum();
                assert!((e1 - expected).abs() < 2e-12 * expected.max(1e-30));
            }
        }
    }
}
