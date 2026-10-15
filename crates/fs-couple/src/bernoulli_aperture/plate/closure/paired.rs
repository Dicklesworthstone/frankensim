//! A common planar contact surface between independently moving source blades.
//!
//! The two original charts must have identical planar nodes, triangles and slit
//! edges: this is an explicit matching-interface image, not inferred proximity
//! between unrelated meshes. Material sections, supports and modal shapes may
//! differ. Local gap is g_i + phi_0i*x_0 + phi_1i*x_1. One interface supplies
//! the slit area and contact work; neither area nor normal force is doubled.

use super::{PlateApertureReduction, PlateClosure, PlateClosureProbe, PlateClosureSpec,
    AcousticRealizeError, Obstacle, invalid, positive_linear_mean};
use crate::unilateral_contact::distributed::ApertureContactStep;
use fs_dcontact::ContactStorage;
use fs_phs::Storage;

/// Same-time contact pullback into each original blade coordinate.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PairedClosureResponse {
    /// Conservative generalized forces on the two opening-positive blades [N].
    pub elastic_forces_n: [f64; 2],
    /// Generalized forces after LOCAL nonadhesive unloading clamps [N].
    pub forces_n: [f64; 2],
    /// Sum of physical normal quadrature forces [N], not the modal force sum.
    pub normal_force_n: f64,
    /// Pointwise irreversible power [W]; each interface point contributes once.
    pub dissipated_power_w: f64,
}

/// Source-bound matching-interface geometry and the original contact owner.
/// Numerical trials allocate nothing; construction allocates one scalar law per
/// contact point so every secant uses the existing stable aperture/contact code.
/// This fixed small-slope image does not model separation-dependent pressure
/// footprints, curved shells, friction between blades or nonmatching meshes.
pub struct PairedPlateClosure {
    reference: PlateClosure,
    shapes: Vec<[f64; 2]>,
    points: Vec<Obstacle>,
    storage: ContactStorage,
}
struct ZeroStorage;
impl Storage for ZeroStorage {
    fn hamiltonian(&self, _: &[f64]) -> f64 { 0.0 }
    fn gradient(&self, _: &[f64], out: &mut [f64]) { out.fill(0.0); }
}
fn contact_error(e: impl core::fmt::Display) -> AcousticRealizeError {
    AcousticRealizeError::Nonlinear(e.to_string())
}

impl PlateApertureReduction {
    /// Bind one supplied clearance/contact map to TWO original blade modes.
    /// The first chart defines interface coordinates; the second uses the same
    /// planar indexing with its own opening-positive out-of-plane direction.
    /// Different sections and supports are allowed; no shape is averaged.
    ///
    /// # Errors
    /// Nonmatching surfaces/slits, already-paired reductions, invalid contact
    /// geometry/law, or any refusal from the existing scalar closure compiler.
    pub fn compile_paired_closure(&self, other: &Self, spec: PlateClosureSpec)
        -> Result<PairedPlateClosure, AcousticRealizeError>
    {
        if self.moving_blade_count()!=1 || other.moving_blade_count()!=1
            || self.chart().mesh.nodes!=other.chart().mesh.nodes
            || self.chart().mesh.tris!=other.chart().mesh.tris {
            return Err(invalid("paired closure requires two single blades on the same explicit planar interface mesh"));
        }
        let edges = |p: &PlateApertureReduction| {
            p.options().slit_edges.iter().map(|&[a,b]| (a.min(b),a.max(b)))
                .collect::<std::collections::BTreeSet<_>>()
        };
        if edges(self)!=edges(other) {
            return Err(invalid("paired closure requires the same physical slit edges on both blades"));
        }
        // Reuse the original triangle areas, contact-node ordering, slit edges,
        // and contact-law admission. This is not a second geometry quadrature.
        let reference=self.compile_closure(spec)?;
        let shapes=reference.shape.iter().zip(other.shape_per_opening())
            .map(|(&a,b)| [a,b[0]]).collect();
        PairedPlateClosure::build(reference,shapes)
    }
}

impl PairedPlateClosure {
    fn build(reference: PlateClosure, shapes: Vec<[f64;2]>) -> Result<Self,AcousticRealizeError> {
        if shapes.len()!=reference.shape.len() || shapes.iter().flatten().any(|x| !x.is_finite()) {
            return Err(invalid("paired closure requires complete finite original blade shapes"));
        }
        let mut points=Vec::with_capacity(reference.nodes.len());
        let mut rows=Vec::with_capacity(2*reference.nodes.len());
        let mut gaps=Vec::with_capacity(reference.nodes.len());
        for (row,&node) in reference.nodes.iter().enumerate() {
            let gap=reference.spec.nodal_rest_gap_m[node];
            let weight=reference.law.weights()[row];
            let point=Obstacle::new(vec![-1.0],1,1,vec![gap],vec![weight],
                reference.law.stiffness(),reference.law.alpha(),reference.spec.provenance.clone())
                .and_then(|p|p.with_internal_loss(reference.law.internal_loss())).map_err(contact_error)?;
            points.push(point);
            rows.extend(shapes[node].map(|x| -x));
            gaps.push(gap);
        }
        let law=Obstacle::new(rows,reference.nodes.len(),2,gaps,reference.law.weights().to_vec(),
            reference.law.stiffness(),reference.law.alpha(),reference.spec.provenance.clone())
            .and_then(|p|p.with_internal_loss(reference.law.internal_loss())).map_err(contact_error)?;
        let storage=ContactStorage::new(Box::new(ZeroStorage),2,vec![law]).map_err(contact_error)?;
        Ok(Self {reference,shapes,points,storage})
    }
    /// Unmodified common clearances, interface law and penetration allowance.
    #[must_use]
    pub fn spec(&self) -> &PlateClosureSpec { self.reference.spec() }
    /// Original mesh node corresponding to each contact quadrature point.
    #[must_use]
    pub fn contact_nodes(&self) -> &[usize] { self.reference.contact_nodes() }
    /// Independent local displacement shapes in the SAME interface indexing.
    #[must_use]
    pub fn shapes(&self) -> &[[f64;2]] { &self.shapes }
    /// Original fs-dcontact potential, counted once, in two displacement coordinates.
    pub fn stored_energy_j(&self, displacement: [f64;2]) -> Result<f64,AcousticRealizeError> {
        if displacement.iter().any(|x| !x.is_finite()) { return Err(invalid("paired contact needs finite blade displacements")); }
        let energy=self.storage.hamiltonian(&[displacement[0],displacement[1],0.0,0.0]);
        if !energy.is_finite() || energy<0.0 { return Err(invalid("paired contact storage overflowed")); }
        Ok(energy)
    }
    /// Local physical clearance; displacements are from each blade's rest shape.
    pub fn nodal_gap_m(&self,node:usize,displacement:[f64;2]) -> Result<f64,AcousticRealizeError> {
        let shape=self.shapes.get(node).ok_or_else(||invalid("paired closure node is outside its interface"))?;
        let gap=self.spec().nodal_rest_gap_m[node]+shape[0]*displacement[0]+shape[1]*displacement[1];
        if displacement.iter().any(|x| !x.is_finite()) || !gap.is_finite() {
            return Err(invalid("paired local gap overflowed"));
        }
        Ok(gap)
    }
    /// Integrate max(local gap,0) along the original P1 slit edges. A partially
    /// closed edge keeps its triangular open portion, not a clipped mean gap.
    pub fn open_area_m2(&self,displacement:[f64;2]) -> Result<f64,AcousticRealizeError> {
        let mut area=0.0;
        for &(edge,length) in &self.reference.edges {
            area+=length*positive_linear_mean(self.nodal_gap_m(edge[0],displacement)?,
                self.nodal_gap_m(edge[1],displacement)?);
        }
        if !area.is_finite() || area<0.0 { return Err(invalid("paired slit area overflowed")); }
        Ok(area)
    }
    /// Read physical slit area, contact activity and worst local penetration.
    pub fn probe(&self,displacement:[f64;2]) -> Result<PlateClosureProbe,AcousticRealizeError> {
        let mut out=PlateClosureProbe {open_area_m2:self.open_area_m2(displacement)?,
            active_lay_points:0,active_lay_area_m2:0.0,max_penetration_m:0.0};
        for (i,&node) in self.reference.nodes.iter().enumerate() {
            let gap=self.nodal_gap_m(node,displacement)?;
            if gap<0.0 {
                out.active_lay_points+=1;
                out.active_lay_area_m2+=self.reference.law.weights()[i];
                out.max_penetration_m=out.max_penetration_m.max(-gap);
            }
        }
        if !out.active_lay_area_m2.is_finite() { return Err(invalid("paired active contact area overflowed")); }
        Ok(out)
    }
    /// State admission, separate from unrestricted finite nonlinear trial queries.
    pub fn validate_displacements(&self,displacement:[f64;2]) -> Result<(),AcousticRealizeError> {
        if self.probe(displacement)?.max_penetration_m>self.spec().max_penetration_m {
            return Err(invalid("paired interface exceeds its local penetration allowance"));
        }
        self.stored_energy_j(displacement)?;
        Ok(())
    }
    /// Use the original stable scalar contact secant at EVERY physical point,
    /// then pull its reaction back into both blade coordinates. Local opening
    /// velocity is phi_0*v_0+phi_1*v_1; clamp unloading BEFORE summing forces.
    /// For discrete energy balance use v_i=(after_i-before_i)/dt. Trials change
    /// no history and allocate nothing; material evolution remains elsewhere.
    pub fn response(&self,before:[f64;2],after:[f64;2],velocity:[f64;2])
        -> Result<PairedClosureResponse,AcousticRealizeError>
    {
        if before.iter().chain(&after).chain(&velocity).any(|x| !x.is_finite()) {
            return Err(invalid("paired contact trial needs finite coordinates and velocities"));
        }
        let mut result=PairedClosureResponse::default();
        for (law,&node) in self.points.iter().zip(&self.reference.nodes) {
            let shape=self.shapes[node];
            let dot=|x:[f64;2]| shape[0]*x[0]+shape[1]*x[1];
            let local=ApertureContactStep::new(law,dot(before))
                .and_then(|p|p.response(dot(after),dot(velocity))).map_err(contact_error)?;
            for i in 0..2 {
                result.elastic_forces_n[i]+=shape[i]*local.elastic_force;
                result.forces_n[i]+=shape[i]*local.force;
            }
            result.normal_force_n+=local.force;
            result.dissipated_power_w+=local.dissipated_power;
        }
        if result.forces_n.iter().chain(&result.elastic_forces_n)
            .chain([&result.normal_force_n,&result.dissipated_power_w]).any(|x| !x.is_finite())
            || result.dissipated_power_w<0.0 {
            return Err(invalid("paired contact force or work overflowed"));
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // Numerical port fixture, not a claimed plate eigenmode or measured reed.
    fn fixture() -> PairedPlateClosure {
        let spec=PlateClosureSpec {nodal_rest_gap_m:vec![0.0001,-0.0001,0.0002],
            lay_triangles:vec![0],stiffness_pa_per_m_alpha:3e9,alpha:2.0,
            internal_loss_s_per_m:2.0,provenance:"authored paired port test".into(),max_penetration_m:0.1};
        let law=Obstacle::new(vec![-1.0,-0.5,0.2],3,1,vec![0.0;3],vec![0.0001;3],
            spec.stiffness_pa_per_m_alpha,spec.alpha,spec.provenance.clone()).unwrap()
            .with_internal_loss(spec.internal_loss_s_per_m).unwrap();
        let reference=PlateClosure {spec,law,rest_coordinate_m:0.0,shape:vec![1.0,0.5,-0.2],
            edges:vec![([0,1],0.01)],nodes:vec![0,1,2]};
        PairedPlateClosure::build(reference,vec![[1.0,0.3],[0.5,1.2],[-0.2,0.7]]).unwrap()
    }
    #[test]
    fn paired_pointwise_contact_balances_each_blades_virtual_work() {
        let c=fixture(); let dt=0.002;
        for a in [[0.0,0.0],[-0.001,0.0002],[0.002,-0.003]] {
            for b in [[0.0,0.0],[-0.002,0.001],[0.003,-0.002]] {
                let v=[(b[0]-a[0])/dt,(b[1]-a[1])/dt];
                let f=c.response(a,b,v).unwrap();
                let work=f.forces_n[0]*(b[0]-a[0])+f.forces_n[1]*(b[1]-a[1]);
                let e0=c.stored_energy_j(a).unwrap(); let e1=c.stored_energy_j(b).unwrap();
                let loss=dt*f.dissipated_power_w;
                assert!((e1-e0+work+loss).abs() < 2e-12*(e0+e1+work.abs()+loss).max(1e-30));
                assert!(loss>=0.0);
            }
        }
    }
    #[test]
    fn different_shapes_change_gap_even_when_total_tip_motion_cancels() {
        let c=fixture();
        let a=c.open_area_m2([0.0,0.0]).unwrap();
        let b=c.open_area_m2([0.0001,-0.0001]).unwrap();
        assert!((a-2.5e-7).abs()<1e-20);
        assert!((a-b).abs()>1e-8);
        let force=c.response([-0.001,-0.001],[-0.001,-0.001],[0.0,0.0]).unwrap();
        assert_ne!(force.forces_n[0],force.forces_n[1]);
        assert!(force.normal_force_n>0.0);
    }
    #[test]
    fn partial_edge_closure_is_not_clipped_endpoint_quadrature() {
        let c=fixture();
        // Rest edge is +0.1 mm to -0.1 mm. Mean gap is zero, area is not.
        assert!((c.open_area_m2([0.0,0.0]).unwrap()-0.01*0.0001/4.0).abs()<1e-20);
        let p=c.probe([0.0,0.0]).unwrap();
        assert_eq!(p.active_lay_points,1);
        assert_eq!(p.active_lay_area_m2,0.0001);
    }
    #[test]
    fn refused_geometry_trials_do_not_change_the_interface() {
        let c=fixture(); let before=c.probe([0.0,0.0]).unwrap();
        assert!(c.response([0.0,0.0],[f64::NAN,0.0],[0.0,0.0]).is_err());
        assert!(c.nodal_gap_m(99,[0.0,0.0]).is_err());
        assert!(c.validate_displacements([-1.0,-1.0]).is_err());
        assert_eq!(before,c.probe([0.0,0.0]).unwrap());
    }
}
