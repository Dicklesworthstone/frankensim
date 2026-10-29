//! Source-bound segment, tensile-link and gravity assembly on one native owner.
//! Full source admission precedes assembly. No derivative probes leave a
//! material's declared strain domain. Material conditions and reference mass
//! remain fixed throughout an equilibrium solve.

use super::SpatialEnergy;
use super::super::SpatialRodError;
use crate::rod::Rod;
use crate::rod::derivatives::analytic::{Assembly, Kinematics, Linearization,
    PairResponse, array, vector, vadd, vsub, join};
use fs_ga::Vec3;

impl SpatialEnergy<'_> {
    pub(super) fn assemble_analytic(&self,rod:&Rod,free:&[usize],tangent:bool,
        checkpoint:&impl Fn()->Result<(),SpatialRodError>) -> Result<Linearization,SpatialRodError> {
        checkpoint()?;
        // Includes every layer and link, including completely prescribed regions.
        // A derivative cannot launder a source-domain refusal into a zero force.
        self.owner.equilibrium_energy_for(rod)?;
        let mut assembly=Assembly::new(rod.positions.len(),free,tangent)?;
        for (segment,s) in self.owner.stiffness.iter().enumerate() {
            checkpoint()?;
            let kinematics=Kinematics::new(rod,segment)?;
            let (g,k)=self.owner.reference_relative_strains(segment,kinematics.strains())?;
            let (force,moment)=self.owner.resultants_from_strains(segment,g,k)?;
            let response=kinematics.linearize(vector(force),vector(moment),
                [s.ea,s.ga,s.ga],[s.gj,s.ei,s.ei],tangent)?;
            assembly.add_pair([segment,segment+1],&response);
        }
        if let Some(forces)=self.owner.gravity_nodal_forces_n() {
            for (node,&force) in forces.iter().enumerate() {
                checkpoint()?;
                assembly.add_node(node,join(vector(force).scale(-1.0),Vec3::new(0.0,0.0,0.0)),[[0.0;6];6]);
            }
        }
        if let Some(fluid)=self.owner.submerged_fluid_load() {
            for (node,&force) in fluid.nodal_forces_n().iter().enumerate() {
                checkpoint()?;
                assembly.add_node(node,join(vector(force).scale(-1.0),Vec3::new(0.0,0.0,0.0)),[[0.0;6];6]);
            }
        }
        for (index,link) in self.owner.tension_links.iter().enumerate() {
            checkpoint()?;
            let state=link.state_for(rod)?;
            let attachments=link.attachments();
            let nodes=(*attachments).map(|a|a.node);
            let offsets=(*attachments).map(|a|vector(a.offset_body_m));
            let rotations=nodes.map(|node|rod.frames[node].matrix());
            let forces=state.force_n;
            let gradient=std::array::from_fn(|side|join(vector(forces[side]).scale(-1.0),
                vector(state.moment_body_nm[side]).scale(-1.0)));
            let matrix=if tangent {
                let mut matrix=[[0.0;12];12];
                // At the C1 energy's taut/slack kink choose the slack-sided
                // force derivative explicitly. Re-evaluate on EVERY state;
                // slack links are never removed from the interaction graph.
                if state.tensile_strain>0.0 {
                    let f=vector(forces[0]);
                    let delta=vadd(vsub(vector(rod.positions[nodes[1]]),vector(rod.positions[nodes[0]])),
                        vsub(rotations[1].apply(offsets[1]),rotations[0].apply(offsets[0])));
                    let direction=vector(array(delta).map(|v| v/state.length_m));
                    let stiffness=link.section().axial_stiffness_n()
                        * f64::from(link.parallel_strands())/link.natural_length_m();
                    let transverse=state.tension_n/state.length_m;
                    if !stiffness.is_finite() || !transverse.is_finite() {
                        return Err(SpatialRodError::Invalid("tension link tangent is unrepresentable"));
                    }
                    for col in 0..12 {
                        let side=col/6;let component=col%6;
                        let mut unit=[0.0;3];unit[component%3]=1.0;let unit=vector(unit);
                        let dx=if component<3 {unit} else {rotations[side].apply(unit.cross(offsets[side]))};
                        let dd=dx.scale(if side==0 {-1.0} else {1.0});
                        let df=vadd(dd.scale(transverse),direction.scale((stiffness-transverse)*direction.dot(dd)));
                        for row_side in 0..2 {
                            let mut body=rotations[row_side].transpose().apply(df);
                            if component>=3 && row_side==side {
                                body=vsub(body,unit.cross(rotations[row_side].transpose().apply(f)));
                            }
                            let sign=if row_side==0 {-1.0} else {1.0};
                            let values=join(df.scale(sign),offsets[row_side].cross(body).scale(sign));
                            for row in 0..6 {matrix[6*row_side+row][col]=values[row];}
                        }
                    }
                }
                Some(matrix)
            } else {None};
            let response=PairResponse {gradient,tangent:matrix};
            response.check()?;
            assembly.add_pair(nodes,&response);
            // Centerline weight and buoyancy include the link inventory. Only
            // rotating-offset moments (and their derivatives) remain here.
            if let Some(force)=self.owner.link_body_load_force_n(index) {
                let force=vector(force);
                for side in 0..2 {
                    let body=rotations[side].transpose().apply(force);
                    let moment=offsets[side].cross(body).scale(-1.0);
                    let mut local=[[0.0;6];6];
                    if tangent {
                        for col in 0..3 {
                            let mut unit=[0.0;3];unit[col]=1.0;
                            let derivative=offsets[side].cross(vector(unit).cross(body));
                            let entries=array(derivative);
                            for row in 0..3 {local[row+3][col+3]=entries[row];}
                        }
                    }
                    assembly.add_node(nodes[side],join(Vec3::new(0.0,0.0,0.0),moment),local);
                }
            }
        }
        checkpoint()?;
        Ok(assembly.finish()?)
    }
}
