//! Finite acoustic boundary from the SAME prepared physical board.
//! No duplicate eigensolve, nearest-facet projection, smoothing, box or lid.
use super::{Boundary, Specification, MotionSurface, SpherePanels, MAX_PANELS,
    closed_components, norm, sub};
use super::super::board_geometry::motion::skin::{Skin, LABEL};

impl Specification {
    /// A section-derived skin already uses the structural metre frame and
    /// contains only one moving part. Refuse all ignored transforms/labels.
    pub fn require_board_skin(&self) -> Result<(), String> {
        if self.scale_m != 1. || self.origin_obj != [0.;3]
            || self.rules.len() != 1 || self.rules.get(LABEL) != Some(&true) {
            return Err("board-skin requires obj-scale-m,1; obj-origin,0,0,0; and only moving,soundboard_skin (no missing rigid components)".into());
        }
        Ok(())
    }
}
impl Boundary {
    /// Uses the actual equilibrium midsurface and the source's per-element
    /// sections, with exact source-site motion. The BEM work limit is unchanged;
    /// known embeddings need no all-pairs facet search or its search budget.
    pub fn from_board_skin(source: &str, spec: &Specification, motion: &MotionSurface, continuous: bool)
        -> Result<(Self, String), String> {
        spec.require_board_skin()?;
        let skin = if continuous { Skin::continuous_from_source(motion, source, spec.offset_m, MAX_PANELS)? }
            else { Skin::from_source(motion, source, spec.offset_m, MAX_PANELS)? };
        let triangles = skin.panel_triangles();
        if triangles.iter().flatten().flatten().any(|v| !v.is_finite() || v.abs() > 100.) {
            return Err("section skin exceeds the finite 100 m acoustic-coordinate limit".into());
        }
        // Keep the existing geometric (not just indexed) closure admission.
        let components = closed_components(&triangles)?;
        let center = std::array::from_fn(|c| {
            let lo = skin.vertices.iter().map(|p| p[c]).fold(f64::INFINITY, f64::min);
            let hi = skin.vertices.iter().map(|p| p[c]).fold(f64::NEG_INFINITY, f64::max);
            f64::midpoint(lo, hi)
        });
        let radius = skin.vertices.iter().map(|&p| norm(sub(p, center))).fold(0., f64::max);
        let surface = SpherePanels::from_triangles(triangles).map_err(|e| e.to_string())?;
        let weights = skin.normal_weights(motion, surface.normals())?;
        let image = if continuous { "explicit volume-preserving continuous" } else { "section-column" };
        let report = format!("{image} board skin: {} panels, {} closed components, section/skin volume {:.12e}/{:.12e} m3, maximum facet mean-thickness change {:.3e} m; 1 nm height grid; actual equilibrium coordinates; no cabinet/lid or extra rib/bridge acoustic solids",
            skin.triangles.len(), components, skin.section_volume_m3, skin.volume_m3, skin.maximum_thickness_change_m);
        Ok((Self {surface, weights, center, radius, components}, report))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn fixture() -> (MotionSurface, String, Specification) {
        let mesh=fs_plate::ShellMesh::new(vec![[0.,0.,0.],[0.1,0.,0.],[0.1,0.1,0.],[0.,0.1,0.]],
            vec![[0,1,2],[0,2,3]]).unwrap();
        let motion=MotionSurface::new(mesh,vec![vec![[0.,0.,1.,0.,0.,0.];4]]).unwrap();
        let source="frankensim-board-geometry-si-v1\ntriangle,0,0,1,2,0.02,450,1e10,8e8,0.3,6e8,0\ntriangle,1,0,2,3,0.02,450,1e10,8e8,0.3,6e8,0\n".into();
        let spec=Specification::read(&super::super::tests::specification().replace("moving,skin","moving,soundboard_skin")).unwrap();
        (motion,source,spec)
    }
    #[test]
    fn section_boundary_uses_existing_geometry_and_bem_with_signed_two_sided_motion() {
        let (motion,source,spec)=fixture();
        let (boundary,report)=Boundary::from_board_skin(&source,&spec,&motion,false).unwrap();
        assert!(report.contains("section-column"));assert_eq!(boundary.components,1);
        let skin=Skin::from_source(&motion,&source,0.02,MAX_PANELS).unwrap();
        let reimported=Boundary::from_obj(&skin.obj(),&spec,&motion).unwrap();
        assert_eq!(boundary.surface.areas(),reimported.surface.areas());
        for (a,b) in boundary.weights.iter().flatten().zip(reimported.weights.iter().flatten()) {
            assert!((a-b).abs()<1e-14);
        }
        let receivers=[[0.05,0.05,1.],[0.05,0.05,-1.]];
        let sample=boundary.sample_grid(&[std::f64::consts::TAU*100.],&receivers,spec.medium,6.).unwrap();
        let a=sample.values[0][0][0];let b=sample.values[1][0][0];
        assert!(a.abs()>1e-10);assert!((a+b).abs()<0.02*a.abs());
        let flux:f64=boundary.weights[0].iter().zip(boundary.surface.areas()).map(|(w,a)|w*a).sum();
        assert!(flux.abs()<1e-14);
    }
    #[test]
    fn generated_skin_cannot_ignore_transforms_parts_or_offset_limits() {
        let (motion,source,mut spec)=fixture();
        spec.scale_m=0.001;assert!(Boundary::from_board_skin(&source,&spec,&motion,false).is_err());
        spec.scale_m=1.;spec.origin_obj=[1.,0.,0.];assert!(spec.require_board_skin().is_err());
        spec.origin_obj=[0.;3];spec.rules.insert("lid".into(),false);assert!(spec.require_board_skin().is_err());
        spec.rules.remove("lid");spec.offset_m=0.001;
        assert!(Boundary::from_board_skin(&source,&spec,&motion,false).is_err());
    }
}
