//! Regional bending forms in the retained ORIGINAL DKT eigenbasis.
use super::*;

impl RetainedPlateBasis {
    /// Original DKT curvature quadrature on the same complete physical modal
    /// coordinates as force, pressure, slit motion and material stiffness.
    pub(crate) fn bending_strain_points(
        &self, source: &ResolvedPlateChart, coordinate_scale: &[f64],
        max_projection_terms: usize, gate: &CancelGate,
    ) -> Result<Vec<crate::thin_plate::RetainedPlateBendingPoint>, AcousticRealizeError> {
        let n = self.mass_normalized_modes.len();
        let mesh = &source.chart().mesh;
        let poll = || if gate.is_requested() { Err(AcousticRealizeError::Cancelled) } else { Ok(()) };
        poll()?;
        if self.source_identity != source.specimen_identity() || n==0 || n>64
            || coordinate_scale.len()!=n || coordinate_scale.iter().any(|v| !v.is_finite() || *v<=0.0)
            || self.dof_map.len()!=3*mesh.nodes.len()
            || mesh.tris.len().checked_mul(90*n+512).is_none_or(|w| w>max_projection_terms) {
            return Err(refuse("bending strain projection requires the original complete bounded material basis"));
        }
        let mut points = Vec::with_capacity(3*mesh.tris.len());
        for (element, tri) in mesh.tris.iter().enumerate() {
            poll()?;
            let region = source.region_of_element(element);
            let section = &source.regions()[region].section;
            let (curvature, twice_area) = fs_plate::dkt_curvature_quadrature(
                &tri.map(|v| mesh.nodes[v].0), &tri.map(|v| mesh.nodes[v].1), element,
            ).map_err(map_plate)?;
            let volume_m3 = twice_area/6.0*section.thickness;
            let depth = section.thickness/fs_math::det::sqrt(12.0);
            if !volume_m3.is_finite() || volume_m3<=0.0 || !depth.is_finite() || depth<=0.0 {
                return Err(refuse("bending material quadrature needs positive original volume and thickness"));
            }
            let dofs: [Option<usize>;9] = core::array::from_fn(|i| self.dof_map[3*tri[i/3]+i%3]);
            let local: Vec<[f64;9]> = self.mass_normalized_modes.iter().zip(coordinate_scale)
                .map(|(mode,scale)| dofs.map(|d| d.map_or(0.0,|i| mode[i]*scale))).collect();
            for b in curvature {
                let columns: Vec<[f64;3]> = local.iter().map(|mode| core::array::from_fn(|i|
                    depth*(0..9).map(|j| b[i][j]*mode[j]).sum::<f64>())).collect();
                if columns.iter().flatten().any(|v| !v.is_finite()) {
                    return Err(refuse("bending material strain columns overflow"));
                }
                points.push(crate::thin_plate::RetainedPlateBendingPoint { region,columns,volume_m3 });
            }
        }
        poll()?;
        Ok(points)
    }

    /// Project each original region separately, including every cross-mode term.
    /// `coordinate_scale` maps the original mass-normalized coordinates to the
    /// consumer's physical coordinates. Both matrix axes use that same map.
    pub(crate) fn regional_bending_forms(
        &self,
        source: &ResolvedPlateChart,
        coordinate_scale: &[f64],
        max_projection_terms: usize,
        gate: &CancelGate,
    ) -> Result<Vec<Vec<f64>>, AcousticRealizeError> {
        let n = self.mass_normalized_modes.len();
        let mesh = &source.chart().mesh;
        let poll = || {
            if gate.is_requested() {
                Err(AcousticRealizeError::Cancelled)
            } else {
                Ok(())
            }
        };
        poll()?;
        if self.source_identity != source.specimen_identity()
            || n == 0
            || n > 64
            || coordinate_scale.len() != n
            || coordinate_scale.iter().any(|v| !v.is_finite() || *v <= 0.0)
            || self.dof_map.len() != 3 * mesh.nodes.len()
        {
            return Err(refuse(
                "regional bending requires the original complete material basis",
            ));
        }
        let work = mesh
            .tris
            .len()
            .checked_mul(81 * n * n + 9 * n)
            .and_then(|w| {
                source
                    .regions()
                    .len()
                    .checked_mul(4 * n * n * n)
                    .and_then(|factor| w.checked_add(factor))
            })
            .ok_or_else(|| refuse("regional bending projection work overflows"))?;
        if work > max_projection_terms {
            return Err(refuse(
                "regional bending projection exceeds its declared work budget",
            ));
        }
        let mut forms = Vec::with_capacity(source.regions().len());
        for region in source.regions() {
            let mut form = vec![0.0; n * n];
            for &element in &region.input.region.triangle_indices {
                poll()?;
                let tri = mesh.tris[element];
                let (stiffness, _) = fs_plate::dkt_stiffness(
                    &tri.map(|v| mesh.nodes[v].0),
                    &tri.map(|v| mesh.nodes[v].1),
                    &region.section.d,
                    element,
                )
                .map_err(map_plate)?;
                let dofs: [Option<usize>; 9] =
                    core::array::from_fn(|i| self.dof_map[3 * tri[i / 3] + i % 3]);
                let local: Vec<[f64; 9]> = self
                    .mass_normalized_modes
                    .iter()
                    .zip(coordinate_scale)
                    .map(|(mode, scale)| dofs.map(|dof| dof.map_or(0.0, |i| mode[i] * scale)))
                    .collect();
                for i in 0..n {
                    poll()?;
                    for j in i..n {
                        let mut value = 0.0;
                        for a in 0..9 {
                            for b in 0..9 {
                                value += local[i][a] * stiffness[9 * a + b] * local[j][b];
                            }
                        }
                        form[n * i + j] += value;
                        form[n * j + i] = form[n * i + j];
                    }
                }
            }
            if form.iter().any(|v| !v.is_finite()) {
                return Err(refuse("regional bending form is not representable"));
            }
            forms.push(form);
        }
        poll()?;
        Ok(forms)
    }
}
