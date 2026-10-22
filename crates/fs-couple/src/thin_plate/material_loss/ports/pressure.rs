//! Signed distributed pressure loads on the existing material plate eigenbasis.
use super::*;

impl ResolvedPlateChart {
    /// Consistent P1 pressure loads on this original mesh, without another
    /// eigensolve. Retained-basis consumers use the same quadrature as the
    /// pressure-port reducer below, then project their existing signed shapes.
    pub(crate) fn nodal_pressure_loads(
        &self, pressure_at_nodes: &[Vec<f64>], max_projection_terms: usize,
    ) -> Result<Vec<Vec<f64>>, AcousticRealizeError> {
        pressure_loads(&self.chart().mesh, pressure_at_nodes, max_projection_terms)
            .map(|(loads, _)| loads)
    }

    /// Reduce the material body once and project nodal pressure basis functions
    /// on EXACTLY those modes. Returns radiators and row-major structural-by-
    /// pressure overlaps C[r,j] = integral(phi_r * psi_j dA)/sqrt(modal_mass).
    /// The same C distributes pressure forces and observes displaced volume.
    ///
    /// Both traces are piecewise linear on the supplied triangles: the consistent
    /// triangle mass matrix integrates their product, rather than multiplying
    /// separately averaged traces. This is a P1 trace of the DKT displacement,
    /// not a claim of exact integration of its higher-order interior field.
    /// Signed and zero-mean pressure modes are retained. No independently solved
    /// eigenbasis, unsigned-area substitution or modal frequency is introduced.
    ///
    /// Internally positive/negative nodal loads enter the existing force-port
    /// reducer as unit-total-force footprints, then recombine with their actual
    /// signed integrals. Never divide by the signed net force: it may be zero.
    /// `max_projection_terms` bounds quadrature/decomposition plus the original
    /// force-port owner's projection work. Material and modal limits still apply.
    ///
    /// # Errors
    /// Empty/mismatched/nonfinite pressure basis, invalid mesh, projection budget
    /// overflow, unrepresentable normalization or ordinary material/modal refusal.
    pub fn certified_radiators_with_pressure_ports(
        &self, options: &PlateChartRadiation, pressure_at_nodes: &[Vec<f64>],
        max_projection_terms: usize,
    ) -> Result<(Vec<CompactBody>, Vec<f64>), AcousticRealizeError> {
        let mesh = &self.chart().mesh;
        let (loads, work) = pressure_loads(mesh, pressure_at_nodes, max_projection_terms)?;
        let mut footprints = Vec::new();
        let mut terms = Vec::with_capacity(loads.len());
        for load in &loads {
            let mut signed = Vec::with_capacity(2);
            for sign in [1.0, -1.0] {
                let magnitude: f64 = load.iter().map(|v| (sign*v).max(0.0)).sum();
                if !magnitude.is_finite() { return Err(refuse("pressure-load integral overflows")); }
                if magnitude == 0.0 { continue; }
                let weights: Vec<_> = load.iter().map(|v| (sign*v).max(0.0)/magnitude).collect();
                if weights.iter().zip(load).any(|(w,v)| sign*v > 0.0 && *w == 0.0) {
                    return Err(refuse("pressure-load normalization underflows"));
                }
                signed.push((footprints.len(), sign*magnitude));
                footprints.push(weights);
            }
            terms.push(signed);
        }
        let (radiators, ports) = self.certified_radiators_with_force_ports(
            options, &footprints, max_projection_terms-work)?;
        let count = radiators.len().checked_mul(loads.len())
            .ok_or_else(|| refuse("pressure overlap dimensions overflow"))?;
        let mut overlaps = vec![0.0; count];
        for (j, signed) in terms.iter().enumerate() {
            for r in 0..radiators.len() {
                let c = signed.iter().map(|(port,scale)| ports[*port][r]*scale).sum::<f64>();
                if !c.is_finite() { return Err(refuse("pressure overlap is unrepresentable")); }
                overlaps[r*loads.len()+j] = c;
            }
        }
        Ok((radiators, overlaps))
    }
}

// Consistent P1 face mass action. Count work before allocation or triangle scans.
fn pressure_loads(mesh: &PlateMesh, pressure: &[Vec<f64>], budget: usize)
    -> Result<(Vec<Vec<f64>>, usize), AcousticRealizeError>
{
    if pressure.is_empty() || mesh.nodes.is_empty() || mesh.tris.is_empty() {
        return Err(refuse("pressure projection needs a nonempty mesh and basis"));
    }
    let work = mesh.nodes.len().checked_mul(8)
        .and_then(|n| mesh.tris.len().checked_mul(12).and_then(|t| n.checked_add(t)))
        .and_then(|n| n.checked_mul(pressure.len()))
        .ok_or_else(|| refuse("pressure projection work overflows"))?;
    if work > budget { return Err(refuse("pressure projection exceeds its work budget")); }
    if pressure.iter().any(|p| p.len() != mesh.nodes.len() || p.iter().any(|v| !v.is_finite())) {
        return Err(refuse("pressure basis needs one finite value per original plate node"));
    }
    let mut loads = vec![vec![0.0;mesh.nodes.len()];pressure.len()];
    for tri in &mesh.tris {
        if tri.iter().any(|&i| i >= mesh.nodes.len()) {
            return Err(refuse("pressure projection triangle has an invalid node"));
        }
        let [i,j,k] = *tri;
        let (a,b,c) = (mesh.nodes[i],mesh.nodes[j],mesh.nodes[k]);
        let area = 0.5*((b.0-a.0)*(c.1-a.1)-(b.1-a.1)*(c.0-a.0));
        if !area.is_finite() || area <= 0.0 { return Err(refuse("pressure projection needs positive oriented triangle area")); }
        for (p,load) in pressure.iter().zip(&mut loads) {
            // Diagonal A/6 and off-diagonal A/12; scale before summing.
            for (node,other,third) in [(i,j,k),(j,k,i),(k,i,j)] {
                load[node] += (area/6.0)*p[node] + (area/12.0)*p[other] + (area/12.0)*p[third];
            }
        }
    }
    if loads.iter().flatten().any(|v| !v.is_finite()) { return Err(refuse("pressure quadrature overflows")); }
    Ok((loads,work))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn constant_and_affine_pressure_preserve_area_and_zero_net_load() {
        let mesh = PlateMesh::rectangle(0.4,0.3,4,4);
        let p: Vec<_> = mesh.nodes.iter().map(|&(x,_)| x-0.2).collect();
        let (load,_) = pressure_loads(&mesh,&[vec![1.0;25],p],10000).unwrap();
        assert!((load[0].iter().sum::<f64>()-0.12).abs()<1e-15);
        assert!(load[1].iter().sum::<f64>().abs()<1e-15);
        let work: f64 = mesh.nodes.iter().zip(&load[1]).map(|(&(x,_),f)| (x-0.2)*f).sum();
        // Integral_0^.3 Integral_0^.4 (x-.2)^2 dx dy.
        assert!((work-0.3*0.4_f64.powi(3)/12.0).abs()<1e-16);
        assert!(load[1].iter().any(|v| *v>0.0) && load[1].iter().any(|v| *v<0.0));
    }
    #[test]
    fn pressure_displacement_pairing_is_reciprocal_and_linear() {
        let mesh = PlateMesh::rectangle(0.4,0.3,3,3);
        let p: Vec<_> = mesh.nodes.iter().map(|&(x,y)| x+2.0*y).collect();
        let q: Vec<_> = mesh.nodes.iter().map(|&(x,y)| 0.3-x+y).collect();
        let combination: Vec<_> = p.iter().zip(&q).map(|(p,q)| 2.0*p-3.0*q).collect();
        let (loads,_) = pressure_loads(&mesh,&[p.clone(),q.clone(),combination],10000).unwrap();
        let a: f64 = p.iter().zip(&loads[1]).map(|(p,f)|p*f).sum();
        let b: f64 = q.iter().zip(&loads[0]).map(|(q,f)|q*f).sum();
        assert!((a-b).abs()<1e-16);
        for i in 0..p.len() { assert!((loads[2][i]-(2.0*loads[0][i]-3.0*loads[1][i])).abs()<1e-16); }
    }
    #[test]
    fn invalid_pressure_and_projection_budget_refuse() {
        let mesh=PlateMesh::rectangle(1.0,1.0,2,2);
        for p in [vec![],vec![1.0;8],vec![f64::NAN;9]] {
            assert!(pressure_loads(&mesh,&[p],10000).is_err());
        }
        assert!(pressure_loads(&mesh,&[],10000).is_err());
        let work=9*8+8*12;
        assert!(pressure_loads(&mesh,&[vec![1.0;9]],work-1).is_err());
        pressure_loads(&mesh,&[vec![1.0;9]],work).unwrap();
    }
}
