//! External heat flux on the original complete wall-temperature trace.

use super::*;

impl<'c> ThermoviscousFe<'c> {
    /// Prepare an inward surface heat-flux shape on original no-slip faces.
    ///
    /// `density_per_m2` contains one real hierarchical P2 coefficient [m^-2]
    /// per GLOBAL scalar DOF: vertices first, then edge bubbles. Edge entries
    /// are coefficients, NOT midpoint samples. Only coordinates on `faces`
    /// may be nonzero. A signed drive amplitude [W] produces inward heat flux
    /// [W/m^2]; neither the shape nor its possibly zero integral is normalized.
    ///
    /// The SAME original boundary mass used by passive conductive walls
    /// projects this flux onto every active temperature test. Curved faces use
    /// the original map and declared surface quadrature. Only exact isothermal
    /// coordinates are eliminated, and fully constrained traces refuse.
    /// Selected faces must be unique original no-slip faces, disjoint from
    /// exact isothermal faces. Full-vector moving walls remain supported on
    /// the same fixed reference surface.
    ///
    /// On a finite conductive wall this is an additional heat input at the
    /// gas/solid interface, shared with the original solid temperature. It is
    /// not a prescribed remote backing temperature or a second wall loss.
    /// Harmonic, loaded, conforming-interface and causal heat-source APIs reuse
    /// the same cached source and signed theta*q work, theta=T'/T0. First-order
    /// heat flux and acoustic availability input are different quantities.
    /// No mean heating, heat-production mechanism or moving geometry is inferred.
    pub fn surface_heat_source<'a>(
        &'a self,
        faces: &[[u32; 3]],
        density_per_m2: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermoviscousHeatSource<'a, 'c>, ThermoviscousError> {
        poll(cancelled)?;
        if faces.is_empty()
            || faces.len() > self.wall_faces.len()
            || density_per_m2.len() != self.space.ndof
            || density_per_m2.iter().any(|x| !x.is_finite())
        {
            return Err(ThermoviscousError::InvalidInput(
                "original nonempty surface and complete finite P2 heat-flux shape required",
            ));
        }
        let mut seen = BTreeSet::new();
        let mut selected = vec![false; self.space.ndof];
        let mut traces = Vec::with_capacity(faces.len());
        let mut active = false;
        // Admit the complete selection before assembling any physical load.
        for &face in faces {
            poll(cancelled)?;
            let face = canonical(face);
            if !self.wall_faces.contains(&face)
                || self.isothermal_faces.contains(&face)
                || !seen.insert(face)
            {
                return Err(ThermoviscousError::InvalidInput(
                    "unique non-isothermal original wall faces required for surface heat",
                ));
            }
            let mut dofs = [0; 6];
            dofs[..3].copy_from_slice(&face.map(|v| v as usize));
            for (offset, (a, b)) in [(0, 1), (0, 2), (1, 2)].into_iter().enumerate() {
                let edge = self
                    .space
                    .complex
                    .edges
                    .binary_search(&[face[a], face[b]])
                    .expect("original admitted no-slip face edge");
                dofs[3 + offset] = self.space.edge_off + edge;
            }
            for &dof in &dofs {
                selected[dof] = true;
                active |= self.temperature_map[dof].is_some();
            }
            traces.push((face, dofs));
        }
        if !active
            || density_per_m2
                .iter()
                .zip(&selected)
                .any(|(&coefficient, &on_surface)| !on_surface && coefficient != 0.0)
        {
            return Err(ThermoviscousError::InvalidInput(
                "active original temperature trace and heat-flux support required",
            ));
        }
        let mut load = vec![Dd::ZERO; self.n - self.temperature_start];
        let mut integral = Dd::ZERO;
        for (face, dofs) in traces {
            poll(cancelled)?;
            // A single original face assembly bounds cancellation latency and
            // temporary sparse storage; no second surface quadrature is copied.
            let mass = self.fe.boundary_mass(&[face])?;
            for global in dofs {
                poll(cancelled)?;
                let (columns, values) = mass.row(global);
                let mut value = Dd::ZERO;
                for (&column, &coefficient) in columns.iter().zip(values) {
                    value =
                        value + Dd::from_f64(coefficient) * Dd::from_f64(density_per_m2[column]);
                }
                // The hierarchical constant is one at vertices and zero at
                // edge bubbles. Sum its FULL trace, before thermal constraints.
                if global < self.space.complex.vertex_count {
                    integral = integral + value;
                }
                if let Some(thermal) = self.temperature_map[global] {
                    load[thermal] = load[thermal] + value;
                }
            }
        }
        let integral = integral.to_f64();
        let load: Vec<_> = load.into_iter().map(Dd::to_f64).collect();
        if !integral.is_finite() || load.iter().any(|x| !x.is_finite()) {
            return Err(ThermoviscousError::NonFinite);
        }
        poll(cancelled)?;
        Ok(ThermoviscousHeatSource {
            model: self,
            load,
            integral,
        })
    }
}
