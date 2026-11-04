//! Exact affine electrical response on the retained material and mass owners.
use super::*;

impl<'a, C: SpatialConductorMaterial> PreparedSpatialConductionKick<'a, C> {
    /// Prepare one held source configuration and its six-dimensional coupling.
    ///
    /// `mass` must cover the SAME original vertices/tetrahedra and material
    /// card/state points as `problem`; its canonical resolved inertia enters
    /// the native body unchanged. Mesh positions, conductivity axes, B and
    /// wrench/velocity coordinates all use that original body frame. The
    /// explicit world pose only locates this fixed snapshot for native momentum.
    ///
    /// The field's symmetric linear Ohmic law, fixed imposed B and native
    /// linear voltage/current traces make Lorentz load affine in rigid velocity.
    /// One source solve and six homogeneous unit-velocity solves recover it;
    /// no differencing of nearly equal driven fields is used. `M-dt*A/2` is
    /// factored without dropping couplings, symmetrizing or clipping drag.
    /// Native field gates admit every probe and the final actual trajectory.
    ///
    /// Contact surface work, inductive storage, changing material inventory,
    /// arbitrary joint constraints and free-configuration motion are excluded.
    /// This stage reports heat for a separate receiver; it does not evolve T.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        problem: &'a SpatialConduction<C>,
        mass: &'a TetrahedralMass<'a>,
        boundary: &SpatialConductionBoundary,
        flux_density_t: &[[f64; 3]],
        world_from_body: Se3,
        duration_s: f64,
        controls: SpatialConductionKickControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || ![
                controls.angular_velocity_tolerance_rad_s,
                controls.linear_velocity_tolerance_m_s,
                controls.energy_tolerance_j,
            ]
            .iter()
            .all(|value| value.is_finite() && *value > 0.0)
            || controls.maximum_field_solves < 7
        {
            return Err(invalid(
                "positive duration/tolerances and seven admitted preparation field solves required",
            ));
        }
        let cells = problem.mesh().tets.len();
        let vertices = problem.positions_m().len();
        if cells == 0
            || cells > controls.loads.max_cells
            || vertices > controls.loads.max_vertices
            || flux_density_t.len() != cells
            || !problem.contacts().is_empty()
            || mass.nodes_m().len() != vertices
            || mass.tetrahedra().len() != cells
        {
            return Err(invalid(
                "complete bounded contact-free field and original tetrahedral mass required",
            ));
        }
        for (original, current) in mass.nodes_m().iter().zip(problem.positions_m()) {
            poll(cancelled)?;
            if original != current {
                return Err(invalid(
                    "mass and field use different original vertex coordinates",
                ));
            }
        }
        for (original, current) in mass.tetrahedra().iter().zip(&problem.mesh().tets) {
            poll(cancelled)?;
            if *original != current.map(|v| v as usize) {
                return Err(invalid("mass and field use different original tetrahedra"));
            }
        }
        for (region, &cell) in mass
            .reference()
            .regions()
            .iter()
            .zip(mass.region_elements())
        {
            poll(cancelled)?;
            let density = region.material();
            let electrical = problem.conductors()[cell].material();
            if density.material() != electrical.material()
                || density.card_identity() != electrical.card_identity()
                || density.query_point().len() != electrical.query_point().len()
                || density
                    .query_point()
                    .iter()
                    .zip(electrical.query_point())
                    .any(|((a, x), (b, y))| a != b || x.to_bits() != y.to_bits())
                || density.axis_quantities() != electrical.axis_quantities()
            {
                return Err(invalid(
                    "electrical and inertial cells must retain the same material card and physical state",
                ));
            }
        }
        let model = ArticulatedModel::new(vec![Link::new(
            "original tetrahedral conductor",
            None,
            Se3::identity(),
            JointModel::FIXED,
            mass.reference().inertia(),
        )])?;
        let selected_cells: Vec<_> = (0..cells).collect();
        let (_, source_wrench) = solve_field(
            problem,
            boundary,
            flux_density_t,
            &selected_cells,
            Twist::zero(),
            controls,
            cancelled,
        )?;
        // Preserve the original trace partition and gauges while suppressing
        // only authored source amplitudes. Every basis response is a native
        // charge-continuity solve on the SAME complete spatial operator.
        let mut homogeneous = boundary.clone();
        for trace in &mut homogeneous.voltage {
            trace.voltage_v = [0.0; 3];
        }
        for trace in &mut homogeneous.current {
            trace.inward_current_density_a_m2 = 0.0;
        }
        let mut response = [0.0; 36];
        for col in 0..6 {
            poll(cancelled)?;
            let mut unit = [0.0; 6];
            unit[col] = 1.0;
            let (_, force) = solve_field(
                problem,
                &homogeneous,
                flux_density_t,
                &selected_cells,
                twist(unit),
                controls,
                cancelled,
            )?;
            for (row, value) in force.to_array().into_iter().enumerate() {
                response[row * 6 + col] = finite(value)?;
            }
        }
        let mut matrix = mass.reference().inertia().matrix()?.m;
        for i in 0..36 {
            matrix[i] = finite((-0.5 * duration_s).mul_add(response[i], matrix[i]))?;
        }
        let mut schur = LuWorkspace::new(6)?;
        schur.factor(&matrix)?;
        let prepared = Self {
            problem,
            mass,
            model,
            world_from_body,
            boundary: boundary.clone(),
            flux_density_t: flux_density_t.to_vec(),
            selected_cells,
            duration_s,
            controls,
            source_wrench: source_wrench.to_array(),
            response,
            schur,
        };
        poll(cancelled)?;
        Ok(prepared)
    }
}
