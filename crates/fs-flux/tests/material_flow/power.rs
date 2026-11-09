//! Original material identity and SI scaling of generalized wall work.

use super::*;

struct Source {
    state: FluidTransportState,
    reads: Rc<Cell<usize>>,
    tag: usize,
}
impl FluidTransport for Source {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.reads.set(self.reads.get() + 1);
        Ok(self.state)
    }
}

#[test]
fn original_homogeneous_and_cell_sources_keep_dynamic_viscosity_power_under_density_swap() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let shear = 0.4;
    let mu = 0.8;
    let velocity = |p: [f64; 2]| [shear * p[1], 0.0];
    let previous: Vec<_> = mesh
        .edges
        .iter()
        .flat_map(|e| {
            let a = velocity(mesh.verts[e.verts.0]);
            let b = velocity(mesh.verts[e.verts.1]);
            [
                0.5 * (a[0] + b[0]) * e.normal[0],
                (b[0] - a[0]) * e.normal[0] / 12.0,
            ]
        })
        .collect();
    let edges: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .collect();
    let left: Vec<_> = edges
        .iter()
        .copied()
        .filter(|&e| mesh.edges[e].normal[0] < -0.5)
        .collect();
    let boundaries: Vec<_> = edges
        .iter()
        .map(|&edge| MaterialStokesBoundary::Velocity {
            edge,
            value: &velocity,
        })
        .collect();
    for density in [2.0, 7.0] {
        let reads = Rc::new(Cell::new(0));
        let source = Source {
            state: transport(density, mu),
            reads: Rc::clone(&reads),
            tag: 99,
        };
        let owner = MaterialFluxSystem::new(
            &system,
            &source,
            MaterialFluxNumerics {
                linear_tolerance: 1e-12,
                ..MaterialFluxNumerics::default()
            },
        )
        .unwrap();
        assert!(std::ptr::eq(owner.material(), &source));
        let result = owner
            .bdf1_stokes_mixed_step(&zero_force, &boundaries, &previous, 0.2)
            .unwrap();
        let before = reads.get();
        let power = result.boundary_power(&edges, &[Some(0.0)]).unwrap();
        near(power.prescribed_velocity_w_per_m, mu * shear * shear, 2e-10);
        assert_eq!(power.imposed_traction_w_per_m, 0.0);
        assert_eq!(power.backflow_resistance_w_per_m, 0.0);
        let shifted = result.boundary_power(&left, &[Some(7.0)]).unwrap();
        near(shifted.prescribed_velocity_w_per_m, 7.0 * shear / 2.0, 2e-9);
        assert_eq!(shifted.pressure_offsets_pa, vec![(0, Some(7.0))]);
        assert!(result.boundary_power(&edges, &[None]).is_err());
        assert_eq!(reads.get(), before);

        // Different ORIGINAL cell densities still share one reference metric;
        // steady shear work depends on mu, not local or reference rho.
        let field = owner
            .bind_stokes_materials(
                (0..2)
                    .map(|cell| Source {
                        state: transport(3.0 + 8.0 * cell as f64, mu),
                        reads: Rc::clone(&reads),
                        tag: cell,
                    })
                    .collect(),
            )
            .unwrap();
        let result = field
            .bdf1_stokes_mixed_step(&zero_force, &boundaries, &previous, 0.2)
            .unwrap();
        let before = reads.get();
        near(
            result
                .boundary_power(&edges, &[Some(0.0)])
                .unwrap()
                .prescribed_velocity_w_per_m,
            mu * shear * shear,
            2e-9,
        );
        near(
            result
                .boundary_power(&left, &[Some(7.0)])
                .unwrap()
                .prescribed_velocity_w_per_m,
            7.0 * shear / 2.0,
            2e-9,
        );
        for (cell, sample) in result.samples().iter().enumerate() {
            assert_eq!(sample.material().tag, cell);
            assert!(Rc::ptr_eq(&sample.material().reads, &reads));
        }
        assert_eq!(
            reads.get(),
            before,
            "reporting must retain original source snapshots"
        );
        assert_eq!(previous.len(), system.n_u);
    }
}
