use super::*;
use crate::conjugate::FlowFace;

struct Caloric;
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> { Ok(300.0 + h) }
    fn maximum_temperature_slope(&self) -> f64 { 1.0 }
}

fn input<'a>(ids: &'a [Option<CellId>], law: &'a Caloric) -> FlowEnthalpyInput<'a> {
    FlowEnthalpyInput {
        cell_ids: ids, material_id: MaterialId(4), caloric: law,
        inlet_specific_h_j_kg: [Some(10.0), Some(11.0), Some(12.0), Some(13.0), Some(14.0), Some(15.0)],
        max_cells: ids.len(), max_connections: 1000,
    }
}

#[test]
fn graded_three_axis_fluxes_keep_direction_area_density_and_unique_face_identity() {
    let domain = VoxelDomain::graded([
        vec![0.5, 1.5, 2.0], vec![0.25, 0.75], vec![0.5, 1.0],
    ]).unwrap();
    let ids: Vec<_> = (0..domain.cell_count()).map(|c| Some(CellId(c as u64))).collect();
    let velocity = [-0.25, 0.125, -0.5];
    let field = FlowField::from_face_velocity(&domain, [FlowFace::Fixed; 6],
        |axis, _| velocity[axis]).unwrap();
    let (flows, boundaries) = mass_graph(&domain, &field, 8.0, input(&ids, &Caloric), &mut || Ok(())).unwrap();
    let mut identities = BTreeSet::new();
    let mut balance = vec![0.0; domain.cell_count()];
    // Each interior face appears exactly once; positive axis is independent of
    // actual flow direction. All three unequal face areas enter in SI units.
    assert_eq!(flows.len(), 20);
    assert_eq!(boundaries.len(), 32);
    for flow in &flows {
        assert!(identities.insert(flow.id));
        let a = flow.a.0 as usize;
        let b = flow.b.0 as usize;
        let ca = domain.coords(a);
        let cb = domain.coords(b);
        let axis = (0..3).find(|&axis| ca[axis] != cb[axis]).unwrap();
        assert_eq!(flow.kg_s, 8.0 * domain.face_area(a, axis) * velocity[axis]);
        assert_eq!(cb[axis], ca[axis] + 1);
        balance[a] -= flow.kg_s;
        balance[b] += flow.kg_s;
    }
    for boundary in &boundaries {
        match *boundary {
            MassBoundary::Inlet { id, cell, kg_s, specific_h_j_kg, .. } => {
                assert!(identities.insert(id));
                let face = (specific_h_j_kg - 10.0) as usize;
                assert!([1, 2, 5].contains(&face));
                assert_eq!(kg_s, 8.0 * domain.face_area(cell.0 as usize, face / 2) * velocity[face / 2].abs());
                balance[cell.0 as usize] += kg_s;
            }
            MassBoundary::Outlet { id, cell, kg_s } => {
                assert!(identities.insert(id));
                balance[cell.0 as usize] -= kg_s;
            }
        }
    }
    assert!(balance.iter().all(|&rate| rate == 0.0));
}

#[test]
fn closed_circulation_has_no_external_enthalpy_and_no_duplicate_faces() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let field = FlowField::from_face_velocity(&domain, [FlowFace::Wall; 6], |axis, centre| {
        if axis == 0 { if centre[1] < 1.0 { 0.25 } else { -0.25 } }
        else if axis == 1 { if centre[0] < 1.0 { -0.25 } else { 0.25 } }
        else { 0.0 }
    }).unwrap();
    let (flows, boundaries) = mass_graph(&domain, &field, 4.0, input(&ids, &Caloric), &mut || Ok(())).unwrap();
    assert_eq!(flows.len(), 4);
    assert!(boundaries.is_empty());
    assert!(flows.iter().all(|f| f.kg_s.abs() == 1.0));
}

#[test]
fn divergence_missing_backflow_and_connection_budget_refuse_before_history_creation() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let ids = [Some(CellId(0)), Some(CellId(1))];
    let mut faces = [FlowFace::Wall; 6];
    faces[0] = FlowFace::Fixed;
    faces[1] = FlowFace::Free;
    let divergent = FlowField::from_face_velocity(&domain, faces,
        |axis, centre| if axis == 0 { 0.125 * (1.0 + centre[0]) } else { 0.0 }).unwrap();
    assert!(matches!(mass_graph(&domain, &divergent, 1.0, input(&ids, &Caloric), &mut || Ok(())),
        Err(DepositionError::InvalidInput(_))));
    let reverse = FlowField::from_face_velocity(&domain, faces,
        |axis, _| if axis == 0 { -0.125 } else { 0.0 }).unwrap();
    let mut missing = input(&ids, &Caloric);
    missing.inlet_specific_h_j_kg[1] = None;
    assert!(matches!(mass_graph(&domain, &reverse, 1.0, missing, &mut || Ok(())),
        Err(DepositionError::InvalidInput(_))));
    let mut bounded = input(&ids, &Caloric);
    bounded.max_connections = 2; // One interior plus two boundary connections.
    assert!(matches!(mass_graph(&domain, &reverse, 1.0, bounded, &mut || Ok(())),
        Err(DepositionError::Budget(_))));
    assert!(matches!(mass_graph(&domain, &reverse, 1.0, input(&ids, &Caloric),
        &mut || Err(DepositionError::Cancelled)), Err(DepositionError::Cancelled)));
}
