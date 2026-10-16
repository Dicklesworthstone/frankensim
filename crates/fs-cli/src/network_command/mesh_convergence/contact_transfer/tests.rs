use super::*;

const BASE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cooling-network/nonmatching-contact-hotspot.json"));
const LIMITS: TetRefinementLimits = TetRefinementLimits {
    max_vertices: 20_000, max_tetrahedra: 100_000,
};

fn contact(root: &J) -> &J {
    &root.path(&["solid", "contacts"]).unwrap().as_array().unwrap()[0]
}
fn side(root: &J, key: &str) -> BTreeSet<[u32; 3]> {
    contact(root).get("nonmatching").unwrap().get(key).unwrap().as_array().unwrap()
        .iter().map(|face| {
            let mut face = indices::<3>(face, key, 20_000).unwrap();
            face.sort_unstable(); face
        }).collect()
}
fn first_solid(root: &J) -> Vec<usize> {
    root.path(&["solid", "element_materials"]).unwrap().as_array().unwrap()
        .iter().enumerate().filter_map(|(i, name)|
            (name.as_str() == Some("spreader")).then_some(i)).collect()
}

#[test]
fn one_sided_refinement_keeps_the_other_trace_and_nonuniform_contact_modes() {
    let gate = CancelGate::new();
    let mut root = J::parse(BASE).unwrap();
    let original = root.clone();
    let untouched = side(&root, "side_b_faces");
    for _ in 0..3 {
        let request = Request::parse(&encode(&root).unwrap()).unwrap();
        let marked = first_solid(&root);
        root = with_context(&request, &gate, |cx|
            refine_request_selected(cx, &root, &request, LIMITS, Some(&marked))).unwrap();
        assert_eq!(side(&root, "side_b_faces"), untouched);
        let refined = Request::parse(&encode(&root).unwrap()).unwrap();
        let mut a = BTreeSet::new();
        let mut b = BTreeSet::new();
        let marked: BTreeSet<_> = first_solid(&root).into_iter().collect();
        for (i, tet) in refined.mesh.complex().tets.iter().enumerate() {
            if marked.contains(&i) { a.extend(tet.iter().copied()); }
            else { b.extend(tet.iter().copied()); }
        }
        assert!(a.is_disjoint(&b), "refinement welded the contact traces");
        let mut temperature = vec![300.0; refined.mesh.vertex_count()];
        let mut lambda = vec![0.0; refined.mesh.vertex_count()];
        for vertex in a {
            let i = vertex as usize;
            let mode = refined.mesh.positions()[i][1] / 0.1 - 0.5;
            temperature[i] += mode;
            lambda[i] = mode;
        }
        let interfaces = &refined.contacts.as_ref().unwrap().interfaces;
        let contraction = with_context(&refined, &gate, |cx|
            interfaces.nonmatching_log_resistance_pullback(cx, "bondline", &temperature, &lambda)
                .map_err(producer)).unwrap().unwrap();
        // Integral on a 0.1 x 0.1 plane of (y/0.1-1/2)^2 / 0.01.
        // Both mean jumps are zero, so a mean-only replacement would fail.
        assert!((contraction - 1.0 / 12.0).abs() < 1e-11);
        let flux = &interfaces.fluxes(&temperature).unwrap()[0];
        assert!((flux.area_m2 - 0.01).abs() < 1e-12);
        assert!(flux.heat_rate_a_to_b_w.abs() < 1e-10);
        let source = refined.solid_data.nodal_source.as_ref().unwrap();
        let power: f64 = refined.mesh.complex().tets.iter().enumerate().map(|(i, tet)|
            refined.mesh.element_volume(i) * tet.iter().map(|&v|source.at(v as usize)).sum::<f64>() / 4.0).sum();
        assert!((power - 1.0).abs() < 1e-11);
        for key in ["resistance_m2_k_w", "source", "side_a_material", "side_b_material"] {
            assert_eq!(contact(&root).get(key), contact(&original).get(key));
        }
    }
    assert!(side(&root, "side_a_faces").len() > side(&original, "side_a_faces").len());
}

#[test]
fn global_refinement_preserves_the_original_contact_budget_on_failure() {
    let gate = CancelGate::new();
    let root = J::parse(BASE).unwrap();
    let request = Request::parse(BASE).unwrap();
    let first = with_context(&request, &gate, |cx|
        refine_request(cx, &root, &request, LIMITS)).unwrap();
    assert_eq!(side(&first, "side_a_faces").len(), 8);
    assert_eq!(side(&first, "side_b_faces").len(), 48);
    let accepted = first.clone();
    let request = Request::parse(&encode(&first).unwrap()).unwrap();
    let failure = with_context(&request, &gate, |cx|
        refine_request(cx, &first, &request, LIMITS)).unwrap_err();
    assert_eq!(failure.code, "cooling-network-mesh-budget");
    assert!(failure.message.contains("24976"));
    assert!(failure.message.contains("max_pair_tests=10000"));
    assert_eq!(first, accepted, "a rejected transfer changed its accepted input");
}
