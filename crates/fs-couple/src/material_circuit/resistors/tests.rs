use super::*;

#[test]
fn only_explicit_native_field_inductors_bypass_numeric_storage_admission() {
    for placeholder in [f64::NAN, 0.0, -1.0] {
        let graph = CircuitGraph {
            node_count: 3,
            branches: vec![
                (
                    1,
                    0,
                    Branch::Inductor {
                        henries: placeholder,
                    },
                ),
                (2, 0, Branch::Capacitor { farads: 1.0 }),
                (2, 0, Branch::Inductor { henries: 2.0 }),
            ],
            transformers: vec![],
        };
        assert!(lower_resistors(&graph, vec![]).is_err());
        let (lowered, sources) =
            lower_resistors_with_field_inductors(&graph, vec![], &[0]).unwrap();
        assert!(sources.is_empty());
        let Branch::Inductor { henries } = lowered.branches[0].2 else {
            panic!("native field branch was replaced")
        };
        assert_eq!(henries.to_bits(), placeholder.to_bits());
        assert_eq!(lowered.branches[1], graph.branches[1]);
        assert_eq!(lowered.branches[2], graph.branches[2]);
        for fields in [&[0, 0][..], &[1][..], &[3][..]] {
            assert!(lower_resistors_with_field_inductors(&graph, vec![], fields).is_err());
        }
        let mut bad_capacitor = graph.clone();
        bad_capacitor.branches[1].2 = Branch::Capacitor { farads: 0.0 };
        assert!(lower_resistors_with_field_inductors(&bad_capacitor, vec![], &[0]).is_err());
        let mut bad_inductor = graph;
        bad_inductor.branches[2].2 = Branch::Inductor { henries: f64::NAN };
        assert!(lower_resistors_with_field_inductors(&bad_inductor, vec![], &[0]).is_err());
    }
}
