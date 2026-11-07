//! Actual-temperature source support for nonisothermal binary transport.

use super::*;

#[test]
fn g1_equal_compositions_keep_each_actual_temperature_and_source_receipt() {
    let (material, ids) = card(vec![claim(surface())]);
    let law = bind(material.clone(), query(300.0, 0.25));
    let fractions = [[0.25, 0.75]; 2];
    let resolved = law
        .resolve_nonisothermal_face(
            [Temperature::new(300.0), Temperature::new(400.0)],
            fractions,
            || false,
        )
        .unwrap();
    // The original triangulation supplies D(350,.25)=4.5e-9. Averaging the
    // endpoint D values would give 4e-9 and lose the source's joint response.
    close(resolved.left.properties()[0].value_si(), 2.5e-9);
    close(resolved.diffusivity_m2_s, 4.5e-9);
    close(resolved.face.properties()[0].value_si(), 4.5e-9);
    close(resolved.right.properties()[0].value_si(), 5.5e-9);
    for (receipt, temperature) in [
        (&resolved.left, 300.0),
        (&resolved.face, 350.0),
        (&resolved.right, 400.0),
    ] {
        assert_eq!(receipt.state_coordinate("temperature"), Some(temperature));
        assert_eq!(receipt.state_coordinate("x_a"), Some(0.25));
        assert_eq!(receipt.state_coordinate("pressure"), Some(100_000.0));
        assert_eq!(receipt.card_identity(), material.content_hash());
        assert_eq!(receipt.properties()[0].answer().receipt.selected, ids[0]);
    }
    // Equal compositions alone must never alias a receipt from a different T.
    assert_ne!(resolved.left, resolved.face);
    assert_ne!(resolved.left, resolved.right);
    let reverse = law
        .resolve_nonisothermal_face(
            [Temperature::new(400.0), Temperature::new(300.0)],
            fractions,
            || false,
        )
        .unwrap();
    assert_eq!(reverse.left, resolved.right);
    assert_eq!(reverse.face, resolved.face);
    assert_eq!(reverse.right, resolved.left);
    assert_eq!(reverse.diffusivity_m2_s, resolved.diffusivity_m2_s);
}

#[test]
fn g1_joint_source_support_qualifies_actual_endpoints_without_fictitious_crossstates() {
    let (material, ids) = card(vec![claim(PropertyValue::Surface2D {
        axes: ["temperature".into(), "x_a".into()],
        axis_dims: [Temperature::DIMS, Dimensionless::DIMS],
        samples: vec![
            ([300.0, 0.0], 1e-9),
            ([400.0, 1.0], 9e-9),
            ([400.0, 0.8], 5e-9),
        ],
        triangles: vec![[0, 1, 2]],
        dims: DIMS,
    })]);
    let law = bind(material, query(300.0, 0.0));
    let before = law.clone();
    let fractions = [[0.0, 1.0], [1.0, 0.0]];
    let resolved = law
        .resolve_nonisothermal_face(
            [Temperature::new(300.0), Temperature::new(400.0)],
            fractions,
            || false,
        )
        .unwrap();
    close(resolved.left.properties()[0].value_si(), 1e-9);
    close(resolved.diffusivity_m2_s, 5e-9);
    close(resolved.right.properties()[0].value_si(), 9e-9);
    assert_eq!(resolved.face.state_coordinate("temperature"), Some(350.0));
    assert_eq!(resolved.face.state_coordinate("x_a"), Some(0.5));
    for receipt in [&resolved.left, &resolved.face, &resolved.right] {
        assert_eq!(receipt.properties()[0].answer().receipt.selected, ids[0]);
    }
    // At 350 K the source only supports .4<=x_A<=.5. Neither crossed endpoint
    // (350,0) nor (350,1) is part of this actual nonisothermal face transaction.
    assert!(matches!(
        law.resolve_face(Temperature::new(350.0), fractions, || false),
        Err(BinaryLiquidError::Source(_))
    ));
    let reverse = law
        .resolve_nonisothermal_face(
            [Temperature::new(400.0), Temperature::new(300.0)],
            [fractions[1], fractions[0]],
            || false,
        )
        .unwrap();
    assert_eq!(reverse.left, resolved.right);
    assert_eq!(reverse.face, resolved.face);
    assert_eq!(reverse.right, resolved.left);
    assert_eq!(law, before);
}
