use super::*;

fn section(q: FreeSurfaceConvectionQuery<'_>, start: f64, end: f64) -> PlanarFreeSurfaceSection {
    PlanarFreeSurfaceSection {
        face: q.surface.face,
        outward_normal: q.outward_normal,
        start: Length::new(start),
        end: Length::new(end),
        downstream_direction: [0.0, 1.0, 0.0],
        relative_velocity_m_s: [0.0, 12.0, 0.0],
        regime: PlanarForcedConvectionRegime::Laminar,
    }
}

#[test]
fn unequal_plate_strips_resolve_local_heat_without_restarting_the_boundary_layer() {
    let model = source::model_transport(0.75, true);
    let bulk = model.evaluate(&source::point(320.0, 100000.0)).unwrap();
    let film = model.evaluate(&source::point(400.0, 100000.0)).unwrap();
    let base = query(480.0, &bulk, &film);
    let bounds = [0.0, 0.01, 0.04, 0.07, 0.1];
    let width = 0.05;
    let bindings: Vec<_> = bounds
        .windows(2)
        .enumerate()
        .map(|(i, pair)| {
            let mut strip = section(base, pair[0], pair[1]);
            strip.face.liquid[1] = i + 2;
            strip.face.gas[1] = i + 2;
            strip.face.area_m2 = width * (pair[1] - pair[0]);
            strip
        })
        .collect();
    let provider = PlanarFreeSurfaceSections {
        sections: &bindings,
    };
    let mut heat = 0.0;
    let mut previous = f64::INFINITY;
    for (i, binding) in bindings.iter().enumerate() {
        let mut q = base;
        q.patch_index = i;
        q.surface.face = binding.face;
        let coefficient = provider.coefficient(q).unwrap().unwrap();
        let h = coefficient.value.value();
        assert!(h < previous, "heat transfer must decline downstream");
        previous = h;
        let gas = film.state();
        let pr = gas.specific_heat_cp * gas.dynamic_viscosity / gas.thermal_conductivity;
        let expected = 0.664
            * gas.thermal_conductivity
            * pr.cbrt()
            * (12.0 * gas.density / gas.dynamic_viscosity).sqrt()
            / (binding.end.value().sqrt() + binding.start.value().sqrt());
        assert!((h / expected - 1.0).abs() < 2e-12);
        heat += h * binding.face.area_m2 * (480.0 - 320.0);
        // The gas-side film thickness stays geometry, never a second resistance.
        let direct = evaluate_flat_plate_section_from_material(
            CorrelationId::FlatPlateLaminarAverage,
            &film,
            12.0,
            binding.start,
            binding.end,
        )
        .unwrap();
        assert_eq!(&coefficient, direct.coefficient());
    }
    let mean = forced().coefficient(base).unwrap().unwrap().value.value();
    let expected = mean * width * 0.1 * 160.0;
    assert!((heat / expected - 1.0).abs() < 2e-12);
    assert!(!provider.admits_uniform_boundary_evolution());
}

#[test]
fn sections_reject_reversed_flow_stale_geometry_and_wrong_film_source() {
    let model = source::model(0.75);
    let bulk = model.evaluate(&source::point(300.0, 100000.0)).unwrap();
    let film = model.evaluate(&source::point(350.0, 100000.0)).unwrap();
    let q = query(400.0, &bulk, &film);
    let original = section(q, 0.02, 0.1);
    let coefficients = |binding| {
        let sections = [binding];
        PlanarFreeSurfaceSections {
            sections: &sections,
        }
        .coefficient(q)
    };
    for binding in [
        PlanarFreeSurfaceSection {
            relative_velocity_m_s: [0.0, -12.0, 0.0],
            ..original
        },
        PlanarFreeSurfaceSection {
            relative_velocity_m_s: [0.01, 12.0, 0.0],
            ..original
        },
        PlanarFreeSurfaceSection {
            downstream_direction: [0.0, 0.0, 1.0],
            ..original
        },
        PlanarFreeSurfaceSection {
            outward_normal: [-1.0, 0.0, 0.0],
            ..original
        },
        PlanarFreeSurfaceSection {
            start: Length::new(0.1),
            ..original
        },
    ] {
        assert!(coefficients(binding).is_err());
    }
    let bindings = [original];
    let provider = PlanarFreeSurfaceSections {
        sections: &bindings,
    };
    let mut changed = q;
    changed.surface.face.area_m2 *= 2.0;
    assert!(provider.coefficient(changed).is_err());
    changed = q;
    changed.patch_index = 1;
    assert!(provider.coefficient(changed).is_err());
    changed = q;
    changed.film = &bulk;
    assert!(provider.coefficient(changed).is_err());
    // At exact equilibrium all geometric and regime refusals remain active.
    changed = query(350.0, &film, &film);
    changed.surface.face.liquid[0] += 1;
    assert!(provider.coefficient(changed).is_err());
}

#[test]
fn section_heat_is_invariant_under_a_world_frame_rotation() {
    let model = source::model(0.75);
    let bulk = model.evaluate(&source::point(300.0, 100000.0)).unwrap();
    let film = model.evaluate(&source::point(350.0, 100000.0)).unwrap();
    let q = query(400.0, &bulk, &film);
    let original = section(q, 0.03, 0.1);
    let bindings = [original];
    let expected = PlanarFreeSurfaceSections {
        sections: &bindings,
    }
    .coefficient(q)
    .unwrap()
    .unwrap()
    .value
    .value();
    let rotate = |v: [f64; 3]| {
        [
            (-10.0 * v[0] + 2.0 * v[1] + 11.0 * v[2]) / 15.0,
            (2.0 * v[0] - v[1] + 2.0 * v[2]) / 3.0,
            (5.0 * v[0] + 14.0 * v[1] + 2.0 * v[2]) / 15.0,
        ]
    };
    let bindings = [PlanarFreeSurfaceSection {
        outward_normal: rotate(original.outward_normal),
        downstream_direction: rotate(original.downstream_direction),
        relative_velocity_m_s: rotate(original.relative_velocity_m_s),
        ..original
    }];
    let mut q = q;
    q.outward_normal = rotate(q.outward_normal);
    let actual = PlanarFreeSurfaceSections {
        sections: &bindings,
    }
    .coefficient(q)
    .unwrap()
    .unwrap()
    .value
    .value();
    assert!((actual / expected - 1.0).abs() < 2e-13);
}
