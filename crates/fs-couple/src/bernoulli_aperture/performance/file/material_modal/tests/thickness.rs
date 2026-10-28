//! Authored taper is geometry, with the original material and loss owners.
//! These meshes and material claims are manufactured integration fixtures.
use super::*;

const TENSOR: &str = "loss tensor eta11 eta22 eta66 eta12 eta16 eta26 1 1000000";

fn mapped(text: &str, mut thickness: impl FnMut(usize, usize, f64) -> f64) -> String {
    let mut output = String::new();
    let mut blade = 0;
    let mut seen_blade = false;
    let mut element = 0;
    let mut defaults = Vec::<f64>::new();
    for line in text.lines() {
        let fields = line.split_ascii_whitespace().collect::<Vec<_>>();
        match fields.first().copied() {
            Some("blade") => {
                if seen_blade {
                    blade += 1;
                }
                seen_blade = true;
                element = 0;
                defaults.clear();
            }
            Some("section") => defaults.push(fields[5].parse().unwrap()),
            Some("triangle") => {
                let section: usize = fields[4].parse().unwrap();
                output.push_str(&format!(
                    "{line} thickness {}\n",
                    thickness(blade, element, defaults[section])
                ));
                element += 1;
                continue;
            }
            _ => {}
        }
        output.push_str(line);
        output.push('\n');
    }
    output
}

#[test]
fn explicit_uniform_thickness_preserves_regions_modes_loss_and_every_pressure_sample() {
    let text = input();
    let original = build(&text, 37);
    let explicit = build(&mapped(&text, |_, _, h| h), 1);
    for blade in 0..2 {
        let before = valve(&original).material_plates()[blade].source();
        let after = valve(&explicit).material_plates()[blade].source();
        assert_eq!(before.specimen_identity(), after.specimen_identity());
        assert_eq!(before.mass_kg(), after.mass_kg());
        assert_eq!(before.regions().len(), after.regions().len());
    }
    assert_eq!(explicit.material_records().len(), 4);
    for record in explicit.material_records() {
        assert_eq!(record.section(), record.region());
    }
    assert_eq!(
        valve(&original).viscosity_matrix_n_s_m(),
        valve(&explicit).viscosity_matrix_n_s_m()
    );
    assert_eq!(render(original).0, render(explicit).0);
}

#[test]
fn local_taper_changes_integrated_mass_bending_loss_and_microphone_without_new_material() {
    let text = input();
    let tapered = mapped(&text, |blade, element, h| {
        if blade == 0 && element == 0 {
            h * 0.5
        } else {
            h
        }
    });
    let original = build(&text, 37);
    let changed = build(&tapered, 37);
    let before = valve(&original).material_plates()[0].source();
    let after = valve(&changed).material_plates()[0].source();
    let mesh = &before.chart().mesh;
    let triangle = mesh.tris[0];
    let (a, b, c) = (
        mesh.nodes[triangle[0]],
        mesh.nodes[triangle[1]],
        mesh.nodes[triangle[2]],
    );
    let area = 0.5 * ((b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0));
    let old = &before.regions()[0];
    let thin = after
        .regions()
        .iter()
        .find(|region| region.input.region.triangle_indices.contains(&0))
        .unwrap();
    let removed_mass = old.section.density * old.section.thickness * 0.5 * area;
    assert!((before.mass_kg() - after.mass_kg() - removed_mass).abs() < 1e-13 * before.mass_kg());
    assert_eq!(thin.section.thickness, old.section.thickness * 0.5);
    for (new, old) in thin.section.d.iter().zip(old.section.d) {
        assert!((new - old * 0.125).abs() <= 1e-13 * old.abs());
    }
    assert_eq!(
        thin.input.material.identity(),
        old.input.material.identity()
    );
    assert_ne!(before.specimen_identity(), after.specimen_identity());
    assert_ne!(
        valve(&original).viscosity_matrix_n_s_m(),
        valve(&changed).viscosity_matrix_n_s_m()
    );
    assert_eq!(
        valve(&original).material_plates()[1]
            .source()
            .specimen_identity(),
        valve(&changed).material_plates()[1]
            .source()
            .specimen_identity()
    );
    let records = changed
        .material_records()
        .iter()
        .filter(|r| r.blade() == 0)
        .collect::<Vec<_>>();
    assert_eq!(
        records.iter().map(|r| r.section()).collect::<Vec<_>>(),
        vec![0, 0, 1]
    );
    assert_eq!(
        records.iter().map(|r| r.region()).collect::<Vec<_>>(),
        vec![0, 1, 2]
    );
    assert!(
        changed
            .material_records()
            .iter()
            .all(|r| std::ptr::eq(r.pack(), changed.material_records()[0].pack()))
    );
    let baseline = render(original).0;
    let tapered_pressure = render(changed).0;
    let replacement = render(build(
        &tapered.replacen("section blade.matpack", "section altered.matpack", 1),
        37,
    ))
    .0;
    assert!(tapered_pressure.iter().all(|p| p.is_finite()));
    assert!(
        baseline
            .iter()
            .zip(&tapered_pressure)
            .any(|(a, b)| (a - b).abs() > 1e-15)
    );
    assert!(
        replacement
            .iter()
            .zip(&tapered_pressure)
            .any(|(a, b)| (a - b).abs() > 1e-15)
    );
}

#[test]
fn more_than_sixty_four_thicknesses_share_one_authored_material_section() {
    let mut text = input();
    let (nx, ny) = (10, 4);
    let triangles = 2 * nx * ny;
    let mut geometry = format!(
        "sections 1\nsection blade.matpack isotropic 293.15 100 0.0003 0\nnodes {}\n",
        (nx + 1) * (ny + 1)
    );
    for y in 0..=ny {
        for x in 0..=nx {
            geometry.push_str(&format!(
                "node {} {} {} {}\n",
                0.025 * x as f64 / nx as f64,
                0.01 * y as f64 / ny as f64,
                usize::from(x == 0),
                usize::from(x == nx && y == ny / 2)
            ));
        }
    }
    geometry.push_str(&format!("triangles {triangles}\n"));
    let mut element = 0;
    for y in 0..ny {
        for x in 0..nx {
            let a = y * (nx + 1) + x;
            let b = a + 1;
            let c = a + nx + 1;
            let d = c + 1;
            for [a, b, c] in [[a, b, d], [a, d, c]] {
                let h = 0.0003 * (1.0 - 0.2 * (element + 1) as f64 / triangles as f64);
                geometry.push_str(&format!("triangle {a} {b} {c} 0 thickness {h}\n"));
                element += 1;
            }
        }
    }
    geometry.push_str(&format!("slit_edges {ny}\n"));
    for y in 0..ny {
        let a = y * (nx + 1) + nx;
        geometry.push_str(&format!("edge {a} {}\n", a + nx + 1));
    }
    let start = text.find("sections 2\n").unwrap();
    let end = start + text[start..].find("initial 3\n").unwrap();
    text.replace_range(start..end, &geometry);
    let mut loads = 0;
    let source = MaterialModalValvePerformance::from_bytes_with_sources(
        text.as_bytes(),
        37,
        &CancelGate::new_clock_free(),
        None,
        |path| {
            assert_eq!(path, "blade.matpack");
            loads += 1;
            Ok(pack(900.0, false))
        },
    )
    .unwrap();
    assert_eq!(loads, 1);
    let regions = valve(&source).material_plates()[0].source().regions();
    assert_eq!(regions.len(), triangles);
    assert!(regions.len() > 64);
    let records = source
        .material_records()
        .iter()
        .filter(|r| r.blade() == 0)
        .collect::<Vec<_>>();
    assert_eq!(records.len(), triangles);
    for (region, record) in records.iter().enumerate() {
        assert_eq!(record.section(), 0);
        assert_eq!(record.region(), region);
        assert!(std::ptr::eq(
            record.pack(),
            source.material_records()[0].pack()
        ));
        assert_eq!(
            regions[record.region()].input.region.triangle_indices.len(),
            1
        );
    }
}

#[test]
fn tapered_regions_keep_original_maxwell_sources_in_the_dg_consumer() {
    use super::super::dg::MaterialModalDgValvePerformance;
    const DG: &str =
        include_str!("../../../../../../examples/material_modal_valve/paired-dg.performance");
    let text = mapped(DG, |blade, element, h| {
        if blade == 0 && element == 0 {
            h * 0.5
        } else {
            h
        }
    })
    .replace(
        TENSOR,
        "loss regional-maxwell relaxed 1 boundary_modulus boundary_time",
    );
    let source = MaterialModalDgValvePerformance::from_bytes_with_sources(
        text.as_bytes(),
        37,
        &CancelGate::new_clock_free(),
        None,
        |_| Ok(pack(900.0, false)),
    )
    .unwrap();
    let LoadedAperture::Modal(valve) = source.renderer().system().aperture() else {
        panic!("authored material DG source must retain the modal material owner");
    };
    let memory = valve.material_relaxation().unwrap();
    assert!(memory.is_regional());
    assert_eq!(memory.material_spectra()[0].len(), 3);
    assert_eq!(memory.material_spectra()[1].len(), 2);
    for spectrum in memory.material_spectra().iter().flatten() {
        assert_eq!(spectrum.terms, vec![(1e8, 0.001)]);
    }
    for record in source.material_records() {
        let region = &valve.material_plates()[record.blade()].source().regions()[record.region()];
        assert!(record.section() < 2);
        assert_eq!(
            region.input.material.card_identity(),
            record.pack().card().content_hash()
        );
    }
}

#[test]
fn invalid_element_thickness_or_geometry_suffix_is_refused() {
    for suffix in [
        "thickness 0",
        "thickness -0.0001",
        "thickness NaN",
        "thickness inf",
        "thickness",
        "thickness 0.0002 extra",
        "height 0.0002",
    ] {
        let text = input().replacen("triangle 0 1 8 0", &format!("triangle 0 1 8 0 {suffix}"), 1);
        assert!(
            MaterialModalValvePerformance::from_bytes_with_sources(
                text.as_bytes(),
                37,
                &CancelGate::new_clock_free(),
                None,
                |_| Ok(pack(900.0, false))
            )
            .is_err(),
            "{suffix}"
        );
    }
}
