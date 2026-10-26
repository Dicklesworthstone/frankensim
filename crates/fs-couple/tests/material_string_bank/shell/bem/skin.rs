//! A declared acoustic skin carries the same Galerkin air load as the shell
//! path, on a prefix of the body modes, with a lumped neck's mouths as caps.
//! These synthetic fixtures exercise the owners; no instrument claim.
use super::*;
use fs_couple::air_path::cylinder::CylinderFlowResolution;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::observation::one_way::{
    ShellAcousticSkin, SkinAnchor, SkinAperture, skin_identity,
};
use fs_couple::bowed_string::runtime::reciprocal::bank::radiation::bem::ShellBemReport;
use fs_couple::bowed_string::runtime::reciprocal::plate::MaterialBodyCavity;
use fs_couple::cavity_modal::{CavityNeckLoss, HelmholtzModalCavity};
use fs_couple::render::plate::impact::radiation::Model;

/// The angular level of the tests that load `resolution()`'s whole band. The
/// owner's chart admits k times the box diagonal up to 4 * 2^level: 630 Hz at
/// level 0 on this 0.2 m icosphere, short of the band's 750 Hz; level 1 admits
/// it, under the 857 Hz slab chart.
const LEVEL: u32 = 1;

/// The skin that is the closed body itself: its nodes, its triangles in the
/// shell path's order, each vertex anchored at one corner of a triangle.
fn coincident_skin(body: &SharedMaterialShellBody) -> ShellAcousticSkin {
    assert!(
        body.specimen.pressure_surfaces().is_none(),
        "the shell path reads every facet, positive side out"
    );
    let mesh = body.specimen.mesh();
    let mut anchors = vec![Vec::new(); mesh.nodes.len()];
    for (t, tri) in mesh.tris.iter().enumerate() {
        for (corner, &v) in tri.iter().enumerate() {
            if anchors[v].is_empty() {
                let mut barycentric = [0.0; 3];
                barycentric[corner] = 1.0;
                anchors[v].push(SkinAnchor {
                    weight: 1.0,
                    component: 0,
                    triangle: t,
                    barycentric,
                });
            }
        }
    }
    ShellAcousticSkin {
        vertices: mesh.nodes.clone(),
        anchors,
        triangles: mesh.tris.clone(),
        apertures: Vec::new(),
    }
}

/// Two caps, the triangle fans around the +y and -y poles, sharing one flow.
fn capped(mut skin: ShellAcousticSkin) -> ShellAcousticSkin {
    let fan = |sign: f64| {
        let pole = (0..skin.vertices.len())
            .max_by(|&a, &b| (sign * skin.vertices[a][1]).total_cmp(&(sign * skin.vertices[b][1])))
            .unwrap();
        (0..skin.triangles.len())
            .filter(|&t| skin.triangles[t].contains(&pole))
            .collect::<Vec<_>>()
    };
    let (upper, lower) = (fan(1.0), fan(-1.0));
    skin.apertures = vec![
        SkinAperture {
            id: "upper".into(),
            triangles: upper,
            flow_share: 0.6,
        },
        SkinAperture {
            id: "lower".into(),
            triangles: lower,
            flow_share: 0.4,
        },
    ];
    skin
}

/// The closed body with a second free nodal translation, at a five-triangle
/// vertex on the far side from the first: two modes, not one.
fn two_mode_body(gate: &CancelGate) -> SharedMaterialShellBody {
    let mut body = closed_body(700.0, gate);
    let mesh = body.specimen.mesh();
    let mut fan = vec![0_usize; mesh.nodes.len()];
    for tri in &mesh.tris {
        for &v in tri {
            fan[v] += 1;
        }
    }
    let second = (0..mesh.nodes.len())
        .filter(|&v| fan[v] == 5)
        .min_by(|&a, &b| mesh.nodes[a][0].total_cmp(&mesh.nodes[b][0]))
        .unwrap();
    assert!(
        body.constraints[second].mask[0],
        "the second translation was free already"
    );
    body.constraints[second].mask[0] = false;
    body.n_modes = 2;
    body.projection_budget.max_modes = 2;
    body
}

fn declared() -> fs_scenario::HelmholtzCavity {
    fs_scenario::HelmholtzCavity {
        volume_m3: 4.0e-3,
        neck_radius_m: 0.004,
        neck_length_m: 0.03,
        inner_end: fs_scenario::NeckEnd::Flanged,
        outer_end: fs_scenario::NeckEnd::Flanged,
    }
}

fn cavity(loss: CavityNeckLoss) -> MaterialBodyCavity {
    MaterialBodyCavity {
        cavity: HelmholtzModalCavity::new(declared(), gas(101_325.0), loss).unwrap(),
        modal_budget: ModalAcousticTimeBudget::audible_reference(),
        listener_m: 1.0,
    }
}

fn bank(
    cx: &Cx<'_>,
    gate: &CancelGate,
    body: SharedMaterialShellBody,
    cavity: Option<MaterialBodyCavity>,
) -> ReciprocalStringBank {
    let mut settings = options(0);
    settings.network.max_modes = 64;
    settings.network.max_connections = 64;
    settings.network.max_setup_terms = 1_048_576;
    ReciprocalStringBank::from_material_shear_shell(
        cx,
        vec![member(7800.0).source],
        gas(101_325.0),
        body,
        cavity,
        law(),
        settings,
        gate,
    )
    .unwrap()
}

fn body_hz(bank: &ReciprocalStringBank) -> Vec<f64> {
    bank.system().components()[1]
        .modes()
        .iter()
        .map(|m| m.angular_frequency_rad_s / core::f64::consts::TAU)
        .collect()
}

fn model_bits(model: &Model) -> (usize, Vec<(u64, u64, Vec<u64>)>) {
    let poles = model
        .poles
        .iter()
        .map(|p| {
            (
                p.omega.to_bits(),
                p.zeta.to_bits(),
                p.coupling.iter().map(|c| c.to_bits()).collect(),
            )
        })
        .collect();
    (model.ports, poles)
}

fn refuses(result: Result<(Model, ShellBemReport), ShellBemError>, naming: &str) {
    match result {
        Err(ShellBemError::Invalid(message)) => {
            assert!(message.contains(naming), "refused, but with: {message}");
        }
        other => panic!(
            "expected the refusal naming {naming:?}, got {:?}",
            other.map(|_| ())
        ),
    }
}

/// 32 released steps against the same bank without the load: the work balance
/// closes, the pressure differs and the load stores energy.
fn closes_against(mut bare: ReciprocalStringBank, loaded: &mut ReciprocalStringBank) {
    let mut force = vec![0.0; loaded.system().mode_count()];
    force[4] = 0.02;
    let released = FrictionDrive {
        speed_m_s: 0.0,
        normal_force_n: 0.0,
    };
    let mut pressure_difference = 0.0;
    for _ in 0..32 {
        let a = *bare.step(&force, released).unwrap();
        let b = *loaded.step(&force, released).unwrap();
        assert!(b.energy_residual_j.abs() <= b.energy_tolerance_j);
        pressure_difference += (a.observer_pressure_pa - b.observer_pressure_pa).abs();
    }
    assert!(pressure_difference > 1e-12);
    assert!(loaded.system().radiation_energy_j().unwrap() > 0.0);
}

#[test]
fn a_skin_on_the_shell_itself_carries_the_shell_load_bit_for_bit() {
    with_cx(|cx, gate| {
        let skin = coincident_skin(&closed_body(700.0, gate));
        let (shell, a) = closed_bank(cx, gate, 700.0, 101_325.0)
            .prepare_material_shell_bem_with_angular_refinement(resolution(), LEVEL, gate)
            .unwrap();
        // An exact binding: every vertex is its anchored node, at distance 0.
        let (on_skin, b) = closed_bank(cx, gate, 700.0, 101_325.0)
            .prepare_material_skin_bem(resolution(), &skin, 0.0, LEVEL, gate)
            .unwrap();
        assert_eq!(a.skin_identity, None);
        assert_eq!(b.skin_identity, Some(skin_identity(&skin)));
        assert_eq!(
            (b.body_modes, b.one_way_body_modes, b.opening_ports),
            (1, 0, 0)
        );
        assert_eq!(b.basis_identity, a.basis_identity);
        assert_eq!(
            (b.surface_panels, b.reciprocal.unknowns, b.sweep_terms),
            (a.surface_panels, a.reciprocal.unknowns, a.sweep_terms)
        );
        // The anchored corner translations and the shell's facet mean are the
        // same sums, so the source, its fit and its gates agree to the bit.
        assert_eq!(model_bits(&on_skin), model_bits(&shell));
        for (x, y) in [
            (
                a.reciprocity_relative_residual,
                b.reciprocity_relative_residual,
            ),
            (a.peak_error, b.peak_error),
            (a.rms_error, b.rms_error),
            (a.resistance_peak_error, b.resistance_peak_error),
        ] {
            assert_eq!(x.to_bits(), y.to_bits());
        }
        let mut moved = skin.clone();
        moved.vertices[0][0] += 1e-9;
        refuses(
            closed_bank(cx, gate, 700.0, 101_325.0).prepare_material_skin_bem(
                resolution(),
                &moved,
                0.0,
                LEVEL,
                gate,
            ),
            "farther from its anchored point",
        );
    });
}

#[test]
fn a_band_below_the_upper_mode_loads_the_prefix_and_leaves_the_rest_one_way() {
    with_cx(|cx, gate| {
        let skin = coincident_skin(&two_mode_body(gate));
        let build = || bank(cx, gate, two_mode_body(gate), None);
        let hz = body_hz(&build());
        assert_eq!(hz.len(), 2);
        assert!(
            resolution().minimum_hz < hz[0] && hz[0] < hz[1] && hz[0] < resolution().maximum_hz,
            "the fixture needs two distinct modes, the lower inside the band: {hz:?}"
        );
        let mut band = resolution();
        band.maximum_hz = resolution().maximum_hz.min(0.5 * (hz[0] + hz[1]));
        let (mut loaded, report) = build()
            .with_material_skin_bem(band, &skin, 0.0, 0, gate)
            .unwrap();
        assert_eq!(
            (
                report.body_modes,
                report.one_way_body_modes,
                report.opening_ports
            ),
            (1, 1, 0)
        );
        assert!(report.reciprocity_relative_residual <= 1e-12);
        let binding = loaded.source().body_radiation().unwrap();
        assert_eq!(binding.source_body_coordinates(), &[0]);
        assert_eq!(binding.model().ports, 1);
        let retained = binding.material_shell_bem().unwrap();
        assert_eq!(retained.loading_identity, report.loading_identity);
        assert_eq!(retained.skin_identity, Some(skin_identity(&skin)));
        assert_eq!(retained.one_way_body_modes, 1);
        closes_against(build(), &mut loaded);
        // No loaded mode may lie below the band, and the prefix cannot be empty.
        let mut late = band;
        late.minimum_hz = 0.5 * (hz[0] + band.maximum_hz);
        refuses(
            build().prepare_material_skin_bem(late, &skin, 0.0, 0, gate),
            "ascending prefix",
        );
        let mut low = band;
        low.maximum_hz = 0.5 * (band.minimum_hz + hz[0]);
        refuses(
            build().prepare_material_skin_bem(low, &skin, 0.0, 0, gate),
            "ascending prefix",
        );
    });
}

#[test]
fn a_skin_prefix_load_survives_both_orders_of_string_gas_loading() {
    with_cx(|cx, gate| {
        let skin = coincident_skin(&two_mode_body(gate));
        let build = || bank(cx, gate, two_mode_body(gate), None);
        let hz = body_hz(&build());
        let mut band = resolution();
        band.maximum_hz = resolution().maximum_hz.min(0.5 * (hz[0] + hz[1]));
        let flow = CylinderFlowResolution {
            minimum_wavenumber_radius: 0.01,
            maximum_wavenumber_radius: 4.0,
            intervals: 2,
        };
        // The first rebuild restores the stored prefix map [0], not the full [0, 1].
        let mut first = build()
            .with_material_skin_bem(band, &skin, 0.0, 0, gate)
            .unwrap()
            .0
            .with_cylinder_gas_loading(flow, gate)
            .unwrap();
        let mut second = build()
            .with_cylinder_gas_loading(flow, gate)
            .unwrap()
            .with_material_skin_bem(band, &skin, 0.0, 0, gate)
            .unwrap()
            .0;
        assert_eq!(first.system().mode_count(), second.system().mode_count());
        for loaded in [&first, &second] {
            let binding = loaded.source().body_radiation().unwrap();
            assert_eq!(binding.source_body_coordinates(), &[0]);
            assert_eq!(
                binding.material_shell_bem().unwrap().skin_identity,
                Some(skin_identity(&skin))
            );
        }
        let mut force = vec![0.0; first.system().mode_count()];
        force[4] = 0.02;
        let released = FrictionDrive {
            speed_m_s: 0.0,
            normal_force_n: 0.0,
        };
        for _ in 0..24 {
            assert_eq!(
                first.step(&force, released).unwrap(),
                second.step(&force, released).unwrap()
            );
        }
    });
}

#[test]
fn an_exterior_owned_mouth_radiates_through_its_caps_and_every_other_pairing_refuses() {
    with_cx(|cx, gate| {
        let plain = coincident_skin(&closed_body(700.0, gate));
        let caps = capped(plain.clone());
        let exterior = || {
            bank(
                cx,
                gate,
                closed_body(700.0, gate),
                Some(cavity(CavityNeckLoss::ExteriorMouth)),
            )
        };
        let (mut loaded, report) = exterior()
            .with_material_skin_bem(resolution(), &caps, 0.0, LEVEL, gate)
            .unwrap();
        assert_eq!(
            (
                report.body_modes,
                report.one_way_body_modes,
                report.opening_ports
            ),
            (1, 0, 1)
        );
        assert!(report.reciprocity_relative_residual <= 1e-12);
        let binding = loaded.source().body_radiation().unwrap();
        // The neck is the first coordinate of the component after the body.
        assert_eq!(binding.source_body_coordinates(), &[0, 1]);
        assert_eq!(binding.model().ports, 2);
        let z = binding
            .model()
            .impedance(core::f64::consts::TAU * 150.0)
            .unwrap();
        assert!((z[1] - z[2]).abs() <= 1e-12 * (1.0 + z[1].abs()));
        assert!(
            z[1].abs() > 0.0,
            "the wall and the caps react on each other"
        );
        assert!(
            z[3].re > 0.0 && z[3].im < 0.0,
            "the caps carry the outer end's radiation and exterior mass: {:?}",
            z[3]
        );
        closes_against(exterior(), &mut loaded);
        // The shell path has no caps for the neck's port.
        refuses(
            exterior().prepare_material_shell_bem(resolution(), gate),
            "caps of a declared skin",
        );
        refuses(
            bank(
                cx,
                gate,
                closed_body(700.0, gate),
                Some(cavity(CavityNeckLoss::Resistance(100.0))),
            )
            .prepare_material_skin_bem(resolution(), &caps, 0.0, LEVEL, gate),
            "declare it exterior-owned",
        );
        refuses(
            exterior().prepare_material_skin_bem(resolution(), &plain, 0.0, LEVEL, gate),
            "radiates through declared skin apertures",
        );
        refuses(
            closed_bank(cx, gate, 700.0, 101_325.0).prepare_material_skin_bem(
                resolution(),
                &caps,
                0.0,
                LEVEL,
                gate,
            ),
            "need a lumped cavity neck",
        );
        // A band that starts above the neck's own resonance cannot carry its mouth.
        let atmosphere = gas(101_325.0);
        let geometry = declared();
        let neck_hz = fs_phs::helmholtz_natural_frequency(
            geometry.volume_m3,
            geometry.neck_radius_m,
            geometry.neck_length_m,
            atmosphere.density,
            atmosphere.sound_speed,
            fs_couple::cavity_modal::neck_ends(&geometry),
        )
        .unwrap()
            / core::f64::consts::TAU;
        let hz = body_hz(&exterior());
        assert!(
            neck_hz < 0.9 * hz[0],
            "the fixture needs the neck below the body mode: {neck_hz} vs {hz:?}"
        );
        let mut late = resolution();
        late.minimum_hz = hz[0];
        refuses(
            exterior().prepare_material_skin_bem(late, &caps, 0.0, LEVEL, gate),
            "reach the cavity neck resonance",
        );
    });
}
