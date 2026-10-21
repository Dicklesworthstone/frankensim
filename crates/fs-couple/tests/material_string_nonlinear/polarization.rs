//! Physical two-plane invariants on the existing nonlinear owner.
//! Synthetic material claims isolate wiring; they are not violin measurements.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::polarization::TransverseDirections;
use fs_couple::bowed_string::runtime::reciprocal::bank::nonlinear::{
    NonlinearMaterialString, supported::SupportedMaterialStrings,
};
use fs_couple::render::plate::impact::supported::TranslatingSupport;
use fs_couple::render::plate::impact::audio::ImpactSource;

fn directions(bow: [f64; 2], bridge: [f64; 2]) -> TransverseDirections {
    TransverseDirections::try_new(bow, bridge).unwrap()
}
fn states(x: [f64; 3], y: [f64; 3]) -> Vec<ModalAcousticState> {
    let norm = (7800.0 * core::f64::consts::PI * 0.0003_f64.powi(2) * 0.5 * 0.5).sqrt();
    x.into_iter().chain(y).map(|q| ModalAcousticState {
        displacement_m_sqrt_kg: q * norm, velocity_m_sqrt_kg_per_s: 0.0,
    }).collect()
}
fn rotate(q: &[f64]) -> Vec<f64> {
    let n = q.len() / 2;
    let mut r = q.to_vec();
    for i in 0..n {
        r[i] = 0.6 * q[i] - 0.8 * q[n+i];
        r[n+i] = 0.8 * q[i] + 0.6 * q[n+i];
    }
    r
}
fn near(a: f64, b: f64, floor: f64) {
    assert!((a-b).abs() <= floor + 2e-9 * a.abs().max(b.abs()), "{a:e} != {b:e}");
}
fn potential(image: &NonlinearMaterialString) -> &fs_couple::render::plate::impact::string::StringPotential {
    let BodyPotential::String(s) = &image.bodies()[0].potential else {
        panic!("original nonlinear string storage was lost")
    };
    s
}

#[test]
fn both_planes_share_one_spatial_elongation_and_material_tension() {
    with_cx(|cx, gate| {
        let physical_x = [0.002, -0.0004, 0.0002];
        let physical_y = [0.0015, 0.0003, -0.0001];
        let state = states(physical_x, physical_y);
        let mut energies = Vec::new();
        for young in [50e9, 200e9] {
            let image = source(7800.0, young, 3)
                .with_transverse_polarizations(directions([0.6, 0.8], [1.0, 0.0])).unwrap()
                .into_nonlinear(cx, gas(), 48_000, 0.2, state.clone(), gate).unwrap();
            assert_eq!(image.mode_count(), 6);
            assert_eq!(image.transverse_mode_count(), 6);
            assert_eq!(image.bodies().len(), 1);
            assert_eq!(image.bodies()[0].initial, state);
            let s = potential(&image);
            assert_eq!(s.transverse_planes(), 2);
            assert_eq!(&s.omegas()[..3], &s.omegas()[3..]);
            assert_eq!(&image.bodies()[0].damping_per_s[..3], &image.bodies()[0].damping_per_s[3..]);
            let q: Vec<_> = state.iter().map(|s| s.displacement_m_sqrt_kg).collect();
            let observed = s.observe(&q).unwrap();
            // Independent integral of the two physical continuous sine fields.
            let mut strain = 0.0;
            let cells = 1024;
            for cell in 0..cells {
                let x = 0.5 * (cell as f64 + 0.5) / cells as f64;
                let mut slope = [0.0; 2];
                for k in 0..3 {
                    let wave = (k+1) as f64 * core::f64::consts::PI / 0.5;
                    slope[0] += physical_x[k] * wave * (wave*x).cos();
                    slope[1] += physical_y[k] * wave * (wave*x).cos();
                }
                strain += 0.5 * (slope[0]*slope[0] + slope[1]*slope[1]) / cells as f64;
            }
            near(observed.additional_strain, strain, 1e-16);
            near(observed.tension_n, 20.0 + image.axial_rigidity_n()*strain, 1e-12);
            near(observed.stretching_energy_j, 0.5*image.axial_rigidity_n()*0.5*strain*strain, 1e-15);
            let mut x_only = q.clone(); x_only[3..].fill(0.0);
            let mut y_only = q.clone(); y_only[..3].fill(0.0);
            let independent = s.observe(&x_only).unwrap().stretching_energy_j
                + s.observe(&y_only).unwrap().stretching_energy_j;
            assert!(observed.stretching_energy_j > 1.2 * independent);
            // Motion in one plane changes the other plane's restoring tangent.
            let mut coupled = [0.0; 6];
            s.hessian_vector(&q, &[0.0, 0.0, 0.0, 1.0, 0.0, 0.0], &mut coupled);
            assert!(coupled[0].abs() > 1.0);
            energies.push(observed.stretching_energy_j);
        }
        near(energies[1], 4.0*energies[0], 1e-15);
    });
}

#[test]
fn potential_force_tangent_and_continuous_slope_admission_are_rotation_covariant() {
    with_cx(|cx, gate| {
        let image = source(7800.0, 200e9, 3)
            .with_transverse_polarizations(directions([1.0, 0.0], [0.6, 0.8])).unwrap()
            .into_nonlinear(cx, gas(), 48_000, 0.2, states([0.0;3], [0.0;3]), gate).unwrap();
        let s = potential(&image);
        let q = [6e-5, -2e-5, 1e-5, 4e-5, 1e-5, -7e-6];
        let d = [0.2, -0.1, 0.3, -0.3, 0.4, 0.1];
        let rq = rotate(&q); let rd = rotate(&d);
        near(s.potential(&q), s.potential(&rq), 1e-15);
        let a = s.observe(&q).unwrap(); let b = s.observe(&rq).unwrap();
        near(a.slope_bound, b.slope_bound, 1e-15);
        near(a.additional_strain, b.additional_strain, 1e-16);
        near(a.tension_n, b.tension_n, 1e-12);
        let (mut g, mut rg, mut h, mut rh) = ([0.0;6], [0.0;6], [0.0;6], [0.0;6]);
        s.gradient(&q, &mut g); s.gradient(&rq, &mut rg);
        s.hessian_vector(&q, &d, &mut h); s.hessian_vector(&rq, &rd, &mut rh);
        for (a,b) in rotate(&g).iter().zip(rg) { near(*a, b, 1e-11); }
        for (a,b) in rotate(&h).iter().zip(rh) { near(*a, b, 1e-7); }
        let eps = 1e-9;
        let plus: Vec<_> = q.iter().zip(d).map(|(q,d)| q+eps*d).collect();
        let minus: Vec<_> = q.iter().zip(d).map(|(q,d)| q-eps*d).collect();
        let (mut gp, mut gm) = ([0.0;6], [0.0;6]);
        s.gradient(&plus, &mut gp); s.gradient(&minus, &mut gm);
        for i in 0..6 {
            assert!((h[i] - (gp[i]-gm[i])/(2.0*eps)).abs() < 2e-6*h[i].abs().max(1.0));
        }
        // This threshold admits the vector slope but would reject a sum of
        // absolute plane amplitudes after a basis rotation.
        let state = states([0.003,0.0,0.0], [0.0;3]);
        let limit = 1.01 * core::f64::consts::PI * 0.003 / 0.5;
        let q: Vec<_> = state.iter().map(|s| s.displacement_m_sqrt_kg).collect();
        let rotated: Vec<_> = rotate(&q).into_iter().map(|q| ModalAcousticState {
            displacement_m_sqrt_kg: q, velocity_m_sqrt_kg_per_s: 0.0,
        }).collect();
        for state in [state, rotated] {
            assert!(source(7800.0, 200e9, 3)
                .with_transverse_polarizations(directions([1.0,0.0],[1.0,0.0])).unwrap()
                .into_nonlinear(cx, gas(), 48_000, limit, state, gate).is_ok());
        }
    });
}

#[test]
fn original_nonlinear_steps_preserve_rotated_motion_and_shared_port_work() {
    with_cx(|cx, gate| {
        let state = states([0.002,-0.0004,0.0002], [0.0015,0.0003,-0.0001]);
        let q: Vec<_> = state.iter().map(|s| s.displacement_m_sqrt_kg).collect();
        let rotated: Vec<_> = rotate(&q).into_iter().map(|q| ModalAcousticState {
            displacement_m_sqrt_kg: q, velocity_m_sqrt_kg_per_s: 0.0,
        }).collect();
        let make = |state, d| source(7800.0, 200e9, 3)
            .with_transverse_polarizations(d).unwrap()
            .into_nonlinear(cx, gas(), 48_000, 0.2, state, gate).unwrap()
            .prepare(vec![], vec![], vec![], vec![], config(), gate).unwrap();
        let mut a = make(state, directions([1.0,0.0],[0.6,0.8]));
        let mut b = make(rotated, directions([0.6,0.8],[-0.28,0.96]));
        let row_a = a.source().transverse_shape(0.27).unwrap();
        let row_b = b.source().transverse_shape(0.27).unwrap();
        let initial = a.system().stored_energy_j();
        let mut net = 0.0;
        for step in 0..192 {
            let force = if step < 48 { 0.005 } else { 0.0 };
            let fa: Vec<_> = row_a.iter().map(|v| force*v).collect();
            let fb: Vec<_> = row_b.iter().map(|v| force*v).collect();
            let af = a.step(&fa, gate).unwrap(); let bf = b.step(&fb, gate).unwrap();
            net += af.supplied_work_j - af.dissipated_energy_j;
            assert!(af.dissipated_energy_j >= -1e-15);
            assert!((af.stored_energy_j-initial-net).abs() < 1e-8);
            near(af.stored_energy_j, bf.stored_energy_j, 1e-10);
            near(af.supplied_work_j, bf.supplied_work_j, 1e-10);
            for offset in 0..2 {
                let qa: Vec<_> = (0..6).map(|i| a.state()[2*i+offset]).collect();
                for (i,expected) in rotate(&qa).iter().enumerate() {
                    near(*expected, b.state()[2*i+offset], 2e-9);
                }
            }
        }
        assert_eq!(a.samples_rendered(), 192);
        assert_eq!(a.source().mode_count(), 6);
        assert_eq!(a.system().string_observation(0).unwrap().tension_n,
            potential(a.source()).observe(&(0..6).map(|i|a.state()[2*i]).collect::<Vec<_>>()).unwrap().tension_n);
    });
}

#[test]
fn original_material_twist_and_directional_force_torque_ports_survive_the_shared_stretching_body() {
    with_cx(|cx, gate| {
        let source = source(7800.0, 200e9, 3).with_torsion("G", "eta_s", 2).unwrap()
            .with_transverse_polarizations(directions([0.6,0.8],[1.0,0.0])).unwrap();
        let image = source.into_nonlinear(cx, gas(), TORSION_RATE_HZ, 0.2,
            vec![ModalAcousticState::default();8], gate).unwrap();
        assert_eq!(image.transverse_mode_count(), 6);
        assert_eq!(image.mode_count(), 8);
        assert_eq!(image.bodies().len(), 2);
        assert_eq!(image.bodies()[1].initial.len(), 2);
        let normal = image.transverse_shape(0.27).unwrap();
        let bridge = image.bridge_shape(0.27).unwrap();
        let surface = image.surface_shape(0.27).unwrap();
        let angle = image.angular_shape(0.27).unwrap();
        assert_eq!(&normal[..6], &surface[..6]);
        assert!(bridge[3..].iter().all(|x| *x == 0.0));
        assert!(normal[6..].iter().all(|x| *x == 0.0));
        assert!(angle[..6].iter().all(|x| *x == 0.0));
        for i in 0..3 {
            near(normal[i], 0.6*bridge[i], 1e-14);
            near(normal[3+i], 0.8*bridge[i], 1e-14);
        }
        for i in 6..8 {
            assert!(angle[i].abs() > 0.0);
            near(surface[i], 0.0003*angle[i], 1e-12);
        }
        let mut run = image.prepare(vec![],vec![],vec![],vec![],config_at(TORSION_RATE_HZ),gate).unwrap();
        let load: Vec<_> = normal.iter().map(|x| 0.01*x).collect();
        for _ in 0..48 { run.step(&load,gate).unwrap(); }
        assert!(run.state()[..12].iter().any(|x| *x != 0.0));
        assert!(run.state()[12..].iter().all(|x| *x == 0.0));
        assert!(run.system().string_observation(0).is_some());
        assert!(run.system().string_observation(1).is_none());
    });
}

#[test]
fn directed_endpoint_support_retains_both_planes_and_wrong_clocks_still_refuse() {
    with_cx(|cx,gate| {
        let image = source(7800.0,200e9,3)
            .with_transverse_polarizations(directions([1.0,0.0],[1.0,0.0])).unwrap()
            .into_nonlinear(cx,gas(),48_000,0.2,vec![ModalAcousticState::default();6],gate).unwrap();
        let support = TranslatingSupport { mass_kg: 0.05, stiffness_n_m: 1000.0,
            damping_n_s_m: 0.1, initial_position_m: 0.0, initial_velocity_m_s: 0.0,
            maximum_travel_m: 0.01, maximum_slope: 0.2 };
        let supported = SupportedMaterialStrings::in_comoving_gas(vec![image.clone()], support, gate).unwrap();
        assert_eq!(supported.mode_count(),7);
        assert_eq!(supported.sources()[0].transverse_mode_count(),6);
        assert_eq!(supported.bodies()[0].initial.len(),7);
        assert!(image.clone().prepare(vec![],vec![],vec![],vec![],config_at(96_000),gate).is_err());
        assert!(image.angular_shape(0.27).is_err());
        for station in [0.0,1.0,-1.0,f64::NAN] {
            assert!(image.transverse_shape(station).is_err());
            assert!(image.bridge_shape(station).is_err());
            assert!(image.surface_shape(station).is_err());
        }
        let mut run = image.prepare(vec![],vec![],vec![],vec![],config(),gate).unwrap();
        let before = run.state().to_vec();
        let stop = fs_exec::CancelGate::new_clock_free(); stop.request();
        assert!(run.step(&[0.0;6],&stop).is_err());
        assert_eq!(run.state(), before);
        assert_eq!(run.samples_rendered(), 0);
    });
}
