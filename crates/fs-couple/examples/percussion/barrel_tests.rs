use super::*;
use crate::{Mechanics, Stroke, cavity, mallets, shaft_playing};
use fs_couple::render::plate::impact::{ImpactSystem, VolumeSpring};
use fs_couple::render::plate::impact::cavity::cylinder::CylinderSpec;
use fs_couple::vibroacoustic::AcousticMedium;

const DT: f64 = 1e-7;
const CARD: &str = "frankensim-drum-barrel-v1\nmaterial,20000000,0.3,1000,0\naxial_intervals,2\nband_hz,0,100000\n";
const ORTHOTROPIC_CARD: &str = "frankensim-drum-barrel-v2\northotropic,20000000,2000000,0.25,800000,1000,0,0\naxial_intervals,2\nband_hz,0,100000\n";

// A deliberately compliant small full-pencilled fixture, not a measured shell.
pub(super) fn drum() -> drum_spec::Spec {
    let mut spec = drum_spec::Spec { radial_intervals: 1, azimuths: 8,
        band_hz: [0.0, 2000.0], ..drum_spec::Spec::reference() };
    for head in &mut spec.heads { head.damping_ratio = 0.0; }
    spec
}

#[test]
fn barrel_requires_declared_material_mesh_and_window_and_refuses_foreign_commands() {
    let spec = Spec::read(CARD).unwrap();
    for family in ["drum", "drum-modal", "drum-stretch", "snare", "snare-off"] {
        for suffix in ["", "-wav", "-mic"] { admit_command(Some(&spec), &format!("{family}{suffix}")).unwrap(); }
    }
    for command in ["hihat", "splash-mic", "unknown"] { assert!(admit_command(Some(&spec), command).is_err()); }
    for bad in [CARD.replace("material,20000000,0.3,1000,0\n", ""),
        CARD.replace("20000000", "-1"), CARD.replace("1000,0", "1000,-0.1"),
        CARD.replace("axial_intervals,2", "axial_intervals,1"),
        CARD.replace("band_hz,0,100000", "band_hz,NaN,100000"),
        format!("{CARD}axial_intervals,2\n")] { assert!(Spec::read(&bad).is_err()); }
    assert!(spec.prepare(&drum(), DT, true).is_err(), "mechanical fixture window is not an audio claim");
    let mut args = vec!["--elastic-barrel".into()]; assert!(option(&mut args).is_err());
    let mut args = vec!["--elastic-barrel".into(), "x".into(), "--elastic-barrel".into(), "y".into()];
    assert!(option(&mut args).is_err());
}

#[test]
fn orthotropic_barrel_requires_all_plane_stress_constants_and_an_explicit_grain_angle() {
    let spec = Spec::read(ORTHOTROPIC_CARD).unwrap();
    assert_eq!(spec.material, Material::Orthotropic { young_1_pa:2e7, young_2_pa:2e6,
        poisson_12:0.25, shear_12_pa:8e5, grain_angle_rad:0.0 });
    assert!(Spec::read(include_str!("illustrative-orthotropic-barrel.fsb")).is_ok());
    for bad in [ORTHOTROPIC_CARD.replace("-v2", "-v1"), CARD.replace("-v1", "-v2"),
        ORTHOTROPIC_CARD.replace("20000000,2000000", "20000000,0"),
        ORTHOTROPIC_CARD.replace("0.25,800000", "4,800000"),
        ORTHOTROPIC_CARD.replace("0.25,800000", "0.25,-1"),
        ORTHOTROPIC_CARD.replace("1000,0,0", "1000,0,NaN"),
        ORTHOTROPIC_CARD.replace("1000,0,0", "1000,0,4"),
        ORTHOTROPIC_CARD.replace("1000,0,0", "1000,0"),
        format!("{ORTHOTROPIC_CARD}material,20000000,0.3,1000,0\n")] {
        assert!(Spec::read(&bad).is_err(), "{bad}");
    }
}

// Alternating triangle edge frames around the entire cylinder. An axial grain
// is one physical direction even when the element's first edge is diagonal.
fn grain_mesh() -> ShellMesh {
    let mut nodes = Vec::new();
    for z in [-0.1,0.1] { for i in 0..8 {
        let angle = core::f64::consts::TAU*i as f64/8.0;
        nodes.push([0.18*angle.cos(),0.18*angle.sin(),z]);
    }}
    let mut triangles = Vec::new();
    for i in 0..8 { let j = (i+1)%8; triangles.extend([[i,j,8+j],[i,8+j,8+i]]); }
    ShellMesh::new(nodes,triangles).unwrap()
}

#[test]
fn one_physical_grain_has_the_same_bending_energy_in_every_actual_facet_frame() {
    let mesh = grain_mesh(); let thickness = 0.006;
    for angle in [0.0,core::f64::consts::FRAC_PI_2,0.41,-0.41] {
        let spec = Spec { material:Material::Orthotropic { young_1_pa:2e7, young_2_pa:2e6,
            poisson_12:0.25, shear_12_pa:8e5, grain_angle_rad:angle },
            ..Spec::read(ORTHOTROPIC_CARD).unwrap() };
        let sections = spec.sections(&mesh,thickness).unwrap();
        let (s,c) = angle.sin_cos();
        // Rank-one curvature v tensor v, v = 0.37*axial + 0.91*hoop.
        // Independent principal-axis expression includes the engineering
        // shear factor four and distinguishes positive/negative grain angles.
        let v1 = 0.37*c+0.91*s; let v2 = 0.37*s-0.91*c;
        let scale = thickness.powi(3)/12.0;
        let denominator = 1.0-0.25*0.25*2e6/2e7;
        let expected = scale*((2e7*v1.powi(4)+2e6*v2.powi(4)
            +2.0*0.25*2e6*v1*v1*v2*v2)/denominator+4.0*8e5*v1*v1*v2*v2);
        for (facet,section) in sections.iter().enumerate() {
            let frame = mesh.facet(facet).unwrap().frame;
            let hoop = cross([0.0,0.0,1.0],frame[2]);
            let v:[f64;3] = core::array::from_fn(|i| 0.91*hoop[i]+if i==2 {0.37}else{0.0});
            let x = dot(v,frame[0]); let y = dot(v,frame[1]);
            let curvature = [x*x,y*y,2.0*x*y];
            let actual:f64 = (0..3).map(|i|(0..3).map(|j|
                curvature[i]*section.d[3*i+j]*curvature[j]).sum::<f64>()).sum();
            assert!((actual-expected).abs()<1e-11*expected, "facet {facet}, angle {angle}: {actual} != {expected}");
        }
    }
    let legacy = Spec::read(CARD).unwrap();
    let isotropic_limit = Spec { material:Material::Orthotropic { young_1_pa:2e7, young_2_pa:2e7,
        poisson_12:0.3, shear_12_pa:2e7/(2.0*1.3), grain_angle_rad:0.41 }, ..legacy };
    for (a,b) in legacy.sections(&mesh,thickness).unwrap().iter()
        .zip(isotropic_limit.sections(&mesh,thickness).unwrap()) {
        assert_eq!(a.d,b.d, "an isotropic material cannot acquire a mesh grain");
        assert_eq!(a.thickness,b.thickness); assert_eq!(a.density,b.density);
    }
}

#[test]
fn rotating_grain_changes_wall_modes_and_real_sidewall_compliance_without_changing_geometry() {
    let drum = drum(); let axial = Spec::read(ORTHOTROPIC_CARD).unwrap();
    let hoop = Spec { material:Material::Orthotropic { young_1_pa:2e7, young_2_pa:2e6,
        poisson_12:0.25, shear_12_pa:8e5, grain_angle_rad:core::f64::consts::FRAC_PI_2 }, ..axial };
    let a = axial.prepare(&drum,DT,false).unwrap(); let b = hoop.prepare(&drum,DT,false).unwrap();
    assert_eq!(a.mode_count(),b.mode_count());
    assert_eq!(a.outer_positions,b.outer_positions); assert_eq!(a.triangles,b.triangles);
    assert_eq!(a.inner_points,b.inner_points); assert_eq!(a.inner_areas,b.inner_areas);
    assert!(a.omegas.iter().zip(&b.omegas).any(|(x,y)|(x-y).abs()>1e-3*x.max(*y)));
    let strike = crate::barrel_strike::Spec { azimuth_rad:0.3,height_below_batter_m:0.4*drum.depth_m };
    let pa = strike.port(&a).unwrap(); let pb = strike.port(&b).unwrap();
    assert_eq!(pa.position_m,pb.position_m);
    // Exact static response of each retained linear image to a 1 N force at
    // the ACTUAL outer-skin point. Modal signs cannot change this compliance.
    let compliance = |wall:&Prepared, weights:&[f64]| -> f64 {
        weights.iter().zip(&wall.omegas).map(|(w,omega)|(w/omega).powi(2)).sum()
    };
    let ca = compliance(&a,&pa.weights); let cb = compliance(&b,&pb.weights);
    assert!(ca.is_finite() && ca>0.0 && cb.is_finite() && cb>0.0);
    assert!((ca-cb).abs()>1e-3*ca.max(cb), "grain must change the physical wall response: {ca}, {cb}");
    for wall in [&a,&b] {
        assert!(wall.compression.iter().any(|c|c.abs()>1e-8));
        assert!(wall.outer_weights.iter().flatten().any(|w|w.abs()>1e-8));
        assert!(wall.body().damping_per_s.iter().all(|d|*d==0.0));
    }
}

#[test]
fn inner_skin_pressure_excites_the_barrel_through_both_existing_cavity_images() {
    // The exact prepared image propagates the whole coupled network in one dense exponential
    // (at most 42 fully forced modes, .168). The card's full window holds 48 shell modes, up to
    // 668 Hz, beside 1 free mass, 2 head modes and 6 air modes. Keeping the shell to 460 Hz
    // retains 32 modes, including its breathing modes, and cuts below the 481.5 Hz degenerate pair.
    let spec = drum(); let card = CARD.replace("band_hz,0,100000", "band_hz,0,460");
    let barrel = Spec::read(&card).unwrap().prepare(&spec, DT, false).unwrap();
    let gate = CancelGate::new_clock_free();
    let air = CylindricalCavity::new(CylinderSpec { radius_m: spec.radius_m, depth_m: spec.depth_m,
        radial_intervals: 32, maximum_azimuthal_order: 4, maximum_axial_order: 2,
        maximum_frequency_hz: 1100.0, maximum_modes: 8, eigen_residual_tolerance: 1e-7 },
        AcousticMedium { rho0: 1.2, c0: 343.0 }, &gate).unwrap();
    let c = barrel.coupling(&air, &gate).unwrap(); let nc = air.modes().len();
    assert!(barrel.compression.iter().any(|a| a.abs() > 1e-8), "retain a genuine breathing mode");
    for (i, compression) in barrel.compression.iter().enumerate() {
        assert!((c[i*nc] + compression).abs() < 1e-12, "constant pressure must equal negative compression");
    }
    let (films, modes) = spec.prepare(DT, false).unwrap();
    let first = 1 + modes.iter().map(Vec::len).sum::<usize>();
    let total = first + barrel.mode_count();
    let mut areas = vec![0.0; total]; let mut bodies = vec![ImpactBody::free_mass(1.0, 0.0, 0.0).unwrap().0];
    let mut offset = 1;
    for head in 0..2 {
        let omega: Vec<_> = modes[head].iter().map(|m| m.lambda.sqrt()).collect();
        let mut body = crate::zero_body(BodyPotential::Linear(omega.clone()), &omega);
        body.damping_per_s.fill(0.0);
        for (i, mode) in modes[head].iter().enumerate() {
            areas[offset+i] = (if head==0 {1.0} else {-1.0}) * films[head].modal_area(&mode.phi).unwrap();
        }
        if head == 0 {
            let i = (0..omega.len()).max_by(|&a,&b| areas[offset+a].abs().total_cmp(&areas[offset+b].abs())).unwrap();
            // A displaced real head stores compression; no external force or
            // initial motion is assigned to ANY barrel coordinate.
            body.initial[i].displacement_m_sqrt_kg = 1e-7 / areas[offset+i];
        }
        offset += omega.len(); bodies.push(body);
    }
    bodies.push(barrel.body()); areas[first..].copy_from_slice(barrel.compression_areas());
    let volume = VolumeSpring { bulk_modulus_pa: 1.2*343.0*343.0, volume_m3: spec.volume_m3(), areas };
    let uniform = ImpactSystem::new(bodies.clone(), vec![], vec![], vec![volume], crate::config(32, DT)).unwrap();
    let (distributed, pressure) = cavity::build_with_pads_and_barrel(&films, &modes, bodies.clone(),
        vec![], vec![], vec![], spec.radius_m, spec.depth_m, 32, DT, None, 0.0, Some((first,&barrel)),
        AcousticMedium { rho0: 1.2, c0: 343.0 }).unwrap();
    let (prepared, same_pressure) = cavity::build_prepared_with_losses_and_barrel(&films, &modes, bodies,
        vec![], vec![], spec.radius_m, spec.depth_m, crate::mechanics::coupled_config(32,DT,false).unwrap(),
        None, 0.0, Some((first,&barrel)), AcousticMedium { rho0: 1.2, c0: 343.0 }).unwrap();
    assert_eq!(distributed.state(), prepared.state());
    assert_eq!(pressure.coupling.structural_modes(), total);
    assert_eq!(same_pressure.coupling.structural_modes(), total);
    assert!(pressure.uniform_pressure(distributed.state()).unwrap() > 0.0);
    let mut systems = [Mechanics::Reference(uniform).into_analytic_nonlinear().unwrap(),
        Mechanics::Reference(distributed).into_analytic_nonlinear().unwrap(), Mechanics::Prepared(prepared)];
    for system in &mut systems {
        let n = system.state().len()/2; let force = vec![0.0; n];
        assert!(system.state()[2*first..2*total].iter().all(|v| *v == 0.0));
        let mut initial = None;
        for tick in 0..32 {
            if tick == 16 {
                let before = system.state().to_vec(); let stop = CancelGate::new_clock_free(); stop.request();
                assert!(system.step(&force,&stop).is_err()); assert_eq!(system.state(), before);
            }
            let frame = system.step(&force,&gate).unwrap();
            let e = *initial.get_or_insert(frame.stored_energy_j);
            assert!(frame.balance_residual_j.abs() < 1e-7);
            assert!((frame.stored_energy_j-e).abs() < 1e-7);
        }
        let outward_flow: f64 = barrel.compression.iter().enumerate()
            .map(|(i,a)| -a*system.state()[2*(first+i)+1]).sum();
        assert!(outward_flow > 0.0, "positive pressure must expand the initially unforced barrel");
        assert!(system.state()[2*first..2*total].iter().any(|v| v.abs() > 1e-16));
    }
}

#[test]
fn full_wire_bank_and_flexible_shaft_keep_their_physical_addresses_around_barrel() {
    use fs_couple::render::plate::impact::striker::{RadiusStation, flexible::FlexibleStriker};
    use fs_plate::shell::stiffened::beam::RoundBeamSpec;
    let material = Spec::read(CARD).unwrap(); let drum = drum();
    let barrel_count = material.prepare(&drum,DT,false).unwrap().mode_count();
    let shafts = shaft_playing::Selection { first: Some(FlexibleStriker::new(
        &[(0.0,0.005),(0.4,0.005)].map(|(position_m,radius_m)| RadiusStation {position_m,radius_m}),
        RoundBeamSpec {young_pa:12e9,density_kg_m3:800.0,pivot_m:0.1,contact_m:0.39,hand_m:0.16,
            subdivisions:8,maximum_hz:3000.0,maximum_modes:17},0.001).unwrap()), second: None };
    let stroke = Stroke {speed_m_s:0.0,position_m:Some([0.06,0.01])};
    let second = Some(Stroke {speed_m_s:0.0,position_m:Some([-0.05,0.01])});
    let make = |material| crate::drum_with_barrel(1,DT,false,false,Some(crate::snare::SnareSet::reference(false)),
        false,stroke,true,None,Some(drum.clone()),second,&[],0.0,None,false,None,&mallets::Selection::default(),&shafts,material,None).unwrap();
    let rigid = make(None); let elastic = make(Some(&material));
    let old_shaft = rigid.flexible_sticks[0].as_ref().unwrap(); let shaft = elastic.flexible_sticks[0].as_ref().unwrap();
    assert_eq!(elastic.second_stick.unwrap().coordinate,rigid.second_stick.unwrap().coordinate);
    assert!(old_shaft.elastic_start() >= 162, "all 160 wire modes precede the shell and shaft");
    assert_eq!(shaft.elastic_start(),old_shaft.elastic_start()+barrel_count);
    assert_eq!(shaft.elastic_modes(),old_shaft.elastic_modes());
    let prefix = 2*old_shaft.elastic_start();
    assert_eq!(&elastic.system.state()[..prefix],&rigid.system.state()[..prefix]);
    assert_eq!(elastic.force.len(),rigid.force.len()+barrel_count);
    assert_eq!(elastic.air.as_ref().unwrap().coupling.structural_modes(),
        rigid.air.as_ref().unwrap().coupling.structural_modes()+barrel_count);
    assert!(elastic.force[old_shaft.elastic_start()..shaft.elastic_start()].iter().all(|v|*v==0.0));
    assert!(elastic.observer_a[old_shaft.elastic_start()..].iter().all(|v|*v==0.0));
    assert!(elastic.observer_b[old_shaft.elastic_start()..].iter().all(|v|*v==0.0));
    let side = crate::drum_with_barrel(1,DT,false,false,Some(crate::snare::SnareSet::reference(false)),
        false,Stroke {position_m:None,..stroke},true,None,Some(drum.clone()),second,&[],0.0,None,false,None,
        &mallets::Selection::default(),&shafts,Some(&material),Some(crate::barrel_strike::Spec {
            azimuth_rad:0.3,height_below_batter_m:0.4*drum.depth_m })).unwrap();
    assert_eq!(side.second_stick.unwrap().coordinate,elastic.second_stick.unwrap().coordinate);
    assert_eq!(side.flexible_sticks[0].as_ref().unwrap().elastic_start(),shaft.elastic_start());
    assert_eq!(side.system.state(),elastic.system.state(), "sidewall selection must not move the launch, wire bank or pressure state");
    assert_eq!(side.pressure.as_ref().unwrap().areas,elastic.pressure.as_ref().unwrap().areas);
}
