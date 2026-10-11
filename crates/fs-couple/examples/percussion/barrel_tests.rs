use super::*;
use crate::{Mechanics, Stroke, cavity, mallets, shaft_playing};
use fs_couple::render::plate::impact::{ImpactSystem, VolumeSpring};
use fs_couple::render::plate::impact::cavity::cylinder::CylinderSpec;
use fs_couple::vibroacoustic::AcousticMedium;

const DT: f64 = 1e-7;
const CARD: &str = "frankensim-drum-barrel-v1\nmaterial,20000000,0.3,1000,0\naxial_intervals,2\nband_hz,0,100000\n";

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
fn inner_skin_pressure_excites_the_barrel_through_both_existing_cavity_images() {
    let spec = drum(); let barrel = Spec::read(CARD).unwrap().prepare(&spec, DT, false).unwrap();
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
        vec![], vec![], vec![], spec.radius_m, spec.depth_m, 32, DT, None, 0.0, Some((first,&barrel))).unwrap();
    let (prepared, same_pressure) = cavity::build_prepared_with_losses_and_barrel(&films, &modes, bodies,
        vec![], vec![], spec.radius_m, spec.depth_m, crate::mechanics::coupled_config(32,DT,false).unwrap(),
        None, 0.0, Some((first,&barrel))).unwrap();
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
        false,stroke,true,None,Some(drum),second,&[],0.0,None,false,None,&mallets::Selection::default(),&shafts,material).unwrap();
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
}
