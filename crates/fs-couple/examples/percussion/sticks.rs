//! Two independent inertial strikers coupled to the SAME batter-head basis.
//! The first stick and both head ranges retain their original addresses. The
//! second stick is appended after the heads, before wires and cavity inertia.
//! No second drum, resampled impact, direct sound source or stick-stick law.
use super::{Error, ImpactBody, ModePair, Obstacle, Stroke, TensionedDisk};
use fs_couple::render::plate::impact::linear::wire::film_shapes;

#[derive(Clone, Copy, Debug)]
pub struct Port {
    pub coordinate: usize,
    /// Work-conjugate physical tip participation, 1/sqrt(kg).
    pub weight: f64,
}

/// An explicit XY station is mandatory; the second launch defaults to rest.
/// Extract only these controls, preserving the first stick's complete input.
pub fn option(args: &mut Vec<String>) -> Result<Option<Stroke>, Error> {
    let mut position = None;
    let mut speed = None;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--second-stick-position-m" => {
                if position.is_some() { return Err("duplicate second-stick position".into()); }
                let x: f64 = args.get(i+1).ok_or("second-stick position needs x and y")?.parse()?;
                let y: f64 = args.get(i+2).ok_or("second-stick position needs x and y")?.parse()?;
                if !x.is_finite() || !y.is_finite() { return Err("second-stick position must be finite metres".into()); }
                position = Some([x,y]);
                args.drain(i..i+3);
            }
            "--second-stick-speed-m-s" => {
                if speed.is_some() { return Err("duplicate second-stick launch speed".into()); }
                let value: f64 = args.get(i+1).ok_or("second-stick speed needs a value")?.parse()?;
                if !value.is_finite() || !(0.0..=20.0).contains(&value) {
                    return Err("second-stick speed must be finite in 0..=20 m/s".into());
                }
                speed = Some(value);
                args.drain(i..i+2);
            }
            _ => i += 1,
        }
    }
    match position {
        Some(p) => Ok(Some(Stroke { position_m: Some(p), speed_m_s: speed.unwrap_or(0.0) })),
        None if speed.is_some() => Err("second-stick speed requires --second-stick-position-m X Y".into()),
        None => Ok(None),
    }
}

pub fn admit_command(enabled: bool, command: &str) -> Result<(), Error> {
    if enabled && !matches!(command, "drum"|"drum-wav"|"drum-mic"|
        "drum-stretch"|"drum-stretch-wav"|"drum-stretch-mic"|
        "drum-modal"|"drum-modal-wav"|"drum-modal-mic"|
        "snare"|"snare-wav"|"snare-mic"|"snare-off"|"snare-off-wav"|"snare-off-mic") {
        return Err("two sticks require a drum/snare command; cymbal and stick-stick contact are not implemented".into());
    }
    Ok(())
}

/// Compile an independent tip against actual interpolated head displacement.
/// Its equal-and-opposite contact reactions enter the existing joint solve.
pub fn build(stroke: Stroke, film: &TensionedDisk, modes: &[ModePair],
    coordinate: usize, total: usize) -> Result<(ImpactBody, Obstacle, Port), Error>
{
    if modes.is_empty() || coordinate <= modes.len() || coordinate >= total {
        return Err("second-stick coordinate overlaps the original batter basis or exceeds the mechanical layout".into());
    }
    let position = stroke.position_m.ok_or("second stick requires an explicit physical XY station")?;
    let shapes = film_shapes(film, modes, &[position])?.remove(0);
    if shapes.iter().all(|b| *b == 0.0) {
        return Err("second-stick station has no retained moving-head participation; rim strikes are not implemented".into());
    }
    let (body, weight) = super::stick_with_speed(stroke.speed_m_s)?;
    let mut column = vec![0.0; total];
    column[coordinate] = weight;
    for (i, shape) in shapes.iter().enumerate() { column[1+i] = -shape; }
    let contact = super::elastic_contact(column)?;
    eprintln!("second physical stick: xy_m={position:?}, launch_m_s={}, coordinate={coordinate}; same estimated stick/contact cards, independent inertia, shared head; no direct gas/radiation participation",stroke.speed_m_s);
    Ok((body, contact, Port { coordinate, weight }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::{drum_with_sticks, drum_with_playing, mechanics::Mechanics};
    use fs_exec::CancelGate;
    fn args(s: &str) -> Vec<String> { s.split_whitespace().map(str::to_owned).collect() }
    fn primary() -> Stroke { Stroke { speed_m_s: 4.0, position_m: Some([0.06,0.01]) } }
    fn secondary() -> Stroke { Stroke { speed_m_s: 2.5, position_m: Some([-0.05,0.02]) } }

    #[test]
    fn independent_si_controls_preserve_the_original_stick_and_reject_ambiguity() {
        let mut a = args("drum-modal 192 --second-stick-speed-m-s 2.5 --strike-speed-m-s 4 --second-stick-position-m -0.05 0.02 --strike-position-m 0.06 0.01");
        assert_eq!(option(&mut a).unwrap(), Some(secondary()));
        let (rest, first) = super::super::playing::parse(a).unwrap();
        assert_eq!(rest, ["drum-modal", "192"]); assert_eq!(first, primary());
        let mut a = args("--second-stick-position-m 0.01 0.02");
        assert_eq!(option(&mut a).unwrap().unwrap().speed_m_s, 0.0);
        for s in ["--second-stick-speed-m-s 2", "--second-stick-position-m 0",
            "--second-stick-position-m NaN 0", "--second-stick-position-m 0 0 --second-stick-speed-m-s -1",
            "--second-stick-position-m 0 0 --second-stick-speed-m-s inf",
            "--second-stick-position-m 0 0 --second-stick-position-m 0 0"] {
            assert!(option(&mut args(s)).is_err(), "{s}");
        }
        assert!(admit_command(true,"splash-wav").is_err());
        assert!(admit_command(true,"snare-off-mic").is_ok());
    }

    #[test]
    fn two_actual_impacts_share_the_head_and_close_the_full_mechanical_energy() {
        for prepared in [false,true] {
            let mut one = drum_with_playing(192,2e-6,false,prepared,None,false,primary()).unwrap();
            let mut two = drum_with_sticks(192,2e-6,false,prepared,None,false,primary(),false,None,None,Some(secondary())).unwrap();
            let port = two.second_stick.unwrap();
            assert_eq!(port.coordinate,one.force.len());
            assert_eq!(&two.system.state()[..2*port.coordinate],one.system.state());
            assert_eq!(&two.observer_a[..port.coordinate],one.observer_a);
            assert_eq!(two.observer_a[port.coordinate],0.0);
            assert_eq!(two.pressure.as_ref().unwrap().areas[port.coordinate],0.0);
            let initial = 0.5*(primary().speed_m_s/two.stick_weight).powi(2)
                +0.5*(secondary().speed_m_s/port.weight).powi(2);
            let gate = CancelGate::new_clock_free();
            let mut loss = 0.0; let mut changed = 0.0_f64;
            for _ in 0..192 {
                one.system.step(&one.force,&gate).unwrap();
                let f = two.system.step(&two.force,&gate).unwrap();
                loss += f.dissipated_energy_j;
                assert!(f.balance_residual_j.abs()<1e-7);
                assert!((f.stored_energy_j+loss-initial).abs()<1e-6);
                for (&a,&b) in one.system.state()[2..].iter().zip(&two.system.state()[2..2*port.coordinate]) {
                    changed=changed.max((a-b).abs());
                }
            }
            assert!(changed>1e-12,"the second contact must change the shared head, not just add a free body");
            assert!((two.system.state()[2*port.coordinate+1]*port.weight-secondary().speed_m_s).abs()>1e-8,
                "the second stick must receive the head's contact reaction");
        }
    }

    #[test]
    fn second_striker_does_not_shift_head_radiation_or_replace_a_snare_coordinate() {
        let dt = super::super::acoustics::MECHANICAL_DT;
        let wires = super::super::snare::SnareSet::reference(false);
        let e = drum_with_sticks(4,dt,true,true,Some(wires),false,primary(),true,None,None,Some(secondary())).unwrap();
        let port = e.second_stick.unwrap();
        let Mechanics::Prepared(system) = &e.system else { panic!("prepared multi-contact host"); };
        assert_eq!(system.contact_count(),2+wires.strands*wires.contact_cells);
        let structural = e.air.as_ref().unwrap().coupling.structural_modes();
        assert_eq!(structural,port.coordinate+1+wires.mode_count().unwrap());
        assert!(e.system.state()[2*(port.coordinate+1)..].iter().all(|v|*v==0.0));
        assert!(e.acoustics.is_some(),"the same two-head acoustic boundary remains available");
        assert_eq!(e.observer_a[port.coordinate],0.0); assert_eq!(e.observer_b[port.coordinate],0.0);
        assert_eq!(e.pressure.as_ref().unwrap().areas[port.coordinate],0.0);
    }

    #[test]
    fn independent_player_files_drive_shared_contacts_and_retry_together() {
        use super::super::mechanics::drive::{Input,Program,StickDrive};
        let make = || drum_with_sticks(256,2e-6,false,true,None,false,primary(),
            false,None,None,Some(secondary())).unwrap();
        let mut driven = make(); let mut manual = make();
        let port = driven.second_stick.unwrap(); let n = driven.force.len();
        let first = "0,0\n0.000128,0\n0.000256,2\n0.000512,0";
        let second = "0,0\n0.000096,0\n0.000192,-1\n0.000384,0";
        let inputs = vec![
            Input {program:Program::parse(first).unwrap(),coordinate:0,tip_weight:driven.stick_weight},
            Input {program:Program::parse(second).unwrap(),coordinate:port.coordinate,tip_weight:port.weight},
        ];
        let mut first_stage = StickDrive::new(Program::parse(first).unwrap(),2e-6,256,
            manual.stick_weight,n).unwrap();
        let mut second_stage = StickDrive::new_inputs(vec![Input {
            program:Program::parse(second).unwrap(),coordinate:port.coordinate,tip_weight:port.weight,
        }],2e-6,256,n).unwrap();
        driven.system = driven.system.with_stick_drives(inputs,2e-6,256,n).unwrap();
        let gate = CancelGate::new_clock_free(); let mut work = 0.0_f64;
        for tick in 0..256 {
            if tick==80 {
                let before = driven.system.state().to_vec();
                let mut invalid = vec![0.0;n];invalid[port.coordinate]=1e7;
                assert!(driven.system.step(&invalid,&gate).is_err());
                assert_eq!(driven.system.state(),before);
            }
            let f = driven.system.step(&driven.force,&gate).unwrap();
            let first_load = first_stage.forces(&manual.force).unwrap();
            let both = second_stage.forces(first_load).unwrap();
            manual.system.step(both,&gate).unwrap();
            first_stage.accept();second_stage.accept();
            assert_eq!(driven.system.state(),manual.system.state());
            assert!(f.balance_residual_j.abs()<1e-7);
            work += f.supplied_work_j.abs();
        }
        assert!(work>0.0);
        assert!((driven.system.state()[2*port.coordinate+1]*port.weight-secondary().speed_m_s).abs()>1e-8);
    }

    #[test]
    fn nonlinear_two_stick_heads_keep_distinct_distributed_air_and_neck_coordinates() {
        let neck = super::super::cavity::NeckOptions {radius_m:0.005,effective_length_m:0.012,
            resistance_pa_s_m3:1000.0,azimuth_rad:0.4,axial_position_m:0.08};
        let mut compact = drum_with_sticks(128,2e-6,false,false,None,true,primary(),
            false,None,None,Some(secondary())).unwrap();
        let solid = compact.force.len();
        let mut distributed = drum_with_sticks(128,2e-6,false,false,None,true,primary(),
            true,Some(neck),None,Some(secondary())).unwrap();
        let port = distributed.second_stick.unwrap();
        let air = distributed.air.as_ref().unwrap();
        assert_eq!(air.coupling.structural_modes(),solid);
        assert_eq!(port.coordinate+1,solid);
        assert_eq!(&distributed.system.state()[..2*solid],compact.system.state());
        assert!(distributed.observer_a[solid..].iter().all(|v|*v==0.0));
        assert!(distributed.force[solid..].iter().all(|v|*v==0.0));
        // Merely moving the external stick must not masquerade as cavity gas
        // motion. Pressure changes only after its reaction moves a real head.
        let mut displaced = distributed.system.state().to_vec();
        displaced[2*port.coordinate] += 0.01/port.weight;
        assert_eq!(air.uniform_pressure(&displaced).unwrap(),0.0);
        assert_eq!(air.points(&displaced).unwrap(),(0.0,0.0));
        compact.system = compact.system.into_prepared_nonlinear().unwrap();
        distributed.system = distributed.system.into_prepared_nonlinear().unwrap();
        let gate = CancelGate::new_clock_free();let mut changed = 0.0_f64;
        let mut nonuniform = 0.0_f64;
        for _ in 0..128 {
            compact.system.step(&compact.force,&gate).unwrap();
            let f = distributed.system.step(&distributed.force,&gate).unwrap();
            assert!(f.balance_residual_j.abs()<1e-7);
            for (&a,&b) in compact.system.state().iter().zip(distributed.system.state()) {
                changed=changed.max((a-b).abs());
            }
            let air = distributed.air.as_ref().unwrap();let x = distributed.system.state();
            let (a,b) = air.points(x).unwrap();nonuniform=nonuniform.max((a-b).abs());
            let slug = air.coupling.neck_observation(x,0).unwrap();
            let volume = distributed.pressure.as_ref().unwrap();
            let expected = super::super::cavity_pressure(volume,x)
                -volume.bulk_modulus_pa/volume.volume_m3*slug.displaced_volume_m3;
            assert!((air.uniform_pressure(x).unwrap()-expected).abs()<1e-7*(1.0+expected.abs()));
        }
        assert!(changed>1e-12 && nonuniform>1e-5);
        assert!(distributed.system.membrane_observation(1).unwrap().stretching_energy_j>0.0);
    }

    #[test]
    fn missing_or_outside_second_station_refuses_instead_of_snapping_to_the_head() {
        for p in [None,Some([1.0,0.0]),Some([f64::NAN,0.0])] {
            let second = Stroke { speed_m_s:1.0,position_m:p };
            assert!(drum_with_sticks(1,2e-6,false,true,None,false,primary(),false,None,None,Some(second)).is_err());
        }
    }
}
