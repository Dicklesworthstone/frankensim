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
    fn missing_or_outside_second_station_refuses_instead_of_snapping_to_the_head() {
        for p in [None,Some([1.0,0.0]),Some([f64::NAN,0.0])] {
            let second = Stroke { speed_m_s:1.0,position_m:p };
            assert!(drum_with_sticks(1,2e-6,false,true,None,false,primary(),false,None,None,Some(second)).is_err());
        }
    }
}
