//! CLI composition; both cymbals use one time owner and one pressure scene.
use super::*;
use std::path::Path;
use acoustics::stereo::radiation_spec;

pub fn is_command(command:Option<&str>)->bool {matches!(command,Some("hihat"|"hihat-wav"|"hihat-mic"))}

pub fn run(mut args:Vec<String>)->Result<(),Error> {
    let analytic=mechanics::analytic_option(&mut args)?;
    let substeps=mechanics::substeps_option(&mut args)?;
    let prepared=mechanics::prepared_option(&mut args)?||analytic||substeps.is_some();
    let right=acoustics::stereo::option(&mut args)?;
    let radiation=radiation_spec::option(&mut args)?;
    let feedback=acoustics::stereo::feedback::option(&mut args)?;
    let second=sticks::option(&mut args)?;
    let first_force=mechanics::drive::option(&mut args)?;
    let second_force=mechanics::drive::second_option(&mut args)?;
    let (args,stroke)=playing::parse(args)?;
    let command=args.first().map(String::as_str);
    if !is_command(command)||args.len()<2 {return Err("usage: hihat INPUT.fshh [steps]; hihat-wav INPUT [frames] [full_scale_pa]; hihat-mic INPUT [frames] [full_scale_pa] [x y z]; see HIHAT.md".into());}
    acoustics::stereo::feedback::admit_command(feedback,command.unwrap(),false)?;
    let mic=command==Some("hihat-mic");let audio=command!=Some("hihat");
    if (!audio&&args.len()>3)||(audio&&args.len()>4&&!(mic&&args.len()==7))
        ||(!mic&&right.is_some())||(!audio&&radiation.is_some())||second_force.is_some()&&second.is_none() {
        return Err("hi-hat received incompatible output, receiver or second-stick controls".into());
    }
    let count=if args.len()>=3 {args[2].parse::<u64>()?}else{if audio{48000}else{4096}};
    if count==0||count>if audio{480000}else{1_000_000}{return Err("hi-hat duration exceeds the existing output budget".into());}
    let scale=if args.len()>=4 {args[3].parse::<f64>()?}else{20.};
    if !scale.is_finite()||scale<=0.{return Err("hi-hat PCM scale must be positive finite Pa".into());}
    let dt=if audio{acoustics::MECHANICAL_DT}else{2e-6};
    let steps=if audio{count.checked_mul(acoustics::SUBSTEPS as u64).ok_or("hi-hat clock overflow")?}else{count};
    let receiver=if mic {acoustics::Receiver::FinitePoint(if args.len()==7{
        [args[4].parse()?,args[5].parse()?,args[6].parse()?]
    }else{[0.08,0.05,0.35]})}else{acoustics::Receiver::FarField([1.5,0.7,1.5])};
    let spec=Spec::load(Path::new(&args[1]))?;let [upper,lower]=spec.shells()?;
    let pair=build(&spec,&upper,&lower,stroke,second,steps,dt,audio)?;
    let mut receivers=vec![receiver];if let Some(p)=right{receivers.push(acoustics::Receiver::FinitePoint(p));}
    let (mut e,loaded)=if feedback {
        let (e,bake)=acoustics::stereo::feedback::prepare(pair.experiment,usize::try_from(count)?,scale,
            &receivers,radiation.unwrap_or_default(),&CancelGate::new_clock_free())?;(e,Some(bake))
    }else{(pair.experiment,None)};
    if prepared {e.system=if analytic{e.system.into_analytic_nonlinear()?}else{e.system.into_prepared_nonlinear()?};}
    if let Some(b)=substeps{e.system=e.system.with_impact_substeps(b)?;}
    let mut inputs=vec![mechanics::drive::Input{program:mechanics::drive::Program::parse(&spec.pedal)?,
        coordinate:pair.pedal.coordinate,tip_weight:pair.pedal.weight}];
    if let Some(program)=first_force{inputs.push(mechanics::drive::Input{program,coordinate:0,tip_weight:e.stick_weight});}
    if let Some(program)=second_force{let p=e.second_stick.ok_or("missing hi-hat second stick")?;
        inputs.push(mechanics::drive::Input{program,coordinate:p.coordinate,tip_weight:p.weight});}
    e.system=e.system.with_stick_drives(inputs,dt,steps,e.force.len())?;
    eprintln!("paired cymbals: upper_modes={}, lower_modes={}, contact_sites={}, one joint mechanics; supplied masses/geometry and authored contact, not a calibrated hi-hat; axial carriage, no rocking or squeeze-film air",
        pair.upper_modes.len(),pair.lower_modes.len(),pair.collision.n_points());
    let gate=CancelGate::new_clock_free();let stdout=std::io::stdout();let mut out=std::io::BufWriter::new(stdout.lock());
    if audio {
        let wav=match loaded {
            Some(bake)=>bake.render(&mut e,usize::try_from(count)?,scale,&gate)?,
            None=>acoustics::stereo::render_receivers_with_spec(&mut e,usize::try_from(count)?,scale,&receivers,
                radiation.unwrap_or_default(),&gate)?,
        };
        out.write_all(&wav)?;out.flush()?;return Ok(());
    }
    writeln!(out,"time_s,upper_down_m,lower_down_m,pedal_down_m,pedal_speed_m_s,min_contact_gap_m,active_sites,total_energy_j,loss_j,player_work_j,balance_j")?;
    for _ in 0..steps {
        let f=e.system.step(&e.force,&gate)?;let x=e.system.state();
        let displacement=|row:&[f64]|row.iter().enumerate().map(|(i,b)|b*x[2*i]).sum::<f64>();
        let mut gap=f64::INFINITY;let mut active=0;
        for (row,&clearance) in pair.collision.collocation().chunks(e.force.len()).zip(pair.collision.gaps()) {
            let g=clearance-displacement(row);gap=gap.min(g);active+=usize::from(g<0.);
        }
        writeln!(out,"{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{},{:.17e},{:.17e},{:.17e},{:.17e}",
            f.time_s,displacement(&e.observer_a),displacement(&e.observer_b),x[2*pair.pedal.coordinate]*pair.pedal.weight,
            x[2*pair.pedal.coordinate+1]*pair.pedal.weight,gap,active,f.stored_energy_j,f.dissipated_energy_j,
            f.supplied_work_j,f.balance_residual_j)?;
    }
    out.flush()?;Ok(())
}
