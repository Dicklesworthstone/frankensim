//! Growing clamped road: actual births -> heat -> material history -> geometry.
mod driver;
mod fixture;

fn main() -> Result<(),Box<dyn std::error::Error>> {
    let law = fixture::material();
    let mut road = fixture::build(&law)?;
    println!("time_s,cells,length_m,mass_kg,tip_z_m,min_temperature_k,max_temperature_k,birth_equilibrium_pending");
    for target in [0.25,0.5,0.75,1.0,1.5,2.0,3.0] {
        let (next,report) = road.advance_to(target,fixture::stepping(0.03125),|| Ok(()))?;
        let tip = next.beam().deformed_position_m(next.beam().length_m(),0.0)?;
        let temperatures: Vec<_> = next.thermal().cells().keys()
            .map(|&id| next.thermal().temperature_k(id)).collect::<Result<_,_>>()?;
        let min = temperatures.iter().copied().fold(f64::INFINITY,f64::min);
        let max = temperatures.iter().copied().fold(f64::NEG_INFINITY,f64::max);
        assert!((next.beam().mass_kg()-next.thermal().mass_kg()).abs() < 1e-14);
        assert!(report.thermal.residual_j.abs() < 1e-9);
        if let Some(mechanics) = &report.last_mechanics {
            assert!((mechanics.vertical_reactions_n[0]-next.beam().mass_kg()*9.81).abs() < 1e-10);
        }
        println!("{target:.4},{},{:.8},{:.10e},{:.10e},{min:.6},{max:.6},{}",
            next.beam().cells().len(),next.beam().length_m(),next.beam().mass_kg(),tip[1],report.birth_equilibrium_pending);
        road = next;
    }
    assert_eq!(road.beam().cells().len(),8);
    assert!(road.beam().deformed_position_m(road.beam().length_m(),0.0)?[1] < 0.0);
    Ok(())
}
