//! A single SI player force acting through a geometry-owned displacement row.
use super::{Error, Input, Program, StickDrive};

pub struct SpatialInput {
    pub program: Program,
    /// Physical point displacement / generalized displacement; signed, not normalized.
    pub weights: Vec<f64>,
}
impl StickDrive {
    /// Admit a complete performance without requiring an artificial scalar
    /// port. Every spatial hand and scalar stick/pedal/jaw shares one clock.
    /// Construction is all-or-nothing; no program or physical tick is consumed.
    pub fn new_mixed(inputs: Vec<Input>, spatial: Vec<SpatialInput>, dt_s: f64,
        steps: u64, modes: usize) -> Result<Self, Error>
    {
        if modes == 0 || inputs.len() > 4 || spatial.len() > 4
            || !(1..=4).contains(&(inputs.len() + spatial.len())) {
            return Err("player drive needs one to four physical inputs and a nonempty force basis".into());
        }
        // Reuse scalar admission whenever scalar ports are present. With
        // only spatial hands, add_spatial_input admits the same clock and
        // force bounds using each row's largest absolute component.
        let mut drive = if inputs.is_empty() {
            Self { inputs, spatial_inputs: Vec::new(), dt_s, steps,
                accepted: 0, force: vec![0.0; modes] }
        } else {
            Self::new_inputs(inputs, dt_s, steps, modes)?
        };
        for input in spatial { drive.add_spatial_input(input)?; }
        Ok(drive)
    }

    /// Cold attachment before the first accepted tick. Rows may overlap because
    /// two physical hands can act on the same body; their forces are summed.
    pub fn add_spatial_input(&mut self, input: SpatialInput) -> Result<(), Error> {
        if self.accepted != 0 || self.inputs.len() + self.spatial_inputs.len() >= 4
            || input.weights.len() != self.force.len()
            || input.weights.iter().any(|x| !x.is_finite()) {
            return Err("spatial player needs an unstarted clock, at most four inputs and a finite full-basis row".into());
        }
        let maximum = input.weights.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
        // Admit the largest component but RETAIN every original signed weight.
        // The mechanical owner separately admits the final sum of all forces.
        input.program.admit(self.dt_s, self.steps, maximum)?;
        self.spatial_inputs.push(input);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn drive() -> StickDrive {
        StickDrive::new_inputs(vec![Input { program: Program::parse("0,0\n1,0").unwrap(),
            coordinate: 1, tip_weight: 1.0 }],0.5,2,4).unwrap()
    }
    fn hand(weights: Vec<f64>) -> SpatialInput {
        SpatialInput { program: Program::parse("0,0\n0.25,4\n0.5,0\n0.75,-4\n1,0").unwrap(), weights }
    }
    #[test]
    fn signed_hand_row_preserves_work_and_one_retryable_clock() {
        let mut d=drive();let row=vec![2.0,0.0,-3.0,0.5];
        d.add_spatial_input(hand(row.clone())).unwrap();
        let external=[1.0,7.0,-2.0,3.0];let velocity=[0.2,-0.4,0.1,0.8];
        let staged=d.forces(&external).unwrap().to_vec();
        assert_eq!(staged,[5.0,7.0,-8.0,4.0]);
        let power=staged.iter().zip(external).zip(velocity).map(|((f,e),v)|(f-e)*v).sum::<f64>();
        let hand_speed=row.iter().zip(velocity).map(|(b,v)|b*v).sum::<f64>();
        assert!((power-2.0*hand_speed).abs()<1e-14);
        assert!(d.forces(&[0.0]).is_err());assert_eq!(d.accepted,0);
        assert_eq!(d.forces(&external).unwrap(),staged);
        d.accept();assert_eq!(d.forces(&external).unwrap(),[-3.0,7.0,4.0,2.0]);
    }
    #[test]
    fn spatial_only_hands_need_no_dummy_port_and_sum_reciprocal_work() {
        let a=vec![2.0,0.0,-3.0,0.5];let b=vec![-1.0,4.0,1.0,0.0];
        let mut d=StickDrive::new_mixed(Vec::new(),vec![hand(a.clone()),hand(b.clone())],0.5,2,4).unwrap();
        assert!(d.inputs.is_empty());assert_eq!(d.spatial_inputs.len(),2);
        let external=[1.0,7.0,-2.0,3.0];let velocity=[0.2,-0.4,0.1,0.8];
        let forces=d.forces(&external).unwrap().to_vec();
        assert_eq!(forces,[3.0,15.0,-6.0,4.0]);
        let power=forces.iter().zip(external).zip(velocity).map(|((f,e),v)|(f-e)*v).sum::<f64>();
        let physical=2.0*a.iter().zip(velocity).map(|(w,v)|w*v).sum::<f64>()
            +2.0*b.iter().zip(velocity).map(|(w,v)|w*v).sum::<f64>();
        assert!((power-physical).abs()<1e-14);
        assert!(d.forces(&[f64::NAN;4]).is_err());assert_eq!(d.accepted,0);
        assert_eq!(d.forces(&external).unwrap(),forces);
        d.accept();assert_eq!(d.forces(&external).unwrap(),[-1.0,-1.0,2.0,2.0]);
        d.accept();assert!(d.forces(&external).is_err());
    }
    #[test]
    fn spatial_only_construction_admits_clock_duration_and_combined_port_budget() {
        for (dt,steps) in [(0.0,2),(f64::NAN,2),(0.5,0),(0.5,1),(0.5,(1_u64<<53)+1)] {
            assert!(StickDrive::new_mixed(Vec::new(),vec![hand(vec![1.0;4])],dt,steps,4).is_err());
        }
        assert!(StickDrive::new_mixed(Vec::new(),Vec::new(),0.5,2,4).is_err());
        for row in [vec![],vec![0.0;4],vec![f64::NAN;4],vec![f64::MAX;4]] {
            assert!(StickDrive::new_mixed(Vec::new(),vec![hand(row)],0.5,2,4).is_err());
        }
        assert!(StickDrive::new_mixed(Vec::new(),(0..4).map(|_|hand(vec![1.0;4])).collect(),0.5,2,4).is_ok());
        let scalar=Input{program:Program::parse("0,0\n1,0").unwrap(),coordinate:0,tip_weight:1.0};
        assert!(StickDrive::new_mixed(vec![scalar],(0..4).map(|_|hand(vec![1.0;4])).collect(),0.5,2,4).is_err());
    }
    #[test]
    fn shape_clock_and_force_refusals_do_not_attach_partial_inputs() {
        let mut d=drive();
        for row in [vec![],vec![0.0;4],vec![f64::NAN;4],vec![f64::MAX;4]] {
            assert!(d.add_spatial_input(hand(row)).is_err());assert!(d.spatial_inputs.is_empty());
        }
        let late=SpatialInput {program:Program::parse("0,0\n2,0").unwrap(),weights:vec![1.0;4]};
        assert!(d.add_spatial_input(late).is_err());assert!(d.spatial_inputs.is_empty());
        d.add_spatial_input(hand(vec![1.0,0.0,0.0,0.0])).unwrap();
        d.accept();assert!(d.add_spatial_input(hand(vec![1.0;4])).is_err());
        assert_eq!(d.spatial_inputs.len(),1);
    }
}
