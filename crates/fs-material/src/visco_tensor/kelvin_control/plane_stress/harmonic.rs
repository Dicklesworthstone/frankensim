//! Coherent drive of every coupled original Kelvin/Maxwell state coordinate.
use super::*;

pub struct TensorKelvinPlaneStressHarmonicWindow<'p> {
    control: &'p PreparedTensorKelvinPlaneStress,
    source: TensorPronyHarmonicSpectrum<'p>,
    temporal: HarmonicEnergySpectrum,
    kernels: Vec<PreparedHarmonicDecayInterval>,
    phases: Vec<HarmonicPhaseInterval>,
    start: Dd,
    duration: f64,
    payload_bound: usize,
}
impl PreparedTensorKelvinPlaneStress {
    /// Prepare positive physical-time moments once for repeated points. Peak
    /// phasors use exp(-i omega t), with a separately retained DD start clock.
    pub fn prepare_harmonic_window(&self, frequencies:&[f64],start:Dd,duration:f64,
        budget:TensorPlaneStressHarmonicBudget,mut cancelled:impl FnMut()->bool)
        ->Result<TensorKelvinPlaneStressHarmonicWindow<'_>>
    {
        poll(&mut cancelled)?;
        let f=frequencies.len();let n=self.state_coordinates();
        let kernel_bytes=PreparedHarmonicDecayInterval::required_workspace_bytes(f).map_err(moment_error)?;
        let kernel_work=PreparedHarmonicDecayInterval::required_work(f).map_err(moment_error)?;
        let bytes=n.checked_mul(kernel_bytes+256).and_then(|v|v.checked_add(
            core::mem::size_of::<TensorKelvinPlaneStressHarmonicWindow<'_>>() + 160*1024
            +128*f*(self.model.branches().len()+1)+512*f+64*f*f))
            .ok_or_else(||invalid("Kelvin harmonic payload extent overflow"))?;
        let work=n.checked_mul(kernel_work).and_then(|v|v.checked_add(128*f*f+2048*f*(n+1)+32*n*n+8192))
            .ok_or_else(||invalid("Kelvin harmonic work extent overflow"))?;
        if bytes>budget.maximum_payload_bytes || work>budget.maximum_scalar_work {
            return Err(invalid("Kelvin harmonic interval exceeds its explicit payload or work budget"));
        }
        let source=self.model.prepare_harmonic_spectrum(frequencies)
            .map_err(|_|invalid("Kelvin harmonic carriers leave the original material domain"))?;
        let start=number(Dd::from_pair(start.hi,start.lo))?;
        let temporal=HarmonicEnergySpectrum::new_at_clock(source.distinct_frequencies(),start,duration)
            .map_err(|_|invalid("invalid physical Kelvin harmonic interval"))?;
        let limits=HarmonicDecayIntervalLimits { maximum_tones:LINES,
            maximum_workspace_bytes:kernel_bytes,maximum_work:kernel_work };
        let mut kernels=Vec::with_capacity(n);
        for &rate in &self.rates {
            poll(&mut cancelled)?;
            kernels.push(PreparedHarmonicDecayInterval::new_cancellable(source.distinct_frequencies(),
                start,duration,1.0/rate,limits,&mut cancelled).map_err(moment_error)?);
        }
        let phases=source.distinct_frequencies().iter().map(|&omega|
            phase_angular_interval_from_clock_dd(omega,start,duration)
                .map_err(|_|invalid("unrepresentable original Kelvin harmonic phase")))
            .collect::<Result<Vec<_>>>()?;
        poll(&mut cancelled)?;
        Ok(TensorKelvinPlaneStressHarmonicWindow {control:self,source,temporal,kernels,phases,
            start,duration,payload_bound:bytes})
    }
}

impl TensorKelvinPlaneStressHarmonicWindow<'_> {
    pub const fn start_time_s_dd(&self)->Dd {self.start}
    pub const fn duration_s(&self)->f64 {self.duration}
    pub fn angular_frequencies_rad_s(&self)->&[f64] {self.source.angular_frequencies_rad_s()}
    /// Conservative retained storage plus one point's scratch and candidate;
    /// separately owned original dense constitutive preparation is excluded.
    pub const fn numerical_payload_bound_bytes(&self)->usize {self.payload_bound}

    /// Drive the actual original state with engineering [epsilon_aa,epsilon_bb,
    /// gamma_ab]. Input strain must join the accepted history. Zero transverse
    /// traction determines its rate without a hidden strain or memory jump.
    /// Singular transverse viscosity additionally requires the original initial
    /// spring traction to satisfy its algebraic null-space constraints.
    /// Equal carriers combine before source admission and recovery. The complete
    /// imposed waveform and a conservative full transient strain envelope retain
    /// the original source domain. History changes only on candidate adoption.
    pub fn evaluate(&self,strains:&[[C64;3]],bias:&[f64;3],history:&TensorPronyHistory,
        mut cancelled:impl FnMut()->bool)->Result<TensorKelvinPlaneStressHarmonicPoint>
    {
        poll(&mut cancelled)?;
        let control=self.control;let model=&control.model;
        if !Arc::ptr_eq(model,&history.model) || bias.iter().any(|v|!v.is_finite()) {
            return Err(invalid("Kelvin plane stress requires its original history and finite bias"));
        }
        if let Some(constraint)=&control.constraint {constraint.admit_history(control,history)?;}
        let grouped=self.source.group_strains(strains).map_err(|_|invalid("invalid coherent Kelvin plane strains"))?;
        let frequencies=self.source.distinct_frequencies();let f=frequencies.len();
        let n=control.state_coordinates();let m=control.rates.len();let tolerance=8192.0*f64::EPSILON*(n+6) as f64;
        let mut mapped=[[[Dd::ZERO;2];6];LINES];
        let mut embedded=mapped;
        let mut mapped_bias=[Dd::ZERO;6];
        for axis in 0..6 {
            let mut sum=Sum::default();
            for j in 0..3 {sum.add([control.map[j][axis],bias[j],1.0,1.0]);}
            mapped_bias[axis]=number(sum.retained_value())?;
            for k in 0..f {for part in 0..2 {
                let mut sum=Sum::default();
                for j in 0..3 {add(&mut sum,grouped[k][j][part],control.map[j][axis]);}
                mapped[k][axis][part]=number(sum.retained_value())?;
            }}
        }
        for k in 0..f {embedded[k][..3].copy_from_slice(&grouped[k]);}
        self.source.admit_grouped_strains_with_bias(&mapped[..f],&mapped_bias)
            .map_err(|_|invalid("full imposed Kelvin waveform leaves the original material domain"))?;
        let mut target=history.strain;let mut delta=[Dd::ZERO;6];
        let mut du=[Dd::ZERO;3];let mut end_rate=[Dd::ZERO;6];
        for j in 0..3 {
            let mut start=Dd::from_f64(bias[j]);let mut end=start;let mut rate=Sum::default();
            let mut magnitude=bias[j].abs()+history.strain[control.inside[j]].abs()/ENGINEERING[j];
            for k in 0..f {
                let [a,b]=grouped[k][j];let phase=self.phases[k];
                start=start+a*phase.start[0]-b*phase.start[1];
                end=end+a*phase.end[0]-b*phase.end[1];
                du[j]=du[j]+a*phase.increment[0]-b*phase.increment[1];
                add(&mut rate,a*phase.end[1]+b*phase.end[0],frequencies[k]);
                magnitude+=l1(grouped[k][j]);
            }
            let axis=control.inside[j];
            if (start*Dd::from_f64(ENGINEERING[j])-Dd::from_f64(history.strain[axis])).abs().to_f64()
                >128.0*f64::EPSILON*ENGINEERING[j]*magnitude {
                return Err(invalid("in-plane Kelvin harmonic strain is discontinuous at the accepted clock"));
            }
            delta[axis]=number(du[j]*Dd::from_f64(ENGINEERING[j]))?;
            target[axis]=number(end*Dd::from_f64(ENGINEERING[j]))?.to_f64();
            end_rate[axis]=number(rate.retained_value()*Dd::from_f64(ENGINEERING[j]))?;
        }
        let mut lag=control.outside.map(|i|Dd::from_f64(history.strain[i]));
        for i in 0..3 {for j in 0..3 {
            lag[i]=lag[i]-Dd::from_f64(control.map[j][control.outside[i]])
                *(Dd::from_f64(history.strain[control.inside[j]])/Dd::from_f64(ENGINEERING[j]));
        }}
        let initial_elastic=control.equilibrium.apply(lag)?;
        let mut initial_state=vec![Dd::ZERO;n];initial_state[..3].copy_from_slice(&initial_elastic);
        for row in 3..n {initial_state[row]=Dd::from_f64(history.q[(row-3)/6][(row-3)%6]);}
        let identity=[[(1.0,0),(0.0,0)];LINES];
        let mut work=Sum::default();
        self.temporal.add_stress_work_with_bias(&model.c_eq,&mapped[..f],&identity[..f],&mapped_bias,&mut work);
        let direct=self.temporal.kelvin_heat(&control.direct_gram,&embedded[..f])
            .map_err(|_|invalid("unresolved original direct Kelvin heat"))?;
        let mut rate_spectrum=[[(0.0,0),(0.0,0)];LINES];
        for k in 0..f {rate_spectrum[k][1]=(frequencies[k],0);}
        self.temporal.add_stress_work(&control.direct_work,&embedded[..f],&rate_spectrum[..f],&mut work);
        let mut heat=Sum::default();heat.add([direct,1.0,1.0,1.0]);
        let mut absolute=mapped_bias.map(|v|v.hi.abs()+v.lo.abs());
        let mut excursion=history.strain.map(f64::abs);
        for k in 0..f {for i in 0..6 {
            absolute[i]+=l1(mapped[k][i]);
            excursion[i]+=(frequencies[k]*self.duration).min(2.0)*l1(mapped[k][i]);
        }}
        let mut changes=Vec::with_capacity(m);
        for mode in 0..m {
            poll(&mut cancelled)?;
            let root=control.rate_roots[mode];let mut initial=Sum::default();
            for row in 0..n {for word in [initial_state[row].hi,initial_state[row].lo] {
                initial.add([control.modes[row*m+mode],word,root,1.0]);
            }}
            let initial=number(initial.retained_value())?;
            let drive=rate_fields(&control.forcing[mode],&grouped[..f],frequencies,root,1.0)?;
            let port=rate_fields(&control.forcing[mode],&grouped[..f],frequencies,1.0,root)?;
            let response=self.kernels[mode].evaluate_cancellable(initial,&drive[..f],&port[..f],&mut cancelled).map_err(moment_error)?;
            add(&mut heat,response.mean_square,self.duration);add(&mut work,response.mean_mixed,self.duration);
            changes.push(response.increment);
            let tau=1.0/control.rates[mode];let forcing:f64=drive[..f].iter().copied().map(l1).sum();
            let driven=self.duration.min(tau)*forcing;
            let decay=if self.duration>=tau {1.0}else{self.duration/tau};
            for i in 0..6 {
                let direction=control.strain_modes[mode][i].abs()/root;
                absolute[i]+=direction*(initial.abs().to_f64()+driven);
                excursion[i]+=direction*(decay*initial.abs().to_f64()+driven);
            }
        }
        let maximum=core::array::from_fn(|i|absolute[i].min(excursion[i]));
        model.validate_mandel_strain(&maximum).map_err(|_|invalid("full transient Kelvin strain envelope exceeds the original material domain"))?;
        let mut state_change=vec![Dd::ZERO;n];
        for row in 0..n {
            poll(&mut cancelled)?;let mut sum=Sum::default();
            for mode in 0..m {for word in [changes[mode].hi,changes[mode].lo] {
                sum.add_ratio(&[control.modes[row*m+mode],word],&[control.rate_roots[mode]]);
            }}
            state_change[row]=number(sum.retained_value())?;
        }
        let free=control.equilibrium.solve(core::array::from_fn(|i|state_change[i]),false)?;
        for i in 0..3 {
            let axis=control.outside[i];delta[axis]=free[i];
            for j in 0..3 {delta[axis]=delta[axis]+Dd::from_f64(control.map[j][axis])*du[j];}
            target[axis]=number(Dd::from_f64(history.strain[axis])+delta[axis])?.to_f64();
        }
        let mut q=history.q.clone();
        for row in 3..n {q[(row-3)/6][(row-3)%6]=number(initial_state[row]+state_change[row])?.to_f64();}
        model.validate_mandel_strain(&target).map_err(|_|invalid("Kelvin endpoint leaves the original strain domain"))?;
        // Reconstruct the actual original stress and compatible endpoint rate.
        let eq1=apply(&model.equilibrium_gram,&target);let mut spring=[Dd::ZERO;6];
        let mut scale=[0.0;6];
        for i in 0..6 {
            for row in 0..6 {
                let term=Dd::from_f64(model.equilibrium_gram[row][i])*Dd::from_f64(eq1[row]);
                spring[i]=spring[i]+term;scale[i]+=term.abs().to_f64();
            }
            for (j,branch) in model.branches().iter().enumerate() {
                let gram=branch.resolved_gram_rows(&model.factors[j]);
                for row in 0..6 {
                    let term=Dd::from_f64(gram[row][i])*Dd::from_f64(q[j][row]);
                    spring[i]=spring[i]+term;scale[i]+=term.abs().to_f64();
                }
            }
        }
        let viscosity=model.kv_pa_s().expect("admitted original Kelvin tensor");
        let rhs=control.outside.map(|axis| {
            let mut value=-spring[axis];for &j in &control.inside {value=value-Dd::from_f64(viscosity[axis][j])*end_rate[j];} value
        });
        let rate=if let Some(constraint)=&control.constraint {
            constraint.endpoint_rate(control,&spring,&end_rate,&q)?
        } else {
            let factor=control.viscosity.as_ref().expect("admitted definite transverse viscosity");
            factor.solve(factor.solve(rhs,true)?,false)?
        };
        for i in 0..3 {end_rate[control.outside[i]]=rate[i];}
        let mut stress=[0.0;6];
        for i in 0..6 {
            let mut value=spring[i];
            for j in 0..6 {let term=Dd::from_f64(viscosity[i][j])*end_rate[j];value=value+term;scale[i]+=term.abs().to_f64();}
            stress[i]=number(value)?.to_f64();
            if control.outside.contains(&i) {
                if !scale[i].is_finite() || stress[i].abs()>tolerance*scale[i] {return Err(invalid("original Kelvin plane-stress constraint is unresolved"));}
                stress[i]=0.0;
            }
        }
        let eq0=apply(&model.equilibrium_gram,&history.strain);
        let eq_delta=model.equilibrium_gram.map(|row|(0..6).fold(Dd::ZERO,|s,i|s+Dd::from_f64(row[i])*delta[i]));
        let mut change=Dd::ZERO;let mut exchange=0.0;let mut stored=0.5*dot(&eq1,&eq1);
        for i in 0..n+3 {
            let (old,difference)=if i<6 {(eq0[i],eq_delta[i])}
                else {(history.q[(i-6)/6][(i-6)%6],state_change[i-3])};
            let term=(Dd::from_f64(old)+Dd::from_f64(0.5)*difference)*difference;
            change=change+term;exchange+=term.abs().to_f64();
        }
        for memory in &q {stored+=0.5*dot(memory,memory);}
        let work=number(work.retained_value())?;let heat=number(heat.retained_value())?;
        let residual=number(work-change-heat)?.to_f64();
        let cumulative_work=(Dd::from_f64(history.work)+work).to_f64();
        let cumulative_heat=(Dd::from_f64(history.heat)+heat).to_f64();
        if q.iter().any(|v|!finite(v)) || ![stored,exchange,cumulative_work,cumulative_heat].iter().all(|v|v.is_finite())
            || stored<0.0 || heat.lt(Dd::ZERO)
            || residual.abs()>tolerance*(exchange+work.abs().to_f64()+heat.to_f64()) {
            return Err(invalid("independent Kelvin plane-stress work/storage/heat is unresolved"));
        }
        let account=TensorKelvinPlaneStressHarmonicAccount {start_stored_energy_j_m3:history.stored,
            end_stored_energy_j_m3:stored,stored_energy_change_j_m3:number(change)?.to_f64(),
            material_heat_j_m3:heat.to_f64(),stress_work_j_m3:work.to_f64(),energy_balance_residual_j_m3:residual};
        let mut next=history.clone();next.strain=target;next.q=q;next.stress=stress;
        next.stored=stored;next.work=cumulative_work;next.heat=cumulative_heat;
        poll(&mut cancelled)?;
        Ok(TensorKelvinPlaneStressHarmonicPoint {history:next,account,strain_increment:tensor(&delta.map(Dd::to_f64)),
            strain_rate:tensor(&end_rate.map(Dd::to_f64)),strain_bounds:[tensor(&maximum.map(|v|-v)),tensor(&maximum)]})
    }
}

fn rate_fields(row:&[f64;3],strains:&[PlaneParts],frequencies:&[f64],numerator:f64,denominator:f64)
    ->Result<[[Dd;2];LINES]>
{
    let mut result=[[Dd::ZERO;2];LINES];
    for (k,strain) in strains.iter().enumerate() {
        let mut sum=[Sum::default();2];
        for j in 0..3 {for [re,im] in [[strain[j][0].hi,strain[j][1].hi],[strain[j][0].lo,strain[j][1].lo]] {
            sum[0].add_ratio(&[row[j],frequencies[k],im,numerator],&[denominator]);
            sum[1].add_ratio(&[-row[j],frequencies[k],re,numerator],&[denominator]);
        }}
        result[k]=[number(sum[0].retained_value())?,number(sum[1].retained_value())?];
    }
    Ok(result)
}
