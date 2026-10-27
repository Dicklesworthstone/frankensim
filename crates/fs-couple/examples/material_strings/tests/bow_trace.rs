//! `--bow-trace`: the same WAV bytes, with each accepted sample's bow port beside them.
use super::*;

fn bowed()->input::Spec {
    // The bow lands at 10 ms, reverses at 20 ms and lifts at 40 ms; 2048 samples cover all three.
    let text=include_str!("../pluck-and-bow.material-strings").replace("48000 4096 0.01 128","48000 2048 0.01 128");
    input::Spec::parse(text.as_bytes()).unwrap()
}

#[test]
fn a_traced_render_writes_the_same_wav_and_one_row_per_sample_and_bow_port() {
    scoped(|cx,gate| {
        let mut plain=compile(cx,gate,bowed(),7800.0);
        let (wav,summary)=plain.write(Cursor::new(Vec::<u8>::new()),gate).unwrap();
        let mut traced=compile(cx,gate,bowed(),7800.0);
        let mut trace=Vec::<u8>::new();
        let (traced_wav,traced_summary)=traced.write_traced(Cursor::new(Vec::<u8>::new()),&mut trace,gate).unwrap();
        assert_eq!(traced_wav.into_inner(),wav.into_inner());
        assert_eq!(traced_summary,summary);
        let text=String::from_utf8(trace).unwrap();
        let mut lines=text.lines();
        assert_eq!(lines.next(),Some("sample\tport\tregime\tslip_velocity_m_s\ttraction_n\tnormal_force_n"));
        let rows:Vec<Vec<&str>>=lines.map(|line|line.split('\t').collect()).collect();
        assert_eq!(rows.len(),2048);
        let mut regimes=std::collections::BTreeSet::new();
        for (k,row) in rows.iter().enumerate() {
            assert_eq!(row.len(),6);
            assert_eq!(row[0].parse::<u64>().unwrap(),k as u64+1);
            assert_eq!(row[1],"0");
            regimes.insert(row[2]);
            for value in &row[3..] { assert!(value.parse::<f64>().unwrap().is_finite()); }
            // A released port carries no load and no traction.
            if row[2]=="Released" { assert_eq!(row[5].parse::<f64>().unwrap(),0.0); }
        }
        // Before the bow lands, while it sticks and slips, and after it lifts.
        assert!(["Released","Sticking","Sliding"].iter().all(|r|regimes.contains(r)),"regimes {regimes:?}");
        assert_eq!(rows[100][2],"Released");
        assert_eq!(rows[2047][2],"Released");
    });
}

#[test]
fn decimated_pressure_retains_every_mechanical_trace_and_the_same_final_state() {
    scoped(|cx,gate| {
        let source=include_str!("../pluck-and-bow.material-strings")
            .replace("48000 4096 0.01 128","96000 4096 0.01 128");
        let build=|block|compile(cx,gate,input::Spec::parse(source
            .replace("96000 4096 0.01 128",&format!("96000 4096 0.01 {block}")).as_bytes()).unwrap(),7800.0);
        let mut raw=build(128);
        let mut decimated=build(37);
        let mut plain=build(128);
        let mut raw_trace=Vec::new();
        let mut low_trace=Vec::new();
        raw.write_traced(Cursor::new(Vec::new()),&mut raw_trace,gate).unwrap();
        let (observed,summary)=decimated.write_observed(Cursor::new(Vec::new()),Some(&mut low_trace),true,gate).unwrap();
        let (without_trace,other)=plain.write_observed(Cursor::new(Vec::new()),None,true,gate).unwrap();
        assert_eq!(summary.sample_rate_hz,48_000);
        assert_eq!(summary.samples,2048);
        assert_eq!(summary,other);
        assert_eq!(observed.into_inner(),without_trace.into_inner());
        assert_eq!(raw_trace,low_trace,"output sampling cannot discard mechanical trace rows");
        assert_eq!(String::from_utf8(low_trace).unwrap().lines().count(),4097);
        assert_eq!(raw.performance.bank().system().last_frame(),decimated.performance.bank().system().last_frame());
        assert_eq!(raw.performance.bank().system().last_frame(),plain.performance.bank().system().last_frame());
        assert_eq!(decimated.performance.samples_rendered(),4096);
    });
}

#[test]
fn a_trace_io_failure_does_not_publish_a_decimated_group_or_restart_missing_history() {
    struct FailedTrace { remaining:usize }
    impl Write for FailedTrace {
        fn write(&mut self,bytes:&[u8])->std::io::Result<usize> {
            if self.remaining==0 { return Err(std::io::Error::other("injected trace I/O failure")); }
            let count=bytes.len().min(self.remaining);self.remaining-=count;Ok(count)
        }
        fn flush(&mut self)->std::io::Result<()> { Ok(()) }
    }
    scoped(|cx,gate| {
        let text=include_str!("../pluck-and-bow.material-strings")
            .replace("48000 4096 0.01 128","96000 4096 0.01 128");
        let mut ready=compile(cx,gate,input::Spec::parse(text.as_bytes()).unwrap(),7800.0);
        let header="sample\tport\tregime\tslip_velocity_m_s\ttraction_n\tnormal_force_n\n";
        let mut trace=FailedTrace { remaining:header.len() };
        let mut output=Cursor::new(Vec::new());
        let error=ready.write_observed(&mut output,Some(&mut trace),true,gate).unwrap_err();
        assert!(error.to_string().contains("injected trace I/O failure"));
        assert_eq!(ready.performance.samples_rendered(),1);
        assert_eq!(output.into_inner().len(),44,"failed group has no emitted pressure payload");
        assert!(ready.write_observed(Cursor::new(Vec::new()),None,true,gate).is_err());
        assert_eq!(ready.performance.samples_rendered(),1,"a new filter must not discard accepted history");
    });
}
