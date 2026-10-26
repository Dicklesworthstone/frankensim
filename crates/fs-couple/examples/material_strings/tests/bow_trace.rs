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
