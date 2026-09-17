use super::*;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64,Ordering};
static NEXT:AtomicU64=AtomicU64::new(0);
struct Scratch(PathBuf);
impl Scratch {
    fn new()->Self {
        let path=std::env::temp_dir().join(format!("fs-enclosure-uq-{}-{}",std::process::id(),NEXT.fetch_add(1,Ordering::SeqCst)));
        std::fs::create_dir(&path).unwrap();Self(path)
    }
    fn uq(&self,plan:&str,args:&[&str])->Output {
        std::fs::write(self.0.join("uq.json"),plan).unwrap();
        Command::new(env!("CARGO_BIN_EXE_frankensim")).current_dir(&self.0)
            .args(["--json","cooling-network-uq","base.json","uq.json"]).args(args).output().unwrap()
    }
}
impl Drop for Scratch {fn drop(&mut self){let _=std::fs::remove_dir_all(&self.0);}}

#[test]
fn uncertain_enclosure_finishes_use_real_samples_and_exact_checkpoint_resume() {
    let dir=Scratch::new();std::fs::write(dir.0.join("base.json"),BASE).unwrap();
    let plan=r#"{"schema":"frankensim.cooling-network-uq.v1","seed":"73","samples":2,"wall_seconds":600,"temperature_limit_k":340,"correlation":{"kind":"independent"},"parameters":[{"target":{"kind":"radiation-emissivity","surface":"emitter"},"distribution":{"kind":"uniform","lo":0.2,"hi":1}}]}"#;
    let zero=plan.replace("\"lo\":0.2,\"hi\":1","\"lo\":0.8,\"hi\":0.8");
    let baseline=run(&J::parse(BASE).unwrap());
    let control=success(&dir.uq(&zero,&[]));
    near(n(&control,"mean_k"),n(baseline.get("objective").unwrap(),"value_k"),1e-10);
    assert_eq!(n(&control,"std_dev_k"),0.0);
    let full=dir.uq(plan,&["--checkpoint","full.uqcp"]);let result=success(&full);
    assert!(n(&result,"std_dev_k")>1e-3);
    let chunk=dir.uq(plan,&["--checkpoint","part.uqcp","--max-new-samples","1"]);
    assert_eq!(chunk.status.code(),Some(6));
    let saved=std::fs::read(dir.0.join("part.uqcp")).unwrap();
    let resumed=dir.uq(plan,&["--resume","part.uqcp","--checkpoint","done.uqcp"]);success(&resumed);
    assert_eq!(full.stdout,resumed.stdout);
    assert_eq!(std::fs::read(dir.0.join("full.uqcp")).unwrap(),std::fs::read(dir.0.join("done.uqcp")).unwrap());
    assert_eq!(saved,std::fs::read(dir.0.join("part.uqcp")).unwrap());
    let invalid=plan.replace("\"hi\":1","\"hi\":1.1");
    let refusal=dir.uq(&invalid,&[]);assert!(!refusal.status.success());assert!(refusal.stdout.is_empty());
    let ambient=plan.replace("radiation-emissivity","radiation-ambient-temperature")
        .replace("\"lo\":0.2,\"hi\":1","\"lo\":300,\"hi\":310");
    let refusal=dir.uq(&ambient,&[]);assert!(!refusal.status.success());assert!(refusal.stdout.is_empty());
}
