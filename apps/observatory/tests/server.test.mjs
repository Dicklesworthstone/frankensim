import {test} from 'node:test';
import assert from 'node:assert/strict';
import {mkdtemp, writeFile, rm} from 'node:fs/promises';
import {tmpdir} from 'node:os';
import {join} from 'node:path';
import http from 'node:http';
import {createObservatory} from '../server.mjs';
const inputs = {voltage:0.2,duration:0.1,resistivity:2,support:'tension',strain:0.0001,limit:0.01};
async function start(engine, options={}) {
  const app=createObservatory({engine,...options});
  await new Promise(r=>app.server.listen(0,'127.0.0.1',r));
  const url=`http://127.0.0.1:${app.server.address().port}`;
  const status=await (await fetch(url+'/api/status')).json();
  return {...app,url,status,run:(data=inputs,token=status.token)=>fetch(url+'/api/run', {method:'POST',
    headers:{'Content-Type':'application/json','x-observatory-token':token},body:JSON.stringify(data)})};
}
test('loopback service is available without a native binary and has no file browsing API', async()=> {
  const app=await start('/nonexistent/frankensim');
  try {
    assert.equal(app.status.ready,false); assert.equal((await app.run()).status,503);
    assert.equal((await app.run(inputs,'bad-token')).status,403);
    assert.equal((await fetch(app.url+'/../../etc/passwd')).status,404);
    assert.equal((await fetch(app.url+'/api/status',{headers:{Origin:'https://example.com'}})).status,403);
    const hostStatus = await new Promise(resolve=>{
      http.get(app.url+'/api/status',{headers:{Host:'attacker.example'}},res=>{res.resume();resolve(res.statusCode);});
    });
    assert.equal(hostStatus,403);
    assert.match((await fetch(app.url+'/api/status')).headers.get('content-security-policy'),/frame-ancestors 'none'/);
  } finally {app.shutdown();}
});
const mock = `#!/usr/bin/env node
import {writeFileSync} from 'node:fs';
const a=process.argv.slice(2);const i=a[0]==='--demo'?5:6;
const csv='time_s,cold_pressure_pa,heated_pressure_pa\\n'+Array.from({length:800},(_,j)=>j/16000+',0.01,0.02').join('\\n');
setTimeout(()=>{
 writeFileSync(a[i]+'.csv',csv);
 console.error('support=None; tension_cold_N=10; tension_heated_N=10; heated_layer_stresses_Pa=None');
 console.error('T_layers_K=[301, 302]; EA_cold_N=150000; EA_heated_N=149000; electrical_J=1; stored_heat_J=1; thermal_residual_J=0');
},300);
`;
test('actual process launch/parse, one-job limit, and timeout cleanup (test double, not Rust physics)', {skip:process.platform==='win32'},async()=> {
  const dir=await mkdtemp(join(tmpdir(),'observatory-test-'));const engine=join(dir,'fixture.mjs');
  await writeFile(engine,mock,{mode:0o700});const app=await start(engine);
  try {
    const first=app.run();
    // Wait for the job to be admitted, not a wall-clock guess.
    for(let n=0;n<100;n++){if((await (await fetch(app.url+'/api/status')).json()).busy)break;await new Promise(r=>setTimeout(r,5));}
    assert.equal((await app.run()).status,409);
    const response=await first;assert.equal(response.status,200);
    assert.equal((await response.json()).snapshot.metrics.hot_tension_n,10);
  } finally {app.shutdown();}
  const slow=await start(engine,{deadline:10});
  try {
    const r=await slow.run();assert.equal(r.status,400);assert.match((await r.json()).error,/execution budget/);
    assert.equal((await (await fetch(slow.url+'/api/status')).json()).busy,false);
  } finally {slow.shutdown();await rm(dir,{recursive:true,force:true});}
});
