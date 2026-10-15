// Convert existing native fs-viz ASCII output; does not run or approximate physics.
import {readFile, writeFile, stat} from 'node:fs/promises';
import {parseVTU} from '../import-vtu.mjs';
import {MESH_LIMITS, joinFrames} from '../mesh-data.mjs';
const usage='node apps/observatory/tools/vtu-to-observatory.mjs --units m --output result.json [--times 0,0.1,0.2] frame0.vtu [frame1.vtu ...]';
async function main() {
  const args=process.argv.slice(2), files=[]; let lengthUnit,output,times;
  for(let i=0;i<args.length;i++) {
    if(args[i]==='--help'){console.log(usage);return;}
    if(['--units','--output','--times'].includes(args[i])) {
      const key=args[i], value=args[++i]; if(!value)throw new Error(usage);
      if(key==='--units')lengthUnit=value;
      if(key==='--output')output=value;
      if(key==='--times')times=value.split(',').map(v=>v.trim()===''?NaN:Number(v));
    } else if(args[i].startsWith('--'))throw new Error(usage); else files.push(args[i]);
  }
  if(!output || !files.length || files.length>MESH_LIMITS.frames || (times && times.length!==files.length) || (files.length>1&&!times))throw new Error(usage);
  const snapshots=[]; let bytes=0;
  for(let i=0;i<files.length;i++) {
    bytes+=(await stat(files[i])).size; if(bytes>MESH_LIMITS.bytes)throw new Error('Sequence exceeds the aggregate 32 MiB input budget.');
    snapshots.push(parseVTU(await readFile(files[i],'utf8'),{lengthUnit}));
  }
  const result=times?joinFrames(snapshots,times):snapshots[0];
  if(times) {
    result.source.notes=result.source.notes.filter(n=>!n.startsWith('Static snapshot:'));
    result.source.notes.push('Recorded frame times explicitly supplied through --times. Frames are replayed without interpolation.');
  }
  result.title=files.length>1?'Recorded VTU field sequence':'VTU field snapshot';
  const json=JSON.stringify(result)+'\n'; if(Buffer.byteLength(json)>MESH_LIMITS.bytes)throw new Error('Converted snapshot exceeds 32 MiB.');
  await writeFile(output,json,{flag:'wx'});
  console.log(`Wrote ${result.mesh.time_s.length} frame(s), ${result.mesh.fields.length} field(s), ${result.mesh.triangles.length} exterior triangles. No physics was computed.`);
}
main().catch(e=>{console.error(e.message);process.exitCode=1;});
