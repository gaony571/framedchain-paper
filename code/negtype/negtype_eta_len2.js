// eta_+ for random two-bond chains (S3 Text 2), using the helpers of framedchain-paper/code/negtype/negtype.js
const path = process.argv[2];
const src = require('fs').readFileSync(path + '/negtype.js', 'utf8').split("const mode =")[0];
eval(src);
function angleBetween(Q1,Q2){ const P = matMul(Q1, transpose(Q2)); const c=(P[0]+P[4]+P[8]-1)/2; return Math.acos(Math.max(-1,Math.min(1,c))); }
const fD = th => 0.5*Math.min(Math.sqrt(6), 4*Math.sin(th/4));
for (const [N, reps] of [[60,3],[100,5],[200,3],[400,2]]) {
  const out = [];
  for (let r = 0; r < reps; r++) {
    const Qs=[]; for(let i=0;i<N;i++) Qs.push(randRot());
    const Dm=Qs.map(a=>Qs.map(b=>fD(angleBetween(a,b))));
    const rm=Dm.map(r=>r.reduce((p,q)=>p+q,0)/N), all=rm.reduce((p,q)=>p+q,0)/N;
    const C=Dm.map((r,i)=>r.map((v,j)=>v-rm[i]-rm[j]+all));
    const ev=jacobiEig(C); const pos=ev.filter(v=>v>1e-9).reduce((p,q)=>p+q,0), tot=ev.reduce((p,q)=>p+Math.abs(q),0);
    out.push((pos/tot).toExponential(2) + ' (max eig ' + Math.max(...ev).toExponential(2) + ')');
  }
  console.log('N=' + N + ': eta_+ =', out.join(', '));
}
