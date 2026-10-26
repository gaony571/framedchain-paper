const src = require('fs').readFileSync('negtype.js','utf8').split("const mode =")[0];
eval(src);
function angleBetween(Q1,Q2){ const P = matMul(Q1, transpose(Q2)); const c=(P[0]+P[4]+P[8]-1)/2; return Math.acos(Math.max(-1,Math.min(1,c))); }
const fD = th => 0.5*Math.min(Math.sqrt(6), 4*Math.sin(th/4));
let maxErr=0;
for (let t=0;t<200;t++){ const Q=randRot(), Qp=randRot(); const I=[1,0,0,0,1,0,0,0,1];
  const x=[scale(I,0.5),scale(Q,0.5)], y=[scale(I,0.5),scale(Qp,0.5)];
  const d=D(x,y), f=fD(angleBetween(Q,Qp)); maxErr=Math.max(maxErr,Math.abs(d-f)); }
console.log('공식과 실제 D 의 최대 차:', maxErr.toExponential(2));
// 무작위 길이 2 구조 N 개에서 J D J 의 최대 고유값
for (const N of [60,120,200]) {
  const Qs=[]; for(let i=0;i<N;i++) Qs.push(randRot());
  const Dm=Qs.map(a=>Qs.map(b=>fD(angleBetween(a,b))));
  const rm=Dm.map(r=>r.reduce((p,q)=>p+q,0)/N), all=rm.reduce((p,q)=>p+q,0)/N;
  const C=Dm.map((r,i)=>r.map((v,j)=>v-rm[i]-rm[j]+all));
  const ev=jacobiEig(C).sort((a,b)=>b-a);
  console.log('N='+N+': J D J 고유값 상위 3 =', ev.slice(0,3).map(v=>v.toExponential(3)).join(' '), ' (음의 유형이면 모두 ≤ 0)');
}
