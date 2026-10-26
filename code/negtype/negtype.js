// 이 공간의 정렬 거리 D 가 음의 유형인지 작은 구조들로 시험한다.
// 구조 = 변마다 3x3 행렬의 열. D² = Σ‖X‖² + Σ‖Y‖² − 2 S,
// S = max(0, max_{순서를 지키는 부분 대응 π} max_{U∈SO(3)} Σ_{(i,j)∈π} ⟨X_i, U Y_j⟩).
// max_U tr(U M) = σ1 + σ2 + sign(det M) σ3,  M = Σ Y_j X_iᵀ.
// 음의 유형 위반 ⇔ J D J 가 양의 고유값을 가짐.

function matMul(A, B) { const C = [0,0,0,0,0,0,0,0,0]; for (let i=0;i<3;i++) for (let j=0;j<3;j++) { let s=0; for (let k=0;k<3;k++) s += A[3*i+k]*B[3*k+j]; C[3*i+j]=s; } return C; }
function transpose(A) { return [A[0],A[3],A[6],A[1],A[4],A[7],A[2],A[5],A[8]]; }
function det3(A) { return A[0]*(A[4]*A[8]-A[5]*A[7]) - A[1]*(A[3]*A[8]-A[5]*A[6]) + A[2]*(A[3]*A[7]-A[4]*A[6]); }
function frob2(A) { let s=0; for (const v of A) s+=v*v; return s; }

// 대칭 행렬의 고유값 (Jacobi)
function jacobiEig(Ain) {
  const n = Ain.length; const A = Ain.map(r => r.slice());
  for (let sweep=0; sweep<100; sweep++) {
    let off=0; for (let p=0;p<n;p++) for (let q=p+1;q<n;q++) off += A[p][q]*A[p][q];
    if (off < 1e-26) break;
    for (let p=0;p<n;p++) for (let q=p+1;q<n;q++) {
      if (Math.abs(A[p][q]) < 1e-300) continue;
      const th = (A[q][q]-A[p][p])/(2*A[p][q]);
      const t = Math.sign(th||1)/(Math.abs(th)+Math.sqrt(th*th+1));
      const c = 1/Math.sqrt(t*t+1), s = t*c;
      for (let k=0;k<n;k++) { const akp=A[k][p], akq=A[k][q]; A[k][p]=c*akp - s*akq; A[k][q]=s*akp + c*akq; }
      for (let k=0;k<n;k++) { const apk=A[p][k], aqk=A[q][k]; A[p][k]=c*apk - s*aqk; A[q][k]=s*apk + c*aqk; }
    }
  }
  return A.map((r,i)=>r[i]);
}
function sym3Eig(M) { const A=[[M[0],M[1],M[2]],[M[3],M[4],M[5]],[M[6],M[7],M[8]]]; return jacobiEig(A); }

function maxTrSO3(M) {
  const ev = sym3Eig(matMul(transpose(M), M)).map(v => Math.sqrt(Math.max(v,0))).sort((a,b)=>b-a);
  const sgn = det3(M) >= 0 ? 1 : -1;
  return ev[0] + ev[1] + sgn*ev[2];
}

// 순서를 지키는 부분 대응을 모두 늘어놓는다
function alignments(n, m) {
  const out = [];
  function rec(i, j, cur) {
    out.push(cur.slice());
    for (let a=i; a<n; a++) for (let b=j; b<m; b++) { cur.push([a,b]); rec(a+1, b+1, cur); cur.pop(); }
  }
  rec(0, 0, []);
  return out;
}
const alignCache = {};
function D(x, y) {
  const key = x.length + ',' + y.length;
  if (!alignCache[key]) alignCache[key] = alignments(x.length, y.length).filter(p => p.length > 0);
  let S = 0;
  for (const pi of alignCache[key]) {
    const M = [0,0,0,0,0,0,0,0,0];
    for (const [i,j] of pi) { const P = matMul(y[j], transpose(x[i])); for (let k=0;k<9;k++) M[k]+=P[k]; }
    const v = maxTrSO3(M); if (v > S) S = v;
  }
  let nx=0, ny=0; for (const X of x) nx+=frob2(X); for (const Y of y) ny+=frob2(Y);
  return Math.sqrt(Math.max(nx + ny - 2*S, 0));
}

// 무작위 회전 (사원수)
let seed = 12345; function rnd() { seed = (seed*1103515245 + 12345) % 2147483648; return seed/2147483648; }
function gauss() { const u=rnd()||1e-12, v=rnd(); return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v); }
function randRot() {
  let q=[gauss(),gauss(),gauss(),gauss()]; const n=Math.hypot(...q); q=q.map(v=>v/n);
  const [w,x,y,z]=q;
  return [1-2*(y*y+z*z), 2*(x*y-w*z), 2*(x*z+w*y), 2*(x*y+w*z), 1-2*(x*x+z*z), 2*(y*z-w*x), 2*(x*z-w*y), 2*(y*z+w*x), 1-2*(x*x+y*y)];
}
function scale(A, s) { return A.map(v=>v*s); }
function randStruct(maxLen) { const n = 2 + Math.floor(rnd()*(maxLen-1)); const out=[]; for (let i=0;i<n;i++) out.push(scale(randRot(), 0.5)); return out; }

function violation(pts) {
  const N = pts.length; const Dm = [];
  for (let i=0;i<N;i++) { Dm.push([]); for (let j=0;j<N;j++) Dm[i].push(i===j?0:(j<i?Dm[j][i]:D(pts[i],pts[j]))); }
  // J D J
  const rowMean = Dm.map(r => r.reduce((a,b)=>a+b,0)/N); const all = rowMean.reduce((a,b)=>a+b,0)/N;
  const C = Dm.map((r,i)=>r.map((v,j)=> v - rowMean[i] - rowMean[j] + all));
  const ev = jacobiEig(C); const lmax = Math.max(...ev); let dmax=0; for (const r of Dm) for (const v of r) dmax=Math.max(dmax,v);
  return { lmax, rel: lmax/dmax, Dm };
}

const mode = process.argv[2] || 'search';
const maxLen = +(process.argv[3] || 3), N = +(process.argv[4] || 8), trials = +(process.argv[5] || 300);
let best = { rel: -Infinity };
for (let t=0; t<trials; t++) {
  let pts = []; for (let i=0;i<N;i++) pts.push(randStruct(maxLen));
  let v = violation(pts);
  // 언덕 오르기: 점 하나를 바꿔 위반을 키운다
  for (let it=0; it<60; it++) {
    const k = Math.floor(rnd()*N); const old = pts[k]; pts[k] = randStruct(maxLen);
    const v2 = violation(pts); if (v2.rel > v.rel) v = v2; else pts[k] = old;
  }
  if (v.rel > best.rel) best = { rel: v.rel, lmax: v.lmax, pts: pts.map(p=>p.length), Dm: v.Dm };
}
console.log('최대 상대 위반 λmax(JDJ)/max D =', best.rel.toExponential(4), ' λmax =', best.lmax.toExponential(4));
console.log('구조 길이:', best.pts.join(','));
if (best.rel > 1e-8) console.log('D 행렬:\n' + best.Dm.map(r=>r.map(v=>v.toFixed(4)).join(' ')).join('\n'));
